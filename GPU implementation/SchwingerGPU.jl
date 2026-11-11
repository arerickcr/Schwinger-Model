module SchwingerGPU

using ITensors, ITensorMPS, Adapt
using LinearAlgebra, Random, TOML, Dates

export Model, Solver, Backend, cpu_backend, build_problem, solve, diagnostics,
       initial_labels, gauge_mpo, entropy_profile, middle_tensors, write_result,
       run_case, compare_case, parse_options, options_config, main

"Parameters use exactly the normalization of the original Julia files."
Base.@kwdef struct Model
    N::Int = 20
    x::Float64 = 10.0
    y::Float64 = 1.0
    mass::Float64 = 0.33
    theta::Float64 = pi
    g::Float64 = 0.0
    spin::Float64 = 0.0 # Initial product-state Sz; NOT an enforced sector.
end

Base.@kwdef struct Solver
    nsweeps::Int = 10
    maxdim::Int = 100
    cutoff::Float64 = 0.0
    noise::Float64 = 1e-12
    energy_tol::Float64 = 1e-10
    min_sweeps::Int = 2
    krylovdim::Int = 10
    eig_maxiter::Int = 1
    eig_tol::Float64 = 1e-14
    num_states::Int = 1
    seed::Int = 1234
    outputlevel::Int = 1
end

"General device interface: upload tensors, download tensors, and synchronize."
struct Backend{U,D,S}
    name::String
    upload::U
    download::D
    synchronize::S
    metadata::Dict{String,Any}
end
Backend(name, up, down, sync; metadata=Dict{String,Any}()) =
    Backend(name, up, down, sync, Dict{String,Any}(metadata))
cpu_backend() = Backend("cpu", deepcopy, deepcopy, () -> nothing;
    metadata=Dict("cpu" => Sys.CPU_NAME, "julia_threads" => Threads.nthreads(),
                  "blas_threads" => BLAS.get_num_threads()))

function validate(p::Model, c::Solver)
    p.N >= 2 || error("N must be at least 2")
    all(isfinite, (p.x,p.y,p.mass,p.theta,p.g,p.spin)) || error("Nonfinite model parameter")
    p.x > 0 || error("x must be positive")
    c.nsweeps > 0 && c.maxdim > 0 && c.num_states > 0 || error("Positive sizes required")
    c.min_sweeps >= 2 || error("min_sweeps must be at least 2")
    c.krylovdim >= 2 && c.eig_maxiter > 0 || error("Invalid eigensolver settings")
    all(t -> isfinite(t) && t >= 0, (c.cutoff,c.noise,c.energy_tol,c.eig_tol)) ||
        error("Tolerances and noise must be finite and nonnegative")
    initial_labels(p.N, p.spin)
    return nothing
end

"Alternating state with requested initial Sz; supports even and odd N."
function initial_labels(N::Int, spin::Real)
    nup = N / 2 + spin
    isinteger(nup) && 0 <= nup <= N ||
        error("N=$N requires integer Sz for even N, half-integer Sz for odd N, |Sz| <= N/2")
    labels = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    delta = Int(nup) - count(==("Up"), labels)
    candidates = delta > 0 ? reverse(findall(==("Dn"), labels)) : findall(==("Up"), labels)
    for j in Iterators.take(candidates, abs(delta))
        labels[j] = delta > 0 ? "Up" : "Dn"
    end
    return labels
end

"Bond-dimension-three representation of sum_{n=1}^{N-1} L_n^2."
function gauge_mpo(sites, theta)
    N = length(sites)
    links = [Index(3, "Link,gauge,l=$j") for j in 0:N]
    tensors = ITensor[]
    identity2 = Matrix{Float64}(I, 2, 2)
    sz = [0.5 0.0; 0.0 -0.5]
    for n in 1:N
        left, right = links[n], links[n+1]
        out, input = prime(sites[n]), dag(sites[n])
        W = ITensor(Float64, out, input, left, right)
        # Explicit index identities avoid depending on filterinds ordering.
        function block!(a, b, matrix)
            for row in 1:2, col in 1:2
                W[left=>a, right=>b, out=>row, input=>col] = matrix[row,col]
            end
        end
        for a in 1:3
            block!(a, a, identity2)
        end
        if n < N
            a = sz + ((-1.0)^n / 2 + (n == 1 ? theta/(2pi) : 0.0))*identity2
            block!(1, 2, a)
            block!(1, 3, (N-n)*(a*a))
            block!(2, 3, 2*(N-n)*a)
        end
        push!(tensors, W)
    end
    left = ITensor(links[1]); left[links[1]=>1] = 1.0
    right = ITensor(links[end]); right[links[end]=>3] = 1.0
    tensors[1] *= left
    tensors[end] *= right
    return MPO(tensors)
end

function build_problem(p::Model)
    initial_labels(p.N, p.spin)
    sites = siteinds("S=1/2", p.N; conserve_qns=false)
    hopping, mass = OpSum(), OpSum()
    for j in 1:p.N-1
        hopping += 1.0, "S+", j, "S-", j+1
        hopping += 1.0, "S-", j, "S+", j+1
        if p.g != 0
            hopping += -4*p.g, "Sz", j, "Sz", j+1
        end
    end
    for j in 1:p.N
        mass += (-1.0)^j, "Sz", j
    end
    Hc, Hm, Hxy = gauge_mpo(sites, p.theta), MPO(mass, sites), MPO(hopping, sites)
    mu = 2sqrt(p.x)*p.mass - 0.25
    # Exact direct sum, with no MPO compression or truncation.
    H = +(p.y*Hc, mu*Hm, p.x*Hxy; alg="directsum")
    psi0 = MPS(Float64, sites, initial_labels(p.N, p.spin))
    return (; sites, H, psi0, Hc, Hm, Hxy)
end

mutable struct EnergyObserver <: AbstractObserver
    tolerance::Float64
    min_sweeps::Int
    energies::Vector{Float64}
    truncation_errors::Vector{Float64}
    converged::Bool
end
EnergyObserver(c::Solver) = EnergyObserver(c.energy_tol, c.min_sweeps, Float64[], Float64[], false)

function ITensorMPS.measure!(o::EnergyObserver; sweep, spec, kwargs...)
    while length(o.truncation_errors) < sweep
        push!(o.truncation_errors, 0.0)
    end
    o.truncation_errors[sweep] = max(o.truncation_errors[sweep], spec.truncerr)
    return nothing
end

function ITensorMPS.checkdone!(o::EnergyObserver; energy, sweep, kwargs...)
    E = Float64(real(energy))
    previous = isempty(o.energies) ? Inf : o.energies[end]
    push!(o.energies, E)
    o.converged = o.tolerance > 0 && sweep >= o.min_sweeps &&
        abs(E-previous) <= o.tolerance*max(1.0, abs(E))
    return o.converged
end

function solve(H, psi0, p::Model, c::Solver, backend::Backend)
    states, physical, objectives, observers = MPS[], Float64[], Float64[], EnergyObserver[]
    Random.seed!(c.seed)
    backend.synchronize()
    start = time_ns()
    for k in 1:c.num_states
        obs = EnergyObserver(c) # Independent history for EACH state.
        kwargs = (; nsweeps=c.nsweeps, maxdim=c.maxdim, cutoff=c.cutoff,
            noise=c.noise, observer=obs, eigsolve_krylovdim=c.krylovdim,
            eigsolve_maxiter=c.eig_maxiter, eigsolve_tol=c.eig_tol,
            ishermitian=true, outputlevel=c.outputlevel)
        objective, psi = if k == 1
            dmrg(H, deepcopy(psi0); kwargs...)
        else
            dmrg(H, states, deepcopy(psi0); weight=100sqrt(p.x), kwargs...)
        end
        push!(states, psi); push!(objectives, real(objective)); push!(observers, obs)
    end
    backend.synchronize()
    dmrg_seconds = (time_ns()-start)/1e9
    # Physical energy excludes the excited-state penalty and is normalized.
    for psi in states
        push!(physical, real(inner(psi', H, psi)/inner(psi, psi)))
    end
    return (; states, energies=physical, objectives, observers, dmrg_seconds)
end

"Diagnostics on host tensors; deliberately excluded from DMRG timing."
function diagnostics(H, result, backend::Backend)
    Hcpu = backend.download(H)
    states = [backend.download(psi) for psi in result.states]
    rows = Dict{String,Any}[]
    for (k, psi) in enumerate(states)
        norm2 = real(inner(psi, psi))
        E = real(inner(psi', Hcpu, psi))/norm2
        H2 = real(inner(Hcpu, psi, Hcpu, psi))/norm2
        raw_variance = H2-E^2
        obs = result.observers[k]
        push!(rows, Dict("state"=>k-1, "energy"=>E, "dmrg_objective"=>result.objectives[k],
            "norm2"=>norm2, "variance_raw"=>raw_variance,
            "residual_norm"=>sqrt(max(0.0,raw_variance)),
            "total_Sz"=>sum(expect(psi,"Sz")),
            "maxlinkdim"=>maxlinkdim(psi), "sweeps"=>length(obs.energies),
            "energy_converged"=>obs.converged, "sweep_energies"=>obs.energies,
            "sweep_max_truncation_error"=>obs.truncation_errors))
    end
    overlaps = [abs(inner(states[i],states[j]))/sqrt(real(inner(states[i],states[i])*inner(states[j],states[j])))
                for i in eachindex(states), j in eachindex(states)]
    return (; rows, overlaps, states)
end

function entropy_profile(psi_input::MPS)
    psi = deepcopy(psi_input) # Host MPS; download before calling on CUDA states.
    values = Float64[]
    for b in 1:length(psi)-1
        orthogonalize!(psi,b)
        left = b == 1 ? (siteind(psi,b),) : (linkind(psi,b-1),siteind(psi,b))
        _, S, _ = svd(psi[b], left)
        probs = [abs2(S[n,n]) for n in 1:dim(S,1)]
        probs ./= sum(probs)
        push!(values, -sum(q*log(q) for q in probs if q > 0))
    end
    L = length(psi)
    return (; entropy=values, chord=[2L/pi*sin(pi*b/L) for b in 1:L-1])
end

function middle_tensors(psi_input::MPS)
    psi = deepcopy(psi_input)
    middle = div(length(psi),2)
    orthogonalize!(psi,middle)
    single = deepcopy(psi[middle])
    pair = psi[middle]*psi[middle+1]
    fusion = combiner(siteind(psi,middle),siteind(psi,middle+1); tags="Combined_Site,n=$middle")
    return (; middle, single, pair, fused=pair*fusion, fusion)
end

fields(x) = Dict(string(n)=>getfield(x,n) for n in fieldnames(typeof(x)))
function write_result(path, p, c, backend, result, diagnostic; build_seconds, transfer_seconds)
    mkpath(dirname(abspath(path)))
    data = Dict("created_at"=>string(Dates.now()), "julia_version"=>string(VERSION),
        "ITensors_version"=>string(pkgversion(ITensors)), "ITensorMPS_version"=>string(pkgversion(ITensorMPS)),
        "backend"=>backend.name, "hardware"=>backend.metadata, "model"=>fields(p), "solver"=>fields(c),
        "build_seconds"=>build_seconds, "transfer_seconds"=>transfer_seconds,
        "dmrg_seconds"=>result.dmrg_seconds, "states"=>diagnostic.rows,
        "absolute_overlaps"=>[collect(row) for row in eachrow(diagnostic.overlaps)])
    open(path,"w") do io
        TOML.print(io,data; sorted=true)
    end
end

function run_case(p::Model, c::Solver, backend::Backend; output=nothing, problem=nothing)
    validate(p,c)
    start = time_ns()
    prob = isnothing(problem) ? build_problem(p) : problem
    build_seconds = (time_ns()-start)/1e9
    backend.synchronize(); start = time_ns()
    H, psi0 = backend.upload(prob.H), backend.upload(prob.psi0)
    backend.synchronize(); transfer_seconds = (time_ns()-start)/1e9
    result = solve(H,psi0,p,c,backend)
    diagnostic = diagnostics(H,result,backend)
    if !isnothing(output)
        write_result(output,p,c,backend,result,diagnostic; build_seconds,transfer_seconds)
    end
    return (; result, diagnostic, build_seconds, transfer_seconds)
end

function compare_case(p, c, backend; directory="results", label="comparison", atol=1e-8, rtol=1e-9)
    # Shared indices and identical host initial state permit meaningful overlaps.
    problem = build_problem(p)
    cpu = run_case(p,c,cpu_backend(); problem,output=joinpath(directory,"$(label)_cpu.toml"))
    gpu = run_case(p,c,backend; problem,output=joinpath(directory,"$(label)_$(backend.name).toml"))
    rows = Dict{String,Any}[]
    for k in 1:c.num_states
        a, b = cpu.diagnostic.rows[k], gpu.diagnostic.rows[k]
        delta = abs(a["energy"]-b["energy"])
        u,v = cpu.diagnostic.states[k],gpu.diagnostic.states[k]
        fidelity = abs2(inner(u,v))/real(inner(u,u)*inner(v,v))
        push!(rows,Dict("state"=>k-1,"energy_abs_difference"=>delta,
            "energy_agrees"=>delta<=atol+rtol*abs(a["energy"]),"fidelity"=>fidelity,
            "both_energy_converged"=>a["energy_converged"] && b["energy_converged"]))
    end
    report = Dict("model"=>fields(p),"solver"=>fields(c),"states"=>rows,
        "atol"=>atol,"rtol"=>rtol,"cpu_dmrg_seconds"=>cpu.result.dmrg_seconds,
        "gpu_dmrg_seconds"=>gpu.result.dmrg_seconds,
        "dmrg_speedup"=>cpu.result.dmrg_seconds/gpu.result.dmrg_seconds,
        "timing_note"=>"DMRG only; excludes setup, transfer, diagnostics. Warm up before interpretation.")
    open(joinpath(directory,"$(label)_comparison.toml"),"w") do io
        TOML.print(io,report; sorted=true)
    end
    return report
end

function parse_options(args)
    opts = Dict{String,String}()
    for arg in args
        startswith(arg,"--") && occursin('=',arg) || error("Use --key=value; got $arg")
        key,value = split(arg[3:end],'=';limit=2)
        haskey(opts,key) && error("Duplicate option $key")
        opts[key]=value
    end
    allowed = Set(["N","x","y","mass","theta","g","spin","nsweeps","maxdim","cutoff",
        "noise","energy_tol","min_sweeps","krylovdim","eig_maxiter","eig_tol","num_states",
        "seed","outputlevel","output","blas_threads","device","sizes","repeats","warmup","directory"])
    isempty(setdiff(Set(keys(opts)),allowed)) || error("Unknown option(s): $(setdiff(Set(keys(opts)),allowed))")
    return opts
end

function options_config(opts; N=parse(Int,get(opts,"N","20")))
    modelargs = Dict{Symbol,Any}(:N=>N)
    for name in fieldnames(Model)
        name == :N && continue
        default = name == :spin ? (iseven(N) ? 0.0 : 0.5) : getfield(Model(),name)
        modelargs[name] = parse(Float64,get(opts,string(name),string(default)))
    end
    defaults = Solver()
    solverargs = Dict(name=>parse(fieldtype(Solver,name),get(opts,string(name),string(getfield(defaults,name))))
                      for name in fieldnames(Solver))
    return Model(;modelargs...),Solver(;solverargs...)
end

function main(args, backend_factory)
    opts = parse_options(args)
    BLAS.set_num_threads(parse(Int,get(opts,"blas_threads","4")))
    p,c = options_config(opts)
    backend = backend_factory(opts)
    path = get(opts,"output",joinpath(@__DIR__,"results","$(backend.name)_N$(p.N).toml"))
    run = run_case(p,c,backend;output=path)
    for row in run.diagnostic.rows
        println("state=$(row["state"]) E=$(row["energy"]) residual=$(row["residual_norm"]) energy_converged=$(row["energy_converged"])")
    end
    println("Saved $path")
    return run
end
end
