using Test, ITensors, ITensorMPS, LinearAlgebra, TOML
include("../SchwingerGPU.jl")
using .SchwingerGPU

# Load the original helper file without editing it or running its N=20 driver.
module Legacy
include("../reference/OBC_Schwinger_dmrg.jl")
end

"Independent spin-basis Hamiltonian, only for tiny validation sizes."
function exact_matrix(p)
    D = 2^p.N
    H = zeros(D,D)
    mu = 2sqrt(p.x)*p.mass-0.25
    for bits in 0:D-1
        sz = [iszero(bits & (1 << (n-1))) ? 0.5 : -0.5 for n in 1:p.N]
        field = p.theta/(2pi)
        energy = mu*sum((-1.0)^n*sz[n] for n in 1:p.N)
        for n in 1:p.N-1
            field += sz[n]+(-1.0)^n/2
            energy += p.y*field^2 - 4p.x*p.g*sz[n]*sz[n+1]
            if sz[n] != sz[n+1]
                flipped = xor(bits,(1 << (n-1)) | (1 << n))
                H[flipped+1,bits+1] += p.x
            end
        end
        H[bits+1,bits+1] += energy
    end
    return H
end

function dense_mpo(H,sites)
    tensor = prod(H)
    return reshape(Array(tensor,prime.(sites)...,dag.(sites)...),2^length(sites),2^length(sites))
end

BLAS.set_num_threads(4)
@testset "Initial states and validation" begin
    @test initial_labels(4,0) == ["Up","Dn","Up","Dn"]
    # Original positive-spin convention flips even sites from the right.
    @test initial_labels(6,2) == ["Up","Dn","Up","Up","Up","Up"]
    @test count(==("Up"),initial_labels(5,0.5)) == 3
    @test count(==("Up"),initial_labels(5,-0.5)) == 2
    @test_throws ErrorException initial_labels(5,0)
    @test_throws ErrorException initial_labels(4,3)
    p,c = options_config(Dict{String,String}();N=11)
    @test p.spin == 0.5
    @test_throws ErrorException parse_options(["--typo=1"])
end

@testset "Hamiltonian: independent matrix and unchanged original" begin
    for N in 2:6, theta in (0.0,pi,0.73)
        p = Model(N=N,theta=theta,g=0.17,y=0.8,mass=0.21,spin=iseven(N) ? 0.0 : 0.5)
        prob = build_problem(p)
        actual = dense_mpo(prob.H,prob.sites)
        @test actual ≈ exact_matrix(p) atol=2e-12 rtol=2e-13
        @test actual ≈ actual' atol=2e-12
        Hc = Legacy.get_MPO_neutral((W,n)->Legacy.schwinger_gauge(W,n,N,theta),prob.sites,3,false)
        @test dense_mpo(Hc,prob.sites) ≈ dense_mpo(prob.Hc,prob.sites) atol=2e-12
        mu = 2sqrt(p.x)*p.mass-0.25
        original = Legacy.directsum_MPOs(Legacy.directsum_MPOs(p.y*Hc,mu*prob.Hm),p.x*prob.Hxy)
        @test dense_mpo(original,prob.sites) ≈ actual atol=2e-12
    end
end

@testset "CPU DMRG: exact energies, residuals, observables, original solver" begin
    for N in (4,5,6)
        p = Model(N=N,spin=iseven(N) ? 0.0 : 0.5)
        c = Solver(nsweeps=16,maxdim=32,noise=0.0,eig_maxiter=8,outputlevel=0)
        run = run_case(p,c,cpu_backend();output=joinpath(@__DIR__,"..","results","validation_N$N.toml"))
        row = run.diagnostic.rows[1]
        Eexact = eigmin(Hermitian(exact_matrix(p)))
        println("N=$N DMRG=$(row["energy"]) exact=$Eexact error=$(abs(row["energy"]-Eexact))")
        @test row["energy"] ≈ Eexact atol=1e-9
        @test row["norm2"] ≈ 1 atol=1e-12
        @test row["residual_norm"] < 5e-6
        @test row["energy_converged"]
        psi = run.diagnostic.states[1]
        entropy = entropy_profile(psi)
        @test length(entropy.entropy) == N-1
        @test all(isfinite,entropy.entropy)
        # Independent Schmidt decomposition of the complete state vector.
        v = vec(Array(prod(psi),siteinds(psi)...))
        for b in 1:N-1
            q = svdvals(reshape(v,2^b,2^(N-b))).^2
            expected = -sum(t*log(t) for t in q if t>0)
            @test entropy.entropy[b] ≈ expected atol=1e-11
        end
        mid = middle_tensors(psi)
        @test norm(mid.fused*dag(mid.fusion)-mid.pair) < 1e-12
        @test dim(combinedind(mid.fusion)) == 4
        prob = build_problem(p)
        @test all(abs.(entropy_profile(prob.psi0).entropy) .< 1e-14)
        if N == 4
            original_Hc = Legacy.get_MPO_neutral((W,n)->Legacy.schwinger_gauge(W,n,N,p.theta),prob.sites,3,false)
            oldE,oldpsi = Legacy.Schwinger(N,p.x,p.y,p.mass,p.theta,16,32,p.spin,
                prob.psi0,original_Hc,prob.Hm,prob.Hxy,1)
            @test oldE[1] ≈ row["energy"] atol=1e-8
        end
    end
end

@testset "Excited states and near-zero convergence" begin
    p = Model(N=4)
    c = Solver(nsweeps=20,maxdim=16,noise=0.0,eig_maxiter=10,num_states=2,outputlevel=0)
    run = run_case(p,c,cpu_backend();output=joinpath(@__DIR__,"..","results","validation_excited_N4.toml"))
    exact = eigvals(Hermitian(exact_matrix(p)))
    @test run.result.energies ≈ exact[1:2] atol=1e-8
    @test run.diagnostic.overlaps[1,2] < 1e-7
    @test run.result.observers[1] !== run.result.observers[2]
    obs = SchwingerGPU.EnergyObserver(c)
    @test !ITensorMPS.checkdone!(obs;energy=0.0,sweep=1)
    @test ITensorMPS.checkdone!(obs;energy=0.0,sweep=2)
end

@testset "Generic backend and report path" begin
    p,c = Model(N=4),Solver(nsweeps=10,maxdim=16,noise=0.0,eig_maxiter=5,outputlevel=0)
    # Exercises the adapter interface and comparison/report code without claiming a GPU run.
    other = Backend("host_adapter",deepcopy,deepcopy,()->nothing)
    report = compare_case(p,c,other;directory=joinpath(@__DIR__,"..","results"),label="host_adapter_N4")
    @test report["states"][1]["energy_agrees"]
    @test report["states"][1]["fidelity"] > 1-1e-12
end
