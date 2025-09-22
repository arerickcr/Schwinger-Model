# This code gives the spectrum of Schwinger model with OBC with preservation (by hand) of the translation by one site symmetry (at theta = pi):
# H = 1/2 * [H(starting at site = 0) + H(starting at site = 1)]

using ITensors, ITensorMPS
using ITensors.Strided
using LinearAlgebra
using Printf
using Random
using HDF5
BLAS.set_num_threads(4)

II = 0 + 1im;
X = [0 1;1 0];
Y = [0 -II; II 0];
sp = [0 1;0 0];
sm = [0 0;1 0];
Z = [1 0;0 -1];
Id = [1 0;0 1];

function assign_itensor_elements(T::ITensor, link_indices::Tuple{Int,Int}, M::AbstractMatrix)
    a, b = link_indices
    nonzero_dict = extract_nonzero_elements(M)
    
    # Get the relevant indices from the ITensor
    l1, l2 = filterinds(T, "Link")
    s1, s2 = filterinds(T, "S=1/2")
    
    # Assign values to the ITensor
    for ((i, j), value) in nonzero_dict
        # println(((a,b,i, j), value))
        T[l1 => a, l2 => b, s1 => i, s2 => j] = value
        # println("done")
    end
end

# Manual building of the MPO for the Coulomb part for site starting at 1
function schwinger_gauge_site1(ham_bulk, n, M, theta)
    N = M/2
    ham_bulk_1 = ham_bulk
    assign_itensor_elements(ham_bulk_1,(1,1),Id)
    assign_itensor_elements(ham_bulk_1,(2,2),Id)
    assign_itensor_elements(ham_bulk_1,(3,3),Id)
    if n == N+1
        assign_itensor_elements(ham_bulk_1,(1,2),1/2*Z + 1/2*(-1)^(n) * Id + theta/(2*pi) * Id)
        assign_itensor_elements(ham_bulk_1,(1,3), (N-1) * (1/2*Z + 1/2*(-1)^(n) * Id + theta/(2*pi) * Id)^2 )
        assign_itensor_elements(ham_bulk_1,(2,3), 2*(N-1) * (1/2*Z + 1/2*(-1)^(n) * Id + theta/(2*pi)* Id)  )
    elseif n > N+1 && n < 2*N
        assign_itensor_elements(ham_bulk_1,(1,2),1/2*Z + 1/2*(-1)^(n) * Id)
        assign_itensor_elements(ham_bulk_1,(1,3), (N-n) * (1/2*Z + 1/2*(-1)^(n) * Id)^2 )
        assign_itensor_elements(ham_bulk_1,(2,3), 2*(N-n) * (1/2*Z + 1/2*(-1)^(n) * Id) )
    end

    return ham_bulk_1
end

# Manual building of the MPO for the Coulomb part for site starting at 0
function schwinger_gauge_site0(ham_bulk, n, M, theta)
    N = M/2
    ham_bulk_1 = ham_bulk
    assign_itensor_elements(ham_bulk_1,(1,1),Id)
    assign_itensor_elements(ham_bulk_1,(2,2),Id)
    assign_itensor_elements(ham_bulk_1,(3,3),Id)
    if n == 1
        assign_itensor_elements(ham_bulk_1,(1,2),1/2*Z + 1/2*(-1)^(n-1) * Id + theta/(2*pi) * Id)
        assign_itensor_elements(ham_bulk_1,(1,3), (N-1) * (1/2*Z + 1/2*(-1)^(n-1) * Id + theta/(2*pi) * Id)^2 )
        assign_itensor_elements(ham_bulk_1,(2,3), 2*(N-1) * (1/2*Z + 1/2*(-1)^(n-1) * Id + theta/(2*pi)* Id)  )
    elseif n < N
        assign_itensor_elements(ham_bulk_1,(1,2),1/2*Z + 1/2*(-1)^(n-1) * Id)
        assign_itensor_elements(ham_bulk_1,(1,3), (N-n) * (1/2*Z + 1/2*(-1)^(n-1) * Id)^2 )
        assign_itensor_elements(ham_bulk_1,(2,3), 2*(N-n) * (1/2*Z + 1/2*(-1)^(n-1) * Id) )
    end

    return ham_bulk_1
end

# Manual building of the MPO for the mass and kinetic terms
# function schwinger_fermion(ham_bulk, n, x, mu)
#     ham_bulk_1 = ham_bulk
#     assign_itensor_elements(ham_bulk_1,(1,1),Id)
#     assign_itensor_elements(ham_bulk_1,(1,2),sp)
#     assign_itensor_elements(ham_bulk_1,(1,3),sm)
#     assign_itensor_elements(ham_bulk_1,(1,4), mu/2 * (-1)^(n) * Z)
#     assign_itensor_elements(ham_bulk_1,(2,4), x * sm)
#     assign_itensor_elements(ham_bulk_1,(3,4), x * sp)
#     assign_itensor_elements(ham_bulk_1,(4,4),Id)

#     return ham_bulk_1
# end

function extract_nonzero_elements(matrix)
    d = size(matrix, 1)
    if size(matrix, 2) != d
        error("Input must be a square matrix")
    end
    
    T = eltype(matrix)
    nonzero_dict = Dict{Tuple{Int, Int}, T}()
    nonzero_list = T[]
    
    for j in 1:d, i in 1:d
        if !iszero(matrix[i, j])
            push!(nonzero_list, matrix[i, j])
            nonzero_dict[(i, j)] = matrix[i, j]
        end
    end
    
    return nonzero_dict
end

function get_MPO_neutral(fn, sites, D_ham)
    N = length(sites)
    links_extend = [Index([QN() => 3],string("Link,l=",l)) for l in 0:N]
    # links = links_extend[2:N]
    tensors = [ITensor() for i in 1:N]
    for i in 1:N
        # Define your operator here. This is just an example.
        tt = ITensor(Float64, sites[i]', dag(sites[i]), dag(links_extend[i]), links_extend[i+1])
        tensors[i] = fn(tt, i)
    end
    v_i = ITensor(links_extend[1])
    v_i[links_extend[1] => 1] = 1
    v_f = ITensor(dag(links_extend[N+1]))
    v_f[links_extend[N+1] => D_ham] = 1
    tensors[1] = tensors[1] * v_i
    tensors[N] = tensors[N] * v_f
    return MPO(tensors)
end

function create_custom_sites(N::Int)
    sites = Vector{Index}(undef, N)
    
    for i in 1:N
        if i % 2 == 0  # Odd sites
            sites[i] = Index([
                QN(("q", -1)) => 1,  # Sz = -1/2, q = -1
                QN(("q",  0)) => 1   # Sz = +1/2, q =  0
            ], "Site,S=1/2,n=$i")
        else  # Even sites
            sites[i] = Index([
                QN(("q",  0)) => 1,  # Sz = -1/2, q = 0
                QN(("q", +1)) => 1   # Sz = +1/2, q = 1
            ], "Site,S=1/2,n=$i")
        end
    end
    
    return sites
end

# Function that implements the direct sum of 2 MPOs
function directsum_MPOs(MPO1, MPO2)
  N = length(MPO1)
  all_tags = tags.(linkinds(MPO1))
  # @show all_tags
  tensors = Vector{ITensor}(undef,N)
  ind_in = undef
  for (T1, T2, i) in zip(MPO1, MPO2, 1:N)
      tags = all_tags[(i==1 ? 1 : i-1):(i==N ? N-1 : i)]
      inds1 = [filterinds(T1, tags=tt)[1] for tt in tags]
      inds2 = [filterinds(T2, tags=tt)[1] for tt in tags]
      T_out, inds_out = directsum(T1 => Tuple(inds1), T2 => Tuple(inds2), tags = tags)
      if i > 1
          replaceind!(T_out, inds_out[1], dag(ind_in))
      end
      ind_in = inds_out[end]
      tensors[i] = T_out
  end
  return MPO(tensors)
end

function chrd_dist(j,L)
  return 2*L/π * sin(π/L * j)
end

#Calculate the entanglement entropy vs subinterval length
function EEcalc(psi_input)
  psi = copy(psi_input)#make a copy to not alter the original
  L = length(psi)
  eevals =  zeros(L-1)#Create a zero vector to store the entanglement entropy values. Initialized to zero
  for r=2:L-1 #from first point till last
      orthogonalize!(psi,r)#shifting orthogonality centre to position r
      #performing singular value decomposition
      U,S,V = svd(psi[r], (linkind(psi, r-1), siteind(psi,r)))
      for n=1:dim(S, 1)
          p = S[n,n]^2 #square of the diagonal elements give p_n
          eevals[r] -= p * log(p)
      end
  end
  return eevals, map(r -> chrd_dist(r, L), 1:L-1)
end

mutable struct DemoObserver <: AbstractObserver
   energy_tol::Float64
   last_energy::Float64

   DemoObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

#Function that defines the Observer (checks the convergence of each sweep and stops once the tolerance is reached)
function ITensorMPS.checkdone!(o::DemoObserver;kwargs...)
  sw = kwargs[:sweep]
  energy = kwargs[:energy]
  if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
    println("Stopping DMRG after sweep $sw, with energy convergence < ", o.energy_tol)
    return true
  end
  # Otherwise, update last_energy and keep going
  o.last_energy = energy
  return false
end

function ITensorMPS.measure!(o::DemoObserver; kwargs...)
  energy = kwargs[:energy]
  sweep = kwargs[:sweep]
  bond = kwargs[:bond]
  outputlevel = kwargs[:outputlevel]

  #if outputlevel > 0
  #  println("Sweep $sweep at bond $bond, the energy is $energy")
  #end
end

function Schwinger(N,x,y,me,theta,nsweeps,maxdim,spin,psi0_init,Hc1,Hm1,Hxy1,Hc0,Hm0,Hxy0,num_states)
  mu = 2*sqrt(abs(x))*me - 1/4
  #mu = 2*sqrt(x)*me

  weight = 100*sqrt(abs(x))       # Weight for excited states, defined as H' = H + w |psi0><psi0|
                            # w has to be at least bigger than the energy gap

  #H = y*Hc + mu*Hm + x*Hxy

  H1 = directsum_MPOs(directsum_MPOs(y*Hc1,mu*Hm1), x*Hxy1)
  H0 = directsum_MPOs(directsum_MPOs(y*Hc0,mu*Hm0), x*Hxy0)
  
  H = directsum_MPOs(1/2*H0,1/2*H1)

  # Set maximum error allowed when adapting bond dimensions
  cutoff = [1E-15]
  #cutoff = [0]
  # If DMRG is far from the global minumum then there is no guarantee that DMRG will be able to find the true ground state. 
  # This problem is exacerbated for quantum number conserving DMRG where the search space is more constrained.
  # If this happens, a way out is to turn on the noise term feature to be a very small number.
  noise = [1E-15]
   
  eigsolve_krylovdim = 10 # Maximum dimension of Krylov space to locally solve problem. Try setting to a higher
                          #     value if convergence is slow or the Hamiltonian is close to a critical point.
  eigsolve_maxiter = 1    # Number of times Krylov space can be rebuild
  
  obs = DemoObserver(1E-10)

  psi_states = Array{Any}(undef,num_states)
  energies = Array{Any}(undef,num_states)

  for i in 1:num_states
    if i == 1
    energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise, observer=obs, eigsolve_krylovdim, eigsolve_maxiter)
    psi_states[i] = psi0
    energies[i] = energy0
    else
    energy1, psi1 = dmrg(H,[psi_states[j] for j in 1:(i-1)],psi0_init; nsweeps, maxdim, cutoff, noise, observer=obs, weight, eigsolve_krylovdim, eigsolve_maxiter) 
    psi_states[i] = psi1
    energies[i] = energy1
    end
  end

  return (energies,psi_states)

end

let
ITensors.enable_threaded_blocksparse()

N = 20    # Number of lattice sites
x, y = 10, 1    # Specify parameters

# Plan to do nsweeps DMRG sweeps:
nsweeps = 10

# Set maximum MPS bond dimensions for each sweep (truncation m)
maxdim = 500
theta = pi
spin = 0    #Spin projection
num_states = 10 

sites = siteinds("S=1/2", 2*N; conserve_qns=true) 

  # Initialize the state with the quantum number (spin) desired to initialize the Lanczos algorithm
  if spin==0
    state = [isodd(n) ? "Up" : "Dn" for n in 1:2*N]
  elseif spin==1
    state = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" elseif n>N+1 && n<2*N+1 if isodd(n) "Up" else "Dn" 
                end elseif (n==N+1) "Dn" end for n in 1:2*N]
  elseif spin==-1
    state = [if (n>1 && n<N+1) if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" elseif (n>N && n<2*N) if isodd(n) "Up" else "Dn" 
                end elseif (n==2*N) "Up" end for n in 1:2*N]
  end

  # Create an MPS for previous states
  psi_init = MPS(sites, state)

psi0_init = psi_init;

me = 0.3

os0 = OpSum();   # Initialize a sum of operators for the Hamiltonian
os1 = OpSum();

 for j in 1:(N - 1)
    os0 += 1, "S+", j, "S-", j + 1
    os0 += 1, "S-", j, "S+", j + 1
    os1 += 1, "S+", N + j, "S-", N + j + 1
    os1 += 1, "S-", N + j, "S+", N + j + 1
 end

 # Create the MPO for the Hamiltonian
 Hxy1 = MPO(os1, sites);
 Hxy0 = MPO(os0, sites);

 # Mass term in the Hamiltonian
 os1 = OpSum();
 os0 = OpSum();

 for j in 1:N
     os0 .+= (-1)^(j-1), "Sz", j   #Sites start at 0
     #os .+= mu/2, "Id", j
 end
 for j in N+1:2*N
    os1 .+= (-1)^(j), "Sz", j     #Sites start at 1
 end

 Hm1 = MPO(os1, sites);
 Hm0 = MPO(os0, sites);
 
 fn_gauge_site1 = (ham1, n1) -> schwinger_gauge_site1(ham1, n1, 2*N, 2*pi*spin + theta);
 Hc1 = get_MPO_neutral(fn_gauge_site1, sites, 3);

 fn_gauge_site0 = (ham0, n0) -> schwinger_gauge_site0(ham0, n0, 2*N, theta);
 Hc0 = get_MPO_neutral(fn_gauge_site0, sites, 3);

 energies,psi_states = Schwinger(N,x,y,me,theta,nsweeps,maxdim,spin,psi0_init,Hc1,Hm1,Hxy1,Hc0,Hm0,Hxy0,num_states);
 psi0_initial = psi_states[1]

#Check normalization of states
  println("\nCheck normalization of states")
  for i in 1:num_states
    println("<psi",string(i-1),"|psi",string(i-1),"> = ", inner(psi_states[i]',psi_states[i]))
  end

  #Check orthogonalization of states
  println("\nCheck orthogonalization of states")
  for i in 1:num_states
    for j in (i+1):num_states
      println("<psi",string(i-1),"|psi",string(j-1),"> = ", inner(psi_states[i]',psi_states[j]))
    end
  end

  println("N = ", string(N), "\nx = ", string(x))
  println("m/e = ", string(me))

  for i in 1:num_states
    println("E",string(i-1)," = ", string(energies[i]))
  end

 end
ITensors.disable_threaded_blocksparse()
end

##############################################################################################
#=
Ground energy of Hamiltonian starting at site = 1
N=50 
m/e = 0.3
E0 = -300.9833783663694
spin = 0
theta = π

Ground energy of Hamiltonian starting at site = 0
N=50 
m/e = 0.3
E0 = -298.7293294637072
spin = 0
theta = π 

av= 1/2*(-300.9833783663694-298.7293294637072) = -299.8563539150383

Ground energy of average Hamiltonian
N=50 
m/e = 0.3
E0 = -299.8563539150383
spin = 0
theta = π 
=#
