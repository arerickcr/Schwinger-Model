# This code gives the spectrum of the Schwinger model with Periodic Boundary Conditions at fixed N and x.
# Note that we include the code "S5(modified).jl" which is a modification of the operators of the spin = 5 sites code. This is needed for the 
# to truncate the bosonic dof L. Here we set Lmax = 5. A similar modification can be done to other spins.

using ITensors
using ITensors.Strided
using LinearAlgebra
using Printf
using Random
using ITensors.HDF5
include("S5(modified).jl")

function Schwinger(N,x,y,e,mu,theta,Lmax,nsweeps,maxdim)

  sites = siteinds(n->isodd(n) ? "S=1/2" : "Lmax="*string(Lmax),2*N; conserve_sz=false)

  weight = 100*sqrt(x)       # Weight for excited states, defined as H' = H + w |psi0><psi0|
                            # w has to be at least bigger than the energy gap

  #---------------------------------Hamiltonian for OBC---------------------------------------------
  
  os = OpSum()     # Initialize a sum of operators for the Hamiltonian
  # Calculate Coulomb Hamiltonian (Sz=\sigma/2)
  for j in 1:2:2*N-3
    os2 = OpSum()
    for k in 1:2:j
       os2 .+= 1,"Sz", k
       os2 .+= 1/2*(-1)^k, "Id", k
    end
    os2 .+= theta/(2*pi),"Id", j
    # Square the previous sum 
    J=Int((j+1)/2)
    for k in 1:2*J+1
       os = y*os2*os2[k] + os
    end
  end
  # Add the OBC last term (corresponding to L(N)=L+theta/(2*pi))
  os += y*(theta/(2*pi))^2, "Id", 2*N-1

  # Mass term in the Hamiltonian
  for j in 1:2:2*N-1
    os .+= mu*(-1)^j, "Sz", j
  end

  # XY Hamiltonian 
  for j in 1:2:(2*N - 3)
    os .+= x, "S+", j, "S-", j + 2
    os .+= x, "S-", j, "S+", j + 2
  end

  #-------------------------------------PBC terms---------------------------------------------
  #Coulomb additional terms
  for j in 1:2:2*N-3
    for k in 1:2:j
      os .+= 2*e*y,"Sz", k, "Sz", 2*N
      os .+= y*e*(-1)^k, "Id", k, "Sz", 2*N
    end
    os .+= 2*y*e*theta/(2*pi),"Id", j, "Sz", 2*N
    os .+= y*e, "Id", j, "Sz", 2*N, "Sz", 2*N
  end
  # Additional terms of (L+theta/(2*pi))^2
  os .+= 2*y*e*theta/(2*pi),"Id", 2*N-1, "Sz", 2*N
  os .+= y*e, "Id", 2*N-1, "Sz", 2*N, "Sz", 2*N
  
  #XY additional terms 
  os .+= x*e,"S+", 2*N-1, "S+", 2*N, "S-", 1
  os .+= x*e,"S-", 2*N-1, "S-", 2*N, "S+", 1
  
  # Create the MPO for the Hamiltonian
  H = MPO(os, sites)    

  # Initialize the state with the quantum number (spin) desired to initialize the Lanczos algorithm
  state = [if isodd(n) if isodd(Int((n+1)/2)) "Up" else "Dn" end else "0" end for n in 1:2*N]
  #state = [if n>N/2 "Up" else "Dn" end for n in 1:N]  

  # Create an MPS for the previous state
  psi0_init = MPS(sites, state)


  # Set maximum error allowed when adapting bond dimensions
  cutoff = [0]

  # If DMRG is far from the global minumum then there is no guarantee that DMRG will be able to find the true ground state. 
  # This problem is exacerbated for quantum number conserving DMRG where the search space is more constrained.
  # If this happens, a way out is to turn on the noise term feature to be a very small number.
  noise = [1E-11]

  # Run the DMRG algorithm, returning energy and optimized MPS of ground state
  energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)
  #energy01, psi01 = dmrg(H1, psi0_init; nsweeps, maxdim, cutoff, noise)
  #energy02, psi02 = dmrg(H2, psi0_init; nsweeps, maxdim, cutoff, noise)

  println()
  
  # Initialize the first excited state for Lanczos
  state1 = [if isodd(n) if isodd(Int((n+1)/2)) "Up" else "Dn" end else "-1" end for n in 1:2*N]
  #state1 = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
  psi1_init = MPS(sites,state1)

  # Run DMRG for energy and optimized MPS for first excited state
  energy1,psi1 = dmrg(H,[psi0],psi1_init; nsweeps, maxdim, cutoff, noise, weight)
  #energy1,psi1 = dmrg(H,psi1_init; nsweeps, maxdim, cutoff, noise)
  println()

  # Check if psi1 is orthogonal to psi0
  @show inner(psi1,psi0)

  #Check the eigenvalues of the states with respect to the Sz operator
  Sop = OpSum();
  for i=1:2:2*N-1
    Sop .+= 1,"Sz", i;
  end
  Sop2 = OpSum();
  for i=2:2:2*N
    Sop2 .+= 1,"Sz", i;
  end
  S = MPO(Sop, sites);
  Z = MPO(Sop2, sites);
  S0 = inner(psi0',S,psi0);
  S1 = inner(psi1',S,psi1);

  @show inner(psi0',Z,psi0)
  @show inner(psi1',Z,psi1)
 
  return (energy0,energy1,S0,S1)
end

let
  N = 4    # Number of lattice sites
  x, y, e, theta = 1, 1, 1, pi    # Specify parameters 
  Lmax = 5
  # Plan to do nsweeps DMRG sweeps:
  nsweeps = 10

  # Set maximum MPS bond dimensions for each sweep (truncation m)
  maxdim = 100

  #mu = 1
  me = 0.1
  mu = 2*sqrt(x)*me - 1/4

  E0,E1,S0,S1 = Schwinger(N,x,y,e,mu,theta,Lmax,nsweeps,maxdim)
  println("E1 - E0 = ", string(E1 - E0))
end
