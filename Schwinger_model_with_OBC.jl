# This code gives the spectrum of the Schwinger model with Open Boundary Conditions at fixed N and x.
# Note that the spectrum depends on the site that we start with, that is, there is a difference if we start at site = 0 or at site = 1 (since we are 
# breaking translation by one site symmetry).

using ITensors, ITensorMPS
using ITensors.Strided
using LinearAlgebra
using Printf
using Random
using HDF5

function Schwinger(N,x,y,mu,theta,nsweeps,maxdim,spin)

  sites = siteinds("S=1/2", N; conserve_qns=true)  # Dimensions of each site (spin 1/2 in this case) 
  #Notice that eigenvalues here are +-2Sz: +1 or -1 in this case. Here we set projection command in the conserved numbers parameter.

  weight = 1000*x    # Weight for excited states, defined as H' = H + w |psi0><psi0|
                            # w has to be at least bigger than the energy gap

  os = OpSum()     # Initialize a sum of operators for the Hamiltonian
  osP = OpSum()    # Initialize a sum of operators for P = ̄ψψ
  osQ = OpSum()    # Initialize a sum of operators for Q = i ψ γ_5 ψ

  # Calculate Coulomb Hamiltonian (Sz=\sigma/2) starting at site = 1
  #for j in 1:N-1
  #  os2 = OpSum()
  #  for k in 1:j
  #     os2 .+= 1,"Sz", k
  #     os2 .+= 1/2*(-1)^(k), "Id", k
  #  end
  #  os2 .+= theta/(2*pi),"Id", j
  #  # Square the previous sum 
  #  for k in 1:(2*j+1)
  #     os = y*os2*os2[k] + os
  #  end
  #end
  # Add the last term (corresponding to L(N)=theta/(2*pi))
  #os += y*(theta/(2*pi))^2, "Id", N

  # Mass term in the Hamiltonian
  #for j in 1:N
  #  os .+= mu*(-1)^(j), "Sz", j
  #  #os .+= mu/2, "Id", j
  #  osP .+= (-1)^j, "Sz", j
  #  osP .+= 1/2, "Id", j
  #end

  # XY Hamiltonian 
  #for j in 1:(N - 1)
  #  os += x, "S+", j, "S-", j + 1
  #  os += x, "S-", j, "S+", j + 1
  #  osQ += (-1)^j, "S+", j, "S-", j + 1
  #  osQ += (-1)^j, "S-", j, "S+", j + 1
  #end

  # Calculate Coulomb Hamiltonian (Sz=\sigma/2) starting at site = 0
  for j in 1:N-1
   os2 = OpSum()
    for k in 1:j
       os2 .+= 1,"Sz", k
       os2 .+= 1/2*(-1)^(k-1), "Id", k
    end
    os2 .+= theta/(2*pi),"Id", j
    # Square the previous sum 
    for k in 1:(2*j+1)
       os = y*os2*os2[k] + os
    end
  end
  
  # Mass term in the Hamiltonian
  for j in 1:N
    os .+= mu*(-1)^(j-1), "Sz", j
    #os .+= mu/2, "Id", j
    osP .+= (-1)^(j-1), "Sz", j
    osP .+= 1/2, "Id", j
  end

  # XY Hamiltonian 
  for j in 1:(N - 1)
    os += x, "S+", j, "S-", j + 1
    os += x, "S-", j, "S+", j + 1
    osQ += (-1)^(j-1), "S+", j, "S-", j + 1
    osQ += (-1)^(j-1), "S-", j, "S+", j + 1
  end
  
  # Create the MPO for the Hamiltonian
  H = MPO(os, sites)    
 
  # Create the MPO for P and Q
  P = MPO(osP, sites)
  Q = MPO(osQ, sites)

  # Initialize the state with the quantum number (spin) desired to initialize the Lanczos algorithm
  if spin==0 
    state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    state1 = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    state2 = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    state3 = [isodd(n) ? "Up" : "Dn" for n in 1:N]
  elseif spin==1
    state = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
    state1 = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
    state2 = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
    state3 = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
  elseif spin==-1
    state = [if n>1 if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" end for n in 1:N]
    state1 = [if n>1 if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" end for n in 1:N]
    state2 = [if n>1 if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" end for n in 1:N]
    state3 = [if n>1 if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" end for n in 1:N]
  end

  # Create an MPS for previous states
  psi0_init = MPS(sites, state)
  psi1_init = MPS(sites,state1)
  psi2_init = MPS(sites,state2)
  psi3_init = MPS(sites,state3)
  # Set maximum error allowed when adapting bond dimensions
  cutoff = [1E-15]

  # If DMRG is far from the global minumum then there is no guarantee that DMRG will be able to find the true ground state. 
  # This problem is exacerbated for quantum number conserving DMRG where the search space is more constrained.
  # If this happens, a way out is to turn on the noise term feature to be a very small number.
  noise = [1E-15]

  eigsolve_krylovdim = 10 # Maximum dimension of Krylov space to locally solve problem. Try setting to a higher
                          #     value if convergence is slow or the Hamiltonian is close to a critical point.
  eigsolve_maxiter = 1    # Number of times Krylov space can be rebuild

  eigsolve_tol = 1e-14    # Krylov eigensolver tolerance

  # Run the DMRG algorithm, returning energy and optimized MPS of ground state
  energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise, eigsolve_krylovdim, eigsolve_tol, eigsolve_maxiter)
  #energy01, psi01 = dmrg(H1, psi0_init; nsweeps, maxdim, cutoff, noise)
  #energy02, psi02 = dmrg(H2, psi0_init; nsweeps, maxdim, cutoff, noise)

  println()

  # Run DMRG for energy and optimized MPS for first excited state
  energy1,psi1 = dmrg(H,[psi0],psi1_init; nsweeps, maxdim, cutoff, noise, weight, eigsolve_krylovdim, eigsolve_tol, eigsolve_maxiter)
  #energy1,psi1 = dmrg(H,psi1_init; nsweeps, maxdim, cutoff, noise)
  println()

  # Run DMRG for energy and optimized MPS for second excited state
  energy2,psi2 = dmrg(H,[psi0,psi1],psi2_init; nsweeps, maxdim, cutoff, noise, weight)
  #energy2,psi2 = dmrg(H,[psi0],psi2_init; nsweeps, maxdim, cutoff, noise, weight)

  println()

  # Run DMRG for energy and optimized MPS for third excited state
  energy3,psi3 = dmrg(H,[psi0,psi1,psi2],psi3_init; nsweeps, maxdim, cutoff, noise, weight)

  # Check if psi1 is orthogonal to psi0
  @show inner(psi1,psi0)
  @show inner(psi2,psi0)
  @show inner(psi3,psi0)
  @show inner(psi2,psi1)
  @show inner(psi3,psi1)
  @show inner(psi3,psi2)

  CC = sqrt(x)*inner(psi0',P,psi0)/N;     #Chiral condensate calculation
  OP = sqrt(x)*inner(psi0',Q,psi0)/N;     #Order parameter calculation 

  #Check the eigenvalues of the states with respect to the Sz operator
  Sop = OpSum();
  for i=1:N
    Sop += 1,"Sz",i;
  end

  #Create the Z2 operator which is a product of all X_i operators 
  Str=OpSum();
  Str=1,"Id",1;
  Sop2 = Str[1];
  for i=1:N
    S=OpSum();
    S+=2,"Sz",i;
    Sop2 = Sop2 * S[1];  
  end
  S = MPO(Sop, sites);
  S0 = inner(psi0',S,psi0);
  S1 = inner(psi1',S,psi1);
  S2 = inner(psi2',S,psi2);
  S3 = inner(psi3',S,psi3);
  Sx = MPO(Sop2, sites);
  @show inner(psi0',Sx,psi0)
  @show inner(psi1',Sx,psi1)
  @show inner(psi2',Sx,psi2)
  @show inner(psi3',Sx,psi3)
  return (energy0,energy1,energy2,energy3)
end
let
  N = 8     # Number of lattice sites
  x, y, theta = 10, 1, 0    # Specify parameters 
  spin = 0

  # Plan to do nsweeps DMRG sweeps ground state:
  nsweeps = 5

  # Set maximum MPS bond dimensions for each sweep (truncation m)
  maxdim = 500

  me=0.1
  #mu=0
  mu = 2*sqrt(x)*me - 1/4
  E0,E1,E2,E3 = Schwinger(N,x,y,mu,theta,nsweeps,maxdim,spin)   
  println("E0 = ", string(E0),"\nE1 = ", string(E1),"\nE2 = ", string(E2),"\nE3 = ", string(E3),"\n")
end
