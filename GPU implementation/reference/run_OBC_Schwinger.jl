using ITensors, ITensorMPS
using ITensors.Strided
using LinearAlgebra
using Random
include("OBC_Schwinger_dmrg.jl")
#using ITensors.HDF5
BLAS.set_num_threads(4)

II = 0 + 1im;
X = [0 1;1 0];
Y = [0 -II; II 0];
sp = [0 1;0 0];
sm = [0 0;1 0];
Z = [1 0;0 -1];
Id = [1 0;0 1];


#let
ITensors.enable_threaded_blocksparse()
x = 10  
y = 1    # Specify parameters
g = 0.0  # Quartic interaction

# Plan to do nsweeps DMRG sweeps:
nsweeps = 10

# Set maximum MPS bond dimensions for each sweep (truncation m)
maxdim = 100
theta = pi
spin = 0  #Spin projection
num_states = 1

mass = 0.33

# file = open("/home/earguell/Julia/Schwinger/OBC/massive_data/precision_data/SpecN_MSM_S" *string(spin) * "_nsweeps" *
# string(nsweeps) * "_m" * string(maxdim) * "_x" * string(x) * "_theta" * string(theta) * "_g" * string(g) * ".txt","w")

# for k in 2:10
#N = 100 + 50*k   # Number of lattice sites
N = 20

conserved = false
sites = siteinds("S=1/2", N; conserve_qns=conserved) 

 #N-even
  # Initialize the state with the quantum number (spin) desired to initialize the Lanczos algorithm
  if spin==0
    state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
  elseif spin==1
    state = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Up" end for n in 1:N]
  elseif spin==-1
    state = [if n>1 if isodd(n) "Up" else "Dn" end elseif (n==1) "Dn" end for n in 1:N]
  elseif spin==2
    state = [if (n!=N && n!=N-2) if isodd(n) "Up" else "Dn" end elseif (n==N || n==N-2) "Up" end for n in 1:N]
  elseif spin==-2
    state = [if (n!=1 && n!=3) if isodd(n) "Up" else "Dn" end elseif (n==1 || n==3) "Dn" end for n in 1:N]
  end

# N-odd
  # Initialize the state with the quantum number (spin) desired to initialize the Lanczos algorithm
  #if spin==1/2
  #  state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
  #elseif spin==-1/2
  #  state = [if n<N if isodd(n) "Up" else "Dn" end elseif (n==N) "Dn" end for n in 1:N]
  #end

  # Create an MPS for previous states
  psi_init = MPS(sites, state)

psi0_init = psi_init;

os = OpSum();   # Initialize a sum of operators for the Hamiltonian
 for j in 1:(N - 1)
    os += 1, "S+", j, "S-", j + 1
    os += 1, "S-", j, "S+", j + 1
    os += -4*g, "Sz", j, "Sz", j + 1
 end
 # Create the MPO for the Hamiltonian
 Hxy1 = MPO(os, sites);

 # Mass term in the Hamiltonian
 os = OpSum();
 for j in 1:N
     os .+= (-1)^(j), "Sz", j
     #os .+= mu/2, "Id", j
 end
 Hm1 = MPO(os, sites);

 fn_gauge = (ham, n) -> schwinger_gauge(ham, n, N, theta);
 Hc1 = get_MPO_neutral(fn_gauge, sites, 3, conserved);

 E,psi = Schwinger(N,x,y,mass,theta,nsweeps,maxdim,spin,psi0_init,Hc1,Hm1,Hxy1,num_states);

    #Check normalization of states
      println("\nCheck normalization of states")
      for i in 1:num_states
        println("<psi",string(i-1),"|psi",string(i-1),"> = ", inner(psi[i]',psi[i]))
      end

      #Check orthogonalization of states
      println("\nCheck orthogonalization of states")
      for i in 1:num_states
        for j in (i+1):num_states
          println("<psi",string(i-1),"|psi",string(j-1),"> = ", inner(psi[i]',psi[j]))
        end
      end

      println("N = ", string(N), "\ny = ", string(y), "\n(m/e) = ", string(mass), "\nValue of g = ", string(g))
      #write(file,"N = ", string(N), "\ny = ", string(y), "\n(m/e) = ", string(mass), "\nValue of g = ", string(g),"\n")
      for i in 1:num_states
        println("E",string(i-1)," = ", string(E[i]))
        # write(file, "E",string(i-1)," = ", string(E[i]), "\n")
      end

#Save MPS
#fpsi0 = h5open("/home/earguell/Julia/Schwinger/OBC/massive_data/precision_data/MPS/psi0_N"*string(N)*
#  "_m"*string(maxdim)* "_x"*string(x)* "_theta" * string(theta)* "_mass"*string(z)*".h5","w")
#  write(fpsi0,"psi0",psi0)
#close(fpsi0)

#end
#close(file)
ITensors.disable_threaded_blocksparse()
#end


### This part returns the tensor T_middle in the middle of the chain, namely extracts the psi0[N/2] part ###

L = length(psi[1])            # L = N
jmid = div(L, 2)              # jmid = N/2

# Put the orthogonality center at site N/2 
orthogonalize!(psi[1], jmid)

T_middle = deepcopy(psi[1][jmid])

@show typeof(T_middle)       # should be a Tensor


### If need physical dimension = 4 (combination of 2 adjacent sites), we combine the Tensor sites at N/2 and N/2 + 1 ###
T_dimer = deepcopy(psi[1][jmid]*psi[1][jmid+1])

s1 = siteind(psi[1],jmid)
s2 = siteind(psi[1],jmid + 1)

C = combiner(s1, s2; tags="Combined_Site,n=$(jmid)")
T_fused = T_dimer * C

@show typeof(T_fused)