from triqs.gf import SemiCircular, inverse, iOmega_n
from h5 import HDFArchive
from triqs.utility.mpi import is_master_node
from numpy import array,zeros
from numpy.random import randint#, rand
from time import time
from ctint_tutorial import Solver, Solver2

# Utils
def pad(u):
   n = max([len(x) for x in u])
   res = zeros((len(u),n),dtype=type(u[0][0]))
   for i in range(len(u)):
      res[i,:len(u[i])] = u[i]
   return res

# Parameters
beta = 6.             # Inverse temperature
U = 4.                # Hubbard interaction
mu = U/2              # Chemical potential
delta = .1
delta0 = .5
half_bandwidth=1.0   # Half bandwidth (energy unit)
n_iw = 128           # Number of Matsubara frequencies
n_cycles = 10000     # Number of MC cycles
length_cycle = 50   # Number of MC cycles
n_warmup_cycles = 5000

N = 10               # Number of solver samples
r = randint(0,100000,size=N)

S = Solver2(beta, n_iw) # Initialize the solver
nobc = 1
err = 0

for k in range(N):

   S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function
   for name, G0 in S.G0_iw:
      G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

   print(f"solve # {k+1} of {N}")
   try:
      S.solve(U, delta, delta0, nobc=nobc, n_cycles=n_cycles,length_cycle=length_cycle, n_warmup_cycles=n_warmup_cycles, seed=r[k]) # Solve the impurity problem
   except RuntimeError:
      err +=1
   print(f"Total errors : {err} of ")