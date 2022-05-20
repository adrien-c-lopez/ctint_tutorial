from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot
import triqs.utility.mpi as mpi
from numpy.random import randint
#from numpy import ones

from ctint_tutorial import Solver, Solver2

# Parameters
U = 10.               # Hubbard interaction
mu = U/2.0           # Chemical potential
half_bandwidth=1.0   # Half bandwidth (energy unit)
beta = 1.           # Inverse temperature
n_iw = 128           # Number of Matsubara frequencies
n_cycles = 10000     # Number of MC cycles
length_cycles = 50   # Number of MC cycles
n_warmup_cycles = 5000
delta = .1           # delta parameter
N = 10               # Number of solver samples
r = randint(0,100000,size=N)

print(f"beta = {beta}")
print(f"U = {U}")

with HDFArchive("data/test_solver.h5",'w') as A:

   S = Solver(beta, n_iw) # Initialize the solver
   
   for k in range(N):

      S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

      for name, G0 in S.G0_iw:
         G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

      S.solve(U, delta, n_cycles, length_cycles, n_warmup_cycles, seed=r[k])

      G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

      if mpi.is_master_node():
         A[f'G_iw_{k}'] = G_sym 
         A[f'G0_iw_{k}'] = S.G0_iw
         A[f'beta'] = beta
         A[f'U'] = U
         A[f'hist_{k}'] = S.Hist
         A[f'hist_sign_{k}'] = S.Hist_sign
         A[f'd_{k}'] = S.D
         A['N'] = k+1

with HDFArchive("data/test_solver2.h5",'w') as A:

   S = Solver2(beta, n_iw) # Initialize the solver

   for k in range(N):
      S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

      for name, G0 in S.G0_iw:
         G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

      S.solve(U, delta, n_cycles,length_cycles//2, n_warmup_cycles, seed=r[k]) # Solve the impurity problem

      G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

      if mpi.is_master_node():
         A[f'G_iw_{k}'] = G_sym
         A[f'G0_iw_{k}'] = S.G0_iw 
         A[f'beta'] = beta
         A[f'U'] = U
         A[f'hist_{k}'] = S.Hist
         A[f'hist_sign_{k}'] = S.Hist_sign
         A[f'd_{k}'] = S.D
         A[f'd0_{k}'] = S.D0
         A['N'] = k+1
