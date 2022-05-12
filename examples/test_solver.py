from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot
import triqs.utility.mpi as mpi
from random import randint

from ctint_tutorial import Solver, Solver2

# Parameters
U = 10.             # Hubbard interaction
mu = U/2.0         # Chemical potential
half_bandwidth=1.0 # Half bandwidth (energy unit)
beta = 8.          # Inverse temperature
n_iw = 128         # Number of Matsubara frequencies
n_cycles = 10000   # Number of MC cycles
delta = .1         # delta parameter
N = 10             # Number of solver samples

print(f"beta = {beta}")
print(f"U = {U}")

with HDFArchive("data/test_solver.h5",'w') as A:

   S = Solver(beta, n_iw) # Initialize the solver

   S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

   for name, G0 in S.G0_iw:
      G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

   S.solve(U, delta, n_cycles) # Solve the impurity problem

   G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

   if mpi.is_master_node():
      A['G'] = G_sym # Save G from every iteration to file
      A['hist'] = S.Hist
      A[f'hist_sign'] = S.Hist_sign
      A['beta'] = beta

with HDFArchive("data/test_solver2.h5",'w') as A:

   S = Solver2(beta, n_iw) # Initialize the solver

   for k in range(N):
      S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

      for name, G0 in S.G0_iw:
         G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

      S.solve(U, delta, n_cycles, seed=10*(k+7)) # Solve the impurity problem

      G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

      if mpi.is_master_node():
         A[f'G_{k}'] = G_sym # Save G from every iteration to file
         A[f'hist_{k}'] = S.Hist
         A[f'hist_sign_{k}'] = S.Hist_sign
         A['N'] = k+1
