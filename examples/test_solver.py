from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot
import triqs.utility.mpi as mpi

from ctint_tutorial import Solver, Solver2

# Parameters
U = 2.5            # Hubbard interaction
mu = U/2.0         # Chemical potential
half_bandwidth=1.0 # Half bandwidth (energy unit)
beta = 100.0        # Inverse temperature
n_iw = 128         # Number of Matsubara frequencies
n_cycles = 10000   # Number of MC cycles
delta = 0.1        # delta parameter
n_iterations = 1   # Number of DMFT iterations

S = Solver(beta, n_iw) # Initialize the solver

S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

with HDFArchive("data/test_solver.h5",'w') as A:

 for it in range(n_iterations): # DMFT loop
  for name, G0 in S.G0_iw:
   G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

  S.solve(U, delta, n_cycles) # Solve the impurity problem

  G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

  if mpi.is_master_node():
   A['G'] = G_sym # Save G from every iteration to file

S = Solver2(beta, n_iw) # Initialize the solver

S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

with HDFArchive("data/test_solver2.h5",'w') as A:

 for it in range(n_iterations): # DMFT loop
  for name, G0 in S.G0_iw:
   G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0

  S.solve(U, delta, n_cycles) # Solve the impurity problem

  G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

  if mpi.is_master_node():
   A['G'] = G_sym # Save G from every iteration to file
