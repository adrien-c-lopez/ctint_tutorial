from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot
import triqs.utility.mpi as mpi
from numpy.random import rand,randint
#from numpy import ones

from ctint_tutorial import Solver, Solver2

# Parameters
beta = 5             # Inverse temperature
U = 6                # Hubbard interaction
mu = U/2            # Chemical potential
half_bandwidth=1.0   # Half bandwidth (energy unit)
n_iw = 128           # Number of Matsubara frequencies
n_cycles = 10000     # Number of MC cycles
length_cycles = 50   # Number of MC cycles
n_warmup_cycles = 5000

N = 10               # Number of solver samples
r = randint(0,100000,size=N)

nid = 1

for id in range(nid):
   delta = .1
   delta0 = .5

   print(f"id: {id+1} of {nid}")
   print(f"beta = {beta}")
   print(f"U = {U}")
   print(f"delta = {delta}")
   print(f"delta0 = {delta0}")

   if True:
      with HDFArchive(f"data/{id}test_solver.h5",'w') as A:

         S = Solver(beta, n_iw) # Initialize the solver
         
         for k in range(N):

            S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

            for name, G0 in S.G0_iw:
               G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
            print(f"solve # {k+1} of {N}")
            S.solve(U, delta, delta0, n_cycles, length_cycles, n_warmup_cycles, seed=r[k])

            G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

            if mpi.is_master_node():
               A['N'] = k+1
               A[f'beta'] = beta
               A[f'U'] = U
               A[f'delta'] = delta
               A[f'delta0'] = delta0
               A[f'G0_iw_{k}'] = S.G0_iw
               A[f'hist_{k}'] = S.Hist
               A[f'hist_sign_{k}'] = S.Hist_sign
               A[f'n_{k}'] = S.N
               A[f'd_{k}'] = S.D
               A[f'G_iw_{k}'] = G_sym 


   if True:
      with HDFArchive(f"data/{id}test_solver2.h5",'w') as A:

         S = Solver2(beta, n_iw) # Initialize the solver

         for k in range(N):
            S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

            for name, G0 in S.G0_iw:
               G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
            print(f"solve # {k+1} of {N}")
            S.solve(U, delta, delta0, n_cycles,length_cycles, n_warmup_cycles, seed=r[k]) # Solve the impurity problem

            G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

            if mpi.is_master_node():
               A['N'] = k+1
               A[f'beta'] = beta
               A[f'U'] = U
               A[f'delta'] = delta
               A[f'delta0'] = delta0
               A[f'G0_iw_{k}'] = S.G0_iw
               A[f'hist_{k}'] = S.Hist
               A[f'hist_sign_{k}'] = S.Hist_sign
               A[f'n_{k}'] = S.N
               A[f'd_{k}'] = S.D
               A[f'G_iw_{k}'] = G_sym 