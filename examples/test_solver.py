from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot
import triqs.utility.mpi as mpi
from numpy import array
from numpy.random import rand,randint
#from numpy import ones
from time import time
from ctint_tutorial import Solver, Solver2

# Parameters
delta = .15
delta0 = .5
half_bandwidth=1.0   # Half bandwidth (energy unit)
n_iw = 128           # Number of Matsubara frequencies
n_cycles = 10000     # Number of MC cycles
length_cycles = 50   # Number of MC cycles
n_warmup_cycles = 5000

N = 10               # Number of solver samples
r = randint(0,100000,size=N)

nid = 5
id0 = 20


for id in range(id0,id0+nid):
   beta = (id-id0+1)/nid*6              # Inverse temperature
   U = 4                # Hubbard interaction
   mu = U/2            # Chemical potential

   print(f"id: {id-id0+1} of {nid}")

   if True:
      with HDFArchive(f"data/{id}test_solver.h5",'w') as A:

         S = Solver(beta, n_iw) # Initialize the solver

         Ns = 0
         hist = []
         hist_sign = []
         n = []
         hist_n = []
         d = []
         hist_d = []
         G0_iw = []
         G_iw = []
         M_iw = []
         pert_k = []
         Mk_iw = []
         t = []

         for k in range(N):
            k0 = (2*k+1)%8

            S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

            for name, G0 in S.G0_iw:
               G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
            print(f"solve # {k+1} of {N}")
            t0 = time()
            S.solve(U, delta, delta0, k0, n_cycles, length_cycles, n_warmup_cycles, seed=r[k])
            dt = time()-t0
            G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

            Ns += 1
            hist.append(S.Hist)
            hist_sign.append(S.Hist_sign)
            n.append(S.N)
            hist_n.append(S.Hist_n)
            d.append(S.D)
            hist_d.append(S.Hist_d)
            G0_iw.append(S.G0_iw)
            G_iw.append(G_sym)
            M_iw.append(S.M_iw)
            pert_k.append(k0) 
            Mk_iw.append(S.Mk_iw) 
            t.append(dt)

         if mpi.is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'delta'] = delta
            A[f'delta0'] = delta0
            A[f'hist'] = hist
            A[f'hist_sign'] = hist_sign
            A[f'n'] = n
            A[f'hist_n'] = hist_n
            A[f'd'] = d
            A[f'hist_d'] = hist_d
            A[f'G0_iw'] = G0_iw
            A[f'G_iw'] = G_iw 
            A[f'M_iw'] = M_iw 
            A[f'k'] = pert_k 
            A[f'Mk_iw'] = Mk_iw
            A[f't'] = t

   if True:
      with HDFArchive(f"data/{id}test_solver2.h5",'w') as A:

         S = Solver2(beta, n_iw) # Initialize the solver
         
         Ns = 0
         hist = []
         hist_sign = []
         n = []
         hist_n = []
         d = []
         hist_d = []
         G0_iw = []
         G_iw = []
         M_iw = []
         pert_k = []
         Mk_iw = []
         t = []

         for k in range(N):
            k0 = (2*k+1)%8

            S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

            for name, G0 in S.G0_iw:
               G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
            print(f"solve # {k+1} of {N}")
            t0 = time()
            try:
               S.solve(U, delta, delta0, k0, n_cycles,length_cycles, n_warmup_cycles, seed=r[k]) # Solve the impurity problem
               dt = time()-t0
               G_sym = (S.G_iw['up']+S.G_iw['down'])/2 # Impose paramagnetic solution

               Ns += 1
               hist.append(S.Hist)
               hist_sign.append(S.Hist_sign)
               n.append(S.N)
               hist_n.append(S.Hist_n)
               d.append(S.D)
               hist_d.append(S.Hist_d)
               G0_iw.append(S.G0_iw)
               G_iw.append(G_sym)
               M_iw.append(S.M_iw)
               pert_k.append(k0) 
               Mk_iw.append(S.Mk_iw) 
               t.append(dt)
            except RuntimeError:
               pass

         
         if mpi.is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'delta'] = delta
            A[f'delta0'] = delta0
            A[f'hist'] = array(hist)
            A[f'hist_sign'] = array(hist_sign)
            A[f'n'] = array(n)
            A[f'hist_n'] = array(hist_n)
            A[f'd'] = array(d)
            A[f'hist_d'] = array(hist_d)
            A[f'G0_iw'] = array(G0_iw)
            A[f'G_iw'] = array(G_iw)
            A[f'M_iw'] = array(M_iw) 
            A[f'k'] = array(pert_k)
            A[f'Mk_iw'] = array(Mk_iw)
            A[f't'] = array(t)