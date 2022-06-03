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
beta = 15.            # Inverse temperature
U = 4.                # Hubbard interaction
mu = U/2              # Chemical potential
delta = .1
delta0 = .5
half_bandwidth=1.0   # Half bandwidth (energy unit)
n_iw = 128           # Number of Matsubara frequencies
n_cycles0 = 10000     # Number of MC cycles
length_cycles = 50   # Number of MC cycles
n_warmup_cycles = 5000

N = 10               # Number of solver samples
r = randint(0,100000,size=N)

nid = 1
id0 = 0

for id in range(id0,id0+nid):
   print(f"id: {id-id0+1} of {nid}")

   if False:
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
         warn = 0

         for k in range(N):
            k0 = (2*k+1)%8

            S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function

            for name, G0 in S.G0_iw:
               G0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
            print(f"solve # {k+1} of {N}")
            t0 = time()
            try:
               S.solve(U, delta, delta0, k0, n_cycles0, length_cycles, n_warmup_cycles, seed=r[k])
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
            except Warning:
               warn += 1

         if is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'delta'] = delta
            A[f'delta0'] = delta0
            A[f'hist'] = pad(hist)
            A[f'hist_sign'] = pad(hist_sign)
            A[f'n'] = array(n)
            A[f'hist_n'] = pad(hist_n)
            A[f'd'] = array(d)
            A[f'hist_d'] = pad(hist_d)
            A[f'G0_iw'] = (G0_iw)
            A[f'G_iw'] = (G_iw)
            A[f'M_iw'] = (M_iw) 
            A[f'k'] = array(pert_k)
            A[f'Mk_iw'] = (Mk_iw)
            A[f't'] = array(t)

   if True:
      n_cycles = n_cycles0
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
         n_meas = []

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
               n_meas.append(n_cycles)
            except RuntimeError:
               pass

         
         if is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'delta'] = delta
            A[f'delta0'] = delta0
            A[f'hist'] = pad(hist)
            A[f'hist_sign'] = pad(hist_sign)
            A[f'n'] = array(n)
            A[f'hist_n'] = pad(hist_n)
            A[f'd'] = array(d)
            A[f'hist_d'] = pad(hist_d)
            A[f'G0_iw'] = (G0_iw)
            A[f'G_iw'] = (G_iw)
            A[f'M_iw'] = (M_iw) 
            A[f'k'] = array(pert_k)
            A[f'Mk_iw'] = (Mk_iw)
            A[f't'] = array(t)
            A[f'n_cycles'] = n_cycles