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
U = 4.                # Hubbard interaction
doping = 1.9
delta = .1
delta0 = .85
half_bandwidth=1.0   # Half bandwidth (energy unit)
n_iw = 128           # Number of Matsubara frequencies
n_cycles = 20000     # Number of MC cycles
length_cycle = 50   # Number of MC cycles
n_warmup_cycles = 5000

N = 3               # Number of solver samples
r = randint(0,100000,size=N)

nid = 5
id0 = 0

for id in range(id0,id0+nid):
   print(f"id: {id-id0+1} of {nid}")
   beta = 2.+(id-id0+1)/nid*7             # Inverse temperature

   if True:
      with HDFArchive(f"data/{id}solver.h5",'w') as A:

         S = Solver(beta, n_iw) # Initialize the solver

         Ns = 0
         hist = []
         hist_sign = []
         n = []
         Density = []
         Mu = []
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
            print(f"solve # {k+1} of {N}")
            k0 = (2*k+1)%8
            density = []
            mu = []
            mu0 = U/2
            dmu = U/2
            t0 = time()
            while abs(dmu) > 1e-2:

               S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function
               for name, G0 in S.G0_iw:
                  G0 << inverse(iOmega_n + mu0 - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
               delta0 = .5+.35*(mu0-U/2)/5
               try:
                  S.solve(U, delta, delta0, k0, n_cycles=n_cycles, length_cycle=length_cycle, n_warmup_cycles=n_warmup_cycles, seed=r[k])
               except RuntimeError:
                  pass
               except KeyboardInterrupt:
                  raise KeyboardInterrupt
               n0 = sum(S.N.real)
               print(f"n0: {n0}\nmu0: {mu0}")
               density.append(n0)
               mu.append(mu0)
               if (n0-doping)*dmu >= 0:
                  dmu /= -10
               mu0 += dmu
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
            Mu.append(mu)
            Density.append(density)

         if is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'mu'] = Mu
            A[f'density'] = Density
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
            A[f'n_warmup_cylces'] = n_warmup_cycles

   if True:
      n_cycles2 = n_cycles//5
      n_warmup_cycles2 = n_warmup_cycles//10
      with HDFArchive(f"data/{id}solver2.h5",'w') as A:

         S = Solver(beta, n_iw) # Initialize the solver

         Ns = 0
         hist = []
         hist_sign = []
         n = []
         Density = []
         Mu = []
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
            print(f"solve # {k+1} of {N}")
            k0 = (2*k+1)%8

            density = []
            mu = []
            mu0 = U/2
            dmu = U/2
            t0 = time()
            while abs(dmu) > 1e-2:

               S.G_iw << SemiCircular(half_bandwidth) # Initialize the Green's function
               for name, G0 in S.G0_iw:
                  G0 << inverse(iOmega_n + mu0 - (half_bandwidth/2.0)**2 * S.G_iw[name] ) # Set G0
               delta0 = .5+.35*(mu0-U/2)/5
               try:
                  S.solve(U, delta, delta0, k0, n_cycles=n_cycles2, length_cycle=length_cycle, n_warmup_cycles=n_warmup_cycles2, seed=r[k])
               except RuntimeError:
                  pass
               except KeyboardInterrupt:
                  raise KeyboardInterrupt
               n0 = sum(S.N.real)
               print(f"n0: {n0}\nmu0: {mu0}")
               density.append(n0)
               mu.append(mu0)
               if (n0-doping)*dmu >= 0:
                  dmu /= -10
               mu0 += dmu

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
            Mu.append(mu)
            Density.append(density)            

         if is_master_node():
            A['N'] = Ns
            A[f'beta'] = beta
            A[f'U'] = U
            A[f'mu'] = Mu
            A[f'density'] = Density
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
            A[f'n_cycles'] = n_cycles2
            A[f'n_warmup_cylces'] = n_warmup_cycles2