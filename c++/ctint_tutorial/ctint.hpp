#include <triqs/gfs.hpp>
#include <triqs/mc_tools.hpp>
#include <triqs/det_manip.hpp>
#include <mpi/mpi.hpp>

// ------------ The main class of the solver -----------------------

using namespace triqs::gfs; 

namespace ctint_tutorial {

  enum spin { up, down };

  class solver {

  double beta;
  triqs::gfs::block_gf<imfreq> g0_iw, g0tilde_iw, g_iw, M_iw; //M_iw is the interacting Green's function kernel
  triqs::gfs::block_gf<imtime> g0tilde_tau;
  std::vector<double> hist;
  std::vector<dcomplex> hist_sign;

  public:
  /// Access non-interacting Matsubara Green function
  triqs::gfs::block_gf_view<imfreq> G0_iw() { return g0_iw; }

  /// Access non-interacting imaginary-time Green function
  triqs::gfs::block_gf_view<imtime> G0_tau() { return g0tilde_tau; }

  /// Access interacting Matsubara Green function
  triqs::gfs::block_gf_view<imfreq> G_iw() { return g_iw; }

  /// Access order histogram
  std::vector<double> Hist() { return hist; };
  std::vector<dcomplex> Hist_sign() { return hist_sign; };

  /// Construct a ctint solver
  solver(double beta_, int n_iw = 1024, int n_tau = 100001);

  /// Method that performs the QMC calculation
  void solve(double U, double delta, int n_cycles, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = "",
              int max_time = -1);

  };

  class solver2 {

  double beta;
  triqs::gfs::block_gf<imfreq> g0_iw, g0tilde_iw, g_iw, M_iw; //M_iw is the interacting Green's function kernel
  triqs::gfs::block_gf<imtime> g0tilde_tau;
  std::vector<double> hist;
  std::vector<dcomplex> hist_sign;

  public:
  /// Access non-interacting Matsubara Green function
  triqs::gfs::block_gf_view<imfreq> G0_iw() { return g0_iw; }

  /// Access non-interacting imaginary-time Green function
  triqs::gfs::block_gf_view<imtime> G0_tau() { return g0tilde_tau; }

  /// Access interacting Matsubara Green function
  triqs::gfs::block_gf_view<imfreq> G_iw() { return g_iw; }

  /// Access order histogram
  std::vector<double> Hist() { return hist; };
  std::vector<dcomplex> Hist_sign() { return hist_sign; };

  /// Construct a ctint solver
  solver2(double beta_, int n_iw = 1024, int n_tau = 100001);

  /// Method that performs the QMC calculation
  void solve(double U, double delta, int n_cycles, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = "",
              int max_time = -1, int seed=34788);
  };
  
} // namespace ctint_tutorial
