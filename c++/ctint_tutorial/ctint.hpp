#include <triqs/gfs.hpp>
#include <triqs/mc_tools.hpp>
#include <triqs/det_manip.hpp>
#include <mpi/mpi.hpp>

// ------------ The main class of the solver -----------------------

using namespace triqs::gfs; 

namespace ctint_tutorial {

  enum spin { up, down };

  struct arg_t {
    double tau; // The imaginary time
    int s;      // The auxiliary spin
  };

// The function that appears in the calculation of the determinant
struct g0bar_tau {
  gf<imtime> const &gt;
  double beta, delta, delta0;
  int s;

  dcomplex operator()(arg_t const &x, arg_t const &y) const {
    if ((x.tau == y.tau)) { // G_\sigma(0^-) - \alpha(\sigma s)
      return 1.0 + gt[0](0, 0) - (delta0 + (2 * (s == x.s ? 1 : 0) - 1) * delta);
    }
    auto x_y = x.tau - y.tau;
    bool b   = (x_y >= 0);
    if (!b) x_y += beta;
    dcomplex res = gt[closest_mesh_pt(x_y)](0, 0);
    return (b ? res : -res); // take into account antiperiodicity
  }
};

  class solver {

  double beta;
  triqs::gfs::block_gf<imfreq> g0_iw, g0tilde_iw, g_iw, M_iw; //M_iw is the interacting Green's function kernel
  triqs::gfs::block_gf<imtime> g0tilde_tau;
  std::vector<double> hist;
  std::vector<dcomplex> hist_sign;
  std::vector<dcomplex> n;
  std::vector<dcomplex> hist_n;
  dcomplex d;
  std::vector<dcomplex> hist_d;

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

  /// Access double occupancy
  //dcomplex D0() {return d0;};

  /// Access double occupancy
  std::vector<dcomplex> N() {return n;};

  /// Access double occupancy histogram
  std::vector<dcomplex> Hist_n() { return hist_n; };

  /// Access double occupancy
  dcomplex D() {return d;};

  /// Access double occupancy histogram
  std::vector<dcomplex> Hist_d() { return hist_d; };

  /// Construct a ctint solver
  solver(double beta_, int n_iw = 1024, int n_tau = 100001);

  /// Method that performs the QMC calculation
  void solve(double U, double delta, double delta0=.5, int n_cycles=10000, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = "",
              int max_time = -1, int seed=34788);

  };

  class solver2 {

  double beta;
  triqs::gfs::block_gf<imfreq> g0_iw, g0tilde_iw, g_iw, M_iw; //M_iw is the interacting Green's function kernel
  triqs::gfs::block_gf<imtime> g0tilde_tau;
  std::vector<double> hist;
  std::vector<dcomplex> hist_sign;
  std::vector<dcomplex> n;
  std::vector<dcomplex> hist_n;
  dcomplex d;
  std::vector<dcomplex> hist_d;

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

  /// Access double occupancy
  std::vector<dcomplex> N() {return n;};
  std::vector<dcomplex> Hist_n() { return hist_n; };

  /// Access double occupancy
  //dcomplex D0() {return d0;};

  /// Access double occupancy
  dcomplex D() {return d;};

  /// Access double occupancy histogram
  std::vector<dcomplex> Hist_d() { return hist_d; };

  /// Construct a ctint solver
  solver2(double beta_, int n_iw = 1024, int n_tau = 100001);

  /// Method that performs the QMC calculation
  void solve(double U, double delta, double delta0=.5, int n_cycles=10000, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = "",
              int max_time = -1, int seed=34788);
  };
  
} // namespace ctint_tutorial
