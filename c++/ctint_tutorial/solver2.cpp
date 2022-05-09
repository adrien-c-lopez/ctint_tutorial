#include "ctint.hpp"

#include <triqs/mc_tools.hpp>
#include <triqs/det_manip.hpp>
#include <mpi/mpi.hpp>

using namespace triqs::gfs; 
using namespace triqs::arrays; 
using namespace ctint_tutorial;

// --------------- The QMC configuration ----------------

// Argument type of g0bar
struct arg_t {
  double tau; // The imaginary time
  int s;      // The auxiliary spin
};

// The function that appears in the calculation of the determinant
struct g0bar_tau {
  gf<imtime> const &gt;
  double beta, delta;
  int s;

  dcomplex operator()(arg_t const &x, arg_t const &y) const {
    if ((x.tau == y.tau)) { // G_\sigma(0^-) - \alpha(\sigma s)
      return 1.0 + gt[0](0, 0) - (0.5 + (2 * (s == x.s ? 1 : 0) - 1) * delta);
    }
    auto x_y = x.tau - y.tau;
    bool b   = (x_y >= 0);
    if (!b) x_y += beta;
    dcomplex res = gt[closest_mesh_pt(x_y)](0, 0);
    return (b ? res : -res); // take into account antiperiodicity
  }
};

// The Monte Carlo configuration
struct configuration {
  // M-matrices for up and down
  std::vector<triqs::det_manip::det_manip<g0bar_tau>> Mmatrices_even;
  std::vector<triqs::det_manip::det_manip<g0bar_tau>> Mmatrices_odd;

  int perturbation_order() const { return Mmatrices_odd[up].size(); } //odd perturbation order

  configuration(block_gf<imtime> &g0tilde_tau, double beta, double delta) {
    // Initialize the M-matrices. 100 is the initial matrix size
    for (auto spin : {up, down})
      Mmatrices_even.emplace_back(g0bar_tau{g0tilde_tau[spin], beta, delta, spin}, 100);
    for (auto spin : {up, down})
      Mmatrices_odd.emplace_back(g0bar_tau{g0tilde_tau[spin], beta, delta, spin}, 100);
    double tau = 0;
    int s = 1;
    for (auto &d : Mmatrices_odd) 
      d.insert(0, 0, {tau,s}, {tau,s}); // Insert first time
  }
};

// ------------ QMC move : inserting a vertex ------------------

/*

struct move_insert {
  configuration *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;

  dcomplex attempt() { // Insert an interaction vertex at time tau with aux spin s
    double tau     = rng(beta);
    int s          = rng(2);
    auto k         = config->perturbation_order();
    auto det_ratio = config->Mmatrices[up].try_insert(k, k, {tau, s}, {tau, s}) * config->Mmatrices[down].try_insert(k, k, {tau, s}, {tau, s});
    return -beta * U / (k + 1) * det_ratio; // The Metropolis ratio
  }

  dcomplex accept() {
    for (auto &d : config->Mmatrices) d.complete_operation(); // Finish insertion
    return 1.0;
  }

  void reject() {
    for (auto &d : config->Mmatrices) d.reject_last_try(); // Finish insertion
  }
};

*/

struct move_insert {
  configuration *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;

  dcomplex attempt() { // Insert an interaction vertex at time tau with aux spin s
    auto k         = config->perturbation_order();
    auto t         = config->Mmatrices_odd[up].get_x(k-1);
    double tau     = t.tau;
    double tau0    = rng(beta);
    double tau1    = rng(beta);
    int s          = t.s;
    int s0         = rng(2);
    int s1         = rng(2);

    auto det_ratio = config->Mmatrices_odd[up].try_remove(k-1, k-1) * config->Mmatrices_odd[down].try_remove(k-1, k-1);
    for (auto &d : config->Mmatrices_odd) d.reject_last_try();

    auto det_ratio_even = config->Mmatrices_even[up].try_insert2(k,k+1, k,k+1, {tau, s},{tau1, s1}, {tau0, s0},{tau0, s0}) 
                  *config->Mmatrices_even[down].try_insert2(k,k+1, k,k+1, {tau, s},{tau, s1}, {tau0, s0},{tau0, s0});
    
    auto det_ratio_odd = config->Mmatrices_odd[up].try_insert2(k,k+1, k,k+1, {tau0, s0},{tau1, s1}, {tau0, s0},{tau1, s1}) 
                  *config->Mmatrices_odd[down].try_insert2(k,k+1, k,k+1, {tau0, s0},{tau1, s1}, {tau0, s0},{tau1, s1});
                  
    return -beta * beta * U * U / (k + 2) /(k + 1)
          * ((k+2) * det_ratio * det_ratio_even / beta / U+det_ratio_odd)/(k * det_ratio / beta / U + 1); // The Metropolis ratio
  }

  dcomplex accept() {
    for (auto &d : config->Mmatrices_even)
      d.complete_operation(); // Finish insertion
    for (auto &d : config->Mmatrices_odd)
      d.complete_operation(); // Finish insertion
    return 1.0;
  }

  void reject() {
    for (auto &d : config->Mmatrices_even)
      d.reject_last_try(); // Finish insertion
    for (auto &d : config->Mmatrices_odd)
      d.reject_last_try(); // Finish insertion
  }
};

// ------------ QMC move : deleting a vertex ------------------

/*

struct move_remove {
  configuration *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;

  dcomplex attempt() {
    auto k = config->perturbation_order();
    if (k <= 0) return 0;    // Config is empty, trying to remove makes no sense
    int p          = rng(k); // Choose one of the operators for removal
    auto det_ratio = config->Mmatrices[up].try_remove(p, p) * config->Mmatrices[down].try_remove(p, p);
    return -k / (beta * U) * det_ratio; // The Metropolis ratio
  }

  dcomplex accept() {
    for (auto &d : config->Mmatrices) d.complete_operation();
    return 1.0;
  }

  void reject() {
    for (auto &d : config->Mmatrices) d.reject_last_try(); // Finish insertion
  }                                                        // Nothing to do
};

*/

struct move_remove {
  configuration *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;

  dcomplex attempt() {
    auto k = config->perturbation_order();
    if (k <= 1) return 0;    // Config is empty, trying to remove makes no sense
    int p0          = rng(k); // Choose one of the operators for removal
    int p1          = rng(k); // Choose one of the operators for removal
    if (p1 == p0) {
      if (p0 == 0) p1 = 1;
      else p1 = p0-1;
    }

    auto det_ratio = config->Mmatrices_odd[up].try_remove(k-1, k-1) * config->Mmatrices_odd[down].try_remove(k-1, k-1);
    for (auto &d : config->Mmatrices_odd) d.reject_last_try();

    auto det_ratio_even = config->Mmatrices_even[up].try_remove2(p0, p1, p0, p1) 
                       *config->Mmatrices_even[down].try_remove2(p0, p1, p0, p1);
    
    auto det_ratio_odd = config->Mmatrices_odd[up].try_remove2(p0, p1, p0, p1)  
                      *config->Mmatrices_odd[down].try_remove2(p0, p1, p0, p1);
                       
    return k * (k-1) / beta / beta / U / U 
        * (k * det_ratio / beta / U + 1)/((k-2) * det_ratio * det_ratio_even / beta / U + det_ratio_odd); // The Metropolis ratio
  }

  dcomplex accept() {
    for (auto &d : config->Mmatrices_even)
      d.complete_operation(); // Finish insertion
    for (auto &d : config->Mmatrices_odd)
      d.complete_operation(); // Finish insertion
    return 1.0;
  }

  void reject() {
    for (auto &d : config->Mmatrices_even)
      d.reject_last_try(); // Finish insertion
    for (auto &d : config->Mmatrices_odd)
      d.reject_last_try(); // Finish insertion  
  }                        // Nothing to do
};

//  -------------- QMC measurement ----------------

struct measure_M {

  configuration *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mw;  // reference to M-matrix
  double beta;
  double U;
  dcomplex Z = 0;
  long count = 0;

  measure_M(configuration *config_, block_gf<imfreq> &Mw_, double beta_, double U_) : config(config_), Mw(Mw_), beta(beta_), U(U_) { Mw() = 0; }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    auto k = config->perturbation_order();
    auto det_ratio = config->Mmatrices_even[up].determinant() * config->Mmatrices_even[down].determinant() 
                    / config->Mmatrices_odd[up].determinant() *  config->Mmatrices_odd[down].determinant();

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda_even = [this, spin, sign, k, det_ratio](arg_t const &x, arg_t const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(2 * phase_step);
        for (auto const &om : mesh) {
          this->Mw[spin][om](0, 0) += sign * M * coeff * coeff * k * det_ratio / this->beta / this->U;
          coeff *= fact;
        }
      };

      auto lambda_odd = [this, spin, sign](arg_t const &x, arg_t const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(2 * phase_step);
        for (auto const &om : mesh) {
          this->Mw[spin][om](0, 0) -=  sign * M * coeff;
          coeff *= fact;
        }
      };

      foreach (config->Mmatrices_even[spin], lambda_even)
        ;

      foreach (config->Mmatrices_odd[spin], lambda_odd)
        ;
      
      this->Mw[spin] /= k * det_ratio / this->beta / this->U - 1;
    }
  }

  void collect_results(mpi::communicator const &c) {
    Mw = mpi::all_reduce(Mw, c);
    Z  = mpi::all_reduce(Z, c);
    Mw = Mw / (-Z * beta);

    // Print the sign
    if (c.rank() == 0) std::cerr << "Average sign " << Z / c.size() / count << std::endl;
  }
};

// ------------ The main class of the solver ------------------------

solver2::solver2(double beta_, int n_iw, int n_tau)
   : beta{beta_},
     g0_iw{make_block_gf({"up", "down"}, gf<imfreq>{{beta, Fermion, n_iw}, {1, 1}})},
     g0tilde_iw{g0_iw},
     g_iw{g0_iw},
     M_iw{g0_iw},
     g0tilde_tau{make_block_gf({"up", "down"}, gf<imtime>{{beta, Fermion, n_tau}, {1, 1}})} {
      std::cout << "--------- /!\\ Using Solver2 /!\\ ---------\n";
     }

// The method that runs the qmc
void solver2::solve(double U, double delta, int n_cycles, int length_cycle, int n_warmup_cycles, std::string random_name, int max_time) {

  mpi::communicator world;
  triqs::clef::placeholder<0> spin_;
  triqs::clef::placeholder<1> om_;

  for (auto spin : {up, down}) { // Apply shift to g0_iw and Fourier transform
    g0tilde_iw[spin](om_) << 1.0 / (1.0 / g0_iw[spin](om_) - U / 2);
    array<dcomplex, 3> mom{{{0}}, {{1}}}; // Fix the moments: 0 + 1/omega
    g0tilde_tau()[spin] = triqs::gfs::fourier(g0tilde_iw[spin], make_const_view(mom));
  }

  // Rank-specific variables
  int verbosity   = (world.rank() == 0 ? 3 : 0);
  int random_seed = 34788 + 928374 * world.rank();

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<dcomplex> CTQMC(random_name, random_seed, verbosity);

  // Prepare the configuration
  auto config = configuration{g0tilde_tau, beta, delta};

  // Register moves and measurements
  CTQMC.add_move(move_insert{&config, CTQMC.get_rng(), beta, U}, "insertion");
  CTQMC.add_move(move_remove{&config, CTQMC.get_rng(), beta, U}, "removal");
  CTQMC.add_measure(measure_M{&config, M_iw, beta, U}, "M measurement");

  // Run and collect results
  CTQMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(max_time));
  CTQMC.collect_results(world);

  // Compute the Green function from Mw
  g_iw[spin_](om_) << g0tilde_iw[spin_](om_) + g0tilde_iw[spin_](om_) * M_iw[spin_](om_) * g0tilde_iw[spin_](om_);

  // Set the tail of g_iw to 1/w
  triqs::arrays::array<dcomplex, 3> mom{{{0}}, {{1}}}; // 0 + 1/omega
  for (auto &g : g_iw) replace_by_tail_in_fit_window(g(), make_const_view(mom));
}
