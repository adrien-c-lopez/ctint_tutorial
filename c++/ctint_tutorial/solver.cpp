#include "ctint.hpp"

using namespace triqs::gfs; 
using namespace triqs::arrays; 
using namespace ctint_tutorial;

// --------------- The QMC configuration ----------------

// The Monte Carlo configuration
struct configuration {
  // M-matrices for up and down
  std::vector<triqs::det_manip::det_manip<g0bar_tau>> Mmatrices;

  int perturbation_order() const { return Mmatrices[up].size(); }

  configuration(block_gf<imtime> &g0tilde_tau, double beta, double delta, double delta0) {
    // Initialize the M-matrices. 100 is the initial matrix size
    for (auto spin : {up, down}) Mmatrices.emplace_back(g0bar_tau{g0tilde_tau[spin], beta, delta, delta0, spin}, 100);
  }
};

// ------------ QMC move : inserting a vertex ------------------

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

// ------------ QMC move : deleting a vertex ------------------

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

//  -------------- QMC measurement ----------------
struct measure_histogram {

  // The Monte-Carlo configuration
  configuration const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<double> &histogram;

  // Accumulation counter
  long N;

  measure_histogram(configuration const *config_, std::vector<double> &histogram_)
      : config(config_), histogram(histogram_) {histogram = std::vector<double>(2); N=0;}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    //std::cout << "accumulating  ";
    int k = config->perturbation_order();
    while (k >= histogram.size()) histogram.resize(2 * histogram.size());
    histogram[k] += 1.;
    N += 1;
    //std::cout << "->  accumulated k:" << k << '\n';
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    std::cout << "collecting  ";
    N = mpi::all_reduce(N, comm);
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram = mpi::all_reduce(histogram, comm);
    for (auto &h_k : histogram) h_k = h_k / N;
  }
};

struct measure_histogram_sign {

  // The Monte-Carlo configuration
  configuration const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_sign;


  measure_histogram_sign(configuration const *config_, std::vector<dcomplex> &histogram_sign_)
      : config(config_), histogram_sign(histogram_sign_) {histogram_sign = std::vector<dcomplex>(2);}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    //std::cout << "accumulating  ";
    int k = config->perturbation_order();
    while (k >= histogram_sign.size()) histogram_sign.resize(2 * histogram_sign.size());
    histogram_sign[k] += sign;
    //std::cout << "->  accumulated k:" << k << '\n';
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    std::cout << "collecting  ";
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_sign.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_sign.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_sign = mpi::all_reduce(histogram_sign, comm);

    //for (auto &h_k : histogram_sign) h_k = h_k / histogram_sign[0];
  }
};

struct measure_n {

  configuration *config; // Pointer to the MC configuration
  std::vector<dcomplex> &n;        // reference to M-matrix
  dcomplex Z;

  measure_n(configuration *config_, std::vector<dcomplex> &n_) : config(config_), n(n_) { 
    for (auto &n_s : n) n_s = 0;
    Z=0;
  }


  void accumulate(dcomplex sign) {
    Z += sign;
    int k = config->perturbation_order();
    arg_t t {0,0};
    for (auto s : {0,1}) {
      t.s = s;
      for (auto spin : {up,down}) {
        n[spin] += sign*config->Mmatrices[spin].try_insert(k,k,t,t)/2;
        config->Mmatrices[spin].reject_last_try();
      }
    }
  }

  void collect_results(mpi::communicator const &c) {
    n = mpi::all_reduce(n, c);
    Z  = mpi::all_reduce(Z, c);
    for (auto &n_s : n) n_s = n_s / Z;
  }

};

struct measure_histogram_n {

  // The Monte-Carlo configuration
  configuration *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_n;
  dcomplex Z;


  measure_histogram_n(configuration *config_, std::vector<dcomplex> &histogram_n_)
      : config(config_), histogram_n(histogram_n_) {histogram_n = std::vector<dcomplex>(2); Z=0;}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    //std::cout << "accumulating  ";
    Z += sign;
    int k = config->perturbation_order();
    while (k >= histogram_n.size()) histogram_n.resize(2 * histogram_n.size());

    arg_t t {0,0};
    for (auto s : {0,1}) {
      t.s = s;
      for (auto spin : {up,down}) {
        histogram_n[k] += sign*config->Mmatrices[spin].try_insert(k,k,t,t)/2;
        config->Mmatrices[spin].reject_last_try();
      }
    }
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    std::cout << "collecting  ";
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_n.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_n.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_n = mpi::all_reduce(histogram_n, comm);
    Z  = mpi::all_reduce(Z, comm);

    for (auto &h_k : histogram_n) h_k = h_k / Z;
  }
};

struct measure_d {

  configuration *config; // Pointer to the MC configuration
  dcomplex &d;        // reference to M-matrix
  dcomplex Z;

  measure_d(configuration *config_, dcomplex &d_) : config(config_), d(d_) { d=0; Z=0;}

  void accumulate(dcomplex sign) {
    Z += sign;
    dcomplex B;
    int k = config->perturbation_order();
    arg_t t {0.,0};

    for (auto s : {0,1}) {
      t.s = s;
      B = 1.;
      for (auto &m : config->Mmatrices) {
          B *= m.try_insert(k,k,t,t);
          m.reject_last_try();
      }
      d += sign*B/2;
    }
  }

  void collect_results(mpi::communicator const &c) {
    d = mpi::all_reduce(d, c);
    Z  = mpi::all_reduce(Z, c);
    d = d / Z;
  }

};

struct measure_histogram_d {

  // The Monte-Carlo configuration
  configuration *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_d;
  dcomplex Z;


  measure_histogram_d(configuration *config_, std::vector<dcomplex> &histogram_d_)
      : config(config_), histogram_d(histogram_d_) {histogram_d = std::vector<dcomplex>(2); Z=0;}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    //std::cout << "accumulating  ";
    Z += sign;
    dcomplex B;
    int k = config->perturbation_order();
    while (k >= histogram_d.size()) histogram_d.resize(2 * histogram_d.size());

    arg_t t {0.,0};
    for (auto s : {0,1}) {
      t.s = s;
      B = 1.;
      for (auto &m : config->Mmatrices) {
          B *= m.try_insert(k,k,t,t);
          m.reject_last_try();
      }
      histogram_d[k] += sign*B/2;
      }
    }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    std::cout << "collecting  ";
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_d.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_d.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_d = mpi::all_reduce(histogram_d, comm);
    Z  = mpi::all_reduce(Z, comm);

    for (auto &h_k : histogram_d) h_k = h_k / Z;
  }
};

struct measure_M {

  configuration const *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mw;        // reference to M-matrix
  double beta;
  dcomplex Z = 0;
  long count = 0;

  measure_M(configuration const *config_, block_gf<imfreq> &Mw_, double beta_) : config(config_), Mw(Mw_), beta(beta_) { Mw() = 0; }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda = [this, spin, sign](arg_t const &x, arg_t const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(2 * phase_step);
        for (auto const &om : mesh) {
          this->Mw[spin][om](0, 0) += sign * M * coeff;
          coeff *= fact;
        }
      };

      foreach (config->Mmatrices[spin], lambda)
        ;
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

struct measure_Mk {

  configuration const *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mkw;        // reference to M-matrix
  int k;
  double beta;
  dcomplex Z;
  long count;

  measure_Mk(configuration const *config_, block_gf<imfreq> &Mkw_, int k_, double beta_) : config(config_), Mkw(Mkw_), k(k_), beta(beta_) { Mkw() = 0; Z=0; count=0;}

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;

    if (config->perturbation_order()==k || config->perturbation_order()==k-1) {

      for (auto spin : {up, down}) {

        // A lambda to measure the M-matrix in frequency
        auto lambda = [this, spin, sign](arg_t const &x, arg_t const &y, dcomplex M) {
          auto const &mesh = this->Mkw[spin].mesh();
          auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
          auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
          auto fact        = std::exp(2 * phase_step);
          for (auto const &om : mesh) {
            this->Mkw[spin][om](0, 0) += sign * M * coeff;
            coeff *= fact;
          }
        };

        foreach (config->Mmatrices[spin], lambda);
      }
    }
  }

  void collect_results(mpi::communicator const &c) {
    Mkw = mpi::all_reduce(Mkw, c);
    Z  = mpi::all_reduce(Z, c);
    Mkw = Mkw / (-Z * beta);

    // Print the sign
    if (c.rank() == 0) std::cerr << "Average sign " << Z / c.size() / count << std::endl;
  }
};
/*
struct measure_M_wn {

  configuration const *config; // Pointer to the MC configuration
  std::vector<dcomplex> &histogram_m;
  int n;
  double beta;
  dcomplex Z;
  long count;

  measure_M_wn(configuration const *config_, std::vector<dcomplex> &histogram_m_, int n_, double beta_) : config(config_), histogram_m(histogram_m_), n(n_), beta(beta_) {
    histogram_m = std::vector<dcomplex>(2); Z=0; count=0;
  }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    int k = config->perturbation_order();
    while (k >= histogram_m.size()) histogram_m.resize(2 * histogram_m.size());

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda = [this, spin, sign, k](arg_t const &x, arg_t const &y, dcomplex M) {
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * this->n + 1) * phase_step);
        this->histogram_m[k] += sign * M * coeff;
      };

      foreach (config->Mmatrices[spin], lambda)
        ;
    }
  }

  void collect_results(mpi::communicator const &comm) {
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_m.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_m.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_m = mpi::all_reduce(histogram_m, comm);
    Z  = mpi::all_reduce(Z, comm);

    for (auto &h_k : histogram_m) h_k = h_k / (-Z * beta);
  }
};
*/
// ------------ The main class of the solver ------------------------

solver::solver(double beta_, int n_iw, int n_tau)
   : beta{beta_},
     g0_iw{make_block_gf({"up", "down"}, gf<imfreq>{{beta, Fermion, n_iw}, {1, 1}})},
     g0tilde_iw{g0_iw},
     g_iw{g0_iw},
     m_iw{g0_iw},
     mk_iw{g0_iw},
     g0tilde_tau{make_block_gf({"up", "down"}, gf<imtime>{{beta, Fermion, n_tau}, {1, 1}})},
     hist{std::vector<double>(2)},
     hist_sign{std::vector<dcomplex>(2)},
     n{std::vector<dcomplex>(2)},
     hist_n{std::vector<dcomplex>(2)},
     d{0},
     hist_d{std::vector<dcomplex>(2)}
     //hist_m{std::vector<dcomplex>(2)}
     {
       std::cout << "--------- /!\\ Using Solver /!\\ ---------\n";
     }

// The method that runs the qmc
void solver::solve(double U, double delta, double delta0, int k, int n_cycles, int length_cycle, int n_warmup_cycles, std::string random_name, int max_time, int seed) {
  std::cout << "--------- /!\\ Using Solver /!\\ ---------\n";

  mpi::communicator world;
  triqs::clef::placeholder<0> spin_;
  triqs::clef::placeholder<1> om_;

  for (auto spin : {up, down}) { // Apply shift to g0_iw and Fourier transform
    g0tilde_iw[spin](om_) << 1.0 / (1.0 / g0_iw[spin](om_) - U * delta0);
    array<dcomplex, 3> mom{{{0}}, {{1}}}; // Fix the moments: 0 + 1/omega
    g0tilde_tau()[spin] = triqs::gfs::fourier(g0tilde_iw[spin], make_const_view(mom));
  }

  // Rank-specific variables
  int verbosity   = (world.rank() == 0 ? 3 : 0);
  int random_seed = seed + 928374 * world.rank();

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<dcomplex> CTQMC(random_name, random_seed, verbosity);

  // Prepare the configuration
  auto config = configuration{g0tilde_tau, beta, delta, delta0};

  // Register moves and measurements
  CTQMC.add_move(move_insert{&config, CTQMC.get_rng(), beta, U}, "insertion");
  CTQMC.add_move(move_remove{&config, CTQMC.get_rng(), beta, U}, "removal");
  CTQMC.add_measure(measure_histogram{&config, hist}, "histogram measurement");
  CTQMC.add_measure(measure_histogram_sign{&config, hist_sign}, "sign histogram measurement");
  CTQMC.add_measure(measure_n{&config, n}, "n measurement");
  CTQMC.add_measure(measure_histogram_n{&config, hist_n}, "n histogram measurement");
  CTQMC.add_measure(measure_d{&config, d}, "double occupancy measurement");
  CTQMC.add_measure(measure_histogram_d{&config, hist_d}, "double occupancy histogram measurement");
  CTQMC.add_measure(measure_M{&config, m_iw, beta}, "M measurement");
  if (k > 0)
    CTQMC.add_measure(measure_Mk{&config, mk_iw, k, beta}, "M kth order measurement");


  // Run and collect results
  CTQMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(max_time));
  CTQMC.collect_results(world);

  // Compute the Green function from Mw
  g_iw[spin_](om_) << g0tilde_iw[spin_](om_) + g0tilde_iw[spin_](om_) * m_iw[spin_](om_) * g0tilde_iw[spin_](om_);

  // Set the tail of g_iw to 1/w
  triqs::arrays::array<dcomplex, 3> mom{{{0}}, {{1}}}; // 0 + 1/omega
  for (auto &g : g_iw) replace_by_tail_in_fit_window(g(), make_const_view(mom));
}
