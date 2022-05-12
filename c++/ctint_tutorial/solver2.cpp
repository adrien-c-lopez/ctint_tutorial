#include "ctint.hpp"

using namespace triqs::gfs; 
using namespace triqs::arrays; 
using namespace ctint_tutorial;

// --------------- The QMC configuration ----------------

// Argument type of g0bar
struct arg_t2 {
  double tau; // The imaginary time
  int s;      // The auxiliary spin
};

// The function that appears in the calculation of the determinant
struct g0bar_tau2 {
  gf<imtime> const &gt;
  double beta, delta;
  int s;

  dcomplex operator()(arg_t2 const &x, arg_t2 const &y) const {
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
struct configuration2 {
  // M-matrices for up and down
  std::vector<triqs::det_manip::det_manip<g0bar_tau2>> Mmatrices_even;
  std::vector<triqs::det_manip::det_manip<g0bar_tau2>> Mmatrices_odd;
  dcomplex det_ratio;

  int perturbation_order() const {
    assert (Mmatrices_odd[up].size()==Mmatrices_odd[down].size());
    assert (Mmatrices_odd[up].size()==Mmatrices_even[down].size()+1); 
    assert (Mmatrices_odd[up].size()==Mmatrices_even[up].size()+1);
    return Mmatrices_odd[up].size(); } //odd perturbation order

  void set_det_ratio() {
    //std::cout << "setting det ratio ->\n";
    det_ratio = 1;
    for (auto &d : Mmatrices_odd) {
      //std::cout << "set_det_ratio, remove";
      det_ratio *= d.try_remove(perturbation_order()-1,perturbation_order()-1);
      //std::cout << "\t -> done \n";

      d.reject_last_try();
    }
    //std::cout << "k: " << perturbation_order() << " det_ratio: " << det_ratio << '\n';
    //std::cout << "-> det ratio set\n";
  }

  void print() {
    for (auto spin : {up, down}) {
      auto lambda = [this, spin](arg_t2 const &x, arg_t2 const &y, dcomplex M) {
        std::cout << M << '\t';
      };
      std::cout << "spin: " << spin << "\teven\n";
      foreach (Mmatrices_even[spin], lambda);
      std::cout <<'\n';
      std::cout << "odd\n";
      foreach (Mmatrices_odd[spin], lambda);
      std::cout <<'\n';
    }
  }

  configuration2(block_gf<imtime> &g0tilde_tau, double beta, double delta) {
    //std::cout << "--------- /!\\ Initializing double config ";
    // Initialize the M-matrices. 100 is the initial matrix size
    det_ratio = 1;
    double tau = 0;
    int s = 1;
    for (auto spin : {up, down}) {
      Mmatrices_even.emplace_back(g0bar_tau2{g0tilde_tau[spin], beta, delta, spin}, 100);
      Mmatrices_odd.emplace_back(g0bar_tau2{g0tilde_tau[spin], beta, delta, spin}, 100);
      det_ratio /= Mmatrices_odd[spin].insert_at_end({tau,s},{tau,s});
      //Mmatrices_even[spin].set_n_operations_before_check(50);
      //Mmatrices_odd[spin].set_n_operations_before_check(50);
      Mmatrices_even[spin].set_n_operations_before_check(1000000000);
      Mmatrices_odd[spin].set_n_operations_before_check(1000000000);

    }
    //std::cout << "-> Initialized double config /!\\ ---------\n";
  }
};

// ------------ QMC move : inserting a vertex ------------------


struct move_insert2 {
  configuration2 *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;
  
  dcomplex attempt() { // Insert an interaction vertex at time tau with aux spin s
    //std::cout << "--------- /!\\ Using double insert attempt \n";
    auto k = config->perturbation_order();
/*    if (k >= 5) {
      std::cout << "--------- /!\\ Using double insert attempt \n";
      std::cout << "odd , up tau: " << config->Mmatrices_odd[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-2).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-1).tau << '\n';
      std::cout << "\teven, up tau: " << config->Mmatrices_even[up].get_x(k-5).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-2).tau << '\n';
    }
*/
    auto t         = config->Mmatrices_odd[up].get_x(k-1);
    double tau     = t.tau;
    double tau0    = rng(beta);
    double tau1    = rng(beta);
    int s          = t.s;
    int s0         = rng(2);
    int s1         = rng(2);

    config->set_det_ratio();
    //std::cout << "insert attempt even";
    auto det_ratio_even = config->Mmatrices_even[up].try_insert2(k-1, k, k-1, k, {tau, s}, {tau0, s0}, {tau, s}, {tau0, s0}) 
                       *config->Mmatrices_even[down].try_insert2(k-1, k, k-1, k, {tau, s}, {tau0, s0}, {tau, s}, {tau0, s0});
    //std::cout << "\t -> done \n";

    //std::cout << "insert attempt odd";
    auto det_ratio_odd = config->Mmatrices_odd[up].try_insert2(k, k+1, k, k+1, {tau0, s0}, {tau1, s1}, {tau0, s0}, {tau1, s1}) 
                      *config->Mmatrices_odd[down].try_insert2(k, k+1, k, k+1, {tau0, s0}, {tau1, s1}, {tau0, s0}, {tau1, s1});
    //std::cout << "\t -> done \n";
/*  std::cout << "-> Used double insert attempt /!\\ ---------\n";
    std::cout << "k: " << k << "\t k-1: " << config->Mmatrices_even[up].size() << '\n';
    std::cout << "det ratio: " << det_ratio << "\t det_ratio_even: " << det_ratio_even << "\t det_ratio_odd: " << det_ratio_odd <<'\n';
*/  return beta * beta * U * U / ((k + 2) * (k + 1))
          * ((k+2) * config->det_ratio * det_ratio_even / (beta * U)  - det_ratio_odd )
          / (  k   * config->det_ratio / (beta * U) - 1); // The Metropolis ratio
  //}
  }

  dcomplex accept() {
    auto k = config->perturbation_order();
    //std::cout << "--------- /!\\ Using double insert accept\n";
    for (auto spin : {up, down}) {

      //std::cout << "insert accept odd, spin: " << spin << " k: " << k <<'\n';
      //config->print();
      config->Mmatrices_odd[spin].complete_operation();
      //std::cout << "\t -> done \n"; // Finish insertion

      //std::cout << "insert accept even";
      config->Mmatrices_even[spin].complete_operation();
      //std::cout << "\t -> done \n"; // Finish insertion
    }
/*    auto k = config->perturbation_order();
    if (k >= 5) {
      std::cout << "odd , up tau: " << config->Mmatrices_odd[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-2).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-1).tau << '\n';
      std::cout << "even, up tau: " << config->Mmatrices_even[up].get_x(k-5).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-2).tau << '\n';
    }
*/    //std::cout << "-> Used double insert accept /!\\ ---------\n";
    return 1.0;
  }

  void reject() {
    //std::cout << "--------- /!\\ Using double insert reject ";
    for (auto spin : {up, down}) {
      config->Mmatrices_even[spin].reject_last_try(); // Finish insertion
      config->Mmatrices_odd[spin].reject_last_try(); // Finish insertion
    }
    //std::cout << "-> Used double insert reject /!\\ ---------\n";
  }
};

// ------------ QMC move : deleting a vertex ------------------

struct move_remove2 {
  configuration2 *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;

  dcomplex attempt() {
    //std::cout << "--------- /!\\ Using double remove attempt \n";
    auto k = config->perturbation_order();
    if (k < 3) return 0;    // Config is empty, trying to remove makes no sense

/*    if (k >= 5) {
      std::cout << "odd , up tau: " << config->Mmatrices_odd[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-2).tau << '\t';
      std::cout << config->Mmatrices_odd[up].get_x(k-1).tau << '\n';
      std::cout << "\teven, up tau: " << config->Mmatrices_even[up].get_x(k-5).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-4).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-3).tau << '\t';
      std::cout << config->Mmatrices_even[up].get_x(k-2).tau << '\n';
    }
*/
    int p0 = 1+rng(k-1); // Choose one of the operators for removal
    int p1 = rng(p0);    // Choose one of the operators for removal
    int p2;
    if (p0 == k-1) {
      if (p1 == k-2)
        p2 = k-3;
      else
        p2 = k-2;
    }
    else
      p2 = p0;

    config->set_det_ratio();
    
    //std::cout << "remove attempt even";
    auto det_ratio_even = config->Mmatrices_even[up].try_remove2(p2, p1, p2, p1) 
                       *config->Mmatrices_even[down].try_remove2(p2, p1, p2, p1);
    //std::cout << "\t -> done \n";

    //std::cout << "remove attempt odd";
    auto det_ratio_odd = config->Mmatrices_odd[up].try_remove2(p0, p1, p0, p1)  
                      *config->Mmatrices_odd[down].try_remove2(p0, p1, p0, p1);
    //std::cout << "\t -> done \n";
/*    std::cout << "-> Used double remove attempt /!\\ ---------\n";
    std::cout << "k: " << k << "\t k-1: " << config->Mmatrices_even[up].size() << '\n';
    std::cout << "det ratio: " << det_ratio << "\t det_ratio_even: " << det_ratio_even << "\t det_ratio_odd: " << det_ratio_odd <<'\n';                           
*/    return k * (k-1) / (beta * U * beta * U) 
          * ((k-2) * config->det_ratio * det_ratio_even / (beta * U) - det_ratio_odd )
          / (  k   * config->det_ratio / (beta * U) - 1); // The Metropolis ratio
  }

  dcomplex accept() {
    //std::cout << "--------- /!\\ Using double remove accept ";
      for (auto spin : {up, down}) {
        //std::cout << "remove accept even";
        config->Mmatrices_even[spin].complete_operation();
        //std::cout << "\t -> done \n"; // Finish insertion

        //std::cout << "remove accept odd";
        config->Mmatrices_odd[spin].complete_operation();
        //std::cout << "\t -> done \n"; // Finish insertion
      }
    //std::cout << "-> Used double remove accept /!\\ ---------\n";
    return 1.0;
  }

  void reject() {
    //std::cout << "--------- /!\\ Using double remove reject ";
    for (auto spin : {up, down}) {
      config->Mmatrices_even[spin].reject_last_try(); // Finish insertion
      config->Mmatrices_odd[spin].reject_last_try(); // Finish insertion
    }
    //std::cout << "-> Used double remove reject /!\\ ---------\n"; 
  }                        // Nothing to do
};

//  -------------- QMC measurement ----------------

struct measure_M2 {

  configuration2 *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mw;  // reference to M-matrix
  double beta;
  double U;
  dcomplex Z = 0;
  long count = 0;

  measure_M2(configuration2 *config_, block_gf<imfreq> &Mw_, double beta_, double U_) : config(config_), Mw(Mw_), beta(beta_), U(U_) { Mw() = 0; }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    auto k = config->perturbation_order();
    config->set_det_ratio();

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda_even = [this, spin, sign, k](arg_t2 const &x, arg_t2 const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -2.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(phase_step);
        for (auto const &om : mesh) {
          this->Mw[spin][om](0, 0) += sign * M * coeff * k * config->det_ratio / (this->beta * this->U);
          coeff *= fact;
        }
      };

      auto lambda_odd = [this, spin, sign](arg_t2 const &x, arg_t2 const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -2.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp(( mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp( phase_step);
        for (auto const &om : mesh) {
          this->Mw[spin][om](0, 0) -=  sign * M * coeff;
          coeff *= fact;
        }
      };

      foreach (config->Mmatrices_even[spin], lambda_even);

      foreach (config->Mmatrices_odd[spin], lambda_odd);
      
      this->Mw[spin] /= std::abs((k * config->det_ratio / (this->beta * this->U - 1)));
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

struct measure_histogram_sign2 {

  // The Monte-Carlo configuration
  configuration2 const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_sign;

  // Accumulation counter
  long N = 0;

  measure_histogram_sign2(configuration2 const *config_, std::vector<dcomplex> &histogram_sign_)
      : config(config_), histogram_sign(histogram_sign_) {}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    int k = config->perturbation_order();
    while (k >= histogram_sign.size()) histogram_sign.resize(2 * histogram_sign.size());
    histogram_sign[k] += sign;
    ++N;
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_sign.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_sign.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_sign = mpi::all_reduce(histogram_sign, comm);

    for (auto &h_k : histogram_sign) h_k = h_k / N;
  }
};

struct measure_histogram2 {

  // The Monte-Carlo configuration
  configuration2 const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<double> &histogram;

  // Accumulation counter
  long N = 0;

  measure_histogram2(configuration2 const *config_, std::vector<double> &histogram_)
      : config(config_), histogram(histogram_) {}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    int k = config->perturbation_order();
    while (k >= histogram.size()) histogram.resize(2 * histogram.size());
    histogram[k] += 1.;
    ++N;
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
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

// ------------ The main class of the solver ------------------------

solver2::solver2(double beta_, int n_iw, int n_tau)
   : beta{beta_},
     g0_iw{make_block_gf({"up", "down"}, gf<imfreq>{{beta, Fermion, n_iw}, {1, 1}})},
     g0tilde_iw{g0_iw},
     g_iw{g0_iw},
     M_iw{g0_iw},
     g0tilde_tau{make_block_gf({"up", "down"}, gf<imtime>{{beta, Fermion, n_tau}, {1, 1}})}, 
     hist{std::vector<double>(2)},
     hist_sign{std::vector<dcomplex>(2)} {
      std::cout << "--------- /!\\ Using Solver2 /!\\ ---------\n";
     }

// The method that runs the qmc
void solver2::solve(double U, double delta, int n_cycles, int length_cycle, int n_warmup_cycles, std::string random_name, int max_time, int seed) {

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
  int random_seed = seed + 928374 * world.rank();

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<dcomplex> CTQMC(random_name, random_seed, verbosity);

  // Prepare the configuration
  auto config = configuration2{g0tilde_tau, beta, delta};

  // Register moves and measurements
  CTQMC.add_move(move_insert2{&config, CTQMC.get_rng(), beta, U}, "insertion");
  CTQMC.add_move(move_remove2{&config, CTQMC.get_rng(), beta, U}, "removal");
  CTQMC.add_measure(measure_M2{&config, M_iw, beta, U}, "M measurement");
  CTQMC.add_measure(measure_histogram2{&config, hist}, "histogram measurement");
  CTQMC.add_measure(measure_histogram_sign2{&config, hist_sign}, "sign histogram measurement");

  // Run and collect results
  CTQMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(max_time));
  CTQMC.collect_results(world);

  // Compute the Green function from Mw
  g_iw[spin_](om_) << g0tilde_iw[spin_](om_) + g0tilde_iw[spin_](om_) * M_iw[spin_](om_) * g0tilde_iw[spin_](om_);

  // Set the tail of g_iw to 1/w
  triqs::arrays::array<dcomplex, 3> mom{{{0}}, {{1}}}; // 0 + 1/omega
  for (auto &g : g_iw) replace_by_tail_in_fit_window(g(), make_const_view(mom));
}
