#include "ctint.hpp"

using namespace triqs::gfs; 
using namespace triqs::arrays; 
using namespace ctint_tutorial;

// --------------- The QMC configuration ----------------

// The Monte Carlo configuration
struct configuration2 {
  // M-matrices for up and down up to even order
  std::vector<triqs::det_manip::det_manip<g0bar_tau>> Mmatrices_even;
  // M-matrices for up and down up to odd order
  std::vector<triqs::det_manip::det_manip<g0bar_tau>> Mmatrices_odd;
  dcomplex det_ratio;
  dcomplex d0;

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
      auto lambda = [this, spin](arg_t const &x, arg_t const &y, dcomplex M) {
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

  void print_taus() {
    for (int i=0; i<perturbation_order();i++)
      std::cout << Mmatrices_odd[up].get_x(i).tau << "  ";
    std::cout << "\n";

    for (int i=0; i<perturbation_order()-1;i++) {
      assert (Mmatrices_odd[up].get_x(i).tau == Mmatrices_even[up].get_x(i).tau 
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_even[up].get_y(i).tau 
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_even[down].get_x(i).tau
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_even[down].get_y(i).tau 
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_odd[up].get_y(i).tau 
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_odd[down].get_x(i).tau
           && Mmatrices_odd[up].get_x(i).tau == Mmatrices_odd[down].get_y(i).tau);
      std::cout << Mmatrices_even[up].get_x(i).tau << "  ";
    }
    std::cout << "\n";
  }

  configuration2(block_gf<imtime> &g0tilde_tau, double beta, double delta, double delta0, int nobc) {
    //std::cout << "--------- /!\\ Initializing double config ";
    // Initialize the M-matrices. 100 is the initial matrix size
    det_ratio = 1;

    double tau = 0;
    int s = 1;
    std::vector<g0bar_tau> g;

    for (auto spin : {up, down}) {
      g.emplace_back(g0bar_tau{g0tilde_tau[spin], beta, delta, delta0, spin});
      Mmatrices_even.emplace_back(g[spin], 100);
      Mmatrices_odd.emplace_back(g[spin], 100);
      det_ratio /= Mmatrices_odd[spin].insert_at_end({tau,s},{tau,s});
      Mmatrices_even[spin].set_n_operations_before_check(nobc);
      Mmatrices_odd[spin].set_n_operations_before_check(nobc);
      //std::cout << "nobc: " << Mmatrices_even[spin].get_n_operations_before_check() << " " << Mmatrices_odd[spin].get_n_operations_before_check() << '\n';
    }

    //print_taus();
    //std::cout << "-> Initialized double config /!\\ ---------\n";
  }
};

// ------------ QMC move : inserting a vertex ------------------


struct move_insert2 {
  configuration2 *config;
  triqs::mc_tools::random_generator &rng;
  double beta, U;
  double tau0, tau1;
  
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
    auto t = config->Mmatrices_odd[up].get_x(k-1);
    arg_t t0 {rng(beta),rng(2)};
    arg_t t1 {rng(beta),rng(2)};
    tau0 = t0.tau;
    tau1 = t1.tau;
    int p0 = 1+rng(k+1);   // Choose one of the operators for insertion
    int p1 = rng(p0);      // Choose one of the operators for insertion
    arg_t t2,t3;
    int p2,p3;
    if (p0<=k) {
      p2 = p0;
      t2 = t0;
      p3 = p1;
      t3 = t1;
    }
    else if (p1<k) {
      p2 = k;
      t2 = t;
      p3 = p1;
      t3 = t1;
    }
    else {
      p2 = p1;
      t2 = t1;
      p3 = k-1;
      t3 = t;
    }

    config->set_det_ratio();
    //std::cout << "insert attempt even";
    auto det_ratio_even = config->Mmatrices_even[up].try_insert2(p2, p3, p2, p3, t2, t3, t2, t3)
                       *config->Mmatrices_even[down].try_insert2(p2, p3, p2, p3, t2, t3, t2, t3);
    //std::cout << "\t -> done \n";

    //std::cout << "insert attempt odd";
    auto det_ratio_odd = config->Mmatrices_odd[up].try_insert2(p0, p1, p0, p1, t0, t1, t0, t1) 
                      *config->Mmatrices_odd[down].try_insert2(p0, p1, p0, p1, t0, t1, t0, t1);
    //std::cout << "\t -> done \n";
/*  std::cout << "-> Used double insert attempt /!\\ ---------\n";
    std::cout << "k: " << k << "\t k-1: " << config->Mmatrices_even[up].size() << '\n';
    std::cout << "det ratio: " << det_ratio << "\t det_ratio_even: " << det_ratio_even << "\t det_ratio_odd: " << det_ratio_odd <<'\n';
*/  return beta * beta * U * U / ((k + 2) * (k + 1))
          * ((k+2) * config->det_ratio * det_ratio_even / ( beta * U)  - det_ratio_odd )
          / ( k   * config->det_ratio / ( beta * U) - 1); // The Metropolis ratio
  //}
  }

  dcomplex accept() {
    //auto k = config->perturbation_order();
    //std::cout << "--------- /!\\ Using double insert accept\n";
    //std::cout << "insert " << tau0 << "  " << tau1 << " in" << '\n';
    //config->print_taus();

    for (auto spin : {up, down}) {

      //std::cout << "insert accept odd, spin: " << spin << " k: " << k <<'\n';
      //config->print();
      config->Mmatrices_odd[spin].complete_operation();
      //std::cout << "\t -> done \n"; // Finish insertion

      //std::cout << "insert accept even";
      config->Mmatrices_even[spin].complete_operation();
      //std::cout << "\t -> done \n"; // Finish insertion
    }
    //config->print_taus();
    //std::cout << "-> Used double insert accept /!\\ ---------\n";
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
  double tau0, tau1;

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
    tau0 = config->Mmatrices_odd[up].get_x(p0).tau;
    tau1 = config->Mmatrices_odd[up].get_x(p1).tau;
    assert (tau0!=tau1);
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
          * ((k-2) * config->det_ratio * det_ratio_even / ( beta * U) - det_ratio_odd )
          / ( k   * config->det_ratio / ( beta * U) - 1); // The Metropolis ratio
  }

  dcomplex accept() {
    //std::cout << "--------- /!\\ Using double remove accept ";
    //std::cout << "remove: " << tau0 << "  " << tau1 << '\n';
      for (auto spin : {up, down}) {
        //std::cout << "remove accept even";
        config->Mmatrices_even[spin].complete_operation();
        //std::cout << "\t -> done \n"; // Finish insertion

        //std::cout << "remove accept odd";
        config->Mmatrices_odd[spin].complete_operation();
        //std::cout << "\t -> done \n"; // Finish insertion
      }
    //std::cout << "-> Used double remove accept /!\\ ---------\n";
    //config->print_taus();
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

struct measure_histogram2 {

  // The Monte-Carlo configuration
  configuration2 const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<double> &histogram;

  // Accumulation counter
  long N;

  measure_histogram2(configuration2 const *config_, std::vector<double> &histogram_)
      : config(config_), histogram(histogram_) {
        histogram = std::vector<double>(2);
        for (auto &h_k : histogram) h_k = 0;
        N=0;
      }

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

struct measure_histogram_sign2 {

  // The Monte-Carlo configuration
  configuration2 const *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_sign;
  
  double beta, U;

  measure_histogram_sign2(configuration2 const *config_, std::vector<dcomplex> &histogram_sign_, double beta_, double U_)
      : config(config_), histogram_sign(histogram_sign_), beta(beta_), U(U_) {
        histogram_sign = std::vector<dcomplex>(2);
        for (auto &h_k : histogram_sign) h_k=0;
      }

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    int k = config->perturbation_order();
    while (k >= histogram_sign.size()) histogram_sign.resize(2 * histogram_sign.size());
    histogram_sign[k] += sign;
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
  
    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_sign.size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_sign.resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_sign = mpi::all_reduce(histogram_sign, comm);

    //for (auto &h_k : histogram_sign) h_k = h_k *(1 - U * beta * config->d0) / histogram_sign[1];
  }
};

struct measure_n2 {

  configuration2 *config; // Pointer to the MC configuration
  std::vector<dcomplex> &n;        // reference to M-matrix
  double beta, U;
  dcomplex Z;

  measure_n2(configuration2 *config_, std::vector<dcomplex> &n_, double beta_, double U_) : config(config_), n(n_), beta(beta_), U(U_) {
    for (auto &n_s : n) n_s = 0;
    Z = 0; 
  }

  void accumulate(dcomplex sign) {
    Z += sign;
    int k = config->perturbation_order();
    config->set_det_ratio();

    arg_t t {0,0};
    std::vector<dcomplex> m{std::vector<dcomplex>(2)};
    for (auto s : {0,1}) {
      t.s = s;
      for (auto spin : {up,down}) {
        m[spin] += config->Mmatrices_even[spin].try_insert(k-1,k-1,t,t)*k* config->det_ratio / (beta * U);
        config->Mmatrices_even[spin].reject_last_try();
      }

      for (auto spin : {up,down}) {
        m[spin] -= config->Mmatrices_odd[spin].try_insert(k,k,t,t);
        config->Mmatrices_odd[spin].reject_last_try();
      }
    }

    for (auto &m_s : m) m_s *= sign / 2 / (k * config->det_ratio / (beta * U) - 1);

    for (auto spin : {up,down}) n[spin] += m[spin];
  }

  void collect_results(mpi::communicator const &c) {
    n = mpi::all_reduce(n, c);
    Z  = mpi::all_reduce(Z, c);
    for (auto &n_s : n) n_s /= Z;
  }

};

struct measure_histogram_n2 {

  // The Monte-Carlo configuration
  configuration2 *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_n;
  double beta, U;
  dcomplex Z;


  measure_histogram_n2(configuration2 *config_, std::vector<dcomplex> &histogram_n_, double beta_, double U_)
      : config(config_), histogram_n(histogram_n_), beta(beta_), U(U_) {
        histogram_n = std::vector<dcomplex>(2);
        for (auto &h_k : histogram_n) h_k = 0;
        Z=0;}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    //std::cout << "accumulating  ";
    Z += sign;
    int k = config->perturbation_order();
    while (k >= histogram_n.size()) histogram_n.resize(2 * histogram_n.size());
    config->set_det_ratio();

    arg_t t {0,0};
    dcomplex n = 0;
    for (auto s : {0,1}) {
      t.s = s;
      for (auto spin : {up,down}) {
        n += config->Mmatrices_even[spin].try_insert(k-1,k-1,t,t)*k* config->det_ratio / (beta * U);
        config->Mmatrices_even[spin].reject_last_try();
      }

      for (auto spin : {up,down}) {
        n -= config->Mmatrices_odd[spin].try_insert(k,k,t,t);
        config->Mmatrices_odd[spin].reject_last_try();
      }
    }

    n *= sign / 2 / (k * config->det_ratio / (beta * U) - 1);

    histogram_n[k] += n;
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    //std::cout << "collecting  ";
  
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

struct measure_d2 {

  configuration2 *config; // Pointer to the MC configuration
  dcomplex &d;        // reference to M-matrix
  double beta,U;
  dcomplex Z;

  measure_d2(configuration2 *config_, dcomplex &d_, double beta_, double U_) : config(config_), d(d_), beta(beta_), U(U_) { d = 0; Z=0;}

  void accumulate(dcomplex sign) {
    Z += sign;
    dcomplex B_even;
    dcomplex B_odd;

    int k = config->perturbation_order();
    config->set_det_ratio();

    arg_t t {0.,0};

    for (auto s : {0,1}) {
      t.s = s;
      B_even = 1.;
      for (auto &m : config->Mmatrices_even) {
        B_even *= m.try_insert(k-1,k-1,t,t);
        m.reject_last_try();
      }

      B_odd = 1.;
      for (auto &m : config->Mmatrices_odd) {
        B_odd *= m.try_insert(k,k,t,t);
        m.reject_last_try();
      }

    d += sign*(B_even * k * config->det_ratio / (beta * U) - B_odd)
        / 2 / (k * config->det_ratio / (beta * U) - 1);
        
    }

  }

  void collect_results(mpi::communicator const &c) {
    d = mpi::all_reduce(d, c);
    Z  = mpi::all_reduce(Z, c);
    d = d / Z;
  }

};

struct measure_histogram_d2 {

  // The Monte-Carlo configuration
  configuration2 *config; // Pointer to the MC configuration

  // Reference to accumulation vector
  std::vector<dcomplex> &histogram_d;  
  double beta, U;
  dcomplex Z;


  measure_histogram_d2(configuration2 *config_, std::vector<dcomplex> &histogram_d_, double beta_, double U_)
      : config(config_), histogram_d(histogram_d_), beta(beta_), U(U_) {
        histogram_d = std::vector<dcomplex>(2);
        for (auto &h_k : histogram_d) h_k = 0;
        Z=0;}

  /// Accumulate perturbation order into histogram
  void accumulate(dcomplex sign) {
    Z += sign;
    dcomplex B_even;
    dcomplex B_odd;

    int k = config->perturbation_order();
    while (k >= histogram_d.size()) histogram_d.resize(2 * histogram_d.size());
    config->set_det_ratio();

    arg_t t {0.,0};
    for (auto s : {0,1}) {
      t.s = s;
      B_even = 1.;
      for (auto &m : config->Mmatrices_even) {
        B_even *= m.try_insert(k-1,k-1,t,t);
        m.reject_last_try();
      }

      B_odd = 1.;
      for (auto &m : config->Mmatrices_odd) {
        B_odd *= m.try_insert(k,k,t,t);
        m.reject_last_try();
      }

      histogram_d[k] += sign*(B_even * k * config->det_ratio / (beta * U) - B_odd)
        / 2 / (k * config->det_ratio / (beta * U) - 1);
    }
  }

  /// Reduce and normalize
  void collect_results(mpi::communicator const &comm) {
    //std::cout << "collecting  ";
  
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

struct measure_M2 {

  configuration2 *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mw;  // reference to M-matrix
  double beta;
  double U;
  dcomplex Z;
  long count;
  block_gf<imfreq> Kw;  // reference to M-matrix


  measure_M2(configuration2 *config_, block_gf<imfreq> &Mw_, double beta_, double U_) : config(config_), Mw(Mw_), beta(beta_), U(U_), count(0), Kw{Mw_} { Mw() = 0; Z=0; count=0;}

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    auto k = config->perturbation_order();
    config->set_det_ratio();
    Kw()=0;

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda_even = [this, spin, k](arg_t const &x, arg_t const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(2 * phase_step);
        for (auto const &om : mesh) {
          this->Kw[spin][om] += M * coeff * k * config->det_ratio / (this->beta * this->U);
          coeff *= fact;
        }
      };

      auto lambda_odd = [this, spin](arg_t const &x, arg_t const &y, dcomplex M) {
        auto const &mesh = this->Mw[spin].mesh();
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
        auto fact        = std::exp(2 * phase_step);
        for (auto const &om : mesh) {
          this->Kw[spin][om] -= M * coeff;
          coeff *= fact;
        }
      };

      foreach (config->Mmatrices_even[spin], lambda_even);

      foreach (config->Mmatrices_odd[spin], lambda_odd);
      
      Kw[spin] *= sign/(k * config->det_ratio / (beta * U) - 1);

      Mw[spin] += Kw[spin];
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

struct measure_Mk2 {

  configuration2 *config; // Pointer to the MC configuration
  block_gf<imfreq> &Mkw;  // reference to M-matrix
  int k;
  double beta;
  double U;
  dcomplex Z;
  long count;
  block_gf<imfreq> Kw;  // reference to M-matrix


  measure_Mk2(configuration2 *config_, block_gf<imfreq> &Mkw_, int k_, double beta_, double U_) : config(config_), Mkw(Mkw_), k(k_), beta(beta_), U(U_), Kw{Mkw_} {
     Mkw() = 0; Z=0; count=0;
     }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    if (k==config->perturbation_order()) {
      config->set_det_ratio();
      Kw()=0;

      for (auto spin : {up, down}) {

        // A lambda to measure the M-matrix in frequency
        auto lambda_even = [this, spin, sign](arg_t const &x, arg_t const &y, dcomplex M) {
          auto const &mesh = this->Mkw[spin].mesh();
          auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
          auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
          auto fact        = std::exp(2 * phase_step);
          for (auto const &om : mesh) {
            //std::cout << om;
            this->Kw[spin][om] += M * coeff * this->k * config->det_ratio / (this->beta * this->U);
            coeff *= fact;
          }
        };

        auto lambda_odd = [this, spin, sign](arg_t const &x, arg_t const &y, dcomplex M) {
          auto const &mesh = this->Mkw[spin].mesh();
          auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
          auto coeff       = std::exp((2 * mesh.first_index() + 1) * phase_step);
          auto fact        = std::exp(2 * phase_step);
          for (auto const &om : mesh) {
            this->Kw[spin][om] -= M * coeff;
            coeff *= fact;
          }
        };

        foreach (config->Mmatrices_even[spin], lambda_even);

        foreach (config->Mmatrices_odd[spin], lambda_odd);
        
        Kw[spin] *= sign/(k * config->det_ratio / (beta * U) - 1);

        Mkw[spin] += Kw[spin];
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

  configuration2 const *config; // Pointer to the MC configuration
  std::vector<dcomplex> &histogram_m;
  int n;
  double beta,U;
  dcomplex Z;
  long count;

  measure_M_wn(configuration2 const *config_, std::vector<dcomplex> &histogram_m_, int n_, double beta_, double U_) : config(config_), histogram_m(histogram_m_), n(n_), beta(beta_), U(U_) {
    histogram_m = std::vector<dcomplex>(2); Z=0; count=0;
    }

  void accumulate(dcomplex sign) {
    Z += sign;
    count++;
    int k = config->perturbation_order();
    while (k >= histogram_m.size()) histogram_m.resize(2 * histogram_m.size());

    config->set_det_ratio();
    dcomplex K = 0;

    for (auto spin : {up, down}) {

      // A lambda to measure the M-matrix in frequency
      auto lambda_even = [this, spin, k, K](arg_t const &x, arg_t const &y, dcomplex M) {
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * this->n + 1) * phase_step);
        K += M * coeff * k * config->det_ratio / (this->beta * this->U);
      };

      auto lambda_odd = [this, spin, K](arg_t const &x, arg_t const &y, dcomplex M) {
        auto phase_step  = -1.0i * M_PI * (x.tau - y.tau) / beta;
        auto coeff       = std::exp((2 * this->n + 1) * phase_step);
        K -= M * coeff;
      };

      foreach (config->Mmatrices_even[spin], lambda_even);

      foreach (config->Mmatrices_odd[spin], lambda_odd);
      
      K *= sign/(k * config->det_ratio / (beta * U) - 1);

      histogram_m[k] += K;
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

solver2::solver2(double beta_, int n_iw, int n_tau)
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
      std::cout << "--------- CTINT2-TUTORIAL ---------\n";
     }

// The method that runs the qmc
void solver2::solve(double U, double delta, double delta0, int k, int nobc, int n_cycles, int length_cycle, int n_warmup_cycles, std::string random_name, int max_time, int seed) {
  std::cout << "--------- CTINT2-TUTORIAL ---------\n";

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
  auto config = configuration2{g0tilde_tau, beta, delta, delta0, nobc};

  // Register moves and measurements
  CTQMC.add_move(move_insert2{&config, CTQMC.get_rng(), beta, U,0,0}, "insertion");
  CTQMC.add_move(move_remove2{&config, CTQMC.get_rng(), beta, U,0,0}, "removal");
  CTQMC.add_measure(measure_histogram2{&config, hist}, "histogram measurement");
  CTQMC.add_measure(measure_histogram_sign2{&config, hist_sign, beta, U}, "sign histogram measurement");
  //CTQMC.add_measure(measure_n2{&config, n, beta, U}, "density measurement");
  //CTQMC.add_measure(measure_histogram_n2{&config, hist_n, beta, U}, "density histogram measurement");
  //CTQMC.add_measure(measure_d2{&config, d, beta, U}, "double occupancy measurement");
  //CTQMC.add_measure(measure_histogram_d2{&config, hist_d, beta, U}, "double occupancy histogram measurement");
  //CTQMC.add_measure(measure_M2{&config, m_iw, beta, U}, "M measurement");
  //if (k>0)
  //  CTQMC.add_measure(measure_Mk2{&config, mk_iw, k, beta, U}, "M kth order measurement");

  // Run and collect results
  CTQMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(max_time));
  CTQMC.collect_results(world);

  for(auto &n_s : n) {
    n_s += delta0;
    d += delta0 * n_s;
  }
  d -= delta0 * delta0 - delta * delta;

  // Compute the Green function from Mw
  g_iw[spin_](om_) << g0tilde_iw[spin_](om_) + g0tilde_iw[spin_](om_) * m_iw[spin_](om_) * g0tilde_iw[spin_](om_);

  // Set the tail of g_iw to 1/w
  triqs::arrays::array<dcomplex, 3> mom{{{0}}, {{1}}}; // 0 + 1/omega
  for (auto &g : g_iw) replace_by_tail_in_fit_window(g(), make_const_view(mom));
}
