# Generated automatically using the command :
# c++2py ../../c++/ctint_tutorial/solver.hpp -p --members_read_only -N ctint_tutorial -a ctint_tutorial -m ctint_tutorial_module -o ctint_tutorial_module --moduledoc="The ctint_tutorial python module" -C pytriqs --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "ctint_tutorial_module", doc = r"The ctint_tutorial python module", app_name = "ctint_tutorial")

# Imports
module.add_imports(*['triqs.gf'])

# Add here all includes
module.add_include("ctint_tutorial/ctint.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/string.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace ctint_tutorial;
""")

module.add_enum("spin", ['spin::up', 'spin::down'], "ctint_tutorial", doc = r"""""")

# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "ctint_tutorial::solver",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(double beta_, int n_iw = 500, int n_tau = 5001)""", doc = r"""Construct a ctint solver""")

c.add_method("""void solve (double U, double delta=0, double delta0=0, int k=-1, int n_cycles=10000, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = \"\", int max_time = -1, int seed = 34788)""",
             doc = r"""Method that performs the QMC calculation""")

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G0_iw ()"),
               doc = r"""Access non-interacting Matsubara Green function""")

c.add_property(name = "G0_tau",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> G0_tau ()"),
               doc = r"""Access non-interacting imaginary-time Green function""")

c.add_property(name = "G_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")

c.add_property(name = "Hist",
               getter = cfunction("std::vector<double> Hist ()"),
               doc = r"""Access order histogram""")

c.add_property(name = "Hist_sign",
                getter = cfunction("std::vector<dcomplex> Hist_sign ()"),
                doc = r"""Access order sign histogram""")

c.add_property(name = "N",
               getter = cfunction("std::vector<dcomplex> N ()"),
               doc = r"""Access density""")

c.add_property(name = "Hist_n",
               getter = cfunction("std::vector<dcomplex> Hist_n ()"),
               doc = r"""Access order density histogram""")

c.add_property(name = "D",
                getter = cfunction("dcomplex D ()"),
                doc = r"""Access double occupancy""")

c.add_property(name = "Hist_d",
               getter = cfunction("std::vector<dcomplex> Hist_d ()"),
               doc = r"""Access order double occupancy histogram""")

c.add_property(name = "M_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> M_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")

c.add_property(name = "Mk_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> Mk_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")              
module.add_class(c)

# The class solver2
c2 = class_(
        py_type = "Solver2",  # name of the python class
        c_type = "ctint_tutorial::solver2",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c2.add_constructor("""(double beta_, int n_iw = 500, int n_tau = 5001)""", doc = r"""Construct a ctint solver""")

c2.add_method("""void solve (double U, double delta=0, double delta0=0, int k=-1, int nobc=50, int n_cycles=10000, int length_cycle = 50, int n_warmup_cycles = 5000, std::string random_name = \"\", int max_time = -1, int seed = 34788)""",
             doc = r"""Method that performs the QMC calculation""")

c2.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G0_iw ()"),
               doc = r"""Access non-interacting Matsubara Green function""")

c2.add_property(name = "G0_tau",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> G0_tau ()"),
               doc = r"""Access non-interacting imaginary-time Green function""")

c2.add_property(name = "G_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")

c2.add_property(name = "Hist",
               getter = cfunction("std::vector<double> Hist ()"),
               doc = r"""Access order histogram""")

c2.add_property(name = "Hist_sign",
               getter = cfunction("std::vector<dcomplex> Hist_sign ()"),
               doc = r"""Access order sign histogram""")

c2.add_property(name = "N",
               getter = cfunction("std::vector<dcomplex> N ()"),
               doc = r"""Access density""")

c2.add_property(name = "Hist_n",
               getter = cfunction("std::vector<dcomplex> Hist_n ()"),
               doc = r"""Access order density histogram""")

c2.add_property(name = "D",
                getter = cfunction("dcomplex D ()"),
                doc = r"""Access double occupancy""")

c2.add_property(name = "Hist_d",
               getter = cfunction("std::vector<dcomplex> Hist_d ()"),
               doc = r"""Access order double occupancy histogram""")

c2.add_property(name = "M_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> M_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")

c2.add_property(name = "Mk_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> Mk_iw ()"),
               doc = r"""Access interacting Matsubara Green function""")

module.add_class(c2)

module.generate_code()
