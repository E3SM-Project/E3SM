#include "shoc_main_wrap.hpp"
#include "shoc_test_data.hpp"
#include "shoc_ic_cases.hpp"

#include "share/eamxx_types.hpp"
#include "share/eamxx_session.hpp"
#include "share/util/eamxx_utils.hpp"

#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <vector>

namespace {
using namespace scream;
using namespace scream::shoc;

/* shoc_run_and_cmp can be run in 2 modes. First, generate_baseline
 * runs the baseline (aka reference, probably git master) version of
 * the code and saves its output as a raw binary file. Then run_and_cmp
 * runs the new/experimental version of the code and compares it against
 * the baseline data you've saved to file. Both baseline and cmp modes
 * start from an initial condition in ../shoc_ic_cases.cpp. Each call to
 * shoc_main loops through nadv=15 steps with dt=5 min. On top of this,
 * shoc_main is called iteratively num_iters=10 steps, performing checks
 * and potentiallywriting output each time. This means that shoc_run_and_cmp
 * is really a single 150-step shoc run.
 */

 /* When called with the below 3 args, compare loops over all variables
  * and calls the above version of "compare" to check for and report
  * large discrepancies.
  */
 Int compare (const double& tol,
             const FortranData::Ptr& ref, const FortranData::Ptr& d) {

  Int nerr = 0;
  FortranDataIterator refi(ref), di(d);
  EKAT_ASSERT(refi.nfield() == di.nfield());
  for (Int i = 0, n = refi.nfield(); i < n; ++i) {
    const auto& fr = refi.getfield(i);
    const auto& fd = di.getfield(i);
    EKAT_ASSERT(fr.size == fd.size);

    // tkh is an input/output of shoc_main() in shoc.F90,
    // but is treated as a local variable in the c++
    // version (tkh values are reset before its used in both).
    // So we just skip the comparison.
    if (fr.name == "tkh") continue;

    nerr += scream::compare(fr.name, fr.data, fd.data, fr.size, tol);
  }
  return nerr;
}

struct Baseline {

  Baseline (const Int nsteps, const Real dt, const Int ncol, const Int nlev, const Int num_qtracers, const Int nadv, const Int repeat)
  {
    params_.push_back({ic::Factory::standard, repeat, nsteps, ncol, nlev, num_qtracers, nadv, dt});
  }

  Int generate_baseline (const std::string& filename) {
    std::ofstream ofile (filename, std::ios::binary);
    EKAT_REQUIRE_MSG (ofile.good(), "generate_baseline can't write '" + filename + "'\n");
    Int nerr = 0;

    // These times are thrown out, I just wanted to be able to use auto
    Int duration = 0;

    for (auto ps : params_) {
      for (Int r = -1; r < ps.repeat; ++r) {
        // Run reference shoc on this set of parameters.
        const auto d = ic::Factory::create(ps.ic, ps.ncol, ps.nlev, ps.num_qtracers);
        set_params(ps, *d);

        if (ps.repeat > 0 && r == -1) {
          std::cout << "Running SHOC with ni=" << d->shcol << ", nk=" << d->nlev
                    << ", dt=" << d->dtime << ", ts=" << ps.nsteps;

          std::cout << ", small_packn=" << SCREAM_SMALL_PACK_SIZE;
          std::cout << std::endl;
        }

        for (int it = 0; it < ps.nsteps; ++it) {
          Int current_microsec = shoc_main(*d);

          if (r != -1 && ps.repeat > 0) { // do not count the "cold" run
            duration += current_microsec;
          }

          if (ps.repeat == 0) {
            write(ofile, d);
          }
        }
      }
      if (ps.repeat > 0) {
        const double report_time = (1e-6*duration) / ps.repeat;

        printf("Time = %1.3e seconds\n", report_time);
      }
    }
    return nerr;
  }

  Int run_and_cmp (const std::string& filename, const double& tol, bool no_baseline) {
    std::ifstream ifile;
    if (!no_baseline) {
      ifile.open(filename,std::ios::binary);
      EKAT_REQUIRE_MSG( ifile.good(), "run_and_cmp can't read '" + filename + "'\n");
    }
    Int nerr = 0, ne;
    int case_num = 0;
    for (auto ps : params_) {
      case_num++;
      if (no_baseline) {
        const auto d = ic::Factory::create(ps.ic, ps.ncol, ps.nlev, ps.num_qtracers);
        set_params(ps, *d);
        for (int it=0; it<ps.nsteps; it++) {
          std::cout << "--- running case # " << case_num << ", timestep # " << it+1 << " of " << ps.nsteps << " ---\n" << std::flush;
          shoc_main(*d);
        }
      }
      else {
        // Read the reference impl's data from the baseline file.
        const auto d_ref = ic::Factory::create(ps.ic, ps.ncol, ps.nlev, ps.num_qtracers);
        set_params(ps, *d_ref);
        // Now run a sequence of other impls. This includes the reference
        // implementation b/c it's likely we'll want to change it as we go.
        {
          const auto d = ic::Factory::create(ps.ic, ps.ncol, ps.nlev, ps.num_qtracers);
          set_params(ps, *d);
          for (int it = 0; it < ps.nsteps; it++) {
            std::cout << "--- checking case # " << case_num << ", timestep # = " << (it+1)*ps.nadv
                      << " ---\n" << std::flush;
            read(ifile, d_ref);
            shoc_main(*d);
            ne = compare(tol, d_ref, d);
            if (ne) std::cout << "Ref impl failed.\n";
            nerr += ne;
          }
        }
      }
    }
    return nerr;
  }

private:
  struct ParamSet {
    ic::Factory::IC ic;
    Int repeat, nsteps, ncol, nlev, num_qtracers, nadv;
    Real dt;
  };

  static void set_params (const ParamSet& ps, FortranData& d) {
    // ncol, nlev, num_qtracers are already set by the Factory.
    d.nadv  = ps.nadv;
    d.dtime = ps.dt;
  }

  std::vector<ParamSet> params_;

  static void write (std::ofstream& ofile, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      impl::write_scalars(ofile,f.dim);
      impl::write_scalars(ofile,f.extent);
      impl::write_scalars(ofile,f.data,f.size);
    }
  }

  static void read (std::ifstream& ifile, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      int dim;
      impl::read_scalars(ifile,dim);
      EKAT_REQUIRE_MSG(dim == f.dim,
                      "For field " << f.name << " read expected dim " <<
                      f.dim << " but got " << dim);
      int extent[3];
      impl::read_scalars(ifile,extent);
      for (int i = 0; i < dim; ++i)
        EKAT_REQUIRE_MSG(extent[i] == f.extent[i],
                        "For field " << f.name << " read expected dim "
                        << i << " to have extent " << f.extent[i] << " but got "
                        << extent[i]);
      impl::read_scalars(ifile,f.data,f.size);
    }
  }
};

void expect_another_arg (int i, int argc) {
  EKAT_REQUIRE_MSG(i != argc-1, "Expected another cmd-line arg.");
}

} // namespace anon

int main (int argc, char** argv) {
  int nerr = 0;

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options] baseline-filename\n"
      "Options:\n"
      "  -g                Generate baseline file.\n"
      "  -c                Compare  baseline file. Default False.\n"
      "  -n                Run without baseline actions. Default True.\n"
      "  -b <baseline_path>  Path to directory containing baselines.\n"
      "  -t <tol>          Tolerance for relative error.\n"
      "  -s <steps>        Number of timesteps. Default=10.\n"
      "  -dt <seconds>     Length of timestep. Default=150.\n"
      "  -i <cols>         Number of columns(ncol). Default=8.\n"
      "  -k <nlev>         Number of vertical levels. Default=72.\n"
      "  -q <num_qtracers> Number of q tracers. Default=3.\n"
      "  -l <nadv>         Number of SHOC loops per timestep. Default=15.\n"
      "  -r <repeat>       Number of repetitions, implies timing run (generate + no I/O). Default=0.\n";

    return 1;
  }

  bool generate = false, no_baseline = true;
  scream::Real tol = SCREAM_BFB_TESTING ? 0 : std::numeric_limits<Real>::infinity();
  Int nsteps = 10;
  Int dt = 150;
  Int ncol = 8;
  Int nlev = 72;
  Int num_qtracers = 3;
  Int nadv = 15;
  Int repeat = 0;
  std::string baseline_fn;
  std::string device;
  for (int i = 1; i < argc; ++i) {
    if (ekat::argv_matches(argv[i], "-g", "--generate")) { generate = true; no_baseline = false; }
    if (ekat::argv_matches(argv[i], "-c", "--compare"))  { no_baseline = false; }
    if (ekat::argv_matches(argv[i], "-b", "--baseline-file")) {
      expect_another_arg(i, argc);
      ++i;
      baseline_fn = argv[i];
    }
    if (ekat::argv_matches(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-s", "--steps")) {
      expect_another_arg(i, argc);
      ++i;
      nsteps = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-dt", "--dt")) {
      expect_another_arg(i, argc);
      ++i;
      dt = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-i", "--ncol")) {
      expect_another_arg(i, argc);
      ++i;
      ncol = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-k", "--nlev")) {
      expect_another_arg(i, argc);
      ++i;
      nlev = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-q", "--num-qtracers")) {
      expect_another_arg(i, argc);
      ++i;
      num_qtracers = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-l", "--nadv")) {
      expect_another_arg(i, argc);
      ++i;
      nadv = std::atoi(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-r", "--repeat")) {
      expect_another_arg(i, argc);
      ++i;
      repeat = std::atoi(argv[i]);
      if (repeat > 0) {
        generate = true;
      }
    }
    if (std::string(argv[i])=="--kokkos-device-id=") {
      auto tokens = ekat::split(argv[i],"=");
      device = tokens[1];
    }
  }

  // Compute full baseline file name with precision.
  baseline_fn += "/shoc_run_and_cmp.baseline" + std::to_string(sizeof(scream::Real));

  scream::initialize_eamxx_session(argc, argv);
  {
    Baseline bln(nsteps, static_cast<Real>(dt), ncol, nlev, num_qtracers, nadv, repeat);
    if (generate) {
      std::cout << "Generating to " << baseline_fn << "\n";
      nerr += bln.generate_baseline(baseline_fn);
    } else if (no_baseline) {
      printf("Running with no baseline actions\n");
      nerr += bln.run_and_cmp(baseline_fn, tol, no_baseline);
    }
    else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol, no_baseline);
    }
  }
  scream::finalize_eamxx_session();

  return nerr != 0 ? 1 : 0;
}
