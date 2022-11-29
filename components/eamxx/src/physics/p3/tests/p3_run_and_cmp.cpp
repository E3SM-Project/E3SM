#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

#include "p3_main_wrap.hpp"
#include "p3_functions_f90.hpp"
#include "p3_ic_cases.hpp"

#include "ekat/util/ekat_file_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <chrono>
#include <vector>

namespace {
using namespace scream;
using namespace scream::p3;

/*
 * p3_run_and_cmp can be run in 2 modes. First, generate_baseline
 * runs the baseline (aka reference, probably git master) version of
 * the code and saves its output as a raw binary file. Then run_and_cmp
 * runs the new/experimental version of the code and compares it against
 * the baseline data you've saved to file. Both baseline and cmp modes
 * start from an initial condition in ../p3_ic_cases.cpp. By default,
 * tests are run with log_PredictNc=true and false and also for 1 step or
 * 6 steps. This creates a total of 4 different cases. When 6 steps are run,
 * output for each step is considered separately.
 */


/* Given a column of data for variable "label" from the reference run
 * (probably master) and from your new exploratory run, loop over all
 * heights and confirm whether or not the relative difference between
 * runs is within tolerance "tol". If not, print debug info. Here, "a"
 * is the value from the reference run and "b" is from the new run.
 */
template <typename Scalar>
static Int compare (const std::string& label, const Scalar* a,
                    const Scalar* b, const Int& n, const Real& tol) {

  Int nerr1 = 0;
  Int nerr2 = 0;
  Real den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Real worst = 0;
  for (Int i = 0; i < n; ++i) {
    if (std::isnan(a[i]) || std::isinf(a[i]) ||
        std::isnan(b[i]) || std::isinf(b[i])) {
      ++nerr1;
      continue;
    }

    // The code below is to force a result difference. This is used by the
    // scream/scripts internal testing to verify that various DIFFs are detected.
#if defined(SCREAM_FORCE_RUN_DIFF)
    const Real num = std::abs(a[i]*Real(1.2) - b[i]);
#else
    const Real num = std::abs(a[i] - b[i]);
#endif
    if (num > tol*den) {
      ++nerr2;
      worst = std::max(worst, num);
    }
  }

  if (nerr1) {
    std::cout << label << " has " << nerr1 << " infs + nans.\n";

  }

  if (nerr2) {
    std::cout << label << " > tol " << nerr2 << " times. Max rel diff= " << (worst/den)
             << " normalized by ref impl val=" << den << ".\n";

  }

  return nerr1 + nerr2;
}

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
    nerr += compare(fr.name, fr.data, fd.data, fr.size, tol);
  }
  return nerr;
}

struct Baseline {
  Baseline (const Int nsteps, const Real dt, const Int ncol, const Int nlev, const Int repeat, const std::string predict_nc, const std::string prescribed_CCN)
  {
    //If predict_nc="both", start looping at i_start=0 (false) and end after i_start=1 (true)
    //otherwise, modify start and end to only loop over case of interest. Test that predict_nc
    //is only yes, no, or both is done below so not bothering here.
    int i_start=0;
    int i_end=2;
    if (predict_nc == "yes") {i_start=1;}
    if (predict_nc == "no" ) {i_end=1;}

    //Like above, get indexing for prescribe_CCN options
    int j_start=0;
    int j_end=2;
    if (prescribed_CCN == "yes") {j_start=1;}
    if (prescribed_CCN == "no" ) {j_end=1;}


    for (int i = i_start; i < i_end; ++i) { // predict_nc is false or true
      for (int j = j_start; j< j_end; ++j) { //prescribed_CCN is false or true
  //                 initial condit,     repeat, nsteps, ncol, nlev, dt, prescribe or predict nc, prescribe CCN or not
  params_.push_back({ic::Factory::mixed, repeat, nsteps, ncol, nlev, dt, i>0,                     j>0 });
      }
    }
  }

  Int generate_baseline (const std::string& filename, bool use_fortran) {
    auto fid = ekat::FILEPtr(fopen(filename.c_str(), "w"));
    EKAT_REQUIRE_MSG( fid, "generate_baseline can't write " << filename);
    Int nerr = 0;

    Int total_duration_microsec = 0;

    for (auto ps : params_) {
      // Run reference p3 on this set of parameters.
      for (Int r = -1; r < ps.repeat; ++r) {
        const auto d = ic::Factory::create(ps.ic, ps.ncol, ps.nlev);
        set_params(ps, *d);
        p3_init();

        if (ps.repeat > 0 && r == -1) {
          std::cout << "Running P3 with ni=" << d->ncol << ", nk=" << d->nlev
                    << ", dt=" << d->dt << ", ts=" << d->it
                    << ", predict_nc=" << d->do_predict_nc
                    << ", prescribed_CCN=" << d->do_prescribed_CCN;

          if (!use_fortran) {
            std::cout << ", small_packn=" << SCREAM_SMALL_PACK_SIZE;
          }
          std::cout << std::endl;
        }

        for (int it=0; it<ps.nsteps; it++) {
          Int current_microsec = p3_main_wrap(*d, use_fortran);

          if (r != -1 && ps.repeat > 0) { // do not count the "cold" run
            total_duration_microsec += current_microsec;
          }

          if (ps.repeat == 0) {
            write(fid, d); // Save the fields to the baseline file.
          }
        }
      }

      if (ps.repeat > 0) {
        const double report_time = (1e-6*total_duration_microsec) / ps.repeat;

        printf("Time = %1.3e seconds\n", report_time);
      }
    }
    return nerr;
  }

  Int run_and_cmp (const std::string& filename, const double& tol, bool use_fortran) {
    auto fid = ekat::FILEPtr(fopen(filename.c_str(), "r"));
    EKAT_REQUIRE_MSG( fid, "generate_baseline can't read " << filename);
    Int nerr = 0, ne;
    int case_num = 0;
    for (auto ps : params_) {
      case_num++;
      // Read the reference impl's data from the baseline file.
      const auto d_ref = ic::Factory::create(ps.ic, ps.ncol, ps.nlev);
      set_params(ps, *d_ref);
      // Now run a sequence of other impls. This includes the reference
      // implementation b/c it's likely we'll want to change it as we go.
      {
        const auto d = ic::Factory::create(ps.ic, ps.ncol, ps.nlev);
        set_params(ps, *d);
        p3_init();
        for (int it=0; it<ps.nsteps; it++) {
          std::cout << "--- checking case # " << case_num << ", timestep # " << it+1 << " of " << ps.nsteps << " ---\n" << std::flush;
          read(fid, d_ref);
          p3_main_wrap(*d, use_fortran);
          ne = compare(tol, d_ref, d);
          if (ne) std::cout << "Ref impl failed.\n";
          nerr += ne;
        }
      }
    }
    return nerr;
  }

private:

  // Full specification for a run
  struct ParamSet {
    ic::Factory::IC ic;
    Int repeat, nsteps, ncol, nlev;
    Real dt;
    bool do_predict_nc, do_prescribed_CCN;
  };

  static void set_params (const ParamSet& ps, FortranData& d) {
    // Items not set by factory
    d.dt                = ps.dt;
    d.it                = ps.nsteps;
    d.do_predict_nc     = ps.do_predict_nc;
    d.do_prescribed_CCN = ps.do_prescribed_CCN;
  }

  std::vector<ParamSet> params_;

  static void write (const ekat::FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      ekat::write(&f.dim, 1, fid);
      ekat::write(f.extent, f.dim, fid);
      ekat::write(f.data, f.size, fid);
    }
  }

  static void read (const ekat::FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      int dim, ds[3];
      ekat::read(&dim, 1, fid);
      EKAT_REQUIRE_MSG(dim == f.dim,
                      "For field " << f.name << " read expected dim " <<
                      f.dim << " but got " << dim);
      ekat::read(ds, dim, fid);
      for (int i = 0; i < dim; ++i)
        EKAT_REQUIRE_MSG(ds[i] == f.extent[i],
                        "For field " << f.name << " read expected dim "
                        << i << " to have extent " << f.extent[i] << " but got "
                        << ds[i]);
      ekat::read(f.data, f.size, fid);
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
      "  -g                  Generate baseline file. Default False.\n"
      "  -f                  Use fortran impls instead of c++. Default False.\n"
      "  -t <tol>            Tolerance for relative error. Default 0.\n"
      "  -s <steps>          Number of timesteps. Default=6.\n"
      "  -dt <seconds>       Length of timestep. Default=300.\n"
      "  -i <cols>           Number of columns. Default=3.\n"
      "  -k <nlev>           Number of vertical levels. Default=72.\n"
      "  -r <repeat>         Number of repetitions, implies timing run (generate + no I/O). Default=0.\n"
      "  -p <predict_nc>     yes|no|both. Default=both.\n"
      "  -c <prescribed_ccn> yes|no|both. Default=both.\n";
    return 1;
  }

  bool generate = false, use_fortran = false;
  scream::Real tol = SCREAM_BFB_TESTING ? 0 : std::numeric_limits<Real>::infinity();
  Int timesteps = 6;
  Int dt = 300;
  Int ncol = 3;
  Int nlev = 72;
  Int repeat = 0;
  std::string device;
  std::string predict_nc = "both";
  std::string prescribed_ccn = "both";
  std::string baseline_fn;
  for (int i = 1; i < argc-1; ++i) {
    if (ekat::argv_matches(argv[i], "-g", "--generate")) generate = true;
    if (ekat::argv_matches(argv[i], "-f", "--fortran")) use_fortran = true;
    if (ekat::argv_matches(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
    if (ekat::argv_matches(argv[i], "-b", "--baseline-file")) {
      expect_another_arg(i, argc);
      ++i;
      baseline_fn = argv[i];
    }
    if (ekat::argv_matches(argv[i], "-s", "--steps")) {
      expect_another_arg(i, argc);
      ++i;
      timesteps = std::atoi(argv[i]);
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
    if (std::string(argv[i])=="--ekat-kokkos-device") {
      expect_another_arg(i, argc);
      ++i;
      device = argv[i];
    }
    if (ekat::argv_matches(argv[i], "-r", "--repeat")) {
      expect_another_arg(i, argc);
      ++i;
      repeat = std::atoi(argv[i]);
      if (repeat > 0) {
        generate = true;
      }
    }
    if (ekat::argv_matches(argv[i], "-p", "--predict-nc")) {
      expect_another_arg(i, argc);
      ++i;
      predict_nc = std::string(argv[i]);
      EKAT_REQUIRE_MSG(predict_nc == "yes" || predict_nc == "no" || predict_nc == "both",
                       "Predict option value must be one of yes|no|both");
    }
    if (ekat::argv_matches(argv[i], "-c", "--prescribed-ccn")) {
      expect_another_arg(i, argc);
      ++i;
      prescribed_ccn = std::string(argv[i]);
      EKAT_REQUIRE_MSG(prescribed_ccn == "yes" || prescribed_ccn == "no" || prescribed_ccn == "both",
                       "Prescribed CCN option value must be one of yes|no|both");
    }
  }

  // Decorate baseline name with precision.
  baseline_fn += std::to_string(sizeof(scream::Real));

  std::vector<char*> args;
  for (int i=0; i<argc; ++i) {
    args.push_back(argv[i]);
  }

  // If "--ekat-kokkos-device <N>" was specified, add kokkos
  // initialization flag to argv
  // Create it outside the if, so its c_str pointer survives
  std::string dev_arg;
  if (device!="") {
    auto is_int = [] (const std::string& s)->bool {
      std::istringstream is(s);
      int d;
      is >> d;
      return !is.fail() && is.eof();
    };

    EKAT_REQUIRE_MSG (is_int(device), "Error! Invalid device specification.\n");

    if (std::stoi(device) != -1) {
      dev_arg = "--kokkos-device-id=" + device;
      args.push_back(const_cast<char*>(dev_arg.c_str()));
    }
  }

  scream::initialize_scream_session(args.size(), args.data()); {
    Baseline bln(timesteps, static_cast<Real>(dt), ncol, nlev, repeat, predict_nc, prescribed_ccn);
    if (generate) {
      std::cout << "Generating to " << baseline_fn << "\n";
      nerr += bln.generate_baseline(baseline_fn, use_fortran);
    } else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol, use_fortran);
    }
    P3GlobalForFortran::deinit();
  } scream::finalize_scream_session();

  return nerr != 0 ? 1 : 0;
}
