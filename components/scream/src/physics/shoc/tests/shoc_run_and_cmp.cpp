#include "ekat/scream_session.hpp"
#include "ekat/util/file_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"

#include "physics/shoc/shoc_f90.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/shoc/shoc_ic_cases.hpp"

#include <vector>

namespace {
using namespace scream;
using namespace scream::util;
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
    
    const auto num = std::abs(a[i] - b[i]);
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
  scream_assert(refi.nfield() == di.nfield());
  for (Int i = 0, n = refi.nfield(); i < n; ++i) {
    const auto& fr = refi.getfield(i);
    const auto& fd = di.getfield(i);
    scream_assert(fr.size == fd.size);
    nerr += compare(fr.name, fr.data, fd.data, fr.size, tol);
  }
  return nerr;
}

struct Baseline {

  // Number of iterations (nadv steps of size dtime per iteration).
  const int num_iters = 10;

  Baseline () {
    //                 ic, shcol, nlev, num_qtracers, nadv, dtime
    params_.push_back({ic::Factory::standard, 8, 72, 3, 15, 150});
  }

  Int generate_baseline (const std::string& filename, bool use_fortran) {
    auto fid = FILEPtr(fopen(filename.c_str(), "w"));
    scream_require_msg( fid, "generate_baseline can't write " << filename);
    Int nerr = 0;
    for (auto ps : params_) {
      // Run reference shoc on this set of parameters.
      const auto d = ic::Factory::create(ps.ic, ps.shcol, ps.nlev, ps.num_qtracers);
      set_params(ps, *d);
      shoc_init(ps.nlev, use_fortran);
      for (int it = 0; it < num_iters; ++it) {
        shoc_main(*d);
        write(fid, d);
      }
    }
    return nerr;
  }

  Int run_and_cmp (const std::string& filename, const double& tol, bool use_fortran) {
    auto fid = FILEPtr(fopen(filename.c_str(), "r"));
    scream_require_msg( fid, "generate_baseline can't read " << filename);
    Int nerr = 0, ne;
    int case_num = 0;
    for (auto ps : params_) {
      case_num++;
      // Read the reference impl's data from the baseline file.
      const auto d_ref = ic::Factory::create(ps.ic, ps.shcol, ps.nlev, ps.num_qtracers);
      set_params(ps, *d_ref);
      // Now run a sequence of other impls. This includes the reference
      // implementation b/c it's likely we'll want to change it as we go.
      {
        const auto d = ic::Factory::create(ps.ic, ps.shcol, ps.nlev, ps.num_qtracers);
        set_params(ps, *d);
        shoc_init(ps.nlev, use_fortran);
        for (int it = 0; it < num_iters; it++) {
          std::cout << "--- checking case # " << case_num << ", timestep # = " << (it+1)*ps.nadv
                     << " ---\n" << std::flush;
          read(fid, d_ref);
          shoc_main(*d);
          ne = compare(tol, d_ref, d);
          if (ne) std::cout << "Ref impl failed.\n";
          nerr += ne;
        }
      }
    }
    return nerr;
  }

private:

  struct ParamSet {
    ic::Factory::IC ic;
    Int shcol, nlev, num_qtracers, nadv;
    Real dtime;
  };

  static void set_params (const ParamSet& ps, FortranData& d) {
    // shcol, nlev, num_qtracers are already set by the Factory.
    d.nadv = ps.nadv;
    d.dtime = ps.dtime;
  }

  std::vector<ParamSet> params_;

  static void write (const FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      util::write(&f.dim, 1, fid);
      util::write(f.extent, f.dim, fid);
      util::write(f.data, f.size, fid);
    }
  }

  static void read (const FILEPtr& fid, const FortranData::Ptr& d) {
    FortranDataIterator fdi(d);
    for (Int i = 0, n = fdi.nfield(); i < n; ++i) {
      const auto& f = fdi.getfield(i);
      int dim, ds[3];
      util::read(&dim, 1, fid);
      scream_require_msg(dim == f.dim,
                      "For field " << f.name << " read expected dim " <<
                      f.dim << " but got " << dim);
      util::read(ds, dim, fid);
      for (int i = 0; i < dim; ++i)
        scream_require_msg(ds[i] == f.extent[i],
                        "For field " << f.name << " read expected dim "
                        << i << " to have extent " << f.extent[i] << " but got "
                        << ds[i]);
      util::read(f.data, f.size, fid);
    }
  }
};

void expect_another_arg (int i, int argc) {
  scream_require_msg(i != argc-1, "Expected another cmd-line arg.");
}

} // namespace anon

int main (int argc, char** argv) {
  int nerr = 0;

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options] baseline-filename\n"
      "Options:\n"
      "  -g        Generate baseline file.\n"
      "  -f        Use fortran impls instead of c++.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return 1;
  }

  bool generate = false, use_fortran = false;
  scream::Real tol = 0;
  for (int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-g", "--generate")) generate = true;
    if (util::eq(argv[i], "-f", "--fortran")) use_fortran = true;
    if (util::eq(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
  }

  // Decorate baseline name with precision.
  std::string baseline_fn(argv[argc-1]);
  baseline_fn += std::to_string(sizeof(scream::Real));

  scream::initialize_scream_session(argc, argv); {
    Baseline bln;
    if (generate) {
      std::cout << "Generating to " << baseline_fn << "\n";
      nerr += bln.generate_baseline(baseline_fn, use_fortran);
    } else {
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      nerr += bln.run_and_cmp(baseline_fn, tol, use_fortran);
    }
  } scream::finalize_scream_session();

  return nerr != 0 ? 1 : 0;
}
