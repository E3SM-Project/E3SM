#include <string>
#include <fstream>
#include <vector>

using Real = double;
using Int  = int;

namespace {

struct FakeBaselineTest
{
  FakeBaselineTest () {}

  Int generate_baseline (const std::string& filename) {
    std::ofstream ofile(filename,std::ios::binary);
    EKAT_REQUIRE_MSG( ofile.good(), "generate_baseline can't write " + filename + "\n");


    Real result = 42;
    ofile << result;

    return 0;
  }

  Int run_and_cmp (const std::string& filename, const Scalar& tol, bool no_baseline) {
    std::ifstream ifile;
    if (!no_baseline) {
      ifile.open(filename,std::ios::binary);
      EKAT_REQUIRE_MSG( ifile.good(), "run_and_cmp can't read " + filename + "\n");
    }

    Real result = 42;

    if (!no_baseline) {
      Real ref;
      ifile >> ref;
      if (ref != result) {
        return 1;
      }
    }

    return 0;
  }
};


void expect_another_arg (int i, int argc) {
  EKAT_REQUIRE_MSG(i != argc-1, "Expected another cmd-line arg.");
}

} // empty namespace

int main (int argc, char** argv) {
  int nerr = 0;

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options]\n"
      "Options:\n"
      "  -g                  Generate baseline file. Default False.\n"
      "  -c                  Compare baseline file. Default False.\n"
      "  -n                  Run without baseline actions. Default True.\n"
      "  -b <baseline_path>  Path to directory containing baselines.\n"
      "  -t <tol>            Tolerance for relative error. Default machine eps (*10000 for Release).\n";
    return 1;
  }

  // Set up options with defaults
  bool generate = false, no_baseline = true;
  Real tol = 0;
  std::string baseline_fn;

  // Parse options
  for (int i = 1; i < argc-1; ++i) {
    if (ekat::argv_matches(argv[i], "-g", "--generate")) {generate = true; no_baseline = false;}
    if (ekat::argv_matches(argv[i], "-c", "--compare"))  {no_baseline = false;}
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
  }

  // Compute full baseline file name with precision.
  baseline_fn += "/baseline_test.baseline" + std::to_string(sizeof(Real));

  std::vector<char*> args;
  for (int i=0; i<argc; ++i) {
    args.push_back(argv[i]);
  }

  FakeBaselineTest bln;
  if (generate) {
    std::cout << "Generating to " << baseline_fn << "\n";
    nerr += bln.generate_baseline(baseline_fn);
  } else if (no_baseline) {
    printf("Running with no baseline actions\n");
    nerr += bln.run_and_cmp(baseline_fn, tol, no_baseline);
  } else {
    printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
    nerr += bln.run_and_cmp(baseline_fn, tol, no_baseline);
  }

  return nerr != 0 ? 1 : 0;
}
