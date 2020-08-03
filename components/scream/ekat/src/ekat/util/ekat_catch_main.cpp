#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include "ekat/scream_session.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_test_utils.hpp"

#include <mpi.h>

int main (int argc, char **argv) {

  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Read possible scream-specific arguments
  auto const readCommaSeparaterParams = [] (const std::string& cmd_line_arg) {
    if (cmd_line_arg=="") {
      return;
    }
    auto& ts = scream::util::TestSession::get();

    std::stringstream input(cmd_line_arg);
    std::string option;
    while (getline(input,option,',')) {
      // Split option according to key=val
      auto pos = option.find('=');
      scream_require_msg(pos!=std::string::npos, "Error! Badly formatted command line options.\n");
      std::string key = option.substr(0,pos);
      std::string val = option.substr(pos+1);
      scream_require_msg(key!="", "Error! Empty key found in the list of command line options.\n");
      scream_require_msg(val!="", "Error! Empty value for key '" + key + "'.\n");
      ts.params[key] = val;
    }
  };
  Catch::Session catch_session;
  auto cli = catch_session.cli();
  cli |= Catch::clara::Opt(readCommaSeparaterParams, "key1=val1[,key2=val2[,...]]")
             ["--scream-test-params"]
             ("list of parameters to forward to the test");
  catch_session.cli(cli);

  scream_require_msg(catch_session.applyCommandLine(argc,argv)==0,
                     "Error! Something went wrong while parsing command line.\n");

  // If we are on a gpu build, we might have a test device id to use
  // so start by creating a copy of argv that we can extend
  std::vector<char*> args;
  for (int i=0; i<argc; ++i) {
    args.push_back(argv[i]);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int dev_id = scream::util::get_test_device(rank);
  // Create it outside the if, so its c_str pointer survives
  std::string new_arg;
  if (dev_id>=0) {
    new_arg = "--kokkos-device=" + std::to_string(dev_id);
    args.push_back(const_cast<char*>(new_arg.c_str()));
  }

  // Initialize scream;
  scream::initialize_scream_session(args.size(),args.data());

  // Run tests
  int num_failed = catch_session.run(argc, argv);

  // Finalize scream
  scream::finalize_scream_session();

  // Finalize MPI
  MPI_Finalize();

  // Return test result
  return num_failed != 0 ? 1 : 0;
}
