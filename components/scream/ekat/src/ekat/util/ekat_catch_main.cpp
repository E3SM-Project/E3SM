#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include "ekat/scream_session.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_test_utils.hpp"

#include <mpi.h>

int main (int argc, char **argv) {

  // Initialize MPI
  MPI_Init(&argc,&argv);

  // Initialize scream;
  scream::initialize_scream_session(argc, argv);

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

  // Run tests
  int num_failed = catch_session.run(argc, argv);

  // Finalize scream
  scream::finalize_scream_session();

  // Finalize MPI
  MPI_Finalize();

  // Return test result
  return num_failed != 0 ? 1 : 0;
}
