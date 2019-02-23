#include <catch2/catch.hpp>
#include <control/atmosphere_driver.hpp>

namespace {

TEST_CASE("process_factory", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a parameter list for inputs
  ParameterList ad_params("Atmosphere Driver");
  auto& params = ad_params.sublist("Atmosphere Processes");

  params.set("Number of Entries",2);
  params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Surface Coupling");

  auto& p1 = params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Dynamics");

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  AtmosphereDriver ad;
  ad.initialize(atm_comm,ad_params);

  // If we get here, the ad initialized correctly
  REQUIRE (true);
}

} // empty namespace
