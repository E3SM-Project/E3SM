#include <catch2/catch.hpp>

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "diagnostics/potential_temperature.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream {


std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm) {

  const int num_local_elems = 4;
  const int np = 4;
  const int nlevs = 32;
  const int num_local_cols = 13;
  const int num_global_cols = num_local_cols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<std::string>("Reference Grid", "Point Grid");
  gm_params.sublist("Mesh Free").set<int>("Number of Global Columns", num_global_cols);
  gm_params.sublist("Mesh Free").set<int>("Number of Local Elements", num_local_elems);
  gm_params.sublist("Mesh Free").set<int>("Number of Vertical Levels", nlevs);
  gm_params.sublist("Mesh Free").set<int>("Number of Gauss Points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_all_grids();

  return gm;
}

TEST_CASE("potential_temperature_diagnostic", "") {

  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager
  auto gm = create_gm(comm);

  // Construct the Diagnostic
  ekat::ParameterList params;
  params.set<std::string>("Diagnostic Name", "Potential Temperature");
  params.set<std::string>("Grid", "Point Grid");
  auto diag = std::make_shared<PotentialTemperatureDiagnostic>(comm,params);
  diag->set_grids(gm);

  // Set the required fields for the diagnostic.
  for (const auto& req : diag->get_required_field_requests()) {
    Field f(req.fid);
    f.allocate_view();
    const auto name = f.name();
    if (name == "T_mid") {
      f.deep_copy(270.0);
    } else if (name == "p_mid") {
      f.deep_copy(10000.0);
    } else {
      REQUIRE(false);
    }
    f.get_header().get_tracking().update_time_stamp(t0);
    diag->set_required_field(f.get_const());
    REQUIRE_THROWS(diag->set_computed_field(f));
  }

  // Initialize the diagnostic
  diag->initialize(t0,RunType::Initial);

  // Run the diagnostic
  const auto& diag_out = diag->get_diagnostic(100.0);
 
  // Finalize the diagnostic
  diag->finalize(); 


}

} //namspace scream
