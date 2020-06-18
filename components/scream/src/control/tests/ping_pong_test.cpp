#include <catch2/catch.hpp>
#include "share/grid/user_provided_grids_manager.hpp"
#include "ekat/scream_parse_yaml_file.hpp"
#include "control/atmosphere_driver.hpp"

#include "dummy_atm_setup.hpp"

#include <numeric>

namespace scream {

TEST_CASE("ping-pong", "") {
  using namespace scream;

  constexpr int num_cols   = 2;

  // Load ad parameter list
  std::string fname = "ping_pong.yaml";
  ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );
  const int vec_length = ad_params.sublist("Atmosphere Processes").sublist("Process 0").get<int>("Number of vector components");

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Setup the atm factories and grid manager
  dummy_atm_init (num_cols);

  // Create the driver
  control::AtmosphereDriver ad;

  // Init and run (to finalize, wait till checks are completed,
  // or you'll clear the field repo!)
  util::TimeStamp init_time(2019,0,0,0);
  util::TimeStamp end_time (2019,0,1,0);
  const Real dt = 3500.0; // This should trigger an adjustment of the last time step
  ad.initialize(atm_comm,ad_params,init_time);

  // Fill the field with initial guess
  std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
  std::vector<int> dims = {num_cols, vec_length};
  FieldLayout layout (tags,dims);
  FieldIdentifier fid("field_0",layout,units::m,"Physics_fwd");
  const auto& repo = ad.get_field_repo();
  const auto& field = repo.get_field(fid);
  auto d_view = field.get_view();
  auto h_view = Kokkos::create_mirror_view(d_view);
  const int size = h_view.size();
  decltype(h_view) answer("",size);

  Kokkos::deep_copy(d_view,0.0);
  Kokkos::deep_copy(h_view,0.0);
#ifdef SCREAM_DEBUG
  Kokkos::deep_copy(ad.get_bkp_field_repo().get_field(fid).get_view(),d_view);
#endif

  // Run
  int iter = 0;
  for (auto time=init_time; time<end_time; time+=dt, ++iter) {
    const auto& ts = ad.get_atm_time_stamp();
    std::cout << " -------------------------------------------------------\n";
    std::cout << "   Start of atmosphere time step\n";
    std::cout << "    - current time: " << ts.to_string() << "\n";
    Real dt_adj = dt;
    bool adjusted = false;
    if (end_time<(ts+dt)) {
      dt_adj = end_time-ts;
      adjusted = true;
    }
    std::cout << "    - time step: " << dt_adj << (adjusted ? "s (adjusted to hit final time)\n" : "s\n");
    std::cout << "    - end time: " << (ts+dt_adj).to_string() << "\n";
    ad.run(dt_adj);

    // Prepare the answer
    // Every atm proc does out(:) = in(:) op 2.0, where op is +,*,-,/
    // in a cyclic fashion, in that order (one op per iter).
    // Since there are two atm procs, doing the very same op,
    // we can exactly compute what the answer should be.
    int rem = iter % 4;
    for (int i=0; i<size; ++i) {
      switch (rem) {
        // Do update twice, since both atm procs do the same update
        case 0: answer[i] += 2.0; answer[i] += 2.0; break;
        case 1: answer[i] *= 2.0; answer[i] *= 2.0; break;
        case 2: answer[i] -= 2.0; answer[i] -= 2.0; break;
        case 3: answer[i] /= 2.0; answer[i] /= 2.0; break;
      }
    }
  }

  // Check the answer
  Kokkos::deep_copy(h_view,d_view);
  for (int i=0; i<h_view.extent_int(0); ++i) {
    REQUIRE (h_view(i) == answer[i]);
  }

  // Finalize and clean up
  ad.finalize();

  // Clean up factories and grids manager
  dummy_atm_cleanup ();
}

} // empty namespace
