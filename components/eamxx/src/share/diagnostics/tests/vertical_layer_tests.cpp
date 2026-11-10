#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/grid/point_grid.hpp"
#include "share/physics/physics_constants.hpp"

namespace scream {


//-----------------------------------------------------------------------------------------------//
template<typename DeviceT, int N>
void run (const std::string& diag_name, const std::string& location)
{
  using PC = scream::physics::Constants<Real>;

  const int packsize = N;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  
  // Construct the Diagnostic
  ekat::ParameterList params;

  params.set<std::string>("diag_name", diag_name);
  params.set<std::string>("vert_location",location);
  auto diag = diag_factory.create("VerticalLayer",comm,params, grid);
  const bool needs_phis = diag_name=="z" or diag_name=="geopotential";

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    auto scalar2d = grid->get_2d_scalar_layout();
    for (const auto& fname : diag->get_input_fields_names()) {
      auto layout = (fname == "phis") ? scalar2d : scalar3d;
      FieldIdentifier fid(fname, layout, Pa, grid->name());
      Field f(fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    }
  }

  // Note: we are not testing the calculate_dz utility. We are testing
  //       the diag class, so use some inputs that make checking results easier
  //       With these inputs, T_virt=T, and dz=8*rd/g
  const Real g = PC::gravit.value;
  const Real rho_val = 4;
  const Real qv_val = 0;
  const Real p_val = 2;
  const Real T_val = 4;
  const Real c1 = -PC::ONE + PC::ONE / PC::ep_2.value;
  const Real Tvirt_val = T_val*(PC::ONE + c1*qv_val);
  const Real dz_val = (PC::RD.value/g) * rho_val*Tvirt_val / p_val;
  const Real phis_val = 3;

  input_fields["T_mid"].deep_copy(T_val);
  input_fields["p_mid"].deep_copy(p_val);
  input_fields["pseudo_density"].deep_copy(rho_val);
  input_fields["qv"].deep_copy(qv_val);
  if (needs_phis) {
    input_fields["phis"].deep_copy(phis_val);
  }

  // Initialize and run the diagnostic
  diag->initialize();
  diag->compute(t0);
  const auto& diag_out = diag->get();
  diag_out.sync_to_host();
  auto d_h = diag_out.get_view<Real**,Host>();

  // Compare against expected value
  const auto last_int = nlevs;
  const auto last_mid = last_int-1;

  // Precompute surface value and increment depending on the diag type
  Real delta, surf_val;
  if (diag_name=="height") {
    surf_val = 0;
    delta = dz_val;
  } else if (diag_name=="z") {
    surf_val = phis_val/g;
    delta = dz_val;
  } else {
    surf_val = phis_val;
    delta = dz_val*g;
  }

  for (int icol=0; icol<ncols; ++icol) {
    Real prev_int_val = surf_val;

    if (location=="interfaces") {
      // Check surface value
      REQUIRE (d_h(icol,nlevs)==prev_int_val);
    }

    for (int ilev=last_mid; ilev>=0; --ilev) {
      if (diag_name=="dz") {
        REQUIRE (d_h(icol,ilev)==dz_val);
      } else {
        // If interface, check value, otherwise perform int->mid averaging and check value
        auto int_val = prev_int_val + delta;
        if (location=="interfaces") {
          REQUIRE_THAT(d_h(icol,ilev), Catch::Matchers::WithinRel(int_val,Real(1e-5)));
        } else {
          auto mid_val = (int_val + prev_int_val) / 2;
          REQUIRE_THAT(d_h(icol,ilev), Catch::Matchers::WithinRel(mid_val,Real(1e-5)));
        }
        prev_int_val = int_val;
      }
    }
  }

} // run()

TEST_CASE("vertical_layer_test", "vertical_layer_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  ekat::Comm comm(MPI_COMM_WORLD);
  auto root_print = [&](const std::string& msg) {
    if (comm.am_i_root()) {
      printf("%s",msg.c_str());
    }
  };
  auto do_run = [&] (auto int_const) {
    constexpr int N = decltype(int_const)::value;
    root_print("\n");
    root_print(" -> Testing diagnostic for pack_size=" + std::to_string(N) + "\n");
    for (std::string loc : {"midpoints","interfaces"}) {
      for (std::string diag : {"geopotential","height","z"}) {
        std::string msg = "    -> Testing diag=" + diag + " at " + loc + " ";
        std::string dots (50-msg.size(),'.');
        root_print (msg + dots + "\n");
        run<Device,N>(diag, loc);
        root_print (msg + dots + " PASS!\n");
      }
    }
    std::string msg = "    -> Testing diag=dz ";
    std::string dots (50-msg.size(),'.');
    root_print (msg + dots + "\n");
    run<Device, N>("dz", "midpoints");
    root_print (msg + dots + " PASS!\n");
  };

  if (SCREAM_PACK_SIZE!=1) {
    do_run(std::integral_constant<int,1>());
  }
  do_run(std::integral_constant<int,SCREAM_PACK_SIZE>());

} // TEST_CASE

} // namespace
