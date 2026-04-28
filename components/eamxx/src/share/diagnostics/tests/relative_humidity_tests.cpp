#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/physics_functions.hpp" 

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_pack.hpp>

namespace scream {

template<typename DeviceT>
void run()
{
  using PF      = scream::physics::Functions<Real, DefaultDevice>;
  using PC      = scream::physics::Constants<Real>;
  using KT      = ekat::KokkosTypes<DeviceT>;
  using Pack    = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // Create a grid
  const int ncols = 1;
  const int nlevs = 2*SCREAM_PACK_SIZE+1;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("RelativeHumidity",comm,params, grid);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    for (const auto& fname : diag->get_input_fields_names()) {
      FieldIdentifier fid(fname, scalar3d, Pa, grid->name());
      Field f(fid);
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    }
  }
  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Randomize inputs (except for dpdry)
    const auto& T_mid_f     = input_fields["T_mid"];
    const auto& p_dry_mid_f = input_fields["p_dry_mid"];
    const auto& dpwet_f = input_fields["pseudo_density"];
    const auto& dpdry_f = input_fields["pseudo_density_dry"];
    const auto& qv_f = input_fields["qv"];

    randomize_uniform(T_mid_f,seed++,200,400);
    randomize_uniform(p_dry_mid_f,seed++,10,PC::P0.value);
    randomize_uniform(qv_f,seed++,0,1e-2);
    randomize_uniform(dpwet_f,seed++,1,100);

    const auto& T_mid_v     = T_mid_f.get_view<Pack**>();
    const auto& p_dry_mid_v = p_dry_mid_f.get_view<Pack**>();
    const auto& dpwet_v = dpwet_f.get_view<Pack**>();
    const auto& dpdry_v = dpdry_f.get_view<Pack**>();
    const auto& qv_v = qv_f.get_view<Pack**>();

    // Run diagnostic and compare with manual calculation
    Field rh_f = T_mid_f.clone();
    rh_f.deep_copy(0);
    const auto& rh_v = rh_f.get_view<Pack**>();
    auto manual = KOKKOS_LAMBDA(const int icol, const int jpack) {
      auto range_pack = ekat::range<Pack>(jpack*Pack::n);
      auto range_mask = range_pack < nlevs;
      dpdry_v(icol,jpack) = dpwet_v(icol,jpack) - dpwet_v(icol,jpack)*qv_v(icol,jpack);
      auto qv_sat_l = PF::qv_sat_dry(T_mid_v(icol,jpack), p_dry_mid_v(icol,jpack), true, range_mask);
      qv_sat_l *=  dpdry_v(icol,jpack) ;
      qv_sat_l /=  dpwet_v(icol,jpack) ;
      rh_v(icol,jpack) = qv_v(icol,jpack)/qv_sat_l;
    };
    const int npacks = ekat::PackInfo<SCREAM_PACK_SIZE>::num_packs(nlevs);
    MDRange policy({0,0},{ncols,npacks});
    Kokkos::parallel_for("", policy, manual);
    Kokkos::fence();

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic(t0);
    const auto& diag_out = diag->get_diagnostic();

    REQUIRE(views_are_equal(diag_out,rh_f));
  }
}

TEST_CASE("relative_humidity_test", "relative_humidity_test]")
{
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for (int irun=0; irun<num_runs; ++irun) {
    printf(" -> Test run %d\n",irun);
    run<Device>();
  }
} // TEST_CASE

} // namespace scream
