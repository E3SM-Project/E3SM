#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/util/eamxx_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_team_policy_utils.hpp>

#include <iomanip>

namespace scream {


//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using TPF        = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using MemberType = typename KT::MemberType;

  constexpr Real gravit = PC::gravit.value;
  constexpr Real macheps = PC::macheps;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_alpha(0.1,0.9);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  auto& diag_factory = DiagnosticFactory::instance();
  register_diagnostics();
  ekat::ParameterList params;

  REQUIRE_THROWS (diag_factory.create("WaterPath",comm,params, grid)); // No 'water_kind'
  params.set<std::string>("water_kind","Foo");
  REQUIRE_THROWS (diag_factory.create("WaterPath",comm,params, grid)); // Invalid 'water_kind'

  std::map<std::string,std::shared_ptr<AbstractDiagnostic>> diags;

  // Vapor
  params.set<std::string>("water_kind","Vap");
  auto diag_vap = diag_factory.create("WaterPath",comm,params, grid);
  diags.emplace("vwp",diag_vap);
  // Liquid
  params.set<std::string>("water_kind","Liq");
  auto diag_liq = diag_factory.create("WaterPath",comm,params, grid);
  diags.emplace("lwp",diag_liq);
  // Ice
  params.set<std::string>("water_kind","Ice");
  auto diag_ice = diag_factory.create("WaterPath",comm,params, grid);
  diags.emplace("iwp",diag_ice);
  // Rime
  params.set<std::string>("water_kind","Rime");
  auto diag_rime = diag_factory.create("WaterPath",comm,params, grid);
  diags.emplace("mwp",diag_rime);
  // Rain
  params.set<std::string>("water_kind","Rain");
  auto diag_rain = diag_factory.create("WaterPath",comm,params, grid);
  diags.emplace("rwp",diag_rain);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    for (const auto& dd : diags) {
      const auto& diag = dd.second;
      for (const auto& fname : diag->get_input_fields_names()) {
        if (input_fields.find(fname)==input_fields.end()) {
          FieldIdentifier fid(fname, scalar3d, Pa, grid->name());
          Field f(fid);
          f.allocate_view();
          f.get_header().get_tracking().update_time_stamp(t0);
          input_fields.emplace(fname, f);
        }
        diag->set_input_field(input_fields.at(fname));
      }
      // Initialize the diagnostic
      diag->initialize();
    }
  }


  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    auto& qv_f = input_fields.at("qv");
    auto& qc_f = input_fields.at("qc");
    auto& qi_f = input_fields.at("qi");
    auto& qm_f = input_fields.at("qm");
    auto& qr_f = input_fields.at("qr");
    auto& dp_f = input_fields.at("pseudo_density");

    randomize_uniform(qv_f,seed++,1e-12,1e-6);
    randomize_uniform(dp_f,seed++,1,100);
    randomize_uniform(qc_f,seed++,0,1e-3);
    randomize_uniform(qi_f,seed++,0,1e-3);
    randomize_uniform(qm_f,seed++,0,1e-3);
    randomize_uniform(qr_f,seed++,0,1e-3);

    // Grab views for each of the water path diagnostics
    const auto& vwp = diags["vwp"]->get();
    const auto& lwp = diags["lwp"]->get();
    const auto& iwp = diags["iwp"]->get();
    const auto& mwp = diags["mwp"]->get();
    const auto& rwp = diags["rwp"]->get();
    const auto& vwp_v = vwp.get_view<const Real*>();
    const auto& lwp_v = lwp.get_view<const Real*>();
    const auto& iwp_v = iwp.get_view<const Real*>();
    const auto& mwp_v = mwp.get_view<const Real*>();
    const auto& rwp_v = rwp.get_view<const Real*>();
    const auto& vwp_h = vwp.get_view<const Real*,Host>();
    const auto& lwp_h = lwp.get_view<const Real*,Host>();
    const auto& iwp_h = iwp.get_view<const Real*,Host>();
    const auto& mwp_h = mwp.get_view<const Real*,Host>();
    const auto& rwp_h = rwp.get_view<const Real*,Host>();

    const auto& qv_v = qv_f.get_view<Real**>();
    const auto& qc_v = qc_f.get_view<Real**>();
    const auto& qi_v = qi_f.get_view<Real**>();
    const auto& qm_v = qm_f.get_view<Real**>();
    const auto& qr_v = qr_f.get_view<Real**>();
    const auto& dp_v = dp_f.get_view<Real**>();


    for (const auto& dd : diags) {
      dd.second->compute(t0);
    }
    // Test 1: The water path should be >= the maximum cell level mass per column
    {
      auto policy = TPF::get_default_team_policy(ncols, nlevs);
      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        auto dp_icol = ekat::subview(dp_v, icol);

        // Test qv mass
        Real qv_mass_max=0.0;
        auto qv_icol = ekat::subview(qv_v, icol);
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,nlevs), [&] (Int idx, Real& lmax) {
            if (qv_icol(idx)*dp_icol(idx)/gravit > lmax) {
              lmax = qv_icol(idx)*dp_icol(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qv_mass_max));
        EKAT_KERNEL_REQUIRE(vwp_v(icol)>qv_mass_max);
        team.team_barrier();
        // Test qc mass
        Real qc_mass_max=0.0;
        auto qc_icol = ekat::subview(qc_v, icol);
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,nlevs), [&] (Int idx, Real& lmax) {
            if (qc_icol(idx)*dp_icol(idx)/gravit > lmax) {
              lmax = qc_icol(idx)*dp_icol(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qc_mass_max));
        EKAT_KERNEL_REQUIRE(lwp_v(icol)>qc_mass_max);
        team.team_barrier();
        // Test qi mass
        Real qi_mass_max=0.0;
        auto qi_icol = ekat::subview(qi_v, icol);
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,nlevs), [&] (Int idx, Real& lmax) {
            if (qi_icol(idx)*dp_icol(idx)/gravit > lmax) {
              lmax = qi_icol(idx)*dp_icol(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qi_mass_max));
        EKAT_KERNEL_REQUIRE(iwp_v(icol)>qi_mass_max);
        team.team_barrier();
        // Test qm mass
        Real qm_mass_max=0.0;
        auto qm_icol = ekat::subview(qm_v, icol);
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,nlevs), [&] (Int idx, Real& lmax) {
            if (qm_icol(idx)*dp_icol(idx)/gravit > lmax) {
              lmax = qm_icol(idx)*dp_icol(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qm_mass_max));
        EKAT_KERNEL_REQUIRE(mwp_v(icol)>qm_mass_max);
        team.team_barrier();
        // Test qr mass
        Real qr_mass_max=0.0;
        auto qr_icol = ekat::subview(qr_v, icol);
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,nlevs), [&] (Int idx, Real& lmax) {
            if (qr_icol(idx)*dp_icol(idx)/gravit > lmax) {
              lmax = qr_icol(idx)*dp_icol(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qr_mass_max));
        EKAT_KERNEL_REQUIRE(rwp_v(icol)>qr_mass_max);
      });
      Kokkos::fence();
    }
    // Test 2: If the cell-wise mass is scaled by constant alpha then the water
    //         path should also be scaled by alpha.
    {
      Field vwp_copy_f = diags["vwp"]->get().clone();
      Field lwp_copy_f = diags["lwp"]->get().clone();
      Field iwp_copy_f = diags["iwp"]->get().clone();
      Field mwp_copy_f = diags["mwp"]->get().clone();
      Field rwp_copy_f = diags["rwp"]->get().clone();
      const auto& vwp_copy_v = vwp_copy_f.get_view<Real*>();
      const auto& lwp_copy_v = lwp_copy_f.get_view<Real*>();
      const auto& iwp_copy_v = iwp_copy_f.get_view<Real*>();
      const auto& mwp_copy_v = mwp_copy_f.get_view<Real*>();
      const auto& rwp_copy_v = rwp_copy_f.get_view<Real*>();

      Real alpha_qv = pdf_alpha(engine);
      Real alpha_qc = pdf_alpha(engine);
      Real alpha_qi = pdf_alpha(engine);
      Real alpha_qm = pdf_alpha(engine);
      Real alpha_qr = pdf_alpha(engine);
      REQUIRE(alpha_qv*alpha_qc*alpha_qi*alpha_qr != 1.0);

      qv_f.scale(alpha_qv);
      qc_f.scale(alpha_qc);
      qi_f.scale(alpha_qi);
      qm_f.scale(alpha_qm);
      qr_f.scale(alpha_qr);

      vwp_copy_f.scale(alpha_qv);
      lwp_copy_f.scale(alpha_qc);
      iwp_copy_f.scale(alpha_qi);
      mwp_copy_f.scale(alpha_qm);
      rwp_copy_f.scale(alpha_qr);

      // Change inputs timestamp, to prevent early return and trigger diag recalculation
      for (auto f : {qv_f, qi_f, qc_f, qr_f, qm_f} ) {
        f.get_header().get_tracking().update_time_stamp(t0+1);
      }
      for (const auto& dd : diags) {
        dd.second->compute(t0+1);
      }
      vwp_copy_f.sync_to_host();
      lwp_copy_f.sync_to_host();
      iwp_copy_f.sync_to_host();
      mwp_copy_f.sync_to_host();
      rwp_copy_f.sync_to_host();
      const auto& vwp_copy_h = vwp_copy_f.get_view<Real*,Host>();
      const auto& lwp_copy_h = lwp_copy_f.get_view<Real*,Host>();
      const auto& iwp_copy_h = iwp_copy_f.get_view<Real*,Host>();
      const auto& mwp_copy_h = mwp_copy_f.get_view<Real*,Host>();
      const auto& rwp_copy_h = rwp_copy_f.get_view<Real*,Host>();
      vwp.sync_to_host();
      lwp.sync_to_host();
      iwp.sync_to_host();
      mwp.sync_to_host();
      rwp.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(std::abs(vwp_copy_h(icol)-vwp_h(icol))<macheps);
        REQUIRE(std::abs(lwp_copy_h(icol)-lwp_h(icol))<macheps);
        REQUIRE(std::abs(iwp_copy_h(icol)-iwp_h(icol))<macheps);
        REQUIRE(std::abs(mwp_copy_h(icol)-mwp_h(icol))<macheps);
        REQUIRE(std::abs(rwp_copy_h(icol)-rwp_h(icol))<macheps);
      }
    }
    // Test 3: If mass moves from one phase to another than the total water path
    //         should remain unchanged.
    {
      auto total_mass_f = diags["vwp"]->get().clone();
      auto total_mass_v = total_mass_f.get_view<Real*>();
      const auto alpha_qv_to_qc = pdf_alpha(engine);
      const auto alpha_qc_to_qi = pdf_alpha(engine);
      const auto alpha_qi_to_qr = pdf_alpha(engine);
      const auto alpha_qr_to_qv = pdf_alpha(engine);
      Kokkos::parallel_for("",ncols*nlevs,
                           KOKKOS_LAMBDA(const int& idx) {
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        total_mass_v(icol) = vwp_v(icol) + lwp_v(icol) + iwp_v(icol) + rwp_v(icol);

        auto qv_to_qc = alpha_qv_to_qc*qv_v(icol,ilev);
        auto qc_to_qi = alpha_qc_to_qi*qc_v(icol,ilev);
        auto qi_to_qr = alpha_qi_to_qr*qi_v(icol,ilev);
        auto qr_to_qv = alpha_qr_to_qv*qr_v(icol,ilev);

        qv_v(icol,ilev) -= qv_to_qc;
        qc_v(icol,ilev) -= qc_to_qi;
        qi_v(icol,ilev) -= qi_to_qr;
        qr_v(icol,ilev) -= qr_to_qv;

        qv_v(icol,ilev) += qr_to_qv;
        qc_v(icol,ilev) += qv_to_qc;
        qi_v(icol,ilev) += qc_to_qi;
        qr_v(icol,ilev) += qi_to_qr;
      });
      Kokkos::fence();
      // Change inputs timestamp, to prevent early return and trigger diag recalculation
      for (auto f : {qv_f, qi_f, qc_f, qr_f, qm_f} ) {
        f.get_header().get_tracking().update_time_stamp(t0+2);
      }
      for (const auto& dd : diags) {
        dd.second->compute(t0+2);
      }
      total_mass_f.sync_to_host();
      auto total_mass_h = total_mass_f.get_view<Real*,Host>();
      vwp.sync_to_host();
      lwp.sync_to_host();
      iwp.sync_to_host();
      rwp.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        const auto new_total_mass = vwp_h(icol) + lwp_h(icol) + iwp_h(icol) + rwp_h(icol);
        REQUIRE(std::abs(total_mass_h(icol)-new_total_mass)<macheps);
      }
    }
    // Test 4: Delta_TWP = Delta_LWP + Delta_IWP + Delta_RWP
    {
      const auto alpha_qc_precip = pdf_alpha(engine);
      const auto alpha_qi_precip = pdf_alpha(engine);
      const auto alpha_qr_precip = pdf_alpha(engine);
      const int surf_lev = nlevs-1;
      auto total_mass_f = diags["vwp"]->get().clone();
      auto delta_mass_f = diags["vwp"]->get().clone();
      auto total_mass_v = total_mass_f.get_view<Real*>();
      auto delta_mass_v = delta_mass_f.get_view<Real*>();
      Kokkos::parallel_for("",ncols,KOKKOS_LAMBDA(const int& icol) {

        total_mass_v(icol) = vwp_v(icol) + lwp_v(icol) + iwp_v(icol) + rwp_v(icol);
        // simulate the loss of mass due to precipitation
        const auto& qc_sub      = ekat::subview(qc_v,icol);
        const auto& qi_sub      = ekat::subview(qi_v,icol);
        const auto& qr_sub      = ekat::subview(qr_v,icol);
        const auto& dp_sub      = ekat::subview(dp_v,icol);

        const Real dp_surf = dp_sub(surf_lev);

        const Real qc_precip = alpha_qc_precip * qc_sub(surf_lev);
        const Real qi_precip = alpha_qi_precip * qi_sub(surf_lev);
        const Real qr_precip = alpha_qr_precip * qr_sub(surf_lev);

        qc_sub(surf_lev) -= qc_precip;
        qi_sub(surf_lev) -= qi_precip;
        qr_sub(surf_lev) -= qr_precip;

        delta_mass_v(icol) = -(qc_precip + qi_precip + qr_precip) * dp_surf/gravit;
      });
      Kokkos::fence();
      // Change inputs timestamp, to prevent early return and trigger diag recalculation
      for (auto f : {qv_f, qi_f, qc_f, qr_f, qm_f} ) {
        f.get_header().get_tracking().update_time_stamp(t0+3);
      }
      for (const auto& dd : diags) {
        dd.second->compute(t0+3);
      }
      total_mass_f.sync_to_host();
      delta_mass_f.sync_to_host();
      auto total_mass_h = total_mass_f.get_view<Real*,Host>();
      auto delta_mass_h = delta_mass_f.get_view<Real*,Host>();
      vwp.sync_to_host();
      lwp.sync_to_host();
      iwp.sync_to_host();
      rwp.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        const auto new_total_mass = vwp_h(icol) + lwp_h(icol) + iwp_h(icol) + rwp_h(icol);
        const auto new_delta_mass = new_total_mass - total_mass_h(icol);
        REQUIRE(std::abs(delta_mass_h(icol)-new_delta_mass)<macheps);
      }
    }
    // Test 5: Verify water path calculation with pseudo_density=g and qx(k)=k+1
    //         X*sum(k=0,k=N-1)[k+1] = X*(N-1)*N/2
    {
      Kokkos::deep_copy(dp_v,gravit);
      Kokkos::parallel_for("",ncols*nlevs,
                           KOKKOS_LAMBDA(const int& idx) {
        const int icol  = idx / nlevs;
        const int ilev  = idx % nlevs;

        qv_v(icol,ilev) = (icol+1) * (idx+1);
        qc_v(icol,ilev) = (icol+1) * (idx+1);
        qi_v(icol,ilev) = (icol+1) * (idx+1);
        qm_v(icol,ilev) = (icol+1) * (idx+1);
        qr_v(icol,ilev) = (icol+1) * (idx+1);
      });
      Kokkos::fence();
      // Change inputs timestamp, to prevent early return and trigger diag recalculation
      for (auto f : {qv_f, qi_f, qc_f, qr_f, qm_f} ) {
        f.get_header().get_tracking().update_time_stamp(t0+4);
      }
      for (const auto& dd : diags) {
        dd.second->compute(t0+4);
      }
      vwp.sync_to_host();
      lwp.sync_to_host();
      iwp.sync_to_host();
      rwp.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(vwp_h(icol) == (icol+1)*nlevs*(nlevs+1)/2);
        REQUIRE(lwp_h(icol) == (icol+1)*nlevs*(nlevs+1)/2);
        REQUIRE(iwp_h(icol) == (icol+1)*nlevs*(nlevs+1)/2);
        REQUIRE(rwp_h(icol) == (icol+1)*nlevs*(nlevs+1)/2);
      }
    }
  }
}

TEST_CASE("water_path_test", "water_path_test]"){
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
}

} // namespace scream
