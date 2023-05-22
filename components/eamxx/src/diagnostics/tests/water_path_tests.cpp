#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/rain_water_path.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_local_elems = 4;
  const int np = 4;
  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_local_elements", num_local_elems);
  gm_params.set<int>("number_of_vertical_levels", nlevs);
  gm_params.set<int>("number_of_gauss_points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

template<typename VT>
typename VT::HostMirror
cmvdc (const VT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC         = scream::physics::Constants<Real>;
  using Pack       = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Pack>;
  using rview_1d   = typename KT::template view_1d<Real>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  const     int num_mid_packs    = ekat::npack<Pack>(num_levs);
  constexpr Real gravit = PC::gravit;
  constexpr Real macheps = PC::macheps;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d
    pseudo_density("pseudo_density",num_mid_packs),
    qv("qv",num_mid_packs),
    qc("qc",num_mid_packs),
    qr("qr",num_mid_packs),
    qi("qi",num_mid_packs),
    qm("qm",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF 
    pdf_qv(1e-12,1e-6),
    pdf_qx(0.0,1e-3),
    pdf_pseudodens(1.0,100.0),
    pdf_alpha(0.1,0.9);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  std::map<std::string,std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();
  ekat::ParameterList params;
  // Vapor
  auto diag_vap = diag_factory.create("VapWaterPath",comm,params);
  diag_vap->set_grids(gm);
  diags.emplace("vwp",diag_vap);
  // Liquid
  auto diag_liq = diag_factory.create("LiqWaterPath",comm,params);
  diag_liq->set_grids(gm);
  diags.emplace("lwp",diag_liq);
  // Ice
  auto diag_ice = diag_factory.create("IceWaterPath",comm,params);
  diag_ice->set_grids(gm);
  diags.emplace("iwp",diag_ice);
  // Rime
  auto diag_rime = diag_factory.create("RimeWaterPath",comm,params);
  diag_rime->set_grids(gm);
  diags.emplace("mwp",diag_rime);
  // Rain
  auto diag_rain = diag_factory.create("RainWaterPath",comm,params);
  diag_rain->set_grids(gm);
  diags.emplace("rwp",diag_rain);


  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& dd : diags) {
    const auto& diag = dd.second;
    for (const auto& req : diag->get_required_field_requests()) {
      if (input_fields.find(req.fid.name())==input_fields.end()) {
        Field f(req.fid);
        auto & f_ap = f.get_header().get_alloc_properties();
        f_ap.request_allocation(packsize);
        f.allocate_view();
        const auto name = f.name();
        f.get_header().get_tracking().update_time_stamp(t0);
        diag->set_required_field(f.get_const());
        REQUIRE_THROWS(diag->set_computed_field(f));
        input_fields.emplace(name,f);
      } else {
        auto& f = input_fields[req.fid.name()];
        const auto name = f.name();
        f.get_header().get_tracking().update_time_stamp(t0);
        diag->set_required_field(f.get_const());
        REQUIRE_THROWS(diag->set_computed_field(f));
        input_fields.emplace(name,f);
      }
    }
    // Initialize the diagnostic
    diag->initialize(t0,RunType::Initial);
  }


  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& qv_f          = input_fields["qv"];
    const auto& qv_v          = qv_f.get_view<Pack**>();
    const auto& qc_f          = input_fields["qc"];
    const auto& qc_v          = qc_f.get_view<Pack**>();
    const auto& qi_f          = input_fields["qi"];
    const auto& qi_v          = qi_f.get_view<Pack**>();
    const auto& qm_f          = input_fields["qm"];
    const auto& qm_v          = qm_f.get_view<Pack**>();
    const auto& qr_f          = input_fields["qr"];
    const auto& qr_v          = qr_f.get_view<Pack**>();
    const auto& pseudo_dens_f = input_fields["pseudo_density"];
    const auto& pseudo_dens_v = pseudo_dens_f.get_view<Pack**>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& qv_sub      = ekat::subview(qv_v,icol);
      const auto& qc_sub      = ekat::subview(qc_v,icol);
      const auto& qi_sub      = ekat::subview(qi_v,icol); 
      const auto& qm_sub      = ekat::subview(qm_v,icol);
      const auto& qr_sub      = ekat::subview(qr_v,icol);
      const auto& dp_sub      = ekat::subview(pseudo_dens_v,icol);
      ekat::genRandArray(dview_as_real(pseudo_density),  engine, pdf_pseudodens);
      Kokkos::deep_copy(dp_sub,pseudo_density);

      ekat::genRandArray(dview_as_real(qv),              engine, pdf_qv);
      ekat::genRandArray(dview_as_real(qc),              engine, pdf_qx);
      ekat::genRandArray(dview_as_real(qi),              engine, pdf_qx);
      ekat::genRandArray(dview_as_real(qm),              engine, pdf_qx);
      ekat::genRandArray(dview_as_real(qr),              engine, pdf_qx);
      Kokkos::deep_copy(qv_sub,qv);
      Kokkos::deep_copy(qc_sub,qc);
      Kokkos::deep_copy(qi_sub,qi);
      Kokkos::deep_copy(qm_sub,qm);
      Kokkos::deep_copy(qr_sub,qr);
    }
    // Grab views for each of the water path diagnostics
    const auto& vwp = diags["vwp"]->get_diagnostic();
    const auto& lwp = diags["lwp"]->get_diagnostic();
    const auto& iwp = diags["iwp"]->get_diagnostic();
    const auto& mwp = diags["mwp"]->get_diagnostic();
    const auto& rwp = diags["rwp"]->get_diagnostic();
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

    for (const auto& dd : diags) {
      dd.second->compute_diagnostic();
    }
    // Test 1: The water path should be >= the maximum cell level mass per column
    {
      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        auto dp_s = ekat::scalarize(ekat::subview(pseudo_dens_v, icol));

        // Test qv mass
        Real qv_mass_max=0.0;
        auto qv_s = ekat::scalarize(ekat::subview(qv_v, icol));
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,num_levs), [&] (Int idx, Real& lmax) {
            if (qv_s(idx)*dp_s(idx)/gravit > lmax) {
              lmax = qv_s(idx)*dp_s(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qv_mass_max));
        EKAT_KERNEL_REQUIRE(vwp_v(icol)>qv_mass_max);
        // Test qc mass
        Real qc_mass_max=0.0;
        auto qc_s = ekat::scalarize(ekat::subview(qc_v, icol));
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,num_levs), [&] (Int idx, Real& lmax) {
            if (qc_s(idx)*dp_s(idx)/gravit > lmax) {
              lmax = qc_s(idx)*dp_s(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qc_mass_max));
        EKAT_KERNEL_REQUIRE(lwp_v(icol)>qc_mass_max);
        // Test qi mass
        Real qi_mass_max=0.0;
        auto qi_s = ekat::scalarize(ekat::subview(qi_v, icol));
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,num_levs), [&] (Int idx, Real& lmax) {
            if (qi_s(idx)*dp_s(idx)/gravit > lmax) {
              lmax = qi_s(idx)*dp_s(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qi_mass_max));
        EKAT_KERNEL_REQUIRE(iwp_v(icol)>qi_mass_max);
        // Test qm mass
        Real qm_mass_max=0.0;
        auto qm_s = ekat::scalarize(ekat::subview(qm_v, icol));
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,num_levs), [&] (Int idx, Real& lmax) {
            if (qm_s(idx)*dp_s(idx)/gravit > lmax) {
              lmax = qm_s(idx)*dp_s(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qm_mass_max));
        EKAT_KERNEL_REQUIRE(mwp_v(icol)>qm_mass_max);

        // Test qr mass
        Real qr_mass_max=0.0;
        auto qr_s = ekat::scalarize(ekat::subview(qr_v, icol));
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team,num_levs), [&] (Int idx, Real& lmax) {
            if (qr_s(idx)*dp_s(idx)/gravit > lmax) {
              lmax = qr_s(idx)*dp_s(idx)/gravit;
            }
          }, Kokkos::Max<Real>(qr_mass_max));
        EKAT_KERNEL_REQUIRE(rwp_v(icol)>qr_mass_max);
      });
      Kokkos::fence();
    }
    // Test 2: If the cell-wise mass is scaled by constant alpha then the water
    //         path should also be scaled by alpha.
    {
      Field vwp_copy_f = diags["vwp"]->get_diagnostic().clone();
      Field lwp_copy_f = diags["lwp"]->get_diagnostic().clone();
      Field iwp_copy_f = diags["iwp"]->get_diagnostic().clone();
      Field mwp_copy_f = diags["mwp"]->get_diagnostic().clone();
      Field rwp_copy_f = diags["rwp"]->get_diagnostic().clone();
      const auto& vwp_copy_v = vwp_copy_f.get_view<Real*>();
      const auto& lwp_copy_v = lwp_copy_f.get_view<Real*>();
      const auto& iwp_copy_v = iwp_copy_f.get_view<Real*>();
      const auto& mwp_copy_v = mwp_copy_f.get_view<Real*>();
      const auto& rwp_copy_v = rwp_copy_f.get_view<Real*>();

      const auto alpha_qv = pdf_alpha(engine);
      const auto alpha_qc = pdf_alpha(engine);
      const auto alpha_qi = pdf_alpha(engine);
      const auto alpha_qm = pdf_alpha(engine);
      const auto alpha_qr = pdf_alpha(engine);
      REQUIRE(alpha_qv*alpha_qc*alpha_qi*alpha_qr != 1.0);

      Kokkos::parallel_for("",ncols*num_mid_packs,KOKKOS_LAMBDA(const int& idx) {
        const int icol  = idx / num_mid_packs;
        const int jpack = idx % num_mid_packs;

        qv_v(icol,jpack) *= alpha_qv;
        qc_v(icol,jpack) *= alpha_qc;
        qi_v(icol,jpack) *= alpha_qi;
        qm_v(icol,jpack) *= alpha_qm;
        qr_v(icol,jpack) *= alpha_qr;

        if (jpack==0) {
          vwp_copy_v(icol) *= alpha_qv;
          lwp_copy_v(icol) *= alpha_qc;
          iwp_copy_v(icol) *= alpha_qi;
          mwp_copy_v(icol) *= alpha_qm;
          rwp_copy_v(icol) *= alpha_qr;
        }
      });
      Kokkos::fence();
      for (const auto& dd : diags) {
        dd.second->compute_diagnostic();
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
      rview_1d total_mass("",ncols);
      const auto alpha_qv_to_qc = pdf_alpha(engine);
      const auto alpha_qc_to_qi = pdf_alpha(engine);
      const auto alpha_qi_to_qr = pdf_alpha(engine);
      const auto alpha_qr_to_qv = pdf_alpha(engine);
      Kokkos::parallel_for("",ncols*num_mid_packs,KOKKOS_LAMBDA(const int& idx) {
        const int icol  = idx / num_mid_packs;
        const int jpack = idx % num_mid_packs;

        total_mass(icol) = vwp_v(icol) + lwp_v(icol) + iwp_v(icol) + rwp_v(icol);

        auto qv_to_qc = alpha_qv_to_qc*qv_v(icol,jpack);
        auto qc_to_qi = alpha_qc_to_qi*qc_v(icol,jpack);
        auto qi_to_qr = alpha_qi_to_qr*qi_v(icol,jpack);
        auto qr_to_qv = alpha_qr_to_qv*qr_v(icol,jpack);

        qv_v(icol,jpack) -= qv_to_qc;
        qc_v(icol,jpack) -= qc_to_qi;
        qi_v(icol,jpack) -= qi_to_qr;
        qr_v(icol,jpack) -= qr_to_qv;

        qv_v(icol,jpack) += qr_to_qv;
        qc_v(icol,jpack) += qv_to_qc;
        qi_v(icol,jpack) += qc_to_qi;
        qr_v(icol,jpack) += qi_to_qr;
      });
      Kokkos::fence();
      for (const auto& dd : diags) {
        dd.second->compute_diagnostic();
      }
      auto total_mass_h = cmvdc(total_mass);
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
      const int surf_pack = num_mid_packs-1;
      const int surf_lev  = std::max(0,num_levs % Pack::n - 1);
      rview_1d total_mass("",ncols);
      rview_1d delta_mass("",ncols);
      Kokkos::parallel_for("",ncols,KOKKOS_LAMBDA(const int& icol) {

        total_mass(icol) = vwp_v(icol) + lwp_v(icol) + iwp_v(icol) + rwp_v(icol);
        // simulate the loss of mass due to precipitation
        const auto& qc_sub      = ekat::subview(qc_v,icol);
        const auto& qi_sub      = ekat::subview(qi_v,icol);
        const auto& qr_sub      = ekat::subview(qr_v,icol);
        const auto& dp_sub      = ekat::subview(pseudo_dens_v,icol);

        const Real dp_surf = dp_sub(surf_pack)[surf_lev];

        const Real qc_precip = alpha_qc_precip * qc_sub(surf_pack)[surf_lev];
        const Real qi_precip = alpha_qi_precip * qi_sub(surf_pack)[surf_lev];
        const Real qr_precip = alpha_qr_precip * qr_sub(surf_pack)[surf_lev];

        qc_sub(surf_pack)[surf_lev] -= qc_precip;
        qi_sub(surf_pack)[surf_lev] -= qi_precip;
        qr_sub(surf_pack)[surf_lev] -= qr_precip;

        delta_mass(icol) = -(qc_precip + qi_precip + qr_precip) * dp_surf/gravit;

      });
      Kokkos::fence();
      for (const auto& dd : diags) {
        dd.second->compute_diagnostic();
      }
      auto total_mass_h = cmvdc(total_mass);
      auto delta_mass_h = cmvdc(delta_mass);
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
      Kokkos::deep_copy(pseudo_dens_v,gravit);
      Kokkos::parallel_for("",ncols*num_levs,KOKKOS_LAMBDA(const int& idx) {
        const int icol  = idx / num_levs;
        const int ilev  = idx % num_levs;
        const int kpack = ilev / Pack::n;
        const int klev  = ilev % Pack::n;

        qv_v(icol,kpack)[klev] = (icol+1) * (idx+1);
        qc_v(icol,kpack)[klev] = (icol+1) * (idx+1);
        qi_v(icol,kpack)[klev] = (icol+1) * (idx+1);
        qm_v(icol,kpack)[klev] = (icol+1) * (idx+1);
        qr_v(icol,kpack)[klev] = (icol+1) * (idx+1);
      });
      Kokkos::fence();
      for (const auto& dd : diags) {
        dd.second->compute_diagnostic();
      }
      vwp.sync_to_host();
      lwp.sync_to_host();
      iwp.sync_to_host();
      rwp.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(vwp_h(icol) == (icol+1)*num_levs*(num_levs+1)/2);
        REQUIRE(lwp_h(icol) == (icol+1)*num_levs*(num_levs+1)/2);
        REQUIRE(iwp_h(icol) == (icol+1)*num_levs*(num_levs+1)/2);
        REQUIRE(rwp_h(icol) == (icol+1)*num_levs*(num_levs+1)/2);
      }
    }

  }
 
  // Finalize the diagnostic
  for (const auto& dd : diags) {
    const auto& diag = dd.second;
    diag->finalize();
  } 

} // run()

TEST_CASE("water_path_test", "water_path_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace
