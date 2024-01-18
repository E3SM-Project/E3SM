#include <physics/mam/eamxx_mam_optics_process_interface.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h" // for SCREAM_CIME_BUILD

#include <ekat/ekat_assert.hpp>

namespace scream
{

MAMOptics::MAMOptics(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params),
    aero_config_() {
}

AtmosphereProcessType MAMOptics::type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAMOptics::name() const {
  return "mam4_optics";
}

void MAMOptics::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // number of levels per column
  nswbands_ = 14;                           // number of shortwave bands
  nlwbands_ = 16;                           // number of longwave bands

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Define aerosol optics fields computed by this process.
  auto nondim = Units::nondimensional();
  FieldLayout scalar3d_swband_layout { {COL, SWBND, LEV}, {ncol_, nswbands_, nlev_} };
  FieldLayout scalar3d_lwband_layout { {COL, LWBND, LEV}, {ncol_, nlwbands_, nlev_} };

  // shortwave aerosol scattering asymmetry parameter [-]
  add_field<Computed>("aero_g_sw",   scalar3d_swband_layout, nondim, grid_name);
  // shortwave aerosol single-scattering albedo [-]
  add_field<Computed>("aero_ssa_sw", scalar3d_swband_layout, nondim, grid_name);
  // shortwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_sw", scalar3d_swband_layout, nondim, grid_name);
  // longwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_lw", scalar3d_lwband_layout, nondim, grid_name);

  // FIXME: this field doesn't belong here, but this is a convenient place to
  // FIXME: put it for now.
  // number mixing ratio for CCN
  using Spack      = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;
  using Pack       = ekat::Pack<Real,Spack::n>;
  constexpr int ps = Pack::n;
  FieldLayout scalar3d_layout_mid { {COL, LEV}, {ncol_, nlev_} };
  add_field<Computed>("nccn", scalar3d_layout_mid, 1/kg, grid_name, ps);
}

void MAMOptics::initialize_impl(const RunType run_type) {
}

void MAMOptics::run_impl(const double dt) {

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // get the aerosol optics fields
  auto aero_g_sw   = get_field_out("aero_g_sw").get_view<Real***>();
  auto aero_ssa_sw = get_field_out("aero_ssa_sw").get_view<Real***>();
  auto aero_tau_sw = get_field_out("aero_tau_sw").get_view<Real***>();
  auto aero_tau_lw = get_field_out("aero_tau_lw").get_view<Real***>();

  auto aero_nccn   = get_field_out("nccn").get_view<Real**>(); // FIXME: get rid of this

  // Compute optical properties on all local columns.
  // (Strictly speaking, we don't need this parallel_for here yet, but we leave
  //  it in anticipation of column-specific aerosol optics to come.)
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    const Int icol = team.league_rank(); // column index

    auto g_sw = ekat::subview(aero_g_sw, icol);
    auto ssa_sw = ekat::subview(aero_ssa_sw, icol);
    auto tau_sw = ekat::subview(aero_tau_sw, icol);
    auto tau_lw = ekat::subview(aero_tau_lw, icol);

    // populate these fields with reasonable representative values
    Kokkos::deep_copy(g_sw, 0.5);
    Kokkos::deep_copy(ssa_sw, 0.7);
    Kokkos::deep_copy(tau_sw, 0.0);
    Kokkos::deep_copy(tau_lw, 0.0);

    // FIXME: Get rid of this
    auto nccn = ekat::subview(aero_nccn, icol);
    Kokkos::deep_copy(nccn, 50.0);
  });
}

void MAMOptics::finalize_impl()
{
}

} // namespace scream
