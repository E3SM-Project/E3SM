#ifndef EAMXX_MAM_MICROPHYSICS_HPP
#define EAMXX_MAM_MICROPHYSICS_HPP

#include <physics/mam/mam_coupling.hpp>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>
#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosol microphysics. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMMicrophysics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // views for single- and multi-column data
  using view_1d_int   = typename KT::template view_1d<int>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  // unmanaged views (for buffer and workspace manager)
  using uview_1d = Unmanaged<typename KT::template view_1d<Real>>;
  using uview_2d = Unmanaged<typename KT::template view_2d<Real>>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

public:

  // Constructor
  MAMMicrophysics(const ekat::Comm& comm, const ekat::ParameterList& params);

protected_except_cuda:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // performs some checks on the tracers group
  void set_computed_group_impl(const FieldGroup& group) override;

private_except_cuda:

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;

    // on host: initializes preprocess functor with necessary state data
    void initialize(const int ncol, const int nlev, const Real z_surf,
                    const view_1d_int& convert_wet_dry_idx_d,
                    const mam_coupling::AtmosphericState& atm,
                    const mam_coupling::AerosolState& aero) {
      ncol_ = ncol;
      nlev_ = nlev;
      z_surf_ = z_surf;
      convert_wet_dry_idx_d_ = convert_wet_dry_idx_d;
      atm_ = atm;
      aero_ = aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank(); // column index

      // compute vertical layer heights
      const auto dz_i = ekat::subview(atm_.dz,    i);
      auto z_iface_i  = ekat::subview(atm_.z_iface, i);
      auto z_mid_i    = ekat::subview(atm_.z_mid, i);
      PF::calculate_z_int(team, nlev_, dz_i, z_surf_, z_iface_i);
      team.team_barrier();  // TODO: is this barrier necessary?
      PF::calculate_z_mid(team, nlev_, z_iface_i, z_mid_i);
      // barrier here allows the kernels that follow to use layer heights
      team.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_),
        [&] (const int k) {
          //--------------------------
          // Vertical velocity from pressure to height
          //--------------------------
          const auto rho = PF::calculate_density(atm_.pdel(i,k), dz_i(k));
          atm_.w_updraft(i,k) = PF::calculate_vertical_velocity(atm_.omega(i,k), rho);

          //--------------------------
          // Wet to dry mixing ratios
          //--------------------------
          //
          // Since tracers from the host model (or AD) are wet mixing ratios, and
          // MAM4 expects these tracers in dry mixing ratios, we convert the wet
          // mixing ratios to dry mixing ratios for all the tracers.
          //
          // The function calculate_drymmr_from_wetmmr takes 2 arguments:
          // 1. wet mmr
          // 2. "wet" water vapor mixing ratio
          //
          // Units of all tracers become [kg/kg(dry-air)] for mass mixing ratios and
          // [#/kg(dry-air)] for number mixing ratios after the following
          // conversion.
          const auto qv_ik = atm_.qv_wet(i,k);

          // calculate dry atmospheric mixing ratios
          atm_.qv_dry(i,k) = PF::calculate_drymmr_from_wetmmr(atm_.qv_wet(i,k), qv_ik);
          atm_.qc_dry(i,k) = PF::calculate_drymmr_from_wetmmr(atm_.qc_wet(i,k), qv_ik);
          atm_.nc_dry(i,k) = PF::calculate_drymmr_from_wetmmr(atm_.nc_wet(i,k), qv_ik);
          atm_.qi_dry(i,k) = PF::calculate_drymmr_from_wetmmr(atm_.qi_wet(i,k), qv_ik);
          atm_.ni_dry(i,k) = PF::calculate_drymmr_from_wetmmr(atm_.ni_wet(i,k), qv_ik);

          // calculate dry aerosol mixing ratios
          for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
            aero_.dry_int_aero_nmr[m](i,k) = PF::calculate_drymmr_from_wetmmr(aero_.wet_int_aero_nmr[m](i,k), qv_ik);
            aero_.dry_cld_aero_nmr[m](i,k) = PF::calculate_drymmr_from_wetmmr(aero_.wet_cld_aero_nmr[m](i,k), qv_ik);
            for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
              if (aero_.dry_int_aero_mmr[m][a].data()) {
                aero_.dry_int_aero_mmr[m][a](i,k) = PF::calculate_drymmr_from_wetmmr(aero_.wet_int_aero_mmr[m][a](i,k), qv_ik);
              }
              if (aero_.dry_cld_aero_mmr[m][a].data()) {
                aero_.dry_cld_aero_mmr[m][a](i,k) = PF::calculate_drymmr_from_wetmmr(aero_.wet_cld_aero_mmr[m][a](i,k), qv_ik);
              }
            }
          }
          for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
            aero_.dry_gas_mmr[g](i,k) = PF::calculate_drymmr_from_wetmmr(aero_.dry_gas_mmr[g](i,k), qv_ik);
          }
        });
      // make sure all operations are done before exiting kernel
      team.team_barrier();
    } // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // height of bottom of atmosphere
    Real z_surf_;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d_;

    // local atmospheric and aerosol state column variables
    mam_coupling::AtmosphericState atm_;
    mam_coupling::AerosolState     aero_;

  }; // MAMMicrophysics::Preprocess

  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    // on host: initializes postprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const view_1d_int& convert_wet_dry_idx_d,
                    const mam_coupling::AtmosphericState& atm,
                    const mam_coupling::AerosolState& aero) {
      ncol_ = ncol;
      nlev_ = nlev;
      convert_wet_dry_idx_d_ = convert_wet_dry_idx_d;
      atm_ = atm;
      aero_ = aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      // after these updates, all non-const tracers are converted from dry mmr to wet mmr
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_),
        [&] (const int k) {
          const auto qv_ik = atm_.qv_dry(i,k);
          for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
            aero_.wet_int_aero_nmr[m](i,k) = PF::calculate_wetmmr_from_drymmr(aero_.dry_int_aero_nmr[m](i,k), qv_ik);
            aero_.wet_cld_aero_nmr[m](i,k) = PF::calculate_wetmmr_from_drymmr(aero_.dry_cld_aero_nmr[m](i,k), qv_ik);
            for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
              if (aero_.wet_int_aero_mmr[m][a].data()) {
                aero_.wet_int_aero_mmr[m][a](i,k) = PF::calculate_wetmmr_from_drymmr(aero_.dry_int_aero_mmr[m][a](i,k), qv_ik);
              }
            }
          }
          for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
            aero_.wet_gas_mmr[g](i,k) = PF::calculate_wetmmr_from_drymmr(aero_.wet_gas_mmr[g](i,k), qv_ik);
          }
      });
      team.team_barrier();
    } // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d_;

    // local atmospheric and aerosol state column variables
    mam_coupling::AtmosphericState atm_;
    mam_coupling::AerosolState     aero_;

  }; // MAMMicrophysics::Postprocess

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // aerosol processes
  std::unique_ptr<mam4::NucleationProcess> nucleation_;

  // pre- and postprocessing scratch pads
  Preprocess preprocess_;
  Postprocess postprocess_;

  // atmospheric and aerosol state variables
  mam_coupling::AtmosphericState atm_;
  mam_coupling::AerosolState     aero_;

  // workspace manager for internal local variables
  //ekat::WorkspaceManager<Real, KT::Device> workspace_mgr_;
  mam_coupling::Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMMicrophysics

} // namespace scream

#endif // EAMXX_MAM_MICROPHYSICS_HPP
