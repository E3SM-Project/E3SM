#ifndef EAMXX_MAM_MICROPHYSICS_HPP
#define EAMXX_MAM_MICROPHYSICS_HPP

#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>
#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected public
#define private public
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosols. The AD stores exactly ONE
// instance of this class in its list of subcomponents.
class MAMMicrophysics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // views for single- and multi-column data
  using view_1d_int   = typename KT::template view_1d<int>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_1d_const = typename KT::template view_1d<const Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using view_2d_const = typename KT::template view_2d<const Real>;

  // unmanaged views (for buffer and workspace manager)
  using uview_1d = Unmanaged<typename KT::template view_1d<Real>>;
  using uview_2d = Unmanaged<typename KT::template view_2d<Real>>;

  using ColumnView = mam4::ColumnView;
  using ThreadTeam = mam4::ThreadTeam;

public:

  // Constructor
  MAMMicrophysics(const ekat::Comm& comm, const ekat::ParameterList& params);

protected:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  /*
  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;
  */

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // MAM4xx updates the 'tracers' group.
  void set_computed_group_impl(const FieldGroup& group) override;

private:

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_), [&](const int k) {
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
        // conversion. qv is converted to dry mmr in the next parallel for.
        q_soag_(i,k)  = PF::calculate_drymmr_from_wetmmr(q_soag_(i,k), qv_(i,k));
        q_h2so4_(i,k) = PF::calculate_drymmr_from_wetmmr(q_h2so4_(i,k), qv_(i,k));
        q_nh3_(i,k)   = PF::calculate_drymmr_from_wetmmr(q_nh3_(i,k), qv_(i,k));

        q_aitken_so4_(i,k) = PF::calculate_drymmr_from_wetmmr(q_aitken_so4_(i,k), qv_(i,k));

        // convert qv to dry mmr
        qv_(i,k) = PF::calculate_drymmr_from_wetmmr(qv_(i,k), qv_(i,k));
      });
      team.team_barrier();
    } // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d_;

    // local atmospheric state column variables
    view_2d_const T_mid_;  // temperature at grid midpoints [K]
    view_2d_const p_mid_;  // total pressure at grid midpoints [Pa]
    view_2d       qv_;     // water vapor mass mixing ratio, not const because it
                           // must be converted from wet to dry [kg vapor/kg dry air]
    view_2d       height_; // height at grid interfaces [m]
    view_2d       pdel_;   // hydrostatic "pressure thickness" at grid
                           // interfaces [Pa]
    view_1d_const pblh_;   // planetary boundary layer height [m]

    // local aerosol-related gases
    view_2d       q_soag_;  // secondary organic aerosol gas [kg gas/kg dry air]
    view_2d       q_h2so4_; // H2SO3 gas [kg/kg dry air]
    view_2d       q_nh3_;   // NH3 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d       q_aitken_so4_; // SO4 aerosol in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol, const int nlev,
                       const view_2d_const& T_mid,
                       const view_2d_const& p_mid,
                       const view_2d&       qv,
                       const view_2d&       height,
                       const view_2d&       pdel,
                       const view_1d_const& pblh,
                       const view_2d&       q_soag,
                       const view_2d&       q_h2so4,
                       const view_2d&       q_nh3,
                       const view_2d&       q_aitken_so4) {
      ncol_ = ncol;
      nlev_ = nlev;
      T_mid_ = T_mid;
      p_mid_ = p_mid;
      qv_ = qv;
      height_ = height;
      pdel_ = pdel;
      pblh_ = pblh;
      q_soag_ = q_soag;
      q_h2so4_ = q_h2so4;
      q_nh3_ = q_nh3;
      q_aitken_so4_ = q_aitken_so4;
    } // set_variables
  }; // MAM4AerosolMicrophysics::Preprocess

  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      // Copy the updated aerosol tracers back into the tracer array.
      // FIXME

      // After these updates, all tracers are converted from dry mmr to wet mmr
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_), [&](const int k) {
        // Here we convert our dry mmrs back to wet mmrs for EAMxx.
        // NOTE: calculate_wetmmr_from_drymmr takes 2 arguments:
        // 1. dry mmr
        // 2. "dry" water vapor mixing ratio
        q_soag_(i,k)  = PF::calculate_wetmmr_from_drymmr(q_soag_(i,k), qv_(i,k));
        q_h2so4_(i,k) = PF::calculate_wetmmr_from_drymmr(q_h2so4_(i,k), qv_(i,k));
        q_nh3_(i,k)   = PF::calculate_wetmmr_from_drymmr(q_nh3_(i,k), qv_(i,k));

        q_aitken_so4_(i,k) = PF::calculate_wetmmr_from_drymmr(q_aitken_so4_(i,k), qv_(i,k));

        qv_(i,k) = PF::calculate_wetmmr_from_drymmr(qv_(i,k), qv_(i,k));
      });
      team.team_barrier();
    } // operator()

    // Local variables
    int ncol_, nlev_;
    view_2d qv_;

    // local aerosol-related gases
    view_2d       q_soag_;  // secondary organic aerosol gas [kg gas/kg dry air]
    view_2d       q_h2so4_; // H2SO3 gas [kg/kg dry air]
    view_2d       q_nh3_;   // NH3 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d       q_aitken_so4_; // SO4 aerosol in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol,
                       const int nlev,
                       const view_2d& qv,
                       const view_2d& q_soag,
                       const view_2d& q_h2so4,
                       const view_2d& q_nh3,
                       const view_2d& q_aitken_so4) {
      ncol_ = ncol;
      nlev_ = nlev;
      qv_ = qv;
      q_soag_ = q_soag;
      q_h2so4_ = q_h2so4;
      q_nh3_ = q_nh3;
      q_aitken_so4_ = q_aitken_so4;
    } // set_variables
  }; // MAM4AerosolMicrophysics::Postprocess

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // aerosol processes
  std::unique_ptr<mam4::NucleationProcess> nucleation_;

  // pre- and postprocessing scratch pads
  Preprocess preprocess_;
  Postprocess postprocess_;

  /*
  // WSM for internal local variables
  ekat::WorkspaceManager<Real, KT::Device> workspace_mgr_;
  Buffer buffer_;
  */

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMMicrophysics

} // namespace scream

#endif // EAMXX_MAM_MICROPHYSICS_HPP
