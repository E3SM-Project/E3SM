#ifndef SCREAM_MAM4_AEROSOL_MICROPHYSICS_HPP
#define SCREAM_MAM4_AEROSOL_MICROPHYSICS_HPP

#include <share/atm_process/atmosphere_process.hpp>
//#include <share/util/scream_common_physics_functions.hpp>
//#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
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
class MAM4AerosolMicrophysics final : public scream::AtmosphereProcess
{
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;

  using ColumnView   = mam4::ColumnView;
  using ThreadTeam   = mam4::ThreadTeam;

public:

  // Constructor
  MAM4AerosolMicrophysics(const ekat::Comm& comm,
                          const ekat::ParameterList& params);

protected:

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
  void run_impl(const int dt) override;
  void finalize_impl() override;

  // MAM4xx updates the 'tracers' group.
  void set_computed_group_impl(const FieldGroup& group) override;

private:

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct MAM4Preprocess {
    MAM4Preprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&](const int k) {
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
        q_soag(i,k)  = PF::calculate_drymmr_from_wetmmr(q_soag(i,k), qv(i,k));
        q_h2so4(i,k) = PF::calculate_drymmr_from_wetmmr(q_h2so4(i,k), qv(i,k));
        q_nh3(i,k)   = PF::calculate_drymmr_from_wetmmr(q_nh3(i,k), qv(i,k));

        q_aitken_so4(i,k = PF::calculate_drymmr_from_wetmmr(q_aitken_so4(i,k), qv(i,k));

        // convert qv to dry mmr
        qv(i,k) = PF::calculate_drymmr_from_wetmmr(qv(i,k), qv(i,k));
      });
      team.team_barrier();
    } // operator

    int ncol, nlev;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d;

    // local atmospheric state column variables
    view_2d_const T_mid;  // temperature at grid midpoints [K]
    view_2d_const p_mid;  // total pressure at grid midpoints [Pa]
    view_2d       qv;     // water vapor mass mixing ratio, not const because it
                          // must be converted from wet to dry [kg vapor/kg dry air]
    view_2d       height; // height at grid interfaces [m]
    view_2d       pdel;   // hydrostatic "pressure thickness" at grid
                          // interfaces [Pa]
    view_1d_const pblh;   // planetary boundary layer height [m]

    // local aerosol-related gases
    view_2d       q_soag;  // secondary organic aerosol gas [kg gas/kg dry air]
    view_2d       q_h2so4; // H2SO3 gas [kg/kg dry air]
    view_2d       q_nh3;   // NH3 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d       q_aitken_so4; // SO4 aerosol in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol_, const int nlev_,
                       const view_2d_const& T_mid_,
                       const view_2d_const& p_mid_,
                       const view_2d&       qv_,
                       const view_2d&       height_,
                       const view_2d&       pdel_,
                       const view_1d_const& pblh_,
                       const view_2d&       q_soag_,
                       const view_2d&       q_h2so4,
                       const view_2d&       q_nh3,
                       const view_2d&       q_aitken_so4) {
      ncol = ncol_;
      nlev = nlev_;
      T_mid = T_mid_;
      p_mid = p_mid_;
      qv = qv_;
      height = height_;
      pdel = pdel_;
      pblh = pblh_;
      q_soag = q_soag_;
      q_h2so4 = q_h2so4_;
      q_nh3 = q_nh3_;
      q_aitken_so4 = q_aitken_so4_;
    } // set_variables
  }; // MAM4Preprocess

  // Postprocessing functor
  struct MAM4Postprocess {
    MAM4Postprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      // Copy the updated aerosol tracers back into the tracer array.
      // FIXME

      // After these updates, all tracers are converted from dry mmr to wet mmr
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&](const int k) {
        // Here we convert our dry mmrs back to wet mmrs for EAMxx.
        // NOTE: calculate_wetmmr_from_drymmr takes 2 arguments:
        // 1. dry mmr
        // 2. "dry" water vapor mixing ratio
        q_soag(i,k)  = PF::calculate_wetmmr_from_drymmr(q_soag(i,k), qv(i,k));
        q_h2so4(i,k) = PF::calculate_wetmmr_from_drymmr(q_h2so4(i,k), qv(i,k));
        q_nh3(i,k)   = PF::calculate_wetmmr_from_drymmr(q_nh3(i,k), qv(i,k));

        q_aitken_so4(i,k = PF::calculate_wetmmr_from_drymmr(q_aitken_so4(i,k), qv(i,k));

        qv(i,k) = PF::calculate_wetmmr_from_drymmr(qv(i,k), qv(i,k));
      });
      team.team_barrier();
    } // operator

    // Local variables
    int ncol, nlev;
    view_2d qv;

    // local aerosol-related gases
    view_2d       q_soag;  // secondary organic aerosol gas [kg gas/kg dry air]
    view_2d       q_h2so4; // H2SO3 gas [kg/kg dry air]
    view_2d       q_nh3;   // NH3 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d       q_aitken_so4; // SO4 aerosol in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol_, const int nlev_,
                       const view_2d& qv_,
                       const view_2d& q_soag_,
                       const view_2d& q_h2so4,
                       const view_2d& q_nh3,
                       const view_2d& q_aitken_so4) {
    {
      ncol = ncol_;
      nlev = nlev_;
      qv = qv_;
      q_soag = q_soag_;
      q_h2so4 = q_h2so4_;
      q_nh3 = q_nh3_;
      q_aitken_so4 = q_aitken_so4_;
    } // set_variables
  }; // MAM4Postprocess

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_1d_scalar_ncol = 18;
    static constexpr int num_1d_scalar_nlev = 1;
    static constexpr int num_2d_vector_mid  = 22;
    static constexpr int num_2d_vector_int  = 13;
    static constexpr int num_2d_vector_tr   = 1;

    uview_1d<Real> cell_length;
    uview_1d<Real> wpthlp_sfc;
    uview_1d<Real> wprtp_sfc;
    uview_1d<Real> upwp_sfc;
    uview_1d<Real> vpwp_sfc;

    uview_1d<Spack> pref_mid;

    uview_2d<Spack> z_mid;
    uview_2d<Spack> z_int;
    uview_2d<Spack> rrho;
    uview_2d<Spack> rrho_i;
    uview_2d<Spack> thv;
    uview_2d<Spack> dz;
    uview_2d<Spack> zt_grid;
    uview_2d<Spack> zi_grid;
    uview_2d<Spack> wtracer_sfc;
    uview_2d<Spack> wm_zt;

    Spack* wsm_data;
  };
  /* --------------------------------------------------------------------------------------------*/

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // Aerosol processes
  mam4::NucleationProcess nucleation_;

  // WSM for internal local variables
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
};

} // namespace scream

#endif // SCREAM_MAM4_AEROSOLS_HPP
