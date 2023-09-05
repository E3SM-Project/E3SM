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
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
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

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank(); // column index

      // Compute vertical layer heights
      const auto dz_i = ekat::subview(dz_,    i);
      auto z_iface_i = ekat::subview(z_iface_, i);
      auto z_mid_i   = ekat::subview(z_mid_, i);
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
          const auto rho = PF::calculate_density(pdel_(i,k), dz_i(k));
          w_updraft_(i,k) = PF::calculate_vertical_velocity(omega_(i,k), rho);

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
          const auto qv_ik = qv_(i,k);
          // const fields need separate storage for "dry" values
          qv_dry_(i,k) = PF::calculate_drymmr_from_wetmmr(qv_(i,k), qv_ik);
          qc_dry_(i,k) = PF::calculate_drymmr_from_wetmmr(qc_(i,k), qv_ik);
          n_qc_dry_(i,k) = PF::calculate_drymmr_from_wetmmr(n_qc_(i,k), qv_ik);
          qi_dry_(i,k) = PF::calculate_drymmr_from_wetmmr(qi_(i,k), qv_ik);
          n_qi_dry_(i,k) = PF::calculate_drymmr_from_wetmmr(n_qi_(i,k), qv_ik);

          // non-const fields can be overwritten; we'll convert back to moist
          // air ratios during postprocess
          q_h2so4_(i,k) = PF::calculate_drymmr_from_wetmmr(q_h2so4_(i,k), qv_ik);
          q_aitken_so4_(i,k) = PF::calculate_drymmr_from_wetmmr(q_aitken_so4_(i,k), qv_ik);
          n_aitken_(i,k) = PF::calculate_drymmr_from_wetmmr(n_aitken_(i,k), qv_ik);
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

    // local atmospheric state column variables
    const_view_2d T_mid_;   // temperature at grid midpoints [K]
    const_view_2d p_mid_;   // total pressure at grid midpoints [Pa]
    const_view_2d qv_;      // water vapor specific humidity [kg vapor / kg moist air]
    view_2d qv_dry_;  // water vapor mixing ratio [kg vapor / kg dry air]
    const_view_2d qc_;      // cloud liquid water mass mixing ratio [kg vapor/kg moist air]
    view_2d qc_dry_;
    const_view_2d n_qc_;      // cloud liquid water number mixing ratio [kg cloud water / kg moist air]
    view_2d n_qc_dry_; // cloud liquid water number mixing ratio (dry air)
    const_view_2d qi_;      // cloud ice water mass mixing ratio
    view_2d qi_dry_;     // [kg vapor/kg dry air]
    const_view_2d n_qi_;      // cloud ice water number mixing ratio
    view_2d n_qi_dry_;
    view_2d z_mid_;   // height at layer midpoints [m]
    view_2d z_iface_; // height at layer interfaces [m]
    view_2d dz_;      // layer thickness [m]
    const_view_2d pdel_;    // hydrostatic "pressure thickness" at grid
                      // interfaces [Pa]
    const_view_2d cloud_f_; // cloud fraction [-]
    const_view_2d omega_; // vertical pressure velocity [Pa/s]
    view_2d w_updraft_;  // updraft velocity [m/s]
    const_view_1d pblh_;    // planetary boundary layer height [m]

    // local aerosol-related gases
    view_2d q_h2so4_; // H2SO4 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d n_aitken_; // aitken mode number mixing ratio [1/kg dry air]
    view_2d q_aitken_so4_; // SO4 mass mixing ratio in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol, const int nlev, const Real z_surf,
                       const view_1d_int& convert_wet_dry_idx_d,
                       const const_view_2d&     T_mid,
                       const const_view_2d&     p_mid,
                       const const_view_2d&     qv,
                       const view_2d&     qv_dry,
                       const view_2d&     qc,
                       const view_2d&     n_qc,
                       const view_2d&     qc_dry,
                       const view_2d&     n_qc_dry,
                       const const_view_2d&     qi,
                       const const_view_2d&     n_qi,
                       const view_2d&     qi_dry,
                       const view_2d&     n_qi_dry,
                       const view_2d&     z_mid,
                       const view_2d&     z_iface,
                       const view_2d&     dz,
                       const const_view_2d&     pdel,
                       const const_view_2d&     cf,
                       const const_view_2d&     omega,
                       const view_2d&           w_updraft,
                       const const_view_1d&     pblh,
                       const view_2d&     q_h2so4,
                       const view_2d&     q_aitken_so4,
                       const view_2d&     n_aitken) {
      ncol_ = ncol;
      nlev_ = nlev;
      z_surf_ = z_surf;
      convert_wet_dry_idx_d_ = convert_wet_dry_idx_d;
      T_mid_ = T_mid;
      p_mid_ = p_mid;
      qv_ = qv;
      qv_dry_ = qv_dry;
      qc_ = qc;
      n_qc_ = n_qc;
      qc_dry_ = qc_dry;
      n_qc_dry_ = n_qc_dry;
      qi_ = qi;
      n_qi_ = n_qi;
      qi_dry_ = qi_dry;
      n_qi_dry_ = n_qi_dry;
      z_mid_ = z_mid;
      z_iface_ = z_iface;
      dz_ = dz;
      pdel_ = pdel;
      cloud_f_ = cf;
      omega_ = omega;
      w_updraft_ = w_updraft;
      pblh_ = pblh;
      q_h2so4_ = q_h2so4;
      q_aitken_so4_ = q_aitken_so4;
      n_aitken_ = n_aitken;
    } // set_variables
  }; // MAMMicrophysics::Preprocess

  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      // After these updates, all non-const tracers are converted from dry mmr to wet mmr
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_),
        [&] (const int k) {
          const auto qv_ik = qv_dry_(i,k);
          q_h2so4_(i,k) = PF::calculate_wetmmr_from_drymmr(q_h2so4_(i,k), qv_ik);
          q_aitken_so4_(i,k) = PF::calculate_wetmmr_from_drymmr(q_aitken_so4_(i,k), qv_ik);
          n_aitken_(i,k) = PF::calculate_wetmmr_from_drymmr(n_aitken_(i,k), qv_ik);
      });
      team.team_barrier();
    } // operator()

    // Local variables
    int ncol_, nlev_;
    view_2d qv_dry_;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d_;

    // local aerosol-related gases
    view_2d q_h2so4_; // H2SO4 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d q_aitken_so4_; // SO4 aerosol in aitken mode [kg/kg dry air]

    // modal quantities
    view_2d n_aitken_;

    // assigns local variables
    void set_variables(const int ncol,
                       const int nlev,
                       const view_1d_int& convert_wet_dry_idx_d,
                       const view_2d& qv_dry,
                       const view_2d& q_h2so4,
                       const view_2d& q_aitken_so4,
                       const view_2d& n_aitken) {
      ncol_ = ncol;
      nlev_ = nlev;
      convert_wet_dry_idx_d_ = convert_wet_dry_idx_d;
      qv_dry_ = qv_dry;
      q_h2so4_ = q_h2so4;
      q_aitken_so4_ = q_aitken_so4;
      n_aitken_ = n_aitken;
    } // set_variables
  }; // MAMMicrophysics::Postprocess

  // storage for local variables, initialized with ATMBufferManager
  struct Buffer {
    // number of local fields stored at column midpoints
    static constexpr int num_2d_mid = 11;

    // local column midpoint fields
    uview_2d z_mid;             // height at midpoints
    uview_2d dz;                // layer thickness
    uview_2d qv_dry;            // water vapor mixing ratio (dry air)
    uview_2d qc_dry;            // cloud water mass mixing ratio
    uview_2d n_qc_dry;          // cloud water number mixing ratio
    uview_2d qi_dry;            // cloud ice mass mixing ratio
    uview_2d n_qi_dry;          // cloud ice number mixing ratio
    uview_2d w_updraft;         // vertical wind velocity
    uview_2d q_h2so4_tend;      // tendency for H2SO4 gas
    uview_2d n_aitken_tend;     // tendency for aitken aerosol mode
    uview_2d q_aitken_so4_tend; // tendency for aitken mode sulfate aerosol

    // number of local fields stored at column interfaces
    static constexpr int num_2d_iface = 1;

    // local column interface fields
    uview_2d z_iface; // height at interfaces

    // storage
    Real* wsm_data;
  };

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // aerosol processes
  std::unique_ptr<mam4::NucleationProcess> nucleation_;

  // pre- and postprocessing scratch pads
  Preprocess preprocess_;
  Postprocess postprocess_;

  // local atmospheric state column variables
  const_view_2d T_mid_;   // temperature at grid midpoints [K]
  const_view_2d p_mid_;   // total pressure at grid midpoints [Pa]
  const_view_2d qv_;      // water vapor specific humidity [kg h2o vapor / kg moist air]
                          // we keep the specific humidity to use in the PostProcess step.
  view_2d qc_;      // cloud liquid mass mixing ratio
                    // must be converted from wet to dry [kg cloud water /kg dry air]
  view_2d n_qc_;      // cloud liquid number mixing ratio
                    // must be converted from wet to dry [1 /kg dry air]
  const_view_2d qi_;      // cloud ice mass mixing ratio
                    // must be converted from wet to dry [kg cloud water /kg dry air]
  const_view_2d n_qi_;      // cloud ice number mixing ratio
                    // must be converted from wet to dry [1 /kg dry air]
  const_view_2d pdel_;    // hydrostatic "pressure thickness" at grid
                          // interfaces [Pa]
  const_view_2d cloud_f_; // cloud fraction [-]
  const_view_1d pblh_;    // planetary boundary layer height [m]

  // local aerosol-related gases
  view_2d q_h2so4_; // H2SO3 gas [kg/kg dry air]

  // local aerosols (more to appear as we improve this atm process)
  view_2d n_aitken_; // aitken mode number mixing ratio [1/kg dry air]
  view_2d q_aitken_so4_; // SO4 mass mixing ratio in aitken mode [kg/kg dry air]

  // workspace manager for internal local variables
  //ekat::WorkspaceManager<Real, KT::Device> workspace_mgr_;
  Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMMicrophysics

} // namespace scream

#endif // EAMXX_MAM_MICROPHYSICS_HPP
