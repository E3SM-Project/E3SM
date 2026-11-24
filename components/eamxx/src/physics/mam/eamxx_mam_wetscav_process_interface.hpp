#ifndef EAMXX_MAM_WETSCAV_HPP
#define EAMXX_MAM_WETSCAV_HPP

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For declaring wetscav class derived from atm process class
#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

// For ZM convective tracer transport (only available when ZM is built)
#ifdef EAMXX_HAS_ZM
#include <physics/zm/zm_functions.hpp>
#endif

// For component name
#include <memory>
#include <string>

namespace scream
{

/*
 * The class responsible to handle the aerosol wetscavenging
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 */

class MAMWetscav : public MAMGenericInterface
{
  using KT           = ekat::KokkosTypes<DefaultDevice>;
  using view_2d      = typename KT::template view_2d<Real>;
  using view_2d_host = typename KT::template view_2d<Real>::host_mirror_type;
  using int_view_1d  = typename KT::template view_1d<int>;
  using int_view_2d  = typename KT::template view_2d<int>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_1d_int = typename KT::template view_1d<const int>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

public:
  // Constructor
  MAMWetscav(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the subcomponent
  std::string
  name() const override
  {
    return "mam4_wetscav";
  }

  // Set the grid and input output variables
  void create_requests() override;

  // management of common atm process memory
  // ON HOST, returns the number of bytes of device memory needed by the above
  // Buffer type given the number of columns and vertical levels
  size_t
  requested_buffer_size_in_bytes() const override
  {
    return mam_coupling::buffer_size(ncol_, nlev_, num_2d_scratch_, len_temporary_views_);
  }
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl() override { /*Do nothing*/ };

private:
  // -----------------------------------------------
  // Local variables
  // ------------------------------------------------

  // aerosol metadata
  mam4::AeroConfig aero_config_;

  // Work arrays
  view_2d work_;
  view_2d work_convproc_;  // Separate work array for convection processing scratch1Dviews
  int_view_2d isprx_;

  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now.
  // NOTE: rprddp, icwmrdp, and evapcdp are no longer member variables;
  // when do_convproc_ is true they alias FM views directly to avoid
  // redundant copies.  When false, the zero-initialized member views below
  // serve as placeholders.

  view_2d sh_frac_;

  // Deep convective cloud fraction [fraction]
  view_2d dp_frac_;

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  view_2d evapcsh_;

  // Rain production, shallow convection [kg/kg/s]
  view_2d rprdsh_;

  // In cloud water mixing ratio, shallow convection
  view_2d icwmrsh_;

  // Deep convection zero-placeholder views (only allocated when !do_convproc_).
  // When do_convproc_ is true, run_impl aliases FM views directly instead.
  view_2d rprddp_zero_;
  view_2d icwmrdp_zero_;
  view_2d evapcdp_zero_;

  // Detraining cld H20 from deep convection [kg/kg/s]
  view_2d dlf_;

  // Temporary view for pressure thickness in mb (converted from Pa)
  view_2d dp_tmp_;
  view_2d dpdry_tmp_;

  int num_2d_scratch_ = 39;
  //
  view_2d scavimptblnum_;

  view_2d scavimptblvol_;

  // Aerosol states
  mam_coupling::AerosolState dry_aero_tends_;

  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;
  // parameters for calcsize
  mam4::modal_aero_opt::CalcsizeData calcsize_data_;
  int get_len_temporary_views();
  void init_temporary_views();
  int len_temporary_views_{0};

  // scav_fraction_in_cloud_strat_ is in-cloud scavenging fraction
  Real scav_fraction_in_cloud_strat_ = 1.00;

  // scav_fraction_in_cloud_conv_ is in-cloud convective scavenging fraction
  Real scav_fraction_in_cloud_conv_ = 0.00;

  // scav_fraction_below_cloud_strat_ is below-cloud scavenging fraction
  Real scav_fraction_below_cloud_strat_ = 0.03;

  // activation_fraction_in_cloud_conv_ is convection activation fraction
  Real activation_fraction_in_cloud_conv_ = 0.40;

  // convproc_do_aer: Enable aerosol processing in convection
  bool convproc_do_aer_ = false;

  // convproc_do_gas: Enable gas processing in convection
  bool convproc_do_gas_ = false;

  // do_convproc: Combined flag for any convection processing
  bool do_convproc_ = false;

  // ConvProc configuration for species classification and resuspension mapping
  mam4::ConvProc::Config convproc_config_;

  // Device views for species classification and resuspension mapping (for GPU kernels)
  int_view_1d species_class_device_;
  int_view_1d mmtoo_prevap_resusp_device_;

  // ---------------------------------------------------------------------------
  // ZM convective tracer transport (EAM: zm_conv_tend_2 -> zm_transport_tracer)
  // ---------------------------------------------------------------------------
  // zm_tracer_transport: transport gases (and aerosols when convproc_do_aer is
  // false) through the ZM plume after aero_model_wetdep.
  bool zm_tracer_transport_ = false;

  // Per-column scratch for the transport call, carved from the MAM buffer:
  //   [ q (nlev*pcnst) | dqdt (nlev*pcnst) | mu (nlev) | md (nlev) | dp (nlev) ]
  view_2d zm_transport_work_;
  // Length (in Reals) of the per-column scratch above; 0 without ZM
  int get_zm_transport_work_len() const;

#ifdef EAMXX_HAS_ZM
public:
  // NOTE: the ZM helpers below are public because CUDA does not allow
  // extended (KOKKOS_LAMBDA) lambdas inside private member functions.
  using ZMF = zm::Functions<Real, DefaultDevice>;

  // Everything the column kernel needs for ZM tracer transport, gathered on
  // the host by get_zm_transport_inputs() and captured by value in the kernel.
  struct ZmTransportInputs {
    bool enabled = false;               // zm_tracer_transport_
    const_view_1d cldtop_klev;          // ZM cloud top level index
    const_view_1d msemax_klev;          // ZM cloud base level index
    const_view_1d activity;             // ZM activity flag (1/0)
    const_view_2d mflx_up, mflx_dn;     // up/downdraft mass flux [kg/m2/s]
    const_view_2d entr_up, detr_up;     // updraft entrainment/detrainment [1/s]
    const_view_2d entr_dn;              // downdraft entrainment [1/s]
    Kokkos::View<const bool *> doconvtran;  // transport mask (pcnst)
    const_view_2d fracis;               // insoluble fraction (nlev, pcnst)
    view_2d work;                       // per-column scratch
    ZMF::WorkspaceManager wsm;          // workspace for zm_transport_tracer
    ZMF::ZmRuntimeOpt runtime_opt;      // ZM options (zm_microp)
    int ktm = 0;                        // min cloud top over active columns
    int kbm = 0;                        // min cloud base over active columns
  };

  // Host: retrieve ZM fields and compute ktm/kbm (enabled=false if off)
  ZmTransportInputs get_zm_transport_inputs() const;

  // Host: domain-wide cloud top/base bounds over active columns
  static void zm_cloud_level_bounds(const int ncol, const int nlev,
                                    const const_view_1d &cldtop_klev,
                                    const const_view_1d &msemax_klev,
                                    const const_view_1d &activity,
                                    int &ktm, int &kbm);

  // Per-column scratch length needed by zm_tracer_transport_column
  static int zm_tracer_transport_work_len(const int nlev);

  // Device (whole team): ZM tracer transport for column icol, applied in
  // place to the MAM dry state in progs
  KOKKOS_INLINE_FUNCTION
  static void zm_tracer_transport_column(const ThreadTeam &team,
                                         const ZmTransportInputs &zm,
                                         const int icol, const int nlev,
                                         const Real dt,
                                         const mam4::Atmosphere &atm,
                                         const mam4::Prognostics &progs,
                                         const const_view_2d &p_del);

private:
  // doconvtran(m): true if constituent m is transported by ZM
  Kokkos::View<bool *> zm_doconvtran_;
  // fraction of each tracer that is insoluble (nlev, pcnst); set to 1
  view_2d zm_fracis_;
  // workspace for the 10 internal (pcnst x nlev) arrays of zm_transport_tracer
  std::shared_ptr<ZMF::WorkspaceManager> zm_wsm_;
#endif

}; // class MAMWetscav

} // namespace scream

#endif // EAMXX_MAM_WETSCAV_HPP
