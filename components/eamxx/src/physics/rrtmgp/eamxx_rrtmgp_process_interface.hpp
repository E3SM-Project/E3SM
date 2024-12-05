#ifndef SCREAM_RRTMGP_RADIATION_HPP
#define SCREAM_RRTMGP_RADIATION_HPP

#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include <string>

namespace scream {
/*
 * Class responsible for atmosphere radiative transfer. The AD should store
 * exactly ONE instance of this class in its list of subcomponents.
 */

// rrtmgp is performance tuned for layout left views but will accept any
// view. We probably want to stick with layout left views for performance
// reasons even though this requires us to make copies of our fields (they
// are layout right).
#define RRTMGP_LAYOUT_LEFT

class RRTMGPRadiation : public AtmosphereProcess {
public:
  using KT         = ekat::KokkosTypes<DefaultDevice>;
#ifdef RRTMGP_LAYOUT_LEFT
  using layout_t   = Kokkos::LayoutLeft;
#else
  using layout_t   = typename ekat::KokkosTypes<DefaultDevice>::Layout;
#endif
  using real1dk    = Kokkos::View<Real*, DefaultDevice>;
  using real2dk    = Kokkos::View<Real**, layout_t, DefaultDevice>;
  using real3dk    = Kokkos::View<Real***, layout_t, DefaultDevice>;
  using creal1dk   = Kokkos::View<const Real*, DefaultDevice>;
  using creal2dk   = Kokkos::View<const Real**, layout_t, DefaultDevice>;
  using creal3dk   = Kokkos::View<const Real***, layout_t, DefaultDevice>;
  using ureal1dk  = Unmanaged<real1dk>;
  using ureal2dk  = Unmanaged<real2dk>;
  using ureal3dk  = Unmanaged<real3dk>;
  using cureal1dk  = Unmanaged<creal1dk>;
  using cureal2dk  = Unmanaged<creal2dk>;
  using cureal3dk  = Unmanaged<creal3dk>;

  using ci_string = ekat::CaseInsensitiveString;

  using lrreal2dk   = typename KT::template view_2d<Real>;
  using ulrreal2dk  = Unmanaged<lrreal2dk>;

#ifdef RRTMGP_ENABLE_KOKKOS
  using interface_t = rrtmgp::rrtmgp_interface<Real, layout_t, DefaultDevice>;
#endif

  // Constructors
  RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of the subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "rrtmgp"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grid_manager);

// NOTE: cannot use lambda functions for CUDA devices if these are protected!
public:
  // The three main interfaces for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Keep track of number of columns and levels
  int m_ncol;
  int m_num_col_chunks;
  int m_col_chunk_size;
  std::vector<int> m_col_chunk_beg;
  int m_nlay;
  Field m_lat;
  Field m_lon;

  // Whether we use aerosol forcing in radiation
  bool m_do_aerosol_rad;
  // Whether we do extra aerosol forcing calls
  bool m_extra_clnsky_diag;
  bool m_extra_clnclrsky_diag;

  // The orbital year, used for zenith angle calculations:
  // If > 0, use constant orbital year for duration of simulation
  // If < 0, use year from timestamp for orbital parameters
  Int m_orbital_year;
  // Orbital parameters, used for zenith angle calculations.
  // If >= 0, bypass computation based on orbital year and use fixed parameters
  // If <  0, compute based on orbital year, specified above
  Real m_orbital_eccen;  // Eccentricity
  Real m_orbital_obliq;  // Obliquity
  Real m_orbital_mvelp;  // Vernal Equinox Mean Longitude of Perihelion

  // Value for prescribing an invariant solar constant (i.e. total solar irradiance
  // at TOA).  Used for idealized experiments such as RCE. This is only used when a
  // positive value is supplied.
  Real m_fixed_total_solar_irradiance;

  // Fixed solar zenith angle to use for shortwave calculations
  // This is only used if a positive value is supplied
  Real m_fixed_solar_zenith_angle;

  // Need to hard-code some dimension sizes for now.
  // TODO: find a better way of configuring this
  const int m_nswbands = 14;
  const int m_nlwbands = 16;
  int m_nswgpts;
  int m_nlwgpts;

  // These are the gases that we keep track of
  int m_ngas;
  std::vector<ci_string>   m_gas_names;
  real1dk                  m_gas_mol_weights;
#ifdef RRTMGP_ENABLE_YAKL
  GasConcs                 m_gas_concs;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
  GasConcsK<Real, layout_t, DefaultDevice> m_gas_concs_k;
#endif

  // Prescribed greenhouse gas surface concentrations in moles / moles air
  Real m_co2vmr;
  Real m_n2ovmr;
  Real m_ch4vmr;
  Real m_f11vmr;
  Real m_f12vmr;
  Real m_n2vmr;
  Real m_covmr;

  // Rad frequency in number of steps
  int m_rad_freq_in_steps;

  // Whether or not to do subcolumn sampling of cloud state for MCICA
  bool m_do_subcol_sampling;

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_1d_ncol        = 10;
    static constexpr int num_2d_nlay        = 16;
    static constexpr int num_2d_nlay_p1     = 23;
    static constexpr int num_2d_nswbands    = 2;
    static constexpr int num_3d_nlev_nswbands = 4;
    static constexpr int num_3d_nlev_nlwbands = 2;
    static constexpr int num_3d_nlay_nswbands = 4;
    static constexpr int num_3d_nlay_nlwbands = 2;
    static constexpr int num_3d_nlay_nswgpts = 1;
    static constexpr int num_3d_nlay_nlwgpts = 1;

    // 1d size (ncol)
    ureal1dk cosine_zenith;
#ifdef RRTMGP_ENABLE_YAKL
    real1d mu0;
    real1d sfc_alb_dir_vis;
    real1d sfc_alb_dir_nir;
    real1d sfc_alb_dif_vis;
    real1d sfc_alb_dif_nir;
    real1d sfc_flux_dir_vis;
    real1d sfc_flux_dir_nir;
    real1d sfc_flux_dif_vis;
    real1d sfc_flux_dif_nir;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal1dk mu0_k;
    ureal1dk sfc_alb_dir_vis_k;
    ureal1dk sfc_alb_dir_nir_k;
    ureal1dk sfc_alb_dif_vis_k;
    ureal1dk sfc_alb_dif_nir_k;
    ureal1dk sfc_flux_dir_vis_k;
    ureal1dk sfc_flux_dir_nir_k;
    ureal1dk sfc_flux_dif_vis_k;
    ureal1dk sfc_flux_dif_nir_k;
#endif

    // 2d size (ncol, nlay)
    ureal2dk d_dz;
#ifdef RRTMGP_ENABLE_YAKL
    real2d p_lay;
    real2d t_lay;
    real2d z_del;
    real2d p_del;
    real2d qc;
    real2d nc;
    real2d qi;
    real2d cldfrac_tot;
    real2d eff_radius_qc;
    real2d eff_radius_qi;
    real2d tmp2d;
    real2d lwp;
    real2d iwp;
    real2d sw_heating;
    real2d lw_heating;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal2dk p_lay_k;
    ureal2dk t_lay_k;
    ureal2dk z_del_k;
    ureal2dk p_del_k;
    ureal2dk qc_k;
    ureal2dk nc_k;
    ureal2dk qi_k;
    ureal2dk cldfrac_tot_k;
    ureal2dk eff_radius_qc_k;
    ureal2dk eff_radius_qi_k;
    ureal2dk tmp2d_k;
    ureal2dk lwp_k;
    ureal2dk iwp_k;
    ureal2dk sw_heating_k;
    ureal2dk lw_heating_k;
#endif

    // 2d size (ncol, nlay+1)
    ureal2dk d_tint;
#ifdef RRTMGP_ENABLE_YAKL
    real2d p_lev;
    real2d t_lev;
    real2d sw_flux_up;
    real2d sw_flux_dn;
    real2d sw_flux_dn_dir;
    real2d lw_flux_up;
    real2d lw_flux_dn;
    real2d sw_clnclrsky_flux_up;
    real2d sw_clnclrsky_flux_dn;
    real2d sw_clnclrsky_flux_dn_dir;
    real2d sw_clrsky_flux_up;
    real2d sw_clrsky_flux_dn;
    real2d sw_clrsky_flux_dn_dir;
    real2d sw_clnsky_flux_up;
    real2d sw_clnsky_flux_dn;
    real2d sw_clnsky_flux_dn_dir;
    real2d lw_clnclrsky_flux_up;
    real2d lw_clnclrsky_flux_dn;
    real2d lw_clrsky_flux_up;
    real2d lw_clrsky_flux_dn;
    real2d lw_clnsky_flux_up;
    real2d lw_clnsky_flux_dn;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal2dk p_lev_k;
    ureal2dk t_lev_k;
    ureal2dk sw_flux_up_k;
    ureal2dk sw_flux_dn_k;
    ureal2dk sw_flux_dn_dir_k;
    ureal2dk lw_flux_up_k;
    ureal2dk lw_flux_dn_k;
    ureal2dk sw_clnclrsky_flux_up_k;
    ureal2dk sw_clnclrsky_flux_dn_k;
    ureal2dk sw_clnclrsky_flux_dn_dir_k;
    ureal2dk sw_clrsky_flux_up_k;
    ureal2dk sw_clrsky_flux_dn_k;
    ureal2dk sw_clrsky_flux_dn_dir_k;
    ureal2dk sw_clnsky_flux_up_k;
    ureal2dk sw_clnsky_flux_dn_k;
    ureal2dk sw_clnsky_flux_dn_dir_k;
    ureal2dk lw_clnclrsky_flux_up_k;
    ureal2dk lw_clnclrsky_flux_dn_k;
    ureal2dk lw_clrsky_flux_up_k;
    ureal2dk lw_clrsky_flux_dn_k;
    ureal2dk lw_clnsky_flux_up_k;
    ureal2dk lw_clnsky_flux_dn_k;
#endif

    // 3d size (ncol, nlay+1, nswbands)
#ifdef RRTMGP_ENABLE_YAKL
    real3d sw_bnd_flux_up;
    real3d sw_bnd_flux_dn;
    real3d sw_bnd_flux_dir;
    real3d sw_bnd_flux_dif;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal3dk sw_bnd_flux_up_k;
    ureal3dk sw_bnd_flux_dn_k;
    ureal3dk sw_bnd_flux_dir_k;
    ureal3dk sw_bnd_flux_dif_k;
#endif

    // 3d size (ncol, nlay+1, nlwbands)
#ifdef RRTMGP_ENABLE_YAKL
    real3d lw_bnd_flux_up;
    real3d lw_bnd_flux_dn;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal3dk lw_bnd_flux_up_k;
    ureal3dk lw_bnd_flux_dn_k;
#endif

    // 2d size (ncol, nswbands)
#ifdef RRTMGP_ENABLE_YAKL
    real2d sfc_alb_dir;
    real2d sfc_alb_dif;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal2dk sfc_alb_dir_k;
    ureal2dk sfc_alb_dif_k;
#endif

    // 3d size (ncol, nlay, n[sw,lw]bands)
#ifdef RRTMGP_ENABLE_YAKL
    real3d aero_tau_sw;
    real3d aero_ssa_sw;
    real3d aero_g_sw;
    real3d aero_tau_lw;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal3dk aero_tau_sw_k;
    ureal3dk aero_ssa_sw_k;
    ureal3dk aero_g_sw_k;
    ureal3dk aero_tau_lw_k;
#endif

    // 3d size (ncol, nlay, n[sw,lw]bnds)
#ifdef RRTMGP_ENABLE_YAKL
    real3d cld_tau_sw_bnd;
    real3d cld_tau_lw_bnd;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal3dk cld_tau_sw_bnd_k;
    ureal3dk cld_tau_lw_bnd_k;
#endif

    // 3d size (ncol, nlay, n[sw,lw]gpts)
#ifdef RRTMGP_ENABLE_YAKL
    real3d cld_tau_sw_gpt;
    real3d cld_tau_lw_gpt;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    ureal3dk cld_tau_sw_gpt_k;
    ureal3dk cld_tau_lw_gpt_k;
#endif

  };

protected:

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  std::shared_ptr<const AbstractGrid>   m_grid;

  // Struct which contains local variables
  Buffer m_buffer;
};  // class RRTMGPRadiation

}  // namespace scream

#endif  // SCREAM_RRTMGP_RADIATION_HPP
