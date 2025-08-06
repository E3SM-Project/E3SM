#ifndef SCREAM_RRTMGP_RADIATION_HPP
#define SCREAM_RRTMGP_RADIATION_HPP

#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_parameter_list.hpp>
#include <ekat_string_utils.hpp>

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

  using interface_t = rrtmgp::rrtmgp_interface<Real, layout_t, DefaultDevice>;

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
  GasConcsK<Real, layout_t, DefaultDevice> m_gas_concs_k;

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
    static constexpr int num_1d_ncol        = 8;
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
    ureal1dk sfc_alb_dir_vis_k;
    ureal1dk sfc_alb_dir_nir_k;
    ureal1dk sfc_alb_dif_vis_k;
    ureal1dk sfc_alb_dif_nir_k;
    ureal1dk sfc_flux_dir_vis_k;
    ureal1dk sfc_flux_dir_nir_k;
    ureal1dk sfc_flux_dif_vis_k;
    ureal1dk sfc_flux_dif_nir_k;

    // 2d size (ncol, nlay)
    ureal2dk d_dz;
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

    // 2d size (ncol, nlay+1)
    ureal2dk d_tint;
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

    // 3d size (ncol, nlay+1, nswbands)
    ureal3dk sw_bnd_flux_up_k;
    ureal3dk sw_bnd_flux_dn_k;
    ureal3dk sw_bnd_flux_dir_k;
    ureal3dk sw_bnd_flux_dif_k;

    // 3d size (ncol, nlay+1, nlwbands)
    ureal3dk lw_bnd_flux_up_k;
    ureal3dk lw_bnd_flux_dn_k;

    // 2d size (ncol, nswbands)
    ureal2dk sfc_alb_dir_k;
    ureal2dk sfc_alb_dif_k;

    // 3d size (ncol, nlay, n[sw,lw]bands)
    ureal3dk aero_tau_sw_k;
    ureal3dk aero_ssa_sw_k;
    ureal3dk aero_g_sw_k;
    ureal3dk aero_tau_lw_k;

    // 3d size (ncol, nlay, n[sw,lw]bnds)
    ureal3dk cld_tau_sw_bnd_k;
    ureal3dk cld_tau_lw_bnd_k;

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    ureal3dk cld_tau_sw_gpt_k;
    ureal3dk cld_tau_lw_gpt_k;
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
