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

class RRTMGPRadiation : public AtmosphereProcess {
public:
  using view_1d_real     = typename ekat::KokkosTypes<DefaultDevice>::template view_1d<Real>;
  using view_2d_real     = typename ekat::KokkosTypes<DefaultDevice>::template view_2d<Real>;
  using view_3d_real     = typename ekat::KokkosTypes<DefaultDevice>::template view_3d<Real>;
  using view_2d_real_const = typename ekat::KokkosTypes<DefaultDevice>::template view_2d<const Real>;
  using ci_string        = ekat::CaseInsensitiveString;

  using KT               = ekat::KokkosTypes<DefaultDevice>;
  template<typename ScalarT>
  using uview_1d         = Unmanaged<typename KT::template view_1d<ScalarT>>;

  // Constructors
  RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of the subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Radiation"; }

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
  view_1d_real             m_gas_mol_weights;
  GasConcs                 m_gas_concs;

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
    static constexpr int num_2d_nlay        = 13;
    static constexpr int num_2d_nlay_p1     = 12;
    static constexpr int num_2d_nswbands    = 2;
    static constexpr int num_3d_nlev_nswbands = 4;
    static constexpr int num_3d_nlev_nlwbands = 2;
    static constexpr int num_3d_nlay_nswbands = 3;
    static constexpr int num_3d_nlay_nlwbands = 1;
    static constexpr int num_3d_nlay_nswgpts = 1;
    static constexpr int num_3d_nlay_nlwgpts = 1;

    // 1d size (ncol)
    real1d mu0;
    real1d sfc_alb_dir_vis;
    real1d sfc_alb_dir_nir;
    real1d sfc_alb_dif_vis;
    real1d sfc_alb_dif_nir;
    uview_1d<Real> cosine_zenith;
    real1d sfc_flux_dir_vis;
    real1d sfc_flux_dir_nir;
    real1d sfc_flux_dif_vis;
    real1d sfc_flux_dif_nir;

    // 2d size (ncol, nlay)
    real2d p_lay;
    real2d t_lay;
    real2d p_del;
    real2d qc;
    real2d qi;
    real2d cldfrac_tot;
    real2d eff_radius_qc;
    real2d eff_radius_qi;
    real2d tmp2d;
    real2d lwp;
    real2d iwp;
    real2d sw_heating;
    real2d lw_heating;

    // 2d size (ncol, nlay+1)
    real2d p_lev;
    real2d t_lev;
    real2d sw_flux_up;
    real2d sw_flux_dn;
    real2d sw_flux_dn_dir;
    real2d lw_flux_up;
    real2d lw_flux_dn;
    real2d sw_clrsky_flux_up;
    real2d sw_clrsky_flux_dn;
    real2d sw_clrsky_flux_dn_dir;
    real2d lw_clrsky_flux_up;
    real2d lw_clrsky_flux_dn;

    // 3d size (ncol, nlay+1, nswbands)
    real3d sw_bnd_flux_up;
    real3d sw_bnd_flux_dn;
    real3d sw_bnd_flux_dir;
    real3d sw_bnd_flux_dif;

    // 3d size (ncol, nlay+1, nlwbands)
    real3d lw_bnd_flux_up;
    real3d lw_bnd_flux_dn;

    // 2d size (ncol, nswbands)
    real2d sfc_alb_dir;
    real2d sfc_alb_dif;

    // 3d size (ncol, nlay, n[sw,lw]bands)
    real3d aero_tau_sw;
    real3d aero_ssa_sw;
    real3d aero_g_sw;
    real3d aero_tau_lw;

    // 3d size (ncol, nlay, n[sw,lw]gpts)
    real3d cld_tau_sw_gpt;
    real3d cld_tau_lw_gpt;
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
