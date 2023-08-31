#ifndef MAM_COUPLING_HPP
#define MAM_COUPLING_HPP

#include <mam4xx/mam4.hpp>

// These data structures and functions are used to move data between EAMxx
// and mam4xx.

namespace scream::mam_coupling {

KOKKOS_INLINE_FUNCTION
constexpr int num_aero_modes() {
  return mam4::AeroConfig::num_modes();
}

KOKKOS_INLINE_FUNCTION
constexpr int num_aero_species() {
  return mam4::AeroConfig::num_aerosol_ids();
}

KOKKOS_INLINE_FUNCTION
constexpr int num_aero_gases() {
  return mam4::AeroConfig::num_gas_ids();
}

// This type stores views related to the atmospheric state used by MAM. It has
// wet and dry representations of tracer variables, as well as a simple height
// coordinate.
struct AtmosphericState {
  using KT            = ekat::KokkosTypes<DefaultDevice>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  const_view_2d T_mid;     // temperature at grid midpoints [K]
  const_view_2d p_mid;     // total pressure at grid midpoints [Pa]
  const_view_2d qv_wet;    // wet water vapor specific humidity [kg vapor / kg moist air]
  view_2d       qv_dry;    // dry water vapor mixing ratio [kg vapor / kg dry air]
  const_view_2d qc_wet;    // wet cloud liquid water mass mixing ratio [kg cloud water/kg moist air]
  view_2d       qc_dry;    // dry cloud liquid water mass mixing ratio [kg cloud water/kg dry air]
  const_view_2d nc_wet;    // wet cloud liquid water number mixing ratio [# / kg moist air]
  view_2d       nc_dry;    // dry cloud liquid water number mixing ratio [# / kg dry air]
  const_view_2d qi_wet;    // wet cloud ice water mass mixing ratio [kg cloud ice water / kg moist air]
  view_2d       qi_dry;    // dry cloud ice water mass mixing ratio [kg cloud ice water / kg dry air]
  const_view_2d ni_wet;    // wet cloud ice water number mixing ratio [# / kg moist air]
  view_2d       ni_dry;    // dry cloud ice water number mixing ratio [# / kg dry air]
  view_2d       z_mid;     // height at layer midpoints [m]
  view_2d       z_iface;   // height at layer interfaces [m]
  view_2d       dz;        // layer thickness [m]
  const_view_2d pdel;      // hydrostatic "pressure thickness" at grid interfaces [Pa]
  const_view_2d cldfrac;   // cloud fraction [-]
  const_view_2d omega;     // vertical pressure velocity [Pa/s]
  view_2d       w_updraft; // updraft velocity [m/s]
  const_view_1d pblh;      // planetary boundary layer height [m]
};

// This type stores aerosol number and mass mixing ratios evolved by MAM. It
// has wet and dry representations of each quantity. These quantities are stored
// by mode (and species, for mass mixing ratio) in the same way as they are in
// mam4xx, and indexed using mam4::AeroConfig.
struct AerosolState {
  using KT      = ekat::KokkosTypes<DefaultDevice>;
  using view_2d = typename KT::template view_2d<Real>;

  view_2d wet_aero_nmr[num_aero_modes]; // wet modal aerosol number mixing ratios [# / kg moist air]
  view_2d dry_aero_nmr[num_aero_modes]; // dry modal aerosol number mixing ratios [# / kg dry air]
  view_2d wet_aero_mmr[num_aero_modes][num_aero_species]; // wet aerosol mass mixing ratios [kg aerosol / kg moist air]
  view_2d dry_aero_mmr[num_aero_modes][num_aero_species]; // dry aerosol mass mixing ratios [kg aerosol / kg dry air]
  view_2d wet_gas_mmr[num_aero_gases]; // wet gas mass mixing ratios [kg gas / kg moist air]
  view_2d dry_gas_mmr[num_aero_gases]; // dry gas mass mixing ratios [kg gas / kg dry air]
};

} // namespace scream::mam_coupling

#endif
