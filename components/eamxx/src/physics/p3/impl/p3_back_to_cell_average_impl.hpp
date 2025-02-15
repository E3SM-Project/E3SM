#ifndef P3_BACK_TO_CELL_AVERAGE_IMPL_HPP
#define P3_BACK_TO_CELL_AVERAGE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cell averaging function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::back_to_cell_average(
  const Spack& cld_frac_l, const Spack& cld_frac_r,
  const Spack& cld_frac_i, Spack& qc2qr_accret_tend, Spack& qr2qv_evap_tend,
  Spack& qc2qr_autoconv_tend, Spack& nc_accret_tend, Spack& nc_selfcollect_tend,
  Spack& nc2nr_autoconv_tend, Spack& nr_selfcollect_tend, Spack& nr_evap_tend,
  Spack& ncautr,
  Spack& qi2qv_sublim_tend, Spack& nr_ice_shed_tend, Spack& qc2qi_hetero_freeze_tend,
  Spack& qr2qi_collect_tend, Spack& qc2qr_ice_shed_tend, Spack& qi2qr_melt_tend,
  Spack& qc2qi_collect_tend, Spack& qr2qi_immers_freeze_tend, Spack& ni2nr_melt_tend,
  Spack& nc_collect_tend, Spack& ncshdc, Spack& nc2ni_immers_freeze_tend,
  Spack& nr_collect_tend, Spack& ni_selfcollect_tend, Spack& qv2qi_vapdep_tend,
  Spack& nr2ni_immers_freeze_tend, Spack& ni_sublim_tend, Spack& qv2qi_nucleat_tend,
  Spack& ni_nucleat_tend, Spack& qc2qi_berg_tend,
  Spack& ncheti_cnt, Spack& qcheti_cnt, Spack& nicnt, Spack& qicnt, Spack& ninuc_cnt,
  Spack& qinuc_cnt, const Smask& context)
{
  Spack ir_cldm, il_cldm, lr_cldm;
  ir_cldm = min(cld_frac_i,cld_frac_r); // Intersection of ICE and RAIN cloud
  il_cldm = min(cld_frac_i,cld_frac_l); // Intersection of ICE and LIQUID cloud
  lr_cldm = min(cld_frac_l,cld_frac_r); // Intersection of LIQUID and RAIN cloud

  // Some process rates take place within the intersection of liquid, rain and
  // ice cloud fractions. We calculate the intersection as the minimum between
  // combinations of cloud fractions and use these values to map back to
  // cell-average quantities where applicable.

  // map warm-phase process rates to cell-avg
  qc2qr_accret_tend.set(context, qc2qr_accret_tend * lr_cldm);  // Accretion of liquid to rain
  qr2qv_evap_tend.set(context, qr2qv_evap_tend * cld_frac_r);    // Evaporation of rain
  qc2qr_autoconv_tend.set(context, qc2qr_autoconv_tend * cld_frac_l);    // Autoconversion of liquid
  nc_accret_tend.set(context, nc_accret_tend * lr_cldm);  // Number change due to accretion
  nc_selfcollect_tend.set(context, nc_selfcollect_tend * cld_frac_l);    // Self collection occurs locally in liq. cloud
  nc2nr_autoconv_tend.set(context, nc2nr_autoconv_tend * cld_frac_l);   // Impact of autoconversion on number
  nr_selfcollect_tend.set(context, nr_selfcollect_tend * cld_frac_r);    // Self collection occurs locally in rain cloud
  nr_evap_tend.set(context, nr_evap_tend * cld_frac_r);    // Change in rain number due to evaporation
  ncautr.set(context, ncautr * lr_cldm); // Autoconversion of rain drops within rain/liq cloud

  // map ice-phase  process rates to cell-avg
  qi2qv_sublim_tend.set(context, qi2qv_sublim_tend * cld_frac_i);    // Sublimation of ice in ice cloud
  nr_ice_shed_tend.set(context, nr_ice_shed_tend * il_cldm); // Rain # increase due to shedding from rain-ice collisions, occurs when ice and liquid interact
  qc2qi_hetero_freeze_tend.set(context, qc2qi_hetero_freeze_tend * il_cldm); // Immersion freezing of cloud drops
  qr2qi_collect_tend.set(context, qr2qi_collect_tend * ir_cldm);  // Collection of rain mass by ice
  qc2qr_ice_shed_tend.set(context, qc2qr_ice_shed_tend * il_cldm);  // Rain mass growth due to shedding of fain drops after collisions with ice, occurs when ice and liquid interact
  qi2qr_melt_tend.set(context, qi2qr_melt_tend * cld_frac_i);    // Melting of ice
  qc2qi_collect_tend.set(context, qc2qi_collect_tend * il_cldm);  // Collection of water by ice
  qr2qi_immers_freeze_tend.set(context, qr2qi_immers_freeze_tend * cld_frac_r);   // Immersion freezing of rain
  ni2nr_melt_tend.set(context, ni2nr_melt_tend * cld_frac_i);    // Change in number due to melting
  nc_collect_tend.set(context, nc_collect_tend * il_cldm);  // Cloud # change due to collection of cld water by ice
  ncshdc.set(context, ncshdc * il_cldm); // Number change due to shedding, occurs when ice and liquid interact
  nc2ni_immers_freeze_tend.set(context, nc2ni_immers_freeze_tend * cld_frac_l);   // Number change associated with freexzing of cld drops
  nr_collect_tend.set(context, nr_collect_tend * ir_cldm);  // Rain number change due to collection from ice
  ni_selfcollect_tend.set(context, ni_selfcollect_tend * cld_frac_i);    // Ice self collection
  qv2qi_vapdep_tend.set(context, qv2qi_vapdep_tend * cld_frac_i);    // Vapor deposition to ice phase
  nr2ni_immers_freeze_tend.set(context, nr2ni_immers_freeze_tend * cld_frac_r);   // Change in number due to immersion freezing of rain
  ni_sublim_tend.set(context, ni_sublim_tend * cld_frac_i);    // Number change due to sublimation of ice
  qc2qi_berg_tend.set(context, qc2qi_berg_tend * il_cldm); // Bergeron process

  ncheti_cnt.set(context,ncheti_cnt*cld_frac_l);
  qcheti_cnt.set(context, qcheti_cnt*cld_frac_l);
  nicnt.set(context, nicnt*cld_frac_l);
  qicnt.set(context, qicnt*cld_frac_l);
  ninuc_cnt.set(context, ninuc_cnt*cld_frac_l);
  qinuc_cnt.set(context, qinuc_cnt*cld_frac_l);

  // AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
  //qv2qi_nucleat_tend = qv2qi_nucleat_tend;           // Deposition and condensation-freezing nucleation, already cell-averaged
  //ni_nucleat_tend = ni_nucleat_tend;           // Number change due to deposition and condensation-freezing, already cell-averaged
}

} // namespace p3
} // namespace scream

#endif
