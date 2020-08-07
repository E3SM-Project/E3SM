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
  const Spack& lcldm, const Spack& rcldm,
  const Spack& icldm, Spack& qcacc, Spack& qrevp,
  Spack& qcaut, Spack& ncacc, Spack& ncslf,
  Spack& ncautc, Spack& nrslf, Spack& nrevp,
  Spack& ncautr,
  Spack& qisub, Spack& nrshdr, Spack& qcheti,
  Spack& qrcol, Spack& qcshd, Spack& qimlt,
  Spack& qccol, Spack& qrheti, Spack& nimlt,
  Spack& nccol, Spack& ncshdc, Spack& ncheti,
  Spack& nrcol, Spack& nislf, Spack& qidep,
  Spack& nrheti, Spack& nisub, Spack& qinuc,
  Spack& ninuc, Spack& qiberg,
  const Smask& context)
{
  Spack ir_cldm, il_cldm, lr_cldm;
  ir_cldm = min(icldm,rcldm); // Intersection of ICE and RAIN cloud
  il_cldm = min(icldm,lcldm); // Intersection of ICE and LIQUID cloud
  lr_cldm = min(lcldm,rcldm); // Intersection of LIQUID and RAIN cloud

  // Some process rates take place within the intersection of liquid, rain and
  // ice cloud fractions. We calculate the intersection as the minimum between
  // combinations of cloud fractions and use these values to map back to
  // cell-average quantities where applicable.

  // map warm-phase process rates to cell-avg
  qcacc.set(context, qcacc * lr_cldm);  // Accretion of liquid to rain
  qrevp.set(context, qrevp * rcldm);    // Evaporation of rain
  qcaut.set(context, qcaut * lcldm);    // Autoconversion of liquid
  ncacc.set(context, ncacc * lr_cldm);  // Number change due to accretion
  ncslf.set(context, ncslf * lcldm);    // Self collection occurs locally in liq. cloud
  ncautc.set(context, ncautc * lcldm);   // Impact of autoconversion on number
  nrslf.set(context, nrslf * rcldm);    // Self collection occurs locally in rain cloud
  nrevp.set(context, nrevp * rcldm);    // Change in rain number due to evaporation
  ncautr.set(context, ncautr * lr_cldm); // Autoconversion of rain drops within rain/liq cloud

  // map ice-phase  process rates to cell-avg
  qisub.set(context, qisub * icldm);    // Sublimation of ice in ice cloud
  nrshdr.set(context, nrshdr * il_cldm); // Rain # increase due to shedding from rain-ice collisions, occurs when ice and liquid interact
  qcheti.set(context, qcheti * il_cldm); // Immersion freezing of cloud drops
  qrcol.set(context, qrcol * ir_cldm);  // Collection of rain mass by ice
  qcshd.set(context, qcshd * il_cldm);  // Rain mass growth due to shedding of fain drops after collisions with ice, occurs when ice and liquid interact
  qimlt.set(context, qimlt * icldm);    // Melting of ice
  qccol.set(context, qccol * il_cldm);  // Collection of water by ice
  qrheti.set(context, qrheti * rcldm);   // Immersion freezing of rain
  nimlt.set(context, nimlt * icldm);    // Change in number due to melting
  nccol.set(context, nccol * il_cldm);  // Cloud # change due to collection of cld water by ice
  ncshdc.set(context, ncshdc * il_cldm); // Number change due to shedding, occurs when ice and liquid interact
  ncheti.set(context, ncheti * lcldm);   // Number change associated with freexzing of cld drops
  nrcol.set(context, nrcol * ir_cldm);  // Rain number change due to collection from ice
  nislf.set(context, nislf * icldm);    // Ice self collection
  qidep.set(context, qidep * icldm);    // Vapor deposition to ice phase
  nrheti.set(context, nrheti * rcldm);   // Change in number due to immersion freezing of rain
  nisub.set(context, nisub * icldm);    // Number change due to sublimation of ice
  qiberg.set(context, qiberg * il_cldm); // Bergeron process

  // AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
  //qinuc = qinuc;           // Deposition and condensation-freezing nucleation, already cell-averaged
  //ninuc = ninuc;           // Number change due to deposition and condensation-freezing, already cell-averaged
}

} // namespace p3
} // namespace scream

#endif
