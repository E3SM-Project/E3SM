#ifndef P3_CONSERVATION_IMPL_HPP
#define P3_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp"
#include "share/physics/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_conservation(const Pack& qc, const Scalar dt,
  Pack& qc2qr_autoconv_tend, Pack& qc2qr_accret_tend, Pack &qc2qi_collect_tend, Pack& qc2qi_hetero_freeze_tend, 
  Pack& qc2qr_ice_shed_tend, Pack& qc2qi_berg_tend, Pack& qi2qv_sublim_tend, Pack& qv2qi_vapdep_tend,
  Pack& qcheti_cnt, Pack& qicnt, const bool& use_hetfrz_classnuc, const Mask& context,
  const Pack& cld_frac_l, const Pack& cld_frac_i, const P3Runtime& runtime_options)
{

  Pack sinks;
  if(use_hetfrz_classnuc){
    sinks = (qc2qr_autoconv_tend+qc2qr_accret_tend+qc2qi_collect_tend+qcheti_cnt+qc2qr_ice_shed_tend+qc2qi_berg_tend)*dt; // Sinks of cloud water
  }
  else{
    sinks = (qc2qr_autoconv_tend+qc2qr_accret_tend+qc2qi_collect_tend+qc2qi_hetero_freeze_tend+qc2qr_ice_shed_tend+qc2qi_berg_tend)*dt; // Sinks of cloud water
  }
  const auto sources = qc; // Source of cloud water
  // il_cldm is the intersection of ice and liquid cloud fractions
  const auto il_cldm = (runtime_options.use_separate_ice_liq_frac)
                           ? min(cld_frac_i, cld_frac_l)
                           : Pack(1);
  const auto cld_frac_glaciated = (runtime_options.use_separate_ice_liq_frac)
                           ? max(cld_frac_i-il_cldm, 0.0001)
                           : Pack(1);
  Pack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Mask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  Mask nothing_todo = !enforce_conservation && context;

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qc2qr_autoconv_tend.set(enforce_conservation, qc2qr_autoconv_tend*ratio);
    qc2qr_accret_tend.set(enforce_conservation, qc2qr_accret_tend*ratio);
    qc2qi_collect_tend.set(enforce_conservation, qc2qi_collect_tend*ratio);
    if(use_hetfrz_classnuc){
         qcheti_cnt.set(enforce_conservation, qcheti_cnt*ratio);
         qicnt.set(enforce_conservation, qicnt*ratio);
    }
    else{
      qc2qi_hetero_freeze_tend.set(enforce_conservation, qc2qi_hetero_freeze_tend*ratio);
    }
    qc2qr_ice_shed_tend.set(enforce_conservation, qc2qr_ice_shed_tend*ratio);
    qc2qi_berg_tend.set(enforce_conservation, qc2qi_berg_tend*ratio);
  }

  if(nothing_todo.any()){
    ratio.set(nothing_todo, 1); // If not limiting sinks on qc then most likely did not run out of qc
  }

  //PMC: ratio is also frac of step w/ liq. thus we apply qc2qi_berg_tend for
  //"ratio" of timestep and vapor deposition and sublimation  for the
  //remaining frac of the timestep.  Only limit if there will be cloud
  // water to begin with.
  // for the case of separate ice and liquid cloud fractions,
  // in instances where ratio < 1 and we have qc, qidep needs to take over
  // but this is an in addition to the qidep we computed outside the mixed
  // phase cloud. qidep*(1._rtype-ratio)*(il_cldm/cld_frac_i) is the additional
  // vapor depositional growth rate that takes place within the mixed phase cloud
  // after qc is depleted
  enforce_conservation = sources > qtendsmall && context;
  if (enforce_conservation.any()){
    if (runtime_options.use_separate_ice_liq_frac) {
      qv2qi_vapdep_tend.set(enforce_conservation, qv2qi_vapdep_tend + qv2qi_vapdep_tend*(1-ratio)*(il_cldm/cld_frac_glaciated));
      qi2qv_sublim_tend.set(enforce_conservation, qi2qv_sublim_tend + qi2qv_sublim_tend*(1-ratio)*(il_cldm/cld_frac_glaciated));
    } else {
      qv2qi_vapdep_tend.set(enforce_conservation, qv2qi_vapdep_tend*(1-ratio));
      qi2qv_sublim_tend.set(enforce_conservation, qi2qv_sublim_tend*(1-ratio));
    }
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_water_conservation(
  const Pack& qr, const Pack& qc2qr_autoconv_tend, const Pack& qc2qr_accret_tend, 
  const Pack& qi2qr_melt_tend, const Pack& qc2qr_ice_shed_tend, const Scalar dt,
  Pack& qr2qv_evap_tend, Pack& qr2qi_collect_tend, Pack& qr2qi_immers_freeze_tend,
  const Mask& context)
{
  const auto sinks   = (qr2qv_evap_tend+qr2qi_collect_tend+qr2qi_immers_freeze_tend)*dt; // Sinks of rain water
  const auto sources = qr + (qc2qr_autoconv_tend+qc2qr_accret_tend+qi2qr_melt_tend+qc2qr_ice_shed_tend)*dt; // Sources of rain water
  Pack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Mask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qr2qv_evap_tend.set(enforce_conservation, qr2qv_evap_tend*ratio);
    qr2qi_collect_tend.set(enforce_conservation, qr2qi_collect_tend*ratio);
    qr2qi_immers_freeze_tend.set(enforce_conservation, qr2qi_immers_freeze_tend*ratio);
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_water_conservation(
  const Pack& qi,const Pack& qv2qi_vapdep_tend,const Pack& qv2qi_nucleat_tend,const Pack& qc2qi_berg_tend, 
  const Pack &qr2qi_collect_tend,const Pack &qc2qi_collect_tend,const Pack& qr2qi_immers_freeze_tend,
  const Pack& qc2qi_hetero_freeze_tend,const Scalar dt, Pack &qinuc_cnt, Pack &qcheti_cnt, Pack &qicnt,
  Pack& qi2qv_sublim_tend, Pack& qi2qr_melt_tend, const bool& use_hetfrz_classnuc,
  const Mask& context)
{
  const auto sinks = (qi2qv_sublim_tend+qi2qr_melt_tend)*dt; // Sinks of ice water

  Pack sources; 
  if(use_hetfrz_classnuc){
    sources = qi + (qv2qi_vapdep_tend+qv2qi_nucleat_tend+qr2qi_collect_tend+qc2qi_collect_tend
                          + qr2qi_immers_freeze_tend+qc2qi_berg_tend+qinuc_cnt+qcheti_cnt+qicnt)*dt; // Sources of ice water
  }
  else{
    sources = qi + (qv2qi_vapdep_tend+qv2qi_nucleat_tend+qr2qi_collect_tend+qc2qi_collect_tend
                          + qr2qi_immers_freeze_tend+qc2qi_hetero_freeze_tend+qc2qi_berg_tend)*dt; // Sources of ice water
  }
  Pack ratio;
  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Mask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  if(enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qi2qv_sublim_tend.set(enforce_conservation, qi2qv_sublim_tend*ratio);
    qi2qr_melt_tend.set(enforce_conservation, qi2qr_melt_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
