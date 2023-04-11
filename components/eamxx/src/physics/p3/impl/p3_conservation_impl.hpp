#ifndef P3_CONSERVATION_IMPL_HPP
#define P3_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_conservation(const Spack& qc, const Scalar dt,
  Spack& qc2qr_autoconv_tend, Spack& qc2qr_accret_tend, Spack &qc2qi_collect_tend, Spack& qc2qi_hetero_freeze_tend, 
  Spack& qc2qr_ice_shed_tend, Spack& qc2qi_berg_tend, Spack& qi2qv_sublim_tend, Spack& qv2qi_vapdep_tend,
  const Smask& context)
{
  const auto sinks = (qc2qr_autoconv_tend+qc2qr_accret_tend+qc2qi_collect_tend+qc2qi_hetero_freeze_tend+qc2qr_ice_shed_tend+qc2qi_berg_tend)*dt; // Sinks of cloud water
  const auto sources = qc; // Source of cloud water
  Spack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  Smask nothing_todo = !enforce_conservation && context;

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qc2qr_autoconv_tend.set(enforce_conservation, qc2qr_autoconv_tend*ratio);
    qc2qr_accret_tend.set(enforce_conservation, qc2qr_accret_tend*ratio);
    qc2qi_collect_tend.set(enforce_conservation, qc2qi_collect_tend*ratio);
    qc2qi_hetero_freeze_tend.set(enforce_conservation, qc2qi_hetero_freeze_tend*ratio);
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
  enforce_conservation = sources > qtendsmall && context;
  if (enforce_conservation.any()){
    qv2qi_vapdep_tend.set(enforce_conservation, qv2qi_vapdep_tend*(1-ratio));
    qi2qv_sublim_tend.set(enforce_conservation, qi2qv_sublim_tend*(1-ratio));
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_water_conservation(
  const Spack& qr, const Spack& qc2qr_autoconv_tend, const Spack& qc2qr_accret_tend, 
  const Spack& qi2qr_melt_tend, const Spack& qc2qr_ice_shed_tend, const Scalar dt,
  Spack& qr2qv_evap_tend, Spack& qr2qi_collect_tend, Spack& qr2qi_immers_freeze_tend,
  const Smask& context)
{
  const auto sinks   = (qr2qv_evap_tend+qr2qi_collect_tend+qr2qi_immers_freeze_tend)*dt; // Sinks of rain water
  const auto sources = qr + (qc2qr_autoconv_tend+qc2qr_accret_tend+qi2qr_melt_tend+qc2qr_ice_shed_tend)*dt; // Sources of rain water
  Spack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary

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
  const Spack& qi,const Spack& qv2qi_vapdep_tend,const Spack& qv2qi_nucleat_tend,const Spack& qc2qi_berg_tend, 
  const Spack &qr2qi_collect_tend,const Spack &qc2qi_collect_tend,const Spack& qr2qi_immers_freeze_tend,
  const Spack& qc2qi_hetero_freeze_tend,const Scalar dt,
  Spack& qi2qv_sublim_tend, Spack& qi2qr_melt_tend,
  const Smask& context)
{
  const auto sinks = (qi2qv_sublim_tend+qi2qr_melt_tend)*dt; // Sinks of ice water
  const auto sources = qi + (qv2qi_vapdep_tend+qv2qi_nucleat_tend+qr2qi_collect_tend+qc2qi_collect_tend
                          + qr2qi_immers_freeze_tend+qc2qi_hetero_freeze_tend+qc2qi_berg_tend)*dt; // Sources of ice water
  Spack ratio;
  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  if(enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qi2qv_sublim_tend.set(enforce_conservation, qi2qv_sublim_tend*ratio);
    qi2qr_melt_tend.set(enforce_conservation, qi2qr_melt_tend*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
