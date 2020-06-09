#ifndef P3_FUNCTIONS_CONSERVATION_IMPL_HPP
#define P3_FUNCTIONS_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp"
#include "physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_conservation(const Spack& qc, const Scalar dt,
  Spack& qcaut, Spack& qcacc, Spack &qccol, Spack& qcheti, Spack& qcshd, Spack& qiberg, Spack& qisub, Spack& qidep,
  const Smask& context)
{
  const auto sinks = (qcaut+qcacc+qccol+qcheti+qcshd+qiberg)*dt; // Sinks of cloud water
  const auto sources = qc; // Source of cloud water
  Spack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  Smask nothing_todo = !enforce_conservation && context;

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qcaut.set(enforce_conservation, qcaut*ratio);
    qcacc.set(enforce_conservation, qcacc*ratio);
    qccol.set(enforce_conservation, qccol*ratio);
    qcheti.set(enforce_conservation, qcheti*ratio);
    qcshd.set(enforce_conservation, qcshd*ratio);
    qiberg.set(enforce_conservation, qiberg*ratio);
  }

  if(nothing_todo.any()){
    ratio.set(nothing_todo, 1); // If not limiting sinks on qc then most likely did not run out of qc
  }

  //PMC: ratio is also frac of step w/ liq. thus we apply qiberg for
  //"ratio" of timestep and vapor deposition and sublimation  for the
  //remaining frac of the timestep.  Only limit if there will be cloud
  // water to begin with.
  enforce_conservation = sources > qtendsmall && context;
  if (enforce_conservation.any()){
    qidep.set(enforce_conservation, qidep*(1-ratio));
    qisub.set(enforce_conservation, qisub*(1-ratio));
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_water_conservation(
  const Spack& qr, const Spack& qcaut, const Spack& qcacc, const Spack& qimlt, const Spack& qcshd, const Scalar dt,
  Spack& qrevp, Spack& qrcol, Spack& qrheti,
  const Smask& context)
{
  const auto sinks   = (qrevp+qrcol+qrheti)*dt; // Sinks of rain water
  const auto sources = qr + (qcaut+qcacc+qimlt+qcshd)*dt; // Sources of rain water
  Spack ratio;

  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qrevp.set(enforce_conservation, qrevp*ratio);
    qrcol.set(enforce_conservation, qrcol*ratio);
    qrheti.set(enforce_conservation, qrheti*ratio);
  }
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_water_conservation(
  const Spack& qitot,const Spack& qidep,const Spack& qinuc,const Spack& qiberg, const Spack &qrcol,const Spack &qccol,const Spack& qrheti,const Spack& qcheti,const Scalar dt,
  Spack& qisub, Spack& qimlt,
  const Smask& context)
{
  const auto sinks = (qisub+qimlt)*dt; // Sinks of ice water
  const auto sources = qitot + (qidep+qinuc+qrcol+qccol+qrheti+qcheti+qiberg)*dt; // Sources of ice water
  Spack ratio;
  constexpr Scalar qtendsmall = C::QTENDSMALL;
  Smask enforce_conservation  = sinks > sources && sinks >= qtendsmall && context;  // determine if  conservation corrction is necessary
  if(enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qisub.set(enforce_conservation, qisub*ratio);
    qimlt.set(enforce_conservation, qimlt*ratio);
  }
}

} // namespace p3
} // namespace scream

#endif
