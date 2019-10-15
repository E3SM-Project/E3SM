#ifndef P3_FUNCTIONS_UPDATE_PROGNOSTICS_IMPL_HPP
#define P3_FUNCTIONS_UPDATE_PROGNOSTICS_IMPL_HPP

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::update_prognostic_ice(const Spack& qcheti, const Spack& qccol,
    const Spack& qcshd,  const Spack& nccol,  const Spack& ncheti, const Spack& ncshdc,
    const Spack& qrcol,  const Spack& nrcol,  const Spack& qrheti, const Spack& nrheti,
    const Spack& nrshdr, const Spack& qimlt,  const Spack& nimlt,  const Spack& qisub,
    const Spack& qidep,  const Spack& qinuc,  const Spack& ninuc,  const Spack& nislf,
    const Spack& nisub,  const Spack& qiberg, const Spack& exner,  const Spack& xxls,
    const Spack& xlf,    const Smask log_predictNc, const Smask log_wetgrowth, const Scalar dt,
    const Spack& nmltratio, const Spack& rhorime_c, Spack& th, Spack& qv, Spack& qitot,
    Spack& nitot, Spack& qirim, Spack& birim, Spack& qc,  Spack& nc, Spack& qr,
    Spack& nr)
{
  /*const auto sinks = (qcaut+qcacc+qccol+qcheti+qcshd+qiberg)*dt; // Sinks of cloud water
  const auto sources = qc + (qcnuc)*dt; // Source of cloud water
  Spack ratio;

  Smask enforce_conservation  = sinks > sources && sinks >= C::QTENDSMALL;  // determine if  conservation corrction is necessary
  Smask nothing_todo = !enforce_conservation;

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
    ratio.set(nothing_todo, 1.0); // If not limiting sinks on qc then most likely did not run out of qc
  }


  //PMC: ratio is also frac of step w/ liq. thus we apply qiberg for
  //"ratio" of timestep and vapor deposition and sublimation  for the
  //remaining frac of the timestep.  Only limit if there will be cloud
  // water to begin with.
  enforce_conservation = sources > C::QTENDSMALL;
  if (enforce_conservation.any()){
    qidep.set(enforce_conservation, qidep*(1.0-ratio));
    qisub.set(enforce_conservation, qisub*(1.0-ratio));
  }
*/
}

} // namespace p3
} // namespace scream

#endif 