#ifndef P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP
#define P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP

#include "p3_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_deposition_sublimation(
  const Spack& qi_incld, const Spack& ni_incld, const Spack& T_atm,   const Spack& qv_sat_l,
  const Spack& qv_sat_i,         const Spack& epsi,        const Spack& abi, const Spack& qv, const Scalar& inv_dt,
  Spack& qv2qi_vapdep_tend, Spack& qi2qv_sublim_tend, Spack& ni_sublim_tend, Spack& qc2qi_berg_tend,
  const Smask& context)
{
  constexpr Scalar QSMALL   = C::QSMALL;
  constexpr Scalar T_zerodegc = C::T_zerodegc;

  Spack qi_tend;   //temporary var for mass tend before splitting into sublim or depos

  //INITIALIZE EVERYTHING TO 0:
  qc2qi_berg_tend=0;
  qv2qi_vapdep_tend=0;
  qi2qv_sublim_tend=0;
  ni_sublim_tend=0;

  //CAN'T HAVE DEPOSITION/SUBLIMATION IF NO ICE MASS
  const auto qi_incld_not_small = qi_incld > QSMALL && context;

  if (qi_incld_not_small.any()) {

    //COMPUTE NET MASS TENDENCY: USING MIN IN THE LINE BELOW TO PREVENT SUBLIM OR DEPOS
    //FROM PUSHING QV BEYOND ICE SAT WITHIN THE GIVEN TIMESTEP. APPLYING MIN HERE IS
    //EQUIVALENT TO LIMITING THE TENDENCY LATER TO ENSURE END-OF-STEP QV ISN'T
    //INAPPROPRIATELY SUPER OR SUBSATURATED.
    qi_tend.set(qi_incld_not_small,min(epsi/abi,inv_dt) * (qv - qv_sat_i));

    //SUBLIMATE WHERE qi_tend<0. MAKE POSITIVE TO MATCH CONVENTION
    const auto neg_qi_tend = (qi_tend < 0 );
    qi2qv_sublim_tend.set(qi_incld_not_small && neg_qi_tend, -qi_tend);
    ni_sublim_tend.set(qi_incld_not_small && neg_qi_tend, qi2qv_sublim_tend*(ni_incld/qi_incld));

    //DEPOSITION (FROM VAPOR OR LIQ) ONLY OCCURS BELOW FREEZING:
    const auto T_lt_frz = (T_atm < T_zerodegc);

    //BERGERON OCCURS WHERE LIQUID IS PRESENT AND DEPOSITION FROM VAPOR OCCURS WHERE IT ISN'T.
    //IF ALL LIQUID IS CONSUMED PARTWAY THROUGH A STEP, BERGERON SHOULD BE ACTIVE FOR THE
    //FRACTION OF THE STEP WHEN LIQUID IS PRESENT AND DEPOSITION FROM VAPOR SHOULD BE ACTIVE FOR
    //THE REST OF THE STEP. THE FRACTION OF THE STEP WITH LIQUID ISN'T KNOWN UNTIL THE 'CONSERVATION
    //CHECKS' AT THE END OF THE STEP, SO WE COMPUTE BERGERON AND VAPOR DEPOSITION HERE ASSUMING
    //LIQUID IS OR ISN'T PRESENT FOR THE WHOLE STEP (RESPECTIVELY).

    //VAPOR DEPOSITION
    qv2qi_vapdep_tend.set(qi_incld_not_small && T_lt_frz && !neg_qi_tend, qi_tend);

    //BERGERON: NOTE THAT AS FORMULATED, BERG DOESN'T HAVE ANYTHING TO DO WITH QV, SO CAN'T
    //PUSH IT BEYOND SATURATION. THUS, NOT LIMITING WITH INV_DT HERE.
    qc2qi_berg_tend.set(qi_incld_not_small && T_lt_frz, max(epsi/abi*(qv_sat_l - qv_sat_i), 0));

  } //end if at least 1 qi is greater than qmall

}
} // namespace p3
} // namespace scream

#endif // P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP
