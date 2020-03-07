#ifndef P3_FUNCTIONS_ICE_NUCLEATION_IMPL_HPP
#define P3_FUNCTIONS_ICE_NUCLEATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_functions_math_impl.hpp"


namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_nucleation(const Spack& temp, const Spack& inv_rho, const Spack& nitot, const Spack& naai,
                 const Spack& supi, const Spack& odt, const Smask& log_predictNc, 
                 Spack& qinuc, Spack& ninuc)
{
   constexpr Scalar nsmall  = C::NSMALL;
   constexpr Scalar tmelt   = C::Tmelt;
   constexpr Scalar icenuct = C::Tmelt-15.0;
   constexpr Scalar zero    = C::ZERO;
   constexpr Scalar piov3   = C::PIOV3;
   constexpr Scalar mi0     = 4.0*piov3*900.0*1.e-18;
 
   const auto t_lt_icenuct = temp < icenuct;
   const auto supi_ge_005 = supi >= 0.05;
   
   const auto any_if_log     = t_lt_icenuct && supi_ge_005 && log_predictNc;
   const auto any_if_not_log = t_lt_icenuct && supi_ge_005 && (!log_predictNc);

   Spack dum{0.0}, N_nuc{0.0}, Q_nuc{0.0};

   dum.set(any_if_not_log,
           sp(0.005)*pack::exp(sp(0.304)*(tmelt-temp))*sp(1.0e3)*inv_rho);

   dum.set(any_if_not_log,
           pack::min(dum, sp(1.0e5)*inv_rho));

   N_nuc.set(any_if_not_log,
             pack::max(zero, (dum-nitot)*odt));

   const auto n_nuc_ge_nsmall = N_nuc >= nsmall;

   Q_nuc.set(any_if_not_log && n_nuc_ge_nsmall,
             pack::max(zero, (dum-nitot)*mi0*odt));

   qinuc.set(any_if_not_log && n_nuc_ge_nsmall,
             Q_nuc);

   ninuc.set(any_if_not_log && n_nuc_ge_nsmall,
             N_nuc);

   ninuc.set(any_if_log,
             pack::max(zero, (naai-nitot)*odt));

   qinuc.set(any_if_log,
             ninuc*mi0);

}

} // namespace p3
} // namespace scream

#endif
