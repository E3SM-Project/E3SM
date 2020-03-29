#ifndef P3_FUNCTIONS_INCLOUD_MIXINGRATIOS_IMPL_HPP
#define P3_FUNCTIONS_INCLOUD_MIXINGRATIOS_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_functions_math_impl.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calculate_incloud_mixingratios(const Spack& qc, const Spack& qr, const Spack& qitot, const Spack& qirim, const Spack& nc,
                                 const Spack& nr, const Spack& nitot, const Spack& birim, const Spack& inv_lcldm, 
                                 const Spack& inv_icldm, const Spack& inv_rcldm, 
                                 Spack& qc_incld, Spack& qr_incld, Spack& qitot_incld, Spack& qirim_incld,
                                 Spack& nc_incld, Spack& nr_incld, Spack& nitot_incld, Spack& birim_incld)
{
   constexpr Scalar qsmall = C::QSMALL;
   constexpr Scalar incloud_limit = C::incloud_limit;
   constexpr Scalar precip_limit  = C::precip_limit;

   const auto qc_ge_qsmall = qc >= qsmall;
   const auto not_qc_ge_qsmall = ! qc_ge_qsmall;

   const auto qitot_ge_qsmall = qitot >= qsmall;
   const auto not_qitot_ge_qsmall = ! qitot_ge_qsmall;

   const auto qitot_and_qirim_ge_qsmall = (qitot >= qsmall) && (qirim >= qsmall);
   const auto not_qitot_and_qirim_ge_qsmall = ! qitot_and_qirim_ge_qsmall;

   const auto qr_ge_qsmall = qr >= qsmall;
   const auto not_qr_ge_qsmall = ! qr_ge_qsmall;

   qc_incld.set(qc_ge_qsmall,
                qc*inv_lcldm);
   nc_incld.set(qc_ge_qsmall,
                pack::max(nc*inv_lcldm,sp(0)));

   qc_incld.set(not_qc_ge_qsmall,
                sp(0));
   nc_incld.set(not_qc_ge_qsmall,
                sp(0));

//       if (qc.ge.qsmall) then
//          qc_incld = qc*inv_lcldm
//          nc_incld = max(nc*inv_lcldm,0._rtype)
//          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
//       else
//          qc_incld = 0._rtype
//          nc_incld = 0._rtype
//       end if 

   qitot_incld.set(qitot_ge_qsmall,
                   qitot*inv_icldm);
   nitot_incld.set(qitot_ge_qsmall,
                   pack::max(nitot*inv_icldm, sp(0)));

   qitot_incld.set(not_qitot_ge_qsmall,
                  sp(0));
   nitot_incld.set(not_qitot_ge_qsmall,
                  sp(0));

//       if (qitot.ge.qsmall) then
//          qitot_incld = qitot*inv_icldm
//          nitot_incld = max(nitot*inv_icldm,0._rtype)
//          !AaronDonahue, kai has something about if nicons then ni=ninst/rho
//       else
//          qitot_incld = 0._rtype
//          nitot_incld = 0._rtype
//       end if 

   qirim_incld.set(qitot_and_qirim_ge_qsmall,
                   qirim*inv_icldm);
   birim_incld.set(qitot_and_qirim_ge_qsmall,
                   pack::max(birim*inv_lcldm, sp(0)));

   qirim_incld.set(not_qitot_and_qirim_ge_qsmall,
                   sp(0));
   birim_incld.set(not_qitot_and_qirim_ge_qsmall,
                   sp(0));

//       if (qirim.ge.qsmall.and.qitot.ge.qsmall) then
//          qirim_incld = qirim*inv_icldm
//          birim_incld = max(birim*inv_lcldm,0._rtype)
//       else
//          qirim_incld = 0._rtype
//          birim_incld = 0._rtype
//       end if 

   qr_incld.set(qr_ge_qsmall,
                qr*inv_rcldm);
   nr_incld.set(qr_ge_qsmall,
                pack::max(nr*inv_rcldm, sp(0)));

   qr_incld.set(not_qr_ge_qsmall,
                sp(0));
   nr_incld.set(not_qr_ge_qsmall,
                sp(0));

//       if (qr.ge.qsmall) then
//          qr_incld = qr*inv_rcldm
//          nr_incld = max(nr*inv_rcldm,0._rtype)
//          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
//       else
//          qr_incld = 0._rtype
//          nr_incld = 0._rtype
//       end if

   const auto any_gt_limit = (qc_incld > incloud_limit) || (qitot_incld > incloud_limit) || (qr_incld > precip_limit)
                          || (birim_incld > incloud_limit);

   qc_incld.set(any_gt_limit,
                pack::max(qc_incld, incloud_limit));
   qitot_incld.set(any_gt_limit,
                   pack::max(qitot_incld, incloud_limit));
   birim_incld.set(any_gt_limit,
                   pack::max(birim_incld, incloud_limit));
   qr_incld.set(any_gt_limit,
                pack::max(qr_incld, precip_limit));


//       if (qc_incld.gt.incloud_limit .or.qitot_incld.gt.incloud_limit .or. qr_incld.gt.precip_limit .or.birim_incld.gt.incloud_limit) then
//!          write(errmsg,'(a3,i4,3(a5,1x,e16.8,1x))') 'k: ', k, ', qc:',qc_incld, ', qi:',qitot_incld,', qr:',qr_incld
//          qc_incld    = max(qc_incld,incloud_limit)
//          qitot_incld = max(qitot_incld,incloud_limit)
//          birim_incld = max(birim_incld,incloud_limit)
//          qr_incld    = max(qr_incld,precip_limit)
//!          if (masterproc) write(iulog,*)  errmsg
//!          call handle_errmsg('Micro-P3 (Init)',subname='In-cloud mixing
//!          ratio too large',extra_msg=errmsg)
//       end if

}

} // namespace p3
} // namespace scream

#endif
