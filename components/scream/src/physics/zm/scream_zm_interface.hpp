#ifndef SCREAM_ZM_INTERFACE_HPP
#define SCREAM_ZM_INTERFACE_HPP

#include "ekat/ekat_assert.hpp"
//#include "ekat/util/scream_utils.hpp"
#include "ekat/util/ekat_file_utils.hpp"

#include "share/scream_types.hpp"

// Put everything into a scream namespace


namespace scream {

extern "C"
{

// Fortran routines to be called from C
void zm_init_f90     (const Real& limcnv_in, const bool& no_deep_pbl_in);
void zm_main_f90(const Real& lchnk, const Real& ncol, Real* t, Real* qh, Real* prec,
			Real* jctop, Real* jcbot, Real* pblh, Real *zm, Real* geos, Real* zi,
			Real* qtnd, Real* heat, Real* pap, Real* paph, Real* dpp, const Real &delt,
			Real* mcon, Real* cme, Real* cape, Real* tpert, Real* dlf, Real* plfx,
			Real* zdu, Real* rprd, Real* mu, Real* md, Real* du, Real* eu, 
			Real* ed, Real* dp, Real* dsubcld, Real* jt, Real* maxg, Real* ideep,
			const Real& lengath, Real* ql, Real* rliq, Real* landfrac, Real* hu_nm1,
			Real* cnv_nm1, Real* tm1, Real* qm1, Real** t_star, Real** q_star, 
			Real* dcape, Real* q, Real** tend_s,
			Real** tend_q, Real** cld, Real* snow, Real* ntprprd, Real* ntsnprd,
			Real** flxprec, Real** flxsnow, const Real& ztodt, Real* pguall, Real* pgdall, 
			Real* icwu, const Real& ncnst, Real*** fracis 
			); 
void zm_finalize_f90 ();

} // extern "C"

} // namespace scream

#endif // SCREAM_SHOC_INTERFACE_HPP
