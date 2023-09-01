#ifndef DP_IOP_SETINITIAL_IMPL_HPP
#define DP_IOP_SETINITIAL_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_setinitial. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::iop_setinitial(
    const Int& plev,
    const Int& pcnst,
    const Int& nelemd,
    const Int& np,
    const Int& nstep,
    const bool& use_replay,
    const bool& dynproc,
    const uview_1d<element_t>& elem,
    const uview_1d<Spack>& tobs,
    const uview_1d<Spack>& qobs)
{
  // Made these up
  constexpr Int inumliq = 0;
  constexpr Int inumice = 1;
  constexpr Int icldliq = 2;
  constexpr Int icldice = 3;

  if (!use_replay && nstep == 0 && dynproc) {
    for (Int ie = 0; ie < nelemd; ++ie) {
      for (Int j = 0; j < np; ++j) {
        for (Int i = 0; i < np; ++i) {

          // Find level where tobs is no longer zero
          Int thelev=0;
          for (Int k = 0; k < plev; ++k) {
            if (tobs(k)[0] != 0) { // TODO
              thelev=k;
              break;
            }
          }

          if (nstep <= 1) {
            for (Int k = 0; k < thelev; ++k) {
              tobs(k)[0]=elem(ie).m_forcing.m_ft(ie,i,j,k)[0]; // TODO
              //qobs(k)=elem(ie).m_state.Q(i,j,k,0);
            }
          }
          else {
            for (Int k = 0; k < plev; ++k) {
              tobs(k)[0]=elem(ie).m_forcing.m_ft(ie,i,j,k)[0]; // TODO
              //qobs(k)=elem(ie).m_state.Q(i,j,k,0);
            }
          }

//           if (get_nstep() .eq. 0) then
//             do cix = 1, pcnst
//                if (scm_zero_non_iop_tracers) elem(ie)%state%Q(i,j,:,cix) = qmin(cix)
//             end do
//             do k=thelev, PLEV
//               if (have_t) elem(ie)%derived%FT(i,j,k)=tobs(k)
//               if (have_q) elem(ie)%state%Q(i,j,k,1)=qobs(k)
//             enddo

//             do k=1,PLEV
//               if (have_ps) elem(ie)%state%ps_v(i,j,1) = psobs
//               if (have_u) elem(ie)%state%v(i,j,1,k,1) = uobs(k)
//               if (have_v) elem(ie)%state%v(i,j,2,k,1) = vobs(k)
//               if (have_numliq) elem(ie)%state%Q(i,j,k,inumliq) = numliqobs(k)
//               if (have_cldliq) elem(ie)%state%Q(i,j,k,icldliq) = cldliqobs(k)
//               if (have_numice) elem(ie)%state%Q(i,j,k,inumice) = numiceobs(k)
//               if (have_cldice) elem(ie)%state%Q(i,j,k,icldice) = cldiceobs(k)
//               !  If DP-CRM mode we do NOT want to write over the dy-core vertical
//               !    velocity with the large-scale one.  wfld is used in forecast.F90
//               !    for the compuation of the large-scale subsidence.
//               if (have_omega .and. .not. dp_crm) elem(ie)%derived%omega_p(i,j,k) = wfld(k)
//               if (dp_crm) elem(ie)%derived%omega_p(i,j,k) = 0.0_real_kind
//             enddo

//           endif

//         enddo
//       enddo
//     enddo
//   endif

        }
      }
    }
  }

//   ! If DP-CRM mode then SHOC/CLUBB needs to know about grid
//   !   length size.  The calculations of this based on a sphere in the
//   !   SHOC and CLUBB interefaces are not valid for a planar grid, thus
//   !   save the grid length from the dycore. Note that planar dycore
//   !   only supports uniform grids, thus we only save one value.
//   ! Set this if it is the first time step or the first restart step
//   if ((get_nstep() .eq. 0 .or. is_first_restart_step()) .and. dp_crm .and. par%dynproc) then
//     do ie=1,nelemd
//         dyn_dx_size = elem(ie)%dx_short * 1000.0_real_kind
//     enddo
//   endif
}

} // namespace dp
} // namespace scream

#endif
