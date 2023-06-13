#ifndef CAM
#include "config.h"

module newphysics

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use control_mod,          only: theta_hydrostatic_mode,&
                   case_planar_bubble, bubble_prec_type, bubble_rj_cpstar_hy, bubble_rj_cpstar_nh,bubble_rj_nosedim
use dimensions_mod,       only: np, nlev, nlevp , qsize, qsize_d, nelemd
use element_mod,          only: element_t
use element_state,        only: nt=>timelevels
use hybrid_mod,           only: hybrid_t
use kinds,                only: rl=>real_kind, iulog
use parallel_mod,         only: abortmp,iam
use reduction_mod,        only: parallelmax, parallelmin

use physical_constants,   only: bubble_const1, bubble_const2, bubble_const3, bubble_const4, &
                                bubble_t0_const, bubble_epsilo, bubble_e0, &
                                latvap, latice, &
                                gravit=>g, p0, &
                                rdry=>rgas, cpdry=>cp, cpv=>cpwater_vapor, cl, cvdry, cvv, &
                                rvapor=>rwater_vapor, rhow

implicit none

real(rl) :: Tforcl  = 300.0_rl

contains

!pi is the result, it is \pi at midlevels
subroutine construct_hydro_pressure(dp,ptop,pi)

  real(rl), dimension(nlev), intent(in)  :: dp
  real(rl),                  intent(in)  :: ptop
  real(rl), dimension(nlev), intent(out) :: pi
  integer :: k
  real(rl) :: pi_upper

  pi_upper = ptop
  do k = 1, nlev
    pi_upper = pi_upper + dp(k)
    pi(k) = pi_upper - dp(k)/2
  enddo

end subroutine construct_hydro_pressure

subroutine convert_to_dry(q1,q2,q3,dp,dpdry,q1dry,q2dry,q3dry)

  real(rl), dimension(nlev), intent(in)  :: q1,q2,q3,dp,dpdry
  real(rl), dimension(nlev), intent(out) :: q1dry,q2dry,q3dry
  q1dry = q1 * dp /dpdry
  q2dry = q2 * dp /dpdry
  q3dry = q3 * dp /dpdry

end subroutine convert_to_dry

subroutine convert_to_wet(q1dry,q2dry,q3dry,dp,dpdry,q1,q2,q3)

  real(rl), dimension(nlev), intent(in)  :: q1dry,q2dry,q3dry,dp,dpdry
  real(rl), dimension(nlev), intent(out) :: q1,q2,q3
  q1 = q1dry * dpdry / dp
  q2 = q2dry * dpdry / dp
  q3 = q3dry * dpdry / dp

end subroutine convert_to_wet

subroutine compute_mass(qvdry,qcdry,qrdry,dpdry,mass)

  real(rl), dimension(nlev), intent(in)    :: qvdry,qcdry,qrdry
  real(rl), dimension(nlev), intent(in)    :: dpdry
  real(rl),                  intent(inout) :: mass

  integer  :: k

  mass = 0.0
  do k=1,nlev
    mass = mass + dpdry(k)*( 1.0 + qvdry(k) + qcdry(k) + qrdry(k) )
  enddo

end subroutine compute_mass




subroutine rj_new(qv_c,ql_c,T_c,dp_c,p_c,zi_c,ptop,massout,energyout,&
                                        en1, en2cp, en2cV, en3, encl, wasiactive)

  real(rl), dimension(nlev), intent(inout) :: p_c
  real(rl), dimension(nlev), intent(inout) :: dp_c
  real(rl), dimension(nlev), intent(inout) :: qv_c,ql_c,T_c
  real(rl), dimension(nlevp),intent(inout) :: zi_c
  real(rl),                  intent(inout) :: massout, energyout, en1, en2cp, en2cV, en3, encl
  real(rl),                  intent(in)    :: ptop
  logical,                   intent(inout) :: wasiactive

  real(rl) :: qsat, dp_loc, qv_loc, ql_loc, dpdry_loc, qvdry_loc, qldry_loc, qsatdry, dq_loc, &
                vapor_mass_change
  real(rl) :: T_loc, p_loc, pi_loc, L_old, L_new, rstar_old, rstar_new, hold, T_new
  real(rl) :: cpstarTerm_new, rstardp, oldQ1mass, oldQ2mass, olddphi, pi_top_int, enb, ena
  real(rl), dimension(nlev)  :: pi, dphi, zero, rain, d_pnh, rain2
  real(rl), dimension(nlevp) :: p_int
  integer  :: k, ll

  massout = 0.0
  energyout = 0.0
  zero = 0.0
  rain = 0.0
  en1 = 0; en2cp = 0; en2cV=0; en3 = 0; encl = 0;

  call construct_hydro_pressure(dp_c,ptop,pi)

  !compute en1, energy before condensation
  if(bubble_rj_cpstar_nh) then
    call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),p_c,en1)
  else
    call energycp_hy_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),    en1)
  endif


  !since one option uses geo, go from bottom to up
  do k=nlev, 1, -1
    !call qsat_rj2(p_c(k), T_c(k), qsat)
    call qsat_rj2(pi(k), T_c(k), qsat)

    if (qv_c(k) > qsat) then

       !condensation stage

       wasiactive = .true. 
       !compute dry values
       dp_loc = dp_c(k)
       qv_loc = qv_c(k)
       ql_loc = ql_c(k)
       dpdry_loc = dp_loc * (1.0 - qv_loc - ql_loc)
       qvdry_loc = qv_loc * dp_loc / dpdry_loc
       qldry_loc = ql_loc * dp_loc / dpdry_loc
       qsatdry = qsat * dp_loc / dpdry_loc     ! qv_new
       dq_loc = qvdry_loc - qsatdry            ! > 0 , is qliq_dry_new
       vapor_mass_change = dpdry_loc*dq_loc
       T_loc = T_c(k)
       p_loc = p_c(k)
       pi_loc = pi(k)

       !new Q will be qsatdry
       L_old = (latvap + latice) * qvdry_loc + latice * qldry_loc
       L_new = (latvap + latice) * qsatdry   + latice * (qldry_loc + dq_loc)

       rstar_old = rdry * 1.0 + rvapor * qvdry_loc
       rstar_new = rdry * 1.0 + rvapor * qsatdry

       cpstarTerm_new = cpdry*1.0 + cpv*qsatdry + cl*(qldry_loc + dq_loc)

       !use extra term in dh=0 rule in case of NH or not
       if(bubble_rj_cpstar_nh) then

       ! L and Rstar are in terms of q, so, enthalpy too
         hold = T_loc*( cpdry*1.0 + cpv*qvdry_loc + cl*qldry_loc) + &
             (pi_loc/p_loc - 1) * rstar_old * T_loc + L_old
         T_new = (hold - L_new)/ &
         (   cpstarTerm_new + ( pi_loc/p_loc - 1)*rstar_new   )

       elseif (bubble_rj_cpstar_hy) then
         hold = T_loc*( cpdry*1.0 + cpv*qvdry_loc + cl*qldry_loc ) + L_old
         T_new = (hold - L_new)/(   cpstarTerm_new   )
       endif

       T_c(k)  = T_new
       rain(k) = vapor_mass_change
       !this is before raining out, so dp did not change
       qv_c(k) = (dp_loc*qv_loc - vapor_mass_change)/dp_loc
       ql_c(k) = (dp_loc*ql_loc + vapor_mass_change)/dp_loc
 
    endif
  enddo

  !recompute geop
  do k=nlev, 1, -1
    rstar_new = rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k)
    olddphi = rstar_new * dp_c(k) * T_c(k) / p_c(k)
    zi_c(k) = zi_c(k+1) + olddphi/gravit
  enddo

!if(wasiactive)then
!print *, 'ql_c', ql_c
!compute en2, energy before sedimentation
!  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),p_c,en2)
!print *, 'en2 via P',en2
!  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c       ,p_c,en2)
!print *, 'en2 via V',en2
!endif

  if(bubble_rj_cpstar_nh) then
    !sanity check is to use the opposite formulation
    call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),p_c,en2cp)
    call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,en2cV)
  else
    !for HY update, there is no version of energy functional with cV
    call energycp_hy_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),    en2cp)
  endif

  ! in HY update, what to do about energyout and how/when to recompute p and phi?
  ! decision: keep both updates the same for cpstar HY and NH

  if(.not.bubble_rj_nosedim)then

  do k=nlev, 1, -1

     !compute current level of energy
     !we can save on perf by reusing 1 of these calls in next k iteration

     !rain2 contains what'd not been rained out yet
     call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,enb)
 
     !above, we only need to save dphi
     !compute dphi above the cell
     do ll=k-1,1,-1
       dphi(ll) = (zi_c(ll) - zi_c(ll+1))*gravit
     enddo

     !k cell, we need to recompute state vars and dphi
     !sedimentation stage
     oldQ1mass         = dp_c(k)*qv_c(k)
     vapor_mass_change = dp_c(k)*ql_c(k)
     dp_c(k) = dp_c(k) - vapor_mass_change
     p_c(k)  = p_c(k)  - vapor_mass_change
     qv_c(k) = oldQ1mass/dp_c(k)
     ql_c(k) = 0.0
     rstardp = dp_c(k) * (rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k))
     dphi(k) = rstardp * T_c(k) / p_c(k)

     !below, we need to recompute dphi from changed pressure
     do ll=k+1,nlev
        p_c(ll) = p_c(ll) - vapor_mass_change
        rstardp = dp_c(ll) * (rdry * (1.0 - qv_c(ll) - ql_c(ll)) + rvapor * qv_c(ll)) 
        dphi(ll) = rstardp * T_c(ll) / p_c(ll)
     enddo

     !recompute zi_c for all levels
     do ll=nlev, 1, -1
        zi_c(ll) = zi_c(ll+1) + dphi(ll)/gravit
     enddo

     call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,ena)
    
     massout = massout + vapor_mass_change
     energyout = energyout + (enb - ena)
     encl = encl + (Tforcl*cl+latice)*vapor_mass_change

  enddo
  endif

  !interchange V and P formulas of energy for debugging
  !call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),p_c,en3)
  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,en3)
 
end subroutine rj_new



subroutine rj_new_volume(qv_c,ql_c,T_c,dp_c,p_c,zi_c,ptop,massout,energyout,&
                         en1, en2, en3, encl, wasiactive)

  real(rl), dimension(nlev), intent(inout) :: p_c
  real(rl), dimension(nlev), intent(inout) :: dp_c
  real(rl), dimension(nlev), intent(inout) :: qv_c,ql_c,T_c
  real(rl), dimension(nlevp),intent(in)    :: zi_c
  real(rl),                  intent(inout) :: massout,energyout, en1, en2, en3, encl
  real(rl),                  intent(in)    :: ptop
  logical,                   intent(inout) :: wasiactive

  real(rl) :: qsat, dp_loc, qv_loc, ql_loc, dpdry_loc, qvdry_loc, qldry_loc, &
              qsatdry, dq_loc, vapor_mass_change
  real(rl) :: T_loc, p_loc, pi_loc, L_old, L_new, Iold, T_new
  real(rl) :: cvstarTerm_new, rstardp, oldQ1mass, oldQ2mass, olddphi
  real(rl), dimension(nlev) :: pi, zero, rain, ppp_c
  integer  :: k

  zero = 0.0
  massout = 0.0
  energyout = 0.0
  rain = 0.0
  en1 = 0; en2 = 0; en3 = 0; encl = 0;

  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en1)
!print *, 'from cV', en1
!  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),p_c,en1)
!print *, 'from cp', en1


  call construct_hydro_pressure(dp_c,ptop,pi)

  do k=1, nlev
    !call qsat_rj2(p_c(k), T_c(k), qsat)
    call qsat_rj2(pi(k), T_c(k), qsat)

    if (qv_c(k) > qsat) then
       !condensation stage

       wasiactive = .true.
       !compute dry values
       dp_loc = dp_c(k)
       qv_loc = qv_c(k)
       ql_loc = ql_c(k)
       dpdry_loc = dp_loc * (1.0 - qv_loc - ql_loc)
       qvdry_loc = qv_loc * dp_loc / dpdry_loc
       qldry_loc = ql_loc * dp_loc / dpdry_loc
       qsatdry = qsat * dp_loc / dpdry_loc     ! qv_new
       dq_loc = qvdry_loc - qsatdry            ! > 0 , is qliq_dry_new
       vapor_mass_change = dpdry_loc*dq_loc
       T_loc = T_c(k)
       p_loc = p_c(k)
       pi_loc = pi(k)

       !new Q will be qsatdry
       L_old = (latvap + latice) * qvdry_loc + latice * qldry_loc
       L_new = (latvap + latice) * qsatdry   + latice * (qldry_loc + dq_loc)

       Iold = T_loc*( cvdry*1.0 + cvv*qvdry_loc + cl*qldry_loc) + L_old

       cvstarTerm_new = cvdry*1.0 + cvv*qsatdry + cl*(qldry_loc + dq_loc)
       T_new = ( Iold - L_new )/(   cvstarTerm_new   )

       T_c(k)  = T_new
       rain(k) = vapor_mass_change
       !this is before raining out, so dp did not change
       qv_c(k) = (dp_loc*qv_loc - vapor_mass_change)/dp_loc
       ql_c(k) = (dp_loc*ql_loc + vapor_mass_change)/dp_loc

       !recompute pressure here in case there is no sedim
       rstardp = dp_c(k)* (rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k))
       p_c(k)  = - rstardp * T_c(k) / (zi_c(k+1) - zi_c(k) ) / gravit
     endif
  enddo

!if(wasiactive)then
!print *, 'pressure after', p_c
!print *, 'p diff', p_c - ppp_c
!stop
!endif

  !call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c),dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en2)
  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c),dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),p_c,en2)

  if(.not.bubble_rj_nosedim)then

  !sedimentation stage
  do k=1, nlev
      
     !sedimentation stage
     oldQ1mass = dp_c(k)*qv_c(k)
     vapor_mass_change = dp_c(k)*ql_c(k)
     dp_c(k) = dp_c(k) - vapor_mass_change
     ql_c(k) = 0.0
     qv_c(k) = oldQ1mass/dp_c(k)

     rstardp = dp_c(k)* (rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k))
     p_c(k)  = - rstardp * T_c(k) / (zi_c(k+1) - zi_c(k) ) / gravit

     massout = massout + vapor_mass_change
     energyout = energyout + vapor_mass_change*(T_c(k)*cl + latice + (zi_c(k) + zi_c(k+1))*gravit/2.0)
     !encl      = encl + vapor_mass_change*(T_c(k)*cl + latice)
     encl = encl + (Tforcl*cl+latice)*vapor_mass_change
  enddo
  endif

  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c),dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en3)

end subroutine rj_new_volume



! ABANDONED in this branch
subroutine rj_old(qv_c,T_c,dp_c,p_c,zi_c,ptop,massout,wasiactive)

  real(rl), dimension(nlev), intent(inout) :: p_c 
  real(rl), dimension(nlev), intent(inout) :: dp_c
  real(rl), dimension(nlev), intent(inout) :: qv_c,T_c
  real(rl), dimension(nlevp),intent(inout) :: zi_c
  real(rl),                  intent(in)    :: ptop
  real(rl),                  intent(inout) :: massout
  logical,                   intent(inout) :: wasiactive

  real(rl) :: qsat, dq_loc, rstardp, oldQ1mass, vapor_mass_lost
  integer :: k
  real(rl), dimension(nlev) :: pi
  real(rl), dimension(nlev) :: dphi

  ! this code is wrt wet mmrs!
  ! this code is not like old RJ implementation in part because sedimentation will not be done
  ! with const pressure assumption

  massout = 0.0
  
  !should be negative
  dphi(1:nlev) = gravit*(zi_c(1:nlev) - zi_c(2:nlevp))

  call construct_hydro_pressure(dp_c,ptop,pi)

  do k=1, nlev
    !call qsat_rj2(p_c(k), T_c(k), qsat)
    call qsat_rj2(pi(k), T_c(k), qsat)

    if (qv_c(k) > qsat) then
       wasiactive = .true.
!original RJ

       ! condensation stage, const dp, const pnh

#if 0
!old version with 1st order qsat substitution
       dq_loc = (qv_c(k) - qsat) &
         / (1.0 + (latvap/cpdry) * bubble_epsilo * latvap * qsat / (rdry*T_c(k)*T_c(k)))
#else
!version that matches the rest of routines, no 1st order approx
       dq_loc  = qv_c(k) - qsat
#endif
       T_c(k)  = T_c(k)  + latvap / cpdry * dq_loc
       qv_c(k) = qv_c(k) - dq_loc
       vapor_mass_lost = dq_loc * dp_c(k)

       !dummy ql_c = dq_loc


       ! sedimentation
       ! we recompute dp here after the fact to avoid calling _cam functions
       ! we also recomopute everything else (in this case, after dp, pnh and dphi change)
       !dummy ql_c = 0
       !mass_lost is dq_loc * dp_c(k)
       !new dp is dp -= mass_lost
     
       oldQ1mass = dp_c(k) * qv_c(k)
       massout = massout + vapor_mass_lost
       dp_c(k) = dp_c(k) - vapor_mass_lost
       p_c(k)  = p_c(k)  - vapor_mass_lost
       qv_c(k) = oldQ1mass / dp_c(k)

       !dphi = - dp * Rstar * T / pnh
       rstardp = dp_c(k) * (rdry * (1.0 - qv_c(k)) + rvapor * qv_c(k))
       dphi(k) = rstardp * T_c(k) / p_c(k)

     endif
  enddo

  !recompute zi
  do k=nlev, 1, -1
     zi_c(k) = zi_c(k+1) + dphi(k)/gravit
  enddo


end subroutine rj_old




subroutine energycp_hy_via_mass(dpdry_c,dpv_c,dpc_c,dpr_c,T_c,ptop,zbottom,energy)

  real(rl), dimension(nlev), intent(in) :: dpdry_c, dpv_c, dpc_c, dpr_c
  real(rl), dimension(nlev), intent(in) :: T_c
  real(rl),                  intent(in) :: ptop, zbottom
  real(rl),                  intent(inout) :: energy

  real(rl) :: pis, cpterm, Lterm
  integer  :: k

  pis = ptop + sum(dpdry_c + dpv_c + dpc_c + dpr_c)

  energy = zbottom * pis * gravit

  do k=1,nlev
    cpterm = cpdry*dpdry_c(k) + cpv * dpv_c(k) + cl * (dpc_c(k) + dpr_c(k))
    Lterm = (latvap+latice) * dpv_c(k) + latice * (dpc_c(k) + dpr_c(k))
    energy = energy + T_c(k)*cpterm + Lterm 
  enddo

end subroutine energycp_hy_via_mass




subroutine energycp_nh_via_mass(dpdry_c,dpv_c,dpc_c,dpr_c,T_c,ptop,zbottom,p_c,energy)

  real(rl), dimension(nlev), intent(in) :: dpdry_c, dpv_c, dpc_c, dpr_c, p_c
  real(rl), dimension(nlev), intent(in) :: T_c
  real(rl),                  intent(in) :: ptop, zbottom
  real(rl),                  intent(inout) :: energy

  real(rl) :: pis, cpterm, Lterm, nhterm, rstar
  integer  :: k
  real(rl), dimension(nlev) :: ppi
  real(rl), dimension(nlev) :: pnh, dpi

  dpi = dpdry_c + dpv_c + dpc_c + dpr_c
  pis = ptop + sum(dpi)

  !derived pressure values, on midlevels
  call construct_hydro_pressure(dpi,ptop,ppi)
  pnh = p_c

  energy = zbottom * pis * gravit

  do k=1,nlev
    cpterm = cpdry*dpdry_c(k) + cpv * dpv_c(k) + cl * (dpc_c(k) + dpr_c(k))

    Lterm  = (latvap+latice) * dpv_c(k) + latice * (dpc_c(k) + dpr_c(k))

    rstar = rdry*dpdry_c(k)/dpi(k) + rvapor*dpv_c(k)/dpi(k)

    nhterm = rstar*T_c(k)*(ppi(k)/pnh(k)-1.0)*dpi(k)

    energy = energy + T_c(k)*cpterm + Lterm + nhterm
  enddo

end subroutine energycp_nh_via_mass


subroutine energycV_nh_via_mass(dpdry_c,dpv_c,dpc_c,dpr_c,T_c,ptop,zi_c,p_c,energy)

  real(rl), dimension(nlev), intent(in) :: dpdry_c, dpv_c, dpc_c, dpr_c, p_c
  real(rl), dimension(nlevp),intent(in) :: zi_c
  real(rl), dimension(nlev), intent(in) :: T_c
  real(rl),                  intent(in) :: ptop
  real(rl),                  intent(inout) :: energy

  real(rl) :: ztop, cvterm, Lterm
  integer  :: k
  real(rl), dimension(nlev) :: zm, dpi

  ztop = zi_c(1)
  dpi = dpdry_c + dpv_c + dpc_c + dpr_c
  zm(1:nlev) = (zi_c(1:nlev) + zi_c(2:nlevp))/2.0

  energy = ztop * ptop * gravit

  do k=1,nlev
    cvterm = cvdry*dpdry_c(k) + cvv * dpv_c(k) + cl * (dpc_c(k) + dpr_c(k))

    Lterm  = (latvap+latice) * dpv_c(k) + latice * (dpc_c(k) + dpr_c(k))

    energy = energy + T_c(k)*cvterm + zm(k)*gravit*dpi(k) + Lterm
  enddo

end subroutine energycV_nh_via_mass


!mu update formula, not used here
subroutine energycp_nh_via_mass_nhupdate(dpdry_c,dpv_c,dpc_c,dpr_c,T_c,ptop,zi_c,p_c,energy)

  real(rl), dimension(nlev), intent(in) :: dpdry_c, dpv_c, dpc_c, dpr_c, p_c
  real(rl), dimension(nlevp),intent(in) :: zi_c
  real(rl), dimension(nlev), intent(in) :: T_c
  real(rl),                  intent(in) :: ptop
  real(rl),                  intent(inout) :: energy

  real(rl) :: zbottom, cpterm, Lterm, factor
  integer  :: k
  real(rl), dimension(nlev) :: zm, dpi
  real(rl), dimension(nlevp):: p_int

  zbottom = zi_c(nlevp)

  dpi = dpdry_c + dpv_c + dpc_c + dpr_c
  zm(1:nlev) = (zi_c(1:nlev) + zi_c(2:nlevp))/2.0

  p_int(2:nlev) = (p_c(1:nlev-1) + p_c(2:nlev))/2.0
  p_int(1) = p_c(1) - dpi(1)/2.0
  p_int(nlevp) = p_c(nlev) + dpi(nlev)/2.0
  
  energy = zbottom * p_int(nlevp) * gravit

  do k=1,nlev
    factor = (1 - (p_int(k) - p_int(k+1)) /dpi(k) )
    cpterm = cpdry*dpdry_c(k) + cpv * dpv_c(k) + cl * (dpc_c(k) + dpr_c(k))
    Lterm  = (latvap+latice) * dpv_c(k) + latice * (dpc_c(k) + dpr_c(k))
    energy = energy + T_c(k)*cpterm + zm(k)*gravit*dpi(k)*factor + Lterm
  enddo

end subroutine energycp_nh_via_mass_nhupdate






!!!!!!!!!!! copy from dcmip16 
subroutine qsat_rj2(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_epsilo * bubble_e0 / p * &
         exp(-(latvap/rvapor) * ((1.0/T)-(1.0/bubble_t0_const)))
end subroutine qsat_rj2


end module newphysics
#endif
