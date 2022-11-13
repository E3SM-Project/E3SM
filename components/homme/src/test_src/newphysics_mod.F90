#ifndef CAM
#include "config.h"

module newphysics

! Implementation of the dcmip2012 dycore tests for the preqx dynamics target

use control_mod,          only: theta_hydrostatic_mode,&
                   case_planar_bubble, bubble_prec_type, bubble_rj_cpdry, bubble_rj_cpstar
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
                                rdry=>rgas, cpdry=>cp, cpv=>cpwater_vapor, cl, &
                                rvapor=>rwater_vapor, rhow

implicit none

!kessler constants, put into physonst mod
real(rl), parameter:: k1=0.001,k2=2.2,k3=0.875,a1=0.001

contains

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

subroutine accrecion_and_accumulation(qcdry,qrdry,dt)

  real(rl), dimension(nlev), intent(inout) :: qcdry,qrdry
  real(rl),                  intent(in)    :: dt
  integer  :: k
  real(rl) :: change

  !no movement between levels, no T change
  !accretion step, collection Cr = max ( k2 qc qr^k3, 0)
  do k=1,nlev
    change = dt * max ( k2*qcdry(k)*qrdry(k)**k3, 0.0)
    qcdry(k) = qcdry(k) - change
    qrdry(k) = qrdry(k) + change
  enddo

  !auto-accumulation step Ar = max ( k1 (qc - a), 0 )
  do k=1,nlev
    change = dt * max (k1*(qcdry(k) - a1), 0.0)
    qcdry(k) = qcdry(k) - change
    qrdry(k) = qrdry(k) + change
  enddo

end subroutine accrecion_and_accumulation

subroutine get_geo_from_drydp(Tempe, dpdry, pidry, zbottom, zi, zm, dz)

  real(rl), dimension(nlevp), intent(out) :: zi
  real(rl), dimension(nlev),  intent(out) :: zm, dz
  real(rl), dimension(nlev),  intent(in)  :: Tempe, dpdry, pidry
  real(rl),                   intent(in)  :: zbottom

  integer  :: k

  zi(nlevp) = zbottom
  do k=nlev,1,-1
    zi(k) = zi(k+1) + rdry*Tempe(k)*dpdry(k)/pidry(k)/gravit
  enddo
  dz(1:nlev) = zi(1:nlev) - zi(2:nlevp)
  zm(1:nlev) = (zi(1:nlev) + zi(2:nlevp))/2

end subroutine get_geo_from_drydp


subroutine sedimentation(qvdry,qcdry,qrdry,tempe,dpdry,pidry,zbottom,massleft,energyleft,dt)

  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,qrdry,tempe
  real(rl), dimension(nlev), intent(in)    :: dpdry, pidry
  real(rl),                  intent(in)    :: zbottom,dt
  real(rl),                  intent(out)   :: massleft, energyleft

  integer  :: k, kk, d_ind
  real(rl) :: change, positt, velqr, T_new, cptermold, cptermnew
  real(rl), dimension(nlev)  :: zm, dz, rhodry
  real(rl), dimension(nlevp) :: zi

  massleft = 0.0; energyleft = 0.0

  call get_geo_from_drydp(tempe, dpdry, pidry, zbottom, zi, zm, dz)

  rhodry = dpdry / dz / gravit

  do k=1,nlev

    if(qrdry(k)>0.0) then
      ! Liquid water terminal velocity, using dry densities instead of wet
      ! in meter/sec
      velqr  = 36.34d0*(qrdry(k)*rhodry(k))**0.1364*sqrt(rhodry(nlev)/rhodry(k))

!print *, 'k, velqr, qr', k, velqr(k), qrdry(k)

      !check where it ends
       d_ind = k
       positt = zm(k) - velqr * dt
       !cell with kk index has boundaries zi(kk) and zi(kk+1)
       do kk=k+1,nlev
         if ( (positt > zi(kk)) .and. (positt <= zi(kk+1)) ) then
           d_ind = kk
           !kk = nlev ! does F allow this
         endif
       enddo

print *, 'positt - zm(k)', positt - zm(k)
if(d_ind > kk) then
print *, k,kk
endif

       if(positt < zi(nlevp)) then
         !hit the bottom
         d_ind = -1
         massleft   = massleft   + qrdry(k) * dpdry(k)
!add nh term later
         energyleft = energyleft + qrdry(k) * dpdry(k) * ( cl*tempe(k) + latice )
         qrdry(k) = 0.0
       else
         !stuck in between
         !arrives to dest cell, recompute ratio for the dest cell
         change = qrdry(k)*dpdry(k)/dpdry(d_ind)
         !compute average temperature with new mass
         cptermold = cpdry + cpv * qvdry(d_ind) + cl * (qcdry(d_ind) + qrdry(d_ind))
         cptermnew = cptermold + cl*change
!add nh term later
         T_new = ( tempe(d_ind) * dpdry(d_ind) * cptermold + tempe(k) * dpdry(k) * qrdry(k) * cl ) / &
                  cptermnew / dpdry(d_ind)

!print *, 'T_new', T_new 
         tempe(d_ind) = T_new

         qrdry(k) = 0.0
         qrdry(d_ind) = qrdry(d_ind) + change
       endif
     endif !there was rain in cell
   end do ! k loop for sedim

end subroutine sedimentation





end module newphysics
#endif
