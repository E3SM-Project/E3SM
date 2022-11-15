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
real(rl), parameter:: k1=0.001,k2=2.2,k3=0.875,a1=0.00001 !a1 should be 0.001

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

!print *, 'max of qc', maxval(qcdry)
!print *, k1,a1
!if(maxval(qcdry)> 0.001)then
!print *, 'hey you!!!!!!!!!!!'
!endif


  !auto-accumulation step Ar = max ( k1 (qc - a), 0 )
  do k=1,nlev
    change = dt * max (k1*(qcdry(k) - a1), 0.0)
    qcdry(k) = qcdry(k) - change
    qrdry(k) = qrdry(k) + change

!if (change > 0.0) then
!print *, 'hey rain happned!'
!stop
!endif
!if(qcdry(k)> 0.001)then
!print *, "change is ", change
!endif


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

!print *, 'positt - zm(k)', positt - zm(k)
!if(d_ind > kk) then
!print *, k,kk
!endif

!hack
!positt = -100.0
!hack
if (k < nlev) then
  positt = zm(k+1)
  d_ind = k+1
else
  positt = -100.0
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


subroutine recompute_pressures(qvdry,qcdry,qrdry,dpdry,pidry,pprime,pi,pnh,dp)

  real(rl), dimension(nlev), intent(in)    :: qvdry,qcdry,qrdry
  real(rl), dimension(nlev), intent(in)    :: dpdry,pidry,pprime
  real(rl), dimension(nlev), intent(out)   :: dp,pi,pnh

  real(rl)  :: ptop
  real(rl), dimension(nlev)  :: dpi_from_water, pi_from_water

  ptop = 0.0

  dpi_from_water = dpdry*( qvdry + qcdry + qrdry )
  dp = dpi_from_water + dpdry

  call construct_hydro_pressure(dpi_from_water,ptop,pi_from_water)

  pi = pidry + pi_from_water
  pnh = pi + pprime

end subroutine recompute_pressures


subroutine condensation_and_back_again(qvdry,qcdry,qrdry,tempe,dpdry,dp,pi,pnh,dt)

  real(rl), dimension(nlev), intent(in)    :: qrdry
  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,tempe
  real(rl), dimension(nlev), intent(in)    :: dpdry, dp, pi, pnh
  real(rl),                  intent(in)    :: dt

  integer  :: k
  real(rl) :: dq, qsat, qsatdry, T_new, cptermold, cptermnew, Lold, Lnew

  do k=1, nlev
    call qsat_kessler2(pnh(k), tempe(k), qsat)
    qsatdry = qsat*dp(k)/dpdry(k)

    !assume condensation
    dq = qvdry(k) - qsatdry
    !except if
    if ( dq < 0.0 ) then
      !evaporation
      dq = max( dq, -qcdry(k) )
    endif
    ! new values: qv = qv - dq, qc = qc + dq

    cptermold = cpdry + cpv *  qvdry(k)       + cl * (qcdry(k) + qrdry(k))
    cptermnew = cpdry + cpv * (qvdry(k) - dq) + cl * (qcdry(k) + qrdry(k) + dq)

    Lold = (latvap+latice) * qvdry(k)     + latice * (qcdry(k) + qrdry(k))
    Lnew = (latvap+latice) *(qvdry(k)-dq) + latice * (qcdry(k) + qrdry(k) + dq)

    !dpdry is cancelled from both sides
    T_new = ( tempe(k) * cptermold + Lold - Lnew ) / cptermnew 

    tempe(k) = T_new

    qvdry(k) = qvdry(k) - dq
    qcdry(k) = qcdry(k) + dq

!if (dq > 0) then
! print *, 'yay, condensation',' qc', qcdry(k)
!endif

  enddo

end subroutine condensation_and_back_again


subroutine compute_energy(qvdry,qcdry,qrdry,tempe,dpdry,pi,zbottom,energy)

  real(rl), dimension(nlev), intent(in)    :: qrdry
  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,tempe
  real(rl), dimension(nlev), intent(in)    :: dpdry, pi
  real(rl),                  intent(inout) :: energy
  real(rl),                  intent(in)    :: zbottom

  integer  :: k
  real(rl) :: cpterm, Lterm

  !dont do 1/g term
  energy = zbottom*pi(nlevp)
  do k=1,nlev
    cpterm = cpdry + cpv * qvdry(k) + cl * (qcdry(k) + qrdry(k))
    Lterm = (latvap+latice) * qvdry(k) + latice * (qcdry(k) + qrdry(k))
    energy = energy + dpdry(k)*( tempe(k)*cpterm + Lterm )
  enddo

end subroutine compute_energy









!!!!!!!!!!! copy from dcmip16 FIX THAT!
subroutine qsat_kessler2(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_const1 / p * exp( bubble_const2 * (T - bubble_const3) / ( T - bubble_const4 ) )
end subroutine qsat_kessler2




end module newphysics
#endif
