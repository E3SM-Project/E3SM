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

!use physical_constants,   only:  pi=>dd_pi

use physical_constants,   only: bubble_const1, bubble_const2, bubble_const3, bubble_const4, &
                                bubble_t0_const, bubble_epsilo, bubble_e0, &
                                latvap, latice, &
                                gravit=>g, p0, &
                                rdry=>rgas, cpdry=>cp, cpv=>cpwater_vapor, cl, &
                                rvapor=>rwater_vapor, rhow

implicit none

!kessler constants, put into phys_const mod
real(rl), parameter:: k1=0.001,k2=2.2,k3=0.875,a1=0.001

contains

subroutine construct_hydro_pressure(dp,ptop,pi)

  real(rl), dimension(nlev), intend(in)  :: dp
  real(rl),                  intend(in)  :: ptop
  real(rl), dimension(nlev), intend(out) :: p
  integer :: k
  real(rl) :: pi_upper

  pi_upper = ptop
  do k = 1, nlev
    pi_upper = pi_upper + dp_c(k)
    pi(k) = pi_upper - dp_c(k)/2
  enddo

end subroutine construct_hydro_pressure

subroutine convert_to_dry(q1,q2,q3,dp,dpdry,q1dry,q2dry,q3dry)

  real(rl), dimension(nlev), intend(in)  :: q1,q2,q3,dp,dpdry
  real(rl), dimension(nlev), intend(out) :: q1dry,q2dry,q3dry
  q1dry = q1 * dp /dpdry
  q2dry = q2 * dp /dpdry
  q3dry = q3 * dp /dpdry

end subroutine convert_to_dry

subroutine convert_to_wet(q1dry,q2dry,q3dry,dp,dpdry,q1,q2,q3)

  real(rl), dimension(nlev), intend(in)  :: q1dry,q2dry,q3dry,dp,dpdry
  real(rl), dimension(nlev), intend(out) :: q1,q2,q3
  q1 = q1dry * dpdry / dp
  q2 = q2dry * dpdry / dp
  q3 = q3dry * dpdry / dp

end subroutine convert_to_dry

subroutine accrecion_and_accumulation(qcdry,qrdry,dt)

  real(rl), dimension(nlev), intend(inout) :: qcdry,qrdry
  real(rl),                  intend(in)    :: dt
  integer  :: k
  real(rl) :: change

  !no movement between levels, no T change
  !accretion step, collection Cr = max ( k2 qc qr^k3, 0)
  do k=1,nlev
    change = dt * max ( k2*qcdry_c(k)*qrdry(k)**k3, 0.0)
    qcdry(k) = qcdry(k) - change
    qrdry(k) = qrdry(k) + change
  enddo

  !auto-accumulation step Ar = max ( k1 (qc - a), 0 )
  do k=1,nlev
    change = dt * max (k1*(qcdry(k) - a1), 0.0)
    qcdry(k) = qcdry(k) - change
    qrdry(k) = qrdry(k) + change
  enddo

end subroutine accrecion

subroutine get_geo_from_drydp(Tempe, dpdry, pidry, zbottom, zi, zm, dz)

  real(rl), dimension(nlevp), intend(out) :: zi
  real(rl), dimension(nlev),  intend(out) :: zm, dz
  real(rl), dimension(nlev),  intend(in)  :: Tempe, dpdry, pidry
  real(rl),                   intend(in)  :: zbottom

  integer  :: k

  zi(nlevp) = zbottom
  do k=nlev,1,-1
    zi(k) = zi(k+1) + rdry*Tempe(k)*dpdry(k)/pidry(k)/gravit
  enddo
  dz(1:nlev) = zi(1:nlev) - zi(2:nlevp)
  zm(1:nlev) = (zi(1:nlev) + zi(2:nlevp))/2

end subroutine get_geo_from_drydp


subroutine sedimentation(qvdry,qcdry,qrdry,tempe,dpdry,pidry,zi,massleft,energyleft,dt)

  real(rl), dimension(nlev), intend(inout) :: qvdry,qcdry,qrdry,tempe
  real(rl), dimension(nlev), intend(in)    :: dpdry, pidry, zi
  real(rl),                  intend(in)    :: dt
  real(rl),                  intend(out)   :: massleft, energyleft

  integer  :: k, kk, d_ind
  real(rl) :: change, positt, velqr, T_new, cptermold, cptermnew
  real(rl), dimension(nlev) :: zm, dz, rhodry

  massleft = 0.0; energyleft = 0.0

  call get_geo_from_drydp(tempe, dpdry, pidry, zbottom, zi, zm, dz)

  rhodry = dpdry / dz / gravit

  do k=1,nlev

    if(qrdry(k)>0.0) then
      ! Liquid water terminal velocity, using dry densities instead of wet
      ! in meter/sec
      velqr  = 36.34d0*(qrdry(k)*rhodry(k))**0.1364*sqrt(rhodry(nlev)/rhodry(k))

!print *, 'k, velqr, qr', k, velqr(k), qrdry_c(k)

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

print *, 'positt - zm_c(k)', positt - zm_c(k)
if(d_ind > kk) then
print *, k,kk
endif

       if(positt < zi(nlevp)) then
         !hit the bottom
         d_ind = -1
         massleft   = massleft   + qrdry(k) * dpdry(k)
         energyleft = energyleft + qrdry(k) * dpdry(k) * ( cl*tempe(k) + latice )
         qrdry_c(k) = 0.0
       else
         !stuck in between
         !arrives to dest cell, recompute ratio for the dest cell
         change = qrdry(k)*dpdry(k)/dpdry(d_ind)
         !compute average temperature with new mass
         cptermold = cpdry + cpwater_vapor * qvdry_c(d_ind) + cl * (qcdry_c(d_ind) + qrdry_c(d_ind))
         cptermnew = cptermold + cl*change
         T_new = ( tempe(d_ind) * dpdry(d_ind) * cptermold + tempe(k) * dpdry(k) * qrdry(k) * cl ) / &
                  cptermnew / dpdry_c(d_ind)

!print *, 'T_new', T_new 
         tempe(d_ind) = T_new

         qrdry(k) = 0.0
         qrdry(d_ind) = qrdry(d_ind) + change
       endif
     endif !there was rain in cell
   end do ! k loop for sedim

end subroutine sedimentation





















!_______________________________________________________________________
subroutine bubble_new_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer :: i,j,k,kk,ie,qind,d_ind
  real(rl), dimension(np,np,nlev) :: u,v,w,T,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: T0,qv0,qc0,qr0
  real(rl), dimension(np,np)      :: ps 
  real(rl), dimension(nlev)       :: p_c,qv_c,qc_c,qr_c,dp_c,T_c,pi,dz_c,rhodry,velqr
  real(rl), dimension(nlev)       :: dpdry_c, pidry, zm_c, qvdry_c, qcdry_c, qrdry_c, change_c
  real(rl) :: max_w, max_precl, min_ps

  real(rl) :: pi_upper,qsat,qsatdry,dp_loc,qv_loc,dpdry_loc,qvdry_loc,dq_loc,vapor_mass_change
  real(rl) :: T_loc,p_loc,pi_loc,L_old,L_new,rstar_old,rstar_new,hold,cpstarTerm_new
  real(rl) :: T_new, change, pidry_upper, mass_prect, energy_prect, positt

  real(rl) :: zi(np,np,nlevp), zi_c(nlevp)
!  integer, parameter :: test = 1

!  real(rl), parameter:: gravit = 9.80616, rair = 287.0, cpair = 1.0045e3, cpv = 1810.0, cl = 4188.0, &
!                        rvapor = 461.5, epsilo = rair/rvapor, &
!                        latvap=2.501e6, latice=3.337e5, e0=610.78, &
!                        T0const = 273.16, rhow = 1000.0, &

  !kessler constants
  real(rl), parameter:: k1=0.001,k2=2.2,k3=0.875,a1=0.001


  if (qsize .ne. 3) call abortmp('ERROR: moist bubble test requires qsize=3')

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    elem(ie)%derived%FM(:,:,1,:) = 0.0 !(u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = 0.0 !(v - v0)/dt

    !precl(:,:,ie) = -1.0d0
    precl(:,:,ie) = 0.0d0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    ! use qind to avoid compiler warnings when qsize_d<5
    qind=1;  qv  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=2;  qc  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp
    qind=3;  qr  = elem(ie)%state%Qdp(:,:,:,qind,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere

    ! save un-forced prognostics
    T0=T; qv0=qv; qc0=qc; qr0=qr;

    ! apply forcing to columns
    do j=1,np; do i=1,np

      !column values
      qv_c = qv(i,j,:); qc_c = qc(i,j,:); qr_c = qr(i,j,:); 
      p_c  = p(i,j,:); dp_c = dp(i,j,:); T_c = T(i,j,:);

      pi_upper = hvcoord%hyai(1) * hvcoord%ps0
      do k = 1, nlev
        pi_upper = pi_upper + dp_c(k)
        pi(k) = pi_upper - dp_c(k)/2
      enddo


      ! RJ part
      if(bubble_prec_type == 1) then

      do k=1, nlev
        call qsat_rj(p_c(k), T_c(k), qsat)

        if (qv_c(k) > qsat) then

!print *, 'HEY RAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!stop
        !condensation happens
        if(bubble_rj_cpstar) then
        !compute dry values
          dp_loc = dp_c(k)
          qv_loc = qv_c(k)
          dpdry_loc = dp_loc * (1.0 - qv_loc)
          qvdry_loc = qv_loc * dp_loc / dpdry_loc
          qsatdry = qsat * dp_loc / dpdry_loc     ! qv_new
          dq_loc = qvdry_loc - qsatdry            ! > 0 , is qliq_dry_new
          vapor_mass_change = dpdry_loc*dq_loc
          T_loc = T_c(k)
          p_loc = p_c(k)
          pi_loc = pi(k)

          !new Q will be qsatdry
          L_old = (latvap + latice) * qvdry_loc
          L_new = (latvap + latice) * qsatdry   + latice * dq_loc

          rstar_old = Rgas * 1.0 + Rwater_vapor * qvdry_loc
          rstar_new = Rgas * 1.0 + Rwater_vapor * qsatdry

!#define HYY
#undef HYY
#ifndef HYY
          ! L and Rstar are in terms of q, so, enthalpy too
          hold = T_loc*( cp*1.0 + cpwater_vapor*qvdry_loc ) + &
               (pi_loc/p_loc - 1) * rstar_old * T_loc + L_old
#else
          hold = T_loc*( cp*1.0 + cpwater_vapor*qvdry_loc ) + &
                                                        L_old
#endif
          cpstarTerm_new = cp*1.0 + cpwater_vapor*qsatdry + cl*dq_loc
          !hnew = T_new * (   cpstarTerm_new + ( (pi_loc - vapor_mass_change)/(p_loc - vapor_mass_change) - 1)*rstar_new   ) + &
          !       L_new
          !     = hold
#ifndef HYY
          T_new = (hold - L_new)/ &
          (   cpstarTerm_new + ( pi_loc/p_loc - 1)*rstar_new   )  
#else
          T_new = (hold - L_new)/ &
          (   cpstarTerm_new   )
#endif
print *, 'T_new - T_c', T_new - T_c(k)

#if 0
if( (T_new - T_c(k)) < 0 )then
print *, 'dq_loc', dq_loc
print *, 'T_loc, T_c', T_loc, T_c(k)
print *, 'qvdry_loc', qvdry_loc
print *, 'dq_loc', dq_loc
print *, 'dq_loc', dq_loc
stop
endif
#endif

          !this does not need conversion wet-dry
          T_c(k)  = T_new
          qv_c(k) = qsat
          precl(i,j,ie) = precl(i,j,ie) + vapor_mass_change / dt / rhow / g
print *, 'precl',  precl(i,j,ie), vapor_mass_change, rhow

        elseif(bubble_rj_cpdry) then
!original RJ
          dq_loc = (qv_c(k) - qsat) &
            / (1.0 + (latvap/cp) * bubble_epsilo * latvap * qsat / (Rgas*T_c(k)*T_c(k)))
          T_c(k) = T_c(k) + latvap / cp * dq_loc
          qv_c(k) = qv_c(k) - dq_loc
          precl(i,j,ie) = precl(i,j,ie) + dq_loc * dp_c(k) / (dt * rhow) / g
        endif

      endif ! if prect happened
      enddo ! k loop

      ! Kessler part
      elseif(bubble_prec_type == 0) then


#if 1

!available quantities
      !column values
!      qv_c, qc_c, qr_c
!      p_c, dp_c, T_c, pi

!unlike RJ, Kessler needs all column values at once (for sedimentation, at least)
!and to split all processes
!and the same temp vars that were used for RJ do not work here

if(theta_hydrostatic_mode) then
print *, 'A! switch to NH'
stop
endif

        if( any(qv_c>0).or.any(qc_c>0).or.any(qr_c>0) ) then

        !dry density, the only quantity that won't change
        dpdry_c = dp_c*(1.0 - qv_c - qc_c - qr_c)
        pidry_upper = hvcoord%hyai(1) * hvcoord%ps0
        do k = 1, nlev
          pidry_upper = pidry_upper + dpdry_c(k)
          pidry(k) = pidry_upper - dpdry_c(k)/2
        enddo

        !convert all tracers to dry ratios
        qvdry_c = qv_c * dp_c / dpdry_c
        qcdry_c = qc_c * dp_c / dpdry_c
        qrdry_c = qr_c * dp_c / dpdry_c

#if 1
        !no movement between levels, no T change
        !accretion step, collection Cr = max ( k2 qc qr^k3, 0)
        do k=1,nlev
        change_c(k) = dt * max ( k2*qcdry_c(k)*qrdry_c(k)**k3, 0.0)
        enddo
        !apply change
        qcdry_c = qcdry_c - change_c
        qrdry_c = qrdry_c + change_c

        !auto-accumulation step Ar = max ( k1 (qc - a), 0 )
        do k=1,nlev
        change_c(k) = dt * max (k1*(qcdry_c(k) - a1), 0.0)
!if(change_c(k)>0) then
!print *, 'change_c', change_c(k)
!endif  

        enddo

!print *, 'max of qcdry', maxval(qcdry_c)
        !apply change
        qcdry_c = qcdry_c - change_c
        qrdry_c = qrdry_c + change_c
#endif

!print *, 'rain', qrdry_c

        !sedimentation needs rho, which needs dz
        zi_c(nlevp) = zi(i,j,nlevp)
        do k=nlev,1,-1
          zi_c(k) = zi_c(k+1) + Rgas*T_c(k)*dpdry_c(k)/pidry(k)/g
        enddo
        dz_c(1:nlev) = zi_c(1:nlev) - zi_c(2:nlevp)
        zm_c(1:nlev) = (zi_c(1:nlev) + zi_c(2:nlevp))/2
        
        rhodry = dpdry_c / dz_c / g
        mass_prect = 0.0; energy_prect = 0.0

        !sedimentation
        !
        !compute velocity of rain
        do k=1,nlev

          if(qrdry_c(k)>0.0) then
          ! Liquid water terminal velocity, using dry densities instead of wet
          ! in meter/sec
          velqr(k)  = 36.34d0*(qrdry_c(k)*rhodry(k))**0.1364*sqrt(rhodry(nlev)/rhodry(k))

!print *, 'k, velqr, qr', k, velqr(k), qrdry_c(k)

          !check where it ends
          d_ind = k
          positt = zm_c(k) - velqr(k) * dt
          !cell with kk index has boundaries zi(kk) and zi(kk+1)
          do kk=k+1,nlev
            if ( (positt > zi_c(kk)) .and. (positt <= zi_c(kk+1)) ) then
               d_ind = kk
               !kk = nlev ! does F allow this
            endif   
          enddo

print *, 'positt - zm_c(k)', positt - zm_c(k)

if(d_ind > kk) then
print *, k,kk
endif

          if(positt < zi_c(nlevp)) then 
            !hit the bottom
            d_ind = -1
            mass_prect = mass_prect + qrdry_c(k) * dpdry_c(k) 
            energy_prect =            qrdry_c(k) * dpdry_c(k) * ( cl*T_c(k) + latice ) 
            qrdry_c(k) = 0.0
          else
            !stuck in between
            !arrives to dest cell
            change = qrdry_c(k)
            !compute average temperature with new mass
            T_new = ( T_c(d_ind)*dpdry_c(d_ind)*(cp + cpwater_vapor * qvdry_c(d_ind) + cl * (qcdry_c(d_ind) + qrdry_c(d_ind)) ) + &
                     T_c(k)    *dpdry_c(k)    *cl*change ) / &
                     ( cp + cpwater_vapor * qvdry_c(d_ind) + cl * (qcdry_c(d_ind) + qrdry_c(d_ind) + change) ) / &
                     dpdry_c(d_ind)

!print *, 'T_new', T_new 

            T_c(d_ind) = T_new

            qrdry_c(k)     = qrdry_c(k)     - change
            qrdry_c(d_ind) = qrdry_c(d_ind) + change
          endif
          endif !there was rain in cell
        end do ! k loop for sedim

!print *, 'old-new', T(i,j,:)-T_c
!print *, 'new', T_c

        ! evaporation of rain
        ! condensation - evaporation

        !update q fields, 
        dp_c = dpdry_c*(1.0 + qvdry_c + qcdry_c + qrdry_c)
        qv_c = qvdry_c * dpdry_c / dp_c 
        qc_c = qcdry_c * dpdry_c / dp_c 
        qr_c = qrdry_c * dpdry_c / dp_c 
        precl(i,j,ie) = precl(i,j,ie) + mass_prect / (dt * rhow) / g

      endif ! any water >0

#endif
      endif ! RJ or Kessler choice

      !now update 3d fields here
      T(i,j,:)  = T_c(:)
      qv(i,j,:) = qv_c(:)
      qc(i,j,:) = qc_c(:)
      qr(i,j,:) = qr_c(:)

      ! set dynamics forcing
      elem(ie)%derived%FT(:,:,:)   = (T - T0)/dt

      ! set tracer-mass forcing. conserve tracer mass
      ! one call say, this is "rain" step, liquid water dissapears from the column
      elem(ie)%derived%FQ(i,j,:,:) = 0.0
      elem(ie)%derived%FQ(i,j,:,1) = (qv(i,j,:) - qv0(i,j,:))/dt !(rho_dry/rho)*dp*(qv-qv0)/dt
      elem(ie)%derived%FQ(i,j,:,2) = (qc(i,j,:) - qc0(i,j,:))/dt !(rho_dry/rho)*dp*(qv-qv0)/dt
      elem(ie)%derived%FQ(i,j,:,3) = (qr(i,j,:) - qr0(i,j,:))/dt !(rho_dry/rho)*dp*(qv-qv0)/dt
 

    enddo; enddo; !j,i loop

    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo !ie loop

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine bubble_new_forcing











subroutine toy_init(rcd)
  real(rl), intent(inout) :: rcd(6)
  rcd(1) = 1; rcd(2) = -1; rcd(3) = 1; rcd(4) = -1; rcd(5) = 1; rcd(6) = -1
end subroutine toy_init

subroutine toy_rcd(q, rcd)
  real(rl), intent(in) :: q(:,:,:,:)
  real(rl), intent(inout) :: rcd(6)

  rcd(1) = min(rcd(1), minval(q(:,:,:,1)))
  rcd(2) = max(rcd(2), maxval(q(:,:,:,1)))
  rcd(3) = min(rcd(3), minval(q(:,:,:,2)))
  rcd(4) = max(rcd(4), maxval(q(:,:,:,2)))
  rcd(5) = min(rcd(5), minval(q(:,:,:,1) + 2*q(:,:,:,2)))
  rcd(6) = max(rcd(6), maxval(q(:,:,:,1) + 2*q(:,:,:,2)))
end subroutine toy_rcd

subroutine toy_print(hybrid, nstep, rcd)
  use reduction_mod, only: ParallelMin, ParallelMax

  type(hybrid_t), intent(in) :: hybrid
  integer, intent(in) :: nstep
  real(rl), intent(inout) :: rcd(6)

  if (modulo(nstep,50) == 0) then
     rcd(1) = ParallelMin(rcd(1), hybrid)
     rcd(2) = ParallelMax(rcd(2), hybrid)
     rcd(3) = ParallelMin(rcd(3), hybrid)
     rcd(4) = ParallelMax(rcd(4), hybrid)
     rcd(5) = ParallelMin(rcd(5), hybrid)
     rcd(6) = ParallelMax(rcd(6), hybrid)
     if (hybrid%masterthread) &
          write(*,'(a,i5,es11.3,es11.3,es11.3,es11.3,es11.3,es11.3)') 'toy>', nstep, rcd
  end if
end subroutine toy_print

subroutine dcmip2016_test1_pg_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)
  use element_ops, only: get_field
  use gllfvremap_mod
  use perf_mod, only: t_startf, t_stopf
  use control_mod, only: ftype

  ! to DSS precl
  use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
  use bndry_mod, only: bndry_exchangev

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer, parameter :: iqv = 1
  real(rl), parameter :: one = 1.0_rl

  integer :: i,j,k,ie,qi
  real(rl), dimension(np,np,nlev) :: u,v,w,T,p,dp,rho,rho_dry,z,exner_kess,theta_kess
  real(rl), dimension(np,np,nlev) :: rho_new,p_pk
  real(rl), dimension(nlev)       :: u_c,v_c,p_c,qv_c,qc_c,qr_c,rho_c,z_c, th_c
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, lon, dz_top(np,np),zi(np,np,nlevp),zi_c(nlevp), ps(np,np), &
       wrk(np,np), rd, wrk3(np,np,nlev), wrk4(np,np,nlev,2), wf(np*np,1)

  integer :: nf, ncol
  real(rl), dimension(np,np,nlev) :: dp_fv, p_fv, u_fv, v_fv, T_fv, exner_kess_fv, &
       theta_kess_fv, Rstar, rho_fv, rho_dry_fv, u0, v0, T0, z_fv, ddt_cl, ddt_cl2
  real(rl), dimension(np,np,nlev,qsize) :: Q_fv, Q0_fv
  real(rl), dimension(np,np,nlevp) :: phi_i, zi_fv
  real(rl), dimension(np,np) :: zs_fv, ps_fv, delta_ps
  real(rl) :: precl_fv(np,np,1), rcd(6)
  real(rl), allocatable :: qmin(:,:,:), qmax(:,:,:)

  integer :: pbl_type, prec_type
  integer, parameter :: test = 1

  nf = pg_data%nphys
  ncol = nf*nf

  prec_type = dcmip16_prec_type
  pbl_type  = dcmip16_pbl_type

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  call t_startf('gfr_dyn_to_fv_phys')
  call gfr_dyn_to_fv_phys(hybrid, nt, hvcoord, elem, nets, nete, &
       pg_data%ps, pg_data%zs, pg_data%T, pg_data%uv, pg_data%omega_p, pg_data%q)
  call t_stopf('gfr_dyn_to_fv_phys')

  do ie = nets,nete
     ! for max_w
     call get_field(elem(ie), 'w', w, hvcoord, nt, ntQ)

     T_fv(:nf,:nf,:) = reshape(pg_data%T(:,:,ie), (/nf,nf,nlev/))
     u_fv(:nf,:nf,:) = reshape(pg_data%uv(:,1,:,ie), (/nf,nf,nlev/))
     v_fv(:nf,:nf,:) = reshape(pg_data%uv(:,2,:,ie), (/nf,nf,nlev/))
     Q_fv(:nf,:nf,:,:) = reshape(pg_data%q(:,:,:,ie), (/nf,nf,nlev,qsize/))
     ps_fv(:nf,:nf) = reshape(pg_data%ps(:,ie), (/nf,nf/))
     zs_fv(:nf,:nf) = reshape(pg_data%zs(:,ie), (/nf,nf/))
     zs_fv(:nf,:nf) = zs_fv(:nf,:nf)/g
     do k = 1,nlev
        p_fv(:nf,:nf,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps_fv(:nf,:nf)
        dp_fv(:nf,:nf,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                           (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps_fv(:nf,:nf)
     end do

     ! Rederive the remaining vars so they are self-consistent; use
     ! hydrostatic assumption.
     if (use_moisture) then
        Rstar(:nf,:nf,:) = Rgas + (Rwater_vapor - Rgas)*Q_fv(:nf,:nf,:,iqv)
     else
        Rstar(:nf,:nf,:) = Rgas
     end if
     rho_fv(:nf,:nf,:) = p_fv(:nf,:nf,:)/(Rstar(:nf,:nf,:)*T_fv(:nf,:nf,:))
     phi_i(:nf,:nf,nlevp) = g*zs_fv(:nf,:nf)
     do k = nlev,1,-1
        phi_i(:nf,:nf,k) = phi_i(:nf,:nf,k+1) + &
             (Rstar(:nf,:nf,k)*(dp_fv(:nf,:nf,k)*T_fv(:nf,:nf,k)))/p_fv(:nf,:nf,k)
     end do
     do k=1,nlev
        z_fv(:nf,:nf,k) = (phi_i(:nf,:nf,k)+phi_i(:nf,:nf,k+1))/(2*g)
     end do
     do k=1,nlevp
        zi_fv(:nf,:nf,k) = phi_i(:nf,:nf,k)/g
     end do

     rho_dry_fv(:nf,:nf,:) = (1 - Q_fv(:nf,:nf,:,iqv))*rho_fv(:nf,:nf,:)

     ! Compute form of exner pressure expected by Kessler physics.
     exner_kess_fv(:nf,:nf,:) = (p_fv(:nf,:nf,:)/p0)**(Rgas/Cp)
     theta_kess_fv(:nf,:nf,:) = T_fv(:nf,:nf,:)/exner_kess_fv(:nf,:nf,:)

     u0(:nf,:nf,:) = u_fv(:nf,:nf,:); v0(:nf,:nf,:) = v_fv(:nf,:nf,:); T0(:nf,:nf,:) = T_fv(:nf,:nf,:)
     Q0_fv(:nf,:nf,:,:) = Q_fv(:nf,:nf,:,:);

     ! Convert to dry mixing ratios.
     do i = 1,3
        Q_fv(:nf,:nf,:,i) = (rho_fv(:nf,:nf,:)/rho_dry_fv(:nf,:nf,:))*Q_fv(:nf,:nf,:,i)
     end do

     do j = 1,nf
        do i = 1,nf
           u_c  = u_fv(i,j,nlev:1:-1)
           v_c  = v_fv(i,j,nlev:1:-1)
           qv_c = Q_fv(i,j,nlev:1:-1,1)
           qc_c = Q_fv(i,j,nlev:1:-1,2)
           qr_c = Q_fv(i,j,nlev:1:-1,3)
           p_c  = p_fv(i,j,nlev:1:-1)
           rho_c= rho_dry_fv(i,j,nlev:1:-1)
           z_c  = z_fv(i,j,nlev:1:-1)
           zi_c = zi_fv(i,j,nlevp:1:-1)
           th_c = theta_kess_fv(i,j,nlev:1:-1)

           call gfr_f_get_latlon(ie, i, j, lat, lon)

           ! Get forced versions of u,v,p,qv,qc,qr. rho is constant.
           call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, &
                z_c, zi_c, lat, nlev, precl_fv(i,j,1), pbl_type, prec_type)

           u_fv(i,j,:)   = u_c(nlev:1:-1)
           v_fv(i,j,:)   = v_c(nlev:1:-1)
           Q_fv(i,j,:,1) = qv_c(nlev:1:-1)
           Q_fv(i,j,:,2) = qc_c(nlev:1:-1)
           Q_fv(i,j,:,3) = qr_c(nlev:1:-1)
           theta_kess_fv(i,j,:) = th_c(nlev:1:-1)

           do k=1,nlev
              call tendency_terminator(lat*rad2dg, lon*rad2dg, Q_fv(i,j,k,4), Q_fv(i,j,k,5), &
                   dt, ddt_cl(i,j,k), ddt_cl2(i,j,k))
           enddo
        enddo
     enddo

     do i = 1,3
        Q_fv(:nf,:nf,:,i) = (rho_dry_fv(:nf,:nf,:)/rho_fv(:nf,:nf,:))*Q_fv(:nf,:nf,:,i)
     end do
     Q_fv(:nf,:nf,:,4) = Q_fv(:nf,:nf,:,4) + dt*ddt_cl(:nf,:nf,:)
     Q_fv(:nf,:nf,:,5) = Q_fv(:nf,:nf,:,5) + dt*ddt_cl2(:nf,:nf,:)

     ! Convert from theta to T w.r.t. new model state.
     ! Assume hydrostatic pressure pi changed by qv forcing.
     ! Assume NH pressure perturbation is unchanged.
     delta_ps(:nf,:nf) = sum(dp_fv(:nf,:nf,:)*(Q_fv(:nf,:nf,:,iqv) - Q0_fv(:nf,:nf,:,iqv)), 3)
     do k = 1,nlev
        p_fv(:nf,:nf,k) = p_fv(:nf,:nf,k) + hvcoord%hybm(k)*delta_ps(:nf,:nf)
     enddo
     exner_kess_fv(:nf,:nf,:) = (p_fv(:nf,:nf,:)/p0)**(Rgas/Cp)
     T_fv(:nf,:nf,:) = exner_kess_fv(:nf,:nf,:)*theta_kess_fv(:nf,:nf,:)

     ! These gfr calls are special to this routine, to handle
     ! DCMIP-specific precl.
     wf(:ncol,1) = reshape(precl_fv(:nf,:nf,1), (/ncol/))
     call gfr_f2g_scalar(ie, elem(ie)%metdet, wf(:,:1), wrk3(:,:,:1))
     call gfr_g_make_nonnegative(elem(ie)%metdet, wrk3(:,:,:1))
     precl(:,:,ie) = wrk3(:,:,1)

     ! T, uv tendencies
     pg_data%T(:,:,ie) = reshape((T_fv(:nf,:nf,:) - T0(:nf,:nf,:))/dt, (/ncol,nlev/))
     pg_data%uv(:,1,:,ie) = reshape((u_fv(:nf,:nf,:) - u0(:nf,:nf,:))/dt, (/ncol,nlev/))
     pg_data%uv(:,2,:,ie) = reshape((v_fv(:nf,:nf,:) - v0(:nf,:nf,:))/dt, (/ncol,nlev/))
     ! q state
     do i = 1,qsize
        pg_data%q(:,:,i,ie) = reshape(Q_fv(:nf,:nf,:,i), (/ncol,nlev/))
     end do
     
     ! Measure max w and max prect. w is not used in the physics, so
     ! just look at the GLL values.
     max_w = max(max_w, maxval(w))
     ! ps isn't updated by the physics, so just look at the GLL values.
     min_ps = min(min_ps, minval(elem(ie)%state%ps_v(:,:,nt)))
  enddo

  call t_startf('gfr_fv_phys_to_dyn')
  call gfr_fv_phys_to_dyn(hybrid, nt, hvcoord, elem, nets, nete, &
       pg_data%T, pg_data%uv, pg_data%q)
  call t_stopf('gfr_fv_phys_to_dyn')
  ! dp_coupling doesn't do the DSS; stepon does. Thus, this DCMIP test
  ! also needs to do its own DSS.
  call gfr_f2g_dss(hybrid, elem, nets, nete)
  call gfr_pg1_reconstruct(hybrid, nt, hvcoord, elem, nets, nete)

  call toy_init(rcd)
  do ie = nets,nete
     do i = 1,2
        wrk4(:,:,:,i) = elem(ie)%state%Q(:,:,:,i+3)
     end do
     call toy_rcd(wrk4, rcd)
  end do
  call toy_print(hybrid, tl%nstep, rcd)

  if (ftype == 0) then
     ! Convert FQ from state to Qdp tendency.
     do ie = nets,nete
        do k = 1,nlev
           dp(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(:,:,nt)
        end do
        do i = 1,qsize
           elem(ie)%derived%FQ(:,:,:,i) = &
                dp*(elem(ie)%derived%FQ(:,:,:,i) - elem(ie)%state%Q(:,:,:,i))/dt
        end do
     end do
  end if

  ! DSS precl
  do ie = nets,nete
     precl(:,:,ie) = precl(:,:,ie)*elem(ie)%spheremp
     call edgeVpack_nlyr(edge_g, elem(ie)%desc, precl(:,:,ie), 1, 0, 1)
  end do
  call bndry_exchangeV(hybrid, edge_g)
  do ie = nets,nete
     call edgeVunpack_nlyr(edge_g, elem(ie)%desc, precl(:,:,ie), 1, 0, 1)
     precl(:,:,ie) = precl(:,:,ie)*elem(ie)%rspheremp
     max_precl = max( max_precl, maxval(precl(:,:,ie)) )
  end do

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)
end subroutine dcmip2016_test1_pg_forcing

!_______________________________________________________________________
subroutine dcmip2016_test2_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl, test)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure
  integer,            intent(in)            :: test                     ! dcmip2016 test number

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl), dimension(np,np,nlev) :: u,v,w,T,exner_kess,theta_kess,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: u0,v0,T0,qv0,qc0,qr0
  real(rl), dimension(np,np,nlev) :: rho_dry,rho_new,Rstar,p_pk
  real(rl), dimension(np,np)      :: delta_ps(np,np)
  real(rl), dimension(nlev)       :: u_c,v_c,p_c,qv_c,qc_c,qr_c,rho_c,z_c, th_c
  real(rl) :: max_w, max_precl, min_ps
  real(rl) :: lat, dz_top(np,np), zi(np,np,nlevp),zi_c(nlevp), ps(np,np)

  integer :: pbl_type, prec_type

  prec_type = dcmip16_prec_type
  pbl_type  = dcmip16_pbl_type

  if(test==3) prec_type=0 ! kessler only for test 3

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    precl(:,:,ie) = -1.0d0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)
    theta_kess = T/exner_kess

    ! get wet mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere

    rho_dry = (1-qv)*rho  ! convert to dry density using wet mixing ratio

    ! convert to dry mixing ratios
    qv  = qv*rho/rho_dry
    qc  = qc*rho/rho_dry
    qr  = qr*rho/rho_dry


    ! save un-forced prognostics (DRY)
    u0=u; v0=v; T0=T; qv0=qv; qc0=qc; qr0=qr


    ! apply forcing to columns
    do j=1,np; do i=1,np

      ! invert column
      u_c  = u  (i,j,nlev:1:-1)
      v_c  = v  (i,j,nlev:1:-1)
      qv_c = qv (i,j,nlev:1:-1)
      qc_c = qc (i,j,nlev:1:-1)
      qr_c = qr (i,j,nlev:1:-1)
      p_c  = p  (i,j,nlev:1:-1)
      rho_c= rho_dry(i,j,nlev:1:-1)
      z_c  = z  (i,j,nlev:1:-1)
      zi_c = zi (i,j,nlevp:1:-1)
      th_c = theta_kess(i,j,nlev:1:-1)

      lat=0.0
      ! get forced versions of u,v,p,qv,qc,qr. rho is constant
      call DCMIP2016_PHYSICS(test, u_c, v_c, p_c, th_c, qv_c, qc_c, qr_c, rho_c, dt, z_c, zi_c, lat, nlev, &
                             precl(i,j,ie), pbl_type, prec_type)

      ! revert column
      u(i,j,:)  = u_c(nlev:1:-1)
      v(i,j,:)  = v_c(nlev:1:-1)
      qv(i,j,:) = qv_c(nlev:1:-1)
      qc(i,j,:) = qc_c(nlev:1:-1)
      qr(i,j,:) = qr_c(nlev:1:-1)
      theta_kess(i,j,:) = th_c(nlev:1:-1)

    enddo; enddo;
    ! convert from theta to T w.r.t. new model state
    ! assume hydrostatic pressure pi changed by qv forcing
    ! assume NH pressure perturbation unchanged
    delta_ps = sum( (rho_dry/rho)*dp*(qv-qv0) , 3 )
    do k=1,nlev
       p(:,:,k) = p(:,:,k) + hvcoord%hybm(k)*delta_ps(:,:)
    enddo
    exner_kess = (p/p0)**(Rgas/Cp)
    T = exner_kess*theta_kess


    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = (u - u0)/dt
    elem(ie)%derived%FM(:,:,2,:) = (v - v0)/dt
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt


    ! set tracer-mass forcing. conserve tracer mass
    ! rho_dry*(qv-qv0)*dz = FQ deta, dz/deta = -dp/(g*rho)
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt


    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine dcmip2016_test2_forcing

!_______________________________________________________________________
subroutine dcmip2016_test3_forcing(elem,hybrid,hvcoord,nets,nete,nt,ntQ,dt,tl)

  type(element_t),    intent(inout), target :: elem(:)                  ! element array
  type(hybrid_t),     intent(in)            :: hybrid                   ! hybrid parallel structure
  type(hvcoord_t),    intent(in)            :: hvcoord                  ! hybrid vertical coordinates
  integer,            intent(in)            :: nets,nete                ! start, end element index
  integer,            intent(in)            :: nt, ntQ                  ! time level index
  real(rl),           intent(in)            :: dt                       ! time-step size
  type(TimeLevel_t),  intent(in)            :: tl                       ! time level structure

  integer :: i,j,k,ie                                                     ! loop indices
  real(rl):: lat
  real(rl), dimension(np,np,nlev) :: u,v,w,T,theta_kess,exner_kess,p,dp,rho,z,qv,qc,qr
  real(rl), dimension(np,np,nlev) :: T0,qv0,qc0,qr0
  real(rl), dimension(np,np,nlev) :: rho_dry,rho_new,Rstar,p_pk
  real(rl), dimension(np,np,nlev) :: theta_inv,qv_inv,qc_inv,qr_inv,rho_inv,exner_inv,z_inv ! inverted columns
  real(rl), dimension(np,np,nlevp):: zi
  real(rl), dimension(np,np)      :: ps,delta_ps(np,np)
  real(rl) :: max_w, max_precl, min_ps

  max_w     = -huge(rl)
  max_precl = -huge(rl)
  min_ps    = +huge(rl)

  do ie = nets,nete

    precl(:,:,ie) = 0.0

    ! get current element state
    call get_state(u,v,w,T,p,dp,ps,rho,z,zi,g,elem(ie),hvcoord,nt,ntQ)

    ! get mixing ratios
    qv  = elem(ie)%state%Qdp(:,:,:,1,ntQ)/dp
    qc  = elem(ie)%state%Qdp(:,:,:,2,ntQ)/dp
    qr  = elem(ie)%state%Qdp(:,:,:,3,ntQ)/dp

    ! ensure positivity
    where(qv<0); qv=0; endwhere
    where(qc<0); qc=0; endwhere
    where(qr<0); qr=0; endwhere

    rho_dry = (1-qv)*rho  ! convert to dry density using wet mixing ratio

    ! convert to dry mixing ratios
    qv  = qv*rho/rho_dry
    qc  = qc*rho/rho_dry
    qr  = qr*rho/rho_dry

    ! compute form of exner pressure expected by Kessler physics
    exner_kess = (p/p0)**(Rgas/Cp)
    theta_kess = T/exner_kess

    ! save un-forced prognostics
    T0=T; qv0=qv; qc0=qc; qr0=qr


    ! invert columns (increasing z)
    theta_inv= theta_kess(:,:,nlev:1:-1)
    qv_inv   = qv   (:,:,nlev:1:-1)
    qc_inv   = qc   (:,:,nlev:1:-1)
    qr_inv   = qr   (:,:,nlev:1:-1)
    rho_inv  = rho_dry  (:,:,nlev:1:-1)
    exner_inv= exner_kess(:,:,nlev:1:-1)
    z_inv    = z    (:,:,nlev:1:-1)

    ! apply forcing to columns
    do j=1,np; do i=1,np

      ! apply physics to (u,v,p, qv,qc,qr). (rho held constant)
      CALL KESSLER(       &
        theta_inv(i,j,:), &
        qv_inv(i,j,:),    &
        qc_inv(i,j,:),    &
        qr_inv(i,j,:),    &
        rho_inv(i,j,:),   &
        exner_inv(i,j,:), &
        dt,               &
        z_inv(i,j,:),     &
        nlev,             &
        precl(i,j,ie))

    enddo; enddo;

    ! revert columns (increasing eta)
    theta_kess = theta_inv(:,:,nlev:1:-1)
    qv    = qv_inv   (:,:,nlev:1:-1)
    qc    = qc_inv   (:,:,nlev:1:-1)
    qr    = qr_inv   (:,:,nlev:1:-1)


    ! convert from theta to T w.r.t. new model state
    ! assume hydrostatic pressure pi changed by qv forcing
    ! assume NH pressure perturbation unchanged
    delta_ps = sum( (rho_dry/rho)*dp*(qv-qv0) , 3 )
    do k=1,nlev
       p(:,:,k) = p(:,:,k) + hvcoord%hybm(k)*delta_ps(:,:)
    enddo
    exner_kess = (p/p0)**(Rgas/Cp)
    T = exner_kess*theta_kess

    ! set dynamics forcing
    elem(ie)%derived%FM(:,:,1,:) = 0
    elem(ie)%derived%FM(:,:,2,:) = 0
    elem(ie)%derived%FT(:,:,:)   = (T-T0)/dt

    ! set tracer-mass forcing. conserve tracer mass.  
    ! rho_dry*(qv-qv0)*dz = FQ deta, dz/deta = -dp/(g*rho)
    elem(ie)%derived%FQ(:,:,:,1) = (rho_dry/rho)*dp*(qv-qv0)/dt
    elem(ie)%derived%FQ(:,:,:,2) = (rho_dry/rho)*dp*(qc-qc0)/dt
    elem(ie)%derived%FQ(:,:,:,3) = (rho_dry/rho)*dp*(qr-qr0)/dt
    ! after forcings above are applied:
    ! Qdp_new = (rho_dry/rho)*dp*qv   


    ! perform measurements of max w, and max prect
    max_w     = max( max_w    , maxval(w    ) )
    max_precl = max( max_precl, maxval(precl(:,:,ie)) )
    min_ps    = min( min_ps,    minval(ps) )

  enddo

  call dcmip2016_append_measurements(max_w,max_precl,min_ps,tl,hybrid)

end subroutine


subroutine qsat_kessler(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_const1 / p * exp( bubble_const2 * (T - bubble_const3) / ( T - bubble_const4 ) )
end subroutine qsat_kessler

subroutine qsat_rj(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_epsilo * bubble_e0 / p * &
         exp(-(latvap/Rwater_vapor) * ((1.0/T)-(1.0/bubble_t0_const)))
end subroutine qsat_rj


end module dcmip16_wrapper
#endif
