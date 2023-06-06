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

!kessler constants, put into physconst mod
real(rl), parameter:: k1=0.001,k2=2.2,k3=0.875,a1=0.00001 !a1 should be 0.001
real(rl), parameter:: r1=0.525,r2=5.4e5,r3=2.55e6
real(rl), parameter:: r4=1.6,r5=124.9,r6=0.2046

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

!!!!!!!!!!!!!!!!!! there should be only one phi, "wet", this sub is not valid
!does not need moisture
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


subroutine sedimentation_nh(qvdry,qcdry,qrdry,tempe,dpdry,ptop,zi,pnh,massleft,energyleft,dt)

  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,qrdry,tempe,pnh
  real(rl), dimension(nlevp),intent(inout) :: zi
  real(rl), dimension(nlev), intent(in)    :: dpdry
  real(rl),                  intent(in)    :: dt,ptop
  real(rl),                  intent(inout) :: massleft, energyleft

  integer  :: k, kk, d_ind
  real(rl) :: change, positt, velqr, T_new, ctermold, ctermnew, zbottom
  real(rl) :: deltamass, deltaenergy
  real(rl), dimension(nlev)  :: zm, dz, rhodry, dpi, pidry, ppi, dphi, rstar

  massleft = 0.0; energyleft = 0.0

  zbottom=zi(nlevp)
  zm = (zi(1:nlev)+zi(2:nlevp))/2.0
  dpi = dpdry*(1.0 + qvdry + qcdry + qrdry)
  !derived pressure values
  !not used
  !call construct_hydro_pressure(dpi,  ptop,ppi  )
  call construct_hydro_pressure(dpdry,ptop,pidry)
  call get_geo_from_drydp(tempe, dpdry, pidry, zbottom, zi, zm, dz)
  rhodry = dpdry / dz / gravit

  do k=1,nlev

    if(qrdry(k)>0.0) then
      ! WHY NOT WET rho? 
      ! Liquid water terminal velocity, using dry densities instead of wet
      ! in meter/sec
      velqr  = 36.34d0*(qrdry(k)*rhodry(k))**0.1364*sqrt(rhodry(nlev)/rhodry(k))

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

if (k < nlev) then
!stick to the level below
!  positt = zm(k+1)
!  d_ind = k+1

!rain it down immediately
positt = zi(nlevp)-100.0
   
else
  positt = zi(nlevp)-100.0 !to ensure cond below
endif

      deltamass = qrdry(k) * dpdry(k)
      !hy
      !deltaenergy =  qrdry(k) * dpdry(k) * ( cl*tempe(k) + latice + zbottom*gravit )
      !nh
      deltaenergy = qrdry(k) * dpdry(k) * ( cl*tempe(k) + latice + zm(k)*gravit )

      if(positt < zi(nlevp)) then
        !hit the bottom
        d_ind = -1
        massleft   = massleft   + deltamass
        energyleft = energyleft + deltaenergy
        qrdry(k) = 0.0
      else
        change = deltamass/dpdry(d_ind)
#if 0
        !hy update in terms of cp
        !compute average temperature with new mass
        !arrives to dest cell, recompute qrdry in terms of dpdry of the dest cell
        ctermold = cpdry + cpv * qvdry(d_ind) + cl * (qcdry(d_ind) + qrdry(d_ind))
        ctermnew = cptermold + cl*change
        T_new = ( tempe(d_ind) * dpdry(d_ind) * cptermold + tempe(k) * dpdry(k) * qrdry(k) * cl ) / &
                 cptermnew / dpdry(d_ind)
#endif

#if 1
        !nh update in terms of cv
        ! Tnew * cvstar * (dp+newmass) + phi(dp+newmass) = Told * cvstar * dp + phi*dp +newenergy
        ! this is CV update, this won't preserve NH pressure
        ctermold = cvdry + cvv*qvdry(d_ind) + cl*(qcdry(d_ind) + qrdry(d_ind))
        ctermnew = ctermold + cl*change
    
        T_new = ( tempe(d_ind)*ctermold*dpdry(d_ind) - gravit*zm(d_ind)*deltamass + deltaenergy ) &
                / ctermnew / dpdry(d_ind)

        tempe(d_ind) = T_new
#endif
        qrdry(k) = 0.0
        qrdry(d_ind) = qrdry(d_ind) + change
      endif
    endif !there was rain in cell
  end do ! k loop for sedim

  !none of the above used NH pressure, and all above assumed const volume
  ! so nh pressure needs to be recomputed
  !we could have recomputed pnh for each iteration for k cell and its destination, 
  !but this is simpler
  dphi = gravit*( zi(1:nlev) - zi(2:nlevp) )

!!!!!!!!!!!!!!do i need to change this Rstar too???
  !homme rstar
  !rstar = rdry*dpdry*(1.0 + qcdry + qrdry) + rvapor*dpdry*qvdry
  !correct Rstar
  !rstar = rdry*dpdry + rvapor*dpdry*qvdry
  !pnh = rstar*tempe/dphi
  
end subroutine sedimentation_nh


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




subroutine rain_evaporation(qvdry,qcdry,qrdry,tempe,dpdry,ptop,zbottom,p_c)

  real(rl), dimension(nlev), intent(in)    :: qcdry
  real(rl), dimension(nlev), intent(inout) :: qvdry,qrdry,tempe,p_c
  real(rl), dimension(nlev), intent(in)    :: dpdry
  real(rl),                  intent(in)    :: zbottom,ptop

  integer  :: k
  real(rl) :: dq, qsat, qsatdry, cval
  real(rl), dimension(nlev)  :: zm,dz,rhodry,dpi,ppi,pnh,pidry
  real(rl), dimension(nlevp) :: zi

  dpi = dpdry*(1.0 + qvdry + qcdry + qrdry)

  !derived pressure values
  call construct_hydro_pressure(dpi,  ptop,ppi  )
  call construct_hydro_pressure(dpdry,ptop,pidry)
  pnh = p_c
  
  !we can use wet rho, but we will ignore thah here
                          !input                        !output
  call get_geo_from_drydp(tempe, dpdry, pidry, zbottom, zi, zm, dz)
  rhodry = dpdry / dz / gravit

  do k=1, nlev

    !r1=0.525,r2=5.4e5,r3=2.55e6 
    !r4=1.6,r5=124.9,r6=0.2046

    !call qsat_kessler2(pnh(k), tempe(k), qsat)
    call qsat_kessler2(pidry(k), tempe(k), qsat)
    qsatdry = qsat*dpi(k)/dpdry(k)

    cval = r4 + r5*( (rhodry(k)*qrdry(k))**r6 )
    dq = (1 - qvdry(k)/qsatdry) * cval**r1
    dq = dq / ( r2 + r3/pnh(k)/qsat ) / rhodry(k)
    !change the sign to conform to the phase_change... routine
    dq = -dq
    dq = max( dq, -qcdry(k) )

    ! new values: qv = qv - dq, qr = qr + dq
    !notice liquid_in_use var is now qrdry
    call phase_change_gas_liquid_level( &
         qvdry(k),qrdry(k),qcdry(k),dq,tempe(k),dpdry(k),dpi(k),ppi(k),pnh(k) )

!if (dq > 0) then
! print *, 'yay, condensation',' qc', qcdry(k)
!endif

  enddo

end subroutine rain_evaporation



subroutine condensation_and_back_again(qvdry,qcdry,qrdry,tempe,dpdry,ptop,zi,pnh,ie)

  real(rl), dimension(nlev), intent(in)    :: qrdry, dpdry
  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,tempe,pnh
  real(rl), dimension(nlevp),intent(inout) :: zi
  real(rl),                  intent(in)    :: ptop

  integer, intent(in), optional :: ie

  integer  :: k, ii, jj
  real(rl) :: dq, qsat, qsatdry  !, dphi, rstar, olddphi
  real(rl), dimension(nlev)  :: dpi,ppi,pidry, dphi, rstar

  dpi = dpdry*(1.0 + qvdry + qcdry + qrdry)

  !derived pressure values
  call construct_hydro_pressure(dpi,  ptop,ppi  )
  call construct_hydro_pressure(dpdry,ptop,pidry)

  do k=1, nlev
    call qsat_kessler2(pnh(k), tempe(k), qsat)
    !call qsat_kessler2(pidry(k), tempe(k), qsat)
    qsatdry = qsat*dpi(k)/dpdry(k)

    !assume condensation
    dq = qvdry(k) - qsatdry
    !except if
    if ( dq < 0.0 ) then
      !evaporation
      dq = max( dq, -qcdry(k) )
    endif
    ! new values: qv = qv - dq, qc = qc + dq

!if(abs(dq)>0.0)then
!print *,'qv,qc,qsat,dq',qvdry(k),qcdry(k),qsatdry,dq
!print *, 'T old', tempe(k)
!endif
!rstar = rdry*dpdry(k) + rvapor*dpdry(k)*qvdry(k)
!print *, 'rstar', rstar
    !pnh = rstar*tempe(k)/dphi*dpi(k)
!print *, 'dphi as from EOS, and computed' ,rstar*tempe(k)/pnh(k),zi(k) - zi(k+1)
!print *, 'dphi as from EOS HY, and NH' ,rstar*tempe(k)/pnh(k)/gravit,rstar*tempe(k)/ppi(k)/gravit


    !this update needs only level values
    !there will be updates that need whole column
    call phase_change_gas_liquid_level_hydroupdate( &
         qvdry(k),qcdry(k),qrdry(k),dq,tempe(k),dpdry(k) )

!if(abs(dq)>0.0)then
!print *, 'after update -------------'
!print *,'qv,qc,qsat,dq',qvdry(k),qcdry(k),qsatdry,dq
!print *, 'T new', tempe(k)
!endif

    !recompute zi for all levels <= k
!    rstar = rdry*dpdry(k) + rvapor*dpdry(k)*qvdry(k)
!    !pnh = rstar*tempe(k)/dphi*dpi(k)
!    dphi = rstar*tempe(k)/pnh(k)
!print *,'AFTER qv,qc,qsat,dq',qvdry(k),qcdry(k),qsatdry,dq
!print *, 'T new', tempe(k)
!print *, 'rstar', rstar
!print *, 'olddphi and new dphi', zi(k) - zi(k+1), dphi
!    olddphi = zi(k) - zi(k+1)
!    zi(1:k) = zi(1:k) + (dphi - olddphi)/gravit

!if (dq > 0) then
! print *, 'yay, condensation',' qc', qcdry(k)
!endif

  enddo


!!!!!!!!!!!!!! this is where init, that is HY phi or NH phi might matter
!updating 1 layer at a time (shift all z above k cell)
!does not seem to work 
!let's do a final zi update here instead
!print *, 'old zi ', zi



!!!!!!!!!!!!!! correct R*
!rstar = rdry*dpdry + rvapor*dpdry*qvdry

!!!!!!!!!!!!!! incorrect, matching homme
rstar = rdry*dpdry*(1.0+qvdry+qcdry+qrdry) + (rvapor-rdry)*dpdry*qvdry
dphi = rstar*tempe/pnh

#if 0

if(present(ie))then
print *, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
ii=110; jj=110;
print *, 'in here old g*zi is', gravit*zi(ii:jj)
endif

do k=nlev,1,-1
zi(k) = zi(k+1) + dphi(k)/gravit
enddo

if(present(ie))then
print *, 'rstar HOMME here', rstar(ii:jj)
!print *, 'rstar Homme pieces', rdry,dpdry(ii:jj)*(1.0+qvdry(ii:jj)+qcdry(ii:jj)+qrdry(ii:jj)),(rvapor-rdry),dpdry(ii:jj)*qvdry(ii:jj)
!print *, 'rstar WL here', rdry*dpdry(ii:jj) + rvapor*dpdry(ii:jj)*qvdry(ii:jj)
print *, 'new vapor mass', dpdry(ii:jj)*qvdry(ii:jj)
print *, 'new qc mass', dpdry(ii:jj)*qcdry(ii:jj)
print *, 'tempe here', tempe(ii:jj)
print *, 'pnh here', pnh(ii:jj)
print *, 'dphi here is',dphi(ii:jj)
print *, 'in here new g*zi is', gravit*zi(ii:jj)
endif
#endif

!print *, 'new zi ', zi(120:128)


end subroutine condensation_and_back_again


! convension: vapor = vapor - dq, liquid_in_use = liquid_in_use + dq
subroutine phase_change_gas_liquid_level_hydroupdate(&
           vapor,liquid_in_use,liquid_not_in_use,dq,tempe,dpdry)

  real(rl),  intent(in)    :: liquid_not_in_use, dpdry, dq
  real(rl),  intent(inout) :: vapor,liquid_in_use,tempe

  real(rl) :: T_new, cptermold, cptermnew, Lold, Lnew

  cptermold = cpdry     + cpv *  vapor     + cl * (liquid_in_use + liquid_not_in_use)
  cptermnew = cptermold + cpv * ( -dq )    + cl * dq

  Lold =        (latvap+latice) * vapor    + latice * (liquid_in_use + liquid_not_in_use)
  Lnew = Lold + (latvap+latice) * ( -dq )  + latice * dq

  !dpdry is cancelled from both sides
  T_new = ( tempe * cptermold + Lold - Lnew ) / cptermnew

  tempe = T_new

  vapor = vapor - dq
  liquid_in_use = liquid_in_use + dq

!if (dq > 0) then
! print *, 'yay, condensation',' qc', qcdry(k)
!endif

end subroutine phase_change_gas_liquid_level_hydroupdate


!!!! use as template
! convension: vapor = vapor - dq, liquid_in_use = liquid_in_use + dq
subroutine phase_change_gas_liquid_level(vapor,liquid_in_use,liquid_not_in_use,dq,tempe,dpdry,dp,pi,pnh)

  real(rl),  intent(in)    :: liquid_not_in_use
  real(rl),  intent(inout) :: vapor,liquid_in_use,tempe
  real(rl),  intent(in)    :: dpdry, dp, pi, pnh
  real(rl),  intent(in)    :: dq

  integer  :: k
  real(rl) :: T_new, cptermold, cptermnew, Lold, Lnew

  cptermold = cpdry     + cpv *  vapor     + cl * (liquid_in_use + liquid_not_in_use)
  cptermnew = cptermold + cpv * ( -dq )    + cl * dq

  Lold =        (latvap+latice) * vapor    + latice * (liquid_in_use + liquid_not_in_use)
  Lnew = Lold + (latvap+latice) * ( -dq )  + latice * dq

  !dpdry is cancelled from both sides
  T_new = ( tempe * cptermold + Lold - Lnew ) / cptermnew

  tempe = T_new

  vapor = vapor - dq
  liquid_in_use = liquid_in_use + dq

!if (dq > 0) then
! print *, 'yay, condensation',' qc', qcdry(k)
!endif

end subroutine phase_change_gas_liquid_level




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
                                        en1, en2, en3, wasiactive)

  real(rl), dimension(nlev), intent(inout) :: p_c
  real(rl), dimension(nlev), intent(inout) :: dp_c
  real(rl), dimension(nlev), intent(inout) :: qv_c,ql_c,T_c
  real(rl), dimension(nlevp),intent(inout) :: zi_c
  real(rl),                  intent(inout) :: massout, energyout, en1, en2, en3
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
  en1 = 0; en2 = 0; en3 = 0; 

  call construct_hydro_pressure(dp_c,ptop,pi)

  !compute en1, energy before condensation
  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),p_c,en1)

#if 0
!experiment
  qv_c = 0.01
  ql_c = qv_c

!recompute pressure
!  do k=1,nlev
!     rstar_new = rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k)
!     olddphi = gravit*(zi_c(k) - zi_c(k+1))
!     p_c(k) = rstar_new * dp_c(k) * T_c(k) / olddphi
!  enddo

!what if we recompute dz instead?

  do k=nlev, 1, -1
    rstar_new = rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k)
    olddphi = rstar_new * dp_c(k) * T_c(k) / p_c(k)
    zi_c(k) = zi_c(k+1) + olddphi/gravit
  enddo

  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c(nlevp),p_c,en1)
print *, 'with liquid P formula', en1
  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c       ,p_c,en1)
print *, 'with liquid V formula', en1
#endif

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
         hold  = T_loc*( cpdry*1.0 + cpv*qvdry_loc + cl*qldry_loc ) + L_old
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

  !sanity check to use the opposite formulation
  !call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),p_c,en2)
  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,en2)

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

  enddo
  endif

  !interchange V and P formulas of energy for debugging
  call energycp_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c(nlevp),p_c,en3)
  !call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c, dp_c*ql_c, zero,T_c,ptop,zi_c,p_c,en3)
 
 
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

  real(rl) :: qsat, dp_loc, qv_loc, ql_loc, dpdry_loc, qvdry_loc, qsatdry, dq_loc, vapor_mass_change
  real(rl) :: T_loc, p_loc, pi_loc, L_old, L_new, Iold, T_new
  real(rl) :: cvstarTerm_new, rstardp, oldQ1mass, oldQ2mass
  real(rl), dimension(nlev) :: pi, zero, rain
  integer  :: k

  zero = 0.0
  massout = 0.0
  energyout = 0.0
  rain = 0.0
  en1 = 0; en2 = 0; en3 = 0; encl = 0;

  call construct_hydro_pressure(dp_c,ptop,pi)

  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c), dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en1)

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
       qsatdry = qsat * dp_loc / dpdry_loc     ! qv_new
       dq_loc = qvdry_loc - qsatdry            ! > 0 , is qliq_dry_new
       vapor_mass_change = dpdry_loc*dq_loc
       T_loc = T_c(k)
       p_loc = p_c(k)
       pi_loc = pi(k)

       !new Q will be qsatdry
       L_old = (latvap + latice) * qvdry_loc
       L_new = (latvap + latice) * qsatdry   + latice * dq_loc


!fix this
       Iold = T_loc*( cvdry*1.0 + cvv*qvdry_loc ) + L_old

!and this
       cvstarTerm_new = cvdry*1.0 + cvv*qsatdry + cl*dq_loc
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

  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c),dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en2)

  if(.not.bubble_rj_nosedim)then

  !sedimentation stage
  do k=1, nlev
      
     !sedimentation stage
     oldQ1mass = dp_c(k)*qv_c(k)
     oldQ2mass = dp_c(k)*ql_c(k)
     dp_c(k) = dp_c(k) - rain(k)
     !set a limiter here?
     ql_c(k) = (oldQ2mass - rain(k))/dp_c(k)
     qv_c(k) = oldQ1mass/dp_c(k)

     rstardp = dp_c(k)* (rdry * (1.0 - qv_c(k) - ql_c(k)) + rvapor * qv_c(k))
     p_c(k)  = - rstardp * T_c(k) / (zi_c(k+1) - zi_c(k) ) / gravit

     massout = massout + rain(k)
     energyout = energyout + rain(k)*(T_c(k)*cl + latice + (zi_c(k) + zi_c(k+1))*gravit/2.0)
     encl      = encl + rain(k)*(T_c(k)*cl + latice)
  enddo
  endif

  call energycV_nh_via_mass(dp_c*(1-qv_c-ql_c),dp_c*qv_c,dp_c*ql_c,zero,T_c,ptop,zi_c,p_c,en3)

end subroutine rj_new_volume





!in NH case, uses NH pressure
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



!operate with state
subroutine kessler_new(qv_c,qc_c,qr_c,T_c,dp_c,dpdry_c,p_c,ptop,zi_c,&
                             massout,energyout,dt,wasiactive, ie)
 
  integer, intent(in), optional :: ie

  real(rl), dimension(nlev), intent(in)    :: dpdry_c
  real(rl), dimension(nlev), intent(inout) :: p_c, dp_c
  real(rl), dimension(nlevp),intent(inout) :: zi_c
  real(rl), dimension(nlev), intent(inout) :: qv_c,qc_c,qr_c,T_c
  real(rl),                  intent(inout) :: massout, energyout
  real(rl),                  intent(in)    :: ptop, dt
  logical,                   intent(inout) :: wasiactive

  real(rl), dimension(nlev) :: ppi, ppidry, ploc_c, dploc_c
  real(rl), dimension(nlev) :: qvdry_c, qcdry_c, qrdry_c

  real(rl) :: zbottom, energy_start_timestep, energyhy_before, energyhy_after, &
              energynh_before, energynh_after, loc_mass_p, loc_energy_p
  integer  :: k

  real(rl), parameter:: tol_energy = 1e-12, tol_mass = 1e-12

  massout = 0.0; energyout = 0.0
  zbottom = zi_c(nlevp)

  !derived pressure values
  call construct_hydro_pressure(dp_c,ptop,ppi)

  !derived dry pressure values
  !dry density, the only quantity that won't change
  call construct_hydro_pressure(dpdry_c,ptop,ppidry)

  !convert all tracers to dry ratios
  call convert_to_dry(qv_c, qc_c, qr_c, dp_c, dpdry_c, qvdry_c, qcdry_c, qrdry_c)

  !if there is any water int he column
  if( any(qv_c>0).or.any(qc_c>0).or.any(qr_c>0) ) then

     !not exactly, as below conditions might not get triggered
     wasiactive = .true.
     ! Cr, Ar stages ----------------------------------------------------------
     ! these can convert masses of qc, qr into each other only
     ! so, no updates of other state variables
!     call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_before)
!     call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_before)
     energy_start_timestep = energynh_before
     call accrecion_and_accumulation(qcdry_c, qrdry_c, dt)
!     call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_after)
!     call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_after)
     !print *, 'ACC-C stage HY: enbef - enafter', (energyhy_before - energyhy_after)/energyhy_after
     !print *, 'ACC-C stage NH: enbef - enafter', (energynh_before - energynh_after)/energynh_after

     ! sedimentation ----------------------------------------------------------
     ! right now nh term is not used, so, no need to recompute wet hydro pressure and total nh pressure
     !so far, it is only part that has fluxes out
!     call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_before)
!     call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_before)
     call sedimentation_nh (qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zi_c,p_c, loc_mass_p,loc_energy_p,dt)
!     call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_after)
!     call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_after)
!     if(loc_energy_p > 300.0)then
!     print *, 'Sedime HY:enb-ena(up to flux)', (energyhy_before - energyhy_after - loc_energy_p)/energyhy_before, loc_energy_p
!     print *, 'Sedime NH:enb-ena(up to flux)', (energynh_before - energynh_after - loc_energy_p)/energynh_before, loc_energy_p
!     print *, 'mass out', loc_mass_p
     !stop
!     endif
     massout = massout + loc_mass_p; energyout = energyout + loc_energy_p;
     !!!!!!!!!!! update p here!

     !do this later
     ! evaporation of rain ----------------------------------------------------
     !call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energy_before)
     !call rain_evaporation(qvdry_c,qcdry_c,qrdry_c, T_c, dpdry_c,ptop,zbottom,p_c)
     !call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energy_after)
     !print *, 'Rain evap: enbefore - enafter(up to flux)', (energy_before - energy_after)/energy_before

     ! condensation <-> evaporation -------------------------------------------
     !call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_before)
     !call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_before)

!if(present(ie))then
!     call condensation_and_back_again(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zi_c,p_c,ie)
!else
     call condensation_and_back_again(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zi_c,p_c)
!endif

     !call energy_nh_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,p_c,energynh_after)
     !call energy_hy_via_dry(qvdry_c,qcdry_c,qrdry_c,T_c,dpdry_c,ptop,zbottom,energyhy_after)
     !print *, 'Condensation: enbefore - enafter(up to flux)', (energy_before - energy_after)/energy_before


     !this works for now
     !if(energyout > tol_energy)then
     !print *, 'Total: en - en(up to flux)', (energy_start_timestep - energy_after - energyout)/energy_start_timestep
     !print *, 'energy flux, total energy after', energyout, energy_after
     !endif

     !if(energyout > tol_energy)then
     !print *, 'Energy flux comparison:', (energyout - cl*T_c(nlev)*massout)/energyout
     !endif

     !update q fields
     ! if WITH RESPECT TO OLD PRESSURE!
     ! then this is wrong to do: dp_c = dpdry_c*(1.0 + qvdry_c + qcdry_c + qrdry_c)

#if 1
     ! update that avoids forcing_tracers code
     ! otherwise, use old pressure below
     dp_c = dpdry_c*(1.0 + qvdry_c + qcdry_c + qrdry_c)
#else
     ! update that uses forcing_tracers code: do nothing
#endif

     call convert_to_wet(qvdry_c, qcdry_c, qrdry_c, dp_c, dpdry_c, qv_c, qc_c, qr_c)

!!!! where is dp3d updated here?


  endif ! any water >0

end subroutine kessler_new






subroutine energycp_hy_via_mass(dpdry_c,dpv_c,dpc_c,dpr_c,T_c,ptop,zi_c,energy)

  real(rl), dimension(nlev), intent(in) :: dpdry_c, dpv_c, dpc_c, dpr_c
  real(rl), dimension(nlevp),intent(in) :: zi_c
  real(rl), dimension(nlev), intent(in) :: T_c
  real(rl),                  intent(in) :: ptop
  real(rl),                  intent(inout) :: energy

  real(rl) :: zbottom, pis, cpterm, Lterm
  integer  :: k

  zbottom = zi_c(nlevp)
  pis = ptop + sum(dpdry_c + dpv_c + dpc_c + dpr_c)

  !derived pressure values
  !call construct_hydro_pressure(dp_c,ptop,ppi)

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







subroutine energy_hy_via_dry(qvdry,qcdry,qrdry,tempe,dpdry,ptop,zbottom,energy)

  real(rl), dimension(nlev), intent(in)    :: qrdry
  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,tempe
  real(rl), dimension(nlev), intent(in)    :: dpdry
  real(rl),                  intent(inout) :: energy
  real(rl),                  intent(in)    :: zbottom, ptop

  integer  :: k
  real(rl) :: cpterm, Lterm, pis
  real(rl), dimension(nlev) :: dpi

  dpi = dpdry*(1.0 + qvdry + qcdry + qrdry)
  pis = ptop + sum(dpi)

  !dont do 1/g term
  energy = zbottom*pis*gravit
  do k=1,nlev
    cpterm = cpdry + cpv * qvdry(k) + cl * (qcdry(k) + qrdry(k))

    Lterm = (latvap+latice) * qvdry(k) + latice * (qcdry(k) + qrdry(k))

    energy = energy + dpdry(k)*( tempe(k)*cpterm + Lterm )
  enddo

end subroutine energy_hy_via_dry




subroutine energy_nh_via_dry(qvdry,qcdry,qrdry,tempe,dpdry,ptop,zbottom,p_c,energy)

  real(rl), dimension(nlev), intent(in)    :: qrdry, p_c
  real(rl), dimension(nlev), intent(inout) :: qvdry,qcdry,tempe
  real(rl), dimension(nlev), intent(in)    :: dpdry
  real(rl),                  intent(inout) :: energy
  real(rl),                  intent(in)    :: zbottom, ptop

  integer  :: k
  real(rl) :: cpterm, Lterm, rstar, nhterm, pis
  real(rl), dimension(nlev) :: pnh, dpi
  real(rl), dimension(nlev) :: ppi

  dpi = dpdry*(1.0 + qvdry + qcdry + qrdry)

  !derived pressure values
  call construct_hydro_pressure(dpi,ptop,ppi)
  pnh = p_c
  pis = sum(dpi)+ptop

  !dont do 1/g term
  energy = zbottom*pis*gravit
  do k=1,nlev
    cpterm = cpdry + cpv * qvdry(k) + cl * (qcdry(k) + qrdry(k))

    Lterm = (latvap+latice) * qvdry(k) + latice * (qcdry(k) + qrdry(k))

    rstar = rdry + rvapor*qvdry(k)
 
    nhterm = rstar*tempe(k)*(ppi(k)/pnh(k)-1.0)

    energy = energy + dpdry(k)*( tempe(k)*cpterm + Lterm + nhterm )
  enddo

end subroutine energy_nh_via_dry



























!!!!!!!!!!! copy from dcmip16 FIX THAT!
subroutine qsat_kessler2(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_const1 / p * exp( bubble_const2 * (T - bubble_const3) / ( T - bubble_const4 ) )
end subroutine qsat_kessler2

subroutine qsat_rj2(p, T, qsat)
  real(rl),         intent(out):: qsat
  real(rl),         intent(in) :: p, T
  qsat = bubble_epsilo * bubble_e0 / p * &
         exp(-(latvap/rvapor) * ((1.0/T)-(1.0/bubble_t0_const)))
end subroutine qsat_rj2


end module newphysics
#endif
