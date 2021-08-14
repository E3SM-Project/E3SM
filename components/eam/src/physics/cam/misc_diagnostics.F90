module misc_diagnostics

use shr_kind_mod,   only: r8 => shr_kind_r8

implicit none
public

contains


!------------------------------------------------
! saturation specific humidity wrt ice.
! The calculation in this subroutine follows 
! what is used in the model for 
!  - history output
!  - ice nucleation parameterization
! 
subroutine qsat_ice( ncol, pver, tair, pair, qsati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)

  real(r8),intent(out) :: qsati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qsati = qsatw*esi/esl

end subroutine qsat_ice
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt water (liquid) given as mixing ratio
!
subroutine supersat_q_water( ncol, pver, tair, pair, qv, qssatw )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssatw(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  qssatw = qv-qsatw

end subroutine supersat_q_water
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt ice given as mixing ratio
!
subroutine supersat_q_ice( ncol, pver, tair, pair, qv, qssati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qssati = qv - qsatw*esi/esl


end subroutine supersat_q_ice
!----------------------------------

subroutine relhum_water_percent( ncol, pver, tair, pair, qv,  rhw_percent )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhw_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  rhw_percent = qv/qsatw*100._r8

end subroutine relhum_water_percent


!----------------------------------
subroutine relhum_ice_percent( ncol, pver, tair, pair, qv,  rhi_percent )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhi_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  rhi_percent = 100._r8* qv/qsatw* esl/esi

end subroutine relhum_ice_percent

subroutine compute_cape( state, pbuf, pcols, pver, cape )
!----------------------------------------------------------------------
! Purpose: compute CAPE (convecitve available potential energy)
!          using subroutine buoyan_dilute from the ZM deep convection
!          parameterization (module file zm_conv.F90)
! History: first version by Hui Wan, 2021-05
!----------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  real(r8),                intent(out) :: cape(pcols)

  ! local variables used for providing input to subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8),pointer ::    qv(:,:)
  real(r8),pointer ::  temp(:,:)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.
  integer  :: zmaxi(pcols)      ! index of level with largest moist static energy.
  logical iclosure              ! switch on sequence of call to buoyan_dilute to derive DCAPE

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  msg = limcnv - 1  ! limcnv is the top interface level limit for convection

  idx = pbuf_get_index('tpert') ; call pbuf_get_field( pbuf, idx, tpert )

  pmid_in_hPa(1:ncol,:) = state%pmid(1:ncol,:) * 0.01_r8
  pint_in_hPa(1:ncol,:) = state%pint(1:ncol,:) * 0.01_r8

  qv   => state%q(:,:,1)
  temp => state%t

  ! Surface elevation (m) is needed to calculate height above sea level (m) 
  ! Note that zm (and zi) stored in state are height above surface. 
  ! The layer midpoint height provided to buoyan_dilute is height above sea level. 

  zs(1:ncol) = state%phis(1:ncol)/gravit

  !------------------------------------------
  ! Height above sea level at layer midpoints
  !------------------------------------------
  do kk = 1,pver
     zmid_above_sealevel(1:ncol,kk) = state%zm(1:ncol,kk)+zs(1:ncol)
  end do

  !--------------------------
  ! layer index for PBL top
  !--------------------------
  idx = pbuf_get_index('pblh')  ; call pbuf_get_field( pbuf, idx, pblh )

  pblt(:) = pver
  do kk = pver-1, msg+1, -1
     where( abs(zmid_above_sealevel(:ncol,kk)-zs(:ncol)-pblh(:ncol))   &
            < (state%zi(:ncol,kk)-state%zi(:ncol,kk+1))*0.5_r8       ) & 
     pblt(:ncol) = kk
  end do

  !----------------
  ! Calculate CAPE
  !----------------
  iclosure = .true. !standard calculation, scanning for launching level up to 600 hPa
  call buoyan_dilute(lchnk ,ncol, qv, temp,            &! in
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap, cape, pblt,               &! in, out, in
                     zlcl, zlel, zlon, zmaxi,          &! out
                     rair, gravit, cpair, msg, tpert, iclosure )! in

 end subroutine compute_cape
!---------------------------

end module misc_diagnostics
