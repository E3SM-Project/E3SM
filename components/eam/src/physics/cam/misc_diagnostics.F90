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

subroutine compute_cape_diags( state, pbuf, pcols, pver, cape_out, dcape_out )
!-------------------------------------------------------------------------------------------
! Purpose: 
! - CAPE, the convecitve available potential energy
! - dCAPE, change in CAPE 
! - dCAPEe, dCAPE caused by environment change
! - dCAPEp, dCAPE caused by parcel property change 
!
! History: 
! First version (2021): algorithm by Xiaoliang Song and Guang Zhang; code by Hui Wan.
! Update to EAMv2 (2022).
!-------------------------------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols

  real(r8),                 intent(out) ::  cape_out(pcols)
  real(r8),optional,        intent(out) :: dcape_out(pcols,3)

  ! local variables used for providing the same input to different calls of subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8) ::   zs(pcols)
  integer  :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  logical :: iclosure = .true.  ! set to .true. to avoid interference with trig_dcape

  ! variables that distinguish different calls of buoyan_dilute 

  real(r8),pointer ::    qv_new(:,:)  ! new qv   from current state
  real(r8),pointer ::  temp_new(:,:)  ! new temp from current state

  logical :: use_old_parcel_tq    ! whether or not to use old launching level and parcel T, q when calculating CAPE

  integer  ::    mx_new(pcols)  ! index of launching level in new environment,  calculated by buoyan_dilute
  real(r8) ::  q_mx_new(pcols)  ! new specific humidity at new launching level, calculated by buoyan_dilute
  real(r8) ::  t_mx_new(pcols)  ! new temperature       at new launching level, calculated by buoyan_dilute

  integer, pointer ::   mx_old(:)  ! old launching level from pbuf
  real(r8),pointer :: q_mx_old(:)  ! old qv   at launching level from pbuf
  real(r8),pointer :: t_mx_old(:)  ! old temp at launching level from pbuf

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.
  integer  ::   zmx(pcols)      ! launching level index 

  ! CAPE calculated using different combinations of environmental profiles and parcel properties

  real(r8)         :: cape_new_pcl_new_env(pcols) ! cape in new environment (assuming new launching level and parcel properties) 
  real(r8)         :: cape_old_pcl_new_env(pcols) ! cape in new environment (assuming old launching level and parcel properties) 

  real(r8),pointer :: cape_old_pcl_old_env(:)     ! cape in old environment (assuming old launching level and parcel proper.
                                                  ! This variable is a pointer because the values are saved in pbuf

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  !-----------------------------------------------------------------------------------------
  ! If diagnosing dCAPE and its decomposition, retrieve old CAPE and old parcel properties 
  !-----------------------------------------------------------------------------------------
  if (PRESENT(dcape_out)) then 

    idx = pbuf_get_index('CAPE_old_dCAPEd')  ; call pbuf_get_field( pbuf, idx, cape_old_pcl_old_env )
    idx = pbuf_get_index('Q_mx_old_dCAPEd')  ; call pbuf_get_field( pbuf, idx, q_mx_old )
    idx = pbuf_get_index('T_mx_old_dCAPEd')  ; call pbuf_get_field( pbuf, idx, t_mx_old )
    idx = pbuf_get_index(  'mx_old_dCAPEd')  ; call pbuf_get_field( pbuf, idx,   mx_old )

  end if 

  !-----------------------------------
  ! Time-independent quantities 
  !-----------------------------------
  msg = limcnv - 1  ! limcnv is the top interface level limit for convection

  idx = pbuf_get_index('tpert') ; call pbuf_get_field( pbuf, idx, tpert )

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

  !-------------------------------------------------------------------
  ! Temperature, specific humidity, and pressure in new environment
  !-------------------------------------------------------------------
  qv_new   => state%q(:,:,1)
  temp_new => state%t

  pmid_in_hPa(1:ncol,:) = state%pmid(1:ncol,:) * 0.01_r8
  pint_in_hPa(1:ncol,:) = state%pint(1:ncol,:) * 0.01_r8

  !------------------------------------------------------------------------
  ! Calculate CAPE using the new state; also return launching level index
  ! and T, qv values at (new) launching level
  !------------------------------------------------------------------------
  iclosure = .true.
  call buoyan_dilute(lchnk ,ncol,                      &! in
                     qv_new, temp_new,                 &! in  !!
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap,                           &! in
                     cape_new_pcl_new_env,             &! out !!
                     pblt,                             &! in
                     zlcl, zlel, zlon,                 &! out
                     mx_new,                           &! out !!
                     rair, gravit, cpair, msg, tpert,  &! in
                     iclosure,                         &! in
                     use_input_parcel_tq_in = .false., &! in  !!
                     q_mx = q_mx_new,                  &! out !!
                     t_mx = t_mx_new                   )! out !!

  cape_out(:ncol) = cape_new_pcl_new_env(:ncol)
 
  !---------------------------------------------------------------------
  ! If requested output is CAPE only, then we are done. Otherwise,
  ! diagnose dCAPE and its decomposition (dCAPEe and dCAPEp) 
  !---------------------------------------------------------------------
  if (PRESENT(dcape_out)) then 

    !-----------------------------------------------------------------
    ! dCAPE is the difference between the new CAPE calculated above 
    ! and the old CAPE retrieved from pbuf

     dcape_out(:ncol,1) = cape_new_pcl_new_env(:ncol) - cape_old_pcl_old_env(:ncol)

    !-----------------------------------------------------------------
    ! Calculate cape_old_pcl_new_env using
    !  - new state (T, qv profiles)
    !  - old launching level and parcel T, qv

    iclosure = .true.
    call buoyan_dilute(lchnk ,ncol,                      &! in
                       qv_new, temp_new,                 &! in  !!!
                       pmid_in_hPa, zmid_above_sealevel, &! in
                       pint_in_hPa,                      &! in
                       ztp, zqstp, ztl,                  &! out
                       latvap,                           &! in
                       cape_old_pcl_new_env,             &! out !!!
                       pblt,                             &! in
                       zlcl, zlel, zlon,                 &! out
                       zmx,                              &! out 
                       rair, gravit, cpair, msg, tpert,  &! in
                       iclosure,                         &! in
                       dcapemx = mx_old,                 &! in  !!!
                       use_input_parcel_tq_in = .true.,  &! in  !!!
                       q_mx = q_mx_old,                  &! in  !!!
                       t_mx = t_mx_old                   )! in  !!!

    ! dCAPEp = CAPE(new parcel, new env) - CAPE( old parcel, new env)
    dcape_out(:ncol,2) = cape_new_pcl_new_env(:ncol) - cape_old_pcl_new_env(:ncol)

    ! dCAPEe = CAPE(old parcel, new env) - CAPE( old parcel, old env)
    dcape_out(:ncol,3) = cape_old_pcl_new_env(:ncol) - cape_old_pcl_old_env(:ncol)

    !-----------------------------------------------------------------
    ! Update "old" CAPE and parcel properties in pbuf for next call

    cape_old_pcl_old_env(:ncol) = cape_new_pcl_new_env(:ncol)

    t_mx_old(:ncol) = t_mx_new(:ncol)
    q_mx_old(:ncol) = q_mx_new(:ncol)
      mx_old(:ncol) =   mx_new(:ncol)

  end if 

 end subroutine compute_cape_diags

 subroutine dcape_diags_register( pcols )

   use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4

   integer,intent(in) :: pcols
   integer :: idx

   call pbuf_add_field('CAPE_old_dCAPEd', 'global', dtype_r8,(/pcols/), idx)
   call pbuf_add_field(  'mx_old_dCAPEd', 'global', dtype_i4,(/pcols/), idx)
   call pbuf_add_field('Q_mx_old_dCAPEd', 'global', dtype_r8,(/pcols/), idx)
   call pbuf_add_field('T_mx_old_dCAPEd', 'global', dtype_r8,(/pcols/), idx)

 end subroutine dcape_diags_register

 subroutine dcape_diags_init( pbuf, pver )

   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

   implicit none

   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: pver

   real(r8), pointer ::  ptr2d(:,:)
   real(r8), pointer ::  ptr1d(:)
   integer,  pointer :: iptr1d(:)

   call pbuf_get_field( pbuf, pbuf_get_index('CAPE_old_dCAPEd'), ptr1d ); ptr1d = 0._r8 
   call pbuf_get_field( pbuf, pbuf_get_index('Q_mx_old_dCAPEd'), ptr1d ); ptr1d = 5.e-3_r8
   call pbuf_get_field( pbuf, pbuf_get_index('T_mx_old_dCAPEd'), ptr1d ); ptr1d = 273._r8
   call pbuf_get_field( pbuf, pbuf_get_index(  'mx_old_dCAPEd'), iptr1d); iptr1d= pver

 end subroutine dcape_diags_init
!---------------------------

end module misc_diagnostics
