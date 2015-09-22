
subroutine engy_te(cwava   ,w       ,t       ,u      ,v        , &
                   phis    ,pdel    ,engy    , nlon  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to total energy
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plon
  use physconst, only: cpair

  implicit none
!
!------------------------------Arguments--------------------------------
!
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: t   (plon,plev)      ! temperature
  real(r8), intent(in)  :: u   (plon,plev)      ! u-component
  real(r8), intent(in)  :: v   (plon,plev)      ! v-component
  real(r8), intent(in)  :: phis(plon)           ! Geopotential
  real(r8), intent(in)  :: pdel(plon,plev)      ! pressure diff between interfaces
  real(r8), intent(out) :: engy                 ! accumulator
!
!---------------------------Local variables-----------------------------
!
  integer i,k               ! longitude, level indices
  real(r8) const            ! temporary constant
!
!-----------------------------------------------------------------------
!
! Integration factor (the 0.5 factor arises because gaussian weights sum to 2)
!
  const = cwava*w*0.5_r8
  engy = 0._r8
!
  do k=1,plev
     do i=1,nlon
        engy = engy + ( cpair*t(i,k)                          + &
                        0.5_r8*( u(i,k)*u(i,k) + v(i,k)*v(i,k) ) + &
                        phis(i) )*pdel(i,k)
     end do
  end do

  engy = engy*const

  return
end subroutine engy_te
