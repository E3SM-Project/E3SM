
subroutine engy_tdif(cwava   ,w       ,t       ,tm1     ,pdel    , &
                     difft   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to del-T integral
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
  implicit none
!
!------------------------------Arguments--------------------------------
!
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: t   (plon,plev)      ! temperature
  real(r8), intent(in)  :: tm1 (plon,plev)      ! temperature (previous timestep)
  real(r8), intent(in)  :: pdel(plon,plev)      ! pressure diff between interfaces
  real(r8), intent(out) :: difft                ! accumulator
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
  difft = 0._r8
!
! Compute mass integral
!
  do k=1,plev
     do i=1,nlon
        difft = difft + pdel(i,k)
     end do
  end do

  difft = difft*const

  return
end subroutine engy_tdif
