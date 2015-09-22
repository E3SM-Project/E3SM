
subroutine flxint (w       ,flx     ,flxlat  ,nlon    )
!----------------------------------------------------------------------- 
! 
! Purpose: Calculate contribution of current latitude to energy flux integral
! 
! Method: 
! 
! Author: Jerry Olson
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
   real(r8), intent(in) :: w          ! gaussian weight this latitude
   real(r8), intent(in) :: flx(plon)  ! energy field

   integer, intent(in) :: nlon        ! number of longitudes

   real(r8), intent(out) :: flxlat    ! accumulator for given latitude
!
! Local variables
!
   integer :: i                       ! longitude index
!
!-----------------------------------------------------------------------
!
   flxlat = 0._r8
!
   do i=1,nlon
      flxlat = flxlat + flx(i)
   end do
!
! Integration factor (the 0.5 factor arises because gaussian weights
! sum to 2)
!
   flxlat = flxlat*w*0.5_r8/real(nlon,r8)
!
   return
end subroutine flxint
