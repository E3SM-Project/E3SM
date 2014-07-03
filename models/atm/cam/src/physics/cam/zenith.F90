subroutine zenith(calday  ,clat    , clon   ,coszrs  ,ncol    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cosine of solar zenith angle for albedo and radiation
!   computations.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_orb_mod
   use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of positions
   real(r8), intent(in) :: calday              ! Calendar day, including fraction
   real(r8), intent(in) :: clat(ncol)          ! Current centered latitude (radians)
   real(r8), intent(in) :: clon(ncol)          ! Centered longitude (radians)
!
! Output arguments
!
   real(r8), intent(out) :: coszrs(ncol)       ! Cosine solar zenith angle
!
!---------------------------Local variables-----------------------------
!
   integer i         ! Position loop index
   real(r8) delta    ! Solar declination angle  in radians
   real(r8) eccf     ! Earth orbit eccentricity factor
!
!-----------------------------------------------------------------------
!
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf      )
!
! Compute local cosine solar zenith angle,
!
   do i=1,ncol
      coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta )
   end do

end subroutine zenith
