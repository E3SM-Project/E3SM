module orbit

contains

subroutine zenith(calday  ,clat    , clon   ,coszrs  ,ncol, dt_avg    )
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
   use shr_const_mod, only: SHR_CONST_PI
   use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp, constant_zenith_deg
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of positions
   real(r8), intent(in) :: calday              ! Calendar day, including fraction
   real(r8), intent(in) :: clat(ncol)          ! Current centered latitude (radians)
   real(r8), intent(in) :: clon(ncol)          ! Centered longitude (radians)
   real(r8), intent(in), optional :: dt_avg    ! if present, time step to use for the shr_orb_cosz calculation
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
   real(r8) constant_cosz ! consine of constant zenith angle
!
!-----------------------------------------------------------------------
!
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf      )
!
! Compute local cosine solar zenith angle,
!
   if (constant_zenith_deg>0) then
      constant_cosz = cos( constant_zenith_deg * SHR_CONST_PI/180. )
      do i=1,ncol
         coszrs(i) = constant_cosz
      end do
   else
      do i=1,ncol
         coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta, dt_avg )
      end do
   end if

end subroutine zenith
end module orbit
