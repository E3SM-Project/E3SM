module orbit

use shr_kind_mod, only: r8 => shr_kind_r8
use seq_infodata_mod, only: seq_infodata_type, seq_infodata_SolarDecl
use seq_infodata_mod, only: seq_infodata_UseEquationOfTime
use shr_sys_mod, only: shr_sys_abort
implicit none

type(seq_infodata_type), private, save :: solar_context
logical, private, save :: solar_context_initialized = .false.

contains

subroutine orbit_set_context(infodata)
   type(seq_infodata_type), intent(in) :: infodata
   solar_context = infodata
   solar_context_initialized = .true.
end subroutine orbit_set_context

subroutine orbital_decl(calday, delta, eccf, offset, eot_seconds)
   use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
   use shr_orb_mod, only: shr_orb_decl
   use time_manager, only: timemgr_get_calendar_cf, get_curr_date

   real(r8), intent(in) :: calday
   real(r8), intent(out) :: delta, eccf
   integer, optional, intent(in) :: offset
   real(r8), optional, intent(out) :: eot_seconds
   integer :: year, month, day, seconds

   if (present(eot_seconds)) eot_seconds = 0.0_r8
   if (timemgr_get_calendar_cf() == 'gregorian' .or. seq_infodata_UseEquationOfTime(solar_context)) then
      if (.not. solar_context_initialized) then
         call shr_sys_abort('EAM Gregorian solar geometry requires initialized orbital context')
      end if
      call get_curr_date(year, month, day, seconds, offset=offset)
      call seq_infodata_SolarDecl(solar_context, year * 10000 + month * 100 + day, &
           real(seconds, r8), timemgr_get_calendar_cf(), delta, eccf, eot_seconds)
   else
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
   end if
end subroutine orbital_decl

subroutine zenith(calday, clat, clon, coszrs, ncol, dt_avg, uniform_angle )
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
   real(r8), intent(in), optional :: dt_avg    ! if present, time step to use for the shr_orb_cosz calculation
   real(r8), intent(in), optional :: uniform_angle ! if present, use a globally uniform zenith angle [radians]
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
   real(r8) eot_seconds
!
!-----------------------------------------------------------------------
!
   call orbital_decl(calday, delta, eccf, eot_seconds=eot_seconds)
!
! Compute local cosine solar zenith angle,
!
   do i=1,ncol
      if (present(uniform_angle)) then
         coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta, dt_avg, &
              uniform_angle=uniform_angle, eot_seconds=eot_seconds )
      else
         coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta, dt_avg, eot_seconds=eot_seconds )
      end if
   end do

end subroutine zenith
end module orbit
