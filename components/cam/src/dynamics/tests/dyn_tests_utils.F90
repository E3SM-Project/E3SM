module dyn_tests_utils
!----------------------------------------------------------------------- 
! 
! Utility data (and code) for dynamics testing
!
! The public items in this module are items used both by internal code
!   (e.g., analytic initial conditions) and by infrastructure which uses
!   the internal code (e.g., read_inidat). They cannot be members of the 
!   internal code because that is conditionally compiled.
!
!-----------------------------------------------------------------------


  implicit none
  private
  save

  integer, parameter :: vc_moist_pressure = 0 ! Moist pressure vertical coord
  integer, parameter :: vc_dry_pressure   = 1 ! Dry pressure vertical coord
  integer, parameter :: vc_height         = 2 ! Height vertical coord
  public :: vc_moist_pressure, vc_dry_pressure, vc_height

end module dyn_tests_utils
