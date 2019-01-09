MODULE TIMEINFO
!  Contains the variables describing the start and current time of the model run.

   USE shr_kind_mod, only: r8 => shr_kind_r8
   
IMPLICIT NONE

! Note: Fractional Julian day is specified as the day of year during iyr (1.0 = 00Z Jan 1)

    integer :: iyr              ! Year of current time step
    
    real (kind=r8) :: rjday0    ! Fractional Julian day at the start of the model simulation
    real (kind=r8) :: rjday     ! Fractional Julian day of the current model time step
    
END MODULE timeinfo