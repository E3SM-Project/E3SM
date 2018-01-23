MODULE TIMEINFO
!  Contains the variables describing the start and current time of the model run.

   USE shr_kind_mod, only: dbl_kind => shr_kind_r8
   
IMPLICIT NONE

! Note: Fractional Julian day is specified as the day of year during iyr (0.0 = 00Z Jan 1)

    integer :: iyr                    ! Year of current time step
    integer :: imonth                 ! Month of current time step
    integer :: iday                   ! Calendar day of current time step
    
    real (kind=dbl_kind) :: rjday0    ! Fractional Julian day at the start of the model simulation
    real (kind=dbl_kind) :: rjday     ! Fractional Julian day of the current model time step
    real (kind=dbl_kind) :: utc_time  ! UTC time of the current model time step 
    
END MODULE timeinfo