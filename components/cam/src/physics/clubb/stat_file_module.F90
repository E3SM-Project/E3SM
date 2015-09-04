!-------------------------------------------------------------------------------
! $Id: stat_file_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module stat_file_module
 

! Description:
!   Contains two derived types for describing the contents and location of
!   either NetCDF or GrADS files.
!-------------------------------------------------------------------------------
   use clubb_precision, only: & 
       stat_rknd,  & ! Variable
       time_precision, &
       core_rknd
 
   implicit none

   public :: variable, stat_file

   private ! Default scope
   
  ! Structure to hold the description of a variable

   type variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:,:,:), pointer :: ptr 

     character(len = 30) :: name        ! Variable name
     character(len = 100) :: description ! Variable description
     character(len = 20) :: units       ! Variable units

     integer :: indx ! NetCDF module Id for var / GrADS index
   end type variable

  ! Structure to hold the description of a NetCDF output file
  ! This makes the new code as compatible as possible with the
  ! GrADS output code

   type stat_file

     ! File information

     character(len = 200) ::  & 
       fname,   & ! File name without suffix
       fdir    ! Path where fname resides

     integer :: iounit  ! This number is used internally by the 
                        ! NetCDF module to track the data set, or by 
                        ! GrADS to track the actual file unit.
     integer :: &
       nrecord, & ! Number of records written
       ntimes     ! Number of times written

     logical :: &
       l_defined,  &  ! Whether nf90_enddef() has been called
       l_byte_swapped ! Is this a file in the opposite byte ordering?

     ! NetCDF datafile dimensions indices
     integer ::  & 
       LatDimId, LongDimId, AltDimId, TimeDimId, & 
       LatVarId, LongVarId, AltVarId, TimeVarId

     ! Grid information

     integer :: ia, iz  ! Vertical extent

     integer :: nlat, nlon ! The number of points in the X and Y

     real( kind = core_rknd ), dimension(:), pointer ::  & 
       z ! Height of vertical levels [m]

     ! Time information

     integer :: day, month, year ! Date of starting time

     real( kind = core_rknd ), dimension(:), pointer :: & 
       rlat, & ! Latitude                   [Degrees N]
       rlon    ! Longitude                  [Degrees E]

     real(kind=time_precision) :: & 
       dtwrite ! Interval between output    [Seconds]

     real(kind=time_precision) ::  & 
       time    ! Start time                 [Seconds]

     ! Statistical Variables

     integer :: nvar  ! Number of variables for this file

     type (variable), dimension(:), pointer ::  & 
       var ! List and variable description

   end type stat_file

 end module stat_file_module
