!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
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

   public :: grid_avg_variable, samples_of_variable, stat_file

   ! These are used in a 2D or 3D host model to output multiple columns
   ! Set clubb_i and clubb_j according to the column within the host model;
   ! The indices must not exceed nlon (for i) or nlat (for j).
   integer, save, public :: clubb_i = 1, clubb_j = 1
!$omp threadprivate(clubb_i, clubb_j)

   private ! Default scope

  ! Structures to hold the description of a variable

   type grid_avg_variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:,:,:), pointer :: ptr

     character(len = 30) :: name        ! Variable name
     character(len = 100) :: description ! Variable description
     character(len = 25) :: units       ! Variable units

     integer :: indx ! NetCDF module Id for var / GrADS index

     logical :: l_silhs ! If true, we sample this variable once for each SILHS
                        ! sample point per timestep, rather than just once per
                        ! timestep.
   end type grid_avg_variable

   type samples_of_variable
     ! Pointer to the array
     real(kind=stat_rknd), dimension(:,:,:,:), pointer :: ptr

     character(len = 30) :: name        ! Variable name
     character(len = 100) :: description ! Variable description
     character(len = 25) :: units       ! Variable units

     integer :: indx ! NetCDF module Id for var / GrADS index

     logical :: l_silhs ! If true, we sample this variable once for each SILHS
                        ! sample point per timestep, rather than just once per
                        ! timestep.
   end type samples_of_variable

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

     ! NetCDF datafile dimensions indices (Samp*Id for SILHS samples)
     integer ::  & 
       SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId, &
       SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId

     ! Grid information

     integer :: ia, iz  ! Vertical extent

     integer :: nlat, nlon ! The number of points in the X and Y

     ! Number of SILHS samples (i.e. subcolumns).  Initialized to zero
     ! to be safe, but will be updated if appropriate
     integer :: nsamp = 0

     real( kind = core_rknd ), dimension(:), allocatable ::  & 
       z ! Height of vertical levels [m]

     ! Time information

     integer :: day, month, year ! Date of starting time

     real( kind = core_rknd ), dimension(:), allocatable :: & 
       lat_vals, & ! Latitude                   [Degrees N]
       lon_vals, & ! Longitude                  [Degrees E]
       samp_idx   ! SILHS subcolumn index

     real( kind = core_rknd ) :: & 
       dtwrite ! Interval between output    [Seconds]

     real( kind = time_precision ) ::  & 
       time    ! Start time                 [Seconds]

     ! Statistical Variables

     integer :: nvar  ! Number of variables for this file

     type (grid_avg_variable), dimension(:), allocatable ::  &
       grid_avg_var ! List and variable description

     type (samples_of_variable), dimension(:), allocatable :: &
       samples_of_var

   end type stat_file

 end module stat_file_module
