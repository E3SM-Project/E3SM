module scream_scorpio_interface

!==============================================================================!
! This module handles the Fortran interface to the PIO library of input/output
! subroutines.  The essential set of steps to enable and use PIO for creating
! output are as follows:
!                                  OUTPUT
! Initialization:
! 1) Gather the pio_subsystem and pio_iotype information for the EAM component
! as assigned by the component coupler.
!    This is accomplished during eam_init_pio_1 by calling
!    'eam_init_pio_subsystem'
! 2) For each output file "create" a file in PIO and record the unique file
! descriptor. 
!    This is accomplished during the eam_init_pio_1 by calling
!    'eam_pio_createfile'
! 3) For each output file define the "header" information, which is essentially
! a set of metadata strings that describe the output file.
!    This is accomplished during eam_init_pio_1 by calling 'eam_pio_createHeader'
! 4) Define all of the dimensions that variables in this file will be defined
! on.  Examples would time, lat, lon, vertical coordinate, # of consituents,
! etc.
!    This is accomplished during 'register_dimension' by which calls 'PIO_def_dim'
! 5) Define all of the variables that will be written to this file.  This
! includes any dimension that should also be defined as a variable, such as lat,
! lon, time.
!    This is accomplished during 'register_variable' which calls to 'PIO_def_var'
! 6) Determine the unique PIO identifiers for each domain decomposition.  This
! is essentially a decomposition of what will be written to the file, multiple
! variables can have the same decomposition.  Every arrangement of dimensions
! requires a pio decomposition.
!    This is accomplished during 'register_variable' by calling
!    'PIO_initdecomp'
! 7) Close the PIO file definition step.  In other words, tell PIO that all of
! the dimensions, variables and decompositions associated with this output have
! been defined and no new ones will be added.
!    This is accomplished during eam_init_pio_2 by calling 'PIO_enddef'
!
! Writing Output:
! 
! Finalization: 
!==============================================================================!

 
  use shr_pio_mod,  only: shr_pio_getrearranger
  use shr_kind_mod,   only: rtype=>shr_kind_r8, rtype_short=>shr_kind_r4  ! Need to change this to use scream real types
  use piolib_mod, only : PIO_init, PIO_finalize, PIO_createfile, PIO_closefile, &
      PIO_initdecomp, PIO_freedecomp, PIO_syncfile, PIO_openfile, PIO_setframe
  use pio_types,  only : iosystem_desc_t, file_desc_t, &
      pio_noerr, PIO_iotype_netcdf, var_desc_t, io_desc_t, PIO_int, &
      pio_clobber, PIO_nowrite, PIO_unlimited, pio_global, PIO_real, &
      PIO_double
  use pio_kinds,  only : PIO_OFFSET_KIND, i4 
  use pio_nf,     only : PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef, PIO_inq_dimid, &
                         PIO_inq_dimlen
  use piodarray,  only : PIO_write_darray, PIO_read_darray 
  use pionfatt_mod, only : PIO_put_att   => put_att
  use pionfput_mod, only : PIO_put_var   => put_var

  implicit none
  save
                     
  public :: eam_init_pio_1,       & ! Get pio subsystem info from main code, create PIO files
            eam_init_pio_2,       & ! Register variables and dimensions with PIO files
            eam_history_finalize, & ! Run any final PIO commands (currently unused)
            eam_history_write,    & ! Write updated data for each variable to file
            register_variable,    & ! Register a variable with a particular pio output file  
            register_dimension      ! Register a dimension with a particular pio output file
 
  private :: errorHandle

  ! Universal PIO variables for the module
  integer               :: pio_mpicom  
  integer               :: pio_iotype
  type(iosystem_desc_t), pointer, public :: pio_subsystem => null()
  integer               :: pio_rearranger
  integer               :: pio_ntasks
  integer               :: pio_myrank

  ! TYPES to handle history coordinates and files
  integer,parameter :: max_hcoordname_len = 16
  integer,parameter :: max_chars = 256
  integer,parameter :: max_hvarname_len = 16
  integer,parameter :: max_hvar_dimlen  = 5

!----------------------------------------------------------------------
  type, public :: hist_coord_t
    character(len=max_hcoordname_len) :: name = ''          ! coordinate name
    integer                  :: dimsize = 0                 ! size of dimension
    integer                  :: dimid                       ! Unique PIO Id for this dimension
    character(len=max_chars) :: long_name = ''              ! 'long_name' attribute
    character(len=max_chars) :: units = ''                  ! 'units' attribute
    character(len=max_chars) :: bounds_name = ''            ! 'bounds' attribute (& name of bounds variable)
    character(len=max_chars) :: standard_name = ''          ! 'standard_name' attribute
    character(len=4)         :: positive = ''               ! 'positive' attribute ('up' or 'down')
    integer,  pointer        :: integer_values(:) => null() ! dim values if integer
    real(rtype), pointer     :: real_values(:) => null()    ! dim values if real
    real(rtype), pointer     :: bounds(:,:) => null()       ! dim bounds
    logical                  :: integer_dim                 ! .true. iff dim has integral values
    logical                  :: vertical_coord              ! .true. iff dim is vertical
  end type hist_coord_t
!----------------------------------------------------------------------
  type, public :: hist_var_t
    character(len=max_hvarname_len) :: name   ! coordinate name
    character(len=max_chars) :: long_name     ! 'long_name' attribute
    character(len=max_chars) :: units         ! 'units' attribute
    type(var_desc_t) :: piovar                ! netCDF variable ID
    integer          :: dtype                 ! data type
    integer          :: numdims               ! Number of dimensions in out field
    type(io_desc_t), pointer  :: iodesc       ! PIO decomp associated with this variable
    integer, allocatable :: compdof(:)        ! Global locations in output array for this process
    integer, allocatable :: dimid(:)          ! array of PIO dimension id's for this variable
    integer, allocatable :: dimlen(:)         ! array of PIO dimension lengths for this variable
  end type hist_var_t
!----------------------------------------------------------------------
  ! The iodesc_list allows us to cache existing PIO decompositions
  ! The tag needs the dim lengths, the dtype and map id (+ optional permutation)
  ! Define a recursive structure because we do not know ahead of time how many
  ! decompositions will be require
  integer, parameter      :: tag_len           = 48
  type iodesc_list
    character(tag_len)          :: tag              ! Unique tag associated with this decomposition
    type(io_desc_t),    pointer :: iodesc => NULL() ! PIO - decomposition
    type(iodesc_list),  pointer :: next => NULL()   ! Needed for recursive definition
  end type iodesc_list
  ! Define the first iodesc_list 
  type(iodesc_list), target :: iodesc_list_top
!----------------------------------------------------------------------
  type hist_coord_list
    type(hist_coord_t),    pointer :: coord => NULL() ! Pointer to a history dimension structure
    type(hist_coord_list), pointer :: next => NULL()  ! Needed for recursive definition
  end type hist_coord_list
!----------------------------------------------------------------------
  type hist_var_list
    type(hist_var_t),    pointer :: var => NULL()  ! Pointer to a history variable structure
    type(hist_var_list), pointer :: next => NULL() ! Needed for recursive definition
  end type hist_var_list
!----------------------------------------------------------------------
  type pio_file_list
    type(pio_atm_output), pointer :: pio_file => NULL() ! Pointer to an atm. pio file
    type(pio_file_list),  pointer :: next => NULL()     ! Needed for recursive definition
  end type pio_file_list
  ! Define the first pio_file_list
  type(pio_file_list), target :: pio_file_list_top
!----------------------------------------------------------------------
  type, public :: pio_atm_output
        !> @brief Output filename.
        character(len=max_chars) :: filename = ""

        !> @brief Contains data identifying the file.
        type(file_desc_t)     :: pioFileDesc

        !> @brief Number of output dimensions, and counter to track them during
        !  registration
        integer               :: DimCounter = 0
        !> @brief Recursive list of variables
        type(hist_coord_list)   :: coord_list_top

        !> @brief Number of output variable and counter to track them during
        !  registration
        integer               :: VarCounter = 0
        !> @brief Recursive list of variables
        type(hist_var_list)   :: var_list_top

        !> @brief Number of history records on this file
        integer               :: numRecs

        !> @brief Coordinate Dimensions Array
        type(hist_coord_t), allocatable :: dimensions(:)

  end type pio_atm_output

!----------------------------------------------------------------------
  interface grid_write_data_array
    module procedure grid_write_darray_1d_int
    module procedure grid_write_darray_2d_int
    module procedure grid_write_darray_3d_int
    module procedure grid_write_darray_1d_double
    module procedure grid_write_darray_2d_double
    module procedure grid_write_darray_3d_double
    module procedure grid_write_darray_1d_real
    module procedure grid_write_darray_2d_real
    module procedure grid_write_darray_3d_real
  end interface
!----------------------------------------------------------------------
contains
!=====================================================================!
  ! Subroutine to initialize the pio interface.  May be deprecated in the future
  ! with the C++ code directly a) calling init_pio_subsystem and b) creating all
  ! new files.
  subroutine eam_init_pio_1(mpicom,atm_id)
    integer, intent(in) :: mpicom  ! MPI communication group for ATM
    integer, intent(in) :: atm_id  ! Unique identifier assigned by component coupler

    type(pio_atm_output), pointer :: current_atm_file => null()
    integer :: ierr
    logical :: found

    ! Gather the pio_subsystem information for the atmosphere component.
    call eam_init_pio_subsystem(mpicom,atm_id)

    ! BEGIN PLACEHOLDER FOR DEVELOPMENT TESTING ONLY
    ! TODO: Delete the following code after the Field Manager assumes control
    ! over pio file creation.
    
    ! FIRST TESTING PIO OUTPUT FILE
    ! The following two lines will be typical, assign a filename and check that
    ! the file is in fact new.
    call get_pio_atm_file("example_pio_structured.nc",current_atm_file,found)
    if (found) call errorHandle("PIO Error: Creating file: "//trim("example_pio_structured.nc")//" already exists",999)
    ! Create a temporary header for this file.  TODO: Have the field manager or
    ! scream-ad take care of adding header metadata when the C++ interface is
    ! developed.
    call eam_pio_createHeader(current_atm_file%pioFileDesc)

    ! SECOND TESTING PIO OUTPUT FILE 
    ! Allocate the dimension and variable arrays for the output file.
    call get_pio_atm_file("example_pio_structured_v2.nc",current_atm_file,found)
    if (found) call errorHandle("PIO Error: Creating file: "//trim("example_pio_structured_v2.nc")//" already exists",999)

  end subroutine eam_init_pio_1
!=====================================================================!
  subroutine eam_init_pio_2()
    ! WARNING: This code will be eventually deprecated when the C++ interface is
    ! developed and all new variables and dimensions for a file is handled by
    ! the scream-AD and not locally.
    ! THIS ROUTINE IS ONLY PRESENT FOR TESTING PURPOSES.  MUST BE DELETED IN THE
    ! FINAL PHASE OF CODE INTEGRATION.
    type(pio_atm_output), pointer :: current_atm_file => null()
    integer :: ierr
    logical :: found

    character(len=10) :: var_dim_list(3) ! Needed for cases with time and space dimensions, need to create the array of strings ahead of time for GCC.

    ! Register all dimensions with the output file
    call register_dimension("example_pio_structured.nc","x","horizontal distance",10)
    call register_dimension("example_pio_structured.nc","y","vertical distance",3)
    call register_dimension("example_pio_structured.nc","time","time",0)

    ! Register all variables with the output file
    call register_variable("example_pio_structured.nc","time","time",1,(/ "time" /), PIO_double,"t")
    call register_variable("example_pio_structured.nc","x","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    call register_variable("example_pio_structured.nc","y","answer to space and time",1,(/ "y" /), PIO_double,"y-real")

    call register_variable("example_pio_structured.nc","bar2","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    var_dim_list(1) = trim("x")
    var_dim_list(2) = trim("time")
    call register_variable("example_pio_structured.nc","real_foo","answer to space and time",2,var_dim_list(:2), PIO_double,"xt-real")
    call register_variable("example_pio_structured.nc","foo","answer to space and time",2,var_dim_list(:2),PIO_int,"xt-int")
    var_dim_list(1) = trim("x")
    var_dim_list(2) = trim("y")
    var_dim_list(3) = trim("time")
    call register_variable("example_pio_structured.nc","bar","answer to space and time",3,var_dim_list, PIO_double,"xyt-real")

    ! Finish the "definition" phase of the PIO file.  This is an essential step
    ! for the netCDF file to be ready for variables to be written to it.
    call get_pio_atm_File(trim("example_pio_structured.nc"),current_atm_file,found)
    ! It is essential that PIO_enddef is called when all dimensions and
    ! variables have been registered with a specific file.
    ierr = PIO_enddef(current_atm_file%pioFileDesc)
    call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),ierr)
   

    ! Register all dimensions with the output file
    call register_dimension("example_pio_structured_v2.nc","x","horizontal distance",10)
    call register_dimension("example_pio_structured_v2.nc","y","vertical distance",3)
    call register_dimension("example_pio_structured_v2.nc","time","time",0)
    ! Register all variables with the output file
    call register_variable("example_pio_structured_v2.nc","time","time",1,(/ "time" /), PIO_double,"t")
    call register_variable("example_pio_structured_v2.nc","x","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    call register_variable("example_pio_structured_v2.nc","y","answer to space and time",1,(/ "y" /), PIO_double,"y-real")
    var_dim_list(1) = trim("x")
    var_dim_list(2) = trim("y")
    var_dim_list(3) = trim("time")
    call register_variable("example_pio_structured_v2.nc","bar","answer to space and time",3,var_dim_list, PIO_double,"xyt-real")
    call register_variable("example_pio_structured_v2.nc","foo","answer to space and time",3,var_dim_list, PIO_double,"xyt-real")
    var_dim_list(1) = trim("y")
    var_dim_list(2) = trim("x")
    var_dim_list(3) = trim("time")
    call register_variable("example_pio_structured_v2.nc","bar_flip","answer to space and time",3,var_dim_list, PIO_double,"yxt-real")
    call register_variable("example_pio_structured_v2.nc","foo_flip","answer to space and time",3,var_dim_list, PIO_double,"yxt-real")

    call get_pio_atm_File(trim("example_pio_structured_v2.nc"),current_atm_file,found)
    ! It is essential that PIO_enddef is called when all dimensions and
    ! variables have been registered with a specific file.
    ierr = PIO_enddef(current_atm_file%pioFileDesc)
    call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),ierr)

  end subroutine eam_init_pio_2
!=====================================================================!
  ! Register a dimension with a specific pio output file 
  subroutine register_dimension(pio_atm_filename,shortname,longname,length)
    character(len=*), intent(in)        :: pio_atm_filename   ! Name of file to register the dimension on.
    character(len=*), intent(in)        :: shortname,longname ! Short- and long- names for this dimension, short: brief identifier and name for netCDF output, long: longer descriptor sentence to be included as meta-data in file.
    integer, intent(in)                 :: length             ! Length of the dimension, 0: unlimited (like time), >0 actual length of dimension

    type(pio_atm_output), pointer       :: pio_atm_file
    type(hist_coord_t), pointer         :: hist_coord
    type(hist_coord_list), pointer      :: curr=>null(), prev=>null()
    integer                             :: ierr
    character(len=150)                  :: errstr
    logical                             :: found

    ! Make sure the dimension length is reasonable
    if (length<0) call errorHandle("PIO Error: dimension "//trim(shortname)//", can't have a negative dimension length",-999)
 
    ! Find the pointer for this file
    call get_pio_atm_File(trim(pio_atm_filename),pio_atm_file,found)
    if (.not.found) call errorHandle("PIO Error: Could not register dimension "//trim(shortname)//", pio file "//trim(pio_atm_filename)//" not found",-999)
    ! Get a new dimension pointer in coord_list
    curr => pio_atm_file%coord_list_top
    do while (associated(curr))
      if (associated(curr%coord)) then 
        if(trim(curr%coord%name)==trim(shortname)) call errorHandle("PIO Error: Could not register dimension"//trim(shortname)//", already exists in file: "//trim(pio_atm_filename),-999)
      end if
      prev => curr
      curr => prev%next
    end do
    allocate(prev%next)
    curr => prev%next
    allocate(curr%coord)
    hist_coord => curr%coord
    pio_atm_file%dimcounter = pio_atm_file%dimcounter + 1
    ! Register this dimension
    hist_coord%name      = trim(shortname)
    hist_coord%long_name = trim(longname)
    hist_coord%dimsize   = length
    if (length.eq.0) then
      ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), pio_unlimited , hist_coord%dimid)
    else
      ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), length , hist_coord%dimid)
    end if
    call errorHandle("PIO ERROR: could not define dimension "//trim(shortname)//" on file: "//trim(pio_atm_filename),ierr)
    
    return
  end subroutine register_dimension
!=====================================================================!
  ! Register a variable with a specific pio output file
  subroutine register_variable(pio_atm_filename,shortname,longname,numdims,var_dimensions,dtype,pio_decomp_tag)
    character(len=*), intent(in) :: pio_atm_filename         ! Name of the file to register this variable with
    character(len=*), intent(in) :: shortname,longname       ! short and long names for the variable.  Short: variable name in file, Long: more descriptive name
    integer, intent(in)          :: numdims                  ! Number of dimensions for this variable, including time dimension
    character(len=*), intent(in) :: var_dimensions(numdims)  ! String array with shortname descriptors for each dimension of variable.
    integer, intent(in)          :: dtype                    ! datatype for this variable, REAL, DOUBLE, INTEGER, etc.
    character(len=*), intent(in) :: pio_decomp_tag           ! Unique tag for this variables decomposition type, to be used to determine if the io-decomp already exists.

    ! Local variables
    type(pio_atm_output),pointer :: pio_atm_file
    integer                      :: loc_len
    logical                      :: found
    type(io_desc_t), pointer     :: iodesc_loc
    type(hist_var_t), pointer    :: hist_var
    integer                      :: dim_ii
    integer                      :: ierr
    integer, allocatable         :: dimlen(:)
    integer                      :: total_dimlen, my_dof_len, extra_procs
    integer                      :: ii, istart, istop
    logical                      :: has_t_dim  ! Logical to flag whether this variable has a time-dimension.  This is important for the decomposition step.

    type(hist_var_list), pointer :: curr => null(), prev => null()
 
    ! Find the pointer for this file
    call get_pio_atm_File(trim(pio_atm_filename),pio_atm_file)
    ! Update the number of variables on file 
    pio_atm_file%varcounter = pio_atm_file%varcounter + 1

    ! Get a new variable pointer in var_list
    curr => pio_atm_file%var_list_top
    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(shortname)) call errorHandle("PIO Error: Could not register variable"//trim(shortname)//", already exists in file: "//trim(pio_atm_filename),-999)
      end if
      prev => curr
      curr => prev%next
    end do
    allocate(prev%next)
    curr => prev%next
    allocate(curr%var)
    hist_var => curr%var
    ! Populate meta-data associated with this variable
    hist_var%name      = trim(shortname)
    hist_var%long_name = trim(longname)
    hist_var%numdims   = numdims
    hist_var%dtype     = dtype
    ! Determine the dimension id's saved in the netCDF file and associated with
    ! this variable, check if variable has a time dimension
    has_t_dim = .false.
    allocate(hist_var%dimid(numdims),hist_var%dimlen(numdims))
    do dim_ii = 1,numdims
      ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(var_dimensions(dim_ii)),hist_var%dimid(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to find dimension id for "//trim(var_dimensions(dim_ii)),ierr)
      ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,hist_var%dimid(dim_ii),hist_var%dimlen(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension "//trim(var_dimensions(dim_ii)),ierr)
      if (hist_var%dimlen(dim_ii).eq.0) has_t_dim = .true.
    end do

    ! Distribute responsibility for writing cores over all PIO ranks
    ! i.e. compute the degrees of freedom that this rank will contribute to PIO
    ! TODO: Possible todo would be to allow for a subset of cores to be assigned
    ! PIO.
    if (has_t_dim) then
      call get_compdof(numdims-1,hist_var%dimlen(:numdims-1),my_dof_len,istart,istop)
    else
      call get_compdof(numdims,hist_var%dimlen(:numdims),my_dof_len,istart,istop)
    end if
    allocate( hist_var%compdof(my_dof_len) )
    hist_var%compdof(:my_dof_len) = (/ (ii, ii=istart,istop, 1) /)

    ! Register Variable with PIO
    ierr = PIO_def_var(pio_atm_file%pioFileDesc, trim(shortname), hist_var%dtype, hist_var%dimid(:numdims), hist_var%pioVar)
    call errorHandle("PIO ERROR: could not define variable "//trim(shortname),ierr)

    ! Gather the pio decomposition for this variable, and assign the pointer.
    if (has_t_dim) then
      loc_len = max(1,numdims-1)
      call get_decomp(pio_decomp_tag,hist_var%dtype,hist_var%dimlen(:loc_len),hist_var%compdof,hist_var%iodesc)
    else
      call get_decomp(pio_decomp_tag,hist_var%dtype,hist_var%dimlen,hist_var%compdof,hist_var%iodesc)
    end if

    return
  end subroutine register_variable
!=====================================================================!
  subroutine eam_history_write()
    ! WARNING: THIS SUBROUTINE IS A PLACEHOLDER TO ROUTINE TO TEST OUTPUT
    ! All of the following data itself is for testing purposes only.  This
    ! routine will either be changed, or more likely removed when the C++
    ! interface has been developed.  It does serve as a guide to how data is
    ! written to the output file.
    type(pio_atm_output),pointer :: pio_atm_file
    type(hist_var_t), pointer :: var 
    integer :: ivar, ierr, jj
    integer :: my_dummy_value, mylen
    integer, allocatable :: fdata_int(:)
    real(rtype) :: my_dummy_r
    real(rtype), allocatable :: fdata_real(:)
    real(rtype) :: time
    character(len=100) :: varname
    integer :: fdata_int_1d(10), fdata_int_2d(10,3)
    real(rtype) :: fdata_real_1d(10), fdata_real_2d(10,3), fdata_real_2d_inv(3,10)


    call get_pio_atm_file("example_pio_structured.nc",pio_atm_file)
    pio_atm_file%numRecs = pio_atm_file%numRecs + 1
    time = pio_atm_file%numRecs * 300.0_rtype

    ! DUMMY Array output for testing
    my_dummy_value = 100*pio_myrank + pio_atm_file%numRecs
    my_dummy_r     = 100.01_rtype * pio_myrank + real(pio_atm_file%numRecs)
    fdata_int_1d(:)   = my_dummy_value
    fdata_int_2d(:,1) = my_dummy_value
    fdata_int_2d(:,2) = my_dummy_value! + 1000
    fdata_real_1d(:)   = my_dummy_r
    fdata_real_2d(:,1) = my_dummy_r
    fdata_real_2d(:,2) = my_dummy_r! + 1000
    ! Record TIME
    call get_var(pio_atm_file,'time',var)
    ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,(/ pio_atm_file%numRecs /), (/ 1 /), (/ time /))
    ! Record the outputs
    call grid_write_data_array("example_pio_structured.nc", fdata_int_1d,  "foo")
    call grid_write_data_array("example_pio_structured.nc", fdata_real_1d, "real_foo")
    call grid_write_data_array("example_pio_structured.nc", fdata_real_1d, "x")
    call grid_write_data_array("example_pio_structured.nc", fdata_real_1d, "y")
    call grid_write_data_array("example_pio_structured.nc", fdata_real_1d, "bar2")
    call grid_write_data_array("example_pio_structured.nc", fdata_real_2d, "bar")

    call PIO_syncfile(pio_atm_file%pioFileDesc)

    call get_pio_atm_file("example_pio_structured_v2.nc",pio_atm_file)
    pio_atm_file%numRecs = pio_atm_file%numRecs + 1
    time = pio_atm_file%numRecs * 300.0_rtype
    my_dummy_r = 100.0_rtype*(pio_myrank+1) + pio_atm_file%numRecs
    do ivar = 1,3
      do jj = 1,10
      fdata_real_1d(jj) = jj
      fdata_real_2d(jj,ivar) = 10*ivar + jj!my_dummy_r + 1000*ivar
      fdata_real_2d_inv(ivar,jj) = 10*ivar + jj!my_dummy_r + 1000*ivar
      end do
    end do
    ! Record TIME
    call get_var(pio_atm_file,'time',var)
    ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,(/ pio_atm_file%numRecs /), (/ 1 /), (/ time /))
    ! Record outputs
    if (pio_atm_file%numRecs == 1) then
      call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_1d, "x")
      call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_1d, "y")
    end if
    call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_2d_inv, "bar")
    call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_2d, "foo")
    call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_2d_inv, "bar_flip")
    call grid_write_data_array("example_pio_structured_v2.nc", fdata_real_2d, "foo_flip")
    call PIO_syncfile(pio_atm_file%pioFileDesc)
    
  end subroutine eam_history_write
!=====================================================================!
  subroutine eam_history_finalize()

  ! For now nothing to do, could be used in the future.  If not, should be
  ! deleted.

  end subroutine eam_history_finalize
!=====================================================================!
  ! Assign header metadata to a specific pio output file.
  subroutine eam_pio_createHeader(File)

    type(file_desc_t), intent(in) :: File             ! Pio file Handle
    integer                       :: retval

    ! TODO change options below to match specific simulation case
    retval=pio_put_att (File, PIO_GLOBAL, 'source', 'SCREAM')
    retval=pio_put_att (File, PIO_GLOBAL, 'case', 'TEST 1')
    retval=pio_put_att (File, PIO_GLOBAL, 'title', 'SCORPIO TEST')
    retval=pio_put_att (File, PIO_GLOBAL, 'logname','THE GIT LOG HASH')
    retval=pio_put_att (File, PIO_GLOBAL, 'host', 'THE HOST')
    retval=pio_put_att (File, PIO_GLOBAL, 'Version', &
           '0')
    retval=pio_put_att (File, PIO_GLOBAL, 'revision_Id', &
           'None')
    retval=pio_put_att (File, PIO_GLOBAL, 'initial_file', 'NONE FOR NOW')
    retval=pio_put_att (File, PIO_GLOBAL, 'topography_file', 'NONE FOR NOW')
    
  end subroutine eam_pio_createHeader
!=====================================================================!
  ! Query the pio subsystem, pio rank and number of pio ranks from the component
  ! coupler.
  subroutine eam_init_pio_subsystem(mpicom,atm_id)
    use shr_pio_mod,   only: shr_pio_getiosys, shr_pio_getiotype
    
    integer, intent(in) :: mpicom
    integer, intent(in) :: atm_id

    integer :: ierr

    pio_mpicom = mpicom

    call MPI_Comm_rank(pio_mpicom, pio_myrank, ierr)
    call MPI_Comm_size(pio_mpicom, pio_ntasks , ierr)
    
    pio_subsystem  => shr_pio_getiosys(atm_id)
    pio_iotype     = shr_pio_getiotype(atm_id)
    pio_rearranger = shr_pio_getrearranger(atm_id)

  end subroutine eam_init_pio_subsystem
!=====================================================================!
  ! Create a file with the appropriate name
  subroutine eam_pio_createfile(File,fname)

    type(file_desc_t), intent(inout) :: File             ! Pio file Handle
    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    integer                          :: retval           ! PIO error return value
    integer                          :: mode             ! Mode for how to handle the new file

    mode = pio_clobber ! Set to CLOBBER for now, TODO: fix to allow for optional mode type like in CAM
    retval = pio_createfile(pio_subsystem,File,pio_iotype,fname,mode) 

  end subroutine eam_pio_createfile
!=====================================================================!
  subroutine eam_pio_finalize()
    ! May not be needed, possibly handled by PIO directly.

    integer :: ierr

    call PIO_finalize(pio_subsystem, ierr)

  end subroutine eam_pio_finalize
!=====================================================================!
  ! Handle any errors that crop up
  subroutine errorHandle(errMsg, retVal)
    use iso_c_binding
    implicit none

    interface
      subroutine finalize_scream_session() bind(C)
        ! No inputs or anything, just interface to C code.
      end subroutine finalize_scream_session
    end interface

    character(len=*),  intent(in)    :: errMsg
    integer,           intent(in)    :: retVal

    type(pio_file_list), pointer :: curr => NULL()

    if (retVal .ne. PIO_NOERR) then
      write(*,*) retVal,errMsg
      ! Close all the PIO Files before aborting run
      curr => pio_file_list_top
      do while (associated(curr))
        if (associated(curr%pio_file)) call PIO_closefile(curr%pio_file%pioFileDesc)
        curr => curr%next
      end do
      ! Kill run
      call finalize_scream_session()
      call mpi_abort(pio_mpicom,0,retVal)
    end if

  end subroutine errorHandle
!=====================================================================!
  ! Algorithm to determine the degrees-of-freedom in the global array that this
  ! PIO rank is responsible for writing.  Needed in the pio_write interface.
  subroutine get_compdof(numdims,dimension_len,dof_len,istart,istop)

    integer, intent(in)  :: numdims                ! Number of dimensions
    integer, intent(in)  :: dimension_len(numdims) ! Array of each dimension length
    integer, intent(out) :: dof_len, istart, istop ! Length of degrees-of-freedom for this rank (dof), start and stop in global array (flattened to 1d)

    integer :: extra_procs, total_dimlen

    ! Get the total number of array elements for this output
    total_dimlen = product(dimension_len)
    dof_len   = total_dimlen/pio_ntasks
    ! If the number of pio tasks does not evenly divide the total number of
    ! array elements we need to assign less degrees of freedom to the final PIO
    ! task.
    extra_procs = mod(total_dimlen,pio_ntasks)
    if (extra_procs > 0) dof_len = dof_len + 1
    ! Determine the starting and finishing array location for the output chunk
    ! handled by this PIO task
    istart = pio_myrank * dof_len + 1
    istop  = istart +  dof_len - 1
    ! Special treatment for the final PIO task, which may have less dof's to
    ! write
    if (pio_myrank == pio_ntasks-1) then
      istop = total_dimlen
      dof_len = istop-istart+1
    end if
    return 
 
  end subroutine get_compdof
!=====================================================================!
 ! Determine the unique pio_decomposition for this output grid, if it hasn't
 ! been defined create a new one.
  subroutine get_decomp(tag,dtype,dimension_len,compdof,iodesc)
    ! TODO: CAM code creates the decomp tag for the user.  Theoretically it is
    ! unique because it is based on dimensions and datatype.  But the tag ends
    ! up not being very descriptive.  The todo item is to revisit how tags are
    ! handled and decide if we want the code to create a tag or let the use
    ! assign a tag.
    character(len=*)          :: tag              ! Unique tag string describing this output grid
    integer, intent(in)       :: dtype            ! Datatype associated with the output
    integer, intent(in)       :: dimension_len(:) ! Array of the dimension lengths for this decomp
    integer, intent(in)       :: compdof(:)       ! The degrees of freedom this rank is responsible for
    type(io_desc_t), pointer  :: iodesc           ! The pio decomposition that has been found or created

    logical                   :: found            ! Keep track if a decomp has been found among the previously defined decompositions 
    type(iodesc_list),pointer :: curr, prev       ! Used to toggle through the recursive list of decompositions
    type(io_desc_t), pointer  :: iodesc_loc       ! Local iodesc structure
    integer                   :: loc_len          ! Used to keep track of how many dimensions there are in decomp 
    
    ! Assign a PIO decomposition to variable, if none exists, create a new one:
    found = .false.
    curr => iodesc_list_top
    ! Cycle through all current iodesc to see if the decomp has already been
    ! created
    do while(associated(curr) .and. (.not.found))
      if (trim(tag) == trim(curr%tag)) then
        found = .true.
        iodesc => curr%iodesc
      else
        prev => curr
        curr => curr%next
      end if
    end do
    ! If we didn't find an iodesc then we need to create one
    if (.not.found) then
      curr => prev ! Go back and allocate the new iodesc in curr%next
      if(associated(curr%iodesc)) then
        allocate(curr%next)
        curr => curr%next
        nullify(curr%iodesc)  ! Extra step to ensure clean iodesc
        nullify(curr%next)  ! Extra step to ensure clean iodesc
      end if
      allocate(curr%iodesc)
      curr%tag = trim(tag)
      loc_len = size(dimension_len)
      if ( loc_len.eq.1 .and. dimension_len(loc_len).eq.0 ) then
        allocate(curr%iodesc)
      else
        call pio_initdecomp(pio_subsystem, dtype, dimension_len, compdof, curr%iodesc, rearr=pio_rearranger)
      end if
      iodesc => curr%iodesc
    end if 

  end subroutine get_decomp
!=====================================================================!
  ! Query the hist_var_t pointer for a specific variable on a specific file.
  subroutine get_var(pio_file,varname,var)

    type(pio_atm_output), pointer :: pio_file ! Pio output file structure
    character(len=*)              :: varname  ! Name of the variable to query
    type(hist_var_t), pointer     :: var      ! Pointer to the variable structure that has been found
    
    type(hist_var_list), pointer  :: curr     ! Used to cycle through recursive list of variables

    curr => pio_file%var_list_top
    do while (associated(curr))
      var => curr%var
      if (associated(var)) then
        if (trim(varname) == trim(var%name)) return
      end if
      curr => curr%next
    end do

    ! If we got this far we didn't find the variable
    call errorHandle("PIO ERROR: unable to find variable: "//trim(varname)//" in file: "//trim(pio_file%filename),999)

  end subroutine get_var
!=====================================================================!
  ! Lookup pointer for pio file based on filename.
  subroutine get_pio_atm_file(filename,pio_atm_file,found)

    character(len=*),intent(in)   :: filename     ! Name of file to be found
    type(pio_atm_output), pointer :: pio_atm_file ! Pointer to pio_atm_output structure associated with this filename
    logical,optional,intent(out)  :: found        ! return TRUE if pio file already exists

    type(pio_file_list), pointer  :: curr => NULL(), prev => NULL() ! Used to cycle through recursive list of pio atm files

    if (present(found)) found = .false.

    ! Starting at the top of the current list of PIO_FILES find the next empty
    ! one in list.
    curr => pio_file_list_top
    do while (associated(curr))
      if (associated(curr%pio_file)) then
        if (trim(filename)==trim(curr%pio_file%filename)) then
          pio_atm_file => curr%pio_file
          if (present(found)) found = .true.
          return
        end if
      end if
      prev => curr
      curr => prev%next
    end do
    allocate(prev%next)
    curr => prev%next
    allocate(curr%pio_file)
    ! Create and initialize the new pio file:
    curr%pio_file%filename = trim(filename)
    curr%pio_file%numRecs = 0
    call eam_pio_createfile(curr%pio_file%pioFileDesc,trim(curr%pio_file%filename))
    pio_atm_file => curr%pio_file  

  end subroutine get_pio_atm_file
!=====================================================================!
  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(in)    :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_1d_int: Error writing variable',ierr)
  end subroutine grid_write_darray_1d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(in)    :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    integer, allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_2d_int: Error writing variable',ierr)
  end subroutine grid_write_darray_2d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(in)    :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    integer, allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_3d_int: Error writing variable',ierr)
  end subroutine grid_write_darray_3d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_double(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_1d_double: Error writing variable',ierr)
  end subroutine grid_write_darray_1d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_double(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    real(rtype), allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_2d_double: Error writing variable',ierr)
  end subroutine grid_write_darray_2d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_double(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    real(rtype), allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_3d_double: Error writing variable',ierr)

  end subroutine grid_write_darray_3d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype_short),         intent(in)    :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_1d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_1d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype_short),         intent(in)    :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    real(rtype_short), allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_2d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_2d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype_short),         intent(in)    :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_output),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    real(rtype_short), allocatable, dimension(:)   :: vbuf                 

    call get_pio_atm_file(trim(filename),pio_atm_file)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    vbuf = pack(hbuf,.true.)
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, vbuf(var%compdof), ierr)
    call errorHandle( 'cam_grid_write_darray_3d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_3d_real

end module scream_scorpio_interface
