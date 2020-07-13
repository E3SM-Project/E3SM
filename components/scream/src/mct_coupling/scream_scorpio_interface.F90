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
!    This is accomplished during eam_pio_enddef by calling 'PIO_enddef' on all
!    defined files.
!
! Writing Output:
! 
! Finalization: 
!==============================================================================!

  ! TODO: have the code rely on shr_pio_mod only when the build is of the full
  ! model. Note, these three variables from shr_pio_mod are only needed in the
  ! case when we want to use the pio_subsystem that has already been defined by
  ! the component coupler. 
  use shr_pio_mod,  only: shr_pio_getrearranger, shr_pio_getiosys, shr_pio_getiotype
  !------------
  use physics_utils, only: rtype, rtype8, itype, btype
  use piolib_mod, only : PIO_init, PIO_finalize, PIO_createfile, PIO_closefile, &
      PIO_initdecomp, PIO_freedecomp, PIO_syncfile, PIO_openfile, PIO_setframe, &
      pio_init
  use pio_types,  only : iosystem_desc_t, file_desc_t, &
      pio_noerr, PIO_iotype_netcdf, var_desc_t, io_desc_t, PIO_int, &
      pio_clobber, PIO_nowrite, PIO_unlimited, pio_global, PIO_real, &
      PIO_double, pio_rearr_subset
  use pio_kinds,  only : PIO_OFFSET_KIND, i4 
  use pio_nf,     only : PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef, PIO_inq_dimid, &
                         PIO_inq_dimlen, PIO_inq_varid
  use piodarray,  only : PIO_write_darray, PIO_read_darray 
  use pionfatt_mod, only : PIO_put_att   => put_att
  use pionfput_mod, only : PIO_put_var   => put_var

  implicit none
  save

#include "scream_config.f"
                     
  public :: & 
            eam_pio_closefile,      & ! Close a specfic pio file.
            eam_pio_enddef,         & ! Register variables and dimensions with PIO files
            eam_init_pio_subsystem, & ! Gather pio specific data from the component coupler
            eam_pio_finalize,       & ! Run any final PIO commands
            register_outfile,       & ! Create a pio output file
            register_infile,        & ! Open a pio input file
            register_variable,      & ! Register a variable with a particular pio output file  
            get_variable,           & ! Register a variable with a particular pio output file  
            register_dimension,     & ! Register a dimension with a particular pio output file
            set_decomp,             & ! Set the pio decomposition for all variables in file.
            set_dof,                & ! Set the pio dof decomposition for specific variable in file.
            grid_write_data_array,  & ! Write gridded data to a pio managed netCDF file
            grid_read_data_array,   & ! Read gridded data from a pio managed netCDF file
            eam_sync_piofile,       & ! Syncronize the piofile, to be done after all output is written during a single timestep
            eam_update_time           ! Update the timestamp (i.e. time variable) for a given pio netCDF file
 
  private :: errorHandle
  ! Universal PIO variables for the module
  integer               :: pio_mpicom  
  integer               :: pio_iotype
  type(iosystem_desc_t), pointer, public :: pio_subsystem
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
    character(len=max_chars) :: pio_decomp_tag ! PIO decomposition label used by this variable.
    character(len=max_chars) :: units         ! 'units' attribute
    type(var_desc_t) :: piovar                ! netCDF variable ID
    integer          :: dtype                 ! data type
    integer          :: numdims               ! Number of dimensions in out field
    type(io_desc_t), pointer  :: iodesc       ! PIO decomp associated with this variable
    integer, allocatable :: compdof(:)        ! Global locations in output array for this process
    integer, allocatable :: dimid(:)          ! array of PIO dimension id's for this variable
    integer, allocatable :: dimlen(:)         ! array of PIO dimension lengths for this variable
    logical              :: has_t_dim         ! true, if variable has a time dimension
    logical              :: compdof_set = .false. ! true after the dof for this rank has been set.
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
    type(pio_atm_file_t), pointer :: pio_file => NULL() ! Pointer to an atm. pio file
    type(pio_file_list),  pointer :: next => NULL()     ! Needed for recursive definition
  end type pio_file_list
  ! Define the first pio_file_list
  type(pio_file_list), target  :: pio_file_list_top
  type(pio_file_list), pointer :: pio_file_list_bottom
!----------------------------------------------------------------------
  type, public :: pio_atm_file_t
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

        !> @brief Whether or not this pio file is still open
        logical                         :: isopen = .false.

        !> @brief Whether or not the dim/var definition phase is still open
        logical                         :: is_enddef = .false.

  end type pio_atm_file_t

!----------------------------------------------------------------------
  interface grid_read_data_array
    module procedure grid_read_darray_1d_real
    module procedure grid_read_darray_2d_real
    module procedure grid_read_darray_3d_real
    module procedure grid_read_darray_4d_real
    module procedure grid_read_darray_1d_int
    module procedure grid_read_darray_2d_int
    module procedure grid_read_darray_3d_int
    module procedure grid_read_darray_4d_int
  end interface grid_read_data_array
!----------------------------------------------------------------------
  interface grid_write_data_array
    module procedure grid_write_darray_1d_int
    module procedure grid_write_darray_2d_int
    module procedure grid_write_darray_3d_int
    module procedure grid_write_darray_4d_int
    module procedure grid_write_darray_1d_real
    module procedure grid_write_darray_2d_real
    module procedure grid_write_darray_3d_real
    module procedure grid_write_darray_4d_real
  end interface
!----------------------------------------------------------------------
contains
!=====================================================================!
  subroutine register_outfile(filename)

    character(len=*), intent(in) :: filename

    type(pio_atm_file_t), pointer :: current_atm_file => null()

    call get_new_pio_atm_file(filename,current_atm_file,1)
    call eam_pio_createHeader(current_atm_file%pioFileDesc)
    
  end subroutine register_outfile
!=====================================================================!
  subroutine register_infile(filename)

    character(len=*), intent(in) :: filename

    type(pio_atm_file_t), pointer :: current_atm_file => null()

    call get_new_pio_atm_file(filename,current_atm_file,2)
    
  end subroutine register_infile
!=====================================================================!
  subroutine eam_pio_enddef(filename)

    character(len=*), intent(in) :: filename

    type(pio_atm_file_t), pointer :: current_atm_file => null()
    integer                       :: ierr
    logical                       :: found

    call lookup_pio_atm_file(filename,current_atm_file,found)

    ! Gather the pio decomposition for all variables in this file, and assign them pointers.
    call set_decomp(trim(filename))
    ! Officially close the definition step for this file.
    ierr = PIO_enddef(current_atm_file%pioFileDesc)
    call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),ierr)
    current_atm_file%is_enddef = .true.

  end subroutine eam_pio_enddef
!=====================================================================!
  ! Register a dimension with a specific pio output file 
  subroutine register_dimension(pio_atm_filename,shortname,longname,length)
    character(len=*), intent(in)        :: pio_atm_filename   ! Name of file to register the dimension on.
    character(len=*), intent(in)        :: shortname,longname ! Short- and long- names for this dimension, short: brief identifier and name for netCDF output, long: longer descriptor sentence to be included as meta-data in file.
    integer, intent(in)                 :: length             ! Length of the dimension, 0: unlimited (like time), >0 actual length of dimension

    type(pio_atm_file_t), pointer       :: pio_atm_file
    type(hist_coord_t), pointer         :: hist_coord
    type(hist_coord_list), pointer      :: curr=>null(), prev=>null()
    integer                             :: ierr
    logical                             :: found

    ! Make sure the dimension length is reasonable
    if (length<0) call errorHandle("PIO Error: dimension "//trim(shortname)//", can't have a negative dimension length",-999)
 
    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(pio_atm_filename),pio_atm_file,found)
    if (.not.found) call errorHandle("PIO Error: can't find pio_atm_file associated with file: "//trim(pio_atm_filename),-999) 
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
  subroutine get_variable(pio_atm_filename,shortname,longname,numdims,var_dimensions,dtype,pio_decomp_tag)
    character(len=*), intent(in) :: pio_atm_filename         ! Name of the file to register this variable with
    character(len=*), intent(in) :: shortname,longname       ! short and long names for the variable.  Short: variable name in file, Long: more descriptive name
    integer, intent(in)          :: numdims                  ! Number of dimensions for this variable, including time dimension
    character(len=*), intent(in) :: var_dimensions(numdims)  ! String array with shortname descriptors for each dimension of variable.
    integer, intent(in)          :: dtype                    ! datatype for this variable, REAL, DOUBLE, INTEGER, etc.
    character(len=*), intent(in) :: pio_decomp_tag           ! Unique tag for this variables decomposition type, to be used to determine if the io-decomp already exists.

    ! Local variables
    type(pio_atm_file_t),pointer :: pio_atm_file
    integer                      :: loc_len
    type(hist_var_t), pointer    :: hist_var
    integer                      :: dim_ii
    integer                      :: ierr
    integer, allocatable         :: dimlen(:)
    integer                      :: my_dof_len
    integer, allocatable         :: compdof(:)
    integer                      :: ii, istart, istop
    logical                      :: found

    type(hist_var_list), pointer :: curr => null(), prev => null()
 
    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(pio_atm_filename),pio_atm_file,found)
    ! Update the number of variables on file 
    pio_atm_file%varcounter = pio_atm_file%varcounter + 1

    ! Get a new variable pointer in var_list
    curr => pio_atm_file%var_list_top
    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(shortname)) call errorHandle("PIO Error: Could not register variable "//trim(shortname)//", already exists in file: "//trim(pio_atm_filename),-999)
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
    hist_var%pio_decomp_tag = trim(pio_decomp_tag) 
    ! Determine the dimension id's saved in the netCDF file and associated with
    ! this variable, check if variable has a time dimension
    hist_var%has_t_dim = .false.
    allocate(hist_var%dimid(numdims),hist_var%dimlen(numdims))
    do dim_ii = 1,numdims
      ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(var_dimensions(dim_ii)),hist_var%dimid(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to find dimension id for "//trim(var_dimensions(dim_ii)),ierr)
      ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,hist_var%dimid(dim_ii),hist_var%dimlen(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension "//trim(var_dimensions(dim_ii)),ierr)
      if (hist_var%dimlen(dim_ii).eq.0) hist_var%has_t_dim = .true.
    end do

    ! Register Variable with PIO
    ! check to see if variable already is defined with file (for use with input)
    ierr = PIO_inq_varid(pio_atm_file%pioFileDesc,trim(shortname),hist_var%piovar)
    call errorHandle("PIO ERROR: could not find variable "//trim(shortname)//" in file "//trim(pio_atm_filename),ierr)
    
    return
  end subroutine get_variable
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
    type(pio_atm_file_t),pointer :: pio_atm_file
    integer                      :: loc_len
    type(hist_var_t), pointer    :: hist_var
    integer                      :: dim_ii
    integer                      :: ierr
    integer, allocatable         :: dimlen(:)
    integer                      :: my_dof_len
    integer, allocatable         :: compdof(:)
    integer                      :: ii, istart, istop
    logical                      :: found

    type(hist_var_list), pointer :: curr => null(), prev => null()
 
    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(pio_atm_filename),pio_atm_file,found)
    ! Update the number of variables on file 
    pio_atm_file%varcounter = pio_atm_file%varcounter + 1

    ! Get a new variable pointer in var_list
    curr => pio_atm_file%var_list_top
    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(shortname)) call errorHandle("PIO Error: Could not register variable "//trim(shortname)//", already exists in file: "//trim(pio_atm_filename),-999)
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
    hist_var%pio_decomp_tag = trim(pio_decomp_tag) 
    ! Determine the dimension id's saved in the netCDF file and associated with
    ! this variable, check if variable has a time dimension
    hist_var%has_t_dim = .false.
    allocate(hist_var%dimid(numdims),hist_var%dimlen(numdims))
    do dim_ii = 1,numdims
      ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(var_dimensions(dim_ii)),hist_var%dimid(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to find dimension id for "//trim(var_dimensions(dim_ii)),ierr)
      ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,hist_var%dimid(dim_ii),hist_var%dimlen(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension "//trim(var_dimensions(dim_ii)),ierr)
      if (hist_var%dimlen(dim_ii).eq.0) hist_var%has_t_dim = .true.
    end do

    ! Register Variable with PIO
    ! First, check to see if variable already is defined with file
    ierr = PIO_inq_varid(pio_atm_file%pioFileDesc,trim(shortname),hist_var%piovar)
    if (ierr == PIO_NOERR) call errorHandle("PIO ERROR: could not define variable "//trim(shortname)//" in file "//trim(pio_atm_filename)//", already exists",-999)
    
    ! if ierr is not pio_noerror then the variable needs to be defined
    if (ierr.ne.pio_noerr) ierr = PIO_def_var(pio_atm_file%pioFileDesc, trim(shortname), hist_var%dtype, hist_var%dimid(:numdims), hist_var%piovar)
    call errorHandle("PIO ERROR: could not define variable "//trim(shortname),ierr)

    return
  end subroutine register_variable
!=====================================================================!
  subroutine eam_update_time(filename,time)
    character(len=*), intent(in) :: filename       ! PIO filename
    real(rtype), intent(in)      :: time 

    type(hist_var_t), pointer    :: var 
    type(pio_atm_file_t),pointer   :: pio_atm_file
    integer                      :: ierr
    logical                      :: found

    call lookup_pio_atm_file(filename,pio_atm_file,found)
    pio_atm_file%numRecs = pio_atm_file%numRecs + 1
    call get_var(pio_atm_file,'time',var)
    ! Only update time on the file if a valid time is provided
    if (time>=0) ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,(/ pio_atm_file%numRecs /), (/ 1 /), (/ time /))
  end subroutine eam_update_time
!=====================================================================!
  subroutine eam_sync_piofile(filename)
    character(len=*),          intent(in)    :: filename       ! PIO filename
    
    type(pio_atm_file_t),pointer             :: pio_atm_file
    logical                      :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call PIO_syncfile(pio_atm_file%pioFileDesc)
  end subroutine eam_sync_piofile
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
  subroutine eam_init_pio_subsystem(mpicom,atm_id,local)
    
    integer, intent(in) :: mpicom
    integer, intent(in) :: atm_id
    logical, intent(in) :: local

    integer :: ierr
    integer :: my_task
    integer :: nprocs
    integer :: stride
    integer :: num_aggregator
    integer :: base

    if (associated(pio_subsystem)) call errorHandle("PIO ERROR: local pio_subsystem pointer has already been established.",-999)

    pio_mpicom = mpicom
    call MPI_Comm_rank(pio_mpicom, pio_myrank, ierr)
    call MPI_Comm_size(pio_mpicom, pio_ntasks , ierr)
   
    if (.not.local) then 
      pio_subsystem  => shr_pio_getiosys(atm_id)
      pio_iotype     = shr_pio_getiotype(atm_id)
      pio_rearranger = shr_pio_getrearranger(atm_id)
    else
      allocate(pio_subsystem)
      stride         = 1
      num_aggregator = 0
      pio_rearranger = pio_rearr_subset
      pio_iotype     = PIO_iotype_netcdf
      base           = 0
      call PIO_init(pio_myrank, pio_mpicom, pio_ntasks, num_aggregator, stride, &
           pio_rearr_subset, pio_subsystem, base=base)
    end if

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
    call errorHandle("PIO ERROR: unable to create file: "//trim(fname),retval)

  end subroutine eam_pio_createfile
!=====================================================================!
  subroutine eam_pio_openfile(File,fname)

    type(file_desc_t), intent(inout) :: File             ! Pio file Handle
    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    integer                          :: retval           ! PIO error return value
    integer                          :: mode             ! Mode for how to handle the new file

    mode = 0 ! TODO: make sure this is correct. 
    retval = pio_openfile(pio_subsystem,File,pio_iotype,fname,mode)
    call errorHandle("PIO ERROR: unable to open file: "//trim(fname),retval)

  end subroutine eam_pio_openfile
!=====================================================================!
  subroutine eam_pio_closefile(fname)

    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    type(pio_atm_file_t),pointer     :: pio_atm_file => null()
    logical                          :: found

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(fname),pio_atm_file,found)
    if (found) then
      call PIO_closefile(pio_atm_file%pioFileDesc)
      pio_atm_file%isopen = .false.
    else
      call errorHandle("PIO ERROR: unable to close file: "//trim(fname)//", was not found",-999)
    end if

  end subroutine eam_pio_closefile
!=====================================================================!
  subroutine eam_pio_finalize()
    ! May not be needed, possibly handled by PIO directly.

    integer :: ierr
    type(pio_file_list), pointer :: curr => NULL()

    ! Close all the PIO Files 
    curr => pio_file_list_top
    do while (associated(curr))
      if (associated(curr%pio_file)) then
        if (curr%pio_file%isopen) call PIO_closefile(curr%pio_file%pioFileDesc)
        curr%pio_file%isopen = .false.
      end if
      curr => curr%next
    end do
    call PIO_finalize(pio_subsystem, ierr)
    nullify(pio_subsystem)

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
      write(*,'(I8,2x,A100)') retVal,trim(errMsg)
      ! Kill run
      call eam_pio_finalize() 
      call finalize_scream_session()
      call mpi_abort(pio_mpicom,0,retVal)
    end if

  end subroutine errorHandle
!=====================================================================!
  ! Algorithm to determine the degrees-of-freedom in the global array that this
  ! PIO rank is responsible for writing.  Needed in the pio_write interface.
  ! TODO: For unit test at least this isn't used.  Should we delete it and
  ! always expect the C++ code to establish the DOF locally, or keep this as a
  ! tool that can be used if needed?
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
  subroutine set_dof(filename,varname,dof_len,dof_vec)
    character(len=*), intent(in)            :: filename
    character(len=*), intent(in)            :: varname
    integer, intent(in)                     :: dof_len
    integer, intent(in), dimension(dof_len) :: dof_vec

    type(pio_atm_file_t),pointer            :: pio_atm_file
    type(hist_var_t), pointer               :: var
    logical                                 :: found
    integer                                 :: ii
    
    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    if (allocated(var%compdof)) deallocate(var%compdof)
    allocate( var%compdof(dof_len) )
    do ii = 1,dof_len
      var%compdof(ii) = dof_vec(ii)
    end do
    var%compdof_set = .true.
 
  end subroutine set_dof
!=====================================================================!
  subroutine set_decomp(filename)

    character(len=*)              :: filename  ! Name of the pio file to set decomp for

    type(pio_atm_file_t), pointer :: current_atm_file => null()
    type(hist_var_list), pointer  :: curr     ! Used to cycle through recursive list of variables
    type(hist_var_t), pointer     :: hist_var ! Pointer to the variable structure that has been found
    integer                       :: loc_len
    logical                       :: found

    call lookup_pio_atm_file(filename,current_atm_file,found)
    if (current_atm_file%is_enddef) call errorHandle("PIO ERROR: unable to set decomposition in file: "//trim(current_atm_file%filename)//", definition phase has ended (pio_enddef)",999) 
    curr => current_atm_file%var_list_top
    do while (associated(curr))
      if (associated(curr%var)) then
        hist_var => curr%var
        if (.not.hist_var%compdof_set) call errorHandle("PIO ERROR: unable to set decomp for file, var: "//trim(current_atm_file%filename)//", "//trim(hist_var%name)//". Set DOF.",999) 
        ! Assign decomp
        if (hist_var%has_t_dim) then
          loc_len = max(1,hist_var%numdims-1)
          call get_decomp(hist_var%pio_decomp_tag,hist_var%dtype,hist_var%dimlen(:loc_len),hist_var%compdof,hist_var%iodesc)
        else
          call get_decomp(hist_var%pio_decomp_tag,hist_var%dtype,hist_var%dimlen,hist_var%compdof,hist_var%iodesc)
        end if
      end if
      curr => curr%next
    end do

  end subroutine set_decomp
!=====================================================================!
  ! Query the hist_var_t pointer for a specific variable on a specific file.
  subroutine get_var(pio_file,varname,var)

    type(pio_atm_file_t), pointer :: pio_file ! Pio output file structure
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
  ! Diagnostic routine to determine how many pio files are currently open:
  subroutine count_pio_atm_file(total_count)
    integer, intent(out) :: total_count

    type(pio_file_list), pointer :: curr => NULL(), prev => NULL() ! Used to cycle through recursive list of pio atm files

    total_count = 0
    curr => pio_file_list_top
    do while (associated(curr))
      if (associated(curr%pio_file)) then
        if (curr%pio_file%isopen) total_count = total_count+1
      end if
      prev => curr
      curr => prev%next
    end do
  end subroutine count_pio_atm_file
!=====================================================================!
  ! Lookup pointer for pio file based on filename.
  subroutine lookup_pio_atm_file(filename,pio_file,found)

    character(len=*),intent(in)   :: filename     ! Name of file to be found
    type(pio_atm_file_t), pointer :: pio_file     ! Pointer to pio_atm_output structure associated with this filename
    logical, intent(out)          :: found        ! whether or not the file was found
!    type(pio_file_list), optional, pointer :: curr_bottom

    type(pio_file_list), pointer :: curr => NULL(), prev => NULL() ! Used to cycle through recursive list of pio atm files
    integer :: cnt
    ! Starting at the top of the current list of PIO_FILES search for this
    ! filename.
    cnt = 0
    found = .false.
    curr => pio_file_list_top
    do while (associated(curr))
      cnt = cnt+1
      if (associated(curr%pio_file)) then
        if (trim(filename)==trim(curr%pio_file%filename).and.curr%pio_file%isopen) then
          pio_file => curr%pio_file
          found = .true.
          return
        end if
      end if
      prev => curr
      curr => prev%next
    end do
    allocate(prev%next)
    pio_file_list_bottom => prev%next

  end subroutine lookup_pio_atm_file
!=====================================================================!
  ! Lookup pointer for pio file based on filename.
  subroutine get_new_pio_atm_file(filename,pio_file,purpose)

    character(len=*),intent(in)   :: filename     ! Name of file to be found
    type(pio_atm_file_t), pointer :: pio_file     ! Pointer to pio_atm_output structure associated with this filename
    integer,intent(in)            :: purpose      ! Purpose for this file lookup, 0 = find already existing, 1 = create new as output, 2 = open new as input

    logical                      :: found
    type(pio_file_list), pointer :: curr => NULL()

    ! Make sure a there isn't a pio_atm_file pointer already estalished for a
    ! file with this filename.
    call lookup_pio_atm_file(trim(filename),pio_file,found)
    if (found) call errorHandle("PIO Error: get_pio_atm_file with filename = "//trim(filename)//", has already been registered with the pio_file list.",-999)
    curr => pio_file_list_bottom 
    allocate(curr%pio_file)
    pio_file => curr%pio_file
    pio_file_list_bottom => curr%next
    ! Create and initialize the new pio file:
    pio_file%filename = trim(filename)
    pio_file%isopen = .true.
    if (purpose == 1) then  ! Will be used for output.  Set numrecs to zero and create the new file.
      pio_file%numRecs = 0
      call eam_pio_createfile(pio_file%pioFileDesc,trim(pio_file%filename))
    elseif (purpose == 2) then ! Will be used for input, just open it
      call eam_pio_openfile(pio_file%pioFileDesc,trim(pio_file%filename))
    else
      call errorHandle("PIO Error: get_pio_atm_file with filename = "//trim(filename)//", purpose (int) assigned to this lookup is not valid" ,-999)
    end if

  end subroutine get_new_pio_atm_file
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
    type(pio_atm_file_t),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_1d_int: Error writing variable',ierr)
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
    type(pio_atm_file_t), pointer              :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_2d_int: Error writing variable',ierr)
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
    type(pio_atm_file_t), pointer              :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_3d_int: Error writing variable',ierr)
  end subroutine grid_write_darray_3d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_4d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_4d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(in)    :: hbuf(:,:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_4d_int: Error writing variable',ierr)
  end subroutine grid_write_darray_4d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t), pointer              :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_1d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_1d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t), pointer              :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_2d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_2d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(in)    :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t), pointer              :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_3d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_3d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_4d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_4d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*), intent(in)           :: filename       ! PIO filename
    real(rtype), intent(in)                :: hbuf(:,:,:,:)
    character(len=*), intent(in)           :: varname

    ! Local variables
    type(pio_atm_file_t),pointer             :: pio_atm_file
    type(hist_var_t), pointer              :: var
    integer                                :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_write_darray_4d_real: Error writing variable',ierr)
  end subroutine grid_write_darray_4d_real
!=====================================================================!
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_1d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_1d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(out)   :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_1d_real: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_1d_real
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_2d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_2d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(out)   :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_2d_real: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_2d_real
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_3d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_3d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(out)   :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_3d_real: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_3d_real
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_4d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_4d_real(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    real(rtype),               intent(out)   :: hbuf(:,:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_4d_real: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_4d_real
!=====================================================================!
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_1d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_1d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(out)   :: hbuf(:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer             :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_1d_int: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_1d_int
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_2d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_2d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(out)   :: hbuf(:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_2d_int: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_2d_int
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_3d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_3d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(out)   :: hbuf(:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_3d_int: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_3d_int
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_4d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_4d_int(filename, hbuf, varname)

    ! Dummy arguments
    character(len=*),          intent(in)    :: filename       ! PIO filename
    integer,                   intent(out)   :: hbuf(:,:,:,:)
    character(len=*),          intent(in)    :: varname

    ! Local variables
    type(pio_atm_file_t),pointer               :: pio_atm_file
    type(hist_var_t), pointer                :: var
    integer                                  :: ierr
    logical                                  :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, hbuf, ierr)
    call errorHandle( 'eam_grid_read_darray_4d_int: Error reading variable '//trim(varname),ierr)

  end subroutine grid_read_darray_4d_int
!=====================================================================!

end module scream_scorpio_interface
