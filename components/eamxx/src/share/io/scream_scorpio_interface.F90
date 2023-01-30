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

  !------------
  use pio_types,    only: iosystem_desc_t, file_desc_t, var_desc_t, io_desc_t, &
                          pio_noerr, pio_global, &
                          PIO_int, PIO_real, PIO_double, PIO_float=>PIO_real
  use pio_kinds,    only: PIO_OFFSET_KIND
  use pio_nf,       only: PIO_enddef, PIO_inq_dimid, PIO_inq_dimlen, PIO_inq_varid, &
                          PIO_inquire, PIO_inquire_variable
  use pionfatt_mod, only: PIO_put_att   => put_att

  use mpi, only: mpi_abort, mpi_comm_size, mpi_comm_rank

  use iso_c_binding, only: c_float, c_double, c_int
  implicit none
  save

  public :: &
            lookup_pio_atm_file,         & ! Checks if a pio file is present
            eam_pio_closefile,           & ! Close a specfic pio file.
            eam_pio_enddef,              & ! Register variables and dimensions with PIO files
            eam_init_pio_subsystem,      & ! Gather pio specific data from the component coupler
            is_eam_pio_subsystem_inited, & ! Query whether the pio subsystem is inited already
            eam_pio_finalize,            & ! Run any final PIO commands
            register_file,               & ! Creates/opens a pio input/output file
            register_variable,           & ! Register a variable with a particular pio output file
            set_variable_metadata,       & ! Sets a variable metadata (always char data)
            get_variable,                & ! Register a variable with a particular pio output file
            register_dimension,          & ! Register a dimension with a particular pio output file
            set_decomp,                  & ! Set the pio decomposition for all variables in file.
            set_dof,                     & ! Set the pio dof decomposition for specific variable in file.
            grid_write_data_array,       & ! Write gridded data to a pio managed netCDF file
            grid_read_data_array,        & ! Read gridded data from a pio managed netCDF file
            eam_update_time,             & ! Update the timestamp (i.e. time variable) for a given pio netCDF file
            get_int_attribute,           & ! Retrieves an integer global attribute from the nc file
            set_int_attribute,           & ! Writes an integer global attribute to the nc file
            get_dimlen,                  & ! Returns the length of a specific dimension in a file
            read_time_at_index,          & ! Returns the time stamp for a specific time index
            has_variable                   ! Checks if given file contains a certain variable

  private :: errorHandle, get_coord

  ! Universal PIO variables for the module
  integer               :: atm_mpicom
  integer               :: pio_iotype
  type(iosystem_desc_t), pointer, public :: pio_subsystem
  integer               :: pio_rearranger
  integer               :: pio_mode
  integer               :: time_dimid = -1

  ! TYPES to handle history coordinates and files
  integer,parameter :: max_hcoordname_len = 16
  integer,parameter :: max_chars = 256
  integer,parameter :: max_hvarname_len = 64
  integer,parameter :: max_hvar_dimlen  = 5
  integer,parameter :: file_purpose_not_set  = 0
  integer,parameter :: file_purpose_in  = 1
  integer,parameter :: file_purpose_out = 2

  type, public :: hist_coord_t
    character(len=max_hcoordname_len) :: name = ''  ! coordinate name
    integer                  :: dimsize = 0         ! size of dimension
    integer                  :: dimid               ! Unique PIO Id for this dimension
    character(len=max_chars) :: long_name = ''      ! 'long_name' attribute
    character(len=max_chars) :: units = ''          ! 'units' attribute
    logical                  :: is_partitioned      ! whether the dimension is partitioned across ranks
  end type hist_coord_t

  type, public :: hist_var_t
    character(len=max_hvarname_len) :: name   ! coordinate name
    character(len=max_chars) :: long_name     ! 'long_name' attribute
    character(len=max_chars) :: pio_decomp_tag ! PIO decomposition label used by this variable.
    character(len=max_chars) :: units         ! 'units' attribute
    type(var_desc_t) :: piovar                ! netCDF variable ID
    integer          :: dtype                 ! data type used to pass data to read/write routines
    integer          :: nc_dtype              ! data type used in the netcdf files
    integer          :: numdims               ! Number of dimensions in out field
    type(io_desc_t), pointer  :: iodesc       ! PIO decomp associated with this variable
    type(iodesc_list_t), pointer  :: iodesc_list ! PIO decomp list with metadata about PIO decomp
    integer(kind=pio_offset_kind), allocatable :: compdof(:)        ! Global locations in output array for this process
    integer, allocatable :: dimid(:)          ! array of PIO dimension id's for this variable
    integer, allocatable :: dimlen(:)         ! array of PIO dimension lengths for this variable
    logical              :: has_t_dim         ! true, if variable has a time dimension
    logical              :: is_set = .false.  ! Safety measure to ensure a deallocated hist_var_t is never used
    logical              :: is_partitioned    ! Whether at least one of the dims is partitioned
  end type hist_var_t

  ! The iodesc_list allows us to cache existing PIO decompositions
  ! The tag needs the dim lengths, the dtype and map id (+ optional permutation)
  ! Define a recursive structure because we do not know ahead of time how many
  ! decompositions will be require
  type iodesc_list_t
    character(max_chars)         :: tag              ! Unique tag associated with this decomposition
    type(io_desc_t),     pointer :: iodesc => NULL() ! PIO - decomposition
    type(iodesc_list_t), pointer :: next => NULL()   ! Needed for recursive definition, the next list
    type(iodesc_list_t), pointer :: prev => NULL()   ! Needed for recursive definition, the list that points to this one
    logical                      :: iodesc_set = .false.
    integer                      :: num_customers = 0 ! Track the number of currently active variables that use this pio decomposition
    integer                      :: location = 0      ! where am in the recursive list
  end type iodesc_list_t

  ! Define the first iodesc_list_t
  type(iodesc_list_t), pointer :: iodesc_list_top
!----------------------------------------------------------------------
  type hist_coord_list_t
    type(hist_coord_t),      pointer :: coord => NULL() ! Pointer to a history dimension structure
    type(hist_coord_list_t), pointer :: next => NULL()  ! Needed for recursive definition
  end type hist_coord_list_t
!----------------------------------------------------------------------
  type hist_var_list_t
    type(hist_var_t),      pointer :: var => NULL()  ! Pointer to a history variable structure
    type(hist_var_list_t), pointer :: next => NULL() ! Needed for recursive definition
  end type hist_var_list_t
!----------------------------------------------------------------------
  type pio_file_list_t
    type(pio_atm_file_t),  pointer :: pio_file => NULL() ! Pointer to an atm. pio file
    type(pio_file_list_t), pointer :: next => NULL()     ! Needed for recursive definition
    type(pio_file_list_t), pointer :: prev => NULL()     ! A doubly-linked list is easier to handle
  end type pio_file_list_t
  ! Define the first pio_file_list
  type(pio_file_list_t), pointer :: pio_file_list_front
  type(pio_file_list_t), pointer :: pio_file_list_back

!----------------------------------------------------------------------
  type, public :: pio_atm_file_t
    character(len=max_chars) :: filename = ""
    integer :: purpose = file_purpose_not_set       ! Input or Output file
    type(file_desc_t)        :: pioFileDesc         ! Contains data identifying the file.
    type(hist_coord_list_t)  :: coord_list_top      ! Recursive list of variables
    type(hist_var_list_t)    :: var_list_top        ! Recursive list of variables
    integer                  :: numRecs             ! Number of history records on file
    logical                  :: is_enddef = .false. ! Whether definition phase is open
    integer                  :: num_customers       ! The number of customer that requested to open the file.
  end type pio_atm_file_t

!----------------------------------------------------------------------
  interface grid_read_data_array
    module procedure grid_read_darray_double
    module procedure grid_read_darray_float
    module procedure grid_read_darray_int
  end interface grid_read_data_array
!----------------------------------------------------------------------
  interface grid_write_data_array
    module procedure grid_write_darray_float
    module procedure grid_write_darray_double
    module procedure grid_write_darray_int
  end interface
!----------------------------------------------------------------------

contains
!=====================================================================!
  ! Register a PIO file to be used for input/output operations.
  ! If file is already open, ensures file_purpose matches the current one
  subroutine register_file(filename,file_purpose)

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: file_purpose

    type(pio_atm_file_t), pointer :: pio_file

    if (.not.associated(pio_subsystem)) then
      call errorHandle("PIO ERROR: local pio_subsystem pointer has not been established yet.",-999)
    endif

    call get_pio_atm_file(filename,pio_file,file_purpose)
  end subroutine register_file
!=====================================================================!
  ! Mandatory call to finish the variable and dimension definition phase
  ! of a new PIO file.  Once this routine is called it is not possible
  ! to add new dimensions or variables to the file.
  subroutine eam_pio_enddef(filename)

    character(len=*), intent(in) :: filename

    type(pio_atm_file_t), pointer :: current_atm_file
    integer                       :: ierr
    logical                       :: found

    call lookup_pio_atm_file(filename,current_atm_file,found)
    if (.not.found) then
      call errorHandle("PIO ERROR: error running enddef on file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    ! It could happen that we are running a test, with an input file opening the
    ! same file that an output stream just wrote. In this case, the def phase ended
    ! during the output setup.
    ! 
    if (.not. current_atm_file%is_enddef) then
      ! Gather the pio decomposition for all variables in this file, and assign them pointers.
      call set_decomp(trim(filename))
      ! Officially close the definition step for this file.
      ierr = PIO_enddef(current_atm_file%pioFileDesc)
      call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),ierr)
      current_atm_file%is_enddef = .true.
    endif

  end subroutine eam_pio_enddef
!=====================================================================!
  ! Register a dimension with a specific pio output file.  Mandatory inputs
  ! include:
  ! filename:  Name of file to add the dimension to.
  ! shortname: Short name descriptor for this dimension.  This will be
  !            name to find the dimension in the netCDF file.
  ! longname:  A longer character string with a more descriptive name of
  !            the dimension.
  ! length:    The dimension length (must be >=0).  Choosing 0 marks the
  !            dimensions as having "unlimited" length which is used for
  !            dimensions such as time.
  subroutine register_dimension(filename,shortname,longname,length,is_partitioned)
    use pio_types, only: pio_unlimited
    use pio_nf,    only: PIO_def_dim

    character(len=*), intent(in)        :: filename   ! Name of file to register the dimension on.
    character(len=*), intent(in)        :: shortname,longname ! Short- and long- names for this dimension, short: brief identifier and name for netCDF output, long: longer descriptor sentence to be included as meta-data in file.
    integer, intent(in)                 :: length             ! Length of the dimension, 0: unlimited (like time), >0 actual length of dimension
    logical, intent(in)                 :: is_partitioned   ! whether this dimension is partitioned across ranks

    type(pio_atm_file_t), pointer       :: pio_atm_file
    type(hist_coord_t), pointer         :: hist_coord
    type(hist_coord_list_t), pointer    :: curr, prev
    logical                             :: found, dim_found
    integer                             :: ierr

    dim_found = .false.

    ! Make sure the dimension length is reasonable
    if (length<0) call errorHandle("PIO Error: dimension "//trim(shortname)//", can't have a negative dimension length",-999)

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    if (.not.found ) then
      call errorHandle("PIO ERROR: error registering dimension "//trim(shortname)//" in file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    ! Get a new dimension pointer in coord_list
    curr => pio_atm_file%coord_list_top
    do while (associated(curr))
      if (associated(curr%coord)) then
        if(trim(curr%coord%name)==trim(shortname)) then
          dim_found = .true.
          exit
        endif
      end if
      prev => curr
      curr => prev%next
    end do
    ! If the dim was not found, create it
    if (.not. dim_found) then
      allocate(prev%next)
      curr => prev%next
      allocate(curr%coord)
      hist_coord => curr%coord

      ! Register this dimension
      hist_coord%name           = trim(shortname)
      hist_coord%long_name      = trim(longname)
      hist_coord%dimsize        = length
      hist_coord%is_partitioned = is_partitioned

      if (length.eq.0) then
        ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), pio_unlimited , hist_coord%dimid)
        time_dimid = hist_coord%dimid
      else
        ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), length , hist_coord%dimid)
      end if
      call errorHandle("PIO ERROR: could not define dimension "//trim(shortname)//" on file: "//trim(filename),ierr)
    else
      ! The dim was already registered by another input/output instance. Check that everything matches
      hist_coord => curr%coord
      if (trim(hist_coord%name) .ne. trim(shortname) .or. &
          trim(hist_coord%long_name) .ne. trim(longname) .or. &
          hist_coord%dimsize .ne. length) then
        call errorHandle("PIO ERROR: dimension "//trim(shortname)//" was already defined on file "//trim(filename)//", but with a different length",ierr)

      endif
    endif
  end subroutine register_dimension
!=====================================================================!
  ! Register a variable with a specific pio input file. Mandatory inputs
  ! include:
  ! filename:       The name of the netCDF file this variable will be
  !                 registered with.
  ! shortname:      A shortname descriptor (tag) for this variable.  This will be
  !                 used to label the variable in the netCDF file as well.
  ! longname:       A longer character string describing the variable.
  ! numdims:        The number of dimensions associated with this variable,
  !                 including time (if applicable).
  ! var_dimensions: An array of character strings with the dimension shortnames
  !                 for each dimension used by this variable.  Should have
  !                 'numdims' entries.
  ! dtype:          The data type for this variable using the proper netCDF
  !                 integer tag.
  ! pio_decomp_tag: A string that describes this particular dimension
  !                 arrangement which will be used to create a unique PIO
  !                 decomposition for reading this variable.  It is ok to reuse
  !                 the pio_decomp_tag for variables that have the same
  !                 dimensionality.  See get_decomp for more details.
  subroutine get_variable(filename,shortname,longname,numdims,var_dimensions,dtype,pio_decomp_tag)
    use pio_nf, only: PIO_inq_vartype
    character(len=*), intent(in) :: filename         ! Name of the file to register this variable with
    character(len=*), intent(in) :: shortname,longname       ! short and long names for the variable.  Short: variable name in file, Long: more descriptive name
    integer, intent(in)          :: numdims                  ! Number of dimensions for this variable, including time dimension
    character(len=*), intent(in) :: var_dimensions(numdims)  ! String array with shortname descriptors for each dimension of variable.
    integer, intent(in)          :: dtype                    ! datatype for this variable, REAL, DOUBLE, INTEGER, etc.
    character(len=*), intent(in) :: pio_decomp_tag           ! Unique tag for this variables decomposition type, to be used to determine if the io-decomp already exists.

    ! Local variables
    type(pio_atm_file_t),pointer :: pio_atm_file
    type(hist_var_t), pointer    :: hist_var
    integer                      :: dim_ii
    integer                      :: ierr
    logical                      :: found,var_found

    type(hist_var_list_t), pointer :: curr, prev

    var_found = .false.

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    if (.not.found ) then
      call errorHandle("PIO ERROR: error registering variable "//trim(shortname)//" in file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    ! Get a new variable pointer in var_list
    if (len_trim(shortname)>max_hvarname_len) call errorHandle("PIO Error: variable shortname "//trim(shortname)//" is too long, consider increasing max_hvarname_len or changing the variable shortname",-999)
    curr => pio_atm_file%var_list_top

    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(shortname) .and. curr%var%is_set) then
          var_found = .true.
          exit
        endif
      endif
      prev => curr
      curr => prev%next
    end do

    if (.not. var_found) then
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
      hist_var%is_partitioned = .false.
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
      call errorHandle("PIO ERROR: could not find variable "//trim(shortname)//" in file "//trim(filename),ierr)

      ! Not really needed, but just in case, store var data type in the nc file
      ierr = PIO_inq_vartype(pio_atm_file%pioFileDesc,hist_var%piovar,hist_var%nc_dtype)
      call errorHandle("EAM PIO ERROR: Unable to retrieve dtype for variable "//shortname,ierr)

      ! Set that the new variable has been set
      hist_var%is_set = .true.
    else
      ! The var was already registered by another input/output instance. Check that everything matches
      hist_var => curr%var
      if ( trim(hist_var%long_name) .ne. trim(longname) ) then
        ! Different long name
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different longname, in file: "//trim(filename),-999)
      elseif (hist_var%dtype .ne. dtype) then
        ! Different data type
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different dtype, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_out .and. & ! Out files must match the decomp tag
              (hist_var%numdims .ne. numdims .or. &
               trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag))) then
        ! Different decomp tag in output file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      elseif (hist_var%numdims .ne. numdims .and. &
              hist_var%numdims .ne. (numdims+1)) then
        ! Invalid dimlen
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different dimlen, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_in .and. &
              hist_var%numdims .eq. numdims .and. &
              trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag)) then
        ! Same dimlen, but different decomp tag in input file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_in .and. & ! In files *may* use a decomp tag
             (hist_var%numdims .eq. (numdims+1) .and. &           ! without "-time" at the end
              trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag)//"-time")) then
        ! Different dimlen, but different decomp tag even if attaching "-time" in input file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      endif
    endif
  end subroutine get_variable
!=====================================================================!
  ! Register a variable with a specific pio output file. Mandatory inputs
  ! include:
  ! pio_atm_filename: The name of the netCDF file this variable will be
  !                   registered with.
  ! shortname:      A shortname descriptor (tag) for this variable.  This will be
  !                 used to label the variable in the netCDF file as well.
  ! longname:       A longer character string describing the variable.
  ! numdims:        The number of dimensions associated with this variable,
  !                 including time (if applicable).
  ! var_dimensions: An array of character strings with the dimension shortnames
  !                 for each dimension used by this variable.  Should have
  !                 'numdims' entries.
  ! dtype:          The data type for this variable using the proper netCDF
  !                 integer tag.
  ! pio_decomp_tag: A string that describes this particular dimension
  !                 arrangement which will be used to create a unique PIO
  !                 decomposition for reading this variable.  It is ok to reuse
  !                 the pio_decomp_tag for variables that have the same
  !                 dimensionality.  See get_decomp for more details.
  subroutine register_variable(filename,shortname,longname,units, &
                               numdims,var_dimensions,            &
                               dtype,nc_dtype,pio_decomp_tag)
    use pio_nf, only: PIO_def_var

    character(len=*), intent(in) :: filename         ! Name of the file to register this variable with
    character(len=*), intent(in) :: shortname,longname       ! short and long names for the variable.  Short: variable name in file, Long: more descriptive name
    character(len=*), intent(in) :: units                    ! units for variable
    integer, intent(in)          :: numdims                  ! Number of dimensions for this variable, including time dimension
    character(len=*), intent(in) :: var_dimensions(numdims)  ! String array with shortname descriptors for each dimension of variable.
    integer, intent(in)          :: dtype                    ! datatype for arrays that will be passed to read/write routines
    integer, intent(in)          :: nc_dtype                 ! datatype for this variable in nc files
    character(len=*), intent(in) :: pio_decomp_tag           ! Unique tag for this variables decomposition type, to be used to determine if the io-decomp already exists.

    ! Local variables
    type(pio_atm_file_t),pointer :: pio_atm_file
    type(hist_var_t), pointer    :: hist_var
    integer                      :: dim_ii
    logical                      :: found,var_found
    integer                      :: ierr
    character(len=256)           :: dimlen_str
    type(hist_coord_t), pointer  :: hist_coord

    type(hist_var_list_t), pointer :: curr, prev
    
    var_found = .false.

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    if (.not.found ) then
      call errorHandle("PIO ERROR: error registering variable "//trim(shortname)//" in file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    ! Get a new variable pointer in var_list
    if (len_trim(shortname)>max_hvarname_len) call errorHandle("PIO Error: variable shortname "//trim(shortname)//" is too long, consider increasing max_hvarname_len or changing the variable shortname",-999)
    curr => pio_atm_file%var_list_top

    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(shortname) .and. curr%var%is_set) then
          var_found = .true.
          exit
        endif
      end if
      prev => curr
      curr => prev%next
    end do

    ! If the var was not found, allocate the new var
    if (.not. var_found) then
      allocate(prev%next)
      curr => prev%next
      allocate(curr%var)
      hist_var => curr%var
      ! Populate meta-data associated with this variable
      hist_var%name      = trim(shortname)
      hist_var%long_name = trim(longname)
      hist_var%units = trim(units)
      hist_var%numdims   = numdims
      hist_var%dtype     = dtype
      hist_var%nc_dtype  = nc_dtype
      hist_var%pio_decomp_tag = trim(pio_decomp_tag)
      ! Determine the dimension id's saved in the netCDF file and associated with
      ! this variable, check if variable has a time dimension
      hist_var%has_t_dim = .false.
      hist_var%is_partitioned = .false.
      allocate(hist_var%dimid(numdims),hist_var%dimlen(numdims))
      do dim_ii = 1,numdims
        ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(var_dimensions(dim_ii)),hist_var%dimid(dim_ii))
        call errorHandle("EAM_PIO ERROR: Unable to find dimension id for "//trim(var_dimensions(dim_ii)),ierr)
        ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,hist_var%dimid(dim_ii),hist_var%dimlen(dim_ii))
        call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension "//trim(var_dimensions(dim_ii)),ierr)
        if (hist_var%dimlen(dim_ii).eq.0) hist_var%has_t_dim = .true.
        call convert_int_2_str(hist_var%dimlen(dim_ii),dimlen_str)
        hist_var%pio_decomp_tag = hist_var%pio_decomp_tag//"_"//trim(dimlen_str)

        call get_coord (filename,var_dimensions(dim_ii),hist_coord)
        if (hist_coord%is_partitioned) then
          hist_var%is_partitioned = .true.
        endif
      end do

      ierr = PIO_def_var(pio_atm_file%pioFileDesc, trim(shortname), hist_var%nc_dtype, hist_var%dimid(:numdims), hist_var%piovar)
      call errorHandle("PIO ERROR: could not define variable "//trim(shortname),ierr)

      !PMC
      ierr=PIO_put_att(pio_atm_file%pioFileDesc, hist_var%piovar, 'units', hist_var%units )
      ierr=PIO_put_att(pio_atm_file%pioFileDesc, hist_var%piovar, 'long_name', hist_var%long_name )

      ! Set that new variable has been created
      hist_var%is_set = .true.
    else
      ! The var was already registered by another input/output instance. Check that everything matches
      hist_var => curr%var
      if ( trim(hist_var%long_name) .ne. trim(longname) ) then
        ! Different long name
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different longname, in file: "//trim(filename),-999)
     elseif ( trim(hist_var%units) .ne. trim(units) ) then
        ! Different units
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different units, in file: "//trim(filename),-999)
      elseif (hist_var%nc_dtype .ne. nc_dtype .or. hist_var%dtype .ne. dtype) then
        ! Different data type
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different dtype, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_out .and. & ! Out files must match the decomp tag
              (hist_var%numdims .ne. numdims .or. &
               trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag))) then
        ! Different decomp tag in output file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      elseif (hist_var%numdims .ne. numdims .and. &
              hist_var%numdims .ne. (numdims+1)) then
        ! Invalid dimlen
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different dimlen, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_in .and. &
              hist_var%numdims .eq. numdims .and. &
              trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag)) then
        ! Same dimlen, but different decomp tag in input file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      elseif (pio_atm_file%purpose .eq. file_purpose_in .and. & ! In files *may* use a decomp tag
             (hist_var%numdims .eq. (numdims+1) .and. &           ! without "-time" at the end
              trim(hist_var%pio_decomp_tag) .ne. trim(pio_decomp_tag)//"-time")) then
        ! Different dimlen, but different decomp tag even if attaching "-time" in input file
        call errorHandle("PIO Error: variable "//trim(shortname)//", already registered with different decomp tag, in file: "//trim(filename),-999)
      endif
    endif

  end subroutine register_variable
!=====================================================================!
  subroutine set_variable_metadata(filename, varname, metaname, metaval)
    use pionfatt_mod, only: PIO_put_att => put_att

    character(len=256), intent(in) :: filename
    character(len=256), intent(in) :: varname
    character(len=256), intent(in) :: metaname
    character(len=256), intent(in) :: metaval

    ! Local variables
    type(pio_atm_file_t),pointer :: pio_file
    type(hist_var_t),    pointer :: var
    integer                      :: ierr
    logical                      :: found

    type(hist_var_list_t), pointer :: curr

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(filename),pio_file,found)
    if (.not.found ) then
      call errorHandle("PIO ERROR: error setting metadata for variable "//trim(varname)//" in file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    ! Find the variable in the file
    curr => pio_file%var_list_top

    found = .false.
    do while (associated(curr))
      if (associated(curr%var)) then
        if (trim(curr%var%name)==trim(varname) .and. curr%var%is_set) then
          found = .true.
          var => curr%var
          exit
        endif
      endif
      curr => curr%next
    end do
    if (.not.found ) then
      call errorHandle("PIO ERROR: error setting metadata for variable "//trim(varname)//" in file "//trim(filename)//".\n Variable not found.",-999)
    endif

    ierr = PIO_put_att(pio_file%pioFileDesc, var%piovar, metaname, metaval)
    if (ierr .ne. 0) then
      call errorHandle("Error setting attribute '" // trim(metaname) &
                       // "' on variable '" // trim(varname) &
                       // "' in pio file " // trim(filename) // ".", -999)
    endif

  end subroutine set_variable_metadata
!=====================================================================!
  ! Update the time dimension for a specific PIO file.  This is needed when
  ! reading or writing multiple time levels.  Unlimited dimensions are treated
  ! differently in netCDF than typical static length variables.  Note, here
  ! "time" is hardcoded as the only unlimited variable.  If, in the future,
  ! scream decides to allow for other "unlimited" dimensions to be used our
  ! input/output than this routine will need to be adjusted.
  subroutine eam_update_time(filename,time)
    use pionfput_mod, only: PIO_put_var   => put_var

    character(len=*), intent(in) :: filename       ! PIO filename
    real(c_double), intent(in)   :: time

    type(hist_var_t), pointer    :: var
    type(pio_atm_file_t),pointer :: pio_atm_file
    integer                      :: ierr
    logical                      :: found

    call lookup_pio_atm_file(filename,pio_atm_file,found)
    pio_atm_file%numRecs = pio_atm_file%numRecs + 1
    call get_var(pio_atm_file,'time',var)
    ! Only update time on the file if a valid time is provided
    if (time>=0) ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,(/ pio_atm_file%numRecs /), (/ 1 /), (/ time /))
  end subroutine eam_update_time
!=====================================================================!
  ! Assign header metadata to a specific pio output file.  TODO: Fix this to be
  ! more general.  Right now it is all dummy boiler plate.  Would make the most
  ! sense to pass a structure with all of the relevant header info contained
  ! within it.
  subroutine eam_pio_createHeader(File)

    type(file_desc_t), intent(in) :: File             ! Pio file Handle
    integer                       :: retval

    ! TODO change options below to match specific simulation case
    retval=pio_put_att (File, PIO_GLOBAL, 'source', 'E3SM Atmosphere Model Version 4')
    retval=pio_put_att (File, PIO_GLOBAL, 'case', 'TEST 1') ! NEED TO FIX THIS!!!
    retval=pio_put_att (File, PIO_GLOBAL, 'title', 'EAMv4 History File')
    retval=pio_put_att (File, PIO_GLOBAL, 'git_hash','THE GIT LOG HASH')  ! NEED TO FIX THIS!!!
    retval=pio_put_att (File, PIO_GLOBAL, 'host', 'THE HOST')  ! NEED TO FIX THIS!!!
    retval=pio_put_att (File, PIO_GLOBAL, 'Version', '1.0')
    retval=pio_put_att (File, PIO_GLOBAL, 'revision_Id', 'None')  !WHAT IS THIS? NOT IN EAM.
    retval=pio_put_att (File, PIO_GLOBAL, 'initial_file', 'NONE FOR NOW')  !NEED TO FIX THIS
    retval=pio_put_att (File, PIO_GLOBAL, 'topography_file', 'NONE FOR NOW')  !NEED TO FIX THIS
    retval=pio_put_att (File, PIO_GLOBAL, 'contact', 'e3sm-data-support@llnl.gov')
    retval=pio_put_att (File, PIO_GLOBAL, 'institution_id', 'E3SM-Project')
    retval=pio_put_att (File, PIO_GLOBAL, 'product', 'model-output')
    retval=pio_put_att (File, PIO_GLOBAL, 'realm','atmos')
    retval=pio_put_att (File, PIO_GLOBAL, 'Conventions','None yet')
    retval=pio_put_att (File, PIO_GLOBAL, 'institution', 'LLNL (Lawrence Livermore National Laboratory, &
    &Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven &
    &National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, &
    &NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak &
    &Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, &
    &Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). Mailing &
    &address: LLNL Climate Program, c/o David C. Bader, Principal Investigator, L-103, 7000 East Avenue, &
    &Livermore, CA 94550, USA')

  end subroutine eam_pio_createHeader
!=====================================================================!
  ! Ensures a pio system is in place, by either creating a new one
  ! or getting the one created by CIME (for CIME builds)
  subroutine eam_init_pio_subsystem(mpicom,atm_id)
#ifdef SCREAM_CIME_BUILD
    use shr_pio_mod,  only: shr_pio_getrearranger, shr_pio_getiosys, &
                            shr_pio_getiotype, shr_pio_getioformat
#else
    use pio_types,  only: pio_rearr_subset, PIO_iotype_netcdf, PIO_64BIT_DATA
    use piolib_mod, only: pio_init

    integer :: ierr, stride, atm_rank, atm_size, num_aggregator
#endif

    integer, intent(in) :: mpicom
    integer, intent(in) :: atm_id

    if (associated(pio_subsystem)) call errorHandle("PIO ERROR: local pio_subsystem pointer has already been established.",-999)

    atm_mpicom = mpicom

#ifdef SCREAM_CIME_BUILD
    pio_subsystem  => shr_pio_getiosys(atm_id)
    pio_iotype     = shr_pio_getiotype(atm_id)
    pio_rearranger = shr_pio_getrearranger(atm_id)
    pio_mode       = shr_pio_getioformat(atm_id)
#else
    ! WARNING: we're assuming *every atm rank* is an I/O rank
    call MPI_Comm_rank(atm_mpicom, atm_rank, ierr)
    call MPI_Comm_size(atm_mpicom, atm_size, ierr)

    ! Just for removing unused dummy warnings
    if (.false.) print *, atm_id

    stride = 1
    num_aggregator = 0

    allocate(pio_subsystem)
    pio_rearranger = pio_rearr_subset
    pio_iotype     = PIO_iotype_netcdf
    pio_mode       = PIO_64BIT_DATA ! Default to 64 bit
    call PIO_init(atm_rank, atm_mpicom, atm_size, num_aggregator, stride, &
                  pio_rearr_subset, pio_subsystem, base=0)
#endif

    ! Init the list of pio files so that begin==end==null
    pio_file_list_back   => null()
    pio_file_list_front => null()

    ! Init the iodecomp 
    iodesc_list_top => null()

  end subroutine eam_init_pio_subsystem
!=====================================================================!
  ! Query whether the pio subsystem is inited already
  ! This can be useful to avoid double-init or double-finalize calls.
  function is_eam_pio_subsystem_inited() result(is_it) bind(c)
    use iso_c_binding, only: c_bool

    logical(kind=c_bool) :: is_it

    is_it = LOGICAL(associated(pio_subsystem),kind=c_bool)
  end function is_eam_pio_subsystem_inited
!=====================================================================!
  ! Create a pio netCDF file with the appropriate name.
  subroutine eam_pio_createfile(File,fname)
    use piolib_mod, only: pio_createfile
    use pio_types,  only: pio_clobber

    type(file_desc_t), intent(inout) :: File             ! Pio file Handle
    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    integer                          :: retval           ! PIO error return value
    integer                          :: mode             ! Mode for how to handle the new file

    mode = ior(pio_mode,pio_clobber) ! Set to CLOBBER for now, TODO: fix to allow for optional mode type like in CAM
    retval = pio_createfile(pio_subsystem,File,pio_iotype,fname,mode)
    call errorHandle("PIO ERROR: unable to create file: "//trim(fname),retval)

  end subroutine eam_pio_createfile
!=====================================================================!
  ! Open an already existing netCDF file.
  subroutine eam_pio_openfile(pio_file,fname)
    use piolib_mod, only: pio_openfile
    use pio_types,  only: pio_write, pio_nowrite

    type(pio_atm_file_t), pointer, intent(in) :: pio_file     ! Pointer to pio file struct associated with this filename
    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    integer                          :: retval           ! PIO error return value
    integer                          :: mode             ! Mode for how to handle the new file

    if (pio_file%purpose .eq. file_purpose_in) then
      mode = pio_nowrite
    else
      mode = pio_write
    endif
    retval = pio_openfile(pio_subsystem,pio_file%pioFileDesc,pio_iotype,fname,mode)
    call errorHandle("PIO ERROR: unable to open file: "//trim(fname),retval)

  end subroutine eam_pio_openfile
!=====================================================================!
  ! Close a netCDF file.  To be done as a last step after all input or output
  ! for that file has been finished.
  subroutine eam_pio_closefile(fname)
    use piolib_mod, only: PIO_syncfile, PIO_closefile

    character(len=*),  intent(in)    :: fname            ! Pio file name
    !--
    type(pio_atm_file_t),pointer     :: pio_atm_file
    type(pio_file_list_t), pointer   :: pio_file_list_ptr
    logical                          :: found
    type(hist_var_list_t), pointer   :: curr_var_list
    type(hist_var_t), pointer        :: var

    ! Find the pointer for this file
    call lookup_pio_atm_file(trim(fname),pio_atm_file,found,pio_file_list_ptr)
    if (found) then
      if (pio_atm_file%num_customers .eq. 1) then
        if (pio_atm_file%purpose .eq. file_purpose_out) then
          call PIO_syncfile(pio_atm_file%pioFileDesc)
        endif
        call PIO_closefile(pio_atm_file%pioFileDesc)
        pio_atm_file%num_customers = pio_atm_file%num_customers - 1

        ! Remove all variables from this file as customers for the stored pio
        ! decompostions
        curr_var_list => pio_atm_file%var_list_top  ! Start with the first variable in the file
        do while (associated(curr_var_list))
          var => curr_var_list%var  ! The actual variable pointer
          if (associated(var)) then
            ! Remove this variable as a customer of the associated iodesc
            var%iodesc_list%num_customers = var%iodesc_list%num_customers - 1
            ! Dellocate select memory from this variable.  Note we can't just
            ! deallocate the whole var structure because this would also
            ! deallocate the iodesc_list.
            call deallocate_hist_var_t(var)
          end if ! associated(var)
          curr_var_list => curr_var_list%next  ! Move on to the next variable
        end do ! associated(curr_var_list)

        ! Adjust pointers in the pio file list
        if (associated(pio_file_list_ptr%prev)) then
          pio_file_list_ptr%prev%next => pio_file_list_ptr%next
        else
          ! We're deleting the first item in the lists. Update pio_file_list_front
          pio_file_list_front => pio_file_list_ptr%next
        endif
        if (associated(pio_file_list_ptr%next)) then
          pio_file_list_ptr%next%prev => pio_file_list_ptr%prev
        else
          ! We're deleting the last item in the lists. Update pio_file_list_back
          pio_file_list_back => pio_file_list_ptr%prev
        endif

        ! Now that we have closed this pio file and purged it from the list we
        ! can deallocate the structure.
        deallocate(pio_atm_file)
      else if (pio_atm_file%num_customers .gt. 1) then
        pio_atm_file%num_customers = pio_atm_file%num_customers - 1
      else
        call errorHandle("PIO ERROR: while closing file: "//trim(fname)//", found num_customers<=0",-999)
      endif
    else
      call errorHandle("PIO ERROR: unable to close file: "//trim(fname)//", was not found",-999)
    end if

    ! Final step, free any pio decompostion memory that is no longer needed.
    !   Update: We are trying to reuse decompostions maximally, so we're
    ! skipping this step.
    !call free_decomp()

  end subroutine eam_pio_closefile
!=====================================================================!
  ! Helper function to debug list of decomps 
  subroutine print_decomp()
    type(iodesc_list_t),   pointer :: iodesc_ptr

    integer :: total
    integer :: cnt 
    logical :: assoc

    if (associated(iodesc_list_top)) then
      total = 0
      cnt   = 0
      write(*,*) "            PRINT DECOMP            "
      write(*,*) "            ------------            "
      write(*,'(8X,A10,A15,A50,A15)') "Location", "Associated?", "IODESC TAG", "# Customers"
    else
      write(*,*) "No DECOMP List to print"
      return
    end if
    iodesc_ptr => iodesc_list_top
    do while(associated(iodesc_ptr))
      total = total + 1
      assoc = .false.
      if (associated(iodesc_ptr%iodesc).and.iodesc_ptr%iodesc_set) then
        cnt = cnt + 1
        assoc = .true.
      end if
      write(*,'(I3,A5,I10,L15,A50,I15)') total, ": ", iodesc_ptr%location, assoc , trim(iodesc_ptr%tag), iodesc_ptr%num_customers
      iodesc_ptr => iodesc_ptr%next
    end do
    write(*,*) "            total: ", total, ", cnt: ",cnt
    write(*,*) "            ------------            "

  end subroutine print_decomp
!=====================================================================!
  ! Free the memory stored in a hist_var_t derived type
  subroutine deallocate_hist_var_t(var)

    type(hist_var_t), pointer        :: var
    
    deallocate(var%compdof)
    deallocate(var%dimid)
    deallocate(var%dimlen)
    var%is_set = .false.

  end subroutine deallocate_hist_var_t
  !=====================================================================!
  ! Free pio decomposition memory in PIO for any decompositions that are no
  ! longer needed.  Previously, we thought that this is an important memory
  ! management step that should be taken whenever a file is closed.  Now we're
  ! trying to keep decomps persistent so they can be reused.  Thus, calling
  ! this routine is optional.
  subroutine free_decomp()
    use piolib_mod, only: PIO_freedecomp
    type(iodesc_list_t),   pointer :: iodesc_ptr, next

    ! Free all decompositions from PIO
    iodesc_ptr => iodesc_list_top
    do while(associated(iodesc_ptr))
      next => iodesc_ptr%next
      if (associated(iodesc_ptr%iodesc).and.iodesc_ptr%iodesc_set) then
        if (iodesc_ptr%num_customers .eq. 0) then
          ! Free decomp
          call pio_freedecomp(pio_subsystem,iodesc_ptr%iodesc)
          ! Nullify this decomp
          ! If we are at iodesc_list_top we need to make iodesc_ptr%next the new
          ! iodesc_list_top:
          if (associated(iodesc_ptr%prev)) then
            iodesc_ptr%prev%next => iodesc_ptr%next
          else
            ! We are deleting the first item in the list, update
            ! iodesc_list_front
            iodesc_list_top => iodesc_ptr%next
          end if
          if (associated(iodesc_ptr%next)) then
            iodesc_ptr%next%prev => iodesc_ptr%prev
          end if
          deallocate(iodesc_ptr)
        end if
      end if
      iodesc_ptr => next
    end do

  end subroutine free_decomp
!=====================================================================!
  ! Finalize a PIO session within scream.  Close all open files and deallocate
  ! the pio_subsystem session.
  subroutine eam_pio_finalize()
    use piolib_mod, only: PIO_finalize, pio_freedecomp
    ! May not be needed, possibly handled by PIO directly.

#if !defined(SCREAM_CIME_BUILD)
    integer :: ierr
#endif
    type(pio_file_list_t), pointer :: curr_file_ptr, prev_file_ptr
    type(iodesc_list_t),   pointer :: iodesc_ptr

    ! Close all the PIO Files
    curr_file_ptr => pio_file_list_front
    do while (associated(curr_file_ptr))
      call eam_pio_closefile(curr_file_ptr%pio_file%filename)
      prev_file_ptr => curr_file_ptr
      curr_file_ptr => curr_file_ptr%next
      deallocate(prev_file_ptr)
    end do
    ! Free all decompositions from PIO
    iodesc_ptr => iodesc_list_top
    do while(associated(iodesc_ptr))
      if (associated(iodesc_ptr%iodesc).and.iodesc_ptr%iodesc_set) then
        call pio_freedecomp(pio_subsystem,iodesc_ptr%iodesc)
      end if
      iodesc_ptr => iodesc_ptr%next
    end do

#if !defined(SCREAM_CIME_BUILD)
    call PIO_finalize(pio_subsystem, ierr)
    nullify(pio_subsystem)
#endif

  end subroutine eam_pio_finalize
!=====================================================================!
  ! Handle any errors that occur in this module and print to screen an error
  ! message.
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

    integer :: ierr

    if (retVal .ne. PIO_NOERR) then
      write(*,'(I8,2x,A200)') retVal,trim(errMsg)
      ! Kill run
      call eam_pio_finalize()
      call finalize_scream_session()
      call mpi_abort(atm_mpicom,retVal,ierr)
    end if

  end subroutine errorHandle
!=====================================================================!
 ! Determine the unique pio_decomposition for this output grid, if it hasn't
 ! been defined create a new one.
  subroutine get_decomp(tag,dtype,dimension_len,compdof,iodesc_list)
    use piolib_mod, only: pio_initdecomp
    ! TODO: CAM code creates the decomp tag for the user.  Theoretically it is
    ! unique because it is based on dimensions and datatype.  But the tag ends
    ! up not being very descriptive.  The todo item is to revisit how tags are
    ! handled and decide if we want the code to create a tag or let the use
    ! assign a tag.
    character(len=*)          :: tag              ! Unique tag string describing this output grid
    integer, intent(in)       :: dtype            ! Datatype associated with the output
    integer, intent(in)       :: dimension_len(:) ! Array of the dimension lengths for this decomp
    integer(kind=pio_offset_kind), intent(in) :: compdof(:)       ! The degrees of freedom this rank is responsible for
    type(iodesc_list_t), pointer :: iodesc_list   ! The pio decomposition list that holds this iodesc 

    logical                     :: found            ! Whether a decomp has been found among the previously defined decompositions
    type(iodesc_list_t),pointer :: curr, prev       ! Used to toggle through the recursive list of decompositions
    integer                     :: loc_len          ! Used to keep track of how many dimensions there are in decomp

    ! Assign a PIO decomposition to variable, if none exists, create a new one:
    found = .false.
    curr => iodesc_list_top
    prev => iodesc_list_top
    ! Cycle through all current iodesc to see if the decomp has already been
    ! created
    do while(associated(curr) .and. (.not.found))
      if (trim(tag) == trim(curr%tag)) then
        found = .true.
      else
        prev => curr
        curr => curr%next
      end if
    end do
    ! If we didn't find an iodesc then we need to create one
    if (.not.found) then
      curr => prev ! Go back and allocate the new iodesc in curr%next
      ! We may have no iodesc to begin with, so we need to associate the
      ! beginning of the list.
      if(.not.associated(curr)) then
        allocate(curr)
        curr%location = 1
        iodesc_list_top => curr
      end if
      if(associated(curr%iodesc)) then
        allocate(curr%next)
        curr%next%prev => curr
        curr => curr%next
        nullify(curr%iodesc)  ! Extra step to ensure clean iodesc
        nullify(curr%next)  ! Extra step to ensure clean iodesc
      end if
      allocate(curr%iodesc)
      curr%tag = trim(tag)
      if (associated(curr%prev)) then
        curr%location = prev%location+1
      end if
      loc_len = size(dimension_len)
      if ( loc_len.eq.1 .and. dimension_len(loc_len).eq.0 ) then
        allocate(curr%iodesc)
      else
        call pio_initdecomp(pio_subsystem, dtype, dimension_len, compdof, curr%iodesc, rearr=pio_rearranger)
        curr%iodesc_set = .true.
      end if
    end if
    iodesc_list => curr

  end subroutine get_decomp
!=====================================================================!
  ! Set the degrees of freedom (dof) this MPI rank is responsible for
  ! reading/writing from/to file.
  ! Briefly, PIO distributes the work of reading/writing the total array
  ! of any variable data over all of the MPI ranks assigned to PIO.
  ! DOF's are assigned considering a 1D flattening of any array, ignoring
  ! unlimited dimensions, which are handled elsewhere.
  ! Once it is known which dof a rank is responsible for this routine is used to
  ! set those locally for use with all read and write statements.
  ! Arguments:
  ! filename: name of PIO input/output file.
  ! varname:  shortname of variable to set the dof for.
  ! dof_len:  number of dof associated for this rank and this variable
  ! dof_vec:  a vector of length dof_len which includes the global indices of
  !           the flattened "varname" array that this MPI rank is responsible
  !           for.
  ! --- Example: If variable F has dimensions (x,y,t) of size (2,5,t), which means a ---
  ! total of 2x5=10 dof, not counting the time dimension.
  ! If 3 MPI ranks are used then a typical breakdown for dof might be:
  ! Rank 1: (1,2,3)
  ! Rank 2: (4,5,6)
  ! Rank 3: (7,8,9,10)
  subroutine set_dof(filename,varname,dof_len,dof_vec)
    character(len=*), intent(in)              :: filename
    character(len=*), intent(in)              :: varname
    integer, intent(in)                       :: dof_len
    integer(kind=pio_offset_kind), intent(in) :: dof_vec(dof_len)

    type(pio_atm_file_t),pointer              :: pio_atm_file
    type(hist_var_t), pointer                 :: var
    logical                                   :: found
    integer                                   :: ii

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)
    if (allocated(var%compdof)) deallocate(var%compdof)
    allocate( var%compdof(dof_len) )
    do ii = 1,dof_len
      var%compdof(ii) = dof_vec(ii)
    end do

  end subroutine set_dof
!=====================================================================!
  ! Get and assign all pio decompositions for a specific PIO file.  This is a
  ! mandatory step to be taken after all dimensions and variables have been
  ! registered with an input or output file.
  subroutine set_decomp(filename)

    character(len=*)               :: filename  ! Name of the pio file to set decomp for

    type(pio_atm_file_t), pointer  :: current_atm_file
    type(hist_var_list_t), pointer :: curr     ! Used to cycle through recursive list of variables
    type(hist_var_t), pointer      :: hist_var ! Pointer to the variable structure that has been found
    integer                        :: loc_len
    logical                        :: found

    call lookup_pio_atm_file(filename,current_atm_file,found)
    if (.not. found) then
      call errorHandle("PIO ERROR: pio file '"//trim(filename)//"' not found.",999)
    endif

    curr => current_atm_file%var_list_top
    do while (associated(curr))
      if (associated(curr%var)) then
        hist_var => curr%var
        if (.not.associated(hist_var)) call errorHandle("PIO ERROR: unable to set decomp for file, var: "//trim(current_atm_file%filename)//", "//trim(hist_var%name)//". Set DOF.",999)
        ! Assign decomp
        if (hist_var%has_t_dim) then
          loc_len = max(1,hist_var%numdims-1)
          call get_decomp(hist_var%pio_decomp_tag,hist_var%dtype,hist_var%dimlen(:loc_len),hist_var%compdof,hist_var%iodesc_list)
        else
          call get_decomp(hist_var%pio_decomp_tag,hist_var%dtype,hist_var%dimlen,hist_var%compdof,hist_var%iodesc_list)
        end if
        hist_var%iodesc => hist_var%iodesc_list%iodesc
        hist_var%iodesc_list%num_customers = hist_var%iodesc_list%num_customers + 1  ! Add this variable as a customer of this pio decomposition
      end if
      curr => curr%next
    end do

  end subroutine set_decomp
!=====================================================================!
  ! Query the hist_var_t pointer for a specific variable on a specific file.
  subroutine get_var(pio_file,varname,var)

    type(pio_atm_file_t), pointer  :: pio_file ! Pio output file structure
    character(len=*)               :: varname  ! Name of the variable to query
    type(hist_var_t), pointer      :: var      ! Pointer to the variable structure that has been found

    type(hist_var_list_t), pointer :: curr     ! Used to cycle through recursive list of variables

    curr => pio_file%var_list_top
    do while (associated(curr))
      var => curr%var
      if (associated(var)) then
        if (trim(varname) == trim(var%name) .and. var%is_set) return
      end if
      curr => curr%next
    end do

    ! If we got this far we didn't find the variable
    call errorHandle("PIO ERROR: unable to find variable: "//trim(varname)//" in file: "//trim(pio_file%filename),999)

  end subroutine get_var
!=====================================================================!
  ! Retrieves an integer global attribute from the nc file
  function get_int_attribute (file_name, attr_name) result(val)
    use pionfatt_mod, only: PIO_get_att => get_att
    character(len=*), intent(in) :: file_name  ! Name of the filename
    character(len=*), intent(in) :: attr_name  ! Name of the attribute
    type(pio_atm_file_t), pointer :: pio_atm_file
    integer :: val, ierr
    logical :: found

    call lookup_pio_atm_file(trim(file_name),pio_atm_file,found)
    if (.not.found) then
      call errorHandle("PIO Error: can't find pio_atm_file associated with file: "//trim(file_name),-999)
    endif
    ierr = PIO_get_att(pio_atm_file%pioFileDesc, PIO_GLOBAL, attr_name, val)
    if (ierr .ne. 0) then
      call errorHandle("Error retrieving global attribute '" // trim(attr_name) &
                       // "' in pio file " // trim(file_name) // ".", -999)
    endif
  end function get_int_attribute

  ! Writes an integer global attribute to the nc file
  subroutine set_int_attribute (file_name, attr_name, val)
    use pionfatt_mod, only: PIO_put_att => put_att
    use pio_nf, only: pio_redef, PIO_inq_att

    character(len=*), intent(in) :: file_name  ! Name of the filename
    character(len=*), intent(in) :: attr_name  ! Name of the attribute
    integer, intent(in) :: val
    type(pio_atm_file_t), pointer :: pio_atm_file
    integer(pio_offset_kind) :: len
    integer :: ierr,xtype
    logical :: found, enddef_needed

    call lookup_pio_atm_file(trim(file_name),pio_atm_file,found)
    if (.not.found) then
      call errorHandle("PIO Error: can't find pio_atm_file associated with file: "//trim(file_name),-999)
    endif

    ! If this attribute does not exist, we need to re-open the nc file for definition,
    ! then re-close it to put it in data mode again.
    ! NOTE: this check step is only for pre-NetCDF4 format, where attributes can
    !       only be defined while in 'define' mode. For NetCDF4/HDF5, attributes
    !       can be defined at any time.
    ! TODO: add check on netcdf format, to see if this inq_att shenanigans is needed.
    ierr = PIO_inq_att(pio_atm_file%pioFileDesc,PIO_GLOBAL,attr_name,xtype,len)
    enddef_needed = .false.
    if (ierr .ne. PIO_NOERR) then
      ! In theory, there are several reason why this could fail. However, pio.F90
      ! does *not* expose all the nc error codes like pio.h does (e.g., no PIO_ENOTATT).
      ! So we just *assume* that the attribute was not found, and try to define it
      ierr = PIO_redef(pio_atm_file%pioFileDesc)
      if (ierr .ne. 0) then
        call errorHandle("Error while re-opening pio file " // trim(file_name) // ".", -999)
      endif
      enddef_needed = .true.
    endif
    ierr = PIO_put_att(pio_atm_file%pioFileDesc, PIO_GLOBAL, attr_name, val)
    if (ierr .ne. 0) then
      call errorHandle("Error setting global attribute '" // trim(attr_name) &
                       // "' in pio file " // trim(file_name) // ".", -999)
    endif
    if (enddef_needed) then
      ierr = PIO_enddef(pio_atm_file%pioFileDesc)
      if (ierr .ne. 0) then
        call errorHandle("Error while re-closing pio file " // trim(file_name) // ".", -999)
      endif
    endif
  end subroutine set_int_attribute
!=====================================================================!
  ! Lookup pointer for pio file based on filename.
  subroutine lookup_pio_atm_file(filename,pio_file,found,pio_file_list_ptr_in)

    character(len=*),intent(in)   :: filename     ! Name of file to be found
    type(pio_atm_file_t), pointer :: pio_file     ! Pointer to pio_atm_output structure associated with this filename
    logical, intent(out)          :: found        ! whether or not the file was found
    type(pio_file_list_t), pointer, optional :: pio_file_list_ptr_in

    type(pio_file_list_t), pointer :: pio_file_list_ptr

    ! Scan pio file list, search for this filename
    found = .false.
    pio_file_list_ptr => pio_file_list_front
    pio_file => null()
    do while (associated(pio_file_list_ptr))
      if (trim(filename)==trim(pio_file_list_ptr%pio_file%filename)) then
        pio_file => pio_file_list_ptr%pio_file
        found = .true.
        if (present(pio_file_list_ptr_in)) then
          pio_file_list_ptr_in => pio_file_list_ptr
        endif
        return
      end if
      pio_file_list_ptr => pio_file_list_ptr%next
    end do

  end subroutine lookup_pio_atm_file
!=====================================================================!
  ! Create a new pio file pointer based on filename.
  subroutine get_pio_atm_file(filename,pio_file,purpose)

    character(len=*),intent(in)   :: filename     ! Name of file to be found
    type(pio_atm_file_t), pointer :: pio_file     ! Pointer to pio_atm_output structure associated with this filename
    integer,intent(in)            :: purpose      ! Purpose for this file lookup, 0 = find already existing, 1 = create new as output, 2 = open new as input

    logical                        :: found
    type(pio_file_list_t), pointer :: new_list_item

    integer                        :: ierr, time_id

    ! Sanity check
    if (purpose .ne. file_purpose_in .and. purpose .ne. file_purpose_out) then
      call errorHandle("PIO Error: unrecognized file purpose for file '"//filename//"'.",-999)
    endif

    ! If the file already exists, return that file
    call lookup_pio_atm_file(trim(filename),pio_file,found)
    if (found) then
      if (purpose .ne. file_purpose_in .or. &
          pio_file%purpose .ne. file_purpose_in ) then
        ! We only allow multiple customers of the file if they all use it in read mode.
        call errorHandle("PIO Error: file '"//trim(filename)//"' was already open for writing.",-999)
      else
        pio_file%purpose = purpose
        call eam_pio_openfile(pio_file,trim(pio_file%filename))
        pio_file%num_customers = pio_file%num_customers + 1
      endif
    else
      allocate(new_list_item)
      allocate(new_list_item%pio_file)
      pio_file => new_list_item%pio_file

      ! Create and initialize the new pio file:
      pio_file%filename = trim(filename)
      pio_file%numRecs = 0
      pio_file%num_customers = 1
      pio_file%purpose = purpose
      if (purpose == file_purpose_out) then  ! Will be used for output.  Set numrecs to zero and create the new file.
        call eam_pio_createfile(pio_file%pioFileDesc,trim(pio_file%filename))
        call eam_pio_createHeader(pio_file%pioFileDesc)
      elseif (purpose == file_purpose_in) then ! Will be used for input, just open it
        call eam_pio_openfile(pio_file,trim(pio_file%filename))
        pio_file%is_enddef = .true. ! Files open in read mode are in data mode already
        ! Update the numRecs to match the number of recs in this file.
        ierr = pio_inq_dimid(pio_file%pioFileDesc,"time",time_id)
        if (ierr.ne.0) then
          ! time is not present in the file, set the number of Recs to 1
          pio_file%numRecs = 1
        else
          ! time is present, set the number of Recs accordingly
          ierr = pio_inq_dimlen(pio_file%pioFileDesc,time_id,pio_file%numRecs)
          call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension time in file "//trim(pio_file%filename),ierr)
        end if
      else
        call errorHandle("PIO Error: get_pio_atm_file with filename = "//trim(filename)//", purpose (int) assigned to this lookup is not valid" ,-999)
      end if

      ! Update the pio file list
      if (associated(pio_file_list_back)) then
        ! 1) Link new file to the new_list_itement back of the list
        new_list_item%prev => pio_file_list_back
        ! 2) Link the current last element of the list to the new one
        pio_file_list_back%next => new_list_item
        ! 3) and update the pointer to the last
        pio_file_list_back => new_list_item
      else
        ! The list was empty. Set both front/back to point to the new item
        pio_file_list_front => new_list_item
        pio_file_list_back  => new_list_item
      endif
    endif
  end subroutine get_pio_atm_file
!=====================================================================!
  ! Retrieve the time value for a specific time_index
  ! If the input arg time_index is <= 0 then it is assumed the user wants
  ! the last time entry.
  function read_time_at_index(filename,time_index) result(val)
    use pio,          only: PIO_get_var
    use pio_nf,       only: PIO_inq_varid
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: time_index
    real(c_double)               :: val
    real(c_double)               :: val_buf(1)
    
    type(pio_atm_file_t), pointer :: pio_atm_file
    logical                       :: found
    integer                       :: dim_id, time_len, ierr
    type(var_desc_t)              :: varid ! netCDF variable ID
    integer                       :: strt(1), cnt(1)

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    if (.not.found) call errorHandle("read_time_at_index ERROR: File "//trim(filename)//" not found",-999)
    ierr = PIO_inq_varid(pio_atm_file%pioFileDesc,"time",varid)
    call errorHandle('read_time_at_index: Error finding variable ID for "time" in file '//trim(filename)//'.',ierr);

    ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim("time"),dim_id)
    call errorHandle("read_time_at_index ERROR: dimension 'time' not found in file "//trim(filename)//".",ierr)
    ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,dim_id,time_len)
    if (time_index .gt. time_len) then
      call errorHandle("read_time_at_index ERROR: time_index arg larger than length of time dimension",-999)
    end if

    if (time_index .gt. 0) then
      strt(1) = time_index
    else
      strt(1) = int(pio_atm_file%numRecs,kind=pio_offset_kind)
    end if

    cnt(1)  = 1
    ierr = PIO_get_var(pio_atm_file%pioFileDesc,varid,strt,cnt,val_buf)
    call errorHandle('read_time_at_index: Error reading variable "time" in file '//trim(filename)//'.',ierr);
    val  = val_buf(1)
  end function read_time_at_index
!=====================================================================!
  ! Retrieve the dimension length for a file.
  function get_dimlen(filename,dimname) result(val)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dimname
    integer                      :: val

    type(pio_atm_file_t), pointer :: pio_atm_file
    integer                       :: dim_id, ierr
    logical                       :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    if (.not.found) call errorHandle("pio_inq_dimlen ERROR: File "//trim(filename)//" not found",-999)
    ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(dimname),dim_id)
    call errorHandle("pio_inq_dimlen ERROR: dimension "//trim(dimname)//" not found in file "//trim(filename)//".",ierr)
    ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,dim_id,val)

  end function get_dimlen

  function has_variable(filename,varname) result(has)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname

    logical                       :: has, found
    type(pio_atm_file_t), pointer :: pio_atm_file
    integer                       :: nvars, ivar, ierr
    character (len=256)           :: vname

    call register_file(filename,file_purpose_in)
    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)

    ierr = PIO_inquire(pio_atm_file%pioFileDesc,nVariables=nvars)
    call errorHandle("pio_inquire on file "//trim(filename)//" returned nonzero error code.",ierr)

    has = .false.
    do ivar=1,nvars
      ierr = pio_inquire_variable(pio_atm_file%pioFileDesc, ivar, name=vname)
      if (trim(vname) == trim(varname)) then
        has = .true.
        exit
      endif
    enddo

    call eam_pio_closefile(filename)
  end function has_variable
!=====================================================================!
  ! Write output to file based on type (int or real)
  ! --Note-- that any dimensionality could be written if it is flattened to 1D
  ! before calling a write routine.
  ! Mandatory inputs are:
  ! filename:     the netCDF filename to be written to.
  ! varname:      The shortname for the variable being written.
  ! var_data_ptr: An array of data that will be used to write the output.  NOTE:
  !               If PIO MPI ranks > 1, var_data_ptr should be the subset of the global array
  !               which includes only those degrees of freedom that have been
  !               assigned to this rank.  See set_dof above.
  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_float(filename, varname, buf, buf_size)
    use pionfput_mod, only: PIO_put_var   => put_var
    use piolib_mod, only: PIO_setframe
    use pio_types, only: PIO_max_var_dims
    use piodarray,  only: PIO_write_darray

    ! Dummy arguments
    character(len=*),    intent(in) :: filename       ! PIO filename
    character(len=*),    intent(in) :: varname
    integer(kind=c_int), intent(in) :: buf_size
    real(kind=c_float),  intent(in) :: buf(buf_size)

    ! Local variables

    type(pio_atm_file_t), pointer :: pio_atm_file
    type(hist_var_t), pointer     :: var
    integer                       :: ierr,jdim
    integer                       :: start(pio_max_var_dims), count(pio_max_var_dims)
    logical                       :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    if (var%has_t_dim) then
      ! Set the time index we are writing
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    endif

    if (var%is_partitioned) then
      call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    else
      if (var%has_t_dim) then
        do jdim=1,var%numdims
          if (var%dimid(jdim) .eq. time_dimid) then
            start (jdim) = int(max(1,pio_atm_file%numRecs))
            count (jdim) = 1
          else
            start (jdim) = 1
            count (jdim) = var%dimlen(jdim)
          endif
        enddo
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,start(:var%numdims),count(:var%numdims),buf)
      else
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,buf)
      endif
    endif

    call errorHandle( 'eam_grid_write_darray_float: Error writing variable '//trim(varname),ierr)
  end subroutine grid_write_darray_float
  subroutine grid_write_darray_double(filename, varname, buf, buf_size)
    use pionfput_mod, only: PIO_put_var   => put_var
    use pio_types, only: PIO_max_var_dims
    use piolib_mod, only: PIO_setframe
    use piodarray,  only: PIO_write_darray

    ! Dummy arguments
    character(len=*),    intent(in) :: filename       ! PIO filename
    character(len=*),    intent(in) :: varname
    integer(kind=c_int), intent(in) :: buf_size
    real(kind=c_double), intent(in) :: buf(buf_size)

    ! Local variables

    type(pio_atm_file_t), pointer :: pio_atm_file
    type(hist_var_t), pointer     :: var
    integer                       :: ierr,jdim
    integer                       :: start(pio_max_var_dims), count(pio_max_var_dims)
    logical                       :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    if (var%has_t_dim) then
      ! Set the time index we are writing
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    endif

    if (var%is_partitioned) then
      call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    else
      if (var%has_t_dim) then
        do jdim=1,var%numdims
          if (var%dimid(jdim) .eq. time_dimid) then
            start (jdim) = int(max(1,pio_atm_file%numRecs))
            count (jdim) = 1
          else
            start (jdim) = 1
            count (jdim) = var%dimlen(jdim)
          endif
        enddo
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,start(:var%numdims),count(:var%numdims),buf)
      else
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,buf)
      endif
    endif

    call errorHandle( 'eam_grid_write_darray_double: Error writing variable '//trim(varname),ierr)
  end subroutine grid_write_darray_double
  subroutine grid_write_darray_int(filename, varname, buf, buf_size)
    use pionfput_mod, only: PIO_put_var   => put_var
    use piolib_mod, only: PIO_setframe
    use pio_types, only: PIO_max_var_dims
    use piodarray,  only: PIO_write_darray

    ! Dummy arguments
    character(len=*),    intent(in) :: filename       ! PIO filename
    character(len=*),    intent(in) :: varname
    integer(kind=c_int), intent(in) :: buf_size
    integer(kind=c_int), intent(in) :: buf(buf_size)

    ! Local variables

    type(pio_atm_file_t), pointer :: pio_atm_file
    type(hist_var_t), pointer     :: var
    integer                       :: ierr,jdim
    integer                       :: start(pio_max_var_dims), count(pio_max_var_dims)
    logical                       :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    if (var%has_t_dim) then
      ! Set the time index we are writing
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(max(1,pio_atm_file%numRecs),kind=pio_offset_kind))
    endif

    if (var%is_partitioned) then
      call pio_write_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    else
      if (var%has_t_dim) then
        do jdim=1,var%numdims
          if (var%dimid(jdim) .eq. time_dimid) then
            start (jdim) = int(max(1,pio_atm_file%numRecs))
            count (jdim) = 1
          else
            start (jdim) = 1
            count (jdim) = var%dimlen(jdim)
          endif
        enddo
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,start(:var%numdims),count(:var%numdims),buf)
      else
        ierr = pio_put_var(pio_atm_file%pioFileDesc,var%piovar,buf)
      endif
    endif

    call errorHandle( 'eam_grid_write_darray_int: Error writing variable '//trim(varname),ierr)
  end subroutine grid_write_darray_int
!=====================================================================!
  ! Read output from file based on type (int or real)
  ! --Note-- that any dimensionality could be read if it is flattened to 1D
  ! before calling a write routine.
  ! Mandatory inputs are:
  ! filename:     the netCDF filename to be read from.
  ! varname:      The shortname for the variable being read.i
  ! var_data_ptr: An c_ptr of data that will be used to read the input.  NOTE:
  !               If PIO MPI ranks > 1, var_data_ptr should be the subset of the global array
  !               which includes only those degrees of freedom that have been
  !               assigned to this rank.  See set_dof above.
  ! time_index:   The 1-based time index for this variable (0 is ignored).
  !---------------------------------------------------------------------------
  !
  !  grid_read_darray_1d: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_read_darray_double(filename, varname, buf, buf_size, time_index)
    use piolib_mod, only: PIO_setframe
    use piodarray,  only: PIO_read_darray

    ! Dummy arguments
    character(len=*),     intent(in) :: filename       ! PIO filename
    character(len=*),     intent(in) :: varname
    integer (kind=c_int), intent(in) :: buf_size
    real(kind=c_double),  intent(out) :: buf(buf_size)
    integer, intent(in)          :: time_index

    ! Local variables
    type(pio_atm_file_t),pointer       :: pio_atm_file
    type(hist_var_t), pointer          :: var
    integer                            :: ierr, var_size
    logical                            :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    ! Set the timesnap we are reading
    if (time_index .gt. 0) then
      ! The user has set a valid time index to read from
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(time_index,kind=pio_offset_kind))
    else
      ! Otherwise default to the last time_index in the file
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(pio_atm_file%numRecs,kind=pio_offset_kind))
    end if
    
    ! We don't want the extent along the 'time' dimension
    var_size = SIZE(var%compdof)

    ! Now we know the exact size of the array, and can shape the f90 pointer
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    call errorHandle( 'eam_grid_read_darray_double: Error reading variable '//trim(varname),ierr)
  end subroutine grid_read_darray_double
  subroutine grid_read_darray_float(filename, varname, buf, buf_size, time_index)
    use piolib_mod, only: PIO_setframe
    use piodarray,  only: PIO_read_darray

    ! Dummy arguments
    character(len=*),     intent(in) :: filename       ! PIO filename
    character(len=*),     intent(in) :: varname
    integer (kind=c_int), intent(in) :: buf_size
    real(kind=c_float),  intent(out) :: buf(buf_size)
    integer, intent(in)          :: time_index

    ! Local variables
    type(pio_atm_file_t),pointer       :: pio_atm_file
    type(hist_var_t), pointer          :: var
    integer                            :: ierr, var_size
    logical                            :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    ! Set the timesnap we are reading
    if (time_index .gt. 0) then
      ! The user has set a valid time index to read from
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(time_index,kind=pio_offset_kind))
    else
      ! Otherwise default to the last time_index in the file
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(pio_atm_file%numRecs,kind=pio_offset_kind))
    end if
    
    ! We don't want the extent along the 'time' dimension
    var_size = SIZE(var%compdof)

    ! Now we know the exact size of the array, and can shape the f90 pointer
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    call errorHandle( 'eam_grid_read_darray_float: Error reading variable '//trim(varname),ierr)
  end subroutine grid_read_darray_float
  subroutine grid_read_darray_int(filename, varname, buf, buf_size, time_index)
    use piolib_mod, only: PIO_setframe
    use piodarray,  only: PIO_read_darray

    ! Dummy arguments
    character(len=*),     intent(in) :: filename       ! PIO filename
    character(len=*),     intent(in) :: varname
    integer (kind=c_int), intent(in) :: buf_size
    integer (kind=c_int), intent(out) :: buf(buf_size)
    integer, intent(in)          :: time_index

    ! Local variables
    type(pio_atm_file_t),pointer       :: pio_atm_file
    type(hist_var_t), pointer          :: var
    integer                            :: ierr, var_size
    logical                            :: found

    call lookup_pio_atm_file(trim(filename),pio_atm_file,found)
    call get_var(pio_atm_file,varname,var)

    ! Set the timesnap we are reading
    if (time_index .gt. 0) then
      ! The user has set a valid time index to read from
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(time_index,kind=pio_offset_kind))
    else
      ! Otherwise default to the last time_index in the file
      call PIO_setframe(pio_atm_file%pioFileDesc,var%piovar,int(pio_atm_file%numRecs,kind=pio_offset_kind))
    end if
    
    ! We don't want the extent along the 'time' dimension
    var_size = SIZE(var%compdof)

    ! Now we know the exact size of the array, and can shape the f90 pointer
    call pio_read_darray(pio_atm_file%pioFileDesc, var%piovar, var%iodesc, buf, ierr)
    call errorHandle( 'eam_grid_read_darray_int: Error reading variable '//trim(varname),ierr)
  end subroutine grid_read_darray_int
!=====================================================================!
  subroutine convert_int_2_str(int_in,str_out)
    integer, intent(in)           :: int_in
    character(len=*), intent(out) :: str_out

    character(len=8) :: fmt_str

    if (int_in < 0) then
      fmt_str = "(A1,"
    else
      fmt_str = "("
    end if

    if (abs(int_in)<10) then
      fmt_str = trim(fmt_str)//"I1)"
    elseif (abs(int_in)<1e2) then
      fmt_str = trim(fmt_str)//"I2)"
    elseif (abs(int_in)<1e3) then
      fmt_str = trim(fmt_str)//"I3)"
    elseif (abs(int_in)<1e4) then
      fmt_str = trim(fmt_str)//"I4)"
    elseif (abs(int_in)<1e5) then
      fmt_str = trim(fmt_str)//"I5)"
    elseif (abs(int_in)<1e6) then
      fmt_str = trim(fmt_str)//"I6)"
    elseif (abs(int_in)<1e7) then
      fmt_str = trim(fmt_str)//"I7)"
    elseif (abs(int_in)<1e8) then
      fmt_str = trim(fmt_str)//"I8)"
    elseif (abs(int_in)<1e9) then
      fmt_str = trim(fmt_str)//"I9)"
    endif
    
    if (int_in < 0) then
      write(str_out,fmt_str) "n", int_in
    else
      write(str_out,fmt_str) int_in
    end if 

  end subroutine convert_int_2_str

  subroutine get_coord (filename,shortname,hist_coord)
    character(len=256), intent(in)           :: filename
    character(len=256), intent(in)           :: shortname
    type(hist_coord_t), pointer, intent(out) :: hist_coord

    type(pio_atm_file_t), pointer    :: pio_atm_file
    type(hist_coord_list_t), pointer :: curr, prev
    logical                          :: file_found, dim_found

    hist_coord => NULL()

    call lookup_pio_atm_file(trim(filename),pio_atm_file,file_found)
    if (.not. file_found ) then
      call errorHandle("PIO ERROR: error retrieving dimension "//trim(shortname)//" in file "//trim(filename)//".\n PIO file not found or not open.",-999)
    endif

    curr => pio_atm_file%coord_list_top
    do while (associated(curr))
      if (associated(curr%coord)) then
        if(trim(curr%coord%name)==trim(shortname)) then
          dim_found = .true.
          exit
        endif
      end if
      prev => curr
      curr => prev%next
    end do

    if (.not. dim_found ) then
      call errorHandle("PIO ERROR: error retrieving dimension "//trim(shortname)//" from file "//trim(filename)//". Dimension not found.",-999)
    endif

    hist_coord => curr%coord
  end subroutine get_coord
!=====================================================================!
end module scream_scorpio_interface
