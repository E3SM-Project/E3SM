module scream_scorpio_interface

!==============================================================================!
! This module handles the Fortran interface to the PIO library of input/output
! subroutines.  The essential set of steps to enable and use PIO for creating
! output are as follows:
!                                  OUTPUT
! Initialization:
! 1) Gather the pio_subsystem and pio_iotype information for the EAM component
! as assigned by the component coupler.
!    This is accomplished during eam_init_pio by calling
!    'eam_init_pio_subsystem'
! 2) For each output file "create" a file in PIO and record the unique file
! descriptor. 
!    This is accomplished during the eam_init_pio by calling
!    'eam_pio_createfile'
! 3) For each output file define the "header" information, which is essentially
! a set of metadata strings that describe the output file.
!    This is accomplished during eam_init_pio by calling 'eam_pio_createHeader'
! 4) Define all of the dimensions that variables in this file will be defined
! on.  Examples would time, lat, lon, vertical coordinate, # of consituents,
! etc.
!    This is accomplished during 'h_define' by repeated calls to 'PIO_def_dim'
! 5) Define all of the variables that will be written to this file.  This
! includes any dimension that should also be defined as a variable, such as lat,
! lon, time.
!    This is accomplished during 'h_define' by repeated calls to 'PIO_def_var'
! 6) Determine the unique PIO identifiers for each domain decomposition.  This
! is essentially a decomposition of what will be written to the file, multiple
! variables can have the same decomposition.  Every arrangement of dimensions
! requires a pio decomposition.
!    This is accomplished during 'h_define' by repeated calls to
!    'PIO_initdecomp'
! 7) Close the PIO file definition step.  In other words, tell PIO that all of
! the dimensions, variables and decompositions associated with this output have
! been defined and no new ones will be added.
!    This is accomplished during eam_init_pio by calling 'PIO_enddef'
!
! Writing Output:
! 
! Finalization: 
!==============================================================================!

 
  use shr_pio_mod,  only: shr_pio_getrearranger
  use shr_kind_mod,   only: rtype=>shr_kind_r8, rtype_short=>shr_kind_r4  ! Need to change this to use scream real types
  use pio_mods, only: io_desc_t, iosystem_desc_t, file_desc_t, var_desc_t, pio_global, &
                      pio_init, var_desc_t, pio_unlimited, pio_int, pio_def_dim,& 
                      pio_def_var, pio_enddef, pio_noerr, pio_closefile, pio_inq_dimid, &
                      pio_real, pio_double, pio_initdecomp, pio_inq_dimlen, pio_setframe, &
                      pio_offset_kind, pio_write_darray, pio_syncfile, pio_put_var

  implicit none
  save
                     
  public :: eam_init_pio,   &  ! Get pio subsystem info from main code
            eam_h_define,   &      ! Create a new NetCDF file for PIO writing
            eam_h_finalize, &
            eam_h_write,    &
            register_variable, &
            register_dimension 
 
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
  type, public :: recur_test
    integer :: A
    type(recur_test), pointer :: next => NULL()
  end type recur_test
!----------------------------------------------------------------------
  type, public :: hist_coord_t
    character(len=max_hcoordname_len) :: name = ''  ! coordinate name
    integer                  :: dimsize = 0       ! size of dimension
    integer                  :: dimid             ! Unique PIO Id for this dimension
    character(len=max_chars) :: long_name = ''    ! 'long_name' attribute
    character(len=max_chars) :: units = ''        ! 'units' attribute
    character(len=max_chars) :: bounds_name = ''  ! 'bounds' attribute (& name of bounds variable)
    character(len=max_chars) :: standard_name = ''! 'standard_name' attribute
    character(len=4)         :: positive = ''     ! 'positive' attribute ('up' or 'down')
    integer,  pointer        :: integer_values(:) => null() ! dim values if integer
    real(rtype), pointer     :: real_values(:) => null() ! dim values if real
    real(rtype), pointer     :: bounds(:,:) => null() ! dim bounds
    logical                  :: integer_dim       ! .true. iff dim has integral values
    logical                  :: vertical_coord    ! .true. iff dim is vertical
  end type hist_coord_t
!----------------------------------------------------------------------
  type, public :: hist_var_t
    character(len=max_hvarname_len) :: name   ! coordinate name
    character(len=max_chars) :: long_name     ! 'long_name' attribute
    character(len=max_chars) :: units         ! 'units' attribute
    type(var_desc_t) :: piovar                    ! netCDF variable ID
    integer          :: dtype                     ! data type
    integer          :: numdims                   ! Number of dimensions in out field
    type(io_desc_t), pointer  :: iodesc           ! PIO decomp associated with this variable
    integer, allocatable :: compdof(:)            ! Global locations in output array for this process
    integer, allocatable :: dimid(:)              ! array of PIO dimension id's for this variable
    integer, allocatable :: dimlen(:)             ! array of PIO dimension lengths for this variable
  end type hist_var_t
!----------------------------------------------------------------------
  ! The iodesc_list allows us to cache existing PIO decompositions
  ! The tag needs the dim lengths, the dtype and map id (+ optional permutation)
  ! Define a recursive structure because we do not know ahead of time how many
  ! decompositions will be require
  integer, parameter      :: tag_len           = 48
  type iodesc_list
    character(tag_len)          :: tag
    type(io_desc_t),    pointer :: iodesc => NULL()
    type(iodesc_list),  pointer :: next => NULL()
  end type iodesc_list
  ! Define the first iodesc_list 
  type(iodesc_list), target :: iodesc_list_top
!----------------------------------------------------------------------
  type hist_var_list
    type(hist_var_t),    pointer :: var => NULL()
    type(hist_var_list), pointer :: next => NULL()
  end type hist_var_list
!----------------------------------------------------------------------
  type, public :: pio_atm_output
        !> @brief Output filename.
        character(len=max_chars) :: filename

        !> @brief Contains data identifying the file.
        type(file_desc_t)     :: pioFileDesc

        !> @brief Number of output dimensions, and counter to track them during
        !  registration
        integer               :: numDim
        integer               :: DimCounter = 0

        !> @brief Number of output variable and counter to track them during
        !  registration
        integer               :: numVar
        integer               :: VarCounter = 0

        !> @brief Number of history records on this file
        integer               :: numRecs

        !> @brief Coordinate Dimensions Array
        type(hist_coord_t), allocatable :: dimensions(:)

        !> @brief Variable Array
        type(hist_var_t), allocatable :: variables(:)

        !> @brief Recursive list of variables
        type(hist_var_list), pointer :: var_list => NULL()

  end type pio_atm_output
!----------------------------------------------------------------------
  
  type(pio_atm_output), target, allocatable :: atm_output_files(:)

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
  subroutine eam_init_pio(mpicom,atm_id,numdim,numvar)

    integer, intent(in) :: mpicom  ! MPI communication group for ATM
    integer, intent(in) :: atm_id  ! Unique identifier assigned by component coupler
    integer, intent(in) :: numdim  ! Total number of possible dimensions for outputs
    integer, intent(in) :: numvar  ! Total number of variables for output

    integer :: numFiles = 2  ! TODO: When switch to more general is complete this should be set by AD
    type(pio_atm_output), pointer :: current_atm_file => null()
    integer :: ierr
    
    type(hist_var_list), pointer :: curr => null()

    ! Gather the pio_subsystem information for the atmosphere component.
    call eam_init_pio_subsystem(mpicom,atm_id)

    ! Allocate one atm_output_files (for testing), TODO: make this more general
    ! to allow for multiple output files.
    allocate( atm_output_files(numFiles) )

    ! Allocate the dimension and variable arrays for the output file.
    current_atm_file => atm_output_files(1)
    current_atm_file%numDim = numdim
    current_atm_file%numVar = numvar
    current_atm_file%numRecs = 0
    current_atm_file%var_list => NULL()
    allocate( current_atm_file%dimensions(current_atm_file%numdim), &
              current_atm_file%variables(current_atm_file%numvar) )

    ! Set the filename TODO: Allow for mutliple files
    current_atm_file%filename = "example_pio_structured.nc"
    write(*,*) 'ASD : ', trim(current_atm_file%filename), pio_subsystem, pio_iotype 
    ! Create the new file via PIO
    call eam_pio_createfile(current_atm_file%pioFileDesc,trim(current_atm_file%filename))
    ! Create the header for the new PIO file
    call eam_pio_createHeader(current_atm_file%pioFileDesc)

    ! Register all dimensions with the output file
    call register_dimension(current_atm_file,"x","horizontal distance",10)
    call register_dimension(current_atm_file,"y","vertical distance",3)
    call register_dimension(current_atm_file,"time","time",0)

    ! Register all variables with the output file
    call register_variable(current_atm_file,"time","time",1,(/ "time" /), PIO_double,"t")
    call register_variable(current_atm_file,"x","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    call register_variable(current_atm_file,"y","answer to space and time",1,(/ "y" /), PIO_double,"y-real")

    call register_variable(current_atm_file,"real_foo","answer to space and time",2,(/ "x", "time" /), PIO_double,"xt-real")
    call register_variable(current_atm_file,"bar2","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    call register_variable(current_atm_file,"foo","answer to space and time",2,(/ "x", "time" /),PIO_int,"xt-int")
    call register_variable(current_atm_file,"bar","answer to space and time",3,(/ "x", "y", "time" /), PIO_double,"xyt-real")

    ! Finish the "definition" phase of the PIO file.  This is an essential step
    ! for the netCDF file to be ready for variables to be written to it.
    ierr = PIO_enddef(current_atm_file%pioFileDesc)
    call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),current_atm_file%pioFileDesc,ierr)
   

    ! SECOND PIO OUTPUT FILE 
    ! Allocate the dimension and variable arrays for the output file.
    current_atm_file => atm_output_files(2)
    current_atm_file%numDim = numdim
    current_atm_file%numVar = numvar
    current_atm_file%numRecs = 0
    allocate( current_atm_file%dimensions(current_atm_file%numdim), &
              current_atm_file%variables(current_atm_file%numvar) )
    current_atm_file%filename = "example_pio_structured_v2.nc"
    ! Create the new file via PIO
    call eam_pio_createfile(current_atm_file%pioFileDesc,trim(current_atm_file%filename))
    ! Register all dimensions with the output file
    call register_dimension(current_atm_file,"x","horizontal distance",10)
    call register_dimension(current_atm_file,"y","vertical distance",3)
    call register_dimension(current_atm_file,"time","time",0)
    ! Register all variables with the output file
    call register_variable(current_atm_file,"time","time",1,(/ "time" /), PIO_double,"t")
    call register_variable(current_atm_file,"x","answer to space and time",1,(/ "x" /), PIO_double,"x-real")
    call register_variable(current_atm_file,"y","answer to space and time",1,(/ "y" /), PIO_double,"y-real")
    call register_variable(current_atm_file,"bar","answer to space and time",3,(/ "x", "y", "time" /), PIO_double,"xyt-real")
    call register_variable(current_atm_file,"foo","answer to space and time",3,(/ "x", "y", "time" /), PIO_double,"xyt-real")
    call register_variable(current_atm_file,"bar_flip","answer to space and time",3,(/ "y", "x", "time" /), PIO_double,"yxt-real")
    call register_variable(current_atm_file,"foo_flip","answer to space and time",3,(/ "y", "x", "time" /), PIO_double,"yxt-real")
    ierr = PIO_enddef(current_atm_file%pioFileDesc)
    call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(current_atm_file%filename),current_atm_file%pioFileDesc,ierr)

    do ierr = 1,2
    current_atm_file => atm_output_files(ierr)
    write(*,*) "File: ", trim(current_atm_file%filename)
    curr => current_atm_file%var_list
    do while (associated(curr))
      curr => curr%next
      write(*,*) "File: ", trim(current_atm_file%filename), ", Var: ", trim(curr%var%name)
    end do
    end do
  end subroutine eam_init_pio
!=====================================================================!
  subroutine register_dimension(pio_atm_file,shortname,longname,length)
    type(pio_atm_output), pointer, intent(inout) :: pio_atm_file
    character(len=*), intent(in)        :: shortname,longname
    integer, intent(in)                 :: length

    type(hist_coord_t),pointer          :: hist_coord
    integer :: ierr
  
    pio_atm_file%dimcounter = pio_atm_file%dimcounter + 1
    if (pio_atm_file%dimcounter.gt.pio_atm_file%numdim) call errorHandle("EAM_PIO ERROR: Attempted to register more dimensions than originally declared, "//trim(shortname),pio_atm_file%pioFileDesc,-999)
    hist_coord => pio_atm_file%dimensions(pio_atm_file%dimcounter)
    hist_coord%name      = trim(shortname)
    hist_coord%long_name = trim(longname)
    hist_coord%dimsize   = length
    if (length.eq.0) then
      ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), pio_unlimited , hist_coord%dimid)
    else
      ierr = PIO_def_dim(pio_atm_file%pioFileDesc, trim(shortname), length , hist_coord%dimid)
    end if
    call errorHandle("PIO ERROR: could not define dimension "//trim(shortname),pio_atm_file%pioFileDesc,ierr)
    
    return
  end subroutine register_dimension
!=====================================================================!
  subroutine register_variable(pio_atm_file,shortname,longname,length,dims,dtype,tag)
    type(pio_atm_output),pointer, intent(inout) :: pio_atm_file
    character(len=*), intent(in) :: shortname,longname
    integer, intent(in)          :: length
    character(len=*), intent(in) :: dims(length)
    integer, intent(in)          :: dtype
    character(len=*), intent(in) :: tag     ! Unique tag for this variables decomposition type

    ! Local variables
    integer :: loc_len
    logical                      :: found
    type(io_desc_t), pointer     :: iodesc_loc
    type(hist_var_t), pointer    :: hist_var
    integer :: dim_ii
    integer :: ierr
    integer, allocatable :: dimlen(:)
    integer :: total_dimlen, my_dof_len, extra_procs
    integer :: ii, istart, istop
    logical :: has_t_dim  ! Logical to flag whether this variable has a time-dimension.  This is important for the decomposition step.

    type(hist_var_list), pointer :: curr => null(), prev => null()
 
    write(*,*) '(', pio_myrank, ')',' ASD Reg. Variable: ', trim(shortname), ' for file: ', trim(pio_atm_file%filename)
 
    has_t_dim = .false.
    pio_atm_file%varcounter = pio_atm_file%varcounter + 1
    if (pio_atm_file%varcounter.gt.pio_atm_file%numvar) call errorHandle("EAM_PIO ERROR: Attempted to register more variables than originally declared, "//trim(shortname),pio_atm_file%pioFileDesc,-999)
    hist_var => pio_atm_file%variables(pio_atm_file%varcounter)
    hist_var%name      = trim(shortname)
    hist_var%long_name = trim(longname)
    hist_var%numdims   = length
    hist_var%dtype    = dtype
    allocate(hist_var%dimid(length),hist_var%dimlen(length))
    do dim_ii = 1,length
      ierr = pio_inq_dimid(pio_atm_file%pioFileDesc,trim(dims(dim_ii)),hist_var%dimid(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to find dimension id for "//trim(dims(dim_ii)),pio_atm_file%pioFileDesc,ierr)
      ierr = pio_inq_dimlen(pio_atm_file%pioFileDesc,hist_var%dimid(dim_ii),hist_var%dimlen(dim_ii))
      call errorHandle("EAM_PIO ERROR: Unable to determine length for dimension "//trim(dims(dim_ii)),pio_atm_file%pioFileDesc,ierr)
      if (hist_var%dimlen(dim_ii).eq.0) has_t_dim = .true.
    end do

    ! Distribute responsibility for writing cores over all PIO ranks
    ! i.e. compute the degrees of freedom that this rank will contribute to PIO
    if (has_t_dim) then
      call get_compdof(length-1,hist_var%dimlen(:length-1),my_dof_len,istart,istop)
    else
      call get_compdof(length,hist_var%dimlen(:length),my_dof_len,istart,istop)
    end if
    allocate( hist_var%compdof(my_dof_len) )
    hist_var%compdof(:my_dof_len) = (/ (ii, ii=istart,istop, 1) /)

    ! Register Variable with PIO
    ierr = PIO_def_var(pio_atm_file%pioFileDesc, trim(shortname), hist_var%dtype, hist_var%dimid(:length), hist_var%pioVar)
    call errorHandle("PIO ERROR: could not define variable "//trim(shortname),pio_atm_file%pioFileDesc,ierr)

    ! Gather the pio decomposition for this variable
    if (has_t_dim) then
      loc_len = max(1,length-1)
      call get_decomp(tag,hist_var%dtype,hist_var%dimlen(:loc_len),hist_var%compdof,hist_var%iodesc)
    else
      call get_decomp(tag,hist_var%dtype,hist_var%dimlen,hist_var%compdof,hist_var%iodesc)
    end if

!    write(*,*) '(', pio_myrank, ')',' ASD Reg. Variable - Rec: ', trim(shortname), ' for file: ', trim(pio_atm_file%filename)
!    curr => pio_atm_file%var_list
!    if (associated(curr)) then
!    do while (associated(curr))
!      prev => curr
!      curr => prev%next
!    end do
!    allocate(curr)
!    curr%var => hist_var
!    else
!      allocate(pio_atm_file%var_list)
!      pio_atm_file%var_list%var => hist_var 
!    end if
!    write(*,*) '(', pio_myrank, ')',' ASD Reg. Variable - Rec DONE: ', trim(shortname), trim(curr%var%name), ' for file: ', trim(pio_atm_file%filename)

    return
  end subroutine register_variable
!=====================================================================!
  subroutine eam_h_write()

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
    
  end subroutine eam_h_write
!=====================================================================!
  ! MAY BE DEPRECATED
  subroutine eam_h_define(numdim,dimlen,dimnames)

   integer, intent(in) :: numdim
   integer, intent(in) :: dimlen(numdim)
   character(len=*), intent(in) :: dimnames(numdim)

   character(len=100) :: fname
   integer :: ii, ierr

   integer :: numvar = 1
   character(len=10) :: varnames(1)
   integer :: vartype(1)

!   varnames(1) = "foo"
!   vartype(1) = PIO_int
!
!   fname = "eam_pio_example.nc"
!
!   ! Create the file
!   call eam_pio_createfile(pioFile,trim(fname)) ! TODO set up File and fname inputs
!   ! Create netCDF Header info (like caseid, title, etc.)
!   call eam_pio_createHeader(pioFile)
!   ! Define dimensions
!   allocate(pioDimId(numdim))
!   do ii = 1,numdim
!     if (dimlen(ii).eq.0) then
!       ierr = PIO_def_dim(pioFile, trim(dimnames(ii)), pio_unlimited , pioDimId(ii))
!     else
!       ierr = PIO_def_dim(pioFile, trim(dimnames(ii)), dimlen(ii) , pioDimId(ii))
!     end if
!     call errorHandle("PIO ERROR: could not define dimension "//trim(dimnames(ii)),ierr)
!   end do
!
!   do ii = 1,numvar
!     ierr = PIO_def_var(pioFile, trim(varnames(ii)), vartype(ii), (/pioDimId/), pioVar)
!     call errorHandle("PIO ERROR: could not define variable "//trim(varnames(ii)),ierr)
!   end do
!
!   ierr = PIO_enddef(pioFile)
!   call errorHandle("PIO ERROR: issue arose with PIO_enddef for file"//trim(fname),ierr)

   ! TODO define vars from field manager
   ! TODO define optional dimensions for nonstandard variables.  Alternatively,
   ! define dimensions based on field manager (second option is probably better).
   ! Define Grid attribute
   ! TODO Define grid attribute routine ala cam_grid_write_attr


  end subroutine eam_h_define
!=====================================================================!
  subroutine eam_h_finalize()
!    type(pio_atm_output),pointer :: pio_atm_file

!    pio_atm_file => atm_output_files(1)
!    call PIO_syncfile(pio_atm_file%pioFileDesc)
!    call eam_pio_finalize()

  end subroutine eam_h_finalize
!=====================================================================!
  subroutine eam_pio_createHeader(File)

    use pio_mods, only : PIO_put_att

    type(file_desc_t), intent(in) :: File             ! Pio file Handle
    integer :: retval

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
  subroutine eam_pio_createfile(File,fname)
    use pio_mods, only:  pio_createfile, pio_clobber

    type(file_desc_t), intent(inout) :: File             ! Pio file Handle
    character(len=*),  intent(in)    :: fname
    !--
    integer :: retval
    integer                                   :: mode

    mode = pio_clobber ! Set to CLOBBER for now, TODO: fix to allow for optional mode type like in CAM
    retval = pio_createfile(pio_subsystem,File,pio_iotype,fname,mode) 

  end subroutine eam_pio_createfile
!=====================================================================!
  subroutine eam_pio_finalize()

    use pio_mods, only: pio_finalize

    integer :: ierr

    call PIO_finalize(pio_subsystem, ierr)

  end subroutine eam_pio_finalize
!=====================================================================!
    subroutine errorHandle(errMsg, pioFile, retVal)

        implicit none

        type(file_desc_t), intent(in)    :: pioFile
        character(len=*),  intent(in)    :: errMsg
        integer,           intent(in)    :: retVal

        if (retVal .ne. PIO_NOERR) then
            write(*,*) retVal,errMsg
            call PIO_closefile(pioFile)
            call mpi_abort(pio_mpicom,0,retVal)
        end if

    end subroutine errorHandle
!=====================================================================!
  subroutine get_compdof(length,dimlen,dof_len,istart,istop)

    integer, intent(in)  :: length
    integer, intent(in)  :: dimlen(length)
    integer, intent(out) :: dof_len, istart, istop

    integer :: extra_procs, total_dimlen

    ! Get the total number of array elements for this output
    total_dimlen = product(dimlen)
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
  subroutine get_decomp(tag,dtype,dimlen,compdof,iodesc)

    character(len=*)            :: tag
    integer, intent(in)         :: dtype
    integer, intent(in)         :: dimlen(:)
    integer, intent(in)         :: compdof(:)
    type(io_desc_t),    pointer :: iodesc

    logical                      :: found
    type(iodesc_list),pointer    :: curr, prev
    type(io_desc_t), pointer     :: iodesc_loc
    integer                      :: loc_len
    
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
      loc_len = size(dimlen)
      if ( loc_len.eq.1 .and. dimlen(loc_len).eq.0 ) then
        allocate(curr%iodesc)
      else
        call pio_initdecomp(pio_subsystem, dtype, dimlen, compdof, curr%iodesc, rearr=pio_rearranger)
      end if
      iodesc => curr%iodesc
    end if 

  end subroutine get_decomp
!=====================================================================!
  subroutine get_var(pio_file,varname,var)

    type(pio_atm_output), pointer :: pio_file
    character(len=*)              :: varname
    type(hist_var_t), pointer     :: var
    
    integer :: ii
    type(hist_var_list), pointer :: curr

!    curr => pio_file%var_list
!    do while (associated(curr))
!      var => curr%var
!      write(*,*) "get_var: ", trim(varname), ", ", trim(var%name)
!      if (trim(varname) == trim(var%name)) return
!      curr => curr%next
!    end do

    do ii = 1,pio_file%varcounter
      var => pio_file%variables(ii)
      if (trim(varname)==trim(var%name)) return
    end do
    ! If we got this far we didn't find the variable
    call errorHandle("PIO ERROR: unable to find variable: "//trim(varname)//" in file: "//trim(pio_file%filename),pio_file%pioFileDesc,999)

  end subroutine get_var
!=====================================================================!
  subroutine get_pio_atm_file(filename,pio_atm_file)

    character(len=*)              :: filename
    type(pio_atm_output), pointer :: pio_atm_file

    integer :: ii

    do ii = 1,size(atm_output_files)
      pio_atm_file => atm_output_files(ii)
      if (trim(filename)==trim(pio_atm_file%filename)) return
    end do  
    ! If we got this far we didn't find the file
    call errorHandle("PIO ERROR: unable to find file: "//trim(filename),pio_atm_file%pioFileDesc,999)

  end subroutine get_pio_atm_file
!=====================================================================!
  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_int(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_INT

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
    call errorHandle( 'cam_grid_write_darray_1d_int: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_1d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_int(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_INT

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
    call errorHandle( 'cam_grid_write_darray_2d_int: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_2d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_int(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_INT

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
    call errorHandle( 'cam_grid_write_darray_3d_int: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_3d_int

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_double(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE

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
    call errorHandle( 'cam_grid_write_darray_1d_double: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_1d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_double(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE

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
    call errorHandle( 'cam_grid_write_darray_2d_double: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_2d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_double(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE

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
    call errorHandle( 'cam_grid_write_darray_3d_double: Error writing variable',pio_atm_file%pioFileDesc,ierr)

  end subroutine grid_write_darray_3d_double

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_1d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_1d_real(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_REAL

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
    call errorHandle( 'cam_grid_write_darray_1d_real: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_1d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_2d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_2d_real(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_REAL

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
    call errorHandle( 'cam_grid_write_darray_2d_real: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_2d_real

  !---------------------------------------------------------------------------
  !
  !  grid_write_darray_3d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine grid_write_darray_3d_real(filename, hbuf, varname)
    use pio,           only: file_desc_t
    use pio,           only: pio_write_darray, PIO_REAL

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
    call errorHandle( 'cam_grid_write_darray_3d_real: Error writing variable',pio_atm_file%pioFileDesc,ierr)
  end subroutine grid_write_darray_3d_real

end module scream_scorpio_interface
