module cam_grid_support
  use shr_kind_mod,        only: r8=>shr_kind_r8, r4=>shr_kind_r4, max_chars=>shr_kind_cl
  use shr_kind_mod,        only: i8=>shr_kind_i8, i4=>shr_kind_i4
  use shr_sys_mod,         only: shr_sys_flush
  use pio,                 only: iMap=>PIO_OFFSET_KIND, var_desc_t
  use cam_abortutils,      only: endrun
  use cam_logfile,         only: iulog
  use spmd_utils,          only: masterproc
  use cam_pio_utils,       only: cam_pio_handle_error
  use cam_map_utils,       only: cam_filemap_t
!!XXgoldyXX:
use spmd_utils, only: iam, npes, mpicom, mpi_integer,  mpi_sum, mpi_max, mpi_min
use shr_sys_mod, only: shr_sys_flush
use cam_map_utils, only: goldy_debug
!!XXgoldyXX:

  implicit none
  private

  public iMap

  integer, parameter, public :: max_hcoordname_len   = 16

  !---------------------------------------------------------------------------
  !
  !  horiz_coord_t: Information for horizontal dimension attributes
  !
  !---------------------------------------------------------------------------
  type, public :: horiz_coord_t
    private
    character(len=max_hcoordname_len) :: name = ''  ! coordinate name
    character(len=max_hcoordname_len) :: dimname = ''  ! dimension name
         ! NB: If dimname is blank, it is assumed to be name
    integer                   :: dimsize = 0       ! global size of dimension
    character(len=max_chars)  :: long_name = ''    ! 'long_name' attribute
    character(len=max_chars)  :: units = ''        ! 'units' attribute
    real(r8),         pointer :: values(:) => NULL() ! dim values (local if map)
    integer(iMap),    pointer :: map(:) => NULL()  ! map (dof) for dist. coord
    logical                   :: latitude          ! .false. means longitude
    type(var_desc_t), pointer :: vardesc => NULL() ! If we are to write coord
  contains
    procedure                 :: get_coord_len  => horiz_coord_len
    procedure                 :: num_elem       => horiz_coord_num_elem
    procedure                 :: global_size    => horiz_coord_find_size
    procedure                 :: get_coord_name => horiz_coord_name
    procedure                 :: get_dim_name   => horiz_coord_dim_name
    procedure                 :: get_long_name  => horiz_coord_long_name
    procedure                 :: get_units      => horiz_coord_units
    procedure                 :: write_attr     => write_horiz_coord_attr
    procedure                 :: write_var      => write_horiz_coord_var
  end type horiz_coord_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_t: Auxiliary quantity for a CAM grid
  !
  !---------------------------------------------------------------------------
  type, abstract :: cam_grid_attribute_t
    character(len=max_hcoordname_len)    :: name = ''      ! attribute name
    character(len=max_chars)             :: long_name = '' ! attribute long_name
    type(var_desc_t), pointer            :: vardesc => NULL()
! We aren't going to use this until we sort out PGI issues
    class(cam_grid_attribute_t), pointer :: next => NULL()
  contains
    procedure                                :: cam_grid_attr_init
    procedure(write_cam_grid_attr), deferred :: write_attr
    procedure(write_cam_grid_attr), deferred :: write_val
    procedure(print_attr_spec),     deferred :: print_attr
    procedure                                :: print_attr_base
  end type cam_grid_attribute_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_0d_int_t: Global integral attribute
  !
  !---------------------------------------------------------------------------
  type, extends(cam_grid_attribute_t) :: cam_grid_attribute_0d_int_t
    integer                             :: ival
  contains
    procedure :: cam_grid_attr_init_0d_int
    procedure :: write_attr => write_cam_grid_attr_0d_int
    procedure :: write_val  => write_cam_grid_val_0d_int
    procedure :: print_attr => print_attr_0d_int
  end type cam_grid_attribute_0d_int_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_0d_char_t: Global string attribute
  !
  !---------------------------------------------------------------------------
  type, extends(cam_grid_attribute_t) :: cam_grid_attribute_0d_char_t
    character(len=max_chars)            :: val
  contains
    procedure :: cam_grid_attr_init_0d_char
    procedure :: write_attr => write_cam_grid_attr_0d_char
    procedure :: write_val  => write_cam_grid_val_0d_char
    procedure :: print_attr => print_attr_0d_char
  end type cam_grid_attribute_0d_char_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_1d_int_t: 1-d integer attribute
  !
  !---------------------------------------------------------------------------
  type, extends(cam_grid_attribute_t) :: cam_grid_attribute_1d_int_t
    character(len=max_hcoordname_len)   :: dimname    ! attribute dimension
    integer                             :: dimsize    ! Global array/map size
    integer,        pointer             :: values(:)   => NULL()
    integer(iMap),  pointer             :: map(:) => NULL() ! map (dof) for I/O
  contains
    procedure :: cam_grid_attr_init_1d_int
    procedure :: write_attr => write_cam_grid_attr_1d_int
    procedure :: write_val  => write_cam_grid_val_1d_int
    procedure :: print_attr => print_attr_1d_int
  end type cam_grid_attribute_1d_int_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_1d_r8_t: 1-d real*8 attribute
  !
  !---------------------------------------------------------------------------
  type, extends(cam_grid_attribute_t) :: cam_grid_attribute_1d_r8_t
    character(len=max_hcoordname_len)   :: dimname    ! attribute dimension
    integer                             :: dimsize    ! Global array/map size
    real(r8),       pointer             :: values(:)   => NULL()
    integer(iMap),  pointer             :: map(:) => NULL() ! map (dof) for I/O
  contains
    procedure :: cam_grid_attr_init_1d_r8
    procedure :: write_attr => write_cam_grid_attr_1d_r8
    procedure :: write_val  => write_cam_grid_val_1d_r8
    procedure :: print_attr => print_attr_1d_r8
  end type cam_grid_attribute_1d_r8_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attr_ptr_t: linked list of CAM grid attributes
  !
  !---------------------------------------------------------------------------
  type :: cam_grid_attr_ptr_t
    private
    class(cam_grid_attribute_t), pointer :: attr => NULL()
    type(cam_grid_attr_ptr_t),   pointer :: next => NULL()
  contains
    private
    procedure, public :: initialize => initializeAttrPtr
    procedure, public :: getAttr => getAttrPtrAttr
    procedure, public :: getNext => getAttrPtrNext
    procedure, public :: setNext => setAttrPtrNext
  end type cam_grid_attr_ptr_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_t: Information for a CAM grid (defined by a dycore)
  !
  !---------------------------------------------------------------------------
  type :: cam_grid_t
    character(len=max_hcoordname_len)  :: name = ''     ! grid name
    integer                            :: id            ! e.g., dyn_decomp
    type(horiz_coord_t), pointer       :: lat_coord => NULL() ! Latitude coord
    type(horiz_coord_t), pointer       :: lon_coord => NULL() ! Longitude coord
    logical                            :: unstructured  ! Is this needed?
    logical                            :: block_indexed ! .false. for lon/lat
    logical                            :: attrs_defined = .false.
    type(cam_filemap_t),       pointer :: map => null() ! global dim map (dof)
    type(cam_grid_attr_ptr_t), pointer :: attributes => NULL()
  contains
    procedure :: print_cam_grid
    procedure :: is_unstructured        => cam_grid_unstructured
    procedure :: is_block_indexed       => cam_grid_block_indexed
    procedure :: coord_lengths          => cam_grid_get_dims
    procedure :: coord_names            => cam_grid_coord_names
    procedure :: dim_names              => cam_grid_dim_names
    procedure :: num_elem               => cam_grid_local_size
    procedure :: set_map                => cam_grid_set_map
    procedure :: get_patch_mask         => cam_grid_get_patch_mask
    procedure :: get_lon_lat            => cam_grid_get_lon_lat
    procedure :: find_src_dims          => cam_grid_find_src_dims
    procedure :: find_dest_dims         => cam_grid_find_dest_dims
    procedure :: find_dimids            => cam_grid_find_dimids
    procedure :: get_decomp             => cam_grid_get_pio_decomp
    procedure :: read_darray_2d_int     => cam_grid_read_darray_2d_int
    procedure :: read_darray_3d_int     => cam_grid_read_darray_3d_int
    procedure :: read_darray_2d_double  => cam_grid_read_darray_2d_double
    procedure :: read_darray_3d_double  => cam_grid_read_darray_3d_double
    procedure :: read_darray_2d_real    => cam_grid_read_darray_2d_real
    procedure :: read_darray_3d_real    => cam_grid_read_darray_3d_real
    procedure :: write_darray_2d_int    => cam_grid_write_darray_2d_int
    procedure :: write_darray_3d_int    => cam_grid_write_darray_3d_int
    procedure :: write_darray_2d_double => cam_grid_write_darray_2d_double
    procedure :: write_darray_3d_double => cam_grid_write_darray_3d_double
    procedure :: write_darray_2d_real   => cam_grid_write_darray_2d_real
    procedure :: write_darray_3d_real   => cam_grid_write_darray_3d_real
  end type cam_grid_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_patch_t: Information for a patch of a CAM grid
  !
  !---------------------------------------------------------------------------
  type, public :: cam_grid_patch_t
    private
    integer                      :: grid_id = -1  ! grid containing patch points
    integer                      :: global_size = 0      ! var patch dim size
    integer                      :: global_lat_size = 0  ! lat patch dim size
    integer                      :: global_lon_size = 0  ! lon patch dim size
    real(r8)                     :: lon_range(2)
    real(r8)                     :: lat_range(2)
    type(cam_filemap_t), pointer :: mask       => null() ! map for active pts
    integer(iMap),       pointer :: latmap(:) => null() ! map for patch coords
    integer(iMap),       pointer :: lonmap(:) => null() ! map for patch coords
  contains
    procedure :: gridid              => cam_grid_patch_get_id
    procedure :: get_axis_names      => cam_grid_patch_get_axis_names
    procedure :: get_coord_long_name => cam_grid_patch_get_coord_long_name
    procedure :: get_coord_units     => cam_grid_patch_get_coord_units
    procedure :: set_patch           => cam_grid_patch_set_patch
    procedure :: get_decomp          => cam_grid_patch_get_decomp
    procedure :: compact             => cam_grid_patch_compact
    procedure :: active_cols         => cam_grid_patch_get_active_cols
    procedure :: write_coord_vals    => cam_grid_patch_write_vals
    procedure :: grid_index          => cam_grid_patch_get_grid_index
    procedure :: deallocate          => cam_grid_patch_deallocate
!!XXgoldyXX: PGI workaround?
! COMPILER_BUG(goldy, 2014-11-28, pgi <= 14.9); Commented code should work
!    procedure :: global_size_map     => cam_grid_patch_get_global_size_map
!    procedure :: global_size_axes    => cam_grid_patch_get_global_size_axes
!    generic   :: get_global_size     => global_size_map, global_size_axes
    procedure :: cam_grid_patch_get_global_size_map
    procedure :: cam_grid_patch_get_global_size_axes
    generic   :: get_global_size     => cam_grid_patch_get_global_size_map, cam_grid_patch_get_global_size_axes
  end type cam_grid_patch_t

  !---------------------------------------------------------------------------
  !
  !  cam_grid_header_info_t: Hold NetCDF dimension information for a CAM grid
  !
  !---------------------------------------------------------------------------
  type, public :: cam_grid_header_info_t
    private
    integer                       :: grid_id = -1 ! e.g., dyn_decomp
    integer,          allocatable :: hdims(:)     ! horizontal dimension ids
    type(var_desc_t), pointer     :: lon_varid => NULL() ! lon coord variable
    type(var_desc_t), pointer     :: lat_varid => NULL() ! lat coord variable
  contains
    procedure  :: get_gridid    => cam_grid_header_info_get_gridid
    procedure  :: set_gridid    => cam_grid_header_info_set_gridid
    procedure  :: set_hdims     => cam_grid_header_info_set_hdims
    procedure  :: num_hdims     => cam_grid_header_info_num_hdims
    procedure  :: get_hdimid    => cam_grid_header_info_hdim
    !!XXgoldyXX: Maybe replace this with horiz_coords for patches?
    procedure  :: set_varids    => cam_grid_header_info_set_varids
    procedure  :: get_lon_varid => cam_grid_header_info_lon_varid
    procedure  :: get_lat_varid => cam_grid_header_info_lat_varid
    procedure  :: deallocate    => cam_grid_header_info_deallocate
  end type cam_grid_header_info_t

  !---------------------------------------------------------------------------
  !
  !  END: types BEGIN: interfaces for types
  !
  !---------------------------------------------------------------------------

  ! Abstract interface for write_attr procedure of cam_grid_attribute_t class
  ! NB: This will not compile on some pre-13 Intel compilers
  !     (fails on 12.1.0.233 on Frankfurt, passes on 13.0.1.117 on Yellowstone)
  abstract interface
    subroutine write_cam_grid_attr(attr, File)
      use pio, only: file_desc_t
      import      :: cam_grid_attribute_t
      ! Dummy arguments
      class(cam_grid_attribute_t), intent(inout) :: attr
      type(file_desc_t),           intent(inout) :: File ! PIO file Handle
    end subroutine write_cam_grid_attr
  end interface

  ! Abstract interface for print_attr procedure of cam_grid_attribute_t class
  abstract interface
    subroutine print_attr_spec(this)
      import      :: cam_grid_attribute_t
      ! Dummy arguments
      class(cam_grid_attribute_t), intent(in)    :: this
    end subroutine print_attr_spec
  end interface

  !! Grid variables
  integer, parameter                  :: maxhgrids =  16   ! arbitrary limit
  integer, save                       :: registeredhgrids = 0
  type(cam_grid_t), save              :: cam_grids(maxhgrids)

  public     :: horiz_coord_create

  ! Setup and I/O functions for grids rely on the grid's ID, not its index.
  public     :: cam_grid_register, cam_grid_attribute_register
  public     :: cam_grid_attribute_copy
  public     :: cam_grid_write_attr, cam_grid_write_var
  public     :: cam_grid_read_dist_array, cam_grid_write_dist_array
  ! Access functions for grids rely on the grid's ID, not its index.
  public     :: cam_grid_dimensions, cam_grid_num_grids
  public     :: cam_grid_check ! T/F if grid ID exists
  public     :: cam_grid_id    ! Grid ID (decomp) or -1 if error
  public     :: cam_grid_get_local_size
  public     :: cam_grid_get_file_dimids
  public     :: cam_grid_get_decomp
  public     :: cam_grid_get_gcid
  public     :: cam_grid_get_array_bounds
  public     :: cam_grid_get_coord_names, cam_grid_get_dim_names
  public     :: cam_grid_has_blocksize, cam_grid_get_block_count
  public     :: cam_grid_get_latvals,   cam_grid_get_lonvals
  public     :: cam_grid_is_unstructured, cam_grid_is_block_indexed
  ! Functions for dealing with patch masks
  public     :: cam_grid_compute_patch

  interface cam_grid_attribute_register
    module procedure add_cam_grid_attribute_0d_int
    module procedure add_cam_grid_attribute_0d_char
    module procedure add_cam_grid_attribute_1d_int
    module procedure add_cam_grid_attribute_1d_r8
  end interface

  interface cam_grid_dimensions
    module procedure cam_grid_dimensions_id
    module procedure cam_grid_dimensions_name
  end interface

  interface cam_grid_read_dist_array
    module procedure cam_grid_read_dist_array_2d_int
    module procedure cam_grid_read_dist_array_3d_int
    module procedure cam_grid_read_dist_array_2d_double
    module procedure cam_grid_read_dist_array_3d_double
    module procedure cam_grid_read_dist_array_2d_real
    module procedure cam_grid_read_dist_array_3d_real
  end interface

  interface cam_grid_write_dist_array
    module procedure cam_grid_write_dist_array_2d_int
    module procedure cam_grid_write_dist_array_3d_int
    module procedure cam_grid_write_dist_array_2d_double
    module procedure cam_grid_write_dist_array_3d_double
    module procedure cam_grid_write_dist_array_2d_real
    module procedure cam_grid_write_dist_array_3d_real
  end interface

  ! Private interfaces
  interface get_cam_grid_index
    module procedure get_cam_grid_index_char ! For lookup by name
    module procedure get_cam_grid_index_int  ! For lookup by ID
  end interface

contains

!!#######################################################################
!!
!! Horizontal coordinate functions
!!
!!#######################################################################

  integer function horiz_coord_find_size(this, dimname) result(dimsize)
    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    character(len=*),     intent(in)    :: dimname

    dimsize = -1
    if (len_trim(this%dimname) == 0) then
      if(trim(dimname) == trim(this%name)) then
        dimsize = this%dimsize
      end if
    else
      if(trim(dimname) == trim(this%dimname)) then
        dimsize = this%dimsize
      end if
    end if

  end function horiz_coord_find_size

  integer function horiz_coord_num_elem(this)
    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this

    if (associated(this%values)) then
      horiz_coord_num_elem = size(this%values)
    else
      horiz_coord_num_elem = 0
    end if

  end function horiz_coord_num_elem

  subroutine horiz_coord_len(this, clen)
    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    integer,              intent(out)   :: clen

    clen = this%dimsize
  end subroutine horiz_coord_len
    
  subroutine horiz_coord_name(this, name)
    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    character(len=*),     intent(out)   :: name

    if (len(name) < len_trim(this%name)) then
      call endrun('horiz_coord_name: input name too short')
    end if
    name = trim(this%name)
  end subroutine horiz_coord_name

  subroutine horiz_coord_dim_name(this, dimname)
    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    character(len=*),     intent(out)   :: dimname

    if (len_trim(this%dimname) > 0) then
      ! We have a separate dimension name (e.g., ncol)
      if (len(dimname) < len_trim(this%dimname)) then
        call endrun('horiz_coord_dimname: input name too short')
      end if
      dimname = trim(this%dimname)
    else
      ! No dimension name so we use the coordinate's name
      ! i.e., The dimension name is the same as the coordinate variable
      if (len(dimname) < len_trim(this%name)) then
        call endrun('horiz_coord_dimname: input name too short')
      end if
      dimname = trim(this%name)
    end if
  end subroutine horiz_coord_dim_name

  subroutine horiz_coord_long_name(this, name)

    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    character(len=*),     intent(out)   :: name

    if (len(name) < len_trim(this%long_name)) then
      call endrun('horiz_coord_long_name: input name too short')
    else
      name = trim(this%long_name)
    end if

  end subroutine horiz_coord_long_name

  subroutine horiz_coord_units(this, units)

    ! Dummy arguments
    class(horiz_coord_t), intent(in)    :: this
    character(len=*),     intent(out)   :: units

    if (len(units) < len_trim(this%units)) then
      call endrun('horiz_coord_units: input units too short')
    else
      units = trim(this%units)
    end if

  end subroutine horiz_coord_units

  function horiz_coord_create(name, dimname, dimsize, long_name, units,       &
       lbound, ubound, values, map) result(newcoord)

    ! Dummy arguments
    character(len=*),      intent(in)                  :: name
    character(len=*),      intent(in)                  :: dimname
    integer,               intent(in)                  :: dimsize
    character(len=*),      intent(in)                  :: long_name
    character(len=*),      intent(in)                  :: units
    ! NB: Sure, pointers would have made sense but . . . PGI
    integer,               intent(in)                  :: lbound
    integer,               intent(in)                  :: ubound
    real(r8),              intent(in)                  :: values(lbound:ubound)
    integer(iMap),         intent(in), optional        :: map(ubound-lbound+1)
    type(horiz_coord_t),               pointer         :: newcoord

    ! Local variables
    integer                                            :: i

    allocate(newcoord)

    newcoord%name      = trim(name)
    newcoord%dimname   = trim(dimname)
    newcoord%dimsize   = dimsize
    newcoord%long_name = trim(long_name)
    newcoord%units     = trim(units)
    ! Figure out if this is a latitude or a longitude using CF standard
    ! http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#latitude-coordinate
    ! http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#longitude-coordinate
    if ( (trim(units) == 'degrees_north')    .or.                             &
         (trim(units) == 'degree_north')     .or.                             &
         (trim(units) == 'degree_N')         .or.                             &
         (trim(units) == 'degrees_N')        .or.                             &
         (trim(units) == 'degreeN')          .or.                             &
         (trim(units) == 'degreesN')) then
      newcoord%latitude  = .true.
    else if ((trim(units) == 'degrees_east') .or.                             &
         (trim(units) == 'degree_east')      .or.                             &
         (trim(units) == 'degree_E')         .or.                             &
         (trim(units) == 'degrees_E')        .or.                             &
         (trim(units) == 'degreeE')          .or.                             &
         (trim(units) == 'degreesE')) then
      newcoord%latitude  = .false.
    else
      call endrun("horiz_coord_create: unsupported units: '"//trim(units)//"'")
    end if
    allocate(newcoord%values(lbound:ubound))
    newcoord%values(:) = values(:)

    if (present(map)) then
      if (ANY(map < 0)) then
        call endrun("horiz_coord_create "//trim(name)//": map vals < 0")
      end if
      allocate(newcoord%map(ubound - lbound + 1))
      newcoord%map(:) = map(:)
    else
      nullify(newcoord%map)
    end if

  end function horiz_coord_create

  !---------------------------------------------------------------------------
  !
  !  write_horiz_coord_attr
  !
  !  Write the dimension and coordinate attributes for a horizontal grid
  !  coordinate.
  !
  !---------------------------------------------------------------------------

  subroutine write_horiz_coord_attr(this, File, dimid_out)
    use pio, only: file_desc_t, pio_put_att, pio_noerr, pio_double
    use pio, only: pio_bcast_error, pio_seterrorhandling, pio_inq_varid
    use cam_pio_utils, only: cam_pio_def_dim, cam_pio_def_var

    ! Dummy arguments
    class(horiz_coord_t), intent(inout) :: this
    type(file_desc_t),    intent(inout) :: File         ! PIO file Handle
    integer,    optional, intent(out)   :: dimid_out

    ! Local variables
    type(var_desc_t)                    :: vardesc
    character(len=max_hcoordname_len)   :: dimname
    integer                             :: dimid        ! PIO dimension ID
    integer                             :: err_handling
    integer                             :: ierr

    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = File%iosystem%error_handling !! Hack
    call pio_seterrorhandling(File, PIO_BCAST_ERROR,err_handling)

    ! Make sure the dimension exists in the file
    call this%get_dim_name(dimname)
    call cam_pio_def_dim(File, trim(dimname), this%dimsize, dimid,       &
         existOK=.true.)
    ! Should we define the variable?
    ierr = pio_inq_varid(File, trim(this%name), vardesc)
    if (ierr /= PIO_NOERR) then
      ! Variable not already defined, it is up to us to define the variable
      if (associated(this%vardesc)) then
        ! This should not happen (i.e., internal error)
        call endrun('write_horiz_coord_attr: vardesc already allocated for '//trim(dimname))
      end if
      allocate(this%vardesc)
      call cam_pio_def_var(File, trim(this%name), pio_double,                 &
           (/ dimid /), this%vardesc, existOK=.false.)
      ! long_name
      ierr=pio_put_att(File, this%vardesc, 'long_name', trim(this%long_name))
      call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_horiz_coord_attr')
      ! units
      ierr=pio_put_att(File, this%vardesc, 'units', trim(this%units))
      call cam_pio_handle_error(ierr, 'Error writing "units" attr in write_horiz_coord_attr')
    end if ! We define the variable

    if (present(dimid_out)) then
      dimid_out = dimid
    end if

    ! Back to old error handling
    call pio_seterrorhandling(File, err_handling)

  end subroutine write_horiz_coord_attr

  !---------------------------------------------------------------------------
  !
  !  write_horiz_coord_var
  !
  !  Write the coordinate values for this coordinate
  !
  !---------------------------------------------------------------------------

  subroutine write_horiz_coord_var(this, File)
    use cam_pio_utils, only: cam_pio_get_decomp
    use pio,           only: file_desc_t, pio_double, iosystem_desc_t
    use pio,           only: pio_put_var, pio_write_darray
    use pio,           only: pio_bcast_error, pio_seterrorhandling
    !!XXgoldyXX: HACK to get around circular dependencies. Fix this!!
    !!XXgoldyXX: The issue is cam_pio_utils depending on stuff in this module
    use pio,          only: pio_initdecomp, io_desc_t, pio_freedecomp
    use cam_instance, only: atm_id
    use shr_pio_mod,  only: shr_pio_getiosys
    !!XXgoldyXX: End of this part of the hack


    ! Dummy arguments
    class(horiz_coord_t),    intent(inout) :: this
    type(file_desc_t),       intent(inout) :: File ! PIO file Handle

    ! Local variables
    character(len=120)                     :: errormsg
    integer                                :: ierr
    integer                                :: ldims(1)
    integer                                :: fdims(1)
    integer                                :: err_handling
    type(io_desc_t)                        :: iodesc
    !!XXgoldyXX: HACK to get around circular dependencies. Fix this!!
    type(iosystem_desc_t), pointer         :: piosys
    !!XXgoldyXX: End of this part of the hack

    ! Check to make sure we are supposed to write this var
    if (associated(this%vardesc)) then
      ! We will handle errors for this routine
      !!XXgoldyXX: This hack should be replaced with the PIO interface
      !err_handling = File%iosystem%error_handling !! Hack
      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

      ! Write out the values for this dimension variable
      if (associated(this%map)) then
        ! This is a distributed variable, use pio_write_darray
#if 0
        ldims(1) = this%num_elem()
        call this%get_coord_len(fdims(1))
        allocate(iodesc)
        call cam_pio_get_decomp(iodesc, ldims, fdims, PIO_DOUBLE, this%map)
        call pio_write_darray(File, this%vardesc, iodesc, this%values, ierr, -900._r8)
        nullify(iodesc) ! CAM PIO system takes over memory management of iodesc
#else
        !!XXgoldyXX: HACK to get around circular dependencies. Fix this!!
        piosys => shr_pio_getiosys(atm_id)
        call pio_initdecomp(piosys, pio_double, (/this%dimsize/), this%map,   &
             iodesc)
        call pio_write_darray(File, this%vardesc, iodesc, this%values,        &
             ierr, -900._r8)
        call pio_freedecomp(File, iodesc)
#endif
        !!XXgoldyXX: End of this part of the hack
      else
        ! This is a local variable, pio_put_var should work fine
        ierr = pio_put_var(File, this%vardesc, this%values)
      end if
      write(errormsg, *) 'Error writing variable values for ',trim(this%name),&
           ' in write_horiz_coord_var'
      call cam_pio_handle_error(ierr, errormsg)

      ! Back to old error handling
      call pio_seterrorhandling(File, err_handling)

      ! We are done with this variable descriptor, reset for next file
      deallocate(this%vardesc)
      nullify(this%vardesc)
    end if ! Do we write the variable?

  end subroutine write_horiz_coord_var

!!#######################################################################
!!
!! CAM grid functions
!!
!!#######################################################################

  integer function get_cam_grid_index_char(gridname)
    ! Dummy arguments
    character(len=*), intent(in)  :: gridname
    ! Local variables
    integer :: i
    
    get_cam_grid_index_char = -1
    do i = 1, registeredhgrids
      if(trim(gridname) == trim(cam_grids(i)%name)) then
        get_cam_grid_index_char = i
        exit
      end if
    end do

  end function get_cam_grid_index_char

  integer function get_cam_grid_index_int(gridid)
    ! Dummy arguments
    integer, intent(in) :: gridid
    ! Local variables
    integer :: i
    
    get_cam_grid_index_int = -1
    do i = 1, registeredhgrids
      if(gridid == cam_grids(i)%id) then
        get_cam_grid_index_int = i
        exit
      end if
    end do

  end function get_cam_grid_index_int

  subroutine find_cam_grid_attr(gridind, name, attr)
    ! Dummy arguments
    integer,                              intent(in)     :: gridind
    character(len=*),                     intent(in)     :: name
    class(cam_grid_attribute_t), pointer, intent(out)    :: attr
    ! Local variable
    type(cam_grid_attr_ptr_t),   pointer                 :: attrPtr

    nullify(attr)
    attrPtr => cam_grids(gridind)%attributes
    do while (associated(attrPtr))
!!XXgoldyXX: Is this not working in PGI?
!      attr => attrPtr%getAttr()
      attr => attrPtr%attr
      if (trim(name) == trim(attr%name)) then
        exit
      else
!!XXgoldyXX: Is this not working in PGI?
!        attrPtr => attrPtr%getNext()
        attrPtr => attrPtr%next
        nullify(attr)
      end if
    end do
    return ! attr should be NULL if not found
  end subroutine find_cam_grid_attr

  integer function num_cam_grid_attrs(gridind)
    ! Dummy arguments
    integer,                             intent(in)     :: gridind

    ! Local variables
    class(cam_grid_attr_ptr_t), pointer                 :: attrPtr

    num_cam_grid_attrs = 0
    attrPtr => cam_grids(gridind)%attributes
    do while (associated(attrPtr))
      num_cam_grid_attrs = num_cam_grid_attrs + 1
!!XXgoldyXX: Is this not working in PGI?
!      attrPtr => attrPtr%getNext()
      attrPtr => attrPtr%next
    end do
  end function num_cam_grid_attrs

  subroutine cam_grid_register(name, id, lat_coord, lon_coord, map,           &
       unstruct, block_indexed, src_in, dest_in)
    ! Dummy arguments
    character(len=*),             intent(in) :: name
    integer,                      intent(in) :: id
    type(horiz_coord_t), pointer, intent(in) :: lat_coord
    type(horiz_coord_t), pointer, intent(in) :: lon_coord
    integer(iMap),       pointer, intent(in) :: map(:,:)
    logical,  optional,           intent(in) :: unstruct
    logical,  optional,           intent(in) :: block_indexed
    integer,  optional,           intent(in) :: src_in(2)
    integer,  optional,           intent(in) :: dest_in(2)

    ! Local variables
    character(len=max_hcoordname_len)       :: latdimname, londimname
    character(len=120)                      :: errormsg
    integer                                 :: i
    integer                                 :: src(2), dest(2)
    character(len=*), parameter             :: subname = 'CAM_GRID_REGISTER'

    ! For a values grid, we do not allow multiple calls
    if (get_cam_grid_index(trim(name)) > 0) then
      call endrun(trim(subname)//': Grid, '//trim(name)//', already exists')
    else if (get_cam_grid_index(id) > 0) then
      i = get_cam_grid_index(id)
      write(errormsg, '(4a,i5,3a)') trim(subname), ': Attempt to add grid, ', &
           trim(name), ' with id = ', id, ', however, grid ',                 &
           trim(cam_grids(i)%name), ' already has that ID'
      call endrun(trim(errormsg))
    else if (registeredhgrids >= maxhgrids) then
      call endrun(trim(subname)//": Too many grids")
    else
      registeredhgrids = registeredhgrids + 1
      cam_grids(registeredhgrids)%name       = trim(name)
      cam_grids(registeredhgrids)%id         = id
      ! Quick sanity checks to make sure these aren't mixed up
      if (.not. lat_coord%latitude) then
        call endrun(subname//': lat_coord is not a latitude coordinate')
      end if
      if (lon_coord%latitude) then
        call endrun(subname//': lon_coord is not a longitude coordinate')
      end if
      cam_grids(registeredhgrids)%lat_coord => lat_coord
      cam_grids(registeredhgrids)%lon_coord => lon_coord
      call lat_coord%get_dim_name(latdimname)
      call lon_coord%get_dim_name(londimname)
      if (present(unstruct)) then
        cam_grids(registeredhgrids)%unstructured = unstruct
      else
        if (trim(latdimname) == trim(londimname)) then
          cam_grids(registeredhgrids)%unstructured  = .true.
        else
          cam_grids(registeredhgrids)%unstructured  = .false.
        end if
      end if
      if (present(block_indexed)) then
        cam_grids(registeredhgrids)%block_indexed = block_indexed
      else
        cam_grids(registeredhgrids)%block_indexed = cam_grids(registeredhgrids)%unstructured
      end if
      if (associated(cam_grids(registeredhgrids)%map)) then
        write(errormsg, *) 
        call endrun(trim(subname)//": new grid map should not be associated")
      end if
      if (present(src_in)) then
        src = src_in
      else
        src(1) = 1
        src(2) = -1
      end if
      if (present(dest_in)) then
        dest = dest_in
      else
        dest(1) = 1
        if (cam_grids(registeredhgrids)%unstructured) then
          dest(2) = 0
        else
          dest(2) = 2
        end if
      end if
      allocate(cam_grids(registeredhgrids)%map)
      call cam_grids(registeredhgrids)%map%init(map,                          &
           cam_grids(registeredhgrids)%unstructured, src, dest)
      call cam_grids(registeredhgrids)%print_cam_grid()
    end if

  end subroutine cam_grid_register

  subroutine print_cam_grid(this)
    class(cam_grid_t)                         :: this

    type(cam_grid_attr_ptr_t),   pointer      :: attrPtr
    class(cam_grid_attribute_t), pointer      :: attr
    if (masterproc) then
      write(iulog, '(3a,i4,4a,2(a,l2))') 'Grid: ', trim(this%name),           &
           ', ID = ', this%id,                                                &
           ', lat coord = ', trim(this%lat_coord%name),                       &
           ', lon coord = ', trim(this%lon_coord%name),                       &
           ', unstruct = ', this%unstructured,                                &
           ', block_ind = ', this%block_indexed
      attrPtr => this%attributes
      do while (associated(attrPtr))
!!XXgoldyXX: Is this not working in PGI?
!      attr => attrPtr%getAttr()
      attr => attrPtr%attr
        call attr%print_attr()
!!XXgoldyXX: Is this not working in PGI?
!      attrPtr => attrPtr%getNext()
      attrPtr => attrPtr%next
      end do
    end if
  end subroutine print_cam_grid

  integer function cam_grid_num_grids()
    cam_grid_num_grids = registeredhgrids
  end function cam_grid_num_grids

  ! Return .true. iff id represents a valid CAM grid
  logical function cam_grid_check(id)
    ! Dummy argument
    integer, intent(in)    :: id

    cam_grid_check = ((get_cam_grid_index(id) > 0) .and.                      &
         (get_cam_grid_index(id) <= cam_grid_num_grids()))
  end function cam_grid_check

  integer function cam_grid_id(name)
    ! Dummy argument
    character(len=*),   intent(in)    :: name

    ! Local variable
    integer                           :: index

    index = get_cam_grid_index(name)
    if (index > 0) then
      cam_grid_id = cam_grids(index)%id
    else
      cam_grid_id = -1
    end if

  end function cam_grid_id

  ! Return the size of a local array for grid, ID.
  ! With no optional argument, return the basic 2D array size
  ! nlev represents levels or the total column size (product(mdims))
  integer function cam_grid_get_local_size(id, nlev)

    ! Dummy arguments
    integer,                    intent(in)    :: id
    integer,          optional, intent(in)    :: nlev

    ! Local variables
    integer                                   :: gridid
    character(len=128)                        :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      cam_grid_get_local_size = cam_grids(gridid)%num_elem()
      if (present(nlev)) then
        cam_grid_get_local_size = cam_grid_get_local_size * nlev
      end if
    else
      write(errormsg, *) 'cam_grid_get_local_size: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end function cam_grid_get_local_size

  ! Given some array information, find the dimension NetCDF IDs on <File> for this grid
  subroutine cam_grid_get_file_dimids(id, File, dimids)
    use pio,           only: file_desc_t

    ! Dummy arguments
    integer,                   intent(in)    :: id
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(out)   :: dimids(:)

    ! Local variables
    integer                                  :: gridid
    character(len=128)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%find_dimids(File, dimids)
    else
      write(errormsg, *) 'cam_grid_get_file_dimids: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_get_file_dimids

  ! Given some array information, find or compute a PIO decomposition
  subroutine cam_grid_get_decomp(id, field_lens, file_lens, dtype, iodesc,    &
       field_dnames, file_dnames)
    use pio,           only: io_desc_t

    ! Dummy arguments
    integer,                    intent(in)    :: id
    integer,                    intent(in)    :: field_lens(:) ! Array dim sizes
    integer,                    intent(in)    :: file_lens(:)  ! File dim sizes
    integer,                    intent(in)    :: dtype
    type(io_desc_t),  pointer,  intent(out)   :: iodesc
    character(len=*), optional, intent(in)    :: field_dnames(:)
    character(len=*), optional, intent(in)    :: file_dnames(:)

    ! Local variables
    integer                                   :: gridid
    character(len=128)                        :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%get_decomp(field_lens, file_lens, dtype, iodesc, &
           field_dnames, file_dnames)
    else
      write(errormsg, *) 'cam_grid_get_decomp: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_get_decomp

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_2d_int
  !
  !  Interface function for the grid%read_darray_2d_int method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_2d_int(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_2d_int(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_2d_int: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_2d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_3d_int
  !
  !  Interface function for the grid%read_darray_2d_ method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_3d_int(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_3d_int(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_3d_int: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_3d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_2d_double
  !
  !  Interface function for the grid%read_darray_2d_double method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_2d_double(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_2d_double(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_2d_double: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_2d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_3d_double
  !
  !  Interface function for the grid%read_darray_3d_double method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_3d_double(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_3d_double(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_3d_double: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_3d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_2d_real
  !
  !  Interface function for the grid%read_darray_2d_real method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_2d_real(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_2d_real(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_2d_real: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_2d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_dist_array_3d_real
  !
  !  Interface function for the grid%read_darray_3d_real method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_dist_array_3d_real(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%read_darray_3d_real(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_read_dist_array_3d_real: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_read_dist_array_3d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_2d_int
  !
  !  Interface function for the grid%write_darray_2d_int method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_2d_int(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_2d_int(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_2d_int: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_2d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_3d_int
  !
  !  Interface function for the grid%write_darray_3d_int method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_3d_int(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_3d_int(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_3d_int: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_3d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_2d_double
  !
  !  Interface function for the grid%write_darray_2d_double method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_2d_double(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_2d_double(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_2d_double: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_2d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_3d_double
  !
  !  Interface function for the grid%write_darray_3d_double method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_3d_double(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_3d_double(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_3d_double: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_3d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_2d_real
  !
  !  Interface function for the grid%write_darray_2d_real method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_2d_real(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_2d_real(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_2d_real: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_2d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_dist_array_3d_real
  !
  !  Interface function for the grid%write_darray_3d_real method
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_dist_array_3d_real(File, id, adims, fdims, hbuf, varid)
    use pio, only: file_desc_t

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: id
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variable
    integer                                  :: gridid
    character(len=120)                       :: errormsg

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%write_darray_3d_real(File, adims, fdims, hbuf, varid)
    else
      write(errormsg, *) 'cam_grid_write_dist_array_3d_real: Bad grid ID, ', id
      call endrun(errormsg)
    end if

  end subroutine cam_grid_write_dist_array_3d_real

  subroutine cam_grid_get_gcid(id, gcid)

    ! Dummy arguments
    integer,                      intent(in)    :: id
    integer(iMap), pointer                      :: gcid(:)

    ! Local variables
    integer                                     :: gridid
    integer                                     :: fieldbounds(2,2)
    integer                                     :: fieldlens(2)
    integer                                     :: filelens(2)
    type(cam_filemap_t), pointer                :: map

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      map => cam_grids(gridid)%map
      call cam_grids(gridid)%coord_lengths(filelens)
      call map%array_bounds(fieldbounds)
      fieldlens(:) = fieldbounds(:,2) - fieldbounds(:,1) + 1
      call map%get_filemap(fieldlens, filelens, gcid)
    else
      call endrun('cam_grid_get_gcid: Bad grid ID')
    end if
 end subroutine cam_grid_get_gcid

  subroutine cam_grid_get_array_bounds(id, dims)

    ! Dummy arguments
    integer,                  intent(in)    :: id
    integer,                  intent(inout) :: dims(:,:)

    ! Local variables
    integer                                 :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      if (.not. associated(cam_grids(gridid)%map)) then
        call endrun('cam_grid_get_array_bounds: Grid, '//trim(cam_grids(gridid)%name)//', has no map')
      else
        call cam_grids(gridid)%map%array_bounds(dims)
      end if
    else
      call endrun('cam_grid_get_array_bounds: Bad grid ID')
    end if

  end subroutine cam_grid_get_array_bounds

  subroutine cam_grid_get_coord_names(id, name1, name2)

    ! Dummy arguments
    integer,                  intent(in)    :: id
    character(len=*),         intent(out)   :: name1
    character(len=*),         intent(out)   :: name2

    ! Local variables
    integer                                 :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%coord_names(name1, name2)
    else
      call endrun('cam_grid_get_coord_names: Bad grid ID')
    end if

  end subroutine cam_grid_get_coord_names

  subroutine cam_grid_get_dim_names(id, name1, name2)

    ! Dummy arguments
    integer,                  intent(in)    :: id
    character(len=*),         intent(out)   :: name1
    character(len=*),         intent(out)   :: name2

    ! Local variables
    integer                                 :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%dim_names(name1, name2)
    else
      call endrun('cam_grid_get_dim_names: Bad grid ID')
    end if

  end subroutine cam_grid_get_dim_names

  logical function cam_grid_has_blocksize(id)

    ! Dummy arguments
    integer,                  intent(in)    :: id

    ! Local variables
    integer                                 :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      if (.not. associated(cam_grids(gridid)%map)) then
        call endrun('cam_grid_has_blocksize: Grid, '//trim(cam_grids(gridid)%name)//', has no map')
      else
        cam_grid_has_blocksize = cam_grids(gridid)%map%has_blocksize()
      end if
    else
      call endrun('cam_grid_has_blocksize: Bad grid ID')
    end if
  end function cam_grid_has_blocksize

  ! Return the number of active columns in the block specified by block_id
  integer function cam_grid_get_block_count(id, block_id) result(ncol)

    ! Dummy arguments
    integer,                  intent(in)    :: id
    integer,                  intent(in)    :: block_id

    ! Local variables
    integer                                 :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      if (.not. associated(cam_grids(gridid)%map)) then
        call endrun('cam_grid_get_block_count: Grid, '//trim(cam_grids(gridid)%name)//', has no map')
      else
        ncol = cam_grids(gridid)%map%blocksize(block_id)
      end if
    else
      call endrun('cam_grid_get_block_count: Bad grid ID')
    end if
  end function cam_grid_get_block_count

  function cam_grid_get_latvals(id) result(latvals)

    ! Dummy argument
    integer,                  intent(in) :: id
    real(r8), pointer                    :: latvals(:)

    ! Local variables
    integer                              :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      if (.not. associated(cam_grids(gridid)%lat_coord%values)) then
        nullify(latvals)
      else
        latvals => cam_grids(gridid)%lat_coord%values
      end if
    else
      call endrun('cam_grid_get_latvals: Bad grid ID')
    end if
  end function cam_grid_get_latvals

  function cam_grid_get_lonvals(id) result(lonvals)

    ! Dummy arguments
    integer,                  intent(in) :: id
    real(r8), pointer                    :: lonvals(:)

    ! Local variables
    integer                              :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      if (.not. associated(cam_grids(gridid)%lon_coord%values)) then
        nullify(lonvals)
      else
        lonvals => cam_grids(gridid)%lon_coord%values
      end if
    else
      call endrun('cam_grid_get_lonvals: Bad grid ID')
    end if
  end function cam_grid_get_lonvals

  logical function cam_grid_is_unstructured(id) result(unstruct)

    ! Dummy arguments
    integer,                  intent(in) :: id

    ! Local variables
    integer                              :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      unstruct = cam_grids(gridid)%is_unstructured()
    else
      call endrun('cam_grid_is_unstructured: Bad grid ID')
    end if
  end function cam_grid_is_unstructured

  logical function cam_grid_is_block_indexed(id) result(block_indexed)

    ! Dummy arguments
    integer,                  intent(in) :: id

    ! Local variables
    integer                              :: gridid
    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      block_indexed = cam_grids(gridid)%is_block_indexed()
    else
      call endrun('s: Bad grid ID')
    end if
  end function cam_grid_is_block_indexed

  ! Compute or update a grid patch mask
  subroutine cam_grid_compute_patch(id, patch, lonl, lonu, latl, latu)

    ! Dummy arguments
    integer,                         intent(in)    :: id
    type(cam_grid_patch_t),          intent(inout) :: patch
    real(r8),                        intent(in)    :: lonl
    real(r8),                        intent(in)    :: lonu
    real(r8),                        intent(in)    :: latl
    real(r8),                        intent(in)    :: latu

    ! Local variables
    integer                                        :: gridid

    gridid = get_cam_grid_index(id)
    if (gridid > 0) then
      call cam_grids(gridid)%get_patch_mask(lonl, lonu, latl, latu, patch)
    else
      call endrun('cam_grid_compute_patch: Bad grid ID')
    end if

  end subroutine cam_grid_compute_patch

!!#######################################################################
!!
!! CAM grid attribute functions
!!
!!#######################################################################

  subroutine cam_grid_attr_init(this, name, long_name, next)
    ! Dummy arguments
    class(cam_grid_attribute_t)                        :: this
    character(len=*),                    intent(in)    :: name
    character(len=*),                    intent(in)    :: long_name
    class(cam_grid_attribute_t), pointer               :: next

    this%name      = trim(name)
    this%long_name = trim(long_name)
    this%next => next
  end subroutine cam_grid_attr_init

  subroutine print_attr_base(this)
    ! Dummy arguments
    class(cam_grid_attribute_t), intent(in)             :: this
    if (masterproc) then
      write(iulog, '(5a)') 'Attribute: ', trim(this%name), ", long name = '", &
           trim(this%long_name), "'"
    end if
  end subroutine print_attr_base

  subroutine cam_grid_attr_init_0d_int(this, name, long_name, val)
    ! Dummy arguments
    class(cam_grid_attribute_0d_int_t)                  :: this
    character(len=*),                    intent(in)     :: name
    character(len=*),                    intent(in)     :: long_name
    integer,                             intent(in)     :: val

!    call this%cam_grid_attr_init(name, '')
    this%name      = trim(name)
    this%long_name = trim(long_name)
    this%ival      = val
  end subroutine cam_grid_attr_init_0d_int

  subroutine print_attr_0d_int(this)
    ! Dummy arguments
    class(cam_grid_attribute_0d_int_t), intent(in)      :: this

    call this%print_attr_base()
    if (masterproc) then
      write(iulog, *) '          value = ', this%ival
    end if
  end subroutine print_attr_0d_int

  subroutine cam_grid_attr_init_0d_char(this, name, long_name, val)
    ! Dummy arguments
    class(cam_grid_attribute_0d_char_t)                 :: this
    character(len=*),                    intent(in)     :: name
    character(len=*),                    intent(in)     :: long_name
    character(len=*),                    intent(in)     :: val

!    call this%cam_grid_attr_init(name, '')
    this%name      = trim(name)
    this%long_name = trim(long_name)
    this%val       = trim(val)
  end subroutine cam_grid_attr_init_0d_char

  subroutine print_attr_0d_char(this)
    ! Dummy arguments
    class(cam_grid_attribute_0d_char_t), intent(in)     :: this

    call this%print_attr_base()
    if (masterproc) then
      write(iulog, *) '          value = ', trim(this%val)
    end if
  end subroutine print_attr_0d_char

  subroutine cam_grid_attr_init_1d_int(this, name, long_name, dimname,        &
       dimsize, values, map)
    ! Dummy arguments
    class(cam_grid_attribute_1d_int_t)                  :: this
    character(len=*),                    intent(in)     :: name
    character(len=*),                    intent(in)     :: long_name
    character(len=*),                    intent(in)     :: dimname
    integer,                             intent(in)     :: dimsize
    integer,                     target, intent(in)     :: values(:)
    integer(iMap),     optional, target, intent(in)     :: map(:)

!    call this%cam_grid_attr_init(trim(name), trim(long_name))
    if (len_trim(name) > max_hcoordname_len) then
       call endrun('cam_grid_attr_1d_int: name too long')
    end if
    this%name      = trim(name)
    if (len_trim(long_name) > max_chars) then
       call endrun('cam_grid_attr_1d_int: long_name too long')
    end if
    this%long_name = trim(long_name)
    
    if (len_trim(dimname) > max_hcoordname_len) then
       call endrun('cam_grid_attr_1d_int: dimname too long')
    end if
    this%dimname =  trim(dimname)
    this%dimsize =  dimsize
    this%values  => values
    ! Fill in the optional map
    if (present(map)) then
      this%map => map
    else
      nullify(this%map)
    end if
  end subroutine cam_grid_attr_init_1d_int

  subroutine cam_grid_attr_init_1d_r8(this, name, long_name, dimname,         &
       dimsize, values, map)
    ! Dummy arguments
    class(cam_grid_attribute_1d_r8_t)                   :: this
    character(len=*),                    intent(in)     :: name
    character(len=*),                    intent(in)     :: long_name
    character(len=*),                    intent(in)     :: dimname
    integer,                             intent(in)     :: dimsize
    real(r8),                    target, intent(in)     :: values(:)
    integer(iMap),     optional, target, intent(in)     :: map(:)

!    call this%cam_grid_attr_init(trim(name), trim(long_name), next)
    this%name      = trim(name)
    this%long_name = trim(long_name)
    
    this%dimname =  trim(dimname)
    this%dimsize =  dimsize
    this%values  => values
    ! Fill in the optional map
    if (present(map)) then
      this%map => map
    else
      nullify(this%map)
    end if
  end subroutine cam_grid_attr_init_1d_r8

  subroutine print_attr_1d_int(this)
    ! Dummy arguments
    class(cam_grid_attribute_1d_int_t), intent(in)      :: this
    call this%print_attr_base()
    if (masterproc) then
      write(iulog, *) '          dimname = ', trim(this%dimname)
    end if
  end subroutine print_attr_1d_int

  subroutine print_attr_1d_r8(this)
    ! Dummy arguments
    class(cam_grid_attribute_1d_r8_t), intent(in)       :: this
    call this%print_attr_base()
    if (masterproc) then
      write(iulog, *) '          dimname = ', trim(this%dimname)
    end if
  end subroutine print_attr_1d_r8

  subroutine insert_grid_attribute(gridind, attr)
    integer,                              intent(in) :: gridind
    class(cam_grid_attribute_t), pointer             :: attr

    ! Push a new attribute onto the grid
    type(cam_grid_attr_ptr_t),  pointer              :: attrPtr => NULL()
    allocate(attrPtr)
    call attrPtr%initialize(attr)
    call attrPtr%setNext(cam_grids(gridind)%attributes)
    cam_grids(gridind)%attributes => attrPtr
    call attrPtr%attr%print_attr()
  end subroutine insert_grid_attribute

  subroutine add_cam_grid_attribute_0d_int(gridname, name, long_name, val)
    ! Dummy arguments
    character(len=*),      intent(in)                   :: gridname
    character(len=*),      intent(in)                   :: name
    character(len=*),      intent(in)                   :: long_name
    integer,               intent(in)                   :: val

    ! Local variables
    type(cam_grid_attribute_0d_int_t), pointer          :: attr
    class(cam_grid_attribute_t),       pointer          :: attptr
    character(len=120)                                  :: errormsg
    integer                                             :: gridind

    gridind = get_cam_grid_index(trim(gridname))
    if (gridind > 0) then
      call find_cam_grid_attr(gridind, trim(name), attptr)
      if (associated(attptr)) then
        ! Attribute found, can't add it again!
        write(errormsg, '(4a)')                                               &
             'add_cam_grid_attribute_0d_int: attribute ', trim(name),         &
           ' already exists for ', cam_grids(gridind)%name
        call endrun(errormsg)
      else
        ! Need a new attribute.
        allocate(attr)
        call attr%cam_grid_attr_init_0d_int(trim(name), trim(long_name), val)
        attptr => attr
        call insert_grid_attribute(gridind, attptr)
      end if
    else
      write(errormsg, '(3a)') 'add_cam_grid_attribute_0d_int: grid ',         &
           trim(gridname), ' was not found'
      call endrun(errormsg)
    end if
!    call cam_grids(gridind)%print_cam_grid()
  end subroutine add_cam_grid_attribute_0d_int

  subroutine add_cam_grid_attribute_0d_char(gridname, name, val)
    ! Dummy arguments
    character(len=*),      intent(in)                   :: gridname
    character(len=*),      intent(in)                   :: name
    character(len=*),      intent(in)                   :: val

    ! Local variables
    type(cam_grid_attribute_0d_char_t), pointer         :: attr
    class(cam_grid_attribute_t),        pointer         :: attptr
    character(len=120)                                  :: errormsg
    integer                                             :: gridind

    gridind = get_cam_grid_index(trim(gridname))
    if (gridind > 0) then
      call find_cam_grid_attr(gridind, trim(name), attptr)
      if (associated(attptr)) then
        ! Attribute found, can't add it again!
        write(errormsg, '(4a)')                                               &
             'add_cam_grid_attribute_0d_char: attribute ', trim(name),        &
           ' already exists for ', cam_grids(gridind)%name
        call endrun(errormsg)
      else
        ! Need a new attribute.
        allocate(attr)
        call attr%cam_grid_attr_init_0d_char(trim(name), '', val)
        attptr => attr
        call insert_grid_attribute(gridind, attptr)
      end if
    else
      write(errormsg, '(3a)') 'add_cam_grid_attribute_0d_char: grid ',        &
           trim(gridname), ' was not found'
      call endrun(errormsg)
    end if
!    call cam_grids(gridind)%print_cam_grid()
  end subroutine add_cam_grid_attribute_0d_char

  subroutine add_cam_grid_attribute_1d_int(gridname, name, long_name,         &
       dimname, values, map)
    ! Dummy arguments
    character(len=*),      intent(in)                   :: gridname
    character(len=*),      intent(in)                   :: name
    character(len=*),      intent(in)                   :: long_name
    character(len=*),      intent(in)                   :: dimname
    integer,               intent(in), target           :: values(:)
    integer(iMap),         intent(in), target, optional :: map(:)

    ! Local variables
    type(cam_grid_attribute_1d_int_t), pointer          :: attr   => NULL()
    class(cam_grid_attribute_t),       pointer          :: attptr => NULL()
    character(len=120)                                  :: errormsg
    integer                                             :: gridind
    integer                                             :: dimsize

    gridind = get_cam_grid_index(trim(gridname))
    if (gridind > 0) then
      call find_cam_grid_attr(gridind, trim(name), attptr)
      if (associated(attptr)) then
        ! Attribute found, can't add it again!
        write(errormsg, '(4a)')                                               &
             'add_cam_grid_attribute_1d_int: attribute ', trim(name),         &
             ' already exists for ', cam_grids(gridind)%name
        call endrun(errormsg)
      else
        ! Need a new attribute.
        dimsize = cam_grids(gridind)%lat_coord%global_size(trim(dimname))
        if (dimsize < 1) then
          dimsize = cam_grids(gridind)%lon_coord%global_size(trim(dimname))
        end if
        if (dimsize < 1) then
          write(errormsg, *) 'add_cam_grid_attribute_1d_int: attribute ',     &
               'dimension ', trim(dimname), ' for ', trim(name), ', not found'
          call endrun(errormsg)
        end if
        allocate(attr)
        call attr%cam_grid_attr_init_1d_int(trim(name), trim(long_name),      &
             trim(dimname), dimsize, values, map)
        attptr => attr
        call insert_grid_attribute(gridind, attptr)
      end if
    else
      write(errormsg, '(3a)') 'add_cam_grid_attribute_1d_int: grid ',         &
           trim(gridname), ' was not found'
      call endrun(errormsg)
    end if
!    call cam_grids(gridind)%print_cam_grid()
  end subroutine add_cam_grid_attribute_1d_int

  subroutine add_cam_grid_attribute_1d_r8(gridname, name, long_name,          &
       dimname, values, map)
    ! Dummy arguments
    character(len=*),      intent(in)                   :: gridname
    character(len=*),      intent(in)                   :: name
    character(len=*),      intent(in)                   :: long_name
    character(len=*),      intent(in)                   :: dimname
    real(r8),              intent(in), target           :: values(:)
    integer(iMap),         intent(in), target, optional :: map(:)

    ! Local variables
    type(cam_grid_attribute_1d_r8_t),  pointer          :: attr
    class(cam_grid_attribute_t),       pointer          :: attptr
    character(len=120)                                  :: errormsg
    integer                                             :: gridind
    integer                                             :: dimsize

    gridind = get_cam_grid_index(trim(gridname))
    if (gridind > 0) then
      call find_cam_grid_attr(gridind, trim(name), attptr)
      if (associated(attptr)) then
        ! Attribute found, can't add it again!
        write(errormsg, '(4a)')                                               &
             'add_cam_grid_attribute_1d_r8: attribute ', trim(name),          &
             ' already exists for ', cam_grids(gridind)%name
        call endrun(errormsg)
      else
        ! Need a new attribute.
        dimsize = cam_grids(gridind)%lat_coord%global_size(trim(dimname))
        if (dimsize < 1) then
          dimsize = cam_grids(gridind)%lon_coord%global_size(trim(dimname))
        end if
        if (dimsize < 1) then
          write(errormsg, *) 'add_cam_grid_attribute_1d_r8: attribute ',      &
               'dimension ', trim(dimname), ' for ', trim(name), ', not found'
          call endrun(errormsg)
        end if
        allocate(attr)
        call attr%cam_grid_attr_init_1d_r8(trim(name), trim(long_name),       &
             trim(dimname), dimsize, values, map)
        attptr => attr
        call insert_grid_attribute(gridind, attptr)
      end if
    else
      write(errormsg, '(3a)') 'add_cam_grid_attribute_1d_r8: grid ',          &
           trim(gridname), ' was not found'
      call endrun(errormsg)
    end if
!    call cam_grids(gridind)%print_cam_grid()
  end subroutine add_cam_grid_attribute_1d_r8

!!#######################################################################
!!
!! CAM grid attribute pointer (list node) functions
!!
!!#######################################################################

  subroutine initializeAttrPtr(this, attr)
    ! Dummy arguments
    class(cam_grid_attr_ptr_t)           :: this
    class(cam_grid_attribute_t), target  :: attr

    if (associated(this%next)) then
      if (masterproc) then
        write(iulog, *) 'WARNING: Overwriting attr pointer for cam_grid_attr_ptr_t'
      end if
    end if
    this%attr => attr
  end subroutine initializeAttrPtr

  function getAttrPtrAttr(this)
    ! Dummy variable
    class(cam_grid_attr_ptr_t)                 :: this
    class(cam_grid_attribute_t), pointer       :: getAttrPtrAttr

    getAttrPtrAttr => this%attr
  end function getAttrPtrAttr

  function getAttrPtrNext(this)
    ! Dummy arguments
    class(cam_grid_attr_ptr_t)                 :: this
    type(cam_grid_attr_ptr_t), pointer         :: getAttrPtrNext

    getAttrPtrNext => this%next
  end function getAttrPtrNext

  subroutine setAttrPtrNext(this, next)
    ! Dummy arguments
    class(cam_grid_attr_ptr_t)                 :: this
    type(cam_grid_attr_ptr_t),  pointer        :: next

    if (associated(this%next)) then
      if (masterproc) then
        write(iulog, *) 'WARNING: Overwriting next pointer for cam_grid_attr_ptr_t'
      end if
    end if
    this%next => next
  end subroutine setAttrPtrNext

  !---------------------------------------------------------------------------
  !
  !  write_cam_grid_attr_0d_int
  !
  !  Write a grid attribute
  !
  !---------------------------------------------------------------------------

  subroutine write_cam_grid_attr_0d_int(attr, File)
    use pio,           only: file_desc_t, pio_put_att, pio_noerr, pio_int,    &
         pio_inq_att, PIO_GLOBAL, PIO_OFFSET_KIND
    use cam_pio_utils, only: cam_pio_def_var

    ! Dummy arguments
    class(cam_grid_attribute_0d_int_t), intent(inout) :: attr
    type(file_desc_t),                  intent(inout) :: File ! PIO file Handle

    ! Local variables
    character(len=120)                  :: errormsg
    integer                             :: attrtype
    integer(PIO_OFFSET_KIND)            :: attrlen
    integer                             :: ierr

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is associated, that grid defined the attribute
    if (.not. associated(attr%vardesc)) then
      if (len_trim(attr%long_name) > 0) then
        ! This 0d attribute is a scalar variable with a long_name attribute
        ! First, define the variable
        allocate(attr%vardesc)
        call cam_pio_def_var(File, trim(attr%name), pio_int, attr%vardesc,    &
             existOK=.false.)
        ierr=pio_put_att(File, attr%vardesc, 'long_name', trim(attr%long_name))
        call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_cam_grid_attr_0d_int')
      else
        ! This 0d attribute is a global attribute
        ! Check to see if the attribute already exists in the file
        ierr = pio_inq_att(File, PIO_GLOBAL, attr%name, attrtype, attrlen)
        if (ierr /= PIO_NOERR) then
          ! Time to define the attribute
          ierr = pio_put_att(File, PIO_GLOBAL, trim(attr%name), attr%ival)
          call cam_pio_handle_error(ierr, 'Unable to define attribute in write_cam_grid_attr_0d_int')
        end if
      end if
    end if

  end subroutine write_cam_grid_attr_0d_int

  !---------------------------------------------------------------------------
  !
  !  write_cam_grid_attr_0d_char
  !
  !  Write a grid attribute
  !
  !---------------------------------------------------------------------------

  subroutine write_cam_grid_attr_0d_char(attr, File)
    use pio, only: file_desc_t, pio_put_att, pio_noerr,                       &
                   pio_inq_att, PIO_GLOBAL, PIO_OFFSET_KIND

    ! Dummy arguments
    class(cam_grid_attribute_0d_char_t), intent(inout) :: attr
    type(file_desc_t),                   intent(inout) :: File ! PIO file Handle

    ! Local variables
    character(len=120)                  :: errormsg
    integer                             :: attrtype
    integer(PIO_OFFSET_KIND)            :: attrlen
    integer                             :: ierr

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is associated, that grid defined the attribute
    if (.not. associated(attr%vardesc)) then
      ! The 0d char attributes are global attribues
      ! Check to see if the attribute already exists in the file
      ierr = pio_inq_att(File, PIO_GLOBAL, attr%name, attrtype, attrlen)
      if (ierr /= PIO_NOERR) then
        ! Time to define the variable
        ierr = pio_put_att(File, PIO_GLOBAL, trim(attr%name), attr%val)
        call cam_pio_handle_error(ierr, 'Unable to define attribute in write_cam_grid_attr_0d_char')
      end if
    end if

  end subroutine write_cam_grid_attr_0d_char

  !---------------------------------------------------------------------------
  !
  !  write_cam_grid_attr_1d_int
  !
  !  Write a grid attribute
  !
  !---------------------------------------------------------------------------

  subroutine write_cam_grid_attr_1d_int(attr, File)
    use pio,           only: file_desc_t, pio_put_att, pio_noerr
    use pio,           only: pio_inq_dimid, pio_int
    use cam_pio_utils, only: cam_pio_def_var

    ! Dummy arguments
    class(cam_grid_attribute_1d_int_t), intent(inout) :: attr
    type(file_desc_t),                  intent(inout) :: File ! PIO file Handle

    ! Local variables
    integer                             :: dimid      ! PIO dimension ID
    character(len=120)                  :: errormsg
    integer                             :: ierr

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is associated, that grid defined the attribute
    if (.not. associated(attr%vardesc)) then
      ! Check to see if the dimension already exists in the file
      ierr = pio_inq_dimid(File, trim(attr%dimname), dimid)
      if (ierr /= PIO_NOERR) then
        ! The dimension has not yet been defined. This is an error
        ! NB: It should have been defined as part of a coordinate
        write(errormsg, *) 'write_cam_grid_attr_1d_int: dimension, ',         &
             trim(attr%dimname), ', does not exist'
        call endrun(errormsg)
      end if
      ! Time to define the variable
      allocate(attr%vardesc)
      call cam_pio_def_var(File, trim(attr%name), pio_int, (/dimid/),         &
           attr%vardesc, existOK=.false.)
      ierr = pio_put_att(File, attr%vardesc, 'long_name', trim(attr%long_name))
      call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_cam_grid_attr_1d_int')
    end if

  end subroutine write_cam_grid_attr_1d_int

  !---------------------------------------------------------------------------
  !
  !  write_cam_grid_attr_1d_r8
  !
  !  Write a grid attribute
  !
  !---------------------------------------------------------------------------

  subroutine write_cam_grid_attr_1d_r8(attr, File)
    use pio,           only: file_desc_t, pio_put_att, pio_noerr, pio_double, &
         pio_inq_dimid
    use cam_pio_utils, only: cam_pio_def_var

    ! Dummy arguments
    class(cam_grid_attribute_1d_r8_t), intent(inout) :: attr
    type(file_desc_t),                 intent(inout) :: File ! PIO file Handle

    ! Local variables
    integer                             :: dimid      ! PIO dimension ID
    character(len=120)                  :: errormsg
    integer                             :: ierr

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is associated, that grid defined the attribute
    if (.not. associated(attr%vardesc)) then
      ! Check to see if the dimension already exists in the file
      ierr = pio_inq_dimid(File, trim(attr%dimname), dimid)
      if (ierr /= PIO_NOERR) then
        ! The dimension has not yet been defined. This is an error
        ! NB: It should have been defined as part of a coordinate
        write(errormsg, *) 'write_cam_grid_attr_1d_r8: dimension, ',          &
             trim(attr%dimname), ', does not exist'
        call endrun(errormsg)
      end if
      ! Time to define the variable
      allocate(attr%vardesc)
      call cam_pio_def_var(File, trim(attr%name), pio_double, (/dimid/),      &
           attr%vardesc, existOK=.false.)
      ! long_name
      ierr = pio_put_att(File, attr%vardesc, 'long_name', trim(attr%long_name))
      call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_cam_grid_attr_1d_r8')
    end if

  end subroutine write_cam_grid_attr_1d_r8

  !---------------------------------------------------------------------------
  !
  !  cam_grid_attribute_copy
  !
  !  Copy an attribute from a source grid to a destination grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_attribute_copy(src_grid, dest_grid, attribute_name)
    ! Dummy arguments
    character(len=*),         intent(in) :: src_grid
    character(len=*),         intent(in) :: dest_grid
    character(len=*),         intent(in) :: attribute_name

    ! Local variables
    character(len=120)                   :: errormsg
    integer                              :: src_ind, dest_ind
    class(cam_grid_attribute_t), pointer :: attr

    ! Find the source and destination grid indices
    src_ind = get_cam_grid_index(trim(src_grid))
    dest_ind = get_cam_grid_index(trim(dest_grid))

    call find_cam_grid_attr(dest_ind, trim(attribute_name), attr)
    if (associated(attr)) then
      ! Attribute found, can't add it again!
      write(errormsg, '(4a)') 'CAM_GRID_ATTRIBUTE_COPY: attribute ',          &
           trim(attribute_name),' already exists for ',cam_grids(dest_ind)%name
      call endrun(errormsg)
    else
      call find_cam_grid_attr(src_ind, trim(attribute_name), attr)
      if (associated(attr)) then
        ! Copy the attribute
        call insert_grid_attribute(dest_ind, attr)
      else
        write(errormsg, '(4a)') ": Did not find attribute, '",                &
             trim(attribute_name), "' in ", cam_grids(src_ind)%name
        call endrun("CAM_GRID_ATTRIBUTE_COPY"//errormsg)
      end if
    end if

  end subroutine cam_grid_attribute_copy

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_attr
  !
  !  Write the dimension and coordinate attributes for the horizontal history
  !  coordinates.
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_attr(File, grid_id, header_info)
    use pio, only: file_desc_t, PIO_BCAST_ERROR, pio_seterrorhandling
    use pio, only: pio_inq_dimid

    ! Dummy arguments
    type(file_desc_t),            intent(inout) :: File       ! PIO file Handle
    integer,                      intent(in)    :: grid_id
    type(cam_grid_header_info_t), intent(inout) :: header_info

    ! Local variables
    integer                                     :: gridind
    class(cam_grid_attribute_t), pointer        :: attr
    type(cam_grid_attr_ptr_t),   pointer        :: attrPtr
    integer                                     :: dimids(2)
    integer                                     :: err_handling

    gridind = get_cam_grid_index(grid_id)
    !! Fill this in to make sure history finds grid
    header_info%grid_id = grid_id

    if (allocated(header_info%hdims)) then
      ! This shouldn't happen but, no harm, no foul
      deallocate(header_info%hdims)
    end if

    if (associated(header_info%lon_varid)) then
      ! This could be a sign of bad memory management
      call endrun('CAM_GRID_WRITE_ATTR: lon_varid should be NULL')
    end if
    if (associated(header_info%lat_varid)) then
      ! This could be a sign of bad memory management
      call endrun('CAM_GRID_WRITE_ATTR: lat_varid should be NULL')
    end if

    ! Only write this grid if not already defined
    if (cam_grids(gridind)%attrs_defined) then
      ! We need to fill out the hdims info for this grid
      call cam_grids(gridind)%find_dimids(File, dimids)
      if (dimids(2) < 0) then
        allocate(header_info%hdims(1))
        header_info%hdims(1) = dimids(1)
      else
        allocate(header_info%hdims(2))
        header_info%hdims(1:2) = dimids(1:2)
      end if
    else
      ! Write the horizontal coord attributes first so that we have the dims
      call cam_grids(gridind)%lat_coord%write_attr(File, dimids(2))
      call cam_grids(gridind)%lon_coord%write_attr(File, dimids(1))

      if (dimids(2) == dimids(1)) then
        allocate(header_info%hdims(1))
      else
        allocate(header_info%hdims(2))
        header_info%hdims(2) = dimids(2)
      end if
      header_info%hdims(1) = dimids(1)

      ! We will handle errors for this routine
      !!XXgoldyXX: This hack should be replaced with the PIO interface
      !err_handling = File%iosystem%error_handling !! Hack
      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

      attrPtr => cam_grids(gridind)%attributes
      do while (associated(attrPtr))
!!XXgoldyXX: Is this not working in PGI?
!      attr => attrPtr%getAttr()
        attr => attrPtr%attr
        call attr%write_attr(File)
!!XXgoldyXX: Is this not working in PGI?
!      attrPtr => attrPtr%getNext()
        attrPtr => attrPtr%next
      end do

      ! Back to previous I/O error handling
      call pio_seterrorhandling(File, err_handling)

      cam_grids(gridind)%attrs_defined = .true.
    end if

  end subroutine cam_grid_write_attr

  subroutine write_cam_grid_val_0d_int(attr, File)
    use pio, only: file_desc_t, pio_inq_varid, pio_put_var

    ! Dummy arguments
    class(cam_grid_attribute_0d_int_t), intent(inout) :: attr
    type(file_desc_t),                  intent(inout) :: File

    ! Local variables
    character(len=120)               :: errormsg
    integer                          :: ierr

    ! We only write this var if it is a variable
    if (associated(attr%vardesc)) then
      ierr = pio_put_var(File, attr%vardesc, attr%ival)
      call cam_pio_handle_error(ierr, 'Error writing value in write_cam_grid_val_0d_int')
      deallocate(attr%vardesc)
      nullify(attr%vardesc)
    end if

  end subroutine write_cam_grid_val_0d_int

  subroutine write_cam_grid_val_0d_char(attr, File)
    use pio, only: file_desc_t

    ! Dummy arguments
    class(cam_grid_attribute_0d_char_t), intent(inout) :: attr
    type(file_desc_t),                   intent(inout) :: File

    ! This subroutine is a stub because global attributes are written
    ! in define mode
    return
  end subroutine write_cam_grid_val_0d_char

  subroutine write_cam_grid_val_1d_int(attr, File)
    use pio,           only: file_desc_t, pio_put_var, pio_int,               &
         pio_inq_varid, pio_write_darray, io_desc_t, pio_freedecomp
    use cam_pio_utils, only: cam_pio_newdecomp

    ! Dummy arguments
    class(cam_grid_attribute_1d_int_t), intent(inout) :: attr
    type(file_desc_t),                  intent(inout) :: File

    ! Local variables
    character(len=120)               :: errormsg
    integer                          :: ierr
    type(io_desc_t), pointer         :: iodesc => NULL()

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is not associated, another grid write the values
    if (associated(attr%vardesc)) then
      ! Write out the values for this dimension variable
      if (associated(attr%map)) then
        ! This is a distributed variable, use pio_write_darray
        allocate(iodesc)
        call cam_pio_newdecomp(iodesc, (/attr%dimsize/), attr%map, pio_int)
        call pio_write_darray(File, attr%vardesc, iodesc, attr%values, ierr, -900)
        call pio_freedecomp(File, iodesc)
        deallocate(iodesc)
        nullify(iodesc)
      else
        ! This is a local variable, pio_put_var should work fine
        ierr = pio_put_var(File, attr%vardesc, attr%values)
      end if
      call cam_pio_handle_error(ierr, 'Error writing variable values in write_cam_grid_val_1d_int')
      deallocate(attr%vardesc)
      nullify(attr%vardesc)
    end if

  end subroutine write_cam_grid_val_1d_int

  subroutine write_cam_grid_val_1d_r8(attr, File)
    use pio,           only: file_desc_t, pio_put_var, pio_double,            &
         pio_inq_varid, pio_write_darray, io_desc_t, pio_freedecomp
    use cam_pio_utils, only: cam_pio_newdecomp

    ! Dummy arguments
    class(cam_grid_attribute_1d_r8_t), intent(inout) :: attr
    type(file_desc_t),                 intent(inout) :: File

    ! Local variables
    character(len=120)               :: errormsg
    integer                          :: ierr
    type(io_desc_t), pointer         :: iodesc => NULL()

    ! Since more than one grid can share an attribute, assume that if the
    ! vardesc is not associated, another grid write the values
    if (associated(attr%vardesc)) then
      ! Write out the values for this dimension variable
      if (associated(attr%map)) then
        ! This is a distributed variable, use pio_write_darray
        allocate(iodesc)
        call cam_pio_newdecomp(iodesc, (/attr%dimsize/), attr%map, pio_double)
        call pio_write_darray(File, attr%vardesc, iodesc, attr%values, ierr, -900._r8)
        call pio_freedecomp(File, iodesc)
        deallocate(iodesc)
        nullify(iodesc)
      else
        ! This is a local variable, pio_put_var should work fine
        ierr = pio_put_var(File, attr%vardesc, attr%values)
      end if
      call cam_pio_handle_error(ierr, 'Error writing variable values in write_cam_grid_val_1d_r8')
      deallocate(attr%vardesc)
      nullify(attr%vardesc)
    end if

  end subroutine write_cam_grid_val_1d_r8

  subroutine cam_grid_write_var(File, grid_id)
   use pio, only: file_desc_t, pio_bcast_error, pio_seterrorhandling

    ! Dummy arguments
    type(file_desc_t), intent(inout)     :: File        ! PIO file Handle
    integer,           intent(in)        :: grid_id

    ! Local variables
    integer                              :: gridind
    integer                              :: err_handling
    class(cam_grid_attribute_t), pointer :: attr
    type(cam_grid_attr_ptr_t),   pointer :: attrPtr

    gridind = get_cam_grid_index(grid_id)
    ! Only write if not already done
    if (cam_grids(gridind)%attrs_defined) then
      ! Write the horizontal coorinate values
      call cam_grids(gridind)%lon_coord%write_var(File)
      call cam_grids(gridind)%lat_coord%write_var(File)

      ! We will handle errors for this routine
      !!XXgoldyXX: This hack should be replaced with the PIO interface
      !err_handling = File%iosystem%error_handling !! Hack
      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

      ! Write out the variable values for each grid attribute
      attrPtr => cam_grids(gridind)%attributes
      do while (associated(attrPtr))
!!XXgoldyXX: Is this not working in PGI?
!      attr => attrPtr%getAttr()
        attr => attrPtr%attr
        call attr%write_val(File)
!!XXgoldyXX: Is this not working in PGI?
!      attrPtr => attrPtr%getNext()
        attrPtr => attrPtr%next
      end do

      ! Back to previous I/O error handling
      call pio_seterrorhandling(File, err_handling)

      cam_grids(gridind)%attrs_defined = .false.
    end if

  end subroutine cam_grid_write_var

  logical function cam_grid_block_indexed(this)
    class(cam_grid_t)                         :: this

    cam_grid_block_indexed = this%block_indexed
  end function cam_grid_block_indexed

  logical function cam_grid_unstructured(this)
    class(cam_grid_t)                         :: this

    cam_grid_unstructured = this%unstructured
  end function cam_grid_unstructured

  !---------------------------------------------------------------------------
  !
  !  cam_grid_get_dims: Return the dimensions of the grid
  !                For lon/lat grids, this is (nlon, nlat)
  !                For unstructured grids, this is (ncols, 1)
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_get_dims(this, dims)
    ! Dummy arguments
    class(cam_grid_t)                :: this
    integer,           intent(inout) :: dims(2)

    if (this%is_unstructured()) then
      call this%lon_coord%get_coord_len(dims(1))
      dims(2) = 1
    else
      call this%lon_coord%get_coord_len(dims(1))
      call this%lat_coord%get_coord_len(dims(2))
    end if

  end subroutine cam_grid_get_dims

  !---------------------------------------------------------------------------
  !
  !  cam_grid_coord_names: Return the names of the grid axes
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_coord_names(this, name1, name2)
    ! Dummy arguments
    class(cam_grid_t)                :: this
    character(len=*),  intent(out)   :: name1
    character(len=*),  intent(out)   :: name2

    call this%lon_coord%get_coord_name(name1)
    call this%lat_coord%get_coord_name(name2)

  end subroutine cam_grid_coord_names

  !---------------------------------------------------------------------------
  !
  !  cam_grid_dim_names: Return the names of the dimensions of the grid axes.
  !        Note that these may be the same
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_dim_names(this, name1, name2)
    ! Dummy arguments
    class(cam_grid_t)                :: this
    character(len=*),  intent(out)   :: name1
    character(len=*),  intent(out)   :: name2

    call this%lon_coord%get_dim_name(name1)
    call this%lat_coord%get_dim_name(name2)

  end subroutine cam_grid_dim_names

  !---------------------------------------------------------------------------
  !
  !  cam_grid_dimensions_id: Return the dimensions of the grid
  !                For lon/lat grids, this is (nlon, nlat)
  !                For unstructured grids, this is (ncols, 1)
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_dimensions_id(gridid, dims, rank)
    ! Dummy arguments
    integer,           intent(in)     :: gridid
    integer,           intent(inout)  :: dims(2)
    integer, optional, intent(out)    :: rank

    ! Local variables
    integer                           :: index
    character(len=max_hcoordname_len) :: dname1, dname2
    character(len=120)                :: errormsg

    index = get_cam_grid_index(gridid)
    if (index < 0) then
      write(errormsg, *) 'No CAM grid with ID =', gridid
      call endrun(errormsg)
    else
      call cam_grids(index)%coord_lengths(dims)
    end if
    if (present(rank)) then
      call cam_grids(index)%dim_names(dname1, dname2)
      if (trim(dname1) == trim(dname2)) then
        rank = 1
      else
        rank = 2
      end if
    end if

  end subroutine cam_grid_dimensions_id

  !---------------------------------------------------------------------------
  !
  !  cam_grid_dimensions_name: Return the dimensions of the grid
  !                For lon/lat grids, this is (nlon, nlat)
  !                For unstructured grids, this is (ncols, 1)
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_dimensions_name(gridname, dims, rank)
    ! Dummy arguments
    character(len=*),  intent(in)     :: gridname
    integer,           intent(inout)  :: dims(2)
    integer, optional, intent(out)    :: rank

    ! Local variables
    integer                           :: gridind
    character(len=max_hcoordname_len) :: dname1, dname2
    character(len=120)                :: errormsg

    gridind = get_cam_grid_index(trim(gridname))
    if (gridind < 0) then
      write(errormsg, *) 'No CAM grid with name =', trim(gridname)
      call endrun(errormsg)
    else
      call cam_grids(gridind)%coord_lengths(dims)
    end if
    if (present(rank)) then
      call cam_grids(gridind)%dim_names(dname1, dname2)
      if (trim(dname1) == trim(dname2)) then
        rank = 1
      else
        rank = 2
      end if
    end if

  end subroutine cam_grid_dimensions_name

  !---------------------------------------------------------------------------
  !
  !  cam_grid_set_map: Set a grid's distribution map
  !             This maps the local grid elements to global file order
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_set_map(this, map, src, dest)
    use spmd_utils,      only: mpi_sum, mpi_integer, mpicom
    ! Dummy arguments
    class(cam_grid_t)                      :: this
    integer(iMap),     pointer             :: map(:,:)
    integer,                    intent(in) :: src(2)   ! decomp info
    integer,                    intent(in) :: dest(2)  ! Standard dim(s) in file

    ! Local variables
    integer                                :: dims(2)
    integer                                :: dstrt, dend
    integer                                :: gridlen, gridloc, ierr

    ! Check to make sure the map meets our needs
    call this%coord_lengths(dims)
    dend = size(map, 1)
    ! We always have to have one source and one destination
    if (dest(2) > 0) then
      dstrt = dend - 1
    else
      dstrt = dend
    end if
    if ((src(2) /= 0) .and. (dstrt < 3)) then
      call endrun('cam_grid_set_map: src & dest too large for map')
    else if (dstrt < 2) then
      call endrun('cam_grid_set_map: dest too large for map')
    ! No else needed
    end if
    if (dstrt == dend) then
      gridloc = count(map(dend,:) /= 0)
    else
      gridloc = count((map(dstrt,:) /= 0) .and. (map(dend,:) /= 0))
    end if
    call MPI_Allreduce(gridloc, gridlen, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (gridlen /= product(dims)) then
      call endrun('cam_grid_set_map: Bad map size for '//trim(this%name))
    else
      if (.not. associated(this%map)) then
        allocate(this%map)
      end if
      call this%map%init(map, this%unstructured, src, dest)
    end if
  end subroutine cam_grid_set_map

  !---------------------------------------------------------------------------
  !
  !  cam_grid_local_size: return the local size of a 2D array on this grid
  !
  !---------------------------------------------------------------------------
  integer function cam_grid_local_size(this)

    ! Dummy argument
    class(cam_grid_t)                         :: this

    ! Local variable
    character(len=128)                        :: errormsg

    if (.not. associated(this%map)) then
      write(errormsg, *) 'Grid, '//trim(this%name)//', has no map'
      call endrun('cam_grid_local_size: '//trim(errormsg))
    else
      cam_grid_local_size = this%map%num_elem()
    end if

  end function cam_grid_local_size

  !---------------------------------------------------------------------------
  !
  !  cam_grid_get_lon_lat: Find the latitude and longitude for a given
  !                        grid map index. Note if point is not mapped
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_get_lon_lat(this, index, lon, lat, isMapped)

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    integer,                   intent(in)    :: index
    real(r8),                  intent(out)   :: lon
    real(r8),                  intent(out)   :: lat
    logical,                   intent(out)   :: isMapped

    ! Local variables
    integer                                  :: latindex, lonindex
    character(len=*), parameter              :: subname = "cam_grid_get_lon_lat"

    if (this%block_indexed) then
      lonindex = index
      latindex = index
      isMapped = this%map%is_mapped(index)
    else
      call this%map%coord_vals(index, lonindex, latindex, isMapped)
    end if

    !!XXgoldyXX: May be able to relax all the checks
    if ( (latindex < LBOUND(this%lat_coord%values, 1)) .or.                   &
         (latindex > UBOUND(this%lat_coord%values, 1))) then
      call endrun(trim(subname)//": index out of range for latvals")
    else
      lat = this%lat_coord%values(latindex)
    end if

    if ( (lonindex < LBOUND(this%lon_coord%values, 1)) .or.                   &
         (lonindex > UBOUND(this%lon_coord%values, 1))) then
      call endrun(trim(subname)//": index out of range for lonvals")
    else
      lon = this%lon_coord%values(lonindex)
    end if

  end subroutine cam_grid_get_lon_lat

  !---------------------------------------------------------------------------
  !
  !  cam_grid_find_src_dims: Find the correct src array dims for this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_find_src_dims(this, field_dnames, src_out)
    ! Dummy arguments
    class(cam_grid_t)                         :: this
    character(len=*),           intent(in)    :: field_dnames(:)
    integer,           pointer                :: src_out(:)

    ! Local variables
    integer                                   :: i, j
    integer                                   :: num_coords
    character(len=max_hcoordname_len)         :: coord_dimnames(2)

    call this%dim_names(coord_dimnames(1), coord_dimnames(2))
    if (associated(src_out)) then
      deallocate(src_out)
      nullify(src_out)
    end if
    if (trim(coord_dimnames(1)) == trim(coord_dimnames(2))) then
      num_coords = 1
    else
      num_coords = 2
    end if
    allocate(src_out(2)) ! Currently, all cases have two source dims
    do i = 1, num_coords
      do j = 1, size(field_dnames)
        if (trim(field_dnames(j)) == trim(coord_dimnames(i))) then
          src_out(i) = j
        end if
      end do
    end do
    if (num_coords < 2) then
      src_out(2) = -1  ! Assume a block structure for unstructured grids
    end if

  end subroutine cam_grid_find_src_dims

  !---------------------------------------------------------------------------
  !
  !  cam_grid_find_dest_dims: Find the correct file array dims for this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_find_dest_dims(this, file_dnames, dest_out)
    ! Dummy arguments
    class(cam_grid_t)                         :: this
    character(len=*),           intent(in)    :: file_dnames(:)
    integer,           pointer                :: dest_out(:)

    ! Local variables
    integer                                   :: i, j
    integer                                   :: num_coords
    character(len=max_hcoordname_len)         :: coord_dimnames(2)

    call this%dim_names(coord_dimnames(1), coord_dimnames(2))
    if (associated(dest_out)) then
      deallocate(dest_out)
      nullify(dest_out)
    end if
    if (trim(coord_dimnames(1)) == trim(coord_dimnames(2))) then
      num_coords = 1
    else
      num_coords = 2
    end if
    allocate(dest_out(num_coords))
    dest_out = 0
    do i = 1, num_coords
      do j = 1, size(file_dnames)
        if (trim(file_dnames(j)) == trim(coord_dimnames(i))) then
          dest_out(i) = j
        end if
      end do
    end do

  end subroutine cam_grid_find_dest_dims

  !---------------------------------------------------------------------------
  !
  !  cam_grid_get_pio_decomp: Find or create a PIO decomp on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_get_pio_decomp(this, field_lens, file_lens, dtype,      &
       iodesc, field_dnames, file_dnames)
    use pio,           only: io_desc_t
    use cam_pio_utils, only: cam_pio_get_decomp, calc_permutation

    ! Dummy arguments
    class(cam_grid_t)                         :: this
    integer,                    intent(in)    :: field_lens(:)
    integer,                    intent(in)    :: file_lens(:)
    integer,                    intent(in)    :: dtype
    type(io_desc_t), pointer,   intent(out)   :: iodesc
    character(len=*), optional, intent(in)    :: field_dnames(:)
    character(len=*), optional, intent(in)    :: file_dnames(:)

    ! Local variables
    integer,              pointer             :: src_in(:)
    integer,              pointer             :: dest_in(:)
    integer, allocatable                      :: permutation(:)
    logical                                   :: is_perm
    character(len=128)                        :: errormsg

    nullify(src_in)
    nullify(dest_in)
    is_perm = .false.
    if (.not. associated(this%map)) then
      write(errormsg, *) 'Grid, '//trim(this%name)//', has no map'
      call endrun('cam_grid_get_pio_decomp: '//trim(errormsg))
    else
      if (present(field_dnames)) then
        call this%find_src_dims(field_dnames, src_in)
      end if
      if (present(file_dnames)) then
        call this%find_dest_dims(file_dnames, dest_in)
      end if
      if (present(file_dnames) .and. present(field_dnames)) then
        ! This only works if the arrays are the same size
        if (size(file_dnames) == size(field_dnames)) then
          allocate(permutation(size(file_dnames)))
          call calc_permutation(file_dnames, field_dnames, permutation, is_perm)
        end if
      end if
      ! Call cam_pio_get_decomp with the appropriate options
      if (present(field_dnames) .and. present(file_dnames)) then
        if (is_perm) then
          call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype,       &
               this%map, field_dist_in=src_in, file_dist_in=dest_in,          &
               permute=permutation)
        else
          call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype,       &
               this%map, field_dist_in=src_in, file_dist_in=dest_in)
        end if
      else if (present(field_dnames)) then
        call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype,         &
             this%map, field_dist_in=src_in)
      else if (present(file_dnames)) then
        call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype,         &
             this%map, file_dist_in=dest_in)
      else
        call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype, this%map)
      end if
    end if
    if (associated(src_in)) then
      deallocate(src_in)
      nullify(src_in)
    end if
    if (associated(dest_in)) then
      deallocate(dest_in)
      nullify(dest_in)
    end if
    if (allocated(permutation)) then
      deallocate(permutation)
    end if

  end subroutine cam_grid_get_pio_decomp

  !-------------------------------------------------------------------------------
  !
  !  cam_grid_find_dimids: Find the dimension NetCDF IDs on <File> for this grid
  !
  !-------------------------------------------------------------------------------
  subroutine cam_grid_find_dimids(this, File, dimids)
    use pio, only: file_desc_t, pio_noerr, pio_inq_dimid
    use pio, only: pio_seterrorhandling, pio_bcast_error

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(out)   :: dimids(:)

    ! Local vaariables
    integer                                  :: dsize, ierr
    integer                                  :: err_handling
    character(len=max_hcoordname_len)        :: dimname1, dimname2

    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = File%iosystem%error_handling !! Hack
    call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

    call this%dim_names(dimname1, dimname2)
    if (size(dimids) < 1) then
      call endrun('CAM_GRID_FIND_DIMIDS: dimids must have positive size')
    end if
    dimids = -1
    ! Check the first dimension
    ierr = pio_inq_dimid(File, trim(dimname1), dimids(1))
    if(ierr /= PIO_NOERR) then
      call endrun('CAM_GRID_FIND_DIMIDS: '//trim(this%name)//' dimension, '//trim(dimname1)//', does not exist on file')
    end if
    if (trim(dimname1) /= trim(dimname2)) then
      ! Structured grid, find second dimid
      if (size(dimids) < 2) then
        call endrun('CAM_GRID_FIND_DIMIDS: dimids too small for '//trim(this%name))
      end if
      ierr = pio_inq_dimid(File, trim(dimname2), dimids(2))
      if(ierr /= PIO_NOERR) then
        call endrun('CAM_GRID_FIND_DIMIDS: '//trim(this%name)//' dimension, '//trim(dimname2)//', does not exist on file')
      end if
    end if
    
    ! Back to whatever error handling was running before this routine
    call pio_seterrorhandling(File, err_handling)

  end subroutine cam_grid_find_dimids

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_2d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_2d_int(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray, PIO_INT
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_INT, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_2d_int: Error reading variable')
  end subroutine cam_grid_read_darray_2d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_3d_int: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_3d_int(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray, PIO_INT
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_INT, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_3d_int: Error reading variable')
  end subroutine cam_grid_read_darray_3d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_2d_double: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_2d_double(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray
    use pio,           only: PIO_DOUBLE
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_DOUBLE, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_2d_double: Error reading variable')
  end subroutine cam_grid_read_darray_2d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_3d_double: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_3d_double(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray
    use pio,           only: PIO_DOUBLE
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_DOUBLE, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_3d_double: Error reading variable')
  end subroutine cam_grid_read_darray_3d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_2d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_2d_real(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray
    use pio,           only: PIO_REAL
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(out)   :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_REAL, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_2d_real: Error reading variable')
  end subroutine cam_grid_read_darray_2d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_read_darray_3d_real: Read a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_read_darray_3d_real(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t, pio_read_darray
    use pio,           only: PIO_REAL
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(out)   :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_REAL, this%map)
    call pio_read_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_read_darray_2d_: Error reading variable')
  end subroutine cam_grid_read_darray_3d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_2d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_2d_int(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_INT

    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_INT, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_2d_int: Error writing variable')
  end subroutine cam_grid_write_darray_2d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_3d_int: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_3d_int(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_INT
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    integer,                   intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_INT, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_3d_int: Error writing variable')
  end subroutine cam_grid_write_darray_3d_int

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_2d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_2d_double(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_DOUBLE, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_2d_double: Error writing variable')
  end subroutine cam_grid_write_darray_2d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_3d_double: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_3d_double(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r8),                  intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_DOUBLE, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_3d_double: Error writing variable')

  end subroutine cam_grid_write_darray_3d_double

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_2d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_2d_real(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_REAL
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(in)    :: hbuf(:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_REAL, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_2d_real: Error writing variable')
  end subroutine cam_grid_write_darray_2d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_write_darray_3d_real: Write a variable defined on this grid
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_write_darray_3d_real(this, File, adims, fdims, hbuf, varid)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_REAL
    use cam_pio_utils, only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: fdims(:)
    real(r4),                  intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),          intent(inout) :: varid

    ! Local variables
    type(io_desc_t), pointer                 :: iodesc => NULL()
    integer                                  :: ierr

    call cam_pio_get_decomp(iodesc, adims, fdims, PIO_REAL, this%map)
    call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_write_darray_3d_real: Error writing variable')
  end subroutine cam_grid_write_darray_3d_real

  !---------------------------------------------------------------------------
  !
  !  cam_grid_get_patch_mask: Compute a map which is defined for locations
  !       within the input patch.
  !
  !---------------------------------------------------------------------------
  subroutine cam_grid_get_patch_mask(this, lonl, lonu, latl, latu, patch)
    use spmd_utils,      only: mpi_min, mpi_real8, mpicom
    use physconst,       only: pi

    ! Dummy arguments
    class(cam_grid_t)                        :: this
    real(r8),                  intent(in)    :: lonl, lonu ! Longitude bounds
    real(r8),                  intent(in)    :: latl, latu ! Latitude bounds
    type(cam_grid_patch_t),    intent(inout) :: patch

    ! Local arguments
    real(r8)                         :: mindist
    real(r8)                         :: minlondist
    real(r8)                         :: dist    ! Test distance
    real(r8)                         :: londeg, latdeg
    real(r8)                         :: lon,    lat
    real(r8)                         :: londeg_min, latdeg_min
    real(r8)                         :: lonmin, lonmax, latmin, latmax
    integer                          :: minind  ! Location of closest point
    integer                          :: mapind  ! Grid map index
    integer                          :: latind, lonind
    integer                          :: ierr    ! For MPI calls
    integer                          :: gridloc ! local size of grid
    logical                          :: unstructured ! grid type
    logical                          :: findClosest  ! .false. == patch output
    logical                          :: isMapped     ! .true. iff point in map
    logical                          :: closer       ! Test for new closest pt

    real(r8), parameter              :: maxangle = pi / 4.0_r8
    real(r8), parameter              :: deg2rad = pi / 180.0_r8
    real(r8), parameter              :: maxlat = 89.5_r8

    if (.not. associated(this%map)) then
      call endrun('cam_grid_get_patch_mask: Grid, '//trim(this%name)//', has no map')
    end if
    gridloc = this%map%num_elem()
    if (associated(patch%mask)) then
      if (patch%mask%num_elem() /= gridloc) then
        ! The mask needs to be the same size as the map
        call endrun('cam_grid_get_patch_mask: mask is incorrect size')
        ! No else, just needed a check
        ! In particular, we are not zeroing the mask since multiple calls with
        ! the same mask can be used for collected-column output
        ! NB: Compacting the mask must be done after all calls (for a
        !     particular mask) to this function.
      end if
    else
      if (associated(patch%latmap)) then
        call endrun('cam_grid_get_patch_mask: unallocated patch has latmap')
      end if
      if (associated(patch%lonmap)) then
        call endrun('cam_grid_get_patch_mask: unallocated patch has lonmap')
      end if
      call patch%set_patch(lonl, lonu, latl, latu, this%id, this%map)
      if (patch%mask%num_elem() /= gridloc) then
        ! Basic check to make sure the copy worked
        call endrun('cam_grid_get_patch_mask: grid map is invalid')
      end if
      call patch%mask%clear()
      ! Set up the lat/lon maps
      if (associated(this%lat_coord%values)) then
        allocate(patch%latmap(LBOUND(this%lat_coord%values, 1):UBOUND(this%lat_coord%values, 1)))
        patch%latmap = 0
      else
        nullify(patch%latmap)
      end if
      if (associated(this%lon_coord%values)) then
        allocate(patch%lonmap(LBOUND(this%lon_coord%values, 1):UBOUND(this%lon_coord%values, 1)))
        patch%lonmap = 0
      else
        nullify(patch%lonmap)
      end if
    end if

    ! We have to iterate through each grid point to check
    ! We have four cases, structured vs. unstructured grid *
    !   patch area vs. closest column
    ! Note that a 1-d patch 'area' is not allowed for unstructured grids
    unstructured = this%is_unstructured()
    findClosest = .false.
    ! Make sure our search items are in order
    lonmin = min(lonl, lonu)
    lonmax = max(lonl, lonu)
    latmin = min(latl, latu)
    latmax = max(latl, latu)
    if (lonl == lonu) then
      if (latl == latu) then
        findClosest = .true.
      else if (unstructured) then
        call endrun('cam_grid_get_patch_mask: 1-D patch (lon) not allowed for unstructured grids')
      else
        ! Find closest lon line to lonu
        ! This is a lat lon grid so it should have coordinate axes
        lonmin = 365.0_r8
        mindist = 365.0_r8
        if (associated(this%lon_coord%values)) then
          do lonind = LBOUND(this%lon_coord%values, 1), UBOUND(this%lon_coord%values, 1)
            dist = abs(this%lon_coord%values(lonind) - lonu)
            if (dist < mindist) then
              lonmin = this%lon_coord%values(lonind)
              mindist = dist
            end if
          end do
        end if
        ! Get the global minimum
        dist = mindist
        call MPI_allreduce(dist, mindist, 1, mpi_real8, mpi_min, mpicom, ierr)
        if (dist == mindist) then
          ! We have a ringer so use only that longitude
          lonmax = lonmin
        else
          ! We don't have a minimum dist so count no points
          lonmax = lonmin - 1.0_r8
        end if
      end if
    else if (latl == latu) then
      if (unstructured) then
        call endrun('cam_grid_get_patch_mask: 1-D patch (lat) not allowed for unstructured grids')
      else
        ! Find closest lat line to latu
        ! This is a lat lon grid so it should have coordinate axes
        latmin = 91.0_r8
        mindist = 181.0_r8
        if (associated(this%lat_coord%values)) then
          do latind = LBOUND(this%lat_coord%values, 1), UBOUND(this%lat_coord%values, 1)
            dist = abs(this%lat_coord%values(latind) - latl)
            if (dist < mindist) then
              latmin = this%lat_coord%values(latind)
              mindist = dist
            end if
          end do
        end if
        ! Get the global minimum
        dist = mindist
        call MPI_allreduce(dist, mindist, 1, mpi_real8, mpi_min, mpicom, ierr)
        if (dist == mindist) then
          ! We have a ringer so use only that latitude
          latmax = latmin
        else
          ! We don't have a minimum dist so count no points
          latmax = latmin - 1.0_r8
        end if
      end if
    end if

    ! Convert to radians
    lonmin = lonmin * deg2rad
    lonmax = lonmax * deg2rad
    latmin = latmin * deg2rad
    latmax = latmax * deg2rad
    ! Loop through all the local grid elements and find the closest match
    ! (or all matches depending on the value of findClosest)
    mapind = 1
    minind = -1
    londeg_min = 361.0_r8
    latdeg_min = 91.0_r8
    mindist = maxangle * 2.0_r8
    minlondist = 361.0_r8
    do mapind = 1, patch%mask%num_elem()
      call this%get_lon_lat(mapind, londeg, latdeg, isMapped)
      lon = londeg * deg2rad
      lat = latdeg * deg2rad
      if (isMapped) then
        if (findClosest) then
          ! Use the Spherical Law of Cosines to find the great-circle distance.
          ! Might as well use the unit sphere since we just want differences
          if ( (abs(lat - latmin) <= maxangle) .and.                          &
               (abs(lon - lonmin) <= maxangle)) then
            ! maxangle could be pi but why waste all those trig functions?
            ! XXgoldyXX: What should we use for maxangle given coarse Eul grids?
            dist = acos((sin(latmin) * sin(lat)) +                            &
                 (cos(latmin) * cos(lat) * cos(lon - lonmin)))
            closer = dist < mindist
            if (abs(latdeg) > maxlat) then
              closer = closer .and. ((lon - lonmin) < minlondist)
            end if
            if (closer) then
              minind = mapind
              mindist = dist
              londeg_min = londeg
              latdeg_min = latdeg
              minlondist = (londeg_min * deg2rad) - lonmin
            end if
          end if
        else
          if ( (latmin <= lat) .and. (lat <= latmax) .and.                    &
               (lonmin <= lon) .and. (lon <= lonmax)) then
            if (patch%mask%num_elem() >= mapind) then
              call patch%mask%copy_elem(this%map, mapind)
              if (this%block_indexed) then
                 call this%map%coord_dests(mapind, lonind, latind)
                 if (latind == 0) then
                   latind = lonind
                 end if
                 if (associated(patch%latmap)) then
                   patch%latmap(mapind) = latind
                 end if
                 if (associated(patch%lonmap)) then
                   patch%lonmap(mapind) = lonind
                 end if
              else
                 call this%map%coord_vals(mapind, lonind, latind)
                 if (associated(patch%latmap)) then
                   patch%latmap(latind) = latind
                 end if
                 if (associated(patch%lonmap)) then
                   patch%lonmap(lonind) = lonind
                 end if
             end if
            else
              call endrun('cam_grid_get_patch_mask: PE has patch points but mask too small')
            end if
          end if
        end if
      end if
    end do
    if (findClosest) then
      ! We need to find the minimum mindist and use only that value
      dist = mindist
      call MPI_allreduce(dist, mindist, 1, mpi_real8, mpi_min, mpicom, ierr)
      lon = minlondist
      call MPI_allreduce(lon, minlondist, 1, mpi_real8, mpi_min, mpicom, ierr)
      ! Now, only task(s) which has real minimum distance should set their mask
      ! minind test allows for no match
      if ( (dist == mindist) .and.                                            &
           ((abs(latdeg_min) <= maxlat) .or. (lon == minlondist))) then
        if (minind < 0) then
          call endrun("cam_grid_get_patch_mask: No closest point found!!")
        else
          if (patch%mask%num_elem() >= minind) then
            call patch%mask%copy_elem(this%map, minind)
            if (this%block_indexed) then
              call this%map%coord_dests(minind, lonind, latind)
              if (latind == 0) then
                latind = lonind
              end if
              patch%latmap(minind) = latind
              patch%lonmap(minind) = lonind
            else
              call this%map%coord_vals(minind, lonind, latind)
              patch%latmap(latind) = latind
              patch%lonmap(lonind) = lonind
            end if
          else
            call endrun('cam_grid_get_patch_mask: PE has patch closest point but mask too small')
          end if
        end if
      end if
    end if

  end subroutine cam_grid_get_patch_mask

  !---------------------------------------------------------------------------
  !
  !  Grid Patch functions
  !
  !---------------------------------------------------------------------------

  integer function cam_grid_patch_get_id(this) result(id)

    ! Dummy argument
    class(cam_grid_patch_t)                  :: this

    id = this%grid_id
  end function cam_grid_patch_get_id

  subroutine cam_grid_patch_get_global_size_map(this, gsize)

    ! Dummy arguments
    class(cam_grid_patch_t),   intent(in)    :: this
    integer,                   intent(out)   :: gsize

    gsize = this%global_size

  end subroutine cam_grid_patch_get_global_size_map

  subroutine cam_grid_patch_get_global_size_axes(this, latsize, lonsize)

    ! Dummy arguments
    class(cam_grid_patch_t),   intent(in)    :: this
    integer,                   intent(out)   :: latsize
    integer,                   intent(out)   :: lonsize

    latsize = this%global_lat_size
    lonsize = this%global_lon_size

  end subroutine cam_grid_patch_get_global_size_axes

  ! cam_grid_patch_get_axis_names
  !   Collect or compute unique names for the latitude and longitude axes
  !   If the grid is unstructured or col_output is .true., the column
  !     dimension name is also generated (e.g., ncol)
  subroutine cam_grid_patch_get_axis_names(this, lat_name, lon_name,          &
       col_name, col_output)

    ! Dummy arguments
    class(cam_grid_patch_t)                  :: this
    character(len=*),          intent(out)   :: lat_name
    character(len=*),          intent(out)   :: lon_name
    character(len=*),          intent(out)   :: col_name
    logical,                   intent(in)    :: col_output

    ! Local variable
    integer                                  :: index
    character(len=120)                       :: errormsg
    character(len=max_hcoordname_len)        :: grid_name
    logical                                  :: unstruct

    if (cam_grid_check(this%grid_id)) then
      index = this%grid_index()
      unstruct = cam_grids(index)%is_unstructured()
      ! Get coordinate and dim names
      call cam_grids(index)%lat_coord%get_coord_name(lat_name)
      call cam_grids(index)%lon_coord%get_coord_name(lon_name)
      if (unstruct) then
        call cam_grids(index)%lon_coord%get_dim_name(col_name)
      else
        col_name = ''
      end if
      grid_name = cam_grids(index)%name
      if (col_output .or. unstruct) then
        ! In this case, we are using collect_column_output on a lat/lon grid
        col_name = 'ncol_'//trim(grid_name)
        lat_name = trim(lat_name)//'_'//trim(grid_name)
        lon_name = trim(lon_name)//'_'//trim(grid_name)
      else
        ! Separate patch output for a lat/lon grid
        col_name = ''
        lat_name = trim(lat_name)//'_'//trim(grid_name)
        lon_name = trim(lon_name)//'_'//trim(grid_name)
      end if
    else
      write(errormsg, *) 'Bad grid ID:', this%grid_id
      call endrun('cam_grid_patch_get_axis_names: '//errormsg)
    end if

  end subroutine cam_grid_patch_get_axis_names

  subroutine cam_grid_patch_get_coord_long_name(this, axis, name)

    ! Dummy arguments
    class(cam_grid_patch_t)                  :: this
    character(len=*),          intent(in)    :: axis
    character(len=*),          intent(out)   :: name

    ! Local variable
    character(len=120)               :: errormsg
    integer                          :: index

    if (cam_grid_check(this%grid_id)) then
      index = this%grid_index()
      if (trim(axis) == 'lat') then
        call cam_grids(index)%lat_coord%get_long_name(name)
      else if (trim(axis) == 'lon') then
        call cam_grids(index)%lon_coord%get_long_name(name)
      else
        write(errormsg, *) 'Bad axis name:', axis
        call endrun('cam_grid_patch_get_coord_long_name: '//errormsg)
      end if
    else
      write(errormsg, *) 'Bad grid ID:', this%grid_id
      call endrun('cam_grid_patch_get_coord_long_name: '//errormsg)
    end if

  end subroutine cam_grid_patch_get_coord_long_name

  subroutine cam_grid_patch_get_coord_units(this, axis, units)

    ! Dummy arguments
    class(cam_grid_patch_t)                  :: this
    character(len=*),          intent(in)    :: axis
    character(len=*),          intent(out)   :: units

    ! Local variable
    character(len=120)               :: errormsg
    integer                          :: index

    if (cam_grid_check(this%grid_id)) then
      index = this%grid_index()
      if (trim(axis) == 'lat') then
        call cam_grids(index)%lat_coord%get_units(units)
      else if (trim(axis) == 'lon') then
        call cam_grids(index)%lon_coord%get_units(units)
      else
        write(errormsg, *) 'Bad axis name:', axis
        call endrun('cam_grid_patch_get_coord_units: '//errormsg)
      end if
    else
      write(errormsg, *) 'Bad grid ID:', this%grid_id
      call endrun('cam_grid_patch_get_coord_units: '//errormsg)
    end if

  end subroutine cam_grid_patch_get_coord_units

  subroutine cam_grid_patch_set_patch(this, lonl, lonu, latl, latu, id, map)

    ! Dummy arguments
    class(cam_grid_patch_t)                  :: this
    real(r8),                  intent(in)    :: lonl, lonu ! Longitude bounds
    real(r8),                  intent(in)    :: latl, latu ! Latitude bounds
    integer,                   intent(in)    :: id
    type(cam_filemap_t),       intent(in)    :: map

    this%grid_id      = id
    this%lon_range(1) = lonl
    this%lon_range(2) = lonu
    this%lat_range(1) = latl
    this%lat_range(2) = latu
    if (.not. associated(this%mask)) then
      allocate(this%mask)
    end if
    call this%mask%copy(map)
    call this%mask%new_index()

  end subroutine cam_grid_patch_set_patch

  subroutine cam_grid_patch_get_decomp(this, field_lens, file_lens, dtype,    &
       iodesc, file_dest_in)
    use pio,                only: io_desc_t
    use cam_pio_utils,      only: cam_pio_get_decomp

    ! Dummy arguments
    class(cam_grid_patch_t)                  :: this
    integer,                   intent(in)    :: field_lens(:)
    integer,                   intent(in)    :: file_lens(:)
    integer,                   intent(in)    :: dtype
    type(io_desc_t), pointer,  intent(out)   :: iodesc
    integer,         optional, intent(in)    :: file_dest_in(:)

    call cam_pio_get_decomp(iodesc, field_lens, file_lens, dtype, this%mask,  &
         file_dist_in=file_dest_in)

  end subroutine cam_grid_patch_get_decomp

  subroutine cam_grid_patch_compact(this, collected_output)
    use spmd_utils,  only: mpi_sum, mpi_integer, mpicom
    use shr_mpi_mod, only: shr_mpi_chkerr

    ! Dummy arguments
    class(cam_grid_patch_t)               :: this
    logical,         optional, intent(in) :: collected_output

    ! Local variables
    integer                               :: index ! Our grid's index
    logical                               :: unstructured

    index = this%grid_index()
    if (index > 0) then
      unstructured = cam_grids(index)%is_unstructured()
    else
      ! This is probably an error condition but someone else will catch it first
      unstructured = .false.
    end if
    call this%mask%compact(this%lonmap, this%latmap,                          &
         num_lons=this%global_lon_size, num_lats=this%global_lat_size,        &
         num_mapped=this%global_size, columnize=collected_output,             &
         dups_ok_in=unstructured)

  end subroutine cam_grid_patch_compact

  subroutine cam_grid_patch_get_active_cols(this, lchnk, active, srcdim_in)

    ! Dummy arguments
    class(cam_grid_patch_t)                    :: this
    integer,                    intent(in)     :: lchnk
    logical,                    intent(out)    :: active(:)
    integer, optional,          intent(in)     :: srcdim_in

    if (.not. associated(this%mask)) then
      call endrun('cam_grid_patch_get_active_cols: No mask')
    else
      call this%mask%active_cols(lchnk, active, srcdim_in)
    end if

  end subroutine cam_grid_patch_get_active_cols

  ! cam_grid_patch_write_vals: Write lat and lon coord values to File
  subroutine cam_grid_patch_write_vals(this, File, header_info)
    use pio,           only: file_desc_t, io_desc_t
    use pio,           only: pio_write_darray, PIO_DOUBLE
    use pio,           only: pio_initdecomp, pio_freedecomp
    use cam_pio_utils, only: cam_pio_handle_error, pio_subsystem

    ! Dummy arguments
    class(cam_grid_patch_t)                     :: this
    type(file_desc_t),            intent(inout) :: File       ! PIO file handle
    type(cam_grid_header_info_t), intent(inout) :: header_info

    ! Local variables
    type(io_desc_t)                             :: iodesc
    type(var_desc_t), pointer                   :: vdesc      => NULL()
    real(r8),         pointer                   :: coord_p(:) => NULL()
    real(r8),         pointer                   :: coord(:)   => NULL()
    integer(iMap),    pointer                   :: map(:)     => null()
    integer                                     :: field_lens(1)
    integer                                     :: file_lens(1)
    integer                                     :: ierr

    if (this%grid_id /= header_info%get_gridid()) then
      call endrun('CAM_GRID_PATCH_WRITE_VALS: Grid id mismatch')
    end if
    ! Write out lon
    if (associated(this%lonmap)) then
      field_lens(1) = size(this%lonmap, 1)
      map => this%lonmap
    else
      field_lens(1) = 0
      allocate(map(0))
    end if
    file_lens(1) = this%global_lon_size
    !! XXgoldyXX: Think about caching these decomps
    call pio_initdecomp(pio_subsystem, pio_double, file_lens, map, iodesc)
    coord_p => cam_grid_get_lonvals(this%grid_id)
    if (associated(coord_p)) then
      coord => coord_p
    else
      allocate(coord(0))
    end if
    vdesc => header_info%get_lon_varid()
    call pio_write_darray(File, vdesc, iodesc, coord, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_patch_write_vals: Error writing longitude')
    if (.not. associated(this%lonmap)) then
      deallocate(map)
      nullify(map)
    end if
    if (.not. associated(coord_p)) then
      deallocate(coord)
      nullify(coord)
    end if
    ! Write out lat
    if (associated(this%latmap)) then
      field_lens(1) = size(this%latmap, 1)
      map => this%latmap
    else
      field_lens(1) = 0
      allocate(map(0))
    end if
    file_lens(1) = this%global_lat_size
    !! XXgoldyXX: Think about caching these decomps
    call pio_initdecomp(pio_subsystem, pio_double, file_lens, map, iodesc)

    coord_p => cam_grid_get_latvals(this%grid_id)
    if (associated(coord_p)) then
      coord => coord_p
    else
      allocate(coord(0))
    end if
    vdesc => header_info%get_lat_varid()
    call pio_write_darray(File, vdesc, iodesc, coord, ierr)
    call cam_pio_handle_error(ierr, 'cam_grid_patch_write_vals: Error writing latitude')
    if (.not. associated(this%latmap)) then
      deallocate(map)
      nullify(map)
    end if
    if (.not. associated(coord_p)) then
      deallocate(coord)
      nullify(coord)
    end if
    call pio_freedecomp(File, iodesc)

  end subroutine cam_grid_patch_write_vals

  integer function cam_grid_patch_get_grid_index(this) result(index)
    ! Dummy argument
    class(cam_grid_patch_t)                  :: this

    ! Local variable
    integer                                  :: i

    index = -1
    ! Find the grid index associated with our grid_id which is a decomp
    do i = 1, cam_grid_num_grids()
      if (cam_grids(i)%id == this%grid_id) then
        index = i
        exit
      end if
    end do

  end function cam_grid_patch_get_grid_index

  subroutine cam_grid_patch_deallocate(this)
    ! Dummy argument
    class(cam_grid_patch_t)                  :: this

    if (associated(this%mask)) then
      deallocate(this%mask)
      nullify(this%mask)
    end if

  end subroutine cam_grid_patch_deallocate

  integer function cam_grid_header_info_get_gridid(this) result(id)
    ! Dummy argument
    class(cam_grid_header_info_t)           :: this

    id = this%grid_id

  end function cam_grid_header_info_get_gridid

  subroutine cam_grid_header_info_set_gridid(this, id)
    ! Dummy argument
    class(cam_grid_header_info_t)            :: this
    integer,                      intent(in) :: id

    this%grid_id = id

  end subroutine cam_grid_header_info_set_gridid

  subroutine cam_grid_header_info_set_hdims(this, hdim1, hdim2)
    ! Dummy arguments
    class(cam_grid_header_info_t)                :: this
    integer,                       intent(in)    :: hdim1
    integer, optional,             intent(in)    :: hdim2

    ! Local variables
    integer                                      :: hdsize

    if (present(hdim2)) then
      hdsize = 2
    else
      hdsize = 1
    end if

    if (allocated(this%hdims)) then
      ! This can happen, for instance on opening a new version of the file
      if (size(this%hdims) /= hdsize) then
        call endrun('cam_grid_header_info_set_hdims: hdims is wrong size')
      end if
    else
      allocate(this%hdims(hdsize))
    end if
    this%hdims(1) = hdim1
    if (present(hdim2)) then
      this%hdims(2) = hdim2
    end if

  end subroutine cam_grid_header_info_set_hdims

  integer function cam_grid_header_info_num_hdims(this) result(num)
    ! Dummy argument
    class(cam_grid_header_info_t)           :: this

    if (allocated(this%hdims)) then
      num = size(this%hdims)
    else
      num = 0
    end if

  end function cam_grid_header_info_num_hdims

  integer function cam_grid_header_info_hdim(this, index) result(id)
    ! Dummy arguments
    class(cam_grid_header_info_t)               :: this
    integer,                      intent(in)    :: index

    ! Local variable
    character(len=120)               :: errormsg

    if (allocated(this%hdims)) then
      if ((index >= 1) .and. (index <= size(this%hdims))) then
        id = this%hdims(index)
      else
        write(errormsg, '(a,i0,a)') 'Index out of range, (',index,')'
        call endrun('cam_grid_header_info_hdim: '//errormsg)
      end if
    else
      write(errormsg, '(a)') 'No hdims allocated'
      call endrun('cam_grid_header_info_hdim: '//errormsg)
    end if

  end function cam_grid_header_info_hdim

  subroutine cam_grid_header_info_set_varids(this, lon_varid, lat_varid)

    ! Dummy arguments
    class(cam_grid_header_info_t)             :: this
    type(var_desc_t),              pointer    :: lon_varid
    type(var_desc_t),              pointer    :: lat_varid

    if (associated(this%lon_varid)) then
      deallocate(this%lon_varid)
      nullify(this%lon_varid)
    end if
    this%lon_varid => lon_varid
    if (associated(this%lat_varid)) then
      deallocate(this%lat_varid)
      nullify(this%lat_varid)
    end if
    this%lat_varid => lat_varid

  end subroutine cam_grid_header_info_set_varids

  function cam_grid_header_info_lon_varid(this) result(id)

    ! Dummy arguments
    class(cam_grid_header_info_t)               :: this
    type(var_desc_t),   pointer                 :: id

    id => this%lon_varid

  end function cam_grid_header_info_lon_varid

  function cam_grid_header_info_lat_varid(this) result(id)

    ! Dummy arguments
    class(cam_grid_header_info_t)               :: this
    type(var_desc_t),   pointer                 :: id

    id => this%lat_varid

  end function cam_grid_header_info_lat_varid

  subroutine cam_grid_header_info_deallocate(this)
    ! Dummy argument
    class(cam_grid_header_info_t)           :: this

    this%grid_id = -1
    if (allocated(this%hdims)) then
      deallocate(this%hdims)
    end if
    if (associated(this%lon_varid)) then
      deallocate(this%lon_varid)
      nullify(this%lon_varid)
    end if
    if (associated(this%lat_varid)) then
      deallocate(this%lat_varid)
      nullify(this%lat_varid)
    end if

  end subroutine cam_grid_header_info_deallocate

end module cam_grid_support
