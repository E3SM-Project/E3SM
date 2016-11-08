module cam_history_support

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  cam_history_support is used by cam_history as well as by the dycores
!!    (for vertical coordinate and "mdim" support). Some parameters are
!!    also referenced by cam_grid_support (although those could be copied
!!    if necessary).
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use shr_kind_mod,     only: r8=>shr_kind_r8, shr_kind_cl
  use shr_sys_mod,      only: shr_sys_flush
  use pio,              only: var_desc_t, file_desc_t
  use cam_abortutils,   only: endrun
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  use cam_grid_support, only: cam_grid_patch_t, cam_grid_header_info_t
  use cam_grid_support, only: max_hcoordname_len
  use cam_pio_utils,    only: cam_pio_handle_error

  implicit none
  private
  save

  integer, parameter, public :: max_string_len = 256   ! Length of strings
  integer, parameter, public :: max_chars = shr_kind_cl         ! max chars for char variables
  integer, parameter, public :: fieldname_len = 24   ! max chars for field name
  integer, parameter, public :: fieldname_suffix_len =  3 ! length of field name suffix ("&IC")
  integer, parameter, public :: fieldname_lenp2      = fieldname_len + 2 ! allow for extra characters
  ! max_fieldname_len = max chars for field name (including suffix)
  integer, parameter, public :: max_fieldname_len    = fieldname_len + fieldname_suffix_len

  integer, parameter, public :: pflds  = 1000      ! max number of fields for namelist entries fincl and fexcl
                                                   ! also used in write restart
  integer, parameter, public :: ptapes = 12        ! max number of tapes

  ! A special symbol for declaring a field which has no vertical or
  ! non-grid dimensions. It is here (rather than cam_history) so that it
  ! be checked by add_hist_coord
  character(len=10), parameter, public :: horiz_only = 'horiz_only'

  type, public :: history_patch_t
    character(len=max_chars)          :: namelist_entry =  ''
    ! lon_axis_name and lat_axis_name are not used if collected_output = .true.
    character(len=max_fieldname_len)  :: lon_axis_name       =  ''
    character(len=max_fieldname_len)  :: lat_axis_name       =  ''
    logical                           :: collected_output
    ! There is one patch for every grid and one header_info for every patch
    type(cam_grid_patch_t),        pointer :: patches(:)     => NULL()
    type (cam_grid_header_info_t), pointer :: header_info(:) => NULL()
  contains
    procedure                   :: write_attrs  => history_patch_write_attrs
    procedure                   :: write_vals   => history_patch_write_vals
    procedure                   :: field_name   => history_patch_field_name
    procedure                   :: num_hdims    => history_patch_num_hdims
    procedure                   :: get_var_data => history_patch_get_var_data
    procedure                   :: write_var    => history_patch_write_var
    procedure                   :: compact      => history_patch_compact
    procedure                   :: active_cols  => history_patch_active_cols
    procedure                   :: deallocate   => history_patch_deallocate
  end type history_patch_t

!
! dim_index_2d, dim_index_3d: 2-D & 3-D dimension index lower & upper bounds
!
  type, public :: dim_index_2d        ! 2-D dimension index
     integer :: beg1, end1            ! lower & upper bounds of 1st dimension
     integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
   contains
     procedure :: dim_sizes_2d  => dim_index_2d_dim_sizes_2d
     procedure :: dim_sizes_arr => dim_index_2d_dim_size_arr
     generic   :: dim_sizes     => dim_sizes_arr, dim_sizes_2d
  end type dim_index_2d
  
  type, public :: dim_index_3d       ! 3-D dimension index
    integer :: beg1, end1            ! lower & upper bounds of 1st dimension
    integer :: beg2, end2            ! lower & upper bounds of 2nd dimension
    integer :: beg3, end3            ! lower & upper bounds of 3rd dimension
   contains
     procedure :: dim_sizes_3d  => dim_index_3d_dim_sizes_3d
     procedure :: dim_sizes_arr => dim_index_3d_dim_size_arr
     generic   :: dim_sizes     => dim_sizes_arr, dim_sizes_3d
  end type dim_index_3d

  !---------------------------------------------------------------------------
  !
  ! field_info: A derived type containing information in an addfld call.
  !
  !---------------------------------------------------------------------------
  type, public :: field_info

    logical :: flag_xyfill                   ! non-applicable xy points flagged with fillvalue
    logical :: is_subcol                     ! .true. iff field output as subcol
    integer, pointer :: mdims(:) => NULL()   ! indicies into hist_coords list
    integer, pointer :: shape(:) => NULL()   ! shape of field on file

    real(r8) :: fillvalue                    ! fillvalue for this variable, set to default if not explicit in addfld

    integer :: numlev                        ! vertical dimension (.nc file and internal arr)

    integer :: begdim1                       ! on-node dim1 start index
    integer :: enddim1                       ! on-node dim1 end index

    integer :: begdim2                       ! on-node dim2 start index
    integer :: enddim2                       ! on-node dim2 end index

    integer :: begdim3                       ! on-node chunk or lat start index
    integer :: enddim3                       ! on-node chunk or lat end index

    logical :: colperchunk                   ! .true. iff ncols /= chunksize

    integer :: decomp_type                   ! type of decomposition (e.g., physics or dynamics)

    ! meridional and zonal complements are for fields defined as part of a
    ! 2-D vector. These IDs are used to facilitate interpolated history output
    ! At most one of these will be a positive field ID.
    integer :: meridional_complement         ! meridional field id or -1
    integer :: zonal_complement              ! zonal field id or -1

    character(len=max_fieldname_len) :: name ! field name
    character(len=max_chars) :: long_name    ! long name
    character(len=max_chars) :: units        ! units
    character(len=max_chars) :: sampling_seq ! sampling sequence - if not every timestep, how often field is sampled
    ! (i.e., how often "outfld" is called):  every other; only during LW/SW
    ! radiation calcs; etc.
  contains
    procedure :: get_shape   => field_info_get_shape
    procedure :: get_bounds  => field_info_get_bounds
    procedure :: get_dims_2d => field_info_get_dims_2d
    procedure :: get_dims_3d => field_info_get_dims_3d
    generic   :: get_dims    => get_dims_2d, get_dims_3d
  end type field_info

  real(r8), parameter, public :: fillvalue = 1.e36_r8     ! fill value for netcdf fields


  !---------------------------------------------------------------------------
  !
  ! hentry: elements of an entry in the list of active fields on a single
  !         history file
  !         nacs is an accumulation counter which normally counts an entire
  !         chunk (physics) or block (dynamics) as accumulated as a single
  !         entity. The per-chunk counting avoids counting multiple outfld
  !         calls as multiple accumulations. Only the value of the first chunk
  !         or block is written to or read from a history restart file.
  !         For certain actions (e.g., only accumulating on
  !         non-fillvalue or accumulating based on local time), nacs has an
  !         entry for every column.
  !         nacs does not keep track of levels
  !
  !---------------------------------------------------------------------------
  type, public:: hentry
    type (field_info)         :: field       ! field information
    character(len=1)          :: avgflag     ! averaging flag
    character(len=max_chars)  :: time_op     ! time operator (e.g. max, min, avg)

    integer                   :: hwrt_prec   ! history output precision
    real(r8),         pointer :: hbuf(:,:,:) => NULL()
    type(var_desc_t), pointer :: varid(:)    => NULL() ! variable ids
    integer,          pointer :: nacs(:,:)   => NULL() ! accumulation counter
    type(var_desc_t), pointer :: nacs_varid  => NULL()
  end type hentry

  !---------------------------------------------------------------------------
  !
  !  active_entry: derived type containing information for a history tape
  !
  !---------------------------------------------------------------------------
  type, public:: active_entry

    type(hentry), pointer :: hlist(:)

    integer,               pointer :: grid_ids(:) => NULL()
    type(history_patch_t), pointer :: patches(:)  => NULL()

    !
    ! PIO ids
    !

    type(file_desc_t) :: File            ! PIO file id

    type(var_desc_t) :: mdtid            ! var id for timestep
    type(var_desc_t) :: ndbaseid         ! var id for base day
    type(var_desc_t) :: nsbaseid         ! var id for base seconds of base day
    type(var_desc_t) :: nbdateid         ! var id for base date
    type(var_desc_t) :: nbsecid          ! var id for base seconds of base date
    type(var_desc_t) :: ndcurid          ! var id for current day
    type(var_desc_t) :: nscurid          ! var id for current seconds of current day
    type(var_desc_t) :: dateid           ! var id for current date
    type(var_desc_t) :: co2vmrid         ! var id for co2 volume mixing ratio
    type(var_desc_t) :: ch4vmrid         ! var id for ch4 volume mixing ratio
    type(var_desc_t) :: n2ovmrid         ! var id for n2o volume mixing ratio
    type(var_desc_t) :: f11vmrid         ! var id for f11 volume mixing ratio
    type(var_desc_t) :: f12vmrid         ! var id for f12 volume mixing ratio
    type(var_desc_t) :: sol_tsiid        ! var id for total solar irradiance (W/m2)
    type(var_desc_t) :: datesecid        ! var id for curent seconds of current date
#if ( defined BFB_CAM_SCAM_IOP )
    type(var_desc_t) :: bdateid         ! var id for base date
    type(var_desc_t) :: tsecid        ! var id for curent seconds of current date
#endif
    type(var_desc_t) :: nstephid         ! var id for current timestep
    type(var_desc_t) :: timeid           ! var id for time
    type(var_desc_t) :: tbndid           ! var id for time_bnds
    type(var_desc_t) :: date_writtenid   ! var id for date time sample written
    type(var_desc_t) :: time_writtenid   ! var id for time time sample written
    type(var_desc_t) :: f107id           ! var id for f107
    type(var_desc_t) :: f107aid          ! var id for f107a
    type(var_desc_t) :: kpid             ! var id for kp
    type(var_desc_t) :: apid             ! var id for ap

  end type active_entry

  !---------------------------------------------------------------------------
  !
  !  formula_terms_t: Information for formula terms (CF convention) variables
  !                   Used to add a formula-terms variable to the history file
  !                   Also adds a string, '<name>: <var_name>' to the parent
  !                   mdim's 'formula_terms' attribute.
  !
  !---------------------------------------------------------------------------
  type, public :: formula_terms_t
    character(len=max_fieldname_len) :: a_name = ''   ! 'A' term variable name
    character(len=max_string_len)    :: a_long_name = '' ! 'A' long name
    real(r8), pointer                :: a_values(:) => null() ! 'A' variable values
    character(len=max_fieldname_len) :: b_name = ''   ! 'B' term variable name
    character(len=max_string_len)    :: b_long_name = '' ! 'B' long name
    real(r8), pointer                :: b_values(:) => null() ! 'B' variable values
    character(len=max_fieldname_len) :: p0_name = ''  ! 'p0' term variable name
    character(len=max_string_len)    :: p0_long_name = '' ! 'p0' long name
    character(len=max_chars)         :: p0_units = '' ! 'p0' variable units
    real(r8)                         :: p0_value = fillvalue ! 'p0' variable values
    character(len=max_fieldname_len) :: ps_name = ''  ! 'ps' term variable name
  end type formula_terms_t

  !---------------------------------------------------------------------------
  !
  !  hist_coord_t: Information for history variable dimension attributes
  !
  !---------------------------------------------------------------------------
  type, public :: hist_coord_t
    character(len=max_hcoordname_len) :: name = ''  ! coordinate name
    integer                  :: dimsize = 0       ! size of dimension
    character(len=max_chars) :: long_name = ''    ! 'long_name' attribute
    character(len=max_chars) :: units = ''        ! 'units' attribute
    character(len=max_chars) :: bounds_name = ''  ! 'bounds' attribute (& name of bounds variable)
    character(len=max_chars) :: standard_name = ''! 'standard_name' attribute
    character(len=4)         :: positive = ''     ! 'positive' attribute ('up' or 'down')
    integer,  pointer        :: integer_values(:) => null() ! dim values if integral
    real(r8), pointer        :: real_values(:) => null() ! dim values if real
    real(r8), pointer        :: bounds(:,:) => null() ! dim bounds
    type(formula_terms_t)    :: formula_terms     ! vars for formula terms
    logical                  :: integer_dim       ! .true. iff dim has integral values
    logical                  :: vertical_coord    ! .true. iff dim is vertical
  end type hist_coord_t

  ! Some parameters for use with interpolated output namelist items
  integer,          parameter, public :: interp_type_native            = 0
  integer,          parameter, public :: interp_type_bilinear          = 1
  integer,          parameter, public :: interp_gridtype_equal_poles   = 1
  integer,          parameter, public :: interp_gridtype_gauss         = 2
  integer,          parameter, public :: interp_gridtype_equal_nopoles = 3

  !---------------------------------------------------------------------------
  !
  !  interp_info_t: Information for lat/lon interpolated history output
  !
  !---------------------------------------------------------------------------
  type, public :: interp_info_t
    ! store the  lat-lon grid information
    character(len=28)     :: gridname = ''
    integer               :: grid_id  = -1
    ! gridtype = 1      equally spaced, including poles (FV scalars output grid)
    ! gridtype = 2      Gauss grid (CAM Eulerian)
    ! gridtype = 3      equally spaced, no poles (FV staggered velocity)
    integer               :: interp_gridtype = interp_gridtype_equal_poles
    ! interpolate_type = 0: native high order interpolation
    ! interpolate_type = 1: bilinear interpolation
    integer               :: interp_type = interp_type_bilinear
    integer               :: interp_nlat = 0
    integer               :: interp_nlon = 0
    real(r8), pointer     :: interp_lat(:) => NULL()
    real(r8), pointer     :: interp_lon(:) => NULL()
    real(r8), pointer     :: interp_gweight(:) => NULL()
  end type interp_info_t

  !! Coordinate variables
  integer,                     public :: registeredmdims = 0
  integer,                     public :: maxvarmdims     = 1
  character(len=9), parameter, public :: mdim_var_name   = 'mdimnames'
  integer,          parameter         :: maxmdims        = 25  ! arbitrary limit
  type(hist_coord_t),          public :: hist_coords(maxmdims)

  public     :: add_hist_coord, add_vert_coord
  public     :: write_hist_coord_attrs, write_hist_coord_vars
  public     :: lookup_hist_coord_indices, hist_coord_find_levels
  public     :: get_hist_coord_index, hist_coord_name, hist_coord_size
  public     :: sec2hms, date2yyyymmdd

  interface add_hist_coord
    module procedure add_hist_coord_regonly
    module procedure add_hist_coord_int
    module procedure add_hist_coord_r8
  end interface

  interface hist_coord_size
    module procedure hist_coord_size_char
    module procedure hist_coord_size_int
  end interface

  interface assignment(=)
    module procedure field_copy
    module procedure formula_terms_copy
  end interface

  interface check_hist_coord
    ! NB: This is supposed to be a private interface
    ! check_hist_coord: returns 0 if <name> is not registered as an mdim
    !                   returns i if <name> is registered with compatible values
    !              calls endrun if <name> is registered with incompatible values
    ! Versions without the <name> argument return .true. or .false.
    module procedure check_hist_coord_char
    module procedure check_hist_coord_int
    module procedure check_hist_coord_int1
    module procedure check_hist_coord_r8
    module procedure check_hist_coord_r81
    module procedure check_hist_coord_r82
    module procedure check_hist_coord_ft
    module procedure check_hist_coord_all
  end interface

!!---------------------------------------------------------------------------

contains

  subroutine dim_index_2d_dim_sizes_2d(this, dim1, dim2)

    ! Dummy arguments
    class(dim_index_2d)                :: this
    integer,             intent(out)   :: dim1
    integer,             intent(out)   :: dim2

    dim1 = MAX(0, this%end1 - this%beg1 + 1)
    dim2 = MAX(0, this%end2 - this%beg2 + 1)
    
  end subroutine dim_index_2d_dim_sizes_2d

  subroutine dim_index_2d_dim_size_arr(this, dims)

    ! Dummy arguments
    class(dim_index_2d)                :: this
    integer,             intent(out)   :: dims(:)

    if (size(dims) < 2) then
      call endrun('dim_index_2d_dim_sizes: dims must have at least two elements')
    end if

    call this%dim_sizes(dims(1), dims(2))
    
  end subroutine dim_index_2d_dim_size_arr

  subroutine dim_index_3d_dim_sizes_3d(this, dim1, dim2, dim3)

    ! Dummy arguments
    class(dim_index_3d)                :: this
    integer,             intent(out)   :: dim1
    integer,             intent(out)   :: dim2
    integer,             intent(out)   :: dim3

    dim1 = MAX(0, this%end1 - this%beg1 + 1)
    dim2 = MAX(0, this%end2 - this%beg2 + 1)
    dim3 = MAX(0, this%end3 - this%beg3 + 1)
    
  end subroutine dim_index_3d_dim_sizes_3d

  subroutine dim_index_3d_dim_size_arr(this, dims)

    ! Dummy arguments
    class(dim_index_3d)                :: this
    integer,             intent(out)   :: dims(:)

    if (size(dims) < 3) then
      call endrun('dim_index_3d_dim_sizes: dims must have at least three elements')
    end if

    call this%dim_sizes(dims(1), dims(2), dims(3))
    
  end subroutine dim_index_3d_dim_size_arr

  ! field_info_get_dims_2d: Retrieve bounds for stepping through a chunk
  type(dim_index_2d) function field_info_get_dims_2d(this, col) result(dims)
    use cam_grid_support, only: cam_grid_get_block_count

    ! Dummy argument
    class(field_info)                :: this
    integer,           intent(in)    :: col

    ! Local variable
    integer              :: endi

    if (this%colperchunk) then
      endi = this%begdim1 + cam_grid_get_block_count(this%decomp_type, col) - 1
      dims = dim_index_2d(this%begdim1, endi, this%begdim2, this%enddim2)
    else
      dims = dim_index_2d(this%begdim1, this%enddim1, this%begdim2, this%enddim2)
    end if
  end function field_info_get_dims_2d

  ! field_info_get_dims_3d: Retrieve grid decomp bounds
  type(dim_index_3d) function field_info_get_dims_3d(this) result(dims)

    ! Dummy argument
    class(field_info)                :: this

    dims = dim_index_3d(this%begdim1, this%enddim1, this%begdim2, this%enddim2,&
         this%begdim3, this%enddim3)

  end function field_info_get_dims_3d

  ! field_info_get_shape: Return a pointer to the field's global shape.
  !                       Calculate it first if necessary
  subroutine field_info_get_shape(this, shape_out, rank_out)
    use cam_grid_support, only: cam_grid_dimensions

    ! Dummy arguments
    class(field_info)                         :: this
    integer,                    intent(out)   :: shape_out(:)
    integer,                    intent(out)   :: rank_out

    ! Local arguments
    integer                                   :: rank, i, pos
    integer                                   :: gdims(2)

    if (.not. associated(this%shape)) then
      ! Calculate field's global shape
      call cam_grid_dimensions(this%decomp_type, gdims, rank)
      pos = rank
      if (associated(this%mdims)) then
        rank = rank + size(this%mdims)
      end if
      allocate(this%shape(rank))
      this%shape(1:pos) = gdims(1:pos)
      if (rank > pos) then
        do i = 1, size(this%mdims)
          pos = pos + 1
          this%shape(pos) = hist_coords(this%mdims(i))%dimsize
        end do
      end if
    end if

    rank_out = size(this%shape)

    if (size(shape_out) < rank_out) then
      call endrun('field_info_get_shape: shape_out too small')
    end if

    shape_out(1:rank_out) = this%shape(1:rank_out)
    if (size(shape_out) > rank_out) then
      shape_out(rank_out+1:) = 1
    end if

  end subroutine field_info_get_shape

  subroutine field_info_get_bounds(this, dim, beg, end)

    ! Dummy arguments
    class(field_info)                :: this
    integer,           intent(in)    :: dim
    integer,           intent(out)   :: beg
    integer,           intent(out)   :: end

    select case(dim)
    case (1)
      beg = this%begdim1
      end = this%enddim1
    case (2)
      beg = this%begdim2
      end = this%enddim2
    case (3)
      beg = this%begdim3
      end = this%enddim3
    case default
      call endrun('field_info_get_bounds: dim must be 1, 2, or 3')
    end select

  end subroutine field_info_get_bounds

  ! history_patch_write_attrs: Define coordinate variables and attributes
  !               for a patch
  subroutine history_patch_write_attrs(this, File)
    use cam_grid_support, only: cam_grid_is_unstructured
    use pio,              only: file_desc_t, var_desc_t, pio_put_att, pio_double
    use cam_pio_utils,    only: cam_pio_def_dim, cam_pio_def_var, cam_pio_handle_error

    ! Dummy arguments
    class(history_patch_t)                  :: this
    type(file_desc_t),        intent(inout) :: File    ! PIO file Handle

    ! Local variable
    type(cam_grid_patch_t), pointer         :: patchptr
    type(var_desc_t), pointer               :: vardesc_lat => NULL()
    type(var_desc_t), pointer               :: vardesc_lon => NULL()
    character(len=128)                      :: errormsg
    character(len=max_chars)                :: lat_name
    character(len=max_chars)                :: lon_name
    character(len=max_chars)                :: col_name
    character(len=max_chars)                :: temp_str
    integer                                 :: dimid1, dimid2 ! PIO dim IDs
    integer                                 :: num_patches
    integer                                 :: temp1, temp2
    integer                                 :: latid, lonid ! Coordinate dims
    integer                                 :: i, ierr
    logical                                 :: col_only
    logical                                 :: unstruct

    num_patches = size(this%patches)
    if (associated(this%header_info)) then
      ! Make sure header_info is the right size
      if (size(this%header_info) /= num_patches) then
        write(errormsg, '(a,2(i0,a))') 'Size mismatch between header_info (', &
             size(this%header_info), ') and patches (', num_patches, ')'
        call endrun('history_patch_write_attrs: '//errormsg)
      end if
    else
      allocate(this%header_info(num_patches))
    end if

    ! Write attributes for each patch
    do i = 1, num_patches
      patchptr => this%patches(i)
      call this%header_info(i)%set_gridid(patchptr%gridid())
      unstruct = cam_grid_is_unstructured(patchptr%gridid())
      ! What are the dimension(s) for this patch?
      col_only = this%collected_output
      if (num_patches == 1) then
        ! Backwards compatibility
        if (unstruct) then
          col_name = 'ncol'
        else
          col_name = ''
        end if
        lat_name = 'lat'
        lon_name = 'lon'
      else
        call patchptr%get_axis_names(lat_name, lon_name, col_name, col_only)
      end if
      ! Define the dimensions (latx/lonx or ncolx)
      ! col_name is set for unstructured output (ncolx)
      if (len_trim(col_name) > 0) then
        call patchptr%get_global_size(gsize=temp1)
        if (temp1 <= 0) then
          call endrun('history_patch_write_attrs: col dimsize must be positive')
        end if
        if (unstruct .and. (.not. col_only)) then
          ! For the case of unstructured output without collected column
          ! output, we need to make sure that the ncolx dimension is unique
          col_name = trim(col_name)//'_'//trim(this%lon_axis_name)//'_'//trim(this%lat_axis_name)
        end if
        call cam_pio_def_dim(File, trim(col_name), temp1, dimid1, existOK=.false.)
        call this%header_info(i)%set_hdims(dimid1)
        latid = dimid1
        lonid = dimid1
      else
        lat_name = trim(lat_name)//'_'//trim(this%lat_axis_name)
        call patchptr%get_global_size(temp1, temp2)
        if (temp1 <= 0) then
          call endrun('history_patch_write_attrs: lat dimsize must be positive')
        end if
        call cam_pio_def_dim(File, trim(lat_name), temp1, dimid1, existOK=.true.)
        latid = dimid1
        lon_name = trim(lon_name)//'_'//trim(this%lon_axis_name)
        if (temp2 <= 0) then
          call endrun('history_patch_write_attrs: lon dimsize must be positive')
        end if
        call cam_pio_def_dim(File, trim(lon_name), temp2, dimid2, existOK=.true.)
        lonid = dimid2
        call this%header_info(i)%set_hdims(lonid, latid)
      end if
      !! Define the latx (coordinate) variable
      if (unstruct .and. (.not. col_only)) then
        ! We need to make sure the latx name is unique
        lat_name = trim(lat_name)//'_'//trim(this%lon_axis_name)//'_'//trim(this%lat_axis_name)
      end if
      allocate(vardesc_lat)
      call cam_pio_def_var(File, trim(lat_name), pio_double, (/latid/),       &
           vardesc_lat, existOK=.true.)
      ! Coordinate attributes
      call patchptr%get_coord_long_name('lat', temp_str)
      if (len_trim(temp_str) > 0) then
        ierr = pio_put_att(File, vardesc_lat, 'long_name', trim(temp_str))
        call cam_pio_handle_error(ierr, 'history_patch_write_attrs: Unable to define long_name')
      end if
      call patchptr%get_coord_units('lat', temp_str)
      if (len_trim(temp_str) > 0) then
        ierr = pio_put_att(File, vardesc_lat, 'units', trim(temp_str))
        call cam_pio_handle_error(ierr, 'history_patch_write_attrs: Unable to define units')
      end if
      !! Define the lonx (coordinate) variable
      if (unstruct .and. (.not. col_only)) then
        ! We need to make sure the lonx name is unique
        lon_name = trim(lon_name)//'_'//trim(this%lon_axis_name)//'_'//trim(this%lat_axis_name)
      end if
      allocate(vardesc_lon)
      call cam_pio_def_var(File, trim(lon_name), pio_double, (/lonid/),       &
           vardesc_lon, existOK=.true.)
      ! Coordinate attributes
      call patchptr%get_coord_long_name('lon', temp_str)
      if (len_trim(temp_str) > 0) then
        ierr = pio_put_att(File, vardesc_lon, 'long_name', trim(temp_str))
        call cam_pio_handle_error(ierr, 'history_patch_write_attrs: Unable to define long_name')
      end if
      call patchptr%get_coord_units('lon', temp_str)
      if (len_trim(temp_str) > 0) then
        ierr = pio_put_att(File, vardesc_lon, 'units', trim(temp_str))
        call cam_pio_handle_error(ierr, 'history_patch_write_attrs: Unable to define units')
      end if
      call this%header_info(i)%set_varids(vardesc_lon, vardesc_lat)
      nullify(vardesc_lat, vardesc_lon) ! They belong to the header_info now
    end do

  end subroutine history_patch_write_attrs

  ! history_patch_write_vals: Write coordinate variable values for a patch
  subroutine history_patch_write_vals(this, File)
    use pio,           only: file_desc_t, var_desc_t

    ! Dummy arguments
    class(history_patch_t)                  :: this
    type(file_desc_t),        intent(inout) :: File    ! PIO file Handle

    ! Local variable
    type(cam_grid_patch_t), pointer         :: patchptr
    type(var_desc_t), pointer               :: vardesc => NULL()  ! PIO var desc
    character(len=128)                      :: errormsg
    character(len=max_chars)                :: lat_name
    character(len=max_chars)                :: lon_name
    character(len=max_chars)                :: col_name
    character(len=max_chars)                :: temp_str
    integer                                 :: dimid    ! PIO dimension ID
    integer                                 :: num_patches
    integer                                 :: temp1, temp2
    integer                                 :: latid, lonid ! Coordinate dims
    integer                                 :: i
    logical                                 :: col_only

    num_patches = size(this%patches)
    if (.not. associated(this%header_info)) then
      ! We need this for dim and coord variable descriptors
      write(errormsg, '(2a)') 'No header info for ', trim(this%namelist_entry)
      call endrun('history_patch_write_vals: '//errormsg)
    end if

    ! Write attributes for each patch
    do i = 1, num_patches
      patchptr => this%patches(i)
      ! Write the coordinate variables (or just lat/lon for column output)
      call patchptr%write_coord_vals(File, this%header_info(i))
    end do

  end subroutine history_patch_write_vals

  ! history_patch_field_name: Add patch description to field name
  subroutine history_patch_field_name(this, name)
    ! Dummy arguments
    class(history_patch_t)                     :: this
    character(len=*),            intent(inout) :: name

    ! Add patch description info to the variable name
    name = trim(name)//'_'//trim(this%lon_axis_name)//'_'//trim(this%lat_axis_name)
  end subroutine history_patch_field_name

  ! history_patch_num_hdims: Find the number of horizontal dimensions for
  !         the indicated grid
  integer function history_patch_num_hdims(this, gridid)
    ! Dummy arguments
    class(history_patch_t)                     :: this
    integer,                     intent(in)    :: gridid      ! The field's grid

    ! Local variables
    type(cam_grid_patch_t),        pointer     :: patchptr
    character(len=128)                         :: errormsg
    integer                                    :: i
    integer                                    :: num_patches
 
    ! Basic sanity checks, is this patch OK?
    num_patches = size(this%patches)
    if (associated(this%header_info)) then
      ! Make sure header_info is the right size
      if (size(this%header_info) /= num_patches) then
        write(errormsg, '(a,2(i0,a))') 'Size mismatch between header_info (', &
             size(this%header_info), ') and patches (', num_patches, ')'
        call endrun('history_patch_num_hdims: '//errormsg)
      end if
    else
      write(errormsg, *) 'No header info for patch, ', trim(this%namelist_entry)
      call endrun('history_patch_num_hdims: '//errormsg)
    end if

    ! Find the correct patch by matching grid ID
    history_patch_num_hdims = -1
    do i = 1, num_patches
      patchptr => this%patches(i)
      if (patchptr%gridid() == gridid) then
        ! This is the right patch, set the return val and quit loop
        history_patch_num_hdims = this%header_info(i)%num_hdims()
        exit
      else if (i >= num_patches) then
        write(errormsg, '(3a,i0)') 'No grid found for patch, ',               &
             trim(this%namelist_entry), '. Was looking for decomp ', gridid
        call endrun('history_patch_num_hdims: '//errormsg)
      ! No else needed
      end if
    end do
    if (history_patch_num_hdims <= 0) then
      write(errormsg, '(2a,2(a,i0))') 'INTERNAL: No grid patch for ',         &
           trim(this%namelist_entry), ', num_patches = ',num_patches,         &
           ', gridid = ', gridid
      call endrun('history_patch_num_hdims: '//errormsg)
    end if

  end function history_patch_num_hdims

  ! history_patch_get_var_data: Calculate data relevant to history variable
  !         on a patch by substituting patch dimension ids for the horiz. ids
  !         and adding patch information to the variable name
  subroutine history_patch_get_var_data(this, name, dimids, gridid)
    ! Dummy arguments
    class(history_patch_t)                     :: this
    character(len=*),            intent(inout) :: name
    integer,                     intent(inout) :: dimids(:)   ! Grid dimids
    integer,                     intent(in)    :: gridid      ! The field's grid

    ! Local variables
    type(cam_grid_patch_t),        pointer     :: patchptr
    type (cam_grid_header_info_t), pointer     :: histptr
    character(len=128)                         :: errormsg
    integer                                    :: num_patches
    integer                                    :: i

    ! Basic sanity checks, is this patch OK?
    num_patches = size(this%patches)
    if (associated(this%header_info)) then
      ! Make sure header_info is the right size
      if (size(this%header_info) /= num_patches) then
        write(errormsg, '(a,2(i0,a))') 'Size mismatch between header_info (', &
             size(this%header_info), ') and patches (', num_patches, ')'
        call endrun('history_patch_get_var_data: '//errormsg)
      end if
    else
      write(errormsg, *) 'No header info for patch, ', trim(this%namelist_entry)
      call endrun('history_patch_get_var_data: '//errormsg)
    end if

    ! Find the correct patch by matching grid ID
    do i = 1, num_patches
      patchptr => this%patches(i)
      if (patchptr%gridid() == gridid) then
        ! This is the right patch, quit loop
        histptr  => this%header_info(i)
        exit
      else if (i >= num_patches) then
        write(errormsg, '(3a,i0)') 'No grid found for patch, ',               &
             trim(this%namelist_entry), '. Was looking for decomp ', gridid
        call endrun('history_patch_get_var_data: '//errormsg)
      ! No else needed
      end if
    end do

    ! We have the correct patch, replace the horizontal dimension ids
    do i = 1, histptr%num_hdims()
      dimids(i) = histptr%get_hdimid(i)
    end do
    ! Re-define the variable name
    call this%field_name(name)

  end subroutine history_patch_get_var_data

  subroutine history_patch_compact(this)

    ! Dummy arguments
    class(history_patch_t)                   :: this

    ! Local variables
    integer                                  :: num_patches
    integer                                  :: i
 
    num_patches = size(this%patches)

    ! Find the correct patch by matching grid ID
    do i = 1, num_patches
      call this%patches(i)%compact(this%collected_output)
    end do

  end subroutine history_patch_compact

  subroutine history_patch_write_var(this, File, gridid, adims, dtype, hbuf, varid)
    use pio,           only: file_desc_t, var_desc_t, io_desc_t
    use pio,           only: pio_write_darray
    use cam_pio_utils, only: cam_pio_handle_error, cam_pio_var_info

    ! Dummy arguments
    class(history_patch_t)                   :: this
    type(file_desc_t),         intent(inout) :: File       ! PIO file handle
    integer,                   intent(in)    :: gridid
    integer,                   intent(in)    :: adims(:)
    integer,                   intent(in)    :: dtype
    real(r8),                  intent(in)    :: hbuf(:,:,:)
    type(var_desc_t),       pointer          :: varid

    ! Local variables
    type(cam_grid_patch_t), pointer          :: patchptr
    character(len=128)                       :: errormsg
    integer                                  :: num_patches
    integer                                  :: i
    type(io_desc_t),        pointer          :: iodesc
    integer                                  :: ierr, nfdims
    integer                                  :: fdimlens(7), dimids(7)
 
    num_patches = size(this%patches)

    ! Find the correct patch by matching grid ID
    do i = 1, num_patches
      patchptr => this%patches(i)
      if (patchptr%gridid() == gridid) then
        ! This is the right patch, quit loop
        exit
      else if (i >= num_patches) then
        write(errormsg, '(3a,i0)') 'No grid found for patch, ',               &
             trim(this%namelist_entry), '. Was looking for decomp ', gridid
        call endrun('history_patch_write_var: '//trim(errormsg))
      ! No else needed
      end if
    end do

    ! We have the right grid, write the hbuf
    call cam_pio_var_info(File, varid, nfdims, dimids, fdimlens)
    call patchptr%get_decomp(adims, fdimlens(1:nfdims), dtype, iodesc)
    if (size(adims) == 2) then
      call pio_write_darray(File, varid, iodesc, hbuf(:,1,:), ierr)
    else if (size(adims) == 3) then
      call pio_write_darray(File, varid, iodesc, hbuf, ierr)
    else
      call endrun("history_patch_write_var: adims must be rank 2 or 3")
    end if
    call cam_pio_handle_error(ierr, 'history_patch_write_var: Error writing variable')

  end subroutine history_patch_write_var

  subroutine history_patch_active_cols(this, gridid, lchnk, active)
    ! Dummy arguments
    class(history_patch_t)                :: this
    integer,                intent(in)    :: gridid ! desired grid
    integer,                intent(in)    :: lchnk  ! chunk or block number
    logical,                intent(out)   :: active(:)

    ! Local variables
    type(cam_grid_patch_t),        pointer     :: patchptr
    character(len=128)                         :: errormsg
    integer                                    :: num_patches
    integer                                    :: i
 
    num_patches = size(this%patches)

    ! Find the correct patch by matching grid ID
    do i = 1, num_patches
      patchptr => this%patches(i)
      if (patchptr%gridid() == gridid) then
        ! This is the right patch, quit loop
        exit
      else if (i >= num_patches) then
        write(errormsg, '(3a,i0)') 'No grid found for patch, ',               &
             trim(this%namelist_entry), '. Was looking for decomp ', gridid
        call endrun('history_patch_active_cols: '//errormsg)
      ! No else needed
      end if
    end do

    ! If we get here, patchptr is the grid patch we want
    call patchptr%active_cols(lchnk, active)

 end subroutine history_patch_active_cols

  subroutine history_patch_deallocate(this)
    ! Dummy argument
    class(history_patch_t) :: this
    ! Local variable
    integer                :: i

    this%lon_axis_name = ''
    this%lat_axis_name = ''

    if (associated(this%patches)) then
      do i = 1, size(this%patches)
        call this%patches(i)%deallocate()
      end do
      deallocate(this%patches)
      nullify(this%patches)
    end if

    if (associated(this%header_info)) then
      do i = 1, size(this%header_info)
        call this%header_info(i)%deallocate()
      end do
      deallocate(this%header_info)
      nullify(this%header_info)
    end if

  end subroutine history_patch_deallocate

  subroutine field_copy(f_out, f_in)
    type(field_info), intent(in) :: f_in
    type(field_info), intent(out) :: f_out

    f_out%flag_xyfill= f_in%flag_xyfill
    f_out%is_subcol = f_in%is_subcol
    f_out%fillvalue= f_in%fillvalue
    f_out%numlev =  f_in%numlev                      ! vertical dimension (.nc file and internal arr)
    f_out%begdim1 = f_in%begdim1                     ! on-node dim1 start index
    f_out%enddim1 = f_in%enddim1                     ! on-node dim1 end index
    f_out%begdim2 = f_in%begdim2                     ! on-node dim2 start index
    f_out%enddim2 = f_in%enddim2                     ! on-node dim2 end index
    f_out%begdim3 = f_in%begdim3                     ! on-node chunk or lat start index
    f_out%enddim3 = f_in%enddim3                     ! on-node chunk or lat end index
    f_out%decomp_type = f_in%decomp_type             ! type of decomposition (physics or dynamics)

    f_out%meridional_complement = f_in%meridional_complement ! id  or -1
    f_out%zonal_complement = f_in%zonal_complement           ! id  or -1

    f_out%name = f_in%name                           ! field name
    f_out%long_name = f_in%long_name                 ! long name
    f_out%units = f_in%units                         ! units
    f_out%sampling_seq =  f_in%sampling_seq          ! sampling sequence - if not every timestep, how often field is sampled

    if(associated(f_in%mdims)) then
      f_out%mdims=>f_in%mdims
    else
      nullify(f_out%mdims)
    end if

  end subroutine field_copy

  subroutine formula_terms_copy(f_out, f_in)
    type(formula_terms_t), intent(in) :: f_in
    type(formula_terms_t), intent(out) :: f_out

    f_out%a_name = f_in%a_name
    f_out%a_long_name = f_in%a_long_name
    f_out%a_values => f_in%a_values
    f_out%b_name = f_in%b_name
    f_out%b_long_name = f_in%b_long_name
    f_out%b_values => f_in%b_values
    f_out%p0_name = f_in%p0_name
    f_out%p0_long_name = f_in%p0_long_name
    f_out%p0_units = f_in%p0_units
    f_out%p0_value = f_in%p0_value
    f_out%ps_name = f_in%ps_name
  end subroutine formula_terms_copy

  integer function get_hist_coord_index(mdimname)
    ! Input variables
    character(len=*), intent(in)            :: mdimname
    ! Local variable
    integer :: i
    
    get_hist_coord_index = -1
    do i = 1, registeredmdims
      if(trim(mdimname) == trim(hist_coords(i)%name)) then
        get_hist_coord_index = i
        exit
      end if
    end do

  end function get_hist_coord_index

  character(len=max_hcoordname_len) function hist_coord_name(index)
    ! Input variables
    integer, intent(in)            :: index
    
    if ((index > 0) .and. (index <= registeredmdims)) then
      hist_coord_name = hist_coords(index)%name
    else
      call endrun('hist_coord_name: index out of range')
    end if

  end function hist_coord_name

  integer function hist_coord_size_int(index)
    ! Input variables
    integer, intent(in)            :: index
    
    if (index > 0) then
      hist_coord_size_int = hist_coords(index)%dimsize
    else
      hist_coord_size_int = -1
    end if

  end function hist_coord_size_int

  integer function hist_coord_size_char(mdimname)
    ! Input variables
    character(len=*), intent(in)            :: mdimname
    ! Local variable
    integer :: i
    
    i = get_hist_coord_index(mdimname)
    hist_coord_size_char = hist_coord_size(i)

  end function hist_coord_size_char

  ! Functions to check consistent term definition for hist coords
  logical function check_hist_coord_char(defined, input)

    ! Input variables
    character(len=*), intent(in)            :: defined
    character(len=*), intent(in), optional  :: input

    if (len_trim(defined) == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_char = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_char = (trim(input) == trim(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_char = .false.
    end if
  end function check_hist_coord_char

  logical function check_hist_coord_int(defined, input)

    ! Input variables
    integer, intent(in)            :: defined
    integer, intent(in), optional  :: input

    if (defined == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_int = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_int = (input == defined)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_int = .false.
    end if
  end function check_hist_coord_int

  logical function check_hist_coord_int1(defined, input)

    ! Input variables
    integer,             pointer            :: defined(:)
    integer, intent(in),          optional  :: input(:)

    ! Local variables
    integer                                 :: i

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_int1 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_int1 = (size(input) == size(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_int1 = .false.
    end if
    if (check_hist_coord_int1 .and. associated(defined)) then
      ! Need to check the values
      do i = 1, size(defined)
        if (defined(i) /= input(i)) then
          check_hist_coord_int1 = .false.
          exit
        end if
      end do
    end if
  end function check_hist_coord_int1

  logical function check_hist_coord_r8(defined, input)

    ! Input variables
    real(r8), intent(in)            :: defined
    real(r8), intent(in), optional  :: input

    if (defined == fillvalue) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r8 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r8 = (input == defined)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r8 = .false.
    end if
  end function check_hist_coord_r8

  logical function check_hist_coord_r81(defined, input)

    ! Input variables
    real(r8),             pointer            :: defined(:)
    real(r8), intent(in),          optional  :: input(:)

    ! Local variables
    integer                                  :: i

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r81 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r81 = (size(input) == size(defined))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r81 = .false.
    end if
    if (check_hist_coord_r81 .and. associated(defined)) then
      ! Need to check the values
      do i = 1, size(defined)
        if (defined(i) /= input(i)) then
          check_hist_coord_r81 = .false.
          exit
        end if
      end do
    end if
  end function check_hist_coord_r81

  logical function check_hist_coord_r82(defined, input)

    ! Input variables
    real(r8),             pointer            :: defined(:,:)
    real(r8), intent(in),          optional  :: input(:,:)

    ! Local variables
    integer                                  :: i, j

    if (.not. associated(defined)) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_r82 = .true.
    else if (present(input)) then
      ! We have to match definitions
      check_hist_coord_r82 = ((size(input, 1) == size(defined, 1)) .and.    &
                              (size(input, 2) == size(defined, 2)))
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_r82 = .false.
    end if
    if (check_hist_coord_r82 .and. associated(defined)) then
      ! Need to check the values
      do j = 1, size(defined, 2)
        do i = 1, size(defined, 1)
          if (defined(i, j) /= input(i, j)) then
            check_hist_coord_r82 = .false.
            exit
          end if
        end do
      end do
    end if
  end function check_hist_coord_r82

  logical function check_hist_coord_ft(defined, input)

    ! Input variables
    type(formula_terms_t), intent(in)           :: defined
    type(formula_terms_t), intent(in), optional :: input

    ! We will assume that if formula_terms has been defined, a_name has a value
    if (len_trim(defined%a_name) == 0) then
      ! In this case, we assume the current value is undefined so any input OK
      check_hist_coord_ft = .true.
    else if (present(input)) then
      ! We have to match definitions
      ! Need to check the values
      check_hist_coord_ft =                                                   &
           check_hist_coord(defined%a_name,       input%a_name)         .and. &
           check_hist_coord(defined%a_long_name,  input%a_long_name)    .and. &
           check_hist_coord(defined%a_values,     input%a_values)       .and. &
           check_hist_coord(defined%b_name,       input%b_name)         .and. &
           check_hist_coord(defined%b_long_name,  input%b_long_name)    .and. &
           check_hist_coord(defined%b_values,     input%b_values)       .and. &
           check_hist_coord(defined%p0_name,      input%p0_name)        .and. &
           check_hist_coord(defined%p0_long_name, input%p0_long_name)   .and. &
           check_hist_coord(defined%p0_units,     input%p0_units)       .and. &
           check_hist_coord(defined%p0_value,     input%p0_value)       .and. &
           check_hist_coord(defined%ps_name,      input%ps_name)
    else
      ! Not sure here. We have a value and are redefining without one?
      check_hist_coord_ft = .false.
    end if
  end function check_hist_coord_ft

  ! check_hist_coord: returns 0 if <name> is not registered as a hist coord
  !                   returns i if <name> is registered with compatible values
  !                   calls endrun if <name> is registered with incompatible values
  integer function check_hist_coord_all(name, vlen, long_name, units, bounds, &
       i_values, r_values, bounds_name, positive, standard_name, formula_terms)

    ! Input variables
    character(len=*),      intent(in)            :: name
    integer,               intent(in)            :: vlen
    character(len=*),      intent(in),  optional :: long_name
    character(len=*),      intent(in),  optional :: units
    character(len=*),      intent(in),  optional :: bounds_name
    integer,               intent(in),  optional :: i_values(:)
    real(r8),              intent(in),  optional :: r_values(:)
    real(r8),              intent(in),  optional :: bounds(:,:)
    character(len=*),      intent(in),  optional :: positive
    character(len=*),      intent(in),  optional :: standard_name
    type(formula_terms_t), intent(in),  optional :: formula_terms

    ! Local variables
    character(len=120)                           :: errormsg
    integer                                      :: i

    i = get_hist_coord_index(trim(name))
    ! If i > 0, this mdim has already been registered
    if (i > 0) then
      check_hist_coord_all = i
      if (.not. check_hist_coord(hist_coords(i)%dimsize, vlen)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, '//trim(name)//' with incompatible size'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%long_name, long_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different long_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%units, units)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different units'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%bounds_name, bounds_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different bounds_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%standard_name, standard_name)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different standard_name'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%positive, positive)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with a different value of positive'
        call endrun(errormsg)
      end if
      ! Since the integer_dim defaults to .true., double check which to check
      if ((.not. hist_coords(i)%integer_dim) .or.                             &
           associated(hist_coords(i)%real_values)) then
        if (.not. check_hist_coord(hist_coords(i)%real_values, r_values)) then
          write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different values'
          call endrun(errormsg)
        else if (present(i_values)) then
          write(errormsg, *) 'ERROR: Attempt to register integer values for real dimension'
          call endrun(errormsg)
        end if
      else
        if (.not. check_hist_coord(hist_coords(i)%integer_values, i_values)) then
          write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different values'
          call endrun(errormsg)
        else if (present(i_values) .and. present(r_values)) then
          write(errormsg, *) 'ERROR: Attempt to register real values for integer dimension'
          call endrun(errormsg)
        end if
      end if
      if (.not. check_hist_coord(hist_coords(i)%bounds, bounds)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different bounds'
        call endrun(errormsg)
      end if
      if (.not. check_hist_coord(hist_coords(i)%formula_terms, formula_terms)) then
        write(errormsg, *) 'ERROR: Attempt to register dimension, ',trim(name),' with different formula_terms'
        call endrun(errormsg)
      end if
    else
      check_hist_coord_all = 0
    end if
  end function check_hist_coord_all

  subroutine add_hist_coord_regonly(name, index)

    ! Input variable
    character(len=*),  intent(in)    :: name
    integer, optional, intent(out)   :: index

    ! Local variables
    character(len=120)               :: errormsg
    integer                          :: i

    if ((trim(name) == trim(horiz_only)) .or. (len_trim(name) == 0)) then
      call endrun('ADD_HIST_COORD: '//trim(name)//' is not a valid coordinate name')
    end if
    i = get_hist_coord_index(trim(name))
    ! If i > 0, this mdim has already been registered
    if (i <= 0) then
      registeredmdims = registeredmdims + 1
      if (registeredmdims > maxmdims) then
        call endrun('Too many dimensions in add_hist_coord.')
      end if
      if (len_trim(name) > max_hcoordname_len) then
        write(errormsg,'(a,i3,a)') 'History coord name exceeds the ',         &
             max_hcoordname_len, ' character length limit'
        call endrun(errormsg)
      end if
      hist_coords(registeredmdims)%name = trim(name)
      hist_coords(registeredmdims)%dimsize = 0
      hist_coords(registeredmdims)%long_name = ''
      hist_coords(registeredmdims)%units = ''
      hist_coords(registeredmdims)%integer_dim = .true.
      hist_coords(registeredmdims)%positive = ''
      hist_coords(registeredmdims)%standard_name = ''
      if (present(index)) then
        index = registeredmdims
      end if
    else
      if (present(index)) then
        index = i
      end if
    end if

  end subroutine add_hist_coord_regonly

  subroutine add_hist_coord_int(name, vlen, long_name, units, values,         &
       positive, standard_name)

    ! Input variables
    character(len=*), intent(in)                    :: name
    integer,          intent(in)                    :: vlen
    character(len=*), intent(in)                    :: long_name
    character(len=*), intent(in),          optional :: units
    integer,          intent(in),  target, optional :: values(:)
    character(len=*), intent(in),          optional :: positive
    character(len=*), intent(in),          optional :: standard_name

    ! Local variables
    integer                                         :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         i_values=values, positive=positive, standard_name=standard_name)
    ! Register the name if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), i)
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    ! Set the coord's values
    hist_coords(i)%dimsize = vlen
    if (len_trim(long_name) > max_chars) then
      if(masterproc) then
        write(iulog,*) 'WARNING: long_name for ',trim(name),' too long'
      end if
    end if
    hist_coords(i)%long_name = trim(long_name)
    if (present(units)) then
      hist_coords(i)%units = trim(units)
    else
      hist_coords(i)%units = ''
    end if
    hist_coords(i)%integer_dim = .true.
    if (present(values)) then
      hist_coords(i)%integer_values => values
    endif
    if (present(positive)) then
      hist_coords(i)%positive = trim(positive)
    end if
    if (present(standard_name)) then
      hist_coords(i)%standard_name = trim(standard_name)
    end if
    hist_coords(i)%vertical_coord = .false.

  end subroutine add_hist_coord_int

  subroutine add_hist_coord_r8(name, vlen, long_name, units, values,         &
       bounds_name, bounds, positive, standard_name, vertical_coord)

    ! Input variables
    character(len=*),      intent(in)                    :: name
    integer,               intent(in)                    :: vlen
    character(len=*),      intent(in)                    :: long_name
    character(len=*),      intent(in)                    :: units
    real(r8),              intent(in), target            :: values(:)
    character(len=*),      intent(in),          optional :: bounds_name
    real(r8),              intent(in), target,  optional :: bounds(:,:)
    character(len=*),      intent(in),          optional :: positive
    character(len=*),      intent(in),          optional :: standard_name
    logical,               intent(in),          optional :: vertical_coord

    ! Local variables
    character(len=120)                                   :: errormsg
    integer                                              :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         r_values=values, positive=positive, standard_name=standard_name,     &
         bounds_name=bounds_name, bounds=bounds)
    ! Register the name if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), i)
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    ! Set the coord's size
    hist_coords(i)%dimsize = vlen
    if (len_trim(long_name) > max_chars) then
      if(masterproc) then
        write(iulog,*) 'WARNING: long_name for ',trim(name),' too long'
      end if
    end if
    hist_coords(i)%long_name = trim(long_name)
    if (len_trim(units) > 0) then
      hist_coords(i)%units = trim(units)
    else
      hist_coords(i)%units = '1'
    end if
    hist_coords(i)%integer_dim = .false.
    hist_coords(i)%real_values => values
    if (present(positive)) then
      hist_coords(i)%positive = trim(positive)
    end if
    if (present(standard_name)) then
      hist_coords(i)%standard_name = trim(standard_name)
    end if
    if (present(bounds_name)) then
      hist_coords(i)%bounds_name = trim(bounds_name)
      if (.not. present(bounds)) then
        write(errormsg,*) 'bounds must be present for ',trim(bounds_name)
        call endrun(errormsg)
      end if
      hist_coords(i)%bounds => bounds
    else if (present(bounds)) then
      write(errormsg,*) 'bounds_name must be present for bounds values'
      call endrun(errormsg)
    else
      hist_coords(i)%bounds_name = ''
    end if
    if (present(vertical_coord)) then
      hist_coords(i)%vertical_coord = vertical_coord
    else
      hist_coords(i)%vertical_coord = .false.
    end if

  end subroutine add_hist_coord_r8

  subroutine add_vert_coord(name, vlen, long_name, units, values,            &
       positive, standard_name, formula_terms)

    ! Input variables
    character(len=*),      intent(in)                    :: name
    integer,               intent(in)                    :: vlen
    character(len=*),      intent(in)                    :: long_name
    character(len=*),      intent(in)                    :: units
    real(r8),              intent(in), target            :: values(:)
    character(len=*),      intent(in),          optional :: positive
    character(len=*),      intent(in),          optional :: standard_name
    type(formula_terms_t), intent(in),          optional :: formula_terms

    ! Local variable
    integer                                              :: i

    ! First, check to see if it is OK to add this coord
    i = check_hist_coord(name, vlen=vlen, long_name=long_name, units=units,   &
         r_values=values, positive=positive, standard_name=standard_name,     &
         formula_terms=formula_terms)
    ! Register the name and hist_coord values if necessary
    if (i == 0) then
      call add_hist_coord(trim(name), vlen, long_name, units, values,         &
           positive=positive, standard_name=standard_name,                    &
           vertical_coord=.true.)
      i = get_hist_coord_index(trim(name))
      !  if(masterproc) write(iulog,*) 'Registering hist coord',name,'(',i,') with length: ',vlen
    end if

    if (present(formula_terms)) then
      hist_coords(i)%formula_terms = formula_terms
    end if

  end subroutine add_vert_coord

  subroutine write_hist_coord_attr(File, mdimind, boundsdim, dimonly, mdimid)
    use pio, only: file_desc_t, var_desc_t, pio_put_att, pio_noerr,           &
                   pio_int, pio_double, pio_inq_varid, pio_def_var
    use cam_pio_utils, only: cam_pio_def_dim, cam_pio_def_var

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: mdimind        ! Internal dim index
    integer,           intent(in)    :: boundsdim      ! Bounds dimension ID
    logical,           intent(in)    :: dimonly        ! No def_var if .true.
    integer, optional, intent(out)   :: mdimid

    ! Local variables
    integer                          :: dimid          ! PIO dimension ID
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    character(len=120)               :: errormsg
    character(len=max_chars)         :: formula_terms  ! Constructed string
    integer                          :: ierr
    integer                          :: dtype
    logical                          :: defvar         ! True if var exists

    ! Check to see if the dimension already exists in the file
    call cam_pio_def_dim(File, trim(hist_coords(mdimind)%name),               &
           hist_coords(mdimind)%dimsize, dimid, existOK=.false.)
    ! If the caller wants to know the NetCDF dimension ID, set it here
    if (present(mdimid)) then
      mdimid = dimid
    end if
    if (.not. dimonly) then
      ! Time to define the variable (only if there are values)
      if (hist_coords(mdimind)%integer_dim) then
        dtype = pio_int
        defvar = associated(hist_coords(mdimind)%integer_values)
      else
        dtype = pio_double
        defvar = associated(hist_coords(mdimind)%real_values)
      end if
      if (defvar) then
        call cam_pio_def_var(File, trim(hist_coords(mdimind)%name), dtype,    &
             (/dimid/), vardesc, existOK=.false.)
        ! long_name
        ierr=pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%long_name))
        call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_hist_coord_attr')
        ! units
        if(len_trim(hist_coords(mdimind)%units) > 0) then
          ierr=pio_put_att(File, vardesc, 'units',                              &
               trim(hist_coords(mdimind)%units))
          call cam_pio_handle_error(ierr, 'Error writing "units" attr in write_hist_coord_attr')
        end if
        ! positive
        if(len_trim(hist_coords(mdimind)%positive) > 0) then
          ierr=pio_put_att(File, vardesc, 'positive',                           &
               trim(hist_coords(mdimind)%positive))
          call cam_pio_handle_error(ierr, 'Error writing "positive" attr in write_hist_coord_attr')
        end if
        ! standard_name
        if(len_trim(hist_coords(mdimind)%standard_name) > 0) then
          ierr=pio_put_att(File, vardesc, 'standard_name',                      &
               trim(hist_coords(mdimind)%standard_name))
          call cam_pio_handle_error(ierr, 'Error writing "standard_name" attr in write_hist_coord_attr')
        end if
        ! formula_terms
        if(len_trim(hist_coords(mdimind)%formula_terms%a_name) > 0) then
          write(formula_terms, '("a: ",a," b: ",a," p0: ",a," ps: ",a)')        &
               trim(hist_coords(mdimind)%formula_terms%a_name),                 &
               trim(hist_coords(mdimind)%formula_terms%b_name),                 &
               trim(hist_coords(mdimind)%formula_terms%p0_name),                &
               trim(hist_coords(mdimind)%formula_terms%ps_name)
          ierr=pio_put_att(File, vardesc, 'formula_terms', trim(formula_terms))
          call cam_pio_handle_error(ierr, 'Error writing "formula_terms" attr in write_hist_coord_attr')
        end if
        ! bounds
        if (associated(hist_coords(mdimind)%bounds)) then
          ! Write name of the bounds variable
          ierr=pio_put_att(File, vardesc, 'bounds', trim(hist_coords(mdimind)%bounds_name))
          call cam_pio_handle_error(ierr, 'Error writing "bounds" attr in write_hist_coord_attr')
        end if
      end if

      ! Now, we need to define and populate the associated bounds variable
      ! NB: Reusing vardesc, no longer assocated with main variable
      if (associated(hist_coords(mdimind)%bounds)) then
        if (size(hist_coords(mdimind)%bounds,2) /= hist_coords(mdimind)%dimsize) then
          ! If anyone hits this check, add a new dimension for this case
          write(errormsg, *) 'The bounds variable, ',                         &
               trim(hist_coords(mdimind)%bounds_name),                        &
               ', needs to have dimension (2,', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        call cam_pio_def_var(File, trim(hist_coords(mdimind)%bounds_name),    &
             pio_double, (/boundsdim,dimid/), vardesc, existOK=.false.)
      end if

      ! See if we have formula_terms variables to define
      ! Define the "a" variable name
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (associated(hist_coords(mdimind)%formula_terms%a_values)) then
        if (size(hist_coords(mdimind)%formula_terms%a_values) /= hist_coords(mdimind)%dimsize) then
          write(errormsg, *) 'The forumla_terms variable, ',                  &
               trim(hist_coords(mdimind)%formula_terms%a_name),               &
               ', needs to have dimension', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        call cam_pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%a_name), &
             pio_double, (/dimid/), vardesc, existOK=.false.)
        ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%a_long_name))
        call cam_pio_handle_error(ierr, 'Error writing "long_name" attr for "a" formula_term in write_hist_coord_attr')
      end if
      ! Define the "b" variable name
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (associated(hist_coords(mdimind)%formula_terms%b_values)) then
        if (size(hist_coords(mdimind)%formula_terms%b_values) /= hist_coords(mdimind)%dimsize) then
          write(errormsg, *) 'The forumla_terms variable, ',                  &
               trim(hist_coords(mdimind)%formula_terms%b_name),               &
               ', needs to have dimension', hist_coords(mdimind)%dimsize
          call endrun(errormsg)
        end if
        call cam_pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%b_name), &
             pio_double, (/dimid/), vardesc, existOK=.false.)
        ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%b_long_name))
        call cam_pio_handle_error(ierr, 'Error writing "long_name" attr for "b" formula_term in write_hist_coord_attr')
      end if
      ! Maybe define the p0 variable (this may be defined already which is OK)
      ! NB: Reusing vardesc, no longer assocated with previous variables
      if (hist_coords(mdimind)%formula_terms%p0_value /= fillvalue) then
        ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%p0_name), vardesc)
        if (ierr /= PIO_NOERR) then
          ierr = pio_def_var(File, trim(hist_coords(mdimind)%formula_terms%p0_name), &
               pio_double, vardesc)
          call cam_pio_handle_error(ierr, 'Unable to define "p0" formula_terms variable in write_hist_coord_attr')
          ierr = pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%formula_terms%p0_long_name))
          call cam_pio_handle_error(ierr, 'Error writing "long_name" attr for "p0" formula_term in write_hist_coord_attr')
          ierr = pio_put_att(File, vardesc, 'units', trim(hist_coords(mdimind)%formula_terms%p0_units))
          call cam_pio_handle_error(ierr, 'Error writing "units" attr for "p0" formula_term in write_hist_coord_attr')
        end if
      end if
      ! PS is not our responsibility
    end if ! (.not. dimonly)

  end subroutine write_hist_coord_attr

  !---------------------------------------------------------------------------
  !
  !  write_hist_coord_attrs
  !
  !  Write the dimension and coordinate attributes for the defined history
  !  coordinates.
  !
  !---------------------------------------------------------------------------
  subroutine write_hist_coord_attrs(File, boundsdim, mdimids, writemdims_in)
    use pio, only: file_desc_t, var_desc_t, pio_put_att,         &
                   pio_bcast_error, pio_internal_error, pio_seterrorhandling, &
                   pio_char
    use cam_pio_utils, only: cam_pio_def_dim, cam_pio_def_var

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: boundsdim      ! Bounds dimension ID
    integer, optional, allocatable, intent(out)   :: mdimids(:) ! NetCDF dim IDs
    logical, optional, intent(in)    :: writemdims_in  ! Write mdim variable

    ! Local variables
    integer                          :: i
    integer                          :: ierr
    integer                          :: dimids(2)      ! PIO dimension IDs
    logical                          :: writemdims     ! Define an mdim variable
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor

    if (present(mdimids)) then
      allocate(mdimids(registeredmdims))
    end if

    ! We will handle errors for this routine
    call pio_seterrorhandling(File, PIO_BCAST_ERROR)

    if (present(writemdims_in)) then
      writemdims = writemdims_in
    else
      writemdims = .false.
    end if

    ! NB: Currently, writemdims is for restart and we don't need to write
    ! these out in a history-restart file. This could change in the future.
    ! which would require a change to the function of the fourth argument
    ! Fill in the attribute information for each mdim
    do i = 1, registeredmdims
      if (present(mdimids)) then
        call write_hist_coord_attr(File, i, boundsdim, writemdims, mdimids(i))
      else
        call write_hist_coord_attr(File, i, boundsdim, writemdims)
      end if
    end do

    if (writemdims) then
      call cam_pio_def_dim(File, 'mdimslen', max_hcoordname_len, dimids(1),   &
           existOK=.true.)
      call cam_pio_def_dim(File, 'num_mdims', registeredmdims, dimids(2),     &
           existOK=.true.)
      call cam_pio_def_var(File, mdim_var_name, pio_char, dimids, vardesc,    &
           existOK=.false.)
      ierr = pio_put_att(File, vardesc, 'long_name', 'mdim dimension names')
      call cam_pio_handle_error(ierr, 'Error writing "long_name" attr for mdimnames in write_hist_coord_attrs')

    end if

    ! Back to I/O or die trying
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
  end subroutine write_hist_coord_attrs

  subroutine write_hist_coord_var(File, mdimind)
    use pio, only: file_desc_t, var_desc_t, pio_put_var, pio_inq_varid

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    integer,           intent(in)    :: mdimind        ! Internal dim index

    ! Local variables
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    integer                          :: ierr

    if ((hist_coords(mdimind)%integer_dim .and.                               &
         associated(hist_coords(mdimind)%integer_values)) .or.                &
         ((.not. hist_coords(mdimind)%integer_dim) .and.                      &
         associated(hist_coords(mdimind)%real_values))) then
      ! Check to make sure the variable already exists in the file
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%name), vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent dimension variable write_hist_coord_var')
      ! Write out the values for this dimension variable
      if (hist_coords(mdimind)%integer_dim) then
        ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%integer_values)
      else
        ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%real_values)
      end if
      call cam_pio_handle_error(ierr, 'Error writing variable values in write_hist_coord_var')
    end if

    ! Now, we need to possibly write values for the associated bounds variable
    if (associated(hist_coords(mdimind)%bounds)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%bounds_name), vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent bounds variable write_hist_coord_var')
    ! Write out the values for this bounds variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%bounds)
      call cam_pio_handle_error(ierr, 'Error writing bounds values in write_hist_coord_var')
    end if

    ! Write values for the "a" variable name
    if (associated(hist_coords(mdimind)%formula_terms%a_values)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%a_name), vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent "a" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "a" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%a_values)
      call cam_pio_handle_error(ierr, 'Error writing "a" formula_terms values in write_hist_coord_var')
    end if
    ! Write values for the "b" variable name
    if (associated(hist_coords(mdimind)%formula_terms%b_values)) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%b_name), vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent "b" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "b" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%b_values)
      call cam_pio_handle_error(ierr, 'Error writing "b" formula_terms values in write_hist_coord_var')
    end if
    ! Write values for the "p0" variable name (this may be an overwrite, too bad)
    if (hist_coords(mdimind)%formula_terms%p0_value /= fillvalue) then
      ! Check to make sure the variable already exists in the file
      ! NB: Reusing vardesc, no longer assocated with previous variables
      ierr = pio_inq_varid(File, trim(hist_coords(mdimind)%formula_terms%p0_name), vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent "p0" formula_terms variable write_hist_coord_var')
    ! Write out the values for this "p0" formula_terms variable
      ierr = pio_put_var(File, vardesc, hist_coords(mdimind)%formula_terms%p0_value)
      call cam_pio_handle_error(ierr, 'Error writing "p0" formula_terms values in write_hist_coord_var')
    end if

  end subroutine write_hist_coord_var

  subroutine write_hist_coord_vars(File, writemdims_in)
   use pio, only: file_desc_t, var_desc_t, pio_put_var,                       &
                  pio_bcast_error, pio_internal_error,                        &
                  pio_seterrorhandling, pio_inq_varid

    ! Input variables
    type(file_desc_t), intent(inout) :: File           ! PIO file Handle
    logical, optional, intent(in)    :: writemdims_in  ! Write mdim variable

    ! Local variables
    integer                          :: i
    integer                          :: ierr
    logical                          :: writemdims     ! Define an mdim variable
    type(var_desc_t)                 :: vardesc        ! PIO variable descriptor
    character(len=max_hcoordname_len), allocatable :: mdimnames(:)

    ! We will handle errors for this routine
    call pio_seterrorhandling(File, PIO_BCAST_ERROR)

    if (present(writemdims_in)) then
      writemdims = writemdims_in
    else
      writemdims = .false.
    end if

    if (writemdims) then
      allocate(mdimnames(registeredmdims))
    end if

    ! Write out the variable values for each mdim
    do i = 1, registeredmdims
      if (.not. writemdims) then
        ! NB: Currently, writemdims is for restart and we don't need to write
        ! these out in a history-restart file. This could change in the future
        ! which is why it is a separate if block
        ! Fill in the attribute information for each mdim
        call write_hist_coord_var(File, i)
      end if
      if (writemdims) then
        mdimnames(i) = trim(hist_coords(i)%name)
      end if
    end do

    if (writemdims) then
      ierr = pio_inq_varid(File, mdim_var_name, vardesc)
      call cam_pio_handle_error(ierr, 'Error writing values for nonexistent mdimnames variable in write_hist_coord_vars')
      ! Write out the values for mdim names
      ierr = pio_put_var(File, vardesc, mdimnames)
      call cam_pio_handle_error(ierr, 'Error writing values for mdimnames variable in write_hist_coord_vars')
      deallocate(mdimnames)
    end if

    ! Back to I/O or die trying
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

  end subroutine write_hist_coord_vars

  subroutine lookup_hist_coord_indices(mdimnames, mdimindicies)
    ! Dummy arguments
    character(len=*), intent(in) :: mdimnames(:)
    integer, intent(out) :: mdimindicies(:)

    ! Local variables
    integer :: i, j
    integer :: cnt
    character(len=120) :: errormsg
    character(len=16) :: name


    cnt = size(mdimnames)
    mdimindicies = -1


    do j=1,cnt
      name = mdimnames(j)
      do i = 1, registeredmdims
        if(name .eq. hist_coords(i)%name) then
          mdimindicies(j)=i
        end if
      end do
    end do
    do j = 1, cnt
      if(mdimindicies(j) < 0) then
        do i = 1, registeredmdims		
          print *,__FILE__,__LINE__,i,hist_coords(i)%name
        end do
        write(errormsg,*) 'Name ',mdimnames(j),' is not a registered history coordinate'
        call endrun(errormsg)
      end if
    end do

  end subroutine lookup_hist_coord_indices

  ! Find the vertical dimension (if present) in dimnames and return its size
  !    (which is the number of levels). Return -1 if not found
  !    If dimnames is not present, search all of the registered history coords
  integer function hist_coord_find_levels(dimnames) result(levels)
    ! Dummy argument
    character(len=*), optional, intent(in) :: dimnames(:)

    ! Local variables
    integer i, index, dimcnt

    levels = -1  ! Error return value

    if (present(dimnames)) then
      dimcnt = size(dimnames)
    else
      dimcnt = registeredmdims
    end if

    do i = 1, dimcnt
      if (present(dimnames)) then
        index = get_hist_coord_index(trim(dimnames(i)))
        if (i < 0) then
          call endrun('hist_coord_find_levels: '//trim(dimnames(i))//' is not a registred history coordinate')
        end if
      else
        index = i  ! Just cycle through all the registered mdims
      end if

      if (hist_coords(index)%vertical_coord) then
        levels = hist_coords(index)%dimsize
        exit
      end if
    end do

  end function hist_coord_find_levels

  !#######################################################################

  character(len=8) function sec2hms (seconds)

    ! Input arguments

    integer, intent(in) :: seconds

    ! Local workspace

    integer :: hours     ! hours of hh:mm:ss
    integer :: minutes   ! minutes of hh:mm:ss
    integer :: secs      ! seconds of hh:mm:ss

    if (seconds < 0 .or. seconds > 86400) then
      write(iulog,*)'SEC2HRS: bad input seconds:', seconds
      call endrun ()
    end if

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    if (minutes < 0 .or. minutes > 60) then
      write(iulog,*)'SEC2HRS: bad minutes = ',minutes
      call endrun ()
    end if

    if (secs < 0 .or. secs > 60) then
      write(iulog,*)'SEC2HRS: bad secs = ',secs
      call endrun ()
    end if

    write(sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)
    return
  end function sec2hms
  character(len=10) function date2yyyymmdd (date)

    ! Input arguments

    integer, intent(in) :: date

    ! Local workspace

    integer :: year    ! year of yyyy-mm-dd
    integer :: month   ! month of yyyy-mm-dd
    integer :: day     ! day of yyyy-mm-dd

    if (date < 0) then
      call endrun ('DATE2YYYYMMDD: negative date not allowed')
    end if

    year  = date / 10000
    month = (date - year*10000) / 100
    day   = date - year*10000 - month*100

    write(date2yyyymmdd,80) year, month, day
80  format(i4.4,'-',i2.2,'-',i2.2)
    return
  end function date2yyyymmdd

  !#######################################################################


end module cam_history_support
