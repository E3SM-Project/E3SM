module BeTR_GridMod

  use bshr_kind_mod    , only : r8 => shr_kind_r8
  use babortutils      , only : endrun
  use bshr_log_mod     , only : errMsg => shr_log_errMsg
  use betr_constants   , only : betr_filename_length
  use betr_constants   , only : betr_string_length, betr_string_length_long
  use betr_constants   , only : betr_namelist_buffer_size
  use betr_constants   , only : stdout

  use betr_varcon, only : bspval


  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  character(len=*), parameter          :: uniform_str = 'uniform'
  character(len=*), parameter          :: clm_str = 'clm'
  character(len=*), parameter          :: dataset_str = 'dataset'

  integer, parameter :: uniform_grid = 1
  integer, parameter :: clm_grid = 2
  integer, parameter :: dataset_grid = 3

  integer, parameter :: num_columns = 1

  type, public :: betr_grid_type
     character(len=betr_filename_length) :: grid_data_filename
     character(len=betr_string_length)   :: grid_data_format ! file format: netcdf, namelist, csv, etc.
     character(len=betr_string_length)   :: grid_type_str ! uniform, clm
     integer :: grid_type

     integer,  public          :: nlevgrnd
     real(r8), public          :: delta_z
     real(r8), public, pointer :: zsoi(:)  => null() !soil depth, node center 1 : nlevsoi
     real(r8), public, pointer :: zisoi(:) => null() !soil depth, interface,  0 : nlevsoi
     real(r8), public, pointer :: dzsoi(:) => null() !soil layer thickness

     real(r8), public, pointer :: bsw(:) => null() ! clap-hornberg parameter
     real(r8), public, pointer :: watsat(:) => null() ! saturated volumetric water content
     real(r8), public, pointer :: sucsat(:)=> null()
     real(r8), public, pointer :: pctsand(:)=> null()
   contains
     procedure, public  :: Init
     procedure, public  :: ReadNamelist
     procedure, public  :: ReadNetCDFData
     procedure, private :: InitAllocate
     procedure, private :: uniform_vertical_grid
     procedure, private :: clm_exponential_vertical_grid
     procedure, private :: set_interface_depths

  end type betr_grid_type


contains

  ! ---------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer)

    implicit none

    class(betr_grid_type),                    intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in)    :: namelist_buffer

    call this%ReadNameList(namelist_buffer)
    call this%InitAllocate()
    select case (trim(this%grid_type_str))
       case (clm_str)
          this%grid_type = clm_grid
          call this%clm_exponential_vertical_grid()
       case (uniform_str)
          this%grid_type = uniform_grid
          call this%uniform_vertical_grid()
       case (dataset_str)
          this%grid_type = dataset_grid
          ! grid will be initialized by reading the dataset
       case default
          this%grid_type = clm_grid
          call this%clm_exponential_vertical_grid()
          write(*, *) 'WARNING: no grid data type specified, using clm.'
       end select

    ! select read routine based on data format.
    call this%ReadNetCDFData()

    write(*, *) 'dump grid: '
    write(*, *) 'dzsoi = ', this%dzsoi
    write(*, *) 'zsoi = ', this%zsoi
    write(*, *) 'zisoi = ', this%zisoi
    write(*, *) 'bsw = ', this%bsw
    write(*, *) 'watsat = ', this%watsat
    write(*, *) 'sucsat = ', this%sucsat
  end subroutine Init

  ! ---------------------------------------------------------------------------
  subroutine InitAllocate(this)

    implicit none

    class(betr_grid_type), intent(inout) :: this

    allocate(this%zsoi(1:this%nlevgrnd))
    this%zsoi = bspval

    allocate(this%zisoi(0:this%nlevgrnd))
    this%zisoi = bspval

    allocate(this%dzsoi(1:this%nlevgrnd))
    this%dzsoi = bspval

    allocate(this%bsw(1:this%nlevgrnd))
    this%bsw = bspval

    allocate(this%watsat(1:this%nlevgrnd))
    this%watsat = bspval

    allocate(this%sucsat(1:this%nlevgrnd))
    this%sucsat = bspval

    allocate(this%pctsand(1:this%nlevgrnd))
    this%pctsand = bspval
  end subroutine InitAllocate

  ! ---------------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer)

    use betr_constants, only : betr_namelist_buffer_size
    use betr_constants, only : betr_filename_length
    use betr_constants, only : betr_string_length_long

    implicit none

    class(betr_grid_type),                    intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in)    :: namelist_buffer

    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'ReadNameList'
    character(len=betr_string_length)      :: grid_data_format, grid_type_str
    character(len=betr_filename_length)    :: grid_data_filename
    character(len=betr_string_length_long) :: ioerror_msg
    integer                                :: nlevgrnd
    real(r8)                               :: delta_z

    !-----------------------------------------------------------------------

    namelist / betr_grid /                                    &
         grid_type_str, grid_data_filename, grid_data_format, &
         nlevgrnd, delta_z

    grid_data_format   = ''
    grid_data_filename = ''
    grid_type_str      = clm_str
    nlevgrnd           = 15 ! default to clm grid
    delta_z            = bspval

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=betr_grid, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          write(*, *) ioerror_msg
          call endrun(msg="ERROR reading betr_grid namelist "//errmsg(mod_filename, __LINE__))
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr grid :'
       write(stdout, *)
       write(stdout, *) ' betr_grid namelist settings :'
       write(stdout, *)
       write(stdout, betr_grid)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%grid_type_str = trim(grid_type_str)
    this%grid_data_format = trim(grid_data_format)
    this%grid_data_filename = trim(grid_data_filename)
    this%nlevgrnd = nlevgrnd
    this%delta_z = delta_z

  end subroutine ReadNamelist

  !------------------------------------------------------------------------
  subroutine ReadNetCDFData(this)
    !DESCRIPTION
    !read netcdf data
    use ncdio_pio, only : file_desc_t
    use ncdio_pio, only : ncd_nowrite
    use ncdio_pio, only : ncd_pio_openfile
    use ncdio_pio, only : ncd_getvar
    use ncdio_pio, only : ncd_pio_closefile
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    character(len=250)    :: ncf_in_filename
    type(file_desc_t)     :: ncf_in
    real(r8), allocatable :: data(:,:)
    integer               :: j


    ncf_in_filename = trim(this%grid_data_filename)
    call ncd_pio_openfile(ncf_in, ncf_in_filename, mode=ncd_nowrite)

    allocate(data(num_columns, this%nlevgrnd))

    call ncd_getvar(ncf_in, 'WATSAT', data)
    do j = 1, this%nlevgrnd
       this%watsat(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'BSW', data)
    do j = 1, this%nlevgrnd
       this%bsw(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'SUCSAT', data)
    do j = 1, this%nlevgrnd
       this%sucsat(j) = data(num_columns, j)
    enddo

    call ncd_getvar(ncf_in, 'PCTSAND', data)
    do j = 1, this%nlevgrnd
      this%pctsand(j) = data(num_columns, j)
    enddo

    if (this%grid_type == dataset_grid) then
       call ncd_getvar(ncf_in, 'DZSOI', data)
       do j = 1, this%nlevgrnd
          this%dzsoi(j) = data(num_columns, j)
       end do

       call ncd_getvar(ncf_in, 'ZSOI', data)
       do j = 1, this%nlevgrnd
          this%zsoi(j) = data(num_columns, j)
       end do

       call this%set_interface_depths()

    end if

    call ncd_pio_closefile(ncf_in)

    deallocate(data)

  end subroutine ReadNetCDFData

  ! ---------------------------------------------------------------------------
  subroutine uniform_vertical_grid(this)
    !DESCRIPTION
    !set uniform thickness grid
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variable
    integer :: j

    if (this%delta_z == bspval) then
       call endrun(msg="ERROR reading betr_grid namelist must specify delta_z. "//errmsg(mod_filename, __LINE__))
    end if

    ! thickness b/n two interfaces
    this%dzsoi(:) = this%delta_z

    ! node depths
    do j = 1, this%nlevgrnd
       this%zsoi(j) = (real(j, r8) - 0.5_r8) * this%dzsoi(j)
    enddo

    call this%set_interface_depths()

  end subroutine uniform_vertical_grid

  ! ---------------------------------------------------------------------------
  subroutine clm_exponential_vertical_grid(this)
    !
    ! DESCRIPTION
    ! initialize the clm exporential vertical grid for computation
    !
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    real(r8)                             :: scalez = 0.025_r8 ! Soil layer thickness discretization (m)
    integer                              :: j

    ! node depths
    do j = 1, this%nlevgrnd
       this%zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)
    enddo

    ! thickness b/n two interfaces
    this%dzsoi(1) = 0.5_r8*(this%zsoi(1) + this%zsoi(2))
    do j = 2,this%nlevgrnd-1
       this%dzsoi(j) = 0.5_r8*(this%zsoi(j+1) - this%zsoi(j-1))
    enddo
    this%dzsoi(this%nlevgrnd) = this%zsoi(this%nlevgrnd) - this%zsoi(this%nlevgrnd-1)

    call this%set_interface_depths()

  end subroutine clm_exponential_vertical_grid

  ! ---------------------------------------------------------------------------
  subroutine set_interface_depths(this)
    !DESCRIPTION
    !set node depth
    !USES
    use bshr_kind_mod , only: r8 => shr_kind_r8
    implicit none
    !argument
    class(betr_grid_type), intent(inout) :: this
    !temporary variables
    integer :: j

    this%zisoi(0) = 0._r8
    do j = 1, this%nlevgrnd-1
       this%zisoi(j) = 0.5_r8*(this%zsoi(j) + this%zsoi(j+1))
    enddo
    this%zisoi(this%nlevgrnd) = this%zsoi(this%nlevgrnd) + 0.5_r8*this%dzsoi(this%nlevgrnd)
  end subroutine set_interface_depths

  ! ---------------------------------------------------------------------------

end module BeTR_GridMod
