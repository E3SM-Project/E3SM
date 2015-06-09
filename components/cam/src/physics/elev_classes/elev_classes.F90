!
!-------------------------------------------------------------------------------
! initialization code for elevation classes
!-------------------------------------------------------------------------------
module elev_classes

  use shr_kind_mod,   only: r8 => shr_kind_r8, SHR_KIND_CL

  implicit none
  private
  save

  ! Public module variables
  logical, public, protected :: ec_active = .false.
  integer, public, protected :: max_elevation_classes = 1

  ! Private module variables
  character(len=SHR_KIND_CL) :: elevation_classes_filename

  ! Public interface functions
  public :: elevation_classes_readnl
  public :: elevation_classes_init

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine elevation_classes_readnl(NLFilename)
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: masterproc, masterprocid, mpicom
    use spmd_utils,     only: mpi_integer, mpi_logical, mpi_character
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog

    ! Dymmy arguments
    character(len=*),    intent(in)  :: NLFileName

    ! Local variables
    integer                          :: ierr
    integer                          :: unitn
    character(len=SHR_KIND_CL)       :: elevation_classes
    character(len=*), parameter      :: subname = "ELEVATION_CLASSES_READNL"

    namelist /ec_nl/ elevation_classes, ec_active

    ! Default is no elevation classes
    ec_active = .false.
    elevation_classes = ''

    ! Read the namelist
    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(NLFilename), status='old')
      call find_group_name(unitn, 'ec_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, ec_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname//': ERROR reading ec_nl namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)

      elevation_classes_filename = elevation_classes

      if (ec_active .and. (len_trim(elevation_classes_filename) == 0)) then
        call endrun(subname//': elevation_classes filename required')
      end if
      if (ec_active) then
        write(iulog, *) 'Elevation Classes will be read from ',trim(elevation_classes_filename)
      else
        write(iulog, *) 'Elevation Classes are disabled'
      end if
    end if

    call mpi_bcast(elevation_classes_filename, len(elevation_classes_filename), mpi_character, masterprocid, mpicom, ierr)
      
  end subroutine elevation_classes_readnl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  elevation_classes_init is called from phys_grid_init
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine elevation_classes_init(ldof, num_subgrids, subgrid_area, subgrid_elev)
    use cam_logfile,   only: iulog
    use spmd_utils,    only: masterproc
    use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile
    use cam_pio_utils, only: cam_pio_handle_error
    use pio,           only: PIO_NOWRITE, PIO_inq_dimid, PIO_inq_dimlen
    use pio,           only: PIO_seterrorhandling, PIO_read_darray
    use pio,           only: PIO_initdecomp, PIO_freedecomp, PIO_inq_varid
    use pio,           only: file_desc_t, var_desc_t, io_desc_t, iosystem_desc_t
    use pio,           only: PIO_INT, PIO_DOUBLE
    use pio,           only: PIO_OFFSET_KIND, PIO_BCAST_ERROR
    use cam_instance,  only: atm_id
    use shr_pio_mod,   only: shr_pio_getiosys

    !! Dummy arguments
    integer(PIO_OFFSET_KIND), intent(in)  :: ldof(:)    ! Map to file order
    integer,     allocatable, intent(out) :: num_subgrids(:)
    real(r8),    allocatable, intent(out) :: subgrid_area(:,:)
    real(r8),    allocatable, intent(out) :: subgrid_elev(:,:)

    !! Local variables
    type(file_desc_t)                  :: fh_ec
    integer                            :: larray_size ! Size of dest arrays
    integer                            :: ierr
    integer                            :: dimid
    integer                            :: grid_len
    integer                            :: err_handling
    type(io_desc_t)                    :: iodesc
    type(var_desc_t)                   :: varid
    type(iosystem_desc_t), pointer     :: piosys
    character(len=*), parameter        :: subname = 'ELEVATION_CLASSES_INIT'

    call cam_pio_openfile(fh_ec, trim(elevation_classes_filename), PIO_NOWRITE)
    ! We will handle errors for this routine
    call pio_seterrorhandling(fh_ec, PIO_BCAST_ERROR, err_handling)

    ! We need to know the maximum number of elevation classes and grid size
    ierr = PIO_inq_dimid(fh_ec, 'grid_size', dimid)
    call cam_pio_handle_error(ierr, subname//': Error finding dimension, grid_size')
    ierr = PIO_inq_dimlen(fh_ec, dimid, grid_len)
    ierr = PIO_inq_dimid(fh_ec, 'MaxNoClass', dimid)
    call cam_pio_handle_error(ierr, subname//': Error finding dimension, MaxNoClass')
    ierr = PIO_inq_dimlen(fh_ec, dimid, max_elevation_classes)

    larray_size = size(ldof)
    if (allocated(num_subgrids)) then
      deallocate(num_subgrids)
    end if
    allocate(num_subgrids(larray_size))
    num_subgrids = 0

    ! Read the relevant variables
    piosys => shr_pio_getiosys(atm_id)
    call PIO_initdecomp(piosys, PIO_INT, (/ grid_len/), ldof, iodesc)
    ierr = PIO_inq_varid(fh_ec, 'NumSubgrids', varid)
    call cam_pio_handle_error(ierr, subname//': Error finding variable, NumSubgrids')
    call pio_read_darray(fh_ec, varid, iodesc, num_subgrids, ierr)
    call cam_pio_handle_error(ierr, subname//': Error reading NumSubgrids')
    call PIO_freedecomp(fh_ec, iodesc)

    if (allocated(subgrid_area)) then
      deallocate(subgrid_area)
    end if
    allocate(subgrid_area(max_elevation_classes, larray_size))
    subgrid_area = 0.0_r8
    call PIO_initdecomp(piosys, PIO_DOUBLE, (/ max_elevation_classes, grid_len/), ldof, iodesc)
    ierr = PIO_inq_varid(fh_ec, 'SubgridAreaFrac', varid)
    call cam_pio_handle_error(ierr, subname//': Error finding variable, SubgridAreaFrac')
    call pio_read_darray(fh_ec, varid, iodesc, subgrid_area, ierr)
    call cam_pio_handle_error(ierr, subname//': Error reading SubgridAreaFrac')

    if (allocated(subgrid_elev)) then
      deallocate(subgrid_elev)
    end if
    allocate(subgrid_elev(max_elevation_classes, larray_size))
    subgrid_elev = 0.0_r8
    ierr = PIO_inq_varid(fh_ec, 'AveSubgridElv', varid)
    call cam_pio_handle_error(ierr, subname//': Error finding variable, AveSubgridElv')
    call pio_read_darray(fh_ec, varid, iodesc, subgrid_elev, ierr)
    call cam_pio_handle_error(ierr, subname//': Error reading AveSubgridElv')
    call PIO_freedecomp(fh_ec, iodesc)

    ! Back to old error handling
    call pio_seterrorhandling(fh_ec, err_handling)
    ! Finally, close the file
    call cam_pio_closefile(fh_ec)
    if (masterproc) then
      write(iulog, *) subname, ': Read elevation classes from ', trim(elevation_classes_filename)
    end if
  end subroutine elevation_classes_init

end module elev_classes
