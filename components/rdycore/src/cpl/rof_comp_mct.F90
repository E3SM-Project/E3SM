module rof_comp_mct

  ! !USES:

  use esmf
  use netcdf
  use pio
  use seq_infodata_mod
  use mct_mod
  use seq_flds_mod
  use seq_cdata_mod   , only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
  use seq_comm_mct    , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_strdata_mod , only: shr_strdata_type
  use shr_file_mod    , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only: shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only: shr_file_freeunit
  use shr_const_mod   , only: SHR_CONST_REARTH
  use rof_cpl_indices , only: rof_cpl_indices_set, &
                              index_r2x_Forr_rofl,    index_r2x_Forr_rofi,   &
                              index_r2x_Flrr_flood,   index_r2x_Flrr_volr,   &
                              index_r2x_Flrr_volrmch, index_r2x_Flrr_supply, &
                              index_r2x_Flrr_deficit

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer(IN)            :: mpicom_rof          ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: npes                ! number of tasks in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  integer                :: lsize    
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit_rof         ! logging unit number
  integer(IN)            :: rofid               ! mct comp id

  ! temporary for reading mosart file and making dummy domain
  integer :: nlatg, nlong, gsize       ! size of runoff data and number of grid cells
  real(r8), allocatable :: latc(:)     ! latitude of 1d grid cell (deg)
  real(r8), allocatable :: lonc(:)     ! longitude of 1d grid cell (deg)
  real(r8), allocatable :: areac(:)    ! area of 1d grid cell (deg)

  integer, allocatable :: start(:)     ! for gsmap initialization
  integer, allocatable :: length(:)    ! for gsmap initialization
  integer, allocatable :: pe_loc(:)    ! for gsmap initialization

  character(*), parameter :: F00   = "('(rof_comp_init) ',8a)"
  integer     , parameter :: master_task=0 ! task number of master task
  character(*), parameter :: subName = "(rof_init_mct) "

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: rof_init_mct
  !
  ! !DESCRIPTION:
  !     stub rof model init
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine rof_init_mct( EClock, cdata, x2r_r, r2x_r, NLFilename )

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2r_r, r2x_r
    character(len=*), optional  , intent(in)    :: NLFilename

    !EOP
    !-------------------------------------------------------------------------------

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap_rof
    type(mct_gGrid),         pointer :: dom_r        ! runoff model domain
    logical           :: rof_present               ! flag
    logical           :: rof_prognostic            ! flag
    logical           :: rofice_present            ! flag
    logical           :: flood_present             ! flag
    integer           :: shrlogunit                ! original log unit
    integer           :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    logical           :: exists                    ! true if file exists
    integer           :: ierr                      ! error code

    character(CL)     :: filePath ! generic file path
    character(CL)     :: fileName ! generic file name
    character(CS)     :: timeName ! domain file: time variable name
    character(CS)     :: lonName  ! domain file: lon  variable name
    character(CS)     :: latName  ! domain file: lat  variable name
    character(CS)     :: hgtName  ! domain file: hgt  variable name
    character(CS)     :: maskName ! domain file: mask variable name
    character(CS)     :: areaName ! domain file: area variable name

    character(*), parameter :: subName = "(rof_init_mct) "
    !-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, &
         id=rofid, &
         mpicom=mpicom_rof, &
         gsMap=gsMap_rof, &
         dom=dom_r, &
         infodata=infodata)

    ! Determine attribute vector indices
    call rof_cpl_indices_set()

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, &
         read_restart=read_restart)

    ! Determine instance information
    inst_name   = seq_comm_name(rofid)
    inst_index  = seq_comm_inst(rofid)
    inst_suffix = seq_comm_suffix(rofid)

    ! Determine communicator group
    call mpi_comm_rank(mpicom_rof, my_task, ierr)

    !--- open log file ---
    call shr_file_getLogUnit (shrlogunit)
    if (my_task == master_task) then
       inquire(file='rof_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          logunit_rof = shr_file_getUnit()
          call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),logunit_rof)
       end if
       write(logunit_rof,*) "RDycore model initialization"
    else
       logunit_rof = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit_rof)

    if (my_task == master_task) then
       write(logunit_rof,*) ' RDycore npes = ', npes
       write(logunit_rof,*) ' RDycore my_task  = ', my_task
       write(logunit_rof,*) ' inst_name = ', trim(inst_name)
    endif

    !----------------------------------------------------------------------------
    ! Read namelist file (TODO) - set values for now
    !----------------------------------------------------------------------------
    rof_present=.true.
    rofice_present=.false.
    rof_prognostic=.true.
    flood_present=.false.

    !----------------------------------------------------------------------------
    ! Initialize RDycore
    !----------------------------------------------------------------------------

!JW        if (rof_prognostic) then

       ! Initialize rof gsMap for ocean rof and land rof
       call rof_SetGSMap_mct( mpicom_rof, rofid, gsMap_rof)

       ! Initialize rof domain
       lsize = mct_gsMap_lsize(gsMap_rof, mpicom_rof)
       call rof_domain_mct( lsize, gsMap_rof, dom_r )

       ! Initialize cpl -> RDycore attribute vector
       call mct_aVect_init(x2r_r, rList=seq_flds_x2r_fields, lsize=lsize)
       call mct_aVect_zero(x2r_r)

       ! Initialize RDycore -> cpl attribute vector
       call mct_aVect_init(r2x_r, rList=seq_flds_r2x_fields, lsize=lsize)
       call mct_aVect_zero(r2x_r)

       ! Create mct river runoff export state
       call rof_export_mct( r2x_r )
!JW    endif

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from RDycore
    !----------------------------------------------------------------------------
    call seq_infodata_PutData( infodata, rof_present=rof_prognostic, &
            rof_nx = nlong, rof_ny = nlatg, rof_prognostic=rof_prognostic, &
            rofice_present=rofice_present, flood_present=flood_present)


    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit_rof,F00) 'rof_comp_init done'
    call shr_sys_flush(logunit_rof)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine rof_init_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: rof_run_mct
  !
  ! !DESCRIPTION:
  !     stub rof model run
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine rof_run_mct( EClock, cdata, x2r, r2x )

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r, r2x

    !EOP
    !-------------------------------------------------------------------------------
    call rof_export_mct( r2x )

  end subroutine rof_run_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: rof_final_mct
  !
  ! !DESCRIPTION:
  !     rof model finalize
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------
  !
  subroutine rof_final_mct( EClock, cdata, x2r, r2x)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r, r2x

    !EOP
    !-------------------------------------------------------------------------------

  end subroutine rof_final_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: rof_SetGSMap_mct
  !
  ! !DESCRIPTION:
  !     This routine sets up the RDYcore grid numbering for MCT
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------
  !
    subroutine rof_SetGSMap_mct( mpicom_rof, rofid, gsMap_rof)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Set the MCT GS map for the runoff model
    !
    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: mpicom_rof    ! MPI communicator for rof model
    integer        , intent(in)    :: rofid         ! Runoff model identifier
    type(mct_gsMap), intent(inout) :: gsMap_rof     ! MCT gsmap for runoff -> land data
    !
    ! LOCAL VARIABLES
    integer :: n                         ! indices
    integer :: ier                       ! error code
    character(len=32), parameter :: sub = 'rof_SetGSMap_mct'
    !-----------------------------------------------------

    ! TODO; get real mesh sizes
    call rof_read_mosart()

    ! get number of pes
    call mpi_comm_size(mpicom_rof, npes, ier)
    allocate(start(npes), length(npes), pe_loc(npes))
    start = 0
    length = 0
    pe_loc = 0

    ! for testing, simple 1D decomp of 2D grid
    do n = 1,npes
       length(n)  = gsize/npes
       if (n <= mod(gsize,npes)) length(n) = length(n) + 1
       if (n == 1) then
           start(n) = 1
       else
           start(n) = start(n-1) + length(n-1)
       endif
       pe_loc(n) = n-1
    enddo

     call mct_gsMap_init( gsMap_rof, rofid, npes, gsize, start, length, pe_loc)
    ! Init the gsMap with gindex
    ! call mct_gsMap_init( gsMap_rof, gindex, mpicom_r, rofid, lsize, gsize )

  end subroutine rof_SetGSMap_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: rof_domain_mct
  !
  ! !DESCRIPTION:
  !     This routine sets up the MCT domain for RDYcore
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------
  !
  subroutine rof_domain_mct( lsize, gsMap_rof, dom_rof )

    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_rof
    type(mct_gGrid), intent(inout) :: dom_rof
    !
    ! LOCAL VARIABLES
    integer :: n, ni              ! index
    integer :: lstart, lstop      ! index
    integer , pointer :: idata(:) ! temporary
    real(r8), pointer :: data(:)  ! temporary
    real(r8) :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)
    character(len=32), parameter :: sub = 'rof_domain_mct'
    !-----------------------------------------------------

    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    call mct_gGrid_init( GGrid=dom_rof, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Allocate memory
    allocate(data(lsize))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_rof, my_task, idata)
    call mct_gGrid_importIAttr(dom_rof,'GlobGridNum',idata,lsize)

    ! Initialize attribute vector with special value
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_rof,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_rof,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_rof,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_rof,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_rof,"mask" ,data,lsize)

    ! TODO - Fill in correct values for domain components
    lstart = start(my_task+1)
    lstop  = start(my_task+1) + length(my_task+1) - 1

    ni = 0
    do n = lstart, lstop
       ni = ni + 1
       data(ni) = lonc(n)
    end do
    call mct_gGrid_importRattr(dom_rof,"lon",data,lsize)

    ni = 0
    do n = lstart, lstop
       ni = ni + 1
       data(ni) = latc(n)
    end do
    call mct_gGrid_importRattr(dom_rof,"lat",data,lsize)

    ni = 0
    do n = lstart, lstop
       ni = ni + 1
       data(ni) = areac(n)*1.0e-6_r8/(re*re)
    end do
    call mct_gGrid_importRattr(dom_rof,"area",data,lsize)

    ni = 0
    do n = lstart, lstop
       ni = ni + 1
       data(ni) = 1.0_r8
    end do
    call mct_gGrid_importRattr(dom_rof,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_rof,"frac",data,lsize)

    deallocate(start,length,pe_loc)
    deallocate(data)
    deallocate(idata)

  end subroutine rof_domain_mct

  !====================================================================================

  subroutine rof_import_mct( x2r_r)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Obtain the runoff input from the coupler
    ! convert from kg/m2s to m3/s
    !
    ! ARGUMENTS:
    implicit none
    type(mct_aVect), intent(inout) :: x2r_r
    !
    ! LOCAL VARIABLES
    character(len=32), parameter :: sub = 'rof_import_mct'
    !---------------------------------------------------------------------------

  end subroutine rof_import_mct

!====================================================================================

  subroutine rof_export_mct( r2x_r )

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Send the runoff model export state to the coupler
    ! convert from m3/s to kg/m2s
    !
    ! ARGUMENTS:
    implicit none
    type(mct_aVect), intent(inout) :: r2x_r  ! Runoff to coupler export state
    !
    ! LOCAL VARIABLES
    integer :: ni, n, nt, nliq, nfrz
    logical,save :: first_time = .true.
    character(len=32), parameter :: sub = 'rof_export_mct'
    real(R8) :: tmp1
    !---------------------------------------------------------------------------

      do n = 1, lsize
        r2x_r%rattr(index_r2x_Forr_rofl,n)    = float(my_task)
        r2x_r%rattr(index_r2x_Forr_rofi,n)    = float(my_task)
        r2x_r%rattr(index_r2x_Flrr_flood,n)   = float(my_task)
        r2x_r%rattr(index_r2x_Flrr_volr,n)    = float(my_task)
        r2x_r%rattr(index_r2x_Flrr_volrmch,n) = float(my_task)
        r2x_r%rattr(index_r2x_Flrr_supply,n)  = float(my_task)
        r2x_r%rattr(index_r2x_Flrr_deficit,n) = float(my_task)
      enddo

  end subroutine rof_export_mct

!====================================================================================

  subroutine  rof_read_mosart()

    use shr_mpi_mod
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Send the runoff model export state to the coupler
    ! convert from m3/s to kg/m2s
    !
    ! ARGUMENTS:
    implicit none
    !
    ! LOCAL VARIABLES
    integer :: dimid
    integer :: ncid
    integer :: varid, areavarid
    integer :: status
    integer :: i, j, count, ier
    integer :: lat, lon
    character(len=128) :: filename_rof
    integer, parameter :: RKIND = selected_real_kind(13)
    real(kind=RKIND), dimension(:),   allocatable :: lat1D, lon1D, area1D
    real(kind=RKIND), dimension(:,:), allocatable :: area

    !---------------------------------------------------------------------------

    ! open mosart file
    filename_rof = '/global/cfs/cdirs/e3sm/inputdata/rof/mosart/MOSART_global_half_20180721a.nc'

    ! read file only from master_task
    if (my_task == master_task) then
       status = nf90_open(trim(filename_rof), NF90_NOWRITE, ncid)

       ! read in data dimensions
       status = nf90_inq_dimid(ncid, "lat", dimid)
       if (status == PIO_NOERR) then
         status = nf90_inquire_dimension(ncid, dimid, len=lat)
         status = nf90_inq_dimid(ncid, "lon", dimid)
         if (status /= PIO_NOERR) then
            call shr_sys_abort("ERROR: lat dimension defined in the file, but lon dimension is missing. file " // trim(filename_rof))
         endif
         status = nf90_inquire_dimension(ncid, dimid, len=lon)
       else
          status = nf90_inq_dimid(ncid, "gridcell", dimid)
          if (status == PIO_NOERR) then
             status = nf90_inquire_dimension(ncid, dimid, len=lon)
             lat = 1
          else
             call shr_sys_abort("ERROR: Neither (lon,lat) nor (gridcell) dimension defined in the following file: " // trim(filename_rof))
          endif
       endif
    end if

    ! broadcast dimensions
    call shr_mpi_bcast (lat, mpicom_rof, 'rof_read_mosart')
    call shr_mpi_bcast (lon, mpicom_rof, 'rof_read_mosart')

    ! set dimensions on all processors
    gsize = lat*lon
    nlatg = lat
    nlong = lon

    ! allocate arrays
    allocate(lat1D(lat))
    allocate(lon1D(lon))
    allocate(area(lat,lon))
    allocate(lonc(gsize))
    allocate(latc(gsize))
    allocate(areac(gsize))

    ! read in lat, lon, and area only on master_task
    if (my_task == master_task) then
       allocate(area1D(lat))
       status = nf90_inq_varid(ncid, "lat", varid)
       status = nf90_get_var(ncid, varid, lat1D)

       status = nf90_inq_varid(ncid, "lon", varid)
       status = nf90_get_var(ncid, varid, lon1D)

       status = nf90_inq_varid(ncid, "area", areavarid)
       do j = 1, nlong
          status = nf90_get_var(ncid, areavarid, area1D, start=(/1,j/), count=(/lat,1/))
          area(:,j) = area1d(:)
       enddo

       call shr_file_freeUnit(ncid)
       deallocate(area1D)
    endif

    ! broadcast lat1D, lon1D and arear
    call shr_mpi_bcast (lat1D, mpicom_rof, 'rof_read_mosart')
    call shr_mpi_bcast (lon1D, mpicom_rof, 'rof_read_mosart')
    call shr_mpi_bcast (area,  mpicom_rof, 'rof_read_mosart')

    ! load 2d data into 1d arrays
    count = 0
    do j = 1, nlong
       do i = 1, nlatg
         count = count + 1
         latc(count)  = lat1D(i)
         lonc(count)  = lon1D(j)
         areac(count) = area(i,j)
       end do
    end do

    if (my_task == master_task) then
       write(logunit_rof,*) 'Read lat1D ', minval(lat1D), maxval(lat1D)
       write(logunit_rof,*) 'Read lon1D ', minval(lon1D), maxval(lon1D)
       write(logunit_rof,*) 'Read area  ', minval(area),  maxval(area)
    end if

    deallocate(lat1D)
    deallocate(lon1D)
    deallocate(area)

    return

  end subroutine rof_read_mosart
!====================================================================================

end module rof_comp_mct
