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
  use perf_mod        , only: t_startf, t_stopf, t_barrierf
  use rof_cpl_indices , only: rof_cpl_indices_set, &
                              index_r2x_Forr_rofl,    index_r2x_Forr_rofi,   &
                              index_r2x_Flrr_flood,   index_r2x_Flrr_volr,   &
                              index_r2x_Flrr_volrmch, index_r2x_Flrr_supply, &
                              index_r2x_Flrr_deficit, index_x2r_So_ssh
  use rdycoreMod      , only: inst_name, inst_suffix, inst_index
  use rdycoreSpmdMod  , only: masterproc, mpicom_rof, iam, npes, rofid, RDycoreSpmdInit, &
                              MPI_REAL8, MPI_INTEGER, MPI_CHARACTER, MPI_LOGICAL, MPI_SUM
  use RDycoreIO

  use rdycoreMod     , only: rdycore_init
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                :: lsize    

  ! temporary for reading rdycore file and making dummy domain
  integer :: nlatg, nlong, gsize         ! size of runoff data and number of grid cells

  logical :: isgrid2d                    ! true if the grid is 2D, else false

  ! global cell = all cells in the mesh
  real(r8), allocatable :: latc_g(:)     ! global latitude of 1d grid cell (deg)
  real(r8), allocatable :: lonc_g(:)     ! global longitude of 1d grid cell (deg)
  real(r8), allocatable :: areac_g(:)    ! global area of 1d grid cell (deg)

  ! local cell = cells owned by each MPI rank
  real(r8), allocatable :: latc_l(:)     ! local latitude of 1d grid cell (deg)
  real(r8), allocatable :: lonc_l(:)     ! local longitude of 1d grid cell (deg)
  real(r8), allocatable :: areac_l(:)    ! local area of 1d grid cell (deg)

  integer, allocatable :: start(:)     ! for gsmap initialization
  integer, allocatable :: length(:)    ! for gsmap initialization
  integer, allocatable :: pe_loc(:)    ! for gsmap initialization

  logical                 :: do_rdycore

  character(len=256)      :: rdycore_yaml_file
  character(len=256)      :: filename_rof
  character(*), parameter :: F00   = "('(rof_comp_init) ',8a)"
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
    integer           :: shrlogunit                ! original log unit
    integer           :: shrloglev                 ! original log level
    integer           :: ierr                      ! error code
    integer           :: mpicom_loc                ! local mpi communicator
    logical           :: read_restart              ! start from restart
    logical           :: exists                    ! true if file exists
    logical           :: rof_present               ! flag
    logical           :: rof_prognostic            ! flag
    logical           :: rofice_present            ! flag
    logical           :: flood_present             ! flag

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
         mpicom=mpicom_loc, &
         gsMap=gsMap_rof, &
         dom=dom_r, &
         infodata=infodata)

    ! Determine attribute vector indices
    call rof_cpl_indices_set()

    ! Initialize RDycore MPI communicator
    call RDycoreSpmdInit(mpicom_loc)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, &
         read_restart=read_restart)

    ! Determine instance information
    inst_name   = seq_comm_name(rofid)
    inst_index  = seq_comm_inst(rofid)
    inst_suffix = seq_comm_suffix(rofid)

    ! Initialize RDycore pio
    call ncd_pio_init()

    !--- open log file ---
    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
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

    if (masterproc) then
       write(logunit_rof,*) ' RDycore npes = ', npes
       write(logunit_rof,*) ' RDycore iam  = ', iam
       write(logunit_rof,*) ' inst_name = ', trim(inst_name)
    endif

    !----------------------------------------------------------------------------
    ! Read namelist file (TODO) - set values for now
    !----------------------------------------------------------------------------
    rof_present=.true.
    rofice_present=.false.
    rof_prognostic=.true.
    flood_present=.false.
    call rof_read_namelist()

    call rdycore_init(logunit_rof)

    !----------------------------------------------------------------------------
    ! Initialize RDycore
    !----------------------------------------------------------------------------

!JW        if (rof_prognostic) then

       ! Initialize rof gsMap for ocean rof and land rof
       !call rof_SetGSMap_mct( mpicom_rof, rofid, gsMap_rof)
       !lsize = mct_gsMap_lsize(gsMap_rof, mpicom_rof)
       call rof_SetGSMap_From_RDycore_mct( rofid, gsMap_rof)

       ! Initialize rof domain
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

    if (masterproc) write(logunit_rof,F00) 'rof_comp_init done'
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

    ! Map MCT datatype to rof data
    call t_startf ('lc_rof_import')
    call rof_import_mct( x2r )
    call t_stopf ('lc_rof_import')

    !TODO: run rdycore code

    ! Map rof data to MCT datatype
    call t_startf ('lc_rof_export')
    call rof_export_mct( r2x )
    call t_stopf ('lc_rof_export')

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
    call rof_read_rdycore()

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
  subroutine rof_SetGSMap_From_RDycore_mct(rofid, gsMap_rof)
   !
   use rdycoreMod, only : num_cells_owned, num_cells_global, natural_id_cells_owned
   !
   implicit none
   !
   integer        , intent(in)    :: rofid         ! Runoff model identifier
   type(mct_gsMap), intent(inout) :: gsMap_rof     ! MCT gsmap for runoff -> land data
    !
    ! LOCAL VARIABLES
   integer :: i, ibeg, iend             ! indices
   integer :: ier                       ! error code
   integer, allocatable :: gindex(:)
   character(len=32), parameter :: sub = 'rof_SetGSMap_From_RDycore_mct'

   gsize = num_cells_global
   lsize = num_cells_owned

    ! TODO; get real mesh sizes
   call rof_read_rdycore()

   ! get number of pes
   call mpi_comm_size(mpicom_rof, npes, ier)

   ! allocate memory
   allocate(start(npes), length(npes), pe_loc(npes))
   start = 0
   length = 0
   pe_loc = 0

   ! set only values corresponding to my rank
   length(iam + 1) = lsize
   pe_loc(iam + 1) = iam

   call MPI_Scan(lsize, start(iam + 1), 1, MPI_INTEGER, MPI_SUM, mpicom_rof, ier)
   start(iam + 1) = start(iam + 1) + 1 - lsize

   ibeg = start(iam + 1)
   iend = ibeg + lsize - 1

   ! set the indices that are locally owned
   allocate(gindex(ibeg:iend))
   do i = ibeg, iend
      gindex(i) = natural_id_cells_owned(i - ibeg + 1) + 1 ! converting 0-based IDs to 1-based
   enddo

   ! create the gsMap
   call mct_gsMap_init( gsMap_rof, gindex, mpicom_rof, rofid, lsize, gsize )

   ! free up memory
   deallocate(gindex)

   end subroutine rof_SetGSMap_From_RDycore_mct


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

    use rdycoreMod, only : natural_id_cells_owned
    !
    ! ARGUMENTS:
   use rdycoreMod, only : natural_id_cells_owned
   !
    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_rof
    type(mct_gGrid), intent(inout) :: dom_rof
    !
    ! LOCAL VARIABLES
    integer :: n, ni              ! index
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
    call mct_gsMap_orderedPoints(gsMap_rof, iam, idata)
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
    if (isgrid2d) then
       ni = 0
       do n = 1, lsize
          ni = natural_id_cells_owned(n) + 1
          data(n) = lonc_g(ni)
       end do
       call mct_gGrid_importRattr(dom_rof,"lon",data,lsize)

       ni = 0
       do n = 1, lsize
          ni = natural_id_cells_owned(n) + 1
          data(n) = latc_g(ni)
       end do
       call mct_gGrid_importRattr(dom_rof,"lat",data,lsize)

       do n = 1, lsize
          ni = natural_id_cells_owned(n) + 1
          data(n) = areac_g(ni)*1.0e-6_r8/(re*re)
       end do
       call mct_gGrid_importRattr(dom_rof,"area",data,lsize)

    else

       data(:) = lonc_l(:)
       call mct_gGrid_importRattr(dom_rof,"lon",data,lsize)

       data(:) = latc_l(:)
       call mct_gGrid_importRattr(dom_rof,"lat",data,lsize)

       do n = 1, lsize
          data(n) = areac_l(n)*1.0e-6_r8/(re*re)
       end do
       call mct_gGrid_importRattr(dom_rof,"area",data,lsize)
    end if

    do n = 1, lsize
      data(n) = 1.0_r8
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
    integer :: n
    real(kind=R8), dimension(:), allocatable :: ssh
    character(len=32), parameter :: sub = 'rof_import_mct'
    !---------------------------------------------------------------------------

    allocate(ssh(lsize))
    if (index_x2r_So_ssh>0) then
      do n = 1, lsize
         ssh(n) = x2r_r%rAttr(index_x2r_So_ssh,n)
      enddo
    end if

    if (masterproc) write(logunit_rof,*) 'Import ssh ',minval(ssh),maxval(ssh)
    deallocate(ssh)

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
        r2x_r%rattr(index_r2x_Forr_rofl,n)    = float(iam)
        r2x_r%rattr(index_r2x_Forr_rofi,n)    = float(iam)
        r2x_r%rattr(index_r2x_Flrr_flood,n)   = float(iam)
        r2x_r%rattr(index_r2x_Flrr_volr,n)    = float(iam)
        r2x_r%rattr(index_r2x_Flrr_volrmch,n) = float(iam)
        r2x_r%rattr(index_r2x_Flrr_supply,n)  = float(iam)
        r2x_r%rattr(index_r2x_Flrr_deficit,n) = float(iam)
      enddo

  end subroutine rof_export_mct

!====================================================================================

  subroutine rof_read_namelist()

    use shr_mpi_mod
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read the rdycore_in file and associated namelist
    !
    ! ARGUMENTS:
    implicit none
    !
    ! LOCAL VARIABLES
    character(len=256):: nlfilename_rof       ! namelist filename
    character(len=*),parameter :: subname = '(rof_read_namelist) '
    integer :: ier, unitn
    logical :: lexist

    !---------------------------------------------------------------------------

    namelist /rdycore_inparm / do_rdycore, rdycore_yaml_file, filename_rof

    ! default values
    do_rdycore        = .true.
    rdycore_yaml_file = ' '
    filename_rof      = ' '

    ! read namelist from expected file
    nlfilename_rof = "rdycore_in" // trim(inst_suffix)
    inquire (file = trim(nlfilename_rof), exist = lexist)
    if ( .not. lexist ) then
       write(logunit_rof,*) subname // ' ERROR: nlfilename_rof does NOT exist:'&
            //trim(nlfilename_rof)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_rof does not exist')
    end if
    if (masterproc) then
       unitn = shr_file_getunit()
       write(logunit_rof,*) 'Read in rdycore_inparm namelist from: ', trim(nlfilename_rof)
       open( unitn, file=trim(nlfilename_rof), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, rdycore_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on rdycore_inparm read' )
          endif
       end do
       call shr_file_freeunit(unitn)
    end if

    call mpi_bcast (do_rdycore,        1,                      MPI_LOGICAL,   0, mpicom_rof, ier)
    call mpi_bcast (rdycore_yaml_file, len(rdycore_yaml_file), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (filename_rof,      len(filename_rof),      MPI_CHARACTER, 0, mpicom_rof, ier)

    ! print out namelist settings to log
    if (masterproc) then
       write(logunit_rof,*) ' '
       write(logunit_rof,*) 'read from namelist:'
       write(logunit_rof,*) '   do_rdycore        = ', do_rdycore
       write(logunit_rof,*) '   rdycore_yaml_file = ', trim(rdycore_yaml_file)
       write(logunit_rof,*) '   filename_rof      = ', trim(filename_rof)
    end if

    return

  end subroutine rof_read_namelist

!====================================================================================

  subroutine  rof_read_rdycore()

    use shr_mpi_mod
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read the specififed mosart file
    !
    ! ARGUMENTS:
    implicit none
    !
    ! LOCAL VARIABLES
    integer :: i, j, count, ier
    logical :: found
    character(len=*),parameter :: subname = '(rof_read_rdycore) '
    integer, parameter :: RKIND = selected_real_kind(13)
    real(kind=RKIND), dimension(:),   allocatable :: lat1D, lon1D, area1D
    real(kind=RKIND), dimension(:,:), allocatable :: area
    type(file_desc_t) :: ncid                 ! netcdf file id

    !---------------------------------------------------------------------------

    ! open rdycore file
!JW    filename_rof = '/global/cfs/cdirs/e3sm/inputdata/rof/rdycore/MOSART_global_half_20180721a.nc'

    if (masterproc) then
       write(logunit_rof,*) 'Read in RDycore file name: ',trim(filename_rof)
       call shr_sys_flush(logunit_rof)
    endif

    call ncd_pio_openfile(ncid, trim(filename_rof), 0)

    call ncd_inqfdims(ncid, isgrid2d, nlong, nlatg, gsize)

    if (masterproc) then
       write(logunit_rof,*) 'Values for lon/lat: ', nlong, nlatg, gsize
       write(logunit_rof,*) 'Successfully read RDycore dimensions'
       if (isgrid2d) then
        write(logunit_rof,*) 'RDycore input is 2d'
       else
        write(logunit_rof,*) 'RDycore input is 1d'
       endif
       call shr_sys_flush(logunit_rof)
    endif

    ! allocate arrays
    if (isgrid2d) then

       allocate(lon1D(nlong))
       allocate(lat1D(nlatg))

       allocate(lonc_g(gsize))
       allocate(latc_g(gsize))
       allocate(areac_g(gsize))

       ! read the mesh data
       call ncd_io(ncid=ncid, varname='lon', flag='read', data=lon1D, readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore longitudes')
       if (masterproc) write(logunit_rof,*) 'Read lon ',minval(lon1D),maxval(lon1D)

       call ncd_io(ncid=ncid, varname='lat', flag='read', data=lat1D, readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore latitudes')
       if (masterproc) write(logunit_rof,*) 'Read lat ',minval(lat1D),maxval(lat1D)

       call ncd_io(ncid=ncid, varname='area', flag='read', data=areac_g, readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore area')
       if (masterproc) write(logunit_rof,*) 'Read area ',minval(areac_g),maxval(areac_g)

       count = 0
       do j = 1, nlong
          do i = 1, nlatg
             count = count + 1
             latc_g(count)  = lat1D(i)
             lonc_g(count)  = lon1D(j)
          end do
       end do

       ! free up memory
       deallocate(lat1D)
       deallocate(lon1D)

    else

       allocate(lonc_l(lsize))
       allocate(latc_l(lsize))
       allocate(areac_l(lsize))

       ! read the mesh data
       call ncd_io(ncid=ncid, varname='lon', flag='read', data=lonc_l, dim1name='gridcell', readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore longitudes')

       call ncd_io(ncid=ncid, varname='lat', flag='read', data=latc_l, dim1name='gridcell', readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore latitudes')

       call ncd_io(ncid=ncid, varname='area', flag='read', data=areac_l, dim1name='gridcell', readvar=found)
       if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read RDycore area')

    endif

    return

  end subroutine rof_read_rdycore
!====================================================================================

end module rof_comp_mct
