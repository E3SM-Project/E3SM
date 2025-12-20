module atm_comp_mct

  ! !USES:

  use esmf
  use netcdf
  use pio
  use mct_mod
  use perf_mod
  use seq_cdata_mod   , only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
  use seq_comm_mct    , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_strdata_mod , only: shr_strdata_type
  use shr_file_mod    , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only: shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only: shr_file_freeunit
  use seq_flds_mod    , only: seq_flds_a2x_fields, seq_flds_x2a_fields

  use atm_cpl_indices
  use eatmMod
  use eatm_comp_mod
  use eatmSpmdMod
  use eatmIO

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                :: inst_index          ! number of current instance (ie. 1)
  integer(IN)            :: compid              ! mct comp id
  integer(IN),parameter  :: master_task=0       ! task number of master task

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

  logical                :: do_eatm
  character(len=256)     :: filename_eatm

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine atm_init_mct( EClock, cdata, x2a, a2x, NLFilename )

    implicit none

    ! !DESCRIPTION: initialize data atm model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local variables ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    logical           :: atm_present               ! flag
    logical           :: atm_prognostic            ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    logical           :: exists                    ! true if file exists
    logical           :: first_time = .true.
    integer(IN)       :: ierr                      ! error code
    integer           :: mpicom_loc                ! local mpi communicator

    !--- formats ---
    character(*), parameter :: F00   = "('(atm_comp_init) ',8a)"
    integer(IN) , parameter :: master_task=0 ! task number of master task
    character(*), parameter :: subName = "(atm_init_mct) "
    !-------------------------------------------------------------------------------

    if (masterproc) write(logunit_atm,*) 'got to atm_init'
    ! Set cdata pointers to derived types (in coupler)
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom_loc, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata,&
         read_restart=read_restart)

    if (first_time) then

       ! Determine attribute vector indices
       call atm_cpl_indices_set()

       ! Initialize eatm MPI communicator
       call eatmSpmdInit(mpicom_loc)

       ! Determine instance information
       inst_name   = seq_comm_name(compid)
       inst_index  = seq_comm_inst(compid)
       inst_suffix = seq_comm_suffix(compid)

       ! Initialize pio
       call ncd_pio_init(inst_name)

       !--- open log file ---
       call shr_file_getLogUnit (shrlogunit)
       if (masterproc) then
          inquire(file='atm_modelio.nml'//trim(inst_suffix),exist=exists)
          if (exists) then
             logunit_atm = shr_file_getUnit()
             call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix),logunit_atm)
          end if
          write(logunit_atm,*) "eatm model initialization"
       else
          logunit_atm = shrlogunit
       endif

       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (logunit_atm)

       !----------------------------------------------------------------------------
       ! Read input namelists and set present and prognostic flags
       !----------------------------------------------------------------------------

       call t_startf('eatm_readnml')
       atm_present=.true.
       atm_prognostic=.true.
       call eatm_read_namelist()

       call seq_infodata_PutData(infodata, &
            atm_present=atm_present, &
            atm_prognostic=atm_prognostic)
       call t_stopf('eatm_readnml')

       !----------------------------------------------------------------------------
       ! Initialize eatm
       !----------------------------------------------------------------------------

       if (masterproc) write(logunit_atm,*) 'got to atm_init 2'
       ! Initialize atm gsMap
       call atm_SetGSMap_mct( mpicom_atm, compid, gsMap)

       ! Initialize atm domain
       if (masterproc) write(logunit_atm,*) 'got to atm_init 3'
       call atm_domain_mct( gsMap, ggrid )

       ! Initialize cpl -> eatm attribute vector
       call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
       call mct_aVect_zero(x2a)

       ! Initialize eatm -> cpl attribute vector
       call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)
       call mct_aVect_zero(a2x)

       if (masterproc) write(logunit_atm,*) 'got to atm_init 4'
       call eatm_comp_init(Eclock, x2a, a2x, &
            seq_flds_x2a_fields, seq_flds_a2x_fields, &
            gsmap, ggrid, read_restart)

       !----------------------------------------------------------------------------
       ! Fill infodata that needs to be returned from eatm
       !----------------------------------------------------------------------------

       if (masterproc) write(logunit_atm,*) 'got to atm_init 5'
       call seq_infodata_PutData(infodata, &
            atm_nx=lsize_x, &
            atm_ny=lsize_y)
  
       !----------------------------------------------------------------------------
       ! Create initial atm export state
       !----------------------------------------------------------------------------
       call atm_export_mct(a2x)

       !----------------------------------------------------------------------------
       ! Reset shr logging to original values
       !----------------------------------------------------------------------------

       if (masterproc) write(logunit_atm,F00) 'eatm_comp_init done'
       call shr_sys_flush(logunit_atm)
       call shr_sys_flush(shrlogunit)

       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

       first_time = .false.

    else ! so here first_time == .false.

       ! Redirect share output to eatm log

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (logunit_atm)

       call atm_import_mct(x2a)

       call t_startf('EATM_run')
       call eatm_comp_run( &
         EClock = EClock, &
         x2a = x2a, &
         a2x = a2x, &
         gsmap = gsmap, &
         ggrid = ggrid)
       call t_stopf('EATM_run1')

       call atm_export_mct(a2x)

       ! End redirection of share output to eatm log

       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

    end if

    call shr_sys_flush(logunit_atm)

    return

  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct( EClock, cdata,  x2a, a2x)

    ! !DESCRIPTION: run method for eatm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2a
    type(mct_aVect)             ,intent(inout) :: a2x

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    character(*), parameter :: subName = "(atm_run_mct) "
    !-------------------------------------------------------------------------------
    if (masterproc) write(logunit_atm,*) 'got to atm_run'
    call shr_sys_flush(logunit_atm)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call seq_infodata_GetData(infodata, &
         case_name=case_name)

    ! Map MCT datatype to eatm data
    call t_startf ('lc_eatm_import')
    call atm_import_mct( x2a )
    call t_stopf ('lc_eatm_import')

    call eatm_comp_run( &
         EClock = EClock, &
         x2a = x2a, &
         a2x = a2x, &
         gsmap = gsmap, &
         ggrid = ggrid)

    ! Map eatm data to MCT datatype
    call t_startf ('lc_eatm_export')
    call atm_export_mct( a2x )
    call t_stopf ('lc_eatm_export')

  end subroutine atm_run_mct

  !===============================================================================

  subroutine atm_final_mct(EClock, cdata, x2a, a2x)

    ! !DESCRIPTION: finalize method for dead atm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2a
    type(mct_aVect)             ,intent(inout) :: a2x

    !--- formats ---
    character(*), parameter :: subName = "(atm_final_mct) "
    !-------------------------------------------------------------------------------

    call  eatm_comp_final()

  end subroutine atm_final_mct

  !====================================================================================

  subroutine atm_SetGSMap_mct( mpicom_atm, compid, gsMap)

    ! DESCRIPTION: Set the MCT GS map for the atm model
    !
    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: mpicom_atm    ! MPI communicator
    integer        , intent(in)    :: compid        ! eatm model identifier
    type(mct_gsMap), intent(inout) :: gsMap         ! MCT gsmap
    !
    ! LOCAL VARIABLES
    integer :: n                         ! indices
    integer :: ier                       ! error code
    character(len=32), parameter :: sub = 'atm_SetGSMap_mct'
    !-----------------------------------------------------

    call atm_read_eatm()

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

     call mct_gsMap_init( gsMap, compid, npes, gsize, start, length, pe_loc)
    ! Init the gsMap with gindex
    ! call mct_gsMap_init( gsMap, gindex, mpicom_atm, compid, lsize, gsize )

    return

  end subroutine atm_SetGSMap_mct

  !====================================================================================

  subroutine atm_domain_mct( gsMap, dom_atm )

    ! !DESCRIPTION: This routine sets up the MCT domain
 
    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    type(mct_gsMap), intent(in)    :: gsMap
    type(mct_gGrid), intent(inout) :: dom_atm
      !
    ! LOCAL VARIABLES
    integer :: n, ni              ! index
    integer , pointer :: idata(:) ! temporary
    real(r8), pointer :: data(:)  ! temporary
    character(len=32), parameter :: sub = 'atm_domain_mct'
    !-----------------------------------------------------

    ! lat/lon in degrees,  area in radians^2
    call mct_gGrid_init( GGrid=dom_atm, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Allocate memory
    allocate(data(lsize))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap, iam, idata)
    call mct_gGrid_importIAttr(dom_atm,'GlobGridNum',idata,lsize)

    ! Initialize attribute vector with special value
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_atm,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_atm,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_atm,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_atm,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_atm,"mask" ,data,lsize)

    data(:) = lonc_g(:)
    call mct_gGrid_importRattr(dom_atm,"lon",data,lsize)
    data(:) = latc_g(:)
    call mct_gGrid_importRattr(dom_atm,"lat",data,lsize)

    data(:) = areac_g(:)
    call mct_gGrid_importRattr(dom_atm,"area",data,lsize)

    do n = 1, lsize
      data(n) = 1.0_r8
    end do
    call mct_gGrid_importRattr(dom_atm,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_atm,"frac",data,lsize)

    deallocate(start,length,pe_loc)
    deallocate(data)
    deallocate(idata)

  end subroutine atm_domain_mct

  !====================================================================================

  subroutine atm_import_mct( x2a_a)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Obtain the atm input from the coupler
    ! convert from kg/m2s to m3/s
    !
    ! ARGUMENTS:
    implicit none
    type(mct_aVect), intent(inout) :: x2a_a
    !
    ! LOCAL VARIABLES
    integer :: i, j, n
    character(len=32), parameter :: sub = 'atm_import_mct'
    !---------------------------------------------------------------------------

    if (masterproc) write(logunit_atm,*) 'got to atm_import'
    n = 0
    do j = 1, lsize_y
       do i = 1, lsize_x
          n = n + 1
          if (index_x2a_Faxx_sen>0) then
             shf(i,j) = -x2a_a%rAttr(index_x2a_Faxx_sen,n)
          end if
          if (masterproc) write(logunit_atm,*) 'Import shf ',minval(shf),maxval(shf)
          if (index_x2a_Faxx_evap>0) then
             cflx(i,j) = -x2a_a%rAttr(index_x2a_Faxx_evap,n)
          end if
          if (index_x2a_Faxx_lat>0) then
             lhf(i,j) = x2a_a%rAttr(index_x2a_Faxx_lat,n)
          end if
          if (index_x2a_Faxx_taux>0) then
             wsx(i,j) = -x2a_a%rAttr(index_x2a_Faxx_taux,n)
          end if
          if (index_x2a_Faxx_tauy>0) then
             wsy(i,j) = -x2a_a%rAttr(index_x2a_Faxx_tauy,n)
          end if
          if (index_x2a_Faxx_lwup>0) then
             lwup(i,j) = -x2a_a%rAttr(index_x2a_Faxx_lwup,n)
          end if
          if (index_x2a_Sx_avsdr>0) then
             asdir(i,j) = x2a_a%rAttr(index_x2a_Sx_avsdr,n)
          end if
          if (index_x2a_Sx_anidr>0) then
             aldir(i,j) = x2a_a%rAttr(index_x2a_Sx_anidr,n)
          end if
          if (index_x2a_Sx_avsdf>0) then
             asdif(i,j) = x2a_a%rAttr(index_x2a_Sx_avsdf,n)
          end if
          if (index_x2a_Sx_anidf>0) then
             aldif(i,j) = x2a_a%rAttr(index_x2a_Sx_anidf,n)
          end if
          if (index_x2a_Sx_t>0) then
             ts(i,j) = x2a_a%rAttr(index_x2a_Sx_t,n)
          end if
          if (index_x2a_So_t>0) then
             sst(i,j) = x2a_a%rAttr(index_x2a_So_t,n)
          end if
          if (index_x2a_Sl_snowh>0) then
             snowhland(i,j) = x2a_a%rAttr(index_x2a_Sl_snowh,n)
          end if
          if (index_x2a_Si_snowh>0) then
             snowhice(i,j) = x2a_a%rAttr(index_x2a_Si_snowh,n)
          end if
          if (index_x2a_Sx_tref>0) then
             tref(i,j) = x2a_a%rAttr(index_x2a_Sx_tref,n)
          end if
          if (index_x2a_Sx_qref>0) then
             qref(i,j) = x2a_a%rAttr(index_x2a_Sx_qref,n)
          end if
          if (index_x2a_Sx_u10>0) then
             u10(i,j) = x2a_a%rAttr(index_x2a_Sx_u10,n)
          end if
          if (index_x2a_Sx_u10withgusts>0) then
             u10withgusts(i,j) = x2a_a%rAttr(index_x2a_Sx_u10withgusts,n)
          end if
          if (index_x2a_Sf_ifrac>0) then
             icefrac(i,j) = x2a_a%rAttr(index_x2a_Sf_ifrac,n)
          end if
          if (index_x2a_Sf_ofrac>0) then
             ocnfrac(i,j) = x2a_a%rAttr(index_x2a_Sf_ofrac,n)
          end if
          if (index_x2a_Sf_lfrac>0) then
             lndfrac(i,j) = x2a_a%rAttr(index_x2a_Sf_lfrac,n)
          end if
       end do
    end do

  end subroutine atm_import_mct

!====================================================================================

  subroutine atm_export_mct( a2x_a )

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Send the atm model export state to the coupler
    ! convert from m3/s to kg/m2s
    !
    ! ARGUMENTS:
    implicit none
    type(mct_aVect), intent(inout) :: a2x_a  ! atm to coupler export state
    !
    ! LOCAL VARIABLES
    integer :: i, j, n
    logical,save :: first_time = .true.
    character(len=32), parameter :: sub = 'atm_export_mct'
    !---------------------------------------------------------------------------

    !JW nothing to return yet and sending 0's will wreak havoc, so just stop for now
    if (masterproc) write(logunit_atm,*) 'got to atm_export, stopping'
    call shr_sys_flush(logunit_atm)
    stop

    !JW list returned could be different for different emulators
    n = 0
    do j = 1, lsize_y
       do i = 1, lsize_x
          n = n + 1
          a2x_a%rattr(index_a2x_Sa_z,n)       = zbot(i,j)
!JW          a2x_a%rattr(index_a2x_Sa_topo,n)    = 0._r8
          a2x_a%rattr(index_a2x_Sa_u,n)       = ubot(i,j)
          a2x_a%rattr(index_a2x_Sa_v,n)       = vbot(i,j)
          a2x_a%rattr(index_a2x_Sa_tbot,n)    = tbot(i,j)
          a2x_a%rattr(index_a2x_Sa_ptem,n)    = ptem(i,j)
          a2x_a%rattr(index_a2x_Sa_shum,n)    = shum(i,j)
          a2x_a%rattr(index_a2x_Sa_dens,n)    = dens(i,j)
          a2x_a%rattr(index_a2x_Sa_pbot,n)    = pbot(i,j)
          a2x_a%rattr(index_a2x_Sa_pslv,n)    = pslv(i,j)
          a2x_a%rattr(index_a2x_Faxa_lwdn,n)  = lwdn(i,j)
          a2x_a%rattr(index_a2x_Faxa_rainc,n) = rainc(i,j)
          a2x_a%rattr(index_a2x_Faxa_rainl,n) = rainl(i,j)
          a2x_a%rattr(index_a2x_Faxa_snowc,n) = snowc(i,j)
          a2x_a%rattr(index_a2x_Faxa_snowl,n) = snowl(i,j)
          a2x_a%rattr(index_a2x_Faxa_swndr,n) = swndr(i,j)
          a2x_a%rattr(index_a2x_Faxa_swvdr,n) = swvdr(i,j)
          a2x_a%rattr(index_a2x_Faxa_swndf,n) = swndf(i,j)
          a2x_a%rattr(index_a2x_Faxa_swvdf,n) = swvdf(i,j)
          a2x_a%rattr(index_a2x_Faxa_swnet,n) = swnet(i,j)
       enddo
    enddo

  end subroutine atm_export_mct

  !====================================================================================

  subroutine eatm_read_namelist()

    use shr_mpi_mod
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read the eatm_in file and associated namelist
    !
    ! ARGUMENTS:
    implicit none
    !
    ! LOCAL VARIABLES
    character(len=256):: nlfilename_atm       ! namelist filename
    character(len=*),parameter :: subname = '(atm_read_namelist) '
    integer :: ier, unitn
    logical :: lexist

    !---------------------------------------------------------------------------

    namelist /eatm_inparm / do_eatm, filename_eatm

    ! default values
    do_eatm        = .true.
    filename_eatm  = ' '

    ! read namelist from expected file
    nlfilename_atm = "eatm_in" // trim(inst_suffix)
    inquire (file = trim(nlfilename_atm), exist = lexist)
    if ( .not. lexist ) then
       write(logunit_atm,*) subname // ' ERROR: nlfilename_atm does NOT exist:'&
            //trim(nlfilename_atm)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_atm does not exist')
    end if
    if (masterproc) then
       unitn = shr_file_getunit()
       write(logunit_atm,*) 'Read in eatm_inparm namelist from: ', trim(nlfilename_atm)
       open( unitn, file=trim(nlfilename_atm), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, eatm_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on eatm_inparm read' )
          endif
       end do
       call shr_file_freeunit(unitn)
    end if
    call mpi_bcast (do_eatm,        1,                      MPI_LOGICAL,   0, mpicom_atm, ier)
    call mpi_bcast (filename_eatm,  len(filename_eatm),     MPI_CHARACTER, 0, mpicom_atm, ier)

    ! print out namelist settings to log
    if (masterproc) then
       write(logunit_atm,*) ' '
       write(logunit_atm,*) 'read from namelist:'
       write(logunit_atm,*) '   do_eatm           = ', do_eatm
       write(logunit_atm,*) '   filename_eatm     = ', trim(filename_eatm)
    end if

    return

  end subroutine eatm_read_namelist

!====================================================================================

  subroutine atm_read_eatm()

    use shr_mpi_mod
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read the specififed eatm file
    !
    ! ARGUMENTS:
    implicit none
    !
    ! LOCAL VARIABLES
    logical :: found
    integer, dimension(:) :: grid_dims(2)
    character(len=*),parameter :: subname = '(atm_read_eatm) '
    integer, parameter :: RKIND = selected_real_kind(13)
    type(file_desc_t) :: ncid                 ! netcdf file id

    !---------------------------------------------------------------------------

    ! open eatm file
    if (masterproc) then
       write(logunit_atm,*) 'Read in eatm file name: ',trim(filename_eatm)
       call shr_sys_flush(logunit_atm)
    endif

    call ncd_pio_openfile(ncid, trim(filename_eatm), 0)

    call ncd_inqfdims(ncid, gsize)

    call ncd_io(varname='grid_dims', data=grid_dims, flag='read', ncid=ncid, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: reading EATM grid_dims')

    !JW for now let lsize = gsize (local vs global)
    lsize = gsize
    lsize_x = grid_dims(1)
    lsize_y = grid_dims(2)

    if (masterproc) then
       write(logunit_atm,*) 'Values for lon/lat: ', lsize_x, lsize_y
       write(logunit_atm,*) 'Successfully read eatm dimensions'
       call shr_sys_flush(logunit_atm)
    endif

    allocate(lonc_g(gsize))
    allocate(latc_g(gsize))
    allocate(areac_g(gsize))

    ! read the mesh data
    call ncd_io(ncid=ncid, varname='grid_center_lon', flag='read', data=lonc_g, dim1name='grid_size', readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read eatm longitudes')

    call ncd_io(ncid=ncid, varname='grid_center_lat', flag='read', data=latc_g, dim1name='grid_size', readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read eatm latitudes')

    call ncd_io(ncid=ncid, varname='grid_area', flag='read', data=areac_g, dim1name='grid_size', readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read eatm area')

    return

  end subroutine atm_read_eatm

  !===============================================================================

end module atm_comp_mct
