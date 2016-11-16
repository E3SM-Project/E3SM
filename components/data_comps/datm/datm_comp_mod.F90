#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module datm_comp_mod

! !USES:

  use shr_const_mod
  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                               CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use shr_cal_mod      , only: shr_cal_date2julian
  use shr_mpi_mod      , only: shr_mpi_bcast
  use mct_mod
  use esmf
  use perf_mod

  use shr_strdata_mod
  use shr_dmodel_mod
  use shr_pcdf_mod
  use datm_shr_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: seq_flds_a2x_fields, &
                               seq_flds_x2a_fields
  use shr_precip_mod   , only: shr_precip_partition_rain_snow_ramp
!
! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: datm_comp_init
  public :: datm_comp_run
  public :: datm_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  !--- other ---
  character(CS) :: myModelName = 'atm'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: COMPID                ! mct comp id
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance
                                         ! (ie. "_0001" or "")
  character(CL) :: atm_mode              ! mode
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: firstcall = .true.    ! first call logical
  logical       :: scmMode = .false.     ! single column mode
  real(R8)      :: scmLat  = shr_const_SPVAL  ! single column lat
  real(R8)      :: scmLon  = shr_const_SPVAL  ! single column lon
  integer       :: phase                 ! phase of method
  logical       :: read_restart          ! start from restart
  real(R8)      :: orbEccen              ! orb eccentricity (unit-less)
  real(R8)      :: orbMvelpp             ! orb moving vernal eq (radians)
  real(R8)      :: orbLambm0             ! orb mean long of perhelion (radians)
  real(R8)      :: orbObliqr             ! orb obliquity (radians)
  real(R8)      :: tbotmax               ! units detector
  real(R8)      :: tdewmax               ! units detector
  real(R8)      :: anidrmax              ! existance detector
  integer(IN)   :: iradsw                ! radiation logical
  character(CL) :: factorFn              ! file containing correction factors

  character(len=*),parameter :: rpfile = 'rpointer.atm'
  character(len=*),parameter :: nullstr = 'undefined'

  real(R8),parameter :: aerodep_spval = 1.e29_r8    ! special aerosol deposition
  real(R8),parameter :: tKFrz  = SHR_CONST_TKFRZ
  real(R8),parameter :: degtorad = SHR_CONST_PI/180.0_R8
  real(R8),parameter :: pstd   = SHR_CONST_PSTD     ! standard pressure ~ Pa
  real(R8),parameter :: stebol = SHR_CONST_STEBOL   ! Stefan-Boltzmann constant ~ W/m^2/K^4
  real(R8),parameter :: rdair  = SHR_CONST_RDAIR    ! dry air gas constant   ~ J/K/kg
  real(R8),parameter :: avg_c0 =  61.846_R8
  real(R8),parameter :: avg_c1 =   1.107_R8
  real(R8),parameter :: amp_c0 = -21.841_R8
  real(R8),parameter :: amp_c1 =  -0.447_R8
  real(R8),parameter :: phs_c0 =   0.298_R8
  real(R8),parameter :: dLWarc =  -5.000_R8
  real(R8)     ,save :: dTarc(12)
  data   dTarc      / 0.49_R8, 0.06_R8,-0.73_R8,  -0.89_R8,-0.77_R8,-1.02_R8, &
  &                  -1.99_R8,-0.91_R8, 1.72_R8,   2.30_R8, 1.81_R8, 1.06_R8/

  integer(IN) :: kz,ktopo,ku,kv,ktbot,kptem,kshum,kdens,kpbot,kpslv,klwdn
  integer(IN) :: krc,krl,ksc,ksl,kswndr,kswndf,kswvdr,kswvdf,kswnet,kco2p,kco2d
  integer(IN) :: kbid,kbod,kbiw,koid,kood,koiw,kdw1,kdw2,kdw3,kdw4,kdd1,kdd2,kdd3,kdd4
  integer(IN) :: kanidr,kanidf,kavsdr,kavsdf
  integer(IN) :: stbot,swind,sz,spbot,sshum,stdew,srh,slwdn,sswdn,sswdndf,sswdndr
  integer(IN) :: sprecc,sprecl,sprecn,sco2p,sco2d,sswup,sprec,starcf
  !
  ! anomaly forcing
  !
  integer(IN) :: sprecsf
  integer(IN) :: sprec_af,su_af,sv_af,stbot_af,sshum_af,spbot_af,slwdn_af,sswdn_af

  type(shr_strdata_type) :: SDATM
  type(mct_rearr) :: rearr
  type(mct_avect) :: avstrm   ! av of data from stream
  integer(IN), pointer :: imask(:)
  real(R8), pointer :: yc(:)
  real(R8), pointer :: windFactor(:)
  real(R8), pointer :: winddFactor(:)
  real(R8), pointer :: qsatFactor(:)
  !
  ! for anomaly forcing
  !
  integer(IN),parameter :: ktrans  = 66

  character(16),parameter  :: avofld(1:ktrans) = &
       (/"Sa_z            ","Sa_topo         ", &
         "Sa_u            ","Sa_v            ","Sa_tbot         ", &
         "Sa_ptem         ","Sa_shum         ","Sa_dens         ","Sa_pbot         ", &
         "Sa_pslv         ","Faxa_lwdn       ","Faxa_rainc      ","Faxa_rainl      ", &
         "Faxa_snowc      ","Faxa_snowl      ","Faxa_swndr      ","Faxa_swvdr      ", &
         "Faxa_swndf      ","Faxa_swvdf      ","Faxa_swnet      ","Sa_co2prog      ", &
         "Sa_co2diag      ","Faxa_bcphidry   ","Faxa_bcphodry   ","Faxa_bcphiwet   ", &
         "Faxa_ocphidry   ","Faxa_ocphodry   ","Faxa_ocphiwet   ","Faxa_dstwet1    ", &
         "Faxa_dstwet2    ","Faxa_dstwet3    ","Faxa_dstwet4    ","Faxa_dstdry1    ", &
         "Faxa_dstdry2    ","Faxa_dstdry3    ","Faxa_dstdry4    ",                    &
         "Sx_tref         ","Sx_qref         ","Sx_avsdr        ","Sx_anidr        ", &
         "Sx_avsdf        ","Sx_anidf        ","Sx_t            ","So_t            ", &
         "Sl_snowh        ","Sf_lfrac        ","Sf_ifrac        ","Sf_ofrac        ", &
         "Faxx_taux       ","Faxx_tauy       ","Faxx_lat        ","Faxx_sen        ", &
         "Faxx_lwup       ","Faxx_evap       ","Fall_fco2_lnd   ","Faoo_fco2_ocn   ", &
         "Faoo_fdms_ocn   ",  &
         !
         ! add values for bias correction / anomaly forcing
         !
         "Sa_precsf       ", &
         "Sa_prec_af      ","Sa_u_af         ","Sa_v_af         ","Sa_tbot_af      ",&
         "Sa_pbot_af      ","Sa_shum_af      ","Sa_swdn_af      ","Sa_lwdn_af      " &
       /)

  character(16),parameter  :: avifld(1:ktrans) = &
       (/"z               ","topo            ", &
         "u               ","v               ","tbot            ", &
         "ptem            ","shum            ","dens            ","pbot            ", &
         "pslv            ","lwdn            ","rainc           ","rainl           ", &
         "snowc           ","snowl           ","swndr           ","swvdr           ", &
         "swndf           ","swvdf           ","swnet           ","co2prog         ", &
         "co2diag         ","bcphidry        ","bcphodry        ","bcphiwet        ", &
         "ocphidry        ","ocphodry        ","ocphiwet        ","dstwet1         ", &
         "dstwet2         ","dstwet3         ","dstwet4         ","dstdry1         ", &
         "dstdry2         ","dstdry3         ","dstdry4         ",                    &
         "tref            ","qref            ","avsdr           ","anidr           ", &
         "avsdf           ","anidf           ","ts              ","to              ", &
         "snowhl          ","lfrac           ","ifrac           ","ofrac           ", &
         "taux            ","tauy            ","lat             ","sen             ", &
         "lwup            ","evap            ","co2lnd          ","co2ocn          ", &
         ! add precsf
         "dms             ","precsf          ", &
         ! add Sa_precsf for precip scale factor
         "prec_af         ","u_af            ","v_af            ","tbot_af         ", &
         "pbot_af         ","shum_af         ","swdn_af         ","lwdn_af         "  &
       /)

  ! add stream for anomaly forcing
  integer(IN),parameter :: ktranss = 28

  ! The stofld and stifld lists are used for fields that are read but not passed to the
  ! coupler (e.g., they are used to compute fields that are passed to the coupler), and
  ! other fields used in calculations. Fields that are simply read and passed directly to
  ! the coupler do not need to be in these lists.
  character(16),parameter  :: stofld(1:ktranss) = &
       (/"strm_tbot       ","strm_wind       ","strm_z          ","strm_pbot       ", &
         "strm_shum       ","strm_tdew       ","strm_rh         ","strm_lwdn       ", &
         "strm_swdn       ","strm_swdndf     ","strm_swdndr     ","strm_precc      ", &
         "strm_precl      ","strm_precn      ","strm_co2prog    ","strm_co2diag    ", &
         "strm_swup       ","strm_prec       ","strm_tarcf      ", &
         ! add bias correction / anomaly forcing streams
         "strm_precsf     ", &
         "strm_prec_af    ","strm_u_af       ","strm_v_af       ","strm_tbot_af    ", &
         "strm_pbot_af    ","strm_shum_af    ","strm_swdn_af    ","strm_lwdn_af    "  &
       /)

  character(16),parameter  :: stifld(1:ktranss) = &
       (/"tbot            ","wind            ","z               ","pbot            ", &
         "shum            ","tdew            ","rh              ","lwdn            ", &
         "swdn            ","swdndf          ","swdndr          ","precc           ", &
         "precl           ","precn           ","co2prog         ","co2diag         ", &
         ! add precsf
         "swup            ","prec            ","tarcf           ","precsf          ", &
         ! add anomaly forcing streams
         "prec_af         ","u_af            ","v_af            ","tbot_af         ", &
         "pbot_af         ","shum_af         ","swdn_af         ","lwdn_af         "  &
       /)

  character(CL), pointer :: ilist_av(:)     ! input list for translation
  character(CL), pointer :: olist_av(:)     ! output list for translation
  character(CL), pointer :: ilist_st(:)     ! input list for translation
  character(CL), pointer :: olist_st(:)     ! output list for translation
  integer(IN)  , pointer :: count_av(:)
  integer(IN)  , pointer :: count_st(:)

  save

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: datm_comp_init
!
! !DESCRIPTION:
!     initialize data atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine datm_comp_init( EClock, cdata, x2a, a2x, NLFilename )

    use pio         , only : iosystem_desc_t
    use shr_pio_mod , only : shr_pio_getiosys, shr_pio_getiotype
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

   !EOP

    !--- local variables ---
    integer(IN)   :: n,k                   ! generic counters
    integer(IN)   :: ierr                  ! error code
    integer(IN)   :: gsize                 ! global size
    integer(IN)   :: lsize                 ! local size
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: nunit                 ! unit number
    integer(IN)   :: kmask                 ! field reference
    integer(IN)   :: klat                  ! field reference
    integer(IN)   :: kfld                  ! fld index
    integer(IN)   :: cnt                   ! counter
    logical       :: atm_present           ! flag
    logical       :: atm_prognostic        ! flag

    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    type(iosystem_desc_t)  , pointer :: iosystem

    character(CL) :: filePath              ! generic file path
    character(CL) :: fileName              ! generic file name
    character(CS) :: timeName              ! domain file: time variable name
    character(CS) ::  lonName              ! domain file: lon  variable name
    character(CS) ::  latName              ! domain file: lat  variable name
    character(CS) :: maskName              ! domain file: mask variable name
    character(CS) :: areaName              ! domain file: area variable name

    integer(IN)   :: yearFirst             ! first year to use in data stream
    integer(IN)   :: yearLast              ! last  year to use in data stream
    integer(IN)   :: yearAlign             ! data year that aligns with yearFirst

    character(CL) :: atm_in                ! dshr atm namelist
    character(CL) :: decomp                ! decomp strategy
    character(CL) :: rest_file             ! restart filename
    character(CL) :: rest_file_strm        ! restart filename for streams
    character(CL) :: restfilm              ! model restart file namelist
    character(CL) :: restfils              ! stream restart file namelist
    logical       :: exists                ! filename existance
    integer(IN)   :: nu                    ! unit number
    integer(IN)   :: idt                   ! integer timestep
    integer(IN)   :: CurrentYMD            ! model date
    integer(IN)   :: CurrentTOD            ! model sec into model date
    integer(IN)   :: stepno                ! step number
    real(R8)      :: nextsw_cday           ! calendar of next atm sw
    character(CL) :: flds_strm
    logical       :: presaero              ! true => send valid prescribe aero fields to coupler
    logical       :: force_prognostic_true ! if true set prognostic true
    character(CL) :: calendar              ! calendar type
    character(CL) :: bias_correct          ! true => send bias correction fields to coupler
    character(CL) :: anomaly_forcing(8)    ! true => send anomaly forcing fields to coupler

    !----- define namelist -----
    namelist / datm_nml / &
        atm_in, decomp, iradsw, factorFn, restfilm, restfils, presaero, bias_correct, &
        anomaly_forcing, force_prognostic_true

    !--- formats ---
    character(*), parameter :: F00   = "('(datm_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(datm_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(datm_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(datm_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(datm_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(datm_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(datm_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(datm_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(datm_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(datm_comp_init) "
!-------------------------------------------------------------------------------


    call t_startf('DATM_INIT')

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsmap, dom=ggrid, infodata=infodata)
    call seq_infodata_getData(infodata,atm_phase=phase)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    if (phase == 1) then
       ! Determine communicator groups and sizes
       call mpi_comm_rank(mpicom, my_task, ierr)
       call mpi_comm_size(mpicom, npes, ierr)
       
       !--- open log file ---
       if (my_task == master_task) then
          logUnit = shr_file_getUnit()
          call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix),logUnit)
       else
          logUnit = 6
       endif
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    if (phase == 1) then

       !----------------------------------------------------------------------------
       ! Set a Few Defaults
       !----------------------------------------------------------------------------

       call seq_infodata_getData(infodata,single_column=scmMode, &
            scmlat=scmlat, scmlon=scmLon)
       call seq_infodata_GetData(infodata,orb_eccen=orbEccen,orb_mvelpp=orbMvelpp, &
            orb_lambm0=orbLambm0,orb_obliqr=orbObliqr )

       atm_present = .false.
       atm_prognostic = .false.
       call seq_infodata_GetData(infodata,read_restart=read_restart)

       !----------------------------------------------------------------------------
       ! Read datm_in
       !----------------------------------------------------------------------------

       call t_startf('datm_readnml')

       filename = "datm_in"//trim(inst_suffix)
       atm_in   = "unset"
       decomp   = "1d"
       iradsw   = 0
       factorFn = 'null'
       restfilm = trim(nullstr)
       restfils = trim(nullstr)
       presaero = .false.
       force_prognostic_true = .false.

       if (my_task == master_task) then
          nunit = shr_file_getUnit() ! get unused unit number
          open (nunit,file=trim(filename),status="old",action="read")
          read (nunit,nml=datm_nml,iostat=ierr)
          close(nunit)

          call shr_file_freeUnit(nunit)
          if (ierr > 0) then
             write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
             call shr_sys_abort(subName//': namelist read error '//trim(filename))
          end if
          write(logunit,F00)' atm_in   = ',trim(atm_in)
          write(logunit,F00)' decomp   = ',trim(decomp)
          write(logunit,F01)' iradsw   = ',iradsw
          write(logunit,F00)' factorFn = ',trim(factorFn)
          write(logunit,F00)' restfilm = ',trim(restfilm)
          write(logunit,F00)' restfils = ',trim(restfils)
          write(logunit,F0L)' presaero = ',presaero
          write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
          write(logunit,F01) 'inst_index  =  ',inst_index
          write(logunit,F00) 'inst_name   =  ',trim(inst_name)
          write(logunit,F00) 'inst_suffix =  ',trim(inst_suffix)
          call shr_sys_flush(logunit)
       endif

       call shr_mpi_bcast(atm_in,mpicom,'atm_in')
       call shr_mpi_bcast(decomp,mpicom,'decomp')
       call shr_mpi_bcast(iradsw,mpicom,'iradsw')
       call shr_mpi_bcast(factorFn,mpicom,'factorFn')
       call shr_mpi_bcast(restfilm,mpicom,'restfilm')
       call shr_mpi_bcast(restfils,mpicom,'restfils')
       call shr_mpi_bcast(presaero,mpicom,'presaero')
       call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')

       rest_file = trim(restfilm)
       rest_file_strm = trim(restfils)
       if (force_prognostic_true) then
          atm_present    = .true.
          atm_prognostic = .true.
       endif

       !----------------------------------------------------------------------------
       ! Read dshr namelist
       !----------------------------------------------------------------------------

       call shr_strdata_readnml(SDATM,trim(atm_in),mpicom=mpicom)
       call shr_sys_flush(shrlogunit)

       !----------------------------------------------------------------------------
       ! Initialize PIO
       !----------------------------------------------------------------------------

       iosystem => shr_pio_getiosys(trim(inst_name))
       call shr_strdata_pioinit(SDATM, iosystem, &
            shr_pio_getiotype(trim(inst_name)))

       !----------------------------------------------------------------------------
       ! Validate mode
       !----------------------------------------------------------------------------

       atm_mode = trim(SDATM%dataMode)

       ! check that we know how to handle the mode

       if (trim(atm_mode) == 'NULL'      .or. &
            trim(atm_mode) == 'CORE2_NYF' .or. &
            trim(atm_mode) == 'CORE2_IAF' .or. &
            trim(atm_mode) == 'WRF'       .or. &
            trim(atm_mode) == 'CLMNCEP'   .or. &
            trim(atm_mode) == 'CPLHIST'   .or. &
            trim(atm_mode) == 'COPYALL'   ) then
          if (my_task == master_task) then
             write(logunit,F00) ' atm mode = ',trim(atm_mode)
             call shr_sys_flush(logunit)
          end if
       else
          write(logunit,F00) ' ERROR illegal atm mode = ',trim(atm_mode)
          call shr_sys_abort()
       endif

       call t_stopf('datm_readnml')

       !----------------------------------------------------------------------------
       ! Initialize datasets
       !----------------------------------------------------------------------------

       call t_startf('datm_strdata_init')

       if (trim(atm_mode) /= 'NULL') then
          atm_present = .true.
          call seq_timemgr_EClockGetData( EClock, dtime=idt, calendar=calendar )
          if (scmmode) then
             if (my_task == master_task) &
                  write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
             call shr_strdata_init(SDATM, mpicom, compid,name='atm', &
                  scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, &
                  calendar=calendar)
          else
             call shr_strdata_init(SDATM, mpicom, compid, name='atm', &
                  calendar=calendar)
          endif

          if (my_task == master_task) call shr_sys_flush(shrlogunit)

          !--- overwrite mask and frac ---
          k = mct_aVect_indexRA(SDATM%grid%data,'mask')
          SDATM%grid%data%rAttr(k,:) = 1.0_R8
          k = mct_aVect_indexRA(SDATM%grid%data,'frac')
          SDATM%grid%data%rAttr(k,:) = 1.0_R8

          !--- set data needed for cosz t-interp method ---

          call shr_strdata_setOrbs(SDATM,orbEccen,orbMvelpp,orbLambm0,orbObliqr,idt)
       endif

       if (my_task == master_task) then
          call shr_strdata_print(SDATM,'ATM data')
          call shr_sys_flush(shrlogunit)
       endif

       call t_stopf('datm_strdata_init')

       !----------------------------------------------------------------------------
       ! Set flag to specify data components
       !----------------------------------------------------------------------------

       call seq_infodata_PutData(infodata, &
            atm_present=atm_present, atm_prognostic=atm_prognostic, &
            atm_nx=SDATM%nxg, atm_ny=SDATM%nyg )
       call seq_infodata_PutData( infodata, atm_aero=presaero)

       !----------------------------------------------------------------------------
       ! Initialize MCT global seg map, 1d decomp
       !----------------------------------------------------------------------------

       call t_startf('datm_initgsmaps')
       if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
       call shr_sys_flush(logunit)

       call shr_dmodel_gsmapcreate(gsmap,SDATM%nxg*SDATM%nyg,compid,mpicom,decomp)
       call shr_sys_flush(shrlogunit)
       lsize = mct_gsmap_lsize(gsmap,mpicom)

       if (atm_present) then
          call mct_rearr_init(SDATM%gsmap,gsmap,mpicom,rearr)
       endif

       call t_stopf('datm_initgsmaps')

       !----------------------------------------------------------------------------
       ! Initialize MCT domain
       !----------------------------------------------------------------------------

       call t_startf('datm_initmctdom')
       if (my_task == master_task) write(logunit,F00) 'copy domains'
       call shr_sys_flush(logunit)

       if (atm_present)then
          call shr_dmodel_rearrGGrid(SDATM%grid, ggrid, gsmap, rearr, mpicom)
          call shr_sys_flush(shrlogunit)
       end if

       call t_stopf('datm_initmctdom')

       !----------------------------------------------------------------------------
       ! Initialize MCT attribute vectors
       !----------------------------------------------------------------------------

       call t_startf('datm_initmctavs')
       if (my_task == master_task) write(logunit,F00) 'allocate AVs'
       call shr_sys_flush(logunit)

       call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)
       call mct_aVect_zero(a2x)

       kz    = mct_aVect_indexRA(a2x,'Sa_z')
       ktopo = mct_aVect_indexRA(a2x,'Sa_topo')
       ku    = mct_aVect_indexRA(a2x,'Sa_u')
       kv    = mct_aVect_indexRA(a2x,'Sa_v')
       ktbot = mct_aVect_indexRA(a2x,'Sa_tbot')
       kptem = mct_aVect_indexRA(a2x,'Sa_ptem')
       kshum = mct_aVect_indexRA(a2x,'Sa_shum')
       kdens = mct_aVect_indexRA(a2x,'Sa_dens')
       kpbot = mct_aVect_indexRA(a2x,'Sa_pbot')
       kpslv = mct_aVect_indexRA(a2x,'Sa_pslv')
       klwdn = mct_aVect_indexRA(a2x,'Faxa_lwdn')
       krc   = mct_aVect_indexRA(a2x,'Faxa_rainc')
       krl   = mct_aVect_indexRA(a2x,'Faxa_rainl')
       ksc   = mct_aVect_indexRA(a2x,'Faxa_snowc')
       ksl   = mct_aVect_indexRA(a2x,'Faxa_snowl')
       kswndr= mct_aVect_indexRA(a2x,'Faxa_swndr')
       kswndf= mct_aVect_indexRA(a2x,'Faxa_swndf')
       kswvdr= mct_aVect_indexRA(a2x,'Faxa_swvdr')
       kswvdf= mct_aVect_indexRA(a2x,'Faxa_swvdf')
       kswnet= mct_aVect_indexRA(a2x,'Faxa_swnet')
       kco2p = mct_aVect_indexRA(a2x,'Sa_co2prog',perrWith='quiet')
       kco2d = mct_aVect_indexRA(a2x,'Sa_co2diag',perrWith='quiet')

       kbid  = mct_aVect_indexRA(a2x,'Faxa_bcphidry')
       kbod  = mct_aVect_indexRA(a2x,'Faxa_bcphodry')
       kbiw  = mct_aVect_indexRA(a2x,'Faxa_bcphiwet')
       koid  = mct_aVect_indexRA(a2x,'Faxa_ocphidry')
       kood  = mct_aVect_indexRA(a2x,'Faxa_ocphodry')
       koiw  = mct_aVect_indexRA(a2x,'Faxa_ocphiwet')
       kdd1  = mct_aVect_indexRA(a2x,'Faxa_dstdry1')
       kdd2  = mct_aVect_indexRA(a2x,'Faxa_dstdry2')
       kdd3  = mct_aVect_indexRA(a2x,'Faxa_dstdry3')
       kdd4  = mct_aVect_indexRA(a2x,'Faxa_dstdry4')
       kdw1  = mct_aVect_indexRA(a2x,'Faxa_dstwet1')
       kdw2  = mct_aVect_indexRA(a2x,'Faxa_dstwet2')
       kdw3  = mct_aVect_indexRA(a2x,'Faxa_dstwet3')
       kdw4  = mct_aVect_indexRA(a2x,'Faxa_dstwet4')

       call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
       call mct_aVect_zero(x2a)

       kanidr = mct_aVect_indexRA(x2a,'Sx_anidr')
       kanidf = mct_aVect_indexRA(x2a,'Sx_anidf')
       kavsdr = mct_aVect_indexRA(x2a,'Sx_avsdr')
       kavsdf = mct_aVect_indexRA(x2a,'Sx_avsdf')

       !--- figure out what's on the standard streams ---
       cnt = 0
       flds_strm = ''
       do n = 1,SDATM%nstreams
          do k = 1,ktranss
             kfld = mct_aVect_indexRA(SDATM%avs(n),trim(stifld(k)),perrWith='quiet')
             if (kfld > 0) then
                cnt = cnt + 1
                if (cnt == 1) then
                   flds_strm = trim(stofld(k))
                else
                   flds_strm = trim(flds_strm)//':'//trim(stofld(k))
                endif
             endif
          enddo
       enddo

       if (my_task == master_task) write(logunit,F00) ' flds_strm = ',trim(flds_strm)
       call shr_sys_flush(logunit)

       call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
       call mct_aVect_zero(avstrm)

       stbot    = mct_aVect_indexRA(avstrm,'strm_tbot'    ,perrWith='quiet')
       swind    = mct_aVect_indexRA(avstrm,'strm_wind'    ,perrWith='quiet')
       sz       = mct_aVect_indexRA(avstrm,'strm_z'       ,perrWith='quiet')
       spbot    = mct_aVect_indexRA(avstrm,'strm_pbot'    ,perrWith='quiet')
       sshum    = mct_aVect_indexRA(avstrm,'strm_shum'    ,perrWith='quiet')
       stdew    = mct_aVect_indexRA(avstrm,'strm_tdew'    ,perrWith='quiet')
       srh      = mct_aVect_indexRA(avstrm,'strm_rh'      ,perrWith='quiet')
       slwdn    = mct_aVect_indexRA(avstrm,'strm_lwdn'    ,perrWith='quiet')
       sswdn    = mct_aVect_indexRA(avstrm,'strm_swdn'    ,perrWith='quiet')
       sswdndf  = mct_aVect_indexRA(avstrm,'strm_swdndf'  ,perrWith='quiet')
       sswdndr  = mct_aVect_indexRA(avstrm,'strm_swdndr'  ,perrWith='quiet')
       sprecc   = mct_aVect_indexRA(avstrm,'strm_precc'   ,perrWith='quiet')
       sprecl   = mct_aVect_indexRA(avstrm,'strm_precl'   ,perrWith='quiet')
       sprecn   = mct_aVect_indexRA(avstrm,'strm_precn'   ,perrWith='quiet')
       sco2p    = mct_aVect_indexRA(avstrm,'strm_co2p'    ,perrWith='quiet')
       sco2d    = mct_aVect_indexRA(avstrm,'strm_co2d'    ,perrWith='quiet')
       sswup    = mct_aVect_indexRA(avstrm,'strm_swup'    ,perrWith='quiet')
       sprec    = mct_aVect_indexRA(avstrm,'strm_prec'    ,perrWith='quiet')
       starcf   = mct_aVect_indexRA(avstrm,'strm_tarcf'   ,perrWith='quiet')

       ! anomaly forcing
       sprecsf  = mct_aVect_indexRA(avstrm,'strm_precsf'  ,perrWith='quiet')
       sprec_af = mct_aVect_indexRA(avstrm,'strm_prec_af' ,perrWith='quiet')
       su_af    = mct_aVect_indexRA(avstrm,'strm_u_af'    ,perrWith='quiet')
       sv_af    = mct_aVect_indexRA(avstrm,'strm_v_af'    ,perrWith='quiet')
       stbot_af = mct_aVect_indexRA(avstrm,'strm_tbot_af' ,perrWith='quiet')
       spbot_af = mct_aVect_indexRA(avstrm,'strm_pbot_af' ,perrWith='quiet')
       sshum_af = mct_aVect_indexRA(avstrm,'strm_shum_af' ,perrWith='quiet')
       sswdn_af = mct_aVect_indexRA(avstrm,'strm_swdn_af' ,perrWith='quiet')
       slwdn_af = mct_aVect_indexRA(avstrm,'strm_lwdn_af' ,perrWith='quiet')

       allocate(imask(lsize))
       allocate(yc(lsize))
       allocate(windFactor(lsize))
       allocate(winddFactor(lsize))
       allocate(qsatFactor(lsize))

       kmask = mct_aVect_indexRA(ggrid%data,'mask')
       imask(:) = nint(ggrid%data%rAttr(kmask,:))
       klat = mct_aVect_indexRA(ggrid%data,'lat')
       yc(:) = ggrid%data%rAttr(klat,:)

       call t_stopf('datm_initmctavs')

       !----------------------------------------------------------------------------
       ! Read restart
       !----------------------------------------------------------------------------

       if (read_restart) then
          if (trim(rest_file)      == trim(nullstr) .and. &
               trim(rest_file_strm) == trim(nullstr)) then
             if (my_task == master_task) then
                write(logunit,F00) ' restart filenames from rpointer'
                call shr_sys_flush(logunit)
                inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
                if (.not.exists) then
                   write(logunit,F00) ' ERROR: rpointer file does not exist'
                   call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
                endif
                nu = shr_file_getUnit()
                open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
                read(nu,'(a)') rest_file
                read(nu,'(a)') rest_file_strm
                close(nu)
                call shr_file_freeUnit(nu)
                inquire(file=trim(rest_file_strm),exist=exists)
             endif
             call shr_mpi_bcast(rest_file,mpicom,'rest_file')
             call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
          else
             ! use namelist already read
             if (my_task == master_task) then
                write(logunit,F00) ' restart filenames from namelist '
                call shr_sys_flush(logunit)
                inquire(file=trim(rest_file_strm),exist=exists)
             endif
          endif

          call shr_mpi_bcast(exists,mpicom,'exists')

          if (exists) then
             if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
             call shr_strdata_restRead(trim(rest_file_strm),SDATM,mpicom)
          else
             if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
          endif
          call shr_sys_flush(logunit)
       endif

       if (read_restart) then
          call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
          call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt )
          call seq_timemgr_EClockGetData( EClock, calendar=calendar )
          nextsw_cday = datm_shr_getNextRadCDay( CurrentYMD, CurrentTOD, stepno, idt, iradsw, calendar )
       else
          call seq_timemgr_EClockGetData( EClock, curr_cday=nextsw_cday, stepno=stepno )
       endif
       call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday )

    else

       call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
       call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt)
       call seq_timemgr_EClockGetData( EClock, calendar=calendar )
       nextsw_cday = datm_shr_getNextRadCDay( CurrentYMD, CurrentTOD, stepno, idt, iradsw, calendar )
       call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday )

    endif

    !----------------------------------------------------------------------------
    ! Set initial atm state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call datm_comp_run( EClock, cdata,  x2a, a2x)
    call t_adj_detailf(-2)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'datm_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DATM_INIT')

  end subroutine datm_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: datm_comp_run
!
! !DESCRIPTION:
!     run method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine datm_comp_run( EClock, cdata,  x2a, a2x)

   implicit none

   ! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2a        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: a2x        ! dead   -> driver

   !--- local ---
   type(mct_gsMap)        , pointer :: gsMap
   type(mct_gGrid)        , pointer :: ggrid

   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: yy,mm,dd          ! year month day
   integer(IN)   :: n                 ! indices
   integer(IN)   :: lsize             ! size of attr vect
   integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
   logical       :: mssrmlf           ! remove old data
   integer(IN)   :: idt               ! integer timestep
   real(R8)      :: dt                ! timestep
   logical       :: write_restart     ! restart now
   character(CL) :: case_name         ! case name
   character(CL) :: rest_file         ! restart_file
   character(CL) :: rest_file_strm    ! restart_file
   integer(IN)   :: nu                ! unit number
   integer(IN)   :: stepno            ! step number
   real(R8)      :: nextsw_cday       ! calendar of next atm sw
   integer(IN)   :: eday              ! elapsed day
   real(R8)      :: rday              ! elapsed day
   real(R8)      :: cosFactor         ! cosine factor
   real(R8)      :: factor            ! generic/temporary correction factor
   real(R8)      :: avg_alb           ! average albedo
   real(R8)      :: tMin              ! minimum temperature
   character(CL) :: calendar          ! calendar type

   !--- temporaries
   real(R8)      :: uprime,vprime,swndr,swndf,swvdr,swvdf,ratio_rvrf
   real(R8)      :: tbot,pbot,rtmp,vp,ea,e,qsat,frac,qsatT

   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(datm_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(datm_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(datm_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DATM_RUN')

   call t_startf('datm_run1')

  !----------------------------------------------------------------------------
  ! Reset shr logging to my log file
  !----------------------------------------------------------------------------
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (logUnit)

   call seq_cdata_setptrs(cdata, gsMap=gsmap, dom=ggrid)

   call seq_cdata_setptrs(cdata, infodata=infodata)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
   call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
   call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=idt)
   call seq_timemgr_EClockGetData( EClock, calendar=calendar)
   dt = idt * 1.0_r8
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   call t_stopf('datm_run1')

   !--------------------
   ! ADVANCE ATM
   !--------------------

   call t_barrierf('datm_BARRIER',mpicom)
   call t_startf('datm')

   nextsw_cday = datm_shr_getNextRadCDay( CurrentYMD, CurrentTOD, stepno, idt, iradsw, calendar )
   call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday )

   !--- copy all fields from streams to a2x as default ---

   if (trim(atm_mode) /= 'NULL') then
      call t_startf('datm_strdata_advance')
      call shr_strdata_advance(SDATM,currentYMD,currentTOD,mpicom,'datm')
      call t_stopf('datm_strdata_advance')
      call t_barrierf('datm_scatter_BARRIER',mpicom)
      call t_startf('datm_scatter')
      if (trim(atm_mode) /= 'COPYALL') then
         lsize = mct_avect_lsize(a2x)
         do n = 1,lsize
            a2x%rAttr(kbid,n) = aerodep_spval
            a2x%rAttr(kbod,n) = aerodep_spval
            a2x%rAttr(kbiw,n) = aerodep_spval
            a2x%rAttr(koid,n) = aerodep_spval
            a2x%rAttr(kood,n) = aerodep_spval
            a2x%rAttr(koiw,n) = aerodep_spval
            a2x%rAttr(kdd1,n) = aerodep_spval
            a2x%rAttr(kdd2,n) = aerodep_spval
            a2x%rAttr(kdd3,n) = aerodep_spval
            a2x%rAttr(kdd4,n) = aerodep_spval
            a2x%rAttr(kdw1,n) = aerodep_spval
            a2x%rAttr(kdw2,n) = aerodep_spval
            a2x%rAttr(kdw3,n) = aerodep_spval
            a2x%rAttr(kdw4,n) = aerodep_spval
         enddo
      endif
      if (firstcall) then
         allocate(ilist_av(SDATM%nstreams))
         allocate(olist_av(SDATM%nstreams))
         allocate(ilist_st(SDATM%nstreams))
         allocate(olist_st(SDATM%nstreams))
         allocate(count_av(SDATM%nstreams))
         allocate(count_st(SDATM%nstreams))
      end if
      do n = 1,SDATM%nstreams
         if (firstcall) then
            call shr_dmodel_translate_list(SDATM%avs(n),a2x,&
                 avifld(1:ktrans),avofld,ilist_av(n),olist_av(n),count_av(n))
         end if
         if (count_av(n) > 0) then
            call shr_dmodel_translateAV_list(SDATM%avs(n),a2x,&
                 ilist_av(n),olist_av(n),rearr)
         end if
      enddo
      do n = 1,SDATM%nstreams
         if (firstcall) then
            call shr_dmodel_translate_list(SDATM%avs(n),avstrm,&
                 stifld(1:ktranss),stofld,ilist_st(n),olist_st(n),count_st(n))
         end if
         if (count_st(n) > 0) then
            call shr_dmodel_translateAV_list(SDATM%avs(n),avstrm,&
                 ilist_st(n),olist_st(n),rearr)
         end if
      enddo

      call t_stopf('datm_scatter')
   else
      call mct_aVect_zero(a2x)
   endif

   call t_startf('datm_mode')

   select case (trim(atm_mode))

   case('COPYALL')
      ! do nothing extra

   case('CPLHIST')
      ! do nothing extra

   case ('WRF')
      lsize = mct_avect_lsize(a2x)
      do n = 1,lsize

         !--- fabricate required swdn components from total swdn ---
         a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdn,n)*(0.28_R8)
         a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdn,n)*(0.31_R8)
         a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdn,n)*(0.24_R8)
         a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdn,n)*(0.17_R8)

         !--- just a diagnostic, not really needed
         a2x%rAttr(kswnet,n) = avstrm%rAttr(sswdn,n)-avstrm%rAttr(sswup,n)

         !--- convert from hPa
         a2x%rAttr(kpslv,n) = a2x%rAttr(kpslv,n)*100._R8
         a2x%rAttr(kpbot,n) = a2x%rAttr(kpbot,n)*100._R8

         !--- tcraig, file has terrain height on it, set to 10m
         a2x%rAttr(kz,n) = 10.0_R8

         !--- convert to degK from degC
         a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + tKFrz

      enddo

   case('CORE2_NYF','CORE2_IAF')
      if (firstcall) then
         if (sprec < 1 .or. sswdn < 1) then
            write(logunit,F00) 'ERROR: prec and swdn must be in streams for CORE2'
            call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE2')
         endif
         if (trim(atm_mode) == 'CORE2_IAF' ) then
            if (starcf < 1 ) then
               write(logunit,F00) 'ERROR: tarcf must be in an input stream for CORE2_IAF'
               call shr_sys_abort(trim(subname)//'tarcf must be in an input stream for CORE2_IAF')
            endif
         endif
         call datm_shr_CORE2getFactors(factorFn,windFactor,winddFactor,qsatFactor, &
              mpicom,compid,gsmap,ggrid,SDATM%nxg,SDATM%nyg)
      endif
      call shr_cal_date2julian(currentYMD,currentTOD,rday,calendar)
      rday = mod((rday - 1.0_R8),365.0_R8)
      cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

      lsize = mct_avect_lsize(a2x)
      do n = 1,lsize
         a2x%rAttr(kz,n) = 10.0_R8

         !--- correction to NCEP winds based on QSCAT ---
         uprime    = a2x%rAttr(ku,n)*windFactor(n)
         vprime    = a2x%rAttr(kv,n)*windFactor(n)
         a2x%rAttr(ku,n) = uprime*cos(winddFactor(n)*degtorad)- &
                           vprime*sin(winddFactor(n)*degtorad)
         a2x%rAttr(kv,n) = uprime*sin(winddFactor(n)*degtorad)+ &
                           vprime*cos(winddFactor(n)*degtorad)

         !--- density, tbot, & pslv taken directly from input stream, set pbot ---
         a2x%rAttr(kpbot,n) = a2x%rAttr(kpslv,n)

         !--- correction to NCEP Arctic & Antarctic air T & potential T ---
         if      ( yc(n) < -60.0_R8 ) then
            tMin = (avg_c0 + avg_c1*yc(n)) + (amp_c0 + amp_c1*yc(n))*cosFactor + tKFrz
            a2x%rAttr(ktbot,n) = max(a2x%rAttr(ktbot,n), tMin)
         else if ( yc(n) > 60.0_R8 ) then
            factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
            a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + factor * dTarc(mm)
         endif
         a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)

         !---  correction to NCEP relative humidity for heat budget balance ---
         a2x%rAttr(kshum,n) = a2x%rAttr(kshum,n) + qsatFactor(n)

         !--- Dupont correction to NCEP Arctic air T  ---
         !--- don't correct during summer months (July-September)
         !--- ONLY correct when forcing year is 1997->2004
         if (trim(atm_mode) == 'CORE2_IAF' ) then
            a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) +  avstrm%rAttr(starcf,n)
            a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)
         end if

         !-------------------------------------------------------------------------
         ! PRECIPITATION DATA
         !-------------------------------------------------------------------------

         avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)/86400.0_R8        ! convert mm/day to kg/m^2/s

         !  only correct satellite products, do not correct Serreze Arctic data
         if ( yc(n) < 58. ) then
            avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)*1.14168_R8
         endif
         if ( yc(n) >= 58. .and. yc(n) < 68. ) then
            factor = MAX(0.0_R8, 1.0_R8 - 0.1_R8*(yc(n)-58.0_R8) )
            avstrm%rAttr(sprec,n) = avstrm%rAttr(sprec,n)*(factor*(1.14168_R8 - 1.0_R8) + 1.0_R8)
         endif

         a2x%rAttr(krc,n) = 0.0_R8                    ! default zero
         a2x%rAttr(ksc,n) = 0.0_R8
         if (a2x%rAttr(ktbot,n) < tKFrz ) then        ! assign precip to rain/snow components
            a2x%rAttr(krl,n) = 0.0_R8
            a2x%rAttr(ksl,n) = avstrm%rAttr(sprec,n)
         else
            a2x%rAttr(krl,n) = avstrm%rAttr(sprec,n)
            a2x%rAttr(ksl,n) = 0.0_R8
         endif

         !-------------------------------------------------------------------------
         ! RADIATION DATA
         !-------------------------------------------------------------------------

         !--- fabricate required swdn components from net swdn ---
         a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdn,n)*(0.28_R8)
         a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdn,n)*(0.31_R8)
         a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdn,n)*(0.24_R8)
         a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdn,n)*(0.17_R8)

         !--- compute net short-wave based on LY08 latitudinally-varying albedo ---
         avg_alb = ( 0.069 - 0.011*cos(2.0_R8*yc(n)*degtorad ) )
         a2x%rAttr(kswnet,n) = avstrm%rAttr(sswdn,n)*(1.0_R8 - avg_alb)

         !--- corrections to GISS sswdn for heat budget balancing ---
         factor = 1.0_R8
         if      ( -60.0_R8 < yc(n) .and. yc(n) < -50.0_R8 ) then
            factor = 1.0_R8 - (yc(n) + 60.0_R8)*(0.05_R8/10.0_R8)
         else if ( -50.0_R8 < yc(n) .and. yc(n) <  30.0_R8 ) then
            factor = 0.95_R8
         else if (  30.0_R8 < yc(n) .and. yc(n) <  40._R8 ) then
            factor = 1.0_R8 - (40.0_R8 - yc(n))*(0.05_R8/10.0_R8)
         endif
         a2x%rAttr(kswnet,n) = a2x%rAttr(kswnet,n)*factor
         a2x%rAttr(kswvdr,n) = a2x%rAttr(kswvdr,n)*factor
         a2x%rAttr(kswndr,n) = a2x%rAttr(kswndr,n)*factor
         a2x%rAttr(kswvdf,n) = a2x%rAttr(kswvdf,n)*factor
         a2x%rAttr(kswndf,n) = a2x%rAttr(kswndf,n)*factor

         !--- correction to GISS lwdn in Arctic ---
         if ( yc(n) > 60._R8 ) then
            factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
            a2x%rAttr(klwdn,n) = a2x%rAttr(klwdn,n) + factor * dLWarc
         endif

      enddo   ! lsize

   case('CLMNCEP')
      if (firstcall) then
         if (swind < 1 .or. stbot < 1) then
            write(logunit,F00) ' ERROR: wind and tbot must be in streams for CLMNCEP'
            call shr_sys_abort(trim(subname)//' ERROR: wind and tbot must be in streams for CLMNCEP')
         endif
         rtmp = maxval(a2x%rAttr(ktbot,:))
         call shr_mpi_max(rtmp,tbotmax,mpicom,'datm_tbot',all=.true.)
         rtmp = maxval(x2a%rAttr(kanidr,:))
         call shr_mpi_max(rtmp,anidrmax,mpicom,'datm_ani',all=.true.)
         if (stdew > 0) then
            rtmp = maxval(avstrm%rAttr(stdew,:))
            call shr_mpi_max(rtmp,tdewmax,mpicom,'datm_tdew',all=.true.)
         endif
         if (my_task == master_task) &
             write(logunit,*) trim(subname),' max values = ',tbotmax,tdewmax,anidrmax
      endif
      lsize = mct_avect_lsize(a2x)
      do n = 1,lsize
         !--- bottom layer height ---
         if (sz < 1) a2x%rAttr(kz,n) = 30.0_R8

         !--- temperature ---
         if (tbotmax < 50.0_R8) a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + tkFrz
         a2x%rAttr(kptem,n) = a2x%rAttr(ktbot,n)

         !--- pressure ---
         if (spbot < 1) a2x%rAttr(kpbot,n) = pstd
         a2x%rAttr(kpslv,n) = a2x%rAttr(kpbot,n)

         !--- u, v wind velocity ---
         a2x%rAttr(ku,n) = avstrm%rAttr(swind,n)/sqrt(2.0_R8)
         a2x%rAttr(kv,n) = a2x%rAttr(ku,n)

         !--- specific humidity ---
         tbot = a2x%rAttr(ktbot,n)
         pbot = a2x%rAttr(kpbot,n)
         if (sshum > 0) then
            e = datm_shr_esat(tbot,tbot)
            qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
            if (qsat < a2x%rAttr(kshum,n)) then
               a2x%rAttr(kshum,n) = qsat
            endif
         else if (srh > 0) then
            e = avstrm%rAttr(srh,n) * 0.01_R8 * datm_shr_esat(tbot,tbot)
            qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
            a2x%rAttr(kshum,n) = qsat
         else if (stdew > 0) then
            if (tdewmax < 50.0_R8) avstrm%rAttr(stdew,n) = avstrm%rAttr(stdew,n) + tkFrz
            e = datm_shr_esat(avstrm%rAttr(stdew,n),tbot)
            qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
            a2x%rAttr(kshum,n) = qsat
         else
           call shr_sys_abort(subname//'ERROR: cannot compute shum')
         endif

         !--- density ---
         vp = (a2x%rAttr(kshum,n)*pbot) / (0.622_R8 + 0.378_R8 * a2x%rAttr(kshum,n))
         a2x%rAttr(kdens,n) = (pbot - 0.378_R8 * vp) / (tbot*rdair)

         !--- downward longwave ---
         if (slwdn < 1) then
            e  = a2x%rAttr(kpslv,n) * a2x%rAttr(kshum,n) / (0.622_R8 + 0.378_R8 * a2x%rAttr(kshum,n))
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            a2x%rAttr(klwdn,n) = ea * stebol * tbot**4
         endif

         !--- shortwave radiation ---
         if (sswdndf > 0 .and. sswdndr > 0) then
            a2x%rAttr(kswndr,n) = avstrm%rAttr(sswdndr,n) * 0.50_R8
            a2x%rAttr(kswvdr,n) = avstrm%rAttr(sswdndr,n) * 0.50_R8
            a2x%rAttr(kswndf,n) = avstrm%rAttr(sswdndf,n) * 0.50_R8
            a2x%rAttr(kswvdf,n) = avstrm%rAttr(sswdndf,n) * 0.50_R8
         elseif (sswdn > 0) then
            ! relationship between incoming NIR or VIS radiation and ratio of
            ! direct to diffuse radiation calculated based on one year's worth of
            ! hourly CAM output from CAM version cam3_5_55
            swndr = avstrm%rAttr(sswdn,n) * 0.50_R8
            ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr  &
                           -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
            a2x%rAttr(kswndr,n) = ratio_rvrf*swndr
            swndf = avstrm%rAttr(sswdn,n) * 0.50_R8
            a2x%rAttr(kswndf,n) = (1._R8 - ratio_rvrf)*swndf

            swvdr = avstrm%rAttr(sswdn,n) * 0.50_R8
            ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                           -9.0039e-06_R8*swvdr**2 + 8.1351e-09_R8*swvdr**3,0.01_R8))
            a2x%rAttr(kswvdr,n) = ratio_rvrf*swvdr
            swvdf = avstrm%rAttr(sswdn,n) * 0.50_R8
            a2x%rAttr(kswvdf,n) = (1._R8 - ratio_rvrf)*swvdf
         else
            call shr_sys_abort(subName//'ERROR: cannot compute short-wave down')
         endif

         !--- swnet: a diagnostic quantity ---
         if (anidrmax < 1.0e-8 .or. anidrmax > SHR_CONST_SPVAL * 0.9_R8) then
            a2x%rAttr(kswnet,n) = 0.0_R8
         else
            a2x%rAttr(kswnet,n) = (1.0_R8-x2a%rAttr(kanidr,n))*a2x%rAttr(kswndr,n) + &
                                  (1.0_R8-x2a%rAttr(kavsdr,n))*a2x%rAttr(kswvdr,n) + &
                                  (1.0_R8-x2a%rAttr(kanidf,n))*a2x%rAttr(kswndf,n) + &
                                  (1.0_R8-x2a%rAttr(kavsdf,n))*a2x%rAttr(kswvdf,n)
         endif

         !--- rain and snow ---
         if (sprecc > 0 .and. sprecl > 0) then
            a2x%rAttr(krc,n) = avstrm%rAttr(sprecc,n)
            a2x%rAttr(krl,n) = avstrm%rAttr(sprecl,n)
         elseif (sprecn > 0) then
            a2x%rAttr(krc,n) = avstrm%rAttr(sprecn,n)*0.1_R8
            a2x%rAttr(krl,n) = avstrm%rAttr(sprecn,n)*0.9_R8
         else
            call shr_sys_abort(subName//'ERROR: cannot compute rain and snow')
         endif

         !--- split precip between rain & snow ---
         call shr_precip_partition_rain_snow_ramp(tbot, frac)
         a2x%rAttr(ksc,n) = max(0.0_R8, a2x%rAttr(krc,n)*(1.0_R8 - frac) )
         a2x%rAttr(ksl,n) = max(0.0_R8, a2x%rAttr(krl,n)*(1.0_R8 - frac) )
         a2x%rAttr(krc,n) = max(0.0_R8, a2x%rAttr(krc,n)*(         frac) )
         a2x%rAttr(krl,n) = max(0.0_R8, a2x%rAttr(krl,n)*(         frac) )

      enddo

   end select

   call t_stopf('datm_mode')

   !
   ! bias correction / anomaly forcing ( start block )
   ! modify atmospheric input fields if streams exist
   !
   lsize = mct_avect_lsize(avstrm)

   ! bias correct precipitation relative to observed
   ! (via bias_correct nameslist option)
   if (sprecsf > 0) then
      do n = 1,lsize
         a2x%rAttr(ksc,n) = a2x%rAttr(ksc,n)*min(1.e2_r8,avstrm%rAttr(sprecsf,n))
         a2x%rAttr(ksl,n) = a2x%rAttr(ksl,n)*min(1.e2_r8,avstrm%rAttr(sprecsf,n))
         a2x%rAttr(krc,n) = a2x%rAttr(krc,n)*min(1.e2_r8,avstrm%rAttr(sprecsf,n))
         a2x%rAttr(krl,n) = a2x%rAttr(krl,n)*min(1.e2_r8,avstrm%rAttr(sprecsf,n))

      end do
   endif

   ! adjust atmospheric input fields if anomaly forcing streams exist
   ! (via anomaly_forcing nameslist option)

   ! wind
   if (su_af > 0 .and. sv_af > 0) then
      do n = 1,lsize
         a2x%rAttr(ku,n) = a2x%rAttr(ku,n) + avstrm%rAttr(su_af,n)
         a2x%rAttr(kv,n) = a2x%rAttr(kv,n) + avstrm%rAttr(sv_af,n)
      end do
   endif

   ! specific humidity
   if (sshum_af > 0) then
      do n = 1,lsize
         a2x%rAttr(kshum,n) = a2x%rAttr(kshum,n) + avstrm%rAttr(sshum_af,n)

         ! avoid possible negative q values
         if(a2x%rAttr(kshum,n) < 0._r8) then
            a2x%rAttr(kshum,n) = 1.e-6_r8
         endif

      end do
   endif

   ! pressure
   if (spbot_af > 0) then
      do n = 1,lsize
         a2x%rAttr(kpbot,n) = a2x%rAttr(kpbot,n) + avstrm%rAttr(spbot_af,n)
      end do
   endif

   ! temperature
   if (stbot_af > 0) then
      do n = 1,lsize
         a2x%rAttr(ktbot,n) = a2x%rAttr(ktbot,n) + avstrm%rAttr(stbot_af,n)
      end do
   endif

   ! longwave
   if (slwdn_af > 0) then
      do n = 1,lsize
         a2x%rAttr(klwdn,n) = a2x%rAttr(klwdn,n)*avstrm%rAttr(slwdn_af,n)
      end do
   endif

   ! precipitation
   if (sprec_af > 0) then
      do n = 1,lsize
         a2x%rAttr(ksc,n) = a2x%rAttr(ksc,n)*avstrm%rAttr(sprec_af,n)
         a2x%rAttr(ksl,n) = a2x%rAttr(ksl,n)*avstrm%rAttr(sprec_af,n)
         a2x%rAttr(krc,n) = a2x%rAttr(krc,n)*avstrm%rAttr(sprec_af,n)
         a2x%rAttr(krl,n) = a2x%rAttr(krl,n)*avstrm%rAttr(sprec_af,n)
      enddo
   endif
   ! shortwave
   if (sswdn_af > 0) then
      do n = 1,lsize
         a2x%rAttr(kswndr,n) = a2x%rAttr(kswndr,n)*avstrm%rAttr(sswdn_af,n)
         a2x%rAttr(kswvdr,n) = a2x%rAttr(kswvdr,n)*avstrm%rAttr(sswdn_af,n)
         a2x%rAttr(kswndf,n) = a2x%rAttr(kswndf,n)*avstrm%rAttr(sswdn_af,n)
         a2x%rAttr(kswvdf,n) = a2x%rAttr(kswvdf,n)*avstrm%rAttr(sswdn_af,n)
      enddo
   endif
   !
   ! bias correction / anomaly forcing ( end block )
   !

   if (write_restart) then
      call t_startf('datm_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.datm'//trim(inst_suffix)//'.r.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.datm'//trim(inst_suffix)//'.rs1.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm
         close(nu)
         call shr_file_freeUnit(nu)
      endif

      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm),SDATM,mpicom,trim(case_name),'SDATM strdata')
      call shr_sys_flush(logunit)
      call t_stopf('datm_restart')

   endif

   call t_stopf('datm')

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('datm_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if

   firstcall = .false.

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('datm_run2')

   call t_stopf('DATM_RUN')

end subroutine datm_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: datm_comp_final
!
! !DESCRIPTION:
!     finalize method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine datm_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(datm_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(datm_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(datm_comp_final) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DATM_FINAL')

   if (my_task == master_task) then
      write(logunit,F91)
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91)
   end if

   call t_stopf('DATM_FINAL')

end subroutine datm_comp_final
!===============================================================================
!===============================================================================


end module datm_comp_mod
