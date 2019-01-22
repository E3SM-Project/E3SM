#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dice_comp_mod

  ! !USES:
  use esmf
  use mct_mod
  use perf_mod
  use shr_pcdf_mod
  use shr_const_mod
  use shr_sys_mod
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod    , only: shr_file_getunit, shr_file_freeunit
  use shr_cal_mod     , only: shr_cal_date2julian
  use shr_mpi_mod     , only: shr_mpi_bcast
  use shr_frz_mod     , only: shr_frz_freezetemp
  use shr_cal_mod     , only: shr_cal_ymd2julian
  use shr_strdata_mod , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod  , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod  , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV
  use seq_timemgr_mod , only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn

  use dice_shr_mod   , only: datamode       ! namelist input
  use dice_shr_mod   , only: decomp         ! namelist input
  use dice_shr_mod   , only: rest_file      ! namelist input
  use dice_shr_mod   , only: rest_file_strm ! namelist input
  use dice_shr_mod   , only: flux_swpf      ! namelist input -short-wave penatration factor
  use dice_shr_mod   , only: flux_Qmin      ! namelist input -bound on melt rate
  use dice_shr_mod   , only: flux_Qacc      ! namelist input -activates water accumulation/melt wrt Q
  use dice_shr_mod   , only: flux_Qacc0     ! namelist input -initial water accumulation value
  use dice_shr_mod   , only: nullstr
  use dice_flux_atmice_mod, only: dice_flux_atmice

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dice_comp_init
  public :: dice_comp_run
  public :: dice_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(CS) :: myModelName = 'ice'   ! user defined model name
  logical       :: firstcall = .true.    ! first call logical
  character(len=*),parameter :: rpfile = 'rpointer.ice'

  real(R8),parameter  :: pi     = shr_const_pi      ! pi
  real(R8),parameter  :: spval  = shr_const_spval   ! flags invalid data
  real(R8),parameter  :: tFrz   = shr_const_tkfrz   ! temp of freezing
  real(R8),parameter  :: latice = shr_const_latice  ! latent heat of fusion
  real(R8),parameter  :: cDay   = shr_const_cDay    ! sec in calendar day
  real(R8),parameter  :: waterMax = 1000.0_R8        ! wrt iFrac comp & frazil ice (kg/m^2)

  !----- surface albedo constants ------
  real(R8),parameter  :: snwfrac = 0.286_R8 ! snow cover fraction ~ [0,1]
  real(R8),parameter  :: as_nidf = 0.950_R8 ! albedo: snow,near-infr,diffuse
  real(R8),parameter  :: as_vsdf = 0.700_R8 ! albedo: snow,visible  ,diffuse
  real(R8),parameter  :: as_nidr = 0.960_R8 ! albedo: snow,near-infr,direct
  real(R8),parameter  :: as_vsdr = 0.800_R8 ! albedo: snow,visible  ,direct
  real(R8),parameter  :: ai_nidf = 0.700_R8 ! albedo: ice, near-infr,diffuse
  real(R8),parameter  :: ai_vsdf = 0.500_R8 ! albedo: ice, visible  ,diffuse
  real(R8),parameter  :: ai_nidr = 0.700_R8 ! albedo: ice, near-infr,direct
  real(R8),parameter  :: ai_vsdr = 0.500_R8 ! albedo: ice, visible  ,direct
  real(R8),parameter  :: ax_nidf = ai_nidf*(1.0_R8-snwfrac) + as_nidf*snwfrac
  real(R8),parameter  :: ax_vsdf = ai_vsdf*(1.0_R8-snwfrac) + as_vsdf*snwfrac
  real(R8),parameter  :: ax_nidr = ai_nidr*(1.0_R8-snwfrac) + as_nidr*snwfrac
  real(R8),parameter  :: ax_vsdr = ai_vsdr*(1.0_R8-snwfrac) + as_vsdr*snwfrac

  integer(IN) :: kswvdr,kswndr,kswvdf,kswndf,kq,kz,kua,kva,kptem,kshum,kdens,ktbot
  integer(IN) :: kiFrac,kt,kavsdr,kanidr,kavsdf,kanidf,kswnet,kmelth,kmeltw
  integer(IN) :: ksen,klat,klwup,kevap,ktauxa,ktauya,ktref,kqref,kswpen,ktauxo,ktauyo,ksalt
  integer(IN) :: ksalinity
  integer(IN) :: kbcpho, kbcphi, kflxdst
  integer(IN) :: kbcphidry, kbcphodry, kbcphiwet, kocphidry, kocphodry, kocphiwet
  integer(IN) :: kdstdry1, kdstdry2, kdstdry3, kdstdry4, kdstwet1, kdstwet2, kdstwet3, kdstwet4

  ! optional per thickness category fields
  integer(IN) :: kiFrac_01,kswpen_iFrac_01

  type(mct_rearr) :: rearr
  !  type(mct_avect) :: avstrm   ! av of data from stream
  integer(IN) , pointer :: imask(:)
  real(R8)    , pointer :: yc(:)
  real(R8)    , pointer :: water(:)
  real(R8)    , pointer :: tfreeze(:)
  !  real(R8)    , pointer :: ifrac0(:)

  !--------------------------------------------------------------------------
  integer(IN),parameter :: ktrans = 1
  character(16),parameter  :: avofld(1:ktrans) = (/"Si_ifrac        "/)
  character(16),parameter  :: avifld(1:ktrans) = (/"ifrac           "/)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine dice_comp_init(Eclock, x2i, i2x, &
       seq_flds_x2i_fields, seq_flds_i2x_fields, seq_flds_i2o_per_cat, &
       SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon)

    ! !DESCRIPTION: initialize dice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2i, i2x             ! input/output attribute vectors
    character(len=*)       , intent(in)    :: seq_flds_x2i_fields  ! fields from mediator
    character(len=*)       , intent(in)    :: seq_flds_i2x_fields  ! fields to mediator
    logical                , intent(in)    :: seq_flds_i2o_per_cat ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE                ! dice shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap                ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid                ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom               ! mpi communicator
    integer(IN)            , intent(in)    :: compid               ! mct comp id
    integer(IN)            , intent(in)    :: my_task              ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task          ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix          ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name            ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit              ! logging unit number
    logical                , intent(in)    :: read_restart         ! start from restart
    logical                , intent(in)    :: scmMode              ! single column mode
    real(R8)               , intent(in)    :: scmLat               ! single column lat
    real(R8)               , intent(in)    :: scmLon               ! single column lon

    !--- local variables ---
    integer(IN)   :: lsize       ! local size
    integer(IN)   :: kfld        ! field reference
    logical       :: exists      ! file existance logical
    integer(IN)   :: nu          ! unit number
    character(CL) :: calendar    ! calendar type

    !--- formats ---
    character(*), parameter :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dice_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dice_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dice_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dice_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(dice_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(dice_comp_init) ',a,2f10.4)"
    character(*), parameter :: F06   = "('(dice_comp_init) ',a,5l3)"
    character(*), parameter :: F90   = "('(dice_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dice_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dice_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DICE_INIT')

    call shr_strdata_pioinit(SDICE, COMPID)

    !----------------------------------------------------------------------------
    ! Initialize SDICE
    !----------------------------------------------------------------------------

    call t_startf('dice_strdata_init')

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDICE%gsmap and SDICE%ggrid. DICE%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    if (scmmode) then
       if (my_task == master_task) then
          write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       end if
       call shr_strdata_init(SDICE,mpicom,compid,name='ice', &
            scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, &
            calendar=calendar)
    else
       call shr_strdata_init(SDICE,mpicom,compid,name='ice', &
            calendar=calendar)
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
    endif

    call t_stopf('dice_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dice_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the dice_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap, SDICE%nxg*SDICE%nyg, compid, mpicom, decomp)
    lsize = mct_gsmap_lsize(gsmap, mpicom)

    ! create a rearranger from the data model SDICE%gsmap to gsmap
    call mct_rearr_init(SDICE%gsmap, gsmap, mpicom, rearr)
    call t_stopf('dice_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    call shr_dmodel_rearrGGrid(SDICE%grid, ggrid, gsmap, rearr, mpicom)
    call t_stopf('dice_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(i2x, rList=seq_flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x)

    kiFrac = mct_aVect_indexRA(i2x,'Si_ifrac')
    kt     = mct_aVect_indexRA(i2x,'Si_t')
    ktref  = mct_aVect_indexRA(i2x,'Si_tref')
    kqref  = mct_aVect_indexRA(i2x,'Si_qref')
    kavsdr = mct_aVect_indexRA(i2x,'Si_avsdr')
    kanidr = mct_aVect_indexRA(i2x,'Si_anidr')
    kavsdf = mct_aVect_indexRA(i2x,'Si_avsdf')
    kanidf = mct_aVect_indexRA(i2x,'Si_anidf')
    kswnet = mct_aVect_indexRA(i2x,'Faii_swnet')
    ksen   = mct_aVect_indexRA(i2x,'Faii_sen')
    klat   = mct_aVect_indexRA(i2x,'Faii_lat')
    klwup  = mct_aVect_indexRA(i2x,'Faii_lwup')
    kevap  = mct_aVect_indexRA(i2x,'Faii_evap')
    ktauxa = mct_aVect_indexRA(i2x,'Faii_taux')
    ktauya = mct_aVect_indexRA(i2x,'Faii_tauy')
    kmelth = mct_aVect_indexRA(i2x,'Fioi_melth')
    kmeltw = mct_aVect_indexRA(i2x,'Fioi_meltw')
    kswpen = mct_aVect_indexRA(i2x,'Fioi_swpen')
    ktauxo = mct_aVect_indexRA(i2x,'Fioi_taux')
    ktauyo = mct_aVect_indexRA(i2x,'Fioi_tauy')
    ksalt  = mct_aVect_indexRA(i2x,'Fioi_salt')
    kbcpho = mct_aVect_indexRA(i2x,'Fioi_bcpho')
    kbcphi = mct_aVect_indexRA(i2x,'Fioi_bcphi')
    kflxdst= mct_aVect_indexRA(i2x,'Fioi_flxdst')

    ! optional per thickness category fields
    if (seq_flds_i2o_per_cat) then
       kiFrac_01       = mct_aVect_indexRA(i2x,'Si_ifrac_01')
       kswpen_iFrac_01 = mct_aVect_indexRA(i2x,'PFioi_swpen_ifrac_01')
    end if

    call mct_aVect_init(x2i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i)

    kswvdr    = mct_aVect_indexRA(x2i,'Faxa_swvdr')
    kswndr    = mct_aVect_indexRA(x2i,'Faxa_swndr')
    kswvdf    = mct_aVect_indexRA(x2i,'Faxa_swvdf')
    kswndf    = mct_aVect_indexRA(x2i,'Faxa_swndf')
    kq        = mct_aVect_indexRA(x2i,'Fioo_q')
    kz        = mct_aVect_indexRA(x2i,'Sa_z')
    kua       = mct_aVect_indexRA(x2i,'Sa_u')
    kva       = mct_aVect_indexRA(x2i,'Sa_v')
    kptem     = mct_aVect_indexRA(x2i,'Sa_ptem')
    kshum     = mct_aVect_indexRA(x2i,'Sa_shum')
    kdens     = mct_aVect_indexRA(x2i,'Sa_dens')
    ktbot     = mct_aVect_indexRA(x2i,'Sa_tbot')
    ksalinity = mct_aVect_indexRA(x2i,'So_s')
    kbcphidry = mct_aVect_indexRA(x2i,'Faxa_bcphidry')
    kbcphodry = mct_aVect_indexRA(x2i,'Faxa_bcphodry')
    kbcphiwet = mct_aVect_indexRA(x2i,'Faxa_bcphiwet')
    kocphidry = mct_aVect_indexRA(x2i,'Faxa_ocphidry')
    kocphodry = mct_aVect_indexRA(x2i,'Faxa_ocphodry')
    kocphiwet = mct_aVect_indexRA(x2i,'Faxa_ocphiwet')
    kdstdry1  = mct_aVect_indexRA(x2i,'Faxa_dstdry1')
    kdstdry2  = mct_aVect_indexRA(x2i,'Faxa_dstdry2')
    kdstdry3  = mct_aVect_indexRA(x2i,'Faxa_dstdry3')
    kdstdry4  = mct_aVect_indexRA(x2i,'Faxa_dstdry4')
    kdstwet1  = mct_aVect_indexRA(x2i,'Faxa_dstwet1')
    kdstwet2  = mct_aVect_indexRA(x2i,'Faxa_dstwet2')
    kdstwet3  = mct_aVect_indexRA(x2i,'Faxa_dstwet3')
    kdstwet4  = mct_aVect_indexRA(x2i,'Faxa_dstwet4')

    ! call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
    ! call mct_aVect_zero(avstrm)

    allocate(imask(lsize))
    allocate(yc(lsize))
    allocate(water(lsize))
    allocate(tfreeze(lsize))
    ! allocate(iFrac0(lsize))

    kfld = mct_aVect_indexRA(ggrid%data,'mask')
    imask(:) = nint(ggrid%data%rAttr(kfld,:))
    kfld = mct_aVect_indexRA(ggrid%data,'lat')
    yc(:) = ggrid%data%rAttr(kfld,:)

    call t_stopf('dice_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file) == trim(nullstr) .and. &
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
       if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       call shr_pcdf_readwrite('read',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file),mpicom,gsmap=gsmap,rf1=water,rf1n='water',io_format=SDICE%io_format)
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDICE,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! On initial call, x2i is unset, so set for use in run method
    !  These values should have no impact on the solution!!
    !----------------------------------------------------------------------------
    x2i%rAttr(kz,:)    = 10.0_R8
    x2i%rAttr(kua,:)   = 5.0_R8
    x2i%rAttr(kva,:)   = 5.0_R8
    x2i%rAttr(kptem,:) = 260.0_R8
    x2i%rAttr(ktbot,:) = 260.0_R8
    x2i%rAttr(kshum,:) = 0.0014_R8
    x2i%rAttr(kdens,:) = 1.3_R8

    !----------------------------------------------------------------------------
    ! Set initial ice state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dice_comp_run(EClock, x2i, i2x, &
         seq_flds_i2o_per_cat, &
         SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart)
    call t_adj_detailf(-2)

    call t_stopf('DICE_INIT')

  end subroutine dice_comp_init

  !===============================================================================
  subroutine dice_comp_run(EClock, x2i, i2x, &
       seq_flds_i2o_per_cat, &
       SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, case_name)
    use shr_cal_mod, only : shr_cal_ymdtod2string
    ! !DESCRIPTION: run method for dice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2i
    type(mct_aVect)        , intent(inout) :: i2x
    logical                , intent(in)    :: seq_flds_i2o_per_cat ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer(IN)            , intent(in)    :: mpicom               ! mpi communicator
    integer(IN)            , intent(in)    :: compid               ! mct comp id
    integer(IN)            , intent(in)    :: my_task              ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task          ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix          ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit              ! logging unit number
    logical                , intent(in)    :: read_restart         ! start from restart
    character(CL)          , intent(in), optional :: case_name     ! case name

    !--- local ---
    integer(IN)   :: CurrentYMD        ! model date
    integer(IN)   :: CurrentTOD        ! model sec into model date
    integer(IN)   :: yy,mm,dd          ! year month day
    integer(IN)   :: n                 ! indices
    integer(IN)   :: lsize             ! size of attr vect
    integer(IN)   :: idt               ! integer timestep
    real(R8)      :: dt                ! timestep
    integer(IN)   :: nu                ! unit number
    real(R8)      :: qmeltall          ! q that would melt all accumulated water
    real(R8)      :: cosarg            ! for setting ice temp pattern
    real(R8)      :: jday, jday0       ! elapsed day counters
    character(CS) :: calendar          ! calendar type
    logical       :: write_restart     ! restart now
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dice_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dice_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DICE_RUN')

    call t_startf('dice_run1')

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
    call seq_timemgr_EClockGetData( EClock, dtime=idt, calendar=calendar)
    dt = idt * 1.0_r8
    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

    call t_stopf('dice_run1')

    !--------------------
    ! ADVANCE ICE
    !--------------------

    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice')

    !--- copy all fields from streams to i2x as default ---

    if (trim(datamode) /= 'NULL') then
       call t_startf('dice_strdata_advance')
       call shr_strdata_advance(SDICE,currentYMD,currentTOD,mpicom,'dice')
       call t_stopf('dice_strdata_advance')
       call t_barrierf('dice_scatter_BARRIER',mpicom)
       call t_startf('dice_scatter')
       do n = 1,SDICE%nstreams
          call shr_dmodel_translateAV(SDICE%avs(n),i2x,avifld,avofld,rearr)
       enddo
       call t_stopf('dice_scatter')
    else
       call mct_aVect_zero(i2x)
    endif

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dice_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       if (firstcall .and. .not. read_restart) then
          ! iFrac0 = iFrac  ! previous step's ice fraction
          water  = 0.0_R8 ! previous step's water accumulation
          where (i2x%rAttr(kiFrac,:) > 0.0_R8) water(:) = flux_Qacc0
       endif

       ! tcraig, feb 10, 2012, ymd2eday no longer exists, use ymd2julian instead
       !   this could be improved for use in gregorian calendar
       !      call shr_cal_ymd2eday(0,mm,dd,eDay ,calendar)    ! model date
       !      call shr_cal_ymd2eday(0,09,01,eDay0,calendar)    ! sept 1st
       !      cosArg = 2.0_R8*pi*(real(eDay,R8) + real(currentTOD,R8)/cDay - real(eDay0,R8))/365.0_R8
       call shr_cal_ymd2julian(0,mm,dd,currentTOD,jDay ,calendar)    ! julian day for model
       call shr_cal_ymd2julian(0, 9, 1,0         ,jDay0,calendar)    ! julian day for Sept 1
       cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

       lsize = mct_avect_lsize(i2x)

       tfreeze = shr_frz_freezetemp(x2i%rAttr(ksalinity,:)) + tFrz ! convert to Kelvin

       do n = 1,lsize

          !--- fix erroneous iFrac ---
          i2x%rAttr(kiFrac,n) = min(1.0_R8,max(0.0_R8,i2x%rAttr(kiFrac,n)))

          !--- fabricate ice surface T, fix erroneous iFrac ---
          if ( yc(n) > 0.0_R8) then
             i2x%rAttr(kt,n) = 260.0_R8 + 10.0_R8*cos(cosArg)
          else
             i2x%rAttr(kt,n) = 260.0_R8 - 10.0_R8*cos(cosArg)
          end if

          !--- set albedos (constant) ---
          i2x%rAttr(kavsdr,n) = ax_vsdr
          i2x%rAttr(kanidr,n) = ax_nidr
          i2x%rAttr(kavsdf,n) = ax_vsdf
          i2x%rAttr(kanidf,n) = ax_nidf

          !--- swnet is sent to cpl as a diagnostic quantity only ---
          !--- newly recv'd swdn goes with previously sent albedo ---
          !--- but albedos are (currently) time invariant         ---
          i2x%rAttr(kswnet,n) = (1.0_R8 - i2x%rAttr(kavsdr,n))*x2i%rAttr(kswvdr,n) &
                              + (1.0_R8 - i2x%rAttr(kanidr,n))*x2i%rAttr(kswndr,n) &
                              + (1.0_R8 - i2x%rAttr(kavsdf,n))*x2i%rAttr(kswvdf,n) &
                              + (1.0_R8 - i2x%rAttr(kanidf,n))*x2i%rAttr(kswndf,n)

          !--- compute melt/freeze water balance, adjust iFrac  -------------
          if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
             i2x%rAttr(kmelth,n) = min(x2i%rAttr(kq,n),0.0_R8 ) ! q<0 => melt potential
             i2x%rAttr(kmelth,n) = max(i2x%rAttr(kmelth,n),Flux_Qmin   ) ! limit the melt rate
             i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice   ! corresponding water flux

          else                                 ! Q accumulation option is ON
             !--------------------------------------------------------------
             ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
             ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
             ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
             ! 2b) Q>0 & iFrac = 0  =>  accumulated water
             !--------------------------------------------------------------
             if ( x2i%rAttr(kq,n) <  0.0_R8 ) then ! Q<0 => melt
                if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                   i2x%rAttr(kmelth,n) = i2x%rAttr(kiFrac,n)*max(x2i%rAttr(kq,n),Flux_Qmin)
                   i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice
                   !  water(n) = < don't change this value >
                else
                   Qmeltall   = -water(n)*latice/dt
                   i2x%rAttr(kmelth,n) = max(x2i%rAttr(kq,n), Qmeltall, Flux_Qmin )
                   i2x%rAttr(kmeltw,n) = -i2x%rAttr(kmelth,n)/latice
                   water(n) =  water(n) - i2x%rAttr(kmeltw,n)*dt
                end if
             else                       ! Q>0 => freeze
                if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                   i2x%rAttr(kmelth,n) = 0.0_R8
                   i2x%rAttr(kmeltw,n) = 0.0_R8
                   water(n) = 0.0_R8
                else
                   i2x%rAttr(kmelth,n) = 0.0_R8
                   i2x%rAttr(kmeltw,n) = 0.0_R8
                   water(n) = water(n) + dt*x2i%rAttr(kq,n)/latice
                end if
             end if

             if (water(n) < 1.0e-16_R8 ) water(n) = 0.0_R8

             !--- non-zero water => non-zero iFrac ---
             if (i2x%rAttr(kiFrac,n) <= 0.0_R8  .and.  water(n) > 0.0_R8) then
                i2x%rAttr(kiFrac,n) = min(1.0_R8,water(n)/waterMax)
                ! i2x%rAttr(kT,n) = tfreeze(n)     ! T can be above freezing?!?
             end if

             !--- cpl multiplies melth & meltw by iFrac ---
             !--- divide by iFrac here => fixed quantity flux (not per area) ---
             if (i2x%rAttr(kiFrac,n) > 0.0_R8) then
                i2x%rAttr(kiFrac,n) = max( 0.01_R8, i2x%rAttr(kiFrac,n)) ! min iFrac
                i2x%rAttr(kmelth,n) = i2x%rAttr(kmelth,n)/i2x%rAttr(kiFrac,n)
                i2x%rAttr(kmeltw,n) = i2x%rAttr(kmeltw,n)/i2x%rAttr(kiFrac,n)
             else
                i2x%rAttr(kmelth,n) = 0.0_R8
                i2x%rAttr(kmeltw,n) = 0.0_R8
             end if
          end if

          !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
          i2x%rAttr(kt,n) = tfreeze(n) + i2x%rAttr(kiFrac,n)*(i2x%rAttr(kt,n)-tfreeze(n))

       end do

       ! compute atm/ice surface fluxes
       call dice_flux_atmice( &
            iMask              ,x2i%rAttr(kz,:)     ,x2i%rAttr(kua,:)    ,x2i%rAttr(kva,:)  , &
            x2i%rAttr(kptem,:) ,x2i%rAttr(kshum,:)  ,x2i%rAttr(kdens,:)  ,x2i%rAttr(ktbot,:), &
            i2x%rAttr(kt,:)    ,i2x%rAttr(ksen,:)   ,i2x%rAttr(klat,:)   ,i2x%rAttr(klwup,:), &
            i2x%rAttr(kevap,:) ,i2x%rAttr(ktauxa,:) ,i2x%rAttr(ktauya,:) ,i2x%rAttr(ktref,:), &
            i2x%rAttr(kqref,:) ,logunit )

       ! compute ice/oce surface fluxes (except melth & meltw, see above)
       do n=1,lsize
          if (iMask(n) == 0) then
             i2x%rAttr(kswpen,n) = spval
             i2x%rAttr(kmelth,n) = spval
             i2x%rAttr(kmeltw,n) = spval
             i2x%rAttr(ksalt ,n) = spval
             i2x%rAttr(ktauxo,n) = spval
             i2x%rAttr(ktauyo,n) = spval
             i2x%rAttr(kiFrac,n) = 0.0_R8
          else
             !--- penetrating short wave ---
             i2x%rAttr(kswpen,n) = max(0.0_R8, flux_swpf*i2x%rAttr(kswnet,n) ) ! must be non-negative

             !--- i/o surface stress ( = atm/ice stress) ---
             i2x%rAttr(ktauxo,n) = i2x%rAttr(ktauxa,n)
             i2x%rAttr(ktauyo,n) = i2x%rAttr(ktauya,n)

             !--- salt flux ---
             i2x%rAttr(ksalt ,n) = 0.0_R8
          end if

          !         !--- save ifrac for next timestep
          !         iFrac0(n) = i2x%rAttr(kiFrac,n)
       end do

       ! Compute outgoing aerosol fluxes
       do n = 1,lsize
          i2x%rAttr(kbcpho ,n) = x2i%rAttr(kbcphodry,n)
          i2x%rAttr(kbcphi ,n) = x2i%rAttr(kbcphidry,n) + x2i%rAttr(kbcphiwet,n)
          i2x%rAttr(kflxdst,n) = x2i%rAttr(kdstdry1,n) + x2i%rAttr(kdstwet1,n) &
                               + x2i%rAttr(kdstdry2,n) + x2i%rAttr(kdstwet2,n) &
                               + x2i%rAttr(kdstdry3,n) + x2i%rAttr(kdstwet3,n) &
                               + x2i%rAttr(kdstdry4,n) + x2i%rAttr(kdstwet4,n)
       end do

    end select

    !-------------------------------------------------
    ! optional per thickness category fields
    !-------------------------------------------------

    if (seq_flds_i2o_per_cat) then
       do n=1,lsize
          i2x%rAttr(kiFrac_01,n)       = i2x%rAttr(kiFrac,n)
          i2x%rAttr(kswpen_iFrac_01,n) = i2x%rAttr(kswpen,n) * i2x%rAttr(kiFrac,n)
       end do
    end if

    call t_stopf('dice_datamode')

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('dice_restart')
       call shr_cal_ymdtod2string(date_str, yy, mm, dd, currentTOD)
       write(rest_file,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(case_name), '.dice',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),currentYMD,currentTOD
       call shr_pcdf_readwrite('write',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file),mpicom,gsmap,clobber=.true.,rf1=water,rf1n='water')
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
       call shr_strdata_restWrite(trim(rest_file_strm),SDICE,mpicom,trim(case_name),'SDICE strdata')
       call shr_sys_flush(logunit)
       call t_stopf('dice_restart')
    endif

    call t_stopf('dice')

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    call t_startf('dice_run2')
    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if

    firstcall = .false.
    call t_stopf('dice_run2')

    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

  !===============================================================================
  subroutine dice_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN) , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN) , intent(in) :: master_task ! task number of master task
    integer(IN) , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dice_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dice_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dice_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('DICE_FINAL')

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(myModelName),': end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DICE_FINAL')

  end subroutine dice_comp_final
  !===============================================================================
end module dice_comp_mod
