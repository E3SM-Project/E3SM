#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif

module dice_comp_mod

  ! !USES:
  use mct_mod
  use perf_mod
  use shr_pcdf_mod
  use shr_const_mod
  use shr_sys_mod     , only: shr_sys_flush
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod    , only: shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod     , only: shr_mpi_bcast
  use shr_frz_mod     , only: shr_frz_freezetemp
  use shr_cal_mod     , only: shr_cal_datetod2string
  use shr_strdata_mod , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod  , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod  , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV

  use dice_shr_mod    , only: datamode       ! namelist input
  use dice_shr_mod    , only: decomp         ! namelist input
  use dice_shr_mod    , only: rest_file      ! namelist input
  use dice_shr_mod    , only: rest_file_strm ! namelist input
  use dice_shr_mod    , only: flux_swpf      ! namelist input -short-wave penatration factor
  use dice_shr_mod    , only: flux_Qmin      ! namelist input -bound on melt rate
  use dice_shr_mod    , only: flux_Qacc      ! namelist input -activates water accumulation/melt wrt Q
  use dice_shr_mod    , only: flux_Qacc0     ! namelist input -initial water accumulation value
  use dice_shr_mod    , only: nullstr
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
  public :: dice_comp_debug

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  logical                    :: firstcall = .true. ! first call logical
  integer(IN)                :: dbug = 2           ! debug level (higher is more)
  character(len=*),parameter :: rpfile = 'rpointer.ice'

  real(R8),parameter  :: pi     = shr_const_pi      ! pi
  real(R8),parameter  :: spval  = shr_const_spval   ! flags invalid data
  real(R8),parameter  :: tFrz   = shr_const_tkfrz   ! temp of freezing
  real(R8),parameter  :: latice = shr_const_latice  ! latent heat of fusion
  real(R8),parameter  :: waterMax = 1000.0_R8       ! wrt iFrac comp & frazil ice (kg/m^2)

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

  integer(IN) :: km
  integer(IN) :: kswvdr,kswndr,kswvdf,kswndf,kq,kz,kua,kva,kptem,kshum,kdens,ktbot
  integer(IN) :: kiFrac,kt,kavsdr,kanidr,kavsdf,kanidf,kswnet,kmelth,kmeltw
  integer(IN) :: ksen,klat,klwup,kevap,ktauxa,ktauya,ktref,kqref,kswpen,ktauxo,ktauyo,ksalt
  integer(IN) :: ksalinity
  integer(IN) :: kbcpho, kbcphi, kflxdst
  integer(IN) :: kbcphidry, kbcphodry, kbcphiwet, kocphidry, kocphodry, kocphiwet
  integer(IN) :: kdstdry1, kdstdry2, kdstdry3, kdstdry4, kdstwet1, kdstwet2, kdstwet3, kdstwet4
  integer(IN) :: kiFrac_01,kswpen_iFrac_01 ! optional per thickness category fields

  type(mct_rearr)       :: rearr
  integer(IN) , pointer :: imask(:)
  real(R8)    , pointer :: yc(:)
  real(R8)    , pointer :: water(:)
  real(R8)    , pointer :: tfreeze(:)
  !real(R8)   , pointer :: ifrac0(:)

!===============================================================================
contains
!===============================================================================

  subroutine dice_comp_init(x2i, i2x, &
       flds_x2i_fields, flds_i2x_fields, flds_i2o_per_cat, &
       SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scmlat, scmlon, calendar)

    ! !DESCRIPTION: initialize dice model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2i, i2x         ! input/output attribute vectors
    character(len=*)       , intent(in)    :: flds_x2i_fields  ! fields from mediator
    character(len=*)       , intent(in)    :: flds_i2x_fields  ! fields to mediator
    logical                , intent(in)    :: flds_i2o_per_cat ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE            ! dice shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap            ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid            ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom           ! mpi communicator
    integer(IN)            , intent(in)    :: compid           ! mct comp id
    integer(IN)            , intent(in)    :: my_task          ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task      ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix      ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name        ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit          ! logging unit number
    logical                , intent(in)    :: read_restart     ! start from restart
    logical                , intent(in)    :: scmMode          ! single column mode
    real(R8)               , intent(in)    :: scmLat           ! single column lat
    real(R8)               , intent(in)    :: scmLon           ! single column lon
    character(len=*)       , intent(in)    :: calendar         ! calendar type

    !--- local variables ---
    integer(IN)   :: n,k            ! generic counters
    integer(IN)   :: ierr           ! error code
    integer(IN)   :: lsize          ! local size
    integer(IN)   :: kfld           ! field reference
    logical       :: exists,exists1 ! file existance logical
    integer(IN)   :: nu             ! unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter :: F01   = "('(dice_comp_init) ',a,2f10.4)"
    character(*), parameter :: subName = "(dice_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DICE_INIT')

    call shr_strdata_pioinit(SDICE, compid)

    !----------------------------------------------------------------------------
    ! Initialize SDICE
    !----------------------------------------------------------------------------

    call t_startf('dice_strdata_init')

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDICE%gsmap and SDICE%ggrid. DICE%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    if (scmmode) then
       if (my_task == master_task) then
          write(logunit,F01) ' scm lon lat = ',scmlon,scmlat
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

    call mct_aVect_init(i2x, rList=flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x)

    km     = mct_aVect_indexRA(i2x,'Si_imask')
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

    if (flds_i2o_per_cat) then
       kiFrac_01       = mct_aVect_indexRA(i2x,'Si_ifrac_01')
       kswpen_iFrac_01 = mct_aVect_indexRA(i2x,'PFioi_swpen_ifrac_01')
    end if

    call mct_aVect_init(x2i, rList=flds_x2i_fields, lsize=lsize)
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

    allocate(imask(lsize))
    allocate(yc(lsize))
    allocate(water(lsize))
    allocate(tfreeze(lsize))
    ! allocate(iFrac0(lsize))

    ! Note that the module array, imask, does not change after initialization
    kfld = mct_aVect_indexRA(ggrid%data,'mask')
    imask(:) = nint(ggrid%data%rAttr(kfld,:))
    kfld = mct_aVect_indexRA(ggrid%data,'lat')
    yc(:) = ggrid%data%rAttr(kfld,:)

    if (km /= 0) then
       i2x%rAttr(km, :) = imask(:)
    end if

    call t_stopf('dice_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       exists = .false.
       exists1 = .false.
       if (trim(rest_file)      == trim(nullstr) .and. &
           trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer = ',trim(rpfile)
             call shr_sys_flush(logunit)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (exists) then
                nu = shr_file_getUnit()
                open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
                read(nu,'(a)') rest_file
                read(nu,'(a)') rest_file_strm
                close(nu)
                call shr_file_freeUnit(nu)
                inquire(file=trim(rest_file_strm),exist=exists)
                inquire(file=trim(rest_file),exist=exists1)
             endif
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
       call shr_mpi_bcast(exists1,mpicom,'exists1')

       if (exists1) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
          call shr_pcdf_readwrite('read',SDICE%pio_subsystem, SDICE%io_type, &
               trim(rest_file),mpicom,gsmap=gsmap,rf1=water,rf1n='water',io_format=SDICE%io_format)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
       endif

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

    call t_stopf('DICE_INIT')

  end subroutine dice_comp_init

  !===============================================================================

  subroutine dice_comp_run(x2i, i2x, flds_i2o_per_cat, &
       SDICE, gsmap, ggrid, mpicom, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       calendar, modeldt, target_ymd, target_tod, cosArg, avifld, avofld, case_name )

    ! !DESCRIPTION: run method for dice model

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_aVect)        , intent(inout) :: x2i
    type(mct_aVect)        , intent(inout) :: i2x
    logical                , intent(in)    :: flds_i2o_per_cat     ! .true. if select per ice thickness fields from ice
    type(shr_strdata_type) , intent(inout) :: SDICE
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer(IN)            , intent(in)    :: mpicom               ! mpi communicator
    integer(IN)            , intent(in)    :: my_task              ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task          ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix          ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit              ! logging unit number
    logical                , intent(in)    :: read_restart         ! start from restart
    logical                , intent(in)    :: write_restart        ! restart now
    character(len=*)       , intent(in)    :: calendar
    integer(IN)            , intent(in)    :: modeldt
    integer(IN)            , intent(in)    :: target_ymd
    integer(IN)            , intent(in)    :: target_tod
    real(R8)               , intent(in)    :: cosarg               ! for setting ice temp pattern
    character(len=*)       , intent(in)    :: avifld(:)
    character(len=*)       , intent(in)    :: avofld(:)
    character(CL)          , intent(in), optional :: case_name     ! case name

    !--- local ---
    integer(IN)   :: n                 ! indices
    integer(IN)   :: lsize             ! size of attr vect
    real(R8)      :: dt                ! timestep
    integer(IN)   :: nu                ! unit number
    real(R8)      :: qmeltall          ! q that would melt all accumulated water
    character(len=18) :: date_str

    character(*), parameter :: F00   = "('(dice_comp_run) ',8a)"
    character(*), parameter :: F04   = "('(dice_comp_run) ',a,2i8,'s')"
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    !--------------------
    ! ADVANCE ICE
    !--------------------

    call t_startf('DICE_RUN')
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice')

    dt = modeldt * 1.0_r8

    !--- copy all stream fields to i2x as default (avifld in streams -> avofld in i2x)

    if (trim(datamode) /= 'NULL') then
       call t_startf('dice_strdata_advance')
       call shr_strdata_advance(SDICE,target_ymd,target_tod,mpicom,'dice')
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
          if (km /= 0) then
             i2x%rAttr(km, n) = imask(n)
          end if
          ! !--- save ifrac for next timestep
          ! iFrac0(n) = i2x%rAttr(kiFrac,n)
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

    if (flds_i2o_per_cat) then
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
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
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
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),target_ymd,target_tod
       call shr_pcdf_readwrite('write',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file),mpicom,gsmap,clobber=.true.,rf1=water,rf1n='water')
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm),SDICE,mpicom,trim(case_name),'SDICE strdata')
       call shr_sys_flush(logunit)
       call t_stopf('dice_restart')
    endif

    call t_stopf('dice')

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,F04) 'dice : model date ', target_ymd,target_tod
       call shr_sys_flush(logunit)
    end if

    firstcall = .false.

    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

  !===============================================================================

  subroutine dice_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for dice model

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
       write(logunit,F00) 'dice: end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DICE_FINAL')

  end subroutine dice_comp_final

  !===============================================================================

  subroutine dice_comp_debug(my_task, master_task, logunit, target_ymd, target_tod, x2i, i2x)
    integer         , intent(in) :: my_task
    integer         , intent(in) :: master_task
    integer         , intent(in) :: logunit 
    integer         , intent(in) :: target_ymd
    integer         , intent(in) :: target_tod
    type(mct_aVect) , intent(in) :: x2i
    type(mct_aVect) , intent(in) :: i2x

    integer :: n
    character(*), parameter :: F01   = "('(dice_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"

    !----------------------------------------------------------
    ! Debug output
    !----------------------------------------------------------

    if (dbug > 1 .and. my_task == master_task) then
       do n = 1, mct_aVect_lsize(x2i)
          write(logunit,F01)'import: ymd,tod,n,Faxa_swvdr    = ', target_ymd, target_tod, n, x2i%rattr(kswvdr,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_swndr    = ', target_ymd, target_tod, n, x2i%rattr(kswndr,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_swvdf    = ', target_ymd, target_tod, n, x2i%rattr(kswvdf,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_swndf    = ', target_ymd, target_tod, n, x2i%rattr(kswndf,n)
          write(logunit,F01)'import: ymd,tod,n,Fioo_q        = ', target_ymd, target_tod, n, x2i%rattr(kq,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_z          = ', target_ymd, target_tod, n, x2i%rattr(kz,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_u          = ', target_ymd, target_tod, n, x2i%rattr(kua,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_v          = ', target_ymd, target_tod, n, x2i%rattr(kva,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_ptem       = ', target_ymd, target_tod, n, x2i%rattr(kptem,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_shum       = ', target_ymd, target_tod, n, x2i%rattr(kshum,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_dens       = ', target_ymd, target_tod, n, x2i%rattr(kdens,n)
          write(logunit,F01)'import: ymd,tod,n,Sa_tbot       = ', target_ymd, target_tod, n, x2i%rattr(ktbot,n)
          write(logunit,F01)'import: ymd,tod,n,So_s          = ', target_ymd, target_tod, n, x2i%rattr(ksalinity,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_bcphidry = ', target_ymd, target_tod, n, x2i%rattr(kbcphidry,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_bcphodry = ', target_ymd, target_tod, n, x2i%rattr(kbcphodry,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_bcphiwet = ', target_ymd, target_tod, n, x2i%rattr(kbcphiwet,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_ocphidry = ', target_ymd, target_tod, n, x2i%rattr(kocphidry,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_ocphodry = ', target_ymd, target_tod, n, x2i%rattr(kocphodry,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_ocphiwet = ', target_ymd, target_tod, n, x2i%rattr(kocphiwet,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstdry1  = ', target_ymd, target_tod, n, x2i%rattr(kdstdry1,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstdry2  = ', target_ymd, target_tod, n, x2i%rattr(kdstdry2,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstdry3  = ', target_ymd, target_tod, n, x2i%rattr(kdstdry3,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstdry4  = ', target_ymd, target_tod, n, x2i%rattr(kdstdry4,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstwet1  = ', target_ymd, target_tod, n, x2i%rattr(kdstwet1,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstwet2  = ', target_ymd, target_tod, n, x2i%rattr(kdstwet2,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstwet3  = ', target_ymd, target_tod, n, x2i%rattr(kdstwet3,n)
          write(logunit,F01)'import: ymd,tod,n,Faxa_dstwet4  = ', target_ymd, target_tod, n, x2i%rattr(kdstwet4,n)

          write(logunit,F01)'export: ymd,tod,n,Si_ifrac      = ', target_ymd, target_tod, n, i2x%rattr(kiFrac ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_t          = ', target_ymd, target_tod, n, i2x%rattr(kt     ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_tref       = ', target_ymd, target_tod, n, i2x%rattr(ktref  ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_qref       = ', target_ymd, target_tod, n, i2x%rattr(kqref  ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_avsdr      = ', target_ymd, target_tod, n, i2x%rattr(kavsdr ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_anidr      = ', target_ymd, target_tod, n, i2x%rattr(kanidr ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_avsdf      = ', target_ymd, target_tod, n, i2x%rattr(kavsdf ,n)
          write(logunit,F01)'export: ymd,tod,n,Si_anidf      = ', target_ymd, target_tod, n, i2x%rattr(kanidf ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_swnet    = ', target_ymd, target_tod, n, i2x%rattr(kswnet ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_sen      = ', target_ymd, target_tod, n, i2x%rattr(ksen   ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_lat      = ', target_ymd, target_tod, n, i2x%rattr(klat   ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_lwup     = ', target_ymd, target_tod, n, i2x%rattr(klwup  ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_evap     = ', target_ymd, target_tod, n, i2x%rattr(kevap  ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_taux     = ', target_ymd, target_tod, n, i2x%rattr(ktauxa ,n)
          write(logunit,F01)'export: ymd,tod,n,Faii_tauy     = ', target_ymd, target_tod, n, i2x%rattr(ktauya ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_melth    = ', target_ymd, target_tod, n, i2x%rattr(kmelth ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_meltw    = ', target_ymd, target_tod, n, i2x%rattr(kmeltw ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_swpen    = ', target_ymd, target_tod, n, i2x%rattr(kswpen ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_taux     = ', target_ymd, target_tod, n, i2x%rattr(ktauxo ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_tauy     = ', target_ymd, target_tod, n, i2x%rattr(ktauyo ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_salt     = ', target_ymd, target_tod, n, i2x%rattr(ksalt  ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_bcpho    = ', target_ymd, target_tod, n, i2x%rattr(kbcpho ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_bcphi    = ', target_ymd, target_tod, n, i2x%rattr(kbcphi ,n)
          write(logunit,F01)'export: ymd,tod,n,Fioi_flxdst   = ', target_ymd, target_tod, n, i2x%rattr(kflxdst,n)
       end do
    end if
  end subroutine dice_comp_debug

end module dice_comp_mod
