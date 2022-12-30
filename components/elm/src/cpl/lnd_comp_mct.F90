module lnd_comp_mct
  
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Interface of the active land model component of CESM the ELM (E3SM Land Model)
  !  with the main E3SM driver. This is a thin interface taking E3SM driver information
  !  in MCT (Model Coupling Toolkit) format and converting it to use by ELM.
  !
  ! !uses:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_sys_mod      , only : shr_sys_flush
  use mct_mod          , only : mct_avect, mct_gsmap
  use decompmod        , only : bounds_type, ldecomp
  use lnd_import_export
  use iso_c_binding

#ifdef HAVE_MOAB
  use seq_comm_mct,       only: mlnid! id of moab land app
#endif
  !
  ! !public member functions:
  implicit none
  save
  private                     ! by default make data private
  !
  ! !public member functions:
  public :: lnd_init_mct      ! elm initialization
  public :: lnd_run_mct       ! elm run phase
  public :: lnd_final_mct     ! elm finalization/cleanup
  !
  ! !private member functions:
  private :: lnd_setgsmap_mct ! set the land model mct gs map
  private :: lnd_domain_mct   ! set the land model domain information

#ifdef HAVE_MOAB
  private :: init_land_moab   ! create moab mesh (cloud of points)
  private :: lnd_export_moab ! it should be part of lnd_import_export, but we will keep it here
  integer , private :: mblsize, totalmbls
  real (r8) , allocatable, private :: l2x_lm(:,:) ! for tags in MOAB
  logical :: sameg_al ! save it for export :)  
#endif
  !---------------------------------------------------------------------------

contains

  !====================================================================================

  subroutine lnd_init_mct( EClock, cdata_l, x2l_l, l2x_l, NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize land surface model and obtain relevant atmospheric model arrays
    ! back from (i.e. albedos, surface temperature and snow cover over land).
    !
    ! !USES:
    use abortutils       , only : endrun
    use shr_kind_mod     , only : SHR_KIND_CL
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, set_nextsw_cday
    use elm_initializeMod, only : initialize1, initialize2, initialize3
    use elm_instMod      , only : lnd2atm_vars, lnd2glc_vars
    use elm_varctl       , only : finidat,single_column, elm_varctl_set, iulog, noland
    use elm_varctl       , only : inst_index, inst_suffix, inst_name
    use elm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds
    use domainMod        , only : ldomain
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod     , only : shr_file_getLogUnit, shr_file_getLogLevel
    use shr_file_mod     , only : shr_file_getUnit, shr_file_setIO, shr_file_freeunit
    use shr_taskmap_mod  , only : shr_taskmap_write
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use seq_comm_mct     , only : info_taskmap_comp
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use seq_flds_mod     , only : seq_flds_x2l_fields, seq_flds_l2x_fields
    use spmdMod          , only : masterproc, npes, spmd_init
    use elm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use elm_cpl_indices  , only : elm_cpl_indices_set
    use perf_mod         , only : t_startf, t_stopf
    use mct_mod
    use ESMF
#ifdef HAVE_MOAB
    use iMOAB            , only : iMOAB_RegisterApplication
#endif
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock           ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_l          ! Input land-model driver data
    type(mct_aVect),            intent(inout) :: x2l_l, l2x_l     ! land model import and export states
    character(len=*), optional, intent(in)    :: NLFilename       ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    integer                          :: LNDID        ! Land identifier
    integer                          :: mpicom_lnd   ! MPI communicator
    type(mct_gsMap),         pointer :: GSMap_lnd    ! Land model MCT GS map
    type(mct_gGrid),         pointer :: dom_l        ! Land model domain
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer  :: lsz                                  ! size of attribute vector
    integer  :: g,i,j                                ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_elm                            ! elm time-step
    logical  :: exists                               ! true if file exists
    logical  :: verbose_taskmap_output               ! true then use verbose task-to-node mapping format
    logical  :: atm_aero                             ! Flag if aerosol data sent from atm model
    logical  :: atm_present                          ! Flag if atmosphere model present
    real(r8) :: scmlat                               ! single-column latitude
    real(r8) :: scmlon                               ! single-column longitude
    real(r8) :: nextsw_cday                          ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=8)           :: c_inst_index       ! instance number           
    character(len=8)           :: c_npes             ! number of pes
    integer :: nsrest                                ! elm restart type
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: nstep
    type(bounds_type) :: bounds                      ! bounds
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"

#ifdef HAVE_MOAB
    integer :: ierr, nsend
    character*32  appname
    logical :: samegrid_al !
    character(len=SHR_KIND_CL) :: atm_gnam          ! atm grid
    character(len=SHR_KIND_CL) :: lnd_gnam          ! lnd grid
    ! debugIuli
    integer   ::        debugGSMapFile, n
#endif
    !-----------------------------------------------------------------------

    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    ! Determine attriute vector indices

    call elm_cpl_indices_set()

    ! Initialize elm MPI communicator 

    call spmd_init( mpicom_lnd, LNDID )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_mct:start::',lbnum)
    endif
#endif                      

    inst_name   = seq_comm_name(LNDID)
    inst_index  = seq_comm_inst(LNDID)
    inst_suffix = seq_comm_suffix(LNDID)

    ! Initialize io log unit

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "ELM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Identify SMP nodes and process/SMP mapping for this instance
    ! (Assume that processor names are SMP node names on SMP clusters.)
    write(c_inst_index,'(i8)') inst_index

    if (info_taskmap_comp > 0) then

       if (info_taskmap_comp == 1) then
          verbose_taskmap_output = .false.
       else
          verbose_taskmap_output = .true.
       endif

       write(c_npes,'(i8)') npes

       if (masterproc) then
          write(iulog,'(/,3A)') &
             trim(adjustl(c_npes)), &
             ' pes participating in computation of ELM instance #', &
             trim(adjustl(c_inst_index))
          call shr_sys_flush(iulog)
       endif

       call t_startf("shr_taskmap_write")
       call shr_taskmap_write(iulog, mpicom_lnd,                    &
                              'LND #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output        )
       call t_stopf("shr_taskmap_write")

    endif

    ! Use infodata to set orbital values

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Consistency check on namelist filename	

    call control_setNL("lnd_in"//trim(inst_suffix))

    ! Initialize elm
    ! initialize1 reads namelist, grid and surface data (need this to initialize gsmap) 
    ! initialize2 performs rest of initialization	

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, case_name=caseid,    &
                              case_desc=ctitle, single_column=single_column,    &
                              scmlat=scmlat, scmlon=scmlon,                     &
                              brnch_retain_casename=brnch_retain_casename,      &
                              start_type=starttype, model_version=version,      &
                              hostname=hostname, username=username )
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod)
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call elm_varctl_set(caseid_in=caseid, ctitle_in=ctitle,                     &
                        brnch_retain_casename_in=brnch_retain_casename,         &
                        single_column_in=single_column, scmlat_in=scmlat,       &
                        scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                        hostname_in=hostname, username_in=username)

    ! Read namelist, grid and surface data

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland ) then
       call seq_infodata_PutData( infodata, lnd_present   =.false.)
       call seq_infodata_PutData( infodata, lnd_prognostic=.false.)
       return
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    call seq_infodata_GetData(infodata, atm_present=atm_present)
    call seq_infodata_GetData(infodata, atm_aero=atm_aero )
    !DMR 6/12/15 - remove this requirement (CPL_BPYASS mode uses SATM)
    if ( .not. atm_aero .and. atm_present )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to ELM' )
    end if

    ! Initialize elm gsMap, elm domain and elm attribute vectors

    call get_proc_bounds( bounds )

    call lnd_SetgsMap_mct( bounds, mpicom_lnd, LNDID, gsMap_lnd )
    lsz = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)
#ifdef HAVE_MOAB
    appname="LNDMB"//C_NULL_CHAR
    ! first land instance, should be 9
    ierr = iMOAB_RegisterApplication(appname, mpicom_lnd, LNDID, mlnid)
    if (ierr > 0 )  &
       call endrun('Error: cannot register moab app')
    if(masterproc) then
       write(iulog,*) " "
       write(iulog,*) "register MOAB app:", trim(appname), "  mlnid=", mlnid
       write(iulog,*) " "
    endif

#if 0
    if (masterproc) then
      debugGSMapFile = shr_file_getUnit()
      open( debugGSMapFile, file='LndGSmapC.txt')
      write(debugGSMapFile,*) gsMap_lnd%comp_id
      write(debugGSMapFile,*) gsMap_lnd%ngseg
      write(debugGSMapFile,*) gsMap_lnd%gsize
      do n=1,gsMap_lnd%ngseg
          write(debugGSMapFile,*) gsMap_lnd%start(n),gsMap_lnd%length(n),gsMap_lnd%pe_loc(n)
      end do
      close(debugGSMapFile)
      call shr_file_freeunit(debugGSMapFile)
    endif
#endif
!  endif HAVE_MOAB
#endif

    call lnd_domain_mct( bounds, lsz, gsMap_lnd, dom_l )
#ifdef HAVE_MOAB
!   find out samegrid_al or not; from infodata
    samegrid_al = .true.
    call seq_infodata_GetData(infodata         , &
                   atm_gnam=atm_gnam           , &
                   lnd_gnam=lnd_gnam           )
    if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
    call init_land_moab(bounds, samegrid_al)
    sameg_al = samegrid_al ! will use it for export too
#endif
    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsz)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsz)
    call mct_aVect_zero(l2x_l)

#ifdef HAVE_MOAB
    mblsize = lsz
    nsend = mct_avect_nRattr(l2x_l)
    totalmbls = mblsize * nsend ! size of the double array
    allocate (l2x_lm(lsz, nsend) )
#endif
    ! Finish initializing elm

    call initialize2()
    call initialize3()

    ! Check that elm internal dtime aligns with elm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_elm = get_step_size()
    if (masterproc) then
       write(iulog,*)'dtime_sync= ',dtime_sync,&
            ' dtime_elm= ',dtime_elm,' mod = ',mod(dtime_sync,dtime_elm)
    end if
    if (mod(dtime_sync,dtime_elm) /= 0) then
       write(iulog,*)'elm dtime ',dtime_elm,' and Eclock dtime ',&
            dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state 

    if (atm_present) then 
      call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, l2x_l%rattr)
#ifdef HAVE_MOAB
      call lnd_export_moab(bounds, lnd2atm_vars, lnd2glc_vars) ! it is private here
#endif
    endif

    ! Fill in infodata settings

    call seq_infodata_PutData(infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData(infodata, lnd_nx=ldomain%ni, lnd_ny=ldomain%nj)

    ! Get infodata info

    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )
    call set_nextsw_cday(nextsw_cday)

    if (.not. atm_present) then 
      !Calculate next radiation calendar day (since atm model did not run to set
      !this)
      !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
      nstep = get_nstep()
      nextsw_cday = mod((nstep/(86400._r8/dtime_elm))*1.0_r8,365._r8)+1._r8
      call set_nextsw_cday( nextsw_cday )
    end if

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine lnd_init_mct

  !====================================================================================

  subroutine lnd_run_mct(EClock, cdata_l, x2l_l, l2x_l)
    !
    ! !DESCRIPTION:
    ! Run elm model
    !
    ! !USES:
    use shr_kind_mod    ,  only : r8 => shr_kind_r8
    use elm_instMod     , only : lnd2atm_vars, atm2lnd_vars, lnd2glc_vars, glc2lnd_vars
    use elm_driver      ,  only : elm_drv
    use clm_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use clm_time_manager,  only : advance_timestep, set_nextsw_cday,update_rad_dtime
    use decompMod       ,  only : get_proc_bounds
    use abortutils      ,  only : endrun
    use elm_varctl      ,  only : iulog
    use elm_varorb      ,  only : eccen, obliqr, lambm0, mvelpp
    use shr_file_mod    ,  only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod    ,  only : shr_file_getLogUnit, shr_file_getLogLevel
    use seq_cdata_mod   ,  only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,  only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
    use seq_timemgr_mod ,  only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use seq_infodata_mod,  only : seq_infodata_type, seq_infodata_GetData
    use spmdMod         ,  only : masterproc, mpicom
    use perf_mod        ,  only : t_startf, t_stopf, t_barrierf
    use shr_orb_mod     ,  only : shr_orb_decl
    use mct_mod
    use ESMF
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    !
    ! !LOCAL VARIABLES:
    integer      :: ymd_sync             ! Sync date (YYYYMMDD)
    integer      :: yr_sync              ! Sync current year
    integer      :: mon_sync             ! Sync current month
    integer      :: day_sync             ! Sync current day
    integer      :: tod_sync             ! Sync current time of day (sec)
    integer      :: ymd                  ! ELM current date (YYYYMMDD)
    integer      :: yr                   ! ELM current year
    integer      :: mon                  ! ELM current month
    integer      :: day                  ! ELM current day
    integer      :: tod                  ! ELM current time of day (sec)
    integer      :: dtime                ! time step increment (sec)
    integer      :: nstep                ! time step index
    logical      :: rstwr_sync           ! .true. ==> write restart file before returning
    logical      :: rstwr                ! .true. ==> write restart file before returning
    logical      :: nlend_sync           ! Flag signaling last time-step
    logical      :: nlend                ! .true. ==> last time-step
    logical      :: dosend               ! true => send data back to driver
    logical      :: doalb                ! .true. ==> do albedo calculation on this time step
    real(r8)     :: nextsw_cday          ! calday from clock of next radiation computation
    real(r8)     :: caldayp1             ! elm calday plus dtime offset
    integer      :: shrlogunit,shrloglev ! old values for share log unit and log level
    integer      :: lbnum                ! input to memory diagnostic
    integer      :: g,i,lsz              ! counters
    real(r8)     :: calday               ! calendar day for nstep
    real(r8)     :: declin               ! solar declination angle in radians for nstep
    real(r8)     :: declinp1             ! solar declination angle in radians for nstep+1
    real(r8)     :: eccf                 ! earth orbit eccentricity factor
    real(r8)     :: recip                ! reciprical
    logical,save :: first_call = .true.  ! first call work
    logical      :: atm_present
    type(seq_infodata_type),pointer :: infodata             ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_l                ! Land model domain data
    type(bounds_type)               :: bounds               ! bounds
    character(len=32)               :: rdate                ! date char string for restart file names
    character(len=32), parameter    :: sub = "lnd_run_mct"
    !---------------------------------------------------------------------------

    ! Determine processor bounds

    call get_proc_bounds(bounds)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    call seq_cdata_setptrs(cdata_l, infodata=infodata, dom=dom_l)
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )
    dtime = get_step_size()

    call seq_infodata_GetData(infodata, atm_present=atm_present)
    if (.not. atm_present) then 
      !Calcualte next radiation calendar day (since atm model did not run to set this)
      !DMR:  NOTE this assumes a no-leap calendar and equal input/model timesteps
      nstep = get_nstep()
      nextsw_cday = mod((nstep/(86400._r8/dtime))*1.0_r8,365._r8)+1._r8 
      call set_nextsw_cday( nextsw_cday )
    end if
 
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    ! Map MCT to land data type
    ! Perform downscaling if appropriate

    
    ! Map to elm (only when state and/or fluxes need to be updated)

    call t_startf ('lc_lnd_import')
    call lnd_import( bounds, x2l_l%rattr, atm2lnd_vars, glc2lnd_vars)
    call t_stopf ('lc_lnd_import')

    ! Use infodata to set orbital values if updated mid-run

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Loop over time steps in coupling interval

    dosend = .false.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       nstep = get_nstep()
       caldayp1 = get_curr_calday(offset=dtime)
       if (nstep == 0) then
	        doalb = .false.
       else if (nstep == 1) then 
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8) 
       else
          doalb = (nextsw_cday >= -0.5_r8) 
       end if
       call update_rad_dtime(doalb)

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run elm 

       call t_barrierf('sync_elm_run1', mpicom)
       call t_startf ('elm_run')
       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')
       call elm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('elm_run')

       ! Create l2x_l export state - add river runoff input to l2x_l if appropriate

#ifndef CPL_BYPASS       
       call t_startf ('lc_lnd_export')
       call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, l2x_l%rattr)
#ifdef HAVE_MOAB
       call lnd_export_moab(bounds, lnd2atm_vars, lnd2glc_vars) ! it is private here
#endif
       call t_stopf ('lc_lnd_export')
#endif

       ! Advance elm time step
       
       call t_startf ('lc_elm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_elm2_adv_timestep')

    end do

    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' elm ymd=',ymd     ,'  elm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: ELM clock not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call  = .false.

  end subroutine lnd_run_mct

  !====================================================================================

  subroutine lnd_final_mct( EClock, cdata_l, x2l_l, l2x_l)
    !
    ! !DESCRIPTION:
    ! Finalize land surface model

    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
    use seq_timemgr_mod ,only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use mct_mod
    use esmf
    use elm_finalizeMod, only : final
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    !---------------------------------------------------------------------------

    ! fill this in
#ifdef HAVE_MOAB
    ! deallocate moab fields array
      deallocate (l2x_lm)
#endif
    call final()

  end subroutine lnd_final_mct

  !====================================================================================

  subroutine lnd_setgsmap_mct( bounds, mpicom_lnd, LNDID, gsMap_lnd )
    !
    ! !DESCRIPTION:
    ! Set the MCT GS map for the land model
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use domainMod    , only : ldomain
    use mct_mod      , only : mct_gsMap, mct_gsMap_init
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds     ! bounds
    integer           , intent(in)  :: mpicom_lnd ! MPI communicator for the elm land model
    integer           , intent(in)  :: LNDID      ! Land model identifyer number
    type(mct_gsMap)   , intent(out) :: gsMap_lnd  ! Resulting MCT GS map for the land model
    !
    ! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: i, j, n, gi            ! Indices
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code
    !---------------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    allocate(gindex(bounds%begg:bounds%endg),stat=ier)

    ! number the local grid

    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1
    gsize = ldomain%ni * ldomain%nj

    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

  !====================================================================================

  subroutine lnd_domain_mct( bounds, lsz, gsMap_l, dom_l )
    !
    ! !DESCRIPTION:
    ! Send the land model domain information to the coupler
    !
    ! !USES:
    use elm_varcon  , only: re
    use domainMod   , only: ldomain
    use spmdMod     , only: iam
    use mct_mod     , only: mct_gsMap, mct_gGrid, mct_gGrid_importIAttr
    use mct_mod     , only: mct_gGrid_importRAttr, mct_gGrid_init, mct_gsMap_orderedPoints
    use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other
    !
    ! !ARGUMENTS: 
    type(bounds_type), intent(in)  :: bounds  ! bounds
    integer        , intent(in)    :: lsz     ! land model domain data size
    type(mct_gsMap), intent(inout) :: gsMap_l ! Output land model MCT GS map
    type(mct_ggrid), intent(out)   :: dom_l   ! Output domain information for land model
    !
    ! Local Variables
    integer :: g,i,j              ! index
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !---------------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_l, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsz )
    !
    ! Allocate memory
    !
    allocate(data(lsz))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsz)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsz)
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsz)
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsz)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsz)

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"frac",data,lsz)

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct

#ifdef HAVE_MOAB
  subroutine init_land_moab(bounds, samegrid_al)
    use seq_flds_mod     , only :  seq_flds_l2x_fields, seq_flds_x2l_fields
    use shr_kind_mod     , only : CXX => SHR_KIND_CXX
    use spmdMod     , only: iam  ! rank on the land communicator
    use domainMod   , only: ldomain ! ldomain is coming from module, not even passed
    use elm_varcon  , only: re
    use shr_const_mod, only: SHR_CONST_PI
    use iMOAB        , only: iMOAB_CreateVertices, iMOAB_WriteMesh, &
    iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
    iMOAB_ResolveSharedEntities, iMOAB_CreateElements, iMOAB_MergeVertices, iMOAB_UpdateMeshInfo

    type(bounds_type) , intent(in)  :: bounds
    logical :: samegrid_al

    integer,allocatable :: gindex(:)  ! Number the local grid points; used for global ID
    integer lsz !  keep local size
    integer gsize ! global size, that we do not need, actually
    integer n
    ! local variables to fill in data
    integer, dimension(:), allocatable :: vgids
    !  retrieve everything we need from land domain mct_ldom
    ! number of vertices is the size of land domain
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    real(r8)   :: latv, lonv
    integer   dims, i, iv, ilat, ilon, igdx, ierr, tagindex
    integer tagtype, numco, ent_type, mbtype, block_ID
    character*100 outfile, wopts, localmeshfile
    character(CXX) ::  tagname ! hold all fields

    integer, allocatable :: moabconn(:) ! will have the connectivity in terms of local index in verts

    dims  =3 ! store as 3d mesh
    ! number the local grid
    lsz = bounds%endg - bounds%begg + 1

    allocate(vgids(lsz)) ! use it for global ids, for elements in full mesh or vertices in point cloud

    do n = 1, lsz
       vgids(n) = ldecomp%gdc2glo(bounds%begg+n-1) ! local to global !
    end do
    gsize = ldomain%ni * ldomain%nj ! size of the total grid
    ! if ldomain%nv > 3 , create mesh
    if (ldomain%nv .ge. 3 .and.  .not.samegrid_al) then
        ! number of vertices is nv * lsz !
        allocate(moab_vert_coords(lsz*dims*ldomain%nv))
        ! loop over ldomain
        allocate(moabconn(ldomain%nv * lsz))
        do n = bounds%begg, bounds%endg
            i = (n - bounds%begg) * ldomain%nv
            do iv = 1, ldomain%nv
               lonv = ldomain%lonv(n, iv) * SHR_CONST_PI/180.
               latv = ldomain%latv(n, iv) * SHR_CONST_PI/180.
               i = i + 1 ! iv-th vertex of cell n; i starts at 1 ! should we repeat previous if nan
               ! print *, i, n, ldomain%lonv(n, iv) , ldomain%latv(n, iv)
               moab_vert_coords(3*i-2)=COS(latv)*COS(lonv)
               moab_vert_coords(3*i-1)=COS(latv)*SIN(lonv)
               moab_vert_coords(3*i  )=SIN(latv)
               moabconn(i) = i!
            enddo
        enddo
        ierr = iMOAB_CreateVertices(mlnid, lsz * 3 * ldomain%nv, dims, moab_vert_coords)
        if (ierr > 0 )  &
            call endrun('Error: fail to create MOAB vertices in land model')


        mbtype = 2 ! triangle
        if (ldomain%nv .eq. 4) mbtype = 3 ! quad
        if (ldomain%nv .gt. 4) mbtype = 4 ! polygon
        block_ID = 100 !some value
        ierr = iMOAB_CreateElements( mlnid, lsz, mbtype, ldomain%nv, moabconn, block_ID );
        ! define some tags on cells now, not on vertices
        tagtype = 0  ! dense, integer
        numco = 1
        tagname='GLOBAL_ID'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to retrieve GLOBAL_ID tag ')

        ent_type = 1 ! element type
        ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsz , ent_type, vgids)
        if (ierr > 0 )  &
          call endrun('Error: fail to set GLOBAL_ID tag ')

!        ierr = iMOAB_ResolveSharedEntities( mlnid, lsz, vgids );
!        if (ierr > 0 )  &
!          call endrun('Error: fail to resolve shared entities')

!        !there are no shared entities, but we will set a special partition tag, in order to see the
!        ! partitions ; it will be visible with a Pseudocolor plot in VisIt
!        tagname='partition'//C_NULL_CHAR
!        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
!        if (ierr > 0 )  &
!          call endrun('Error: fail to create new partition tag ')
!
!        vgids = iam
!        ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsz , ent_type, vgids)
!        if (ierr > 0 )  &
!          call endrun('Error: fail to set partition tag ')

        ! use moab_vert_coords as a data holder for a frac tag and area tag that we will create
        !   on the vertices; do not allocate other data array
        tagname='frac'//C_NULL_CHAR
        tagtype = 1 ! dense, double
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create frac tag ')

        do i = 1, lsz
          n = i-1 + bounds%begg
          moab_vert_coords(i) = ldomain%frac(n)
        enddo
        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords)
        if (ierr > 0 )  &
          call endrun('Error: fail to set frac tag ')

        tagname='area'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create area tag ')
        do i = 1, lsz
          n = i-1 + bounds%begg
          moab_vert_coords(i) = ldomain%area(n)/(re*re) ! use the same doubles for second tag :)
        enddo

        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        if (ierr > 0 )  &
          call endrun('Error: fail to set area tag ')

        ! aream needed in cime_init for now.
        tagname='aream'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create aream tag ')
        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        if (ierr > 0 )  &
          call endrun('Error: fail to set aream tag ')

        deallocate(moabconn)
        ! use merge vertices new imoab method to fix cells
        deallocate(vgids) ! use it for global ids, for elements in full mesh or vertices in point cloud
        allocate(vgids(lsz*ldomain%nv)) ! 
        do n = 1, lsz
          do i=1,ldomain%nv
            vgids( (n-1)*ldomain%nv+i ) = (ldecomp%gdc2glo(bounds%begg+n-1)-1)*ldomain%nv+i ! local to global !
          end do
        end do
        ent_type = 0 ! vertices now
        tagname = 'GLOBAL_ID'//C_NULL_CHAR
        ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsz , ent_type, vgids )
        if (ierr > 0 )  &
          call endrun('Error: fail to set global ID tag on vertices in land mesh ')
        ierr = iMOAB_UpdateMeshInfo( mlnid )
        if (ierr > 0 )  &
          call endrun('Error: fail to update mesh info ')
        !ierr = iMOAB_MergeVertices(mlnid)
        !if (ierr > 0 )  &
        !  call endrun('Error: fail to fix vertices in land mesh ')

    else ! old point cloud mesh
        allocate(moab_vert_coords(lsz*dims))
        do i = 1, lsz
          n = i-1 + bounds%begg
          lonv = ldomain%lonc(n) *SHR_CONST_PI/180.
          latv = ldomain%latc(n) *SHR_CONST_PI/180.
          moab_vert_coords(3*i-2)=COS(latv)*COS(lonv)
          moab_vert_coords(3*i-1)=COS(latv)*SIN(lonv)
          moab_vert_coords(3*i  )=SIN(latv)
        enddo
        ierr = iMOAB_CreateVertices(mlnid, lsz*3, dims, moab_vert_coords)
        if (ierr > 0 )  &
          call endrun('Error: fail to create MOAB vertices in land model')

        tagtype = 0  ! dense, integer
        numco = 1
        tagname='GLOBAL_ID'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to retrieve GLOBAL_ID tag ')

        ent_type = 0 ! vertex type
        ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsz , ent_type, vgids)
        if (ierr > 0 )  &
          call endrun('Error: fail to set GLOBAL_ID tag ')

        ierr = iMOAB_ResolveSharedEntities( mlnid, lsz, vgids );
        if (ierr > 0 )  &
          call endrun('Error: fail to resolve shared entities')

        !there are no shared entities, but we will set a special partition tag, in order to see the
        ! partitions ; it will be visible with a Pseudocolor plot in VisIt
        tagname='partition'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create new partition tag ')

        vgids = iam
        ierr = iMOAB_SetIntTagStorage ( mlnid, tagname, lsz , ent_type, vgids)
        if (ierr > 0 )  &
          call endrun('Error: fail to set partition tag ')

        ! use moab_vert_coords as a data holder for a frac tag and area tag that we will create
        !   on the vertices; do not allocate other data array
        tagname='frac'//C_NULL_CHAR
        tagtype = 1 ! dense, double
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create frac tag ')

        do i = 1, lsz
          n = i-1 + bounds%begg
          moab_vert_coords(i) = ldomain%frac(n)
        enddo
        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords)
        if (ierr > 0 )  &
          call endrun('Error: fail to set frac tag ')

        tagname='area'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create area tag ')
        do i = 1, lsz
          n = i-1 + bounds%begg
          moab_vert_coords(i) = ldomain%area(n)/(re*re) ! use the same doubles for second tag :)
        enddo

        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        if (ierr > 0 )  &
          call endrun('Error: fail to set area tag ')

        ! aream needed in cime_init for now.
        tagname='aream'//C_NULL_CHAR
        ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
        if (ierr > 0 )  &
          call endrun('Error: fail to create aream tag ')
        ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        if (ierr > 0 )  &
          call endrun('Error: fail to set aream tag ')
    endif
    deallocate(moab_vert_coords)
    deallocate(vgids)
#ifdef MOABDEBUG
    !     write out the mesh file to disk, in parallel
    outfile = 'wholeLnd.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mlnid, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the land mesh file')
#endif
  ! define all tags from seq_flds_l2x_fields
  ! define tags according to the seq_flds_l2x_fields 
    tagtype = 1  ! dense, double
    numco = 1 !  one value per cell / entity
    tagname = trim(seq_flds_l2x_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
    if ( ierr > 0) then
        call endrun('Error: fail to define seq_flds_l2x_fields for land moab mesh')
    endif
    tagname = trim(seq_flds_x2l_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
    if ( ierr > 0) then
        call endrun('Error: fail to define seq_flds_x2l_fields for land moab mesh')
    endif

  end subroutine init_land_moab

  subroutine lnd_export_moab( bounds, lnd2atm_vars, lnd2glc_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the elm model to the moab coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use elm_varctl         , only : iulog, create_glacier_mec_landunit
    use clm_time_manager   , only : get_nstep, get_step_size  
    use domainMod          , only : ldomain
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use iMOAB,  only       : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    use seq_flds_mod, only : seq_flds_l2x_fields
    use seq_comm_mct, only : num_moab_exports
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_vars ! clm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_vars ! clm land to atmosphere exchange data type
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i   ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step   
    integer  :: num   ! counter
    character(len=*), parameter :: sub = 'lnd_export_moab'

    integer :: ent_type, ierr
    character(len=100) :: outfile, wopts, lnum
   character(len=400) :: tagname
    !---------------------------------------------------------------------------

    ! cesm sign convention is that fluxes are positive downward

    l2x_lm(:,:) = 0.0_r8

    do g = bounds%begg,bounds%endg
       i = 1 + (g-bounds%begg)
       l2x_lm(i,index_l2x_Sl_t)        =  lnd2atm_vars%t_rad_grc(g)
       l2x_lm(i,index_l2x_Sl_snowh)    =  lnd2atm_vars%h2osno_grc(g)
       l2x_lm(i,index_l2x_Sl_avsdr)    =  lnd2atm_vars%albd_grc(g,1)
       l2x_lm(i,index_l2x_Sl_anidr)    =  lnd2atm_vars%albd_grc(g,2)
       l2x_lm(i,index_l2x_Sl_avsdf)    =  lnd2atm_vars%albi_grc(g,1)
       l2x_lm(i,index_l2x_Sl_anidf)    =  lnd2atm_vars%albi_grc(g,2)
       l2x_lm(i,index_l2x_Sl_tref)     =  lnd2atm_vars%t_ref2m_grc(g)
       l2x_lm(i,index_l2x_Sl_qref)     =  lnd2atm_vars%q_ref2m_grc(g)
       l2x_lm(i,index_l2x_Sl_u10)      =  lnd2atm_vars%u_ref10m_grc(g)
       l2x_lm(i,index_l2x_Fall_taux)   = -lnd2atm_vars%taux_grc(g)
       l2x_lm(i,index_l2x_Fall_tauy)   = -lnd2atm_vars%tauy_grc(g)
       l2x_lm(i,index_l2x_Fall_lat)    = -lnd2atm_vars%eflx_lh_tot_grc(g)
       l2x_lm(i,index_l2x_Fall_sen)    = -lnd2atm_vars%eflx_sh_tot_grc(g)
       l2x_lm(i,index_l2x_Fall_lwup)   = -lnd2atm_vars%eflx_lwrad_out_grc(g)
       l2x_lm(i,index_l2x_Fall_evap)   = -lnd2atm_vars%qflx_evap_tot_grc(g)
       l2x_lm(i,index_l2x_Fall_swnet)  =  lnd2atm_vars%fsa_grc(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x_lm(i,index_l2x_Fall_fco2_lnd) = -lnd2atm_vars%nee_grc(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x_lm(i,index_l2x_Sl_ram1)     =  lnd2atm_vars%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x_lm(i,index_l2x_Sl_fv)       =  lnd2atm_vars%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x_lm(i,index_l2x_Sl_soilw)    =  lnd2atm_vars%h2osoi_vol_grc(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x_lm(i,index_l2x_Fall_flxdst1)= -lnd2atm_vars%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x_lm(i,index_l2x_Fall_flxdst2)= -lnd2atm_vars%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x_lm(i,index_l2x_Fall_flxdst3)= -lnd2atm_vars%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x_lm(i,index_l2x_Fall_flxdst4)= -lnd2atm_vars%flxdst_grc(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x_lm(i,index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1) = &
               lnd2atm_vars%ddvel_grc(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x_lm(i,index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1) = &
               -lnd2atm_vars%flxvoc_grc(g,:shr_megan_mechcomps_n)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x_lm(i,index_l2x_Fall_methane) = -lnd2atm_vars%flux_ch4_grc(g) 
       endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x_lm(i,index_l2x_Flrl_rofi) = lnd2atm_vars%qflx_rofice_grc(g)
       l2x_lm(i,index_l2x_Flrl_rofsur) = lnd2atm_vars%qflx_rofliq_qsur_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsurp_grc(g)   !  surface ponding
       l2x_lm(i,index_l2x_Flrl_rofsub) = lnd2atm_vars%qflx_rofliq_qsub_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsubp_grc(g)   !  perched drainiage
       l2x_lm(i,index_l2x_Flrl_rofgwl) = lnd2atm_vars%qflx_rofliq_qgwl_grc(g)
  
       l2x_lm(i,index_l2x_Flrl_demand) =  lnd2atm_vars%qflx_irr_demand_grc(g)   ! needs to be filled in
       if (l2x_lm(i,index_l2x_Flrl_demand) > 0.0_r8) then
           write(iulog,*)'lnd2atm_vars%qflx_irr_demand_grc is',lnd2atm_vars%qflx_irr_demand_grc(g)
           write(iulog,*)'l2x_lm(i,index_l2x_Flrl_demand) is',l2x_lm(i,index_l2x_Flrl_demand)
           call endrun( sub//' ERROR: demand must be <= 0.')
       endif
       l2x_lm(i,index_l2x_Flrl_Tqsur)  = lnd2atm_vars%Tqsur_grc(g)
       l2x_lm(i,index_l2x_Flrl_Tqsub)  = lnd2atm_vars%Tqsub_grc(g)
       l2x_lm(i,index_l2x_coszen_str) = lnd2atm_vars%coszen_str(g)
       ! glc coupling

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             l2x_lm(i,index_l2x_Sl_tsrf(num))   = lnd2glc_vars%tsrf_grc(g,num)
             l2x_lm(i,index_l2x_Sl_topo(num))   = lnd2glc_vars%topo_grc(g,num)
             l2x_lm(i,index_l2x_Flgl_qice(num)) = lnd2glc_vars%qice_grc(g,num)
          end do
       end if

    end do
    tagname=trim(seq_flds_l2x_fields)//C_NULL_CHAR
    if (sameg_al) then
      ent_type = 0 ! vertices, cells only if sameg_al false
    else
      ent_type = 1
    endif
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, totalmbls , ent_type, l2x_lm(1,1) )
    if (ierr > 0 )  &
       call shr_sys_abort( sub//' Error: fail to set moab '// trim(seq_flds_l2x_fields) )
 
#ifdef MOABDEBUG
       write(lnum,"(I0.2)")num_moab_exports
       outfile = 'lnd_export_'//trim(lnum)//'.h5m'//C_NULL_CHAR
       wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
       ierr = iMOAB_WriteMesh(mlnid, outfile, wopts)
       if (ierr > 0 )  &
         call shr_sys_abort( sub//' fail to write the land mesh file with data')
#endif

  end subroutine lnd_export_moab
! endif for ifdef HAVE_MOAB
#endif

end module lnd_comp_mct
