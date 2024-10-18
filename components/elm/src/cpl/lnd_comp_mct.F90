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
  use elm_cpl_indices
  use esmf, only: ESMF_clock

#ifdef HAVE_MOAB
  use seq_comm_mct,       only: mlnid! id of moab land app
  use seq_comm_mct,       only: num_moab_exports
#ifdef MOABCOMP 
  use seq_comm_mct , only: seq_comm_compare_mb_mct
#endif
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
  private :: init_moab_land   ! create moab mesh (cloud of points)
  private :: lnd_export_moab ! it could be part of lnd_import_export, but we will keep it here
  private :: lnd_import_moab ! it could be part of lnd_import_export, but we will keep it here
  integer , private :: mblsize, totalmbls
  real (r8) , allocatable, private :: l2x_lm(:,:) ! for tags to be set in MOAB

  integer :: nrecv, totalmblsimp
  real (r8) , allocatable, private :: x2l_lm(:,:) ! for tags from MOAB

  integer  :: mpicom_lnd_moab ! used also for mpi-reducing the difference between moab tags and mct avs
  integer :: rank2

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
    use elm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, set_nextsw_cday
    use elm_initializeMod, only : initialize1, initialize2, initialize3
    use elm_instMod      , only : lnd2atm_vars, lnd2glc_vars
    use elm_instance     , only : elm_instance_init
    use elm_varctl       , only : finidat,single_column, elm_varctl_set, iulog, noland, fatmlndfrc, &
                                  scm_multcols, scm_nx, scm_ny
    use elm_varctl       , only : inst_index, inst_suffix, inst_name, precip_downscaling_method
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
    use seq_flds_mod     , only : seq_flds_x2l_fields, seq_flds_l2x_fields, lnd_rof_two_way
    use spmdMod          , only : masterproc, npes, spmd_init
    use elm_varctl       , only : nsrStartup, nsrContinue, nsrBranch, use_lnd_rof_two_way
    use elm_cpl_indices  , only : elm_cpl_indices_set
    use perf_mod         , only : t_startf, t_stopf
    use mct_mod
    use ESMF

    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock           ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_l          ! Input land-model driver data
    type(mct_aVect),            intent(inout) :: x2l_l, l2x_l     ! land model import and export states
    character(len=*), optional, intent(in)    :: NLFilename       ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    integer                          :: LNDID        ! Land identifyer
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
    integer :: ierr, nsend,n
    character(len=SHR_KIND_CL) :: atm_gnam          ! atm grid
    character(len=SHR_KIND_CL) :: lnd_gnam          ! lnd grid
#endif
    !-----------------------------------------------------------------------

    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    ! Set and save LNDID for easy access by other modules

    call elm_instance_init( LNDID )

    ! Determine attriute vector indices
#ifdef HAVE_MOAB
    mpicom_lnd_moab = mpicom_lnd ! just store it now, for later use
    call shr_mpi_commrank( mpicom_lnd_moab, rank2 ) ! this will be used for differences between mct and moab tags
#endif 

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
                              scm_multcols=scm_multcols,scm_nx=scm_nx,scm_ny=scm_ny,    &
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

    ! If SCM domain mode, force single_column flag to be false for this
    !  block of code, as special treatment is needed
    if (scm_multcols) single_column = .false.
    call elm_varctl_set(caseid_in=caseid, ctitle_in=ctitle,                      &
                        brnch_retain_casename_in=brnch_retain_casename,          &
                        scm_multcols_in=scm_multcols, scm_nx_in=scm_nx, scm_ny_in=scm_ny,&
                        single_column_in=single_column, scmlat_in=scmlat,        &
                        scmlon_in=scmlon, nsrest_in=nsrest, version_in=version,  &
                        hostname_in=hostname, username_in=username)

    use_lnd_rof_two_way = lnd_rof_two_way
    
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

    call lnd_domain_mct( bounds, lsz, gsMap_lnd, dom_l )
#ifdef HAVE_MOAB
    call init_moab_land(bounds, LNDID)
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

    nrecv = mct_avect_nRattr(x2l_l) ! number of fields retrived from MOAB tags, based on names from seq_flds_x2l_fields
    totalmblsimp = mblsize * nrecv ! size of the double array to fill with data from MOAB
    allocate (x2l_lm(lsz, nrecv) )
    if (masterproc) then
       write(iulog,*) sub, 'mblsize= ',mblsize,' nsend, nrecv for moab:', nsend, nrecv
    end if
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
!     Also send data through the MOAB path in driver-moab
      call lnd_export_moab(EClock, bounds, lnd2atm_vars, lnd2glc_vars) ! it is private here
#endif
    endif

    ! Fill in infodata settings

    call seq_infodata_PutData(infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData(infodata, lnd_nx=ldomain%ni, lnd_ny=ldomain%nj, precip_downscaling_method = precip_downscaling_method)

#ifdef HAVE_MOAB
    call seq_infodata_PutData(infodata, lnd_domain= fatmlndfrc)
#endif

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
    use elm_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use elm_time_manager,  only : advance_timestep, set_nextsw_cday,update_rad_dtime
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
#ifdef MOABCOMP
    use seq_flds_mod     , only :   seq_flds_x2l_fields
#endif
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
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list, ent_type
    type(mct_string)    :: mctOStr  !
    character(100) ::tagname, mct_field, modelStr
#endif 
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
    call lnd_import( bounds, x2l_l%rattr, atm2lnd_vars, glc2lnd_vars, lnd2atm_vars)
    
#ifdef HAVE_MOAB
    ! first call moab import 
#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2l_fields
    call mct_list_init(temp_list ,seq_flds_x2l_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 0 ! entity type is vertex for land, always
    if (rank2 .eq. 0) print *, num_moab_exports, trim(seq_flds_x2l_fields), ' lnd import check'
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      modelStr = 'lnd run'
      call seq_comm_compare_mb_mct(modelStr, mpicom_lnd_moab, x2l_l, mct_field,  mlnid, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)

#endif
! calling MOAB's import last means this is what the model will use.
    call lnd_import_moab( EClock, bounds, atm2lnd_vars, glc2lnd_vars)
#endif

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
       call lnd_export_moab(EClock, bounds, lnd2atm_vars, lnd2glc_vars) ! it is private here
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
      deallocate (x2l_lm)
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
  subroutine init_moab_land(bounds, LNDID)
    use seq_flds_mod     , only :  seq_flds_l2x_fields, seq_flds_x2l_fields
    use shr_kind_mod     , only : CXX => SHR_KIND_CXX
    use spmdMod     , only: iam  ! rank on the land communicator
    use domainMod   , only: ldomain ! ldomain is coming from module, not even passed
    use elm_varcon  , only: re
    use shr_const_mod, only: SHR_CONST_PI
    use elm_varctl  ,  only : iulog  ! for messages 
     use spmdmod          , only: masterproc
    use iMOAB        , only: iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_RegisterApplication, &
    iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
    iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo

    type(bounds_type) , intent(in)  :: bounds
    integer , intent(in) :: LNDID ! id of the land app

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
    character*32  appname

    integer, allocatable :: moabconn(:) ! will have the connectivity in terms of local index in verts

    ! define a MOAB app for ELM
    appname="LNDMB"//C_NULL_CHAR
    ! first land instance, should be 9
    ierr = iMOAB_RegisterApplication(appname, mpicom_lnd_moab, LNDID, mlnid)
    if (ierr > 0 )  &
       call endrun('Error: cannot register moab app')
    if(masterproc) then
       write(iulog,*) " "
       write(iulog,*) "register MOAB app:", trim(appname), "  mlnid=", mlnid
       write(iulog,*) " "
    endif

    ! start describing the mesh to MOAB

    dims  =3 ! store as 3d mesh
    ! number the local grid
    lsz = bounds%endg - bounds%begg + 1

    allocate(vgids(lsz)) ! use it for global ids, for elements in full mesh or vertices in point cloud

    do n = 1, lsz
       vgids(n) = ldecomp%gdc2glo(bounds%begg+n-1) ! local to global !
    end do
    gsize = ldomain%ni * ldomain%nj ! size of the total grid

    
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
    ! ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
    ! if (ierr > 0 )  &
    !   call endrun('Error: fail to set aream tag ')
    ierr = iMOAB_UpdateMeshInfo( mlnid )
    if (ierr > 0 )  &
      call endrun('Error: fail to update mesh info ')

  ! add more domain fields that are missing from domain fields: lat, lon, mask, hgt
  tagname = 'lat:lon:mask:hgt'//C_NULL_CHAR
  tagtype = 1 ! dense, double
  numco = 1
  ierr = iMOAB_DefineTagStorage(mlnid, tagname, tagtype, numco,  tagindex )
  if (ierr > 0 )  &
    call endrun('Error: fail to create lat:lon:mask:hgt tags ')

 ! moab_vert_coords is big enough in both case to hold enough data for us: lat, lon, mask
    do i = 1, lsz
      n = i-1 + bounds%begg
      moab_vert_coords(i) = ldomain%latc(n)  ! lat
      moab_vert_coords(lsz + i) = ldomain%lonc(n) ! lon
      moab_vert_coords(2*lsz + i) = ldomain%mask(n) ! mask
    enddo
    tagname = 'lat:lon:mask'//C_NULL_CHAR

    ent_type = 0 ! point cloud usually

    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz*3 , ent_type, moab_vert_coords)
    if (ierr > 0 )  &
      call endrun('Error: fail to set lat lon mask tag ')
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

  end subroutine init_moab_land

  subroutine lnd_export_moab(EClock, bounds, lnd2atm_vars, lnd2glc_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the elm model to the moab coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use shr_kind_mod       , only : CXX => SHR_KIND_CXX
    use elm_varctl         , only : iulog, create_glacier_mec_landunit
    use elm_time_manager   , only : get_nstep, get_step_size  
    use domainMod          , only : ldomain
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use iMOAB,  only       : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
    use seq_flds_mod, only : seq_flds_l2x_fields
    use seq_timemgr_mod, only : seq_timemgr_eclockgetdata
    !
    ! !ARGUMENTS:
    !implicit none
    type(ESMF_Clock),   intent(inout) :: EClock
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

    integer :: ent_type, ierr, cur_lnd_stepno
    character(len=100) :: outfile, wopts, lnum
    character(CXX) :: tagname
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
    ent_type = 0 ! vertices only, from now on
    ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, totalmbls , ent_type, l2x_lm(1,1) )
    if (ierr > 0 )  &
       call shr_sys_abort( sub//' Error: fail to set moab l2x '// trim(seq_flds_l2x_fields) )

    call seq_timemgr_EClockGetData( EClock, stepno=cur_lnd_stepno )
#ifdef MOABDEBUG
       write(lnum,"(I0.2)")cur_lnd_stepno
       outfile = 'lnd_export_'//trim(lnum)//'.h5m'//C_NULL_CHAR
       wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
       ierr = iMOAB_WriteMesh(mlnid, outfile, wopts)
       if (ierr > 0 )  &
         call shr_sys_abort( sub//' fail to write the land mesh file with data')
#endif

  end subroutine lnd_export_moab

  ! lnd_import_moab will be a copy of lnd_import
  ! data will come from moab tags, from mlnid iMOAB app 
  ! the role of x2l_l AV is taken by the local array x2l_lm, allocated at init stage, with
  ! the order of tags given by seq_flds_x2l_fields 

  !===============================================================================
  subroutine lnd_import_moab(EClock, bounds, atm2lnd_vars, glc2lnd_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the moab coupler to the land model 
    use seq_flds_mod     , only :  seq_flds_l2x_fields, seq_flds_x2l_fields
    use iMOAB,  only       : iMOAB_GetDoubleTagStorage, iMOAB_WriteMesh
    use shr_kind_mod     , only : CXX => SHR_KIND_CXX
    ! !USES:
    use elm_varctl       , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit, &
                                 metdata_type, metdata_bypass, metdata_biases, co2_file, aero_file
    use elm_varctl       , only: const_climate_hist, add_temperature, add_co2, use_cn, use_fates
    use elm_varctl       , only: startdate_add_temperature, startdate_add_co2
    use elm_varcon       , only: rair, o2_molar_const, c13ratio
    use elm_time_manager , only: get_nstep, get_step_size, get_curr_calday, get_curr_date 
    use controlMod       , only: NLFilename
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use domainMod        , only: ldomain
    use shr_kind_mod     , only: r8 => shr_kind_r8, CL => shr_kind_CL
    use fileutils        , only: getavu, relavu
    use spmdmod          , only: masterproc, mpicom, iam, npes, MPI_REAL8, MPI_INTEGER, MPI_STATUS_SIZE
    use elm_nlUtilsMod   , only : find_nlgroup_name
    use seq_timemgr_mod, only : seq_timemgr_eclockgetdata
    use netcdf
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),    intent(inout) :: EClock
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars      ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars      ! clm internal input data type
    !
    ! !LOCAL VARIABLES:
    integer  :: g,topo,i,m,thism,nstep,ier  ! indices, number of steps, and error code
    integer status(MPI_STATUS_SIZE)
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e, ea                ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t               ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q               ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    real(r8) :: vp                   ! water vapor pressure (Pa)
    integer  :: thisng, np, num, nu_nml, nml_error                 
    integer  :: ng_all(100000)
    real(r8) :: swndf, swndr, swvdf, swvdr, ratio_rvrf, frac, q
    real(r8) :: thiscosz, avgcosz, szenith
    integer  :: swrad_period_len, swrad_period_start, thishr, thismin
    real(r8) :: timetemp(2)
    real(r8) :: latixy(500000), longxy(500000)
    integer ::  ierr, varid, dimid, yr, mon, day, tod, nindex(2), caldaym(13)
    integer ::  ncid, met_ncids(14), mask_ncid, thisncid, ng, tm
    integer ::  aindex(2), tindex(14,2), starti(3), counti(3)
    integer ::  grid_map(500000), zone_map(500000)
    integer ::  met_nvars, nyears_spinup, nyears_trans, starti_site, endi_site
    real(r8) :: smap05_lat(360), smap05_lon(720)
    real(r8) :: smapt62_lat(94), smapt62_lon(192)
    real(r8) :: smap2_lat(96), smap2_lon(144)
    real(r8) :: thisdist, mindist, thislon
    real(r8) :: tbot, tempndep(1,1,158), thiscalday, wt1(14), wt2(14), thisdoy
    real(r8) :: site_metdata(14,12)
    real(r8) :: var_month_mean(12)
    integer  :: lnfmind(2)
    integer  :: var_month_count(12)
    integer*2 :: temp(1,500000)
    integer :: xtoget, ytoget, thisx, thisy, calday_start
    integer :: sdate_addt, sy_addt, sm_addt, sd_addt
    integer :: sdate_addco2, sy_addco2, sm_addco2, sd_addco2
    character(len=200) metsource_str, thisline
    character(len=*), parameter :: sub = 'lnd_import_moab'
    integer :: av, v, n, nummetdims, g3, gtoget, ztoget, line, mystart, tod_start, thistimelen  
    character(len=20) aerovars(14), metvars(14)
    character(len=3) zst
    integer :: stream_year_first_lightng, stream_year_last_lightng, model_year_align_lightng
    integer :: stream_year_first_popdens, stream_year_last_popdens, model_year_align_popdens
    integer :: stream_year_first_ndep,    stream_year_last_ndep,    model_year_align_ndep
    character(len=CL)  :: metdata_fname  
    character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
    character(len=CL)  :: popdensmapalgo = 'bilinear' 
    character(len=CL)  :: ndepmapalgo    = 'bilinear' 
    character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
    character(len=CL)  :: stream_fldFileName_popdens ! poplulation density stream filename
    character(len=CL)  :: stream_fldFileName_ndep    ! nitrogen deposition stream filename
    logical :: use_sitedata, has_zonefile, use_daymet, use_livneh

! moab extra stuff 
    character(CXX) ::  tagname ! hold all fields names
    integer        :: ent_type  ! for setting data 
    integer        :: cur_lnd_stepno
#ifdef MOABDEBUG
    character*100 outfile, wopts, lnum
#endif

    data caldaym / 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /    

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)
    !
    ! function declarations
    !
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

    namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

    namelist /ndepdyn_nml/        &
        stream_year_first_ndep,  &
    stream_year_last_ndep,   &
        model_year_align_ndep,   &
        ndepmapalgo,             &
        stream_fldFileName_ndep

    stream_fldFileName_lightng = ' '
    stream_fldFileName_popdens = ' '
   
    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if


    call seq_timemgr_EClockGetData( EClock, stepno=cur_lnd_stepno )
#ifdef MOABDEBUG
    write(lnum,"(I0.2)")cur_lnd_stepno
    outfile = 'LndImp_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mlnid, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the moab lnd mesh before import ')
#endif
    tagname=trim(seq_flds_x2l_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices 
    ierr = iMOAB_GetDoubleTagStorage ( mlnid, tagname, totalmblsimp , ent_type, x2l_lm(1,1) )
    if ( ierr > 0) then
      call endrun('Error: fail to get seq_flds_x2l_fields for land moab instance on component')
    endif

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    thisng = bounds%endg - bounds%begg + 1
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       
       ! Determine flooding input, sign convention is positive downward and
       ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
       ! change the sign to indicate addition of water to system.

       atm2lnd_vars%forc_flood_grc(g)   = -x2l_lm(i,index_x2l_Flrr_flood)  

       atm2lnd_vars%volr_grc(g)   = x2l_lm(i,index_x2l_Flrr_volr) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%volrmch_grc(g)= x2l_lm(i,index_x2l_Flrr_volrmch) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%supply_grc(g) = x2l_lm(i,index_x2l_Flrr_supply)
       atm2lnd_vars%deficit_grc(g) = x2l_lm(i,index_x2l_Flrr_deficit)

       atm2lnd_vars%forc_hgt_grc(g)     = x2l_lm(i,index_x2l_Sa_z)         ! zgcmxy  Atm state m
       atm2lnd_vars%forc_u_grc(g)       = x2l_lm(i,index_x2l_Sa_u)         ! forc_uxy  Atm state m/s
       atm2lnd_vars%forc_v_grc(g)       = x2l_lm(i,index_x2l_Sa_v)         ! forc_vxy  Atm state m/s
       atm2lnd_vars%forc_solad_grc(g,2) = x2l_lm(i,index_x2l_Faxa_swndr)   ! forc_sollxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solad_grc(g,1) = x2l_lm(i,index_x2l_Faxa_swvdr)   ! forc_solsxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,2) = x2l_lm(i,index_x2l_Faxa_swndf)   ! forc_solldxy Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,1) = x2l_lm(i,index_x2l_Faxa_swvdf)   ! forc_solsdxy Atm flux  W/m^2

       atm2lnd_vars%forc_th_not_downscaled_grc(g)    = x2l_lm(i,index_x2l_Sa_ptem)      ! forc_thxy Atm state K
       atm2lnd_vars%forc_q_not_downscaled_grc(g)     = x2l_lm(i,index_x2l_Sa_shum)      ! forc_qxy  Atm state kg/kg
       atm2lnd_vars%forc_pbot_not_downscaled_grc(g)  = x2l_lm(i,index_x2l_Sa_pbot)      ! ptcmxy  Atm state Pa
       atm2lnd_vars%forc_t_not_downscaled_grc(g)     = x2l_lm(i,index_x2l_Sa_tbot)      ! forc_txy  Atm state K
       atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = x2l_lm(i,index_x2l_Faxa_lwdn)    ! flwdsxy Atm flux  W/m^2

       forc_rainc                                    = x2l_lm(i,index_x2l_Faxa_rainc)   ! mm/s
       forc_rainl                                    = x2l_lm(i,index_x2l_Faxa_rainl)   ! mm/s
       forc_snowc                                    = x2l_lm(i,index_x2l_Faxa_snowc)   ! mm/s
       forc_snowl                                    = x2l_lm(i,index_x2l_Faxa_snowl)   ! mm/s

       ! atmosphere coupling, for prognostic/prescribed aerosols
       atm2lnd_vars%forc_aer_grc(g,1)  =  x2l_lm(i,index_x2l_Faxa_bcphidry)
       atm2lnd_vars%forc_aer_grc(g,2)  =  x2l_lm(i,index_x2l_Faxa_bcphodry)
       atm2lnd_vars%forc_aer_grc(g,3)  =  x2l_lm(i,index_x2l_Faxa_bcphiwet)
       atm2lnd_vars%forc_aer_grc(g,4)  =  x2l_lm(i,index_x2l_Faxa_ocphidry)
       atm2lnd_vars%forc_aer_grc(g,5)  =  x2l_lm(i,index_x2l_Faxa_ocphodry)
       atm2lnd_vars%forc_aer_grc(g,6)  =  x2l_lm(i,index_x2l_Faxa_ocphiwet)
       atm2lnd_vars%forc_aer_grc(g,7)  =  x2l_lm(i,index_x2l_Faxa_dstwet1)
       atm2lnd_vars%forc_aer_grc(g,8)  =  x2l_lm(i,index_x2l_Faxa_dstdry1)
       atm2lnd_vars%forc_aer_grc(g,9)  =  x2l_lm(i,index_x2l_Faxa_dstwet2)
       atm2lnd_vars%forc_aer_grc(g,10) =  x2l_lm(i,index_x2l_Faxa_dstdry2)
       atm2lnd_vars%forc_aer_grc(g,11) =  x2l_lm(i,index_x2l_Faxa_dstwet3)
       atm2lnd_vars%forc_aer_grc(g,12) =  x2l_lm(i,index_x2l_Faxa_dstdry3)
       atm2lnd_vars%forc_aer_grc(g,13) =  x2l_lm(i,index_x2l_Faxa_dstwet4)
       atm2lnd_vars%forc_aer_grc(g,14) =  x2l_lm(i,index_x2l_Faxa_dstdry4)
       
       !set the topounit-level atmospheric state and flux forcings
       do topo = grc_pp%topi(g), grc_pp%topf(g)
         ! first, all the state forcings
         top_as%tbot(topo)    = x2l_lm(i,index_x2l_Sa_tbot)      ! forc_txy  Atm state K
         top_as%thbot(topo)   = x2l_lm(i,index_x2l_Sa_ptem)      ! forc_thxy Atm state K
         top_as%pbot(topo)    = x2l_lm(i,index_x2l_Sa_pbot)      ! ptcmxy    Atm state Pa
         top_as%qbot(topo)    = x2l_lm(i,index_x2l_Sa_shum)      ! forc_qxy  Atm state kg/kg
         top_as%ubot(topo)    = x2l_lm(i,index_x2l_Sa_u)         ! forc_uxy  Atm state m/s
         top_as%vbot(topo)    = x2l_lm(i,index_x2l_Sa_v)         ! forc_vxy  Atm state m/s
         top_as%zbot(topo)    = x2l_lm(i,index_x2l_Sa_z)         ! zgcmxy    Atm state m
         ! assign the state forcing fields derived from other inputs
         ! Horizontal windspeed (m/s)
         top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
         ! Relative humidity (percent)
         if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
            e = esatw(tdc(top_as%tbot(topo)))
         else
            e = esati(tdc(top_as%tbot(topo)))
         end if
         qsat = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
         top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
         ! partial pressure of oxygen (Pa)
         top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
         ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
         vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
         top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))
         
         ! second, all the flux forcings
         top_af%rain(topo)    = forc_rainc + forc_rainl       ! sum of convective and large-scale rain
         top_af%snow(topo)    = forc_snowc + forc_snowl       ! sum of convective and large-scale snow
         top_af%solad(topo,2) = x2l_lm(i,index_x2l_Faxa_swndr)   ! forc_sollxy  Atm flux  W/m^2
         top_af%solad(topo,1) = x2l_lm(i,index_x2l_Faxa_swvdr)   ! forc_solsxy  Atm flux  W/m^2
         top_af%solai(topo,2) = x2l_lm(i,index_x2l_Faxa_swndf)   ! forc_solldxy Atm flux  W/m^2
         top_af%solai(topo,1) = x2l_lm(i,index_x2l_Faxa_swvdf)   ! forc_solsdxy Atm flux  W/m^2
         top_af%lwrad(topo)   = x2l_lm(i,index_x2l_Faxa_lwdn)    ! flwdsxy Atm flux  W/m^2
         ! derived flux forcings
         top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
       end do
         

       ! Determine optional receive fields
       ! CO2 (and C13O2) concentration: constant, prognostic, or diagnostic
       if (co2_type_idx == 0) then                    ! CO2 constant, value from namelist
         co2_ppmv_val = co2_ppmv
       else if (co2_type_idx == 1) then               ! CO2 prognostic, value from coupler field
         co2_ppmv_val = x2l_lm(i,index_x2l_Sa_co2prog)
       else if (co2_type_idx == 2) then               ! CO2 diagnostic, value from coupler field
         co2_ppmv_val = x2l_lm(i,index_x2l_Sa_co2diag)
       else
         call endrun( sub//' ERROR: Invalid co2_type_idx, must be 0, 1, or 2 (constant, prognostic, or diagnostic)' )
       end if
       ! Assign to topounits, with conversion from ppmv to partial pressure (Pa)
       ! If using C13, then get the c13ratio from elm_varcon (constant value for pre-industrial atmosphere)

       do topo = grc_pp%topi(g), grc_pp%topf(g)
         top_as%pco2bot(topo) = co2_ppmv_val * 1.e-6_r8 * top_as%pbot(topo)
         if (use_c13) then
            top_as%pc13o2bot(topo) = top_as%pco2bot(topo) * c13ratio;
         end if
       end do
       ! CH4
       if (index_x2l_Sa_methane /= 0) then
          do topo = grc_pp%topi(g), grc_pp%topf(g)
            top_as%pch4bot(topo) = x2l_lm(i,index_x2l_Sa_methane)
          end do
       endif

       if (index_x2l_Sa_co2prog /= 0) then
          co2_ppmv_prog = x2l_lm(i,index_x2l_Sa_co2prog)   ! co2 atm state prognostic
       else
          co2_ppmv_prog = co2_ppmv
       end if

       if (index_x2l_Sa_co2diag /= 0) then
          co2_ppmv_diag = x2l_lm(i,index_x2l_Sa_co2diag)   ! co2 atm state diagnostic
       else
          co2_ppmv_diag = co2_ppmv
       end if

       if (index_x2l_Sa_methane /= 0) then
          atm2lnd_vars%forc_pch4_grc(g) = x2l_lm(i,index_x2l_Sa_methane)
       endif

       ! Determine derived quantities for required fields

       forc_t = atm2lnd_vars%forc_t_not_downscaled_grc(g)
       forc_q = atm2lnd_vars%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)
       
       atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of humidity [m]
       atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
       atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
       atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
       atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
       atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                        atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)
       
       atm2lnd_vars%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
       atm2lnd_vars%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl
       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
       atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
       ! Make sure relative humidity is properly bounded
       ! atm2lnd_vars%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_vars%forc_rh_grc(g) )
       ! atm2lnd_vars%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_vars%forc_rh_grc(g) )

       ! Determine derived quantities for optional fields
       ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
       ! Note that forc_pbot is in Pa

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then

          co2_ppmv_val = co2_ppmv_diag 
           if (use_c13) then
             atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
           end if
       else
          co2_ppmv_val = co2_ppmv
          if (use_c13) then
            atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
          end if
       end if
       atm2lnd_vars%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot 

       ! glc coupling 

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             glc2lnd_vars%frac_grc(g,num)  = x2l_lm(i,index_x2l_Sg_frac(num))
             glc2lnd_vars%topo_grc(g,num)  = x2l_lm(i,index_x2l_Sg_topo(num))
             glc2lnd_vars%hflx_grc(g,num)  = x2l_lm(i,index_x2l_Flgg_hflx(num))
          end do
          glc2lnd_vars%icemask_grc(g)  = x2l_lm(i,index_x2l_Sg_icemask)
          glc2lnd_vars%icemask_coupled_fluxes_grc(g)  = x2l_lm(i,index_x2l_Sg_icemask_coupled_fluxes)
       end if

    end do     
  end subroutine lnd_import_moab

! endif for ifdef HAVE_MOAB
#endif

end module lnd_comp_mct
