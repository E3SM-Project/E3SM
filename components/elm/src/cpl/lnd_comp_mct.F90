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

#ifdef HAVE_MOAB
  use seq_comm_mct,       only: mlnid! id of moab land app
  use seq_comm_mct,       only: mb_land_mesh! true if land is full mesh
  use seq_comm_mct,       only: num_moab_exports
  use seq_flds_mod     , only : seq_flds_x2l_fields, seq_flds_l2x_fields

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
  logical :: sameg_al ! save it for export :)  

#ifdef HAVE_MOAB
  integer  :: mpicom_lnd_moab ! used also for mpi-reducing the difference between moab tags and mct avs
  integer :: rank2
#endif

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
    use elm_instance     , only : elm_instance_init
    use elm_varctl       , only : finidat,single_column, elm_varctl_set, iulog, noland, fatmlndfrc
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
    integer :: ierr, nsend
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

    ! Set and save LNDID for easy access by other modules

    call elm_instance_init( LNDID )

    ! Determine attriute vector indices
#ifdef HAVE_MOAB
    mpicom_lnd_moab = mpicom_lnd ! just store it now, for later use; maybe it is the same as mpicom from spmdMod  (or a copy)
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
!   find out samegrid_al or not; from infodata
    samegrid_al = .true.
    call seq_infodata_GetData(infodata         , &
                   atm_gnam=atm_gnam           , &
                   lnd_gnam=lnd_gnam           )
    if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
    mb_land_mesh = .not. samegrid_al ! global variable, saved in seq_comm
    call init_moab_land(bounds, samegrid_al, LNDID)
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
      call lnd_export_moab(bounds, lnd2atm_vars, lnd2glc_vars) ! it is private here
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
#ifdef HAVE_MOAB
    ! first call moab import 
#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2l_fields
    call mct_list_init(temp_list ,seq_flds_x2l_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 0 ! entity type is vertex for land, usually (bigrid case)
    if (mb_land_mesh) ent_type = 1 
    if (rank2 .eq. 0) print *, num_moab_exports, trim(seq_flds_x2l_fields), ' lnd import check'
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      modelStr = 'lnd run'
      !call compare_to_moab_tag_lnd(mpicom_lnd_moab, x2l_l, mct_field,  mlnid, tagname, ent_type, difference)
      call seq_comm_compare_mb_mct(modelStr, mpicom_lnd_moab, x2l_l, mct_field,  mlnid, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)

#endif
    call lnd_import_moab( bounds, atm2lnd_vars, glc2lnd_vars)
#endif

    call lnd_import( bounds, x2l_l%rattr, atm2lnd_vars, glc2lnd_vars, lnd2atm_vars)
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
  subroutine init_moab_land(bounds, samegrid_al, LNDID)
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
    iMOAB_ResolveSharedEntities, iMOAB_CreateElements, iMOAB_MergeVertices, iMOAB_UpdateMeshInfo

    type(bounds_type) , intent(in)  :: bounds
    logical , intent(in) :: samegrid_al
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
        ! ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        ! if (ierr > 0 )  &
        !   call endrun('Error: fail to set aream tag ')

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
        ! ierr = iMOAB_SetDoubleTagStorage ( mlnid, tagname, lsz , ent_type, moab_vert_coords )
        ! if (ierr > 0 )  &
        !   call endrun('Error: fail to set aream tag ')
        ierr = iMOAB_UpdateMeshInfo( mlnid )
        if (ierr > 0 )  &
          call endrun('Error: fail to update mesh info ')
    endif
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
    if (ldomain%nv .ge. 3 .and.  .not.samegrid_al) then
      ent_type = 1 ! cell in tri-grid case
    endif
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

  ! lnd_import_moab will be a copy of lnd_import
  ! data will come from moab tags, from mlnid iMOAB app 
  ! the role of x2l_l AV is taken by the local array x2l_lm, allocated at init stage, with
  ! the order of tags given by seq_flds_x2l_fields 

  !===============================================================================
  subroutine lnd_import_moab( bounds, atm2lnd_vars, glc2lnd_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the moab coupler to the land model 
    use seq_flds_mod     , only :  seq_flds_l2x_fields, seq_flds_x2l_fields
    use iMOAB,  only       : iMOAB_GetDoubleTagStorage
    use shr_kind_mod     , only : CXX => SHR_KIND_CXX
    ! !USES:
    use elm_varctl       , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit, &
                                 metdata_type, metdata_bypass, metdata_biases, co2_file, aero_file
    use elm_varctl       , only: const_climate_hist, add_temperature, add_co2, use_cn, use_fates
    use elm_varctl       , only: startdate_add_temperature, startdate_add_co2
    use elm_varcon       , only: rair, o2_molar_const, c13ratio
    use clm_time_manager , only: get_nstep, get_step_size, get_curr_calday, get_curr_date 
    use controlMod       , only: NLFilename
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use domainMod        , only: ldomain
    use shr_kind_mod     , only: r8 => shr_kind_r8, CL => shr_kind_CL
    use fileutils        , only: getavu, relavu
    use spmdmod          , only: masterproc, mpicom, iam, npes, MPI_REAL8, MPI_INTEGER, MPI_STATUS_SIZE
    use elm_nlUtilsMod   , only : find_nlgroup_name
    use netcdf
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    ! real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    ! this is moab version, will be replaced with x2l_lm from mlnid  
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
    !real(r8) :: hdm1(720,360,1), hdm2(720,360,1) 
    !real(r8) :: lnfm1(192,94,2920)
    !real(r8) :: ndep1(144,96,1), ndep2(144,96,1)
    !real(r8) :: aerodata(14,144,96,14)
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

    tagname=trim(seq_flds_x2l_fields)//C_NULL_CHAR
    if (sameg_al) then
      ent_type = 0 ! vertices, cells only if sameg_al false
    else
      ent_type = 1
    endif
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

       ! Determine required receive fields

#ifdef CPL_BYPASS
        !read forcing data directly, bypass coupler
        atm2lnd_vars%forc_flood_grc(g)   = 0._r8
        atm2lnd_vars%volr_grc(g)   = 0._r8 

        !Get meteorological data, concatenated to include whole record
        !Note we only do this at the first timestep and keep the whole forcing dataset in the memory
       
  !-----------------------------------Meteorological forcing  -----------------------------------

        call get_curr_date( yr, mon, day, tod )
        thiscalday = get_curr_calday()
        nstep = get_nstep()

        !on first timestep, read all the met data for relevant gridcell(s) and store in array.
        !   Met data are held in short integer format to save memory.
        !   Each node must have enough memory to hold these data.
        met_nvars=7
        if (metdata_type(1:3) == 'cpl') met_nvars=14

        if (atm2lnd_vars%loaded_bypassdata == 0) then
          !meteorological forcing
          if (index(metdata_type, 'qian') .gt. 0) then 
            atm2lnd_vars%metsource = 0   
          else if (index(metdata_type,'cru') .gt. 0) then
            atm2lnd_vars%metsource = 1  
          else if (index(metdata_type,'site') .gt. 0) then 
            atm2lnd_vars%metsource = 2
          else if (index(metdata_type,'princeton') .gt. 0) then 
            atm2lnd_vars%metsource = 3
          else if (index(metdata_type,'gswp3') .gt. 0) then
            atm2lnd_vars%metsource = 4
          else if (index(metdata_type,'cpl') .gt. 0) then 
            atm2lnd_vars%metsource = 5
          else
            call endrun( sub//' ERROR: Invalid met data source for cpl_bypass' )
          end if

          use_livneh = .false.
          use_daymet = .false.
          if(index(metdata_type, 'livneh') .gt. 0) then 
              use_livneh = .true.
          else if (index(metdata_type, 'daymet') .gt. 0) then 
              use_daymet = .true.
          end if
 
          metvars(1) = 'TBOT'
          metvars(2) = 'PSRF'
          metvars(3) = 'QBOT'
          if (atm2lnd_vars%metsource .eq. 2) metvars(3) = 'RH'
          if (atm2lnd_vars%metsource .ne. 5) metvars(4) = 'FSDS'
          if (atm2lnd_vars%metsource .ne. 5) metvars(5) = 'PRECTmms'
          if (atm2lnd_vars%metsource .ne. 5) metvars(6) = 'WIND'
          metvars(4) = 'FSDS'
          metvars(5) = 'PRECTmms'
          metvars(6) = 'WIND'
          metvars(7) = 'FLDS'
          if (atm2lnd_vars%metsource .eq. 5) then 
              metvars(4) = 'SWNDF'
              metvars(5) = 'RAINC'
              metvars(6) = 'U'
              metvars(8) = 'SWNDR'
              metvars(9) = 'SWVDF'
              metvars(10) = 'SWVDR'
              metvars(11) = 'RAINL'
              metvars(12) = 'SNOWC'
              metvars(13) = 'SNOWL'
              metvars(14) = 'V'
          else
              metvars(4) = 'FSDS'
              metvars(5) = 'PRECTmms'
              metvars(6) = 'WIND'
          end if

          !set defaults
          atm2lnd_vars%startyear_met       = 1901
          atm2lnd_vars%endyear_met_spinup  = 1920
          if (atm2lnd_vars%metsource == 0) then 
            metsource_str = 'qian'
            atm2lnd_vars%startyear_met       = 1948
            atm2lnd_vars%endyear_met_spinup  = 1972
            atm2lnd_vars%endyear_met_trans   = 2004
          else if (atm2lnd_vars%metsource == 1) then 
            metsource_str = 'cruncep'
            atm2lnd_vars%endyear_met_trans  = 2016
          else if (atm2lnd_vars%metsource == 2) then
            metsource_str = 'site'
            !get year information from file
            ierr = nf90_open(trim(metdata_bypass) // '/all_hourly.nc', nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'start_year', varid) 
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%startyear_met)
            ierr = nf90_inq_varid(ncid, 'end_year', varid)
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%endyear_met_spinup)
            ierr = nf90_close(ncid)
            atm2lnd_vars%endyear_met_trans = atm2lnd_vars%endyear_met_spinup
          else if (atm2lnd_vars%metsource == 3) then 
            metsource_str = 'princeton'
            atm2lnd_vars%endyear_met_trans = 2012 
          else if (atm2lnd_vars%metsource == 4) then 
            atm2lnd_vars%endyear_met_trans  = 2014
          else if (atm2lnd_vars%metsource == 5) then
            atm2lnd_vars%startyear_met      = 566 !76
            atm2lnd_vars%endyear_met_spinup = 590 !100
            atm2lnd_vars%endyear_met_trans  = 590 !100
          end if

          if (use_livneh) then 
              atm2lnd_vars%startyear_met      = 1950
              atm2lnd_vars%endyear_met_spinup = 1969
          else if (use_daymet) then 
              atm2lnd_vars%startyear_met      = 1980
              atm2lnd_vars%endyear_met_spinup = atm2lnd_vars%endyear_met_trans
          end if

          nyears_spinup = atm2lnd_vars%endyear_met_spinup - &
                             atm2lnd_vars%startyear_met + 1
          nyears_trans  = atm2lnd_vars%endyear_met_trans - &
                             atm2lnd_vars%startyear_met  + 1

          !check for site data in run directory (monthly mean T, precip)
          inquire(file=trim(metdata_biases), exist=use_sitedata)

          !get grid lat/lon information, zone mappings
          inquire(file=trim(metdata_bypass) // '/zone_mappings.txt', exist=has_zonefile)
          if (has_zonefile) then
            open(unit=13, file=trim(metdata_bypass) // '/zone_mappings.txt')
          else if (atm2lnd_vars%metsource .ne. 2) then
            call endrun( sub//' ERROR: Zone mapping file does not exist for cpl_bypass' )
          end if

          if (atm2lnd_vars%metsource .ne. 2) then 
            ng = 0     !number of points
            do v=1,500000
              read(13,*, end=10), longxy(v), latixy(v), zone_map(v), grid_map(v)
              ng = ng + 1
            end do
10          continue
            close(unit=13)

            !Figure out the closest point and which zone file to open
            mindist=99999
            do g3 = 1,ng
              thisdist = 100*((latixy(g3) - ldomain%latc(g))**2 + &
                              (longxy(g3) - ldomain%lonc(g))**2)**0.5
              if (thisdist .lt. mindist) then 
                mindist = thisdist
                ztoget = zone_map(g3)
                gtoget = grid_map(g3)
              end if
            end do
          else
            gtoget = 1
          end if

          !get the site metdata for bias correction if they exist (lat/lons must match domain file)
          if (use_sitedata) then 
            open(unit=9, file=trim(metdata_biases),status='old')
            read(9,*) thisline
            site_metdata(:,:)=-999._r8
            do while ((site_metdata(1,1) .lt. ldomain%lonc(g) - 0.01 .or. &
                    site_metdata(1,1) .gt. ldomain%lonc(g) + 0.01) .and. &
                      (site_metdata(2,1) .lt. ldomain%latc(g) - 0.01 .or. &
                       site_metdata(2,1) .gt. ldomain%latc(g) + 0.01))
              read(9,*) site_metdata(1:7,1)
              if (site_metdata(1,1) .lt. 0) site_metdata(1,1) = site_metdata(1,1)+360._r8
            end do
            do line=2,12
              read(9,*) site_metdata(1:7,line)
            end do
            close(unit=9)
          end if

          do v=1,met_nvars
            write(zst, '(I3)') 100+ztoget
            if (atm2lnd_vars%metsource == 0) then 
                metdata_fname =  trim(metsource_str) // '_' // trim(metvars(v)) // '_z' // zst(2:3) // '.nc'
            else if (atm2lnd_vars%metsource == 1) then 
                metdata_fname = 'CRUNCEP.v5_' // trim(metvars(v)) // '_1901-2013_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'CRUNCEP5_Livneh_' // trim(metvars(v)) // '_1950-2013_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'CRUNCEP5_Daymet3_' // trim(metvars(v)) // '_1980-2013_z' // zst(2:3) // '.nc'
                end if
            else if (atm2lnd_vars%metsource == 2) then
                metdata_fname = 'all_hourly.nc'
            else if (atm2lnd_vars%metsource == 3) then 
               metdata_fname = 'Princeton_' // trim(metvars(v)) // '_1901-2012_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Livneh_' // trim(metvars(v)) // '_1950-2012_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Daymet3_' // trim(metvars(v)) // '_1980-2012_z' // zst(2:3) // '.nc'
                end if
            else if (atm2lnd_vars%metsource == 4) then 
                metdata_fname = 'GSWP3_' // trim(metvars(v)) // '_1901-2014_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'GSWP3_Livneh_' // trim(metvars(v)) // '_1950-2010_z' // zst(2:3) // '.nc'                
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'GSWP3_Daymet3_' // trim(metvars(v)) // '_1980-2010_z' // zst(2:3) // '.nc' 
                end if
            else if (atm2lnd_vars%metsource == 5) then 
                    !metdata_fname = 'WCYCL1850S.ne30_' // trim(metvars(v)) // '_0076-0100_z' // zst(2:3) // '.nc'
                    metdata_fname = 'CBGC1850S.ne30_' // trim(metvars(v)) // '_0566-0590_z' // zst(2:3) // '.nc'
            end if
  
            ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, met_ncids(v))
            if (ierr .ne. 0) call endrun(msg=' ERROR: Failed to open cpl_bypass input meteorology file' )
       
            !get timestep information
            ierr = nf90_inq_dimid(met_ncids(v), 'DTIME', dimid)
            ierr = nf90_Inquire_Dimension(met_ncids(v), dimid, len = atm2lnd_vars%timelen(v))

            starti(1) = 1
            counti(1) = 2
            ierr = nf90_inq_varid(met_ncids(v), 'DTIME', varid)
            ierr = nf90_get_var(met_ncids(v), varid, timetemp, starti(1:1), counti(1:1))   
            atm2lnd_vars%timeres(v)        = (timetemp(2)-timetemp(1))*24._r8
            atm2lnd_vars%npf(v)            = 86400d0*(timetemp(2)-timetemp(1))/get_step_size()  
            atm2lnd_vars%timelen_spinup(v) = nyears_spinup*(365*nint(24./atm2lnd_vars%timeres(v)))
    
            ierr = nf90_inq_varid(met_ncids(v), trim(metvars(v)), varid)
            !get the conversion factors
            ierr = nf90_get_att(met_ncids(v), varid, 'scale_factor', atm2lnd_vars%scale_factors(v))
            ierr = nf90_get_att(met_ncids(v), varid, 'add_offset', atm2lnd_vars%add_offsets(v))
            !get the met data         
            starti(1) = 1
            starti(2) = gtoget
            counti(1) = atm2lnd_vars%timelen_spinup(v)
            counti(2) = 1
            if (.not. const_climate_hist .and. (yr .ge. 1850 .or. use_sitedata)) counti(1) = atm2lnd_vars%timelen(v)

            if (i == 1 .and. v == 1)  then 
              allocate(atm2lnd_vars%atm_input       (met_nvars,bounds%begg:bounds%endg,1,1:counti(1)))
            end if 

            ierr = nf90_get_var(met_ncids(v), varid, atm2lnd_vars%atm_input(v,g:g,1,1:counti(1)), starti(1:2), counti(1:2))
            ierr = nf90_close(met_ncids(v))
    
            if (use_sitedata .and. v == 1) then 
                starti_site = max((nint(site_metdata(4,1))-atm2lnd_vars%startyear_met) * &
                                     365*nint(24./atm2lnd_vars%timeres(v))+1,1)
                endi_site   = (min(atm2lnd_vars%endyear_met_trans,nint(site_metdata(5,1))) - &
                                     atm2lnd_vars%startyear_met+1)*(365*nint(24./atm2lnd_vars%timeres(v)))
            end if
             
            atm2lnd_vars%var_offset(v,g,:) = 0._r8
            atm2lnd_vars%var_mult(v,g,:)   = 1._r8

            if (use_sitedata) then 
              !Compute monthly biases for site vs. reanalysis
              var_month_mean(:)  = 0._r8
              var_month_count(:) = 0
              do i=starti_site, endi_site
                thisdoy = mod(i,365*nint(24./atm2lnd_vars%timeres(v)))/(nint(24./atm2lnd_vars%timeres(v)))+1
                do m=1,12
                  if (thisdoy .ge. caldaym(m) .and. thisdoy .lt. caldaym(m+1)) thism = m
                enddo
                var_month_mean(thism) = var_month_mean(thism) + (atm2lnd_vars%atm_input(v,g,1,i)* &
                                          atm2lnd_vars%scale_factors(v) + atm2lnd_vars%add_offsets(v))
                var_month_count(thism) = var_month_count(thism)+1
              end do
     
              do m = 1,12
                var_month_mean(m) = var_month_mean(m)/var_month_count(m)
                !calculate offset and linear bias factors for temperature and precipitation
                if (v .eq. 1) atm2lnd_vars%var_offset(v,g,m) = (site_metdata(6,m)+SHR_CONST_TKFRZ) - var_month_mean(m)
                if (v .eq. 5 .and. var_month_mean(m) .gt. 0) &     
                      atm2lnd_vars%var_mult(v,g,m) = (site_metdata(7,m))/(caldaym(m+1)-caldaym(m))/24._r8/ &
                                                      3600._r8 / var_month_mean(m)
              end do
            end if
        
            !Align spinups and transient simulations
            !figure out which year to start with (assuming spinups always use integer multiple of met cycles)
            mystart = atm2lnd_vars%startyear_met
            do while (mystart > 1850)
              mystart = mystart - nyears_spinup
            end do
            if (atm2lnd_vars%metsource == 5) mystart=1850

            if (yr .lt. 1850) then 
              atm2lnd_vars%tindex(g,v,1) = (mod(yr-1,nyears_spinup) + (1850-mystart)) * 365 * nint(24./atm2lnd_vars%timeres(v))
            else if (yr .le. atm2lnd_vars%endyear_met_spinup) then
              atm2lnd_vars%tindex(g,v,1) = (mod(yr-1850,nyears_spinup) + (1850-mystart)) * 365 * nint(24./atm2lnd_vars%timeres(v))
            else
              atm2lnd_vars%tindex(g,v,1) = (yr - atm2lnd_vars%startyear_met) * 365 * nint(24./atm2lnd_vars%timeres(v))
            end if
            !adjust for starts not at beginning of year (but currently MUST begin at hour 0)
            atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1) + (caldaym(mon)+day-2)* &
                                         nint(24./atm2lnd_vars%timeres(v))
            
            atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,1) + 1
            if (atm2lnd_vars%tindex(g,v,1) == 0) then 
              atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen(v)
              if (yr .le. atm2lnd_vars%endyear_met_spinup) atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen_spinup(v)
             end if
          end do    !end variable loop        
        else
          do v=1,met_nvars
            if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then 
              if (v .eq. 4 .or. v .eq. 5 .or. (v .ge. 8 .and. v .le. 13)) then    !rad/Precipitation
                if (mod(tod/get_step_size(),nint(atm2lnd_vars%npf(v))) == 1 .and. nstep .gt. 3) then
                  atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+1
                  atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+1
                end if
              else  
                if (mod(tod/get_step_size()-1,nint(atm2lnd_vars%npf(v))) <= atm2lnd_vars%npf(v)/2._r8 .and. &
                    mod(tod/get_step_size(),nint(atm2lnd_vars%npf(v))) > atm2lnd_vars%npf(v)/2._r8) then 
                  atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+1
                  atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+1
                end if
              end if
            else
              atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+nint(1/atm2lnd_vars%npf(v))
              atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+nint(1/atm2lnd_vars%npf(v))  
            end if

            if (const_climate_hist .or. yr .le. atm2lnd_vars%startyear_met) then
              if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,1) = 1
              if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,2) = 1
            else if (yr .gt. atm2lnd_vars%endyear_met_trans) then
              if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen(v)) then
                 atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
              end if
              if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen(v)) then
                 atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
              end if
            end if

            !if (yr .gt. atm2lnd_vars%startyear_met) then 
            !  if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen(v)) atm2lnd_vars%tindex(g,v,1) = 1
            !  if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen(v)) atm2lnd_vars%tindex(g,v,2) = 1
            !else
            !  if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,1) = 1
            !  if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,2) = 1
            !end if
          end do
        end if

        tindex = atm2lnd_vars%tindex(g,:,:)

        !get weights for linear interpolation 
        do v=1,met_nvars
          if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then
               wt1(v) = 1._r8 - (mod((tod+86400)/get_step_size()-atm2lnd_vars%npf(v)/2._r8, &
                   atm2lnd_vars%npf(v))*1._r8)/atm2lnd_vars%npf(v)
               wt2(v) = 1._r8 - wt1(v)
          else
             wt1(v) = 0._r8    
             wt2(v) = 1._r8
          end if
        end do

        !Air temperature
        atm2lnd_vars%forc_t_not_downscaled_grc(g)  = min(((atm2lnd_vars%atm_input(1,g,1,tindex(1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(1,2))* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)             
        atm2lnd_vars%forc_th_not_downscaled_grc(g) = min(((atm2lnd_vars%atm_input(1,g,1,tindex(1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(1,2))* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)
       
        tbot = atm2lnd_vars%forc_t_not_downscaled_grc(g)

        !Air pressure
        atm2lnd_vars%forc_pbot_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(2,g,1,tindex(2,1))*atm2lnd_vars%scale_factors(2)+ &
                                                        atm2lnd_vars%add_offsets(2))*wt1(2) + (atm2lnd_vars%atm_input(2,g,1,tindex(2,2)) &
                                                        *atm2lnd_vars%scale_factors(2)+atm2lnd_vars%add_offsets(2))*wt2(2)) * &
                                                        atm2lnd_vars%var_mult(2,g,mon) + atm2lnd_vars%var_offset(2,g,mon), 4e4_r8)       
        !Specific humidity
        atm2lnd_vars%forc_q_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(3,g,1,tindex(3,1))*atm2lnd_vars%scale_factors(3)+ &
                                                     atm2lnd_vars%add_offsets(3))*wt1(3) + (atm2lnd_vars%atm_input(3,g,1,tindex(3,2)) &
                                                     *atm2lnd_vars%scale_factors(3)+atm2lnd_vars%add_offsets(3))*wt2(3)) * &
                                                     atm2lnd_vars%var_mult(3,g,mon) + atm2lnd_vars%var_offset(3,g,mon), 1e-9_r8)

        if (atm2lnd_vars%metsource == 2) then  !convert RH to qbot                             
          if (tbot > SHR_CONST_TKFRZ) then
            e = esatw(tdc(tbot))
          else
            e = esati(tdc(tbot))
          end if
          qsat           = 0.622_r8*e / (atm2lnd_vars%forc_pbot_not_downscaled_grc(g) - 0.378_r8*e)
          atm2lnd_vars%forc_q_not_downscaled_grc(g) = qsat * atm2lnd_vars%forc_q_not_downscaled_grc(g) / 100.0_r8
        end if

        !use longwave from file if provided
        atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ((atm2lnd_vars%atm_input(7,g,1,tindex(7,1))*atm2lnd_vars%scale_factors(7)+ &
                                                        atm2lnd_vars%add_offsets(7))*wt1(7) + (atm2lnd_vars%atm_input(7,g,1,tindex(7,2)) &
                                                        *atm2lnd_vars%scale_factors(7)+atm2lnd_vars%add_offsets(7))*wt2(7)) * &
                                                        atm2lnd_vars%var_mult(7,g,mon) + atm2lnd_vars%var_offset(7,g,mon)  
        if (atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .le. 50 .or. atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .ge. 600) then 
        !Longwave radiation (calculated from air temperature, humidity)
            e =  atm2lnd_vars%forc_pbot_not_downscaled_grc(g) * atm2lnd_vars%forc_q_not_downscaled_grc(g) / &
                 (0.622_R8 + 0.378_R8 * atm2lnd_vars%forc_q_not_downscaled_grc(g) )
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ea * SHR_CONST_STEBOL * tbot**4
        end if 

        !Shortwave radiation (cosine zenith angle interpolation)
        thishr = (tod-get_step_size()/2)/3600
        if (thishr < 0) thishr=thishr+24
        thismin = mod((tod-get_step_size()/2)/60, 60)
        thiscosz = max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr,thismin,0)* &
                        3.14159265358979/180.0d0), 0.001d0)
        avgcosz = 0d0
        if (atm2lnd_vars%npf(4) - 1._r8 .gt. 1e-3) then 
          swrad_period_len   = get_step_size()*nint(atm2lnd_vars%npf(4))
          swrad_period_start = ((tod-get_step_size()/2)/swrad_period_len) * swrad_period_len
          !set to last period if first model timestep of the day
          if (tod-get_step_size()/2 < 0) swrad_period_start = ((86400-get_step_size()/2)/swrad_period_len) * swrad_period_len   

          do tm=1,nint(atm2lnd_vars%npf(4))  
            !Get the average cosine zenith angle over the time resolution of the input data
            thishr  = (swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/3600
            if (thishr > 23) thishr=thishr-24  
            thismin = mod((swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/60, 60) 
            avgcosz  = avgcosz + max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr, thismin, 0) &
                       *3.14159265358979/180.0d0), 0.001d0)/atm2lnd_vars%npf(4)
          end do
        else
          avgcosz = thiscosz
        end if
        if (thiscosz > 0.001d0) then 
          wt2(4) = min(thiscosz/avgcosz, 10.0_r8)
        else
          wt2(4) = 0d0
        end if
        
        if (atm2lnd_vars%metsource == 5) then 
            wt2(4)=1.0   !cosz interp not working 
            wt2(8:10)=1.0
            swndf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                     atm2lnd_vars%add_offsets(4))*wt2(4)), 0.0_r8)
            swndr = max(((atm2lnd_vars%atm_input(8,g,1,tindex(8,2))*atm2lnd_vars%scale_factors(8)+ &
                                     atm2lnd_vars%add_offsets(8))*wt2(8)), 0.0_r8)
            swvdf = max(((atm2lnd_vars%atm_input(9,g,1,tindex(9,2))*atm2lnd_vars%scale_factors(9)+ &
                                     atm2lnd_vars%add_offsets(9))*wt2(9)), 0.0_r8)
            swvdr = max(((atm2lnd_vars%atm_input(10,g,1,tindex(10,2))*atm2lnd_vars%scale_factors(10)+ &
                                     atm2lnd_vars%add_offsets(10))*wt2(10)), 0.0_r8)
            atm2lnd_vars%forc_solad_grc(g,2) = swndr
            atm2lnd_vars%forc_solad_grc(g,1) = swvdr
            atm2lnd_vars%forc_solai_grc(g,2) = swndf
            atm2lnd_vars%forc_solai_grc(g,1) = swvdf
        else
            swndr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                     atm2lnd_vars%add_offsets(4))*wt2(4)) * 0.50_R8, 0.0_r8)
            swndf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr &
                           -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
            atm2lnd_vars%forc_solad_grc(g,2) = ratio_rvrf*swndr
            atm2lnd_vars%forc_solai_grc(g,2) = (1._R8 - ratio_rvrf)*swndf
            ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                               -9.0039e-06_R8*swvdr**2 +8.1351e-09_R8*swvdr**3,0.01_R8))
            atm2lnd_vars%forc_solad_grc(g,1) = ratio_rvrf*swvdr
            atm2lnd_vars%forc_solai_grc(g,1) = (1._R8 - ratio_rvrf)*swvdf
        end if
        !Rain and snow
        if (atm2lnd_vars%metsource == 5) then 
          forc_rainc = max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = max((((atm2lnd_vars%atm_input(11,g,1,tindex(11,2))*atm2lnd_vars%scale_factors(11)+ &
                                        atm2lnd_vars%add_offsets(11)))*atm2lnd_vars%var_mult(11,g,mon) + &
                                        atm2lnd_vars%var_offset(11,g,mon)), 0.0_r8)
          forc_snowc = max((((atm2lnd_vars%atm_input(12,g,1,tindex(12,2))*atm2lnd_vars%scale_factors(12)+ &
                                        atm2lnd_vars%add_offsets(12)))*atm2lnd_vars%var_mult(12,g,mon) + &
                                        atm2lnd_vars%var_offset(12,g,mon)), 0.0_r8)
          forc_snowl = max((((atm2lnd_vars%atm_input(13,g,1,tindex(13,2))*atm2lnd_vars%scale_factors(13)+ &
                                        atm2lnd_vars%add_offsets(13)))*atm2lnd_vars%var_mult(13,g,mon) + &
                                          atm2lnd_vars%var_offset(13,g,mon)), 0.0_r8)
        else
          frac = (atm2lnd_vars%forc_t_not_downscaled_grc(g) - SHR_CONST_TKFRZ)*0.5_R8       ! ramp near freezing
          frac = min(1.0_R8,max(0.0_R8,frac))           ! bound in [0,1]
          !Don't interpolate rainfall data
          forc_rainc = 0.1_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = 0.9_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                         atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                         atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8) 
          forc_snowc = 0.1_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)  
          forc_snowl = 0.9_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5))) * atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8) 
        end if
        !Wind
        atm2lnd_vars%forc_u_grc(g) = (atm2lnd_vars%atm_input(6,g,1,tindex(6,1))*atm2lnd_vars%scale_factors(6)+ &
                                     atm2lnd_vars%add_offsets(6))*wt1(6) + (atm2lnd_vars%atm_input(6,g,1,tindex(6,2))* &
                                     atm2lnd_vars%scale_factors(6)+atm2lnd_vars%add_offsets(6))*wt2(6)
        if (atm2lnd_vars%metsource == 5) then 
          atm2lnd_vars%forc_v_grc(g) = (atm2lnd_vars%atm_input(14,g,1,tindex(14,1))*atm2lnd_vars%scale_factors(14)+ &
                                     atm2lnd_vars%add_offsets(14))*wt1(14) + (atm2lnd_vars%atm_input(14,g,1,tindex(14,2))* &
                                     atm2lnd_vars%scale_factors(14)+atm2lnd_vars%add_offsets(14))*wt2(14)
        else
            atm2lnd_vars%forc_v_grc(g) = 0.0_R8 
        end if
        atm2lnd_vars%forc_hgt_grc(g) = 30.0_R8 !(atm2lnd_vars%atm_input(8,g,1,tindex(1))*wt1 + &
                                             !atm2lnd_vars%atm_input(8,g,1,tindex(2))*wt2)    ! zgcmxy  Atm state, default=30m

  !------------------------------------Fire data -------------------------------------------------------
 
        nindex(1) = yr-1848
        nindex(2) = nindex(1)+1
        if (yr .lt. 1850 .or. const_climate_hist) nindex(1:2) = 2
        if (yr .ge. 2010 .and. .not. const_climate_hist) nindex(1:2) = 161
      
        model_filter: if (use_cn .or. use_fates) then 
          if (atm2lnd_vars%loaded_bypassdata == 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then  
            if (masterproc .and. i .eq. 1) then 
              ! Read pop_dens streams namelist to get filename
              nu_nml = getavu()
              open(nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
              if (nml_error == 0) then
                  read(nu_nml, nml=popd_streams,iostat=nml_error)
                  if (nml_error /= 0) then
                      call endrun(msg='ERROR reading popdens namelist')
                  end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_popdens), NF90_NOWRITE, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lon)
              ierr = nf90_inq_varid(ncid, 'hdm', varid)
              starti(1:2) = 1 
              starti(3)   = nindex(1)
              counti(1) = 720
              counti(2) = 360
              counti(3) = 1       
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm1, starti, counti)
              starti(3) = nindex(2)
              if (nindex(1) .ne. nindex(2)) then 
                  ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm2, starti, counti)
              else
                  atm2lnd_vars%hdm2 = atm2lnd_vars%hdm1 
              end if
              ierr = nf90_close(ncid)
            end if

            if (i .eq. 1) then 
              call mpi_bcast (atm2lnd_vars%hdm1, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (atm2lnd_vars%hdm2, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lon, 720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lat, 360, MPI_REAL8, 0, mpicom, ier)
            end if
          end if

          !figure out which point to get
          if (atm2lnd_vars%loaded_bypassdata == 0) then 
            mindist=99999
            do thisx = 1,720
              do thisy = 1,360
                  if (ldomain%lonc(g) .lt. 0) then
                      if (smap05_lon(thisx) >= 180) smap05_lon(thisx) = smap05_lon(thisx)-360._r8
                  else if (ldomain%lonc(g) .ge. 180) then
                      if (smap05_lon(thisx) < 0) smap05_lon(thisx) = smap05_lon(thisx) + 360._r8
                  end if
                  thisdist = 100*((smap05_lat(thisy) - ldomain%latc(g))**2 + &
                          (smap05_lon(thisx) - ldomain%lonc(g))**2)**0.5
                  if (thisdist .lt. mindist) then
                      mindist = thisdist
                      atm2lnd_vars%hdmind(g,1) = thisx
                      atm2lnd_vars%hdmind(g,2) = thisy
                  end if
              end do
            end do
          end if
          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
          atm2lnd_vars%forc_hdm(g) = atm2lnd_vars%hdm1(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt1(1) + &
                                     atm2lnd_vars%hdm2(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt2(1)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. masterproc .and. i .eq. 1) then 
            ! Read light_streams namelist to get filename
            nu_nml = getavu()
            open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
            call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
            if (nml_error == 0) then
              read(nu_nml, nml=light_streams,iostat=nml_error)
              if (nml_error /= 0) then
                call endrun(msg='ERROR reading light_streams namelist')
              end if
            end if
            close(nu_nml)
            call relavu( nu_nml )

            !Get all of the data (master processor only)
            allocate(atm2lnd_vars%lnfm_all       (192,94,2920))
            ierr = nf90_open(trim(stream_fldFileName_lightng), NF90_NOWRITE, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lon)
            ierr = nf90_inq_varid(ncid, 'lnfm', varid)
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%lnfm_all)
            ierr = nf90_close(ncid)
          end if
          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. i .eq. 1) then
            call mpi_bcast (smapt62_lon, 192, MPI_REAL8, 0, mpicom, ier)
            call mpi_bcast (smapt62_lat, 94, MPI_REAL8, 0, mpicom, ier)
          end if
          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,192
              do thisy = 1,94
                if (ldomain%lonc(g) .lt. 0) then 
                  if (smapt62_lon(thisx) >= 180) smapt62_lon(thisx) = smapt62_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then 
                  if (smapt62_lon(thisx) < 0) smapt62_lon(thisx) = smapt62_lon(thisx) + 360._r8
                end if
                thisdist = 100*((smapt62_lat(thisy) - ldomain%latc(g))**2 + &
                            (smapt62_lon(thisx) - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  lnfmind(1) = thisx
                  lnfmind(2) = thisy
                end if
              end do
            end do
            if (masterproc) then
              atm2lnd_vars%lnfm(g,:) = atm2lnd_vars%lnfm_all(lnfmind(1),lnfmind(2),:)
              do np = 1,npes-1
                if (i == 1) then 
                  call mpi_recv(thisng,  1, MPI_INTEGER, np, 100000+np, mpicom, status, ier)
                  ng_all(np) = thisng
                end if
                if (i <= ng_all(np)) then 
                  call mpi_recv(lnfmind, 2, MPI_INTEGER, np, 200000+np, mpicom, status, ier)
                  call mpi_send(atm2lnd_vars%lnfm_all(lnfmind(1),lnfmind(2),:), 2920, &
                            MPI_REAL8, np, 300000+np, mpicom, ier)
                end if
              end do
            else
              if (i == 1)  call mpi_send(thisng,  1, MPI_INTEGER, 0, 100000+iam, mpicom, ier)
              call mpi_send(lnfmind, 2, MPI_INTEGER, 0, 200000+iam, mpicom, ier) 
              call mpi_recv(atm2lnd_vars%lnfm(g,:), 2920, MPI_REAL8, 0, 300000+iam, mpicom, status, ier)
            end if
          end if

          !Lightning data is 3-hourly.  Does not currently interpolate.
          atm2lnd_vars%forc_lnfm(g) = atm2lnd_vars%lnfm(g, ((int(thiscalday)-1)*8+tod/(3600*3))+1)

   !------------------------------------Nitrogen deposition----------------------------------------------

          !DMR note - ndep will NOT be correct if more than 1850 years of model
          !spinup (model year > 1850)
          nindex(1) = min(max(yr-1848,2), 168)
          nindex(2) = min(nindex(1)+1, 168)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then 
            if (masterproc .and. i .eq. 1) then 
              nu_nml = getavu()
              open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'ndepdyn_nml', status=nml_error)
              if (nml_error == 0) then
                read(nu_nml, nml=ndepdyn_nml,iostat=nml_error)
                if (nml_error /= 0) then
                  call endrun(msg='ERROR reading ndep namelist')
                end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_ndep), nf90_nowrite, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)      
              ierr = nf90_get_var(ncid, varid, smap2_lon)
              ierr = nf90_inq_varid(ncid, 'NDEP_year', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1)   = 144
              counti(2)   = 96
              counti(3)   = 1
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep1, starti, counti)
              if (nindex(1) .ne. nindex(2)) then 
                starti(3) = nindex(2)
                ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep2, starti, counti)
              else
                atm2lnd_vars%ndep2 = atm2lnd_vars%ndep1
              end if
              ierr = nf90_close(ncid)
             end if
             if (i .eq. 1) then
               call mpi_bcast (atm2lnd_vars%ndep1, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (atm2lnd_vars%ndep2, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lon, 144, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lat, 96, MPI_REAL8, 0, mpicom, ier)
             end if
          end if

          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then 
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then 
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then 
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
          end if

          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
  
          atm2lnd_vars%forc_ndep_grc(g)    = (atm2lnd_vars%ndep1(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt1(1) + &
                                              atm2lnd_vars%ndep2(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt2(1)) / (365._r8 * 86400._r8)
       end if model_filter

   !------------------------------------Aerosol forcing--------------------------------------------------
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then 
          if (masterproc .and. i .eq. 1) then 
            aerovars(1) = 'BCDEPWET'
            aerovars(2) = 'BCPHODRY'
            aerovars(3) = 'BCPHIDRY'
            aerovars(4) = 'OCDEPWET'
            aerovars(5) = 'OCPHODRY'
            aerovars(6) = 'OCPHIDRY'
            aerovars(7) = 'DSTX01DD'
            aerovars(8) = 'DSTX02DD'
            aerovars(9) = 'DSTX03DD'
            aerovars(10) = 'DSTX04DD'
            aerovars(11) = 'DSTX01WD'
            aerovars(12) = 'DSTX02WD'
            aerovars(13) = 'DSTX03WD'
            aerovars(14) = 'DSTX04WD'
            ierr = nf90_open(trim(aero_file), nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)      
            ierr = nf90_get_var(ncid, varid, smap2_lon)
            starti(1:2) = 1
            starti(3)   = max((min(yr,2100)-1849)*12+1, 13)-1
            counti(1)   = 144
            counti(2)   = 96
            counti(3)   = 14
            do av=1,14
              ierr = nf90_inq_varid(ncid, trim(aerovars(av)), varid)
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%aerodata(av,:,:,:), starti, counti)
            end do
            ierr = nf90_close(ncid)
          end if
          if (i .eq. 1) then 
             call mpi_bcast (atm2lnd_vars%aerodata, 14*144*96*14, MPI_REAL8, 0, mpicom, ier)
          end if
        end if

        !Use ndep grid indices since they're on the same grid
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. (.not. (use_fates .or. use_cn) )   ) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
        end if

        !get weights for interpolation (note this method doesn't get the month boundaries quite right..)
        aindex(1) = mon+1
        if (thiscalday .le. (caldaym(mon+1)+caldaym(mon))/2._r8) then 
           wt1(1) = 0.5_r8 + (thiscalday-caldaym(mon))/(caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)-1
        else
           wt1(1) = 1.0_r8 - (thiscalday-(caldaym(mon+1)+caldaym(mon))/2._r8)/   &
                          (caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)+1
        end if
        wt2(1) = 1._r8 - wt1(1)

        do av = 1,14
          atm2lnd_vars%forc_aer_grc(g,av)  =  atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(1))*wt1(1)+atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(2))*wt2(1)
        end do    

       !Parse startdate for adding temperature
       if (startdate_add_temperature .ne. '') then 
         call get_curr_date( yr, mon, day, tod )
         read(startdate_add_temperature,*) sdate_addt
         sy_addt     = sdate_addt/10000
         sm_addt     = (sdate_addt-sy_addt*10000)/100
         sd_addt     = sdate_addt-sy_addt*10000-sm_addt*100
         read(startdate_add_co2,*) sdate_addco2
         sy_addco2     = sdate_addco2/10000
         sm_addco2     = (sdate_addco2-sy_addco2*10000)/100
         sd_addco2     = sdate_addco2-sy_addco2*10000-sm_addt*100
       end if 
       if (startdate_add_temperature .ne. '') then
         if ((yr == sy_addt .and. mon == sm_addt .and. day >= sd_addt) .or. &
             (yr == sy_addt .and. mon > sm_addt) .or. (yr > sy_addt)) then
           atm2lnd_vars%forc_t_not_downscaled_grc(g) = atm2lnd_vars%forc_t_not_downscaled_grc(g) + add_temperature
           atm2lnd_vars%forc_th_not_downscaled_grc(g) = atm2lnd_vars%forc_th_not_downscaled_grc(g) + add_temperature
         end if
       end if

       !set the topounit-level atmospheric state and flux forcings (bypass mode)
       do topo = grc_pp%topi(g), grc_pp%topf(g)
         ! first, all the state forcings
         top_as%tbot(topo)    = atm2lnd_vars%forc_t_not_downscaled_grc(g)      ! forc_txy  Atm state K
         top_as%thbot(topo)   = atm2lnd_vars%forc_th_not_downscaled_grc(g)     ! forc_thxy Atm state K
         top_as%pbot(topo)    = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)   ! ptcmxy    Atm state Pa
         top_as%qbot(topo)    = atm2lnd_vars%forc_q_not_downscaled_grc(g)      ! forc_qxy  Atm state kg/kg
         top_as%ubot(topo)    = atm2lnd_vars%forc_u_grc(g)                     ! forc_uxy  Atm state m/s
         top_as%vbot(topo)    = atm2lnd_vars%forc_v_grc(g)                     ! forc_vxy  Atm state m/s
         top_as%zbot(topo)    = atm2lnd_vars%forc_hgt_grc(g)                   ! zgcmxy    Atm state m
         ! assign the state forcing fields derived from other inputs
         ! Horizontal windspeed (m/s)
         top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
         ! Relative humidity (percent)
         if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
            e = esatw(tdc(top_as%tbot(topo)))
         else
            e = esati(tdc(top_as%tbot(topo)))
         end if
         qsat           = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
         top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
         ! partial pressure of oxygen (Pa)
         top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
         ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
         vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
         top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))

         ! second, all the flux forcings
         top_af%rain(topo)    = forc_rainc + forc_rainl            ! sum of convective and large-scale rain
         top_af%snow(topo)    = forc_snowc + forc_snowl            ! sum of convective and large-scale snow
         top_af%solad(topo,2) = atm2lnd_vars%forc_solad_grc(g,2)   ! forc_sollxy  Atm flux  W/m^2
         top_af%solad(topo,1) = atm2lnd_vars%forc_solad_grc(g,1)   ! forc_solsxy  Atm flux  W/m^2
         top_af%solai(topo,2) = atm2lnd_vars%forc_solai_grc(g,2)   ! forc_solldxy Atm flux  W/m^2
         top_af%solai(topo,1) = atm2lnd_vars%forc_solai_grc(g,1)   ! forc_solsdxy Atm flux  W/m^2
         top_af%lwrad(topo)   = atm2lnd_vars%forc_lwrad_not_downscaled_grc(g)     ! flwdsxy Atm flux  W/m^2
         ! derived flux forcings
         top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
       end do
     
  !-----------------------------------------------------------------------------------------------------
#else

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
         
#endif

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

#ifdef CPL_BYPASS
       co2_type_idx = 2
#endif

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then
#ifdef CPL_BYPASS
        !atmospheric CO2 (to be used for transient simulations only)
        if (atm2lnd_vars%loaded_bypassdata .eq. 0) then 
          ierr = nf90_open(trim(co2_file), nf90_nowrite, ncid)
          ierr = nf90_inq_dimid(ncid, 'time', dimid)
          ierr = nf90_Inquire_Dimension(ncid, dimid, len = thistimelen)
          ierr = nf90_inq_varid(ncid, 'CO2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%co2_input(:,:,1:thistimelen))
          ierr = nf90_inq_varid(ncid, 'C13O2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%c13o2_input(:,:,1:thistimelen))
          ierr = nf90_close(ncid)
        end if

        !get weights/indices for interpolation (assume values represent annual averages)
        nindex(1) = min(max(yr,1850),2100)-1764
        if (thiscalday .le. 182.5) then 
          nindex(2) = nindex(1)-1  
        else
          nindex(2) = nindex(1)+1
        end if
        wt1(1) = 1._r8 - abs((182.5 - (thiscalday -1._r8))/365._r8)
        wt2(1) = 1._r8 - wt1(1)

        co2_ppmv_val = atm2lnd_vars%co2_input(1,1,nindex(1))*wt1(1) + atm2lnd_vars%co2_input(1,1,nindex(2))*wt2(1)
        if (startdate_add_co2 .ne. '') then
          if ((yr == sy_addco2 .and. mon == sm_addco2 .and. day >= sd_addco2) .or. &
              (yr == sy_addco2 .and. mon > sm_addco2) .or. (yr > sy_addco2)) then
            co2_ppmv_val=co2_ppmv_val + add_co2
          end if
        end if

        if (use_c13) then 
          atm2lnd_vars%forc_pc13o2_grc(g) = (atm2lnd_vars%c13o2_input(1,1,nindex(1))*wt1(1) + &
               atm2lnd_vars%c13o2_input(1,1,nindex(2))*wt2(1)) * 1.e-6_r8 * forc_pbot
        end if
        co2_type_idx = 1
#else
          co2_ppmv_val = co2_ppmv_diag 
           if (use_c13) then
             atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
           end if
#endif
       else
          co2_ppmv_val = co2_ppmv
          if (use_c13) then
            atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
          end if
       end if
       atm2lnd_vars%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot 

#ifdef CPL_BYPASS
       do topo = grc_pp%topi(g), grc_pp%topf(g)
         top_as%pco2bot(topo) = atm2lnd_vars%forc_pco2_grc(g)
         if (use_c13) then
            top_as%pc13o2bot(topo) = atm2lnd_vars%forc_pc13o2_grc(g)
         end if
       end do
#endif
      
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
#ifdef CPL_BYPASS
    atm2lnd_vars%loaded_bypassdata = 1
#endif

  end subroutine lnd_import_moab

! endif for ifdef HAVE_MOAB
#endif

end module lnd_comp_mct
