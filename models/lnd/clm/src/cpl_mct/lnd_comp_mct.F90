module lnd_comp_mct
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_mct
!
!  Interface of the active land model component of CESM the CLM (Community Land Model)
!  with the main CESM driver. This is a thin interface taking CESM driver information
!  in MCT (Model Coupling Toolkit) format and converting it to use by CLM.
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_sys_mod      , only : shr_sys_flush
  use mct_mod          , only : mct_aVect, mct_gsmap
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: lnd_init_mct               ! clm initialization
  public :: lnd_run_mct                ! clm run phase
  public :: lnd_final_mct              ! clm finalization/cleanup
!
! !PUBLIC DATA MEMBERS: None
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! Dec/18/2009 Make sure subroutines have documentation. Erik Kluzek
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_SetgsMap_mct         ! Set the land model MCT GS map
  private :: lnd_domain_mct           ! Set the land model domain information
  private :: lnd_export_mct           ! export land data to CESM coupler
  private :: lnd_import_mct           ! import data from the CESM coupler to the land model
  private :: sno_export_mct
  private :: sno_import_mct
!
! !PRIVATE DATA MEMBERS:
!
! Time averaged flux fields
!  
  type(mct_aVect)   :: l2x_l_SNAP     ! Snapshot of land to coupler data on the land grid
  type(mct_aVect)   :: l2x_l_SUM      ! Summation of land to coupler data on the land grid

  type(mct_aVect)   :: s2x_s_SNAP     ! Snapshot of sno to coupler data on the land grid
  type(mct_aVect)   :: s2x_s_SUM      ! Summation of sno to coupler data on the land grid 
!
! Time averaged counter for flux fields
!
  integer :: avg_count                ! Number of times snapshots of above flux data summed together
  integer :: avg_count_sno
!
! Atmospheric mode  
!
  logical :: atm_prognostic           ! Flag if active atmosphere component or not

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_mct
!
! !INTERFACE:
  subroutine lnd_init_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                   cdata_s, x2s_s, s2x_s, &
                                   NLFilename )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use abortutils       , only : endrun
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, &
                                  set_nextsw_cday
    use clm_atmlnd       , only : clm_l2a
    use clm_glclnd       , only : clm_s2x
    use clm_initializeMod, only : initialize1, initialize2
    use clm_varctl       , only : finidat,single_column, set_clmvarctl, iulog, noland, &
                                  inst_index, inst_suffix, inst_name, &
                                  create_glacier_mec_landunit 
    use clm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use controlMod       , only : control_setNL
    use decompMod        , only : get_proc_bounds
    use domainMod        , only : ldomain
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                  shr_file_getLogUnit, shr_file_getLogLevel, &
                                  shr_file_getUnit, shr_file_setIO
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use spmdMod          , only : masterproc, spmd_init
    use clm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use clm_cpl_indices  , only : clm_cpl_indices_set, nflds_l2x
    use seq_flds_mod
    use mct_mod
    use ESMF
    implicit none
!
! !ARGUMENTS:
    type(ESMF_Clock),           intent(in)    :: EClock           ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_l          ! Input land-model driver data
    type(mct_aVect),            intent(inout) :: x2l_l, l2x_l     ! land model import and export states
    type(seq_cdata),            intent(inout) :: cdata_s          ! Input snow-model (land-ice) driver data
    type(mct_aVect),            intent(inout) :: x2s_s, s2x_s     ! Snow-model import and export states
    character(len=*), optional, intent(in)    :: NLFilename       ! Namelist filename to read
!
! !LOCAL VARIABLES:
    integer                          :: LNDID	     ! Land identifyer
    integer                          :: mpicom_lnd   ! MPI communicator
    type(mct_gsMap),         pointer :: GSMap_lnd    ! Land model MCT GS map
    type(mct_gGrid),         pointer :: dom_l        ! Land model domain
    type(mct_gsMap),         pointer :: GSMap_sno
    type(mct_gGrid),         pointer :: dom_s
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer  :: lsize                                ! size of attribute vector
    integer  :: g,i,j                                ! indices
    integer  :: dtime_sync                           ! coupling time-step from the input synchronization clock
    integer  :: dtime_clm                            ! clm time-step
    logical  :: exists                               ! true if file exists
    logical  :: atm_aero                             ! Flag if aerosol data sent from atm model
    logical  :: samegrid_al                          ! true if atmosphere and land are on the same grid
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
    integer :: nsrest                                ! clm restart type
    integer :: perpetual_ymd                         ! perpetual date
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    logical :: perpetual_run                         ! flag if should cycle over a perpetual date or not
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: begg, endg
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    ! Determine attriute vector indices

    call clm_cpl_indices_set()

    ! Initialize clm MPI communicator 

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
       write(iulog,format) "CLM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Use infodata to set orbital values

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Consistency check on namelist filename	

    call control_setNL("lnd_in"//trim(inst_suffix))

    ! Initialize clm
    ! initialize1 reads namelist, grid and surface data (need this to initialize gsmap) 
    ! initialize2 performs rest of initialization	

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, perpetual=perpetual_run,                &
                              perpetual_ymd=perpetual_ymd, case_name=caseid,    &
                              case_desc=ctitle, single_column=single_column,    &
                              scmlat=scmlat, scmlon=scmlon,                     &
                              brnch_retain_casename=brnch_retain_casename,      &
                              start_type=starttype, model_version=version,      &
                              hostname=hostname, username=username,             &
                              samegrid_al=samegrid_al                           &
                                )
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod,  perpetual_run_in=perpetual_run,                &
                           perpetual_ymd_in=perpetual_ymd )
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call set_clmvarctl(    caseid_in=caseid, ctitle_in=ctitle,                     &
                           brnch_retain_casename_in=brnch_retain_casename,         &
                           single_column_in=single_column, scmlat_in=scmlat,       &
                           scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                           hostname_in=hostname, username_in=username)

    ! Read namelist, grid and surface data

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland ) then
           call seq_infodata_PutData( infodata, sno_present   =.false.)
           call seq_infodata_PutData( infodata, lnd_present   =.false.)
           call seq_infodata_PutData( infodata, lnd_prognostic=.false.)
       return
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    call seq_infodata_GetData(infodata, atm_aero=atm_aero )
    if ( .not. atm_aero )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to CLM' )
    end if

    ! Initialize lnd gsMap and domain

    call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd ) 	
    lsize = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)
    call lnd_domain_mct( lsize, gsMap_lnd, dom_l )

    ! Initialize lnd attribute vectors coming from driver

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_l)

    call mct_aVect_init(l2x_l_SNAP, rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SNAP)

    call mct_aVect_init(l2x_l_SUM , rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SUM )

    if (masterproc) then
       write(iulog,format)'time averaging the following flux fields over the coupling interval'
       write(iulog,format) trim(seq_flds_l2x_fluxes)
    end if

    ! Finish initializing clm

    call initialize2()

    ! Check that clm internal dtime aligns with clm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_clm = get_step_size()
    if (masterproc) write(iulog,*)'dtime_sync= ',dtime_sync,&
         ' dtime_clm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'clm dtime ',dtime_clm,' and Eclock dtime ',&
            dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state 

    call get_proc_bounds(begg, endg) 
    call lnd_export_mct( clm_l2a, l2x_l, begg, endg )

    if (create_glacier_mec_landunit) then
       call seq_cdata_setptrs(cdata_s, gsMap=gsMap_sno, dom=dom_s)

       ! Initialize sno gsMap (same as gsMap_lnd)
       call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_sno )
       lsize = mct_gsMap_lsize(gsMap_sno, mpicom_lnd)

       ! Initialize sno domain (same as lnd domain)
       call lnd_domain_mct( lsize, gsMap_sno, dom_s )

       ! Initialize sno attribute vectors
       call mct_aVect_init(x2s_s, rList=seq_flds_x2s_fields, lsize=lsize)
       call mct_aVect_zero(x2s_s)

       call mct_aVect_init(s2x_s, rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s)

       ! In contrast to l2x_l_SNAP / l2x_l_SUM, for s2x we accumulate/average all fields,
       ! not just fluxes. This is because glc wants the time-averaged tsrf field (and the
       ! other state field, topo, is not time-varying, so it doesn't matter what we do
       ! with that field)
       call mct_aVect_init(s2x_s_SUM , rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s_SUM )

       call mct_aVect_init(s2x_s_SNAP , rList=seq_flds_s2x_fields, lsize=lsize)
       call mct_aVect_zero(s2x_s_SNAP )

       ! Create mct sno export state
       call sno_export_mct(clm_s2x, s2x_s)
    endif   ! create_glacier_mec_landunit

    ! Initialize averaging counter

    avg_count = 0

    ! Fill in infodata

    call seq_infodata_PutData( infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData( infodata, lnd_nx=ldomain%ni, lnd_ny=ldomain%nj)
    if (create_glacier_mec_landunit) then
       call seq_infodata_PutData( infodata, sno_present=.true.)
       call seq_infodata_PutData( infodata, sno_prognostic=.false.)
       call seq_infodata_PutData( infodata, sno_nx=ldomain%ni, sno_ny=ldomain%nj)
    else
       call seq_infodata_PutData( infodata, sno_present=.false.)
       call seq_infodata_PutData( infodata, sno_prognostic=.false.)
    endif

    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )

    ! Determine atmosphere modes

    call seq_infodata_GetData(infodata, atm_prognostic=atm_prognostic)
    if (masterproc) then
       if ( atm_prognostic )then
          write(iulog,format) 'Atmospheric input is from a prognostic model'
       else
          write(iulog,format) 'Atmospheric input is from a data model'
       end if
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

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_mct
!
! !INTERFACE:
  subroutine lnd_run_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                  cdata_s, x2s_s, s2x_s)
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use shr_kind_mod    ,only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd      ,only : clm_l2a, clm_a2l
    use clm_driver      ,only : clm_drv
    use clm_time_manager,only : get_curr_date, get_nstep, get_curr_calday, get_step_size, &
                                advance_timestep, set_nextsw_cday,update_rad_dtime
    use decompMod       ,only : get_proc_bounds
    use abortutils      ,only : endrun
    use clm_varctl      ,only : iulog, create_glacier_mec_landunit 
    use clm_varorb      ,only : eccen, obliqr, lambm0, mvelpp
    use shr_file_mod    ,only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel
    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use seq_infodata_mod,only : seq_infodata_type, seq_infodata_GetData
    use spmdMod         ,only : masterproc, mpicom
    use perf_mod        ,only : t_startf, t_stopf, t_barrierf
    use clm_glclnd      ,only : clm_s2x, clm_x2s, unpack_clm_x2s
    use shr_orb_mod     ,only : shr_orb_decl
    use clm_varorb      ,only : eccen, mvelpp, lambm0, obliqr
    use clm_cpl_indices ,only : nflds_l2x, nflds_x2l
    use mct_mod
    use ESMF
    implicit none
!
! !ARGUMENTS:
    type(ESMF_Clock) , intent(in)    :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    type(seq_cdata)  , intent(in)    :: cdata_s   ! Input driver data for snow model (land-ice)
    type(mct_aVect)  , intent(inout) :: x2s_s     ! Import state for snow model
    type(mct_aVect)  , intent(inout) :: s2x_s     ! Export state for snow model
!
! !LOCAL VARIABLES:
    integer :: ymd_sync                   ! Sync date (YYYYMMDD)
    integer :: yr_sync                    ! Sync current year
    integer :: mon_sync                   ! Sync current month
    integer :: day_sync                   ! Sync current day
    integer :: tod_sync                   ! Sync current time of day (sec)
    integer :: ymd                        ! CLM current date (YYYYMMDD)
    integer :: yr                         ! CLM current year
    integer :: mon                        ! CLM current month
    integer :: day                        ! CLM current day
    integer :: tod                        ! CLM current time of day (sec)
    integer :: dtime                      ! time step increment (sec)
    integer :: nstep                      ! time step index
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    logical :: dosend                     ! true => send data back to driver
    logical :: doalb                      ! .true. ==> do albedo calculation on this time step
    real(r8):: nextsw_cday                ! calday from clock of next radiation computation
    real(r8):: caldayp1                   ! clm calday plus dtime offset
    integer :: shrlogunit,shrloglev       ! old values for share log unit and log level
    integer :: begg, endg                 ! Beginning and ending gridcell index numbers
    integer :: lbnum                      ! input to memory diagnostic
    type(seq_infodata_type),pointer :: infodata ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_l    ! Land model domain data
    integer  :: g,i,lsize                       ! counters
    logical,save :: first_call = .true.         ! first call work
    logical  :: glcrun_alarm          ! if true, sno data is averaged and sent to glc this step
    logical  :: update_glc2sno_fields ! if true, update glacier_mec fields
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinp1              ! solar declination angle in radians for nstep+1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    real(r8) :: recip                 ! reciprical
    character(len=32)            :: rdate       ! date char string for restart file names
    character(len=32), parameter :: sub = "lnd_run_mct"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

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

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    call get_proc_bounds(begg, endg)

    ! Map MCT to land data type
    ! Perform downscaling if appropriate

    call t_startf ('lc_lnd_import')
    call lnd_import_mct( x2l_l, clm_a2l, begg, endg )
    
    ! Map to clm (only when state and/or fluxes need to be updated)

    if (create_glacier_mec_landunit) then
       update_glc2sno_fields  = .false.
       call seq_infodata_GetData(infodata, glc_g2supdate = update_glc2sno_fields)
       if (update_glc2sno_fields) then
          call sno_import_mct( x2s_s, clm_x2s )
          call unpack_clm_x2s(clm_x2s)
       endif ! update_glc2sno
    endif ! create_glacier_mec_landunit
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

       ! Run clm 

       call t_barrierf('sync_clm_run1', mpicom)
       call t_startf ('clm_run')
       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')
       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('clm_run')

       ! Create l2x_l export state - add river runoff input to l2x_l if appropriate
       
       call t_startf ('lc_lnd_export')
       call lnd_export_mct( clm_l2a, l2x_l, begg, endg )
       call t_stopf ('lc_lnd_export')

       ! Do not accumulate on first coupling freq - consistency with ccsm3

       nstep = get_nstep()
       if (nstep <= 1) then
          call mct_aVect_copy( l2x_l, l2x_l_SUM )
          avg_count = 1
       else
          call mct_aVect_copy( l2x_l, l2x_l_SNAP )
          call mct_aVect_accum( aVin=l2x_l_SNAP, aVout=l2x_l_SUM )
          avg_count = avg_count + 1
       endif
       
       ! Map sno data type to MCT

       if (create_glacier_mec_landunit) then
          call sno_export_mct(clm_s2x, s2x_s)
          if (nstep <= 1) then
             call mct_aVect_copy( s2x_s, s2x_s_SUM )
             avg_count_sno = 1
          else
             call mct_aVect_copy( s2x_s, s2x_s_SNAP )
             call mct_aVect_accum( aVin=s2x_s_SNAP, aVout=s2x_s_SUM )
             avg_count_sno = avg_count_sno + 1
          endif
       endif    ! create_glacier_mec_landunit

       ! Advance clm time step
       
       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')

    end do

    ! Finish accumulation of attribute vector and average and zero out partial sum and counter
    
    if (avg_count /= 0) then
       recip = 1.0_r8/(real(avg_count,r8))
       l2x_l_SUM%rAttr(:,:) = l2x_l_SUM%rAttr(:,:) * recip
    endif
    call mct_aVect_copy( l2x_l_SUM, l2x_l )
    call mct_aVect_zero( l2x_l_SUM) 
    avg_count = 0                   

    if (create_glacier_mec_landunit) then
       call seq_infodata_GetData(infodata, glcrun_alarm = glcrun_alarm )
       if (glcrun_alarm) then
          if (avg_count_sno /= 0) then
             recip = 1.0_r8/(real(avg_count_sno,r8))
             s2x_s_SUM%rAttr(:,:) = s2x_s_SUM%rAttr(:,:) * recip
          endif
          call mct_aVect_copy( s2x_s_SUM, s2x_s )
          call mct_aVect_zero( s2x_s_SUM)
          avg_count_sno = 0
       endif
    endif

    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' clm ymd=',ymd     ,'  clm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: CLM clock not in sync with Master Sync clock" )
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

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_mct
!
! !INTERFACE:
  subroutine lnd_final_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                    cdata_s, x2s_s, s2x_s )
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!
    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use mct_mod
    use esmf
   implicit none
! !ARGUMENTS:
    type(ESMF_Clock) , intent(in)    :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(in)    :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    type(seq_cdata)  , intent(in)    :: cdata_s   ! Input driver data for snow model (land-ice)
    type(mct_aVect)  , intent(inout) :: x2s_s     ! Import state for snow model
    type(mct_aVect)  , intent(inout) :: s2x_s     ! Export state for snow model
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
  end subroutine lnd_final_mct

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_SetgsMap_mct
!
! !INTERFACE:
  subroutine lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd )
!-------------------------------------------------------------------
!
! !DESCRIPTION:
!
! Set the MCT GS map for the land model
!
!-------------------------------------------------------------------
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use decompMod    , only : get_proc_bounds, ldecomp
    use domainMod    , only : ldomain
    use mct_mod      , only : mct_gsMap, mct_gsMap_init
    implicit none
! !ARGUMENTS:
    integer        , intent(in)  :: mpicom_lnd    ! MPI communicator for the clm land model
    integer        , intent(in)  :: LNDID         ! Land model identifyer number
    type(mct_gsMap), intent(out) :: gsMap_lnd     ! Resulting MCT GS map for the land model
!
! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: i, j, n, gi            ! Indices
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code
    integer :: begg, endg             ! Beginning/Ending grid cell index
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    call get_proc_bounds(begg, endg)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = endg-begg+1
    gsize = ldomain%ni * ldomain%nj

    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_export_mct
!
! !INTERFACE:
  subroutine lnd_export_mct( clm_l2a, l2x_l, begg, endg )   
!
! !DESCRIPTION:
!
! Convert the data to be sent from the clm model to the coupler from clm data types
! to MCT data types.
! 
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use clm_varctl         , only : iulog
    use clm_time_manager   , only : get_nstep, get_step_size  
    use clm_atmlnd         , only : lnd2atm_type
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use clm_cpl_indices
    use clmtype
    implicit none
! !ARGUMENTS:
    type(lnd2atm_type), intent(inout) :: clm_l2a    ! clm land to atmosphere exchange data type
    type(mct_aVect)   , intent(inout) :: l2x_l      ! Land to coupler export state on land grid
    integer           , intent(in)    :: begg       ! beginning grid cell index
    integer           , intent(in)    :: endg       ! ending grid cell index
!
! !LOCAL VARIABLES:
    integer  :: g,i                           ! indices
    integer  :: ier                           ! error status
    integer  :: nstep                         ! time step index
    integer  :: dtime                         ! time step   
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    
    ! cesm sign convention is that fluxes are positive downward

    l2x_l%rAttr(:,:) = 0.0_r8

    do g = begg,endg
       i = 1 + (g-begg)
       l2x_l%rAttr(index_l2x_Sl_t,i)        =  clm_l2a%t_rad(g)
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    =  clm_l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    =  clm_l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    =  clm_l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    =  clm_l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    =  clm_l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Sl_tref,i)     =  clm_l2a%t_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_qref,i)     =  clm_l2a%q_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_u10,i)      =  clm_l2a%u_ref10m(g)
       l2x_l%rAttr(index_l2x_Fall_taux,i)   = -clm_l2a%taux(g)
       l2x_l%rAttr(index_l2x_Fall_tauy,i)   = -clm_l2a%tauy(g)
       l2x_l%rAttr(index_l2x_Fall_lat,i)    = -clm_l2a%eflx_lh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_sen,i)    = -clm_l2a%eflx_sh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   = -clm_l2a%eflx_lwrad_out(g)
       l2x_l%rAttr(index_l2x_Fall_evap,i)   = -clm_l2a%qflx_evap_tot(g)
       l2x_l%rAttr(index_l2x_Fall_swnet,i)  =  clm_l2a%fsa(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x_l%rAttr(index_l2x_Fall_fco2_lnd,i) = -clm_l2a%nee(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x_l%rAttr(index_l2x_Sl_ram1,i) = clm_l2a%ram1(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x_l%rAttr(index_l2x_Sl_fv,i)   = clm_l2a%fv(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x_l%rAttr(index_l2x_Sl_soilw,i)   = clm_l2a%h2osoi_vol(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst1,i)= -clm_l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst2,i)= -clm_l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst3,i)= -clm_l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst4,i)= -clm_l2a%flxdst(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x_l%rAttr(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               clm_l2a%ddvel(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x_l%rAttr(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -clm_l2a%flxvoc(g,:shr_megan_mechcomps_n)
       end if

#ifdef LCH4
       if (index_l2x_Fall_methane /= 0) then
          l2x_l%rAttr(index_l2x_Fall_methane,i) = -clm_l2a%flux_ch4(g) 
       endif
#endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x_l%rattr(index_l2x_Flrl_rofliq,i) = clm_l2a%rofliq(g)
       l2x_l%rattr(index_l2x_Flrl_rofice,i) = clm_l2a%rofice(g)

    end do

  end subroutine lnd_export_mct

!====================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_import_mct
!
! !INTERFACE:
  subroutine lnd_import_mct( x2l_l, a2l, begg, endg )
!
! !DESCRIPTION:
!
! Convert the input data from the coupler to the land model from MCT import state
! into internal clm data types.
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13
    use clm_varcon      , only: rair, o2_molar_const
    use clm_varcon      , only: c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use abortutils      , only: endrun
    use mct_mod         , only: mct_aVect
    use clm_cpl_indices
    use clmtype
    implicit none
! !ARGUMENTS:
    type(mct_aVect)   , intent(inout) :: x2l_l   ! Driver MCT import state to land model
    type(atm2lnd_type), intent(inout) :: a2l     ! clm internal input data type
    integer           , intent(in)    :: begg	
    integer           , intent(in)    :: endg
!
! !LOCAL VARIABLES:
    integer  :: g,i,nstep,ier        ! indices, number of steps, and error code
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
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
    character(len=32), parameter :: sub = 'lnd_import_mct'

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
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 27 February 2008: Keith Oleson; Forcing height change
!
!EOP
!---------------------------------------------------------------------------

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
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    

    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine flooding input, sign convention is positive downward and
        ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
        ! change the sign to indicate addition of water to system.

        a2l%forc_flood(g)   = -x2l_l%rattr(index_x2l_Flrr_flood,i)  

        a2l%volr(g)   = x2l_l%rattr(index_x2l_Slrr_volr,i)

        ! Determine required receive fields

        a2l%forc_hgt(g)     = x2l_l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = x2l_l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = x2l_l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = x2l_l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = x2l_l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = x2l_l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = x2l_l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = x2l_l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l_l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l_l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l_l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l_l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! atmosphere coupling, for prognostic/prescribed aerosols
        a2l%forc_aer(g,1)  =  x2l_l%rAttr(index_x2l_Faxa_bcphidry,i)
        a2l%forc_aer(g,2)  =  x2l_l%rAttr(index_x2l_Faxa_bcphodry,i)
        a2l%forc_aer(g,3)  =  x2l_l%rAttr(index_x2l_Faxa_bcphiwet,i)
        a2l%forc_aer(g,4)  =  x2l_l%rAttr(index_x2l_Faxa_ocphidry,i)
        a2l%forc_aer(g,5)  =  x2l_l%rAttr(index_x2l_Faxa_ocphodry,i)
        a2l%forc_aer(g,6)  =  x2l_l%rAttr(index_x2l_Faxa_ocphiwet,i)
        a2l%forc_aer(g,7)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet1,i)
        a2l%forc_aer(g,8)  =  x2l_l%rAttr(index_x2l_Faxa_dstdry1,i)
        a2l%forc_aer(g,9)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet2,i)
        a2l%forc_aer(g,10) =  x2l_l%rAttr(index_x2l_Faxa_dstdry2,i)
        a2l%forc_aer(g,11) =  x2l_l%rAttr(index_x2l_Faxa_dstwet3,i)
        a2l%forc_aer(g,12) =  x2l_l%rAttr(index_x2l_Faxa_dstdry3,i)
        a2l%forc_aer(g,13) =  x2l_l%rAttr(index_x2l_Faxa_dstwet4,i)
        a2l%forc_aer(g,14) =  x2l_l%rAttr(index_x2l_Faxa_dstdry4,i)

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l_l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l_l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv
        end if

#ifdef LCH4
        if (index_x2l_Sa_methane /= 0) then
           a2l%forc_pch4(g) = x2l_l%rAttr(index_x2l_Sa_methane,i)
        endif
#endif

        ! Determine derived quantities for required fields
        a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
        a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
        a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
        a2l%forc_vp(g)    = a2l%forc_q(g) * a2l%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * a2l%forc_q(g))
        a2l%forc_rho(g)   = (a2l%forc_pbot(g) - 0.378_r8 * a2l%forc_vp(g)) &
                            / (rair * a2l%forc_t(g))
        a2l%forc_po2(g)   = o2_molar_const * a2l%forc_pbot(g)
        a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
        a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
                            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)
        a2l%forc_rain(g)  = forc_rainc + forc_rainl
        a2l%forc_snow(g)  = forc_snowc + forc_snowl
        a2l%rainf    (g)  = a2l%forc_rain(g) + a2l%forc_snow(g)

        if (a2l%forc_t(g) > SHR_CONST_TKFRZ) then
           e = esatw(tdc(a2l%forc_t(g)))
        else
           e = esati(tdc(a2l%forc_t(g)))
        end if
        qsat           = 0.622_r8*e / (a2l%forc_pbot(g) - 0.378_r8*e)
        a2l%forc_rh(g) = 100.0_r8*(a2l%forc_q(g) / qsat)
        ! Make sure relative humidity is properly bounded
        ! a2l%forc_rh(g) = min( 100.0_r8, a2l%forc_rh(g) )
        ! a2l%forc_rh(g) = max(   0.0_r8, a2l%forc_rh(g) )
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv_val = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv_val = co2_ppmv_diag 
        else
           co2_ppmv_val = co2_ppmv
        end if
        a2l%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * a2l%forc_pbot(g) 
        if (use_c13) then
           a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
        end if
     end do

   end subroutine lnd_import_mct

!===============================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_domain_mct
!
! !INTERFACE:
  subroutine lnd_domain_mct( lsize, gsMap_l, dom_l )
!
! !DESCRIPTION:
!
! Send the land model domain information to the coupler
!
!---------------------------------------------------------------------------
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use domainMod   , only : ldomain
    use decompMod   , only : get_proc_bounds
    use spmdMod     , only : iam
    use mct_mod     , only : mct_gsMap, mct_gGrid, mct_gGrid_importIAttr, &
                             mct_gGrid_importRAttr, mct_gGrid_init,       &
                             mct_gsMap_orderedPoints
    use seq_flds_mod
    implicit none
! !ARGUMENTS:
    integer        , intent(in)    :: lsize     ! land model domain data size
    type(mct_gsMap), intent(inout) :: gsMap_l   ! Output land model MCT GS map
    type(mct_ggrid), intent(out)   :: dom_l     ! Output domain information for land model
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local Variables
    !
    integer :: g,i,j              ! index
    integer :: begg, endg         ! beginning and ending gridcell indices
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_l, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    call get_proc_bounds(begg, endg)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct
    
!===============================================================================

  subroutine sno_export_mct( s2x, s2x_s )   

    use clm_glclnd      , only : lnd2glc_type
    use decompMod       , only : get_proc_bounds
    use clm_cpl_indices 
    use clm_varctl       , only : iulog

    type(lnd2glc_type), intent(inout) :: s2x
    type(mct_aVect)   , intent(inout) :: s2x_s

    integer :: g,i,num
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------

    call get_proc_bounds(begg, endg)

    ! qice is positive if ice is growing, negative if melting

    s2x_s%rAttr(:,:) = 0.0_r8
    do g = begg,endg
       i = 1 + (g-begg)
       do num = 1,glc_nec
          s2x_s%rAttr(index_s2x_Ss_tsrf(num),i)   = s2x%tsrf(g,num)
          s2x_s%rAttr(index_s2x_Ss_topo(num),i)   = s2x%topo(g,num)
          s2x_s%rAttr(index_s2x_Fgss_qice(num),i) = s2x%qice(g,num)
       end do
    end do  

  end subroutine sno_export_mct

!====================================================================================

  subroutine sno_import_mct( x2s_s, x2s )

    use clm_glclnd      , only: glc2lnd_type
    use decompMod       , only: get_proc_bounds
    use mct_mod         , only: mct_aVect
    use clm_cpl_indices
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2s_s
    type(glc2lnd_type), intent(inout) :: x2s
    !
    ! Local Variables
    !
    integer  :: g,i,num
    integer  :: begg, endg   ! beginning and ending gridcell indices
    !-----------------------------------------------------

    call get_proc_bounds(begg, endg)

    do g = begg,endg
       i = 1 + (g - begg)
       do num = 1,glc_nec
          x2s%frac(g,num)  = x2s_s%rAttr(index_x2s_Sg_frac(num),i)
          x2s%topo(g,num)  = x2s_s%rAttr(index_x2s_Sg_topo(num),i)
          x2s%hflx(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_hflx(num),i)
          x2s%rofi(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_rofi(num),i)
          x2s%rofl(g,num)  = x2s_s%rAttr(index_x2s_Fsgg_rofl(num),i)
       end do
     end do  

   end subroutine sno_import_mct

!====================================================================================

end module lnd_comp_mct
