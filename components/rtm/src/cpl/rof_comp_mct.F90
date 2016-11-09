module rof_comp_mct
  
!========================================================================
! DESCRIPTION:
! Interface of the active runoff component of CESM 
! with the main CESM driver. This is a thin interface taking CESM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by RTM

  use seq_flds_mod
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
  use seq_timemgr_mod  , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
  use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                seq_infodata_start_type_brnch
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use RunoffMod        , only : runoff
  use RtmVar           , only : rtmlon, rtmlat, ice_runoff, iulog, &
                                nsrStartup, nsrContinue, nsrBranch, & 
                                inst_index, inst_suffix, inst_name, RtmVarSet, &
                                rtm_active, flood_active
  use RtmSpmd          , only : masterproc, mpicom_rof, iam, RtmSpmdInit
  use RtmMod           , only : Rtmini, Rtmrun
  use RtmTimeManager   , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep 
  use perf_mod         , only : t_startf, t_stopf, t_barrierf
  use rtm_cpl_indices  , only : rtm_cpl_indices_set, nt_rtm, rtm_tracers
  use rof_import_export

  use mct_mod
  use ESMF
!
! PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! PUBLIC MEMBER FUNCTIONS:

  public :: rof_init_mct               ! rof initialization
  public :: rof_run_mct                ! rof run phase
  public :: rof_final_mct              ! rof finalization/cleanup
!
! PRIVATE MEMBER FUNCTIONS:

  private :: rof_SetgsMap_mct         ! Set the river runoff model MCT GS map
  private :: rof_domain_mct           ! Set the river runoff model domain information
!
! PRIVATE DATA MEMBERS:
  real(r8), pointer :: totrunin(:,:)   ! cell tracer lnd forcing on rtm grid (mm/s)

! REVISION HISTORY:
! Author: Mariana Vertenstein
!===============================================================
contains
!===============================================================

!========================================================================

  subroutine rof_init_mct( EClock, cdata_r, x2r_r, r2x_r, NLFilename)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Initialize runoff model and obtain relevant atmospheric model arrays
    ! back from (i.e. albedos, surface temperature and snow cover over land).
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock     ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_r    ! Input runoff-model driver data
    type(mct_aVect) ,           intent(inout) :: x2r_r      ! River import state
    type(mct_aVect),            intent(inout) :: r2x_r      ! River export state
    character(len=*), optional, intent(in)    :: NLFilename ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    integer :: ROFID	                             ! rof identifyer
    integer :: mpicom_rof                            ! mpi communicator
    type(mct_gsMap),         pointer :: gsMap_rof    ! runoff model MCT GS map
    type(mct_gGrid),         pointer :: dom_r        ! runoff model domain
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer :: lsize                                 ! size of attribute vector
    integer :: g,i,j,n                               ! indices
    logical :: exists                                ! true if file exists
    integer :: nsrest                                ! restart type
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: begr, endr
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=32), parameter :: sub = 'rof_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    !---------------------------------------------------------------------------

    ! Obtain cdata_r (initalized in ccsm_comp_mod.F90 in the call to 
    ! seq_cdata_init for cdata_rr)
    call seq_cdata_setptrs(cdata_r, ID=ROFID, mpicom=mpicom_rof, &
         gsMap=gsMap_rof, dom=dom_r, infodata=infodata)

    ! Determine attriute vector indices
    call rtm_cpl_indices_set()

    ! Initialize rtm MPI communicator 
    call RtmSpmdInit(mpicom_rof)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_init_mct:start::',lbnum)
    endif
#endif                      

    ! Initialize io log unit
    inst_name   = seq_comm_name(ROFID)
    inst_index  = seq_comm_inst(ROFID)
    inst_suffix = seq_comm_suffix(ROFID)

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='rof_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "RTM land model initialization"
    else
       iulog = shrlogunit
    end if
    
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Initialize rtm
    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )

    call seq_infodata_GetData(infodata, case_name=caseid,                  &
                              case_desc=ctitle, start_type=starttype,      &
                              brnch_retain_casename=brnch_retain_casename, &
                              model_version=version,                       &
                              hostname=hostname, username=username)

    call timemgr_setup(calendar_in=calendar,                           &
                       start_ymd_in=start_ymd, start_tod_in=start_tod, &
                       ref_ymd_in=ref_ymd, ref_tod_in=ref_tod,         &
                       stop_ymd_in=stop_ymd, stop_tod_in=stop_tod)  

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( sub//' ERROR: unknown starttype' )
    end if

    call RtmVarSet(caseid_in=caseid, ctitle_in=ctitle,             &
                   brnch_retain_casename_in=brnch_retain_casename, &
                   nsrest_in=nsrest, version_in=version,           &
                   hostname_in=hostname, username_in=username)

    ! Read namelist, grid and surface data
    call Rtmini()

    if (rtm_active) then
       ! Initialize memory for input state
       begr = runoff%begr
       endr = runoff%endr
       allocate (totrunin(begr:endr,nt_rtm))
       
       ! Initialize rof gsMap for ocean rof and land rof
       call rof_SetgsMap_mct( mpicom_rof, ROFID, gsMap_rof)
       
       ! Initialize rof domain
       lsize = mct_gsMap_lsize(gsMap_rof, mpicom_rof)
       call rof_domain_mct( lsize, gsMap_rof, dom_r )
       
       ! Initialize lnd -> rtm attribute vector		
       call mct_aVect_init(x2r_r, rList=seq_flds_x2r_fields, lsize=lsize)
       call mct_aVect_zero(x2r_r)
       
       ! Initialize rtm -> ocn attribute vector		
       call mct_aVect_init(r2x_r, rList=seq_flds_r2x_fields, lsize=lsize)
       call mct_aVect_zero(r2x_r) 
       
       ! Create mct river runoff export state
       call rof_export(r2x_r%rattr)
    end if

    ! Fill in infodata
    call seq_infodata_PutData( infodata, rof_present=rtm_active, &
         rofice_present=.false., rof_prognostic=rtm_active, &
         rof_nx = rtmlon, rof_ny = rtmlat)
    call seq_infodata_PutData( infodata, flood_present=flood_active)

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine rof_init_mct

!---------------------------------------------------------------------------

  subroutine rof_run_mct( EClock, cdata_r, x2r_r, r2x_r)

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Run runoff model

    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_r   ! Input driver data for runoff model
    type(mct_aVect)  , intent(inout) :: x2r_r     ! Import state from runoff model
    type(mct_aVect)  , intent(inout) :: r2x_r     ! Export state from runoff model

    ! LOCAL VARIABLES:
    integer :: ymd_sync, ymd              ! current date (YYYYMMDD)
    integer :: yr_sync, yr                ! current year
    integer :: mon_sync, mon              ! current month
    integer :: day_sync, day              ! current day
    integer :: tod_sync, tod              ! current time of day (sec)
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend                      ! .true. ==> signaling last time-step
    integer :: shrlogunit,shrloglev       ! old values for share log unit and log level
    integer :: lsize                      ! local size
    integer :: lbnum                      ! input to memory diagnostic
    integer :: g,i                        ! indices
    type(mct_gGrid),        pointer :: dom_r    ! runoff model domain
    type(seq_infodata_type),pointer :: infodata ! CESM information from the driver
    real(r8),               pointer :: data(:)  ! temporary
    character(len=32)               :: rdate    ! date char string for restart file names
    character(len=32), parameter    :: sub = "rof_run_mct"
    !-------------------------------------------------------

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_run_mct:start::',lbnum)
    endif
#endif

    if (.not.rtm_active) return

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)

    ! Map MCT to land data type (output is totrunin - module variable) 
    call t_startf ('lc_rof_import')
    call rof_import(x2r_r%rattr, totrunin=totrunin)
    call t_stopf ('lc_rof_import')

    ! First advance rtm time step
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr = seq_timemgr_RestartAlarmIsOn( EClock )
    call advance_timestep()

    ! Run rtm (input is totrunin, output is runoff%runoff)
    call Rtmrun(totrunin, rstwr, nlend, rdate)

    ! Map roff data to coupler datatype (input is runoff%runoff, output is r2x_r%rattr)
    call t_startf ('lc_rof_export')
    call rof_export(r2x_r%rattr)
    call t_stopf ('lc_rof_export')

    ! Check that internal clock is in sync with master clock
    call get_curr_date( yr, mon, day, tod )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' rtm ymd=',ymd     ,'  rtm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: RTM clock is not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_run_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine rof_run_mct

!===============================================================================

  subroutine rof_final_mct( EClock, cdata_r, x2r_r, r2x_r)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Finalize rof surface model
    !
    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_r   ! Input driver data for runoff model
    type(mct_aVect)  , intent(inout) :: x2r_r     ! Import state from runoff model
    type(mct_aVect)  , intent(inout) :: r2x_r     ! Export state from runoff model
    !-----------------------------------------------------

   ! fill this in
  end subroutine rof_final_mct

!===============================================================================

  subroutine rof_SetgsMap_mct( mpicom_r, ROFID, gsMap_rof)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Set the MCT GS map for the runoff model
    !
    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: mpicom_r      ! MPI communicator for rof model
    integer        , intent(in)    :: ROFID         ! Land model identifier
    type(mct_gsMap), intent(inout) :: gsMap_rof     ! MCT gsmap for runoff -> land data
    !
    ! LOCAL VARIABLES
    integer,allocatable :: gindex(:)         ! indexing for runoff grid cells
    integer :: n, ni                         ! indices
    integer :: lsize,gsize                   ! size of runoff data and number of grid cells
    integer :: begr, endr                    ! beg, end runoff indices
    integer :: ier                           ! error code
    character(len=32), parameter :: sub = 'rof_SetgsMap_mct'
    !-----------------------------------------------------

    begr  = runoff%begr
    endr  = runoff%endr
    lsize = runoff%lnumr
    gsize = rtmlon*rtmlat

    ! Check 
    ni = 0
    do n = begr,endr
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Determine gsmap_rof
    allocate(gindex(lsize),stat=ier)
    ni = 0
    do n = begr,endr
       ni = ni + 1
       gindex(ni) = runoff%gindex(n)
    end do
    call mct_gsMap_init( gsMap_rof, gindex, mpicom_r, ROFID, lsize, gsize )
    deallocate(gindex)

  end subroutine rof_SetgsMap_mct

!===============================================================================

  subroutine rof_domain_mct( lsize, gsMap_r, dom_r )

    !-----------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Send the runoff model domain information to the coupler
    !
    ! !ARGUMENTS:
    implicit none
    integer        , intent(in)    :: lsize       ! Size of runoff domain information
    type(mct_gsMap), intent(inout) :: gsMap_r     ! Output MCT GS map for runoff model
    type(mct_ggrid), intent(out)   :: dom_r       ! Domain information from the runoff model
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
    call mct_gGrid_init( GGrid=dom_r, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Allocate memory
    allocate(data(lsize))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_r, iam, idata)
    call mct_gGrid_importIAttr(dom_r,'GlobGridNum',idata,lsize)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_r,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_r,"mask" ,data,lsize) 

    ! Determine bounds numbering consistency
    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       if (ni > runoff%lnumr) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       end if
    end do
    if (ni /= runoff%lnumr) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Fill in correct values for domain components
    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       data(ni) = runoff%lonc(n)
    end do
    call mct_gGrid_importRattr(dom_r,"lon",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       data(ni) = runoff%latc(n)
    end do
    call mct_gGrid_importRattr(dom_r,"lat",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       data(ni) = runoff%area(n)*1.0e-6_r8/(re*re)
    end do
    call mct_gGrid_importRattr(dom_r,"area",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       data(ni) = 1.0_r8
    end do
    call mct_gGrid_importRattr(dom_r,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_r,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine rof_domain_mct

!====================================================================================
 

end module rof_comp_mct
