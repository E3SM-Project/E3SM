module rof_comp_mct
  
!========================================================================
! DESCRIPTION:
! Interface of the active runoff component of CESM 
! with the main CESM driver. This is a thin interface taking CESM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by MOSART

  use seq_flds_mod
  use shr_kind_mod     , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use shr_taskmap_mod  , only : shr_taskmap_write
  use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
  use seq_comm_mct     , only : info_taskmap_comp
  use seq_timemgr_mod  , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
  use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                seq_infodata_start_type_brnch
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use RunoffMod        , only : rtmCTL, TRunoff, THeat, TUnit, Tctl
  use RtmVar           , only : rtmlon, rtmlat, ice_runoff, iulog, &
                                nsrStartup, nsrContinue, nsrBranch, & 
                                inst_index, inst_suffix, inst_name, RtmVarSet, &
                                wrmflag, heatflag, data_bgc_fluxes_to_ocean_flag, &
                                inundflag, use_lnd_rof_two_way, use_ocn_rof_two_way, &
                                sediflag
  use RtmSpmd          , only : masterproc, mpicom_rof, npes, iam, RtmSpmdInit, ROFID
  use RtmMod           , only : Rtmini, Rtmrun
  use RtmTimeManager   , only : timemgr_setup, get_curr_date, get_step_size
  use perf_mod         , only : t_startf, t_stopf, t_barrierf

  use WRM_type_mod     , only : StorWater

  use rof_cpl_indices  , only : rof_cpl_indices_set, nt_rtm, rtm_tracers, &
                                index_x2r_Flrl_rofsur, index_x2r_Flrl_rofi, &
                                index_x2r_Flrl_rofgwl, index_x2r_Flrl_rofsub, &
                                index_x2r_Flrl_rofdto, index_x2r_Flrl_demand, &
                                index_x2r_Flrl_Tqsur, index_x2r_Flrl_Tqsub, &
                                index_x2r_Sa_tbot, index_x2r_Sa_pbot, &
                                index_x2r_Sa_u   , index_x2r_Sa_v   , &
                                index_x2r_Sa_shum, &
                                index_x2r_So_ssh,  &
                                index_x2r_Faxa_lwdn , &
                                index_x2r_Faxa_swvdr, index_x2r_Faxa_swvdf, &
                                index_x2r_Faxa_swndr, index_x2r_Faxa_swndf, &
                                index_x2r_Flrl_rofdto, index_x2r_Flrl_rofmud, &
                                index_r2x_Forr_rofl, index_r2x_Forr_rofi, &
                                index_r2x_Flrr_flood, &
                                index_r2x_Forr_rofDIN, index_r2x_Forr_rofDIP, &
                                index_r2x_Forr_rofDON, index_r2x_Forr_rofDOP, &
                                index_r2x_Forr_rofDOC, index_r2x_Forr_rofPP , &
                                index_r2x_Forr_rofDSi, index_r2x_Forr_rofPOC, &
                                index_r2x_Forr_rofPN , index_r2x_Forr_rofDIC, &
                                index_r2x_Forr_rofFe , &
                                index_r2x_Flrr_volr, index_r2x_Flrr_volrmch, &
                                index_x2r_coszen_str, &
                                index_r2x_Flrr_supply, index_r2x_Flrr_deficit, &
                                index_r2x_Sr_h2orof, index_r2x_Sr_frac_h2orof, &
                                index_x2r_Flrl_inundinf

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
! PUBLIC DATA MEMBERS:
! None
!
! PRIVATE MEMBER FUNCTIONS:
  private :: rof_SetgsMap_mct         ! Set the river runoff model MCT GS map
  private :: rof_domain_mct           ! Set the river runoff model domain information
  private :: rof_export_mct           ! Export the river runoff model data to the CESM coupler
!
! PRIVATE DATA MEMBERS:

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
    logical :: rof_prognostic                        ! flag
    logical :: flood_present                         ! flag
    logical :: rofocn_prognostic                     ! ocn rof two way coupling flag
    integer :: mpicom_loc                            ! mpi communicator
    type(mct_gsMap),         pointer :: gsMap_rof    ! runoff model MCT GS map
    type(mct_gGrid),         pointer :: dom_r        ! runoff model domain
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer :: lsize                                 ! size of attribute vector
    integer :: g,i,j,n                               ! indices
    logical :: exists                                ! true if file exists
    logical :: verbose_taskmap_output                ! true then use verbose task-to-node mapping format
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
    character(len=8)           :: c_inst_index       ! instance number
    character(len=8)           :: c_npes             ! number of pes
    character(len=32), parameter :: sub = 'rof_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    !---------------------------------------------------------------------------

    ! Obtain cdata_r (initalized in ccsm_comp_mod.F90 in the call to 
    ! seq_cdata_init for cdata_rr)
    call seq_cdata_setptrs(cdata_r, ID=ROFID, mpicom=mpicom_loc, &
         gsMap=gsMap_rof, dom=dom_r, infodata=infodata)

    ! Determine attriute vector indices
    call rof_cpl_indices_set()

    ! Initialize mosart MPI communicator 
    call RtmSpmdInit(mpicom_loc)

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
       write(iulog,format) "MOSART model initialization"
    else
       iulog = shrlogunit
    end if
    
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,*) ' mosart npes = ',npes
       write(iulog,*) ' mosart iam  = ',iam
       write(iulog,*) ' inst_name = ',trim(inst_name)
    endif

    ! Identify SMP nodes and process/SMP mapping for this instance.
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
             ' pes participating in computation of MOSART instance #', &
             trim(adjustl(c_inst_index))
          call shr_sys_flush(iulog)
       endif

       call t_startf("shr_taskmap_write")
       call shr_taskmap_write(iulog, mpicom_rof,                    &
                              'ROF #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output        )
       call t_stopf("shr_taskmap_write")

    endif

    ! Initialize mosart
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

    use_lnd_rof_two_way = lnd_rof_two_way
    use_ocn_rof_two_way = ocn_rof_two_way

    ! Read namelist, grid and surface data
    call Rtmini(rtm_active=rof_prognostic,flood_active=flood_present)

    if (rof_prognostic) then
       ! Initialize memory for input state
       begr = rtmCTL%begr
       endr = rtmCTL%endr
       
       ! Initialize rof gsMap for ocean rof and land rof
       call rof_SetgsMap_mct( mpicom_rof, ROFID, gsMap_rof)
       
       ! Initialize rof domain
       lsize = mct_gsMap_lsize(gsMap_rof, mpicom_rof)
       call rof_domain_mct( lsize, gsMap_rof, dom_r )
       
       ! Initialize lnd -> mosart attribute vector
       call mct_aVect_init(x2r_r, rList=seq_flds_x2r_fields, lsize=lsize)
       call mct_aVect_zero(x2r_r)

       ! Initialize mosart -> ocn attribute vector        
       call mct_aVect_init(r2x_r, rList=seq_flds_r2x_fields, lsize=lsize)
       call mct_aVect_zero(r2x_r) 
       
       ! Create mct river runoff export state
       call rof_export_mct( r2x_r )
    else
       call seq_infodata_PutData(infodata, rofice_present=.false.)
    end if

    ! Fill in infodata
    call seq_infodata_PutData( infodata, rof_present=rof_prognostic, rof_nx = rtmlon, rof_ny = rtmlat, &
         rof_prognostic=rof_prognostic, rofocn_prognostic=use_ocn_rof_two_way)
    call seq_infodata_PutData( infodata, flood_present=flood_present)

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

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)

    ! Map MCT to land data type (output is totrunin, subrunin)
    call t_startf ('lc_rof_import')
    call rof_import_mct( x2r_r)
    call t_stopf ('lc_rof_import')

    ! Run mosart (input is *runin, output is rtmCTL%runoff)
    ! First advance mosart time step
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr = seq_timemgr_RestartAlarmIsOn( EClock )
    call Rtmrun(rstwr,nlend,rdate)

    ! Map roff data to MCT datatype (input is rtmCTL%runoff, output is r2x_r)
    call t_startf ('lc_rof_export')
    call rof_export_mct( r2x_r )
    call t_stopf ('lc_rof_export')

    ! Check that internal clock is in sync with master clock
    call get_curr_date( yr, mon, day, tod )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' mosart ymd=',ymd     ,'  mosart tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: MOSART clock is not in sync with Master Sync clock" )
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

    begr  = rtmCTL%begr
    endr  = rtmCTL%endr
    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat

    ! Check 
    ni = 0
    do n = begr,endr
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,rtmCTL%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,rtmCTL%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Determine gsmap_rof
    allocate(gindex(lsize),stat=ier)
    ni = 0
    do n = begr,endr
       ni = ni + 1
       gindex(ni) = rtmCTL%gindex(n)
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
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       if (ni > rtmCTL%lnumr) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,rtmCTL%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       end if
    end do
    if (ni /= rtmCTL%lnumr) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,rtmCTL%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Fill in correct values for domain components
    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       data(ni) = rtmCTL%lonc(n)
    end do
    call mct_gGrid_importRattr(dom_r,"lon",data,lsize) 

    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       data(ni) = rtmCTL%latc(n)
    end do
    call mct_gGrid_importRattr(dom_r,"lat",data,lsize) 

    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       data(ni) = rtmCTL%area(n)*1.0e-6_r8/(re*re)
    end do
    call mct_gGrid_importRattr(dom_r,"area",data,lsize) 

    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       data(ni) = 1.0_r8
    end do
    call mct_gGrid_importRattr(dom_r,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_r,"frac",data,lsize) 

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
    integer :: n2, n, nt, begr, endr, nliq, nfrz, nmud, nsan
    real(R8) :: tmp1, tmp2
    real(R8) :: shum
    character(len=32), parameter :: sub = 'rof_import_mct'
    !---------------------------------------------------------------------------
    
    ! Note that ***runin are fluxes

    nliq = 0
    nfrz = 0
    nmud = 0
    nsan = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
       if (trim(rtm_tracers(nt)) == 'MUD') then
          nmud = nt
       endif
       if (trim(rtm_tracers(nt)) == 'SAN') then
          nsan = nt
       endif
    enddo
    if (nliq == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers LIQ',nliq,rtm_tracers
       call shr_sys_abort()
    endif
    if (nfrz == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers ICE',nfrz,rtm_tracers
       call shr_sys_abort()
    endif
    if (nmud == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers MUD',nmud,rtm_tracers
       call shr_sys_abort()
    endif
    if (nsan == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers SAN',nsan,rtm_tracers
       call shr_sys_abort()
    endif

    begr = rtmCTL%begr
    endr = rtmCTL%endr
    do n = begr,endr
       n2 = n - begr + 1

       rtmCTL%qsur(n,nliq) = x2r_r%rAttr(index_x2r_Flrl_rofsur,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qsub(n,nliq) = x2r_r%rAttr(index_x2r_Flrl_rofsub,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qgwl(n,nliq) = x2r_r%rAttr(index_x2r_Flrl_rofgwl,n2) * (rtmCTL%area(n)*0.001_r8)
       if (index_x2r_Flrl_rofdto > 0) then
          rtmCTL%qdto(n,nliq) = x2r_r%rAttr(index_x2r_Flrl_rofdto,n2) * (rtmCTL%area(n)*0.001_r8)
       else
          rtmCTL%qdto(n,nliq) = 0.0_r8
       endif
       if (wrmflag) then
          rtmCTL%qdem(n,nliq) = x2r_r%rAttr(index_x2r_Flrl_demand,n2) / TUnit%domainfrac(n) * (rtmCTL%area(n)*0.001_r8)
       else
          rtmCTL%qdem(n,nliq) = 0.0_r8
       endif
       rtmCTL%qsur(n,nfrz) = x2r_r%rAttr(index_x2r_Flrl_rofi,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qsub(n,nfrz) = 0.0_r8
       rtmCTL%qgwl(n,nfrz) = 0.0_r8
       rtmCTL%qdto(n,nfrz) = 0.0_r8
       rtmCTL%qdem(n,nfrz) = 0.0_r8

       if (index_x2r_So_ssh>0) then
          rtmCTL%ssh(n)       = x2r_r%rAttr(index_x2r_So_ssh,n2)
       end if

       if(heatflag) then
          rtmCTL%Tqsur(n) = x2r_r%rAttr(index_x2r_Flrl_Tqsur,n2)
          rtmCTL%Tqsub(n) = x2r_r%rAttr(index_x2r_Flrl_Tqsub,n2)
          THeat%Tqsur(n) = rtmCTL%Tqsur(n)
          THeat%Tqsub(n) = rtmCTL%Tqsub(n)
       
          THeat%forc_t(n) = x2r_r%rAttr(index_x2r_Sa_tbot,n2)
          THeat%forc_pbot(n) = x2r_r%rAttr(index_x2r_Sa_pbot,n2)
          tmp1 = x2r_r%rAttr(index_x2r_Sa_u   ,n2)
          tmp2 = x2r_r%rAttr(index_x2r_Sa_v   ,n2)
          THeat%forc_wind(n) = sqrt(tmp1*tmp1 + tmp2*tmp2)
          THeat%forc_lwrad(n)= x2r_r%rAttr(index_x2r_Faxa_lwdn ,n2)
          THeat%forc_solar(n)= x2r_r%rAttr(index_x2r_Faxa_swvdr,n2) + x2r_r%rAttr(index_x2r_Faxa_swvdf,n2) + &
                               x2r_r%rAttr(index_x2r_Faxa_swndr,n2) + x2r_r%rAttr(index_x2r_Faxa_swndf,n2)
          shum = x2r_r%rAttr(index_x2r_Sa_shum,n2)
          THeat%forc_vp(n)   = shum * THeat%forc_pbot(n)  / (0.622_r8 + 0.378_r8 * shum)
          THeat%coszen(n)    = x2r_r%rAttr(index_x2r_coszen_str,n2)
       end if

       rtmCTL%qsur(n,nmud) = 0.0_r8
       rtmCTL%qsur(n,nsan) = 0.0_r8

       if (index_x2r_Flrl_inundinf > 0) then
          rtmCTL%inundinf(n) = x2r_r%rAttr(index_x2r_Flrl_inundinf,n2) * (rtmCTL%area(n)*0.001_r8)
       endif

    enddo

    if(sediflag) then
        do n = begr,endr
           n2 = n - begr + 1
           rtmCTL%qsur(n,nmud) = x2r_r%rAttr(index_x2r_Flrl_rofmud,n2) * (rtmCTL%area(n)) ! kg/m2/s --> kg/s for sediment
           rtmCTL%qsur(n,nsan) = 0.0_r8
        enddo
    end if

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
    
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    r2x_r%rattr(:,:) = 0._r8

    if (first_time) then
       if (masterproc) then
       if ( ice_runoff )then
          write(iulog,*)'Snow capping will flow out in frozen river runoff'
       else
          write(iulog,*)'Snow capping will flow out in liquid river runoff'
       endif
       endif
       first_time = .false.
    end if

    ni = 0
    if ( ice_runoff )then
       ! separate liquid and ice runoff
       do n = rtmCTL%begr,rtmCTL%endr
          ni = ni + 1
          r2x_r%rAttr(index_r2x_Forr_rofl,ni) =  rtmCTL%direct(n,nliq) / (rtmCTL%area(n)*0.001_r8)
          r2x_r%rAttr(index_r2x_Forr_rofi,ni) =  rtmCTL%direct(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             ! liquid and ice runoff are treated separately - this is what goes to the ocean
             r2x_r%rAttr(index_r2x_Forr_rofl,ni) = r2x_r%rAttr(index_r2x_Forr_rofl,ni) + &
                rtmCTL%runoff(n,nliq) / (rtmCTL%area(n)*0.001_r8)
             r2x_r%rAttr(index_r2x_Forr_rofi,ni) = r2x_r%rAttr(index_r2x_Forr_rofi,ni) + &
                rtmCTL%runoff(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
             if (ni > rtmCTL%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
! note runoff has already been divided by area so do not need to do it again for nutrient flux
             if (data_bgc_fluxes_to_ocean_flag) then
               tmp1 = r2x_r%rAttr(index_r2x_Forr_rofl,ni)
               r2x_r%rAttr(index_r2x_Forr_rofDIN,ni) =  tmp1*rtmCTL%concDIN(n)
               r2x_r%rAttr(index_r2x_Forr_rofDIP,ni) =  tmp1*rtmCTL%concDIP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDON,ni) =  tmp1*rtmCTL%concDON(n)
               r2x_r%rAttr(index_r2x_Forr_rofDOP,ni) =  tmp1*rtmCTL%concDOP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDOC,ni) =  tmp1*rtmCTL%concDOC(n)
               r2x_r%rAttr(index_r2x_Forr_rofPP ,ni) =  tmp1*rtmCTL%concPP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDSi,ni) =  tmp1*rtmCTL%concDSi(n)
               r2x_r%rAttr(index_r2x_Forr_rofPOC,ni) =  tmp1*rtmCTL%concPOC(n)
               r2x_r%rAttr(index_r2x_Forr_rofPN ,ni) =  tmp1*rtmCTL%concPN(n)
               r2x_r%rAttr(index_r2x_Forr_rofDIC,ni) =  tmp1*rtmCTL%concDIC(n)
               r2x_r%rAttr(index_r2x_Forr_rofFe,ni)  =  tmp1*rtmCTL%concFe(n)
             end if

          endif
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = rtmCTL%begr,rtmCTL%endr
          ni = ni + 1
          r2x_r%rAttr(index_r2x_Forr_rofl,ni) =  &
             (rtmCTL%direct(n,nfrz)+rtmCTL%direct(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             r2x_r%rAttr(index_r2x_Forr_rofl,ni) = r2x_r%rAttr(index_r2x_Forr_rofl,ni) + &
                (rtmCTL%runoff(n,nfrz)+rtmCTL%runoff(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
             if (ni > rtmCTL%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
! note runoff has already been divided by area so do not need to do it again for nutrient flux
             if (data_bgc_fluxes_to_ocean_flag) then
               tmp1 = r2x_r%rAttr(index_r2x_Forr_rofl,ni)
               r2x_r%rAttr(index_r2x_Forr_rofDIN,ni) =  tmp1*rtmCTL%concDIN(n)
               r2x_r%rAttr(index_r2x_Forr_rofDIP,ni) =  tmp1*rtmCTL%concDIP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDON,ni) =  tmp1*rtmCTL%concDON(n)
               r2x_r%rAttr(index_r2x_Forr_rofDOP,ni) =  tmp1*rtmCTL%concDOP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDOC,ni) =  tmp1*rtmCTL%concDOC(n)
               r2x_r%rAttr(index_r2x_Forr_rofPP ,ni) =  tmp1*rtmCTL%concPP(n)
               r2x_r%rAttr(index_r2x_Forr_rofDSi,ni) =  tmp1*rtmCTL%concDSi(n)
               r2x_r%rAttr(index_r2x_Forr_rofPOC,ni) =  tmp1*rtmCTL%concPOC(n)
               r2x_r%rAttr(index_r2x_Forr_rofPN ,ni) =  tmp1*rtmCTL%concPN(n)
               r2x_r%rAttr(index_r2x_Forr_rofDIC,ni) =  tmp1*rtmCTL%concDIC(n)
               r2x_r%rAttr(index_r2x_Forr_rofFe,ni)  =  tmp1*rtmCTL%concFe(n)
             end if
          endif
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ni = 0
    do n = rtmCTL%begr, rtmCTL%endr
       ni = ni + 1
       r2x_r%rattr(index_r2x_Flrr_flood,ni)   = -rtmCTL%flood(n) / (rtmCTL%area(n)*0.001_r8)
       r2x_r%rattr(index_r2x_Flrr_volr,ni)    = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / rtmCTL%area(n)
       r2x_r%rattr(index_r2x_Flrr_volrmch,ni) = Trunoff%wr(n,nliq) / rtmCTL%area(n)
       r2x_r%rattr(index_r2x_Flrr_supply,ni)  = 0._r8
       r2x_r%rattr(index_r2x_Flrr_deficit,ni)  = 0._r8
       if (wrmflag) then
          r2x_r%rattr(index_r2x_Flrr_supply,ni)  = StorWater%Supply(n) / (rtmCTL%area(n)*0.001_r8)   !converted to mm/s
          r2x_r%rattr(index_r2x_Flrr_deficit,ni)  = (abs(rtmCTL%qdem(n,nliq)) - abs(StorWater%Supply(n))) / (rtmCTL%area(n)*0.001_r8)   !send deficit back to ELM
       endif
    end do

    if ( index_r2x_Sr_h2orof > 0 ) then
      ni = 0
      do n = rtmCTL%begr, rtmCTL%endr
        ni = ni + 1
        r2x_r%rattr(index_r2x_Sr_h2orof,ni)      = rtmCTL%inundwf(n) / (rtmCTL%area(n)*0.001_r8) ! m^3 to mm
        r2x_r%rattr(index_r2x_Sr_frac_h2orof,ni) = rtmCTL%inundff(n)
      enddo
    endif

  end subroutine rof_export_mct

end module rof_comp_mct
