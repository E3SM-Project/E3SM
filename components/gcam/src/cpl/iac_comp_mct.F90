module iac_comp_mct
  
!========================================================================
! DESCRIPTION:
! Interface of integrated assessment component with the main E3SM driver. 
! This is a thin interface taking E3SM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by GCAM
! Shamelessly copied and stolen from other components (mostly
! RTM/ROF) in the ACME development tree, so if something looks wrong
! it's probably because I didn't understand correctly what the
! original code did.

! Note: Iac single letter descriptor is 'z', because all the other logical ones were
! already taken.

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
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst,seq_comm_name

  use iac_mod          , only : iac
  use gcam_var           , only : gcam_lon, gcam_lat,  iulog, &
                                nsrStartup, nsrContinue, nsrBranch, & 
                                inst_index, inst_suffix, inst_name, &
                                gcam_active, gcam_var_set

! Stub in other moduals for running iac
! use IacMod
! use IacVar

!
! PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mct               ! iac initialization
  public :: iac_run_mct                ! iac run phase
  public :: iac_final_mct              ! iac finalization/cleanup
!
! PRIVATE MEMBER FUNCTIONS:

! PRIVATE MEMBER FUNCTIONS:

  private :: iac_SetgsMap_mct         ! Set the iac model MCT GS map
  private :: iac_domain_mct           ! Set the iac model domain information

! REVISION HISTORY:
! Author: Tim Shippert
!===============================================================
contains
!===============================================================

!========================================================================

  subroutine iac_init_mct( EClock, cdata_z, x2z_z, z2x_z, NLFilename)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Initialize iac model and hook in arrays from lnd module
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock     ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_z    ! Input iac driver data
    type(mct_aVect) ,           intent(inout) :: x2z_z      ! Iac import state
    type(mct_aVect),            intent(inout) :: z2x_z      ! Iac export state
    character(len=*), optional, intent(in)    :: NLFilename ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    ! 
    integer :: IACID	                             ! iac identifyer
    integer :: mpicom_iac                            ! mpi communicator
    type(mct_gsMap),         pointer :: gsMap_iac    ! iac model MCT GS map
    type(mct_gGrid),         pointer :: dom_z        ! iac model domain
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

    ! Possibly not useful, as we expect gcam to run with just one proc
    integer :: begr, endr                            ! Region indeces

    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=32), parameter :: sub = 'iac_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    !---------------------------------------------------------------------------
   
    ! Pull out our cdata and communication stuff, from wherever it is set
    ! these days
    call seq_cdata_setptrs(cdata_z, ID=IACID, mpicom=mpicom_iac, &
         gsMap=gsMap_iac, dom=dom_z, infodata=infodata)

    ! Determine attriute vector indices
    call gcam_cpl_indices_set()

    ! Initialize gcam MPI communicator 
    call gcam_mpi_init(mpicom_iac)

    ! I see this kind of thing everywhere, so why not
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','iac_init_mct:start::',lbnum)
    endif
#endif                      

    ! Initialize io logs
    inst_name   = seq_comm_name(IACID)
    inst_index  = seq_comm_inst(IACID)
    inst_suffix = seq_comm_suffix(IACID)

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       ! I'm copying from other comps, but we may want to change teh
       ! namelist file here
       inquire(file='iac_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('iac_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "IAC initialization"
    else
       iulog = shrlogunit
    end if
    
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Generic timer and infodata stuff
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

    ! Startup type - we'll be generic here too for now
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( sub//' ERROR: unknown starttype' )
    end if

    ! Modeled on RTM, we'll modify later - passing stuff back to gcam
    ! from our namelists and configs
    call gcam_var_set(caseid_in=caseid, ctitle_in=ctitle,             &
         brnch_retain_casename_in=brnch_retain_casename, &
         nsrest_in=nsrest, version_in=version,           &
         hostname_in=hostname, username_in=username)

    ! Do whatever init gcam needs
    call gcam_init()

    ! Here it gets tricky, as I copy from rof and other components.
    ! I'm going to go ahead and use the whole begr,endr regional
    ! indeces for now - I'm not sure if there is any meaningful
    ! overhead if we do this with only one proc, but if we ever go
    ! multiproc it will be better to do this now, and it lets me cut
    ! and paste for now.

    if (gcam_active) then
       ! Initialize memory for input state
       begr = iac%begr
       endr = iac%endr
       allocate (totrunin(begr:endr,nt_gcam))
       
       ! Initialize iac gsMap 
       call iac_SetgsMap_mct( mpicom_iac, IACID, gsMap_iac)
       
       ! Initialize iac domain - I feel like this should be mostly null
       lsize = mct_gsMap_lsize(gsMap_iac, mpicom_iac)
       call iac_domain_mct( lsize, gsMap_iac, dom_z )
       
       ! Initialize input attribute vectors
       call mct_aVect_init(x2z_z, rList=seq_flds_x2z_fields, lsize=lsize)
       call mct_aVect_zero(x2z_z)
       
       ! Initialize output attribute vectors
       call mct_aVect_init(z2x_z, rList=seq_flds_z2x_fields, lsize=lsize)
       call mct_aVect_zero(z2x_z) 
       
       ! Create mct iac export state - try to review what this is.
       call iac_export(r2x_r%rattr)
    end if

    ! Fill in infodata - of course, have to review all this
    call seq_infodata_PutData( infodata, iac_present=gcam_active, &
         iacice_present=.false., iac_prognostic=gcam_active, &
         iac_nx = gcamlon, iac_ny = gcamlat)
    call seq_infodata_PutData( infodata, flood_present=flood_active)

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','iac_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif


  end subroutine iac_init_mct

!---------------------------------------------------------------------------

  subroutine iac_run_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Run IAC model

    ! !USES:
    use shr_kind_mod    ,  only : r8 => shr_kind_r8
    use gcam_instMod    , only : iac2atm_vars, atm2iac_vars, iac2lnd_vars, lnd2iac_vars
    use gcam_driver     ,  only : gcam_drv
    use gcam_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use gcam_time_manager,  only : advance_timestep, set_nextsw_cday,update_rad_dtime
    use decompMod       ,  only : get_proc_bounds
    use abortutils      ,  only : endrun
    use gcam_varctl      ,  only : iulog
    use gcam_varorb      ,  only : eccen, obliqr, lambm0, mvelpp
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

    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for iac model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from iac model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from iac model

    ! LOCAL VARIABLES:
    integer      :: ymd_sync             ! Sync date (YYYYMMDD)
    integer      :: yr_sync              ! Sync current year
    integer      :: mon_sync             ! Sync current month
    integer      :: day_sync             ! Sync current day
    integer      :: tod_sync             ! Sync current time of day (sec)
    integer      :: ymd                  ! CLM current date (YYYYMMDD)
    integer      :: yr                   ! CLM current year
    integer      :: mon                  ! CLM current month
    integer      :: day                  ! CLM current day
    integer      :: tod                  ! CLM current time of day (sec)
    integer      :: dtime                ! time step increment (sec)
    integer      :: nstep                ! time step index
    logical      :: rstwr_sync           ! .true. ==> write restart file before returning
    logical      :: rstwr                ! .true. ==> write restart file before returning
    logical      :: nlend_sync           ! Flag signaling last time-step
    logical      :: nlend                ! .true. ==> last time-step
    logical      :: dosend               ! true => send data back to driver
    logical      :: doalb                ! .true. ==> do albedo calculation on this time step
    real(r8)     :: nextsw_cday          ! calday from clock of next radiation computation
    real(r8)     :: caldayp1             ! gcam calday plus dtime offset
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
    logical      :: lnd_present
    type(seq_infodata_type),pointer :: infodata             ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_z                ! iac model domain data
    type(bounds_type)               :: bounds               ! bounds
    character(len=32)               :: rdate                ! date char string for restart file names
    character(len=32), parameter    :: sub = "iac_run_mct"

    ! I feel this is probably important
    if (.not.iac_active) return

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation
    ! This appears to be just for tagging purposes
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)

    ! We might do some mapping here, taking MCT to our GCAM regions.
    ! But I'm not sure yet - see line 284 in rof_comp_mct.F90, for
    ! example.   

    ! Or, this is where we extract our MCT Avects into something gcam
    ! will use?  Maybe generate those iaci, iaco AVects?

    ! First advance iac time step
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr = seq_timemgr_RestartAlarmIsOn( EClock )
    call advance_timestep()

    ! Run gcam.  At this point, iaci (input) and iaco (output) should
    ! be MCT AVects.  Instead, it might make sense to write external
    ! functions to go back and forth - but since in iESM it was all
    ! handled via file I/O I have to rewrite the mapping code
    ! internally anyway, so I might as well simply use the MCT
    ! method, since it is already there.
    call iac_run_mod( Eclock, cdata_z, iaci, iaco)

    ! Some logging and timing checvks
    ! Check that internal clock is in sync with master clock
    call get_curr_date( yr, mon, day, tod )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' gcam ymd=',ymd     ,'  gcam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: RTM clock is not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine iac_run_mct

!===============================================================================

  subroutine iac_final_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Finalize iac surface model
    !
    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for iac model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from iac model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from iac model
    !-----------------------------------------------------

    ! fill this in
    ! For now, no cleanup - we'll see later on
  end subroutine iac_final_mct

!================================================================================
! Internal functions to deal with gsmap and domain

  subroutine iac_SetgsMap_mct( mpicom_z, IACID, gsMap_iac)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Set the MCT GS map for the iac model
    !
    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: mpicom_z      ! MPI communicator for iac model
    integer        , intent(in)    :: IACID         ! iac model identifier
    type(mct_gsMap), intent(inout) :: gsMap_iac     ! MCT gsmap for iac -> land data
    !
    ! LOCAL VARIABLES
    integer,allocatable :: gindex(:)         ! indexing for iac grid cells
    integer :: n, ni                         ! indices
    integer :: lsize,gsize                   ! size of iac data and number of grid cells
    integer :: begr, endr                    ! beg, end iac indices
    integer :: ier                           ! error code
    character(len=32), parameter :: sub = 'iac_SetgsMap_mct'
    !-----------------------------------------------------

    begr  = iac%begr
    endr  = iac%endr
    lsize = iac%lnumr
    gsize = gcam_lon*gcam_lat

    ! Check 
    ni = 0
    do n = begr,endr
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR iac count',n,ni,iac%lnumr
          call shr_sys_abort( sub//' ERROR: iac > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR iac total count',ni,iac%lnumr
       call shr_sys_abort( sub//' ERROR: iac not equal to expected' )
    endif

    ! Determine gsmap_iac
    allocate(gindex(lsize),stat=ier)
    ni = 0
    do n = begr,endr
       ni = ni + 1
       gindex(ni) = iac%gindex(n)
    end do
    call mct_gsMap_init( gsMap_iac, gindex, mpicom_z, IACID, lsize, gsize )
    deallocate(gindex)

  end subroutine iac_SetgsMap_mct

!===============================================================================

  subroutine iac_domain_mct( lsize, gsMap_z, dom_z )

    !-----------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Send the iac model domain information to the coupler
    ! 
    ! Note: GCAM is currently only single proc, so domain information may not be
    ! necessary or useful.  But I want the infrastructure in case
    ! that changes, and so that the calling structure is the same for
    ! this componant is the same as for the others.  Also, its
    ! completely possible that I don't really understand what
    ! "domain" means in the context of E3SM.
    !
    ! !ARGUMENTS:
    implicit none
    integer        , intent(in)    :: lsize       ! Size of iac domain information
    type(mct_gsMap), intent(inout) :: gsMap_z     ! Output MCT GS map for iac model
    type(mct_ggrid), intent(out)   :: dom_z       ! Domain information from the iac model
    !
    ! LOCAL VARIABLES
    integer :: n, ni              ! index
    integer , pointer :: idata(:) ! temporary
    real(r8), pointer :: data(:)  ! temporary
    real(r8) :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)
    character(len=32), parameter :: sub = 'iac_domain_mct'
    !-----------------------------------------------------

    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! TRS - this is something I need to review
    call mct_gGrid_init( GGrid=dom_z, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Allocate memory
    allocate(data(lsize))

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_z, iam, idata)
    call mct_gGrid_importIAttr(dom_z,'GlobGridNum',idata,lsize)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_z,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_z,"mask" ,data,lsize) 

    ! Determine bounds numbering consistency
    ni = 0
    do n = iac%begr,iac%endr
       ni = ni + 1
       if (ni > iac%lnumr) then
          write(iulog,*) sub, ' : ERROR iac count',n,ni,iac%lnumr
          call shr_sys_abort( sub//' ERROR: iac > expected' )
       end if
    end do
    if (ni /= iac%lnumr) then
       write(iulog,*) sub, ' : ERROR iac total count',ni,iac%lnumr
       call shr_sys_abort( sub//' ERROR: iac not equal to expected' )
    endif

    ! Fill in correct values for domain components
    ni = 0
    do n = iac%begr,iac%endr
       ni = ni + 1
       data(ni) = iac%lonc(n)
    end do
    call mct_gGrid_importRattr(dom_z,"lon",data,lsize) 

    ni = 0
    do n = iac%begr,iac%endr
       ni = ni + 1
       data(ni) = iac%latc(n)
    end do
    call mct_gGrid_importRattr(dom_z,"lat",data,lsize) 

    ni = 0
    do n = iac%begr,iac%endr
       ni = ni + 1
       data(ni) = iac%area(n)*1.0e-6_r8/(re*re)
    end do
    call mct_gGrid_importRattr(dom_z,"area",data,lsize) 

    ni = 0
    do n = iac%begr,iac%endr
       ni = ni + 1
       data(ni) = 1.0_r8
    end do
    call mct_gGrid_importRattr(dom_z,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_z,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine iac_domain_mct

!====================================================================================


!===============================================================================

!==============================================================
! Local functions
!==============================================================

end module iac_comp_mct
