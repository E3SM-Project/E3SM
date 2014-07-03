module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Parallel Ocean Program (POP).
!
! !REVISION HISTORY:
!  SVN:$Id:
!
! !USES:
   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_IOUnitsMod
   use POP_MCT_vars_mod

   use mct_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name
   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use perf_mod
   use ocn_communicator,  only: mpi_communicator_ocn

   use kinds_mod,         only: int_kind, r8
   use POP_CplIndices
   use POP_KindsMod
   use POP_ErrorMod
   use POP_InitMod,       only: POP_Initialize1, POP_Initialize2, &
                                timer_total, cpl_ts 
   use communicate,       only: my_task, master_task
   use constants
   use blocks
   use domain,            only: distrb_clinic, POP_haloClinic
   use exit_mod
   use forcing_shf,       only: SHF_QSW
   use forcing_sfwf,      only: lsend_precip_fact, precip_fact
   use forcing_fields
   use forcing_coupled,   only: ncouple_per_day,  &
                                update_ghost_cells_coupler_fluxes, &
                                rotate_wind_stress, pop_set_coupled_forcing, &
                                pop_init_coupled,  &
                                orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
   use ice,               only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, &
                                tlast_ice
   use grid,              only: TLAT, TLON, KMT
   use global_reductions, only: global_sum_prod
   use io_tools,          only: document
   use named_field_mod,   only: named_field_register, named_field_get_index, &
                                named_field_set, named_field_get
   use prognostic
   use timers,            only: get_timer, timer_start, timer_stop
   use diagnostics,       only: check_KE
   use output,            only: output_driver
   use step_mod,          only: step
   use time_management
   use registry
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ocn_export_mct
  private :: ocn_import_mct
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct
!
! !PRIVATE MODULE VARIABLES

  logical (log_kind) ::   &
       ldiag_cpl = .false.

  integer (int_kind), private ::   &
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg      
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp       ! flag id for computing diagnostics

  real (r8),   &
      dimension(:,:,:,:), allocatable ::  &
      SBUFF_SUM           ! accumulated sum of send buffer quantities
                          ! for averaging before being sent
   real (r8) ::  &
     tlast_coupled

   integer (int_kind)  ::   &
      nsend, nrecv

   character(char_len) :: &
      runtype         

   type(seq_infodata_type), pointer :: &
      infodata   

!=======================================================================

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_mct
!
! !INTERFACE:
  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
!
! !DESCRIPTION:
! Initialize POP 
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o, o2x_o
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer(int_kind) ::  &
       OCNID,       &
       mpicom_o,    &
       lsize,       &
       start_ymd,   &
       start_tod,   &
       start_year,  &
       start_day,   &
       start_month, &
       start_hour,  &
       iyear,       &
       ocn_cpl_dt,  &
       pop_cpl_dt,  &
       shrlogunit,  &  ! old values
       shrloglev       ! old values

    type(mct_gsMap), pointer :: &
       gsMap_o

    type(mct_gGrid), pointer :: &
       dom_o

    integer (POP_i4) :: &
       errorCode         ! error code

    integer (int_kind) :: &
       nThreads

    integer (int_kind) :: iam,ierr 
    character(len=32)  :: starttype          ! infodata start type

#ifdef _OPENMP
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif

     integer :: lbnum

!-----------------------------------------------------------------------
!
!  set cdata pointers
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

#ifdef _OPENMP
   nThreads = omp_get_max_threads()
#endif
   call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpicom_o, &
        gsMap=gsMap_o, dom=dom_o, infodata=infodata)

   POP_MCT_OCNID   =  OCNID
   POP_MCT_gsMap_o => gsMap_o
   POP_MCT_dom_o   => dom_o

#if (defined _MEMTRACE)
    call MPI_comm_rank(mpicom_o,iam,ierr)
    if(iam == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:start::',lbnum) 
    endif
#endif


    ! The following communicator module variable will be utilize in init_communicate that
    ! is called by initial - this is done to make the code backwards compatible

    mpi_communicator_ocn = mpicom_o

!-----------------------------------------------------------------------
!
!  initialize the model run 
!
!-----------------------------------------------------------------------

    call POP_CplIndicesSet()

    call seq_infodata_GetData( infodata, case_name=runid )
   
    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(stdout,*) 'ocn_comp_mct ERROR: unknown starttype'
       call exit_POP(sigAbort,' ocn_comp_mct ERROR: unknown starttype')
    end if

    !TODO: check for consistency of pop runid and runtype with seq_infodata
   
!-----------------------------------------------------------------------
!
!  first initializaiton phase of pop2
!  initialize pop2 because grid information is needed for
!  creation of GSMap_ocn.
!  call pop initialization routines in two stages (needed for backwards
!  compatiblity with cpl6 concurrent system
!
!-----------------------------------------------------------------------

   inst_name   = seq_comm_name(OCNID)
   inst_index  = seq_comm_inst(OCNID)
   inst_suffix = seq_comm_suffix(OCNID)

   call t_startf ('pop_init')
   call POP_Initialize1(errorCode)

!-----------------------------------------------------------------------
!
!  register non-standard incoming fields
!
!-----------------------------------------------------------------------

   if (index_x2o_Sa_co2prog > 0) then
      call named_field_register('ATM_CO2_PROG', ATM_CO2_PROG_nf_ind)
   endif
   if (index_x2o_Sa_co2diag > 0) then
      call named_field_register('ATM_CO2_DIAG', ATM_CO2_DIAG_nf_ind)
   endif
   call register_string('pop_init_coupled')
   call flushm (stdout)

!-----------------------------------------------------------------------
!
!  second initialization phase of pop2
!
!-----------------------------------------------------------------------

   call POP_Initialize2(errorCode)

!-----------------------------------------------------------------------
!
!  initialize time-stamp information
!
!-----------------------------------------------------------------------

   call ccsm_char_date_and_time

   call t_stopf ('pop_init')

!----------------------------------------------------------------------------
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)
   
!-----------------------------------------------------------------------
!
!  check for consistency of pop and sync clock initial time
!
!-----------------------------------------------------------------------

   if (runtype == 'initial') then
      call seq_timemgr_EClockGetData(EClock, &
           start_ymd=start_ymd, start_tod=start_tod)
      call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)

      if (iyear0 /= start_year) then
	 if(master_task == my_task)   then
            call document ('ocn_init_mct', 'iyear0     ', iyear0)
            call document ('ocn_init_mct', 'start_year ', start_year)
         endif
         call exit_POP(sigAbort,' iyear0 does not match start_year')
      end if
      if (imonth0 /= start_month) then
	 if(master_task == my_task)   then
            call document ('ocn_init_mct', 'imonth0     ', imonth0)
            call document ('ocn_init_mct', 'start_month ', start_month)
         endif
         call exit_POP(sigAbort,' imonth0 does not match start_year')
      end if
      if (iday0 /= start_day) then
	 if(master_task == my_task)   then
            call document ('ocn_init_mct', 'iday0     ', iday0)
            call document ('ocn_init_mct', 'start_day ', start_day)
         endif
      end if
!#ifndef _HIRES 
!      if (seconds_this_day /= start_tod) then
!         call document ('ocn_init_mct', 'sec0     ', seconds_this_day)
!         call document ('ocn_init_mct', 'start_tod ', start_tod)
!         call exit_POP(sigAbort,' sec0 does not start_tod')
!      end if
!#endif
   end if

!-----------------------------------------------------------------------
!
!  initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

   call t_startf ('pop_mct_init')

   call ocn_SetGSMap_mct( mpicom_o, OCNID, GSMap_o ) 	
   lsize = mct_gsMap_lsize(gsMap_o, mpicom_o)

   ! Initialize mct ocn domain (needs ocn initialization info)
   
   call ocn_domain_mct( lsize, gsMap_o, dom_o )
   
   ! Inialize mct attribute vectors
   
   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
   call mct_aVect_zero(x2o_o)
   
   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize) 
   call mct_aVect_zero(o2x_o)
   
   nsend = mct_avect_nRattr(o2x_o)
   nrecv = mct_avect_nRattr(x2o_o)
   allocate (SBUFF_SUM(nx_block,ny_block,max_blocks_clinic,nsend))

!-----------------------------------------------------------------------
!
!  Initialize flags and shortwave absorption profile
!  Note that these cpl_write_xxx flags have no freqency options
!  set; therefore, they will retain a default value of .false.
!  unless they are explicitly set .true.  at the appropriate times
!
!-----------------------------------------------------------------------

   call init_time_flag('cpl_write_restart',cpl_write_restart, owner = 'ocn_init_mct')
   call init_time_flag('cpl_write_history',cpl_write_history, owner = 'ocn_init_mct')
   call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg,    owner = 'ocn_init_mct')
   call init_time_flag('cpl_diag_global'  ,cpl_diag_global,   owner = 'ocn_init_mct')
   call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp,   owner = 'ocn_init_mct')

   lsmft_avail = .true.
   tlast_coupled = c0

!-----------------------------------------------------------------------
!
!   initialize necessary  coupling info
!
!-----------------------------------------------------------------------

    call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
    pop_cpl_dt = seconds_in_day / ncouple_per_day
    if (pop_cpl_dt /= ocn_cpl_dt) then
       write(stdout,*)'pop_cpl_dt= ',pop_cpl_dt, &
                     ' ocn_cpl_dt= ',ocn_cpl_dt   
       call exit_POP(sigAbort,'ERROR pop_cpl_dt and ocn_cpl_dt must be identical')
    end if

!-----------------------------------------------------------------------
!
!  send intial state to driver
!
!-----------------------------------------------------------------------

   if ( lsend_precip_fact )  then
      call seq_infodata_PutData( infodata, precip_fact=precip_fact)
   end if
   call pop_sum_buffer

   call ocn_export_mct(o2x_o, errorCode)  
   if (errorCode /= POP_Success) then
      call POP_ErrorPrint(errorCode)
      call exit_POP(sigAbort, 'ERROR in ocn_export_mct')
   endif

   call t_stopf ('pop_mct_init')

   call seq_infodata_PutData( infodata, &
        ocn_nx = nx_global , ocn_ny = ny_global)
   call seq_infodata_PutData( infodata, &
	ocn_prognostic=.true., ocnrof_prognostic=.true.)

!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(iam  == 0) then
!        write(6,*) 'ocn_init_mct:end::'
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:end::',lbnum) 
        call memmon_reset_addr()
    endif
#endif

!-----------------------------------------------------------------------
!
!  document orbital parameters
!
!-----------------------------------------------------------------------

   if (registry_match('qsw_distrb_iopt_cosz')) then
   call seq_infodata_GetData(infodata, &
      orb_eccen=orb_eccen, orb_mvelpp=orb_mvelpp, orb_lambm0=orb_lambm0, orb_obliqr=orb_obliqr)

     write(stdout,*) ' '
     call document ('ocn_import_mct', 'orb_eccen   ',  orb_eccen)
     call document ('ocn_import_mct', 'orb_mvelpp  ',  orb_mvelpp)
     call document ('ocn_import_mct', 'orb_lambm0  ',  orb_lambm0)
     call document ('ocn_import_mct', 'orb_obliqr  ',  orb_obliqr)
    endif

!-----------------------------------------------------------------------
!
!  Now document all time flags, because this is the last step of pop2 
!    initialization
!
!-----------------------------------------------------------------------

   call document_time_flags

!-----------------------------------------------------------------------
!
!  output delimiter to log file
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(" End of initialization")')
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      call POP_IOUnitsFlush(POP_stdout)
#ifdef CCSMCOUPLED
      call POP_IOUnitsFlush(stdout)
#endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Run POP for a coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer(int_kind) :: & 
         errorCode           ! error flag

    integer(int_kind) :: &
         ymd, &          ! POP2 current date (YYYYMMDD)
         tod, &          ! POP2 current time of day (sec)
         ymd_sync, &     ! Sync clock current date (YYYYMMDD)
         tod_sync, &     ! Sync clcok current time of day (sec)
         shrlogunit,  &  ! old values
         shrloglev       ! old values

    character(len=char_len_long) :: &
         fname

    character(len=*), parameter  :: &
         SubName = "ocn_run_mct"

    type(seq_infodata_type), pointer :: &
       infodata   ! Input init object

    logical :: &
         lcoupled,   &  ! temporary
         rstwr,      &  ! true => write restart at end of day
         first_time = .true.

    character (char_len)  :: message

    integer(int_kind) :: info_debug

    integer :: lbnum

#if (defined _MEMTRACE)
    if(my_task == 0 ) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':start::',lbnum) 
    endif
#endif

!-----------------------------------------------------------------------
!
!  start up the main timer
!
!-----------------------------------------------------------------------

   call timer_start(timer_total)


!-----------------------------------------------------------------------
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------

    errorCode = POP_Success

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    call seq_cdata_setptrs(cdata_o, infodata=infodata)

!----------------------------------------------------------------------------
!
! restart flag (rstwr) will assume only an eod restart for now
!
!----------------------------------------------------------------------------

    call seq_infodata_GetData( infodata, info_debug=info_debug)
    if (info_debug >= 2) then
       ldiag_cpl = .true. 
       call register_string ('info_debug_ge2')
    else
       ldiag_cpl = .false.
    endif

    rstwr = seq_timemgr_RestartAlarmIsOn(EClock)
    if (rstwr) then
       call override_time_flag(cpl_write_restart,value=.true.)
       call ccsm_char_date_and_time ! set time_management module vars cyear, cmonth, ..
       write(message,'(6a)') 'driver requests restart file at eod  ',  &
            cyear,'/',cmonth,'/',cday
       call document ('ocn_comp_mct(run):', message)
    endif

!-----------------------------------------------------------------------
!
!  advance the model in time over coupling interval
!  write restart dumps and archiving
!
!-----------------------------------------------------------------------

    ! Note that all ocean time flags are evaluated each timestep in time_manager
    ! ocn_import_mct is analogous to pop_unpack_fluxes_from_coupler in cpl6 
    ! ocn_export_mct is analogous to prepare_send_to_coupler in cpl6
    ! tlast_coupled is set to zero at the end of ocn_export_mct

    advance: do 

       ! obtain import state from driver
       if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then
          call ocn_import_mct(x2o_o, errorCode)   

          if (errorCode /= POP_Success) then
             call POP_ErrorPrint(errorCode)
             call exit_POP(sigAbort, 'ERROR in step')
          endif

          call pop_set_coupled_forcing 
       end if
       
       call step(errorCode)

       if (errorCode /= POP_Success) then
          call POP_ErrorPrint(errorCode)
          call exit_POP(sigAbort, 'ERROR in step')
       endif

       if (check_KE(100.0_r8)) then
          !*** exit if energy is blowing
          call output_driver
          call exit_POP(sigAbort,'ERROR: k.e. > 100 ')
       endif
       call output_driver
       
       ! return export state to driver
       call pop_sum_buffer()  
       if (check_time_flag(cpl_ts)) then
          call ocn_export_mct(o2x_o, errorCode)
          if (errorCode /= POP_Success) then
             call POP_ErrorPrint(errorCode)
             call exit_POP(sigAbort, 'ERROR in ocn_export_mct')
          endif

          exit advance
       end if
       
    enddo advance

    if ( lsend_precip_fact ) then
       call seq_infodata_PutData( infodata, precip_fact=precip_fact )
    end if
    
!--------------------------------------------------------------------
!
! check that internal clock is in sync with master clock
!
!--------------------------------------------------------------------

    ymd = iyear*10000 + imonth*100 + iday
    tod = ihour*seconds_in_hour + iminute*seconds_in_minute + isecond
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync )
       write(stdout,*)' pop2 ymd=',ymd     ,'  pop2 tod= ',tod
       write(stdout,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       write(stdout,*)' Internal pop2 clock not in sync with Sync Clock'
       call shr_sys_abort( SubName// &
          ":: Internal pop2 clock not in sync with Sync Clock")
    end if
   
!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)

   call timer_stop(timer_total)

#if (defined _MEMTRACE)
    if(my_task == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',SubName//':end::',lbnum) 
       call memmon_reset_addr()
    endif
#endif
!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_run_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Finalize POP
!
! !USES:
    use POP_FinalMod
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
! !REVISION HISTORY:
! Author: Fei Liu
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode         ! error code

!-----------------------------------------------------------------------

  call POP_Final(errorCode)

  end subroutine ocn_final_mct

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

 subroutine ocn_SetGSMap_mct( mpicom_ocn, OCNID, gsMap_ocn )

! !DESCRIPTION:
!  This routine mct global seg maps for the pop decomposition
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer,allocatable :: &
      gindex(:)

    integer (int_kind) ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier

    type (block) ::       &
      this_block          ! block information for current block

!-----------------------------------------------------------------------
!  Build the POP grid numbering for MCT
!  NOTE:  Numbering scheme is: West to East and South to North starting
!  at the south pole.  Should be the same as what's used in SCRIP
!-----------------------------------------------------------------------

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
       enddo
       enddo
    enddo
    lsize = n

! not correct for padding, use "n" above
!    lsize = block_size_x*block_size_y*nblocks_clinic
    gsize = nx_global*ny_global
    allocate(gindex(lsize),stat=ier)

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          gindex(n) = (this_block%j_glob(j)-1)*(nx_global) + this_block%i_glob(i) 
       enddo
       enddo
    enddo

    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_SetGSMap_mct

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_mct
! !INTERFACE:

 subroutine ocn_domain_mct( lsize, gsMap_o, dom_o )

! !DESCRIPTION:
!  This routine mct global seg maps for the pop decomposition
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o     

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer, pointer :: &
      idata(:)

    real(r8), pointer :: &
      data(:)

    integer (int_kind) ::   &
      i,j, k, n, iblock, &
      ier

    type (block) ::       &
      this_block          ! block information for current block

!-------------------------------------------------------------------
!
!  initialize mct domain type, lat/lon in degrees,
!  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_o, my_task, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain (numbering scheme is: West to East and South to North to South pole)
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize) 

!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          data(n) = TLOND(i,j,iblock)
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize) 

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          data(n) = TLATD(i,j,iblock)
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize) 

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          data(n) = TAREA(i,j,iblock)/(radius*radius)
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"area",data,lsize) 

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          data(n) = float(KMT(i,j,iblock)) 
          if (data(n) > 1.0_r8) data(n) = 1.0_r8
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_mct


!***********************************************************************
!BOP
! !IROUTINE: ocn_import_mct
! !INTERFACE:

 subroutine ocn_import_mct(x2o_o, errorCode)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  roff   -- river runoff flux                        (kg/m2/s)
!    o  ioff   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect)   , intent(inout) :: x2o_o

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::   &
      label,                 &
      message
 
   integer (int_kind) ::  &
      i,j,k,n,iblock

   real (r8), dimension(nx_block,ny_block) ::  &
      WORKB

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2        ! local work space

   real (r8) ::  &
      m2percm2,  &
      gsum

   type (block) :: this_block ! local block info

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   WORK1 = c0
   WORK2 = c0

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         WORK1(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_taux,n)
         WORK2(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_tauy,n)
      enddo
      enddo
   enddo ! iblock

   !***
   !*** do NOT perform halo updates now, because vector updates must
   !***   be done after the rotation is completed.
   !***

!-----------------------------------------------------------------------
!
!  rotate true zonal/meridional wind stress into local coordinates,
!  convert to dyne/cm**2, and shift SMFT to U grid
!
!  halo updates are performed in subroutine rotate_wind_stress, 
!  following the rotation
!
!-----------------------------------------------------------------------

      call rotate_wind_stress(WORK1, WORK2)

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
!
!  unpack and distribute fresh water flux and salt flux
!
!  NOTE: if there are code changes associated with changing the names or
!        the number of fluxes received from the coupler, then subroutine
!        update_ghost_cells_coupler_fluxes will need to be modified also
!
!-----------------------------------------------------------------------


      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         SNOW_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Faxa_snow,n)
         WORKB (i,j       ) = x2o_o%rAttr(index_x2o_Faxa_rain,n)
         EVAP_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_evap,n)
         MELT_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Fioi_meltw,n)
         ROFF_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Forr_roff,n)
         IOFF_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Forr_ioff,n)
         SALT_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Fioi_salt,n)

         PREC_F(i,j,iblock) = WORKB(i,j) + SNOW_F(i,j,iblock)    ! rain + snow

         WORKB(i,j        ) = x2o_o%rAttr(index_x2o_Foxx_swnet,n)
         SHF_QSW(i,j,iblock) = WORKB(i,j)*  &
            RCALCT(i,j,iblock)*hflux_factor  !  convert from W/m**2
         SENH_F(i,j,iblock)  = x2o_o%rAttr(index_x2o_Foxx_sen,n)
         LWUP_F(i,j,iblock)  = x2o_o%rAttr(index_x2o_Foxx_lwup,n)
         LWDN_F(i,j,iblock)  = x2o_o%rAttr(index_x2o_Faxa_lwdn,n)
         MELTH_F(i,j,iblock) = x2o_o%rAttr(index_x2o_Fioi_melth,n)

         WORKB(i,j       ) = x2o_o%rAttr(index_x2o_Si_ifrac,n)
         IFRAC(i,j,iblock) = WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from Pa to dynes/cm**2
         WORKB(i,j       ) = x2o_o%rAttr(index_x2o_Sa_pslv,n)
         ATM_PRESS(i,j,iblock) = c10 * WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from m**2/s**2 to cm**2/s**2
         WORKB(i,j       ) = x2o_o%rAttr(index_x2o_So_duu10n,n)
         U10_SQR(i,j,iblock) = cmperm * cmperm * WORKB(i,j) * RCALCT(i,j,iblock)

      enddo
      enddo

   enddo

!-----------------------------------------------------------------------
!
!  incoming data quality control
!
!-----------------------------------------------------------------------
#ifdef CCSMCOUPLED
      if ( any(IOFF_F < c0) ) then
        write(message, "(A,1x,e10.3,A)") 'Error: incoming IOFF_F has min value', &
                                      minval(IOFF_F), '; value can not be negative.'
        ! call shr_sys_abort ('Error: incoming IOFF_F is negative')
        call shr_sys_abort (trim(message))
      endif
#endif


!-----------------------------------------------------------------------
!
!  update ghost cells for fluxes received from the coupler
!
!-----------------------------------------------------------------------

   call update_ghost_cells_coupler_fluxes(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ocn_import_mct: error in update_ghost_cells_coupler_fluxes')
      return
   endif

!-----------------------------------------------------------------------
!
!  unpack atmospheric CO2
!
!-----------------------------------------------------------------------

   if (index_x2o_Sa_co2prog > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o_o%rAttr(index_x2o_Sa_co2prog,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating PROG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_PROG_nf_ind, WORK1)
   endif

   if (index_x2o_Sa_co2diag > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o_o%rAttr(index_x2o_Sa_co2diag,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating DIAG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_DIAG_nf_ind, WORK1)
   endif
 
!-----------------------------------------------------------------------
!
!  receive orbital parameters
!
!-----------------------------------------------------------------------

   call seq_infodata_GetData(infodata, &
      orb_eccen=orb_eccen, orb_mvelpp=orb_mvelpp, orb_lambm0=orb_lambm0, orb_obliqr=orb_obliqr)


!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then

     write(message,'(6a,1x,5a)')  &
         ' Global averages of fluxes received from cpl at ',  &
           cyear,'/',cmonth ,'/',cday,  chour,':',cminute,':',csecond
     call document ('pop_recv_from_coupler', trim(message))
 
     m2percm2  = mpercm*mpercm
     do k = 1,nrecv

         n = 0
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)

            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORK1(i,j,iblock) = x2o_o%rAttr(k,n)  ! mult. by TAREA in global_sum_prod
            enddo
            enddo
         enddo

         gsum = global_sum_prod(WORK1 , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
         if (my_task == master_task) then
            call seq_flds_getField(label,k,seq_flds_x2o_fields)
            write(stdout,1100)'ocn','recv', label ,gsum
            call shr_sys_flush(stdout)
         endif
      enddo
   endif


1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_mct
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export_mct(o2x_o, errorCode)   

! !DESCRIPTION:
!  This routine calls the routines necessary to send pop fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)   , intent(inout) :: o2x_o

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n, iblock
           
   character (char_len)    :: label
 
   integer (int_kind) ::  &
      i,j,k

   real (r8), dimension(nx_block,ny_block) ::   &
      WORK1, WORK2,      &! local work space
      WORK3, WORK4

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
        WORKA               ! local work space with full block dimension

   real (r8) ::   &
      m2percm2,   &
      gsum

   type (block) :: this_block ! local block info

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)

     call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_u),iblock)
     call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_v),iblock)

     WORK1 = (WORK3*cos(ANGLET(:,:,iblock))+WORK4*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled
     WORK2 = (WORK4*cos(ANGLET(:,:,iblock))-WORK3*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled

     do j=this_block%jb,this_block%je
     do i=this_block%ib,this_block%ie
        n = n + 1
        o2x_o%rAttr(index_o2x_So_u,n) = WORK1(i,j)
        o2x_o%rAttr(index_o2x_So_v,n) = WORK2(i,j)
     enddo
     enddo
  enddo

!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x_o%rAttr(index_o2x_So_t,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_t)/tlast_coupled + T0_Kelvin
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x_o%rAttr(index_o2x_So_s,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_s)*salt_to_ppt/tlast_coupled
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx),iblock)
      call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy),iblock)
 
      WORK1 = (WORK3*cos(ANGLET(:,:,iblock)) + WORK4*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET(:,:,iblock)) - WORK3*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x_o%rAttr(index_o2x_So_dhdx,n) = WORK1(i,j)
         o2x_o%rAttr(index_o2x_So_dhdy,n) = WORK2(i,j)
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in ice.F
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x_o%rAttr(index_o2x_Fioo_q,n) = QFLUX(i,j,iblock)
      enddo
      enddo
   enddo

   tlast_ice = c0
   AQICE     = c0
   QICE      = c0

!-----------------------------------------------------------------------
!
!     pack co2 flux, if requested (kg CO2/m^2/s)
!     units conversion occurs where co2 flux is computed
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x_o%rAttr(index_o2x_Faoo_fco2_ocn,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_Faoo_fco2_ocn)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then
      call ccsm_char_date_and_time
!DEBUG      write(message,'(6a,1x,5a)')' Global averages of fluxes sent to cpl at ', &
!DEBUG           cyear,'/',cmonth, '/',cday,  chour,':',cminute,':',csecond
!DEBUG      call document ('pop_send_to_coupler', message)
      write(stdout,*)'pop_send_to_coupler'

      m2percm2  = mpercm*mpercm
      do k = 1,nsend
        n = 0
        do iblock = 1, nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           do j=this_block%jb,this_block%je
           do i=this_block%ib,this_block%ie
              n = n + 1
              WORKA(i,j,iblock) = o2x_o%rAttr(k,n)
           enddo
           enddo
        enddo

        call POP_HaloUpdate(WORKA,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
       
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'ocn_export_mct: error updating halo for state')
            return
         endif

        gsum = global_sum_prod(WORKA , TAREA, distrb_clinic, &
                                   field_loc_center, RCALCT)*m2percm2
        if (my_task == master_task) then
           call seq_flds_getField(label,k,seq_flds_o2x_fields)
           write(stdout,1100)'ocn','send', label ,gsum
        endif
      enddo ! k
      if (my_task == master_task) call shr_sys_flush(stdout)
   endif

1100 format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)
    
    tlast_coupled = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_mct

!***********************************************************************

!BOP
! !IROUTINE: pop_sum_buffer
! !INTERFACE:

 subroutine pop_sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
!
! !REVISION HISTORY:
!  same as module
! 
!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                ! local work arrays

   real (r8) ::   &
      delt,             & ! time interval since last step
      delt_last           ! time interval for previous step

   integer (int_kind) :: &
      iblock,           & ! block index
      sflux_co2_nf_ind = 0! named field index of fco2

   logical (log_kind) :: &
      first = .true.      ! only true for first call

   save first

!-----------------------------------------------------------------------
!
!  zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

   if (tlast_coupled == c0) SBUFF_SUM = c0
   WORK = c0

!-----------------------------------------------------------------------
!
!  update time since last coupling
!
!-----------------------------------------------------------------------

   if (avg_ts .or. back_to_back) then
      delt = p5*dtt
   else
      delt =    dtt
   endif
   tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!
!  allow for fco2 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      if (sflux_co2_nf_ind == 0) then
         call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!
!  accumulate sums of U,V,T,S and GRADP
!  accumulate sum of co2 flux, if requested
!     implicitly use zero flux if fco2 field not registered yet
!  ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
   SBUFF_SUM(:,:,iblock,index_o2x_So_u) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_u) + delt*  &
                                   UVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_v) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_v) + delt*  &
                                   VVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_t ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_t ) + delt*  &
                                   TRACER(:,:,1,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_s ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_s ) + delt*  &
                                   TRACER(:,:,1,2,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) + delt*  &
                                   GRADPX(:,:,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) + delt*  &
                                   GRADPY(:,:,curtime,iblock)

   if (index_o2x_Faoo_fco2_ocn > 0 .and. sflux_co2_nf_ind > 0) then
      call named_field_get(sflux_co2_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) = &
         SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) + delt_last*WORK(:,:,iblock)
   endif

   enddo
   !$OMP END PARALLEL DO

   first = .false.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_sum_buffer
 
!***********************************************************************

end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
