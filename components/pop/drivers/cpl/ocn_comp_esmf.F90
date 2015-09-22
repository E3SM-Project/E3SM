module ocn_comp_esmf

#ifdef ESMF_INTERFACE
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_esmf
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
   use POP_CplIndices
   use POP_KindsMod
   use POP_ErrorMod
   use POP_InitMod,       only: POP_Initialize1, POP_Initialize2
   use POP_InitMod,       only: timer_total, cpl_ts 

   use esmf
   use esmfshr_util_mod, only : esmfshr_util_StateArrayDestroy
   use esmfshr_util_mod, only : esmfshr_util_ArrayGetIndex
   use esmfshr_util_mod, only : esmfshr_util_ArrayGetSize
   use esmf2mct_mod    , only : esmf2mct_init

   use seq_flds_mod
   use seq_timemgr_mod
   use seq_infodata_mod,only : seq_infodata_start_type_cont, &
                               seq_infodata_start_type_brnch, seq_infodata_start_type_start
   use seq_comm_mct,    only : seq_comm_suffix, seq_comm_inst, seq_comm_name

   use shr_file_mod 
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod
   use perf_mod

   use kinds_mod,         only: int_kind, r8
   use ocn_communicator,  only: mpi_communicator_ocn
   use ocn_import_export, only: ocn_import, ocn_export, POP_sum_buffer
   use ocn_import_export, only: SBUFF_SUM, tlast_coupled
   use communicate,       only: my_task, master_task
   use constants
   use blocks
   use domain,            only: distrb_clinic, POP_haloClinic
   use exit_mod
   use forcing_shf,       only: SHF_QSW
   use forcing_sfwf,      only: lsend_precip_fact, precip_fact
   use forcing_fields
   use forcing_coupled,   only: ncouple_per_day, pop_set_coupled_forcing, pop_init_coupled
   use forcing_coupled,   only: orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
   use grid,              only: TLAT, TLON, KMT
   use global_reductions, only: global_sum_prod
   use io_tools,          only: document
   use named_field_mod,   only: named_field_register, named_field_get_index
   use named_field_mod,   only: named_field_set, named_field_get
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
  public :: ocn_register_esmf
  public :: ocn_init_esmf
  public :: ocn_run_esmf
  public :: ocn_final_esmf
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ocn_DistGrid_esmf
  private :: ocn_domain_esmf
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

   character(char_len) :: &
      runtype         

   integer (int_kind)  ::   &
      nsend, nrecv

   ! Following is needed for call to shr_str_data in ecosys_mod
   ! These variables are pointed to in POP_MCT_vars_mod
   type(mct_gsMap), target  :: gsmap_o
   type(mct_gGrid), target  :: dom_o
   integer                  :: OCNID      

!=======================================================================

contains

!=======================================================================

  subroutine ocn_register_esmf(comp, rc)
    implicit none
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In ocn register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         ocn_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         ocn_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         ocn_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine ocn_register_esmf

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_esmf
!
! !INTERFACE:
  subroutine ocn_init_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Initialize POP 
!
! !INPUT/OUTPUT PARAMETERS:
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)                        :: rc
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer(int_kind) ::  &
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
       
    integer (POP_i4) :: &
       errorCode         ! error code

    integer (int_kind) :: &
       nThreads

    character(len=32)  :: starttype          ! infodata start type

#ifdef _OPENMP
    integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                              ! concurrently in a single parallel region
#endif

    integer                :: mpicom_ocn, mpicom_vm, lsize, gsize
    type(ESMF_ArraySpec)   :: arrayspec
    type(ESMF_DistGrid)    :: distgrid
    type(ESMF_Array)       :: o2x, x2o, dom
    type(ESMF_VM)          :: vm
    character(ESMF_MAXSTR) :: convCIM, purpComp
    integer                :: nfields
    real(R8), pointer      :: fptr (:,:)          ! data pointer into ESMF array

!-----------------------------------------------------------------------
!
!  set cdata pointers
!
!-----------------------------------------------------------------------

    call POP_CplIndicesSet()

    rc = ESMF_SUCCESS
 
  
    ! duplicate the mpi communicator from the current VM 
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_dup(mpicom_vm, mpicom_ocn, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    errorCode = POP_Success
    print * , 'after mpicom setup'
    
    ! Initialize pop id
    call ESMF_AttributeGet(export_state, name="ID", value=OCNID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

#ifdef _OPENMP
    nThreads = omp_get_max_threads()
#endif

#if (defined _MEMTRACE)
    call MPI_comm_rank(mpicom_ocn,iam,ierr)
    if(iam == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','ocn_init_esmf:start::',lbnum) 
    endif
#endif
    
    ! The following communicator module variable will be utilize in init_communicate that
    ! is called by initial - this is done to make the code backwards compatible
    
    mpi_communicator_ocn = mpicom_ocn

!-----------------------------------------------------------------------
!
!  initialize the model run 
!
!-----------------------------------------------------------------------

    call ESMF_AttributeGet(export_state, name="case_name", value=runid, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(stdout,*) 'ocn_comp_esmf ERROR: unknown starttype'
       call exit_POP(sigAbort,' ocn_comp_esmf ERROR: unknown starttype')
    end if

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

   print * , 'begin pop init1'
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

   print * , 'begin pop init2'
   call POP_Initialize2(errorCode)
   print * , 'end pop init2'

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

      ! Check for consistency
      if (iyear0 /= start_year) then
         if(master_task == my_task)   then
             call document ('ocn_init_esmf', 'iyear0     ', iyear0)
             call document ('ocn_init_esmf', 'start_year ', start_year)
         endif
         call exit_POP(sigAbort,' iyear0 does not match start_year')
      end if
      if (imonth0 /= start_month) then
    	 if(master_task == my_task)   then
            call document ('ocn_init_esmf', 'imonth0     ', imonth0)
            call document ('ocn_init_esmf', 'start_month ', start_month)
         endif
         call exit_POP(sigAbort,' imonth0 does not match start_year')
      end if
      if (iday0 /= start_day) then
    	 if(master_task == my_task)   then
            call document ('ocn_init_esmf', 'iday0     ', iday0)
            call document ('ocn_init_esmf', 'start_day ', start_day)
         endif
      end if
   end if

!-----------------------------------------------------------------------
!
!  initialize distgrid, domain, and arrays
!
!-----------------------------------------------------------------------

    call t_startf ('pop_esmf_init')

    !-----------------------------------------
    !  Set arrayspec for dom, o2x and x2o
    !-----------------------------------------
    
    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    distgrid = ocn_distgrid_esmf(gsize)

    call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------
    ! Create dom 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_dom_fields))

    dom = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set values of dom (needs ocn initialization info)

    call ocn_domain_esmf(dom)
   
    !----------------------------------------- 
    !  Create o2x 
    !-----------------------------------------

    ! 1d undistributed index of fields, 2d is packed data

    nfields = shr_string_listGetNum(trim(seq_flds_o2x_fields))

    o2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(o2x, name="mct_names", value=trim(seq_flds_o2x_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    !  Create x2o 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_x2o_fields))

    x2o = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2o, name="mct_names", value=trim(seq_flds_x2o_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    ! Add esmf arrays to import and export state 
    !-----------------------------------------

    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/o2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2o/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   
    call esmf2mct_init(distgrid, OCNID, gsmap_o, mpicom_ocn, gsize=gsize, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom, dom_o, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    POP_MCT_gsMap_o => gsMap_o
    POP_MCT_dom_o   => dom_o
    POP_MCT_OCNID   =  OCNID

    call esmfshr_util_ArrayGetSize(o2x, lsize1=nsend, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate (SBUFF_SUM(nx_block,ny_block,max_blocks_clinic, nsend))

!-----------------------------------------------------------------------
!
!  Initialize flags and shortwave absorption profile
!  Note that these cpl_write_xxx flags have no freqency options
!  set; therefore, they will retain a default value of .false.
!  unless they are explicitly set .true.  at the appropriate times
!
!-----------------------------------------------------------------------

   call init_time_flag('cpl_write_restart',cpl_write_restart, owner = 'ocn_init_esmf')
   call init_time_flag('cpl_write_history',cpl_write_history, owner = 'ocn_init_esmf')
   call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg,    owner = 'ocn_init_esmf')
   call init_time_flag('cpl_diag_global'  ,cpl_diag_global,   owner = 'ocn_init_esmf')
   call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp,   owner = 'ocn_init_esmf')

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

   print * , 'begin esmf export'
   if ( lsend_precip_fact )  then
      call ESMF_AttributeSet(export_state, name="precip_fact", value=precip_fact, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   end if

   call pop_sum_buffer()

   call ESMF_ArrayGet(o2x, localDe=0, farrayPtr=fptr, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   fptr(:,:) = 0._r8
   call ocn_export(fptr, ldiag_cpl, rc)

   errorCode = rc
   if (errorCode /= POP_Success) then
      call POP_ErrorPrint(errorCode)
      call exit_POP(sigAbort, 'ERROR in ocn_export')
   endif

   call t_stopf ('pop_esmf_init')

   call ESMF_AttributeSet(export_state, name="ocn_prognostic", value=.true., rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(export_state, name="ocnrof_prognostic", value=.true., rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(export_state, name="ocn_nx", value=nx_global, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(export_state, name="ocn_ny", value=ny_global, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

!----------------------------------------------------------------------------
!
! Reset shr logging to original values
!
!----------------------------------------------------------------------------

   print * , 'done esmf export'
   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
   if(iam  == 0) then
      !        write(6,*) 'ocn_init_esmf:end::'
      lbnum=1
      call memmon_dump_fort('memmon.out','ocn_init_esmf:end::',lbnum) 
      call memmon_reset_addr()
   endif
#endif

!-----------------------------------------------------------------------
!
!  document orbital parameters
!
!-----------------------------------------------------------------------

   if (registry_match('qsw_distrb_iopt_cosz')) then
     call ESMF_AttributeGet(export_state, name="orb_eccen", value=orb_eccen, rc=rc)
     if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

     call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=orb_mvelpp, rc=rc)
     if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

     call ESMF_AttributeGet(export_state, name="orb_lambm0", value=orb_lambm0, rc=rc)
     if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

     call ESMF_AttributeGet(export_state, name="orb_obliqr", value=orb_obliqr, rc=rc)
     if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

     write(stdout,*) ' '
     call document ('ocn_import_esmf', 'orb_eccen   ',  orb_eccen)
     call document ('ocn_import_esmf', 'orb_mvelpp  ',  orb_mvelpp)
     call document ('ocn_import_esmf', 'orb_lambm0  ',  orb_lambm0)
     call document ('ocn_import_esmf', 'orb_obliqr  ',  orb_obliqr)
    endif

!-----------------------------------------------------------------------
!
!  Now document all time flags, because this is the last step of pop2 
!    initialization
!
!-----------------------------------------------------------------------

   call document_time_flags

   print * , 'done esmf init'

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

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "POP", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "LongName", &
                           "Parallel Ocean Program", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "Description", &
                  "The ocean component of the CESM1.0 is the Parallel " // &
                  "Ocean Program version 2 (POP2).  This model is based " // &
                  "on the POP version 2.1 of the Los Alamos National " // &
                  "Laboratory; however, it includes many physical and " // &
                  "software developments incorporated by the members " // &
                  "of the Ocean Model Working Group (see the notable " // &
                  "improvements page for these developments).", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "Release Date", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ModelType", "Ocean", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "Susan Bates", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "bates@ucar.edu", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_esmf

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_esmf
!
! !INTERFACE:
subroutine ocn_run_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Run POP for a coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)                        :: rc

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
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

    character(len=*), parameter  :: &
         SubName = "ocn_run_esmf"

    logical :: &
         rstwr           ! true => write restart at end of day

    character (char_len)  :: message

    integer(int_kind) :: info_debug

    type(ESMF_Array) :: o2x, x2o

    real(R8), pointer ::  &
         fptr (:,:)          ! data pointer into ESMF array
!-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

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
!
! reset shr logging to my log file
!
!----------------------------------------------------------------------------

    errorCode = POP_Success

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

!----------------------------------------------------------------------------
!
! restart flag (rstwr) will assume only an eod restart for now
!
!----------------------------------------------------------------------------

    call ESMF_AttributeGet(export_state, name="info_debug", value=info_debug, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

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
       call document ('ocn_comp_esmf(run):', message)
    endif

!-----------------------------------------------------------------------
!
!  advance the model in time over coupling interval
!  write restart dumps and archiving
!
!-----------------------------------------------------------------------

    ! Note that all ocean time flags are evaluated each timestep in time_manager
    ! tlast_coupled is set to zero at the end of ocn_export

    advance: do 

       if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then

          ! Obtain input from driver from import state
          call ESMF_StateGet(import_state, itemName="x2d", array=x2o, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          call ESMF_ArrayGet(x2o, localDe=0, farrayPtr=fptr, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          call ocn_import(fptr, ldiag_cpl, rc)

          ! Get orbital values from export state
          call ESMF_AttributeGet(export_state, name="orb_eccen", value=orb_eccen, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          
          call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=orb_mvelpp, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          
          call ESMF_AttributeGet(export_state, name="orb_lambm0", value=orb_lambm0, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          
          call ESMF_AttributeGet(export_state, name="orb_obliqr", value=orb_obliqr, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          errorCode = rc
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

          call ESMF_StateGet(export_state, itemName="d2x", array=o2x, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          call ESMF_ArrayGet(o2x, localDe=0, farrayPtr=fptr, rc=rc)
          if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          call ocn_export(fptr, ldiag_cpl, rc)

          errorCode = rc
          if (errorCode /= POP_Success) then
             call POP_ErrorPrint(errorCode)
             call exit_POP(sigAbort, 'ERROR in ocn_export')
          endif

          exit advance
       end if
       
    enddo advance

    if ( lsend_precip_fact ) then
       call ESMF_AttributeSet(export_state, name="precip_fact", value=precip_fact, rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
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

  end subroutine ocn_run_esmf

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_esmf
!
! !INTERFACE:
subroutine ocn_final_esmf(comp, import_state, export_state, Eclock, rc)

!
! !DESCRIPTION:
! Finalize POP
!
! !USES:
    use POP_FinalMod
!
! !ARGUMENTS:
!
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)                        :: rc
!
! !LOCAL VARIABLES:
!
    type(ESMF_Array)                 :: d2x, x2d
    type(ESMF_DistGrid)              :: distgrid_ref

! Author: Fei Liu
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer (POP_i4) ::  errorCode         ! error code

!-----------------------------------------------------------------------

    call POP_Final(errorCode)
    rc = ESMF_SUCCESS

    ! Destroy ESMF objects
    ! Destroy ESMF objects

    call esmfshr_util_StateArrayDestroy(export_state,"d2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(import_state,"x2d",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine ocn_final_esmf

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_esmf
! !INTERFACE:

  type(ESMF_DistGrid) function ocn_DistGrid_esmf(gsize)

! !DESCRIPTION:
!  This routine creates the ocean distgrid
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer, intent(out)            :: gsize

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
      i,j, n, iblock, &
      lsize,   &
      rc,      &
      ier

    type (block) ::       &
      this_block          ! block information for current block

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
    allocate(gindex(lsize),stat=ier)

    ! not correct for padding, use "n" above
    !    lsize = block_size_x*block_size_y*nblocks_clinic
    gsize = nx_global*ny_global

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

    ocn_distgrid_esmf = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end function ocn_DistGrid_esmf

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_esmf
! !INTERFACE:

 subroutine ocn_domain_esmf( dom )

! !DESCRIPTION:
!  This routine creates the ocean domain
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    type(ESMF_Array), intent(inout)     :: dom

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) :: rc

    integer (int_kind) ::   &
      i,j, n, iblock

    integer (int_kind) ::   &
      klon,klat,karea,kmask,kfrac ! domain fields

    type (block) ::       &
      this_block          ! block information for current block

    real(R8),    pointer ::  &
      fptr (:,:)          ! data pointer into ESMF array

    real(R8)  :: &
      frac                ! temporary var to compute frac/mask from KMT

!-----------------------------------------------------------------------

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

!-------------------------------------------------------------------
!
!  initialize domain type, lat/lon in degrees,
!  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
!  Fill in correct values for domain components
!
!-------------------------------------------------------------------

    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    fptr(:,:) = -9999.0_R8
    n=0

    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          fptr(klon , n)          = TLON(i,j,iblock)*radian
          fptr(klat , n)          = TLAT(i,j,iblock)*radian 
          fptr(karea, n)          = TAREA(i,j,iblock)/(radius*radius)
          frac                    = float(KMT(i,j,iblock)) 
          if (frac > 1.0_r8) frac = 1.0_r8
          fptr(kfrac, n)          = frac
          fptr(kmask, n)          = frac
       enddo
       enddo
    enddo

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_esmf

#endif

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
end module ocn_comp_esmf
