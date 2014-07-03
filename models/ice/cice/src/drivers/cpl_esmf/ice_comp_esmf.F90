module ice_comp_esmf

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: ice_comp_esmf
!
! !DESCRIPTION:
! CICE interface routine for the ccsm cpl7 esmf system
!
! !USES:

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod,  only : shr_sys_abort, shr_sys_flush
  use shr_file_mod, only : shr_file_getlogunit, shr_file_getloglevel,  &
		           shr_file_setloglevel, shr_file_setlogunit
  use esmf

  use seq_flds_mod
  use seq_cdata_mod,   only : seq_cdata, seq_cdata_setptrs
  use seq_timemgr_mod, only : seq_timemgr_eclockgetdata, &
                              seq_timemgr_restartalarmison, &
		              seq_timemgr_eclockdateinsync, &
                              seq_timemgr_stopalarmison
  use seq_infodata_mod,only : seq_infodata_start_type_cont, &
                              seq_infodata_start_type_brnch, seq_infodata_start_type_start
  use seq_comm_mct,    only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use perf_mod,        only : t_startf, t_stopf

  use ice_cpl_indices
  use ice_flux,        only : strairxt, strairyt, strocnxt, strocnyt,    &
    	                      alvdr, alidr, alvdf, alidf, Tref, Qref, Uref, &
                              flat, fsens, flwout, evap, fswabs, fhocn, &
                              fswthru,     &
		              fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa,  &
		              rhoa, swvdr, swvdf, swidr, swidf, flw, frain,    &
		              fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt,&
		              sss, tf, wind, fsw, init_flux_atm, init_flux_ocn,&
                              faero
  use ice_state,       only : vice, vsno, aice, trcr, filename_aero, filename_iage, &
                              filename_volpn, filename_FY, filename_lvl, &
                              tr_aero, tr_iage, tr_FY, tr_pond, tr_lvl
  use ice_domain_size, only : nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_domain,      only : nblocks, blocks_ice, halo_info, distrb_info, profile_barrier
  use ice_blocks,      only : block, get_block, nx_block, ny_block
  use ice_grid,        only : tlon, tlat, tarea, tmask, anglet, hm, ocn_gridcell_frac, &
                              grid_type, t2ugrid_vector
  use ice_constants,   only : c0, c1, puny, tffresh, spval_dbl, rad_to_deg, radius, &
                              field_loc_center, field_type_scalar, field_type_vector, c100
  use ice_communicate, only : my_task, master_task, lprint_stats, MPI_COMM_ICE
  use ice_calendar,    only : idate, mday, time, month, daycal, secday, &
		              sec, dt, dt_dyn, xndt_dyn, calendar,      &
                              calendar_type, nextsw_cday, days_per_year,&
                              get_daycal, leap_year_count, nyr, new_year
  use ice_orbital,     only : eccen, obliqr, lambm0, mvelpp
  use ice_timers
  use ice_probability, only : init_numIceCells, print_numIceCells,  &
                              write_numIceCells, accum_numIceCells2

  use ice_kinds_mod,   only : int_kind, dbl_kind, char_len_long, log_kind
  use ice_boundary,    only : ice_HaloUpdate 
  use ice_scam,        only : scmlat, scmlon, single_column
  use ice_fileunits,   only : nu_diag, inst_index, inst_name, inst_suffix
  use ice_dyn_evp,     only : kdyn
  use ice_prescribed_mod
  use ice_step_mod
  use ice_global_reductions
  use ice_broadcast
  use CICE_RunMod

  use esmfshr_mod

! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ice_register_esmf
  public :: ice_init_esmf
  public :: ice_run_esmf
  public :: ice_final_esmf
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein, Fei Liu
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: ice_export_esmf
  private :: ice_import_esmf
  private :: ice_DistGrid_esmf
  private :: ice_domain_esmf

!
! !PRIVATE VARIABLES

  type(mct_gsmap),save :: gsmap_i
  type(mct_ggrid),save :: dom_i
 
  logical (kind=log_kind),save :: atm_aero

!=======================================================================

contains

!=======================================================================
subroutine ice_register_esmf(comp, rc)
    implicit none
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    print *, "In ice register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      ice_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      ice_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      ice_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

!=======================================================================
!BOP
!
! !IROUTINE: ice_init_esmf
!
! !INTERFACE:
  subroutine ice_init_esmf(comp, import_state, export_state, EClock, rc)
   !----------------------------------------------------------
!
! !DESCRIPTION:
! Initialize thermodynamic ice model and obtain relevant atmospheric model
! arrays back from driver 
!
! !USES:

    use CICE_InitMod
    use ice_restart, only: runid, runtype, restart_dir, restart_format
    use ice_history, only: history_dir, history_file
    implicit none
!
! !ARGUMENTS:
!
    !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
   
!
! !LOCAL VARIABLES:
!
    !----- local -----
    type(ESMF_DistGrid)                   :: distgrid
    type(ESMF_Array)                      :: d2x, x2d, dom
    type(ESMF_VM)                         :: vm

    character(len=32)  :: starttype          ! infodata start type
    integer            :: start_ymd          ! Start date (YYYYMMDD)
    integer            :: start_tod          ! start time of day (s)
    integer            :: curr_ymd           ! Current date (YYYYMMDD)
    integer            :: curr_tod           ! Current time of day (s)
    integer            :: ref_ymd            ! Reference date (YYYYMMDD)
    integer            :: ref_tod            ! reference time of day (s)
    integer            :: iyear              ! yyyy
    integer            :: nyrp               ! yyyy
    integer            :: dtime              ! time step
    integer            :: shrlogunit,shrloglev ! old values
    integer            :: iam,ierr
    integer            :: lbnum
    integer            :: daycal(13)  !number of cumulative days per month
    integer            :: nleaps      ! number of leap days before current year
    integer            :: mpicom_loc, mpicom_vm, gsize
    integer            :: ICEID       ! cesm ID value

    character(ESMF_MAXSTR) :: convCIM, purpComp

! !REVISION HISTORY:
! Author: Fei Liu
!EOP
!-----------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! Determine attribute vector indices
   !---------------------------------------------------------------------------
 
   call ice_cpl_indices_set()

   rc = ESMF_SUCCESS

   ! duplicate the mpi communicator from the current VM 
   call ESMF_VMGetCurrent(vm, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call MPI_Comm_dup(mpicom_vm, mpicom_loc, rc)
   if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Initialize cice id
   
   call ESMF_AttributeGet(export_state, name="ID", value=ICEID, rc=rc)
   if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Determine time of next atmospheric shortwave calculation
   call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Determine if aerosols are coming through the coupler
   call ESMF_AttributeGet(export_state, name="atm_aero", value=atm_aero, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   call ESMF_AttributeGet(export_state, name="ID", value=ICEID, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Determine orbital parameters
   call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !---------------------------------------------------------------------------
    ! use infodata to determine type of run
    !---------------------------------------------------------------------------

    ! Preset single column values

    single_column = .false.
    scmlat = -999.
    scmlon = -999.

    call ESMF_AttributeGet(export_state, name="case_name", value=runid, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="single_column", value=single_column, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="scmlat", value=scmlat, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="scmlon", value=scmlon, rc=rc)
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
       write(nu_diag,*) 'ice_comp_esmf ERROR: unknown starttype'
       call shr_sys_abort()
    end if

    ! Set nextsw_cday to -1 for continue and branch runs.
    if (trim(runtype) /= 'initial') nextsw_cday = -1

    !=============================================================
    ! Set ice dtime to ice coupling frequency
    !=============================================================

    call seq_timemgr_EClockGetData(EClock, dtime=dtime, calendar=calendar_type)
    dt = real(dtime)

    !=============================================================
    ! Initialize cice because grid information is needed for
    ! creation of GSMap_ice.  cice_init also sets time manager info
    !=============================================================

    inst_name   = seq_comm_name(ICEID)
    inst_index  = seq_comm_inst(ICEID)
    inst_suffix = seq_comm_suffix(ICEID)

    call t_startf ('cice_init')
    call cice_init( mpicom_loc )
    call t_stopf ('cice_init')
    call init_numIceCells

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)
   
    !----------------------------------------------------------------------------
    ! use EClock to reset calendar information on initial start
    !----------------------------------------------------------------------------

    ! - the following logic duplicates the logic for the concurrent system - 
    ! cice_init is called then init_cpl is called where the start date is received
    ! from the flux coupler
    ! - in the following calculation for the variable time, iyear is used
    ! rather than (iyear-1) (as in init_cpl) since the sequential system permits you 
    ! to start with year "0" not year "1"
    ! - on restart run 
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)
    ! - on initial run 
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from iyear, month and mday
    !   - istep0 and istep1 are set to 0 

    call seq_timemgr_EClockGetData(EClock, &
         start_ymd=start_ymd, start_tod=start_tod,       &
         curr_ymd=curr_ymd,   curr_tod=curr_tod,         &
         ref_ymd=ref_ymd, ref_tod=ref_tod)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) 'ice_comp_esmf: ref_ymd ',ref_ymd, &
                  ' must equal start_ymd ',start_ymd
             write(nu_diag,*) 'ice_comp_esmf: ref_ymd ',ref_tod, &
                  ' must equal start_ymd ',start_tod
          end if
       end if

       if (my_task == master_task) then
          write(nu_diag,*) '(ice_init_esmf) idate from sync clock = ', &
               start_ymd
          write(nu_diag,*) '(ice_init_esmf)   tod from sync clock = ', &
               start_tod
          write(nu_diag,*) &
               '(ice_init_esmf) resetting idate to match sync clock'
       end if

       idate = curr_ymd
       iyear = (idate/10000)                     ! integer year of basedate
       if (calendar_type .eq. "GREGORIAN") then 	
          nyr = iyear+1
          new_year = .false.
       end if
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100-1    ! day of month of basedate
                                                 ! (starts at 0)

       call get_daycal(year=iyear,days_per_year_in=days_per_year, &
                       daycal_out=daycal)

       nleaps = leap_year_count(iyear)   ! this sets nleaps in ice_calendar
       time  = (((iyear)*days_per_year  + nleaps + daycal(month)+mday)*secday) &
             + start_tod    

       call shr_sys_flush(nu_diag)
    end if
    call calendar(time)     ! update calendar info
 
    !---------------------------------------------------------------------------
    ! Initialize distgrids, domains, and arrays
    !---------------------------------------------------------------------------

    call t_startf ('cice_esmf_init')

    ! Initialize ice distgrid

    distgrid = ice_DistGrid_esmf(gsize,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize ice domain (needs ice initialization info)

    dom = mct2esmf_init(distgrid, attname=seq_flds_dom_fields, name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ice_domain_esmf(dom, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Inialize input/output arrays
    d2x = mct2esmf_init(distgrid, attname=seq_flds_i2x_fields, name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    x2d = mct2esmf_init(distgrid, attname=seq_flds_x2i_fields, name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/d2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2d/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(distgrid,ICEID,gsmap_i,MPI_COMM_ICE,gsize=gsize,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom,dom_i,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmf2mct_copy(dom, dom_i%data, rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------------------------------
    ! Second phase of prescribed ice initialization
    !-----------------------------------------------------------------

    call ice_prescribed_init(ICEID, gsmap_i, dom_i)

    !-----------------------------------------------------------------
    ! get ready for coupling
    !-----------------------------------------------------------------

    call coupling_prep

    !---------------------------------------------------------------------------
    ! send initial state to driver
    !---------------------------------------------------------------------------

    call ice_export_esmf (d2x, rc=rc)  !Send initial state to driver
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="ice_prognostic", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ice_nx", value=nx_global, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ice_ny", value=ny_global, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "CICE", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Community Ice CodE", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "CICE4 is the latest version of the Los Alamos Sea Ice " // &
                  "Model, sometimes referred to as the Community Ice " // &
                  "CodE.  It is the result of a community effort to  " // &
                  "develop a portable, efficient sea ice model that can " // &
                  "be run coupled in a global climate model or uncoupled " // &
                  "as a stand-alone ice model.  It has been released as " // &
                  "the sea ice component of the Community Earth System " // &
                  "Model (CESM), a fully-coupled global climate model " // &
                  "that provides simulations of the earths past, present " // &
                  "and future climate states.  CICE4 is supported on " // &
                  "high- and low-resolution Greenland Pole and tripole " // &
                  "grids, which are identical to those used by the " // &
                  "Parallel Ocean Program (POP) ocean model.  The high " // &
                  "resolution version is best suited for simulating " // &
                  "present-day and future climate scenarios while the low " // &
                  "resolution option is used for paleoclimate simulations " // &
                  "and debugging.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Ice", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "someone", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "someone@someplace", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)

    call t_stopf ('cice_esmf_init')

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    call ice_timer_stop(timer_total) ! time entire run

end subroutine ice_init_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_run_esmf
!
! !INTERFACE:
subroutine ice_run_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Run thermodynamic CICE
!
! !USES:
    use ice_history
    use ice_restart
    use ice_diagnostics
    use ice_aerosol, only: write_restart_aero
    use ice_age, only: write_restart_age
    use ice_meltpond, only: write_restart_pond
    use ice_FY, only: write_restart_FY
    use ice_lvl, only: write_restart_lvl
    use ice_restoring, only: restore_ice, ice_HaloRestore
    use ice_shortwave, only: init_shortwave
    implicit none

! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

! !LOCAL VARIABLES:
    integer :: k             ! index
    logical :: rstwr         ! .true. ==> write a restart file
    logical :: stop_now      ! .true. ==> stop at the end of this run phase
    integer :: ymd           ! Current date (YYYYMMDD)
    integer :: tod           ! Current time of day (sec)
    integer :: curr_ymd      ! Current date (YYYYMMDD)
    integer :: curr_tod      ! Current time of day (s)
    integer :: yr_sync       ! Sync current year
    integer :: mon_sync      ! Sync current month
    integer :: day_sync      ! Sync current day
    integer :: tod_sync      ! Sync current time of day (sec)
    integer :: ymd_sync      ! Current year of sync clock
    integer :: shrlogunit,shrloglev ! old values
    integer :: lbnum
    integer :: n, nyrp
    character(len=char_len_long) :: string1, string2
    character(len=char_len_long) :: fname
    character(len=*), parameter  :: SubName = "ice_run_esmf"
    type(ESMF_Array) :: d2x, x2d
    logical, save :: first_time = .true.
!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein, Fei Liu
!
!EOP
!---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ice_timer_start(timer_total) ! time entire run

    !---------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !---------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)
   
    ! Determine time of next atmospheric shortwave calculation

    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

   ! Determine orbital parameters
   call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=curr_ymd,   curr_tod=curr_tod)

    if (calendar_type .eq. "GREGORIAN") then 	
       nyrp = nyr
       nyr = (curr_ymd/10000)+1           ! integer year of basedate
       if (nyr /= nyrp) then
          new_year = .true.
       else
          new_year = .false.
       end if
    end if

    !-------------------------------------------------------------------
    ! get import state
    !-------------------------------------------------------------------
    
    call t_startf ('cice_import')
    call ice_timer_start(timer_cplrecv)

    call ESMF_StateGet(import_state, itemName="x2d", array=x2d, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ice_import_esmf(x2d, rc=rc )
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_import')
 
    !--------------------------------------------------------------------
    ! timestep update
    !--------------------------------------------------------------------

    call ice_timer_start(timer_step)

    istep  = istep  + 1    ! update time step counters
    istep1 = istep1 + 1
    time = time + dt       ! determine the time and date
    call calendar(time)    ! at the end of the timestep

    if (resttype == 'old' .and. istep == 1) &
       call init_shortwave    ! initialize radiative transfer using current swdn

    !-----------------------------------------------------------------
    ! restoring on grid boundaries
    !-----------------------------------------------------------------

    if (restore_ice) call ice_HaloRestore

    call t_startf ('cice_initmd')
    call init_mass_diags   ! diagnostics per timestep
    call t_stopf ('cice_initmd')

    if(prescribed_ice) then  ! read prescribed ice
       call t_startf ('cice_presc')
       call ice_prescribed_run(idate, sec)
       call t_stopf ('cice_presc')
    endif
    
    call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
    call init_flux_ocn        ! initialize ocean fluxes sent to coupler

    !-----------------------------------------------------------------
    ! Scale radiation fields
    !-----------------------------------------------------------------

    call t_startf ('cice_prep_radiation')
    call ice_timer_start(timer_sw)
    call prep_radiation(dt)
    call ice_timer_stop(timer_sw)
    call t_stopf ('cice_prep_radiation')
    
    !-----------------------------------------------------------------
    ! thermodynamics1
    !-----------------------------------------------------------------

    call t_startf ('cice_therm1')
    call step_therm1(dt)
    call t_stopf ('cice_therm1')
    
    !-----------------------------------------------------------------
    ! thermodynamics2
    !-----------------------------------------------------------------

    if (.not.prescribed_ice) then
       call t_startf ('cice_therm2')
       call step_therm2 (dt)  ! post-coupler thermodynamics
       call t_stopf ('cice_therm2')
    end if

   !-----------------------------------------------------------------
   ! dynamics, transport, ridging
   !-----------------------------------------------------------------

    if (.not.prescribed_ice .and. kdyn>0) then
       if (xndt_dyn > c1) then
          call t_startf ('cice_dyn')
          do k = 1, nint(xndt_dyn)
             call step_dynamics(dt_dyn,dt) ! dynamics, transport, ridging
          enddo
          call t_stopf ('cice_dyn')
       else
          if (mod(time, dt_dyn) == c0) then
             call t_startf ('cice_dyn')
             call step_dynamics(dt_dyn,dt) ! dynamics, transport, ridging
             call t_stopf ('cice_dyn')
          endif
       endif
    endif ! not prescribed_ice
    
    !-----------------------------------------------------------------
    ! radiation
    !-----------------------------------------------------------------

    call t_startf ('cice_radiation')
    call ice_timer_start(timer_sw)
    call step_radiation(dt)
    call ice_timer_stop(timer_sw)
    call t_stopf ('cice_radiation')
    
    !-----------------------------------------------------------------
    ! get ready for coupling
    !-----------------------------------------------------------------

    call coupling_prep

    call ice_timer_stop(timer_step)

    !-----------------------------------------------------------------
    ! write data
    !-----------------------------------------------------------------
    
    call t_startf ('cice_diag')
    call ice_timer_start(timer_diags)
    if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
    call ice_timer_stop(timer_diags)
    call t_stopf ('cice_diag')
    
    call t_startf ('cice_hist')
    call ice_timer_start(timer_hist)
#if (defined _NOIO)
!  Not enought memory on BGL to write a history file yet! 
!    call ice_write_hist (dt)    ! history file
#else
    call ice_write_hist (dt)    ! history file
#endif
    call ice_timer_stop(timer_hist)
    call t_stopf ('cice_hist')

    !--------------------------------------------
    ! Accumualate the number of active ice cells
    !--------------------------------------------
    call accum_numIceCells2(aice)

    rstwr = seq_timemgr_RestartAlarmIsOn(EClock)
    if (rstwr) then
       call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd_sync, curr_tod=tod_sync, &
          curr_yr=yr_sync,curr_mon=mon_sync,curr_day=day_sync)
       fname = restart_filename(yr_sync, mon_sync, day_sync, tod_sync)

       if (my_task == master_task) then
          write(nu_diag,*) &
            'ice_comp_mct: calling dumpfile for restart filename= ', fname
       endif

       if (restart_format /= 'nc') then

          if (tr_pond) then
               n = index(fname,'cice.r') + 6
               string1 = trim(fname(1:n-1))
               string2 = trim(fname(n:lenstr(fname)))
               write(filename_volpn,'(a,a,a,a)') &
                  string1(1:lenstr(string1)),'.volpn', &
                  string2(1:lenstr(string2))
          endif
          if (tr_aero) then
               n = index(fname,'cice.r') + 6
               string1 = trim(fname(1:n-1))
               string2 = trim(fname(n:lenstr(fname)))
               write(filename_aero,'(a,a,a,a)') &
                  string1(1:lenstr(string1)),'.aero', &
                  string2(1:lenstr(string2))
          endif
          if (tr_iage) then
               n = index(fname,'cice.r') + 6
               string1 = trim(fname(1:n-1))
               string2 = trim(fname(n:lenstr(fname)))
               write(filename_iage,'(a,a,a,a)') &
                  string1(1:lenstr(string1)),'.age', &
                  string2(1:lenstr(string2))
          endif
          if (tr_FY) then
               n = index(fname,'cice.r') + 6
               string1 = trim(fname(1:n-1))
               string2 = trim(fname(n:lenstr(fname)))
               write(filename_FY,'(a,a,a,a)') &
                  string1(1:lenstr(string1)),'.FY', &
                  string2(1:lenstr(string2))
          endif
          if (tr_lvl) then
               n = index(fname,'cice.r') + 6
               string1 = trim(fname(1:n-1))
               string2 = trim(fname(n:lenstr(fname)))
               write(filename_lvl,'(a,a,a,a)') &
                  string1(1:lenstr(string1)),'.lvl', &
                  string2(1:lenstr(string2))
          endif

       endif ! restart_format

#if (defined _NOIO)
!  Not enought memory on BGL to call dumpfile  file yet! 
!       call dumpfile(fname)
#else
       call ice_timer_start(timer_readwrite)
       call dumpfile(fname)
       if (restart_format /= 'nc') then
          if (tr_aero) call write_restart_aero(filename_aero)
          if (tr_iage) call write_restart_age(filename_iage)
          if (tr_FY)   call write_restart_FY(filename_FY)
          if (tr_lvl)  call write_restart_lvl(filename_lvl)
          if (tr_pond) call write_restart_pond(filename_volpn)
       endif
       call ice_timer_stop(timer_readwrite)
#endif
    end if

    !-----------------------------------------------------------------
    ! send export state to driver 
    !-----------------------------------------------------------------
    
    call t_startf ('cice_export')
    call ice_timer_start(timer_cplsend)

    call ESMF_StateGet(export_state, itemName="d2x", array=d2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ice_export_esmf (d2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_export')
    
    !--------------------------------------------------------------------
    ! check that internal clock is in sync with master clock
    !--------------------------------------------------------------------

    tod = sec
    ymd = idate
    if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod )) then
       call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync )
       write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
       write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       call shr_sys_abort( SubName// &
          ":: Internal sea-ice clock not in sync with Sync Clock")
    end if
   
    ! reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    !-------------------------------------------------------------------
    ! stop timers and print timer info
    !-------------------------------------------------------------------
    ! Need to have this logic here instead of in ice_final_mct since 
    ! the ice_final_mct.F90 will still be called even in aqua-planet mode
    ! Could put this logic in the driver - but it seems easier here 

    ! Need to stop this at the end of every run phase in a coupled run.
    call ice_timer_stop(timer_total)        ! stop timing

    stop_now = seq_timemgr_StopAlarmIsOn( EClock )
    if (stop_now) then
       call ice_timer_print_all(stats=.true.) ! print timing information
       if(lprint_stats) then 
          call write_numIceCells()
       endif
       call release_all_fileunits
    end if
    
end subroutine ice_run_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_final_esmf
!
! !INTERFACE:
subroutine ice_final_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Finalize CICE
!
! !USES:
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
   implicit none
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
!
! !REVISION HISTORY:
! Author: Fei Liu
!
!EOP
!---------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Finalize routine 
  !----------------------------------------------------------------------------
  
  ! Note that restart for final timestep was written in run phase.
    rc = ESMF_SUCCESS

    ! Destroy ESMF objects

    call esmfshr_util_StateArrayDestroy(export_state,"d2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(import_state,"x2d",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine ice_final_esmf

!=================================================================================
function ice_DistGrid_esmf(gsize,rc)

    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(out)    :: gsize
    integer, intent(out)    :: rc
    ! Return
    type(ESMF_DistGrid)     :: ice_DistGrid_esmf
    !
    ! Local variables
    !
    integer,allocatable :: gindex(:)
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize
    integer     :: ier
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    !-------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Build the CICE grid numbering for distgrid
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP

    ! number the local grid

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
          enddo !i
       enddo    !j
    enddo        !iblk
    lsize = n

! not valid for padded decomps
!    lsize = block_size_x*block_size_y*nblocks
    gsize = nx_global*ny_global

    allocate(gindex(lsize),stat=ier)
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon = this_block%i_glob(i)
             lat = this_block%j_glob(j)
             gi = (lat-1)*nx_global + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo        !iblk
   
    ice_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function ice_DistGrid_esmf


!====================================================================================

subroutine ice_export_esmf(array, rc )   

    implicit none
    !-----------------------------------------------------
    type(ESMF_Array), intent(inout)     :: array
    integer, intent(out)                :: rc

    integer :: i, j, iblk, n
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain

    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
        Tsrf  &      ! surface temperature
     ,  tauxa &      ! atmo/ice stress
     ,  tauya &
     ,  tauxo &      ! ice/ocean stress
     ,  tauyo &
     ,  ailohi       ! fractional ice area

    real (kind=dbl_kind) :: workx, worky           ! tmps for converting grid

    type(block)        :: this_block               ! block information for current block
    logical :: flag
    real(R8), pointer   :: fptr(:,:)
    !-----------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    flag=.false.

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
             
          ! ice fraction
          ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

          ! surface temperature
          Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)             !Kelvin (original ???)
          
          ! wind stress  (on POP T-grid:  convert to lat-lon)
          workx = strairxT(i,j,iblk)                             ! N/m^2
          worky = strairyT(i,j,iblk)                             ! N/m^2
          tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                          - worky*sin(ANGLET(i,j,iblk))
          tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                          + workx*sin(ANGLET(i,j,iblk))
          
          ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
          workx = -strocnxT(i,j,iblk)                            ! N/m^2
          worky = -strocnyT(i,j,iblk)                            ! N/m^2
          tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                          - worky*sin(ANGLET(i,j,iblk))
          tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                          + workx*sin(ANGLET(i,j,iblk))

       enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
          if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
             flag = .true.
          endif
       end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          do j = 1, ny_block
          do i = 1, nx_block
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                write(nu_diag,*) &
                     ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                call shr_sys_flush(nu_diag)
             endif
          end do
          end do
       end do
    endif

    ! Fill export state i2x_i

     fptr(:,:) = spval_dbl

     n=0
     do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            n = n+1

            !-------states-------------------- 
            fptr(index_i2x_Si_ifrac ,n)    = ailohi(i,j,iblk)   

            if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
               !-------states-------------------- 
               fptr(index_i2x_Si_t     ,n)    = Tsrf(i,j,iblk)
               fptr(index_i2x_Si_avsdr ,n)    = alvdr(i,j,iblk)
               fptr(index_i2x_Si_anidr ,n)    = alidr(i,j,iblk)
               fptr(index_i2x_Si_avsdf ,n)    = alvdf(i,j,iblk)
               fptr(index_i2x_Si_anidf ,n)    = alidf(i,j,iblk)
               fptr(index_i2x_Si_tref  ,n)    = Tref(i,j,iblk)
               fptr(index_i2x_Si_qref  ,n)    = Qref(i,j,iblk)
               fptr(index_i2x_Si_u10  ,n)     = Uref(i,j,iblk)
               fptr(index_i2x_Si_snowh ,n)    = vsno(i,j,iblk) &
                                              / ailohi(i,j,iblk)
            
               !--- a/i fluxes computed by ice
               fptr(index_i2x_Faii_taux ,n)   = tauxa(i,j,iblk)    
               fptr(index_i2x_Faii_tauy ,n)   = tauya(i,j,iblk)    
               fptr(index_i2x_Faii_lat  ,n)   = flat(i,j,iblk)     
               fptr(index_i2x_Faii_sen  ,n)   = fsens(i,j,iblk)    
               fptr(index_i2x_Faii_lwup ,n)   = flwout(i,j,iblk)   
               fptr(index_i2x_Faii_evap ,n)   = evap(i,j,iblk)     
               fptr(index_i2x_Faii_swnet,n)   = fswabs(i,j,iblk)
            
               !--- i/o fluxes computed by ice
               fptr(index_i2x_Fioi_melth,n)   = fhocn(i,j,iblk)
               fptr(index_i2x_Fioi_swpen,n)   = fswthru(i,j,iblk) ! hf from melting          
               fptr(index_i2x_Fioi_meltw,n)   = fresh(i,j,iblk)   ! h2o flux from melting    ???
               fptr(index_i2x_Fioi_salt ,n)   = fsalt(i,j,iblk)   ! salt flux from melting   ???
               fptr(index_i2x_Fioi_taux ,n)   = tauxo(i,j,iblk)   ! stress : i/o zonal       ???
               fptr(index_i2x_Fioi_tauy ,n)   = tauyo(i,j,iblk)   ! stress : i/o meridional  ???
            end if
         enddo    !i
         enddo    !j
     enddo        !iblk

end subroutine ice_export_esmf

!====================================================================================

subroutine ice_import_esmf(array, rc)

    !-----------------------------------------------------

    implicit none
    type(ESMF_Array), intent(inout) :: array
    integer, intent(out)            :: rc

    integer :: i, j, iblk, n
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block) :: this_block      ! block information for current block
    real(R8), pointer   :: fptr(:,:)

    integer,parameter :: nflds=15,nfldv=6
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)

    real (kind=dbl_kind) :: &
         workx, worky
    logical (kind=log_kind) :: &
         first_call = .true.

    !-----------------------------------------------------

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    !-----------------------------------------------------------------
    ! Zero stuff while waiting, only filling in active cells.
    !-----------------------------------------------------------------

    rc = ESMF_SUCCESS
    
    call ESMF_ArrayGet(array, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Use aflds to gather the halo updates of multiple fields
    ! Need to separate the scalar from the vector halo updates

    allocate(aflds(nx_block,ny_block,nflds,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
       do i = ilo, ihi

          n = n+1
          aflds(i,j, 1,iblk)   = fptr(index_x2i_So_t,n)
          aflds(i,j, 2,iblk)   = fptr(index_x2i_So_s,n)
          aflds(i,j, 3,iblk)   = fptr(index_x2i_Sa_z,n)
          aflds(i,j, 4,iblk)   = fptr(index_x2i_Sa_ptem,n)
          aflds(i,j, 5,iblk)   = fptr(index_x2i_Sa_tbot,n)
          aflds(i,j, 6,iblk)   = fptr(index_x2i_Sa_shum,n)
          aflds(i,j, 7,iblk)   = fptr(index_x2i_Sa_dens,n)
          aflds(i,j, 8,iblk)   = fptr(index_x2i_Fioo_q,n)
          aflds(i,j, 9,iblk)   = fptr(index_x2i_Faxa_swvdr,n)
          aflds(i,j,10,iblk)   = fptr(index_x2i_Faxa_swndr,n)
          aflds(i,j,11,iblk)   = fptr(index_x2i_Faxa_swvdf,n)
          aflds(i,j,12,iblk)   = fptr(index_x2i_Faxa_swndf,n)
          aflds(i,j,13,iblk)   = fptr(index_x2i_Faxa_lwdn,n)
          aflds(i,j,14,iblk)   = fptr(index_x2i_Faxa_rain,n)
          aflds(i,j,15,iblk)   = fptr(index_x2i_Faxa_snow,n)

         enddo    !i
         enddo    !j

     enddo        !iblk

     if (.not.prescribed_ice) then
        call t_startf ('cice_imp_halo')
        call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
                                              field_type_scalar)
        call t_stopf ('cice_imp_halo')
     endif
 
     !$OMP PARALLEL DO PRIVATE(iblk,i,j)
     do iblk = 1, nblocks
        do j = 1,ny_block
        do i = 1,nx_block
           sst  (i,j,iblk)   = aflds(i,j, 1,iblk)
           sss  (i,j,iblk)   = aflds(i,j, 2,iblk)
           zlvl (i,j,iblk)   = aflds(i,j, 3,iblk)
           potT (i,j,iblk)   = aflds(i,j, 4,iblk)
           Tair (i,j,iblk)   = aflds(i,j, 5,iblk)
           Qa   (i,j,iblk)   = aflds(i,j, 6,iblk)
           rhoa (i,j,iblk)   = aflds(i,j, 7,iblk)
           frzmlt (i,j,iblk) = aflds(i,j, 8,iblk)
           swvdr(i,j,iblk)   = aflds(i,j, 9,iblk)
           swidr(i,j,iblk)   = aflds(i,j,10,iblk)
           swvdf(i,j,iblk)   = aflds(i,j,11,iblk)
           swidf(i,j,iblk)   = aflds(i,j,12,iblk)
           flw  (i,j,iblk)   = aflds(i,j,13,iblk)
           frain(i,j,iblk)   = aflds(i,j,14,iblk)
           fsnow(i,j,iblk)   = aflds(i,j,15,iblk)
        enddo    !i
        enddo    !j
     enddo        !iblk
     !$OMP END PARALLEL DO
 
     deallocate(aflds)
     allocate(aflds(nx_block,ny_block,nfldv,nblocks))
     aflds = c0

     n=0
     do iblk = 1, nblocks
        this_block = get_block(blocks_ice(iblk),iblk)
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi
 
        do j = jlo, jhi
        do i = ilo, ihi
           n = n+1
           aflds(i,j, 1,iblk)   = fptr(index_x2i_So_u,n)
           aflds(i,j, 2,iblk)   = fptr(index_x2i_So_v,n)
           aflds(i,j, 3,iblk)   = fptr(index_x2i_Sa_u,n)
           aflds(i,j, 4,iblk)   = fptr(index_x2i_Sa_v,n)
           aflds(i,j, 5,iblk)   = fptr(index_x2i_So_dhdx,n)
           aflds(i,j, 6,iblk)   = fptr(index_x2i_So_dhdy,n)
        enddo
        enddo
     enddo

     if (.not.prescribed_ice) then
        call t_startf ('cice_imp_halo')
        call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
                                             field_type_vector)
        call t_stopf ('cice_imp_halo')
     endif

     !$OMP PARALLEL DO PRIVATE(iblk,i,j)
     do iblk = 1, nblocks
        do j = 1,ny_block
        do i = 1,nx_block
           uocn (i,j,iblk)   = aflds(i,j, 1,iblk)
           vocn (i,j,iblk)   = aflds(i,j, 2,iblk)
           uatm (i,j,iblk)   = aflds(i,j, 3,iblk)
           vatm (i,j,iblk)   = aflds(i,j, 4,iblk)
           ss_tltx(i,j,iblk) = aflds(i,j, 5,iblk)
           ss_tlty(i,j,iblk) = aflds(i,j, 6,iblk)
        enddo    !i
        enddo    !j
     enddo        !iblk
     !$OMP END PARALLEL DO

     deallocate(aflds)

     !-------------------------------------------------------
     ! Set aerosols from coupler 
     !-------------------------------------------------------
     
      if (first_call) then
         if (tr_aero .and. .not. atm_aero) then
            write(nu_diag,*) 'ice_comp_mct ERROR: atm_aero must be set for tr_aero' 
            call shr_sys_abort()
         end if
        first_call = .false.
      end if

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            n = n+1
            faero(i,j,1,iblk) = fptr(index_x2i_Faxa_bcphodry,n)

            faero(i,j,2,iblk) = fptr(index_x2i_Faxa_bcphidry,n) &
                              + fptr(index_x2i_Faxa_bcphiwet,n)
         ! Combine all of the dust into one category
            faero(i,j,3,iblk) = fptr(index_x2i_Faxa_dstwet1,n) &
                              + fptr(index_x2i_Faxa_dstdry1,n) &
                              + fptr(index_x2i_Faxa_dstwet2,n) &
                              + fptr(index_x2i_Faxa_dstdry2,n) &
                              + fptr(index_x2i_Faxa_dstwet3,n) &
                              + fptr(index_x2i_Faxa_dstdry3,n) &
                              + fptr(index_x2i_Faxa_dstwet4,n) &
                              + fptr(index_x2i_Faxa_dstdry4,n)

         enddo    !i
         enddo    !j

      enddo        !iblk

     !-----------------------------------------------------------------
     ! rotate zonal/meridional vectors to local coordinates
     ! compute data derived quantities
     !-----------------------------------------------------------------
     
     ! Vector fields come in on T grid, but are oriented geographically
     ! need to rotate to pop-grid FIRST using ANGLET
     ! then interpolate to the U-cell centers  (otherwise we
     ! interpolate across the pole)
     ! use ANGLET which is on the T grid !
     
     call t_startf ('cice_imp_ocn')
     !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
     do iblk = 1, nblocks

       do j = 1,ny_block
       do i = 1,nx_block

          ! ocean
          workx      = uocn  (i,j,iblk) ! currents, m/s 
          worky      = vocn  (i,j,iblk)
          uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                         + worky*sin(ANGLET(i,j,iblk))
          vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                         - workx*sin(ANGLET(i,j,iblk))

          workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
          worky      = ss_tlty  (i,j,iblk)
          ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                            + worky*sin(ANGLET(i,j,iblk))
          ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

          sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)
          Tf (i,j,iblk) = -1.8_dbl_kind                 ! hardwired for NCOM
!         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
!         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

       enddo
       enddo
      enddo
      !$OMP END PARALLEL DO
      call t_stopf ('cice_imp_ocn')

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

      if (.not.prescribed_ice) then
         call t_startf ('cice_imp_t2u')
         call t2ugrid_vector(uocn)
         call t2ugrid_vector(vocn)
         call t2ugrid_vector(ss_tltx)
         call t2ugrid_vector(ss_tlty)
         call t_stopf ('cice_imp_t2u')
      end if

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      call t_startf ('cice_imp_atm')
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
  
         ! atmosphere
         workx      = uatm(i,j,iblk) ! wind velocity, m/s
         worky      = vatm(i,j,iblk) 
         uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                         + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
         vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                         - workx*sin(ANGLET(i,j,iblk))

         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call t_stopf ('cice_imp_atm')

end subroutine ice_import_esmf

!=======================================================================

subroutine ice_domain_esmf( dom, rc )

    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Array), intent(inout)     :: dom
    integer, intent(out)                :: rc
    !
    ! Local Variables
    !
    integer :: i, j, iblk, n               ! indices
    integer :: ilo, ihi, jlo, jhi          ! beginning and end of physical domain
    integer :: klon,klat,karea,kmask,kfrac ! domain fields
    type(block)             :: this_block  ! block information for current block

    real(R8),    pointer    :: fptr (:,:)
    !-------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components
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
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          fptr(klon, n) = TLON(i,j,iblk)*rad_to_deg 
          fptr(klat, n) = TLAT(i,j,iblk)*rad_to_deg 
          fptr(karea, n) = tarea(i,j,iblk)/(radius*radius)
          fptr(kmask, n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          if (trim(grid_type) == 'latlon') then
             fptr(kfrac, n) = ocn_gridcell_frac(i,j,iblk)
          else
             fptr(kfrac, n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          end if
       enddo    !i
       enddo    !j
    enddo       !iblk

end subroutine ice_domain_esmf

!=======================================================================
! BOP
!
! !ROUTINE: restart_filename
!
! !INTERFACE:
  character(len=char_len_long) function restart_filename( yr_spec, mon_spec, day_spec, sec_spec )
!
! !DESCRIPTION: 
! Create a restart filename 
!
! !USES:
  use ice_restart, only : restart_file
!
! !INPUT/OUTPUT PARAMETERS:
  integer         , intent(in)  :: yr_spec         ! Simulation year
  integer         , intent(in)  :: mon_spec        ! Simulation month
  integer         , intent(in)  :: day_spec        ! Simulation day
  integer         , intent(in)  :: sec_spec        ! Seconds into current simulation day
!
! EOP
!
  integer             :: i, n      ! Loop variables
  integer             :: year      ! Simulation year
  integer             :: month     ! Simulation month
  integer             :: day       ! Simulation day
  integer             :: ncsec     ! Seconds into current simulation day
  character(len=char_len_long) :: rdate ! char date for restart filename

  !-----------------------------------------------------------------
  ! Determine year, month, day and sec to put in filename
  !-----------------------------------------------------------------

  year  = yr_spec
  month = mon_spec
  day   = day_spec
  ncsec = sec_spec

  if ( year > 99999   ) then
     write(rdate,'(i6.6,"-",i2.2,"-",i2.2,"-",i5.5)') year,month,mday,ncsec
  else if ( year > 9999    ) then
     write(rdate,'(i5.5,"-",i2.2,"-",i2.2,"-",i5.5)') year,month,mday,ncsec
  else
     write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') year,month,mday,ncsec
  end if
  restart_filename = trim(restart_file) // "." // trim(rdate) 

end function restart_filename

end module ice_comp_esmf

