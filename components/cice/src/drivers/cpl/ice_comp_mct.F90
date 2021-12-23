module ice_comp_mct

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: ice_comp_mct
!
! !DESCRIPTION:
! CICE interface routine for the ccsm cpl7 mct system
!
! !USES:

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod,  only : shr_sys_abort, shr_sys_flush
! use shr_mem_mod,  only : shr_get_memusage, shr_init_memusage
  use shr_file_mod, only : shr_file_getlogunit, shr_file_getloglevel,  &
		           shr_file_setloglevel, shr_file_setlogunit
  use mct_mod
#ifdef USE_ESMF_LIB
  use esmf
#else
  use esmf, only: ESMF_clock
#endif

  use seq_flds_mod
  use seq_cdata_mod,   only : seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod,only : seq_infodata_type, seq_infodata_getdata,       &
		              seq_infodata_putdata, seq_infodata_start_type_cont, &
		              seq_infodata_start_type_brnch, seq_infodata_start_type_start
  use seq_timemgr_mod, only : seq_timemgr_eclockgetdata, &
                              seq_timemgr_restartalarmison, &
		              seq_timemgr_eclockdateinsync, &
                              seq_timemgr_stopalarmison
  use seq_comm_mct,    only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use perf_mod,        only : t_startf, t_stopf, t_barrierf

  use ice_cpl_indices
  use ice_import_export
  use ice_state,       only : aice, filename_aero, filename_iage, &
                              filename_volpn, filename_FY, filename_lvl, &
                              tr_aero, tr_iage, tr_FY, tr_pond, tr_lvl
  use ice_domain_size, only : nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_domain,      only : nblocks, blocks_ice, halo_info, distrb_info, profile_barrier
  use ice_blocks,      only : block, get_block, nx_block, ny_block
  use ice_grid,        only : tlon, tlat, tarea, tmask, anglet, hm, ocn_gridcell_frac, &
 		              grid_type, t2ugrid_vector, gridcpl_file
  use ice_constants,   only : c0, c1, spval_dbl, rad_to_deg, radius
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
  use ice_scam,        only : scmlat, scmlon, single_column, scm_multcols
  use ice_fileunits,   only : nu_diag, inst_index, inst_name, inst_suffix
  use ice_dyn_evp,     only : kdyn
  use ice_prescribed_mod
  use ice_step_mod
  use ice_global_reductions
  use ice_broadcast
  use CICE_RunMod

! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: ice_SetGSMap_mct
  private :: ice_domain_mct
  private :: ice_setdef_mct
  private :: ice_coffset_mct
  private :: ice_setcoupling_mct

!
! !PRIVATE VARIABLES

  integer (kind=int_kind) :: ICEID       

  !--- for coupling on other grid from gridcpl_file ---
  type(mct_gsMap) :: gsMap_iloc  ! local gsmaps
  type(mct_gGrid) :: dom_iloc                 ! local domain
  type(mct_aVect) :: x2i_iloc, i2x_iloc
  type(mct_rearr) :: rearr_ice2iloc
  type(mct_rearr) :: rearr_iloc2ice
  integer         :: nxcpl, nycpl  ! size of coupling grid
  logical         :: other_cplgrid    ! using different coupling grid

!=======================================================================

contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_init_mct
!
! !INTERFACE:
  subroutine ice_init_mct( EClock, cdata_i, x2i_i, i2x_i, NLFilename )
!
! !DESCRIPTION:
! Initialize thermodynamic ice model and obtain relevant atmospheric model
! arrays back from driver 
!
! !USES:

    use CICE_InitMod
    use ice_restart, only: runid, runtype, restart_dir, restart_format
    use ice_history, only: history_dir, history_file
!
! !ARGUMENTS:
    type(ESMF_Clock)         , intent(inout) :: EClock
    type(seq_cdata)          , intent(inout) :: cdata_i
    type(mct_aVect)          , intent(inout) :: x2i_i, i2x_i
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename
!
! !LOCAL VARIABLES:
!
    type(mct_gsMap)             , pointer :: gsMap_ice
    type(mct_gGrid)             , pointer :: dom_i
    type(seq_infodata_type)     , pointer :: infodata   ! Input init object
    integer                               :: lsize,lsize_loc
    integer                               :: xoff,yoff
    integer                               :: nxg,nyg
    integer                               :: k
 
    type(mct_gsMap) :: gsmap_extend  ! local gsmaps

    character(len=256) :: drvarchdir         ! driver archive directory
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
    integer            :: mpicom_loc  ! temporary mpicom
    logical (kind=log_kind) :: atm_aero
    real(r8) :: mrss, mrss0,msize,msize0
    character(len=*), parameter  :: SubName = "ice_init_mct"
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Determine attribute vector indices
    !--------------------------------------------------------------------------

    call ice_cpl_indices_set()

    !---------------------------------------------------------------------------
    ! Set cdata pointers
    !---------------------------------------------------------------------------

    call seq_cdata_setptrs(cdata_i, ID=ICEID, mpicom=mpicom_loc, &
         gsMap=gsMap_ice, dom=dom_i, infodata=infodata)

    ! Determine time of next atmospheric shortwave calculation
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    ! Determine if aerosols are coming from the coupler
    call seq_infodata_GetData(infodata, atm_aero=atm_aero )

    ! Determine orbital parameters
    call seq_infodata_GetData(infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr)

    !   call shr_init_memusage()

    !---------------------------------------------------------------------------
    ! use infodata to determine type of run
    !---------------------------------------------------------------------------

    ! Preset single column values

    single_column = .false.
    scm_multcols = .false.
    scmlat = -999.
    scmlon = -999.

    call seq_infodata_GetData( infodata, case_name=runid   ,  &  
       single_column=single_column,scm_multcols=scm_multcols, &
       scmlat=scmlat,scmlon=scmlon)
    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(nu_diag,*) 'ice_comp_mct ERROR: unknown starttype'
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

    !---------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !---------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)
   
    !---------------------------------------------------------------------------
    ! use EClock to reset calendar information on initial start
    !---------------------------------------------------------------------------

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

    call seq_timemgr_EClockGetData(EClock,               &
         start_ymd=start_ymd, start_tod=start_tod,       &
         curr_ymd=curr_ymd,   curr_tod=curr_tod,         &
         ref_ymd=ref_ymd,     ref_tod=ref_tod)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) 'ice_comp_mct: ref_ymd ',ref_ymd, &
                  ' must equal start_ymd ',start_ymd
             write(nu_diag,*) 'ice_comp_mct: ref_ymd ',ref_tod, &
                  ' must equal start_ymd ',start_tod
          end if
       end if

       if (my_task == master_task) then
          write(nu_diag,*) '(ice_init_mct) idate from sync clock = ', &
               start_ymd
          write(nu_diag,*) '(ice_init_mct)   tod from sync clock = ', &
               start_tod
          write(nu_diag,*) &
               '(ice_init_mct) resetting idate to match sync clock'
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
    ! Initialize MCT attribute vectors and indices
    !---------------------------------------------------------------------------

    call t_startf ('cice_mct_init')

    ! Initialize ice gsMap

    if (trim(gridcpl_file) == 'unknown_gridcpl_file') then
       call ice_SetGSMap_mct( MPI_COMM_ICE, ICEID, GSMap_ice ) 
       lsize = mct_gsMap_lsize(gsMap_ice, MPI_COMM_ICE)
       call ice_domain_mct( lsize, gsMap_ice, dom_i )
       other_cplgrid = .false.
       nxg = nx_global
       nyg = ny_global
    else
       call ice_SetGSMap_mct( MPI_COMM_ICE, ICEID, GSMap_iloc ) 
       lsize_loc = mct_gsMap_lsize(gsMap_iloc, MPI_COMM_ICE)
       call ice_domain_mct( lsize_loc, gsMap_iloc, dom_iloc )
 
       call ice_setcoupling_mct(MPI_COMM_ICE, ICEID, gsmap_ice, dom_i)
       lsize = mct_gsMap_lsize(gsMap_ice, MPI_COMM_ICE)
 
       call ice_coffset_mct(xoff,yoff,gsmap_iloc,dom_iloc,gsmap_ice,dom_i,MPI_COMM_ICE)
 
       call ice_SetGSMap_mct( MPI_COMM_ICE, ICEID, gsmap_extend, xoff, yoff, nxcpl, nycpl)
       if (lsize_loc /= mct_gsmap_lsize(gsmap_extend,MPI_COMM_ICE)) then
          write(nu_diag,*) subname,' :: gsmap_extend extended ',lsize_loc, &
             mct_gsmap_lsize(gsmap_extend,MPI_COMM_ICE)
          call shr_sys_abort(subname//' :: error in gsmap_extend extended')
       endif

       call mct_rearr_init(gsmap_ice, gsmap_extend, MPI_COMM_ICE, rearr_ice2iloc)
       call mct_rearr_init(gsmap_extend, gsmap_ice, MPI_COMM_ICE, rearr_iloc2ice)
       call mct_aVect_init(x2i_iloc, rList=seq_flds_x2i_fields, lsize=lsize_loc)
       call mct_aVect_zero(x2i_iloc)
       call mct_aVect_init(i2x_iloc, rList=seq_flds_i2x_fields, lsize=lsize_loc)
       call mct_aVect_zero(i2x_iloc)
       call mct_gsmap_clean(gsmap_extend)
 
       other_cplgrid = .true.
       nxg = nxcpl
       nyg = nycpl
    endif

    ! Inialize mct attribute vectors

    call mct_aVect_init(x2i_i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i_i)

    call mct_aVect_init(i2x_i, rList=seq_flds_i2x_fields, lsize=lsize) 
    call mct_aVect_zero(i2x_i)

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    if (other_cplgrid) then
       call ice_prescribed_init(ICEID, gsmap_iloc, dom_iloc)
    else
       call ice_prescribed_init(ICEID, gsmap_ice, dom_i)
    endif

    !-----------------------------------------------------------------
    ! Get ready for coupling
    !-----------------------------------------------------------------

    call coupling_prep

    !---------------------------------------------------------------------------
    ! Fill in export state for driver
    !---------------------------------------------------------------------------

    if (other_cplgrid) then
       call ice_export (i2x_iloc%rattr)  !Send initial state to driver
       call ice_setdef_mct ( i2x_i )
       call mct_rearr_rearrange(i2x_iloc, i2x_i, rearr_iloc2ice)
    else
       call ice_export (i2x_i%rattr)  !Send initial state to driver
    endif
    call seq_infodata_PutData( infodata, ice_prognostic=.true., &
      iceberg_prognostic=.false., ice_nx = nxg, ice_ny = nyg )
    call t_stopf ('cice_mct_init')

    ! Error check
    if (tr_aero .and. .not. atm_aero) then
       write(nu_diag,*) 'ice_import ERROR: atm_aero must be set for tr_aero' 
       call shr_sys_abort()
    end if

    !---------------------------------------------------------------------------
    ! Reset shr logging to original values
    !---------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    !   call ice_timer_stop(timer_total) ! time entire run
    !   call shr_get_memusage(msize,mrss)
    !   call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,'ice_init_mct mrss0')
    !   call shr_mpi_max(msize,msize0,MPI_COMM_ICE,'ice_init_mct msize0')
    !   if(my_task == 0) then
    !   write(shrlogunit,105) 'ice_init_mct: memory_write: model date = ',start_ymd,start_tod, &
    !           ' memory = ',msize0,' MB (highwater)    ',mrss0,' MB (usage)'
    !   endif
 
  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ice_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_run_mct
!
! !INTERFACE:
  subroutine ice_run_mct( EClock, cdata_i, x2i_i, i2x_i )
!
! !DESCRIPTION:
! Run thermodynamic CICE
!
! !USES:
    use ice_history
    use ice_restart
    use ice_diagnostics
    use ice_aerosol   , only : write_restart_aero
    use ice_age       , only : write_restart_age
    use ice_meltpond  , only : write_restart_pond
    use ice_FY        , only : write_restart_FY
    use ice_lvl       , only : write_restart_lvl
    use ice_restoring , only : restore_ice, ice_HaloRestore
    use ice_shortwave , only : init_shortwave

! !ARGUMENTS:
    type(ESMF_Clock),intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_i
    type(mct_aVect), intent(inout) :: x2i_i
    type(mct_aVect), intent(inout) :: i2x_i

! !LOCAL VARIABLES:
    integer :: k             ! index
    logical :: rstwr         ! .true. ==> write a restart file
    logical :: stop_now      ! .true. ==> stop at the end of this run phase
    integer :: ymd           ! Current date (YYYYMMDD)
    integer :: tod           ! Current time of day (sec)
    integer :: yr_sync       ! Sync current year
    integer :: mon_sync      ! Sync current month
    integer :: day_sync      ! Sync current day
    integer :: tod_sync      ! Sync current time of day (sec)
    integer :: ymd_sync      ! Current year of sync clock
    integer :: curr_ymd           ! Current date (YYYYMMDD)
    integer :: curr_tod           ! Current time of day (s)
    integer :: shrlogunit,shrloglev ! old values
    integer :: lbnum
    integer :: n, nyrp
    type(mct_gGrid)        , pointer :: dom_i
    type(seq_infodata_type), pointer :: infodata   
    type(mct_gsMap)        , pointer :: gsMap_i
    character(len=char_len_long) :: fname
    character(len=char_len_long) :: string1, string2
    character(len=*), parameter  :: SubName = "ice_run_mct"

    real(r8) :: mrss, mrss0,msize,msize0
    logical, save :: first_time = .true.

!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    call ice_timer_start(timer_total) ! time entire run
    call t_barrierf('cice_run_total_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_total')

    !---------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !---------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)
   
    call seq_cdata_setptrs(cdata_i, infodata=infodata, dom=dom_i, &
         gsMap=gsMap_i)

    ! Determine time of next atmospheric shortwave calculation
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    ! Determine orbital parameters
    call seq_infodata_GetData(infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
                              orb_lambm0=lambm0, orb_obliqr=obliqr)

    ! Get clock information
    call seq_timemgr_EClockGetData(EClock,               &
         curr_ymd=curr_ymd, curr_tod=curr_tod)

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
    
    call t_barrierf('cice_run_import_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_import')
    call ice_timer_start(timer_cplrecv)
    if (other_cplgrid) then
       call mct_rearr_rearrange(x2i_i, x2i_iloc, rearr_ice2iloc)
       call ice_import( x2i_iloc%rattr )
    else
       call ice_import( x2i_i%rattr )
    endif
    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')
 
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

    call t_barrierf('cice_run_initmd_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_initmd')
    call init_mass_diags   ! diagnostics per timestep
    call t_stopf ('cice_run_initmd')

    if(prescribed_ice) then  ! read prescribed ice
       call t_barrierf('cice_run_presc_BARRIER',MPI_COMM_ICE)
       call t_startf ('cice_run_presc')
       call ice_prescribed_run(idate, sec)
       call t_stopf ('cice_run_presc')
    endif
    
    call t_barrierf('cice_run_initflux_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_initflux')
    call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
    call init_flux_ocn        ! initialize ocean fluxes sent to coupler
    call t_stopf ('cice_run_initflux')

    !-----------------------------------------------------------------
    ! Scale radiation fields
    !-----------------------------------------------------------------

    call t_barrierf('cice_run_prepradiation_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_prepradiation')
    call ice_timer_start(timer_sw)
    call prep_radiation(dt)
    call ice_timer_stop(timer_sw)
    call t_stopf ('cice_run_prepradiation')
    
    !-----------------------------------------------------------------
    ! thermodynamics1
    !-----------------------------------------------------------------

    call t_barrierf('cice_run_therm1_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_therm1')
    call step_therm1(dt)
    call t_stopf ('cice_run_therm1')
    
    !-----------------------------------------------------------------
    ! thermodynamics2
    !-----------------------------------------------------------------

    if (.not.prescribed_ice) then
       call t_barrierf('cice_run_therm2_BARRIER',MPI_COMM_ICE)
       call t_startf ('cice_run_therm2')
       call step_therm2 (dt)  ! post-coupler thermodynamics
       call t_stopf ('cice_run_therm2')
    end if

   !-----------------------------------------------------------------
   ! dynamics, transport, ridging
   !-----------------------------------------------------------------

    call t_barrierf('cice_run_dyn_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_dyn')
    if (.not.prescribed_ice .and. kdyn>0) then
       if (xndt_dyn > c1) then
          do k = 1, nint(xndt_dyn)
             call step_dynamics(dt_dyn,dt) ! dynamics, transport, ridging
          enddo
       else
          if (mod(time, dt_dyn) == c0) then
             call step_dynamics(dt_dyn,dt) ! dynamics, transport, ridging
          endif
       endif
    endif ! not prescribed_ice
    call t_stopf ('cice_run_dyn')
    
    !-----------------------------------------------------------------
    ! radiation
    !-----------------------------------------------------------------

    call t_barrierf('cice_run_radiation_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_radiation')
    call ice_timer_start(timer_sw)
    call step_radiation(dt)
    call ice_timer_stop(timer_sw)
    call t_stopf ('cice_run_radiation')
    
    !-----------------------------------------------------------------
    ! get ready for coupling
    !-----------------------------------------------------------------

    call t_barrierf('cice_run_couplingprep_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_couplingprep')
    call coupling_prep
    call t_stopf ('cice_run_couplingprep')

    call ice_timer_stop(timer_step)

    !-----------------------------------------------------------------
    ! write data
    !-----------------------------------------------------------------
    
    call t_barrierf('cice_run_diag_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_diag')
    call ice_timer_start(timer_diags)
    if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
    call ice_timer_stop(timer_diags)
    call t_stopf ('cice_run_diag')
    
    call t_barrierf('cice_run_hist_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_hist')
    call ice_timer_start(timer_hist)
#if (defined _NOIO)
!  Not enought memory on BGL to write a history file yet! 
!    call ice_write_hist (dt)    ! history file
#else
    call ice_write_hist (dt)    ! history file
#endif
    call ice_timer_stop(timer_hist)
    call t_stopf ('cice_run_hist')
 
    !--------------------------------------------
    ! Accumualate the number of active ice cells
    !--------------------------------------------
    call t_barrierf('cice_run_accum_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_accum')
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
    call t_stopf ('cice_run_accum')

    !-----------------------------------------------------------------
    ! send export state to driver 
    !-----------------------------------------------------------------
    
    call t_barrierf('cice_run_export_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_export')
    call ice_timer_start(timer_cplsend)
    if (other_cplgrid) then
       call ice_export ( i2x_iloc%rattr )
       call ice_setdef_mct ( i2x_i )
       call mct_rearr_rearrange(i2x_iloc, i2x_i, rearr_iloc2ice)
    else
       call ice_export ( i2x_i%rattr )
    endif
    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_run_export')
    
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
    
!   if(tod == 0) then
!      call shr_get_memusage(msize,mrss)
!      call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,'ice_run_mct mrss0')
!      call shr_mpi_max(msize,msize0,MPI_COMM_ICE,'ice_run_mct msize0')
!      if(my_task == 0 ) then
!          write(shrlogunit,105) 'ice_run_mct: memory_write: model date = ',ymd,tod, &
!               ' memory = ',msize0,' MB (highwater)    ',mrss0,' MB (usage)'
!      endif
!   endif
    call t_stopf ('cice_run_total')
 
  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ice_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_final_mct
!
! !INTERFACE:
  subroutine ice_final_mct( EClock, cdata_i, x2i_i, i2x_i )
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

    type(ESMF_Clock),intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_i
    type(mct_aVect), intent(inout) :: x2i_i
    type(mct_aVect), intent(inout) :: i2x_i
!
! !REVISION HISTORY:
!
!EOP
!---------------------------------------------------------------------------

  end subroutine ice_final_mct

!===============================================================================

  subroutine ice_SetGSMap_mct( mpicom, ID, gsMap_ice, xoff, yoff, nxgin, nygin )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    implicit none
    integer        , intent(in)    :: mpicom
    integer        , intent(in)    :: ID
    type(mct_gsMap), intent(inout) :: gsMap_ice
    integer,optional, intent(in)   :: xoff   ! x offset
    integer,optional, intent(in)   :: yoff   ! y offset
    integer,optional, intent(in)   :: nxgin ! global size
    integer,optional, intent(in)   :: nygin ! global size
    !
    ! Local variables
    !
    integer,allocatable :: gindex(:)
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize,gsize
    integer     :: lxoff,lyoff,nxg,nyg
    integer     :: ier
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    !-------------------------------------------------------------------

    ! Build the CICE grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP

    lxoff = 1
    lyoff = 1
    if (present(xoff)) then
       lxoff = xoff
    endif
    if (present(yoff)) then
       lyoff = yoff
    endif

    nxg = nx_global
    nyg = ny_global
    if (present(nxgin)) then
       nxg = nxgin
    endif
    if (present(nygin)) then
       nyg = nygin
    endif
    gsize = nxg*nyg

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
             lon = this_block%i_glob(i) + lxoff - 1
             lat = this_block%j_glob(j) + lyoff - 1
             gi = (lat-1)*nxg + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo        !iblk
    
    call mct_gsMap_init( gsMap_ice, gindex, mpicom, ID, lsize, gsize )

    deallocate(gindex)

  end subroutine ice_SetGSMap_mct

  subroutine ice_domain_mct( lsize, gsMap_i, dom_i )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i     
    !
    ! Local Variables
    !
    integer :: i, j, iblk, n, gi           ! indices
    integer :: ilo, ihi, jlo, jhi          ! beginning and end of physical domain
    real(dbl_kind), pointer :: work_dom(:) ! temporary
    real(dbl_kind), pointer :: data(:)     ! temporary
    integer       , pointer :: idata(:)    ! temporary
    type(block)             :: this_block  ! block information for current block
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
    !
    call mct_gGrid_init( GGrid=dom_i, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)
    !  
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize) 
    !
    ! Fill in correct values for domain components
    !
    allocate(work_dom(lsize)) 
    work_dom(:) = 0.0_dbl_kind

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
          data(n) = TLON(i,j,iblk)*rad_to_deg 
       enddo    !i
       enddo    !j
    enddo       !iblk
    call mct_gGrid_importRattr(dom_i,"lon",data,lsize) 

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
          data(n) = TLAT(i,j,iblk)*rad_to_deg 
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"lat",data,lsize) 

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
          data(n) = tarea(i,j,iblk)/(radius*radius)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"area",data,lsize) 

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
          data(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"mask",data,lsize) 

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
	  if (trim(grid_type) == 'latlon') then
             data(n) = ocn_gridcell_frac(i,j,iblk)
          else
             data(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          end if
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)
    deallocate(work_dom)

  end subroutine ice_domain_mct

  !=======================================================================

  subroutine ice_setdef_mct( i2x_i )   

    implicit none

    !-----------------------------------------------------
    type(mct_aVect)   , intent(inout) :: i2x_i

    !-----------------------------------------------------

    call mct_aVect_zero(i2x_i)

    ! tcraig : this is where observations could be read in

  end subroutine ice_setdef_mct

  !=======================================================================

  subroutine ice_coffset_mct(xoff,yoff,gsmap_a,dom_a,gsmap_b,dom_b,mpicom_i)
    implicit none

    integer        , intent(out)   :: xoff
    integer        , intent(out)   :: yoff
    type(mct_gsmap), intent(in)    :: gsmap_a
    type(mct_ggrid), intent(in)    :: dom_a
    type(mct_gsmap), intent(in)    :: gsmap_b
    type(mct_ggrid), intent(in)    :: dom_b
    integer        , intent(in)    :: mpicom_i

    type(mct_aVect) :: ava
    type(mct_aVect) :: avag
    integer :: k1,k2,k
    integer :: npt
    integer :: noff,noffg
    real(dbl_kind) :: x1,y1,x2,y2
    real(dbl_kind) :: dist,distmin,distming
    integer :: lsizea,lsizeb
    integer :: iam,ierr
    integer, pointer :: ipoints(:)
    character(len=*),parameter :: subname = "ice_coffset_mct"

    call mpi_comm_rank(mpicom_i,iam,ierr)

    lsizea = mct_aVect_lsize(dom_a%data)
    lsizeb = mct_aVect_lsize(dom_b%data)

    !--- compute lon/lat at dom_a (local) point (1,1)

    call mct_aVect_init(ava,rList='lon:lat',lsize=lsizea)
    call mct_aVect_copy(dom_a%data,ava,'lon:lat')
    call mct_aVect_gather(ava,avag,gsmap_a,0,mpicom_i)

    if (iam == 0) then
       k1 = mct_aVect_indexRA(avag,'lon',dieWith=subname//'_avag')
       k2 = mct_aVect_indexRA(avag,'lat',dieWith=subname//'_avag')
       npt = 1   ! actual corner points screwed up by U average/wraparound
       npt = nx_global + 2  ! use global point (2,2)
       x1 = mod(avag%rAttr(k1,npt)+360.0_r8,360.0_r8)
       y1 = avag%rAttr(k2,npt)
    endif

    call mct_aVect_clean(avag)
    call mct_aVect_clean(ava)

    call shr_mpi_bcast(x1,mpicom_i)
    call shr_mpi_bcast(y1,mpicom_i)

    !--- find x1,y1 point in dom_b (extended grid)

    noff = -1
    noffg = -1

    call mct_gsMap_orderedPoints(gsMap_b, iam, ipoints)
    if (size(ipoints) /= lsizeb) then
       write(nu_diag,*) subname,' size ipoints = ',size(ipoints),lsizeb
       call shr_sys_abort(subname//' :: error size of ipoints')
    endif

    k1 = mct_aVect_indexRA(dom_b%data,'lon',dieWith=subname//'_domb')
    k2 = mct_aVect_indexRA(dom_b%data,'lat',dieWith=subname//'_domb')
    distmin = 1.0e36
    do k = 1,lsizeb
       x2 = mod(dom_b%data%rAttr(k1,k)+360.0_r8,360.0_r8)
       y2 = dom_b%data%rAttr(k2,k)
       dist = abs((x1-x2)*(x1-x2))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
       dist = abs((x1-x2-360.0_r8)*(x1-x2-360.0_r8))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
       dist = abs((x1-x2+360.0_r8)*(x1-x2+360.0_r8))+abs((y1-y2)*(y1-y2))
       if (dist < distmin) then
          distmin = dist
          noff = ipoints(k)
       endif
    enddo

    deallocate(ipoints)

    call shr_mpi_min(distmin,distming,mpicom_i,'distmin',all=.true.)

    if (distming /= distmin) then
       noff = -1
    endif

    call shr_mpi_max(noff,noffg,mpicom_i,'noffg',all=.true.)

    ! subtract extra -1 and -nxcpl for point (2,2)
    xoff = mod(noffg-1-1,nxcpl) + 1
    yoff = (noffg-1-nxcpl)/nxcpl + 1

    if (iam == 0) then
       write(nu_diag,*) subname,' :: x1,y1  = ',x1,y1
       write(nu_diag,*) subname,' :: offset = ',noffg,xoff,yoff
       call shr_sys_flush(nu_diag)
    endif

    if (noffg < 1) then
       call shr_sys_abort(subname//' :: noffg lt 1')
    endif

  end subroutine ice_coffset_mct

  !=======================================================================

  subroutine ice_setcoupling_mct(mpicom_i, ICEID, gsmap_i, dom_i)

    implicit none
    include 'netcdf.inc'

    integer        , intent(in)    :: mpicom_i
    integer        , intent(in)    :: ICEID
    type(mct_gsmap), intent(inout) :: gsmap_i
    type(mct_ggrid), intent(inout) :: dom_i

    integer :: n     ! counter
    integer :: iam   ! pe rank
    integer :: npes  ! number of pes
    integer :: ierr  ! error code
    integer :: rcode ! error code
    integer :: nx,ny ! grid size
    integer :: gsize ! global size
    integer :: lsize ! local size
    integer, pointer :: start(:),length(:),pe_loc(:)
    integer, pointer :: idata(:)
    real(dbl_kind),pointer :: data(:)
    type(mct_avect) :: avg, av1
    integer :: fid,did,vid
    character(len=8) :: avfld,dofld
    character(len=*), parameter  :: SubName = "ice_setcoupling_mct"

    call MPI_comm_rank(mpicom_i,iam,ierr)
    call MPI_comm_size(mpicom_i,npes,ierr)

    allocate(start(npes),length(npes),pe_loc(npes))

    if (iam == 0) then
       rcode = nf_open(gridcpl_file(1:len_trim(gridcpl_file)),NF_NOWRITE,fid)
       rcode = nf_inq_dimid (fid, 'ni', did)
       rcode = nf_inq_dimlen(fid, did, nx)
       rcode = nf_inq_dimid (fid, 'nj', did)
       rcode = nf_inq_dimlen(fid, did, ny)
       gsize = nx*ny
       nxcpl = nx
       nycpl = ny

       length = gsize / npes
       do n = 1,npes
          if (n <= mod(gsize,npes)) length(n) = length(n) + 1
       enddo

       start(1) = 1
       pe_loc(1) = 0
       do n = 2,npes    
          pe_loc(n) = n-1
          start(n) = start(n-1) + length(n-1)
       enddo
       if ((start(npes) + length(npes) - 1) /= gsize) then
          write(nu_diag,*) &
            subname,' gsize, start, length = ',gsize,start(npes),length(npes)
          call shr_sys_flush(nu_diag)
          call shr_sys_abort( SubName//":: decomp inconsistent")
       endif

       write(nu_diag,*) subname,' read ',trim(gridcpl_file)
       write(nu_diag,*) subname,' size ',nx,ny,gsize
    endif

    call shr_mpi_bcast(nxcpl,mpicom_i)
    call shr_mpi_bcast(nycpl,mpicom_i)
    call shr_mpi_bcast(gsize,mpicom_i)
    call mct_gsmap_init(gsmap_i,npes,start,length,pe_loc,0,mpicom_i,ICEID,gsize)
    deallocate(start,length,pe_loc)

    lsize = mct_gsmap_lsize(gsmap_i,mpicom_i)
    call mct_gGrid_init( GGrid=dom_i, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT

    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    deallocate(idata)

    ! Initialize attribute vector with special value

    allocate(data(lsize))
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize) 
    deallocate(data)

    ! Read domain arrays

    if (iam == 0) then
       call mct_avect_init(avg,rList='fld',lsize=gsize)
    endif

    do n = 1,5

       if (n == 1) avfld = 'lat'
       if (n == 1) dofld = 'yc'
       if (n == 2) avfld = 'lon'
       if (n == 2) dofld = 'xc'
       if (n == 3) avfld = 'area'
       if (n == 3) dofld = 'area'
       if (n == 4) avfld = 'frac'
       if (n == 4) dofld = 'frac'
       if (n == 5) avfld = 'mask'
       if (n == 5) dofld = 'mask'
       if (iam == 0) then
          rcode = nf_inq_varid(fid,trim(dofld),vid)
          if (n == 5) then
             allocate(idata(gsize))
             rcode = nf_get_var_int(fid,vid,idata)
             avg%rAttr(1,:) = idata
             deallocate(idata)
          else
             rcode = nf_get_var_double(fid,vid,avg%rAttr(1,:))
          endif
       endif

       call mct_aVect_scatter(avg,av1,gsmap_i,0,mpicom_i)
       call mct_aVect_copy(av1,dom_i%data,'fld',avfld)

       if (iam == 0) then
          call mct_avect_clean(av1)
       endif

    enddo

    if (iam == 0) then
       call mct_avect_clean(avg)
    endif

  end subroutine ice_setcoupling_mct
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

end module ice_comp_mct

