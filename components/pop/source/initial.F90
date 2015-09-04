!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module initial

!BOP
! !MODULE: initial
! !DESCRIPTION:
!  This module contains routines for initializing a POP simulation,
!  mostly by calling individual initialization routines for each
!  POP module.
!
! !REVISION HISTORY:
!  SVN:$Id: initial.F90 44694 2013-03-12 19:58:14Z mlevy@ucar.edu $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod
   use POP_SolversMod
   use POP_ReductionsMod

   use kinds_mod, only: i4, i8, r8, int_kind, log_kind, char_len
   use blocks, only: block, nx_block, ny_block, get_block
   use domain_size
   use domain, only: nblocks_clinic, blocks_clinic, init_domain_blocks,    &
       init_domain_distribution, distrb_clinic
   use constants, only: radian, delim_fmt, blank_fmt, field_loc_center, blank_fmt, &
       c0, ppt_to_salt, mpercm, c1, field_type_scalar, init_constants,     &
       stefan_boltzmann, latent_heat_vapor_mks, vonkar, emissivity, &
       latent_heat_fusion, t0_kelvin, pi, ocn_ref_salinity,     &
       sea_ice_salinity, radius, cp_sw, grav, omega,cp_air,     &
       rho_fw, sound, rho_air, rho_sw, ndelim_fmt
   use communicate, only: my_task, master_task, init_communicate
   use budget_diagnostics, only: init_budget_diagnostics
   use broadcast, only: broadcast_array, broadcast_scalar
   use prognostic, only: init_prognostic, max_blocks_clinic, nx_global,    &
       ny_global, km, nt, TRACER, curtime, RHO, newtime, oldtime
   use grid, only: init_grid1, init_grid2, kmt, kmt_g, n_topo_smooth, zt,    &
       fill_points, sfc_layer_varthick, sfc_layer_type, TLON, TLAT, partial_bottom_cells
   use io
   use io_tools
   use baroclinic, only: init_baroclinic
   use barotropic, only: init_barotropic
   use pressure_grad, only: init_pressure_grad
   use surface_hgt, only: init_surface_hgt
   use vertical_mix, only: init_vertical_mix, vmix_itype, vmix_type_kpp
   use vmix_kpp, only: bckgrnd_vdc2, linertial
   use horizontal_mix, only: init_horizontal_mix
   use advection, only: init_advection
   use diagnostics, only: init_diagnostics
   use state_mod, only: init_state, state, state_itype, state_type_mwjf, state_range_iopt, &
      state_range_enforce 
   use time_management, only: first_step, init_time1, init_time2, &
                              dttxcel, dtuxcel, check_time_flag_int,  &
                              get_time_flag_id, freq_opt_nhour
   use topostress, only: init_topostress
   use ice
   use output, only: init_output
   use tavg, only: ltavg_restart, tavg_id, set_in_tavg_contents,n_tavg_streams, tavg_streams
   !use hydro_sections
   !use current_meters
   !use drifters
   use forcing, only: init_forcing
   use forcing_sfwf, only: sfwf_formulation, lms_balance, sfwf_data_type, lfw_as_salt_flx
   use forcing_shf, only: luse_cpl_ifrac, OCN_WGT, shf_formulation, shf_data_type
   use forcing_ws, only: ws_data_type
   use sw_absorption, only: init_sw_absorption
   use passive_tracers, only: init_passive_tracers, ecosys_on
   use ecosys_mod, only: ecosys_qsw_distrb_const
   use exit_mod, only: sigAbort, exit_pop, flushm
   use restart, only: read_restart, restart_fmt, read_restart_filename
   use ms_balance, only: init_ms_balance
   use forcing_coupled, only: pop_init_coupled, pop_init_partially_coupled, &
       qsw_distrb_iopt, qsw_distrb_iopt_const, ncouple_per_day, coupled_freq_iopt
   use global_reductions, only: init_global_reductions, global_sum
   use timers, only: init_timers
   use registry
   use qflux_mod, only: init_qflux
   use niw_mixing
   use tidal_mixing
   use step_mod, only: init_step
   use gather_scatter
#ifdef CCSMCOUPLED
   use shr_ncread_mod
   use shr_map_mod
#endif
   use overflows
   use overflow_type

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: pop_init_phase1, pop_init_phase2

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      init_ts_file_fmt,    &! format (bin or nc) for input file
      exit_string           ! exit_POP message string

   logical (log_kind), public ::  &! context variables
      lcoupled,                   &! T ==> pop is coupled to another system
      lccsm,                      &! T ==> pop is being run in the ccsm context
      b4b_flag,                   &! T ==> pop is being run in the "bit-for-bit" mode
      lccsm_control_compatible     ! T ==> pop is being run with code that is b4b with the ccsm4 control run
                                   !       this is a temporary flag that will be removed in ccsm4_0_1

!EOC
!***********************************************************************

 contains
!***********************************************************************
!BOP
! !IROUTINE: pop_init_phase1
! !INTERFACE:

 subroutine pop_init_phase1(errorCode)

! !DESCRIPTION:
!  This routine is the first of a two-phase initialization process for
!  a POP run. It calls various module initialization routines and sets up 
!  the initial temperature and salinity
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                      &! dummy vertical level index
      ier                      ! error flag

!-----------------------------------------------------------------------
!
!  initialize message-passing or other communication protocol
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_communicate

!-----------------------------------------------------------------------
!
!  initialize registry, which keeps track of which initialization
!  routines have been called.  This feature is used for error checking
!  in routines whose calling order is important
!
!-----------------------------------------------------------------------

   call init_registry

!-----------------------------------------------------------------------
!
!  initialize constants and i/o stuff
!
!-----------------------------------------------------------------------

   call init_io

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  temporary synching of old and new pop2 infrastructure for CCSM
!
!-----------------------------------------------------------------------
   POP_stdout = stdout
   POP_stderr = stderr
   POP_stdin  = stdin
#endif

!-----------------------------------------------------------------------
!
!  initialize context in which pop is being run
!
!-----------------------------------------------------------------------

   call init_context

!-----------------------------------------------------------------------
!
!  write version information to output log after output redirection
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a)') ' Parallel Ocean Program (POP) '
      write(stdout,'(a)') ' Based on Version 2.1alpha Jan 2005'
      write(stdout,'(a)') ' Modified for CESM  2005-2010'
      write(stdout,blank_fmt)
      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif

   call init_constants

!-----------------------------------------------------------------------
!
!  initialize timers 
!
!-----------------------------------------------------------------------

   call init_timers

!-----------------------------------------------------------------------
!
!  initialize additional communication routines
!
!-----------------------------------------------------------------------

   call init_global_reductions
   call POP_initReductions

!-----------------------------------------------------------------------
!
!  initialize overflows, part I
!
!-----------------------------------------------------------------------

   call init_overflows1

!-----------------------------------------------------------------------
!
!  initialize domain and grid
!
!-----------------------------------------------------------------------

   call init_domain_blocks
   call init_grid1
   call init_domain_distribution(KMT_G)

!-----------------------------------------------------------------------
!
!  initialize overflows, part II. placed here so KMT_G scatter to
!  KMT can be done (and NOT in init_grid2) for possible KMT mods; then
!  finish with domain and grid initialization
!
!-----------------------------------------------------------------------

   call init_overflows2

   call init_grid2(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_phase1: error initializing grid 2')
      return
   endif

!-----------------------------------------------------------------------
!
!  compute time step and initialize time-related quantities
!
!-----------------------------------------------------------------------

   call init_time1

!-----------------------------------------------------------------------
!
!  initialize equation of state
!
!-----------------------------------------------------------------------

   call init_state

!-----------------------------------------------------------------------
!
!  calculate topographic stress (maximum entropy) velocities
!
!-----------------------------------------------------------------------

   call init_topostress(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_phase1: error in init_topostress')
      return
   endif

!-----------------------------------------------------------------------
!
!  initialize niw driven mixing
!
!-----------------------------------------------------------------------
   
   call init_niw_mixing

!-----------------------------------------------------------------------
!
!  initialize tidally driven mixing
!
!-----------------------------------------------------------------------

   call init_tidal_mixing

   if ( overflows_interactive  .and.  .not.ltidal_mixing ) then
     exit_string = 'FATAL ERROR: overflow code is validated only with tidal mixing'
     call document ('pop_init_phase1', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  initialize barotropic elliptic solver
!
!-----------------------------------------------------------------------

   call POP_SolversInit(errorCode)
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_Init: error initializing solvers')
      return
   endif

!-----------------------------------------------------------------------
!
!  modify 9pt coefficients for barotropic solver for overflow use
!
!-----------------------------------------------------------------------

   call init_overflows3

!-----------------------------------------------------------------------
!
!  initialize pressure gradient (pressure averaging)
!  initialize baroclinic (reset to freezing)
!  initialize barotropic (barotropic-related diagnostics)
!  initialize surface_hgt (ssh-related diagnostics)
!
!-----------------------------------------------------------------------

   call init_pressure_grad
   call init_baroclinic
   call init_barotropic
   call init_surface_hgt

!-----------------------------------------------------------------------
!
!  initialize prognostic fields
!
!-----------------------------------------------------------------------

   call init_prognostic

!-----------------------------------------------------------------------
!
!  initialize ice module
!
!-----------------------------------------------------------------------

   call init_ice

!-----------------------------------------------------------------------
!
!  set initial temperature and salinity profiles (includes read of
!  restart file
!
!-----------------------------------------------------------------------

   call init_ts(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_phase1: error in init_ts')
      return
   endif

!-----------------------------------------------------------------------
!
!  finish computing time-related quantities after restart info
!  available
!
!-----------------------------------------------------------------------

   call init_time2


!-----------------------------------------------------------------------
!
!  initialize fields for surface forcing
!       o init_ws
!       o init_shf 
!       o init_sfwf
!       o init_pt_interior
!       o init_s_interior
!       o init_ap
!
!-----------------------------------------------------------------------

   call init_forcing

!-----------------------------------------------------------------------
!
!  initialize generic aspects of coupled forcing (no coupling-specific 
!  references)
!
!-----------------------------------------------------------------------

   call pop_init_coupled

!-----------------------------------------------------------------------
!EOC

 end subroutine pop_init_phase1


!***********************************************************************
!BOP
! !IROUTINE: pop_init_phase2
! !INTERFACE:

 subroutine pop_init_phase2(errorCode)

! !DESCRIPTION:
!  This routine completes the two-phase initialization process for
!  a POP run. 
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode      ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k,                      &! dummy vertical level index
      ier                      ! error flag


!-----------------------------------------------------------------------
!
!     initialize passive tracer modules -- after call init_forcing_coupled
!     do this independently of nt so that
!     1) consistency of nt and selected passive tracer modules
!        can always be checked
!     2) passive_tavg_nonstd_vars gets allocated
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_passive_tracers(init_ts_file_fmt, read_restart_filename, &
                             errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_phase2: error in init_passive_tracers')
      return
   endif

!-----------------------------------------------------------------------
!
!  initialize vertical mixing variables
!  initialize horizontal mixing variables
!
!-----------------------------------------------------------------------

   call init_vertical_mix

   call init_horizontal_mix(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_phase1: error in init_horizontal_mix')
      return
   endif

!-----------------------------------------------------------------------
!
!  initialize overflow regional values
!
!-----------------------------------------------------------------------

  call init_overflows5

!-----------------------------------------------------------------------
!
!  initialize advection variables
!
!-----------------------------------------------------------------------

  call init_advection

!-----------------------------------------------------------------------
!
!  initialize shortwave absorption
!
!-----------------------------------------------------------------------

   call init_sw_absorption

!-----------------------------------------------------------------------
!
!  partial coupling forcing initialization
!
!-----------------------------------------------------------------------

   call pop_init_partially_coupled

!-----------------------------------------------------------------------
!
!  initialize time-averaged qflux information
!
!-----------------------------------------------------------------------

   call init_qflux

!-----------------------------------------------------------------------
!
!  initialize ms_balance
!
!-----------------------------------------------------------------------
 
   if (lcoupled .and. lms_balance) call init_ms_balance

!-----------------------------------------------------------------------
!
!  initialize diagnostics
!
!-----------------------------------------------------------------------

   call init_diagnostics

!-----------------------------------------------------------------------
!
!  initialize overflows output diagnostics filename
!
!-----------------------------------------------------------------------

   call init_overflows4

!-----------------------------------------------------------------------
!
!  initialize output; subroutine init_output calls 
!       o init_restart
!       o init_history
!       o init_movie
!       o init_tavg
!
!-----------------------------------------------------------------------

   call init_output

!-----------------------------------------------------------------------
!
!  initialize drifters, hydrographic sections and current meters
!
!-----------------------------------------------------------------------

   !call init_drifters
   !call init_hydro_sections
   !call init_current_meters

!-----------------------------------------------------------------------
!
!  initialize global budget diagnostics
!
!-----------------------------------------------------------------------

   call init_budget_diagnostics

!-----------------------------------------------------------------------
!
!  initialize step timers
!
!-----------------------------------------------------------------------

   call init_step

!-----------------------------------------------------------------------
!
!  check registry -- have any errors occured?
!
!-----------------------------------------------------------------------

   call trap_registry_failure

!-----------------------------------------------------------------------
!
!  check consistency of model options 
!
!-----------------------------------------------------------------------

   call POP_check

!-----------------------------------------------------------------------
!
!  write model information into log file
!
!-----------------------------------------------------------------------

   call document_constants


!-----------------------------------------------------------------------
!EOC

 end subroutine pop_init_phase2



!***********************************************************************
!BOP
! !IROUTINE: init_context
! !INTERFACE:

 subroutine init_context

! !DESCRIPTION:
!  This routine initializes the context in which POP is being run,
!    including information about coupling and CCSM
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


   integer (int_kind)   ::  &
      nml_error,            &! namelist i/o error flag
      number_of_fatal_errors

   namelist /context_nml/ lcoupled, lccsm, b4b_flag, lccsm_control_compatible

!-----------------------------------------------------------------------
!
!  read context_nml namelist to determine the context in which pop
!  being run. check for errors and broadcast info to all processors
!
!-----------------------------------------------------------------------

   lcoupled                 = .false.
   lccsm                    = .false.
   b4b_flag                 = .false.
   lccsm_control_compatible = .true.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=context_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     exit_string = 'FATAL ERROR: reading context_nml'
     call document ('init_context', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   call broadcast_scalar(lcoupled,                 master_task)
   call broadcast_scalar(lccsm,                    master_task)
   call broadcast_scalar(b4b_flag,                 master_task)
   call broadcast_scalar(lccsm_control_compatible, master_task)

!-----------------------------------------------------------------------
!
!  register information with the registry function, allowing other
!  modules to access this information (avoids circular dependencies)
!
!-----------------------------------------------------------------------
      if (lcoupled)                 call register_string('lcoupled')
      if (lccsm)                    call register_string('lccsm')
      if (b4b_flag)                 call register_string('b4b_flag')
      if (lccsm_control_compatible) call register_string('lccsm_control_compatible')

!-----------------------------------------------------------------------
!
!  document the namelist information
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' Context:'
       write(stdout,blank_fmt)
       write(stdout,*) ' context_nml namelist settings:'
       write(stdout,blank_fmt)
       write(stdout, context_nml)
       write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   number_of_fatal_errors = 0

   if (.not. (lcoupled .eqv. lccsm)) then
     exit_string = 'FATAL ERROR: presently, lcoupled and lccsm must have the same value'
     call document ('init_context', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

#ifdef CCSMCOUPLED
   if (.not. lcoupled) then
     exit_string = 'FATAL ERROR: inconsistent options.' &
             // ' Cpp option coupled is defined, but lcoupled = .false.'
     call document ('init_context', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif
#else
   if (lcoupled) then
     exit_string = 'FATAL ERROR: inconsistent options.' &
             // ' Cpp option coupled is not defined, but lcoupled = .true.'
     call document ('init_context', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif
#endif

   if (number_of_fatal_errors > 0) &
     call exit_POP(sigAbort,'ERROR: subroutine init_context -- see preceeding message')

!-----------------------------------------------------------------------
!EOC

 end subroutine init_context

!***********************************************************************
!BOP
! !IROUTINE: init_ts
! !INTERFACE:

 subroutine init_ts(errorCode)

! !DESCRIPTION:
!  Initializes temperature and salinity and
!  initializes prognostic variables from restart if required
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  namelist input
!
!-----------------------------------------------------------------------

   integer (int_kind) :: nml_error ! namelist i/o error flag

   character (char_len) :: &
      init_ts_option,      &! option for initializing t,s
      init_ts_suboption,   &! suboption for initializing t,s (rest or spunup)
      init_ts_file,        &! filename for input T,S file
      init_ts_outfile,     &! filename for output T,S file
      init_ts_outfile_fmt   ! format for output T,S file (bin or nc)

   namelist /init_ts_nml/ init_ts_option, init_ts_file, init_ts_file_fmt, &
                          init_ts_suboption, init_ts_outfile, &
		  	  init_ts_outfile_fmt

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,icnt,        &! vertical level index
      n,                 &! tracer index
      kk,                &! indices for interpolating levitus data
      nu,                &! i/o unit for mean profile file
      iblock              ! local block address

   integer (int_kind) :: &
      m,                 &! overflows dummy loop index
      ib,ie,jb,je         ! local domain index boundaries

   integer (i4) :: &
      PHCnx,PHCny,PHCnz

   integer (i4), dimension(:,:), allocatable :: &
      PHC_msk,MASK_G

   logical (log_kind) :: &
      lccsm_branch      ,&! flag for ccsm 'ccsm_branch' restart
      lccsm_hybrid        ! flag for ccsm 'ccsm_hybrid' restart

   type (block) ::       &
      this_block          ! block information for current block

   real (r8) ::          &
      sinterp,           &! factor for interpolating levitus data
      dpth_meters         ! depth of level in meters

   real (r8), dimension(km) :: &
      tinit, sinit        ! mean initial state as function of depth

   type (datafile) ::    &
      in_file             ! data file type for init ts file

   type (io_field_desc) :: &
      io_temp, io_salt    ! io field descriptors for input T,S

   type (io_dim) :: &
      i_dim, j_dim, k_dim ! dimension descriptors

   real (r8), dimension(:,:,:,:), allocatable :: &
      TEMP_DATA           ! temp array for reading T,S data

   real (r8), dimension(:,:,:,:), allocatable :: &
      PHC_ktop,PHC_kbot     ! temp array for ncreading 3D PHC T,S data

   real (r8), dimension(:,:), allocatable :: &
      dataSrc,dataDst,     &! temp arrays for remapping PHC T,S data
      PHC_x,               &! lon array for PHC data
      PHC_y,               &! lat array for PHC data
      PHC_kmod_T,          &! PHC data on model zgrid
      PHC_kmod_S,          &! PHC data on model zgrid
      MOD_T,MOD_S,         &! 2D model init T,S data on model zgrid
      tmpkt,tmpkb,         &! 2D model init T,S data on model zgrid
      TLON_G,              &! global tlon array
      TLAT_G                ! global tlat array

   real (r8), dimension(:), allocatable :: &
      tmp1,tmp2,PHC_z             ! temp arrays

#ifdef CCSMCOUPLED
   type(shr_map_mapType)    :: PHC_map           ! used to map PHC data
#endif


   !***
   !*** 1992 Levitus mean climatology for internal generation of t,s
   !***

   real (r8), dimension(33) ::                   &
      depth_levitus = (/                         &
                   0.0_r8,   10.0_r8,   20.0_r8, &
                  30.0_r8,   50.0_r8,   75.0_r8, &
                 100.0_r8,  125.0_r8,  150.0_r8, &
                 200.0_r8,  250.0_r8,  300.0_r8, &
                 400.0_r8,  500.0_r8,  600.0_r8, &
                 700.0_r8,  800.0_r8,  900.0_r8, &
                1000.0_r8, 1100.0_r8, 1200.0_r8, &
                1300.0_r8, 1400.0_r8, 1500.0_r8, &
                1750.0_r8, 2000.0_r8, 2500.0_r8, &
                3000.0_r8, 3500.0_r8, 4000.0_r8, &
                4500.0_r8, 5000.0_r8, 5500.0_r8 /)

   real (r8), dimension(33) ::                   &
      tmean_levitus = (/                         &
                 18.27_r8, 18.22_r8, 18.09_r8,   &
                 17.87_r8, 17.17_r8, 16.11_r8,   &
                 15.07_r8, 14.12_r8, 13.29_r8,   &
                 11.87_r8, 10.78_r8,  9.94_r8,   &
                  8.53_r8,  7.35_r8,  6.38_r8,   &
                  5.65_r8,  5.06_r8,  4.57_r8,   &
                  4.13_r8,  3.80_r8,  3.51_r8,   &
                  3.26_r8,  3.05_r8,  2.86_r8,   &
                  2.47_r8,  2.19_r8,  1.78_r8,   &
                  1.49_r8,  1.26_r8,  1.05_r8,   &
                  0.91_r8,  0.87_r8,  1.00_r8 /)

   real (r8), dimension(33) ::                   &
      smean_levitus = (/                         &
                 34.57_r8, 34.67_r8, 34.73_r8,   &
                 34.79_r8, 34.89_r8, 34.97_r8,   &
                 35.01_r8, 35.03_r8, 35.03_r8,   &
                 34.98_r8, 34.92_r8, 34.86_r8,   &
                 34.76_r8, 34.68_r8, 34.63_r8,   &
                 34.60_r8, 34.59_r8, 34.60_r8,   &
                 34.61_r8, 34.63_r8, 34.65_r8,   &
                 34.66_r8, 34.68_r8, 34.70_r8,   &
                 34.72_r8, 34.74_r8, 34.75_r8,   &
                 34.74_r8, 34.74_r8, 34.73_r8,   &
                 34.73_r8, 34.72_r8, 34.72_r8 /)

!-----------------------------------------------------------------------
!
!  read input namelist and broadcast 
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   init_ts_suboption = 'rest'
   init_ts_outfile   = 'unknown_init_ts_outfile'
   init_ts_outfile_fmt = 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=init_ts_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading init_ts_nml')
   endif

   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' Initial T,S:'
       write(stdout,blank_fmt)
       write(stdout,*) ' init_ts_nml namelist settings:'
       write(stdout,blank_fmt)
       write(stdout, init_ts_nml)
       write(stdout,blank_fmt)
       if (trim(init_ts_option)    == 'ccsm_startup' .and.  &
           trim(init_ts_suboption) == 'spunup') then
                init_ts_option = 'ccsm_startup_spunup'
                luse_pointer_files = .false.
       endif
       select case (init_ts_option)
         case ('ccsm_continue', 'restart', 'ccsm_branch', 'ccsm_hybrid') 
            if (luse_pointer_files) then
               write(stdout,*) ' In this case, the init_ts_file' /&
               &/ ' name will be read from the pointer file.'   
               write(stdout,*) ' '
            endif
       end select
       call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif

   call broadcast_scalar(init_ts_option    , master_task)
   call broadcast_scalar(init_ts_suboption , master_task)
   call broadcast_scalar(luse_pointer_files, master_task)
   call broadcast_scalar(init_ts_file      , master_task)
   call broadcast_scalar(init_ts_file_fmt  , master_task)
   call broadcast_scalar(init_ts_outfile      , master_task)
   call broadcast_scalar(init_ts_outfile_fmt  , master_task)

!-----------------------------------------------------------------------
!
!  initialize t,s or call restart based on init_ts_option
!
!-----------------------------------------------------------------------


   select case (init_ts_option)

!-----------------------------------------------------------------------
!
!  set initial state from restart file
!
!-----------------------------------------------------------------------

   case ('ccsm_continue', 'restart')
      first_step   = .false.
      lccsm_branch = .false.
      lccsm_hybrid = .false.
      if (my_task == master_task .and. .not. luse_pointer_files) then
         write(stdout,'(a35,a)') 'Initial T,S read from restart file:',&
                                  trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
      call read_restart(init_ts_file,lccsm_branch,lccsm_hybrid, &
                        init_ts_file_fmt, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_ts: error in read_restart')
         return
      endif

      ltavg_restart = .true.

!-----------------------------------------------------------------------
!
!  set initial state from restart file
!
!-----------------------------------------------------------------------

   case ('ccsm_branch')
      first_step   = .false.
      lccsm_branch = .true.
      lccsm_hybrid = .false.
      if (my_task == master_task .and. .not. luse_pointer_files) then
         write(stdout,'(a40,a)') &
            'Initial T,S is a ccsm branch starting from the restart file:', &
            trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

      call read_restart(init_ts_file,lccsm_branch,lccsm_hybrid, &
                        init_ts_file_fmt, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_ts: error in read_restart')
         return
      endif

      ltavg_restart = .false.

   case ('ccsm_hybrid', 'branch')  ! ccsm hybrid start or LANL branch start
      first_step   = .false.
      lccsm_branch = .false.
      lccsm_hybrid = .true.
      if (my_task == master_task .and. .not. luse_pointer_files) then
         write(stdout,'(a80,a)') &
            'Initial T,S ccsm_hybrid start from restart file:', &
            trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
      call read_restart(init_ts_file,lccsm_branch,lccsm_hybrid, &
                        init_ts_file_fmt, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_ts: error in read_restart')
         return
      endif

      ltavg_restart = .false.

   case ('ccsm_startup_spunup')
      if(my_task == master_task )  then
         write(stdout,*) ' ccsm_startup_spunup option'
         write(stdout,*) ' init_ts_option = ', init_ts_option
      endif
      first_step   = .false.
      lccsm_branch = .false.
      lccsm_hybrid = .true.
      if (my_task == master_task .and. .not. luse_pointer_files) then
         write(stdout,'(a80,a)') &
            'Initial T,S ccsm_startup run from spun-up restart file:', &
            trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
      call read_restart(init_ts_file,lccsm_branch,lccsm_hybrid, &
                        init_ts_file_fmt, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_ts: error in read_restart')
         return
      endif

      ltavg_restart = .false.
      !*** turn pointer file-creation back on
      luse_pointer_files = .true.

!-----------------------------------------------------------------------
!
!  read full 3-d t,s from input file
!
!-----------------------------------------------------------------------

   case ('ccsm_startup', 'file')
      first_step = .true.

      if (my_task == master_task) then
         write(stdout,'(a31,a)') 'Initial 3-d T,S read from file:', &
                                 trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

      allocate(TEMP_DATA(nx_block,ny_block,km,max_blocks_clinic))

      in_file = construct_file(init_ts_file_fmt,             &
                               full_name=trim(init_ts_file), &
                               record_length = rec_type_dbl, &
                               recl_words=nx_global*ny_global)
      call data_set(in_file,'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)
      k_dim = construct_io_dim('k',km)

      io_temp = construct_io_field('TEMPERATURE', &
                 dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                 field_loc = field_loc_center,    &
                 field_type = field_type_scalar,  &
                 d3d_array=TEMP_DATA)
      io_salt = construct_io_field('SALINITY',    &
                 dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                 field_loc = field_loc_center,    &
                 field_type = field_type_scalar,  &
                 d3d_array=TEMP_DATA)

      call data_set(in_file,'define',io_temp)
      call data_set(in_file,'define',io_salt)

      call data_set(in_file,'read'  ,io_temp)
      do iblock=1,nblocks_clinic
         TRACER(:,:,:,1,curtime,iblock) = TEMP_DATA(:,:,:,iblock)
      end do
      call data_set(in_file,'read'  ,io_salt)
      do iblock=1,nblocks_clinic
         TRACER(:,:,:,2,curtime,iblock) = TEMP_DATA(:,:,:,iblock)
      end do

      call destroy_io_field(io_temp)
      call destroy_io_field(io_salt)

      deallocate(TEMP_DATA)

      call data_set(in_file,'close')
      call destroy_file(in_file)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a12,a)') ' file read: ', trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

!#################### temporary kludge for overflows ####################
!-----------------------------------------------------------------------
!   fill any overflow-deepened points with T,S values from above
!-----------------------------------------------------------------------


      if (overflows_on) then
         ! fill any overflow-deepened points with T,S values from above
         ! fill entire TRACER array for ghost (or halo) points
         do iblock = 1,nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)
            do j=1,ny_block
            do i=1,nx_block
                do n=1,num_ovf
                   do m=1,ovf(n)%num_kmt
                      if( ovf(n)%loc_kmt(m)%i.eq.this_block%i_glob(i).and.&
                          ovf(n)%loc_kmt(m)%j.eq.this_block%j_glob(j) ) then
                         if(ovf(n)%loc_kmt(m)%knew .gt. ovf(n)%loc_kmt(m)%korg) then
                            do k=ovf(n)%loc_kmt(m)%korg+1,ovf(n)%loc_kmt(m)%knew
                               ! use T,S from level above with slight increase in S
                               TRACER(i,j,k,1,curtime,iblock) = &
                                  TRACER(i,j,k-1,1,curtime,iblock)
                               TRACER(i,j,k,2,curtime,iblock) = &
                                  TRACER(i,j,k-1,2,curtime,iblock) * 1.001
                               write(stdout,100) ovf(n)%loc_kmt(m)%i, &
                                  ovf(n)%loc_kmt(m)%j,ovf(n)%loc_kmt(m)%korg, &
                                  k,ovf(n)%loc_kmt(m)%knew
                               100 format(' init_ts: T,S extended from ijKMT = ', &
                                  3(i4,1x),' to k=',i3,' until KMT_new=',i3)
                            enddo
                         endif
                      endif
                   end do
                end do
            enddo
            enddo
         enddo
      endif

      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
!################ end temporary kludge for overflows ####################


      !$OMP PARALLEL DO PRIVATE(iblock, k, n)
      do iblock = 1,nblocks_clinic
         do n=1,nt
         do k=1,km
             where (k > KMT(:,:,iblock)) &
                TRACER(:,:,k,n,curtime,iblock) = c0
         end do
         end do
 
         !*** convert salinity to model units
         TRACER(:,:,:,2,curtime,iblock) = &
         TRACER(:,:,:,2,curtime,iblock)*ppt_to_salt
      end do
      !$OMP END PARALLEL DO

      if (n_topo_smooth > 0) then
         do k=1,km
            call fill_points(k,TRACER(:,:,k,1,curtime,:),errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'init_ts: error in fill_points for temp')
               return
            endif

            call fill_points(k,TRACER(:,:,k,2,curtime,:),errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'init_ts: error in fill_points for salt')
               return
            endif

         enddo
      endif

      do iblock=1,nblocks_clinic
         TRACER(:,:,:,:,newtime,iblock) = TRACER(:,:,:,:,curtime,iblock)
         TRACER(:,:,:,:,oldtime,iblock) = TRACER(:,:,:,:,curtime,iblock)
      end do

!-----------------------------------------------------------------------
!
!  set up t,s from input mean state as function of depth
!
!-----------------------------------------------------------------------

   case ('mean')
      first_step = .true.

      !***
      !*** open input file and read t,s profile
      !***

      call get_unit(nu)
      if (my_task == master_task) then 
         write(stdout,'(a40,a)') &
            'Initial mean T,S profile read from file:', &
            trim(init_ts_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
         open(nu, file=init_ts_file, status='old')
         do k = 1,km
            read(nu,*) tinit(k),sinit(k)
         enddo
         close (nu)
      endif
      call release_unit(nu)

      call broadcast_array(tinit, master_task)
      call broadcast_array(sinit, master_task)

      !***
      !*** fill tracer array with appropriate values
      !***

      !$OMP PARALLEL DO PRIVATE(iblock, k)
      do iblock = 1,nblocks_clinic
         do k=1,km
            where (k <= KMT(:,:,iblock))
               TRACER(:,:,k,1,curtime,iblock) = tinit(k)
               TRACER(:,:,k,2,curtime,iblock) = sinit(k)*ppt_to_salt
            endwhere
         enddo

         TRACER(:,:,:,:,newtime,iblock)=TRACER(:,:,:,:,curtime,iblock)
         TRACER(:,:,:,:,oldtime,iblock)=TRACER(:,:,:,:,curtime,iblock)
      end do
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  set up initial profile from 1992 Levitus mean ocean data
!
!-----------------------------------------------------------------------

   case ('internal')
      first_step = .true.
      if (my_task == master_task) then
         write(stdout,'(a63)') & 
        'Initial T,S profile computed internally from 1992 Levitus data'
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

      !$OMP PARALLEL DO PRIVATE(iblock, k, kk, &
      !$OMP                     dpth_meters, sinterp, tinit, sinit)

      do iblock = 1,nblocks_clinic
         do k=1,km

            dpth_meters = zt(k)*mpercm

            intrp_loop: do kk=1,32
               if (dpth_meters >= depth_levitus(kk) .and. &
                   dpth_meters <  depth_levitus(kk+1)) exit intrp_loop
            end do intrp_loop

            sinterp = (dpth_meters         - depth_levitus(kk))/ &
                      (depth_levitus(kk+1) - depth_levitus(kk))
 
            tinit(k) = (c1 - sinterp)*tmean_levitus(kk) + &
                             sinterp *tmean_levitus(kk+1)
            sinit(k) = (c1 - sinterp)*smean_levitus(kk) + &
                             sinterp *smean_levitus(kk+1)

            where (k <= KMT(:,:,iblock))
               TRACER(:,:,k,1,curtime,iblock) = tinit(k)
               TRACER(:,:,k,2,curtime,iblock) = sinit(k)*ppt_to_salt
            endwhere

         enddo

         TRACER(:,:,:,:,newtime,iblock)=TRACER(:,:,:,:,curtime,iblock)
         TRACER(:,:,:,:,oldtime,iblock)=TRACER(:,:,:,:,curtime,iblock)
      enddo
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  remap PHC Levitus data to POP grid
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   case ('PHC')
      first_step = .true.

      if (my_task == master_task) then
         write(stdout,'(a63)') &
           'Initial T,S profile generated by 3D remapping of filled Levitus-PHC data'
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)

         write(stdout,*) ' init_ts_option = PHC'
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif


      allocate (TLON_G(nx_global,ny_global),TLAT_G(nx_global,ny_global))
      call gather_global(TLON_G, TLON, master_task,distrb_clinic)
      call gather_global(TLAT_G, TLAT, master_task,distrb_clinic)

      if (my_task == master_task) then

       call shr_ncread_varDimSizes(trim(init_ts_file),'TEMP',PHCnx,PHCny,PHCnz)

       allocate(MOD_T(nx_global,ny_global),MOD_S(nx_global,ny_global))
       allocate(PHC_kmod_T(PHCnx,PHCny),PHC_kmod_S(PHCnx,PHCny))
       allocate(PHC_ktop(PHCnx,PHCny,1,1),PHC_kbot(PHCnx,PHCny,1,1))
       allocate(tmpkt(PHCnx,PHCny),tmpkb(PHCnx,PHCny))
       allocate(PHC_x(PHCnx,PHCny),PHC_y(PHCnx,PHCny),PHC_msk(PHCnx,PHCny))
       allocate(tmp1(PHCnx),tmp2(PHCny),PHC_z(PHCnz))
       allocate (MASK_G(nx_global,ny_global))
       allocate(dataSrc(2,PHCnx*PHCny))
       allocate(dataDst(2,nx_global*ny_global))

       PHC_msk(:,:) = 1
       MASK_G(:,:) = 1

       call shr_ncread_tField(trim(init_ts_file),1,'lon',tmp1)
       call shr_ncread_tField(trim(init_ts_file),1,'lat',tmp2)
       call shr_ncread_tField(trim(init_ts_file),1,'depth',PHC_z)

       do j=1,PHCny 
	  PHC_x(:,j) = tmp1/radian
       enddo
       do i=1,PHCnx 
	  PHC_y(i,:) = tmp2/radian
       enddo

       call shr_map_mapSet(PHC_map, PHC_x, PHC_y, PHC_msk, &
         &                 TLON_G,  TLAT_G, MASK_G, &
         &                   name='phc_map',type='remap',algo='bilinear', &
         &                   mask='dstmask',vect='scalar')

       !-------------------------------------------------
       ! copy input data to arrays ordered for mapping
       !-------------------------------------------------

      endif

      do k=1,km
        if (my_task == master_task) then
           dpth_meters = zt(k)*mpercm

           PHC_z_loop: do kk=1,PHCnz-1
            if (dpth_meters >= PHC_z(kk) .and. &
               dpth_meters <  PHC_z(kk+1)) exit PHC_z_loop
           end do PHC_z_loop

           sinterp = (dpth_meters - depth_levitus(kk))/ &
                  (depth_levitus(kk+1) - depth_levitus(kk))

       !-------------------------------------------------
       ! do vertical remap of T
       !-------------------------------------------------
           call shr_ncread_field4dG(trim(init_ts_file),'TEMP', &
		rfld=PHC_ktop, dim3='depth',dim3i=kk)
           call shr_ncread_field4dG(trim(init_ts_file),'TEMP', &
		rfld=PHC_kbot, dim3='depth',dim3i=kk+1)
           tmpkt = reshape(PHC_ktop,(/PHCnx,PHCny/))
           tmpkb = reshape(PHC_kbot,(/PHCnx,PHCny/))
           PHC_kmod_T(:,:) = (c1 - sinterp)*tmpkt(:,:) + &
                sinterp *tmpkb(:,:)

       !-------------------------------------------------
       ! do vertical remap of S
       !-------------------------------------------------
           call shr_ncread_field4dG(trim(init_ts_file),'SALT', &
                rfld=PHC_ktop, dim3='depth',dim3i=kk)
           call shr_ncread_field4dG(trim(init_ts_file),'SALT', &
                rfld=PHC_kbot, dim3='depth',dim3i=kk+1)
           tmpkt = reshape(PHC_ktop,(/PHCnx,PHCny/))
           tmpkb = reshape(PHC_kbot,(/PHCnx,PHCny/))

           PHC_kmod_S(:,:) = (c1 - sinterp)*tmpkt(:,:) + &
                sinterp *tmpkb(:,:)

       !-------------------------------------------------
       ! do horizontal remap of T & S
       !-------------------------------------------------
           icnt = 0
           do j=1,PHCny
           do i=1,PHCnx
              icnt = icnt + 1
              dataSrc(1,icnt) = PHC_kmod_T(i,j)
              dataSrc(2,icnt) = PHC_kmod_S(i,j)
           enddo
           enddo

           call shr_map_mapData(dataSrc, dataDst, PHC_map)

           icnt = 0
           do j=1,ny_global
           do i=1,nx_global
              icnt = icnt + 1
              MOD_T(i,j) = dataDst(1,icnt)
              MOD_S(i,j) = dataDst(2,icnt)
           enddo
           enddo

        endif

        call scatter_global(TRACER(:,:,k,1,curtime,:), MOD_T, &
                master_task, distrb_clinic, field_loc_center, field_type_scalar)
        call scatter_global(TRACER(:,:,k,2,curtime,:), MOD_S, &
                master_task, distrb_clinic, field_loc_center, field_type_scalar)

      enddo

      deallocate(TLON_G,TLAT_G)
      if (my_task == master_task) then
        deallocate(MOD_T, MOD_S)
        deallocate(PHC_kmod_T,PHC_kmod_S)
        deallocate(PHC_ktop,PHC_kbot)
        deallocate(tmpkt,tmpkb,PHC_z)
        deallocate(PHC_x, PHC_y, PHC_msk, tmp1, tmp2)
        deallocate(MASK_G, dataSrc, dataDst)
      endif

      !$OMP PARALLEL DO PRIVATE(iblock, k, n)
      do iblock = 1,nblocks_clinic
         do n=1,nt
         do k=1,km
             where (k > KMT(:,:,iblock)) &
                TRACER(:,:,k,n,curtime,iblock) = c0
         end do
         end do

         !*** convert salinity to model units
         TRACER(:,:,:,2,curtime,iblock) = &
         TRACER(:,:,:,2,curtime,iblock)*ppt_to_salt
      end do
      !$OMP END PARALLEL DO

      do iblock=1,nblocks_clinic
         TRACER(:,:,:,:,newtime,iblock) = TRACER(:,:,:,:,curtime,iblock)
         TRACER(:,:,:,:,oldtime,iblock) = TRACER(:,:,:,:,curtime,iblock)
      end do

      if (trim(init_ts_outfile) /= 'unknown_init_ts_outfile') then
        if (my_task == master_task) then
          write(stdout,*) 'remapped initial T & S written to',trim(init_ts_outfile)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
        endif
        call write_init_ts(trim(init_ts_outfile),trim(init_ts_outfile_fmt))
      endif

#endif
!-----------------------------------------------------------------------
!
!  bad initialization option
!
!-----------------------------------------------------------------------

   case default
      call exit_POP(sigAbort,'Unknown t,s initialization option')
   end select

!-----------------------------------------------------------------------
!
!  check for appropriate initialization when overflows on and interactive
!
!-----------------------------------------------------------------------

   select case (init_ts_option)
   case ('mean', 'internal', 'PHC')
      if( overflows_on .and. overflows_interactive ) then
         write(stdout,*) &
           'init_ts: ERROR initializing for interactive overflows'
         write(stdout,*) &
           'initialization must be either ccsm_startup or file'
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
         call exit_POP(sigAbort,'ERROR wrong initialization with overflows')
      endif
   end select


!-----------------------------------------------------------------------
!
!  calculate RHO from TRACER at time levels curtime and oldtime
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock, k, this_block)
   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      do k=1,km
         call state(k,k,TRACER(:,:,k,1,curtime,iblock), &
                        TRACER(:,:,k,2,curtime,iblock), &
                        this_block,                     &
                        RHOOUT=RHO(:,:,k,curtime,iblock))
         call state(k,k,TRACER(:,:,k,1,oldtime,iblock), &
                        TRACER(:,:,k,2,oldtime,iblock), &
                        this_block,                     &
                        RHOOUT=RHO(:,:,k,oldtime,iblock))
      enddo

   enddo ! block loop
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  register init_ts
!
!-----------------------------------------------------------------------
   call register_string('init_ts')

   call flushm (stdout)
!-----------------------------------------------------------------------
!EOC


 end subroutine init_ts

!***********************************************************************
!BOP
! !IROUTINE: document_constants
! !INTERFACE:

 subroutine document_constants

! !DESCRIPTION:
! This routine writes the values of POP model constants to the output log file

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

   if (my_task == master_task) then
     write(stdout,ndelim_fmt)
     write(stdout,blank_fmt)
     write(stdout,*) ' Constants used in this run:'
     write(stdout,blank_fmt)

     write(stdout,1020) 'grav',    grav,      'cm/s^2'  
     write(stdout,1020) 'omega',   omega,     'rad/s'
     write(stdout,1020) 'radius',  radius,    'cm'
     write(stdout,1020) 'cp_sw',   cp_sw,     'erg/g/K'
     write(stdout,1020) 'cp_air',  cp_air,    'J/kg/K'
     write(stdout,1020) 'rho_air', rho_air,   'kg/m^3'
     write(stdout,1020) 'rho_sw',  rho_sw,    'g/cm^3'
     write(stdout,1020) 'rho_fw',  rho_fw,    'g/cm^3'
     write(stdout,1020) 'sound',   sound,     'cm/s' 
     write(stdout,1020) 'vonkar',  vonkar,    ' '        
     write(stdout,1020) 'emissivity',emissivity, ' '
     write(stdout,1020) 'stefan_boltzmann', stefan_boltzmann,  &
                           'W/m^2/K^4'
     write(stdout,1020) 'latent_heat_vapor_mks',latent_heat_vapor_mks,  &
                           'J/kg'
     write(stdout,1020) 'latent_heat_fusion',latent_heat_fusion, &
                           'erg/g'
     write(stdout,1020) 'ocn_ref_salinity', ocn_ref_salinity, 'psu' 
     write(stdout,1020) 'sea_ice_salinity', sea_ice_salinity, 'psu' 
     write(stdout,1020) 'T0_Kelvin',        T0_Kelvin,        'K' 
     write(stdout,1020) 'pi',               pi,               ' ' 

     write(stdout,blank_fmt)
     write(stdout,ndelim_fmt)
     write(stdout,blank_fmt)

     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)

   endif

1020 format (5x, a20, ' = ', 1pe25.15, 2x, a)

!-----------------------------------------------------------------------
!EOC

 end subroutine document_constants
 

!***********************************************************************
!BOP
! !IROUTINE: POP_check
! !INTERFACE:

 subroutine POP_check

! !DESCRIPTION:
!  This routine tests for consistency between model options, usually involving
!  two or more modules, then writes warning and error messages to the output log file.
!  If one or more error conditions are detected, the pop model will be shut down 
!  after all warning and error messages are printed.

! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len)      ::  &
      string                      ! temporary string

   integer (int_kind)        ::  &
      number_of_fatal_errors,    &! counter for fatal error conditions
      number_of_warnings,        &! counter for warning messages
      n,                         &! loop index 
      ns,                        &! streams loop index
      temp_tavg_id,              &! temporary tavg_id holder
      coupled_flag                ! flag for coupled_ts 

   logical (log_kind)        ::  &
      test_condition,            &! logical test condition
      lref_val,                  &! are any tracers specifying a non-zero ref_val
      ISOP_test,                 &! temporary logical associated with ISOP
      ISOP_on                     ! are any ISOP tavg fields selected?

   character (char_len), dimension(7) ::    &! var names for diag_gm_bolus test
      strings = (/'UISOP    ' , 'VISOP    ' ,  &
                  'WISOP    ' ,                &
                  'ADVT_ISOP' , 'ADVS_ISOP' ,  &
                  'VNT_ISOP ' , 'VNS_ISOP '   /)


   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*)' POP_check: Check for Option Inconsistencies'
      write(stdout,blank_fmt)
   endif


   !====================!
   ! warning conditions !
   !====================!

   number_of_warnings = 0

!-----------------------------------------------------------------------
!
!  'varthick' and dtuxcel /= dttxcel(1)
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
     if (sfc_layer_type == sfc_layer_varthick .and.  &
         dtuxcel /= dttxcel(1) ) then
       exit_string = 'WARNING:  Surface tracer and momentum timesteps are unequal; ' /&     
       &/'may cause instability when using variable-thickness surface layer.'
       call document ('POP_check', exit_string)
       number_of_warnings = number_of_warnings + 1
     endif
   endif
 
!-----------------------------------------------------------------------
!
!  bulk-NCEP and marginal-seas balancing
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
     if (sfwf_formulation == 'bulk-NCEP' .and. lms_balance) then
       exit_string = 'WARNING:  runoff and marginal seas balancing cannot ' /&
       &/ 'be used with the bulk-NCEP option'
       call document ('POP_check', exit_string)
       number_of_warnings = number_of_warnings + 1
     endif
   endif 

!-----------------------------------------------------------------------
!
!  are time-averaging and coupling frequencies compatible?
!
!-----------------------------------------------------------------------

   coupled_flag = get_time_flag_id('coupled_ts')

   if (my_task == master_task) then
     do ns=1,n_tavg_streams
     if (check_time_flag_int(tavg_streams(ns)%field_flag,   freq_opt=.true.) > 0 .and. &
         check_time_flag_int(coupled_flag,freq_opt=.true.) > 0) then

       if (check_time_flag_int(tavg_streams(ns)%field_flag,   freq_opt=.true.) /=  &
           check_time_flag_int(coupled_flag,freq_opt=.true.)) then
           exit_string = 'WARNING:  time-averaging and coupling frequency ' /&
           &/ 'may be incompatible; tavg must be integer multiple of coupling freq'
           call document ('POP_check', exit_string)
           number_of_warnings = number_of_warnings + 1
       else
         if ( mod(check_time_flag_int(tavg_streams(ns)%field_flag,   freq=.true.),  &
                  check_time_flag_int(coupled_flag,freq=.true.)) .ne. 0) then
           exit_string = 'WARNING: time-averaging frequency is incompatible with ' /&
           &/ ' the coupling frequency'
           call document ('POP_check', exit_string)
           number_of_warnings = number_of_warnings + 1
         endif
       endif
     endif
     enddo ! ns
   endif


!-----------------------------------------------------------------------
!
!  Wrap up warning section with message
!
!-----------------------------------------------------------------------

   call broadcast_scalar(number_of_warnings, master_task)
 
   if (number_of_warnings == 0 ) then
     if (my_task == master_task) then
       exit_string = 'No warning messages generated'
       call document ('POP_check', exit_string)
     endif
   endif


   !========================!
   ! fatal error conditions !
   !========================!


   number_of_fatal_errors = 0

!-----------------------------------------------------------------------
!
!  tidal mixing without KPP mixing
!
!-----------------------------------------------------------------------

   if (check_all(ltidal_mixing .and. vmix_itype /= vmix_type_kpp)) then
     exit_string =   &
     'FATAL ERROR:  Tidally driven mixing is only allowed when KPP mixing is enabled' 
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  tidal mixing without bckgrnd_vdc2 = 0.0
!
!-----------------------------------------------------------------------

   if (check_all(ltidal_mixing .and. bckgrnd_vdc2 /= c0)) then
     exit_string =   &
    'FATAL ERROR:  bckgrnd_vdc2 must be zero when tidal_mixing option is enabled'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  diag_gm_bolus = .true., but ISOP variables not activated in tavg_contents file
!
!-----------------------------------------------------------------------

   if (registry_match('diag_gm_bolus') .and. my_task == master_task) then
    ISOP_on = .true.
    exit_string = 'FATAL ERROR: '

    do n=1,7
      ISOP_test = .false. 
      string = trim(strings(n))
      ISOP_test = set_in_tavg_contents (tavg_id(trim(string),quiet=.true.))
      if (.not. ISOP_test) then
        exit_string =  trim(exit_string) // ' ' // trim(string)
        ISOP_on = .false.
      endif
    enddo
    
    if (.not. ISOP_on) then
     exit_string =   trim(exit_string) /&
     &/' must be activated in tavg_contents file when diag_gm_bolus = .T.'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
    endif
  endif ! diag_gm_bolus

!-----------------------------------------------------------------------
!
!  luse_cpl_ifrac is true, but OCN_WGT is not allocated
!
!-----------------------------------------------------------------------

   if (check_all(luse_cpl_ifrac .and. .not. allocated(OCN_WGT))) then
     exit_string =   &
     'FATAL ERROR:  cannot set luse_cpl_ifrac .true. without allocating OCN_WGT'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  ecosystem requests diurnal cycle, but base model is not using it
!
!-----------------------------------------------------------------------

   if (ecosys_on) then

     if ((.not. ecosys_qsw_distrb_const) .and. &
         (qsw_distrb_iopt == qsw_distrb_iopt_const)) then
       exit_string = &
       'FATAL ERROR: cannot set ecosys_qsw_distrb_const=.false. unless qsw_distrb_opt/=const'
       call document ('POP_check', exit_string)
       number_of_fatal_errors = number_of_fatal_errors + 1
     endif

   endif

!-----------------------------------------------------------------------
!
!  untested forcing_coupled option
!
!-----------------------------------------------------------------------

   if (sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx) then
     exit_string =  'FATAL ERROR:  untested/unsupported combination of options'
     exit_string =   trim(exit_string) /&
     &/' (sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx)'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  inertial mixing on
!
!-----------------------------------------------------------------------

  if (.not. lniw_mixing .and. linertial) then
     exit_string =  'FATAL ERROR:  inertial mixing option. '
     exit_string =   trim(exit_string) /&
     &/' This option is untested. DO NOT USE!'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  inertial mixing inconsistencies
!
!-----------------------------------------------------------------------

  if (linertial .and. (.not. registry_match('diag_gm_bolus') .or. partial_bottom_cells)) then
     exit_string =  'FATAL ERROR:  inertial mixing option inconsistency. '
     exit_string =   trim(exit_string) /&
     &/' diag_gm_bolus must be on and partial_bottom_cells must not be on'
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  overflow check: if overflow active, horiz_grid_opt must be 'file'
!
!-----------------------------------------------------------------------

   if ( overflows_on ) then
     if ( overflows_interactive .and. .not. registry_match('topography_opt_file') ) then
        exit_string = 'FATAL ERROR: interactive overflows without topography option = file'
        call document ('POP_check', exit_string)
        number_of_fatal_errors = number_of_fatal_errors + 1
     endif
   endif

!-----------------------------------------------------------------------
!
!  overflow check: if overflow active, must set state_range_iopt  = state_range_enforce
!     and state_itype = state_type_mwjf for consistency
!
!-----------------------------------------------------------------------

   if ( overflows_on ) then
     if ( .not. (state_range_iopt == state_range_enforce .and. state_itype == state_type_mwjf) ) then
        exit_string = 'FATAL ERROR: if overflows are active, must have state_range_opt = enforce '/&
         &/' and state_choice = mwjf for consistency. You can uncomment this and procede at your own risk.'
        call document ('POP_check', exit_string)
        number_of_fatal_errors = number_of_fatal_errors + 1
     endif
   endif

!-----------------------------------------------------------------------
!
!  near-inertial wave mixing without KPP mixing
!
!-----------------------------------------------------------------------

   if (check_all(lniw_mixing .and. vmix_itype /= vmix_type_kpp)) then
     exit_string =   &
     'FATAL ERROR:  Near-inertial wave mixing is only allowed when KPP mixing is enabled' 
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  near-inertial wave mixing and not 2-hour coupling
!
!-----------------------------------------------------------------------

   test_condition = (coupled_freq_iopt == freq_opt_nhour .and. ncouple_per_day == 12)
   if (check_all(lniw_mixing .and. .not. test_condition)  ) then
     call document ('POP_check', 'coupled_freq_iopt                   ', coupled_freq_iopt )
     call document ('POP_check', 'freq_opt_nhour                      ', freq_opt_nhour )
     call document ('POP_check', 'ncouple_per_day                     ', ncouple_per_day )
     call document ('POP_check', '(coupled_freq_iopt == freq_opt_nhour .and. ncouple_per_day == 2) ',  &
                        (coupled_freq_iopt == freq_opt_nhour .and. ncouple_per_day == 2) )
     exit_string =   &
     'FATAL ERROR:  Near-inertial wave mixing is only allowed when coupling every two hours' 
     call document ('POP_check', exit_string)
     number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  Now that error messages have been written, stop if there are fatal errors
!
!-----------------------------------------------------------------------

   call broadcast_scalar(number_of_fatal_errors, master_task)
 
   if (number_of_fatal_errors > 0 ) then
     call exit_POP (sigAbort,  &
         'correct the error condition(s) listed above before continuing')
   else
     if (my_task == master_task) then
       exit_string = 'No fatal error conditions detected'
       call document ('POP_check', exit_string)
     endif
   endif
 

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_check


!***********************************************************************
!BOP
! !IROUTINE: write_init_ts
! !INTERFACE:

 subroutine write_init_ts(outfile, outfile_fmt)

! !DESCRIPTION:
!  This routine writes out initial TEMP and SALT mapped to
!  POP grid for topography_opt='bathymetry'
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      outfile,                  &! input file name (with path)
      outfile_fmt                ! input file format (bin or nc)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (datafile) ::     &
      ts_file           ! io file type for viscosity file

   type (io_field_desc) ::     &
      TEMP_d, SALT_d      ! descriptors for temp and salt fields

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for vertical levels

!-----------------------------------------------------------------------
!
!  construct io file type and open for writing
!
!-----------------------------------------------------------------------

   ts_file =  construct_file(outfile_fmt, root_name=outfile,    &
                               record_length=rec_type_dbl,       &
                               recl_words=nx_global*ny_global)

   call data_set(ts_file, 'open')

!-----------------------------------------------------------------------
!
!  define variables to be written
!
!-----------------------------------------------------------------------

   !*** define dimensions

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', km)

   TEMP_d = construct_io_field('TEMP', dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                   long_name='Potential Temperature',           &
                   units    ='degC',                                 &
                   grid_loc ='3111',                                  &
                   field_loc = field_loc_center,                    &
                   field_type = field_type_scalar,                    &
                   d3d_array = TRACER(:,:,:,1,curtime,:))

   SALT_d = construct_io_field('SALT', dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                   long_name='Salinity',      &
                   units    ='gram/kilogram',                                 &
                   grid_loc ='3111',                                  &
                   field_loc = field_loc_center,                    &
                   field_type = field_type_scalar,                    &
                   d3d_array = TRACER(:,:,:,2,curtime,:))

   call data_set (ts_file, 'define', TEMP_d)
   call data_set (ts_file, 'define', SALT_d)

!-----------------------------------------------------------------------
!
!  write arrays then clean up
!
!-----------------------------------------------------------------------

   call data_set (ts_file, 'write', TEMP_d)
   call data_set (ts_file, 'write', SALT_d)

   call destroy_io_field (TEMP_d)
   call destroy_io_field (SALT_d)

   call data_set (ts_file, 'close')
   call destroy_file(ts_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine write_init_ts


!***********************************************************************
!BOP
! !IROUTINE: check_all
! !INTERFACE:

 function check_all(condition)

! !DESCRIPTION:
!   Tests input logical condition on all processors; if any element is 
!     .true., check_all is set to .true.
! 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in) :: &
      condition               ! logical condition to be checked

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      check_all               ! true if condition is true on any processor

                     
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      int_condition              ! integer form of logical input condition


   if (condition) then
     int_condition = 1
   else
     int_condition = 0
   endif

   check_all =  (global_sum(int_condition,distrb_clinic) > 0)

!-----------------------------------------------------------------------
!EOC

 end function check_all

!***********************************************************************

 end module initial

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
