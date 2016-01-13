!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_front_end.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module cism_front_end
  ! The CISM front-end is used to connect both the standalone driver
  ! (cism_driver) or the CISM interface to CESM (cism_cesm_interface),
  ! to the internal and external dycore interface programs.  These are
  !* cism_internal_dycore_interface and cism_external_dycore_interface.

contains

subroutine cism_init_dycore(model)

  use parallel
  use glimmer_global
  use glide
  use glissade
  use eismint_forcing
  use glimmer_log
  use glimmer_config
  use glide_nc_custom, only: glide_nc_fillall
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glide_io, only: glide_io_writeall

  use cism_external_dycore_interface

!  use glimmer_to_dycore

  use glide_stop, only: glide_finalise
  use glide_diagnostics

  implicit none


  type(glide_global_type) :: model        ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) :: time                   ! model time in years
  integer :: clock,clock_rate

  integer*4 external_dycore_model_index

  integer :: wd
  logical :: do_glide_init

  integer :: tstep_count

  !  print *,'Entering cism_init_dycore'


  !TODO - call this only for parallel runs?
  ! call parallel_initialise     

  call glimmer_GetCommandline()

  ! DMR -- open_log call commented out, since called in gci_init_interface()
  ! start logging
  ! call open_log(unit=50, fname=logname(commandline_configname))
  
  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  call ConfigRead(commandline_configname,config)

#if (! defined CCSMCOUPLED && ! defined CESMTIMERS)
  ! start timing
  call system_clock(clock,clock_rate)
  wall_start_time = real(clock,kind=dp)/real(clock_rate,kind=dp)
#else
  wall_start_time = 0.0
  wall_stop_time = 0.0
#endif

  ! initialise profiling
  call profile_init(model%profile,'glide.profile')

  call t_startf('cism')
 
  ! initialise GLIDE
  call t_startf('initialization')

  call glide_config(model,config)

  ! This call is needed only if running the EISMINT test cases
  call eismint_initialise(model%eismint_climate,config)

  wd = model%options%whichdycore 
!  do_glide_init = (wd == DYCORE_GLIDE) .OR. (wd == DYCORE_BISICLES) .OR. (wd == DYCORE_ALBANYFELIX)
  do_glide_init = (wd == DYCORE_GLIDE)

  if (do_glide_init) then
     call glide_initialise(model)
  else       ! glam/glissade dycore	
     call glissade_initialise(model)
  endif

  call CheckSections(config)
 
 ! fill dimension variables on output files
  call glide_nc_fillall(model)

  time = model%numerics%tstart
  tstep_count = 0
  model%numerics%time = time    ! MJH added 1/10/13 - the initial diagnostic glissade solve won't know 
                                !                     the correct time on a restart unless we set it here.

  ! Set EISMINT forcing for initial time
  call eismint_massbalance(model%eismint_climate,model,time)
  call eismint_surftemp(model%eismint_climate,model,time)

  ! read forcing time slice if needed - this will overwrite values from IC file if there is a conflict.
  call glide_read_forcing(model, model)

  call spinup_lithot(model)

  call t_stopf('initialization')

  if (model%options%whichdycore == DYCORE_BISICLES) then
    call t_startf('initial_diag_var_solve')

    if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
      ! disable further profiling in normal usage
      call t_adj_detailf(+10)
    endif

    call cism_init_external_dycore(model%options%external_dycore_type,model)

    if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
      ! restore profiling to normal settings
      call t_adj_detailf(-10)
    endif

    call t_stopf('initial_diag_var_solve')
  endif

  if (model%options%whichdycore .ne. DYCORE_BISICLES) then
  !MJH Created this block here to fill out initial state without needing to enter time stepping loop.  This allows
  ! a run with tend=tstart to be run without time-stepping at all.  It requires solving all diagnostic (i.e. not
  ! time depdendent) variables (most important of which is velocity) for the initial state and then writing the 
  ! initial state as time 0 (or more accurately, as time=tstart).  Also, halo updates need to occur after the 
  ! diagnostic variables are calculated.

  ! ------------- Calculate initial state and output it -----------------

    call t_startf('initial_diag_var_solve')

    select case (model%options%whichdycore)
      case (DYCORE_GLIDE)

        if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
          ! disable further profiling in normal usage
          call t_adj_detailf(+10)
        endif

        ! Don't call glide_init_state_diagnostic when running old glide
        ! Instead, start with zero velocity
        if (.not. oldglide) then
          call glide_init_state_diagnostic(model)
        endif

        if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
          ! restore profiling to normal settings
          call t_adj_detailf(-10)
        endif

      case (DYCORE_GLAM, DYCORE_GLISSADE, DYCORE_ALBANYFELIX)

        if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
          ! disable further profiling in normal usage
          call t_adj_detailf(+10)
        endif

        ! solve the remaining diagnostic variables for the initial state
        call glissade_diagnostic_variable_solve(model)  ! velocity, usrf, etc.

        if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
          ! restore profiling to normal settings
          call t_adj_detailf(-10)
        endif

        case default

    end select

    call t_stopf('initial_diag_var_solve')

    ! Write initial diagnostic output to log file

    call t_startf('initial_write_diagnostics')
    call glide_write_diagnostics(model,        time,       &
                                 tstep_count = tstep_count)
    call t_stopf('initial_write_diagnostics')

  end if ! whichdycore .ne. DYCORE_BISICLES


  ! --- Output the initial state -------------

  call t_startf('initial_io_writeall')                                                          
  call glide_io_writeall(model, model, time=time)          ! MJH The optional time argument needs to be supplied 
                                                           !     since we have not yet set model%numerics%time
                                                           !WHL - model%numerics%time is now set above
  call t_stopf('initial_io_writeall')

end subroutine cism_init_dycore


subroutine cism_run_dycore(model)

  use parallel
  use glimmer_global
  use glide
  use glissade
  use eismint_forcing
  use glimmer_log
  use glimmer_config
  use glide_nc_custom, only: glide_nc_fillall
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glide_io, only: glide_io_writeall, glide_io_writeall

  use cism_external_dycore_interface
  
  use glide_stop, only: glide_finalise
  use glide_diagnostics

  implicit none


  type(glide_global_type) :: model        ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) :: time                   ! model time in years
  real(kind=dp) :: dt                     ! current time step to use
  real(kind=dp) :: time_eps               ! tolerance within which times are equal 
  integer :: clock,clock_rate
  integer :: tstep_count

  integer*4 :: external_dycore_model_index

!  external_dycore_model_index = this_rank + 1
  external_dycore_model_index = 1

  time = model%numerics%tstart
  tstep_count = 0
  time_eps = model%numerics%tinc/1000.0d0

  ! ------------- Begin time step loop -----------------
 
  ! run an internal or external dycore, depending on setting external_dycore_type

  ! check if we're doing any evolution
  if (time < model%numerics%tend) then
    do while(time + time_eps < model%numerics%tend)

      !!! SFP moved block of code for applying time dependent forcing read in from netCDF here,
      !!! as opposed to at the end of the time step (commented it out in it's original location for now)
      !!! This is a short-term fix. See additioanl discussion as part of issue #19 (in cism-piscees github repo).

      ! Forcing from a 'forcing' data file - will read time slice if needed
      call t_startf('read_forcing')
      call glide_read_forcing(model, model)
      call t_stopf('read_forcing')

      ! Increment time step
      if (model%options%whichdycore /= DYCORE_BISICLES) then
        time = time + model%numerics%tinc
        tstep_count = tstep_count + 1
        model%numerics%time = time  ! TODO This is redundant with what is happening in glide/glissade, but this is needed for forcing to work properly.
      endif
! print *,"external_dycore_type: ",model%options%external_dycore_type

      !if (model%options%external_dycore_type .EQ. 0) then      ! NO_EXTERNAL_DYCORE) then
      !  if (model%options%whichdycore == DYCORE_GLIDE) then
      call t_startf('tstep')

      select case (model%options%whichdycore)
        case (DYCORE_GLIDE)

          call t_startf('glide_tstep_p1')
          call glide_tstep_p1(model,time)
          call t_stopf('glide_tstep_p1')

          call t_startf('glide_tstep_p2')
          call glide_tstep_p2(model)
          call t_stopf('glide_tstep_p2')

          call t_startf('glide_tstep_p3')
          call glide_tstep_p3(model)
          call t_stopf('glide_tstep_p3')

        case (DYCORE_GLAM, DYCORE_GLISSADE, DYCORE_ALBANYFELIX)
          ! glam/glissade dycore

          call glissade_tstep(model,time)

        case (DYCORE_BISICLES)
          ! print *,'Using External Dycore'
          ! The time variable gets incremented within this call:
          dt = model%numerics%tinc
        
          if (time + dt + time_eps > model%numerics%tend) then
             dt = model%numerics%tend - time
          endif
          call cism_run_external_dycore(model%options%external_dycore_model_index, &
                                        time,dt)
          ! time = time + model%numerics%tinc
        case default
      end select

      call t_stopf('tstep')
      !endif

      ! write ice sheet diagnostics to log file at desired interval (model%numerics%dt_diag)

      call t_startf('write_diagnostics')
      call glide_write_diagnostics(model,        time,       &
                                  tstep_count = tstep_count)
      call t_stopf('write_diagnostics')

      ! update time from dycore advance
      model%numerics%time = time

      ! --- Set forcing ---
      ! Setting forcing at the end of the time step maintains consistency
      ! with a forward Euler time step and ensures consistency of the time stamp
      ! to fields in input and output files.  
      ! For forward Euler time stepping we want S^n+1 = g(S^n, F^n)
      ! where S is the model state, F is forcing, and n, n+1 are time levels
      ! We also want a forcing field in the output file to have a time stamp
      ! that matches its time stamp in the input file or the EISMINT analytic function.
      ! The simplest way to ensure both of these criteria is to set forcing at the
      ! end of each time step.
      ! EISMINT forcing
      ! NOTE: these only do something when an EISMINT case is run
      call t_startf('set_forcing')
      call eismint_massbalance(model%eismint_climate,model,time)
      call eismint_surftemp(model%eismint_climate,model,time)
      call t_stopf('set_forcing')

      !!! SFP moved this next block of code to the start of the time step. See additional notes there.

  !    ! Forcing from a 'forcing' data file - will read time slice if needed
  !    call t_startf('read_forcing')
  !    call glide_read_forcing(model, model)
  !    call t_stopf('read_forcing')

      ! Write to output netCDF files at desired intervals
      call t_startf('io_writeall')
      call glide_io_writeall(model,model)
      call t_stopf('io_writeall')
    end do   ! time < model%numerics%tend
  else ! no evolution -- diagnostic run, still want to do IO
      ! (DFM) uncomment this if we want to do an I/O step even if no evoloution 
      !call t_startf('glide_io_writeall')
      !call glide_io_writeall(model,model)
      !call t_stopf('glide_io_writeall')  
  endif    
  
end subroutine cism_run_dycore

subroutine cism_finalize_dycore(model)

  use parallel
  use glimmer_global
  use glide
  use glissade
  use glimmer_log
  use glimmer_config
  use glide_nc_custom, only: glide_nc_fillall
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glide_io, only: glide_io_writeall

  use cism_external_dycore_interface
  
  use glide_stop, only: glide_finalise
  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  integer :: clock,clock_rate

  call t_stopf('cism')

  ! finalise GLIDE
  call glide_finalise(model)

  ! (DFM) -- finalize external dycores 
  if (model%options%whichdycore == DYCORE_BISICLES) then
    call t_startf('finalize_external_dycore')
    call cism_finalize_external_dycore(model%options%external_dycore_type,model)
    call t_stopf('finalize_external_dycore')
  endif



  !TODO - Do we need to call glimmer_write_stats?
#if (! defined CCSMCOUPLED && ! defined CESMTIMERS)
  call system_clock(clock,clock_rate)
  wall_stop_time = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,wall_stop_time-wall_start_time)
#endif

  call close_log

  !TODO - call this only for parallel runs?
  ! call parallel_finalise
end subroutine cism_finalize_dycore

end module cism_front_end
