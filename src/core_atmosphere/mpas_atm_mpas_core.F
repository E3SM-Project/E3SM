! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework

   type (MPAS_Clock_type), pointer :: clock


   contains


   subroutine mpas_core_init(domain, stream_manager, startTimeStamp)

      use mpas_kind_types
      use mpas_stream_manager
      use mpas_grid_types

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND), pointer :: dt
      type (block_type), pointer :: block

      character(len=StrKIND) :: timeStamp
      integer :: i
      integer :: ierr
      logical, pointer :: config_do_restart
     
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: diag
      type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime


      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock


      call mpas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)

      
      !
      ! If this is a restart run, read the restart stream, else read the input
      ! stream.
      ! Regardless of which stream we read for initial conditions, reset the
      ! input alarms for both input and restart before reading any remaining
      ! input streams.
      !
      if (config_do_restart) then
         call MPAS_stream_mgr_read(stream_manager, streamID='restart', ierr=ierr)
      else
         call MPAS_stream_mgr_read(stream_manager, streamID='input', ierr=ierr)
      end if
      if (ierr /= MPAS_STREAM_MGR_NOERR) then
         call mpas_dmpar_global_abort('Error reading initial conditions.') 
      end if
      call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='input', direction=MPAS_STREAM_INPUT, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='restart', direction=MPAS_STREAM_INPUT, ierr=ierr)

      !
      ! Read all other inputs
      ! For now we don't do this here to match results with previous code; to match requires 
      ! that we read in SST and seaice fields after the call to atm_mpas_init_block()
      !
!      call MPAS_stream_mgr_read(stream_manager, ierr=ierr)
!      call MPAS_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_INPUT, ierr=ierr) 

      if (.not. config_do_restart) then
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_initialize_time_levels(state)
            block => block % next
         end do
      end if


      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 


      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_field(state, 'u', u_field, 1)
      call mpas_dmpar_exch_halo_field(u_field)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)

         call atm_mpas_init_block(domain % dminfo, stream_manager, block, mesh, dt)

         call mpas_pool_get_array(state, 'xtime', xtime, 1)
         xtime = startTimeStamp

         block => block % next
      end do

      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_field(diag, 'pv_edge', pv_edge_field)
      call mpas_dmpar_exch_halo_field(pv_edge_field)

      call mpas_pool_get_field(diag, 'ru', ru_field)
      call mpas_dmpar_exch_halo_field(ru_field)

      call mpas_pool_get_field(diag, 'rw', rw_field)
      call mpas_dmpar_exch_halo_field(rw_field)

   end subroutine mpas_core_init


   subroutine atm_simulation_clock_init(core_clock, configs, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: local_err
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_start_time
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_run_duration
      character (len=StrKIND), pointer :: config_stop_time
      character (len=StrKIND) :: startTimeStamp


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)

      if(trim(config_start_time) == 'file') then
         open(22,file=trim(config_restart_timestamp_name),form='formatted',status='old')
         read(22,*) startTimeStamp
         close(22)
      else
        startTimeStamp = config_start_time
      end if
      call mpas_set_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=local_err)
      call mpas_set_timeInterval(timeStep, dt=config_dt, ierr=local_err)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=local_err)
      else
          write(stderrUnit,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if

      !TODO: set phyics alarms here...
      !....
      !....

   end subroutine atm_simulation_clock_init


   subroutine atm_mpas_init_block(dminfo, stream_manager, block, mesh, dt)
   
      use mpas_grid_types
      use atm_time_integration
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
      use mpas_stream_manager
#ifdef DO_PHYSICS
!     use mpas_atmphys_aquaplanet
      use mpas_atmphys_control
      use mpas_atmphys_init
      use mpas_atmphys_manager
#endif
   
      implicit none
   
      type (dm_info), intent(in) :: dminfo
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      type (block_type), intent(inout) :: block
      type (mpas_pool_type), intent(inout) :: mesh     !MGD does this need to be a pointer?
      real (kind=RKIND), intent(in) :: dt

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: tend
      type (mpas_pool_type), pointer :: sfc_input
      type (mpas_pool_type), pointer :: diag_physics
      type (mpas_pool_type), pointer :: atm_input
      
      real (kind=RKIND), dimension(:,:), pointer :: u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, meshScalingDel4
      character(len=StrKIND), pointer :: mminlu

      integer, pointer :: nEdgesSolve
   
      logical, pointer :: config_do_restart, config_do_DAcycling

   
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'state', state)

      call mpas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(block % configs, 'config_do_DAcycling', config_do_DAcycling)

      if (.not. config_do_restart .or. (config_do_restart .and. config_do_DAcycling)) then
         call atm_init_coupled_diagnostics( state, 1, diag, mesh, block % configs)
      end if
      call atm_compute_solve_diagnostics(dt, state, 1, diag, mesh, block % configs)

      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)

      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_reconstruct(mesh, u,                   &
                            uReconstructX,             &
                            uReconstructY,             &
                            uReconstructZ,             &
                            uReconstructZonal,         &
                            uReconstructMeridional     &
                           )
   
#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at least one physics
      !parameterization is called (using the logical moist_physics):
      call physics_namelist_check(mesh, block % configs)

      !proceed with initialization of physics parameterization if moist_physics is set to true:
      call mpas_pool_get_subpool(block % structs, 'sfc_input', sfc_input)

      ! Before calling physics_init, ensure that mminlu contains the name of the land use dataset
      call mpas_pool_get_array(sfc_input, 'mminlu', mminlu)
      if (len_trim(mminlu) == 0) then
            write(0,*) '****************************************************************'
            write(0,*) 'No information on land use dataset is available.'
            write(0,*) 'Assume that we are using ''USGS''.'
            write(0,*) '****************************************************************'
            write(mminlu,'(a)') 'USGS'
      end if


      if (moist_physics) then
         !initialization of some input variables in registry:
         call mpas_pool_get_subpool(block % structs, 'tend', tend)
         call mpas_pool_get_subpool(block % structs, 'diag_physics', diag_physics)
         call mpas_pool_get_subpool(block % structs, 'atm_input', atm_input)
         call physics_registry_init(mesh, block % configs, sfc_input)
         call physics_run_init(block % configs, mesh, state, clock, stream_manager)

         !initialization of all physics:
         call physics_init(dminfo, clock, block % configs, mesh, diag, tend, state, 1, diag_physics, &
                           atm_input, sfc_input)
      endif
#endif
   
      call atm_compute_mesh_scaling(mesh, block % configs)

      call atm_compute_damping_coefs(mesh, block % configs)

      call atm_compute_pgf_coefs(mesh, block % configs)

      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)

      write(0,*) 'min/max of meshScalingDel2 = ', minval(meshScalingDel2(1:nEdgesSolve)), &
                                                  maxval(meshScalingDel2(1:nEdgesSolve))
      write(0,*) 'min/max of meshScalingDel4 = ', minval(meshScalingDel4(1:nEdgesSolve)), &
                                                  maxval(meshScalingDel4(1:nEdgesSolve))

      call atm_adv_coef_compression(mesh)

   end subroutine atm_mpas_init_block
   
   
   subroutine mpas_core_run(domain, stream_manager)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_io_streams, only : MPAS_STREAM_LATEST_BEFORE
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
   
      real (kind=RKIND), pointer :: dt
      logical, pointer :: config_do_restart
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: itimestep
      integer :: ierr

      type (mpas_pool_type), pointer :: state, diag, diag_physics, mesh

      ! For high-frequency diagnostics output
      character (len=StrKIND) :: tempfilename


      ! Eventually, dt should be domain specific
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)

      ! Avoid writing a restart file at the initial time
      call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      ! Also, for restart runs, avoid writing the initial history fields to avoid overwriting those from the preceding run
      if (config_do_restart) then
         call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='output', direction=MPAS_STREAM_OUTPUT, ierr=ierr)
      end if

      if (MPAS_stream_mgr_ringing_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call atm_compute_output_diagnostics(state, 1, diag, diag_physics, mesh)

            block_ptr => block_ptr % next
         end do
      end if
      call mpas_stream_mgr_write(stream_manager, ierr=ierr)
      call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)


      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 1
      do while (.not. mpas_is_clock_stop_time(clock))

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         write(0,*) ' '
         write(0,*) 'Begin timestep ', trim(timeStamp)

         !
         ! Read external field updates
         !
         call MPAS_stream_mgr_read(stream_manager, whence=MPAS_STREAM_LATEST_BEFORE, ierr=ierr)
         if (ierr /= MPAS_STREAM_MGR_NOERR) then
            call mpas_dmpar_global_abort('Error reading one or more input files.') 
         end if
         call MPAS_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_INPUT, ierr=ierr)

         call mpas_timer_start("time integration")
         call atm_do_timestep(domain, dt, itimestep)
         call mpas_timer_stop("time integration")   

         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
         call mpas_pool_shift_time_levels(state)
         
         ! Advance clock before writing output
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)

         !
         ! Write any output streams that have alarms ringing, after computing diagnostics fields
         !
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         if (MPAS_stream_mgr_ringing_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            block_ptr => domain % blocklist
            do while (associated(block_ptr))

               call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
               call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
               call atm_compute_output_diagnostics(state, 1, diag, diag_physics, mesh)

               block_ptr => block_ptr % next
            end do
         end if
         if (MPAS_stream_mgr_ringing_alarms(stream_manager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            block_ptr => domain % blocklist
            do while (associated(block_ptr))

               call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
               call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
               call atm_compute_restart_diagnostics(state, 1, diag, mesh)

               block_ptr => block_ptr % next
            end do
         end if
         call mpas_stream_mgr_write(stream_manager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      end do
   
   end subroutine mpas_core_run
   
   
   subroutine atm_compute_output_diagnostics(state, time_lev, diag, diag_physics, mesh)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to history files
   !
   ! Input: state - contains model prognostic fields
   !        mesh  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
      use mpas_atm_interp_diagnostics
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: diag_physics
      type (mpas_pool_type), intent(in) :: mesh
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND), dimension(:,:), pointer :: theta, rho, theta_m, rho_zz, zz
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base, pressure_p, pressure
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
      call mpas_pool_get_array(diag, 'pressure', pressure)

      call mpas_pool_get_array(mesh, 'zz', zz)

      do iCell=1,nCells
         do k=1,nVertLevels
            theta(k,iCell) = theta_m(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            pressure(k,iCell) = pressure_base(k,iCell) + pressure_p(k,iCell)
         end do
      end do

      call interp_diagnostics(mesh, state, time_lev, diag, diag_physics)
   
   end subroutine atm_compute_output_diagnostics
   
   
   subroutine atm_compute_restart_diagnostics(state, time_lev, diag, mesh)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to restart files
   !
   ! Input: state - contains model prognostic fields
   !        mesh  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                 ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND), dimension(:,:), pointer :: theta, rho, theta_m, rho_zz, zz
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)

      call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, time_lev)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)

      call mpas_pool_get_array(mesh, 'zz', zz)

      do iCell=1,nCells
         do k=1,nVertLevels
            theta(k,iCell) = theta_m(k,iCell) / (1.0_RKIND + rvord * scalars(index_qv,k,iCell))
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
         end do
      end do
   
   end subroutine atm_compute_restart_diagnostics


   subroutine atm_do_timestep(domain, dt, itimestep)
   
      use mpas_grid_types
      use mpas_kind_types
      use atm_time_integration
#ifdef DO_PHYSICS
      use mpas_atmphys_control
      use mpas_atmphys_driver
      use mpas_atmphys_manager
      use mpas_atmphys_update
#endif
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: itimestep
      
      type (MPAS_Time_Type) :: startTime, currTime
      type (MPAS_TimeInterval_Type) :: xtimeTime
      character(len=StrKIND) :: timeStamp
      integer :: s, s_n, s_d
      real (kind=RKIND) :: xtime_s
      integer :: ierr

      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         
      xtimeTime = currTime - startTime
      call mpas_get_timeInterval(interval=xtimeTime, S=s, S_n=s_n, S_d=s_d, ierr=ierr)         
      xtime_s = (s + s_n / s_d)

      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         


#ifdef DO_PHYSICS
      !proceed with physics if moist_physics is set to true:
      if(moist_physics) then
         call physics_timetracker(domain,dt,clock,itimestep,xtime_s)
         call physics_driver(domain,itimestep,xtime_s)
      endif
#endif

      call atm_timestep(domain, dt, timeStamp, itimestep)

   end subroutine atm_do_timestep
   
   
   subroutine mpas_core_finalize(domain, stream_manager)
   
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      integer :: ierr

      call mpas_destroy_clock(clock, ierr)
   
   end subroutine mpas_core_finalize


   subroutine atm_compute_mesh_scaling(mesh, configs)

      use mpas_grid_types

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      integer :: iEdge, cell1, cell2
      integer, pointer :: nEdges
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4
      logical, pointer :: config_h_ScaleWithMesh

      call mpas_pool_get_array(mesh, 'meshDensity', meshDensity)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

      call mpas_pool_get_config(configs, 'config_h_ScaleWithMesh', config_h_ScaleWithMesh)

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.5
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)
         end do
      end if

   end subroutine atm_compute_mesh_scaling


   subroutine atm_compute_damping_coefs(mesh, configs)

      use mpas_grid_types
!      use mpas_configure

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      real (kind=RKIND), pointer :: config_xnutr, config_zd
      real (kind=RKIND) :: z, zt, m1, pii
      real (kind=RKIND), dimension(:,:), pointer :: dss, zgrid

      m1 = -1.0
      pii = acos(m1)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(mesh, 'dss', dss)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)

      call mpas_pool_get_config(configs, 'config_zd', config_zd)
      call mpas_pool_get_config(configs, 'config_xnutr', config_xnutr)

      dss(:,:) = 0.0
      do iCell=1,nCells
         zt = zgrid(nVertLevels+1,iCell)
         do k=1,nVertLevels
            z = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
            if (z > config_zd) then
               dss(k,iCell) = config_xnutr*sin(0.5*pii*(z-config_zd)/(zt-config_zd))**2.0
            end if
         end do
      end do

   end subroutine atm_compute_damping_coefs


   subroutine atm_compute_pgf_coefs(mesh, configs)

      use mpas_grid_types
!      use mpas_configure

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      integer :: iEdge, iCell1, iCell2, k, iCell, nz, nz1
      real (kind=RKIND) :: d1, d2, d3
      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, zgrid, pzp, pzm
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nCells, nEdges, nVertLevels
      logical, pointer :: config_newpx

      call mpas_pool_get_array(mesh, 'cpr', cpr)
      call mpas_pool_get_array(mesh, 'cpl', cpl)
      call mpas_pool_get_array(mesh, 'pzp', pzp)
      call mpas_pool_get_array(mesh, 'pzm', pzm)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_config(configs, 'config_newpx', config_newpx)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

!**** coefficient arrays for new pressure gradient calculation

      cpr(:,:) = 0.0
      cpl(:,:) = 0.0

      if (config_newpx) then
         do iEdge=1,nEdges

            iCell1 = cellsOnEdge(1,iEdge)
            iCell2 = cellsOnEdge(2,iEdge)

            d1       = .25*(zgrid(1,iCell2)+zgrid(2,iCell2)-zgrid(1,iCell1)-zgrid(2,iCell1))
            d2       = d1+.5*(zgrid(3,iCell2)-zgrid(1,iCell2))
            d3       = d2+.5*(zgrid(4,iCell2)-zgrid(2,iCell2))
!            cpr(1,iEdge) = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpr(2,iEdge) = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpr(3,iEdge) = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

            cpr(1,iEdge) =  d2/(d2-d1)
            cpr(2,iEdge) = -d1/(d2-d1)
            cpr(3,iEdge) =  0.

            d1       = .25*(zgrid(1,iCell1)+zgrid(2,iCell1)-zgrid(1,iCell2)-zgrid(2,iCell2))
            d2       = d1+.5*(zgrid(3,iCell1)-zgrid(1,iCell1))
            d3       = d2+.5*(zgrid(4,iCell1)-zgrid(2,iCell1))
!            cpl(1,iEdge) = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpl(2,iEdge) = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!            cpl(3,iEdge) = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

            cpl(1,iEdge) =  d2/(d2-d1)
            cpl(2,iEdge) = -d1/(d2-d1)
            cpl(3,iEdge) =  0.

         end do

!         write(6,*) 'cpr1 = ',cpr(1,1),'  cpl1 = ',cpl(1,1)
!         write(6,*) 'cpr2 = ',cpr(2,1),'  cpl2 = ',cpl(2,1)
!         write(6,*) 'cpr3 = ',cpr(3,1),'  cpl3 = ',cpl(3,1)

      else

!        Coefficients for computing vertical pressure gradient dp/dz
!        dp/dz (k,iCell) = pzp(k,iCell) * (p(k+1,iCell) - p(k,iCell)) +pzm(k,iCell) * (p(k,iCell) - p(k-1,iCell))

         nz1 = nVertLevels
         nz = nz1 + 1

         do iCell=1, nCells

            d1 = zgrid(3,iCell)-zgrid(1,iCell)
            d2 = zgrid(4,iCell)-zgrid(2,iCell)
            d3 = d1+d2
            pzm(1,iCell) =  2.*d3/(d1*d2)
            pzp(1,iCell) = -2.*d1/(d2*d3)

            do k=2,nz1-1
               pzp(k,iCell) = 2.*(zgrid(k+1,iCell)-zgrid(k-1,iCell))/     &
     &                      ((zgrid(k+2,iCell)-zgrid(k  ,iCell))*     &
     &                       (zgrid(k+2,iCell)-zgrid(k  ,iCell)       &
     &                       +zgrid(k+1,iCell)-zgrid(k-1,iCell)))
               pzm(k,iCell) = 2.*(zgrid(k+2,iCell)-zgrid(k  ,iCell))/     &
     &                      ((zgrid(k+1,iCell)-zgrid(k-1,iCell))*     &
     &                       (zgrid(k+2,iCell)-zgrid(k  ,iCell)       &
     &                       +zgrid(k+1,iCell)-zgrid(k-1,iCell)))
            end do

            pzp(nz1,iCell) = 0.
            pzm(nz1,iCell) = 2./(zgrid(nz,iCell)-zgrid(nz1-1,iCell))

         end do

      end if

   end subroutine atm_compute_pgf_coefs


   subroutine atm_adv_coef_compression( mesh )

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh


      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge

      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell
      integer, pointer :: nCells, nEdges
      integer :: cell_list(20), ordered_cell_list(20)
      logical :: addcell


      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

      do iEdge = 1, nEdges
         nAdvCellsForEdge(iEdge) = 0
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         !
         ! do only if this edge flux is needed to update owned cells
         !
         if (cell1 <= nCells .or. cell2 <= nCells) then
 
            cell_list(1) = cell1
            cell_list(2) = cell2
            n = 2 
  
          !  add cells surrounding cell 1.  n is number of cells currently in list
            do i = 1, nEdgesOnCell(cell1)
               if (cellsOnCell(i,cell1) /= cell2) then
                  n = n + 1
                  cell_list(n) = cellsOnCell(i,cell1)
               end if
            end do
  
          !  add cells surrounding cell 2 (brute force approach)
            do iCell = 1, nEdgesOnCell(cell2)
               addcell = .true.
               do i=1,n
                  if (cell_list(i) == cellsOnCell(iCell,cell2)) addcell = .false.
               end do
               if (addcell) then
                  n = n+1
                  cell_list(n) = cellsOnCell(iCell,cell2)
               end if
            end do
  
          ! order the list by increasing cell number (brute force approach)
  
            do i=1,n
               ordered_cell_list(i) = nCells + 2
               j_in = 1
               do j=1,n
                  if (ordered_cell_list(i) > cell_list(j) ) then
                     j_in = j
                     ordered_cell_list(i) = cell_list(j)
                  end if
               end do
!               ordered_cell_list(i) = cell_list(j_in)
               cell_list(j_in) = nCells + 3
            end do
  
            nAdvCellsForEdge(iEdge) = n
            do iCell = 1, nAdvCellsForEdge(iEdge)
               advCellsForEdge(iCell,iEdge) = ordered_cell_list(iCell)
            end do
  
          ! we have the ordered list, now construct coefficients
  
            adv_coefs(:,iEdge) = 0.
            adv_coefs_3rd(:,iEdge) = 0.
          
          ! pull together third and fourth order contributions to the flux
          ! first from cell1
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,1,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(1,1,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell1)
               j_in = 0
               do j=1, n
                 if( ordered_cell_list(j) == cellsOnCell(iCell,cell1) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
            end do
  
          ! pull together third and fourth order contributions to the flux
          ! now from cell2
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,2,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(1,2,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell2)
               j_in = 0
               do j=1, n
                  if( ordered_cell_list(j) == cellsOnCell(iCell,cell2) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,2,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(iCell+1,2,iEdge)
            end do
  
            do j = 1,n
               adv_coefs    (j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs    (j,iEdge) / 12.
               adv_coefs_3rd(j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end do
  
          ! 2nd order centered contribution - place this in the main flux weights
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
            j_in = 0
            do j=1, n
               if( ordered_cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
          !  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply
  
            do j=1,n
               adv_coefs    (j,iEdge) = dvEdge(iEdge) * adv_coefs    (j,iEdge)
               adv_coefs_3rd(j,iEdge) = dvEdge(iEdge) * adv_coefs_3rd(j,iEdge)
            end do
 
         end if  ! only do for edges of owned-cells
         
      end do ! end loop over edges

   end subroutine atm_adv_coef_compression


   !***********************************************************************
   !
   !  routine mpas_core_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    September 2011
   !> \details 
   !>  This routine is intended to correctly configure the packages for this MPAS
   !>   core. It can use any Fortran logic to properly configure packages, and it
   !>   can also make use of any namelist options. All variables in the model are
   !>   *not* allocated until after this routine is called.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_setup_packages(configs, packages, ierr)!{{{

      use mpas_packages

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      integer, intent(out) :: ierr

      ierr = 0

   end subroutine mpas_core_setup_packages!}}}


   !***********************************************************************
   !
   !  routine mpas_core_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_setup_clock(core_clock, configs, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      call atm_simulation_clock_init(core_clock, configs, ierr)

   end subroutine mpas_core_setup_clock


   !***********************************************************************
   !
   !  routine mpas_core_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_get_mesh_stream(configs, stream, ierr)

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      character(len=*), intent(out) :: stream
      integer, intent(out) :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         ierr = 1
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end subroutine mpas_core_get_mesh_stream

end module mpas_core
