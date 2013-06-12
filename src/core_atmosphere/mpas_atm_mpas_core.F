! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework

   type (io_output_object), save :: restart_obj
   type (io_input_object), save :: sfc_update_obj
   integer :: current_outfile_frames
   
   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1 
   integer, parameter :: restartAlarmID = 2 
   integer, parameter :: sfcAlarmID = 3 
   integer, parameter :: hifreqAlarmID = 4

   contains


   subroutine mpas_core_init(domain, startTimeStamp)

      use mpas_configure
      use mpas_kind_types
      use mpas_grid_types

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block

      character(len=StrKIND) :: timeStamp
      integer :: i
      integer :: ierr

      if (.not. config_do_restart) then

         ! Code that was previously in atm_setup_test_case()

         block => domain % blocklist
         do while (associated(block))
            do i=2,nTimeLevs
               call mpas_copy_state(block % state % time_levs(i) % state, block % state % time_levs(1) % state)
            end do
            block => block % next
         end do

      end if


      !
      ! Initialize core
      !
      dt = config_dt

      call atm_simulation_clock_init(domain, dt, startTimeStamp)

      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(1) % state % u)

      block => domain % blocklist
      do while (associated(block))
         call atm_mpas_init_block(domain % dminfo, block, block % mesh, dt)
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
         block => block % next
      end do

      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % pv_edge)
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % ru)
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rw)

      current_outfile_frames = 0

      if (config_sfc_update_interval /= "none") then

         sfc_update_obj % filename = trim(config_sfc_update_name)
         sfc_update_obj % stream = STREAM_SFC

         call mpas_io_input_init(sfc_update_obj, domain % blocklist, domain % dminfo)

         !     
         ! We need to decide which time slice to read from the surface file - read the most recent time slice that falls before or on the start time
         !
         sfc_update_obj % time = MPAS_seekStream(sfc_update_obj % io_stream, trim(startTimeStamp), MPAS_STREAM_LATEST_BEFORE, timeStamp, ierr)
         if (ierr == MPAS_IO_ERR) then
            write(0,*) 'Error: surface update file '//trim(sfc_update_obj % filename)//' did not contain any times at or before '//trim(startTimeStamp)
            call mpas_dmpar_abort(domain % dminfo)
         end if

         write(0,*) 'Starting model with surface time ', trim(timeStamp)

      end if

   end subroutine mpas_core_init


   subroutine atm_simulation_clock_init(domain, dt, startTimeStamp)

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: ierr

      if(trim(config_start_time) == 'file') then
         open(22,file='restart_timestamp',form='formatted',status='old')
         read(22,*) startTimeStamp
         close(22)
      else
        startTimeStamp = config_start_time
      end if
      call mpas_set_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)
      call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_abort(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      ! set high-frequency history alarm, if necessary
      if (trim(config_hifreq_output_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_hifreq_output_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, hifreqAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      ! set sfc alarm, if necessary
      if (trim(config_sfc_update_interval) /= "none") then      
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_sfc_update_interval, ierr=ierr)
         alarmStartTime = startTime
         call mpas_add_clock_alarm(clock, sfcAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if
      
      !TODO: set phyics alarms here...
      !....
      !....

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

   end subroutine atm_simulation_clock_init


   subroutine atm_mpas_init_block(dminfo, block, mesh, dt)
   
      use mpas_grid_types
   !   use atm_advection
      use atm_time_integration
      use mpas_configure
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
#ifdef DO_PHYSICS
!     use mpas_atmphys_aquaplanet
      use mpas_atmphys_control
      use mpas_atmphys_init
      use mpas_atmphys_manager
#endif
   
      implicit none
   
      type (dm_info), intent(in) :: dminfo
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
   
      if (.not. config_do_restart .or. (config_do_restart .and. config_do_DAcycling)) then
         call atm_init_coupled_diagnostics( block % state % time_levs(1) % state, block % diag, mesh)
      end if
      call atm_compute_solve_diagnostics(dt, block % state % time_levs(1) % state, block % diag, mesh)

      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)
      call mpas_reconstruct(mesh, block % state % time_levs(1) % state % u % array, &
                            block % diag % uReconstructX % array,                   &
                            block % diag % uReconstructY % array,                   &
                            block % diag % uReconstructZ % array,                   &
                            block % diag % uReconstructZonal % array,               &
                            block % diag % uReconstructMeridional % array           &
                           )
   
   !
   ! Note: The following initialization calls have been moved to mpas_setup_test_case()
   !       since values computed by these routines are needed to produce initial fields
   !
   !   call atm_initialize_advection_rk(mesh)
   !   call atm_initialize_deformation_weights(mesh)

#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at least one physics
      !parameterization is called (using the logical moist_physics):
      call physics_namelist_check

      !proceed with initialization of physics parameterization if moist_physics is set to true:
      if(moist_physics) then
         !initialization of seom input variables in registry:
         call physics_registry_init(config_do_restart, mesh, block % sfc_input)
         call physics_run_init(mesh,block % state % time_levs(1) % state,clock)

         !initialization of all physics:
         call physics_init(dminfo, clock, config_do_restart, mesh, block % state % time_levs(1) % state, &
                           block % state % time_levs(1) % state, block % diag_physics, block % sfc_input)

      endif
#endif
   
      call atm_compute_mesh_scaling(mesh)

      call atm_compute_damping_coefs(mesh)

      call atm_compute_pgf_coefs(mesh)

      write(0,*) 'min/max of meshScalingDel2 = ', minval(mesh % meshScalingDel2 % array(1:mesh%nEdges)), &
                                                  maxval(mesh % meshScalingDel2 % array(1:mesh%nEdges))
      write(0,*) 'min/max of meshScalingDel4 = ', minval(mesh % meshScalingDel4 % array(1:mesh%nEdges)), &
                                                  maxval(mesh % meshScalingDel4 % array(1:mesh%nEdges))

      call atm_adv_coef_compression(mesh)
   
   end subroutine atm_mpas_init_block
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_io_output
      use mpas_io_input, only: insert_string_suffix => mpas_insert_string_suffix
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: itimestep
      integer :: ierr

      ! For high-frequency diagnostics output
      character (len=StrKIND) :: tempfilename

      ! Eventually, dt should be domain specific
      dt = config_dt

      call atm_write_output_frame(output_obj, output_frame, domain)

      if (trim(config_hifreq_output_interval) /= 'none') then
         block_ptr => domain % blocklist
         call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

         call insert_string_suffix('diagnostics.nc', trim(timeStamp), tempfilename)
         call write_hifreq_output(block_ptr, tempfilename)
      end if

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 1
      do while (.not. mpas_is_clock_stop_time(clock))

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         write(0,*) 'Begin timestep ', trim(timeStamp)

         ! Input external updates (i.e. surface)
         if (mpas_is_alarm_ringing(clock, sfcAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, sfcAlarmID, ierr=ierr)

            call mpas_read_and_distribute_fields(sfc_update_obj)
            sfc_update_obj % time = sfc_update_obj % time + 1
         end if

         call mpas_timer_start("time integration")
         call atm_do_timestep(domain, dt, itimestep)
         call mpas_timer_stop("time integration")   

         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_shift_time_levels_state(domain % blocklist % state)
         

         ! Advance clock before writing output
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         !TODO: MPAS_getClockRingingAlarms is probably faster than multiple MPAS_isAlarmRinging...
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if
            call atm_write_output_frame(output_obj, output_frame, domain)
         end if

         if (mpas_is_alarm_ringing(clock, hifreqAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, hifreqAlarmID, ierr=ierr)

            block_ptr => domain % blocklist
            call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)

            call insert_string_suffix('diagnostics.nc', trim(timeStamp), tempfilename)
            call write_hifreq_output(block_ptr, tempfilename)
         end if

         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
               call atm_compute_restart_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)
               block_ptr => block_ptr % next
            end do

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
         end if

      end do
   
   end subroutine mpas_core_run
   
   
   subroutine atm_write_output_frame(output_obj, output_frame, domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_io_output
   
      implicit none
   
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
   
      type (block_type), pointer :: block_ptr
   
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call atm_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)
         block_ptr => block_ptr % next
      end do
   
      call mpas_output_state_for_domain(output_obj, domain, output_frame)
      output_frame = output_frame + 1
   
      ! reset frame if the maximum number of frames per outfile has been reached
      if (config_frames_per_outfile > 0) then
         current_outfile_frames = current_outfile_frames + 1            
         if(current_outfile_frames >= config_frames_per_outfile) then
            current_outfile_frames = 0
            output_frame = 1
         end if
      end if

   end subroutine atm_write_output_frame
   
   
   subroutine atm_compute_output_diagnostics(state, diag, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to history files
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
      use mpas_atm_interp_diagnostics
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(in) :: grid
   
      integer :: iCell, k

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * grid % zz % array(k,iCell)
         end do
      end do

      call interp_diagnostics(grid,state,diag)
   
   end subroutine atm_compute_output_diagnostics
   
   
   subroutine atm_compute_restart_diagnostics(state, diag, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to restart files
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_constants
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(in) :: grid
   
      integer :: iCell, k

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * grid % zz % array(k,iCell)
         end do
      end do
   
   end subroutine atm_compute_restart_diagnostics


   subroutine write_hifreq_output(block_ptr, fname)

      implicit none

      type (block_type), pointer :: block_ptr
      character (len=*), intent(in) :: fname

      integer :: ierr
      type (MPAS_Stream_type) :: hifreq_stream

      call MPAS_createStream(hifreq_stream, trim(fname), MPAS_IO_PNETCDF, MPAS_IO_WRITE, 1, ierr)

      if (block_ptr % mesh % on_a_sphere) then
         call MPAS_writeStreamAtt(hifreq_stream, 'on_a_sphere', 'YES             ', ierr)
      else
         call MPAS_writeStreamAtt(hifreq_stream, 'on_a_sphere', 'NO              ', ierr)
      end if
      call MPAS_writeStreamAtt(hifreq_stream, 'sphere_radius', block_ptr % mesh % sphere_radius, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'model_name', block_ptr % domain % modelName, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'core_name', block_ptr % domain % coreName, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'model_version', block_ptr % domain % modelVersion, ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'source', 'MPAS', ierr)
      call MPAS_writeStreamAtt(hifreq_stream, 'Conventions', 'MPAS', ierr)

      call MPAS_streamAddField(hifreq_stream, block_ptr % state % time_levs(1) % state % xtime, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % olrtoa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % rainc, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % rainnc, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % refl10cm_max, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag_physics % precipw, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % temperature_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % height_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % uzonal_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % umeridional_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % w_850hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_200hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_500hPa, ierr)
      call MPAS_streamAddField(hifreq_stream, block_ptr % diag % vorticity_850hPa, ierr)
      call MPAS_writeStream(hifreq_stream, 1, ierr)
      call MPAS_closeStream(hifreq_stream, ierr)

   end subroutine write_hifreq_output
   
   
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
   
   
   subroutine mpas_core_finalize(domain)
   
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer :: ierr

      if (config_sfc_update_interval /= "none") call mpas_io_input_finalize(sfc_update_obj, domain % dminfo)

      call mpas_destroy_clock(clock, ierr)
   
   end subroutine mpas_core_finalize


   subroutine atm_compute_mesh_scaling(mesh)

      use mpas_grid_types

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, cell1, cell2
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4

      meshDensity => mesh % meshDensity % array
      meshScalingDel2 => mesh % meshScalingDel2 % array
      meshScalingDel4 => mesh % meshScalingDel4 % array

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,mesh%nEdges
            cell1 = mesh % cellsOnEdge % array(1,iEdge)
            cell2 = mesh % cellsOnEdge % array(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.5
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)
         end do
      end if

   end subroutine atm_compute_mesh_scaling


   subroutine atm_compute_damping_coefs(mesh)

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iCell, k
      real (kind=RKIND) :: z, zt, m1, pii
      real (kind=RKIND), dimension(:,:), pointer :: dss, zgrid

      m1 = -1.0
      pii = acos(m1)

      dss => mesh % dss % array
      zgrid => mesh % zgrid % array

      dss(:,:) = 0.0
      do iCell=1,mesh%nCells
         zt = zgrid(mesh%nVertLevels+1,iCell)
         do k=1,mesh%nVertLevels
            z = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
            if (z > config_zd) then
               dss(k,iCell) = config_xnutr*sin(0.5*pii*(z-config_zd)/(zt-config_zd))**2.0
            end if
         end do
      end do

   end subroutine atm_compute_damping_coefs


   subroutine atm_compute_pgf_coefs(mesh)

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, iCell1, iCell2, k, iCell, nz, nz1
      real (kind=RKIND) :: d1, d2, d3
      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, zgrid, pzp, pzm

      cpr   => mesh % cpr % array
      cpl   => mesh % cpl % array
      pzp   => mesh % pzp % array
      pzm   => mesh % pzm % array
      zgrid => mesh % zgrid % array

!**** coefficient arrays for new pressure gradient calculation

      cpr(:,:) = 0.0
      cpl(:,:) = 0.0

      if (config_newpx) then
         do iEdge=1,mesh%nEdges

            iCell1 = mesh % cellsOnEdge % array(1,iEdge)
            iCell2 = mesh % cellsOnEdge % array(2,iEdge)

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

         nz1 = mesh % nVertLevels
         nz = nz1 + 1

         do iCell=1, mesh % nCells

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


   subroutine atm_adv_coef_compression( grid )

      implicit none

      type (mesh_type), intent(inout) :: grid


      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge

      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell
      integer :: cell_list(20), ordered_cell_list(20)
      logical :: addcell

      deriv_two => grid % deriv_two % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array

      do iEdge = 1, grid % nEdges
         nAdvCellsForEdge(iEdge) = 0
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         !
         ! do only if this edge flux is needed to update owned cells
         !
         if (cell1 <= grid%nCells .or. cell2 <= grid%nCells) then
 
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
               ordered_cell_list(i) = grid % nCells + 2
               j_in = 1
               do j=1,n
                  if (ordered_cell_list(i) > cell_list(j) ) then
                     j_in = j
                     ordered_cell_list(i) = cell_list(j)
                  end if
               end do
!               ordered_cell_list(i) = cell_list(j_in)
               cell_list(j_in) = grid % nCells + 3
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
               adv_coefs    (j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs    (j,iEdge) / 12.
               adv_coefs_3rd(j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
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
               adv_coefs    (j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs    (j,iEdge)
               adv_coefs_3rd(j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_3rd(j,iEdge)
            end do
 
         end if  ! only do for edges of owned-cells
         
      end do ! end loop over edges

   end subroutine atm_adv_coef_compression

end module mpas_core
