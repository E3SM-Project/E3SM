! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_core

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_add_anal_incr

   type (MPAS_Clock_type), pointer :: clock


   contains


   function atm_core_init(domain, startTimeStamp) result(ierr)

      use mpas_timekeeping
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_performance, only : mpas_performance_init
      use mpas_atm_dimensions, only : mpas_atm_set_dims
      use mpas_atm_soundings, only : mpas_atm_soundings_init

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
      integer :: ierr

      real (kind=RKIND), pointer :: dt
      type (block_type), pointer :: block

      character(len=StrKIND) :: timeStamp
      integer :: i
      logical, pointer :: config_do_restart

      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: diag
      type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime

      integer, pointer :: nVertLevels, maxEdges, maxEdges2, num_scalars

      ierr = 0

      !
      ! Set up inner dimensions used by arrays in optimized dynamics routines
      !
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'maxEdges', maxEdges)
      call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
      call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
      call mpas_atm_set_dims(nVertLevels, maxEdges, maxEdges2, num_scalars)

      !
      ! Initialize performance module
      !
      call mpas_performance_init()

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
         call MPAS_stream_mgr_read(domain % streamManager, streamID='restart', ierr=ierr)
      else
         call MPAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=ierr)
      end if
      if (ierr /= MPAS_STREAM_MGR_NOERR) then
         write(0,*) ' '
         write(0,*) '********************************************************************************'
         write(0,*) 'Error reading initial conditions'
         call mpas_dmpar_global_abort('********************************************************************************')
      end if
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', direction=MPAS_STREAM_INPUT, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_INPUT, ierr=ierr)

      !
      ! Read all other inputs
      ! For now we don't do this here to match results with previous code; to match requires 
      ! that we read in SST and seaice fields after the call to atm_mpas_init_block()
      !
!      call MPAS_stream_mgr_read(domain % streamManager, ierr=ierr)
!      call MPAS_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr) 

      if (.not. config_do_restart) then
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', state)
            call mpas_pool_initialize_time_levels(state)
            block => block % next
         end do
      end if

      !
      ! Read a new data stream for Incremental Analysis Update (IAU), if config_IAU_option > 0    : HA (June-15-2016)
      ! FIXME: should I check xtime for the IAU fields? Maybe not.
      ! Note: Because the 'iau' stream has the 'iau' package attached to it, the MPAS_stream_mgr_read( )
      !       call here will actually try to read the stream only if IAU is being used in the run.
      !
      call MPAS_stream_mgr_read(domain % streamManager, streamID='iau', whence=MPAS_STREAM_NEAREST, ierr=ierr)
      if (ierr /= MPAS_STREAM_MGR_NOERR) then
         write(0,*) ' '
         write(0,*) '********************************************************************************'
         write(0,*) 'Error reading IAU files'
         call mpas_dmpar_global_abort('********************************************************************************')
      end if
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='iau', ierr=ierr)

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

         call atm_mpas_init_block(domain % dminfo, domain % streamManager, block, mesh, dt)

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

      !
      ! Set up sounding locations
      !
      call mpas_atm_soundings_init(clock, mesh, '3:00:00')

   end function atm_core_init


   subroutine atm_simulation_clock_init(core_clock, configs, ierr)

      use mpas_timekeeping

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
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsistent: using config_run_duration.'
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

      integer :: iCell,iEdge,iVertex
      
      real (kind=RKIND), dimension(:,:), pointer :: u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:), pointer :: areaCell, invAreaCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invDvEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge, invDcEdge
      real (kind=RKIND), dimension(:), pointer :: areaTriangle, invAreaTriangle
      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      integer :: thread
      character(len=StrKIND), pointer :: mminlu

      logical, pointer :: config_do_restart, config_do_DAcycling

      call atm_compute_signs(mesh)
   
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'state', state)

      call mpas_pool_get_subpool(block % structs, 'state', state)

      call mpas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(block % configs, 'config_do_DAcycling', config_do_DAcycling)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!! Compute inverses to avoid divides later
      
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)      
      call mpas_pool_get_array(mesh, 'invDvEdge', invDvEdge)
      call mpas_pool_get_array(mesh, 'invDcEdge', invDcEdge)      
      
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'invAreaTriangle', invAreaTriangle)
      
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)      
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
 
      do iCell=1,nCells 
         invAreaCell(iCell) = 1.0_RKIND / areaCell(iCell)
      end do
        
      do iEdge=1,nEdges 
         invDvEdge(iEdge) = 1.0_RKIND / dvEdge(iEdge)
      end do
      
      do iEdge=1,nEdges 
         invDcEdge(iEdge) = 1.0_RKIND / dcEdge(iEdge)
      end do
      
      do iVertex=1,nVertices
         invAreaTriangle(iVertex) = 1.0_RKIND / areaTriangle(iVertex)
      end do
      
      !!!!! End compute inverses
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call atm_adv_coef_compression(mesh)

      call atm_couple_coef_3rd_order(mesh, block % configs)
      
      call mpas_pool_get_dimension(state, 'nVertices', nVertices)
      call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)

      allocate(ke_vertex(nVertLevels,nVertices+1))  ! ke_vertex is a module variable defined in mpas_atm_time_integration.F
      allocate(ke_edge(nVertLevels,nEdges+1))       ! ke_edge is a module variable defined in mpas_atm_time_integration.F

!$OMP PARALLEL DO
      do thread=1,block % nThreads
         if (.not. config_do_restart .or. (config_do_restart .and. config_do_DAcycling)) then
            call atm_init_coupled_diagnostics(state, 1, diag, mesh, block % configs, &
                                          block % cellThreadStart(thread), block % cellThreadEnd(thread), &
                                          block % vertexThreadStart(thread), block % vertexThreadEnd(thread), &
                                          block % edgeThreadStart(thread), block % edgeThreadEnd(thread), &
                                          block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread), &
                                          block % vertexSolveThreadStart(thread), block % vertexSolveThreadEnd(thread), &
                                          block % edgeSolveThreadStart(thread), block % edgeSolveThreadEnd(thread))
         end if

         call atm_compute_solve_diagnostics(dt, state, 1, diag, mesh, block % configs, &
                                            block % cellThreadStart(thread), block % cellThreadEnd(thread), &
                                            block % vertexThreadStart(thread), block % vertexThreadEnd(thread), &
                                            block % edgeThreadStart(thread), block % edgeThreadEnd(thread))
      end do
!$OMP END PARALLEL DO

      deallocate(ke_vertex)
      deallocate(ke_edge)

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

      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)

      write(0,*) 'min/max of meshScalingDel2 = ', minval(meshScalingDel2(1:nEdgesSolve)), &
                                                  maxval(meshScalingDel2(1:nEdgesSolve))
      write(0,*) 'min/max of meshScalingDel4 = ', minval(meshScalingDel4(1:nEdgesSolve)), &
                                                  maxval(meshScalingDel4(1:nEdgesSolve))

   end subroutine atm_mpas_init_block
   
   
   function atm_core_run(domain) result(ierr)
   
      use mpas_timekeeping
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_derived_types, only : MPAS_STREAM_LATEST_BEFORE, MPAS_STREAM_INPUT, MPAS_STREAM_INPUT_OUTPUT
      use mpas_timer
      use mpas_atm_soundings, only : mpas_atm_soundings_write
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: ierr
   
      real (kind=RKIND), pointer :: dt
      logical, pointer :: config_do_restart
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      integer :: itimestep

      integer :: stream_dir
      character(len=StrKIND) :: input_stream, read_time

      type (mpas_pool_type), pointer :: state, diag, mesh, diag_physics, tend, tend_physics

      ! For high-frequency diagnostics output
      character (len=StrKIND) :: tempfilename

      ! For timing information
      real (kind=R8KIND) :: integ_start_time, integ_stop_time
      real (kind=R8KIND) :: diag_start_time, diag_stop_time
      real (kind=R8KIND) :: input_start_time, input_stop_time
      real (kind=R8KIND) :: output_start_time, output_stop_time
      
      ierr = 0
      
      ! Eventually, dt should be domain specific
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Avoid writing a restart file at the initial time
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      ! Also, for restart runs, avoid writing the initial history or diagnostics fields to avoid overwriting those from the preceding run
      if (config_do_restart) then
         call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='output', direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='diagnostics', direction=MPAS_STREAM_OUTPUT, ierr=ierr)
      end if

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
         call mpas_dmpar_get_time(diag_start_time)
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call atm_compute_output_diagnostics(state, 1, diag, mesh)
            call atm_compute_pv_diagnostics(state, 1, diag, mesh)

            block_ptr => block_ptr % next
         end do
         call mpas_dmpar_get_time(diag_stop_time)
      end if
      call mpas_dmpar_get_time(output_start_time)
      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call mpas_dmpar_get_time(output_stop_time)
      if (ierr /= MPAS_STREAM_MGR_NOERR .and. &
          ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
          ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_REC) then
         write(0,*) ' '
         write(0,*) '********************************************************************************'
         write(0,*) 'Error writing one or more output streams'
         call mpas_dmpar_global_abort('********************************************************************************')
      end if
      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
         write(0,'(a,f9.4,a)') ' Timing for diagnostic computation: ', (diag_stop_time - diag_start_time), ' s'
         write(0,'(a,f9.4,a)') ' Timing for stream output: ', (output_stop_time - output_start_time), ' s'
      end if
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      block_ptr => domain % blocklist
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
      call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
      call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
      call mpas_timer_start('write_soundings')
      call mpas_atm_soundings_write(mesh, state, diag, diag_physics)
      call mpas_timer_stop('write_soundings')

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
         call MPAS_stream_mgr_begin_iteration(domain % streamManager, ierr=ierr)
         do while (MPAS_stream_mgr_get_next_stream(domain % streamManager, streamID=input_stream, directionProperty=stream_dir))
            if (stream_dir == MPAS_STREAM_INPUT .or. stream_dir == MPAS_STREAM_INPUT_OUTPUT) then
               if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID=input_stream, &
                                                  direction=MPAS_STREAM_INPUT, ierr=ierr)) then
                  call mpas_dmpar_get_time(input_start_time)
                  call MPAS_stream_mgr_read(domain % streamManager, streamID=input_stream, whence=MPAS_STREAM_LATEST_BEFORE, &
                                            actualWhen=read_time, ierr=ierr)
                  call mpas_dmpar_get_time(input_stop_time)
                  if (ierr /= MPAS_STREAM_MGR_NOERR) then
                     write(0,*) ' '
                     write(0,*) '********************************************************************************'
                     write(0,*) 'Error reading input stream '//trim(input_stream)
                     call mpas_dmpar_global_abort('********************************************************************************')
                  end if

                  write(0,*) '----------------------------------------------------------------------'
                  write(0,*) '  Read '''//trim(input_stream)//''' input stream valid at '//trim(read_time)
                  write(0,'(a,f9.4,a)') '   Timing for stream input: ', (input_stop_time - input_start_time), ' s'
                  write(0,*) '----------------------------------------------------------------------'

                  call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID=input_stream, direction=MPAS_STREAM_INPUT, ierr=ierr)
               end if
            end if
         end do

         call mpas_timer_start("time integration")
         call mpas_dmpar_get_time(integ_start_time)
         call atm_do_timestep(domain, dt, itimestep)
         call mpas_dmpar_get_time(integ_stop_time)
         call mpas_timer_stop("time integration")   
         write(0,'(a,f9.4,a)') ' Timing for integration step: ', (integ_stop_time - integ_start_time), ' s'

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
         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
           call mpas_dmpar_get_time(diag_start_time)
           block_ptr => domain % blocklist
           do while (associated(block_ptr))
              call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
              call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
              call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
              call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
              call mpas_pool_get_subpool(block_ptr % structs, 'tend', tend)
              call mpas_pool_get_subpool(block_ptr % structs, 'tend_physics', tend_physics)
              call atm_compute_output_diagnostics(state, 1, diag, mesh)
              call atm_compute_pv_diagnostics(state, 1, diag, mesh)
              call atm_compute_pvBudget_diagnostics(state, 1, diag, mesh, tend, tend_physics)

              block_ptr => block_ptr % next
           end do
           call mpas_dmpar_get_time(diag_stop_time)
         end if
         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
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

         call mpas_dmpar_get_time(output_start_time)
         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_dmpar_get_time(output_stop_time)
         if (ierr /= MPAS_STREAM_MGR_NOERR .and. &
             ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
             ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_REC) then
            write(0,*) ' '
            write(0,*) '********************************************************************************'
            write(0,*) 'Error writing one or more output streams'
            call mpas_dmpar_global_abort('********************************************************************************')
         end if
         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            write(0,'(a,f9.4,a)') ' Timing for diagnostic computation: ', (diag_stop_time - diag_start_time), ' s'
            write(0,'(a,f9.4,a)') ' Timing for stream output: ', (output_stop_time - output_start_time), ' s'
         end if

         ! reset any diagnostics here

         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='diagnostics', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            block_ptr => domain % blocklist
            do while (associated(block_ptr))

               call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
               call atm_reset_diagnostics(diag, diag_physics)

               block_ptr => block_ptr % next
            end do
         end if


         ! Only after we've successfully written the restart file should we we
         !    write the restart_timestamp file
         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            if (domain % dminfo % my_proc_id == 0) then
               open(22,file=trim(config_restart_timestamp_name),form='formatted',status='replace')
               write(22,*) trim(timeStamp)
               close(22)
            end if
         end if

         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

         block_ptr => domain % blocklist
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
         call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
         call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
         call mpas_timer_start('write_soundings')
         call mpas_atm_soundings_write(mesh, state, diag, diag_physics)
         call mpas_timer_stop('write_soundings')

      end do
   
   end function atm_core_run
   
   subroutine atm_compute_output_diagnostics(state, time_lev, diag, mesh)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain to be written to history files
   !
   ! Input: state - contains model prognostic fields
   !        mesh  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_constants
      use mpas_atm_interp_diagnostics
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
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

      call interp_diagnostics(mesh, state, time_lev, diag)
   
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

   subroutine atm_reset_diagnostics(diag, diag_physics)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! reset some diagnostics after output
   !
   ! Input: diag          - contains dynamics diagnostic fields
   !        daig_physics  - contains physics diagnostic fields
   !
   ! Output: whatever diagnostics need resetting after output
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: diag_physics
   
      real (kind=RKIND), dimension(:), pointer :: updraft_helicity_max
      real (kind=RKIND), dimension(:), pointer :: w_velocity_max
      real (kind=RKIND), dimension(:), pointer :: wind_speed_level1_max
      real (kind=RKIND), dimension(:), pointer :: refl10cm_1km_max

      call mpas_pool_get_array(diag, 'updraft_helicity_max', updraft_helicity_max)
      call mpas_pool_get_array(diag, 'w_velocity_max', w_velocity_max)
      call mpas_pool_get_array(diag, 'wind_speed_level1_max', wind_speed_level1_max)
      call mpas_pool_get_array(diag_physics, 'refl10cm_1km_max', refl10cm_1km_max)

      updraft_helicity_max(:) = 0.
      w_velocity_max(:) = 0.
      wind_speed_level1_max(:) = 0.
      refl10cm_1km_max(:) = 0.
   
   end subroutine atm_reset_diagnostics


   subroutine atm_do_timestep(domain, dt, itimestep)
   
      use mpas_timekeeping
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
   
   
   function atm_core_finalize(domain) result(ierr)
   
      use mpas_decomp
      use mpas_timekeeping
      use mpas_performance, only : mpas_performance_finalize
      use mpas_atm_soundings, only : mpas_atm_soundings_finalize
   
#ifdef DO_PHYSICS
      use mpas_atmphys_finalize
#endif

      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr 
      integer :: ierr

      ierr = 0
      call mpas_atm_soundings_finalize()

      call mpas_destroy_clock(clock, ierr)
      call mpas_decomp_destroy_decomp_list(domain % decompositions)

      !
      ! Finalize performance module
      !
      call mpas_performance_finalize()

#ifdef DO_PHYSICS
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call atmphys_finalize(block_ptr%configs)

            block_ptr => block_ptr%next
         end do
#endif
   
   end function atm_core_finalize


   subroutine atm_compute_mesh_scaling(mesh, configs)

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
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.25
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**0.75
         end do
      end if

   end subroutine atm_compute_mesh_scaling


   subroutine atm_compute_signs(mesh)

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh

      integer :: i, j, iCell, iVtx
      integer, pointer :: nCells, nVertices, nEdges, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnVertex, verticesOnEdge, cellsOnEdge, edgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, kiteForCell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnVertex_sign, edgesOnCell_sign
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell


      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'kiteForCell', kiteForCell)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)


      do iVtx=1,nVertices
         do i=1,vertexDegree
            if (edgesOnVertex(i,iVtx) <= nEdges) then
               if (iVtx == verticesOnEdge(2,edgesOnVertex(i,iVtx))) then
                  edgesOnVertex_sign(i,iVtx) = 1.0
               else
                  edgesOnVertex_sign(i,iVtx) = -1.0
               end if
            else
               edgesOnVertex_sign(i,iVtx) = 0.0
            end if
         end do
      end do

      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            if (edgesOnCell(i,iCell) <= nEdges) then
               if (iCell == cellsOnEdge(1,edgesOnCell(i,iCell))) then
                  edgesOnCell_sign(i,iCell) = 1.0
                  zb_cell(:,i,iCell) = zb(:,1,edgesOnCell(i,iCell))
                  zb3_cell(:,i,iCell) = zb3(:,1,edgesOnCell(i,iCell))
               else
                  edgesOnCell_sign(i,iCell) = -1.0
                  zb_cell(:,i,iCell) = zb(:,2,edgesOnCell(i,iCell))
                  zb3_cell(:,i,iCell) = zb3(:,2,edgesOnCell(i,iCell))
               end if
            else
               edgesOnCell_sign(i,iCell) = 0.0
            end if
         end do
      end do

      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iVtx = verticesOnCell(i,iCell)
            if (iVtx <= nVertices) then
               do j=1,vertexDegree
                  if (iCell == cellsOnVertex(j,iVtx)) then
                     kiteForCell(i,iCell) = j
                     exit
                  end if
               end do
               if (j > vertexDegree) then
                  write(stderrUnit,*) 'Unexpected error while identifying kiteForCell'
               end if
            else
               kiteForCell(i,iCell) = 1
            end if
         end do
      end do

   end subroutine atm_compute_signs


   subroutine atm_compute_damping_coefs(mesh, configs)

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
      integer :: cell_list(20)
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
  
  
            nAdvCellsForEdge(iEdge) = n
            do iCell = 1, nAdvCellsForEdge(iEdge)
               advCellsForEdge(iCell,iEdge) = cell_list(iCell)
            end do
  
          ! we have the ordered list, now construct coefficients
  
            adv_coefs(:,iEdge) = 0.
            adv_coefs_3rd(:,iEdge) = 0.
          
          ! pull together third and fourth order contributions to the flux
          ! first from cell1
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,1,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(1,1,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell1)
               j_in = 0
               do j=1, n
                 if( cell_list(j) == cellsOnCell(iCell,cell1) ) j_in = j
               end do
               adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
               adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
            end do
  
          ! pull together third and fourth order contributions to the flux
          ! now from cell2
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell2 ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,2,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(1,2,iEdge)
  
            do iCell = 1, nEdgesOnCell(cell2)
               j_in = 0
               do j=1, n
                  if( cell_list(j) == cellsOnCell(iCell,cell2) ) j_in = j
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
               if( cell_list(j) == cell1 ) j_in = j
            end do
            adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
  
            j_in = 0
            do j=1, n
               if( cell_list(j) == cell2 ) j_in = j
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


   subroutine atm_couple_coef_3rd_order(mesh, configs)

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs

      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs_3rd
      real (kind=RKIND), dimension(:,:,:), pointer :: zb3_cell
      real (kind=RKIND), pointer :: config_coef_3rd_order

      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)

      adv_coefs_3rd(:,:) = config_coef_3rd_order * adv_coefs_3rd(:,:)
      zb3_cell(:,:,:) = config_coef_3rd_order * zb3_cell(:,:,:)

   end subroutine atm_couple_coef_3rd_order

   
   real(kind=RKIND) function dotProduct(a, b, sz)

      implicit none

      real(kind=RKIND), dimension(:), intent(in) :: a, b
      integer, intent(in) :: sz

      integer :: i
      real(kind=RKIND) :: rsum

      rsum = 0.0_RKIND

      do i=1,sz
         rsum = rsum + a(i)*b(i)
      end do

      dotProduct = rsum
   end function dotProduct

   integer function elementIndexInArray(val, array, sz)

      implicit none

      integer, intent(in) :: val
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: sz

      integer :: i, ind
      ind = -1
      do i=1,sz
         if (array(i)==val) then
            ind = i
            elementIndexInArray = ind !This returns, right?
            exit !just in case :)
         end if
      end do
      elementIndexInArray = ind
   end function elementIndexInArray
   
   real(kind=RKIND) function formErtelPV(gradxu, gradtheta, density, unitX, unitY, unitZ)

      use mpas_constants, only : omega_e => omega

      implicit none

      real(kind=RKIND), dimension(3), intent(inout) :: gradxu
      real(kind=RKIND), dimension(3), intent(in) :: gradtheta
      real(kind=RKIND), intent(in) :: density
      real(kind=RKIND), dimension(3), intent(in) :: unitX, unitY, unitZ

      real(kind=RKIND) :: epv, eVort
      real(kind=RKIND), dimension(3) :: eVortDir, eVortComponents

      !earth vorticity is in +z-direction in global Cartesian space
      eVort = 2.0 * omega_e
      eVortDir(1) = 0.0_RKIND
      eVortDir(2) = 0.0_RKIND
      eVortDir(3) = eVort

      eVortComponents(1) = dotProduct(eVortDir, unitX,3)
      eVortComponents(2) = dotProduct(eVortDir, unitY,3)
      eVortComponents(3) = dotProduct(eVortDir, unitZ,3)

      gradxu(:) = gradxu(:) + eVortComponents(:)

      epv = dotProduct(gradxu, gradtheta,3) / density

      epv = epv * 1.0e6 !SI to PVUs
    
      formErtelPV = epv
   end function formErtelPV
   
   subroutine local2FullVorticity(gradxu, unitX, unitY, unitZ)
      !given gradxu, return gradxu+earthVort
      
      use mpas_constants, only : omega_e => omega

      implicit none

      real(kind=RKIND), dimension(3), intent(inout) :: gradxu
      real(kind=RKIND), dimension(3), intent(in) :: unitX, unitY, unitZ
      
      real(kind=RKIND) :: eVort
      real(kind=RKIND), dimension(3) :: eVortDir, eVortComponents

      !earth vorticity is in z-direction in global Cartesian space
      eVort = 2.0 * omega_e
      eVortDir(1) = 0.0_RKIND
      eVortDir(2) = 0.0_RKIND
      eVortDir(3) = eVort

      eVortComponents(1) = dotProduct(eVortDir, unitX,3)
      eVortComponents(2) = dotProduct(eVortDir, unitY,3)
      eVortComponents(3) = dotProduct(eVortDir, unitZ,3)

      gradxu(:) = gradxu(:) + eVortComponents(:)
   end subroutine local2FullVorticity
   
   real(kind=RKIND) function calc_verticalVorticity_cell(c0, level, nVerticesOnCell, verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell, vVortVertex)
      !area weighted average of vorticity at vertices to cell center for the specified cell
      !
      implicit none

      real(kind=RKIND), intent(in) :: areaCell
      integer, intent(in) :: c0, level, nVerticesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:,:), intent(in) :: kiteAreasOnVertex, vVortVertex

      real(kind=RKIND) :: vVortCell
      integer :: i, iVertex, cellIndOnVertex

      vVortCell = 0.0_RKIND
      do i = 1,nVerticesOnCell
         iVertex = verticesOnCell(i,c0)
         cellIndOnVertex = elementIndexInArray(c0, cellsOnVertex(:,iVertex), 3)
         vVortCell = vVortCell + kiteAreasOnVertex(cellIndOnVertex, iVertex)*vVortVertex(level, iVertex)/areaCell
      end do

      calc_verticalVorticity_cell = vVortCell
   end function calc_verticalVorticity_cell

   subroutine coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, c0, xyz)

      implicit none

      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors
      integer, intent(in) :: c0
      real(kind=RKIND), dimension(3,3), intent(out) :: xyz

      integer :: i

      xyz(:,1) = cellTangentPlane(:,1,c0) !are these guaranteed unit vectors?
      xyz(:,2) = cellTangentPlane(:,2,c0)
      xyz(:,3) = localVerticalUnitVectors(:,c0)
      do i=1,2
         call normalizeVector(xyz(:,i), 3)
      end do
   end subroutine coordinateSystem_cell

   real(kind=RKIND) function fluxSign(c0, iEdge, cellsOnEdge)
      
      !For finite volume computations, we'll use a normal pointing out of the cell
      implicit none

      integer, intent(in) :: c0
      integer, intent(in) :: iEdge
      integer, dimension(:,:), intent(in) :: cellsOnEdge

      if (c0 == cellsOnEdge(1,iEdge)) then
         fluxSign = 1.0_RKIND
      else
         fluxSign = -1.0_RKIND
      end if
   end function fluxSign

   real(kind=RKIND) function calc_heightCellCenter(c0, level, zgrid)

      implicit none

      integer, intent(in) :: c0, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      calc_heightCellCenter = 0.5*(zgrid(level,c0)+zgrid(level+1,c0))
   end function calc_heightCellCenter

   real(kind=RKIND) function calc_heightVerticalEdge(c0, c1, level, zgrid)

      implicit none

      integer, intent(in) :: c0, c1, level
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid

      real(kind=RKIND) :: hTop, hBottom

      hTop = .5*(zgrid(level+1,c0)+zgrid(level+1,c1))
      hBottom = .5*(zgrid(level,c0)+zgrid(level,c1))

      calc_heightVerticalEdge = hTop-hBottom
   end function calc_heightVerticalEdge

   subroutine normalizeVector(vals, sz)
      !normalize a vector to unit magnitude
      implicit none

      real (kind=RKIND), dimension(:), intent(inout) :: vals
      integer, intent(in) :: sz

      integer :: i
      real (kind=RKIND) :: mag

      mag = 0.0_RKIND !sqrt(sum(squares))
      do i=1,sz
         mag = mag+vals(i)*vals(i)
      end do
      mag = sqrt(mag)
      vals(:) = vals(:)/mag
   end subroutine normalizeVector

   real(kind=RKIND) function calcVolumeCell(areaCell, nEdges, hEdge)

      implicit none

      integer, intent(in) :: nEdges
      real(kind=RKIND), intent(in) :: areaCell
      real(kind=RKIND), dimension(nEdges), intent(in) :: hEdge

      integer :: i
      real(kind=RKIND) :: avgHt, vol

      avgHt = 0.0_RKIND
      do i=1,nEdges
         avgHt = avgHt + hEdge(i)
      end do
      avgHt = avgHt/nEdges

      vol = areaCell*avgHt
      calcVolumeCell = vol
   end function calcVolumeCell

   real(kind=RKIND) function calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
                                                normalEdge, unitDeriv, volumeCell)
      !normals to edges point out of cell
      implicit none

      integer, intent(in) :: nNbrs
      real(kind=RKIND), dimension(:), intent(in) :: valEdges, dvEdge, dhEdge
      real(kind=RKIND), dimension(3,nNbrs), intent(in) :: normalEdge
      real(kind=RKIND), dimension(3), intent(in) :: unitDeriv
      real(kind=RKIND), intent(in) :: volumeCell

      integer :: i
      real(kind=RKIND) :: vale, rsum, areaFace
      real(kind=RKIND), dimension(3) :: unitNormalEdge

      rsum = 0.0_RKIND
      do i=1,nNbrs
         vale = valEdges(i) !0.5 * (val0 + valNbrs(i))
         areaFace = dvEdge(i) * dhEdge(i)
         unitNormalEdge(:) = normalEdge(:,i)
         call normalizeVector(unitNormalEdge,3)
         areaFace = areaFace*dotProduct(unitNormalEdge, unitDeriv,3)  !* abs(dotProduct(unitNormalEdge, unitDeriv,3))
         rsum = rsum + vale * areaFace
      end do
      rsum = rsum / volumeCell

      calc_horizDeriv_fv = rsum
   end function calc_horizDeriv_fv

   !cell centers are halfway between w faces
   real(kind=RKIND) function calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)

      implicit none

      real(kind=RKIND), intent(in) :: val0, valp, valm, z0,zp,zm !center, plus, minus
      
      real(kind=RKIND) :: dval_dzp, dval_dzm

      !Average 1 sided differences to below and above since not equally spaced pts
      dval_dzp = calc_vertDeriv_one(valp, val0, zp-z0)
      dval_dzm = calc_vertDeriv_one(val0, valm, z0-zm)
      calc_vertDeriv_center = 0.5*(dval_dzp+dval_dzm)

   end function calc_vertDeriv_center

   real(kind=RKIND) function calc_vertDeriv_one(valp, valm, dz)
      !1 sided finite difference

      implicit none

      real(kind=RKIND), intent(in) :: valp, valm, dz

      calc_vertDeriv_one = (valp - valm) / dz

   end function calc_vertDeriv_one
   
   subroutine floodFill_strato(mesh, diag, pvuVal, stratoPV)
      !Searching down each column from TOA to find 2pvu surface is buggy with stratospheric wave breaking,
      !since will find 2 pvu at a higher level than "tropopause". This looks to be worse as mesh gets finer and vertical vorticity jumps.
      !Note that stratospheric blobs may persist for long times w/ slow mixing downstream of mountains or deep convection.
      !A few quicker fixes (make sure <2pvu for a number of layers; search down from 10PVU instead of TOA) are hacky and not robust.
      
      !To alleviate the (hopefully) pockets of wave breaking, we can flood fill from a known
      !stratosphere region (e.g., model top > 2pvu) and hopefully filter down around any trouble regions.
      !The problem w/ using only the flood fill is that strong surface PV anomalies can connect to 2pvu, 
      !and the resulting "flood-filled 2 pvu" can have sizeable areas that are just at the surface while there is clearly a tropopause above (e.g., in a cross-section).
      !To address large surface blobs, take the flood fill mask and try to go up from the surface to 10 pvu w/in column. If can, all stratosphere. Else, disconnect "surface blob".
      
      !The "output" is iLev_DT, which is the vertical index for the level >= pvuVal. If >nVertLevels, pvuVal above column. If <2, pvuVal below column.
      !Communication between blocks during the flood fill may be needed to treat some edge cases appropriately.
     
      implicit none
      
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND), intent(in) :: pvuVal, stratoPV
      
      integer :: iCell, k, nChanged, iNbr, iCellNbr
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
      integer, dimension(:,:), pointer :: cellsOnCell
      
      real(kind=RKIND) :: sgnHemi, sgn
      real(kind=RKIND),dimension(:),pointer:: latCell
      real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
      
      integer, dimension(:,:), allocatable :: candInStrato, inStrato
      
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)

      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      allocate(candInStrato(nVertLevels, nCells))
      allocate(inStrato(nVertLevels, nCells))
      candInStrato(:,:) = 0
      inStrato(:,:) = 0
      !store whether each level above DT to avoid repeating logic. we'll use candInStrato as a isVisited marker further below.
      do iCell=1,nCells
         sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
         if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
         do k=1,nVertLevels
            sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
            if (sgn .GE. 0) candInStrato(k,iCell) = 1
         end do
      end do
      
      !seed flood fill with model top that's above DT.
      !can have model top below 2pvu (eg, tropics)
      nChanged = 0
      do iCell=1,nCells
         do k=nVertLevels-5,nVertLevels
            if (candInStrato(k,iCell) .GT. 0) then
               inStrato(k,iCell) = 1
               candInStrato(k,iCell) = 0
               nChanged = nChanged+1
            end if
         end do
      end do
      
      !flood fill from the given seeds. since I don't know enough fortran,
      !we'll just brute force a continuing loop rather than queue.
      do while(nChanged .GT. 0)
        nChanged = 0
        do iCell=1,nCells
          do k=nVertLevels,1,-1
             !update if candidate and neighbor in strato
             if (candInStrato(k,iCell) .GT. 0) then
                !nbr above
                if (k .LT. nVertLevels) then
                  if (inStrato(k+1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
                
                !side nbrs
                do iNbr = 1, nEdgesOnCell(iCell)
                  iCellNbr = cellsOnCell(iNbr,iCell)
                  if (inStrato(k,iCellNbr) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end do
                
                !nbr below
                if (k .GT. 1) then
                  if (inStrato(k-1,iCell) .GT. 0) then
                    inStrato(k,iCell) = 1
                    candInStrato(k,iCell) = 0
                    nChanged = nChanged+1
                    cycle
                  end if
                end if
             end if !candInStrato
          end do !levels
        end do !cells
      end do !while
      
      !Detach high surface PV blobs w/o vertical connection to "stratosphere"
      do iCell=1,nCells
        if (inStrato(1,iCell) .GT. 0) then
          !see how high up we can walk in the column
          do k=2,nVertLevels
            if (inStrato(k,iCell) .LT. 1) then
              exit
            end if !k is highest connected level to sfc
            sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
            if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
            sgn = ertel_pv(k,iCell)*sgnHemi-stratoPV
            if (sgn .LT. 0) then !not actually connected to "stratosphere"
              inStrato(1:k,iCell) = 0
            end if
          end do !k
        end if !inStrato at sfc
      end do !iCell
      
      !Fill iLev_DT with the lowest level above the tropopause (If DT above column, iLev>nVertLevels. If DT below column, iLev=0.
      nChanged = 0
      do iCell=1,nCells
        do k=1,nVertLevels
          if (inStrato(k,iCell) .GT. 0) then
            nChanged = 1
            exit
          end if
        end do !k
        if (nChanged .GT. 0) then !found lowest level
          if (k .EQ. 1) then 
            sgnHemi = sign(1.0_RKIND, latCell(iCell))
            sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
            if (sgn .GT. 0) then !whole column above DT
              iLev_DT(iCell) = 0
            end if
          else
            iLev_DT(iCell) = k
          end if
        else !whole column below DT
          iLev_DT(iCell) = nVertLevels+2
        end if
      end do !iCell
     
   end subroutine floodFill_strato
   
   subroutine interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
      !compute various fields on 2pvu surface using calculated PVU field
      !potential temperature, uZonal, uMeridional, vertical vorticity
      
      implicit none
      
      type (mpas_pool_type), intent(in)  :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
                                          
      real(kind=RKIND),dimension(:),pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv
      real(kind=RKIND),dimension(:,:),pointer:: uZonal, uMeridional, vorticity, theta, ertel_pv, &
                                                kiteAreasOnVertex
      real(kind=RKIND), dimension(:,:), allocatable :: vVort
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCells)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uMeridional)
      call mpas_pool_get_array(diag, 'u_pv', u_pv)
      call mpas_pool_get_array(diag, 'v_pv', v_pv)
      call mpas_pool_get_array(diag, 'theta_pv', theta_pv)
      call mpas_pool_get_array(diag, 'vort_pv', vort_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      !write(0,*) 'Interpolating u,v,theta,vort to pv '
      
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uZonal, u_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, uMeridional, v_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, theta, theta_pv, missingVal, iLev_DT)
                     
      allocate(vVort(nVertLevels, nCells))
      do iCell=1,nCells
         do k=1,nVertLevels
            vVort(k,iCell) = calc_verticalVorticity_cell(iCell, k, nEdgesOnCell(iCell), verticesOnCell, cellsOnVertex, &
                                                         kiteAreasOnVertex, areaCell(iCell), vorticity)
         end do
      end do
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, vVort, vort_pv, missingVal, iLev_DT)
      deallocate(vVort)
      !write(0,*) 'Done interpolating '
   end subroutine interp_pv_diagnostics     
   
   subroutine interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
      !compute various fields on 2pvu surface using calculated PVU field
      !tend_diab, tend_fric
      
      implicit none
      
      type (mpas_pool_type), intent(in)  :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      real(kind=RKIND) ::  pvuVal, missingVal
      
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: iLev_DT
                                          
      real(kind=RKIND),dimension(:),pointer:: latCell, depv_dt_diab_pv, depv_dt_fric_pv
      real(kind=RKIND),dimension(:,:),pointer:: depv_dt_diab, depv_dt_fric, ertel_pv
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCells)
      
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      
      !write(0,*) 'Interpolating u,v,theta,vort to pv '
      
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, depv_dt_diab, depv_dt_diab_pv, missingVal, iLev_DT)
      call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                     ertel_pv, depv_dt_fric, depv_dt_fric_pv, missingVal, iLev_DT)
      !write(0,*) 'Done interpolating '
   end subroutine interp_pvBudget_diagnostics
   
   subroutine interp_pv( nCells, nLevels, interpVal, &
                         latCell, field0, field1,field, &
                         missingVal, iLev_DT)

      implicit none
      !linear-in-PV interpolate columns of field1 to where field0 is interpVal*sign(lat)
      !using level above tropopause already diagnosed
      
      ! input

      integer :: nCells, nLevels
      integer, intent(in) :: iLev_DT(nCells)
      real(kind=RKIND) ::  interpVal, missingVal
      real(kind=RKIND), intent(in) ::latCell(nCells)
      real(kind=RKIND), intent(in) :: field0(nLevels,nCells), field1(nLevels,nCells)
      real(kind=RKIND), intent(out) :: field(nCells)

      !  local
      
      integer :: iCell, iLev, levInd, indlNbr
      real(kind=RKIND) :: valh, vall, vallNbr, sgnh, sgnl, sgnlNbr
      real(kind=RKIND) :: dv_dl, levFrac, valInterpCell, sgnHemi

      do iCell = 1, nCells
        !starting from top, trap val if values on opposite side
        levInd = -1 !what should happen with missing values?
        levFrac = 0.0
        sgnHemi = sign(1.0_RKIND, latCell(iCell)) !problem at the equator...is sign(0)=0?
        if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
        valInterpCell = interpVal*sgnHemi
        
        iLev = iLev_DT(iCell)
        if (iLev .GT. nLevels) then
          levInd = -1
          sgnl = -1.0
        else if (iLev .LT. 1) then
          levInd = -1
          sgnl = 1.0
        else
          valh = field0(iLev,iCell)
          vall = field0(iLev-1,iCell)
          !sandwiched value. equal in case val0 is a vals[l].
          !get linear interpolation: val0 = vals[l]+dvals/dl * dl
          !Avoid divide by 0 by just assuming value is 
          !halfway between...
   
          dv_dl = valh-vall;
          if (abs(dv_dl)<1.e-6) then
            levFrac = 0.5;
          else
            levFrac = (valInterpCell-vall)/dv_dl
          end if
          
          levInd = iLev-1
        end if !iLev in column

        !find value of field using index we just found
        if (levInd<0) then !didn't trap value
          if (sgnl>0.0) then !column above value, take surface
            field(iCell) = field1(1,iCell)
          else !column below value, take top
            !field(iCell) = missingVal
            field(iCell) = field1(nLevels,iCell)
          end if
        else
          valh = field1(levInd+1,iCell)
          vall = field1(levInd,iCell)
        
          dv_dl = valh-vall
          field(iCell) = vall+dv_dl*levFrac
        end if
      end do
      
   end subroutine interp_pv
   
   subroutine calc_gradxu_cell(gradxu, addEarthVort, &
                             iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
      implicit none
      
      real(kind=RKIND), dimension(3), intent(out) :: gradxu
      integer, intent(in) :: addEarthVort, iCell, level, nVertLevels, nEdgesCell0
      real(kind=RKIND), intent(in) :: areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid,uReconstructX, uReconstructY, uReconstructZ, &
                                                      w, vorticity, kiteAreasOnVertex
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm, dw_dx, dw_dy, du_dz, dv_dz
      real(kind=RKIND), dimension(3) :: unitDeriv, velCell0, velCellp, velCellm
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell
     
     !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      !normal vectors at voronoi polygon edges pointing out of cell
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         !dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid) !vertical thickness of that face
         !if don't want to consider 3d cell since we haven't calculated the cell
         !volume well, set all thicknesses to be the same
         dhEdge(i) = 100.0_RKIND

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)
      
      !w
      val0 = .5*(w(level+1, iCell)+w(level, iCell))
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = .5*(w(level+1, iNbr)+w(level, iNbr))
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      dw_dx = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      dw_dy = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !du/dz and dv/dz
      velCell0(1) = uReconstructX(level,iCell)
      velCell0(2) = uReconstructY(level,iCell)
      velCell0(3) = uReconstructZ(level,iCell)
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      if (level>1) then
         !have cell beneath
         velCellm(1) = uReconstructX(level-1,iCell)
         velCellm(2) = uReconstructY(level-1,iCell)
         velCellm(3) = uReconstructZ(level-1,iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         velCellp(1) = uReconstructX(level+1,iCell)
         velCellp(2) = uReconstructY(level+1,iCell)
         velCellp(3) = uReconstructZ(level+1,iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(valp, val0, zp-z0)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_one(val0, valm, z0-zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valm = dotProduct(velCellp, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         !u
         val0 = dotProduct(velCell0, xyzLocal(:,1),3)
         valp = dotProduct(velCellp, xyzLocal(:,1),3)
         valm = dotProduct(velCellm, xyzLocal(:,1),3)
         du_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
         !v
         val0 = dotProduct(velCell0, xyzLocal(:,2),3)
         valp = dotProduct(velCellp, xyzLocal(:,2),3)
         valm = dotProduct(velCellm, xyzLocal(:,2),3)
         dv_dz = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if

      gradxu(3) = calc_verticalVorticity_cell(iCell, level, nEdgesCell0, verticesOnCell, cellsOnVertex, &
                                              kiteAreasOnVertex, areaCell0, vorticity)

      gradxu(1) = dw_dy-dv_dz
      gradxu(2) = du_dz-dw_dx
      
      if (addEarthVort>0) then
        call local2FullVorticity(gradxu, xyzLocal(:,1), xyzLocal(:,2), xyzLocal(:,3))
      end if
     
   end subroutine calc_gradxu_cell
   
   subroutine calc_grad_cell(gradtheta, &
                             iCell, level, nVertLevels, nEdgesCell0, verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell0, &
                             theta)
      !
      implicit none
      
      real(kind=RKIND), dimension(3), intent(out) :: gradtheta
      real(kind=RKIND), intent(in) :: areaCell0
      real(kind=RKIND), dimension(:), intent(in) :: dvEdge
      real(kind=RKIND), dimension(3,2,*), intent(in) :: cellTangentPlane
      real(kind=RKIND), dimension(3,*), intent(in) :: localVerticalUnitVectors, edgeNormalVectors
      real(kind=RKIND), dimension(:,:), intent(in) :: zgrid, theta, kiteAreasOnVertex
      integer, intent(in) :: iCell, level, nVertLevels, nEdgesCell0
      integer, dimension(:,:), intent(in) :: cellsOnCell, edgesOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      
      integer :: i, iNbr, iEdge
      real(kind=RKIND) :: val0, valNbr, volumeCell, areaFactor, z0, zp, zm, valp, valm
      real(kind=RKIND), dimension(3) :: unitDeriv, velCell0, velCellp, velCellm
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(nEdgesCell0) :: valEdges, dvEdgeCell, dhEdge
      real(kind=RKIND), dimension(3,nEdgesCell0) :: normalEdgeCell

      !local coordinate system
      call coordinateSystem_cell(cellTangentPlane, localVerticalUnitVectors, iCell, xyzLocal)
      !normal vectors at voronoi polygon edges pointing out of cell
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         !dhEdge(i) = calc_heightVerticalEdge(iCell, iNbr, level, zgrid) !vertical thickness of that face
         !if don't want to consider 3d cell since we haven't calculated the cell
         !volume well, set all thicknesses to be the same
         dhEdge(i) = 100.0_RKIND

         iEdge = edgesOnCell(i,iCell)
         dvEdgeCell(i) = dvEdge(iEdge)
         val0 = fluxSign(iCell, iEdge, cellsOnEdge)
         normalEdgeCell(:,i) = edgeNormalVectors(:,iEdge)
         call normalizeVector(normalEdgeCell(:,i),3)
         normalEdgeCell(:,i) = normalEdgeCell(:,i)*val0
      end do

      volumeCell = calcVolumeCell(areaCell0, nEdgesCell0, dhEdge)

      !Need to get 3d curl and grad theta
      !horizontal derivatives
      !calc_horizDeriv_fv(valEdges, nNbrs, dvEdge, dhEdge, &
      !                                         normalEdge, unitDeriv, volumeCell)
      !theta
      val0 = theta(level, iCell)
      do i=1,nEdgesCell0
         iNbr = cellsOnCell(i, iCell)
         valNbr = theta(level,iNbr)
         valEdges(i) = 0.5*(valNbr+val0)
      end do
      unitDeriv(:) = xyzLocal(:,1)
      gradtheta(1) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)
      unitDeriv(:) = xyzLocal(:,2)
      gradtheta(2) = calc_horizDeriv_fv(valEdges, nEdgesCell0, dvEdgeCell, dhEdge, normalEdgeCell, unitDeriv, volumeCell)

      !vertical derivatives
      !calc_heightCellCenter(c0, level, zgrid) calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      !theta
      gradtheta(3) = 0.0_RKIND
      z0 = calc_heightCellCenter(iCell, level, zgrid)
      val0 = theta(level, iCell)
      if (level>1) then
         !have cell beneath
         valm = theta(level-1, iCell)
         zm = calc_heightCellCenter(iCell, level-1, zgrid)
      end if
      if (level<nVertLevels) then
         !have cell above
         valp = theta(level+1, iCell)
         zp = calc_heightCellCenter(iCell, level+1, zgrid)
      end if

      if (level==1) then
         !calc_vertDeriv_one(valp, valm, dz)
         gradtheta(3) = calc_vertDeriv_one(valp, val0, zp-z0)
      else if (level==nVertLevels) then
         gradtheta(3) = calc_vertDeriv_one(val0, valm, z0-zm)
      else
         gradtheta(3) = calc_vertDeriv_center(val0, valp, valm, z0,zp,zm)
      end if
   
   end subroutine calc_grad_cell
   
   subroutine calc_vertical_curl(vorticity, u, dcEdge, areaTriangle, verticesOnEdge, nEdges, nVertices)
      ! Adapted from computation of circulation and relative vorticity at each vertex in atm_compute_solve_diagnostics()
      !This takes scvt face values and computes finite volume curl at scvt vertices (triangle cell centers), but
      !only works on 1 horizontal level at a time
      
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: vorticity
      integer, intent(in) :: nEdges, nVertices
      integer, dimension(:,:), intent(in) :: verticesOnEdge
      real (kind=RKIND), dimension(:), intent(in) :: dcEdge, areaTriangle
      real (kind=RKIND), dimension(:), intent(in) :: u
      
      integer :: iEdge, iVertex
      
      vorticity(:) = 0.0
      do iEdge=1,nEdges
            vorticity(verticesOnEdge(1,iEdge)) = vorticity(verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(iEdge)
            vorticity(verticesOnEdge(2,iEdge)) = vorticity(verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(iEdge)
      end do
      do iVertex=1,nVertices
            vorticity(iVertex) = vorticity(iVertex) / areaTriangle(iVertex)
      end do

   end subroutine calc_vertical_curl
   
   subroutine calc_epv(mesh, time_lev, state, diag)
      
      !EPV= absoluteVorticity/density . grad(theta)
      
      implicit none
      
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh

      integer :: iCell, k
      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
      !real(kind=RKIND) :: rvord
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   theta, uReconstructX, uReconstructY, uReconstructZ, &
                                                   ertel_pv
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellTangentPlane', cellTangentPlane)
      call mpas_pool_get_array(mesh, 'localVerticalUnitVectors', localVerticalUnitVectors)
      call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      
      call mpas_pool_get_array(state, 'w', w, time_lev)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      
      !epv and diabatic component ----------------------
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !vort1/rho1
            call calc_gradxu_cell(gradxu, 1, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
            gradxu(:) = gradxu(:)/rho(k,iCell)
            
            !epv1 -------------
            call calc_grad_cell(gradtheta, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             theta)

            ertel_pv(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6 !SI to PVUs
         end do
      end do
   end subroutine calc_epv
   
   subroutine atm_compute_pv_diagnostics(state, time_lev, diag, mesh)
   ! diagnose epv and some fields on horizontal surfaces
   
      use mpas_constants
   
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal, stratoPV
      real (kind=RKIND), dimension(:,:), pointer :: theta, rho, theta_m, rho_zz, zz, dtheta_dt_mix, tend_theta_euler
      type (field2DReal), pointer :: theta_f, uReconstructX_f, uReconstructY_f, uReconstructZ_f, w_f
      type (field2DReal), pointer :: rthratenlw_f, rthratensw_f, rthcuten_f, rthblten_f, dtheta_dt_mp_f, theta_euler_f, dtheta_dt_mix_f
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
            theta(k,iCell) = theta_m(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
         end do
      end do
      
      !nick szapiro
      write(0,*) "Calculating epv"
      
      !need halo cells for everything w/ horizontal derivative
      call mpas_pool_get_field(state, 'w', w_f, time_lev)
      call mpas_pool_get_field(diag, 'uReconstructX', uReconstructX_f)
      call mpas_pool_get_field(diag, 'uReconstructY', uReconstructY_f)
      call mpas_pool_get_field(diag, 'uReconstructZ', uReconstructZ_f)
      call mpas_pool_get_field(diag, 'theta', theta_f)

      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(uReconstructX_f)
      call mpas_dmpar_exch_halo_field(uReconstructY_f)
      call mpas_dmpar_exch_halo_field(uReconstructZ_f)
      call mpas_dmpar_exch_halo_field(w_f)
      
      call calc_epv(mesh, time_lev, state, diag)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      stratoPV = 10.0_RKIND
      call floodFill_strato(mesh, diag, pvuVal, stratoPV)
      call interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
   
   end subroutine atm_compute_pv_diagnostics
   
   subroutine calc_pvBudget(state, time_lev, diag, mesh, tend, tend_physics)
      
      use mpas_vector_reconstruction
      
      implicit none
      
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: tend_physics
      type (mpas_pool_type), intent(inout) :: tend !modify tend_w_euler to uncouple with density

      integer :: iCell, k, iEdge
      integer, pointer :: nCellsSolve, nVertLevels, nVertices, nCells, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex, verticesOnEdge
      !real(kind=RKIND) :: rvord
      real(kind=RKIND), dimension(3) :: gradxu, gradtheta, gradxf
      real(kind=RKIND), dimension(3,3) :: xyzLocal
      real(kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, areaTriangle, dcEdge
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, vorticity, zgrid, &
                                                   localVerticalUnitVectors, edgeNormalVectors, kiteAreasOnVertex, &
                                                   theta, uReconstructX, uReconstructY, uReconstructZ
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab, depv_dt_fric, depv_dt_phys
      real(kind=RKIND), dimension(:,:), pointer :: tend_u_phys, tend_u_euler, rho_edge, tend_w_euler
      real(kind=RKIND), dimension(:,:), pointer :: rthblten, rthcuten, rthratenlw, rthratensw, &
                                                   dtheta_dt_mp, dtheta_dt_mix
      real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
      
      real(kind=RKIND), dimension(:,:), allocatable :: varVerts, tenduX, tenduY, tenduZ, tenduZonal,tendUMerid
      
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellTangentPlane', cellTangentPlane)
      call mpas_pool_get_array(mesh, 'localVerticalUnitVectors', localVerticalUnitVectors)
      call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      
      call mpas_pool_get_array(state, 'w', w, time_lev)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag, 'uReconstructZ', uReconstructZ)
      
      call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)
      call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
      call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
      call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)
      
      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
      
      call mpas_pool_get_array(diag, 'depv_dt_phys', depv_dt_phys)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      
      call mpas_pool_get_array(diag , 'tend_u_phys', tend_u_phys)
      call mpas_pool_get_array(diag , 'rho_edge', rho_edge)
      call mpas_pool_get_array(tend, 'u_euler', tend_u_euler)
      call mpas_pool_get_array(tend, 'w_euler', tend_w_euler)
      
      allocate(varVerts(nVertLevels,nVertices))
      allocate(tenduX(nVertLevels,nCells))
      allocate(tenduY(nVertLevels,nCells))
      allocate(tenduZ(nVertLevels,nCells))
      allocate(tenduZonal(nVertLevels,nCells))
      allocate(tenduMerid(nVertLevels,nCells))
      
      !diabatic component ----------------------
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !vort1/rho1
            call calc_gradxu_cell(gradxu, 1, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             uReconstructX, uReconstructY, uReconstructZ, w,vorticity)
            gradxu(:) = gradxu(:)/rho(k,iCell)
            
            !depv_dt_lw/sw/mp/ -------------
            !absolute vorticity here should maybe be from before taking timestep (from field that generated that tendency...)
            if (associated(rthratenlw)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthratenlw)
               depv_dt_lw(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6 !SI to PVUs
            else
               depv_dt_lw(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthratensw)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthratensw)
               depv_dt_sw(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_sw(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthblten)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthblten)
               depv_dt_bl(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_bl(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(rthcuten)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                rthcuten)
               depv_dt_cu(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_cu(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(dtheta_dt_mp)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                dtheta_dt_mp)
               depv_dt_mp(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_mp(k,iCell) = 0.0_RKIND
            end if
            
            if (associated(dtheta_dt_mp)) then
               call calc_grad_cell(gradtheta, &
                                iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                                cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                                cellsOnVertex, &
                                cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                                dtheta_dt_mix)
               depv_dt_mix(k,iCell) = dotProduct(gradxu, gradtheta,3)* 1.0e6
            else
               depv_dt_mix(k,iCell) = 0.0_RKIND
            end if
         end do 
      end do
      depv_dt_diab = depv_dt_lw + depv_dt_sw + depv_dt_bl + depv_dt_cu + depv_dt_mp + depv_dt_mix
      
      !frictional component ----------------------
      !vertical curl at vertices ( like SAT analogies tend_u:varVerts :: u:vorticity)
      do iEdge=1,nEdges
         do k=1,nVertLevels
            tend_u_phys(k,iEdge) = tend_u_phys(k,iEdge)+tend_u_euler(k,iEdge)/rho_edge(k,iEdge)
         end do
      end do
      !tend_u_phys = tend_u_phys + tend_u_euler/rho_edge
      do k=1,nVertLevels
         call calc_vertical_curl(varVerts(k,:), tend_u_phys(k,:), dcEdge, areaTriangle, verticesOnEdge, nEdges, nVertices)
      end do
      
      !tend_u at cell centers
      call mpas_reconstruct(mesh, tend_u_phys,         &
                               tenduX,tenduY,tenduZ,   &
                               tenduZonal,tenduMerid)
      !uncouple tend_w_euler
      do iCell=1,nCells
         do k=2,nVertLevels
            !average density to vertical interfaces between cells
            !top of lowest cell is interface 2
            tend_w_euler(k,iCell) = tend_w_euler(k,iCell)/( .5*( rho(k-1,iCell)+rho(k,iCell) ) )
         end do
      end do
      !constant extrapolation for density above and below cell centers
      tend_w_euler(1,1:nCells) = tend_w_euler(1,1:nCells)/rho(1,1:nCells)
      tend_w_euler(nVertLevels+1,1:nCells) = tend_w_euler(nVertLevels+1,1:nCells)/rho(nVertLevels,1:nCells)
      
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            !calculating grad(theta)/rho . (grad x F/rho)
            
            !gradtheta term
            call calc_grad_cell(gradtheta, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             theta)
            !
            gradtheta(:) = gradtheta(:)/rho(k,iCell)
            
            !we can call calc_gradxu_cell where:
            !w: tend_w     uReconstruct{X,Y,Z}: tend_u to cell centers     vorticity: tend_u at vertices
            !
            call calc_gradxu_cell(gradxf, 0, &
                             iCell, k, nVertLevels, nEdgesOnCell(iCell), verticesOnCell, kiteAreasOnVertex, &
                             cellsOnCell, edgesOnCell, cellsOnEdge, dvEdge, edgeNormalVectors, &
                             cellsOnVertex, &
                             cellTangentPlane, localVerticalUnitVectors, zgrid, areaCell(iCell), &
                             tenduX, tenduY,tenduZ, tend_w_euler,varVerts)
            
            depv_dt_fric(k,iCell) = dotProduct(gradxf, gradtheta,3)* 1.0e6
         end do
      end do
      
      deallocate(varVerts)
      deallocate(tenduX)
      deallocate(tenduY)
      deallocate(tenduZ)
      deallocate(tenduZonal)
      deallocate(tenduMerid)
      
   end subroutine calc_pvBudget
   
   subroutine atm_compute_pvBudget_diagnostics(state, time_lev, diag, mesh, tend, tend_physics)
      ! after calculating epv,
      !pv budget in the "classic" formulation: (e.g., Pedlosky) 
      !Depv_Dt = Thermo+Friction = vort3d/rho . grad(Dtheta/Dt) + gradTheta/rho . grad x F/rho
      ! The thermo term gets calculated just like epv but theta replaced w/ Dtheta/Dt
      ! F/rho is tend_{u,v,w} and we'll calculate a cell's vertical and horizontal curl separately.
      
      use mpas_constants
      
      implicit none
      
      type (mpas_pool_type), intent(inout) :: diag, tend
      type (mpas_pool_type), intent(in) :: state, mesh, tend_physics
      integer, intent(in) :: time_lev            ! which time level to use from state
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal
      real (kind=RKIND), dimension(:,:), pointer :: dtheta_dt_mix, tend_theta_euler
      type (field2DReal), pointer :: rthratenlw_f, rthratensw_f, rthcuten_f, rthblten_f, dtheta_dt_mp_f, theta_euler_f, dtheta_dt_mix_f
      type (field2DReal), pointer :: tend_u_phys_f, tend_u_euler_f, tend_w_euler_f
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
      
      !nick szapiro
      write(0,*) "Calculating pvBudget"
      
      !need halo cells for everything w/ horizontal derivative
      !Dtheta/Dt
      call mpas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)
      do iCell=1,nCells
         do k=1,nVertLevels
            !with modified moist potential temperature being the model state variable being mixed,
            ! assume qv field is not mixed and so there's no tend_qv to consider
            dtheta_dt_mix(k,iCell) = tend_theta_euler(k,iCell)/( 1._RKIND + rvord*scalars(index_qv,k,iCell) )
         end do
      end do
      call mpas_pool_get_field(tend_physics, 'rthratenlw', rthratenlw_f)
      call mpas_pool_get_field(tend_physics, 'rthratensw', rthratensw_f)
      call mpas_pool_get_field(tend_physics, 'rthcuten', rthcuten_f)
      call mpas_pool_get_field(tend_physics, 'rthblten', rthblten_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mp', dtheta_dt_mp_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mix', dtheta_dt_mix_f)
      
      call mpas_dmpar_exch_halo_field(rthratenlw_f)
      call mpas_dmpar_exch_halo_field(rthratensw_f)
      call mpas_dmpar_exch_halo_field(rthcuten_f)
      call mpas_dmpar_exch_halo_field(rthblten_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mp_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mix_f)
      
      !friction
      call mpas_pool_get_field(diag , 'tend_u_phys', tend_u_phys_f)
      call mpas_pool_get_field(tend, 'u_euler', tend_u_euler_f)
      call mpas_pool_get_field(tend, 'w_euler', tend_w_euler_f)
      call mpas_dmpar_exch_halo_field(tend_u_phys_f)
      call mpas_dmpar_exch_halo_field(tend_u_euler_f)
      call mpas_dmpar_exch_halo_field(tend_w_euler_f)
      
      call calc_pvBudget(state, time_lev, diag, mesh, tend, tend_physics)
      
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      call interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
   
   end subroutine atm_compute_pvBudget_diagnostics

end module atm_core


