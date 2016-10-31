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

   type (MPAS_Clock_type), pointer :: clock


   contains


   function atm_core_init(domain, startTimeStamp) result(ierr)

      use mpas_timekeeping
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_performance, only : mpas_performance_init
      use mpas_atm_dimensions, only : mpas_atm_set_dims
      use mpas_atm_soundings, only : mpas_atm_soundings_init
      use mpas_atm_diagnostics_manager, only : mpas_atm_diag_setup

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

      call mpas_atm_diag_setup(domain % streamManager, domain % blocklist % structs, domain % clock)

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
      use mpas_atm_diagnostics_manager, only : mpas_atm_diag_update, mpas_atm_diag_compute, mpas_atm_diag_reset
   
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

      call mpas_dmpar_get_time(diag_start_time)

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call atm_compute_output_diagnostics(state, 1, diag, mesh)

            block_ptr => block_ptr % next
         end do
      end if
      call mpas_atm_diag_reset()
      call mpas_atm_diag_update()
      call mpas_atm_diag_compute()

      call mpas_dmpar_get_time(diag_stop_time)

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

      call mpas_atm_diag_reset()

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


         call mpas_dmpar_get_time(diag_start_time)

         !
         ! Write any output streams that have alarms ringing, after computing diagnostics fields
         !
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
           block_ptr => domain % blocklist
           do while (associated(block_ptr))
              call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
              call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
              call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
              call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
              call mpas_pool_get_subpool(block_ptr % structs, 'tend', tend)
              call mpas_pool_get_subpool(block_ptr % structs, 'tend_physics', tend_physics)
              call atm_compute_output_diagnostics(state, 1, diag, mesh)

              block_ptr => block_ptr % next
           end do
         end if

         call mpas_atm_diag_update()
         call mpas_atm_diag_compute()

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
         call mpas_dmpar_get_time(diag_stop_time)

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

         call mpas_atm_diag_reset()

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
   
      real (kind=RKIND), dimension(:), pointer :: refl10cm_1km_max

      call mpas_pool_get_array(diag_physics, 'refl10cm_1km_max', refl10cm_1km_max)

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
      use mpas_atm_diagnostics_manager, only : mpas_atm_diag_cleanup
   
#ifdef DO_PHYSICS
      use mpas_atmphys_finalize
#endif

      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr 
      integer :: ierr

      ierr = 0
      call mpas_atm_soundings_finalize()

      call mpas_atm_diag_cleanup()

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


end module atm_core

