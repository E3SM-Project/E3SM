! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework
   use mpas_timekeeping

   private
   public :: &
        mpas_core_init, &
        mpas_core_run, &
        mpas_core_finalize, &
        mpas_core_setup_packages, &
        mpas_core_setup_decompositions, &
        mpas_core_setup_clock, &
        mpas_core_get_mesh_stream

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  routine mpas_core_init
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine mpas_core_init(domain, startTimeStamp)
   
      use mpas_derived_types
      use mpas_stream_manager
      use cice_forcing, only: cice_forcing_init
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: diagnosticsPool

      logical, pointer :: config_do_restart
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      integer :: ierr

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)


      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock


      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

      !
      ! If this is a restart run, read the restart stream, else read the input stream.
      ! Regardless of which stream we read for initial conditions, reset the
      ! input alarms for both input and restart before reading any remaining input streams.
      !
      if (config_do_restart) then
         call MPAS_stream_mgr_read(domain % streamManager, streamID='restart', ierr=ierr)
      else
         call MPAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=ierr)
      end if
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', direction=MPAS_STREAM_INPUT, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_INPUT, ierr=ierr)

      ! Read all other inputs
      call MPAS_stream_mgr_read(domain % streamManager, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

      ! initialize forcing
      call cice_forcing_init(domain, domain % streamManager)

      !
      ! Initialize core
      !
      dt = config_dt

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
         call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)
         xtime = startTimeStamp

         block => block % next
      end do

      call mpas_init_block(domain, dt)

   end subroutine mpas_core_init

   !***********************************************************************
   !
   !  routine simulation_clock_init
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine simulation_clock_init(core_clock, configs, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: local_err

      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time, config_restart_timestamp_name
      character (len=StrKIND) :: restartTimeStamp
      real (kind=RKIND), pointer :: config_dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)


      ! Set time to the user-specified start time OR use a restart time from file
      if ( trim(config_start_time) == "file" ) then
         open(22, file=config_restart_timestamp_name, form='formatted', status='old')
         read(22,*) restartTimeStamp
         close(22)
         call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=local_err)
      else
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=local_err)
      end if

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
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if

   end subroutine simulation_clock_init
   
   !***********************************************************************
   !
   !  routine mpas_init_block
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine mpas_init_block(domain, dt)
   
      use mpas_derived_types
      
      use cice_initialize, only: cice_init

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      ! initialize general stuff
      call cice_init(domain, clock, dt)

   end subroutine mpas_init_block

   !***********************************************************************
   !
   !  routine mpas_core_run
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------
   
   subroutine mpas_core_run(domain)
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_timer
      use cice_forcing, only: cice_forcing_get
   
      implicit none
   
      type (domain_type), intent(inout) :: domain

      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr
      
      real (kind=RKIND), pointer :: config_dt
      character(len=StrKIND), pointer :: config_restart_timestamp_name

      ! configurations
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Eventually, dt should be domain specific
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
      write(0,*) 'Initial timestep ', trim(timeStamp)

      ! Avoid writing a restart file at the initial time
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      itimestep = 1
      call mpas_advance_clock(clock)

      ! initial forcing
      call cice_forcing_get(domain % streamManager, domain, clock, .true.)

      timeLoop: do

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         write(0,*) 'Doing timestep ', trim(timeStamp)

         ! do the time integration
         call mpas_timer_start("time integration")
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration")

         ! update the restart_timestamp file with the new time if needed
         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr) ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) timeStamp
            close(22)
         end if

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

         ! check if the simulation is finished
         if (mpas_is_clock_stop_time(clock)) exit timeLoop

         ! advance the clock
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         ! forcing
         call cice_forcing_get(domain % streamManager, domain, clock, .false.)

      end do timeLoop

   end subroutine mpas_core_run
   
   !***********************************************************************
   !
   !  routine mpas_timestep
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine mpas_timestep(domain, itimestep, dt, timeStamp)
   
      use mpas_derived_types
      use cice_time_integration
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block

      type (MPAS_pool_type), pointer :: diagnostics
      character(len=strKIND), pointer :: xtime

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)
         call MPAS_pool_get_array(diagnostics, "xtime", xtime)
         xtime = timeStamp

         block => block % next
      end do

      call cice_timestep(domain, clock, dt, itimestep)

   end subroutine mpas_timestep
   
   !***********************************************************************
   !
   !  routine mpas_core_finalize
   !
   !> \brief   
   !> \author  
   !> \date    
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------
   
   subroutine mpas_core_finalize(domain)
   
      use mpas_derived_types
      use mpas_decomp
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr
 
      call mpas_destroy_clock(clock, ierr)
      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end subroutine mpas_core_finalize

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

   subroutine mpas_core_setup_packages(configPool, packagePool, ierr)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool
      integer, intent(out) :: ierr

      ierr = 0

      ! dynamics 
      call setup_packages_dynamics(configPool, packagePool, ierr)

      ! column physics
      call setup_packages_column_physics(configPool, packagePool, ierr)

   end subroutine mpas_core_setup_packages!}}}

   !***********************************************************************
   !
   !  routine setup_packages_dynamics
   !
   !> \brief   
   !> \author  Adrian K. Turner
   !> \date    9th Feburary 2015
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_dynamics(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr

     logical, pointer :: &
          config_use_velocity_solver

     character(len=strKIND), pointer :: &
          config_stress_divergence_scheme, &
          config_variational_basis

     logical, pointer :: &
          pkgWeakActive, &
          pkgVariationalActive, &
          pkgWachspressActive, &
          pkgPieceWiseLinearActive

     !pkgWeak
     !pkgVariational
     !pkgWachspress
     !pkgPieceWiseLinear
     
     call MPAS_pool_get_config(configPool, "config_use_velocity_solver", config_use_velocity_solver)
     call MPAS_pool_get_config(configPool, "config_stress_divergence_scheme", config_stress_divergence_scheme)
     call MPAS_pool_get_config(configPool, "config_variational_basis", config_variational_basis)

     call MPAS_pool_get_package(packagePool, "pkgWeakActive", pkgWeakActive)
     call MPAS_pool_get_package(packagePool, "pkgVariationalActive", pkgVariationalActive)
     call MPAS_pool_get_package(packagePool, "pkgWachspressActive", pkgWachspressActive)
     call MPAS_pool_get_package(packagePool, "pkgPieceWiseLinearActive", pkgPieceWiseLinearActive)

     ! is the velocity solver on
     if (config_use_velocity_solver) then

        ! stress divergence scheme
        if (trim(config_stress_divergence_scheme) == "weak") then

           ! weak stress divergence scheme
           pkgWeakActive = .true.

        else if (trim(config_stress_divergence_scheme) == "variational") then

           ! variational stress divergence scheme
           pkgVariationalActive = .true.

           ! variational basis function scheme
           if (trim(config_variational_basis) == "wachspress") then

              ! wachspress basis functions
              pkgWachspressActive = .true.

           else if (trim(config_variational_basis) == "pwl") then

              ! piecewise linear basis functions
              pkgPieceWiseLinearActive = .true.

           else

              ! unknown value of config_variational_basis
              write(*,*) "Unknown variational scheme: ", trim(config_variational_basis)
              ierr = 1

           endif ! config_variational_basis

        else

           ! unknown value of config_stress_divergence_scheme
           write(*,*) "Unknown stress divergence scheme: ", trim(config_stress_divergence_scheme)
           ierr = 1

        endif ! config_stress_divergence_scheme

     endif ! config_use_velocity_solver

   end subroutine setup_packages_dynamics!}}}

   !***********************************************************************
   !
   !  routine setup_packages_dynamics
   !
   !> \brief   
   !> \author  Adrian K. Turner
   !> \date    9th Feburary 2015
   !> \details 
   !>  
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_column_physics(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr
     
     ! column physics package packages
     logical, pointer :: &
          config_use_column_physics_package, &
          config_use_column_shortwave, &
          config_use_column_vertical_thermodynamics, &
          config_use_column_biogeochemistry, &
          config_use_column_itd_thermodynamics, &
          config_use_column_ridging

     logical, pointer :: &
          pkgColumnShortwaveActive, &
          pkgColumnVerticalThermodynamicsActive, &
          pkgColumnBiogeochemistryActive, &
          pkgColumnITDThermodynamicsActive, &
          pkgColumnRidgingActive, &
          pkgColumnAerosolsActive

     ! column tracers
     logical, pointer :: &
          config_use_ice_age, &
          config_use_first_year_ice, &
          config_use_level_ice, &
          config_use_meltponds, &
          config_use_cesm_meltponds, &
          config_use_level_meltponds, &
          config_use_topo_meltponds, &
          config_use_aerosols, &
          config_use_brine

     logical, pointer :: &
          pkgColumnTracerIceAgeActive, &
          pkgColumnTracerFirstYearIceActive, &
          pkgColumnTracerLevelIceActive, &
          pkgColumnTracerPondsActive, &
          pkgColumnTracerLidThicknessActive, &
          pkgColumnTracerAerosolsActive

     ! other packages
     logical, pointer :: &
          config_use_form_drag

     logical, pointer :: &
          pkgColumnFormDragActive

     !-----------------------------------------------------------------------
     ! column main routines
     !-----------------------------------------------------------------------

     !pkgColumnShortwave
     !pkgColumnVerticalThermodynamics
     !pkgColumnBiogeochemistry
     !pkgColumnITDThermodynamics
     !pkgColumnRidging

     call MPAS_pool_get_config(configPool, "config_use_column_physics_package", config_use_column_physics_package)
     call MPAS_pool_get_config(configPool, "config_use_column_shortwave", config_use_column_shortwave)
     call MPAS_pool_get_config(configPool, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)
     !call MPAS_pool_get_config(configPool, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
     call MPAS_pool_get_config(configPool, "config_use_column_itd_thermodynamics", config_use_column_itd_thermodynamics)
     call MPAS_pool_get_config(configPool, "config_use_column_ridging", config_use_column_ridging)

     call MPAS_pool_get_package(packagePool, "pkgColumnShortwaveActive", pkgColumnShortwaveActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnVerticalThermodynamicsActive", pkgColumnVerticalThermodynamicsActive)
     !call MPAS_pool_get_package(packagePool, "pkgColumnBiogeochemistryActive", pkgColumnBiogeochemistryActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnITDThermodynamicsActive", pkgColumnITDThermodynamicsActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnRidgingActive", pkgColumnRidgingActive)

     pkgColumnShortwaveActive              = config_use_column_shortwave
     pkgColumnVerticalThermodynamicsActive = config_use_column_vertical_thermodynamics
     !pkgColumnBiogeochemistryActive        = config_use_column_biogeochemistry
     pkgColumnITDThermodynamicsActive      = config_use_column_itd_thermodynamics
     pkgColumnRidgingActive                = config_use_column_ridging

     if (.not. config_use_column_physics_package) then
        pkgColumnShortwaveActive              = .false.
        pkgColumnVerticalThermodynamicsActive = .false.
        !pkgColumnBiogeochemistryActive        = .false.
        pkgColumnITDThermodynamicsActive      = .false.
        pkgColumnRidgingActive                = .false.
     endif

     pkgColumnShortwaveActive              = .true.
     pkgColumnVerticalThermodynamicsActive = .true.
     !pkgColumnBiogeochemistryActive        = .true.
     pkgColumnITDThermodynamicsActive      = .true.
     pkgColumnRidgingActive                = .true.

     !-----------------------------------------------------------------------
     ! tracers
     !-----------------------------------------------------------------------

     !pkgColumnTracerIceAge
     !pkgColumnTracerFirstYearIce
     !pkgColumnTracerLevelIce
     !pkgColumnTracerPonds
     !pkgColumnTracerLidThickness
     !pkgColumnTracerAerosols

     call MPAS_pool_get_config(configPool, "config_use_ice_age", config_use_ice_age)
     call MPAS_pool_get_config(configPool, "config_use_first_year_ice", config_use_first_year_ice)
     call MPAS_pool_get_config(configPool, "config_use_level_ice", config_use_level_ice)
     call MPAS_pool_get_config(configPool, "config_use_meltponds", config_use_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_cesm_meltponds", config_use_cesm_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_level_meltponds", config_use_level_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_topo_meltponds", config_use_topo_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_aerosols", config_use_aerosols)
     call MPAS_pool_get_config(configPool, "config_use_brine", config_use_brine)

     call MPAS_pool_get_package(packagePool, "pkgColumnTracerIceAgeActive", pkgColumnTracerIceAgeActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerFirstYearIceActive", pkgColumnTracerFirstYearIceActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerLevelIceActive", pkgColumnTracerLevelIceActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerPondsActive", pkgColumnTracerPondsActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerLidThicknessActive", pkgColumnTracerLidThicknessActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerAerosolsActive", pkgColumnTracerAerosolsActive)

     pkgColumnTracerIceAgeActive       = config_use_ice_age
     pkgColumnTracerFirstYearIceActive = config_use_first_year_ice
     pkgColumnTracerLevelIceActive     = config_use_level_ice
     pkgColumnTracerPondsActive        = config_use_meltponds
     pkgColumnTracerLidThicknessActive = (config_use_level_meltponds .or. config_use_topo_meltponds)
     pkgColumnTracerAerosolsActive     = config_use_aerosols

     if (.not. config_use_column_physics_package) then
        pkgColumnTracerIceAgeActive       = .false.
        pkgColumnTracerFirstYearIceActive = .false.
        pkgColumnTracerLevelIceActive     = .false.
        pkgColumnTracerPondsActive        = .false.
        pkgColumnTracerLidThicknessActive = .false.
        pkgColumnTracerAerosolsActive     = .false.
     endif

     pkgColumnTracerIceAgeActive       = .true.
     pkgColumnTracerFirstYearIceActive = .true.
     pkgColumnTracerLevelIceActive     = .true.
     pkgColumnTracerPondsActive        = .true.
     pkgColumnTracerLidThicknessActive = .true.
     pkgColumnTracerAerosolsActive     = .true.

     !-----------------------------------------------------------------------
     ! other column packages
     !-----------------------------------------------------------------------

     !pkgColumnFormDrag

     ! form drag
     call MPAS_pool_get_config(configPool, "config_use_form_drag", config_use_form_drag)
     call MPAS_pool_get_package(packagePool, "pkgColumnFormDragActive", pkgColumnFormDragActive)
     pkgColumnFormDragActive = (config_use_column_physics_package .and. config_use_column_vertical_thermodynamics .and.config_use_form_drag)

     pkgColumnFormDragActive = .true.

   end subroutine setup_packages_column_physics!}}}


   !***********************************************************************
   !
   !  routine mpas_core_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    September 2011
   !> \details
   !>  This routine is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_setup_decompositions(domain, ierr)!{{{

      use mpas_derived_types
      use mpas_decomp

      implicit none

      type (domain_type), intent(in) :: domain

      integer, intent(out) :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(domain % decompositions)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(domain % decompositions, 'uniform', decompFunc, iErr)

   end subroutine mpas_core_setup_decompositions!}}}

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

      call simulation_clock_init(core_clock, configs, ierr)

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

   !-----------------------------------------------------------------------

end module mpas_core

