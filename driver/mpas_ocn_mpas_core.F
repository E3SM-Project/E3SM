! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_core
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains initialization and timestep drivers for 
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module mpas_core

   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units

   use ocn_forward_mode
   use ocn_analysis_mode
   use ocn_init_mode

   contains

!***********************************************************************
!
!  routine mpas_core_init
!
!> \brief   Initialize MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls all initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine mpas_core_init(domain, stream_manager, startTimeStamp)!{{{

      use mpas_grid_types
      use mpas_stream_manager

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp

      type (dm_info) :: dminfo

      integer :: err

      character (len=StrKIND), pointer :: config_ocean_run_mode

      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_config(domain % configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call ocn_forward_mode_init(domain, stream_manager, startTimeStamp)
      else if ( trim(config_ocean_run_mode) == 'analysis' ) then
         call ocn_analysis_mode_init(domain, stream_manager, startTimeStamp)
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         call ocn_init_mode_init(domain, stream_manager, startTimeStamp)
      end if

   end subroutine mpas_core_init!}}}

!***********************************************************************
!
!  routine mpas_core_run
!
!> \brief   Main driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine includes the time-stepping loop, and calls timer
!>  routines to write output and restart files.
!
!-----------------------------------------------------------------------

   subroutine mpas_core_run(domain, stream_manager)!{{{
   
      use mpas_kind_types
      use mpas_grid_types
      use mpas_stream_manager
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager

      character(len=StrKIND), pointer :: config_ocean_run_mode
   
      call mpas_pool_get_config(domain % configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call ocn_forward_mode_run(domain, stream_manager)
      else if ( trim(config_ocean_run_mode) == 'analysis' ) then
         call ocn_analysis_mode_run(domain, stream_manager)
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         call ocn_init_mode_run(domain, stream_manager)
      end if

   end subroutine mpas_core_run!}}}

   subroutine mpas_core_finalize(domain, stream_manager)!{{{
   
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      integer :: ierr

      character(len=StrKIND), pointer :: config_ocean_run_mode

      call mpas_pool_get_config(domain % configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call ocn_forward_mode_finalize(domain, stream_manager)
      else if (trim(config_ocean_run_mode) == 'analysis' ) then
         call ocn_analysis_mode_finalize(domain, stream_manager)
      else if (trim(config_ocean_run_mode) == 'init' ) then
         call ocn_init_mode_finalize(domain, stream_manager)
      end if

   end subroutine mpas_core_finalize!}}}

!***********************************************************************
!
!  routine mpas_core_setup_packages
!
!> \brief   Package setup routine
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

      use ocn_analysis_driver

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      integer, intent(out) :: ierr

      integer :: err_tmp

      logical, pointer :: forwardModeActive, analysisModeActive, initModeActive
      logical, pointer :: thicknessFilterActive
      logical, pointer :: splitTimeIntegratorActive
      logical, pointer :: bulkForcingActive
      logical, pointer :: frazilIceActive
      logical, pointer :: inSituEOSActive

      logical, pointer :: config_use_freq_filtered_thickness
      logical, pointer :: config_frazil_ice_formation
      character (len=StrKIND), pointer :: config_time_integrator, config_forcing_type
      character (len=StrKIND), pointer :: config_ocean_run_mode, config_pressure_gradient_type

      ! Get Packages
      call mpas_pool_get_package(packagePool, 'forwardModeActive', forwardModeActive)
      call mpas_pool_get_package(packagePool, 'analysisModeActive', analysisModeActive)
      call mpas_pool_get_package(packagePool, 'initModeActive', initModeActive)
      call mpas_pool_get_package(packagePool, 'thicknessFilterActive', thicknessFilterActive)
      call mpas_pool_get_package(packagePool, 'splitTimeIntegratorActive', splitTimeIntegratorActive)
      call mpas_pool_get_package(packagePool, 'bulkForcingActive', bulkForcingActive)
      call mpas_pool_get_package(packagePool, 'frazilIceActive', frazilIceActive)
      call mpas_pool_get_package(packagePool, 'inSituEOSActive', inSituEOSActive)

      call mpas_pool_get_config(configPool, 'config_ocean_run_mode', config_ocean_run_mode)

      ierr = 0

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         forwardModeActive = .true.

         call mpas_pool_get_config(configPool, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
         call mpas_pool_get_config(configPool, 'config_time_integrator', config_time_integrator)
         call mpas_pool_get_config(configPool, 'config_forcing_type', config_forcing_type)
         call mpas_pool_get_config(configPool, 'config_frazil_ice_formation', config_frazil_ice_formation)
         call mpas_pool_get_config(configPool, 'config_pressure_gradient_type', config_pressure_gradient_type)

         if (config_use_freq_filtered_thickness) then
            thicknessFilterActive = .true.
         end if

         if (config_time_integrator == trim('split_explicit') &
            .or. config_time_integrator == trim('unsplit_explicit') ) then

            splitTimeIntegratorActive = .true.
         end if

         if (config_frazil_ice_formation) then
            frazilIceActive = .true.
         end if

         if (config_pressure_gradient_type.eq.'Jacobian_from_TS') then
            inSituEOSActive = .true.
         end if

         call ocn_analysis_setup_packages(configPool, packagePool, err_tmp)
         ierr = ior(ierr, err_tmp)
      else if (trim(config_ocean_run_mode) == 'analysis' ) then
         analysisModeActive = .true.
         call ocn_analysis_setup_packages(configPool, packagePool, ierr)
      else if (trim(config_ocean_run_mode) == 'init' ) then
         initModeActive = .true.
         call ocn_init_validate_configuration(configPool, packagePool, ierr)
      end if

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
   subroutine mpas_core_setup_clock(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      character(len=StrKIND), pointer :: config_ocean_run_mode

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call ocn_forward_mode_simulation_clock_init(core_clock, configs, ierr)
      else if ( trim(config_ocean_run_mode) == 'analysis' ) then
         call ocn_analysis_mode_simulation_clock_init(core_clock, configs, ierr)
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         call ocn_init_mode_simulation_clock_init(core_clock, configs, ierr)
      end if

   end subroutine mpas_core_setup_clock!}}}

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
   subroutine mpas_core_get_mesh_stream(configs, stream, ierr)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      character(len=*), intent(out) :: stream
      integer, intent(out) :: ierr

      logical, pointer :: config_do_restart
      character(len=StrKIND), pointer :: config_ocean_run_mode

      ierr = 0

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' .or. trim(config_ocean_run_mode) == 'analysis' ) then
         call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

         if (.not. associated(config_do_restart)) then
            ierr = 1
         else if (config_do_restart) then
            write(stream,'(a)') 'restart'
         else
            write(stream,'(a)') 'input'
         end if
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         write(stream, '(a)') 'input_init'
      end if

   end subroutine mpas_core_get_mesh_stream!}}}


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
   subroutine mpas_core_setup_decompositions(ierr)!{{{

      use mpas_decomp

      implicit none

      integer, intent(out) :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompositions)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompositions, 'uniform', decompFunc, iErr)

   end subroutine mpas_core_setup_decompositions!}}}

end module mpas_core

! vim: foldmethod=marker
