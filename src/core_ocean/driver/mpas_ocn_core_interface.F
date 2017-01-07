! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_core_interface

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_constants
   use mpas_io_units
   use ocn_core
   use mpas_attlist

   use ocn_forward_mode
   use ocn_analysis_mode

   private

   public :: ocn_setup_core, ocn_setup_domain

   contains
   !***********************************************************************
   !
   !  routine ocn_setup_core
   !
   !> \brief   Ocean core setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a core_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_setup_core(core)!{{{
      type (core_type), pointer :: core

      core % core_init => ocn_core_init
      core % core_run => ocn_core_run
      core % core_finalize => ocn_core_finalize
      core % define_packages => ocn_define_packages
      core % setup_packages => ocn_setup_packages
      core % setup_decompositions => ocn_setup_decompositions
      core % setup_clock => ocn_setup_clock
      core % get_mesh_stream => ocn_get_mesh_stream
      core % setup_immutable_streams => ocn_setup_immutable_streams
      core % setup_derived_dimensions => ocn_setup_derived_dimensions
      core % setup_decomposed_dimensions => ocn_setup_decomposed_dimensions
      core % setup_block => ocn_setup_block
      core % setup_namelist => ocn_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "../inc/core_variables.inc"

   end subroutine ocn_setup_core!}}}


   !***********************************************************************
   !
   !  routine ocn_setup_domain
   !
   !> \brief   Ocean domain setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a domain_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_setup_domain(domain)!{{{
      type (domain_type), pointer :: domain

#include "../inc/domain_variables.inc"

   end subroutine ocn_setup_domain!}}}


   !***********************************************************************
   !
   !  function ocn_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    03/12/2015
   !> \details 
   !>  This function is intended to correctly configure the packages for this MPAS
   !>   core. It can use any Fortran logic to properly configure packages, and it
   !>   can also make use of any namelist options. All variables in the model are
   !>   *not* allocated until after this routine is called.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_packages(configPool, packagePool, iocontext) result(ierr)!{{{

      use ocn_analysis_driver

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer :: ierr

      integer :: err_tmp

      logical, pointer :: forwardModeActive, analysisModeActive
      logical, pointer :: thicknessFilterActive
      logical, pointer :: splitTimeIntegratorActive
      logical, pointer :: surfaceRestoringActive
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
      call mpas_pool_get_package(packagePool, 'thicknessFilterActive', thicknessFilterActive)
      call mpas_pool_get_package(packagePool, 'splitTimeIntegratorActive', splitTimeIntegratorActive)
      call mpas_pool_get_package(packagePool, 'surfaceRestoringActive', surfaceRestoringActive)
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

         if (config_forcing_type == trim('restoring')) then
            surfaceRestoringActive = .true.
         else if (config_forcing_type == trim('bulk')) then
            bulkForcingActive = .true.
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
      end if

   end function ocn_setup_packages!}}}


   !***********************************************************************
   !
   !  routine ocn_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    04/08/2015
   !> \details 
   !>  This routine is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_decompositions(decompList) result(ierr)!{{{

      use mpas_derived_types
      use mpas_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList

      integer :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, iErr)

      if ( iErr == MPAS_DECOMP_NOERR ) then
         iErr = 0
      end if

   end function ocn_setup_decompositions!}}}


   !***********************************************************************
   !
   !  function ocn_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this function is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This function is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_clock(core_clock, configs) result(ierr)!{{{

      use mpas_derived_types

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_ocean_run_mode

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         ierr = ocn_forward_mode_setup_clock(core_clock, configs)
      else if ( trim(config_ocean_run_mode) == 'analysis' ) then
         ierr = ocn_analysis_mode_setup_clock(core_clock, configs)
      end if


   end function ocn_setup_clock!}}}


   !***********************************************************************
   !
   !  function ocn_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This function returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function ocn_get_mesh_stream(configs, stream) result(ierr)!{{{

      use mpas_derived_types
      use mpas_pool_routines

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart
      character(len=StrKIND), pointer :: config_ocean_run_mode

      ierr = 0

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' .or. trim(config_ocean_run_mode) == 'analysis' ) then
         write(stream,'(a)') 'mesh'
      end if

   end function ocn_get_mesh_stream!}}}


   !***********************************************************************
   !
   !  function ocn_setup_block
   !
   !> \brief   Ocean block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  ocean core block.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_block(block) result(iErr)!{{{
      use mpas_derived_types
      type (block_type), pointer :: block
      integer :: iErr

      iErr = 0
      call ocn_generate_structs(block, block % structs, block % dimensions, block % packages)
   end function ocn_setup_block!}}}

#include "../inc/setup_immutable_streams.inc"

#include "../inc/block_dimension_routines.inc"

#include "../inc/define_packages.inc"

#include "../inc/structs_and_variables.inc"

#include "../inc/namelist_call.inc"

#include "../inc/namelist_defines.inc"

end module ocn_core_interface

