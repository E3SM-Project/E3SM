! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module li_core_interface

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_constants
   use mpas_io_units
   use mpas_attlist
   use mpas_abort, only : mpas_dmpar_global_abort
   use li_core
   use li_analysis_driver

   public

   contains
   !***********************************************************************
   !
   !  routine li_setup_core
   !
   !> \brief   Land ice core setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a core_type
   !>  for the land ice core.
   !
   !-----------------------------------------------------------------------
   subroutine li_setup_core(core)!{{{
      type (core_type), pointer :: core

      core % core_init => li_core_init
      core % core_run => li_core_run
      core % core_finalize => li_core_finalize
      core % define_packages => li_define_packages
      core % setup_packages => li_setup_packages
      core % setup_decompositions => li_setup_decompositions
      core % setup_clock => li_setup_clock
      core % get_mesh_stream => li_get_mesh_stream
      core % setup_immutable_streams => li_setup_immutable_streams
      core % setup_derived_dimensions => li_setup_derived_dimensions
      core % setup_decomposed_dimensions => li_setup_decomposed_dimensions
      core % setup_block => li_setup_block
      core % setup_namelist => li_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "../inc/core_variables.inc"

   end subroutine li_setup_core!}}}


   !***********************************************************************
   !
   !  routine li_setup_domain
   !
   !> \brief   Land ice domain setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a domain_type
   !>  for the land ice core.
   !
   !-----------------------------------------------------------------------
   subroutine li_setup_domain(domain)!{{{
      type (domain_type), pointer :: domain

#include "../inc/domain_variables.inc"

   end subroutine li_setup_domain!}}}

!***********************************************************************
!
!  function li_setup_packages
!
!> \brief   Pacakge setup function
!> \author  Doug Jacobsen
!> \date    September 2011
!> \details 
!>  This function is intended to correctly configure the packages for this MPAS
!>   core. It can use any Fortran logic to properly configure packages, and it
!>   can also make use of any namelist options. All variables in the model are
!>   *not* allocated until after this routine is called.
!
!-----------------------------------------------------------------------
   function li_setup_packages(configPool, packagePool, iocontext) result(ierr)

      implicit none
      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      ! Local variables
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_SGH
      logical, pointer :: config_adaptive_timestep_include_DCFL

      logical, pointer :: higherOrderVelocityActive
      logical, pointer :: SIAvelocityActive
      logical, pointer :: hydroActive

      ierr = 0

      call mpas_pool_get_config(configPool, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(configPool, 'config_SGH', config_SGH)

      call mpas_pool_get_package(packagePool, 'SIAvelocityActive', SIAvelocityActive)
      call mpas_pool_get_package(packagePool, 'higherOrderVelocityActive', higherOrderVelocityActive)
      call mpas_pool_get_package(packagePool, 'hydroActive', hydroActive)

      if (trim(config_velocity_solver) == 'sia') then
         SIAvelocityActive = .true.
         write (stdoutUnit,*) 'The SIAVelocity package and associated variables and streams has been enabled because the SIA velocity solver is selected.'
      else
         higherOrderVelocityActive = .true.
         write (stdoutUnit,*) "The 'higherOrderVelocity' package and associated variables have been enabled because a higher-order velocity solver is selected."
      end if

      if (config_SGH) then
         hydroActive = .true.
         write(stdOutUnit,*) "The 'hydro' package and assocated variables have been enabled because 'config_SGH' is set to .true."
      endif

      ! call setup packages in analysis driver
      call li_analysis_setup_packages(configPool, packagePool, ierr) 


   end function li_setup_packages


!***********************************************************************
!
!  function li_setup_clock
!
!> \brief   Pacakge setup function
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
   function li_setup_clock(core_clock, configs) result(ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      call li_simulation_clock_init(core_clock, configs, ierr)

   end function li_setup_clock


   !***********************************************************************
   !
   !  routine li_get_mesh_stream
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
   function li_get_mesh_stream(configs, stream) result(ierr)

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         call mpas_dmpar_global_abort('MPAS-landice: ERROR: config_do_restart is not associated.')
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end function li_get_mesh_stream


   !***********************************************************************
   !
   !  function li_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen
   !> \date    04/08/2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function li_setup_decompositions(decompList) result(ierr)

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

   end function li_setup_decompositions


   !***********************************************************************
   !
   !  function li_setup_block
   !
   !> \brief   Land ice block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  land ice core block.
   !
   !-----------------------------------------------------------------------
   function li_setup_block(block) result(iErr)!{{{
      use mpas_derived_types
      type (block_type), pointer :: block
      integer :: iErr

      iErr = 0
      call li_generate_structs(block, block % structs, block % dimensions, block % packages)
   end function li_setup_block!}}}

#include "../inc/setup_immutable_streams.inc"

#include "../inc/block_dimension_routines.inc"

#include "../inc/define_packages.inc"

#include "../inc/structs_and_variables.inc"

#include "../inc/namelist_call.inc"

#include "../inc/namelist_defines.inc"

end module li_core_interface

