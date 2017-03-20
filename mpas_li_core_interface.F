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
   use mpas_attlist
   use mpas_io_units
   use mpas_log
   use li_core

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
      core % setup_log => li_setup_log
      core % get_mesh_stream => li_get_mesh_stream
      core % setup_immutable_streams => li_setup_immutable_streams
      core % setup_derived_dimensions => li_setup_derived_dimensions
      core % setup_decomposed_dimensions => li_setup_decomposed_dimensions
      core % setup_block => li_setup_block
      core % setup_namelist => li_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "inc/core_variables.inc"

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

#include "inc/domain_variables.inc"

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

      ierr = 0

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
   !  function li_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function li_setup_log(logInfo, domain) result(iErr)!{{{

      use mpas_derived_types
      use mpas_log

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo  !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain  !< domain object to provide info for setting up log manager
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! ATTENTION: DISABLING LOG MANAGER FOR NOW
      logInfo % outputLog % isActive = .false.
      logInfo % errorLog % isActive = .false.

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

   end function li_setup_log!}}}


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
         call mpas_log_write('ERROR: config_do_restart is not associated.', MPAS_LOG_CRIT)
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

#include "inc/setup_immutable_streams.inc"

#include "inc/block_dimension_routines.inc"

#include "inc/define_packages.inc"

#include "inc/structs_and_variables.inc"

#include "inc/namelist_call.inc"

#include "inc/namelist_defines.inc"

end module li_core_interface

