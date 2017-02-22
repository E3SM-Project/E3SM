! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core_interface

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_constants
   use mpas_io_units
   use mpas_attlist
   use test_core

   public

   contains
   !***********************************************************************
   !
   !  routine test_setup_core
   !
   !> \brief   Test core setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a core_type
   !>  for the test core.
   !
   !-----------------------------------------------------------------------
   subroutine test_setup_core(core)!{{{
      type (core_type), pointer :: core

      core % core_init => test_core_init
      core % core_run => test_core_run
      core % core_finalize => test_core_finalize
      core % define_packages => test_define_packages
      core % setup_packages => test_setup_packages
      core % setup_decompositions => test_setup_decompositions
      core % setup_clock => test_setup_clock
      core % setup_log => test_setup_log
      core % get_mesh_stream => test_get_mesh_stream
      core % setup_immutable_streams => test_setup_immutable_streams
      core % setup_derived_dimensions => test_setup_derived_dimensions
      core % setup_decomposed_dimensions => test_setup_decomposed_dimensions
      core % setup_block => test_setup_block
      core % setup_namelist => test_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "inc/core_variables.inc"

   end subroutine test_setup_core!}}}


   !***********************************************************************
   !
   !  routine test_setup_domain
   !
   !> \brief   Test domain setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within a domain_type
   !>  for the test core.
   !
   !-----------------------------------------------------------------------
   subroutine test_setup_domain(domain)!{{{
      type (domain_type), pointer :: domain

#include "inc/domain_variables.inc"

   end subroutine test_setup_domain!}}}


   !***********************************************************************
   !
   !  function test_setup_packages
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
   function test_setup_packages(configPool, packagePool, iocontext) result(ierr)!{{{

      use mpas_derived_types

      implicit none

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      ierr = 0

   end function test_setup_packages!}}}


   !***********************************************************************
   !
   !  routine test_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    04/08/2015
   !> \details 
   !>  This routine is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function test_setup_decompositions(decompList) result(ierr)!{{{

      use mpas_derived_types
      use mpas_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList

      integer :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function test_setup_decompositions!}}}


   !***********************************************************************
   !
   !  function test_setup_clock
   !
   !> \brief   Clock setup routine
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
   function test_setup_clock(core_clock, configs) result(ierr)!{{{

      use mpas_derived_types
      use mpas_io_units

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep
      integer :: local_err

      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time

      ierr = 0

      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=local_err)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
         call mpas_set_timeInterval(timeStep, timeString=config_run_duration, ierr=local_err)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
            if(startTime + runduration /= stopTime) then
               write(stderrUnit,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
         timeStep = stopTime - startTime
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=local_err)
      else
          write(stderrUnit,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if


   end function test_setup_clock!}}}


   !***********************************************************************
   !
   !  function test_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function test_setup_log(logInfo, domain) result(iErr)!{{{

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

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

   end function test_setup_log!}}}


   !***********************************************************************
   !
   !  function test_get_mesh_stream
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
   function test_get_mesh_stream(configs, stream) result(ierr)!{{{

      use mpas_derived_types
      use mpas_pool_routines

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      ierr = 0

      write(stream,'(a)') 'input'

   end function test_get_mesh_stream!}}}


   !***********************************************************************
   !
   !  function test_setup_block
   !
   !> \brief   Test block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  test core block.
   !
   !-----------------------------------------------------------------------
   function test_setup_block(block) result(iErr)!{{{
      use mpas_derived_types
      type (block_type), pointer :: block
      integer :: iErr

      iErr = 0
      call test_generate_structs(block, block % structs, block % dimensions, block % packages)
   end function test_setup_block!}}}

#include "inc/setup_immutable_streams.inc"

#include "inc/block_dimension_routines.inc"

#include "inc/define_packages.inc"

#include "inc/structs_and_variables.inc"

#include "inc/namelist_call.inc"

#include "inc/namelist_defines.inc"

end module test_core_interface

