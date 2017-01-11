! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_core_interface

   use mpas_attlist

   contains


   !***********************************************************************
   !
   !  routine init_atm_setup_core
   !
   !> \brief   Init atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use init_atm_core, only : init_atm_core_init, init_atm_core_run, init_atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => init_atm_core_init
      core % core_run => init_atm_core_run
      core % core_finalize => init_atm_core_finalize
      core % define_packages => init_atm_define_packages
      core % setup_packages => init_atm_setup_packages
      core % setup_decompositions => init_atm_setup_decompositions
      core % setup_clock => init_atm_setup_clock
      core % get_mesh_stream => init_atm_get_mesh_stream
      core % setup_immutable_streams => init_atm_setup_immutable_streams
      core % setup_derived_dimensions => init_atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => init_atm_setup_decomposed_dimensions
      core % setup_block => init_atm_setup_block
      core % setup_namelist => init_atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "inc/core_variables.inc"

write(0,*) ''
#ifdef SINGLE_PRECISION
write(0,'(a)') 'Using default single-precision reals'
#else
write(0,'(a)') 'Using default double-precision reals'
#endif
write(0,*) ''

   end subroutine init_atm_setup_core


   !***********************************************************************
   !
   !  routine init_atm_setup_domain
   !
   !> \brief   Init atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "inc/domain_variables.inc"

   end subroutine init_atm_setup_domain


   !***********************************************************************
   !
   !  function init_atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  init_atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_packages(configs, packages, iocontext) result(ierr)

      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package
      use mpas_io_units, only : stderrUnit

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      logical, pointer :: initial_conds, sfc_update
      logical, pointer :: gwd_stage_in, vertical_stage_in, vertical_stage_out, met_stage_in, met_stage_out
      logical, pointer :: config_native_gwd_static, config_static_interp, config_vertical_grid, config_met_interp
      integer, pointer :: config_init_case


      ierr = 0

      call mpas_pool_get_config(configs, 'config_init_case', config_init_case)
      call mpas_pool_get_config(configs, 'config_static_interp', config_static_interp)
      call mpas_pool_get_config(configs, 'config_native_gwd_static', config_native_gwd_static)
      call mpas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(configs, 'config_met_interp', config_met_interp)

      nullify(initial_conds)
      call mpas_pool_get_package(packages, 'initial_condsActive', initial_conds)

      nullify(sfc_update)
      call mpas_pool_get_package(packages, 'sfc_updateActive', sfc_update)

      nullify(gwd_stage_in)
      call mpas_pool_get_package(packages, 'gwd_stage_inActive', gwd_stage_in)

      nullify(vertical_stage_in)
      call mpas_pool_get_package(packages, 'vertical_stage_inActive', vertical_stage_in)

      nullify(vertical_stage_out)
      call mpas_pool_get_package(packages, 'vertical_stage_outActive', vertical_stage_out)

      nullify(met_stage_in)
      call mpas_pool_get_package(packages, 'met_stage_inActive', met_stage_in)

      nullify(met_stage_out)
      call mpas_pool_get_package(packages, 'met_stage_outActive', met_stage_out)

      if (.not. associated(initial_conds) .or. &
          .not. associated(sfc_update) .or. &
          .not. associated(gwd_stage_in) .or. &
          .not. associated(vertical_stage_in) .or. &
          .not. associated(vertical_stage_out) .or. &
          .not. associated(met_stage_in) .or. &
          .not. associated(met_stage_out)) then
         write(stderrUnit,*) '********************************************************************************'
         write(stderrUnit,*) '* Error while setting up packages for init_atmosphere core.'
         write(stderrUnit,*) '********************************************************************************'
         ierr = 1
         return
      end if

      if (config_init_case == 8) then
         initial_conds = .false.
         sfc_update = .true.
      else
         initial_conds = .true.
         sfc_update = .false.
      end if

      if (config_init_case == 7) then
         gwd_stage_in = (config_native_gwd_static .and. .not. config_static_interp)
         vertical_stage_in = (config_vertical_grid .and. .not. config_static_interp)
         vertical_stage_out = (config_vertical_grid .and. .not. config_met_interp)
         met_stage_in = (config_met_interp .and. .not. config_vertical_grid)
         met_stage_out = config_met_interp
      else if (config_init_case == 8) then
         gwd_stage_in = .false.
         vertical_stage_in = .true.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.
      else
         gwd_stage_in = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .true.
      end if

   end function init_atm_setup_packages


   !***********************************************************************
   !
   !  function init_atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
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
   function init_atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_timekeeping, only : mpas_set_time, mpas_set_timeInterval, mpas_create_clock
      use mpas_derived_types, only : MPAS_Clock_type, MPAS_Time_type, MPAS_Timeinterval_type, mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config
      use mpas_kind_types, only : StrKIND

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_start_time, config_stop_time
      integer, pointer :: config_fg_interval

      type (MPAS_Time_type)         :: start_time, stop_time
      type (MPAS_TimeInterval_type) :: dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_fg_interval', config_fg_interval)

      call mpas_set_time(start_time, dateTimeString=trim(config_start_time))
      call mpas_set_time(stop_time, dateTimeString=trim(config_stop_time))
      call mpas_set_timeInterval(dt, S=config_fg_interval)
      call mpas_create_clock(core_clock, start_time, dt, stopTime=stop_time)

   end function init_atm_setup_clock


   !***********************************************************************
   !
   !  function init_atm_get_mesh_stream
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
   function init_atm_get_mesh_stream(configs, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      ierr = 0

      write(stream,'(a)') 'input'

   end function init_atm_get_mesh_stream


   !***********************************************************************
   !
   !  function init_atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp 

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

   end function init_atm_setup_decompositions


   !***********************************************************************
   !
   !  function init_atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  init atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call init_atm_generate_structs(block, block % structs, block % dimensions, block % packages)

   end function init_atm_setup_block


#include "inc/setup_immutable_streams.inc"

#include "inc/block_dimension_routines.inc"

#include "inc/define_packages.inc"

#include "inc/structs_and_variables.inc"

#include "inc/namelist_call.inc"

#include "inc/namelist_defines.inc"

end module init_atm_core_interface
