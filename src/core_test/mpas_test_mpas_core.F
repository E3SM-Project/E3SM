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

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  routine mpas_core_init
   !
   !> \brief   MPAS Core Initialization Routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This routine performs the necessary initialization of a core. This can
   !>   involve things like setting up coeffiecients for remapping and advection,
   !>   but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_init(domain, startTimeStamp)!{{{
   
      use mpas_configure
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      integer :: ierr

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock


      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

   end subroutine mpas_core_init!}}}


   !***********************************************************************
   !
   !  routine mpas_core_run
   !
   !> \brief   MPAS Core Run Routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This routine performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_run(domain)!{{{
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain


   end subroutine mpas_core_run!}}}

   !***********************************************************************
   !
   !  routine mpas_core_finalize
   !
   !> \brief   MPAS Core Finalization Routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This routine finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_finalize(domain)!{{{
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr
 
     call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize!}}}


   !***********************************************************************
   !
   !  routine mpas_core_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
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

   end subroutine mpas_core_setup_packages!}}}


   !***********************************************************************
   !
   !  routine mpas_core_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
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

      type (MPAS_Time_Type) :: startTime, stopTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep
      integer :: local_err

      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time

      ierr = 0

      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=local_err)
      call mpas_set_timeInterval(timeStep, timeString=config_run_duration, ierr=local_err)

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

   end subroutine mpas_core_setup_clock!}}}


   !***********************************************************************
   !
   !  routine mpas_core_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
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

      ierr = 0

      write(stream,'(a)') 'input'

   end subroutine mpas_core_get_mesh_stream!}}}

   !***********************************************************************
   !
   !  routine mpas_core_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
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



end module mpas_core
