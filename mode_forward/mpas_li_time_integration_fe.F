! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_integration_fe
!
!> \brief MPAS land ice Forward Euler time integration scheme
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains the Forward Euler time integration scheme
!
!-----------------------------------------------------------------------

module li_time_integration_fe

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_timer
   use mpas_vector_reconstruction
   use li_advection
   use li_calving, only: li_calve_ice, li_restore_calving_front
   use li_thermal, only: li_thermal_solver
   use li_diagnostic_vars
   use li_setup

   implicit none
   private

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_time_integrator_forwardeuler

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains
!***********************************************************************


!***********************************************************************
!
!  routine li_time_integrator_forwardeuler
!
!> \brief   Forward Euler time integration scheme 
!> \author  Matthew Hoffman
!> \date    10 January 2012
!> \details 
!>  This routine performs Forward Euler time integration.
!
!-----------------------------------------------------------------------
   subroutine li_time_integrator_forwardeuler(domain, err)

      use li_subglacial_hydro

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object 

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      integer :: err_tmp

      logical, pointer :: config_restore_calving_front

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_restore_calving_front', config_restore_calving_front)

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! (time level 1 should not be modified.)

!!!      procVertexMaskChanged = 0

! === Update subglacial hydrology  ===========
! It's not clear where the best place to call this should be, so putting it first for now
      call mpas_timer_start("subglacial hydro")
      call li_SGH_solve(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("subglacial hydro")

! === Prepare for advection (including CFL checks) ===========

      call mpas_timer_start("advection prep")
      call prepare_advection(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("advection prep")

! === Implicit column physics (vertical temperature diffusion) ===========
      call mpas_timer_start("calculate vertical therm")
      call li_thermal_solver(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calculate vertical therm")

! === Compute new state for prognostic variables ==================================
      call mpas_timer_start("advect thickness and tracers")
      call advection_solver(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("advect thickness and tracers")

! === Calve ice ========================
      call mpas_timer_start("calve_ice")

      if (config_restore_calving_front) then

         ! restore the calving front to its initial position; calving options are ignored
         call li_restore_calving_front(domain, err_tmp)
         err = ior(err, err_tmp)

      else

         ! ice calving
         call li_calve_ice(domain, err_tmp)
         err = ior(err, err_tmp)

      endif

      call mpas_timer_stop("calve_ice")

! === Calculate diagnostic variables for new state =====================

      ! Some velocity solvers require an initial guess.
      ! If needed, that guess should be inserted into normalVelocity
      ! in time level 2 before calling li_calculate_diagnostic_vars.

      call li_calculate_diagnostic_vars(domain, solveVelo=.true., err=err_tmp)
      err = ior(err, err_tmp)


! === Cleanup & Misc. =============================

      ! === error check
      if (err == 1) then
          write (stderrUnit,*) "An error has occurred in li_time_integrator_forwardeuler."
      endif

   !--------------------------------------------------------------------
   end subroutine li_time_integrator_forwardeuler



!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

!***********************************************************************
!
!  routine prepare_advection
!
!> \brief   Preparation for advection, including CFL diagnostics
!> \author  Matthew Hoffman and William Lipscomb
!> \date    January 2016
!> \details
!>  This routine does preparatory calculations for advection of thickness
!>  and tracers:
!>  (1) Compute layer normal velocities.
!>  (2) Compute the advective CFL limit (and optionally, the diffusive CFL limit).
!>  (3) If config_adaptive_timestep = .true., then set deltat based on CFL info.
!>   These calculations were previously done at the same time as advection, but
!>   now are done at the start of the timestep to support an adaptive time step.
!-----------------------------------------------------------------------

   subroutine prepare_advection(domain, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: scratchPool

      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: layerNormalVelocity

      logical, pointer :: config_print_thickness_advection_info
      logical, pointer :: config_adaptive_timestep
      logical, pointer :: config_adaptive_timestep_include_DCFL

      integer :: &
           allowableAdvecDtProcNumberHere, &
           allowableAdvecDtProcNumber

      real (kind=RKIND) :: &
           allowableAdvecDt, &
           allowableAdvecDtOnProc, &
           allowableAdvecDtAllProcs

      type (MPAS_TimeInterval_type) :: &
           allowableAdvecDtOnProcInterval, &
           allowableAdvecDtAllProcsInterval

      character (len=StrKIND) :: &
           allowableAdvecDtOnProcString, &
           allowableAdvecDtAllProcsString

      integer :: &
           allowableDiffDtProcNumberHere, &
           allowableDiffDtProcNumber

      real (kind=RKIND) :: &
           allowableDiffDt, &
           allowableDiffDtOnProc, &
           allowableDiffDtAllProcs

      type (MPAS_TimeInterval_type) :: &
           allowableDiffDtOnProcInterval, &
           allowableDiffDtAllProcsInterval

      character (len=StrKIND) :: &
           allowableDiffDtOnProcString, &
           allowableDiffDtAllProcsString

      real (kind=RKIND), pointer :: &
           allowableDtACFL, &
           allowableDtDCFL

      real (kind=RKIND), pointer :: deltat ! variable in blocks

      real (kind=RKIND) :: dtSeconds ! local variable

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      allowableAdvecDtAllProcs = 0.0_RKIND
      allowableDiffDtAllProcs = 0.0_RKIND

      dminfo => domain % dminfo

      ! Initialize

      err = 0

      allowableAdvecDtOnProc = 1.0e36 ! set to large number
      allowableDiffDtOnProc = 1.0e36 ! set to large number

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)

         ! compute normal velocities and advective CFL limit for this block

         call li_layer_normal_velocity( &
              meshPool,               &
              normalVelocity,         &
              layerNormalVelocity,    &
              allowableAdvecDt,       &
              err_tmp)

         err = ior(err, err_tmp)

         allowableAdvecDtOnProc = min(allowableAdvecDtOnProc, allowableAdvecDt)

         ! Calculate diffusive CFL timestep, if needed
         ! This used to be only calculated if (config_adaptive_timestep_include_DCFL) but for simplicity,
         ! now it is always calculated.  That allows assessment of the DCFL even when it is not being obeyed
         ! by the adaptive time stepper.  Timer is added here to monitor if this unnecessary calculation becomes significant.
         call mpas_timer_start("calculate apparent diffusivity")
         call li_calculate_apparent_diffusivity(meshPool, velocityPool, scratchPool, geometryPool, allowableDiffDt)
         allowableDiffDtOnProc = min(allowableDiffDtOnProc, allowableDiffDt)
         call mpas_timer_stop("calculate apparent diffusivity")

         ! Note: The ACFL and DCFL timesteps could be calculated in
         ! diagnostic_variables_solve_after_velocity.  In that case, we could also add
         ! variables to store their values, rather than just relying on the values
         ! written to the log files.  However, the current logic only calculates these
         ! values if certain config options are set, so that would need to be dealt with.
         ! If the ACFL and DCFL timesteps are moved to diagnostic_variables_solve_after_velocity,
         ! Then the setting of the timestep value could happen at the beginning of the timestep,
         ! probably in li_timestep rather than here.  That might be cleaner, but I have
         ! persisted with doing it here, because since the calculation of ACFL/DCFL have a lot
         ! in common with the advection calculation, and it seems kind of silly to do those calculations
         ! on the previous time step.  That said, the calculations are pretty cheap.

         block => block % next
      end do


      ! Local advective CFL info
      call mpas_set_timeInterval(allowableAdvecDtOnProcInterval, dt=allowableAdvecDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableAdvecDtOnProcInterval, timeString=allowableAdvecDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)

      if (config_print_thickness_advection_info) then
         write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on advective CFL is (Days_hhh:mmm:sss):   ' // trim(allowableAdvecDtOnProcString)
      endif


      ! Local diffusive CFL info
      ! This used to be only calculated if (config_adaptive_timestep_include_DCFL) but for simplicity,
      ! now it is always calculated.  That allows assessment of the DCFL even when it is not being obeyed
      ! by the adaptive time stepper.
      call mpas_set_timeInterval(allowableDiffDtOnProcInterval, dt=allowableDiffDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableDiffDtOnProcInterval, timeString=allowableDiffDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)

      if (config_print_thickness_advection_info) then
         write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on diffusive CFL is (Days_hhh:mmm:sss):   ' // trim(allowableDiffDtOnProcString)
      endif


      ! These calculation could be (and used to be) restricted to only if we are
      ! printing advection debug information or adaptive timestepping, because
      ! it requires 2 unnecessary MPI communications for each of the ACFL and DCFL
      ! However, those communications are probably small relative to other costs in
      ! the model, so now they always happen.  The timer has been added to allow
      ! assessment of that assumption.
      call mpas_timer_start("calculate global CFL limits")

      ! Determine ACFL limit on all procs
      call mpas_dmpar_min_real(dminfo, allowableAdvecDtOnProc, allowableAdvecDtAllProcs)

      ! Determine which processor has the limiting CFL
      if (allowableAdvecDtOnProc == allowableAdvecDtAllProcs) then
        allowableAdvecDtProcNumberHere = dminfo % my_proc_id
      else
        allowableAdvecDtProcNumberHere = -1
      endif

      call mpas_dmpar_max_int(dminfo, allowableAdvecDtProcNumberHere, allowableAdvecDtProcNumber)
      call mpas_set_timeInterval(allowableAdvecDtAllProcsInterval, dt=allowableAdvecDtAllProcs, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableAdvecDtAllProcsInterval, timeString=allowableAdvecDtAllProcsString, ierr=err_tmp)
      err = ior(err,err_tmp)

      ! Repeat for diffusive CFL
      ! This used to be only calculated if (config_adaptive_timestep_include_DCFL) but for simplicity,
      ! now it is always calculated.  That allows assessment of the DCFL even when it is not being obeyed
      ! by the adaptive time stepper.

      ! Determine DCFL limit on all procs
      call mpas_dmpar_min_real(dminfo, allowableDiffDtOnProc, allowableDiffDtAllProcs)

      ! Determine which processor has the limiting CFL
      if (allowableDiffDtOnProc == allowableDiffDtAllProcs) then
        allowableDiffDtProcNumberHere = dminfo % my_proc_id
      else
        allowableDiffDtProcNumberHere = -1
      endif

      call mpas_dmpar_max_int(dminfo, allowableDiffDtProcNumberHere, allowableDiffDtProcNumber)
      call mpas_set_timeInterval(allowableDiffDtAllProcsInterval, dt=allowableDiffDtAllProcs, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableDiffDtAllProcsInterval, timeString=allowableDiffDtAllProcsString, ierr=err_tmp)
      err = ior(err,err_tmp)

      call mpas_timer_stop("calculate global CFL limits")


      ! Write messages if they are turned on
      if (config_print_thickness_advection_info) then
          write(stdoutUnit,*) '  Maximum allowable time step for all processors based on advective CFL is (Days_hhh:mmm:sss): ' // trim(allowableAdvecDtAllProcsString) //  '  Time step is limited by processor number ', allowableAdvecDtProcNumber
          if (config_adaptive_timestep_include_DCFL) then
             write(stdoutUnit,*) '  Maximum allowable time step for all processors based on diffusive CFL is (Days_hhh:mmm:sss): ' // trim(allowableDiffDtAllProcsString) //  '  Time step is limited by processor number ', allowableDiffDtProcNumber
          endif
      endif


      ! Set adaptive timestep if needed
      if (config_adaptive_timestep) then
         call set_timestep(allowableAdvecDtAllProcs, allowableDiffDtAllProcs, domain % clock, dtSeconds, err_tmp)
         err = ior(err,err_tmp)
         ! Set new value on all blocks
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_array(meshPool, 'deltat', deltat)

            deltat = dtSeconds

            block => block % next
         end do
      end if


      ! Check for CFL error before finishing
      if (dtSeconds > allowableAdvecDtOnProc) then
         write(stderrUnit,*) 'ERROR: Advective CFL violation on this processor.  Maximum allowable time step for this processor is (Days_hhh:mmm:sss): ' // trim(allowableAdvecDtOnProcString)
         err = ior(err,1)
      endif

      ! Local diffusive CFL info
      if ( (config_adaptive_timestep_include_DCFL) .and. (dtSeconds > allowableDiffDtOnProc) ) then
         write(stderrUnit,*) 'WARNING: Diffusive CFL violation on this processor.  Maximum allowable time step for this processor is (Days_hhh:mmm:sss): ' // trim(allowableDiffDtOnProcString)
      endif

      if (err > 0) then
           write(stderrUnit,*) 'Error in calculating thickness advection  (possibly CFL violation)'
      endif


      ! set CFL variables if they have been calculated - every block should be set to the same value!
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'allowableDtACFL', allowableDtACFL)
         allowableDtACFL = allowableAdvecDtAllProcs
         call mpas_pool_get_array(meshPool, 'allowableDtDCFL', allowableDtDCFL)
         allowableDtDCFL = allowableDiffDtAllProcs

         block => block % next
      end do


      ! Halo updates
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'layerNormalVelocity')
      call mpas_timer_stop("halo updates")

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in prepare_advection."
      endif

   !--------------------------------------------------------------------
    end subroutine prepare_advection

!***********************************************************************
!
!  routine advection_solver
!
!> \brief   Advect thickness and tracers 
!> \author  Matthew Hoffman and William Lipscomb
!> \date    September 2013; revised December 2015
!> \details
!>  This routine advects thickness and tracers as part of forward Euler
!>  time integration.
!>  Note: This routine replaces much of the old subroutines calculate_tendencies
!>        and update_prognostics. The CFL diagnostics that were previously
!>        in calculate tendencies are now in prepare_advection.
!-----------------------------------------------------------------------

   subroutine advection_solver(domain, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: scratchPool

      real (kind=RKIND), pointer :: deltat ! variable in blocks
      real (kind=RKIND), dimension(:), pointer :: thicknessOld
      real (kind=RKIND), dimension(:), pointer :: thicknessNew
      real (kind=RKIND), dimension(:), pointer :: thickness

      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), dimension(:,:), pointer :: waterfrac
      real (kind=RKIND), dimension(:,:), pointer :: enthalpy

      integer, pointer :: nCells

      character (len=StrKIND), pointer :: config_thickness_advection
      character (len=StrKIND), pointer :: config_tracer_advection

      logical, pointer :: config_print_thickness_advection_info

      !TODO - Replace masktmp with a scratch field?
      integer, dimension(:), allocatable :: masktmp  ! Temporary mask for assessing new thickness field

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      dminfo => domain % dminfo

      ! Halo updates
      ! Note: The layer thickness and tracers must be up to date in halos before calling the advection subroutines.
      !       Currently, layerThicknessEdge is updated at the end of subroutine li_calculate_diagnostic_vars.
      !       The thermal tracers (temperature, waterfrac, enthalpy) are updated at the end of li_thermal_solver.
      !       But thickness (which is used by subroutine li_advection_thickness_tracers) needs an update here.
      !        We use timeLevel=2 to get the value from the previous timestep.

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_timer_stop("halo updates")

      ! ===
      ! === Advect thickness and tracers
      ! ===

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)

         if (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'fo') then

            ! Note: This subroutine requires that thickness and tracers are correct in halos

            if (config_print_thickness_advection_info) then
               write(stdoutUnit,*) 'Advect thickness and tracers, dt =', deltat
            endif

            call li_advection_thickness_tracers(&
                 deltat,                 &
                 meshPool,               &
                 velocityPool,           &
                 geometryPool,           &
                 thermalPool,            &
                 scratchPool,            &
                 err_tmp,                &
                 advectTracersIn = .true.)

            err = ior(err,err_tmp)

         elseif (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'none') then

            if (config_print_thickness_advection_info) then
               write(stdoutUnit,*) 'Advect thickness (but not tracers), dt =', deltat
            endif

            call li_advection_thickness_tracers(&
                 deltat,                 &
                 meshPool,               &
                 velocityPool,           &
                 geometryPool,           &
                 thermalPool,            &
                 scratchPool,            &
                 err_tmp,                &
                 advectTracersIn = .false.)

            err = ior(err,err_tmp)

         endif

         block => block % next
      end do

      ! Reset negative thicknesses to zero if needed

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel=1)

         
         allocate( masktmp(nCells + 1) )
         masktmp = 0

         where (thickness < 0.0_RKIND)
            masktmp = 1
            thickness = 0.0_RKIND
         end where

         if (config_print_thickness_advection_info) then

            if (sum(masktmp) > 0) then
               write(stdoutUnit,*) '  Cells with negative thickness (set to 0):', sum(masktmp)
            endif

            ! Count how many cells have ice.
            masktmp = 0
            where (thickness > 0.0_RKIND)
               masktmp = 1
            end where
            write(stdoutUnit,*) '  Cells with nonzero thickness:', sum(masktmp)

         endif

         deallocate(masktmp)

         block => block % next
      end do

      ! Halo updates
      call mpas_timer_start("halo updates")

      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_dmpar_field_halo_exch(domain, 'temperature')
      call mpas_dmpar_field_halo_exch(domain, 'waterfrac')
      call mpas_dmpar_field_halo_exch(domain, 'enthalpy')

      call mpas_timer_stop("halo updates")

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in advection_solver."
      endif

   !--------------------------------------------------------------------
    end subroutine advection_solver

!***********************************************************************
!
!  routine set_timestep
!
!> \brief   Adjusts the time step based on the CFL condition.
!> \author  Matthew Hoffman
!> \date    23 Jan 2014
!> \details
!>  This routine sdjusts the time step based on the CFL condition.
!
!-----------------------------------------------------------------------
   subroutine set_timestep(allowableAdvecDt, allowableDiffDt, clock, dtSeconds, err)
      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: allowableAdvecDt
      real (kind=RKIND), intent(in) :: allowableDiffDt
      type (MPAS_Clock_type), intent(in) :: clock

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(out) :: dtSeconds
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_adaptive_timestep_include_DCFL
      real (kind=RKIND), pointer :: config_adaptive_timestep_CFL_fraction
      real (kind=RKIND), pointer :: config_max_adaptive_timestep
      real (kind=RKIND), pointer :: config_min_adaptive_timestep
      type (MPAS_Time_type) :: nextForceTime, currTime
      type (MPAS_TimeInterval_type) :: intervalToNextForceTime
      real (kind=RKIND) :: secondsToNextForceTime
      real (kind=RKIND) :: allowableDt
      real (kind=RKIND) :: proposedDt
      integer :: err_tmp

      err = 0
      err_tmp = 0


      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_CFL_fraction', config_adaptive_timestep_CFL_fraction)
      call mpas_pool_get_config(liConfigs, 'config_max_adaptive_timestep', config_max_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_min_adaptive_timestep', config_min_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      if (config_adaptive_timestep_include_DCFL) then
         allowableDt = min(allowableAdvecDt, allowableDiffDt)
      else
         allowableDt = allowableAdvecDt
      endif
      proposedDt = min(allowableDt * config_adaptive_timestep_CFL_fraction, config_max_adaptive_timestep)

      ! Check if we need to force a timestep length to hit the target interval
      currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
      !print *, 'curr', currTime % t % YR, currTime % t % basetime % S, currTime % t % basetime % Sn, currTime % t % basetime % Sd
      err = ior(err,err_tmp)
      nextForceTime = mpas_alarm_get_next_ring_time(clock, 'adaptiveTimestepForceInterval')
      !print *, 'ring', nextForceTime % t % YR, nextForceTime % t % basetime % S, nextForceTime % t % basetime % Sn, nextForceTime % t % basetime % Sd
      intervalToNextForceTime = nextForceTime - currTime
      !print *, 'int', intervalToNextForceTime % ti % YR, intervalToNextForceTime % ti % MM, intervalToNextForceTime % ti % basetime % S, intervalToNextForceTime % ti % basetime % Sn, intervalToNextForceTime % ti % basetime % Sd
      ! Due to roundoff errors, we might be just shy of the desired time.
      ! To avoid this, add one to the numerator of the fractional seconds to
      ! make sure we get pushed over the edge.  The way ESMF does fractional
      ! seconds, this means we get the desired interval to better than 1 part per 100 million seconds
      ! Note that even though this is a *very* tiny fudge factor, it does not
      ! affect conservation within MPAS-LI, but it could have a very, very tiny
      ! effect on a climate model that thinks we ran for, say, 10 years, but we
      ! actually ran for 10 years +/- 1e-8 seconds.
      intervalToNextForceTime % ti % basetime % Sn = intervalToNextForceTime % ti % basetime % Sn + 1
      call mpas_get_timeInterval(intervalToNextForceTime, dt=secondsToNextForceTime, ierr=err_tmp)
      err = ior(err,err_tmp)
      !print *, proposedDt, secondsToNextForceTime
      dtSeconds = min(proposedDt, secondsToNextForceTime)

      write(stdOutUnit,*) '  Setting time step (days) to:', dtSeconds / (86400.0)
      if (dtSeconds < config_min_adaptive_timestep) then
         write(stdErrUnit,*) 'ERROR: New deltat is less than config_min_adaptive_timestep.'
         err = 1
      endif

   !--------------------------------------------------------------------
   end subroutine set_timestep



end module li_time_integration_fe


