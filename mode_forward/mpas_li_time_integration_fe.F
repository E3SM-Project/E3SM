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
   use li_velocity, only: li_velocity_solve
   use li_tendency   !TODO - li_tendency to be deprecated
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
   subroutine li_time_integrator_forwardeuler(domain, deltat, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: deltat  !< Input: time step

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


! === Implicit column physics (vertical temperature diffusion) ===========
      call mpas_timer_start("calculate vertical therm")
      call li_thermal_solver(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calculate vertical therm")

      !WHL - temporary option to choose between old and new advection schemes
      !TODO - Remove the old thickness advection scheme

   if (new_thickness_advection) then

! === Compute new state for prognostic variables ==================================
      call mpas_timer_start("advect thickness and tracers")
      call advection_solver(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("advect thickness and tracers")

   else   ! old thickness advection

! === Calculate tendencies ========================
      call mpas_timer_start("calculate tendencies")
      call calculate_tendencies(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calculate tendencies")

! === Compute new state for prognostic variables ==================================
! (once implicit column physics are added (i.e. temp diffusion), these calculations will need to be adjusted to apply to the new values as needed)
      call mpas_timer_start("calc. new prognostic vars")
      call update_prognostics(domain, deltat, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calc. new prognostic vars")

   endif   ! new_thickness_advection

! === Calve ice ========================                                                                                                                                             
      call mpas_timer_start("calve_ice")

      if (config_restore_calving_front) then

         ! restore the calving front to its initial position; calving options are ignored
         call li_restore_calving_front(domain, err_tmp)
         err = ior(err, err_tmp)

      else

         ! ice calving
         call li_calve_ice(domain, deltat, err_tmp)
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
!  routine calculate_tendencies
!
!> \brief   Calculates tendencies 
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine calculates tendencies of thickness
!>  to use with the Forward Euler time integrator.
!
!-----------------------------------------------------------------------

   subroutine calculate_tendencies(domain, dtSeconds, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: dtSeconds

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
      type (mpas_pool_type), pointer :: tendencyPool
      type (mpas_pool_type), pointer :: scratchPool

      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessTendency
      type (field2DReal), pointer :: layerThicknessTendencyField

      type (field2DReal), pointer :: layerThicknessEdgeField

      logical, pointer :: config_print_thickness_advection_info
      logical, pointer :: config_adaptive_timestep
      logical, pointer :: config_adaptive_timestep_include_DCFL

      integer :: allowableAdvecDtProcNumberHere, allowableAdvecDtProcNumber
      real (kind=RKIND) :: allowableAdvecDt, allowableAdvecDtOnProc, allowableAdvecDtAllProcs
      type (MPAS_TimeInterval_type) :: allowableAdvecDtOnProcInterval, allowableAdvecDtAllProcsInterval
      character (len=StrKIND) :: allowableAdvecDtOnProcString, allowableAdvecDtAllProcsString
      integer :: allowableDiffDtProcNumberHere,  allowableDiffDtProcNumber  
      real (kind=RKIND) :: allowableDiffDt, allowableDiffDtOnProc, allowableDiffDtAllProcs
      type (MPAS_TimeInterval_type) :: allowableDiffDtOnProcInterval, allowableDiffDtAllProcsInterval
      character (len=StrKIND) :: allowableDiffDtOnProcString, allowableDiffDtAllProcsString
      real (kind=RKIND), pointer :: allowableDtACFL, allowableDtDCFL

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      allowableAdvecDtAllProcs = 0.0_RKIND
      allowableDiffDtAllProcs = 0.0_RKIND

      dminfo => domain % dminfo

      ! Halo updates: none currently needed
      ! Note: Currently (Nov. 2015), layerThicknessEdge (which is used by subroutine li_tendency_thickness) is updated 
      !       at the end of subroutine li_calculate_diagnostic_vars.
      !   

      ! ===
      ! === Thickness and tracer tendencies
      ! ===
      allowableAdvecDtOnProc = 1.0e36 ! set to large number
      allowableDiffDtOnProc = 1.0e36 ! set to large number

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'tendency', tendencyPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_array(tendencyPool, 'layerThicknessTendency', layerThicknessTendency)

         ! Calculate thickness tendency using state at time n =========

         call li_tendency_thickness(&
              meshPool,               &
              velocityPool,           &
              geometryPool,           & 
              layerThicknessTendency, &
              dtSeconds,              &
              dminfo,                 &
              allowableAdvecDt,       &
              err_tmp)

         err = ior(err,err_tmp)
         allowableAdvecDtOnProc = min(allowableAdvecDtOnProc, allowableAdvecDt)

         ! Calculate diffusive CFL timestep, if needed
         if (config_adaptive_timestep_include_DCFL) then
            call li_calculate_apparent_diffusivity(meshPool, velocityPool, scratchPool, geometryPool, allowableDiffDt)
            allowableDiffDtOnProc = min(allowableDiffDtOnProc, allowableDiffDt)
         endif

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

      ! Now that we have exited the block loop, do any needed halo updates.
      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'tendency', tendencyPool)
      call mpas_pool_get_field(tendencyPool, 'layerThicknessTendency', layerThicknessTendencyField)
      call mpas_dmpar_exch_halo_field(layerThicknessTendencyField)
      call mpas_timer_stop("halo updates")


      ! Local advective CFL info
      call mpas_set_timeInterval(allowableAdvecDtOnProcInterval, dt=allowableAdvecDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableAdvecDtOnProcInterval, timeString=allowableAdvecDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)

      if (config_print_thickness_advection_info) then
         write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on advective CFL is (Days_hhh:mmm:sss):   ' // trim(allowableAdvecDtOnProcString)
      endif

      ! Local diffusive CFL info
      if (config_adaptive_timestep_include_DCFL) then
          call mpas_set_timeInterval(allowableDiffDtOnProcInterval, dt=allowableDiffDtOnProc, ierr=err_tmp)
          err = ior(err,err_tmp)
          call mpas_get_timeInterval(allowableDiffDtOnProcInterval, timeString=allowableDiffDtOnProcString, ierr=err_tmp)
          err = ior(err,err_tmp)

          if (config_print_thickness_advection_info) then
             write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on diffusive CFL is (Days_hhh:mmm:sss):   ' // trim(allowableDiffDtOnProcString)
          endif
      endif


      ! If we are printing advection debug information or adaptive timestepping,
      ! then find out what the global CFL limit is.  Don't do this otherwise because
      ! it requires 2 unnecessary MPI communications.
      if (config_print_thickness_advection_info .or. config_adaptive_timestep) then
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
          if (config_adaptive_timestep_include_DCFL) then
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
          endif
      endif

      ! Write messages if they are turned on
      if (config_print_thickness_advection_info) then
          write(stdoutUnit,*) '  Maximum allowable time step for all processors based on advective CFL is (Days_hhh:mmm:sss): ' // trim(allowableAdvecDtAllProcsString) //  '  Time step is limited by processor number ', allowableAdvecDtProcNumber
          if (config_adaptive_timestep_include_DCFL) then
             write(stdoutUnit,*) '  Maximum allowable time step for all processors based on diffusive CFL is (Days_hhh:mmm:sss): ' // trim(allowableDiffDtAllProcsString) //  '  Time step is limited by processor number ', allowableDiffDtProcNumber
          endif
      endif


      ! Set adaptive timestep
      call set_timestep(allowableAdvecDtAllProcs, allowableDiffDtAllProcs, domain % clock, dtSeconds, err_tmp)
      err = ior(err,err_tmp)



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
           write(stderrUnit,*) 'Error in calculating thickness tendency  (possibly CFL violation)'
      endif


      ! set CFL variables if they have been calculated - every block should be set to the same value!
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         if (config_print_thickness_advection_info .or. config_adaptive_timestep) then
            call mpas_pool_get_array(meshPool, 'allowableDtACFL', allowableDtACFL)
            allowableDtACFL = allowableAdvecDtAllProcs
            if (config_adaptive_timestep_include_DCFL) then
               call mpas_pool_get_array(meshPool, 'allowableDtDCFL', allowableDtDCFL)
               allowableDtDCFL = allowableDiffDtAllProcs
            endif
         endif

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in calculate_tendencies."
      endif

   !--------------------------------------------------------------------
   end subroutine calculate_tendencies



!***********************************************************************
!
!  routine update_prognostics
!
!> \brief   Updates prognostic variables 
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine updates prognostic variables using a Forward Euler time integrator.
!>  Note: The thickness field is updated here, but tracers (which have just
!>        one time level, are not updated. It is assumed that new tracer values
!>        have already been computed.
!-----------------------------------------------------------------------

   subroutine update_prognostics(domain, deltat, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: deltat

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
      type (mpas_pool_type), pointer :: meshPool, geometryPool, tendencyPool
      integer, pointer :: nCells
      logical, pointer :: config_print_thickness_advection_info

      real (kind=RKIND), dimension(:), pointer :: thickness, layerThicknessFractions
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld, layerThicknessNew, layerThicknessTendency 
      integer, dimension(:), allocatable :: masktmp  ! Temporary mask for assessing new thickness field

      err = 0

      dminfo => domain % dminfo
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'tendency', tendencyPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)

         ! State at time n
         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessOld, timeLevel=2)

         ! State at time n+1 (advanced by dt by Forward Euler)
         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessNew, timeLevel=1)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! Tendencies
         call mpas_pool_get_array(tendencyPool, 'layerThicknessTendency', layerThicknessTendency)
         
         ! Update thickness ======================
         
         layerThicknessNew = layerThicknessOld + layerThicknessTendency * deltat
         thickness = sum(layerThicknessNew, 1)

         !Optionally print some information about the new thickness
         !print *, 'thickness_tend maxval:', maxval(thickness_tend(1:mesh % nCellsSolve))       
         !print *, 'thicknessOld maxval:', maxval(thicknessOld(1:mesh % nCellsSolve))
         !print *, '  thicknessNew maxval:', maxval(thicknessNew(1:mesh % nCellsSolve))
         allocate( masktmp(nCells + 1) )
         masktmp = 0

!!!         ! if holding advance within initial extent of ice, set thickness to 0 anywhere it has expanded beyond initial extent
!!!         if (config_allow_additional_advance .eqv. .false.) then 
!!!             where ( MASK_WAS_INITIALLY_NOT_ICE(cellMaskOld) )
!!!                 thicknessNew = 0.0_RKIND
!!!             end where
!!!         endif

         ! reset negative thickness to 0.  This should not happen unless negative MB is larger than entire ice column.
         where (thickness < 0.0_RKIND)
            masktmp = 1
            thickness = 0.0_RKIND
         end where

         if (config_print_thickness_advection_info) then
            if (sum(masktmp) > 0) then
               write(stdoutUnit,*) '  Cells with negative thickness (set to 0):',sum(masktmp)
            endif

            ! Note how many cells have ice.
            masktmp = 0
            where (thickness > 0.0_RKIND)
               masktmp = 1
            end where
            write(stdoutUnit,*) '  Cells with nonzero thickness:', sum(masktmp)
         endif
         deallocate(masktmp)

         block => block % next
      end do


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in update_prognostics."
      endif


   !--------------------------------------------------------------------
   end subroutine update_prognostics


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
!>  Note: This routine replaces the old subroutines calculate_tendencies
!>        and update_prognostics.
!   TODO - Remove the old subroutines after testing
!-----------------------------------------------------------------------

   subroutine advection_solver(domain, dtSeconds, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: dtSeconds

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

      type (field2DReal), pointer :: layerThicknessField
      type (field2DReal), pointer :: layerThicknessEdgeField

      real (kind=RKIND), dimension(:), pointer :: thickness

      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessNew

      integer, pointer :: nCells

      character (len=StrKIND), pointer :: config_thickness_advection
      character (len=StrKIND), pointer :: config_tracer_advection

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

      !TODO - Replace masktmp with a scratch field?
      integer, dimension(:), allocatable :: masktmp  ! Temporary mask for assessing new thickness field

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      allowableAdvecDtAllProcs = 0.0_RKIND
      allowableDiffDtAllProcs = 0.0_RKIND

      dminfo => domain % dminfo

      ! Halo updates
      ! Note: The layer thickness and tracers must be up to date in halos before calling the advection subroutines.
      !       Currently (Nov. 2015), layerThicknessEdge (which is used by subroutine li_tendency_thickness) is updated 
      !        at the end of subroutine li_calculate_diagnostic_vars.
      !       The thermal tracers (temperature, waterfrac, enthalpy) are updated at the end of li_thermal_solver.
      !       But layerThickness (which is used by subroutine li_advection_thickness_tracers) needs an update here.
      !   
      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_field(geometryPool, 'layerThickness', layerThicknessField)
      call mpas_dmpar_exch_halo_field(layerThicknessField)
      call mpas_timer_stop("halo updates")

      ! ===
      ! === Thickness and tracer tendencies
      ! ===
      allowableAdvecDtOnProc = 1.0e36 ! set to large number
      allowableDiffDtOnProc = 1.0e36 ! set to large number

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! Note: The old and new layer thickness are input and output arguments, respectively, for the main advection subroutines
         ! TODO - Remove the time levels, and simply fetch the layer thickness from the previous timestep?
         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessOld, timeLevel=2)
         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessNew, timeLevel=1)

         if (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'fo') then

            ! Note: This subroutine requires that layerThickness and tracers are correct in halos.

            if (config_print_thickness_advection_info) then
               write(stderrUnit,*) 'Advect thickness and tracers using new scheme...'
            endif

            call li_advection_thickness_tracers(&
                 dtSeconds,              &
                 meshPool,               &
                 velocityPool,           &
                 geometryPool,           &
                 thermalPool,            &
                 scratchPool,            &
                 layerThicknessOld,      &
                 layerThicknessNew,      &
                 allowableAdvecDt,       &
                 err_tmp)

         elseif (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'none') then

            if (config_print_thickness_advection_info) then
               write(stderrUnit,*) 'Advect thickness (but not tracers) using new scheme...'
            endif

            call li_advection_thickness_tracers(&
                 dtSeconds,              &
                 meshPool,               &
                 velocityPool,           &
                 geometryPool,           &
                 thermalPool,            &
                 scratchPool,            &
                 layerThicknessOld,      &
                 layerThicknessNew,      &
                 allowableAdvecDt,       &
                 err_tmp,   advectTracersIn = .false.)

         endif

         err = ior(err,err_tmp)
         allowableAdvecDtOnProc = min(allowableAdvecDtOnProc, allowableAdvecDt)

         ! Calculate diffusive CFL timestep, if needed
         if (config_adaptive_timestep_include_DCFL) then
            call li_calculate_apparent_diffusivity(meshPool, velocityPool, scratchPool, geometryPool, allowableDiffDt)
            allowableDiffDtOnProc = min(allowableDiffDtOnProc, allowableDiffDt)
         endif

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

      ! Now that we have exited the block loop, do any needed halo updates.
      ! TODO - Update new layer thickness or tracers?


      ! Local advective CFL info
      call mpas_set_timeInterval(allowableAdvecDtOnProcInterval, dt=allowableAdvecDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableAdvecDtOnProcInterval, timeString=allowableAdvecDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)

      if (config_print_thickness_advection_info) then
         write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on advective CFL is (Days_hhh:mmm:sss):   ' // trim(allowableAdvecDtOnProcString)
      endif

      ! Local diffusive CFL info
      if (config_adaptive_timestep_include_DCFL) then
          call mpas_set_timeInterval(allowableDiffDtOnProcInterval, dt=allowableDiffDtOnProc, ierr=err_tmp)
          err = ior(err,err_tmp)
          call mpas_get_timeInterval(allowableDiffDtOnProcInterval, timeString=allowableDiffDtOnProcString, ierr=err_tmp)
          err = ior(err,err_tmp)

          if (config_print_thickness_advection_info) then
             write(stdoutUnit,*) '  Maximum allowable time step on THIS processor based on diffusive CFL is (Days_hhh:mmm:sss):   ' // trim(allowableDiffDtOnProcString)
          endif
      endif


      ! If we are printing advection debug information or adaptive timestepping,
      ! then find out what the global CFL limit is.  Don't do this otherwise because
      ! it requires 2 unnecessary MPI communications.
      if (config_print_thickness_advection_info .or. config_adaptive_timestep) then

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
          if (config_adaptive_timestep_include_DCFL) then

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

          endif

      endif  ! config_print_thickness_advection_info or config_adaptive_timestep

      ! Write messages if they are turned on
      if (config_print_thickness_advection_info) then
          write(stdoutUnit,*) '  Maximum allowable time step for all processors based on advective CFL is (Days_hhh:mmm:sss): ' // trim(allowableAdvecDtAllProcsString) //  '  Time step is limited by processor number ', allowableAdvecDtProcNumber
          if (config_adaptive_timestep_include_DCFL) then
             write(stdoutUnit,*) '  Maximum allowable time step for all processors based on diffusive CFL is (Days_hhh:mmm:sss): ' // trim(allowableDiffDtAllProcsString) //  '  Time step is limited by processor number ', allowableDiffDtProcNumber
          endif
      endif


      ! Set adaptive timestep
      call set_timestep(allowableAdvecDtAllProcs, allowableDiffDtAllProcs, domain % clock, dtSeconds, err_tmp)
      err = ior(err,err_tmp)


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
           write(stderrUnit,*) 'Error in calculating thickness tendency  (possibly CFL violation)'
      endif


      ! set CFL variables if they have been calculated - every block should be set to the same value!
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         if (config_print_thickness_advection_info .or. config_adaptive_timestep) then
            call mpas_pool_get_array(meshPool, 'allowableDtACFL', allowableDtACFL)
            allowableDtACFL = allowableAdvecDtAllProcs
            if (config_adaptive_timestep_include_DCFL) then
               call mpas_pool_get_array(meshPool, 'allowableDtDCFL', allowableDtDCFL)
               allowableDtDCFL = allowableDiffDtAllProcs
            endif
         endif

         block => block % next
      end do

      ! Reset negative thicknesses to zero if needed

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! State at time n
!!         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessOld, timeLevel=2)

         ! State at time n+1 (advanced by dt by Forward Euler)
!!         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessNew, timeLevel=1)
         
         allocate( masktmp(nCells + 1) )
         masktmp = 0

         ! reset negative thickness to 0.  This should not happen unless negative MB is larger than entire ice column.
         !TODO - Should layerThickness also be set to 0?
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
      logical, pointer :: config_adaptive_timestep
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

      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)

      if (config_adaptive_timestep) then
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
      endif

   !--------------------------------------------------------------------
   end subroutine set_timestep



end module li_time_integration_fe


