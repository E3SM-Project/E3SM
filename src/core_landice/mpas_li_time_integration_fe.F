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
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_timer
   use mpas_vector_reconstruction
   use mpas_configure
   use li_velocity, only: li_velocity_solve
   use li_tendency
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


      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! (time level 1 should not be modified.)

!!!      procVertexMaskChanged = 0

! === Implicit column physics (vertical temperature diffusion) ===========
      !call ()

! === Calculate Tendencies ========================
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
!>  This routine calculates tendencies of thickness (and eventually tracers)
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
      type (mpas_pool_type), pointer :: tendencyPool

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness_tend
      type (field2DReal), pointer :: layerThickness_tend_field
      logical, pointer :: config_print_thickness_advection_info
      logical, pointer :: config_adaptive_timestep
      integer :: allowableDtProcNumber
      real (kind=RKIND) :: allowableDt, allowableDtOnProc, allowableDtAllProcs
      type (MPAS_TimeInterval_type) :: allowableDtOnProcInterval, allowableDtAllProcsStringInterval
      character (len=StrKIND) :: allowableDtOnProcString, allowableDtAllProcsString
      integer :: err_tmp

      integer :: y, m, d, hh, mm, ss

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)

      dminfo => domain % dminfo

      ! ===
      ! === Thickness tendencies
      ! ===
      allowableDtOnProc = 1.0e36 ! set to large number
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'tendency', tendencyPool)
         call mpas_pool_get_array(tendencyPool, 'layerThickness', layerThickness_tend)

         ! Calculate thickness tendency using state at time n =========
         call li_tendency_thickness(meshPool, velocityPool, geometryPool, layerThickness_tend, dtSeconds, dminfo, allowableDt, err_tmp)
         err = ior(err,err_tmp)
         allowableDtOnProc = min(allowableDtOnProc, allowableDt)

         block => block % next
      end do

      ! Now that we have exited the block loop, do any needed halo updates.
      ! update halos on thickness tend
      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'tendency', tendencyPool)
      call mpas_pool_get_field(tendencyPool, 'layerThickness', layerThickness_tend_field)
      call mpas_dmpar_exch_halo_field(layerThickness_tend_field)
      call mpas_timer_stop("halo updates")


      ! Build a time string of the maximum allowable dt calculated
      ! (We only need this if a CFL violation occurred or config_print_thickness_advection_info is true)
      call mpas_set_timeInterval(allowableDtOnProcInterval, dt=allowableDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableDtOnProcInterval, timeString=allowableDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)
      if (dtSeconds > allowableDtOnProc) then
         write(stderrUnit,*) 'ERROR: CFL violation on this processor.  Maximum allowable time step for this processor is (Days_hhh:mmm:sss): ' // trim(allowableDtOnProcString)
                  err = ior(err,1)
      endif
      if (config_print_thickness_advection_info) then
         write(stdoutUnit,*) '  Maximum allowable time step on THIS processor is (Days_hhh:mmm:sss):   ' // trim(allowableDtOnProcString)
      endif


      ! If we are printing advection debug information,
      ! then find out what the global CFL limit is.  Don't do this otherwise because
      ! it requires 2 unnecessary MPI communications.
      if (config_print_thickness_advection_info .or. config_adaptive_timestep) then
          ! Determine CFL limit on all procs
          call mpas_dmpar_min_real(dminfo, allowableDtOnProc, allowableDtAllProcs)
          ! Determine which processor has the limiting CFL
          if (allowableDt == allowableDtAllProcs) then
            allowableDtProcNumber = dminfo % my_proc_id
          else 
            allowableDtProcNumber = -1
          endif
          call mpas_dmpar_max_int(dminfo, allowableDtProcNumber, allowableDtProcNumber)
          call mpas_set_timeInterval(allowableDtAllProcsStringInterval, dt=allowableDtOnProc, ierr=err_tmp)
          err = ior(err,err_tmp)
          call mpas_get_timeInterval(allowableDtAllProcsStringInterval, timeString=allowableDtAllProcsString, ierr=err_tmp)
          err = ior(err,err_tmp)
      endif
      if (config_print_thickness_advection_info) then
          write(stdoutUnit,*) '  Maximum allowable time step for all processors is (Days_hhh:mmm:sss): ' // trim(allowableDtAllProcsString) //  '  Time step is limited by processor number ', allowableDtProcNumber
      endif
      if (config_adaptive_timestep) then
         call set_timestep(allowableDtAllProcs, dtSeconds, err_tmp)
         err = ior(err,err_tmp)
      elseif (dtSeconds > allowableDtAllProcs) then
          write(stdErrUnit,*) 'Error: CFL violation has occurred on some processor(s)!'
          err = ior(err,1)
      endif


      if (err > 0) then
           write(stderrUnit,*) 'Error in calculating thickness tendency  (possibly CFL violation)'
      endif


      ! ===
      ! === Tracer tendencies
      ! ===
! Commenting out tracer tendencies
!!!      ! Tracer tendencies
!!!      block => domain % blocklist
!!!      do while (associated(block))
!!!         ! Mesh information
!!!         mesh => block % mesh
!!!         ! State at time n
!!!         stateOld => block % state % time_levs(1) % state
!!!         ! Tendencies
!!!         layerThickness_tend => block % tend % layerThickness % array
!!!         tracer_tendency => block % tend % tracers % array

!!!         ! Calculate tracer tendencies ==========
!!!         ! There could be a negative layer thickness with SMB turned on!
!!!         call li_tendency_tracers(mesh, stateOld, layerThickness_tend, tracer_tendency, dt, dminfo, err)

!!!         block => block % next
!!!      end do 

!!!      if (err == 1) then
!!!          call mpas_dmpar_global_abort("An error has occurred in li_tendency_tracers. Aborting...")
!!!      endif

!!!      ! update halos on tracer tend 
!!!      select case (config_tracer_advection)
!!!      case ('None')  !===================================================
!!!          ! Do nothing - no need to waste time doing a halo update if not advecting tracers!  The tendency will be 0 everywhere
!!!      case default
!!!          call mpas_dmpar_exch_halo_field3d_real(dminfo, tracer_tendency, &
!!!                                         size(tracer_tendency,dim=1), mesh % nVertLevels, mesh % nCells, &
!!!                                         domain % blocklist % parinfo % cellsToSend, domain % blocklist % parinfo % cellsToRecv)
!!!      end select


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
!>  This routine updates prognostic variables using a 
!>  Forward Euler time integrator.
!
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

      real (kind=RKIND), dimension(:), pointer :: thicknessNew, layerThicknessFractions
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld, layerThicknessNew, layerThickness_tend
!!!      real (kind=RKIND), dimension(:,:,:), pointer :: tracer_tendency, tracersNew, tracersOld
!!!      integer, dimension(:), pointer :: cellMaskOld
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
!!!         tracersOld => stateOld % tracers % array
!!!         cellMaskOld => stateOld % cellMask % array

         ! State at time n+1 (advanced by dt by Forward Euler)
         call mpas_pool_get_array(geometryPool, 'layerThickness', layerThicknessNew, timeLevel=1)
         call mpas_pool_get_array(geometryPool, 'thickness', thicknessNew, timeLevel=1)
!!!         tracersNew => stateNew % tracers % array

         ! Tendencies
         call mpas_pool_get_array(tendencyPool, 'layerThickness', layerThickness_tend)
!!!         tracer_tendency => block % tend % tracers % array


         ! Update thickness ======================
         
         ! Commented out usage for advecting thickness as a column
         !!!thicknessNew = thicknessOld + thickness_tend * dt / SecondsInYear   
         ! Commented out usage for using FCT for thickness 
         !!!stateNew % sup_thickness % array(1,:,:) = (stateOld % tracers % array(stateOld % index_temperature,:,:) * layerThicknessOld  + tracer_tendency(stateOld % index_temperature, :, :) * dt / SecondsInYear) / (layerThicknessNew+1.0e-12)


         layerThicknessNew = layerThicknessOld + layerThickness_tend * deltat
         thicknessNew = sum(layerThicknessNew, 1)


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
         where (thicknessNew < 0.0_RKIND)
            masktmp = 1
            thicknessNew = 0.0_RKIND
!!!            stateNew % iceArea % array = 0.0_RKIND
         end where

         if (config_print_thickness_advection_info) then
            if (sum(masktmp) > 0) then
               write(stdoutUnit,*) '  Cells with negative thickness (set to 0):',sum(masktmp)
            endif

            ! Note how many cells have ice.
            masktmp = 0
            where (thicknessNew > 0.0_RKIND)
               masktmp = 1
            end where
            write(stdoutUnit,*) '  Cells with nonzero thickness:', sum(masktmp)
         endif
         deallocate(masktmp)



!!!         ! Calculate new tracer values =================
!!!         if (config_tracer_advection .ne. 'None') then
!!!           do iTracer = 1, size(tracersNew, 1)
!!!               where (layerThicknessNew > 0.0_RKIND)
!!!                   tracersNew(iTracer,:,:) = (tracersOld(iTracer,:,:) * layerThicknessOld &
!!!                       + tracer_tendency(iTracer,:,:) * dt / SecondsInYear)  /  (layerThicknessNew)
!!!               elsewhere
!!!                   ! May or may not want to assign tracer values to non-ice cells
!!!                   tracersNew(iTracer,:,:) = 0.0_RKIND
!!!               end where
!!!           end do
!!!         endif

!!!         ! Apply calving after we have updated the new state  - TODO Is this the right place?
!!!         call li_apply_calving(mesh, stateNew, err)


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
!  routine set_timestep
!
!> \brief   Adjusts the time step based on the CFL condition.
!> \author  Matthew Hoffman
!> \date    23 Jan 2014
!> \details
!>  This routine sdjusts the time step based on the CFL condition.
!
!-----------------------------------------------------------------------
   subroutine set_timestep(allowableDt, dtSeconds, err)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: allowableDt

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: dtSeconds
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_adaptive_timestep
      real (kind=RKIND), pointer :: config_adaptive_timestep_CFL_fraction
      real (kind=RKIND), pointer :: config_max_adaptive_timestep
      real (kind=RKIND), pointer :: config_min_adaptive_timestep

      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_CFL_fraction', config_adaptive_timestep_CFL_fraction)
      call mpas_pool_get_config(liConfigs, 'config_max_adaptive_timestep', config_max_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_min_adaptive_timestep', config_min_adaptive_timestep)

      if (config_adaptive_timestep) then
         dtSeconds = min(allowableDt * config_adaptive_timestep_CFL_fraction, config_max_adaptive_timestep)
         write(stdOutUnit,*) '  Setting time step (days) to:', dtSeconds / (86400.0)
         if (dtSeconds < config_min_adaptive_timestep) then
            write(stdErrUnit,*) 'ERROR: New deltat is less than config_min_adaptive_timestep.'
            err = 1
         endif
      endif

   !--------------------------------------------------------------------
   end subroutine set_timestep



end module li_time_integration_fe


