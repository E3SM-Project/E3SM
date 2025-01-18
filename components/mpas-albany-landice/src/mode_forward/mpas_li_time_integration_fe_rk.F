! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_integration_fe_rk
!
!> \brief MPAS land ice Forward Euler and Runge-Kutta time integration schemes
!> \author Matt Hoffman, Trevor Hillebrand
!> \date   17 April 2011, Runge-Kutta added Sept 2023
!> \details
!>  This module contains the Forward Euler and Runge-Kutta time integration schemes
!
!-----------------------------------------------------------------------

module li_time_integration_fe_rk

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timer
   use mpas_vector_reconstruction
   use mpas_log

   use li_advection
   use li_calving, only: li_calve_ice, li_restore_calving_front, &
                         li_remove_icebergs, li_calculate_damage, &
                         li_finalize_damage_after_advection
   use li_thermal, only: li_thermal_solver, li_enthalpy_to_temperature_kelvin
   use li_iceshelf_melt
   use li_diagnostic_vars
   use li_setup
   use li_constants
   use li_mesh
   use li_mask
   use li_tracer_advection_fct_shared

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

   public :: li_time_integrator_forwardeuler_rungekutta

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
!  routine li_time_integrator_forwardeuler_rungekutta
!
!> \brief   Forward Euler and Runge-Kutta time integration schemes
!> \author  Matthew Hoffman, Trevor Hillebrand
!> \date    10 January 2012, updated Sept 2023 for Runge-Kutta
!> \details
!>  This routine performs Forward Euler and Runge-Kutta time integration.
!
!-----------------------------------------------------------------------
   subroutine li_time_integrator_forwardeuler_rungekutta(domain, err)

      use li_subglacial_hydro
      use li_velocity
      use li_bedtopo
      use li_mask
      use li_advection, only: li_grounded_to_floating
      use li_ocean_extrap

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
      type (mpas_pool_type), pointer :: geometryPool, thermalPool, meshPool, velocityPool
      
      logical, pointer :: config_restore_calving_front, &
                          config_restore_calving_front_prevent_retreat
      logical, pointer :: config_calculate_damage
      logical, pointer :: config_finalize_damage_after_advection
      logical, pointer :: config_update_velocity_before_calving
      character (len=StrKIND), pointer :: config_thickness_advection, &
                                          config_tracer_advection
      character (len=StrKIND), pointer :: config_thermal_solver
      character (len=StrKIND), pointer :: config_time_integration
      integer, pointer :: config_rk_order, config_rk3_stages
      integer :: rkStage, iCell, iTracer, k
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, &
                                                    temperature, &
                                                    waterFrac
      real (kind=RKIND), dimension(:), pointer :: thickness, damage, passiveTracer2d
      real (kind=RKIND), dimension(:), pointer :: sfcMassBalApplied, &
                                                  basalMassBalApplied, &
                                                  groundedSfcMassBalApplied, &
                                                  groundedBasalMassBalApplied, &
                                                  floatingBasalMassBalApplied, &
                                                  fluxAcrossGroundingLine, &
                                                  fluxAcrossGroundingLineOnCells, &
                                                  groundedToFloatingThickness
      real (kind=RKIND), dimension(:), allocatable :: sfcMassBalAccum, &
                                                      basalMassBalAccum, &
                                                      groundedSfcMassBalAccum, &
                                                      groundedBasalMassBalAccum, &
                                                      floatingBasalMassBalAccum, &
                                                      fluxAcrossGroundingLineAccum, &
                                                      fluxAcrossGroundingLineOnCellsAccum
      real (kind=RKIND), dimension(:), pointer :: layerCenterSigma, layerThicknessFractions
      integer, dimension(:), pointer :: cellMaskPrev  ! cell mask before advection
      real (kind=RKIND), dimension(:,:), allocatable :: layerThicknessPrev, &
                                                        layerThicknessTmp, &
                                                        temperaturePrev, &
                                                        waterFracPrev, &
                                                        passiveTracer3d, &
                                                        passiveTracer3dPrev, &
                                                        damage3d, &
                                                        damage3dPrev
      integer, pointer :: nVertLevels
      integer, pointer :: nCells, nEdges
      integer :: iCell1, iCell2, iEdge, theGroundedCell
      integer, dimension(:), pointer :: edgeMask, cellMask
      real (kind=RKIND), pointer :: deltat, config_ice_density
      real (kind=RKIND) :: deltatFull
      real (kind=RKIND), dimension(4) :: rkSubstepWeights
      real (kind=RKIND), dimension(4) :: rkSSPweights
      real (kind=RKIND), dimension(4) :: rkTendWeights ! Weights used for calculating budget terms
      integer :: nRKstages
      logical :: solveVeloAfterCalving

      err = 0
      err_tmp = 0
      solveVeloAfterCalving = .false.

      call mpas_pool_get_config(liConfigs, 'config_restore_calving_front', config_restore_calving_front)
      call mpas_pool_get_config(liConfigs, 'config_restore_calving_front_prevent_retreat', config_restore_calving_front_prevent_retreat)
      call mpas_pool_get_config(liConfigs, 'config_calculate_damage', config_calculate_damage)
      call mpas_pool_get_config(liConfigs, 'config_finalize_damage_after_advection', config_finalize_damage_after_advection)
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_rk_order', config_rk_order)
      call mpas_pool_get_config(liConfigs, 'config_rk3_stages', config_rk3_stages)
      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_update_velocity_before_calving', config_update_velocity_before_calving)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'thermal', thermalPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)

      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(velocityPool, 'fluxAcrossGroundingLineOnCells', fluxAcrossGroundingLineOnCells)
      call mpas_pool_get_array(geometryPool, 'groundedToFloatingThickness', groundedToFloatingThickness)

      allocate(temperaturePrev(nVertLevels, nCells+1))
      allocate(waterFracPrev(nVertLevels, nCells+1))
      allocate(layerThicknessPrev(nVertLevels, nCells+1))
      allocate(layerThicknessTmp(nVertLevels, nCells+1))
      allocate(damage3dPrev(nVertLevels, nCells+1))
      allocate(damage3d(nVertLevels, nCells+1))
      allocate(passiveTracer3dPrev(nVertLevels, nCells+1))
      allocate(passiveTracer3d(nVertLevels, nCells+1))
      allocate(cellMaskPrev(nCells+1))

      allocate(sfcMassBalAccum(nCells+1))
      allocate(basalMassBalAccum(nCells+1))
      allocate(groundedSfcMassBalAccum(nCells+1))
      allocate(groundedBasalMassBalAccum(nCells+1))
      allocate(floatingBasalMassBalAccum(nCells+1))
      allocate(fluxAcrossGroundingLineAccum(nEdges+1))
      allocate(fluxAcrossGroundingLineOnCellsAccum(nCells+1))

      temperaturePrev(:,:) = 0.0_RKIND
      waterFracPrev(:,:) = 0.0_RKIND
      layerThicknessPrev(:,:) = 0.0_RKIND
      layerThicknessTmp(:,:) = 0.0_RKIND
      damage3dPrev(:,:) = 0.0_RKIND
      damage3d(:,:) = 0.0_RKIND
      passiveTracer3dPrev(:,:) = 0.0_RKIND
      cellMaskPrev(:) = 0
      passiveTracer3d(:,:) = 0.0_RKIND

      sfcMassBalAccum(:) = 0.0_RKIND
      basalMassBalAccum(:) = 0.0_RKIND
      groundedSfcMassBalAccum(:) = 0.0_RKIND
      groundedBasalMassBalAccum(:) = 0.0_RKIND
      floatingBasalMassBalAccum(:) = 0.0_RKIND
      fluxAcrossGroundingLineAccum(:) = 0.0_RKIND
      fluxAcrossGroundingLineOnCellsAccum(:) = 0.0_RKIND

! === Prepare for advection (including CFL checks) ===========
! This has to come first currently, because it sets the time step!
      call mpas_timer_start("advection prep")
      call prepare_advection(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("advection prep")

! === Advance the clock before all other physics happen ===========
      call mpas_timer_start("advancing clock")
      call advance_clock(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("advancing clock")

!TODO: Determine whether grounded melting should in fact be called first
! === Ocean forcing extrapolation into ice-shelf cavities ===========
      call mpas_timer_start("ocean forcing extrapolation")
      call li_ocean_extrap_solve(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("ocean forcing extrapolation")

! === Face melting for grounded ice ===========
      call mpas_timer_start("face melting for grounded ice")
      call li_face_melt_grounded_ice(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("face melting for grounded ice")

! === Basal melting for floating ice ===========
      call mpas_timer_start("basal melting for floating ice")
      call li_basal_melt_floating_ice(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("basal melting for floating ice")

! === Implicit column physics (vertical temperature diffusion) ===========
      call mpas_timer_start("vertical therm")
      call li_thermal_solver(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("vertical therm")

! *** TODO: Should basal melt rate calculation, column physics, and hydrology go inside RK loop? ***
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
      call mpas_pool_get_array(geometryPool, 'passiveTracer2d', passiveTracer2d)
      call mpas_pool_get_array(geometryPool, 'damage', damage)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
      ! Save relevant fields before RK loop, to be used in update at the end
      temperaturePrev(:,:) = temperature(:,:)
      waterFracPrev(:,:) = waterFrac(:,:)
      layerThicknessPrev(:,:) = layerThickness(:,:)
      cellMaskPrev(:) = cellMask(:)
      do k = 1, nVertLevels
         damage3dPrev(k,:) = damage(:)
         passiveTracer3dPrev(k,:) = passiveTracer2d(:)
      enddo
      deltatFull = deltat ! Save deltat in order to reset it at end of RK loop

      if ( (trim(config_tracer_advection) == 'fct') .or. &
           (trim(config_thickness_advection) == 'fct') ) then
         call li_tracer_advection_fct_shared_init(geometryPool, err_tmp)
         if (err_tmp /= 0) then
            err = 1
            call mpas_log_write(                                 &
               'Error encountered during fct tracer advection shared init', &
               MPAS_LOG_ERR)
         endif
      endif

      ! Set RK weights based on desired time integration method. Note
      ! that rkSubstepWeights are used to update at each sub-step, and
      ! are thus offset from the typical writing of the coefficients
      ! by one index. e.g., the coefficients for RK4 are usually written
      ! (0, 1/2, 1/2, 1), while we use (1/2, 1/2, 1). The last entry of 1.0
      ! is simply for ease of implementation.
      rkSSPweights(:) = -9999.0_RKIND ! Initialized to this value to make it obvious if
                                      ! a weight is used that should not be. Appropriate
                                      ! weights are updated for each case below
      rkTendWeights(:) = -9999.0_RKIND
      rkSubstepWeights(:) = -9999.0_RKIND
      if (trim(config_time_integration) == 'forward_euler') then
         nRKstages = 1
         rkSSPweights(1) = 1.0_RKIND
         rkTendWeights(1) = 1.0_RKIND
         rkSubstepWeights(1) = 1.0_RKIND
      elseif ( (trim(config_time_integration) == 'runge_kutta') .and. &
               (config_rk_order == 2) ) then
         ! use Strong Stability Preserving RK2. Could also
         ! add config option to use standard endpoint or midpoint methods, but
         ! these are algorithmically more complex and less suitable for our domains 
         nRKstages = 2
         rkSubstepWeights(1:2) = 1.0_RKIND

         rkSSPweights(1) = 1.0_RKIND
         rkSSPweights(2) = 0.5_RKIND

         rkTendWeights(1) = 0.5_RKIND
         rkTendWeights(2) = 0.5_RKIND
      elseif ( (trim(config_time_integration) == 'runge_kutta') .and. &
               (config_rk_order == 3) ) then
         if (config_rk3_stages == 3) then
            ! use three-stage Strong Stability Preserving RK3
            nRKstages = 3
            rkSubstepWeights(1:3) = 1.0_RKIND

            rkSSPweights(1) = 1.0_RKIND
            rkSSPweights(2) = 3.0_RKIND / 4.0_RKIND
            rkSSPweights(3) = 1.0_RKIND / 3.0_RKIND

            rkTendWeights(1) = 1.0_RKIND / 6.0_RKIND
            rkTendWeights(2) = 1.0_RKIND / 6.0_RKIND
            rkTendWeights(3) = 2.0_RKIND / 3.0_RKIND
        elseif (config_rk3_stages == 4) then
            ! use four-stage Strong Stability Preserving RK3
            nRKstages = 4
            rkSubstepWeights(:) = 0.5_RKIND

            rkSSPweights(1) = 1.0_RKIND
            rkSSPweights(2) = 0.0_RKIND
            rkSSPweights(3) = 2.0_RKIND / 3.0_RKIND
            rkSSPweights(4) = 0.0_RKIND

            rkTendWeights(1) = 1.0_RKIND / 6.0_RKIND
            rkTendWeights(2) = 1.0_RKIND / 6.0_RKIND
            rkTendWeights(3) = 1.0_RKIND / 6.0_RKIND
            rkTendWeights(4) = 1.0_RKIND / 2.0_RKIND
        else
            err = 1
            call mpas_log_write('config_rk3_stages must 3 or 4', &
                                messageType=MPAS_LOG_ERR)
            return
        endif
      else
         err = 1
         call mpas_log_write('config_time_integration = ' // trim(config_time_integration) &
                             // ' is not supported with config_rk_order = $i', &
                             intArgs=(/config_rk_order/), messageType=MPAS_LOG_ERR)
         return
      endif

      ! Calculate masks prior to RK loop, but do not update masks within the loop
      ! to preserve the accuracy of time integration.
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)

! *** Start RK loop ***
      do rkStage = 1, nRKstages
         call mpas_log_write('beginning rk stage $i of $i', &
                             intArgs=(/rkStage, nRKstages/))
         deltat = deltatFull * rkSubstepWeights(rkStage)

         ! === calculate damage ===========
         if (config_calculate_damage) then
             call mpas_timer_start("damage")
             call li_calculate_damage(domain, err_tmp)
             err = ior(err, err_tmp)
             call mpas_timer_stop("damage")
         endif

         ! === Compute new state for prognostic variables ===
         call mpas_timer_start("advect thickness and tracers")
         call advection_solver(domain, err_tmp)
         err = ior(err, err_tmp)
         call mpas_timer_stop("advect thickness and tracers")

         ! If using SSP RK, then update thickness and tracers incrementally.
         ! For first RK stage, thickness and tracer updates above are sufficient.
         ! Likewise, for the 4-stage SSP RK3 the last stage is just a forward euler update.
         if ( (rkStage > 1) .and. (rkStage < 4) )  then
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
            call mpas_pool_get_array(geometryPool, 'passiveTracer2d', passiveTracer2d)
            call mpas_pool_get_array(geometryPool, 'damage', damage)

            call mpas_pool_get_array(thermalPool, 'temperature', temperature)
            call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)

            layerThicknessTmp(:,:) = layerThickness(:,:)
            layerThickness(:,:) = rkSSPweights(rkStage) * layerThicknessPrev(:,:) + &
                                     (1.0_RKIND - rkSSPweights(rkStage)) * layerThickness(:,:)
            thickness = sum(layerThickness, 1)
            ! Do not calculate masks after updating thickness! We need to keep masks
            ! constant for now to preserve accuracy of time integration

            if (trim(config_thermal_solver) .ne. 'none') then
               do iCell = 1, nCells
                  do k = 1, nVertLevels
                     if (layerThickness(k,iCell) > 0.0_RKIND) then
                        temperature(k,iCell) = ( rkSSPweights(rkStage) * temperaturePrev(k,iCell) * layerThicknessPrev(k,iCell) + &
                                             (1.0_RKIND - rkSSPweights(rkStage)) * temperature(k,iCell) * &
                                             layerThicknessTmp(k,iCell) ) / layerThickness(k,iCell)
                        waterFrac(k,iCell) = ( rkSSPweights(rkStage) * waterFracPrev(k,iCell) * layerThicknessPrev(k,iCell) + &
                                             (1.0_RKIND - rkSSPweights(rkStage)) * waterFrac(k,iCell) * &
                                             layerThicknessTmp(k,iCell) ) / layerThickness(k,iCell)
                     endif
                  enddo
               enddo
            endif

            if (config_calculate_damage) then
               do iCell = 1, nCells
                  do k = 1, nVertLevels
                     if (layerThickness(k,iCell) > 0.0_RKIND) then
                        damage3d(k,iCell) = ( rkSSPweights(rkStage) * damage3dPrev(k,iCell) * layerThicknessPrev(k,iCell) + &
                                              (1.0_RKIND - rkSSPweights(rkStage)) * damage(iCell) * &
                                              layerThicknessTmp(k,iCell) ) / layerThickness(k,iCell)
                     else
                        damage3d(k,iCell) = 0.0_RKIND
                     endif
                  enddo
                  damage(iCell) = sum(damage3d(:, iCell) * layerThicknessFractions)
               enddo
            endif

            do iCell = 1, nCells
               do k = 1, nVertLevels
                  if (layerThickness(k,iCell) > 0.0_RKIND) then
                     passiveTracer3d(k,iCell) = ( rkSSPweights(rkStage) * passiveTracer3dPrev(k,iCell) * layerThicknessPrev(k,iCell) + &
                                                  (1.0_RKIND - rkSSPweights(rkStage)) * passiveTracer2d(iCell) * &
                                                  layerThicknessTmp(k,iCell) ) / layerThickness(k,iCell)
                  else
                     passiveTracer3d(k,iCell) = 0.0_RKIND
                  endif
               enddo
               passiveTracer2d(iCell) = sum(passiveTracer3d(:,iCell) * layerThicknessFractions)
            enddo

         endif
       
         ! === Ensure damage is within bounds before velocity solve ===
         if ( config_finalize_damage_after_advection ) then
             call mpas_timer_start("finalize damage")
             call li_finalize_damage_after_advection(domain, err_tmp)
             err = ior(err, err_tmp)
             call mpas_timer_stop("finalize damage")
         endif

         call mpas_pool_get_array(geometryPool, 'sfcMassBalApplied', sfcMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'groundedSfcMassBalApplied', groundedSfcMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'basalMassBalApplied', basalMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBalApplied', floatingBasalMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'groundedBasalMassBalApplied', groundedBasalMassBalApplied)
         call mpas_pooL_get_array(velocityPool, 'fluxAcrossGroundingLine', fluxAcrossGroundingLine)
         call mpas_pooL_get_array(velocityPool, 'fluxAcrossGroundingLineOnCells', fluxAcrossGroundingLineOnCells)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! update budgets
         sfcMassBalAccum = sfcMassBalAccum + rkTendWeights(rkStage) * sfcMassBalApplied
         groundedSfcMassBalAccum = groundedSfcMassBalAccum + rkTendWeights(rkStage) * groundedSfcMassBalApplied
         basalMassBalAccum = basalMassBalAccum + rkTendWeights(rkStage) * basalMassBalApplied
         groundedBasalMassBalAccum = groundedBasalMassBalAccum + rkTendWeights(rkStage) * groundedBasalMassBalApplied
         floatingBasalMassBalAccum = floatingBasalMassBalAccum + rkTendWeights(rkStage) * floatingBasalMassBalApplied
         fluxAcrossGroundingLineAccum = fluxAcrossGroundingLineAccum + rkTendWeights(rkStage) * fluxAcrossGroundingLine
         fluxAcrossGroundingLineOnCellsAccum = fluxAcrossGroundingLineOnCellsAccum + rkTendWeights(rkStage) * fluxAcrossGroundingLineOnCells

         ! Halo updates
         call mpas_timer_start("halo updates")

         call mpas_dmpar_field_halo_exch(domain, 'thickness')
         call mpas_dmpar_field_halo_exch(domain, 'temperature')
         call mpas_dmpar_field_halo_exch(domain, 'waterFrac')
         call mpas_dmpar_field_halo_exch(domain, 'damage')
         call mpas_dmpar_field_halo_exch(domain, 'passiveTracer2d')

         call mpas_timer_stop("halo updates")

         ! Update velocity for each RK step
         ! === Solve Velocity =====================
         if ( config_update_velocity_before_calving .or. ( (.not. config_update_velocity_before_calving) &
              .and. (rkStage < nRKstages) ) ) then

            if (config_restore_calving_front) then
               ! restore the calving front to its initial position before velocity solve.
               call li_restore_calving_front(domain, err_tmp)
               err = ior(err, err_tmp)
            endif
            ! We need to remove icebergs between RK stages because the
            ! main calving routine is not called until after the RK loop.
            ! This frequently results in icebergs that causes intermediate
            ! RK stage velocity solves to fail.
            call li_remove_icebergs(domain)

            call li_velocity_solve(domain, solveVelo=.true., err=err_tmp)
            err = ior(err, err_tmp)
         endif
! *** end RK loop ***
      enddo

! Finalize budget updates
      ! Update masks after RK integration
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'cellMask')
      call mpas_dmpar_field_halo_exch(domain, 'edgeMask')
      call mpas_dmpar_field_halo_exch(domain, 'vertexMask')
      call mpas_timer_stop("halo updates")

      ! Calculate volume converted from grounded to floating
      ! This needs to be determined after SMB/BMB are applied because those can change floating/grounded state
      call li_grounded_to_floating(cellMaskPrev, cellMask, thickness, groundedToFloatingThickness, nCells)

      sfcMassBalApplied(:) = sfcMassBalAccum(:)
      groundedSfcMassBalApplied(:) = groundedSfcMassBalAccum(:)
      basalMassBalApplied(:) = basalMassBalAccum(:)
      groundedBasalMassBalApplied(:) = groundedBasalMassBalAccum(:)
      floatingBasalMassBalApplied(:) = floatingBasalMassBalAccum(:)
      fluxAcrossGroundingLine(:) = fluxAcrossGroundingLineAccum(:)
      fluxAcrossGroundingLineOnCells(:) = fluxAcrossGroundingLineOnCellsAccum(:)

      ! Deallocate arrays for fct
      if ( (trim(config_thickness_advection) .eq. 'fct') .or. &
           (trim(config_tracer_advection) .eq. 'fct') ) then
         deallocate( nAdvCellsForEdge, &
                     advCellsForEdge, &
                     advCoefs, &
                     advCoefs3rd, &
                     advMaskHighOrder, &
                     advMask2ndOrder)
      endif

! Reset time step to full length after RK loop
      deltat = deltatFull

! === Update subglacial hydrology  ===========
! It's not clear where the best place to call this should be.
! Seems sensible to put it after thermal evolution is complete to get updated basal melting source term.
! Also seems (might be?) better to put it after geometry evolution.
! We want it before the velocity solve since the hydro model can control  the velo basal b.c.
      call mpas_timer_start("subglacial hydro")
      call li_SGH_solve(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("subglacial hydro")

! === Calve ice ========================
      call mpas_timer_start("calve_ice")

      ! ice calving
      if ( config_update_velocity_before_calving ) then
         call li_calve_ice(domain, err_tmp, solveVeloAfterCalving)
      else
         call li_calve_ice(domain, err_tmp)
         solveVeloAfterCalving = .true.
      endif

      err = ior(err, err_tmp)
      if (config_restore_calving_front) then
         ! restore the calving front to its initial position before velocity solve.
         call li_restore_calving_front(domain, err_tmp)
         err = ior(err, err_tmp)
      endif
      call mpas_timer_stop("calve_ice")

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'cellMask')
      call mpas_dmpar_field_halo_exch(domain, 'edgeMask')
      call mpas_dmpar_field_halo_exch(domain, 'vertexMask')
      call mpas_timer_stop("halo updates")

! === Update bed topo =====================
! It's not clear when the best time to do this is.
! Seems cleaner to do it either before or after all of the time evolution of the ice
! is complete.  Putting it after.
      call li_bedtopo_solve(domain, err=err_tmp)
      err = ior(err, err_tmp)

! === Solve velocity for final state =====================
      if (solveVeloAfterCalving) then
         call li_velocity_solve(domain, solveVelo=.true., err=err_tmp)
         err = ior(err, err_tmp)
      endif

! === Calculate diagnostic variables for new state =====================

      call li_calculate_diagnostic_vars(domain, err=err_tmp)
      err = ior(err, err_tmp)


! === Cleanup & Misc. =============================

      ! === error check
      if (err == 1) then
          call mpas_log_write("An error has occurred in li_time_integrator_forwardeuler_rungekutta.", MPAS_LOG_ERR)
      endif

      deallocate(temperaturePrev)
      deallocate(waterFracPrev)
      deallocate(layerThicknessPrev)
      deallocate(layerThicknessTmp)
      deallocate(damage3dPrev)
      deallocate(damage3d)
      deallocate(passiveTracer3dPrev)
      deallocate(passiveTracer3d)
      deallocate(cellMaskPrev)
      deallocate(sfcMassBalAccum)
      deallocate(basalMassBalAccum)
      deallocate(groundedSfcMassBalAccum)
      deallocate(groundedBasalMassBalAccum)
      deallocate(floatingBasalMassBalAccum)
      deallocate(fluxAcrossGroundingLineAccum)
   !--------------------------------------------------------------------
   end subroutine li_time_integrator_forwardeuler_rungekutta



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
      real (kind=RKIND), pointer :: calvingCFLdt, faceMeltingCFLdt
      integer, pointer :: processLimitingTimestep, config_rk_order, config_rk3_stages
      integer, dimension(:), pointer :: edgeMask

      logical, pointer :: config_print_thickness_advection_info
      logical, pointer :: config_adaptive_timestep
      logical, pointer :: config_adaptive_timestep_include_DCFL

      character (len=StrKIND), pointer :: &
           config_thickness_advection, &  ! method for advecting thickness and tracers
           config_time_integration

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
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      call mpas_pool_get_config(liConfigs, 'config_rk_order', config_rk_order)
      call mpas_pool_get_config(liConfigs, 'config_rk3_stages', config_rk3_stages)

      if (trim(config_thickness_advection) == 'none') then
         if (config_adaptive_timestep) then
            call mpas_log_write("Adaptive time stepper cannot be used when config_thickness_advection is set to 'none'", &
                    MPAS_LOG_ERR)
            err = ior(err,1)
         endif
         return ! skip this routine
      endif

      allowableAdvecDtAllProcs = 0.0_RKIND
      allowableDiffDtAllProcs = 0.0_RKIND

      dminfo => domain % dminfo

      ! Initialize

      allowableAdvecDtOnProc = 1.0e36_RKIND ! set to large number
      allowableDiffDtOnProc = 1.0e36_RKIND ! set to large number

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)

         ! compute normal velocities and advective CFL limit for this block

         call li_layer_normal_velocity( &
              meshPool,               &
              normalVelocity,         &
              edgeMask,               &
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

      ! If using 4-stage SSPRK3, CFL number of 2 is theoretically allowed
      if ( (trim(config_time_integration) == 'runge_kutta') .and. &
         (config_rk_order == 3) .and. (config_rk3_stages == 4) ) then
          allowableAdvecDtOnProc = allowableAdvecDtOnProc * 2.0_RKIND
          allowableDiffDtOnProc = allowableDiffDtOnProc * 2.0_RKIND
      endif

      ! Local advective CFL info
      call mpas_set_timeInterval(allowableAdvecDtOnProcInterval, dt=allowableAdvecDtOnProc, ierr=err_tmp)
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(allowableAdvecDtOnProcInterval, timeString=allowableAdvecDtOnProcString, ierr=err_tmp)
      err = ior(err,err_tmp)

      if (config_print_thickness_advection_info) then
         call mpas_log_write('  Maximum allowable time step on THIS processor based on advective CFL is (Days_hhh:mmm:sss):   ' &
            // trim(allowableAdvecDtOnProcString))
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
         call mpas_log_write('  Maximum allowable time step on THIS processor based on diffusive CFL is (Days_hhh:mmm:sss):   ' &
            // trim(allowableDiffDtOnProcString))
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
          call mpas_log_write('  Maximum allowable time step for all processors based on advective CFL is (Days_hhh:mmm:sss): ' &
             // trim(allowableAdvecDtAllProcsString) // '  Time step is limited by processor number $i', &
             intArgs=(/allowableAdvecDtProcNumber/))
          if (config_adaptive_timestep_include_DCFL) then
             call mpas_log_write('  Maximum allowable time step for all processors based on diffusive CFL is (Days_hhh:mmm:sss): '&
                // trim(allowableDiffDtAllProcsString) // '  Time step is limited by processor number $i', &
                intArgs=(/allowableDiffDtProcNumber/))
          endif
      endif


      ! Set adaptive timestep if needed
      if (config_adaptive_timestep) then
         call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
         call mpas_pool_get_array(geometryPool, 'faceMeltingCFLdt', faceMeltingCFLdt)
         call mpas_pool_get_array(meshPool, 'processLimitingTimestep', processLimitingTimestep)
         call set_timestep(allowableAdvecDtAllProcs, allowableDiffDtAllProcs, calvingCFLdt, faceMeltingCFLdt, domain % clock, dtSeconds, processLimitingTimestep, err_tmp)
         err = ior(err,err_tmp)
         ! Set new value on all blocks
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_array(meshPool, 'deltat', deltat)

            deltat = dtSeconds

            block => block % next
         end do
      else
         ! If not adaptive, Get dt from any block to check for CFL violation below
         block => domain % blocklist
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
      end if


      ! Check for CFL error before finishing
      if (deltat > allowableAdvecDtOnProc) then
         call mpas_log_write('Advective CFL violation on this processor.  ' // &
            'Maximum allowable time step for this processor is (Days_hhh:mmm:sss): ' // trim(allowableAdvecDtOnProcString), &
             MPAS_LOG_ERR)
         err = ior(err,1)
      endif

      ! Local diffusive CFL info
      if ( (config_adaptive_timestep_include_DCFL) .and. (deltat > allowableDiffDtOnProc) ) then
         call mpas_log_write('Diffusive CFL violation on this processor.  ' // &
            'Maximum allowable time step for this processor is (Days_hhh:mmm:sss): ' // trim(allowableDiffDtOnProcString), &
            MPAS_LOG_WARN)
      endif

      if (err > 0) then
           call mpas_log_write('Error in calculating thickness advection  (possibly CFL violation)', MPAS_LOG_ERR)
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
          call mpas_log_write("An error has occurred in prepare_advection.", MPAS_LOG_ERR)
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
      use li_mask

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
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: layerNormalVelocity

      integer, dimension(:), pointer :: edgeMask
      real (kind=RKIND) :: allowableDtACFL

      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), dimension(:,:), pointer :: waterFrac
      real (kind=RKIND), dimension(:,:), pointer :: enthalpy

      integer, pointer :: nCells

      character (len=StrKIND), pointer :: config_thickness_advection
      logical, pointer :: config_restore_thickness_after_advection
      character (len=StrKIND), pointer :: config_tracer_advection
      character (len=StrKIND), pointer :: config_time_integration

      logical, pointer :: config_print_thickness_advection_info

      !TODO - Replace masktmp with a scratch field?
      integer, dimension(:), allocatable :: masktmp  ! Temporary mask for assessing new thickness field

      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_restore_thickness_after_advection', config_restore_thickness_after_advection)

      dminfo => domain % dminfo

      ! Halo updates
      ! Note: The layer thickness and tracers must be up to date in halos before calling the advection subroutines.
      !   The thermal tracers (temperature, waterFrac, enthalpy) are updated at the end of li_thermal_solver.
      !   But thickness (which is used by subroutine li_advection_thickness_tracers) needs an update here. TODO: confirm this

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_timer_stop("halo updates")


      ! ===
      ! === Calculate layerThicknessEdge, which is needed for advection
      ! ===
      if (trim(config_thickness_advection) == 'fo' .or. trim(config_tracer_advection) == 'fo' .or. &
          trim(config_thickness_advection) == 'fct' .or. trim(config_tracer_advection) == 'fct') then
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            if (trim(config_time_integration) == "runge_kutta") then
               call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
               call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
               call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)

               call li_layer_normal_velocity( &
                    meshPool,               &
                    normalVelocity,         &
                    edgeMask,               &
                    layerNormalVelocity,    &
                    allowableDtACFL,        &
                    err_tmp)
               err = ior(err,err_tmp)
            endif
            call calculate_layerThicknessEdge(meshPool, geometryPool, velocityPool, err_tmp)
            err = ior(err,err_tmp)

            block => block % next
         end do

         ! Halo update layerThicknessEdge - outer halo may be wrong due to requiring velocity
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'layerThicknessEdge')
         call mpas_timer_stop("halo updates")
      endif ! fo advection option



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

         if ( (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'fo') .or. &
              (trim(config_thickness_advection) == 'fct') .or. &
              (trim(config_thickness_advection) == 'fo' .and. trim(config_tracer_advection) == 'fct') ) then

            ! Note: This subroutine requires that thickness and tracers are correct in halos

            if (config_print_thickness_advection_info) then
               call mpas_log_write('Advect thickness and tracers, dt=$r', realArgs=(/deltat/))
            endif

            call li_advection_thickness_tracers(&
                 domain,                 &
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
               call mpas_log_write('Advect thickness (but not tracers), dt=$r', realArgs=(/deltat/))
            endif

            call li_advection_thickness_tracers(&
                 domain,                 &
                 deltat,                 &
                 meshPool,               &
                 velocityPool,           &
                 geometryPool,           &
                 thermalPool,            &
                 scratchPool,            &
                 err_tmp,                &
                 advectTracersIn = .false.)

                 !call mpas_log_write("errtmp=$i",intArgs=(/err_tmp/))
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
               call mpas_log_write('  Cells with negative thickness (set to 0): $i', intArgs=(/sum(masktmp)/))
            endif

            ! Count how many cells have ice.
            masktmp = 0
            where (thickness > 0.0_RKIND)
               masktmp = 1
            end where
            call mpas_log_write('  Cells with nonzero thickness: $i', intArgs=(/sum(masktmp)/))

         endif

         deallocate(masktmp)

         block => block % next
      end do

      ! Halo updates
      call mpas_timer_start("halo updates")

      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_dmpar_field_halo_exch(domain, 'temperature')
      call mpas_dmpar_field_halo_exch(domain, 'waterFrac')
      call mpas_dmpar_field_halo_exch(domain, 'enthalpy')
      call mpas_dmpar_field_halo_exch(domain, 'damage')
      call mpas_dmpar_field_halo_exch(domain, 'passiveTracer2d')

      call mpas_timer_stop("halo updates")

      ! Update mask and geometry
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         
         ! restore thickness if using config_restore_thickness_after_advection
         if ( config_restore_thickness_after_advection ) then
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            call mpas_pool_get_array(geometryPool, 'thicknessOld', thicknessOld)
            call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
            thickness(:) = thicknessOld(:)

            call li_calculate_layerThickness(meshPool, thickness, layerThickness)
         endif

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         call li_update_geometry(geometryPool)

         block => block % next
      end do


      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in advection_solver.", MPAS_LOG_ERR)
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
   subroutine set_timestep(allowableAdvecDt, allowableDiffDt, calvingCFLdt, faceMeltingCFLdt, clock, dtSeconds, CFLprocess, err)
      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: allowableAdvecDt
      real (kind=RKIND), intent(in) :: allowableDiffDt
      real (kind=RKIND), intent(in) :: calvingCFLdt, faceMeltingCFLdt
      type (MPAS_Clock_type), intent(in) :: clock

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(out) :: dtSeconds  !< Output: time step in seconds determined by this routine
      integer, intent(out) :: CFLprocess !< Output: flag for which process limits the CFL: 1=advective, 2=diffusive, 3=calving, 4=face-melting
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_adaptive_timestep_include_DCFL
      logical, pointer :: config_adaptive_timestep_include_calving, &
                          config_adaptive_timestep_include_face_melting
      real (kind=RKIND), pointer :: config_adaptive_timestep_CFL_fraction
      real (kind=RKIND), pointer :: config_adaptive_timestep_calvingCFL_fraction, &
                                    config_adaptive_timestep_faceMeltingCFL_fraction
      real (kind=RKIND), pointer :: config_max_adaptive_timestep
      real (kind=RKIND), pointer :: config_min_adaptive_timestep
      character (len=StrKIND), pointer :: config_calving, config_front_mass_bal_grounded
      character (len=StrKIND), pointer :: config_damage_calving_method
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
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_calving', config_adaptive_timestep_include_calving)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_face_melting', config_adaptive_timestep_include_face_melting)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_calvingCFL_fraction', config_adaptive_timestep_calvingCFL_fraction)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_faceMeltingCFL_fraction', config_adaptive_timestep_faceMeltingCFL_fraction)
      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)
      call mpas_pool_get_config(liConfigs, 'config_front_mass_bal_grounded', config_front_mass_bal_grounded)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_method', config_damage_calving_method)


      allowableDt = allowableAdvecDt * config_adaptive_timestep_CFL_fraction
      CFLprocess = 1
      if (config_adaptive_timestep_include_DCFL) then
         ! todo: We currently do not have a namelist option for the diffusive CFL fraction
         ! But that's ok - we rarely use it.  Could be added later if needed.
         if (allowableDiffDt < allowableDt) then
            allowableDt = allowableDiffDt
            CFLprocess = 2
         endif
      endif

      ! Only include calving CFL if requested and we are using a calving option that calculates it
      if (config_adaptive_timestep_include_calving .and. ( &
             trim(config_calving) == 'specified_calving_velocity' .or. &
             trim(config_calving) == 'eigencalving' .or. &
             trim(config_calving) == 'von_Mises_stress' .or. &
             trim(config_calving) == 'ismip6_retreat' .or. &
             (trim(config_calving) == 'damagecalving' .and. &
              trim(config_damage_calving_method) == 'calving_rate'))) then
         if (config_adaptive_timestep_calvingCFL_fraction * calvingCFLDt < allowableDt) then
            allowableDt = config_adaptive_timestep_calvingCFL_fraction * calvingCFLDt
            CFLprocess = 3
         endif
      endif

      ! Only include face-melting CFL if requested and we are using a calving option that calculates it.
      if (config_adaptive_timestep_include_face_melting .and. ( &
             trim(config_front_mass_bal_grounded) .ne. 'none' ) ) then
         if (config_adaptive_timestep_faceMeltingCFL_fraction * faceMeltingCFLdt < allowableDt) then
            allowableDt = config_adaptive_timestep_faceMeltingCFL_fraction * faceMeltingCFLdt
            CFLprocess = 4
         endif
      endif

      call mpas_log_write("CFL dt adjusted for fractions (days): advective=$r, diffusive=$r, calving=$r, face-melting=$r", &
         realArgs=(/allowableAdvecDt*config_adaptive_timestep_CFL_fraction/86400.0_RKIND, &
                    allowableDiffDt/86400.0_RKIND, &
                    calvingCFLdt*config_adaptive_timestep_calvingCFL_fraction/86400.0_RKIND, &
                    faceMeltingCFLdt*config_adaptive_timestep_faceMeltingCFL_fraction/86400.0_RKIND/) )
      if (CFLprocess == 1) then
         call mpas_log_write("Timestep limited by advective CFL condition.")
      elseif (CFLprocess == 2) then
         call mpas_log_write("Timestep limited by diffusive CFL condition.")
      elseif (CFLprocess == 3) then
         call mpas_log_write("Timestep limited by calving CFL condition.")
      elseif (CFLprocess == 4) then
         call mpas_log_write("Timestep limited by face-melting CFL condition.")
      endif

      ! Take minimum of the max adaptive timestep setting and the allowable dt from the CFL condition
      proposedDt = min(allowableDt, config_max_adaptive_timestep)

      ! Round down the proposed dt to avoid complications with fractional seconds
      ! (some timekeeping-related functionality, like restarts, don't support them)
      ! (need to perform floor with an 8-bit integer to allow up 293 billion years in seconds)
      proposedDt = real(floor(proposedDt, KIND=8), RKIND)

      ! Check if the proposed timestep is smaller than specified limit.
      ! Do this prior to limiting the timestep for the force interval
      if (proposedDt < config_min_adaptive_timestep) then
         call mpas_log_write('New deltat is less than config_min_adaptive_timestep.', MPAS_LOG_ERR)
         err = ior(err, 1)
      endif

      ! Check if we need to force a timestep length to hit the target interval
      currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
      !print *, 'curr', currTime % t % YR, currTime % t % basetime % S, currTime % t % basetime % Sn, currTime % t % basetime % Sd
      err = ior(err,err_tmp)
      nextForceTime = mpas_alarm_get_next_ring_time(clock, 'adaptiveTimestepForceInterval')
      !print *, 'ring', nextForceTime % t % YR, nextForceTime % t % basetime % S, nextForceTime % t % basetime % Sn, &
      !   nextForceTime % t % basetime % Sd
      intervalToNextForceTime = nextForceTime - currTime
      !print *, 'int', intervalToNextForceTime % ti % YR, intervalToNextForceTime % ti % MM, &
      !   intervalToNextForceTime % ti % basetime % S, intervalToNextForceTime % ti % basetime % Sn, &
      !   intervalToNextForceTime % ti % basetime % Sd
      ! Due to roundoff errors, we might be just shy of the desired time.
      ! To avoid this, add one to the numerator of the fractional seconds to
      ! make sure we get pushed over the edge.  The way ESMF does fractional
      ! seconds, this means we get the desired interval to better than 1 part per 100 million seconds
      ! Note that even though this is a *very* tiny fudge factor, it does not
      ! affect conservation within MPAS-LI, but it could have a very, very tiny
      ! effect on a climate model that thinks we ran for, say, 10 years, but we
      ! actually ran for 10 years +/- 1e-8 seconds.
      !intervalToNextForceTime % ti % basetime % Sn = intervalToNextForceTime % ti % basetime % Sn + 1
      ! Note: commenting above line because it should no longer be relevant after rounding down to full second on line 815,
      !       but we still include a check for fractional seconds.
      call mpas_get_timeInterval(intervalToNextForceTime, dt=secondsToNextForceTime, ierr=err_tmp)
      err = ior(err,err_tmp)
      if (secondsToNextForceTime - real(floor(secondsToNextForceTime, KIND=8), RKIND) /= 0.0_RKIND) then
         call mpas_log_write("set_timestep found secondsToNextForceTime not equal to 0.0: $r, decimal part=$r", &
            MPAS_LOG_ERR, realArgs=(/secondsToNextForceTime, secondsToNextForceTime - &
            real(floor(secondsToNextForceTime, KIND=8), RKIND)/))
         err = ior(err, 1)
      endif
      !print *, proposedDt, secondsToNextForceTime

      ! --- Actually set the dt here ---
      dtSeconds = min(proposedDt, secondsToNextForceTime)

      call mpas_log_write('  Setting time step (days) to: $r', realArgs=(/dtSeconds / (86400.0_RKIND)/))

   !--------------------------------------------------------------------
   end subroutine set_timestep


!***********************************************************************
!
!  routine calculate_layerThicknessEdge
!
!> \brief   Calculates layerThicknessEdge
!> \author  Matthew Hoffman
!> \date    23 May 2017
!> \details
!>  This routine calculates FO upwind thickness on 3d layer edges.
!
!-----------------------------------------------------------------------
   subroutine calculate_layerThicknessEdge(meshPool, geometryPool, velocityPool, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(inout) :: meshPool !< Input/Output: mesh pool
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input/Output: geometryPool
      type (mpas_pool_type), pointer, intent(inout) :: velocityPool !< Input/Output: velocityPool

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessEdge, normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nEdges, nVertLevels
      integer :: iEdge, cell1, cell2, k
      real (kind=RKIND) :: VelSign

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)

      ! Note: SIA velocity solver uses its own local calculation of h_edge that is always 2nd order.
      ! Note: ocn_diagnostic_solve in mpas_ocn_tendency.F has 2, 3, & 4th order calculations for h_edge that can be used.

      ! given thickness, compute layerThickness
      call li_calculate_layerThickness(meshPool, thickness, layerThickness)

      ! If using FO-Upwind then h_edge must be FO.
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1, nVertLevels
            ! Calculate h on edges using first order
            VelSign = sign(1.0_RKIND, normalVelocity(k, iEdge))
            layerThicknessEdge(k,iEdge) = max(VelSign * layerThickness(k, cell1), &
               VelSign * (-1.0_RKIND) * layerThickness(k, cell2))
            ! + velocity goes from index 1 to 2 in the cellsOnEdge array.
            !  Doug does the calculation as: h_edge = max(VelSign, 0.0) * h1 - min(VelSign, 0.0) * h2
            !!! ! Calculate h on edges using second order
            !!! layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k, cell1) + layerThickness(k, cell2))
         end do
      end do

      ! Note: The outmost layerThicknessEdge may be wrong if its upwind cell is off this block -
      !       halo update should be done if this variable will be used.

      !WHL - debug - Commented out, but might be useful for future debugging
!      if (config_print_thickness_advection_info) then
!         call mpas_log_write(' ')
!         call mpas_log_write('End of timestep: iCell (global), new thickness:'
!         do iCell = 1, nCells
!            if (thickness(iCell) > 0.0_RKIND) then
!               call mpas_log_write('$i $r', intArgs=(/indexToCellID(iCell)/), realArgs=(/thickness(iCell)/))
!            endif
!         enddo
!      endif

   !--------------------------------------------------------------------
   end subroutine calculate_layerThicknessEdge

!***********************************************************************
!
!  routine advance_clock
!
!> \brief   Update the clock
!> \author  Holly Han
!> \date    July 2023
!> \details
!>  This routine advances the clock and is called right after the
!>  subroutine 'prepare_advection' such that the clock is updated
!>  earlier in the timestep before other physics are calculated
!-----------------------------------------------------------------------

   subroutine advance_clock(domain, err)

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
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (MPAS_Time_Type) :: currTime  !< current time as time type
      type (MPAS_TimeInterval_type) :: timeStepInterval !< the current time step as an interval
      type (MPAS_Time_type) :: simulationStartTime_timeType

      logical, pointer :: config_adaptive_timestep

      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: simulationStartTime
      character(len=StrKIND) :: timeStamp !< current time as a string

      real (kind=RKIND), pointer :: daysSinceStart
      real (kind=RKIND), pointer :: deltat ! variable in blocks

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)

      ! ===
      ! === Adaptive timestep: update clock information
      ! ===
      ! Set time step in clock object since the time step could have changed
      ! Need to get value out of a block, but all blocks should have the same value, so just use first block
      if (config_adaptive_timestep) then
         block => domain % blocklist
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         ! convert dt in seconds to timeInterval type
         call mpas_set_timeInterval(timeStepInterval, dt=deltat, ierr=err_tmp)
         err = ior(err,err_tmp)
         ! update the clock with the timeInterval
         call mpas_set_clock_timestep(domain % clock, timeStepInterval, err_tmp)
         err = ior(err,err_tmp)
      endif

      ! ===
      ! === Update clock information
      ! ===
      ! Advance clock - needed to wait until after time step is completed in case the dt has changed!
      call mpas_advance_clock(domain % clock)
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_log_write('  Completed timestep.  New time is: ' // trim(timeStamp), flushNow=.true.)

      block => domain % blocklist
      do while (associated(block))
         ! Assign the time stamp for this time step
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'xtime', xtime)
         xtime = timeStamp

         ! compute time since start of simulation, in days
         call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
         call mpas_pool_get_array(meshPool, 'daysSinceStart',daysSinceStart)
         call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
         call mpas_get_timeInterval(currTime - simulationStartTime_timeType, dt=daysSinceStart)
         daysSinceStart = daysSinceStart / seconds_per_day

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in advance_clock.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
    end subroutine advance_clock


end module li_time_integration_fe_rk

