! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mixed_layer_heat_budget
!
!> \brief MPAS ocean analysis mode member: mixed_layer_heat_budget
!> \author Luke Van Roekel
!> \date   5/20/2020
!> \details
!>  MPAS ocean analysis mode member: mixed_layer_heat_budget
!>      Computes the heat budget integrated over the mixed layer
!-----------------------------------------------------------------------

module ocn_mixed_layer_heat_budget

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_log

   use ocn_constants
   use ocn_diagnostics_variables

   implicit none
   private
   save

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

   public :: ocn_init_mixed_layer_heat_budget, &
             ocn_compute_mixed_layer_heat_budget, &
             ocn_restart_mixed_layer_heat_budget, &
             ocn_finalize_mixed_layer_heat_budget

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_mixed_layer_heat_budget
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author Luke Van Roekel
!> \date    5/20/2020
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_mixed_layer_heat_budget(domain, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_mixed_layer_heat_budget !}}}

!***********************************************************************
!
!  routine ocn_compute_mixed_layer_heat_budget
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    5/20/2020
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_mixed_layer_heat_budget(domain, timeLevel, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type(mpas_pool_type), pointer :: mixedLayerHeatBudgetAMPool
      type(mpas_pool_type), pointer :: mixedLayerDepthsAMPool
      type(dm_info) :: dminfo
      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool
      type(mpas_pool_type), pointer :: tracerTendPool
      type(mpas_pool_type), pointer :: tendPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve
      integer :: iTracer, k, iCell, num_tracers
      integer, dimension(:), pointer :: maxLevelCell
      integer, pointer :: index_temperature

      real(kind=RKIND), dimension(:, :, :), pointer :: &
          tracers, tracerTend
      real(kind=RKIND), dimension(:, :), pointer :: activeTracerForcingMLTend, &
          activeTracerNonLocalMLTend, activeTracerHorMixMLTend, activeTracerVertMixMLTend, &
          activeTracerHorAdvectionMLTend, activeTracerVertAdvectionMLTend, activeTracersML, &
          layerThickness, activeTracersTendML
      real(kind=RKIND), dimension(:), pointer :: bruntVaisalaFreqML
      real(kind=RKIND) :: depth, difference
      err = 0

      dminfo = domain%dminfo

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block%structs, 'mixedLayerHeatBudgetAM', mixedLayerHeatBudgetAMPool)
         call mpas_pool_get_subpool(block%structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
         call mpas_pool_get_subpool(block%structs, 'tend', tendPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracerTendPool)

         call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block%dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerForcingMLTend', &
                                  activeTracerForcingMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerNonLocalMLTend', &
                                  activeTracerNonLocalMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerHorMixMLTend', &
                                  activeTracerHorMixMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerVertMixMLTend', &
                                  activeTracerVertMixMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerHorAdvectionMLTend', &
                                  activeTracerHorAdvectionMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracerVertAdvectionMLTend', &
                                  activeTracerVertAdvectionMLTend)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracersML', activeTracersML)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'activeTracersTendML', activeTracersTendML)
         call mpas_pool_get_array(mixedLayerHeatBudgetAMPool, 'bruntVaisalaFreqML', bruntVaisalaFreqML)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(tracerTendPool, 'activeTracersTend', tracerTend)

         num_tracers = size(activeTracerVertMixTendency, dim=1)

         do iCell = 1, nCellsSolve
            activeTracerHorAdvectionMLTend(:, iCell) = 0.0_RKIND
            activeTracerVertAdvectionMLTend(:, iCell) = 0.0_RKIND
            activeTracerHorMixMLTend(:, iCell) = 0.0_RKIND
            activeTracerVertMixMLTend(:, iCell) = 0.0_RKIND
            activeTracerNonLocalMLTend(:, iCell) = 0.0_RKIND
            activeTracerForcingMLTend(:, iCell) = 0.0_RKIND
            activeTracersML(:, iCell) = 0.0_RKIND
            activeTracersTendML(:,iCell) = 0.0_RKIND
            k = 1
            depth = layerThickness(k, iCell)
            do while (depth <= dThreshMLD(iCell) .and. k < maxLevelCell(iCell))
               do iTracer = 1, num_tracers
                  activeTracerHorAdvectionMLTend(iTracer, iCell) = activeTracerHorAdvectionMLTend(iTracer, iCell) + &
                                                 activeTracerHorizontalAdvectionTendency(iTracer, k, iCell)*layerThickness(k, iCell)
                  activeTracerVertAdvectionMLTend(iTracer, iCell) = activeTracerVertAdvectionMLTend(iTracer, iCell) + &
                                                   activeTracerVerticalAdvectionTendency(iTracer, k, iCell)*layerThickness(k, iCell)

                  activeTracerHorMixMLTend(iTracer, iCell) = activeTracerHorMixMLTend(iTracer, iCell) + &
                                                             activeTracerHorMixTendency(iTracer, k, iCell)*layerThickness(k, iCell)
                  activeTracerVertMixMLTend(iTracer, iCell) = activeTracerVertMixMLTend(iTracer, iCell) + &
                                                             activeTracerVertMixTendency(iTracer, k, iCell)*layerThickness(k, iCell)
                  activeTracerNonLocalMLTend(iTracer, iCell) = activeTracerNonLocalMLTend(iTracer, iCell) + &
                                                            activeTracerNonLocalTendency(iTracer, k, iCell)*layerThickness(k, iCell)

                  activeTracerForcingMLTend(iTracer, iCell) = activeTracerForcingMLTend(iTracer, iCell) + &
                                                         activeTracerSurfaceFluxTendency(iTracer, k, iCell)*layerThickness(k, iCell)
                  if (iTracer == index_temperature) THEN
                     activeTracerForcingMLTend(iTracer, iCell) = activeTracerForcingMLTend(iTracer, iCell) + &
                                                                 temperatureShortWaveTendency(k, iCell)*layerThickness(k, iCell)
                  endif
                  activeTracersML(iTracer, iCell) = activeTracersML(iTracer, iCell) + &
                                                    tracers(iTracer, k, iCell)*layerThickness(k, iCell)
                  !Note we don't multiply this term by layerThickness because that term is in 
                  !Tendency already
                  activeTracersTendML(iTracer, iCell) = activeTracersTendML(iTracer, iCell) + &
                                                    tracerTend(iTracer, k, iCell)
               end do
               bruntVaisalaFreqML(iCell) = bruntVaisalaFreqML(iCell) + BruntVaisalaFreqTop(k, iCell)*layerThickness(k, iCell)
               k = k + 1
               depth = depth + layerThickness(k, iCell)
            end do
            !add partial thickness of last layer to get to MLD
            depth = depth - layerThickness(k, iCell)
            difference = dThreshMLD(iCell) - depth
            do iTracer = 1, num_tracers
               activeTracerHorAdvectionMLTend(iTracer, iCell) = activeTracerHorAdvectionMLTend(iTracer, iCell) + &
                                                               activeTracerHorizontalAdvectionTendency(iTracer, k, iCell)*difference
               activeTracerVertAdvectionMLTend(iTracer, iCell) = activeTracerVertAdvectionMLTend(iTracer, iCell) + &
                                                                 activeTracerVerticalAdvectionTendency(iTracer, k, iCell)*difference
               activeTracerHorMixMLTend(iTracer, iCell) = activeTracerHorMixMLTend(iTracer, iCell) + &
                                                          activeTracerHorMixTendency(iTracer, k, iCell)*difference
               activeTracerVertMixMLTend(iTracer, iCell) = activeTracerVertMixMLTend(iTracer, iCell) + &
                                                           activeTracerVertMixTendency(iTracer, k, iCell)*difference
               activeTracerNonLocalMLTend(iTracer, iCell) = activeTracerNonLocalMLTend(iTracer, iCell) + &
                                                            activeTracerNonLocalTendency(iTracer, k, iCell)*difference
               activeTracerForcingMLTend(iTracer, iCell) = activeTracerForcingMLTend(iTracer, iCell) + &
                                                           activeTracerSurfaceFluxTendency(iTracer, k, iCell)*difference
               if (iTracer == index_temperature) THEN
                  activeTracerForcingMLTend(iTracer, iCell) = activeTracerForcingMLTend(iTracer, iCell) + &
                                                              temperatureShortWaveTendency(k, iCell)*difference
               endif
               activeTracersML(iTracer, iCell) = activeTracersML(iTracer, iCell) + tracers(iTracer, k, iCell)* &
                                                 difference 
               activeTracersTendML(iTracer, iCell) = activeTracersTendML(iTracer, iCell) + tracerTend(iTracer,k,iCell)* &
                                                     difference / layerThickness(k,iCell) 

               activeTracersML(iTracer, iCell) = activeTracersML(iTracer, iCell)/dThreshMLD(iCell)
               activeTracersTendML(iTracer, iCell) = activeTracersTendML(iTracer, iCell)/dThreshMLD(iCell)
               activeTracerHorAdvectionMLTend(iTracer, iCell) = activeTracerHorAdvectionMLTend(iTracer, iCell)/dThreshMLD(iCell)
               activeTracerVertAdvectionMLTend(iTracer, iCell) = activeTracerVertAdvectionMLTend(iTracer, iCell) &
                                                                 / dThreshMLD(iCell)
               activeTracerHorMixMLTend(iTracer, iCell) = activeTracerHorMixMLTend(iTracer, iCell)/dThreshMLD(iCell)
               activeTracerVertMixMLTend(iTracer, iCell) = activeTracerVertMixMLTend(iTracer, iCell)/dThreshMLD(iCell)
               activeTracerNonLocalMLTend(iTracer, iCell) = activeTracerNonLocalMLTend(iTracer, iCell)/dThreshMLD(iCell)
               activeTracerForcingMLTend(iTracer, iCell) = activeTracerForcingMLTend(iTracer, iCell)/dThreshMLD(iCell)
            end do
            bruntVaisalaFreqML(iCell) = bruntVaisalaFreqML(iCell) + BruntVaisalaFreqTop(k, iCell)*difference
            bruntVaisalaFreqML(iCell) = bruntVaisalaFreqML(iCell)/dThreshMLD(iCell)
         end do

         block => block%next
      end do

   end subroutine ocn_compute_mixed_layer_heat_budget !}}}

!***********************************************************************
!
!  routine ocn_restart_mixed_layer_heat_budget
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    5/20/2020
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_mixed_layer_heat_budget(domain, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_mixed_layer_heat_budget !}}}

!***********************************************************************
!
!  routine ocn_finalize_mixed_layer_heat_budget
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    5/20/2020
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_mixed_layer_heat_budget(domain, err) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_mixed_layer_heat_budget !}}}

end module ocn_mixed_layer_heat_budget

! vim: foldmethod=marker
