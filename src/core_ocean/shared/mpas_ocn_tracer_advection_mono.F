! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_mono
!
!> \brief MPAS monotonic tracer advection with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen
!> \date   October 2017
!> \details
!>  This module contains routines for monotonic advection of tracers using a FCT
!
!-----------------------------------------------------------------------
module ocn_tracer_advection_mono

#ifdef _ADV_TIMERS
   use mpas_timer
#endif
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use mpas_tracer_advection_helpers

   use ocn_constants

   implicit none
   private
   save

   real (kind=RKIND) :: coef_3rd_order
   integer :: horizOrder
   logical :: vert2ndOrder, vert3rdOrder, vert4thOrder
   logical :: positiveDzDk, monotonicityCheck

   public :: ocn_tracer_advection_mono_tend, &
             ocn_tracer_advection_mono_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_mono_tend
!
!> \brief MPAS monotonic tracer horizontal advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen
!> \date   October 2017
!> \details
!>  This routine computes the monotonic tracer horizontal advection tendencity using a FCT.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_mono_tend(tracers, adv_coefs, adv_coefs_3rd, nAdvCellsForEdge, advCellsForEdge, &!{{{
                                              normalThicknessFlux, w, layerThickness, dt, meshPool, &
                                              scratchPool, diagnosticsPool, tend, maxLevelCell, maxLevelEdgeTop, &
                                              highOrderAdvectionMask, edgeSignOnCell, tracerGroupName)

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs !< Input: Advection coefficients for 2nd order advection
      real (kind=RKIND), dimension(:,:), intent(in) :: adv_coefs_3rd !< Input: Advection coeffs for mising 3rd/4th order advection
      integer, dimension(:), intent(in) :: nAdvCellsForEdge !< Input: Number of advection cells for each edge
      integer, dimension(:,:), intent(in) :: advCellsForEdge !< Input: List of advection cells for each edge
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: Scratch fields
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: pool for traceradvection budget term
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency

      ! These variables are passed in to match the version of the interface in operators.
      integer, dimension(:), pointer :: maxLevelCell !< Input: Index to max level at cell center
      integer, dimension(:), pointer :: maxLevelEdgeTop !< Input: Index to max level at edge with non-land cells on both sides
      integer, dimension(:,:), pointer :: highOrderAdvectionMask !< Input: Mask for high order advection
      integer, dimension(:, :), pointer :: edgeSignOnCell !< Input: Sign for flux from edge on each cell.

      character (len=*), intent(in) :: tracerGroupName ! variable to check for tracer budget

      logical, pointer :: config_compute_active_tracer_budgets

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2, nVertLevels, num_tracers, nCells, nEdges
      integer, pointer :: maxEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell

      real (kind=RKIND) :: signedFactor, tracer_new 
      real (kind=RKIND) :: flux_upwind, tracer_min_new, tracer_max_new, tracer_upwind_new, scale_factor
      real (kind=RKIND) :: flux, tracer_weight, invAreaCell1, invAreaCell2
      real (kind=RKIND) :: verticalWeightK, verticalWeightKm1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: tracer_cur, work_tend, h_new_inv
      real (kind=RKIND), dimension(:,:), pointer :: tracer_max, tracer_min
      real (kind=RKIND), dimension(:,:), pointer :: flux_incoming, flux_outgoing
      real (kind=RKIND), dimension(:,:), pointer :: high_order_flux
      real (kind=RKIND), dimension(:,:), pointer :: h_prov, h_prov_inv

      type (field2DReal), pointer :: tracerCurField, workTendencyField, hNewInvField, tracerMinField, &
                                     tracerMaxField, fluxIncomingField, fluxOutgoingField, &
                                     hProvInvField, hProvField, highOrderFluxField

      real (kind=RKIND), dimension(:,:,:), pointer :: &
              activeTracerHorizontalAdvectionTendency, &
              activeTracerVerticalAdvectionTendency

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

#ifdef _ADV_TIMERS
      call mpas_timer_start('startup')
#endif
      ! Get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      nVertLevels = size(tracers,dim=2)
      num_tracers = size(tracers,dim=1)

      ! Initialize pointers
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

      call mpas_pool_get_config(ocnConfigs, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)
      if (config_compute_active_tracer_budgets) then
         call mpas_pool_get_array(diagnosticsPool,'activeTracerHorizontalAdvectionTendency', &
                 activeTracerHorizontalAdvectionTendency)
         call mpas_pool_get_array(diagnosticsPool,'activeTracerVerticalAdvectionTendency', &
                 activeTracerVerticalAdvectionTendency)
      end if

      call mpas_pool_get_field(scratchPool, 'tracerCur', tracerCurField)
      call mpas_pool_get_field(scratchPool, 'workTendency', workTendencyField)
      call mpas_pool_get_field(scratchPool, 'hNewInv', hNewInvField)
      call mpas_pool_get_field(scratchPool, 'tracerMin', tracerMinField)
      call mpas_pool_get_field(scratchPool, 'tracerMax', tracerMaxField)
      call mpas_pool_get_field(scratchPool, 'fluxIncoming', fluxIncomingField)
      call mpas_pool_get_field(scratchPool, 'fluxOutgoing', fluxOutgoingField)
      call mpas_pool_get_field(scratchPool, 'hProvInv', hProvInvField)
      call mpas_pool_get_field(scratchPool, 'hProv', hProvField)
      call mpas_pool_get_field(scratchPool, 'highOrderFlux', highOrderFluxField)

      call mpas_allocate_scratch_field(tracerCurField, .true., .false.)
      call mpas_allocate_scratch_field(workTendencyField, .true., .false.)
      call mpas_allocate_scratch_field(hNewInvField, .true., .false.)
      call mpas_allocate_scratch_field(tracerMinField, .true., .false.)
      call mpas_allocate_scratch_field(tracerMaxField, .true., .false.)
      call mpas_allocate_scratch_field(fluxIncomingField, .true., .false.)
      call mpas_allocate_scratch_field(fluxOutgoingField, .true., .false.)
      call mpas_allocate_scratch_field(hProvInvField, .true., .false.)
      call mpas_allocate_scratch_field(hProvField, .true., .false.)
      call mpas_allocate_scratch_field(highOrderFluxField, .true., .false.)
      call mpas_threading_barrier()

      ! allocate nCells arrays
      h_prov_inv => hProvInvField % array
      h_prov => hProvField % array
      h_new_inv => hNewInvField % array
      tracer_cur => tracerCurField % array
      tracer_max => tracerMaxField % array
      tracer_min => tracerMinField % array
      work_tend => workTendencyField % array
      flux_incoming => fluxIncomingField % array
      flux_outgoing => fluxOutgoingField % array
      high_order_flux => highOrderFluxField % array

      nCells = nCellsArray( size(nCellsArray) )

      ! Note: This assumes we are in the first part of the horizontal/
      ! vertical operator splitting, which is true because currently
      ! we dont flip order and horizontal is always first.
      ! See notes in commit 2cd4a89d.
      !$omp do schedule(runtime) private(k, i, iEdge, invAreaCell1)
      do iCell = 1, nCells
        invAreaCell1 = 1.0_RKIND / areaCell(iCell)
        do k = 1, maxLevelCell(iCell)
          h_prov(k, iCell) = layerThickness(k, iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            ! Provisional layer thickness is after horizontal thickness flux only
            h_prov(k, iCell) = h_prov(k, iCell) + dt * invAreaCell1 * dvEdge(iEdge) * &
                                            edgeSignOnCell(i, iCell) * normalThicknessFlux(k, iEdge)
          end do
          h_prov_inv(k, iCell) = 1.0_RKIND / h_prov(k, iCell)
          ! New layer thickness is after horizontal and vertical thickness flux
          h_new_inv(k, iCell) = 1.0_RKIND / (h_prov(k, iCell) - dt * w(k, iCell) + dt * w(k+1, iCell))
        end do
      end do
      !$omp end do

#ifdef _ADV_TIMERS
      call mpas_timer_stop('startup')
#endif

      ! Loop over tracers. One tracer is advected at a time. It is copied into a temporary array in order to improve locality
      do iTracer = 1, num_tracers
#ifdef _ADV_TIMERS
        call mpas_timer_start('cell init')
#endif
        nCells = nCellsArray( size(nCellsArray) )
        ! Initialize variables for use in this iTracer iteration
        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCells
          do k=1, nVertLevels
            tracer_cur(k,iCell) = tracers(iTracer,k,iCell)
            work_tend(k, iCell) = 0.0_RKIND
            flux_incoming(k, iCell) = 0.0_RKIND
            flux_outgoing(k, iCell) = 0.0_RKIND
          end do ! k loop
        end do ! iCell loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('cell init')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('horiz flux')
#endif
        nCells = nCellsArray( 2 )
        !  Compute the high and low order vertical fluxes. Also determine bounds on tracer_cur.
        !$omp do schedule(runtime) private(k, i)
        do iCell = 1, nCells
          do k=1, maxLevelCell(iCell)
            tracer_min(k,iCell) = tracer_cur(k,iCell)
            tracer_max(k,iCell) = tracer_cur(k,iCell)
          end do
          ! pull tracer_min and tracer_max from the (horizontal) surrounding cells
          do i = 1, nEdgesOnCell(iCell)
            do k=1, min(maxLevelCell(iCell), maxLevelCell(cellsOnCell(i,iCell)))
              tracer_max(k,iCell) = max(tracer_max(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
              tracer_min(k,iCell) = min(tracer_min(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
            end do ! k loop
          end do ! i loop over nEdgesOnCell
        end do
        !$omp end do

        ! Need all the edges around the 1 halo cells and owned cells
        nEdges = nEdgesArray( 3 )
        !  Compute the high order horizontal flux
        !$omp do schedule(runtime) private(cell1, cell2, k, tracer_weight, i, iCell, flux_upwind)
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          do k = 1, nVertLevels
            high_order_flux(k, iEdge) = 0.0_RKIND
          end do

          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 1, maxLevelCell(iCell)
              tracer_weight = highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) &
                            + coef_3rd_order*sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))

              tracer_weight = normalThicknessFlux(k,iEdge)*tracer_weight
              high_order_flux(k,iEdge) = high_order_flux(k,iEdge) + tracer_weight * tracer_cur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge

          ! Compute 2nd order fluxes where needed.
          ! Also compute low order upwind horizontal flux (monotonic and diffused)
          ! Remove low order flux from the high order flux
          ! Store left over high order flux in high_order_flux array
          do k = 1, maxLevelEdgeTop(iEdge)
            tracer_weight = iand(highOrderAdvectionMask(k, iEdge)+1, 1) * (dvEdge(iEdge) * 0.5_RKIND) &
                          * normalThicknessFlux(k, iEdge)

            high_order_flux(k, iEdge) = high_order_flux(k, iedge) + tracer_weight * (tracer_cur(k, cell1) &
                                            + tracer_cur(k, cell2))

            flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) &
                        + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))
            high_order_flux(k,iEdge) = high_order_flux(k,iEdge) - flux_upwind
          end do ! k loop
        end do ! iEdge loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('horiz flux')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('scale factor build')
#endif
        ! Need one halo of cells around owned cells
        nCells = nCellsArray( 2 )
        !$omp do schedule(runtime) private(k, tracer_max_new, tracer_min_new, tracer_upwind_new, scale_factor, invAreaCell1, i, &
        !$omp                              iEdge, cell1, cell2, flux_upwind, signedFactor)
        do iCell = 1, nCells
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)

          ! Finish computing the low order horizontal fluxes
          ! Upwind fluxes are accumulated in work_tend
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            signedFactor = edgeSignOnCell(i, iCell) * invAreaCell1

            do k = 1, maxLevelEdgeTop(iEdge)
              flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) &
                                          + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))

              ! Here work_tend is the upwind tendency
              work_tend(k,iCell) = work_tend(k,iCell) + signedFactor * flux_upwind

              ! Accumulate remaining high order fluxes
              flux_outgoing(k,iCell) = flux_outgoing(k,iCell) + min(0.0_RKIND, signedFactor &
                                     * high_order_flux(k, iEdge))
              flux_incoming(k,iCell) = flux_incoming(k,iCell) + max(0.0_RKIND, signedFactor &
                                     * high_order_flux(k, iEdge))
            end do
          end do

          ! Build the factors for the FCT
          ! Computed using the bounds that were computed previously, and the bounds on the newly updated value
          ! Factors are placed in the flux_incoming and flux_outgoing arrays
          do k = 1, maxLevelCell(iCell)
            ! Here work_tend is the upwind tendency
            tracer_min_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(work_tend(k,iCell)+flux_outgoing(k,iCell))) &
                           * h_prov_inv(k,iCell)
            tracer_max_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(work_tend(k,iCell)+flux_incoming(k,iCell))) &
                           * h_prov_inv(k,iCell)
            tracer_upwind_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*work_tend(k,iCell)) * h_prov_inv(k,iCell)

            scale_factor = (tracer_max(k,iCell)-tracer_upwind_new)/(tracer_max_new-tracer_upwind_new+eps)
            flux_incoming(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

            scale_factor = (tracer_upwind_new-tracer_min(k,iCell))/(tracer_upwind_new-tracer_min_new+eps)
            flux_outgoing(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          end do ! k loop
        end do ! iCell loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('scale factor build')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('rescale horiz fluxes')
#endif
        ! Need all of the edges around owned cells
        nEdges = nEdgesArray( 2 )
        !  rescale the high order horizontal fluxes
        !$omp do schedule(runtime) private(cell1, cell2, k, flux)
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k = 1, maxLevelEdgeTop(iEdge)
            flux = high_order_flux(k,iEdge)
            flux = max(0.0_RKIND,flux) * min(flux_outgoing(k,cell1), flux_incoming(k,cell2)) &
                 + min(0.0_RKIND,flux) * min(flux_incoming(k,cell1), flux_outgoing(k,cell2))
            high_order_flux(k,iEdge) = flux
          end do ! k loop
        end do ! iEdge loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('rescale horiz fluxes')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('flux accumulate')
#endif

        nCells = nCellsArray( 1 )
        ! Accumulate the scaled high order vertical tendencies, and the upwind tendencies
        !$omp do schedule(runtime) private(invAreaCell1, signedFactor, i, iEdge, flux, k)
        do iCell = 1, nCells
          invAreaCell1 = 1.0_RKIND / areaCell(iCell)

          ! Accumulate the scaled high order horizontal tendencies
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            signedFactor = invAreaCell1 * edgeSignOnCell(i, iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
              ! work_tend on the RHS is the upwind tendency
              ! work_tend on the LHS is the total horizontal advection tendency
              work_tend(k, iCell) = work_tend(k, iCell) + signedFactor * high_order_flux(k, iEdge)
            end do
          end do

          do k = 1, maxLevelCell(iCell)
            ! work_tend on the RHS is the total horizontal advection tendency
            ! tracer_cur on LHS is the  provisional tracer after horizontal fluxes only.
            tracer_cur(k,iCell) = (tracer_cur(k, iCell)*layerThickness(k, iCell) + dt*work_tend(k, iCell)) &
                                  * h_prov_inv(k, iCell)
            tend(iTracer,k,iCell) = tend(iTracer,k,iCell) + work_tend(k,iCell)
          end do

        end do ! iCell loop
        !$omp end do

        if (config_compute_active_tracer_budgets) then
           if (tracerGroupName == 'activeTracers') then
              !$omp do schedule(runtime) private(k)
              do iCell = 1, nCells
                do k = 1, maxLevelCell(iCell)
                    activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = work_tend(k,iCell)
                end do
              end do ! iCell loop
              !$omp end do
           end if
        end if

#ifdef _ADV_TIMERS
        call mpas_timer_stop('flux accumulate')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('monotonic check')
#endif
        if (monotonicityCheck) then
          nCells = nCellsArray( 1 )
          !build min and max bounds on old and new tracer for check on monotonicity.
          !$omp do schedule(runtime) private(k)
          do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
              if(tracer_cur(k,iCell) < tracer_min(k, iCell)-eps) then
                 call mpas_log_write( &
                    'Horizontal minimum out of bounds on tracer: $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer/), realArgs=(/ tracer_min(k, iCell), tracer_cur(k,iCell) /) )
              end if

              if(tracer_cur(k,iCell) > tracer_max(k,iCell)+eps) then
                 call mpas_log_write( &
                    'Horizontal maximum out of bounds on tracer: $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer/), realArgs=(/ tracer_max(k, iCell), tracer_cur(k,iCell) /) )
              end if
            end do
          end do
          !$omp end do
        end if
#ifdef _ADV_TIMERS
        call mpas_timer_stop('monotonic check')
#endif

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  Vertical advection
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#ifdef _ADV_TIMERS
        call mpas_timer_start('cell init')
#endif
        nCells = nCellsArray( size(nCellsArray) )
        ! Initialize variables for use in this iTracer iteration
        !$omp do schedule(runtime) private(k)
        do iCell = 1, nCells
          do k=1, nVertLevels
            high_order_flux(k, iCell) = 0.0_RKIND
            work_tend(k, iCell) = 0.0_RKIND
          end do ! k loop
          high_order_flux(nVertLevels+1, iCell) = 0.0_RKIND

        end do ! iCell loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('cell init')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('vertical flux')
#endif

        ! Need all owned and 1 halo cells
        nCells = nCellsArray( 2 )
        !  Compute the high and low order vertical fluxes. Also determine bounds on tracer_cur.
        !$omp do schedule(runtime) private(k, verticalWeightK, verticalWeightKm1, i, flux_upwind)
        do iCell = 1, nCells

          ! Operate on top cell in column
          k = 1
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k+1,iCell))

          ! Operate on next-to-top cell in column
          k = max(1, min(maxLevelCell(iCell), 2))
          if ( k >= 2 ) then
            verticalWeightK = h_prov(k-1, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
            verticalWeightKm1 = h_prov(k, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
            high_order_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1 &
                *tracer_cur(k-1,iCell))
            tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
            tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          end if

          ! Operate on internal cells in column
          if ( vert4thOrder ) then
             do k=3,maxLevelCell(iCell)-1
                high_order_flux(k, iCell) = mpas_tracer_advection_vflux4( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                       tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell))

                tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
                tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             end do
          else if ( vert3rdOrder ) then
             do k=3,maxLevelCell(iCell)-1
                high_order_flux(k, iCell) = mpas_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                       tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell), coef_3rd_order )

                tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
                tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             end do
          else if ( vert2ndOrder ) then
             do k=3,maxLevelCell(iCell)-1
                verticalWeightK = h_prov(k-1, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
                verticalWeightKm1 = h_prov(k, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
                high_order_flux(k,iCell) = w(k,iCell) * (verticalWeightK * tracer_cur(k,iCell) + verticalWeightKm1 &
                                              * tracer_cur(k-1,iCell))

                tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
                tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             end do
          end if

          ! Operate on deepest vertical cell in column.
          k = max(1, maxLevelCell(iCell))
          verticalWeightK = h_prov(k-1, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
          verticalWeightKm1 = h_prov(k, iCell) / (h_prov(k, iCell) + h_prov(k-1, iCell))
          high_order_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k-1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k-1,iCell))

          !  low order upwind vertical flux (monotonic and diffused)
          !  Remove low order flux from the high order flux.
          !  Store left over high order flux in high_order_flux array.
          !  Upwind fluxes are accumulated in work_tend
          do k = 2, maxLevelCell(iCell)
            flux_upwind = min(0.0_RKIND,w(k,iCell))*tracer_cur(k-1,iCell) + max(0.0_RKIND,w(k,iCell))*tracer_cur(k,iCell)
            work_tend(k-1,iCell) = work_tend(k-1,iCell) + flux_upwind
            work_tend(k  ,iCell) = work_tend(k  ,iCell) - flux_upwind
            high_order_flux(k,iCell) = high_order_flux(k,iCell) - flux_upwind
          end do ! k loop

          ! flux_incoming contains the total remaining high order flux into iCell
          !          it is positive.
          ! flux_outgoing contains the total remaining high order flux out of iCell
          !           it is negative
          do k = 1, maxLevelCell(iCell)
            flux_incoming(k, iCell) = max(0.0_RKIND, high_order_flux(k+1, iCell)) &
                                    - min(0.0_RKIND, high_order_flux(k, iCell))
            flux_outgoing(k, iCell) = min(0.0_RKIND, high_order_flux(k+1, iCell)) &
                                    - max(0.0_RKIND, high_order_flux(k, iCell))
          end do ! k Loop
        end do ! iCell Loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('vertical flux')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('scale factor build')
#endif
        ! Need one halo of cells around owned cells
        nCells = nCellsArray( 2 )
        !$omp do schedule(runtime) private(k, tracer_max_new, tracer_min_new, tracer_upwind_new, scale_factor, i)
        do iCell = 1, nCells

          ! Build the factors for the FCT
          ! Computed using the bounds that were computed previously, and the bounds on the newly updated value
          ! Factors are placed in the flux_incoming and flux_outgoing arrays
          do k = 1, maxLevelCell(iCell)
            ! work_tend on the RHS is the upwind tendency
            tracer_min_new = (tracer_cur(k,iCell)*h_prov(k,iCell) + dt*(work_tend(k,iCell)+flux_outgoing(k,iCell))) &
                           * h_new_inv(k,iCell)
            tracer_max_new = (tracer_cur(k,iCell)*h_prov(k,iCell) + dt*(work_tend(k,iCell)+flux_incoming(k,iCell))) &
                           * h_new_inv(k,iCell)
            tracer_upwind_new = (tracer_cur(k,iCell)*h_prov(k,iCell) + dt*work_tend(k,iCell)) * h_new_inv(k,iCell)

            scale_factor = (tracer_max(k,iCell)-tracer_upwind_new)/(tracer_max_new-tracer_upwind_new+eps)
            flux_incoming(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

            scale_factor = (tracer_upwind_new-tracer_min(k,iCell))/(tracer_upwind_new-tracer_min_new+eps)
            flux_outgoing(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          end do ! k loop
        end do ! iCell loop
        !$omp end do
#ifdef _ADV_TIMERS
        call mpas_timer_stop('scale factor build')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('flux accumulate')
#endif

        nCells = nCellsArray( 1 )
        ! Accumulate the scaled high order vertical tendencies, and the upwind tendencies
        !$omp do schedule(runtime) private(flux, k)
        do iCell = 1, nCells
          ! rescale the high order vertical flux
          do k = 2, maxLevelCell(iCell)
            flux =  high_order_flux(k,iCell)
            flux = max(0.0_RKIND,flux) * min(flux_outgoing(k  ,iCell), flux_incoming(k-1,iCell)) &
                 + min(0.0_RKIND,flux) * min(flux_outgoing(k-1,iCell), flux_incoming(k  ,iCell))
            high_order_flux(k,iCell) = flux
          end do ! k loop

          do k = 1,maxLevelCell(iCell)
            ! work_tend on the RHS is the upwind tendency
            ! work_tend on the LHS is the total vertical advection tendency
            work_tend(k, iCell) = work_tend(k, iCell) + (high_order_flux(k+1, iCell) &
                                  - high_order_flux(k, iCell))
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + work_tend(k, iCell) 
          end do ! k loop

        end do ! iCell loop
        !$omp end do

        if (config_compute_active_tracer_budgets) then
           if (tracerGroupName == 'activeTracers') then
              !$omp do schedule(runtime) private(k)
              do iCell = 1, nCells
                do k = 1, maxLevelCell(iCell)
                    activeTracerVerticalAdvectionTendency(iTracer,k,iCell) = work_tend(k,iCell)
                end do
              end do ! iCell loop
              !$omp end do
           end if
        end if
#ifdef _ADV_TIMERS
        call mpas_timer_stop('flux accumulate')
#endif

#ifdef _ADV_TIMERS
        call mpas_timer_start('monotonic check')
#endif
        if (monotonicityCheck) then
          nCells = nCellsArray( 1 )
          ! Check for monotonicity of new tracer value
          !$omp do schedule(runtime) private(k)
          do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
              ! work_tend on the RHS is the total vertical advection tendency
              tracer_new = (tracer_cur(k, iCell)*h_prov(k, iCell) + dt * work_tend(k, iCell)) * h_new_inv(k, iCell)

              if(tracer_new < tracer_min(k, iCell)-eps) then
                 call mpas_log_write( &
                    'Vertical minimum out of bounds on tracer: $i $i $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer, k, iCell/), realArgs=(/ tracer_min(k, iCell), tracer_new /) )
              end if

              if(tracer_new > tracer_max(k,iCell)+eps) then
                 call mpas_log_write( &
                    'Vertical maximum out of bounds on tracer: $i $i $i $r $r ', &
                    MPAS_LOG_WARN, intArgs=(/iTracer, k, iCell/), realArgs=(/ tracer_max(k, iCell), tracer_new /) )
              end if
            end do
          end do
          !$omp end do
        end if
#ifdef _ADV_TIMERS
        call mpas_timer_stop('monotonic check')
#endif
      end do ! iTracer loop

#ifdef _ADV_TIMERS
      call mpas_timer_start('half step')
#endif
#ifdef _ADV_TIMERS
      call mpas_timer_stop('half step')
#endif

#ifdef _ADV_TIMERS
      call mpas_timer_start('deallocates')
#endif
      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(tracerCurField, .true.)
      call mpas_deallocate_scratch_field(workTendencyField, .true.)
      call mpas_deallocate_scratch_field(hNewInvField, .true.)
      call mpas_deallocate_scratch_field(fluxIncomingField, .true.)
      call mpas_deallocate_scratch_field(fluxOutgoingField, .true.)
      call mpas_deallocate_scratch_field(hProvInvField, .true.)
      call mpas_deallocate_scratch_field(hProvField, .true.)
      call mpas_deallocate_scratch_field(highOrderFluxField, .true.)

#ifdef _ADV_TIMERS
      call mpas_timer_stop('deallocates')
#endif

   end subroutine ocn_tracer_advection_mono_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_mono_init
!
!> \brief MPAS initialize monotonic tracer advection tendency with FCT
!> \author Mark Petersen, David Lee, Doug Jacobsen
!> \date   October 2017
!> \details
!>  This routine initializes the monotonic tracer advection tendencity using a FCT.
!
!-----------------------------------------------------------------------
   subroutine ocn_tracer_advection_mono_init(nHalos, horiz_adv_order, vert_adv_order, coef_3rd_order_in, dzdk_positive, & !{{{
                                             check_monotonicity, err)

      use mpas_dmpar
      integer, intent(in) :: nHalos !< Input: number of halos in current simulation
      integer, intent(in) :: horiz_adv_order !< Input: Order for horizontal advection
      integer, intent(in) :: vert_adv_order !< Input: Order for vertical advection
      real (kind=RKIND), intent(in) :: coef_3rd_order_in !< Input: coefficient for blending advection orders.
      logical, intent(in) :: dzdk_positive !< Input: Logical flag determining if dzdk is positive or negative.
      logical, intent(in) :: check_monotonicity !< Input: Logical flag determining check on monotonicity of tracers
      integer, intent(inout) :: err !< Input/Output: Error Flag

      err = 0

      vert2ndOrder = .false.
      vert3rdOrder = .false.
      vert4thOrder = .false.

      if ( horiz_adv_order == 3) then
          coef_3rd_order = coef_3rd_order_in
      else if(horiz_adv_order == 2 .or. horiz_adv_order == 4) then
          coef_3rd_order = 0.0_RKIND
      end if

      horizOrder = horiz_adv_order

      if (vert_adv_order == 3) then
          vert3rdOrder = .true.
      else if (vert_adv_order == 4) then
          vert4thOrder = .true.
      else
          vert2ndOrder = .true.
          if(vert_adv_order /= 2) then
             call mpas_log_write( &
                'Invalid value for vert_adv_order, defaulting to 2nd order', &
                MPAS_LOG_WARN)
          end if
      end if

      if (nHalos < 3) then
         call mpas_log_write( &
            'Monotonic advection cannot be used with less than 3 halos.', &
            MPAS_LOG_CRIT)
      end if

      positiveDzDk = dzdk_positive
      monotonicityCheck = check_monotonicity

   end subroutine ocn_tracer_advection_mono_init!}}}

end module ocn_tracer_advection_mono

