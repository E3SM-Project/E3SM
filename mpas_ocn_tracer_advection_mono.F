! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
!> \brief MPAS ocean monotonic tracer advection with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for monotonic advection of tracers using a FCT
!
!-----------------------------------------------------------------------
module mpas_ocn_tracer_advection_mono

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_ocn_tracer_advection_helpers

   implicit none
   private
   save 

   real (kind=RKIND) :: coef_3rd_order

   public :: mpas_ocn_tracer_advection_mono_tend, &
             mpas_ocn_tracer_advection_mono_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_mono_tend
!
!> \brief MPAS ocean monotonic tracer advection tendency with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the monotonic tracer advection tendencity using a FCT.
!>  Both horizontal and vertical.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_mono_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, grid, tend_layerThickness, tend)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Tendency for thickness field
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mesh_type), intent(in) :: grid !< Input: Grid information
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2
      integer :: nVertLevels, num_tracers, nCells, nEdges, nCellsSolve
      integer, dimension(:), pointer :: nAdvCellsForEdge, maxLevelCell, maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, advCellsForEdge, highOrderAdvectionMask, lowOrderAdvectionMask, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: flux_upwind, tracer_min_new, tracer_max_new, tracer_upwind_new, scale_factor
      real (kind=RKIND) :: flux, tracer_weight, invDvEdge, invAreaCell1, invAreaCell2
      real (kind=RKIND) :: cur_max, cur_min, new_max, new_min
      real (kind=RKIND) :: verticalWeightK, verticalWeightKm1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_2nd, adv_coefs_3rd
      real (kind=RKIND), dimension(:,:), pointer :: tracer_cur, tracer_new, upwind_tendency, inv_h_new, tracer_max, tracer_min
      real (kind=RKIND), dimension(:,:), pointer :: flux_incoming, flux_outgoing, high_order_horiz_flux, high_order_vert_flux

      real (kind=RKIND), parameter :: eps = 1.e-10

      type (field2dReal), pointer :: high_order_horiz_flux_field

      ! Initialize pointers
      dvEdge      => grid % dvEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      cellsOnCell => grid % cellsOnCell % array
      areaCell    => grid % areaCell % array

      nEdgesOnCell => grid % nEdgesOnCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_2nd => grid % adv_coefs_2nd % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      maxLevelCell => grid % maxLevelCell % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      highOrderAdvectionMask => grid % highOrderAdvectionMask % array
      lowOrderAdvectionMask => grid % lowOrderAdvectionMask % array
      edgeSignOnCell => grid % edgeSignOnCell % array
      edgesOnCell => grid % edgesOnCell % array

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers,dim=1)

      allocate(high_order_horiz_flux_field)
      nullify(high_order_horiz_flux_field % next)
      high_order_horiz_flux_field % block => grid % block
      high_order_horiz_flux_field % sendList => grid % xEdge % sendList
      high_order_horiz_flux_field % recvList => grid % xEdge % recvList
      high_order_horiz_flux_field % copyList => grid % xEdge % copyList
      high_order_horiz_flux_field % dimSizes(1) = nVertLevels
      high_order_horiz_flux_field % dimSizes(2) = nEdges+1
      allocate(high_order_horiz_flux_field % array(high_order_horiz_flux_field % dimSizes(1), high_order_horiz_flux_field % dimSizes(2)))
      high_order_horiz_flux => high_order_horiz_flux_field % array

      ! allocate nCells arrays

      allocate(tracer_new(nVertLevels, nCells+1))
      allocate(tracer_cur(nVertLevels, nCells+1))
      allocate(upwind_tendency(nVertLevels, nCells+1))
      allocate(inv_h_new(nVertLevels, nCells+1))
      allocate(tracer_max(nVertLevels, nCells+1))
      allocate(tracer_min(nVertLevels, nCells+1))
      allocate(flux_incoming(nVertLevels, nCells+1))
      allocate(flux_outgoing(nVertLevels, nCells+1))

      ! allocate nEdges arrays
!     allocate(high_order_horiz_flux(nVertLevels, nEdges))

      ! allocate nVertLevels+1 and nCells arrays
      allocate(high_order_vert_flux(nVertLevels+1, nCells))

      do iCell = 1, nCells
        do k=1, maxLevelCell(iCell)
          inv_h_new(k, iCell) = 1.0 / (layerThickness(k, iCell) + dt * tend_layerThickness(k, iCell))
        end do
      end do

      ! Loop over tracers. One tracer is advected at a time. It is copied into a temporary array in order to improve locality
      do iTracer = 1, num_tracers
        ! Initialize variables for use in this iTracer iteration
        do iCell = 1, nCells
          do k=1, maxLevelCell(iCell)
            tracer_cur(k,iCell) = tracers(iTracer,k,iCell)
            upwind_tendency(k, iCell) = 0.0

            !tracer_new is supposed to be the "new" tracer state. This allows bounds checks.
            if (config_check_tracer_monotonicity) then
              tracer_new(k,iCell) = 0.0
            end if
          end do ! k loop
        end do ! iCell loop

        high_order_vert_flux = 0.0
        high_order_horiz_flux = 0.0

        !  Compute the high order vertical flux. Also determine bounds on tracer_cur. 
        do iCell = 1, nCells
          k = 1
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k+1,iCell))

          k = 2
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
          tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             
          do k=3,maxLevelCell(iCell)-1
             high_order_vert_flux(k,iCell) = mpas_ocn_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                    tracer_cur(k  ,iCell),tracer_cur(k+1,iCell),  &
                                    w(k,iCell), coef_3rd_order )
             tracer_max(k,iCell) = max(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
             tracer_min(k,iCell) = min(tracer_cur(k-1,iCell),tracer_cur(k,iCell),tracer_cur(k+1,iCell))
          end do 
 
          k = maxLevelCell(iCell)
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
          tracer_max(k,iCell) = max(tracer_cur(k,iCell),tracer_cur(k-1,iCell))
          tracer_min(k,iCell) = min(tracer_cur(k,iCell),tracer_cur(k-1,iCell))

          ! pull tracer_min and tracer_max from the (horizontal) surrounding cells
          do i = 1, nEdgesOnCell(iCell)
            do k=1, maxLevelCell(cellsOnCell(i, iCell))
              tracer_max(k,iCell) = max(tracer_max(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
              tracer_min(k,iCell) = min(tracer_min(k,iCell),tracer_cur(k, cellsOnCell(i,iCell)))
            end do ! k loop
          end do ! i loop over nEdgesOnCell
        end do ! iCell Loop

        !  Compute the high order horizontal flux
        do iEdge = 1, nEdges
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 1, maxLevelCell(iCell)
              tracer_weight = lowOrderAdvectionMask(k, iEdge) * adv_coefs_2nd(i,iEdge) & 
                            + highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))

              tracer_weight = normalThicknessFlux(k,iEdge)*tracer_weight
              high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) + tracer_weight* tracer_cur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        end do ! iEdge loop

        !  low order upwind vertical flux (monotonic and diffused)
        !  Remove low order flux from the high order flux.
        !  Store left over high order flux in high_order_vert_flux array.
        !  Upwind fluxes are accumulated in upwind_tendency
        do iCell = 1, nCells
          do k = 2, maxLevelCell(iCell)
            ! dwj 02/03/12 Ocean and Atmosphere are different in vertical
!           flux_upwind = max(0.0_RKIND,w(k,iCell))*tracer_cur(k-1,iCell) + min(0.0_RKIND,w(k,iCell))*tracer_cur(k,iCell)
            flux_upwind = min(0.0_RKIND,w(k,iCell))*tracer_cur(k-1,iCell) + max(0.0_RKIND,w(k,iCell))*tracer_cur(k,iCell)
            upwind_tendency(k-1,iCell) = upwind_tendency(k-1,iCell) + flux_upwind
            upwind_tendency(k  ,iCell) = upwind_tendency(k  ,iCell) - flux_upwind
            high_order_vert_flux(k,iCell) = high_order_vert_flux(k,iCell) - flux_upwind
          end do ! k loop

          ! flux_incoming contains the total remaining high order flux into iCell
          !          it is positive.
          ! flux_outgoing contains the total remaining high order flux out of iCell
          !           it is negative
          do k = 1, maxLevelCell(iCell)
            ! dwj 02/03/12 Ocean and Atmosphere are different in vertical
!           flux_incoming (k,iCell) = -(min(0.0_RKIND,high_order_vert_flux(k+1,iCell))-max(0.0_RKIND,high_order_vert_flux(k,iCell)))
!           flux_outgoing(k,iCell) = -(max(0.0_RKIND,high_order_vert_flux(k+1,iCell))-min(0.0_RKIND,high_order_vert_flux(k,iCell)))

            flux_incoming (k, iCell) = max(0.0_RKIND, high_order_vert_flux(k+1, iCell)) - min(0.0_RKIND, high_order_vert_flux(k, iCell))
            flux_outgoing(k, iCell) = min(0.0_RKIND, high_order_vert_flux(k+1, iCell)) - max(0.0_RKIND, high_order_vert_flux(k, iCell))
          end do ! k Loop
        end do ! iCell Loop

        !  low order upwind horizontal flux (monotinc and diffused)
        !  Remove low order flux from the high order flux
        !  Store left over high order flux in high_order_horiz_flux array
        !  Upwind fluxes are accumulated in upwind_tendency
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          invAreaCell1 = 1.0 / areaCell(cell1)
          invAreaCell2 = 1.0 / areaCell(cell2)

          do k = 1, maxLevelEdgeTop(iEdge)
            flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))
            high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) - flux_upwind
          end do ! k loop
        end do ! iEdge loop

        do iCell = 1, nCells
          invAreaCell1 = 1.0 / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1, maxLevelEdgeTop(iEdge)
              flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell1) + min(0.0_RKIND,normalThicknessFlux(k,iEdge))*tracer_cur(k,cell2))

              upwind_tendency(k,iCell) = upwind_tendency(k,iCell) + edgeSignOncell(i, iCell) * flux_upwind * invAreaCell1

              ! Accumulate remaining high order fluxes
              flux_outgoing(k,iCell) = flux_outgoing(k,iCell) + min(0.0_RKIND, edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge)) * invAreaCell1
              flux_incoming(k,iCell) = flux_incoming(k,iCell) + max(0.0_RKIND, edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge)) * invAreaCell1
            end do
          end do
        end do

        ! Build the factors for the FCT
        ! Computed using the bounds that were computed previously, and the bounds on the newly updated value
        ! Factors are placed in the flux_incoming and flux_outgoing arrays
        do iCell = 1, nCells
          do k = 1, maxLevelCell(iCell)
            tracer_min_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(upwind_tendency(k,iCell)+flux_outgoing(k,iCell))) * inv_h_new(k,iCell)
            tracer_max_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*(upwind_tendency(k,iCell)+flux_incoming(k,iCell))) * inv_h_new(k,iCell)
            tracer_upwind_new = (tracer_cur(k,iCell)*layerThickness(k,iCell) + dt*upwind_tendency(k,iCell)) * inv_h_new(k,iCell)
           
            scale_factor = (tracer_max(k,iCell)-tracer_upwind_new)/(tracer_max_new-tracer_upwind_new+eps)
            flux_incoming(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

            scale_factor = (tracer_upwind_new-tracer_min(k,iCell))/(tracer_upwind_new-tracer_min_new+eps)
            flux_outgoing(k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          end do ! k loop
        end do ! iCell loop

        !  rescale the high order horizontal fluxes
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k = 1, maxLevelEdgeTop(iEdge)
            flux = high_order_horiz_flux(k,iEdge)
            flux = max(0.0_RKIND,flux) * min(flux_outgoing(k,cell1), flux_incoming(k,cell2)) &
                 + min(0.0_RKIND,flux) * min(flux_incoming(k,cell1), flux_outgoing(k,cell2))
            high_order_horiz_flux(k,iEdge) = flux
          end do ! k loop
        end do ! iEdge loop

        ! rescale the high order vertical flux
        do iCell = 1, nCellsSolve
          do k = 2, maxLevelCell(iCell)
            flux =  high_order_vert_flux(k,iCell)
            ! dwj 02/03/12 Ocean and Atmosphere are different in vertical.
!           flux = max(0.0_RKIND,flux) * min(flux_outgoing(k-1,iCell), flux_incoming(k  ,iCell)) &
!                + min(0.0_RKIND,flux) * min(flux_outgoing(k  ,iCell), flux_incoming(k-1,iCell))
            flux = max(0.0_RKIND,flux) * min(flux_outgoing(k  ,iCell), flux_incoming(k-1,iCell)) &
                 + min(0.0_RKIND,flux) * min(flux_outgoing(k-1,iCell), flux_incoming(k  ,iCell))
            high_order_vert_flux(k,iCell) = flux
          end do ! k loop
        end do ! iCell loop

        ! Accumulate the scaled high order horizontal tendencies
        do iCell = 1, nCells
          invAreaCell1 = 1.0 / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) * invAreaCell1

              if(config_check_tracer_monotonicity) then
                tracer_new(k, iCell) = tracer_new(k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) * invAreaCell1
              end if
            end do
          end do
        end do

        ! Accumulate the scaled high order vertical tendencies, and the upwind tendencies
        do iCell = 1, nCellsSolve
          do k = 1,maxLevelCell(iCell)
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + (high_order_vert_flux(k+1, iCell) - high_order_vert_flux(k, iCell)) + upwind_tendency(k,iCell)

            if (config_check_tracer_monotonicity) then
              !tracer_new holds a tendency for now. Only for a check on monotonicity
              tracer_new(k, iCell) = tracer_new(k, iCell) + (high_order_vert_flux(k+1, iCell) - high_order_vert_flux(k, iCell)) + upwind_tendency(k,iCell)

              !tracer_new is now the new state of the tracer. Only for a check on monotonicity
              tracer_new(k, iCell) = (tracer_cur(k, iCell)*layerThickness(k, iCell) + dt * tracer_new(k, iCell)) * inv_h_new(k, iCell)
            end if
          end do ! k loop
        end do ! iCell loop

        if (config_check_tracer_monotonicity) then
          !build min and max bounds on old and new tracer for check on monotonicity.
          do iCell = 1, nCellsSolve
            do k = 1, maxLevelCell(iCell)
              if(tracer_new(k,iCell) < tracer_min(k, iCell)-eps) then
                write(*,*) 'Minimum out of bounds on tracer ', iTracer, tracer_min(k, iCell), tracer_new(k,iCell)
              end if

              if(tracer_new(k,iCell) > tracer_max(k,iCell)+eps) then
                write(*,*) 'Maximum out of bounds on tracer ', iTracer, tracer_max(k, iCell), tracer_new(k,iCell)
              end if
            end do
          end do
        end if
      end do ! iTracer loop

      deallocate(tracer_new)
      deallocate(tracer_cur)
      deallocate(upwind_tendency)
      deallocate(inv_h_new)
      deallocate(tracer_max)
      deallocate(tracer_min)
      deallocate(flux_incoming)
      deallocate(flux_outgoing)
      deallocate(high_order_horiz_flux)
      deallocate(high_order_vert_flux)
      deallocate(high_order_horiz_flux_field)
   end subroutine mpas_ocn_tracer_advection_mono_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_mono_init
!
!> \brief MPAS ocean initialize monotonic tracer advection tendency with FCT
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the monotonic tracer advection tendencity using a FCT.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_mono_init(err)!{{{
      integer, intent(inout) :: err !< Input: Error Flags

      integer :: err_tmp

      err = 0

      if ( config_horiz_tracer_adv_order == 3) then
          coef_3rd_order = config_coef_3rd_order
      else if(config_horiz_tracer_adv_order == 2 .or. config_horiz_tracer_adv_order == 4) then
          coef_3rd_order = 0.0
      end if

   end subroutine mpas_ocn_tracer_advection_mono_init!}}}

end module mpas_ocn_tracer_advection_mono
