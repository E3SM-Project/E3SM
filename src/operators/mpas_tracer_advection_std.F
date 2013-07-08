! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tracer_advection_std
!
!> \brief MPAS standard tracer advection
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains routines for standard advection of tracers
!
!-----------------------------------------------------------------------
module mpas_tracer_advection_std

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_tracer_advection_helpers

   implicit none
   private
   save 

   real (kind=RKIND) :: coef_3rd_order
   logical :: vert2ndOrder, vert3rdOrder, vert4thOrder

   public :: mpas_tracer_advection_std_tend, &
             mpas_tracer_advection_std_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_tend
!
!> \brief MPAS standard tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine computes the standard tracer advection tendencity.
!>  Both horizontal and vertical.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_tend(tracers, normalThicknessFlux, w, layerThickness, verticalCellSize, dt, grid, tend_layerThickness, tend, maxLevelCell_in, maxLevelEdgeTop_in, highOrderAdvectionMask_in, verticalDivergenceFactor_in)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_layerThickness !< Input: Tendency for thickness field
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      type (mesh_type), intent(in) :: grid !< Input: Grid information
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency
      integer, dimension(:), pointer, optional :: maxLevelCell_in !< Input - optional: Index to max level at cell center
      integer, dimension(:), pointer, optional :: maxLevelEdgeTop_in !< Input - optional: Index to max level at edge with non-land cells on both sides
      integer, dimension(:,:), pointer, optional :: highOrderAdvectionMask_in !< Input - optional: Mask for high order advection
      real (kind=RKIND), dimension(:), pointer, optional :: verticalDivergenceFactor_in !< Input - optional: Denominator for vertical divergence. Given as an inverse which can be multiplied.

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2
      integer :: nVertLevels, num_tracers, nCells, nEdges, nCellsSolve
      integer, dimension(:), pointer :: nAdvCellsForEdge, maxLevelCell, maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, advCellsForEdge, highOrderAdvectionMask, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: flux, tracer_weight, invDvEdge, invAreaCell1, invAreaCell2
      real (kind=RKIND) :: verticalWeightK, verticalWeightKm1
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, verticalDivergenceFactor
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension(:,:), pointer :: tracer_cur, high_order_horiz_flux, high_order_vert_flux
      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers,dim=1)

      ! Initialize pointers
      dvEdge      => grid % dvEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array
      edgeSignOnCell => grid % edgeSignOnCell % array
      cellsOnCell => grid % cellsOnCell % array
      areaCell    => grid % areaCell % array

      nEdgesOnCell => grid % nEdgesOnCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      
      ! Setup arrays from optional arguments
      if(present(maxLevelCell_in)) then
        maxLevelCell => maxLevelCell_in
      else
        allocate(maxLevelCell(nCells+1))
        maxLevelCell(:) = nVertLevels
      end if

      if(present(maxLevelEdgeTop_in)) then
        maxLevelEdgeTop => maxLevelEdgeTop_in
      else
        allocate(maxLevelEdgeTop(nEdges+1))
        maxLevelEdgeTop(:) = nVertLevels
      end if

      if(present(highOrderAdvectionMask_in)) then
        highOrderAdvectionMask => highOrderAdvectionMask_in
      else
        allocate(highOrderAdvectionMask(nVertLevels, nEdges+1))
        if(config_horiz_tracer_adv_order == 2) then
          highOrderAdvectionMask = 0
        else
          highOrderAdvectionMask = 1
        end if
      end if

      if(present(verticalDivergenceFactor_in)) then
        verticalDivergenceFactor => verticalDivergenceFactor_in
      else
        allocate(verticalDivergenceFactor(nVertLevels))
        verticalDivergenceFactor = 1.0_RKIND
      end if

      ! Setup high order horizontal flux
      allocate(high_order_horiz_flux(nVertLevels, nEdges+1))

      ! allocate nCells arrays
      allocate(tracer_cur(nVertLevels, nCells+1))

      ! allocate nVertLevels+1 and nCells arrays
      allocate(high_order_vert_flux(nVertLevels+1, nCells))

      ! Loop over tracers. One tracer is advected at a time. It is copied into a temporary array in order to improve locality
      do iTracer = 1, num_tracers
        ! Initialize variables for use in this iTracer iteration
        tracer_cur(:, :) = tracers(iTracer, :, :)

        high_order_vert_flux = 0.0_RKIND
        high_order_horiz_flux = 0.0_RKIND

        !  Compute the high order vertical flux. Also determine bounds on tracer_cur. 
        do iCell = 1, nCells
          k = max(1, min(maxLevelCell(iCell), 2))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
             
          do k=3,maxLevelCell(iCell)-1
             if(vert4thOrder) then
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux4( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell))
             else if(vert3rdOrder) then
               high_order_vert_flux(k, iCell) = mpas_tracer_advection_vflux3( tracer_cur(k-2,iCell),tracer_cur(k-1,iCell),  &
                                      tracer_cur(k  ,iCell),tracer_cur(k+1,iCell), w(k,iCell), coef_3rd_order )
             else if (vert2ndOrder) then
               verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
               high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
             end if
          end do 
 
          k = max(1, maxLevelCell(iCell))
          verticalWeightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          verticalWeightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k, iCell) + verticalCellSize(k-1, iCell))
          high_order_vert_flux(k,iCell) = w(k,iCell)*(verticalWeightK*tracer_cur(k,iCell)+verticalWeightKm1*tracer_cur(k-1,iCell))
        end do ! iCell Loop

        !  Compute the high order horizontal flux
        do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Compute 2nd order fluxes where needed.
          do k = 1, maxLevelEdgeTop(iEdge)
            tracer_weight = iand(highOrderAdvectionMask(k, iEdge)+1, 1) * (dvEdge(iEdge) * 0.5_RKIND) * normalThicknessFlux(k, iEdge)

            high_order_horiz_flux(k, iEdge) = high_order_horiz_flux(k, iedge) + tracer_weight * (tracer_cur(k, cell1) + tracer_cur(k, cell2))
          end do ! k loop

          ! Compute 3rd or 4th fluxes where requested.
          do i = 1, nAdvCellsForEdge(iEdge)
            iCell = advCellsForEdge(i,iEdge)
            do k = 1, maxLevelCell(iCell)
              tracer_weight = highOrderAdvectionMask(k, iEdge) * (adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,normalThicknessFlux(k,iEdge))*adv_coefs_3rd(i,iEdge))

              tracer_weight = normalThicknessFlux(k,iEdge)*tracer_weight
              high_order_horiz_flux(k,iEdge) = high_order_horiz_flux(k,iEdge) + tracer_weight * tracer_cur(k,iCell)
            end do ! k loop
          end do ! i loop over nAdvCellsForEdge
        end do ! iEdge loop

        ! Accumulate the scaled high order horizontal tendencies
        do iCell = 1, nCells
          invAreaCell1 = 1.0 / areaCell(iCell)
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + edgeSignOnCell(i, iCell) * high_order_horiz_flux(k, iEdge) * invAreaCell1
            end do
          end do
        end do

        ! Accumulate the scaled high order vertical tendencies.
        do iCell = 1, nCellsSolve
          do k = 1,maxLevelCell(iCell)
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + verticalDivergenceFactor(k) * (high_order_vert_flux(k+1, iCell) - high_order_vert_flux(k, iCell))
          end do ! k loop
        end do ! iCell loop
      end do ! iTracer loop

      deallocate(tracer_cur)
      deallocate(high_order_horiz_flux)
      deallocate(high_order_vert_flux)

      if(.not.present(maxLevelCell_in)) then
        deallocate(maxLevelCell)
      end if

      if(.not.present(maxLevelEdgeTop_in)) then
        deallocate(maxLevelEdgeTop)
      end if

      if(.not.present(highOrderAdvectionMask_in)) then
        deallocate(highOrderAdvectionMask)
      end if

      if(.not.present(verticalDivergenceFactor_in)) then
        deallocate(verticalDivergenceFactor)
      end if
   end subroutine mpas_tracer_advection_std_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_init
!
!> \brief MPAS initialize standard tracer advection tendency.
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the standard tracer advection tendencity.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_init(err)!{{{
      integer, intent(inout) :: err !< Input: Error Flags

      integer :: err_tmp

      err = 0

      vert2ndOrder = .false.
      vert3rdOrder = .false.
      vert4thOrder = .false.

      if ( config_horiz_tracer_adv_order == 3) then
          coef_3rd_order = config_coef_3rd_order
      else if(config_horiz_tracer_adv_order == 2 .or. config_horiz_tracer_adv_order == 4) then
          coef_3rd_order = 0.0_RKIND
      end if

      if (config_vert_tracer_adv_order == 3) then
          vert3rdOrder = .true.
      else if (config_vert_tracer_adv_order == 4) then
          vert4thOrder = .true.
      else
          vert2ndOrder = .true.
          if(config_vert_tracer_adv_order /= 2) then
            write(6,*) 'Invalid value for config_vert_tracer_adv_order, defaulting to 2nd order'
          end if
      end if

   end subroutine mpas_tracer_advection_std_init!}}}

end module mpas_tracer_advection_std
