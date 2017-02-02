! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_gm

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_threading

   use ocn_constants

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

   public :: ocn_gm_compute_Bolus_velocity, &
             ocn_gm_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   private :: tridiagonal_solve

   ! Config options
   real (kind=RKIND), pointer :: config_gravWaveSpeed_trunc
   real (kind=RKIND), pointer :: config_standardGM_tracer_kappa
   real (kind=RKIND), pointer :: config_max_relative_slope
   real (kind=RKIND), pointer :: config_Redi_kappa
   logical, pointer :: config_disable_redi_k33
   logical, pointer :: config_use_Redi_surface_layer_tapering
   logical, pointer :: config_use_Redi_bottom_layer_tapering
   real (kind=RKIND), pointer :: config_Redi_surface_layer_tapering_extent
   real (kind=RKIND), pointer :: config_Redi_bottom_layer_tapering_depth

   real (kind=RKIND), parameter :: epsGM = 1.0e-12_RKIND

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_gm_compute_Bolus_velocity
!
!> \brief   Computes GM Bolus velocity
!> \author  Qingshan Chen, Mark Petersen, Todd Ringler
!> \date    January 2013
!> \details
!>  This routine is the main driver for the Gent-McWilliams (GM) parameterization.
!>  It computes horizontal and vertical density gradients, the slope
!>  of isopycnal surfaces, and solves a boundary value problem in each column
!>  for the stream function, which is used to compute the Bolus velocity.
!
!-----------------------------------------------------------------------

   subroutine ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input/Output: Scratch structure

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:,:), pointer :: density, displacedDensity, zMid, normalGMBolusVelocity, hEddyFlux, &
         layerThicknessEdge, gradDensityEdge, gradDensityTopOfEdge, gradDensityConstZTopOfEdge, gradZMidEdge, &
         gradZMidTopOfEdge, relativeSlopeTopOfEdge, relativeSlopeTopOfCell, k33, gmStreamFuncTopOfEdge, BruntVaisalaFreqTop, &
         gmStreamFuncTopOfCell, dDensityDzTopOfEdge, dDensityDzTopOfCell, relativeSlopeTapering, relativeSlopeTaperingCell, &
         areaCellSum
      real(kind=RKIND), dimension(:), pointer   :: boundaryLayerDepth
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: maxLevelEdgeTop, maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      integer                          :: i, k, iEdge, cell1, cell2, iCell, N, iter
      real(kind=RKIND)                 :: h1, h2, areaEdge, c, BruntVaisalaFreqTopEdge, rtmp, stmp, maxSlopeK33

      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      type (field2DReal), pointer :: gradDensityEdgeField, gradDensityTopOfEdgeField, gradDensityConstZTopOfEdgeField, &
         gradZMidEdgeField, gradZMidTopOfEdgeField, dDensityDzTopOfCellField, dDensityDzTopOfEdgeField,areaCellSumField

      call mpas_timer_start('gm bolus velocity')

      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfEdge', relativeSlopeTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfCell', relativeSlopeTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTapering', relativeSlopeTapering)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTaperingCell', relativeSlopeTaperingCell)
      call mpas_pool_get_array(diagnosticsPool, 'k33', k33)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'hEddyFlux', hEddyFlux)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', gmStreamFuncTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfCell', gmStreamFuncTopOfCell)

      if(config_use_Redi_surface_layer_tapering) call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',  maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell',  areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge',  dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge',  dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell',  nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell',  edgesOnCell)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_field(scratchPool, 'gradDensityEdge', gradDensityEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradDensityTopOfEdge', gradDensityTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradDensityConstZTopOfEdge', gradDensityConstZTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'dDensityDzTopOfCell', dDensityDzTopOfCellField)
      call mpas_pool_get_field(scratchPool, 'dDensityDzTopOfEdge', dDensityDzTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradZMidEdge', gradZMidEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradZMidTopOfEdge', gradZMidTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'areaCellSum', areaCellSumField)

      call mpas_allocate_scratch_field(gradDensityEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(gradDensityTopOfEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(gradDensityConstZTopOfEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(dDensityDzTopOfCellField, .True., .false.)
      call mpas_allocate_scratch_field(dDensityDzTopOfEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(gradZMidEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(gradZMidTopOfEdgeField, .True., .false.)
      call mpas_allocate_scratch_field(areaCellSumField, .True., .false.)

      call mpas_threading_barrier()

      gradDensityEdge => gradDensityEdgeField % array
      gradDensityTopOfEdge => gradDensityTopOfEdgeField % array
      gradDensityConstZTopOfEdge => gradDensityConstZTopOfEdgeField % array
      dDensityDzTopOfCell => dDensityDzTopOfCellField % array
      dDensityDzTopOfEdge => dDensityDzTopOfEdgeField % array
      gradZMidEdge => gradZMidEdgeField % array
      gradZMidTopOfEdge => gradZMidTopOfEdgeField % array
      areaCellSum => areaCellSumField % array

      allocate(rightHandSide(nVertLevels))
      allocate(tridiagA(nVertLevels))
      allocate(tridiagB(nVertLevels))
      allocate(tridiagC(nVertLevels))

      nCells = nCellsArray( size(nCellsArray) )
      nEdges = nEdgesArray( size(nEdgesArray) )

      ! Assign a huge value to the scratch variables which may manifest itself when
      ! there is a bug.
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges + 1
         do k = 1, nVertLevels
            gradDensityEdge(k, iEdge) = huge(0D0)
            gradDensityTopOfEdge(k, iEdge) = huge(0D0)
            dDensityDzTopOfEdge(k, iEdge) = huge(0D0)
            gradZMidEdge(k, iEdge) = huge(0D0)
            gradZMidTopOfEdge(k, iEdge) = huge(0D0)
            relativeSlopeTopOfEdge(k, iEdge) = 0.0_RKIND
            relativeSlopeTapering(k, iEdge) = 0.0_RKIND
            normalGMBolusVelocity(k, iEdge) = 0.0_RKIND
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells + 1
         do k = 1, nVertLevels
            dDensityDzTopOfCell(k,  iCell) = huge(0D0)
            k33(k, iCell) = 0.0_RKIND
            relativeSlopeTopOfCell(k, iCell) = 0.0_RKIND
            relativeSlopeTaperingCell(k, iCell) = 0.0_RKIND
         end do
      end do
      !$omp end do

      !--------------------------------------------------------------------
      !
      ! Compute vertical derivative of density at top of cell, interpolate to top of edge
      ! This is required for Redi and Bolus parts.
      !
      !--------------------------------------------------------------------

      nCells = nCellsArray( 3 )
      ! Compute vertical derivative of density (dDensityDzTopOfCell) at cell center and layer interface
      ! Note that displacedDensity is used from the upper cell, so that the EOS reference level for
      ! pressure is the same for both displacedDensity(k-1,iCell) and density(k,iCell).
      !$omp do schedule(runtime) private(k, rtmp)
      do iCell = 1, nCells
         do k = 2, maxLevelCell(iCell)
            rtmp = (displacedDensity(k-1,iCell) - density(k,iCell)) / (zMid(k-1,iCell) - zMid(k,iCell))
            dDensityDzTopOfCell(k,iCell) = min(rtmp, -epsGM)
         end do

         ! Approximation of dDensityDzTopOfCell on the top and bottom interfaces through the idea of having
         ! ghost cells above the top and below the bottom layers of the same depths and density.
         ! Essentially, this enforces the boundary condition (d density)/dz = 0 at the top and bottom.
         dDensityDzTopOfCell(1,iCell) = 0.0_RKIND
         dDensityDzTopOfCell(maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
      end do
      !$omp end do

      nEdges = nEdgesArray( 3 )

      ! Interpolate dDensityDzTopOfCell to edge and layer interface
      !$omp do schedule(runtime) private(k, cell1, cell2)
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)+1
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            dDensityDzTopOfEdge(k,iEdge) = 0.5_RKIND * (dDensityDzTopOfCell(k,cell1) + dDensityDzTopOfCell(k,cell2))
         end do
      end do
      !$omp end do

      !--------------------------------------------------------------------
      !
      ! Compute horizontal gradient and mid-layer of edge, interpolate to top of edge
      ! This is required for Redi and Bolus parts.
      !
      !--------------------------------------------------------------------

      nEdges = nEdgesArray( 3 )

      ! Compute density gradient (gradDensityEdge) and gradient of zMid (gradZMidEdge)
      ! along the constant coordinate surface.
      ! The computed variables lives at edge and mid-layer depth
      !$omp do schedule(runtime) private(cell1, cell2, k)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            gradDensityEdge(k,iEdge) = (density(k,cell2) - density(k,cell1)) / dcEdge(iEdge)
            gradZMidEdge(k,iEdge) = (zMid(k,cell2) - zMid(k,cell1)) / dcEdge(iEdge)
         end do
      end do
      !$omp end do

      nEdges = nEdgesArray( 3 )

      ! Interpolate gradDensityEdge and gradZMidEdge to layer interface
      !$omp do schedule(runtime) private(k, h1, h2)
      do iEdge = 1, nEdges
         ! The interpolation can only be carried out on non-boundary edges
         if (maxLevelEdgeTop(iEdge) .GE. 1) then
            do k = 2, maxLevelEdgeTop(iEdge)
               h1 = layerThicknessEdge(k-1,iEdge)
               h2 = layerThicknessEdge(k,iEdge)
               ! Using second-order interpolation below
               gradDensityTopOfEdge(k,iEdge) = (h2 * gradDensityEdge(k-1,iEdge) + h1 * gradDensityEdge(k,iEdge)) / (h1 + h2)
               gradZMidTopOfEdge(k,iEdge) = (h2 * gradZMidEdge(k-1,iEdge) + h1 * gradZMidEdge(k,iEdge)) / (h1 + h2)

            end do

            ! Approximation of values on the top and bottom interfaces through the idea of having ghost cells above
            ! the top and below the bottom layers of the same depths and density.
            gradDensityTopOfEdge(1,iEdge) = gradDensityEdge(1,iEdge)
            gradDensityTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = gradDensityEdge(maxLevelEdgeTop(iEdge),iEdge)
            gradZMidTopOfEdge(1,iEdge) = gradZMidEdge(1,iEdge)
            gradZMidTopOfEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = gradZMidEdge(maxLevelEdgeTop(iEdge),iEdge)
         end if
      end do
      !$omp end do

      !--------------------------------------------------------------------
      !
      ! Compute horizontal gradient required for Bolus part (along constant z)
      !
      !--------------------------------------------------------------------

      nEdges = nEdgesArray( 3 )

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         if (maxLevelEdgeTop(iEdge) .GE. 1) then
            do k = 1, maxLevelEdgeTop(iEdge)+1
               gradDensityConstZTopOfEdge(k,iEdge) = gradDensityTopOfEdge(k,iEdge) - dDensityDzTopOfEdge(k,iEdge) &
                                                   * gradZMidTopOfEdge(k,iEdge)
            end do
         end if
      end do
      !$omp end do

      !--------------------------------------------------------------------
      !
      ! Compute relative slope and k33 for Redi part of GM.
      ! These variables are used in del2 velocity tendency routines.
      !
      !--------------------------------------------------------------------

      nEdges = nEdgesArray( 3 )

      ! Compute relativeSlopeTopOfEdge at edge and layer interface
      ! set relativeSlopeTopOfEdge to zero for horizontal land/water edges.
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
        relativeSlopeTopOfEdge(:, iEdge) = 0.0_RKIND

         ! Beside a full land cell (e.g. missing cell) maxLevelEdgeTop=0, so relativeSlopeTopOfEdge at that edge will remain zero.
         do k = 2, maxLevelEdgeTop(iEdge)
            relativeSlopeTopOfEdge(k,iEdge) = - gradDensityTopOfEdge(k,iEdge) / min(dDensityDzTopOfEdge(k,iEdge),-epsGM)
         end do

         ! Since dDensityDzTopOfEdge is guaranteed to be zero on the top surface, relativeSlopeTopOfEdge on the top
         ! surface is identified with its value on the second interface.
         relativeSlopeTopOfEdge(1,iEdge) = relativeSlopeTopOfEdge(2,iEdge)

         ! dDensityDzTopOfEdge may or may not equal zero on the bottom surface, depending on whether
         ! maxLevelEdgeTop(iEdge) = maxLevelEdgeBottom(iEdge). But here we
         ! take a simplistic approach and identify relativeSlopeTopOfEdge on the bottom surface with its value on
         ! the interface just above.
         relativeSlopeTopOfEdge( maxLevelEdgeTop(iEdge)+1, iEdge ) = relativeSlopeTopOfEdge( max(1,maxLevelEdgeTop(iEdge)), iEdge )

      end do
      !$omp end do

      nEdges = nEdgesArray( 3 )

      ! slope can be unbounded in regions of neutral stability, reset to the large, but bounded, value
      ! values is hardwrite to 1.0, this is equivalent to a slope of 45 degrees
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            relativeSlopeTopOfEdge(k, iEdge) = max( min( relativeSlopeTopOfEdge(k, iEdge), 1.0_RKIND), -1.0_RKIND)
         end do
      end do
      !$omp end do

      ! average relative slope to cell centers
      ! do this by computing (relative slope)^2, then taking sqrt

      nCells = nCellsArray( 2 )

      !$omp do schedule(runtime) private(i, iEdge, areaEdge, rtmp, k)
      do iCell = 1, nCells
         areaCellSum(:, iCell) = 1.0e-34_RKIND
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)

            !contribution of cell area from this edge * 2.0
            areaEdge = 0.5_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
            do k = 1, maxLevelEdgeTop(iEdge)
               rtmp = areaEdge * relativeSlopeTopOfEdge(k, iEdge)**2
               relativeSlopeTopOfCell(k, iCell) = relativeSlopeTopOfCell(k, iCell) + rtmp
               areaCellSum(k, iCell) = areaCellSum(k, iCell) + areaEdge
            end do
         end do
      end do
      !$omp end do

      nCells = nCellsArray( 2 )

      !$omp do schedule(runtime) private(k)
      do iCell=1,nCells
        do k = 1, maxLevelCell(iCell)
           relativeSlopeTopOfCell(k,iCell) = sqrt( relativeSlopeTopOfCell(k,iCell)/areaCellSum(k,iCell) )
        end do
      end do
      !$omp end do

      ! Compute tapering function
      ! Compute k33 at cell center and layer interface

      nCells = nCellsArray( size(nCellsArray) )

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         k33(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      ! use relativeSlopeTaperingCell as a temporary space for smoothing of relativeSlopeTopOfCell
      relativeSlopeTaperingCell = relativeSlopeTopOfCell
      do iter = 1, 5
         !$omp barrier

         nCells = nCellsArray( 2 )

         !$omp do schedule(runtime) private(k, rtmp)
         do iCell=1,nCells
           relativeSlopeTaperingCell(1, iCell) = 0.0_RKIND
           relativeSlopeTaperingCell(maxLevelCell(iCell):nVertLevels, iCell) = 0.0_RKIND
           do k = 2, maxLevelCell(iCell)-1
             rtmp = relativeSlopeTopOfCell(k-1,iCell) + relativeSlopeTopOfCell(k+1,iCell)
             stmp = 2.0_RKIND*relativeSlopeTopOfCell(k,iCell)
             relativeSlopeTaperingCell(k,iCell) = (rtmp+stmp)/4.0_RKIND
           end do
           relativeSlopeTopOfCell(:, iCell) = relativeSlopeTaperingCell(:, iCell)
         end do
         !$omp end do
      end do  ! iter

      nCells = nCellsArray ( 2 )
      ! first, compute tapering across full domain based on a maximum allowable slope
      !$omp do schedule(runtime) private(k)
      do iCell=1,nCells
        do k = 1, maxLevelCell(iCell)
          relativeSlopeTaperingCell(k,iCell) = min(1.0_RKIND, config_max_relative_slope / (relativeSlopeTopOfCell(k,iCell)+epsGM))
        end do
      end do
      !$omp end do

      ! now further taper in the boundary layer
      ! vertical (k33) tapering starts at 2*OBL, increases linearly to OBL and is held uniform across OBL
      ! rtmp = 1 @ zMid = -2.0*OBL, rtmp = 0 @ zMid = -OBL
      if(config_use_Redi_surface_layer_tapering) then
         nCells = nCellsArray ( 2 )
         !$omp do schedule(runtime) private(k, rtmp)
         do iCell=1,nCells
           do k = 1, maxLevelCell(iCell)
             rtmp = -zMid(k,iCell)/max(config_Redi_surface_layer_tapering_extent,boundaryLayerDepth(iCell)+epsGM)
             rtmp = max(0.0_RKIND,rtmp)
             rtmp = min(1.0_RKIND,rtmp)
             relativeSlopeTaperingCell(k,iCell) = rtmp*relativeSlopeTaperingCell(k,iCell)
           end do
         end do
         !$omp end do
      endif ! config_use_Redi_surface_layer_tapering

      ! now further taper in the boundary layer
      ! vertical (k33) tapering starts at 2*OBL, increases linearly to OBL and is held uniform across OBL
      ! rtmp = 1 @ zMid = zMid(maxLevelCell) + config_Redi_bottom_layer_tapering_depth, rtmp = 0 @ zMid = zMid(maxLevelCell)
      if(config_use_Redi_bottom_layer_tapering) then
         nCells = nCellsArray ( 2 )
         !$omp do schedule(runtime) private(k, rtmp)
         do iCell=1,nCells
           do k = 1, maxLevelCell(iCell)
             rtmp = (zMid(k,iCell)-zMid(maxLevelCell(iCell),iCell))/(config_Redi_bottom_layer_tapering_depth+epsGM)
             rtmp = max(0.0_RKIND,rtmp)
             rtmp = min(1.0_RKIND,rtmp)
             relativeSlopeTaperingCell(k,iCell) = rtmp*relativeSlopeTaperingCell(k,iCell)
           end do
         end do
         !$omp end do
      endif ! config_use_Redi_bottom_layer_tapering

      nCells = nCellsArray( 2 )
      !$omp do schedule(runtime) private(k)
      do iCell=1,nCells
        k33(:, iCell) = 0.0_RKIND
        do k = 2, maxLevelCell(iCell)
          k33(k,iCell) = ( relativeSlopeTaperingCell(k,iCell) * relativeSlopeTopOfCell(k,iCell) )**2
        end do
      end do
      !$omp end do

      nEdges = nEdgesArray( 3 )

      ! average tapering function to layer edges
      !$omp do schedule(runtime) private(cell1, cell2, k)
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        do k = 1, maxLevelEdgeTop(iEdge)
          relativeSlopeTapering(k,iEdge) = 0.5_RKIND * (relativeSlopeTaperingCell(k,cell1) + relativeSlopeTaperingCell(k,cell2))
        enddo
      enddo
      !$omp end do

      ! allow disabling of K33 for testing
      if(config_disable_redi_k33) then
        nCells = nCellsArray( size(nCellsArray) )
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           k33(:, iCell) = 0.0_RKIND
        end do
        !$omp end do
      end if

      !--------------------------------------------------------------------
      !
      ! Compute stream function and Bolus velocity for Bolus part of GM
      !
      !--------------------------------------------------------------------

      c = config_gravWaveSpeed_trunc**2

      nEdges = nEdgesArray( 3 )

      !$omp do schedule(runtime) private(cell1, cell2, k, BruntVaisalaFreqTopEdge, N)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         gmStreamFuncTopOfEdge(:, iEdge) = 0.0_RKIND

         ! Construct the tridiagonal matrix
         if (maxLevelEdgeTop(iEdge) .GE. 3) then
            ! First row
            k = 2
            BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
            BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
            tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2 / (layerThicknessEdge(k-1,iEdge) &
                          * layerThicknessEdge(k,iEdge)) - BruntVaisalaFreqTopEdge
            tridiagC(k-1) = 2.0_RKIND * config_gravWaveSpeed_trunc**2 / layerThicknessEdge(k, iEdge) &
                          / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
            rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)

            ! Second to next to the last rows
            do k = 3, maxLevelEdgeTop(iEdge)-1
               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagA(k-2) = 2.0_RKIND * config_gravWaveSpeed_trunc**2 / layerThicknessEdge(k-1, iEdge) &
                             / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
               tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2 / (layerThicknessEdge(k-1, iEdge) &
                             * layerThicknessEdge(k, iEdge) ) - BruntVaisalaFreqTopEdge
               tridiagC(k-1) = 2.0_RKIND * config_gravWaveSpeed_trunc**2 / layerThicknessEdge(k, iEdge) &
                             / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
               rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)
            end do

            ! Last row
            k = maxLevelEdgeTop(iEdge)
            BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
            BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
            tridiagA(k-2) = 2.0_RKIND * config_gravWaveSpeed_trunc**2 / layerThicknessEdge(k-1,iEdge) &
                          / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge))
            tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2 / (layerThicknessEdge(k-1, iEdge) &
                          * layerThicknessEdge(k, iEdge)) - BruntVaisalaFreqTopEdge
            rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)

            ! Total number of rows
            N = maxLevelEdgeTop(iEdge) - 1

            ! Call the tridiagonal solver
            call tridiagonal_solve(tridiagA, tridiagB, tridiagC, rightHandSide, &
                                   gmStreamFuncTopOfEdge(2:maxLevelEdgeTop(iEdge), iEdge), N)
         end if
      end do
      !$omp end do

      nEdges = nEdgesArray( 3 )
      ! Compute normalGMBolusVelocity from the stream function
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)
            normalGMBolusVelocity(k,iEdge) = (gmStreamFuncTopOfEdge(k,iEdge) - gmStreamFuncTopOfEdge(k+1,iEdge)) &
                                           / layerThicknessEdge(k,iEdge)
         end do
      end do
      !$omp end do

      nCells = nCellsArray( 1 )

      ! Interpolate gmStreamFuncTopOfEdge to cell centers for visualization
      !$omp do schedule(runtime) private(i, iEdge, areaEdge, k, rtmp)
      do iCell = 1, nCells
         gmStreamFuncTopOfCell(:, iCell) = 0.0_RKIND
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)

            areaEdge = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)
               rtmp = 0.5_RKIND * ( gmStreamFuncTopOfEdge(k, iEdge) + gmStreamFuncTopOfEdge(k+1, iEdge) ) * areaEdge
               gmStreamFuncTopOfCell(k, iCell) = gmStreamFuncTopOfCell(k, iCell) + rtmp
            end do
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         gmStreamFuncTopOfCell(:, iCell) = gmStreamFuncTopOfCell(:,iCell) / areaCell(iCell)
      end do
      !$omp end do

      deallocate(rightHandSide)
      deallocate(tridiagA)
      deallocate(tridiagB)
      deallocate(tridiagC)

      call mpas_threading_barrier()

      ! Deallocate scratch variables
      call mpas_deallocate_scratch_field(gradDensityEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradDensityTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradDensityConstZTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(dDensityDzTopOfCellField, .true.)
      call mpas_deallocate_scratch_field(dDensityDzTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradZMidEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradZMidTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(areaCellSumField, .true.)

      call mpas_timer_stop('gm bolus velocity')

   end subroutine ocn_gm_compute_Bolus_velocity!}}}

!***********************************************************************
!
!  routine tridiagonal_solve
!
!> \brief   Solve the matrix equation Ax=r for x, where A is tridiagonal.
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
! mrp note:  This subroutine also appears in vmix and should really be put in the framework.
   subroutine tridiagonal_solve(a,b,c,r,x,n) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer,intent(in) :: n
      real (KIND=RKIND), dimension(n), intent(in) :: a,b,c,r

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n), intent(out) :: x

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (KIND=RKIND), dimension(n) :: bTemp,rTemp
      real (KIND=RKIND) :: m
      integer i

      ! Use work variables for b and r
      bTemp(1) = b(1)
      rTemp(1) = r(1)

      ! First pass: set the coefficients
      do i = 2,n
         m = a(i-1)/bTemp(i-1)
         bTemp(i) = b(i) - m*c(i-1)
         rTemp(i) = r(i) - m*rTemp(i-1)
      end do

      x(n) = rTemp(n)/bTemp(n)
       ! Second pass: back-substition
      do i = n-1, 1, -1
         x(i) = (rTemp(i) - c(i)*x(i+1))/bTemp(i)
      end do

   end subroutine tridiagonal_solve !}}}

!***********************************************************************
!
!  routine ocn_gm_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_gm_init(err)!{{{

      !-----------------------------------------------------------------
      !
      ! Output Variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_gravWaveSpeed_trunc',config_gravWaveSpeed_trunc)
      call mpas_pool_get_config(ocnConfigs, 'config_standardGM_tracer_kappa',config_standardGM_tracer_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_max_relative_slope',config_max_relative_slope)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_kappa', config_Redi_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_k33',config_disable_redi_k33)
      call mpas_pool_get_config(ocnConfigs, 'config_use_Redi_surface_layer_tapering',config_use_Redi_surface_layer_tapering)
      call mpas_pool_get_config(ocnConfigs, 'config_use_Redi_bottom_layer_tapering',config_use_Redi_bottom_layer_tapering)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_surface_layer_tapering_extent',config_Redi_surface_layer_tapering_extent)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_bottom_layer_tapering_depth',config_Redi_bottom_layer_tapering_depth)

   end subroutine ocn_gm_init!}}}

!***********************************************************************

end module ocn_gm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
