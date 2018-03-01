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
   real (kind=RKIND), pointer :: config_gravWaveSpeed_trunc, config_standardGM_tracer_kappa, config_density0, &
     config_max_relative_slope, config_Redi_kappa
   logical, pointer :: config_use_standardGM
   logical, pointer :: config_disable_redi_k33

   real (kind=RKIND), parameter :: epsGM = 1.0e-12

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

      real(kind=RKIND), dimension(:,:), pointer :: density, displacedDensity, zMid, normalGMBolusVelocity, hEddyFlux, layerThicknessEdge,  &
         gradDensityEdge, gradDensityTopOfEdge, gradDensityConstZTopOfEdge, gradZMidEdge, &
         gradZMidTopOfEdge, relativeSlopeTopOfEdge, relativeSlopeTopOfCell, k33, gmStreamFuncTopOfEdge, BruntVaisalaFreqTop, gmStreamFuncTopOfCell, &
         dDensityDzTopOfEdge, dDensityDzTopOfCell, relativeSlopeTapering, relativeSlopeTaperingCell, areaCellSum
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer                          :: k, iEdge, cell1, cell2, iCell, N
      real(kind=RKIND)                 :: h1, h2, areaEdge, c, BruntVaisalaFreqTopEdge, rtmp, maxSlopeK33

      ! Dimensions
      integer, pointer :: nCells, nEdges

      type (field2DReal), pointer :: gradDensityEdgeField, gradDensityTopOfEdgeField, gradDensityConstZTopOfEdgeField, &
         gradZMidEdgeField, gradZMidTopOfEdgeField, dDensityDzTopOfCellField, dDensityDzTopOfEdgeField,areaCellSumField

      type (field1DReal), pointer :: rightHandSideField, tridiagAField, tridiagBField, tridiagCField

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

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',  maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell',  areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge',  dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge',  dvEdge)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_field(scratchPool, 'gradDensityEdge', gradDensityEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradDensityTopOfEdge', gradDensityTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradDensityConstZTopOfEdge', gradDensityConstZTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'dDensityDzTopOfCell', dDensityDzTopOfCellField)
      call mpas_pool_get_field(scratchPool, 'dDensityDzTopOfEdge', dDensityDzTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradZMidEdge', gradZMidEdgeField)
      call mpas_pool_get_field(scratchPool, 'gradZMidTopOfEdge', gradZMidTopOfEdgeField)
      call mpas_pool_get_field(scratchPool, 'rightHandSide', rightHandSideField)
      call mpas_pool_get_field(scratchPool, 'tridiagA', tridiagAField)
      call mpas_pool_get_field(scratchPool, 'tridiagB', tridiagBField)
      call mpas_pool_get_field(scratchPool, 'tridiagC', tridiagCField)
      call mpas_pool_get_field(scratchPool, 'areaCellSum', areaCellSumField)

      call mpas_allocate_scratch_field(gradDensityEdgeField, .True.)
      call mpas_allocate_scratch_field(gradDensityTopOfEdgeField, .True.)
      call mpas_allocate_scratch_field(gradDensityConstZTopOfEdgeField, .True.)
      call mpas_allocate_scratch_field(dDensityDzTopOfCellField, .True.)
      call mpas_allocate_scratch_field(dDensityDzTopOfEdgeField, .True.)
      call mpas_allocate_scratch_field(gradZMidEdgeField, .True.)
      call mpas_allocate_scratch_field(gradZMidTopOfEdgeField, .True.)
      call mpas_allocate_scratch_field(rightHandSideField, .True.)
      call mpas_allocate_scratch_field(tridiagAField, .True.)
      call mpas_allocate_scratch_field(tridiagBField, .True.)
      call mpas_allocate_scratch_field(tridiagCField, .True.)
      call mpas_allocate_scratch_field(areaCellSumField, .True.)

      gradDensityEdge => gradDensityEdgeField % array
      gradDensityTopOfEdge => gradDensityTopOfEdgeField % array
      gradDensityConstZTopOfEdge => gradDensityConstZTopOfEdgeField % array
      dDensityDzTopOfCell => dDensityDzTopOfCellField % array
      dDensityDzTopOfEdge => dDensityDzTopOfEdgeField % array
      gradZMidEdge => gradZMidEdgeField % array
      gradZMidTopOfEdge => gradZMidTopOfEdgeField % array
      rightHandSide => rightHandSideField % array
      tridiagA => tridiagAField % array
      tridiagB => tridiagBField % array
      tridiagC => tridiagCField % array
      areaCellSum => areaCellSumField % array

      ! Assign a huge value to the scratch variables which may manifest itself when
      ! there is a bug.
      gradDensityEdge(:,:) = huge(0D0)
      gradDensityTopOfEdge(:,:) = huge(0D0)
      dDensityDzTopOfCell(:,:) = huge(0D0)
      dDensityDzTopOfEdge(:,:) = huge(0D0)
      gradZMidEdge(:,:) = huge(0D0)
      gradZMidTopOfEdge(:,:) = huge(0D0)

      relativeSlopeTopOfEdge(:,:) = 0.0_RKIND
      relativeSlopeTopOfCell(:,:) = 0.0_RKIND
      relativeSlopeTapering(:,:) = 0.0_RKIND
      relativeSlopeTaperingCell(:,:) = 0.0_RKIND
      k33(:,:) = 0.0_RKIND
      normalGMBolusVelocity(:,:) = 0.0_RKIND
      
      !--------------------------------------------------------------------
      !
      ! Compute vertical derivative of density at top of cell, interpolate to top of edge
      ! This is required for Redi and Bolus parts.
      !
      !--------------------------------------------------------------------
      
      ! Compute vertical derivative of density (dDensityDzTopOfCell) at cell center and layer interface
      ! Note that displacedDensity is used from the upper cell, so that the EOS reference level for 
      ! pressure is the same for both displacedDensity(k-1,iCell) and density(k,iCell).
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

      ! Interpolate dDensityDzTopOfCell to edge and layer interface
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)+1
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            dDensityDzTopOfEdge(k,iEdge) = 0.5_RKIND * (dDensityDzTopOfCell(k,cell1) + dDensityDzTopOfCell(k,cell2))
         end do
      end do

      !--------------------------------------------------------------------
      !
      ! Compute horizontal gradient and mid-layer of edge, interpolate to top of edge
      ! This is required for Redi and Bolus parts.
      !
      !--------------------------------------------------------------------

      ! Compute density gradient (gradDensityEdge) and gradient of zMid (gradZMidEdge) 
      ! along the constant coordinate surface.
      ! The computed variables lives at edge and mid-layer depth
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            gradDensityEdge(k,iEdge) = (density(k,cell2) - density(k,cell1)) / dcEdge(iEdge)
            gradZMidEdge(k,iEdge) = (zMid(k,cell2) - zMid(k,cell1)) / dcEdge(iEdge)
         end do
      end do

      ! Interpolate gradDensityEdge and gradZMidEdge to layer interface
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

      !--------------------------------------------------------------------
      !
      ! Compute horizontal gradient required for Bolus part (along constant z)
      !
      !--------------------------------------------------------------------

      do iEdge = 1, nEdges
         if (maxLevelEdgeTop(iEdge) .GE. 1) then
            do k = 1, maxLevelEdgeTop(iEdge)+1
               gradDensityConstZTopOfEdge(k,iEdge) = gradDensityTopOfEdge(k,iEdge) - dDensityDzTopOfEdge(k,iEdge) * gradZMidTopOfEdge(k,iEdge)
            end do
         end if
      end do

      !--------------------------------------------------------------------
      !
      ! Compute relative slope and k33 for Redi part of GM.
      ! These variables are used in del2 velocity tendency routines.
      !
      !--------------------------------------------------------------------

      ! Compute relativeSlopeTopOfEdge at edge and layer interface
      ! set relativeSlopeTopOfEdge to zero for horizontal land/water edges.
      relativeSlopeTopOfEdge = 0.0_RKIND
      do iEdge = 1, nEdges

         ! Beside a full land cell (e.g. missing cell) maxLevelEdgeTop=0, so relativeSlopeTopOfEdge at that edge will remain zero.
         do k = 2, maxLevelEdgeTop(iEdge)
            relativeSlopeTopOfEdge(k,iEdge) = - gradDensityTopOfEdge(k,iEdge) / min(dDensityDzTopOfEdge(k,iEdge),-epsGM)
         end do

         ! Since dDensityDzTopOfEdge is guaranteed to be zero on the top surface, relativeSlopeTopOfEdge on the top surface is identified with its value on the second interface.
         relativeSlopeTopOfEdge(1,iEdge) = relativeSlopeTopOfEdge(2,iEdge)

         ! dDensityDzTopOfEdge may or may not equal zero on the bottom surface, depending on whether maxLevelEdgeTop(iEdge) = maxLevelEdgeBottom(iEdge). But here we
         ! take a simplistic approach and identify relativeSlopeTopOfEdge on the bottom surface with its value on the interface just above.
         relativeSlopeTopOfEdge( maxLevelEdgeTop(iEdge)+1, iEdge ) = relativeSlopeTopOfEdge( max(1,maxLevelEdgeTop(iEdge)), iEdge )

      end do

      ! slope can be unbounded in regions of neutral stability, reset to the large, but bounded, value
      ! values is hardwrite to 1.0, this is equivalent to a slope of 45 degrees
      where(relativeSlopeTopOfEdge < -1.0_RKIND) relativeSlopeTopOfEdge = -1.0_RKIND 
      where(relativeSlopeTopOfEdge >  1.0_RKIND) relativeSlopeTopOfEdge =  1.0_RKIND

      ! average relative slope to cell centers
      ! do this by computing (relative slope)^2, then taking sqrt
      areaCellSum = 1.0e-34_RKIND
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        ! contribution of cell area from this edge:
        areaEdge = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)

        do k = 1, maxLevelEdgeTop(iEdge)

           ! only one component is summed (thus the weighting by a factor of 2.0)
           rtmp = 2.0_RKIND * areaEdge * relativeSlopeTopOfEdge(k,iEdge)**2 
           relativeSlopeTopOfCell(k,cell1) = relativeSlopeTopOfCell(k,cell1) + rtmp 
           relativeSlopeTopOfCell(k,cell2) = relativeSlopeTopOfCell(k,cell2) + rtmp 

           areaCellSum(k,cell1) = areaCellSum(k,cell1) + areaEdge
           areaCellSum(k,cell2) = areaCellSum(k,cell2) + areaEdge

        end do
      end do
      do iCell=1,nCells
        do k = 1, maxLevelCell(iCell)
           relativeSlopeTopOfCell(k,iCell) = sqrt(relativeSlopeTopOfCell(k,iCell)/areaCellSum(k,iCell))
        end do
      end do

      ! Compute tapering function
      ! Compute k33 at cell center and layer interface
      k33(:,:) = 0.0_RKIND
      do iCell=1,nCells
        do k = 2, maxLevelCell(iCell)
          relativeSlopeTaperingCell(k,iCell) = min(1.0_RKIND, config_max_relative_slope**2 / (relativeSlopeTopOfCell(k,iCell)**2+epsGM))
          k33(k,iCell) = relativeSlopeTaperingCell(k,iCell) * (relativeSlopeTopOfCell(k,iCell))**2
        end do
      end do

      ! average tapering function to layer edges
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        do k = 2, maxLevelEdgeTop(iEdge)
          relativeSlopeTapering(k,iEdge) = 0.5_RKIND * (relativeSlopeTaperingCell(k,cell1) + relativeSlopeTaperingCell(k,cell2))
        enddo
      enddo

      ! k33 is still non-dimensional measuring the limited (relative slope)^2 of neutral surfaces.
      ! scale k33 by config_Redi_kappa so it has units of diffusivity
      k33 = config_Redi_kappa * k33

      ! allow disabling of K33 for testing
      if(config_disable_redi_k33) k33=0.0_RKIND

      !--------------------------------------------------------------------
      !
      ! Compute stream function and Bolus velocity for Bolus part of GM
      !
      !--------------------------------------------------------------------

      gmStreamFuncTopOfEdge(:,:) = 0.0_RKIND
      c = config_gravWaveSpeed_trunc**2
      do iEdge = 1, nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! Construct the tridiagonal matrix
         if (maxLevelEdgeTop(iEdge) .GE. 3) then
            ! First row
            k = 2       
            BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
            BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
            tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2/(layerThicknessEdge(k-1,iEdge)*layerThicknessEdge(k,iEdge)) - BruntVaisalaFreqTopEdge
            tridiagC(k-1) = 2.0_RKIND * config_gravWaveSpeed_trunc**2/layerThicknessEdge(k,iEdge)/(layerThicknessEdge(k-1,iEdge)+layerThicknessEdge(k,iEdge))
            rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / config_density0 * gradDensityConstZTopOfEdge(k,iEdge)

            ! Second to next to the last rows
            do k = 3, maxLevelEdgeTop(iEdge)-1
               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagA(k-2) = 2.0_RKIND * config_gravWaveSpeed_trunc**2/layerThicknessEdge(k-1,iEdge)/(layerThicknessEdge(k-1,iEdge)+layerThicknessEdge(k,iEdge))
               tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2/(layerThicknessEdge(k-1,iEdge)*layerThicknessEdge(k,iEdge)) - BruntVaisalaFreqTopEdge
               tridiagC(k-1) = 2.0_RKIND * config_gravWaveSpeed_trunc**2/layerThicknessEdge(k,iEdge)/(layerThicknessEdge(k-1,iEdge)+layerThicknessEdge(k,iEdge))
               rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / config_density0 * gradDensityConstZTopOfEdge(k,iEdge)
            end do

            ! Last row
            k = maxLevelEdgeTop(iEdge)                
            BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
            BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
            tridiagA(k-2) = 2.0_RKIND * config_gravWaveSpeed_trunc**2/layerThicknessEdge(k-1,iEdge)/(layerThicknessEdge(k-1,iEdge)+layerThicknessEdge(k,iEdge))
            tridiagB(k-1) = - 2.0_RKIND * config_gravWaveSpeed_trunc**2/(layerThicknessEdge(k-1,iEdge)*layerThicknessEdge(k,iEdge)) - BruntVaisalaFreqTopEdge
            rightHandSide(k-1) = config_standardGM_tracer_kappa * gravity / config_density0 * gradDensityConstZTopOfEdge(k,iEdge)

            ! Total number of rows
            N = maxLevelEdgeTop(iEdge) - 1

            ! Call the tridiagonal solver
            call tridiagonal_solve(tridiagA, tridiagB, tridiagC, rightHandSide, gmStreamFuncTopOfEdge(2:maxLevelEdgeTop(iEdge),iEdge), N)
         end if

      end do

      ! Compute normalGMBolusVelocity from the stream function
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)
            normalGMBolusVelocity(k,iEdge) = (gmStreamFuncTopOfEdge(k,iEdge) - gmStreamFuncTopOfEdge(k+1,iEdge)) / layerThicknessEdge(k,iEdge)
         end do
      end do

      ! Interpolate gmStreamFuncTopOfEdge to cell centers for visualization
      gmStreamFuncTopOfCell(:,:) = 0.0_RKIND
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         areaEdge = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)

         do k = 1, maxLevelEdgeTop(iEdge)
            rtmp = 0.5_RKIND * ( gmStreamFuncTopOfEdge(k,iEdge) + gmStreamFuncTopOfEdge(k+1,iEdge) ) * areaEdge
            gmStreamFuncTopOfCell(k,cell1) = gmStreamFuncTopOfCell(k,cell1) + rtmp
            gmStreamFuncTopOfCell(k,cell2) = gmStreamFuncTopOfCell(k,cell2) + rtmp
         end do

      end do
      do iCell = 1, nCells
         gmStreamFuncTopOfCell(:, iCell) = gmStreamFuncTopOfCell(:,iCell) / areaCell(iCell)
      end do

      ! Deallocate scratch variables
      call mpas_deallocate_scratch_field(gradDensityEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradDensityTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradDensityConstZTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(dDensityDzTopOfCellField, .true.)
      call mpas_deallocate_scratch_field(dDensityDzTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradZMidEdgeField, .true.)
      call mpas_deallocate_scratch_field(gradZMidTopOfEdgeField, .true.)
      call mpas_deallocate_scratch_field(rightHandSideField, .true.)
      call mpas_deallocate_scratch_field(tridiagAField, .true.)
      call mpas_deallocate_scratch_field(tridiagBField, .true.)
      call mpas_deallocate_scratch_field(tridiagCField, .true.)

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

      call mpas_timer_start("tridiagonal_solve")
      
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

      call mpas_timer_stop("tridiagonal_solve")
      
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
      call mpas_pool_get_config(ocnConfigs, 'config_density0',config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_Redi_kappa', config_Redi_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_use_standardGM',config_use_standardGM)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_redi_k33',config_disable_redi_k33)

   end subroutine ocn_gm_init!}}}

!***********************************************************************

end module ocn_gm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
