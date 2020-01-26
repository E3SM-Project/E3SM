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
   use ocn_config

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

   public :: ocn_GM_compute_Bolus_velocity, &
             ocn_GM_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   private :: tridiagonal_solve

   ! Config options
   real (kind=RKIND), parameter :: epsGM = 1.0e-12_RKIND

   logical :: local_config_GM_lat_variable_c2
   logical :: local_config_GM_kappa_lat_depth_variable
   real (kind=RKIND), parameter :: local_config_GM_min_phase_speed = 0.1_RKIND

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_GM_compute_Bolus_velocity
!
!> \brief   Computes GM Bolus velocity
!> \details
!>  This routine is the main driver for the Gent-McWilliams (GM) parameterization.
!>  It computes horizontal and vertical density gradients, the slope
!>  of isopycnal surfaces, and solves a boundary value problem in each column
!>  for the stream function, which is used to compute the Bolus velocity.
!
!-----------------------------------------------------------------------

   subroutine ocn_GM_compute_Bolus_velocity(statePool, diagnosticsPool, &
      meshPool, scratchPool, timeLevelIn)
      !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool,           &! pool containing state variables
         meshPool              ! pool containing mesh quantities

      integer, intent(in), optional :: &
         timeLevelIn          ! time level for state variables

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool,     &! pool containing some diagnostics
         scratchPool           ! pool containing some scratch space

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:,:), pointer :: density, displacedDensity, zMid, normalGMBolusVelocity, &
         layerThicknessEdge, gradDensityEdge, gradDensityTopOfEdge, gradDensityConstZTopOfEdge, gradZMidEdge, &
         gradZMidTopOfEdge, k33, gmStreamFuncTopOfEdge, BruntVaisalaFreqTop, &
         gmStreamFuncTopOfCell, dDensityDzTopOfEdge, dDensityDzTopOfCell, &
         layerThickness, &
         areaCellSum, inSituThermalExpansionCoeff, inSituSalineContractionCoeff, kappaGM3D

      real(kind=RKIND), dimension(:), pointer   :: gmBolusKappa,  cGMphaseSpeed, bottomDepth, indexBoundaryLayerDepth
      real(kind=RKIND), dimension(:,:,:), pointer :: slopeTriadUp, slopeTriadDown
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: maxLevelEdgeTop, maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell
      integer                          :: i, k, iEdge, cell1, cell2, iCell, N, iter, iCellSelf
      real(kind=RKIND)                 :: h1, h2, areaEdge, c, BruntVaisalaFreqTopEdge, rtmp, stmp, maxSlopeK33
      real(kind=RKIND)                 :: bottomAv, sumN2, countN2, maxN, kappaSum, ltSum
      real(kind=RKIND) :: dcEdgeInv, drhoDx, drhoDT, drhoDS, dTdx, dSdx
      real(kind=RKIND) :: surfaceTaperCoef, bottomTaperCoef
      real(kind=RKIND), dimension(:), allocatable :: dzTop, dTdzTop, dSdzTop, k33Norm
      real(kind=RKIND), dimension(:), allocatable :: slopeTaper, slopeTaperSurface

      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      type (field2DReal), pointer :: gradDensityEdgeField, gradDensityTopOfEdgeField, gradDensityConstZTopOfEdgeField, &
         gradZMidEdgeField, gradZMidTopOfEdgeField, dDensityDzTopOfCellField, dDensityDzTopOfEdgeField,areaCellSumField

      type (mpas_pool_type), pointer :: tracersPool
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, pointer :: indexTemperature, indexSalinity
      integer :: timeLevel
      integer :: iBLD

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_timer_start('gm bolus velocity')

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_subpool(statePool,   'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff',inSituThermalExpansionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      call mpas_pool_get_array(diagnosticsPool, 'cGMphaseSpeed', cGMphaseSpeed)
      call mpas_pool_get_array(diagnosticsPool, 'kappaGM3D', kappaGM3D)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadUp', slopeTriadUp)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadDown', slopeTriadDown)
      call mpas_pool_get_array(diagnosticsPool, 'k33', k33)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', gmStreamFuncTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfCell', gmStreamFuncTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa', gmBolusKappa)
      call mpas_pool_get_array(diagnosticsPool, 'indexBoundaryLayerDepth', indexBoundaryLayerDepth)

      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',  maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell',  maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell',  areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge',  dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge',  dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell',  nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell',  edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

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

! mrp performance note: are these actually needed? check if they are re-initialized below
      ! Assign a huge value to the scratch variables which may manifest itself when
      ! there is a bug.
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            gradDensityEdge(k, iEdge) = huge(0D0)
            gradZMidEdge(k, iEdge) = huge(0D0)
            normalGMBolusVelocity(k, iEdge) = 0.0_RKIND
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels + 1
            gradDensityTopOfEdge(k, iEdge) = huge(0D0)
            dDensityDzTopOfEdge(k, iEdge) = huge(0D0)
            gradZMidTopOfEdge(k, iEdge) = huge(0D0)
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells + 1
         do k = 1, nVertLevels
            dDensityDzTopOfCell(k,  iCell) = huge(0D0)
            k33(k, iCell) = 0.0_RKIND
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
! mrp performance note: are all these variables, like gradDensityEdge and gradZMidEdge actually
! needed?
! mrp performance note: we shouldn't have so many fine grained loops. Merge loops, if they are
! kept.
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
      ! using triads and linear expansion of rho(T,S).
      ! See Griffies et al. 1998 JPO eqn 31 and 32.
      !
      !--------------------------------------------------------------------

      ! Slope calculation:
      ! S = - drho/dx / drho/dz
      ! Expand in T and S using equation of state:
      ! S = -  ( drho/dT * dT/dx + drho/dS * dS/dx ) &
      !      / ( drho/dT * dT/dz + drho/dS * dS/dz )

      ! From Griffies 1998 P 820:
      ! "Should any one of these subfluxes lie adjacent to a solid boundary, its
      ! contribution to the average over the four triads is eliminated."
      ! So we do not use the triad slope beside boundaries or sea surface.

      ! mrp note, we may want to take slope calculation out of this Bolus
      ! routine.
      ! mrp performance note: check which halo is needed here.
      nCells = nCellsArray( 2 )

      allocate(dzTop(nVertLevels+1))
      allocate(dTdzTop(nVertLevels+1))
      allocate(dSdzTop(nVertLevels+1))
      allocate(k33Norm(nVertLevels+1))
      allocate(slopeTaper(nVertLevels))
      allocate(slopeTaperSurface(nVertLevels))
      ! Prepare tapering for this cell. This is the slope of the bottom taper line.
      surfaceTaperCoef = 1.0_RKIND/config_Redi_surface_taper_layers
      bottomTaperCoef = 1.0_RKIND/config_Redi_bottom_taper_layers
      slopeTaperSurface(:) = 1.0_RKIND
      slopeTaperSurface(1) = 0.0_RKIND
      do k = 2, config_Redi_surface_taper_layers
          slopeTaperSurface(k) = slopeTaperSurface(k)*(k-1)*surfaceTaperCoef
      end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ! prep dz, dTdz and dSdz for this column
         do k = 2, maxLevelCell(iCell)
            dzTop(k) = 0.5_RKIND * ( layerThickness(k-1,iCell) + layerThickness(k,iCell) )
            dTdzTop(k) = (  activeTracers(indexTemperature,k-1,iCell) &
                          - activeTracers(indexTemperature,k  ,iCell) ) &
                         / dzTop(k)
            dSdzTop(k) = (  activeTracers(indexSalinity,k-1,iCell) &
                          - activeTracers(indexSalinity,k  ,iCell) ) &
                         / dzTop(k)
         end do
         ! mrp performance: could move above cell loop later
         dzTop(1) = -1e34_RKIND
         dTdzTop(1) = -1e34_RKIND
         dSdzTop(1) = -1e34_RKIND
         dzTop(maxLevelCell(iCell)+1) = -1e34_RKIND
         dTdzTop(maxLevelCell(iCell)+1) = -1e34_RKIND
         dSdzTop(maxLevelCell(iCell)+1) = -1e34_RKIND

         ! Add bottom taper slope for this column.
         ! Top taper slope is computed before the cell loop, because it is the same
         ! for all cells.
         slopeTaper(:) = slopeTaperSurface(:)
         do k = max(maxLevelCell(iCell) - config_Redi_bottom_taper_layers, 2), maxLevelCell(iCell)
            slopeTaper(k) = slopeTaper(k) &
                * (maxLevelCell(iCell)-k)*bottomTaperCoef
         end do

         ! Add taper for boundary layer
         if (config_Redi_zero_in_boundary_layer) then
            iBLD = int(indexBoundaryLayerDepth(iCell))
            do k = 2, iBLD
               slopeTaper(k) = 0.0_RKIND
            end do
            do k = iBLD + 1, &
                   min(iBLD + config_Redi_BLD_taper_layers, maxLevelCell(iCell))
               slopeTaper(k) = slopeTaper(k)*(k-iBLD) / config_Redi_BLD_taper_layers
            end do
         end if

         k33(1:maxLevelCell(iCell)+1,iCell) = 0.0_RKIND
         k33Norm(1:maxLevelCell(iCell)+1) = epsGM
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            ! mrp Performance Note: Replace with a logical array on init later
            if (cell1==iCell) then
                iCellSelf = 1
            elseif (cell2==iCell) then
                iCellSelf = 2
            else
                print *, 'iCellSelf, something went wrong',cell1,cell2,iCell,nCells
            endif
            dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
            areaEdge = dcEdge(iEdge) * dvEdge(iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)
               drhoDT = -inSituThermalExpansionCoeff(k,iCell)
               drhoDS = inSituSalineContractionCoeff(k,iCell)
               dTdx = (  activeTracers(indexTemperature,k,cell2) &
                       - activeTracers(indexTemperature,k,cell1) ) &
                      * dcEdgeInv
               dSdx = (  activeTracers(indexSalinity,k,cell2) &
                       - activeTracers(indexSalinity,k,cell1) ) &
                      * dcEdgeInv
               drhoDx = drhoDT * dTdx + drhoDS * dSdx

               ! Always compute *Up on the top cell and *Down on the bottom
               ! cell, even though they are never used. This avoids an if
               ! statement or separate computation for top and bottom.
               slopeTriadUp(k,iCellSelf,iEdge) = &
                  - drhoDx / &
                  (  drhoDT * dTdzTop(k) &
                   + drhoDS * dSdzTop(k) )
               slopeTriadDown(k,iCellSelf,iEdge) = &
                  - drhoDx / &
                  (  drhoDT * dTdzTop(k+1) &
                   + drhoDS * dSdzTop(k+1) )

               ! Slope can be unbounded in regions of neutral stability.
               ! Reset slope to maximum slope limit.
               slopeTriadUp(k,iCellSelf,iEdge) = &
                  slopeTaper(k) * &
                  max( min( slopeTriadUp(k,iCellSelf,iEdge), &
                  config_mesoscale_eddy_isopycnal_slope_limit), -config_mesoscale_eddy_isopycnal_slope_limit)
               slopeTriadDown(k,iCellSelf,iEdge) = &
                  slopeTaper(k) * &
                  max( min( slopeTriadDown(k,iCellSelf,iEdge), &
                  config_mesoscale_eddy_isopycnal_slope_limit), -config_mesoscale_eddy_isopycnal_slope_limit)

               ! Griffies 1998 eqn 34
               k33(k  ,iCell) = k33(k  ,iCell) + areaEdge*dzTop(k  )*slopeTriadUp(k,iCellSelf,iEdge)**2
               k33(k+1,iCell) = k33(k+1,iCell) + areaEdge*dzTop(k+1)*slopeTriadDown(k,iCellSelf,iEdge)**2
               k33Norm(k  ) = k33Norm(k  ) + areaEdge*dzTop(k  )
               k33Norm(k+1) = k33Norm(k+1) + areaEdge*dzTop(k+1)

               ! Could put in vertical smoothing of slopes here if we wanted.
               ! Oct 2019 mrp, we decided not to smooth the slope variables for
               ! now.

            end do
         end do

         ! Normalize k33
         do k = 2, maxLevelCell(iCell)
            k33(k,iCell) = k33(k,iCell) / k33Norm(k)
         end do
      end do ! iCell
      !$omp end do

      deallocate(dzTop)
      deallocate(dTdzTop)
      deallocate(dSdzTop)
      deallocate(k33Norm)
      deallocate(slopeTaper)
      deallocate(slopeTaperSurface)

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

      if (local_config_GM_lat_variable_c2) then
         !$omp do schedule(runtime) private(cell1, cell2, sumN2, ltSum, countN2, BruntVaisalaFreqTopEdge)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            sumN2 = 0.0
            ltSum = 0.0
            countN2 = 0
            
            do k=2,maxLevelEdgeTop(iEdge)

               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
                
               sumN2 = sumN2 + BruntVaisalaFreqTopEdge*layerThicknessEdge(k,iEdge)
               ltSum = ltSum + layerThicknessEdge(k,iEdge)
               countN2 = countN2 +1

            enddo

            if(countN2 > 0) cGMphaseSpeed(iEdge) = max(local_config_GM_min_phase_speed ,sqrt(sumN2/ltSum)*ltSum / 3.141592_RKIND)

         enddo
         !$omp end do

      else
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            cGMphaseSpeed(iEdge) = config_GM_constant_gravWaveSpeed
         enddo
         !$omp end do
      endif

      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         kappaGM3D(:,iEdge) = gmBolusKappa(iEdge)
      enddo 
      !$omp end do

      if (local_config_GM_kappa_lat_depth_variable) then

         !$omp do schedule(runtime) private(cell1, cell2, k, BruntVaisalaFreqTopEdge, maxN)
         do iEdge = 1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            maxN = -1.0_RKIND
            do k=2,maxLevelEdgeTop(iEdge)
               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

               maxN = max(maxN,BruntVaisalaFreqTopEdge)

            enddo

            do k=2,maxLevelEdgeTop(iEdge)
               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

               kappaGM3D(k,iEdge) = gmBolusKappa(iEdge)*max(config_GM_min_stratification_ratio, &
                       BruntVaisalaFreqTopEdge / (maxN + 1.0E-10_RKIND))
            enddo
         enddo
         !$omp end do
      endif

      nEdges = nEdgesArray( 3 )

      !$omp do schedule(runtime)
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
            tridiagB(k-1) = - 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / (layerThicknessEdge(k-1,iEdge) &
                          * layerThicknessEdge(k,iEdge)) - BruntVaisalaFreqTopEdge
            tridiagC(k-1) = 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / layerThicknessEdge(k, iEdge) &
                          / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
            rightHandSide(k-1) = kappaGM3D(k-1,iEdge) * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)

            ! Second to next to the last rows
            do k = 3, maxLevelEdgeTop(iEdge)-1
               BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagA(k-2) = 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / layerThicknessEdge(k-1, iEdge) &
                             / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
               tridiagB(k-1) = - 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / (layerThicknessEdge(k-1, iEdge) &
                             * layerThicknessEdge(k, iEdge) ) - BruntVaisalaFreqTopEdge
               tridiagC(k-1) = 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / layerThicknessEdge(k, iEdge) &
                             / (layerThicknessEdge(k-1, iEdge) + layerThicknessEdge(k, iEdge))
               rightHandSide(k-1) = kappaGM3D(k-1,iEdge) * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)
            end do

            ! Last row
            k = maxLevelEdgeTop(iEdge)
            BruntVaisalaFreqTopEdge = 0.5_RKIND * (BruntVaisalaFreqTop(k,cell1) + BruntVaisalaFreqTop(k,cell2))
            BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
            tridiagA(k-2) = 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / layerThicknessEdge(k-1,iEdge) &
                          / (layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge))
            tridiagB(k-1) = - 2.0_RKIND * cGMphaseSpeed(iEdge)**2 / (layerThicknessEdge(k-1, iEdge) &
                          * layerThicknessEdge(k, iEdge)) - BruntVaisalaFreqTopEdge
            rightHandSide(k-1) = kappaGM3D(k-1,iEdge) * gravity / rho_sw * gradDensityConstZTopOfEdge(k,iEdge)

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

   end subroutine ocn_GM_compute_Bolus_velocity!}}}

!***********************************************************************
!
!  routine tridiagonal_solve
!
!> \brief   Solve the matrix equation Ax=r for x, where A is tridiagonal.
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
!  routine ocn_GM_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \details
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_GM_init(domain, err)!{{{

      type (domain_type), intent(in) :: domain

      integer, intent(out) :: err !< Output: error flag

      real (kind=RKIND) :: avgCellDiameter, sqrtPiInv

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      real (kind=RKIND), dimension(:), pointer :: areaCell

      integer :: iEdge
      integer, pointer :: nEdges
      integer, pointer :: nVertLevels
      real(kind=RKIND), dimension(:), pointer   :: gmBolusKappa
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOncell

      err = 0

      if ((.not.config_use_GM).and.(.not.config_use_Redi)) return

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa', gmBolusKappa)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

         if ((config_Redi_surface_taper_layers <= 0).or. &
             (config_Redi_surface_taper_layers > nVertLevels )) then
            call mpas_log_write( 'config_Redi_surface_taper_layers must be between 1 and nVertLevels. ', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if

         if ((config_Redi_bottom_taper_layers <= 0).or. &
             (config_Redi_bottom_taper_layers > nVertLevels )) then
            call mpas_log_write( 'config_Redi_bottom_taper_layers must be between 1 and nVertLevels. ', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if

         if (config_GM_closure == 'constant') then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_kappa
            end do
            !$omp end do
         else if (config_GM_closure == 'N2_dependent') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .true.
 
            ! for N2 dependence, we still assign Kappa as a constant.
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_kappa
            end do
            !$omp end do
         else
            call mpas_log_write( 'Invalid choice of config_GM_closure.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if

         ! Add resolution taper
         if (config_eddying_resolution_taper == 'none') then
            ! Nothing to do, as we just keep the same assignment as above.
         else if (config_eddying_resolution_taper == 'ramp') then
            sqrtPiInv = 1.0_RKIND / sqrt(pii)
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               avgCellDiameter = sqrtPiInv*(sqrt(areaCell(cellsOnEdge(1,iEdge))) + sqrt(areaCell(cellsOnEdge(2,iEdge))))
               if (avgCellDiameter <= config_eddying_resolution_ramp_min) then
                  gmBolusKappa(iEdge) = 0.0_RKIND
               else if (avgCellDiameter >= config_eddying_resolution_ramp_max) then
                  ! do nothing, i.e. gmBolusKappa(iEdge) remains unchanged.
               else
                  gmBolusKappa(iEdge) = gmBolusKappa(iEdge) &
                     /(config_eddying_resolution_ramp_max - config_eddying_resolution_ramp_min) &
                     *(avgCellDiameter - config_eddying_resolution_ramp_min)
               end if
            end do
            !$omp end do
         else
            call mpas_log_write( 'Invalid choice of config_eddying_resolution_taper.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain % dminfo)
         end if


         block => block % next
      end do
   end subroutine ocn_GM_init!}}}

!***********************************************************************

end module ocn_gm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
