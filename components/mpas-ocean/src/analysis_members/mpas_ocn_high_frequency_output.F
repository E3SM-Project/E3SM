!Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_high_frequency_output
!
!> \brief MPAS ocean analysis mode member: high_frequency_output
!> \author Todd Ringler and Phillip Wolfram
!> \date   2015/06/12
!> \details
!>  MPAS ocean analysis mode member: high_frequency_output
!>
!-----------------------------------------------------------------------

module ocn_high_frequency_output

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_config
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

   public :: ocn_init_high_frequency_output, &
             ocn_compute_high_frequency_output, &
             ocn_restart_high_frequency_output, &
             ocn_finalize_high_frequency_output

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Replace calls to mpas_allocate_scratch_Field with straight allocates.
   ! These variables were moved from the subroutines, but will be moved back once the
   ! threading implementation is changed to more localized threading
   !
   !--------------------------------------------------------------------
   ! normalThicknessFluxSum: Flux of thickness through an edge over full column, sum(h*u)
   !                  units: m^{2} s^{-1}
   real (kind=RKIND), dimension(:), allocatable :: normalThicknessFluxSum
   ! layerThicknessSumEdge: total thickness at edge, sum(h)
   !                 units: m
   real (kind=RKIND), dimension(:), allocatable :: layerThicknessSumEdge

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_high_frequency_output
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_high_frequency_output(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_init_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_compute_high_frequency_output
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_high_frequency_output(domain, timeLevel, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: highFrequencyOutputAMPool
      type (mpas_pool_type), pointer :: tracersPool

      integer :: iLevel, iCell, iEdge, i, cell1, cell2, k, eoe
      integer :: iLevel0100, iLevel0250, iLevel0700, iLevel2000
      real (kind=RKIND) :: sumLayerThickness
      integer, pointer :: nVertLevels, nCells, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelEdgeTop, maxLevelEdgeBot, maxLevelCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge, edgesOnEdge

      real (kind=RKIND) :: invAreaCell1, layerThicknessEdge1, coeff, weightedNormalVel, cellArea
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, kineticEnergyAt250m, kineticEnergyAtSurface, relativeVorticityAt250m, relativeVorticityAtSurface
      real (kind=RKIND), dimension(:), pointer :: divergenceAt250m, relativeVorticityVertexAt250m
      real (kind=RKIND), dimension(:), pointer :: divergenceAtBottom,relativeVorticityAtBottom,kineticEnergyAtBottom
      real (kind=RKIND), dimension(:), pointer :: vertVelAt250m
      real (kind=RKIND), dimension(:), pointer :: BruntVaisalaFreqTopAtSFC, BruntVaisalaFreqTopAt250m, BruntVaisalaFreqTopAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalVelAtSFC, normalVelAt250m, normalVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialVelAtSFC, tangentialVelAt250m, tangentialVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalVelAtSFC, zonalVelAt250m, zonalVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalVelAtSFC, meridionalVelAt250m, meridionalVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: zonalAreaWeightedCellVelAtSFC, zonalAreaWeightedCellVelAt250m, zonalAreaWeightedCellVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalAreaWeightedCellVelAtSFC, meridionalAreaWeightedCellVelAt250m, meridionalAreaWeightedCellVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: tangentialBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: zonalBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: meridionalBarotropicVel

      real (kind=RKIND), dimension(:), pointer :: normalBaroclinicVelAtSFC, normalBaroclinicVelAt250m, normalBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialBaroclinicVelAtSFC, tangentialBaroclinicVelAt250m, tangentialBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalBaroclinicVelAtSFC, zonalBaroclinicVelAt250m, zonalBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalBaroclinicVelAtSFC, meridionalBaroclinicVelAt250m, meridionalBaroclinicVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalGMBolusVelAtSFC, normalGMBolusVelAt250m, normalGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialGMBolusVelAtSFC, tangentialGMBolusVelAt250m, tangentialGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalGMBolusVelAtSFC, zonalGMBolusVelAt250m, zonalGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalGMBolusVelAtSFC, meridionalGMBolusVelAt250m, meridionalGMBolusVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
      real (kind=RKIND), dimension(:), pointer :: vertTransportVelocityAt250m, vertGMvelocityAt250m
      real (kind=RKIND), dimension(:), pointer :: vertVelSFC, vertTransportVelocitySFC, vertGMvelocitySFC
      real (kind=RKIND), dimension(:), pointer :: barotropicSpeed, columnIntegratedSpeed, dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell, relativeVorticityCell
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAt250m, activeTracersAvgTopto0100, activeTracersAvg0100to0250
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAtSurface, activeTracersAtBottom
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAvg0250to0700, activeTracersAvg0700to2000, activeTracersAvg2000toBottom
      real (kind=RKIND), dimension(:,:), pointer :: relativeVorticity, divergence, layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocity ! baroclinic/barotropic pressume split explicit
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocity ! baroclinic/barotropic pressume split explicit
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         allocate(normalThicknessFluxSum(nEdges))
         allocate(layerThicknessSumEdge(nEdges))

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'highFrequencyOutputAM', highFrequencyOutputAMPool)

         ! get static data from mesh pool
         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

         ! get arrays that will be 'sliced' and put into high frequency output
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence, timeLevel)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

         ! get arrays that can be written to output at high freqency
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAt250m', kineticEnergyAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAtSurface', kineticEnergyAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAt250m', relativeVorticityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAtSurface', relativeVorticityAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'divergenceAt250m', divergenceAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAtBottom', relativeVorticityAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'divergenceAtBottom', divergenceAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAtBottom', kineticEnergyAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityVertexAt250m', relativeVorticityVertexAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAtSurface', activeTracersAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt250m', activeTracersAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAtBottom', activeTracersAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvgTopto0100', activeTracersAvgTopto0100)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0100to0250', activeTracersAvg0100to0250)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0250to0700', activeTracersAvg0250to0700)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0700to2000', activeTracersAvg0700to2000)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg2000toBottom', activeTracersAvg2000toBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertGMvelocitySFC', vertGMvelocitySFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertTransportVelocitySFC', vertTransportVelocitySFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertVelAt250m', vertVelAt250m)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAtSFC', normalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAt250m', normalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAtBottom', normalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAtSFC', tangentialVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAt250m', tangentialVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAtBottom', tangentialVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAtSFC', zonalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAt250m', zonalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAtBottom', zonalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAtSFC', meridionalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAt250m', meridionalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAtBottom', meridionalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAtSFC', zonalAreaWeightedCellVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAt250m', zonalAreaWeightedCellVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAtBottom', zonalAreaWeightedCellVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAtSFC', meridionalAreaWeightedCellVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAt250m', meridionalAreaWeightedCellVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAtBottom', meridionalAreaWeightedCellVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBarotropicVel', normalBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBarotropicVel', tangentialBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBarotropicVel', zonalBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBarotropicVel', meridionalBarotropicVel)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAtSFC', normalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAt250m', normalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAtBottom', normalBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAtSFC', tangentialBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAt250m', tangentialBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAtBottom', tangentialBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAtSFC', zonalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAt250m', zonalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAtBottom', zonalBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAtSFC', meridionalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAt250m', meridionalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAtBottom', meridionalBaroclinicVelAtBottom)

         if (config_use_gm) then
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAtSFC', normalGMBolusVelAtSFC)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAt250m', normalGMBolusVelAt250m)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAtBottom', normalGMBolusVelAtBottom)

            call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAtSFC', tangentialGMBolusVelAtSFC)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAt250m', tangentialGMBolusVelAt250m)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAtBottom', tangentialGMBolusVelAtBottom)

            call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAtSFC', zonalGMBolusVelAtSFC)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAt250m', zonalGMBolusVelAt250m)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAtBottom', zonalGMBolusVelAtBottom)

            call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAtSFC', meridionalGMBolusVelAtSFC)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAt250m', meridionalGMBolusVelAt250m)
            call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAtBottom', meridionalGMBolusVelAtBottom)
         end if

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAtSFC', BruntVaisalaFreqTopAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAt250m', BruntVaisalaFreqTopAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAtBottom', BruntVaisalaFreqTopAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertGMvelocityAt250m', vertGMvelocityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertTransportVelocityAt250m', vertTransportVelocityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertVelSFC', vertVelSFC)

         ! split explicit specific arrays
         if ( config_time_integrator == trim('split_explicit') .or.  config_time_integrator == trim('split_implicit') &
                                                               .or.  config_time_integrator == trim('split_explicit_ab2') ) then
           call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, 1)
           call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity, 1)
         endif

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
         call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
         call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'barotropicSpeed', barotropicSpeed)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'columnIntegratedSpeed', columnIntegratedSpeed)

         ! find vertical level that is just above the 100 m reference level
         ! if even the bottom level isn't deep enough, we still default to the bottom level
         iLevel0100 = nVertLevels
         do iLevel=1,nVertLevels
           if(refBottomDepth(iLevel) > 100.0_RKIND) then
              iLevel0100 = iLevel
              exit
           endif
         enddo

         ! find vertical level that is just above the 250 m reference level
         ! if even the bottom level isn't deep enough, we still default to the bottom level
         iLevel0250 = nVertLevels
         do iLevel=iLevel0100,nVertLevels
           if(refBottomDepth(iLevel) > 250.0_RKIND) then
              iLevel0250 = iLevel
              exit
           endif
         enddo

         ! find vertical level that is just above the 700 m reference level
         ! if even the bottom level isn't deep enough, we still default to the bottom level
         iLevel0700 = nVertLevels
         do iLevel=iLevel0250,nVertLevels
           if(refBottomDepth(iLevel) > 700.0_RKIND) then
              iLevel0700 = iLevel
              exit
           endif
         enddo

         ! find vertical level that is just above the 2000 m reference level
         ! if even the bottom level isn't deep enough, we still default to the bottom level
         iLevel2000 = nVertLevels
         do iLevel=iLevel0700,nVertLevels
           if(refBottomDepth(iLevel) > 2000.0_RKIND) then
              iLevel2000 = iLevel
              exit
           endif
         enddo

         ! copy data into high frequency output fields
         kineticEnergyAt250m(:) = kineticEnergyCell(iLevel0250,:)
         relativeVorticityAt250m(:) = relativeVorticityCell(iLevel0250,:)
         divergenceAt250m(:) = divergence(iLevel0250,:)
         relativeVorticityVertexAt250m(:) = relativeVorticity(iLevel0250,:)
         kineticEnergyAtSurface(:) = kineticEnergyCell(1,:)
         relativeVorticityAtSurface(:) = relativeVorticityCell(1,:)
         activeTracersAtSurface(1,:) = activeTracers(1,1,:)
         activeTracersAtSurface(2,:) = activeTracers(2,1,:)
         activeTracersAt250m(1,:) = activeTracers(1,iLevel0250,:)
         activeTracersAt250m(2,:) = activeTracers(2,iLevel0250,:)
         vertTransportVelocitySFC(:) = vertTransportVelocityTop(1,:)
         vertVelSFC(:) = vertVelocityTop(1,:)
         vertTransportVelocityAt250m(:) = vertTransportVelocityTop(iLevel0250,:)
         vertVelAt250m(:) = vertVelocityTop(iLevel0250,:)

         if (config_use_gm) then
            vertGMvelocitySFC(:) = vertGMBolusVelocityTop(1,:)
            vertGMvelocityAt250m(:) = vertGMBolusVelocityTop(iLevel0250,:)
            do iEdge = 1, nEdges
               normalGMBolusVelAtSFC(iEdge) = normalGMBolusVelocity(1,iEdge)
               normalGMBolusVelAt250m(iEdge) = normalGMBolusVelocity(iLevel0250,iEdge)
               normalGMBolusVelAtBottom(iEdge) = normalGMBolusVelocity(maxLevelEdgeBot(iEdge),iEdge)
            end do

            do iEdge = 1, nEdges
              tangentialGMBolusVelAtSFC(iEdge) = 0.0_RKIND
              tangentialGMBolusVelAt250m(iEdge) = 0.0_RKIND
              tangentialGMBolusVelAtBottom(iEdge) = 0.0_RKIND
              do i = 1, nEdgesOnEdge(iEdge)
                eoe = edgesOnEdge(i,iEdge)
                tangentialGMBolusVelAtSFC(iEdge) = tangentialGMBolusVelAtSFC(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAtSFC(eoe)
                tangentialGMBolusVelAt250m(iEdge) = tangentialGMBolusVelAt250m(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAt250m(eoe)
                tangentialGMBolusVelAtBottom(iEdge) = tangentialGMBolusVelAtBottom(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAtBottom(eoe)
              end do

              zonalGMBolusVelAtSFC(iEdge) = normalGMBolusVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
                - tangentialGMBolusVelAtSFC(iEdge)*sin(angleEdge(iEdge))
              meridionalGMBolusVelAtSFC(iEdge) = normalGMBolusVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
                + tangentialGMBolusVelAtSFC(iEdge)*cos(angleEdge(iEdge))

              zonalGMBolusVelAt250m(iEdge) = normalGMBolusVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
                - tangentialGMBolusVelAt250m(iEdge)*sin(angleEdge(iEdge))
              meridionalGMBolusVelAt250m(iEdge) = normalGMBolusVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
                + tangentialGMBolusVelAt250m(iEdge)*cos(angleEdge(iEdge))

              zonalGMBolusVelAtBottom(iEdge) = normalGMBolusVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
                - tangentialGMBolusVelAtBottom(iEdge)*sin(angleEdge(iEdge))
              meridionalGMBolusVelAtBottom(iEdge) = normalGMBolusVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
                + tangentialGMBolusVelAtBottom(iEdge)*cos(angleEdge(iEdge))
            end do
         end if

         do iEdge = 1, nEdges

           tangentialVelAtSFC(iEdge) = tangentialVelocity(1,iEdge)
           tangentialVelAt250m(iEdge) = tangentialVelocity(iLevel0250,iEdge)
           tangentialVelAtBottom(iEdge) = tangentialVelocity(maxLevelEdgeBot(iEdge),iEdge)

           normalVelAtSFC(iEdge) = normalVelocity(1,iEdge)
           normalVelAt250m(iEdge) = normalVelocity(iLevel0250,iEdge)
           normalVelAtBottom(iEdge) = normalVelocity(maxLevelEdgeBot(iEdge),iEdge)

           zonalVelAtSFC(iEdge) = normalVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAtSFC(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAtSFC(iEdge) = normalVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAtSFC(iEdge)*cos(angleEdge(iEdge))

           zonalVelAt250m(iEdge) = normalVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAt250m(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAt250m(iEdge) = normalVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAt250m(iEdge)*cos(angleEdge(iEdge))

           zonalVelAtBottom(iEdge) = normalVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAtBottom(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAtBottom(iEdge) = normalVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAtBottom(iEdge)*cos(angleEdge(iEdge))
         end do

         do iCell=1,nCells
           relativeVorticityAtBottom(iCell) = relativeVorticityCell(maxLevelCell(iCell),iCell)
           divergenceAtBottom(iCell) = divergence(maxLevelCell(iCell),iCell)
           kineticEnergyAtBottom(iCell) = kineticEnergyCell(maxLevelCell(iCell),iCell)
           activeTracersAtBottom(1,iCell) = activeTracers(1,maxLevelCell(iCell),iCell)
           activeTracersAtBottom(2,iCell) = activeTracers(2,maxLevelCell(iCell),iCell)

           BruntVaisalaFreqTopAtSFC(iCell) = BruntVaisalaFreqTop(1,iCell)
           BruntVaisalaFreqTopAt250m(iCell) = BruntVaisalaFreqTop(iLevel0250,iCell)
           BruntVaisalaFreqTopAtBottom(iCell) = BruntVaisalaFreqTop(maxLevelCell(iCell),iCell)
         end do

         do iCell = 1, nCells
           zonalAreaWeightedCellVelAtSFC(iCell) = 0.0_RKIND
           zonalAreaWeightedCellVelAt250m(iCell) = 0.0_RKIND
           zonalAreaWeightedCellVelAtBottom(iCell) = 0.0_RKIND

           meridionalAreaWeightedCellVelAtSFC(iCell) = 0.0_RKIND
           meridionalAreaWeightedCellVelAt250m(iCell) = 0.0_RKIND
           meridionalAreaWeightedCellVelAtBottom(iCell) = 0.0_RKIND

           cellArea = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               ! note, we don't necessarily have d_{i,j} computed (Fringer, 2006), so we allow averaging between cells via 0.5_RKIND*dcEdge(iEdge) to compute
               ! segment of cell area d_{i,j}-- the key is that the edge-normal velocities are a convex combination
               cellArea = cellArea + 0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)

               weightedNormalVel = normalVelAtSFC(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAtSFC(iCell) = zonalAreaWeightedCellVelAtSFC(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAtSFC(iCell) = meridionalAreaWeightedCellVelAtSFC(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

               weightedNormalVel = normalVelAt250m(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAt250m(iCell) = zonalAreaWeightedCellVelAt250m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAt250m(iCell) = meridionalAreaWeightedCellVelAt250m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

               weightedNormalVel = normalVelAtBottom(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAtBottom(iCell) = zonalAreaWeightedCellVelAtBottom(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAtBottom(iCell) = meridionalAreaWeightedCellVelAtBottom(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
           end do
           zonalAreaWeightedCellVelAtSFC(iCell) = zonalAreaWeightedCellVelAtSFC(iCell)/cellArea
           meridionalAreaWeightedCellVelAtSFC(iCell) = meridionalAreaWeightedCellVelAtSFC(iCell)/cellArea

           zonalAreaWeightedCellVelAt250m(iCell) = zonalAreaWeightedCellVelAt250m(iCell)/cellArea
           meridionalAreaWeightedCellVelAt250m(iCell) = meridionalAreaWeightedCellVelAt250m(iCell)/cellArea

           zonalAreaWeightedCellVelAtBottom(iCell) = zonalAreaWeightedCellVelAtBottom(iCell)/cellArea
           meridionalAreaWeightedCellVelAtBottom(iCell) = meridionalAreaWeightedCellVelAtBottom(iCell)/cellArea
         end do

         if ( config_time_integrator == trim('split_explicit') .or.  config_time_integrator == trim('split_implicit') &
                                                               .or.  config_time_integrator == trim('split_explicit_ab2') ) then

           do iEdge = 1, nEdges
             normalBarotropicVel(iEdge) = normalBarotropicVelocity(iEdge)
             normalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelocity(1,iEdge)
             normalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelocity(iLevel0250,iEdge)
             normalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelocity(maxLevelEdgeBot(iEdge),iEdge)
           end do

           do iEdge = 1, nEdges
             tangentialBarotropicVel(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAtSFC(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAt250m(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAtBottom(iEdge) = 0.0_RKIND
             do i = 1, nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(i,iEdge)
               tangentialBarotropicVel(iEdge) = tangentialBarotropicVel(iEdge) + weightsOnEdge(i,iEdge) * normalBarotropicVel(eoe)
               tangentialBaroclinicVelAtSFC(iEdge) = tangentialBaroclinicVelAtSFC(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAtSFC(eoe)
               tangentialBaroclinicVelAt250m(iEdge) = tangentialBaroclinicVelAt250m(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAt250m(eoe)
               tangentialBaroclinicVelAtBottom(iEdge) = tangentialBaroclinicVelAtBottom(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAtBottom(eoe)
             end do

             zonalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBarotropicVel(iEdge)*sin(angleEdge(iEdge))
             meridionalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBarotropicVel(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAtSFC(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAtSFC(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAt250m(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAt250m(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAtBottom(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAtBottom(iEdge)*cos(angleEdge(iEdge))
           end do

         endif

         !
         ! compute layer averaged tracers
         !
         activeTracersAvgTopto0100(:,:) = -1.0e34_RKIND
         activeTracersAvg0100to0250(:,:) = -1.0e34_RKIND
         activeTracersAvg0250to0700(:,:) = -1.0e34_RKIND
         activeTracersAvg0700to2000(:,:) = -1.0e34_RKIND
         activeTracersAvg2000toBottom(:,:) = -1.0e34_RKIND
         do iCell = 1, nCells
           sumLayerThickness = layerThickness(1,iCell)
           activeTracersAvgTopto0100(:,iCell) = activeTracers(:,1,iCell)*layerThickness(1,iCell)
           do k=2, min(maxLevelCell(iCell),iLevel0100)
              sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
              activeTracersAvgTopto0100(:,iCell) = activeTracersAvgTopto0100(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
           enddo
           activeTracersAvgTopto0100(:,iCell) = activeTracersAvgTopto0100(:,iCell) / max(sumLayerThickness,1.0_RKIND)

           if (iLevel0100+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0100+1,iCell)
              activeTracersAvg0100to0250(:,iCell) = activeTracers(:,iLevel0100+1,iCell)*layerThickness(iLevel0100,iCell)
              do k=iLevel0100+2, min(maxLevelCell(iCell),iLevel0250)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0100to0250(:,iCell) = activeTracersAvg0100to0250(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0100to0250(:,iCell) = activeTracersAvg0100to0250(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel0250+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0250+1,iCell)
              activeTracersAvg0250to0700(:,iCell) = activeTracers(:,iLevel0250+1,iCell)*layerThickness(iLevel0250+1,iCell)
              do k=iLevel0250+2, min(maxLevelCell(iCell),iLevel0700)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0250to0700(:,iCell) = activeTracersAvg0250to0700(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0250to0700(:,iCell) = activeTracersAvg0250to0700(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel0700+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0700+1,iCell)
              activeTracersAvg0700to2000(:,iCell) = activeTracers(:,iLevel0700+1,iCell)*layerThickness(iLevel0700+1,iCell)
              do k=iLevel0700+2, min(maxLevelCell(iCell),iLevel2000)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0700to2000(:,iCell) = activeTracersAvg0700to2000(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0700to2000(:,iCell) = activeTracersAvg0700to2000(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel2000+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel2000+1,iCell)
              activeTracersAvg2000toBottom(:,iCell) = activeTracers(:,iLevel2000+1,iCell)*layerThickness(iLevel2000+1,iCell)
              do k=iLevel2000+2, maxLevelCell(iCell)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg2000toBottom(:,iCell) = activeTracersAvg2000toBottom(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg2000toBottom(:,iCell) = activeTracersAvg2000toBottom(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif
         enddo

         !
         ! compute barotropic kinetic energy
         !

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            ! normal Barotropic Velocity = sum(h*u)/sum(h) on each edge
            layerThicknessEdge1 = 0.5_RKIND*( layerThickness(1,cell1) + layerThickness(1,cell2) )
            normalThicknessFluxSum(iEdge) = layerThicknessEdge1 * normalVelocity(1,iEdge)
            layerThicknessSumEdge(iEdge) = layerThicknessEdge1

            do k=2, maxLevelEdgeTop(iEdge)
               layerThicknessEdge1 = 0.5_RKIND*( layerThickness(k,cell1) + layerThickness(k,cell2) )

               normalThicknessFluxSum(iEdge) = normalThicknessFluxSum(iEdge) + layerThicknessEdge1 * normalVelocity(k,iEdge)
               layerThicknessSumEdge(iEdge) = layerThicknessSumEdge(iEdge) + layerThicknessEdge1
            enddo
         enddo

         do iCell = 1, nCells
            invAreaCell1 = 1.0_RKIND / areaCell(iCell)
            barotropicSpeed(iCell) = 0.0_RKIND
            !columnIntegratedSpeed(iCell) = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)

               coeff = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell1
               ! this is kinetic energy, in units of m^2/sec^2
               barotropicSpeed(iCell)    = barotropicSpeed(iCell) &
                    + coeff * (normalThicknessFluxSum(iEdge) / layerThicknessSumEdge(iEdge))**2

            end do
            barotropicSpeed(iCell)    = sqrt(2.0_RKIND*barotropicSpeed(iCell))

            ! columnIntegratedSpeed = sum(h*sqrt(2*ke)), where ke is kineticEnergyCell
            !   and the sum is over the full column at cell centers.
            columnIntegratedSpeed(iCell) = layerThickness(1,iCell)*sqrt( 2.0_RKIND * kineticEnergyCell(1,iCell) )
            do k=2, maxLevelCell(iCell)
               columnIntegratedSpeed(iCell) = columnIntegratedSpeed(iCell) &
                    + layerThickness(k,iCell)*sqrt( 2.0_RKIND * kineticEnergyCell(k,iCell) )
            enddo

         end do

         ! scratch arrays
         deallocate(normalThicknessFluxSum)
         deallocate(layerThicknessSumEdge)

         block => block % next
      end do

   end subroutine ocn_compute_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_restart_high_frequency_output
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_high_frequency_output(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_restart_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_finalize_high_frequency_output
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_high_frequency_output(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

   end subroutine ocn_finalize_high_frequency_output!}}}

end module ocn_high_frequency_output

! vim: foldmethod=marker
