! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_high_frequency_output
!
!> \brief MPAS ocean analysis mode member: high_frequency_output
!> \author Todd Ringler
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
   use ocn_diagnostics_routines

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

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_high_frequency_output
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
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
!> \author  Todd Ringler
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
      type (mpas_pool_type), pointer :: scratchPool

      integer :: iLevel, iLevelTarget, iCell, iEdge, i, cell1, cell2, k
      integer, pointer :: nVertLevels, nCells, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge

      real (kind=RKIND) :: invAreaCell1, layerThicknessEdge1, coeff
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, kineticEnergyAt100m, kineticEnergyAtSurface, relativeVorticityAt100m
      real (kind=RKIND), dimension(:), pointer :: divergenceAt100m, relativeVorticityVertexAt100m
      real (kind=RKIND), dimension(:), pointer :: barotropicSpeed, columnIntegratedSpeed, dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell, relativeVorticityCell, activeTracersAtSurface
      real (kind=RKIND), dimension(:,:), pointer :: relativeVorticity, divergence, layerThickness, normalVelocity
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      ! scratch space
      type (field1DReal), pointer :: normalThicknessFluxSumField, layerThicknessSumEdgeField
      real (kind=RKIND), dimension(:), pointer :: normalThicknessFluxSum, layerThicknessSumEdge

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
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
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAt100m', kineticEnergyAt100m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAtSurface', kineticEnergyAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAt100m', relativeVorticityAt100m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAtSurface', activeTracersAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'divergenceAt100m', divergenceAt100m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityVertexAt100m', relativeVorticityVertexAt100m)

         !
         ! note for future build out
         ! HERE: interpolate vertically to target z-level or isopycnal or some other surface
         !

         ! for now, just get close enough
         iLevelTarget = 1
         do iLevel=2,nVertLevels
           if(refBottomDepth(iLevel) > 100.0_RKIND) then
              iLevelTarget = iLevel-1
              exit
           endif
         enddo

         ! copy data into high frequency output fields
         ! tracer data will be converted to new tracer infrastrcture (and this line removed) before June 23 2015.
         kineticEnergyAt100m(:) = kineticEnergyCell(iLevelTarget,:)
         relativeVorticityAt100m(:) = relativeVorticityCell(iLevelTarget,:)
         divergenceAt100m(:) = divergence(iLevelTarget,:)
         relativeVorticityVertexAt100m = relativeVorticity(iLevelTarget,:)
         kineticEnergyAtSurface(:) = kineticEnergyCell(1,:)
         activeTracersAtSurface(1,:) = activeTracers(1,1,:)
         activeTracersAtSurface(2,:) = activeTracers(2,1,:)

         !
         ! compute barotropic kinetic energy
         !
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'barotropicSpeed', barotropicSpeed)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'columnIntegratedSpeed', columnIntegratedSpeed)

         ! scratch arrays
         call mpas_pool_get_field(scratchPool, 'normalThicknessFluxSum', normalThicknessFluxSumField)
         call mpas_allocate_scratch_field(normalThicknessFluxSumField, .true.)
         normalThicknessFluxSum => normalThicknessFluxSumField % array

         call mpas_pool_get_field(scratchPool, 'layerThicknessSumEdge', layerThicknessSumEdgeField)
         call mpas_allocate_scratch_field(layerThicknessSumEdgeField, .true.)
         layerThicknessSumEdge => layerThicknessSumEdgeField % array

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

         block => block % next
      end do

   end subroutine ocn_compute_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_restart_high_frequency_output
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
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
!> \author  Todd Ringler
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
