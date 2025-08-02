! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_li_flux_gates
!
!> \brief MPAS land ice analysis mode member: mpas_li_flux_fates
!> \author Trevor Hillebrand
!> \date   03-28-2025
!> \details Calculates ice flux through gates defined in input file.
!>
!>
!-----------------------------------------------------------------------
module li_flux_gates

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use li_mask
   use li_constants

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

   public :: li_compute_flux_across_gates

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
!  routine li_compute_flux_across_gates
!
!> \brief   Compute MPAS-Land Ice analysis member
!> \author  Trevor Hillebrand
!> \date    03/28/2025
!> \details
!>  This routine conducts all computation required for flux through specified gates
!
!-----------------------------------------------------------------------

   subroutine li_compute_flux_across_gates(domain, memberName, timeLevel, nFluxGates, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel
      integer, intent(in) :: nFluxGates  ! need as an input instead of pointer for global reduce
      character (len=*), intent(in) :: memberName

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

      type (mpas_pool_type), pointer :: fluxGatesAMPool
      type (dm_info) :: dminfo
      type (mpas_pool_type), pointer :: fluxGatesPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: fluxGatesAM
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool

      real (kind=RKIND), dimension(:), pointer ::  dvEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, &
                                                    layerNormalVelocity
      integer, dimension(:,:), pointer :: fluxGateEdgeMasks, fluxGateEdgeMasksSigns
      real (kind=RKIND), dimension(:), pointer :: iceFluxThroughGates
      integer, pointer :: nEdgesSolve
      real (kind=RKIND), pointer ::  rhoi
      real (kind=RKIND), dimension(nFluxGates) :: localIceFluxThroughGates
      integer :: iGate

      dminfo = domain % dminfo
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'fluxGatesAM', fluxGatesAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'fluxGates', fluxGatesPool)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(fluxGatesAMPool, 'iceFluxThroughGates', iceFluxThroughGates)
      call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)
      call mpas_pool_get_array(fluxGatesPool, 'fluxGateEdgeMasks', fluxGateEdgeMasks)
      call mpas_pool_get_array(fluxGatesPool, 'fluxGateEdgeMasksSigns', fluxGateEdgeMasksSigns)

      localIceFluxThroughGates(:) = 0.0_RKIND
      iceFluxThroughGates(:) = 0.0_RKIND

      do iGate = 1, nFluxGates
         localIceFluxThroughGates(iGate) = sum( sum(layerNormalVelocity(:,1:nEdgesSolve) * layerThicknessEdge(:,1:nEdgesSolve), 1) * &
                                                dvEdge(1:nEdgesSolve) * rhoi * scyr * fluxGateEdgeMasks(iGate, 1:nEdgesSolve) * &
                                                fluxGateEdgeMasksSigns(iGate, 1:nEdgesSolve) )
      enddo

      call mpas_dmpar_sum_real_array(dminfo, nFluxGates, localIceFluxThroughGates, iceFluxThroughGates)

   end subroutine li_compute_flux_across_gates

end module li_flux_gates
