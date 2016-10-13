! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_transect_transport
!
!> \brief MPAS ocean analysis mode member: transect_transport
!> \author Mark Petersen
!> \date   April 2016
!> \details
!>  MPAS ocean analysis mode member: transect_transport
!>
!-----------------------------------------------------------------------

module ocn_transect_transport

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

   public :: ocn_init_transect_transport, &
             ocn_compute_transect_transport, &
             ocn_restart_transect_transport, &
             ocn_finalize_transect_transport

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_transect_transport
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    April 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_transect_transport(domain, err)!{{{

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

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: transectTransportAMPool
      type (mpas_pool_type), pointer :: meshPool

      type (mpas_pool_type), pointer :: transectPool

      type (mpas_pool_type), pointer :: transectTransportAM

      integer, pointer :: nTransects
      integer :: iTransect

      integer, dimension(:), pointer :: transectEdgeMasksMax
      integer, dimension(:,:), pointer :: transectEdgeMasks

      err = 0

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransects', nTransects)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'transectTransportAM', transectTransportAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'transects', transectPool)

      call mpas_pool_get_array(transectTransportAMPool,'transectEdgeMasksMax',transectEdgeMasksMax)

      transectEdgeMasksMax(:) = 0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'transectTransportAM', transectTransportAMPool)

         call mpas_pool_get_array(transectPool,'transectEdgeMasks',transectEdgeMasks)
         do iTransect = 1,nTransects
            !currentTransect = transectsInGroup(iTransect, transectGroupNumber)
            transectEdgeMasksMax(iTransect) = max(transectEdgeMasksMax(iTransect),maxval(transectEdgeMasks(iTransect,:)))
         end do

         block => block % next
      end do

   end subroutine ocn_init_transect_transport!}}}

!***********************************************************************
!
!  routine ocn_compute_transect_transport
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    April 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_transect_transport(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: transectTransportAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: transectTransportAM

      type (mpas_pool_type), pointer :: transectPool

      integer :: currentTransect, transectGroupNumber, transectsInAddGroup, i
      integer, dimension(:,:), pointer :: transectsInGroup
      character (len=STRKIND), dimension(:), pointer :: transectNames, transectGroupNames
      integer, dimension(:), pointer ::  nTransectsInGroup
      integer, pointer :: nTransectGroups, maxTransectsInGroup
      character (len=STRKIND), pointer :: additionalTransect

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nEdgesSolve, num_tracers, nTransects
      integer :: k, iEdge, iTransect, nTransportVariables, c1,c2
      integer, dimension(:), pointer :: maxLevelEdgeTop, transectEdgeMasksMax
      integer, dimension(:,:), pointer :: transectEdgeMasks, transectEdgeMaskSigns, cellsOnEdge

      real (kind=RKIND) :: m3ps_to_Sv
      real (kind=RKIND), dimension(:), pointer ::  dvEdge, transectVolumeTransport,refLayerThickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, transectVolumeTransportZ
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      real (kind=RKIND), dimension(:,:,:), allocatable ::  sumTransport, totalSumTransport

      err = 0

      dminfo = domain % dminfo

      ! Only computing volume transport right now.  Could add heat, tracer transport later.
      nTransportVariables = 1

      m3ps_to_Sv = 1e-6; ! m^3/sec flux to Sverdrups

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransects',nTransects)

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransectGroups', nTransectGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxTransectsInGroup', maxTransectsInGroup)
      call mpas_pool_get_config(domain % configs, 'config_AM_transectTransport_transect_group', &
                                additionalTransect)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'transects', transectPool)
      call mpas_pool_get_array(transectPool, 'transectsInGroup', transectsInGroup)
      call mpas_pool_get_array(transectPool, 'nTransectsInGroup', nTransectsInGroup)
      call mpas_pool_get_array(transectPool, 'transectNames', transectNames)
      call mpas_pool_get_array(transectPool, 'transectGroupNames', transectGroupNames)
      do i = 1, nTransectGroups
         if (transectGroupNames(i) .eq. additionalTransect) then
            transectGroupNumber = i
         end if
      end do

      transectsInAddGroup = nTransectsInGroup(transectGroupNumber)

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      allocate(sumTransport(nVertLevels,nTransects,nTransportVariables))
      allocate(totalSumTransport(nVertLevels,nTransects,nTransportVariables))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         !call mpas_pool_get_subpool(block % structs, 'tracersPool', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'transectTransportAM', transectTransportAMPool)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'transects', transectPool)

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nTransects', nTransects)

         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness',refLayerThickness);

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         !call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

         call mpas_pool_get_array(transectPool,'transectEdgeMasks',transectEdgeMasks)
         call mpas_pool_get_array(transectPool,'transectEdgeMaskSigns',transectEdgeMaskSigns)

         call mpas_pool_get_array(transectTransportAMPool,'transectEdgeMasksMax',transectEdgeMasksMax)

         sumTransport = 0.0_RKIND

         do iTransect = 1,transectsInAddGroup
            currentTransect = transectsInGroup(iTransect, transectGroupNumber)
            if (transectEdgeMasksMax(currentTransect)==0) cycle
            do iEdge = 1,nEdgesSolve
               if (transectEdgeMasks(currentTransect,iEdge)==0) cycle
               c1 = cellsOnEdge(1,iEdge)
               c2 = cellsOnEdge(2,iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  sumTransport(k,iTransect,1) = sumTransport(k,iTransect,1) + &
                       transectEdgeMaskSigns(currentTransect,iEdge) &
                       * normalVelocity(k,iEdge)*dvEdge(iEdge) &
                       * 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))*m3ps_to_Sv
               end do
            end do
         end do

         block => block % next
      end do

      ! mpi summation over all processors
      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nTransects*nTransportVariables, &
           sumTransport, totalSumTransport)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'transectTransportAM', transectTransportAMPool)
      call mpas_pool_get_array(transectTransportAMPool,'transectVolumeTransport',transectVolumeTransport)
      call mpas_pool_get_array(transectTransportAMPool,'transectVolumeTransportZ',transectVolumeTransportZ)

      do iTransect = 1,nTransects
         transectVolumeTransportZ(:,iTransect) = totalSumTransport(:,iTransect,1)
         transectVolumeTransport(iTransect) = sum(transectVolumeTransportZ(:,iTransect))
      end do

      deallocate(sumTransport, totalSumTransport)

   end subroutine ocn_compute_transect_transport!}}}

!***********************************************************************
!
!  routine ocn_restart_transect_transport
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    April 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_transect_transport(domain, err)!{{{

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

   end subroutine ocn_restart_transect_transport!}}}

!***********************************************************************
!
!  routine ocn_finalize_transect_transport
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    April 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_transect_transport(domain, err)!{{{

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

   end subroutine ocn_finalize_transect_transport!}}}

end module ocn_transect_transport

! vim: foldmethod=marker
