! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_ocean_heat_content
!
!> \brief MPAS ocean analysis mode member: Computes Ocean Heat Content at
!>          standard layers
!> \author Luke Van Roekel
!> \date   February 25, 2020
!> \details
!>   Computes ocean heat content in 0-700, 700-2000, 2000-Bottom and 0-Bottom
!>   layers
!>   Returns a nCells array for each layer
!>
!-----------------------------------------------------------------------

module ocn_ocean_heat_content

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

   public :: ocn_init_ocean_heat_content, &
             ocn_compute_ocean_heat_content, &
             ocn_restart_ocean_heat_content, &
             ocn_finalize_ocean_heat_content

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_ocean_heat_content
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    February 25, 2020
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_ocean_heat_content(domain, err)!{{{

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

      type(domain_type), intent(inout) :: domain

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

   end subroutine ocn_init_ocean_heat_content!}}}

!***********************************************************************
!
!  routine ocn_compute_ocean_heat_content
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    February 25, 2020
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_ocean_heat_content(domain, timeLevel, err)!{{{

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

      type(domain_type), intent(inout) :: domain

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

      type(mpas_pool_type), pointer :: oceanHeatContentAMPool
      type(dm_info) :: dminfo
      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: scratchPool
      type(mpas_pool_type), pointer :: diagnosticsPool
      type(mpas_pool_type), pointer :: tracersPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, index_temperature
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:, :), pointer :: zMid, layerThickness
      real(kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, &
                                                 oceanHeatContentSfcToBot, oceanHeatContentSfcTo700m, &
                                                 oceanHeatContent700mTo2000m, oceanHeatContent2000mToBot

      real(kind=RKIND) :: heatLayer
      err = 0

      dminfo = domain%dminfo

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block%structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block%structs, 'layeredOceanHeatContent', oceanHeatContentAMPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block%dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block%dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block%dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

         call mpas_pool_get_array(oceanHeatContentAMPool, 'oceanHeatContentSfcToBot', oceanHeatContentSfcToBot)
         call mpas_pool_get_array(oceanHeatContentAMPool, 'oceanHeatContentSfcTo700m', oceanHeatContentSfcTo700m)
         call mpas_pool_get_array(oceanHeatContentAMPool, 'oceanHeatContent700mTo2000m', oceanHeatContent700mTo2000m)
         call mpas_pool_get_array(oceanHeatContentAMPool, 'oceanHeatContent2000mToBot', oceanHeatContent2000mToBot)

         do iCell = 1, nCellsSolve
            k = 1
            oceanHeatContentSfcTo700m(iCell) = 0.0_RKIND
            oceanHeatContentSfcToBot(iCell) = 0.0_RKIND
            oceanHeatContent700mTo2000m(iCell) = 0.0_RKIND
            oceanHeatContent2000mToBot(iCell) = 0.0_RKIND
            do while (zMid(k, iCell) .ge. -700.0_RKIND .and. k .le. maxLevelCell(iCell))
               heatLayer = layerThickness(k, iCell)*areaCell(iCell)*cp_sw*rho_sw* &
                           (tracers(index_temperature, k, iCell) + 273.15_RKIND)
               oceanHeatContentSfcTo700m(iCell) = oceanHeatContentSfcTo700m(iCell) + heatLayer
               oceanHeatContentSfcToBot(iCell) = oceanHeatContentSfcToBot(iCell) + heatLayer
               k = k + 1
            enddo

            do while (zMid(k, iCell) .ge. -2000.0_RKIND .and. k .le. maxLevelCell(iCell))
               heatLayer = layerThickness(k, iCell)*areaCell(iCell)*cp_sw*rho_sw* &
                           (tracers(index_temperature, k, iCell) + 273.15_RKIND)
               oceanHeatContent700mTo2000m(iCell) = oceanHeatContent700mTo2000m(iCell) + heatLayer
               oceanHeatContentSfcToBot(iCell) = oceanHeatContentSfcToBot(iCell) + heatLayer
               k = k + 1
            enddo

            do while (k .le. maxLevelCell(iCell))
               heatLayer = layerThickness(k, iCell)*areaCell(iCell)*cp_sw*rho_sw* &
                           (tracers(index_temperature, k, iCell) + 273.15_RKIND)
               oceanHeatContent2000mToBot(iCell) = oceanHeatContent2000mToBot(iCell) + heatLayer
               oceanHeatContentSfcToBot(iCell) = oceanHeatContentSfcToBot(iCell) + heatLayer
               k = k + 1
            end do
         end do

         block => block%next
      end do

   end subroutine ocn_compute_ocean_heat_content!}}}

!***********************************************************************
!
!  routine ocn_restart_ocean_heat_content
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    February 25, 2020
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_ocean_heat_content(domain, err)!{{{

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

      type(domain_type), intent(inout) :: domain

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

   end subroutine ocn_restart_ocean_heat_content!}}}

!***********************************************************************
!
!  routine ocn_finalize_ocean_heat_content
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    February 25, 2020
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_ocean_heat_content(domain, err)!{{{

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

      type(domain_type), intent(inout) :: domain

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

   end subroutine ocn_finalize_ocean_heat_content!}}}

end module ocn_ocean_heat_content

! vim: foldmethod=marker
