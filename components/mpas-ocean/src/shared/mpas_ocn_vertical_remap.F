! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_remap
!
!> \brief MPAS ocean vertical Lagrangian remapping
!> \author Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis
!> \date   July 2021
!> \details
!>  This module contains the vertical remapping routine.
!
!-----------------------------------------------------------------------

module ocn_vertical_remap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_config

   use ppr_1d

   ! TODO add use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   
   type(rmap_opts) :: opts
   
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_remap_vert_state
   public :: ocn_vertical_remap_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_remap_vert_state
!
!> \brief MPAS ocean vertical Lagrangian remapping
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine remaps state variables from current layerThickness to 
!>  a new layerThickness
!
!-----------------------------------------------------------------------

   subroutine ocn_remap_vert_state(block, layerThicknessNew, err)

      type (block_type), intent(in) :: block
      real (kind=RKIND), dimension(:, :), intent(in) :: layerThicknessNew
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges, nVertLevels
      integer :: nLevels, nVars, nDoFs, nTracers
      integer :: iCell, jCell, iEdge, iTrac, k

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_iterator_type) :: groupItr

      type(rmap_work) :: work
      type(rcon_ends), dimension(:), allocatable :: bcUpper, bcLower

      real (kind=RKIND), dimension(:, :, :), allocatable :: dataNow, dataNew
      real (kind=RKIND), dimension(:, :), allocatable :: heightCellNow, heightEdgeNow
      real (kind=RKIND), dimension(:, :), allocatable :: heightCellNew, heightEdgeNew
      
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot
      integer, dimension(:, :), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocity
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThickness,lowFreqDivergence
      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroup

      err = 0

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
      
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 2)
      if (config_use_freq_filtered_thickness) then
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence, 2)
      end if
      
   !----------------------------------------------------------------------------
   !- SETUP VERTICAL GRIDS:
   !- heightCellNow, heightEdgeNow: depth at layer interfaces (lagrangian grid)
   !- heightCellNew, heighyEdgeNew: depth at layer interfaces (new target grid)      

   !- this should all go into a separate "regrid" routine
   !- (eventually) this would handle regridding to "other" coordinates, rather
   !- than just the z* approach here.

      allocate(heightCellNow(nVertLevels + 1, nCells))
      allocate(heightCellNew(nVertLevels + 1, nCells))
      allocate(heightEdgeNow(nVertLevels + 1, nEdges))
      allocate(heightEdgeNew(nVertLevels + 1, nEdges))
      
      do iCell = 1, nCells
         heightCellNow(:, iCell) = -bottomDepth(iCell)
         heightCellNew(:, iCell) = -bottomDepth(iCell)

      !- reconstruct "now" heights

         do k = maxLevelCell(iCell), 1, -1
            heightCellNow(k, iCell) = &
               heightCellNow(k + 1, iCell) + layerThickness(k, iCell)
         end do

      !- reconstruct "new" heights
      !- Assign new layerThicknesses

         heightCellNew(1, iCell) = heightCellNow(1, iCell)
         do k = maxLevelCell(iCell), 2, -1
            heightCellNew(k, iCell) = &
               heightCellNew(k + 1, iCell) + layerThicknessNew(k, iCell)
         end do
      end do
      layerThickness = layerThicknessNew

   !- There's currently a silly issue in PPR: it assumes the "coordinates" for
   !- the data being remapped increases with index. Here, we have the opposite 
   !- (z becomes more -ve with k)... 
   !- I'll make some updates to PPR to handle this, but for now, just flip the
   !- signs to make PPR happy. This shouldn't have any actual impact.
 
      heightCellNow = -1.0_RKIND * heightCellNow
      heightCellNew = -1.0_RKIND * heightCellNew


      do iEdge = 1, nEdges
         iCell = cellsOnEdge(1, iEdge)
         jCell = cellsOnEdge(2, iEdge)

         heightEdgeNow(:, iEdge) = -bottomDepth(iCell)
         heightEdgeNew(:, iEdge) = -bottomDepth(iCell)

      !- this is just a quick hack for the edge heights
      !- really should use layerThicknessEdge here, but it's unclear if this has
      !- been properly initiallised at this point...
      !- layerThicknessEdge could involve upwinding, so need to be careful!

         do k = maxLevelEdgeTop(iEdge), 1, -1
            heightEdgeNow(k, iEdge) = &
                           0.5_RKIND * (heightCellNow(k, iCell) + heightCellNow(k, jCell))
            heightEdgeNew(k, iEdge) = &
                           0.5_RKIND * (heightCellNew(k, iCell) + heightCellNew(k, jCell))
         end do
      end do

   !----------------------------------------------------------------------------
   !- ACTUAL REMAPPING FROM HERE

   !- just simple hard-coded options for PPR at this stage...

      if (associated(highFreqThickness) .and. associated(lowFreqDivergence)) then
         nVars = 3 ! how many tracers to be remapped concurrently per column
      else
         nVars = 1 ! how many tracers to be remapped concurrently per column
      endif
      nDoFs = 1 ! always 1 for MPAS, it'd be >1 for DG methods (with more DoF per layer)

      allocate(dataNow(nDoFs, nVars, nVertLevels))
      allocate(dataNew(nDoFs, nVars, nVertLevels))
      allocate(bcUpper(nVars))
      allocate(bcLower(nVars))
      
      dataNow = 0.0_RKIND
      dataNew = 0.0_RKIND

      bcUpper%bcopt = bcon_loose ! "loose" => extrapolate
      bcLower%bcopt = bcon_loose

      call work%init(nVertLevels + 1, nVars, opts) ! this is internal workspace for PPR

   !- Remap all edge-centred variables from heightEdgeNow to heightEdgeNew
      if (associated(highFreqThickness) .and. associated(lowFreqDivergence)) then
         do iEdge = 1, nEdges
            nLevels = maxLevelEdgeTop(iEdge) + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, :) = normalVelocity(:, iEdge)
            dataNow(1, 2, :) = highFreqThickness(:, iEdge)
            dataNow(1, 3, :) = lowFreqDivergence(:, iEdge)

            call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
               heightEdgeNow(:, iEdge), &
               heightEdgeNew(:, iEdge), &
               dataNow, dataNew, bcUpper, bcLower, work, opts)

            normalVelocity   (1:maxLevelEdgeTop(iEdge), iEdge) = dataNew(1, 1, 1:maxLevelEdgeTop(iEdge))
            highFreqThickness(1:maxLevelEdgeTop(iEdge), iEdge) = dataNew(1, 2, 1:maxLevelEdgeTop(iEdge))
            lowFreqDivergence(1:maxLevelEdgeTop(iEdge), iEdge) = dataNew(1, 3, 1:maxLevelEdgeTop(iEdge))
         
         end do
      else
         do iEdge = 1, nEdges
            nLevels = maxLevelEdgeTop(iEdge) + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, :) = normalVelocity(:, iEdge)

            call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
               heightEdgeNow(:, iEdge), &
               heightEdgeNew(:, iEdge), &
               dataNow, dataNew, bcUpper, bcLower, work, opts)

            normalVelocity(1:maxLevelEdgeTop(iEdge), iEdge) = dataNew(1, 1, 1:maxLevelEdgeTop(iEdge))
         
         end do
      end if
      deallocate(dataNow,dataNew)
      deallocate(bcUpper,bcLower)
      call work%free()
   
   !- Remap all cell-centred variables from heightCellNow to heightCellNew

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 2)
            if ( associated(tracersGroup) ) then
               
            !- we can remap all of the tracers in one go via PPR, but for now
            !- just a simple loop over each one at a time for clarity...

               nTracers = size(tracersGroup, dim=1)
               
               allocate(dataNow(nDoFs, nTracers, nVertLevels))
               allocate(dataNew(nDoFs, nTracers, nVertLevels))
               
               call work%init(nVertLevels + 1, nTracers, opts) ! this is internal workspace for PPR

               allocate(bcUpper(nTracers))
               allocate(bcLower(nTracers))
               ! TODO replace hardcoded bcon_loose with config option 
               bcUpper%bcopt = bcon_loose ! "loose" => extrapolate
               bcLower%bcopt = bcon_loose
               
               do iCell = 1, nCells
                  do iTrac = 1, nTracers
                     nLevels = maxLevelCell(iCell) + 1

                     if (nLevels .lt. 2) cycle

                     dataNow(1, iTrac, :) = tracersGroup(iTrac, :, iCell)
                  end do
                  call rmap1d(nLevels, nLevels, nTracers,nDoFs, &
                     heightCellNow(:, iCell), &
                     heightCellNew(:, iCell), &
                     dataNow, dataNew, bcUpper, bcLower, work, opts)

                  do iTrac = 1, nTracers
                     tracersGroup(iTrac, 1:maxLevelCell(iCell), iCell) = dataNew(1, iTrac, 1:maxLevelCell(iCell))
                  end do
               
               end do

               deallocate(dataNow,dataNew)
               deallocate(bcUpper,bcLower)
               call work%free()
               
            end if
         end if
      end do

   !- dealloc any local data, etc      

      deallocate(heightCellNow)
      deallocate(heightEdgeNow)
      deallocate(heightCellNew)
      deallocate(heightEdgeNew)

   end subroutine ocn_remap_vert_state
  
   subroutine ocn_vertical_remap_init(err)
      
      integer, intent(out) :: err !< Output: Error flag
      
      err = 0

      if ( config_time_integrator == 'RK4' ) then
         CALL mpas_log_write('Vertical remap not supported for RK4', &
            MPAS_LOG_CRIT)
      endif

      ! Options for remapping
      if ( config_vert_remap_order == 1 ) then
         opts%cell_meth = pcm_method ! PCM method
      elseif ( config_vert_remap_order == 2 ) then
         opts%cell_meth = plm_method ! PLM method
      elseif ( config_vert_remap_order == 3 ) then
         opts%cell_meth = ppm_method ! PPM method
         opts%edge_meth = p3e_method ! 3rd-order edge interp.
      elseif ( config_vert_remap_order == 5 ) then
         opts%cell_meth = pqm_method ! PPM method
         opts%edge_meth = p5e_method ! 3rd-order edge interp.
      endif
      opts%cell_lims = mono_limit ! monotone slope limits TODO: consider config_monotonic

   end subroutine ocn_vertical_remap_init

end module ocn_vertical_remap
! vim: foldmethod=marker
