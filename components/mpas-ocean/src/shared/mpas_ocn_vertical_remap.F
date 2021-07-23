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
   use ocn_mesh

   use ppr_1d

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   
   type(rmap_opts) :: opts
   integer :: bc_upper, bc_lower
   
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

      integer :: nCells, nEdges
      integer :: nLayers, nLevels, nVars, nDoFs, nTracers
      integer :: iCell, jCell, iEdge, iTrac, k, kTop, kBot

      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_iterator_type) :: groupItr

      type(rmap_work) :: work
      type(rcon_ends), dimension(:), allocatable :: bcUpper, bcLower

      real (kind=RKIND), dimension(:, :, :), allocatable :: dataNow, dataNew
      real (kind=RKIND), dimension(:, :), allocatable :: heightCellNow, heightEdgeNow
      real (kind=RKIND), dimension(:, :), allocatable :: heightCellNew, heightEdgeNew
      
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocity
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThickness,lowFreqDivergence
      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroup

      err = 0

      itimestepLastRemap = itimestepLastRemap + 1
      if (itimestepLastRemap < config_vert_remap_interval) return
      itimestepLastRemap = 0

      ! Remapping currently only supports one block for the whole domain rather
      ! than sub-blocks
      
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
      
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 2)
      if (config_use_freq_filtered_thickness) then
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence, 2)
      end if
      
      nCells = nCellsAll
      nEdges = nEdgesAll

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
      
   !- There's currently a silly issue in PPR: it assumes the "coordinates" for
   !- the data being remapped increases with index. Here, we have the opposite 
   !- (z becomes more -ve with k)... 
   !- I'll make some updates to PPR to handle this, but for now, just flip the
   !- signs to make PPR happy. This shouldn't have any actual impact.
      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         heightCellNow(:, iCell) = bottomDepth(iCell)
         heightCellNew(:, iCell) = bottomDepth(iCell)

      !- reconstruct "now" heights

         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            heightCellNow(k, iCell) = &
               heightCellNow(k + 1, iCell) - layerThickness(k, iCell)
         end do

      !- reconstruct "new" heights
      !- Assign new layerThicknesses

         heightCellNew(minLevelCell(iCell), iCell) = heightCellNow(minLevelCell(iCell), iCell)
         do k = maxLevelCell(iCell), minLevelCell(iCell)+1, -1
            heightCellNew(k, iCell) = &
               heightCellNew(k + 1, iCell) - layerThicknessNew(k, iCell)
         end do
      end do
      !$omp end do
      !$omp end parallel
      layerThickness = layerThicknessNew

      !$omp parallel
      !$omp do schedule(runtime) private(k, iCell, jCell)
      do iEdge = 1, nEdges
         iCell = cellsOnEdge(1, iEdge)
         jCell = cellsOnEdge(2, iEdge)

         heightEdgeNow(:, iEdge) = bottomDepth(iCell)
         heightEdgeNew(:, iEdge) = bottomDepth(iCell)

      !- this is just a quick hack for the edge heights
      !- really should use layerThicknessEdge here, but it's unclear if this has
      !- been properly initiallised at this point...
      !- layerThicknessEdge could involve upwinding, so need to be careful!

         do k = maxLevelEdgeTop(iEdge), minLevelEdgeBot(iEdge), -1
            heightEdgeNow(k, iEdge) = &
                           0.5_RKIND * (heightCellNow(k, iCell) + heightCellNow(k, jCell))
            heightEdgeNew(k, iEdge) = &
                           0.5_RKIND * (heightCellNew(k, iCell) + heightCellNew(k, jCell))
         end do
      end do
      !$omp end do
      !$omp end parallel

      !----------------------------------------------------------------------------
      ! ACTUAL REMAPPING FROM HERE

      nDoFs = 1 ! always 1 for MPAS, it'd be >1 for DG methods (with more DoF per layer)

      ! how many variables to be remapped concurrently per column
      if (associated(highFreqThickness) .and. associated(lowFreqDivergence)) then
         nVars = 3
      else
         nVars = 1
      endif
      
      allocate(dataNow(nDoFs, nVars, nVertLevels+1))
      allocate(dataNew(nDoFs, nVars, nVertLevels+1))
      allocate(bcUpper(nVars))
      allocate(bcLower(nVars))
      
      ! set boundary conditions
      bcUpper%bcopt = bc_upper
      bcLower%bcopt = bc_lower

      dataNow = 0.0_RKIND
      dataNew = 0.0_RKIND

      call work%init(nVertLevels + 1, nVars, opts) ! this is internal workspace for PPR

      ! Remap all edge-centred variables from heightEdgeNow to heightEdgeNew
      if (associated(highFreqThickness) .and. associated(lowFreqDivergence)) then
         !$omp parallel
         !$omp do schedule(runtime) private(iEdge, kTop, kBot, nLayers, nLevels, dataNow, dataNew)
         do iEdge = 1, nEdges

            kTop = minLevelEdgeBot(iEdge)
            kBot = maxLevelEdgeTop(iEdge)

            nLayers = kBot - kTop + 1
            nLevels = nLayers + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, 1:nLayers) = normalVelocity   (kTop:kBot, iEdge)
            dataNow(1, 2, 1:nLayers) = highFreqThickness(kTop:kBot, iEdge)
            dataNow(1, 3, 1:nLayers) = lowFreqDivergence(kTop:kBot, iEdge)

            call mpas_timer_start("ppr")
            call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
               heightEdgeNow(kTop:kBot+1, iEdge), &
               heightEdgeNew(kTop:kBot+1, iEdge), &
               dataNow, dataNew, &
               bcUpper, bcLower, work, opts)
            call mpas_timer_stop("ppr")

            normalVelocity   (kTop:kBot, iEdge) = dataNew(1, 1, 1:nLayers)
            highFreqThickness(kTop:kBot, iEdge) = dataNew(1, 2, 1:nLayers)
            lowFreqDivergence(kTop:kBot, iEdge) = dataNew(1, 3, 1:nLayers)
         end do
         !$omp end do
         !$omp end parallel
      else
         !$omp parallel
         !$omp do schedule(runtime) private(iEdge, kTop, kBot, nLayers, nLevels, dataNow, dataNew)
         do iEdge = 1, nEdges

            kTop = minLevelEdgeBot(iEdge)
            kBot = maxLevelEdgeTop(iEdge)

            nLayers = kBot - kTop + 1
            nLevels = nLayers + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, 1:nLayers) = normalVelocity(kTop:kBot, iEdge)

            call mpas_timer_start("ppr")
            call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
               heightEdgeNow(kTop:kBot+1, iEdge), &
               heightEdgeNew(kTop:kBot+1, iEdge), &
               dataNow, dataNew, &
               bcUpper, bcLower, work, opts)
            call mpas_timer_stop("ppr")

            normalVelocity   (kTop:kBot, iEdge) = dataNew(1, 1, 1:nLayers)
         
         end do
         !$omp end do
         !$omp end parallel
      end if
      deallocate(dataNow,dataNew)
      deallocate(bcUpper,bcLower)
      call work%free()
   
      ! Remap all cell-centred variables from heightCellNow to heightCellNew

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 2)
            if ( associated(tracersGroup) ) then
               
               nTracers = size(tracersGroup, dim=1)
               
               allocate(dataNow(nDoFs, nTracers, nVertLevels))
               allocate(dataNew(nDoFs, nTracers, nVertLevels))
               
               call work%init(nVertLevels + 1, nTracers, opts) ! this is internal workspace for PPR

               allocate(bcUpper(nTracers))
               allocate(bcLower(nTracers))
               bcUpper%bcopt = bc_upper
               bcLower%bcopt = bc_lower
               
               !$omp parallel
               !$omp do schedule(runtime) private(iCell, iTrac, kTop, kBot, nLayers, nLevels, dataNow, dataNew)
               do iCell = 1, nCells

                  kTop = minLevelCell(iCell)
                  kBot = maxLevelCell(iCell)

                  nLayers = kBot - kTop + 1
                  nLevels = nLayers + 1

                  if (nLevels .lt. 2) cycle

                  do iTrac = 1, nTracers

                     dataNow(1, iTrac, 1:nLayers) = tracersGroup(iTrac, kTop:kBot, iCell)

                  end do
                  call mpas_timer_start("ppr")
                  call rmap1d(nLevels, nLevels, nTracers,nDoFs, &
                     heightCellNow(kTop:kBot+1, iCell), &
                     heightCellNew(kTop:kBot+1, iCell), &
                     dataNow, dataNew, &
                     bcUpper, bcLower, work, opts)
                  call mpas_timer_stop("ppr")

                  do iTrac = 1, nTracers
                     tracersGroup(iTrac, kTop:kBot, iCell) = dataNew(1, iTrac, 1:nLayers)
                  end do
                     
               end do
               !$omp end do
               !$omp end parallel

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
  
!***********************************************************************
!
!  routine ocn_vertical_remap_init
!
!> \brief   Initializes ocean vertical remapping
!> \author  Carolyn Begeman
!> \date    July 2021
!> \details
!>  This routine initializes parameters required for vertical Lagrangian
!>  remapping
!
!-----------------------------------------------------------------------

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
      ! TODO replace hardcoded bcon_loose with config option 
      !select case (config_vert_remap_bc_upper)
      !case ('extrapolate')
      bc_upper = bcon_loose
      !case ('dirichlet')
      !bc_upper = bcon_value
      bc_lower = bcon_loose
      !bc_lower = bcon_value

   end subroutine ocn_vertical_remap_init

end module ocn_vertical_remap
! vim: foldmethod=marker
