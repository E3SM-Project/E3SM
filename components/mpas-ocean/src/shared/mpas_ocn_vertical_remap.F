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
   use ocn_diagnostics_variables
   use ocn_mesh
   use ocn_vertical_regrid

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
   integer :: itimestepLastRemap

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
!> \brief  MPAS ocean vertical Lagrangian remapping
!> \author Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis
!> \date   July 2021
!> \details
!>  This routine remaps state variables from current layerThickness to
!>  a new layerThickness
!
!-----------------------------------------------------------------------

   subroutine ocn_remap_vert_state(block, err)

      type (block_type), intent(in) :: block
      integer, intent(out) :: err

      integer :: nCells, nEdges
      integer :: nLayers, nLevels, nVars, nDoFs, nTracers
      integer :: cell1, cell2, iCell, jCell, iEdge, iTrac
      integer :: k, kmin, kmax, kTop, kBot

      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_iterator_type) :: groupItr

      type(rmap_work) :: work ! Internal workspace for PPR
      type(rcon_ends), dimension(:), allocatable :: bcUpper, bcLower

      real (kind=RKIND), dimension(:, :, :), allocatable :: dataNow, dataNew
      real (kind=RKIND), dimension(:, :), allocatable :: &
         layerThicknessNew,            & ! layerThickness (new target grid)
         layerThickEdgeNew,            & ! layerThickness at edges (new target grid)
         heightCellNow, heightEdgeNow, & ! depth at layer interfaces (lagrangian grid)
         heightCellNew, heightEdgeNew    ! depth at layer interfaces (new target grid)

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

      ! SETUP VERTICAL GRIDS:

      allocate(heightCellNow(nVertLevels + 1, nCells))
      allocate(heightCellNew(nVertLevels + 1, nCells))
      allocate(heightEdgeNow(nVertLevels + 1, nEdges))
      allocate(heightEdgeNew(nVertLevels + 1, nEdges))
      allocate(layerThicknessNew(nVertLevels, nCells + 1))
      allocate(layerThickEdgeNew(nVertLevels, nEdges))

      ! Compute new layer thicknesses based on vertical coordinate choice
      ! Currently, zstar coordinate is hard-coded
      call ocn_vert_regrid(restingThickness, layerThickness, &
                           layerThicknessNew, err)

      ! Compute the layer interface locations between bottomDepth and -ssh
      ! (rather than between -bottomDepth and ssh due to PPR limitations)

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         heightCellNow(:, iCell) = bottomDepth(iCell)
         heightCellNew(:, iCell) = bottomDepth(iCell)

         ! reconstruct "now" heights
         do k = maxLevelCell(iCell), minLevelCell(iCell), -1
            heightCellNow(k, iCell) = &
               heightCellNow(k + 1, iCell) - layerThickness(k, iCell)
         end do

         ! reconstruct "new" heights
         heightCellNew(minLevelCell(iCell), iCell) = heightCellNow(minLevelCell(iCell), iCell)
         do k = maxLevelCell(iCell), minLevelCell(iCell)+1, -1
            heightCellNew(k, iCell) = &
               heightCellNew(k + 1, iCell) - layerThicknessNew(k, iCell)
         end do
      end do
      !$omp end do
      !$omp end parallel

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(layerThickness, layerThickEdgeMean, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop, cellsOnEdge) &
      !$acc    private(k, kmin, kmax, cell1, cell2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, kmin, kmax, cell1, cell2)
#endif
      do iEdge = 1, nEdgesAll
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, nVertLevels
            ! initialize layerThicknessEdgeMean to avoid divide by
            ! zero and NaN problems.
            layerThickEdgeNew(k,iEdge) = -1.0e34_RKIND
         end do
         do k = kmin, kmax
            ! central differenced
            layerThickEdgeNew(k,iEdge) = 0.5_RKIND * &
                                         (layerThicknessNew(k,cell1) + &
                                          layerThicknessNew(k,cell2))
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$omp parallel
      !$omp do schedule(runtime) private(k, kTop, kBot, cell1, cell2)
      do iEdge = 1, nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         kTop = minLevelEdgeBot(iEdge)
         kBot = maxLevelEdgeTop(iEdge)

         ! We only need to define the kBot+1 depth for one active cell
         ! as the absolute kBot+1 depth has no effect on remapping
         if (cell1 <= nCellsAll .and. maxLevelCell(cell1) > 0) then
            heightEdgeNow(:, iEdge) = heightCellNow(kBot+1, cell1)
            heightEdgeNew(:, iEdge) = heightCellNew(kBot+1, cell1)
         elseif (cell2 <= nCellsAll .and. maxLevelCell(cell2) > 0) then
            heightEdgeNow(:, iEdge) = heightCellNow(kBot+1, cell2)
            heightEdgeNew(:, iEdge) = heightCellNew(kBot+1, cell2)
         else
            call mpas_log_write('No valid cell for edge remapping', MPAS_LOG_WARN)
            heightEdgeNow(:,iEdge) = 0.0_RKIND
            cycle
         endif

         do k = kBot, kTop, -1
            heightEdgeNow(k, iEdge) = heightEdgeNow(k+1, iEdge) - layerThickEdgeMean(k, iEdge)
            heightEdgeNew(k, iEdge) = heightEdgeNew(k+1, iEdge) - layerThickEdgeNew(k, iEdge)
         end do

      end do
      !$omp end do
      !$omp end parallel

      ! Assign new layerThicknesses
      layerThickness = layerThicknessNew

      !----------------------------------------------------------------------------
      ! ACTUAL REMAPPING FROM HERE

      nDoFs = 1 ! always 1 for MPAS, it'd be >1 for DG methods (with more DoF per layer)

      ! how many variables to be remapped concurrently per column
      if (config_use_freq_filtered_thickness) then
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

      call work%init(nVertLevels + 1, nVars, opts)

      ! Remap all edge-centred variables from heightEdgeNow to heightEdgeNew
      if (config_use_freq_filtered_thickness) then
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(kTop, kBot, nLayers, nLevels, dataNow, dataNew, work)
         do iEdge = 1, nEdges

            kTop = minLevelEdgeBot(iEdge)
            kBot = maxLevelEdgeTop(iEdge)

            nLayers = kBot - kTop + 1
            nLevels = nLayers + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, 1:nLayers) = normalVelocity   (kTop:kBot, iEdge)
            dataNow(1, 2, 1:nLayers) = highFreqThickness(kTop:kBot, iEdge)
            dataNow(1, 3, 1:nLayers) = lowFreqDivergence(kTop:kBot, iEdge)

            call rmap1d(nLevels, nLevels, nVars ,nDoFs, &
               heightEdgeNow(kTop:kBot+1, iEdge), &
               heightEdgeNew(kTop:kBot+1, iEdge), &
               dataNow, dataNew, &
               bcUpper, bcLower, work, opts)

            normalVelocity   (kTop:kBot, iEdge) = dataNew(1, 1, 1:nLayers)
            highFreqThickness(kTop:kBot, iEdge) = dataNew(1, 2, 1:nLayers)
            lowFreqDivergence(kTop:kBot, iEdge) = dataNew(1, 3, 1:nLayers)
         end do
         !$omp end do
         !$omp end parallel
      else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(kTop, kBot, nLayers, nLevels, dataNow, dataNew, work)
         do iEdge = 1, nEdges

            kTop = minLevelEdgeBot(iEdge)
            kBot = maxLevelEdgeTop(iEdge)

            nLayers = kBot - kTop + 1
            nLevels = nLayers + 1

            if (nLevels .lt. 2) cycle

            dataNow(1, 1, 1:nLayers) = normalVelocity(kTop:kBot, iEdge)

            call rmap1d(nLevels, nLevels, nVars, nDoFs, &
               heightEdgeNow(kTop:kBot+1, iEdge), &
               heightEdgeNew(kTop:kBot+1, iEdge), &
               dataNow, dataNew, &
               bcUpper, bcLower, work, opts)

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

               call work%init(nVertLevels + 1, nTracers, opts)

               allocate(bcUpper(nTracers))
               allocate(bcLower(nTracers))
               bcUpper%bcopt = bc_upper
               bcLower%bcopt = bc_lower

               !$omp parallel
               !$omp do schedule(runtime) &
               !$omp private(iTrac, kTop, kBot, nLayers, nLevels, dataNow, dataNew, work)
               do iCell = 1, nCells

                  kTop = minLevelCell(iCell)
                  kBot = maxLevelCell(iCell)

                  nLayers = kBot - kTop + 1
                  nLevels = nLayers + 1

                  if (nLevels .lt. 2) cycle

                  do iTrac = 1, nTracers

                     dataNow(1, iTrac, 1:nLayers) = tracersGroup(iTrac, kTop:kBot, iCell)

                  end do
                  call rmap1d(nLevels, nLevels, nTracers,nDoFs, &
                     heightCellNow(kTop:kBot+1, iCell), &
                     heightCellNew(kTop:kBot+1, iCell), &
                     dataNow, dataNew, &
                     bcUpper, bcLower, work, opts)

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

      call ocn_vert_regrid_init(err)

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
      else
         CALL mpas_log_write('Vertical remap order $i not supported', &
            MPAS_LOG_CRIT, intArgs=(/config_vert_remap_order/))
      endif

      if ( config_remap_limiter == 'monotonic' ) then
         opts%cell_lims = mono_limit ! monotone slope limits
      elseif ( config_remap_limiter == 'weno' ) then
         opts%cell_lims = weno_limit ! WENO slope limits
         opts%wall_lims = mono_limit ! monotone slope limits
      elseif ( config_remap_limiter == 'none' ) then
         opts%cell_lims = null_limit ! no slope limits
      else
         CALL mpas_log_write('Remap limiter not supported', &
            MPAS_LOG_CRIT)
      endif

      bc_upper = bcon_loose
      bc_lower = bcon_loose

      itimestepLastRemap = 0

   end subroutine ocn_vertical_remap_init

end module ocn_vertical_remap
! vim: foldmethod=marker
