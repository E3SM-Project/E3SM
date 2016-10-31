! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_mask
!
!> \MPAS land-ice mask calculations
!> \author Matt Hoffman
!> \date   10 May 2012
!> \details
!>  This module contains the routines for calculating masks for land ice
!>
!
!-----------------------------------------------------------------------

module li_mask

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timer
   use li_setup

   implicit none

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------
   integer, parameter :: li_mask_ValueIce                   =  32  ! Giving this the highest current value so it is obvious during visualization
   integer, parameter :: li_mask_ValueDynamicIce            =   2
   integer, parameter :: li_mask_ValueFloating              =   4
   integer, parameter :: li_mask_ValueMargin                =   8  ! This is the last cell with ice.
   integer, parameter :: li_mask_ValueDynamicMargin         =  16  ! This is the last dynamically active cell with ice
   integer, parameter :: li_mask_ValueInitialIceExtent      =   1
   integer, parameter :: li_mask_ValueAlbanyActive          =  64  ! These are locations that Albany includes in its solution
   integer, parameter :: li_mask_ValueAlbanyMarginNeighbor  = 128  ! This the first cell beyond the last active albany cell
   integer, parameter :: li_mask_ValueGroundingLine         = 256  ! This is grounded cell that has a floating neighbor, or vertex/edge on that boundary

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   ! all subroutines and functions in this module are public!


   ! interfaces without a suffix return logicals
   ! interfaces with names that end with '_int' return 0/1
   ! TODO Eventually we may decide to only keep and maintain one of these return types.

   interface li_mask_is_ice
      module procedure li_mask_is_ice_logout_1d
      module procedure li_mask_is_ice_logout_0d
   end interface


   interface li_mask_is_ice_int
      module procedure li_mask_is_ice_intout_1d
      module procedure li_mask_is_ice_intout_0d
   end interface


   interface li_mask_is_dynamic_ice
      module procedure li_mask_is_dynamic_ice_logout_1d
      module procedure li_mask_is_dynamic_ice_logout_0d
   end interface


   interface li_mask_is_albany_active
      module procedure li_mask_is_albany_active_logout_1d
      module procedure li_mask_is_albany_active_logout_0d
   end interface


   interface li_mask_is_dynamic_ice_int
      module procedure li_mask_is_dynamic_ice_intout_1d
      module procedure li_mask_is_dynamic_ice_intout_0d
   end interface


   interface li_mask_is_dynamic_margin
      module procedure li_mask_is_dynamic_margin_logout_1d
      module procedure li_mask_is_dynamic_margin_logout_0d
   end interface


   interface li_mask_is_dynamic_margin_int
      module procedure li_mask_is_dynamic_margin_logout_1d
      module procedure li_mask_is_dynamic_margin_logout_0d
   end interface


   interface li_mask_is_albany_margin_neighbor
      module procedure li_mask_is_albany_margin_neighbor_logout_1d
      module procedure li_mask_is_albany_margin_neighbor_logout_0d
   end interface


   interface li_mask_is_floating_ice
      module procedure li_mask_is_floating_ice_logout_1d
      module procedure li_mask_is_floating_ice_logout_0d
   end interface


   interface li_mask_is_floating_ice_int
      module procedure li_mask_is_floating_ice_intout_1d
      module procedure li_mask_is_floating_ice_intout_0d
   end interface

   interface li_mask_is_grounded_ice
      module procedure li_mask_is_grounded_ice_logout_1d
      module procedure li_mask_is_grounded_ice_logout_0d
   end interface

   interface li_mask_is_initial_ice
      module procedure li_mask_is_initial_ice_logout_1d
      module procedure li_mask_is_initial_ice_logout_0d
   end interface

!! SFP: the 1d versions below have not been tested yet - activate and use at your own risk! 
!! Also check the corresponding function below.
   interface li_mask_is_grounded_ice_int
      module procedure li_mask_is_grounded_ice_intout_1d
      module procedure li_mask_is_grounded_ice_intout_0d
   end interface


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



!***********************************************************************

contains



!***********************************************************************
!
!  routine li_calculate_mask_init
!
!> \brief   Calculates masks for land ice for info needed from initial condition only
!> \author  Matt Hoffman
!> \date    25 June 2012
!> \details
!>  This routine Calculates masks for land ice for info needed from initial condition only.
!
!-----------------------------------------------------------------------

   subroutine li_calculate_mask_init(geometryPool, err)

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
      type (mpas_pool_type), intent(inout) :: &
         geometryPool          !< Input/Output: geometry information

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
      integer, dimension(:), pointer :: cellMask
      real(KIND=RKIND), dimension(:), pointer :: thickness
      logical, pointer :: config_do_restart

      err = 0

      ! Assign pointers and variables
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)

      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)

      if (config_do_restart .eqv. .false.) then  ! We only want to set this bit of the mask when a new simulation starts, but not during a restart.
         ! Initialize cell mask to 0 everywhere before we assign anything to it.
         cellMask = 0
         where (thickness > 0.0)
            cellMask = ior(cellMask, li_mask_ValueInitialIceExtent)
         end where
      endif

   !--------------------------------------------------------------------

   end subroutine li_calculate_mask_init



!***********************************************************************
!
!  routine land_ice_calculate_mask
!
!> \brief   Calculates masks for land ice
!> \author  Matt Hoffman
!> \date    10 May 2012
!> \details
!>  This routine Calculates masks for land ice.
!
!-----------------------------------------------------------------------

   subroutine li_calculate_mask(meshPool, velocityPool, geometryPool, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(inout) :: &
         velocityPool          !< Input: velocity information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: &
         geometryPool          !< Input/Output: geometry information

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
      integer, pointer :: nCells, nVertices, nEdges, vertexDegree
      integer, pointer :: nVertInterfaces
      real(KIND=RKIND), dimension(:), pointer :: thickness, bedTopography
      integer, dimension(:), pointer :: nEdgesOnCell, cellMask, vertexMask, edgeMask
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, cellsOnEdge, dirichletVelocityMask
      real (kind=RKIND), pointer :: config_ice_density, config_ocean_density, &
            config_sea_level, config_dynamic_thickness
      character (len=StrKIND), pointer :: config_velocity_solver

      integer :: i, j, iCell
      logical :: isMargin
      logical :: isAlbanyMarginNeighbor
      logical :: aCellOnVertexHasIce, aCellOnVertexHasNoIce, aCellOnVertexHasDynamicIce, aCellOnVertexHasNoDynamicIce, aCellOnVertexIsFloating, aCellOnVertexIsAlbanyActive
      logical :: aCellOnVertexIsGrounded
      logical :: aCellOnEdgeHasIce, aCellOnEdgeHasNoIce, aCellOnEdgeHasDynamicIce, aCellOnEdgeHasNoDynamicIce, aCellOnEdgeIsFloating 
      logical :: aCellOnEdgeIsGrounded
      integer :: numCellsOnVertex
      integer :: numDiriDynamicCells, numDiriNondynamicCells, numExtendedCells
      logical :: validVertex

      call mpas_timer_start('calculate mask')

      err = 0

      ! Assign pointers and variables
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel=1)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)

      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      if ( .not. ((trim(config_velocity_solver) == 'sia') .or. (trim(config_velocity_solver) == 'none')) ) then
         call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMask, timeLevel = 1)
      endif

      ! ====
      ! Calculate cellMask values===========================
      ! ====

      !call mpas_timer_start('calculate mask cell')
      ! Set mask to 0 everywhere, but need to preserve bits the initial ice extent bit
      do i=1, nCells
        cellMask(i) = iand(cellMask(i), li_mask_ValueInitialIceExtent)
      enddo
      
      ! Identify cells with ice
      where (thickness > 0)
          cellMask = ior(cellMask, li_mask_ValueIce)
      end where

      ! Identify cells where the ice is above the ice dynamics thickness limit
      where ( thickness > config_dynamic_thickness )
         cellMask = ior(cellMask, li_mask_ValueDynamicIce)
      end where

      ! Identify cells that Albany would consider as active
      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         ! HO external FEM dycore
         ! Identify cells where the ice is above the ice dynamics thickness limit but not with a dirichletVelocity condition set
         where ( (thickness > config_dynamic_thickness) .and.    &  ! same as for SIA case
                 (maxval(dirichletVelocityMask(1:nVertInterfaces-1, :), dim=1) == 0) )   ! but exclude dirichletVelocityMask locations set as lateral b.c.  We don't want to consider dirichlet b.c. on the basal boundary, so we ignore the basal level
            cellMask = ior(cellMask, li_mask_ValueAlbanyActive)
         end where
      endif

      ! Is it floating? (ice thickness equal to floatation is considered floating)
      ! For now floating ice and grounded ice are mutually exclusive.  
      ! This may change if a grounding line parameterization is added.
      where (  li_mask_is_ice(cellMask) .and. (config_ice_density / config_ocean_density * thickness) <= (config_sea_level - bedTopography) )
          cellMask = ior(cellMask, li_mask_ValueFloating)
      end where

      ! Identify the margin
      ! For a cell, we define the margin as the last cell with ice (the cell has ice and at least one neighbor is a non-ice cell)
      do i=1,nCells      
          if (li_mask_is_ice(cellMask(i))) then
              isMargin = .false.
              do j=1,nEdgesOnCell(i) ! Check if any neighbors are non-ice
                  isMargin = ( isMargin .or. (.not. li_mask_is_ice(cellMask(cellsOnCell(j,i)))) )
              enddo
              if (isMargin) then
                 cellMask(i) = ior(cellMask(i), li_mask_ValueMargin)
              endif
          endif
      enddo

      ! Identify the dynamic margin
      ! For a cell, we define the dynamic margin as the last cell with dynamic ice (the cell is dynamic and at least one neighboring cell is not dynamic)
      do i=1,nCells      
          if (li_mask_is_dynamic_ice(cellMask(i))) then
              isMargin = .false.
              do j=1,nEdgesOnCell(i) ! Check if any neighbors are not dynamic
                  isMargin = ( isMargin .or. (.not. li_mask_is_dynamic_ice(cellMask(cellsOnCell(j,i)))) )
              enddo
              if (isMargin) then
                 cellMask(i) = ior(cellMask(i), li_mask_ValueDynamicMargin)
              endif
          endif
      enddo

      ! Identify first cell outside of active albany extent
      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         do i=1,nCells
            if ( (.not. li_mask_is_albany_active(cellMask(i)))) then  ! check non albany cells only
              isAlbanyMarginNeighbor = .false.
              do j=1,nEdgesOnCell(i) ! Check if any neighbors are dynamic-ice
                  isAlbanyMarginNeighbor = ( isAlbanyMarginNeighbor .or. (li_mask_is_albany_active(cellMask(cellsOnCell(j,i)))) )
              enddo
              if (isAlbanyMarginNeighbor) then
                 cellMask(i) = ior(cellMask(i), li_mask_ValueAlbanyMarginNeighbor)
              endif
          endif
         enddo
      endif

      ! Identify the grounding line
      ! For a cell, we define the GL as a grounded cell with ice with at least one neighbor with floating ice
      do i=1,nCells
          if (li_mask_is_grounded_ice(cellMask(i))) then  ! only need to check grounded cells
              do j=1,nEdgesOnCell(i) ! Check if any neighbors are floating
                 if (li_mask_is_floating_ice(cellMask(cellsOnCell(j,i)))) &
                    cellMask(i) = ior(cellMask(i), li_mask_ValueGroundingLine)
                 cycle  ! no need to look at additional neighbors
              enddo
          endif
      enddo

      !call mpas_timer_stop('calculate mask cell')

      ! ====
      ! Calculate vertexMask values based on cellMask values===========================
      ! ====
      ! Bit: Vertices with ice are ones with at least one adjacent cell with ice
      ! Bit: Vertices with dynamic ice are ones with at least one adjacent cell with dynamic ice
      ! Bit: Floating vertices have at least one neighboring cell floating
      ! Bit: Vertices on margin are vertices with at least one neighboring cell with ice and at least one neighboring cell without ice
      ! Bit: Vertices on dynamic margin are vertices with at least one neighboring cell with dynamic ice and at least one neighboring cell without dynamic ice
      ! NOTE: Vertices are only considered 'valid' if they have three valid neighboring 
      !       cells (i.e., the cells exist in the mesh).  This allows external dycores to use the 
      !       vertexMask to get information about triangles in the Delaunay triangulation.
      !       (This is done in a way which does not assume vertexMask==3.)
      ! Bit: GL is a vertex with at least one neighboring cell grounded ice and at least one neighboring cell floating ice

      !call mpas_timer_start('calculate mask vertex')

      vertexMask = 0
      do i = 1,nVertices
          aCellOnVertexHasIce = .false.
          aCellOnVertexHasNoIce = .false.
          aCellOnVertexHasDynamicIce = .false.
          aCellOnVertexHasNoDynamicIce = .false.
          aCellOnVertexIsFloating = .false.
          aCellOnVertexIsGrounded = .false.
          do j = 1, vertexDegree  ! vertexDegree is usually 3 (e.g. CVT mesh) but could be something else (e.g. 4 for quad mesh)
              iCell = cellsOnVertex(j,i)
              aCellOnVertexHasIce = (aCellOnVertexHasIce .or. li_mask_is_ice(cellMask(iCell)))
              aCellOnVertexHasNoIce = (aCellOnVertexHasNoIce .or. (.not. li_mask_is_ice(cellMask(iCell))))
              aCellOnVertexHasDynamicIce = (aCellOnVertexHasDynamicIce .or. li_mask_is_dynamic_ice(cellMask(iCell)))
              aCellOnVertexHasNoDynamicIce = (aCellOnVertexHasNoDynamicIce .or. (.not. (li_mask_is_dynamic_ice(cellMask(iCell)))))
              aCellOnVertexIsFloating = (aCellOnVertexIsFloating .or. li_mask_is_floating_ice(cellMask(iCell)))
              aCellOnVertexIsGrounded = (aCellOnVertexIsGrounded .or. li_mask_is_grounded_ice(cellMask(iCell)))
          end do
          if (aCellOnVertexHasIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueIce)
          endif
          if (aCellOnVertexHasDynamicIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueDynamicIce)
          endif
          if (aCellOnVertexIsFloating) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueFloating)
          endif
          if (aCellOnVertexIsFloating .and. aCellOnVertexIsGrounded) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueGroundingLine)
          endif
          if (aCellOnVertexHasIce .and. aCellOnVertexHasNoIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueMargin)     ! vertex with both 1+ ice cell and 1+ non-ice cell as neighbors
          endif
          if (aCellOnVertexHasDynamicIce .and. aCellOnVertexHasNoDynamicIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueDynamicMargin)     ! vertex with both 1+ dynamic ice cell(s) and 1+ non-dynamic cell(s) as neighbors
          endif
      end do


      ! Re-loop over vertices for information only needed by Albany
      ! This is 10x more expensive to conditionally include in the above loop.
      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         do i = 1,nVertices
            aCellOnVertexIsAlbanyActive = .false.
            numCellsOnVertex = 0
            validVertex = .false.
            numDiriDynamicCells = 0
            numDiriNondynamicCells = 0
            numExtendedCells = 0
            do j = 1, vertexDegree  ! vertexDegree is usually 3 (e.g. CVT mesh) but could be something else (e.g. 4 for quad mesh)
               iCell = cellsOnVertex(j,i)
               if (iCell < nCells+1) then
                  numCellsOnVertex = numCellsOnVertex + 1
               endif
               aCellOnVertexIsAlbanyActive = (aCellOnVertexIsAlbanyActive .or. li_mask_is_albany_active(cellMask(iCell)))

               !if (li_mask_is_dynamic_ice(cellMask(iCell)) .and. .not. li_mask_is_albany_active(cellMask(iCell))) then ! this finds diri cells
               if ( (maxval(dirichletVelocityMask(1:nVertInterfaces-1, iCell)) > 0) .and. &
                       (li_mask_is_dynamic_ice(cellMask(iCell)) ) ) then
                       numDiriDynamicCells = numDiriDynamicCells + 1
               elseif ( (maxval(dirichletVelocityMask(1:nVertInterfaces-1, iCell)) > 0) .and. &
                       (.not. li_mask_is_dynamic_ice(cellMask(iCell)) ) ) then
                       numDiriNondynamicCells = numDiriNondynamicCells + 1
               elseif (li_mask_is_albany_margin_neighbor(cellMask(iCell))) then
                       numExtendedCells = numExtendedCells + 1
               endif
            end do
            if (numCellsOnVertex == vertexDegree) validVertex = .true.
            if (aCellOnVertexIsAlbanyActive .and. validVertex)  vertexMask(i) = ior(vertexMask(i), li_mask_ValueAlbanyActive)
            if ( (numDiriDynamicCells == 1) .and. (numDiriNondynamicCells == 1) .and. &
                 (numExtendedCells == 1) .and. validVertex) then
               ! This is a special case needed for MISMIP
               vertexMask(i) = ior(vertexMask(i), li_mask_ValueAlbanyActive)
               vertexMask(i) = ior(vertexMask(i), li_mask_ValueAlbanyMarginNeighbor)
            endif
         end do  ! vertices loop
      endif  ! HO velocity solver

      !call mpas_timer_stop('calculate mask vertex')

      ! ====
      ! Calculate edgeMask values based on cellMask values===========================
      ! ====
      ! Bit: Edges with ice are ones with at least one adjacent cell with ice
      ! Bit: Edges with dynamic ice are ones with at least one adjacent cell with dynamic ice
      ! Bit: Floating Edges have at least one neighboring cell floating
      ! Bit: Edges on margin are edges with one neighboring cell with ice and one neighboring cell without ice
      ! Bit: Edges on dynamic margin are edges with one neighboring cell with dynamic ice and one neighboring cell without dynamic ice
      ! Bit: GL is an edge with one cell grounded ice and one cell floating ice

      !call mpas_timer_start('calculate mask edge')

      edgeMask = 0
      do i = 1,nEdges
          aCellOnEdgeHasIce = .false.
          aCellOnEdgeHasNoIce = .false.
          aCellOnEdgeHasDynamicIce = .false.
          aCellOnEdgeHasNoDynamicIce = .false.
          aCellOnEdgeIsFloating = .false.
          aCellOnEdgeIsGrounded = .false.
          do j = 1, 2
              iCell = cellsOnEdge(j,i)
              aCellOnEdgeHasIce = (aCellOnEdgeHasIce .or. li_mask_is_ice(cellMask(iCell)))
              aCellOnEdgeHasNoIce = (aCellOnEdgeHasNoIce .or. (.not. li_mask_is_ice(cellMask(iCell))))
              aCellOnEdgeHasDynamicIce = (aCellOnEdgeHasDynamicIce .or. li_mask_is_dynamic_ice(cellMask(iCell)))
              aCellOnEdgeHasNoDynamicIce = (aCellOnEdgeHasNoDynamicIce .or. (.not. (li_mask_is_dynamic_ice(cellMask(iCell)))))
              aCellOnEdgeIsFloating = (aCellOnEdgeIsFloating .or. li_mask_is_floating_ice(cellMask(iCell)))
              aCellOnEdgeIsGrounded = (aCellOnEdgeIsGrounded .or. li_mask_is_grounded_ice(cellMask(iCell)))
          end do
          if (aCellOnEdgeHasIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueIce)
          endif
          if (aCellOnEdgeHasDynamicIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueDynamicIce)
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueAlbanyActive)  ! Note: Albany does not use edgeMask, but setting this anyway.
          endif
          if (aCellOnEdgeIsFloating) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueFloating)
          endif
          if (aCellOnEdgeIsFloating .and. aCellOnEdgeIsGrounded) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueGroundingLine)
          endif
          if (aCellOnEdgeHasIce .and. aCellOnEdgeHasNoIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueMargin)
          endif
          if (aCellOnEdgeHasDynamicIce .and. aCellOnEdgeHasNoDynamicIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueDynamicMargin)
          endif

      end do
      !call mpas_timer_stop('calculate mask edge')

      ! vertexMask and edgeMask needs halo updates before they can be used.  Halo updates need to occur outside of block loops.  

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_calculate_mask."
      endif

      call mpas_timer_stop('calculate mask')

   !--------------------------------------------------------------------
   end subroutine li_calculate_mask


!***********************************************************************
!
!  routine li_calculate_extrapolate_floating_edgemask
!
!> \brief   Extrapolates floating edges forward as needed by external FEM dycores
!> \author  Matt Hoffman
!> \date    29 January 2015
!> \details
!>  External FEM dycores include the first non-ice cells in their mesh.  They
!>  also use a mask to apply floating lateral boundary conditions on edges.
!>  Because they include extra cell center locations in their meshes, the triangle
!>  edges connecting these extra nodes will not be covered by the standard
!>  MPAS edge mask.  This routine deals with this problem by 'extrapolating'
!>  the floating edge mask forward to cover the edges connecting these extra nodes.
!>  It does so by looping over edges, and setting as floating any edge that has
!>  at least one neighboring vertex that is 'floating'.  This makes use of the
!>  convention that "Floating vertices have at least one neighboring cell floating".
!
!-----------------------------------------------------------------------

   subroutine li_calculate_extrapolate_floating_edgemask(meshPool, vertexMask, floatingEdges)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      integer, dimension(:), intent(inout) :: &
         vertexMask          !< Input/Output: vertexMask

      integer, dimension(:), intent(inout) :: &
         floatingEdges          !< Input/Output: 0/1 mask of floating edges

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer, dimension(:,:), pointer :: verticesOnEdge
      integer, pointer :: nEdges
      integer :: iEdge

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)

      ! Build floatingEdges mask that is extended forward one extra edge
      do iEdge = 1, nEdges
         floatingEdges(iEdge) = maxval(li_mask_is_floating_ice_int(vertexMask(verticesOnEdge(:, iEdge))))
      enddo

      ! Now includes vertices that include 2 cells with ice and one extended floating cell
      ! This will be when the two cells with ice were Dirichlet cells
   end subroutine li_calculate_extrapolate_floating_edgemask


   ! ===================================
   ! Functions for decoding bitmasks - will work with cellMask, edgeMask, or vertexMask
   ! ===================================
   ! Only adding the minimum needed for now.  These should be added as needed.
   ! functions with names that include '_logout' return logical types
   !    -- these should be used with 'if' and 'where' statements
   ! functions with names that include '_intout' return integers types with 0 for false, 1 for true.
   !    -- these should be used when multiplying against numeric arrays


   ! -- Functions that check for presence of ice --
   function li_mask_is_ice_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_ice_logout_1d

      li_mask_is_ice_logout_1d = (iand(mask, li_mask_ValueIce) == li_mask_ValueIce)
   end function li_mask_is_ice_logout_1d

   function li_mask_is_ice_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_ice_logout_0d

      li_mask_is_ice_logout_0d = (iand(mask, li_mask_ValueIce) == li_mask_ValueIce)
   end function li_mask_is_ice_logout_0d


   function li_mask_is_ice_intout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      integer, dimension(size(mask)) :: li_mask_is_ice_intout_1d

      li_mask_is_ice_intout_1d = iand(mask, li_mask_ValueIce) / li_mask_ValueIce
   end function li_mask_is_ice_intout_1d

   function li_mask_is_ice_intout_0d(mask)
      integer, intent(in) :: mask
      integer :: li_mask_is_ice_intout_0d

      li_mask_is_ice_intout_0d = iand(mask, li_mask_ValueIce) / li_mask_ValueIce
   end function li_mask_is_ice_intout_0d


   ! -- Functions that check for presence of dynamic ice --
   function li_mask_is_dynamic_ice_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_dynamic_ice_logout_1d

      li_mask_is_dynamic_ice_logout_1d = (iand(mask, li_mask_ValueDynamicIce) == li_mask_ValueDynamicIce)
   end function li_mask_is_dynamic_ice_logout_1d

   function li_mask_is_dynamic_ice_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_dynamic_ice_logout_0d

      li_mask_is_dynamic_ice_logout_0d = (iand(mask, li_mask_ValueDynamicIce) == li_mask_ValueDynamicIce)
   end function li_mask_is_dynamic_ice_logout_0d

   function li_mask_is_dynamic_ice_intout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      integer, dimension(size(mask)) :: li_mask_is_dynamic_ice_intout_1d

      li_mask_is_dynamic_ice_intout_1d = iand(mask, li_mask_ValueDynamicIce) / li_mask_ValueDynamicIce
   end function li_mask_is_dynamic_ice_intout_1d

   function li_mask_is_dynamic_ice_intout_0d(mask)
      integer, intent(in) :: mask
      integer :: li_mask_is_dynamic_ice_intout_0d

      li_mask_is_dynamic_ice_intout_0d = iand(mask, li_mask_ValueDynamicIce) / li_mask_ValueDynamicIce
   end function li_mask_is_dynamic_ice_intout_0d


   ! -- Functions that check for presence of active albany --
   function li_mask_is_albany_active_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_albany_active_logout_1d

      li_mask_is_albany_active_logout_1d = (iand(mask, li_mask_ValueAlbanyActive) == li_mask_ValueAlbanyActive)
   end function li_mask_is_albany_active_logout_1d

   function li_mask_is_albany_active_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_albany_active_logout_0d

      li_mask_is_albany_active_logout_0d = (iand(mask, li_mask_ValueAlbanyActive) == li_mask_ValueAlbanyActive)
   end function li_mask_is_albany_active_logout_0d


   ! -- Functions that check for presence of dynamic margin --
   function li_mask_is_dynamic_margin_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_dynamic_margin_logout_1d

      li_mask_is_dynamic_margin_logout_1d = (iand(mask, li_mask_ValueDynamicMargin) == li_mask_ValueDynamicMargin)
   end function li_mask_is_dynamic_margin_logout_1d

   function li_mask_is_dynamic_margin_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_dynamic_margin_logout_0d

      li_mask_is_dynamic_margin_logout_0d = (iand(mask, li_mask_ValueDynamicMargin) == li_mask_ValueDynamicMargin)
   end function li_mask_is_dynamic_margin_logout_0d

   function li_mask_is_dynamic_margin_intout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      integer, dimension(size(mask)) :: li_mask_is_dynamic_margin_intout_1d

      li_mask_is_dynamic_margin_intout_1d = iand(mask, li_mask_ValueDynamicMargin) / li_mask_ValueDynamicMargin
   end function li_mask_is_dynamic_margin_intout_1d

   function li_mask_is_dynamic_margin_intout_0d(mask)
      integer, intent(in) :: mask
      integer :: li_mask_is_dynamic_margin_intout_0d

      li_mask_is_dynamic_margin_intout_0d = iand(mask, li_mask_ValueDynamicMargin) / li_mask_ValueDynamicMargin
   end function li_mask_is_dynamic_margin_intout_0d


   ! -- Functions that check for presence of albany margin neighbor --
   function li_mask_is_albany_margin_neighbor_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_albany_margin_neighbor_logout_1d

      li_mask_is_albany_margin_neighbor_logout_1d = (iand(mask, li_mask_ValueAlbanyMarginNeighbor) == li_mask_ValueAlbanyMarginNeighbor)
   end function li_mask_is_albany_margin_neighbor_logout_1d

   function li_mask_is_albany_margin_neighbor_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_albany_margin_neighbor_logout_0d

      li_mask_is_albany_margin_neighbor_logout_0d = (iand(mask, li_mask_ValueAlbanyMarginNeighbor) == li_mask_ValueAlbanyMarginNeighbor)
   end function li_mask_is_albany_margin_neighbor_logout_0d


   ! -- Functions that check for presence of floating ice --
   function li_mask_is_floating_ice_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_floating_ice_logout_1d

      li_mask_is_floating_ice_logout_1d = (iand(mask, li_mask_ValueFloating) == li_mask_ValueFloating)
   end function li_mask_is_floating_ice_logout_1d

   function li_mask_is_floating_ice_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_floating_ice_logout_0d

      li_mask_is_floating_ice_logout_0d = (iand(mask, li_mask_ValueFloating) == li_mask_ValueFloating)
   end function li_mask_is_floating_ice_logout_0d

   function li_mask_is_floating_ice_intout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      integer, dimension(size(mask)) :: li_mask_is_floating_ice_intout_1d

      li_mask_is_floating_ice_intout_1d = iand(mask, li_mask_ValueFloating) / li_mask_ValueFloating
   end function li_mask_is_floating_ice_intout_1d

   function li_mask_is_floating_ice_intout_0d(mask)
      integer, intent(in) :: mask
      integer :: li_mask_is_floating_ice_intout_0d

      li_mask_is_floating_ice_intout_0d = iand(mask, li_mask_ValueFloating) / li_mask_ValueFloating
   end function li_mask_is_floating_ice_intout_0d

   ! -- Functions that check for presence of grounded ice --
   function li_mask_is_grounded_ice_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_grounded_ice_logout_1d

      li_mask_is_grounded_ice_logout_1d = ( (iand(mask, li_mask_ValueFloating) /= li_mask_ValueFloating) &
                                .and. (li_mask_is_ice(mask)) )
   end function li_mask_is_grounded_ice_logout_1d

   function li_mask_is_grounded_ice_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_grounded_ice_logout_0d

      li_mask_is_grounded_ice_logout_0d = ( (iand(mask, li_mask_ValueFloating) /= li_mask_ValueFloating) &
                                .and. (li_mask_is_ice(mask)) )
   end function li_mask_is_grounded_ice_logout_0d

   function li_mask_is_grounded_ice_intout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      integer, dimension(size(mask)) :: li_mask_is_grounded_ice_intout_1d

      where( li_mask_is_ice(mask) )
        li_mask_is_grounded_ice_intout_1d = int(mask*0 + 1) - ( iand(mask, li_mask_ValueFloating) / li_mask_ValueFloating )
      elsewhere
        li_mask_is_grounded_ice_intout_1d = 0
      endwhere
   end function li_mask_is_grounded_ice_intout_1d

   function li_mask_is_grounded_ice_intout_0d(mask)
      integer, intent(in) :: mask
      integer :: li_mask_is_grounded_ice_intout_0d

      if( li_mask_is_ice(mask) )then
        li_mask_is_grounded_ice_intout_0d = 1 - ( iand(mask, li_mask_ValueFloating) / li_mask_ValueFloating ) 
      else
        li_mask_is_grounded_ice_intout_0d = 0
      endif
   end function li_mask_is_grounded_ice_intout_0d

   function li_mask_is_initial_ice_logout_1d(mask)
      integer, dimension(:), intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_initial_ice_logout_1d

      li_mask_is_initial_ice_logout_1d = (iand(mask, li_mask_ValueInitialIceExtent) == li_mask_ValueInitialIceExtent)
   end function li_mask_is_initial_ice_logout_1d

   function li_mask_is_initial_ice_logout_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_initial_ice_logout_0d

      li_mask_is_initial_ice_logout_0d = (iand(mask, li_mask_ValueInitialIceExtent) == li_mask_ValueInitialIceExtent)
   end function li_mask_is_initial_ice_logout_0d



!***********************************************************************
! Private subroutines:
!***********************************************************************

! - no private subroutines - (module is not declared private)


end module li_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

