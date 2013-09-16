! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for calculating masks for land ice
!>
!
!-----------------------------------------------------------------------

module li_mask

   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

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

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   ! all subroutines and functions in this module are public!

   interface li_mask_is_ice
      module procedure li_mask_is_ice_1d
      module procedure li_mask_is_ice_0d
   end interface


   interface li_mask_is_dynamic_ice
      module procedure li_mask_is_dynamic_ice_1d
      module procedure li_mask_is_dynamic_ice_0d
   end interface


   interface li_mask_is_floating_ice
      module procedure li_mask_is_floating_ice_1d
      module procedure li_mask_is_floating_ice_0d
   end interface


   interface li_mask_is_grounded_ice
      module procedure li_mask_is_grounded_ice_1d
      module procedure li_mask_is_grounded_ice_0d
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
!> \version SVN:$Id$
!> \details
!>  This routine Calculates masks for land ice for info needed from initial condition only.
!
!-----------------------------------------------------------------------

   subroutine li_calculate_mask_init(mesh, state, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state          !< Input/Output: state information 

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

      ! Assign pointers and variables
      cellMask => state % cellMask % array
      thickness => state % thickness % array

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
!> \version SVN:$Id$
!> \details
!>  This routine Calculates masks for land ice.
!
!-----------------------------------------------------------------------

   subroutine li_calculate_mask(mesh, state, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state          !< Input/Output: state information 

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


      integer :: nCells, nVertices, nEdges

      real(KIND=RKIND), dimension(:), pointer :: thickness, bedTopography

      integer, dimension(:), pointer :: nEdgesOnCell, cellMask, vertexMask, edgeMask

      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, cellsOnEdge

      integer :: i, j, iCell

      logical :: isMargin

      logical :: aCellOnVertexHasIce, aCellOnVertexHasNoIce, aCellOnVertexHasDynamicIce, aCellOnVertexHasNoDynamicIce, aCellOnVertexIsFloating

      logical :: aCellOnEdgeHasIce, aCellOnEdgeHasNoIce, aCellOnEdgeHasDynamicIce, aCellOnEdgeHasNoDynamicIce, aCellOnEdgeIsFloating 


      err = 0


      ! Assign pointers and variables
      cellMask => state % cellMask % array
      thickness => state % thickness % array
      bedTopography => mesh  % bedTopography % array
      vertexMask => state % vertexMask % array
      edgeMask => state % edgeMask % array
      nEdgesOnCell => mesh % nEdgesOnCell % array
      cellsOnCell => mesh % cellsOnCell % array
      cellsOnVertex => mesh % cellsOnVertex % array
      cellsOnEdge => mesh % cellsOnEdge % array

      nCells = mesh % nCells
      nVertices = mesh % nVertices
      nEdges = mesh % nEdges

      ! ====
      ! Calculate cellMask values===========================
      ! ====

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
      ! see below for an additional check

      ! Is it floating? (ice thickness equal to floatation is considered floating)
      ! For now floating ice and grounded ice are mutually exclusive.  
      ! This may change if a ground line parameterization is added.
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


      ! ====
      ! Calculate vertexMask values based on cellMask values===========================
      ! ====
      ! Bit: Vertices with ice are ones with at least one adjacent cell with ice
      ! Bit: Vertices with dynamic ice are ones with at least one adjacent cell with dynamic ice
      ! Bit: Floating vertices have at least one neighboring cell floating
      ! Bit: Vertices on margin are vertices with at least one neighboring cell with ice and at least one neighboring cell without ice
      ! Bit: Vertices on dynamic margin are vertices with at least one neighboring cell with dynamic ice and at least one neighboring cell without dynamic ice
      vertexMask = 0
      do i = 1,nVertices
          aCellOnVertexHasIce = .false.
          aCellOnVertexHasNoIce = .false.
          aCellOnVertexHasDynamicIce = .false.
          aCellOnVertexHasNoDynamicIce = .false.
          aCellOnVertexIsFloating = .false.
          do j = 1, mesh % vertexDegree  ! vertexDegree is usually 3 (e.g. CVT mesh) but could be something else (e.g. 4 for quad mesh)
              iCell = cellsOnVertex(j,i)
              aCellOnVertexHasIce = (aCellOnVertexHasIce .or. li_mask_is_ice(cellMask(iCell)))
              aCellOnVertexHasNoIce = (aCellOnVertexHasNoIce .or. (.not. li_mask_is_ice(cellMask(iCell))))
              aCellOnVertexHasDynamicIce = (aCellOnVertexHasDynamicIce .or. li_mask_is_dynamic_ice(cellMask(iCell)))
              aCellOnVertexHasNoDynamicIce = (aCellOnVertexHasNoDynamicIce .or. (.not. (li_mask_is_dynamic_ice(cellMask(iCell)))))
              aCellOnVertexIsFloating = (aCellOnVertexIsFloating .or. li_mask_is_floating_ice(cellMask(iCell)))
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
          if (aCellOnVertexHasIce .and. aCellOnVertexHasNoIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueMargin)     ! vertex with both 1+ ice cell and 1+ non-ice cell as neighbors
          endif
          if (aCellOnVertexHasDynamicIce .and. aCellOnVertexHasNoDynamicIce) then
             vertexMask(i) = ior(vertexMask(i), li_mask_ValueDynamicMargin)     ! vertex with both 1+ dynamic ice cell(s) and 1+ non-dynamic cell(s) as neighbors
          endif
      end do       


      ! ====
      ! Calculate edgeMask values based on cellMask values===========================
      ! ====
      ! Bit: Edges with ice are ones with at least one adjacent cell with ice
      ! Bit: Edges with dynamic ice are ones with at least one adjacent cell with dynamic ice
      ! Bit: Floating Edges have at least one neighboring cell floating
      ! Bit: Edges on margin are vertices with one neighboring cell with ice and one neighboring cell without ice
      ! Bit: Edges on dynamic margin are vertices with at least one neighboring cell with dynamic ice and at least one neighboring cell without dynamic ice
      edgeMask = 0
      do i = 1,nEdges
          aCellOnEdgeHasIce = .false.
          aCellOnEdgeHasNoIce = .false.
          aCellOnEdgeHasDynamicIce = .false.
          aCellOnEdgeHasNoDynamicIce = .false.
          aCellOnEdgeIsFloating = .false.
          do j = 1, 2
              iCell = cellsOnEdge(j,i)
              aCellOnEdgeHasIce = (aCellOnEdgeHasIce .or. li_mask_is_ice(cellMask(iCell)))
              aCellOnEdgeHasNoIce = (aCellOnEdgeHasNoIce .or. (.not. li_mask_is_ice(cellMask(iCell))))
              aCellOnEdgeHasDynamicIce = (aCellOnEdgeHasDynamicIce .or. li_mask_is_dynamic_ice(cellMask(iCell)))
              aCellOnEdgeHasNoDynamicIce = (aCellOnEdgeHasNoDynamicIce .or. (.not. (li_mask_is_dynamic_ice(cellMask(iCell)))))
              aCellOnEdgeIsFloating = (aCellOnEdgeIsFloating .or. li_mask_is_floating_ice(cellMask(iCell)))
          end do
          if (aCellOnEdgeHasIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueIce)
          endif
          if (aCellOnEdgeHasDynamicIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueDynamicIce)
          endif
          if (aCellOnEdgeIsFloating) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueFloating)
          endif
          if (aCellOnEdgeHasIce .and. aCellOnEdgeHasNoIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueMargin)
          endif
          if (aCellOnEdgeHasDynamicIce .and. aCellOnEdgeHasNoDynamicIce) then
             edgeMask(i) = ior(edgeMask(i), li_mask_ValueDynamicMargin)
          endif

      end do

      ! vertexMask and edgeMask needs halo updates before they can be used.  Halo updates need to occur outside of block loops.  

   !--------------------------------------------------------------------
   end subroutine li_calculate_mask



   ! ===================================
   ! Functions for decoding bitmasks - will work with cellMask, edgeMask, or vertexMask
   ! ===================================
   ! Only adding the minimum needed for now.  These should be added as needed.


   function li_mask_is_ice_1d(mask)
      integer, dimension(:), pointer, intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_ice_1d

      li_mask_is_ice_1d = (iand(mask, li_mask_ValueIce) == li_mask_ValueIce)
   end function li_mask_is_ice_1d

   function li_mask_is_ice_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_ice_0d

      li_mask_is_ice_0d = (iand(mask, li_mask_ValueIce) == li_mask_ValueIce)
   end function li_mask_is_ice_0d



   function li_mask_is_dynamic_ice_1d(mask)
      integer, dimension(:), pointer, intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_dynamic_ice_1d

      li_mask_is_dynamic_ice_1d = (iand(mask, li_mask_ValueDynamicIce) == li_mask_ValueDynamicIce)
   end function li_mask_is_dynamic_ice_1d

   function li_mask_is_dynamic_ice_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_dynamic_ice_0d

      li_mask_is_dynamic_ice_0d = (iand(mask, li_mask_ValueDynamicIce) == li_mask_ValueDynamicIce)
   end function li_mask_is_dynamic_ice_0d



   function li_mask_is_floating_ice_1d(mask)
      integer, dimension(:), pointer, intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_floating_ice_1d

      li_mask_is_floating_ice_1d = (iand(mask, li_mask_ValueFloating) == li_mask_ValueFloating)
   end function li_mask_is_floating_ice_1d

   function li_mask_is_floating_ice_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_floating_ice_0d

      li_mask_is_floating_ice_0d = (iand(mask, li_mask_ValueFloating) == li_mask_ValueFloating)
   end function li_mask_is_floating_ice_0d



   function li_mask_is_grounded_ice_1d(mask)
      integer, dimension(:), pointer, intent(in) :: mask
      logical, dimension(size(mask)) :: li_mask_is_grounded_ice_1d

      li_mask_is_grounded_ice_1d = ( (iand(mask, li_mask_ValueFloating) /= li_mask_ValueFloating) &
                                .and. (li_mask_is_ice(mask)) )
   end function li_mask_is_grounded_ice_1d

   function li_mask_is_grounded_ice_0d(mask)
      integer, intent(in) :: mask
      logical :: li_mask_is_grounded_ice_0d

      li_mask_is_grounded_ice_0d = ( (iand(mask, li_mask_ValueFloating) /= li_mask_ValueFloating) &
                                .and. (li_mask_is_ice(mask)) )
   end function li_mask_is_grounded_ice_0d



!***********************************************************************
! Private subroutines:
!***********************************************************************

! - no private subroutines - (module is not declared private)


end module li_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

