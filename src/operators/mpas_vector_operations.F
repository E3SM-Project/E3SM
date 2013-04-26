! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_vector_operations
!
!> \brief MPAS vector operations
!> \author Mark Petersen
!> \date    25 September 2012
!> \details
!>  This module contains the routines involving vector operations
!
!-----------------------------------------------------------------------

module mpas_vector_operations

   use mpas_grid_types
   use mpas_constants

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

   public :: mpas_initialize_vectors, &
             mpas_unit_vec_in_r3, &
             mpas_cross_product_in_r3, &
             mpas_tangential_velocity, &
             mpas_vector_3DCell_to_2DEdge

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine mpas_unit_vec_in_r3
!
!> \brief   MPAS 3D unit vector routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine creates a unit vector out of an input point.
!-----------------------------------------------------------------------
  subroutine mpas_unit_vec_in_r3(xin)!{{{
    implicit none
    real (kind=RKIND), intent(inout) :: xin(3) !< Input/Output: Vector and unit vector
    real (kind=RKIND) :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
  end subroutine mpas_unit_vec_in_r3!}}}

!***********************************************************************
!
!  routine mpas_cross_product_in_r3
!
!> \brief   MPAS 3D cross product routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes the cross product of two input vectors.
!-----------------------------------------------------------------------
  subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)!{{{
    real (kind=RKIND), intent(in) :: p_1 (3) !< Input: Vector 1
    real (kind=RKIND), intent(in) :: p_2 (3) !< Input: Vector 2
    real (kind=RKIND), intent(out) :: p_out (3) !< Output: Cross product of vector 1 and vector 2

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine mpas_cross_product_in_r3!}}}

!***********************************************************************
!
!  routine mpas_vector_3DCell_to_2DEdge
!
!> \brief   Convert a 3D cell-centered vector field to 2D vectors at edges
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Convert a 3D cell-centered vector field to 2D vectors at edges, where 
!>  the local coordinate system at the edge is normal and tangent to that edge.
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_3DCell_to_2DEdge(grid, vector3DCell, vector2DEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vector3DCell  !< Input: 3-vector located at cell centers, in x,y,z coordinates

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
         vector2DEdge   !< Output: vector at edge, where 
                        !< index 1=1 is for normal component to edge
                        !< index 1=2 is for tangential component to edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, cell1, cell2, k, nVertLevels

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(3) :: vector3DEdge
      real(kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors, edgeTangentVectors

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      edgeNormalVectors => grid % edgeNormalVectors % array
      edgeTangentVectors => grid % edgeTangentVectors % array
      cellsOnEdge       => grid % cellsOnEdge % array

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels

           ! average neighboring cell-centered vectors to the edge
           vector3DEdge(:) = 0.5*(vector3DCell(:,k,cell1) + vector3DCell(:,k,cell2))

           ! components at edge: take dot products with unit vectors at edge
           vector2DEdge(1,k,iEdge) = sum(vector3DEdge(:)*edgeNormalVectors(:,iEdge))
           vector2DEdge(2,k,iEdge) = sum(vector3DEdge(:)*edgeTangentVectors(:,iEdge))      
 
         enddo
      enddo

   end subroutine mpas_vector_3DCell_to_2DEdge!}}}

!***********************************************************************
!
!  routine mpas_tangential_velocity
!
!> \brief   compute tangential velocity at an edge from the normal velocities
!> \author  Mark Petersen
!> \date    17 April 2013
!> \details 
!>  Compute tangential velocity at an edge from the normal velocities on the 
!>  edges of the two neighboring cells.
!
!-----------------------------------------------------------------------

   subroutine mpas_tangential_velocity(grid, normalVelocity, tangentialVelocity)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity   !< Input: Horizontal velocity normal to edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tangentialVelocity   !< Output: Horizontal velocity tangent to edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges,i,k,nVertLevels, eoe
      integer, dimension(:), pointer :: nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnEdge

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels

      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      weightsOnEdge     => grid % weightsOnEdge % array

      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k=1,nVertLevels
               tangentialVelocity(k,iEdge) = tangentialVelocity(k,iEdge) + weightsOnEdge(i,iEdge) * normalVelocity(k, eoe)
            end do
         end do
      end do

   end subroutine mpas_tangential_velocity!}}}

!***********************************************************************
!
!  routine mpas_initialize_vectors
!
!> \brief   MPAS RBF interpolation initialization routine
!> \author  Xylar Asay-Davis
!> \date    03/28/13
!> \details 
!> This routine computes geometric fields that will be potentially useful for calling
!> the interpolation routines.
!> Input: the grid
!> Output:
!>      edgeNormalVector - the unit vector normal to the edge and tangent to the sphere
!>      edgeTangentVector - the unit vector tangent to the edge and tangent to the sphere
!>      cellTangentPlane - 2 orthogonal unit vectors in the tangent plane of each cell
!>                         The first unit vector is chosen to point toward the center of the first
!>                         edge on the cell.
!>      localVerticalUnitVectors - the unit normal vector of the tangent plane at the center 
!>                         of each cell
!-----------------------------------------------------------------------
  subroutine mpas_initialize_vectors(grid)!{{{

    implicit none

    type (mesh_type), intent(inout) :: grid  !< Input/Output: Grid information

    integer :: nCells, nEdges
    integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell
    integer :: iEdge, iCell, cell1, cell2, vertex1, vertex2
    real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex
    real(kind=RKIND), dimension(:,:), pointer :: localVerticalUnitVectors, edgeNormalVectors, edgeTangentVectors
    real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
    real(kind=RKIND), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(kind=RKIND) :: normalDotRHat
    logical :: on_a_sphere

    xCell       => grid % xCell % array
    yCell       => grid % yCell % array
    zCell       => grid % zCell % array
    xEdge       => grid % xEdge % array
    yEdge       => grid % yEdge % array
    zEdge       => grid % zEdge % array
    xVertex     => grid % xVertex % array
    yVertex     => grid % yVertex % array
    zVertex     => grid % zVertex % array
    cellsOnEdge => grid % cellsOnEdge % array
    verticesOnEdge => grid % verticesOnEdge % array
    edgesOnCell => grid % edgesOnCell % array
    nCells      = grid % nCells
    nEdges      = grid % nEdges
    on_a_sphere = grid % on_a_sphere

    localVerticalUnitVectors => grid % localVerticalUnitVectors % array
    edgeNormalVectors => grid % edgeNormalVectors % array
    edgeTangentVectors => grid % edgeTangentVectors % array
    cellTangentPlane => grid % cellTangentPlane % array

    ! init arrays
    edgeNormalVectors = 0
    edgeTangentVectors = 0
    localVerticalUnitVectors = 0

    ! loop over all cells to be solved on this block
    do iCell=1,nCells
      if(on_a_sphere) then
        localVerticalUnitVectors(1,iCell) = xCell(iCell)
        localVerticalUnitVectors(2,iCell) = yCell(iCell)
        localVerticalUnitVectors(3,iCell) = zCell(iCell)
        call mpas_unit_vec_in_r3(localVerticalUnitVectors(:,iCell))
      else ! on a plane
        localVerticalUnitVectors(:,iCell) = (/ 0., 0., 1. /)
      end if
    end do

    ! Initialize normal unit vectors at each edge
    ! These vectors point from cell to cell.
    ! At boundaries, one cell does not exist, so it points from cell to edge.
    do iEdge = 1,nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)

      if (cell1 == nCells+1) then ! this is a boundary edge
        ! the normal points from the cell location to the edge location
        edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell2)
        edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell2)
        edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell2)

      elseif (cell2 == nCells+1) then ! this is a boundary edge
        ! the normal points from the cell location to the edge location
        edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell1)
        edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell1)
        edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)

      else ! this is not a boundary cell
        ! the normal points from the cell 1 to cell2
        ! mrp problem: on periodic domains, vectors on edges of domain point the wrong way.
        edgeNormalVectors(1,iEdge) = xCell(cell2) - xCell(cell1)
        edgeNormalVectors(2,iEdge) = yCell(cell2) - yCell(cell1)
        edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)

      endif
      call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
    end do

    ! Initialize tangent unit vectors at each edge
    do iEdge = 1,nEdges
      vertex1 = verticesOnEdge(1,iEdge) 
      vertex2 = verticesOnEdge(2,iEdge) 
      ! the tangent vector points from the vertex 1 to vertex 2
      edgeTangentVectors(1,iEdge) = xVertex(vertex2) - xVertex(vertex1)
      edgeTangentVectors(2,iEdge) = yVertex(vertex2) - yVertex(vertex1)
      edgeTangentVectors(3,iEdge) = zVertex(vertex2) - zVertex(vertex1)
      call mpas_unit_vec_in_r3(edgeTangentVectors(:,iEdge))
    end do

    do iCell=1,nCells
      iEdge = edgesOnCell(1,iCell)
      ! xHat and yHat are a local basis in the plane of the horizontal cell
      ! we arbitrarily choose xHat to point toward the first edge
      rHat = localVerticalUnitVectors(:,iCell)
      normalDotRHat = sum(edgeNormalVectors(:,iEdge)*rHat)
      xHatPlane = edgeNormalVectors(:,iEdge) - normalDotRHat*rHat
      call mpas_unit_vec_in_r3(xHatPlane)

      call mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane)
      call mpas_unit_vec_in_r3(yHatPlane) ! just to be sure...
      cellTangentPlane(:,1,iCell) = xHatPlane
      cellTangentPlane(:,2,iCell) = yHatPlane
    end do

  end subroutine mpas_initialize_vectors!}}}


end module mpas_vector_operations


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

