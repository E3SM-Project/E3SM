! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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
!> \date   April 2013
!> \details
!>  This module contains the routines involving vector operations
!
!-----------------------------------------------------------------------
module mpas_vector_operations

   use mpas_derived_types
   use mpas_pool_routines
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

   public :: mpas_vec_mag_in_r3, &
             mpas_initialize_vectors, &
             mpas_initialize_tangent_vectors, &
             mpas_unit_vec_in_r3, &
             mpas_cross_product_in_r3, &
             mpas_tangential_velocity, &
             mpas_tangential_vector_1d, &
             mpas_vector_R3Cell_to_2DEdge, &
             mpas_vector_R3Cell_to_normalVectorEdge, &
             mpas_vector_R3_to_LonLatR, &
             mpas_vector_LonLatR_to_R3, &
             mpas_zonal_meridional_vectors, &
             mpas_fix_periodicity, &
             mpas_unit_test_fix_periodicity


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine mpas_vec_mag_in_r3
!
!> \brief   MPAS 3D vector magnitude routine
!> \author  Matt Hoffman
!> \date    13 Jan. 2015
!> \details
!> This routine calculates the magnitude of a 3d vector.
!-----------------------------------------------------------------------
  real (kind=RKIND) function mpas_vec_mag_in_r3(xin)!{{{
    implicit none
    real (kind=RKIND), dimension(3), intent(in) :: xin !< Input: Vector
    mpas_vec_mag_in_r3 = sqrt(xin(1)**2 + xin(2)**2 + xin(3)**2)
  end function mpas_vec_mag_in_r3!}}}

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
    real (kind=RKIND), dimension(3), intent(inout) :: xin !< Input/Output: Vector and unit vector
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
!  routine mpas_vector_R3Cell_to_2DEdge
!
!> \brief   Convert an R3 cell-centered vector field to 2D vectors at edges
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Convert an R3 cell-centered vector field to 2D vectors at edges, where
!>  the local coordinate system at the edge is normal and tangent to that edge.
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3Cell_to_2DEdge(vectorR3Cell, &
      meshPool, edgeTangentVectors, includeHalo, normalVectorEdge, tangentialVectorEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vectorR3Cell  !< Input: 3-vector located at cell centers, in x,y,z coordinates

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: Mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         edgeTangentVectors   !< Input: unit vector tangent to an edge

      logical, intent(in) :: &
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         normalVectorEdge   !< Output: normal component of vector at edge

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tangentialVectorEdge   !< Output: tangential component of vector at edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2, k
      integer, pointer :: nEdgesCompute, nVertLevels

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(3) :: vectorR3Edge
      real(kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      if (includeHalo) then
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesCompute)
      else
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesCompute)
      endif
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      do iEdge=1,nEdgesCompute
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels

           ! average neighboring cell-centered vectors to the edge
           vectorR3Edge(:) = 0.5_RKIND*(vectorR3Cell(:,k,cell1) + vectorR3Cell(:,k,cell2))

           ! normal component at edge: take dot products with unit vectors at edge
           normalVectorEdge(k,iEdge) = sum(vectorR3Edge(:)*edgeNormalVectors(:,iEdge))
           tangentialVectorEdge(k,iEdge) = sum(vectorR3Edge(:)*edgeTangentVectors(:,iEdge))

         enddo
      enddo

   end subroutine mpas_vector_R3Cell_to_2DEdge!}}}

!***********************************************************************
!
!  routine mpas_vector_R3Cell_to_normalVectorEdge
!
!> \brief   Interpolate from R3 cell-centered vector field vector normal to edge
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Convert an R3 cell-centered vector field to normal vectors at edges.
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3Cell_to_normalVectorEdge(vectorR3Cell, &
      meshPool, includeHalo, normalVectorEdge)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         vectorR3Cell  !< Input: 3-vector located at cell centers, in x,y,z coordinates

      type (mpas_pool_type), intent(in) :: &
         meshPool         !< Input: Mesh information

      logical, intent(in) :: &
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         normalVectorEdge   !< Output: normal component of vector at edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2, k
      integer, pointer :: nEdgesCompute, nVertLevels

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(3) :: vectorR3Edge
      real(kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors

      if (includeHalo) then
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesCompute)
      else
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesCompute)
      endif
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      do iEdge=1,nEdgesCompute
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels

           ! average neighboring cell-centered vectors to the edge
           vectorR3Edge(:) = 0.5_RKIND*(vectorR3Cell(:,k,cell1) + vectorR3Cell(:,k,cell2))

           ! normal component at edge: take dot products with unit vectors at edge
           normalVectorEdge(k,iEdge) = sum(vectorR3Edge(:)*edgeNormalVectors(:,iEdge))

         enddo
      enddo

   end subroutine mpas_vector_R3Cell_to_normalVectorEdge!}}}

!***********************************************************************
!
!  routine mpas_tangential_velocity
!
!> \brief   compute tangential velocity at an edge from the normal velocities
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Compute tangential velocity at an edge from the normal velocities on the
!>  edges of the two neighboring cells.
!
!-----------------------------------------------------------------------

   subroutine mpas_tangential_velocity(normalVelocity, meshPool, includeHalo, tangentialVelocity)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity   !< Input: Horizontal velocity normal to edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: Mesh information

      logical, intent(in) :: &
         includeHalo !< Input: If true, halo cells and edges are included in computation

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

      integer :: iEdge, i, k, eoe
      integer, pointer :: nEdgesCompute, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnEdge

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      if (includeHalo) then
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesCompute)
      else
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesCompute)
      endif
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)

      do iEdge = 1,nEdgesCompute
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
!  routine mpas_tangential_vector_1d
!
!> \brief   compute tangential vector magnitude at an edge from the normal vector magnitudes
!> \author  Matt Hoffman/Mark Petersen
!> \date    December 2014
!> \details
!>  Compute tangential vector magnitude at an edge from the normal vector magnitude on the
!>  edges of the two neighboring cells.
!>  Based off of mpas_tangential_velocity().
!
!-----------------------------------------------------------------------

   subroutine mpas_tangential_vector_1d(normalVector, meshPool, includeHalo, tangentialVector)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         normalVector   !< Input: Horizontal vector magnitude normal to edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: Mesh information

      logical, intent(in) :: &
         includeHalo !< Input: If true, halo cells and edges are included in computation

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         tangentialVector   !< Output: Horizontal vector magnitude tangent to edge

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, i, k, eoe
      integer, pointer :: nEdgesCompute
      integer, dimension(:), pointer :: nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnEdge

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      if (includeHalo) then
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesCompute)
      else
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesCompute)
      endif

      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)

      tangentialVector = 0.0_RKIND
      do iEdge = 1,nEdgesCompute
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            tangentialVector(iEdge) = tangentialVector(iEdge) + weightsOnEdge(i,iEdge) * normalVector(eoe)
         end do
      end do

   end subroutine mpas_tangential_vector_1d!}}}


!***********************************************************************
!
!  routine mpas_zonal_meridional_vectors
!
!> \brief   Computes zonal, meridional, and vertical unit vectors
!> \author  Mark Petersen
!> \date    1 May 2013
!> \details
!>  Given a latitude and longitude location, compute unit vectors pointing
!>  in the zonal, meridional, and vertical directions.
!
!-----------------------------------------------------------------------

   subroutine mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         zonalUnitVector,     &!< Output: zonal unit vector
         meridionalUnitVector,&!< Output: meridional unit vector
         verticalUnitVector    !< Output: vertical unit vector

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND) :: sin_lat, sin_lon, cos_lat, cos_lon


      sin_lat = sin(lat)
      cos_lat = cos(lat)
      sin_lon = sin(lon)
      cos_lon = cos(lon)

      zonalUnitVector(1) = - sin_lon
      zonalUnitVector(2) =   cos_lon
      zonalUnitVector(3) =   0

      meridionalUnitVector(1) = - sin_lat * cos_lon
      meridionalUnitVector(2) = - sin_lat * sin_lon
      meridionalUnitVector(3) =   cos_lat

      verticalUnitVector(1) = cos_lat * cos_lon
      verticalUnitVector(2) = cos_lat * sin_lon
      verticalUnitVector(3) = sin_lat

   end subroutine mpas_zonal_meridional_vectors!}}}


!***********************************************************************
!
!  routine mpas_vector_R3_to_LonLatR
!
!> \brief   Convert an R3 vector to a vector in spherical coordinates
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Given an R3 vector, rotate it to spherical coordinates
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_R3_to_LonLatR(vectorR3, lon, lat, vectorLonLatR)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(in) :: &
         vectorR3  !< Input: vector in R3

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         vectorLonLatR   !< Output: vector in spherical coordinates

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      vectorLonLatR = 0.0_RKIND
      do i=1,3
         do j=1,3
            ! xi = R^T x
            ! where ^T is the transpose.  Note indices on R are (j,i)
            vectorLonLatR(i) = vectorLonLatR(i) + rotationMatrix(j,i)*vectorR3(j)
         enddo
      enddo

   end subroutine mpas_vector_R3_to_LonLatR!}}}


!***********************************************************************
!
!  routine mpas_vector_LonLatR_to_R3
!
!> \brief   Convert a vector in spherical coordinates to an R3 vector
!> \author  Mark Petersen
!> \date    April 2013
!> \details
!>  Given a vector in spherical coordinates, rotate it R3
!
!-----------------------------------------------------------------------

   subroutine mpas_vector_LonLatR_to_R3(vectorLonLatR, lon, lat, vectorR3)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(in) :: &
         vectorLonLatR   !< Input: vector in spherical coordinates

      real (kind=RKIND), intent(in) :: &
         lon, &!< Input: longitude, in radians, ranging [0,2*pi]
         lat   !< Input: latitude,  in radians, ranging [-pi,pi]

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(3), intent(out) :: &
         vectorR3  !< Output: vector in R3

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i,j

      real (kind=RKIND), dimension(3) :: zonalUnitVector, meridionalUnitVector, verticalUnitVector
      real (kind=RKIND), dimension(3,3) :: rotationMatrix

      call mpas_zonal_meridional_vectors(lon, lat, zonalUnitVector, meridionalUnitVector, verticalUnitVector)

      rotationMatrix(:,1) = zonalUnitVector
      rotationMatrix(:,2) = meridionalUnitVector
      rotationMatrix(:,3) = verticalUnitVector

      vectorR3 = 0.0_RKIND
      do i=1,3
         do j=1,3
            vectorR3(i) = vectorR3(i) + rotationMatrix(i,j)*vectorLonLatR(j)
         enddo
      enddo

   end subroutine mpas_vector_LonLatR_to_R3!}}}

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
!> Input: the mesh
!> Output:
!>      edgeNormalVectors - the unit vector normal to the edge and tangent to the sphere
!>      cellTangentPlane - 2 orthogonal unit vectors in the tangent plane of each cell
!>                         The first unit vector is chosen to point toward the center of the first
!>                         edge on the cell.
!>      localVerticalUnitVectors - the unit normal vector of the tangent plane at the center
!>                         of each cell
!-----------------------------------------------------------------------
  subroutine mpas_initialize_vectors(meshPool)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: meshPool  !< Input: Mesh information

    integer :: iEdge, iCell, cell1, cell2
    integer, pointer :: nCells, nEdges
    integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell
    real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge
    real(kind=RKIND), dimension(:,:), pointer :: localVerticalUnitVectors, edgeNormalVectors
    real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
    real(kind=RKIND), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(kind=RKIND) :: normalDotRHat
    logical, pointer :: on_a_sphere, is_periodic
    real(kind=RKIND), pointer :: x_period, y_period

    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)
    call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
    call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
    call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

    call mpas_pool_get_array(meshPool, 'localVerticalUnitVectors', localVerticalUnitVectors)
    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
    call mpas_pool_get_array(meshPool, 'cellTangentPlane', cellTangentPlane)

    ! init arrays
    edgeNormalVectors = 0
    localVerticalUnitVectors = 0

    ! loop over all cells to be solved on this block
    do iCell = 1, nCells
      if (on_a_sphere) then
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
        if (is_periodic) then
          edgeNormalVectors(1,iEdge) = mpas_fix_periodicity(xEdge(iEdge), xCell(cell2), x_period) - xCell(cell2)
          edgeNormalVectors(2,iEdge) = mpas_fix_periodicity(yEdge(iEdge), yCell(cell2), y_period) - yCell(cell2)
          edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell2)
        else
          edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell2)
          edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell2)
          edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell2)
        end if

      else if (cell2 == nCells+1) then ! this is a boundary edge
        ! the normal points from the cell location to the edge location
        if (is_periodic) then
          edgeNormalVectors(1,iEdge) = mpas_fix_periodicity(xEdge(iEdge), xCell(cell1), x_period) - xCell(cell1)
          edgeNormalVectors(2,iEdge) = mpas_fix_periodicity(yEdge(iEdge), yCell(cell1), y_period) - yCell(cell1)
          edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
        else
          edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell1)
          edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell1)
          edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
        end if

      else ! this is not a boundary cell
        ! the normal points from the cell 1 to cell2
        ! mrp problem: on periodic domains, vectors on edges of domain point the wrong way.
        if (is_periodic) then
          edgeNormalVectors(1,iEdge) = mpas_fix_periodicity(xCell(cell2), xCell(cell1), x_period) - xCell(cell1)
          edgeNormalVectors(2,iEdge) = mpas_fix_periodicity(yCell(cell2), yCell(cell1), y_period) - yCell(cell1)
          edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
        else
          edgeNormalVectors(1,iEdge) = xCell(cell2) - xCell(cell1)
          edgeNormalVectors(2,iEdge) = yCell(cell2) - yCell(cell1)
          edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
        end if

      end if
      call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
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

!***********************************************************************
!
!  routine mpas_initialize_tangent_vectors
!
!> \brief   Initialize edgeTangentVectors
!> \author  Mark Petersen
!> \date    June 2013
!> \details
!> This routine initializes the array edgeTangentVectors.  It is not
!> included in the mpas_initialize_vectors subroutine because the
!> array edgeTangentVectors is not included in all cores.
!> Input: the mesh
!> Output:
!>      edgeTangentVectors - the unit vector tangent to the edge and tangent to the sphere
!-----------------------------------------------------------------------
   subroutine mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         edgeTangentVectors   !< Output: unit vector tangent to an edge

      integer :: iEdge, vertex1, vertex2
      integer, pointer :: nEdges
      integer, dimension(:,:), pointer :: verticesOnEdge
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex

      logical, pointer :: is_periodic
      real(kind=RKIND), pointer :: x_period, y_period

      call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
      call mpas_pool_get_config(meshPool, 'x_period', x_period)
      call mpas_pool_get_config(meshPool, 'y_period', y_period)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

      do iEdge = 1,nEdges
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         ! the tangent vector points from the vertex 1 to vertex 2
         if (is_periodic) then
           edgeTangentVectors(1,iEdge) = mpas_fix_periodicity(xVertex(vertex2), xVertex(vertex1), x_period) - xVertex(vertex1)
           edgeTangentVectors(2,iEdge) = mpas_fix_periodicity(yVertex(vertex2), yVertex(vertex1), y_period) - yVertex(vertex1)
           edgeTangentVectors(3,iEdge) = zVertex(vertex2) - zVertex(vertex1)
         else
           edgeTangentVectors(1,iEdge) = xVertex(vertex2) - xVertex(vertex1)
           edgeTangentVectors(2,iEdge) = yVertex(vertex2) - yVertex(vertex1)
           edgeTangentVectors(3,iEdge) = zVertex(vertex2) - zVertex(vertex1)
         end if
         call mpas_unit_vec_in_r3(edgeTangentVectors(:,iEdge))
      end do

  end subroutine mpas_initialize_tangent_vectors!}}}

!***********************************************************************
!
!  subroutine mpas_fix_periodicity
!
!> \brief   Fixes periodicity of point pxi relative to a point xci and xiRef periodicity.
!> \author  Phillip Wolfram & Doug Jacobsen
!> \date    06/29/2015
!> \details
!>  This routine recomputes the location of a point pxi relative to a
!>  point xci with xiRef periodicity. The calculation ensures that the recomputed
!>  point and xci are spatially nearby, permitting reasonable calculation of
!>  relative geometry.  This function operates only on a single dimension "i" for
!>  a point pxi relative to a specific location xci for dimension xiRef
!>  in the periodic direction.  Note, pxi can only be adjusted by at most
!>  a single period, e.g., pxi = 260 + 360 = 620, xci = 0, and xiRef = 360
!>  returns 260, not -100.
!-----------------------------------------------------------------------
    real(KIND=RKIND) function mpas_fix_periodicity(pxi, xci, xiRef) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: pxi, xci, xiRef
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      ! real (kind=RKIND), intent(out) :: mpas_fix_periodicity

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: dist

      dist = pxi - xci

      if (abs(dist) > xiRef * 0.5_RKIND) then
        mpas_fix_periodicity = pxi - (dist/abs(dist)) * xiRef
      else
        mpas_fix_periodicity = pxi
      end if

    end function mpas_fix_periodicity !}}}

!***********************************************************************
!
!  subroutine mpas_unit_test_fix_periodicity
!
!> \brief   Simple unit test for testing the mpas_fix_periodicity routine.
!> \author  Phillip Wolfram
!> \date    06/29/2015
!> \details
!>  This routine tests the mpas_fix_periodicity routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_fix_periodicity(ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: ierr
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: x1, x2, xc, xLen, xnew
      real (kind=RKIND) :: eps = 1.0e-12_RKIND

      ierr = 0

      xLen = 2_RKIND*pii
      x1 = pii
      x2 = 3_RKIND*pii
      xc = pii

      xnew = mpas_fix_periodicity(x1, xc, xLen)
      if (abs(xnew - pii ) > eps) then
         write(stderrUnit,*) "Error in mpas_unit_test_fix_periodicity: ", &
           "x1's periodicity fix has error greater than tolerance -- FAILED."
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_fix_periodicity: x1 test - SUCCESS'
      endif

      xnew = mpas_fix_periodicity(x2, xc, xLen)
      if (abs(xnew - pii ) > eps) then
         write(stderrUnit,*) "Error in mpas_unit_test_fix_periodicity: ", &
           "x2's periodicity fix has error greater than tolerance -- FAILED."
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_fix_periodicity: x2 test - SUCCESS'
      endif

   end subroutine mpas_unit_test_fix_periodicity!}}}

end module mpas_vector_operations


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

