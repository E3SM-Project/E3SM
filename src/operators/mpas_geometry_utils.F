! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_geometry_utils

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_io_units
   use mpas_vector_operations
   use mpas_abort, only : mpas_dmpar_global_abort

   implicit none

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz

      real (kind=RKIND) :: a, b, c          ! Side lengths of spherical triangle ABC

      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC

      real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC

      real (kind=RKIND) :: s                ! Semiperimeter of the triangle
      real (kind=RKIND) :: sin_angle

      a = acos(max(min(bx*cx + by*cy + bz*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (3)
      b = acos(max(min(ax*cx + ay*cy + az*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (2)
      c = acos(max(min(ax*bx + ay*by + az*bz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (1)

      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az

      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az

      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)

      s = 0.5*(a + b + c)
!      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   ! Eqn. (28)
      sin_angle = sqrt(min(1.0_RKIND,max(0.0_RKIND,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)

      if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
         mpas_sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      else
         mpas_sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      end if

   end function mpas_sphere_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_PLANE_ANGLE
   !
   ! Computes the angle between vectors AB and AC, given points A, B, and C, and
   !   a vector (u,v,w) normal to the plane.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_plane_angle(ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w

      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: mAB              ! The magnitude of AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=RKIND) :: mAC              ! The magnitude of AC

      real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC

      real (kind=RKIND) :: cos_angle

      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
      mAB = sqrt(ABx**2.0 + ABy**2.0 + ABz**2.0)

      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
      mAC = sqrt(ACx**2.0 + ACy**2.0 + ACz**2.0)


      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)

      cos_angle = (ABx*ACx + ABy*ACy + ABz*ACz) / (mAB * mAC)

      if ((Dx*u + Dy*v + Dz*w) >= 0.0) then
         mpas_plane_angle =  acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      else
         mpas_plane_angle = -acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      end if

   end function mpas_plane_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_ARC_LENGTH
   !
   ! Returns the length of the great circle arc from A=(ax, ay, az) to
   !    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   !    same sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_arc_length(ax, ay, az, bx, by, bz)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz

      real (kind=RKIND) :: r, c
      real (kind=RKIND) :: cx, cy, cz

      cx = bx - ax
      cy = by - ay
      cz = bz - az

!      r = ax*ax + ay*ay + az*az
!      c = cx*cx + cy*cy + cz*cz
!
!      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

      r = sqrt(ax*ax + ay*ay + az*az)
      c = sqrt(cx*cx + cy*cy + cz*cz)
!      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
      mpas_arc_length = r * 2.0 * asin(c/(2.0*r))

   end function mpas_arc_length!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTine mpas_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpas_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz
      real (kind=RKIND), intent(out) :: cx, cy, cz

      real (kind=RKIND) :: r           ! Radius of the sphere
      real (kind=RKIND) :: d

      r = sqrt(ax*ax + ay*ay + az*az)

      cx = 0.5*(ax + bx)
      cy = 0.5*(ay + by)
      cz = 0.5*(az + bz)

      if (cx == 0. .and. cy == 0. .and. cz == 0.) then
         write(stderrUnit,*) 'Error: arc_bisect: A and B are diametrically opposite'
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if

   end subroutine mpas_arc_bisect!}}}


!***********************************************************************
!
!  routine mpas_distance_plane
!
!> \brief   Calculates distance between two points on a plane.
!> \author  Matthew Hoffman
!> \date    5 February 2015
!> \details
!>  This routine calculates distance between two points on a plane.
!>  It does not work for periodic meshes.
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_distance_plane(a, b)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b  !< Input: 3d (x,y,z) points between which to calculate distance

      mpas_distance_plane = sqrt(sum((a - b)**2))
   end function mpas_distance_plane !}}}


!***********************************************************************
!
!  routine mpas_distance
!
!> \brief   Calculates distance between two points on a plane or sphere
!> \author  Matthew Hoffman
!> \date    5 February 2015
!> \details
!>  This routine calculates distance between two points on a plane or sphere.
!>  It does not work for periodic meshes because mpas_distance_plane does not!
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_distance(a, b, on_a_sphere)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b  !< Input: 3d (x,y,z) points between which to calculate distance
      logical, intent(in) :: on_a_sphere  !< Input: If on a sphere

      if (on_a_sphere) then
         mpas_distance = mpas_arc_length(a(1), a(2), a(3), b(1), b(2), b(3))
      else
         mpas_distance = mpas_distance_plane(a, b)
      endif
   end function mpas_distance !}}}


   subroutine mpas_poly_fit_2(a_in,b_out,weights_in,m,n,ne)!{{{

      use mpas_matrix_operations

      implicit none

      integer, intent(in) :: m,n,ne
      real (kind=RKIND), dimension(ne,ne), intent(in) :: a_in, weights_in
      real (kind=RKIND), dimension(ne,ne), intent(out) :: b_out

      ! local storage

      real (kind=RKIND), dimension(m,n)  :: a
      real (kind=RKIND), dimension(n,m)  :: b
      real (kind=RKIND), dimension(m,m)  :: w,wt,h
      real (kind=RKIND), dimension(n,m)  :: at, ath
!      real (kind=RKIND), dimension(n,n)  :: ata, ata_inv, atha, atha_inv
      real (kind=RKIND), dimension(n,n)  :: ata, atha, atha_inv
      integer, dimension(n) :: indx
!      integer :: i,j

      if ( (ne < n) .or. (ne < m) ) then
         write(stderrUnit,*) ' error in poly_fit_2 inversion ',m,n,ne
         call mpas_dmpar_global_abort('ERROR: in subroutine poly_fit_2()')
      end if

!      a(1:m,1:n) = a_in(1:n,1:m)
      a(1:m,1:n) = a_in(1:m,1:n)
      w(1:m,1:m) = weights_in(1:m,1:m)
      b_out(:,:) = 0.

      wt = transpose(w)
      h = matmul(wt,w)
      at = transpose(a)
      ath = matmul(at,h)
      atha = matmul(ath,a)

      ata = matmul(at,a)

!      if (m == n) then
!         call mpas_migs(a,n,b,indx)
!      else

         call mpas_migs(atha,n,atha_inv,indx)

         b = matmul(atha_inv,ath)

!         call mpas_migs(ata,n,ata_inv,indx)
!         b = matmul(ata_inv,at)
!      end if
      b_out(1:n,1:m) = b(1:n,1:m)

!     do i=1,n
!        write(stderrUnit,*) ' i, indx ',i,indx(i)
!     end do
!
!     write(stderrUnit,*) ' '

   end subroutine mpas_poly_fit_2!}}}


!***********************************************************************
!
!  routine mpas_calculate_barycentric_weights
!
!> \brief   Calculates barycentric weights for a point relative to a triangle
!> \author  Matthew Hoffman
!> \date    8 January 2015
!> \details
!>  This routine calculates calculates barycentric weights (coordinates)
!>  for a point relative to three points provided that form a triangle.
!>  Note this does not handle triangles spanning planar periodic meshes because mpas_triangle_signed_area_plane does not!
!-----------------------------------------------------------------------
   subroutine mpas_calculate_barycentric_weights(point, a, b, c, meshPool, weight_a, weight_b, weight_c, ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: point  !< Input: 3d (x,y,z) point
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
      type (mpas_pool_type), intent(in) :: meshPool  !< Input: Mesh information
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), intent(out) :: weight_a, weight_b, weight_c  !< Output: The 3 bary weights
      integer, intent(out) :: ierr !< Output: error flag
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND) :: triangleArea, error
      integer, pointer :: vertexDegree

      ierr = 0

      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      if (vertexDegree /= 3) then
         write (stderrUnit,*) 'Error: Barycentric weights can only be calculated if vertexDegree is 3'
         ierr = 1
         return
      endif

      triangleArea = mpas_triangle_signed_area(a, b, c, meshPool)
      if (triangleArea > 0.0_RKIND) then
         weight_a = mpas_triangle_signed_area(point, b, c, meshPool) / triangleArea
         weight_b = mpas_triangle_signed_area(point, c, a, meshPool) / triangleArea
         weight_c = mpas_triangle_signed_area(point, a, b, meshPool) / triangleArea
      else
         weight_a = 0.0_RKIND
         weight_b = 0.0_RKIND
         weight_c = 0.0_RKIND
      endif

      error = abs(weight_a + weight_b + weight_c - 1.0_RKIND)
      if (error > 1.0e-12_RKIND) then
          write (stderrUnit,*) 'Error: Barycentric weights sum to substantially different from 1.0. error:', error  !, triangleArea, weight_a, weight_b, weight_c
         ierr = 1
      endif

! An alternate implementation:
! Reference: http://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
!      real(kind=RKIND), dimension(3), u_vec, v_vec, w_vec, n_vec, uxw, wxv
!      u_vec = b - a
!      v_vec = c - a
!      w_vec = point - a
!      call mpas_cross_product_in_r3(u_vec, v_vec, n_vec)
!      call mpas_cross_product_in_r3(u_vec, w_vec, uxw)
!      call mpas_cross_product_in_r3(w_vec, v_vec, wxv)
!      weight_c = sum(uxw * n_vec) / sum(n_vec * n_vec)
!      weight_b = sum(wxv * n_vec) / sum(n_vec * n_vec)
!      weight_a = 1.0_RKIND - weight_c - weight_b

   end subroutine mpas_calculate_barycentric_weights!}}}


!***********************************************************************
!
!  routine mpas_calculate_barycentric_weights_for_points
!
!> \brief   Calculates barycentric weights for a set of points
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates calculates barycentric weights (coordinates)
!>  for a set of points provided.  Requires a mesh with a triangular dual mesh.
!>  The routine will attempt to find which triangle (made of cell centers)
!>  in which each point lies.  Those cell center indices will be stored in the
!>  baryCellsOnPoints array.
!>  The triangle search is limited to a triangle location provided to the routine
!>  and its surrounding triangles.  The triangle location is identified by a vertex
!>  id, which is stored in the input array searchVertex.
!>  If the pointset is vertices, the list of vertices can be passed in.  If the
!>  pointset is something else, then the callling routine must
!>  first determine which vertices each point is closest (or at least close) to.
!>  If no 'owning' triangle is identified in this routine, then
!>  the closest valid triangle that was searched will be used, and therefore
!>  barycentric extrapolation will be set up.  The barycentric weights associated with
!>  the triangle nodes identified in baryCellsOnPoints are stored in the array
!>  baryWeightsOnPoints.
!>
!>  Note this does not handle triangles spanning planar periodic meshes because
!>  mpas_triangle_signed_area_plane does not!
!>
!>  Note that this routine will store local indices in baryCellsOnVertex so this
!>  field will differ across decompositions.  For this reason, baryCellsOnVertex
!>  should not be used as a restart variable.  Calculating baryCellsOnVertex during
!>  init for each simulation is the only robust solution.  Similarly, if you choose
!>  to output baryCellsOnVertex, be aware that you will be looking at local indices
!>  which may be of limited use.
!>
!>  To use this routine to build weights for interpolating from cells to vertices
!>  you can call this routine like:
!>     do iVertex = 1, nVertices
!>        vertexIndicesField % array(iVertex) = iVertex
!>     enddo
!>     call mpas_calculate_barycentric_weights_for_points(meshPool, &
!>            xVertex(1:nVertices), yVertex(1:nVertices), zVertex(1:nVertices), &
!>            vertexIndicesField % array(1:nVertices), &
!>            baryCellsOnVertex(:, 1:nVertices), baryWeightsOnVertex(:, 1:nVertices), err_tmp)
!>  Note that you should only pass in the indices up to nVertices in this case -
!>  this routine will not perform special treatment of the extra 'garbage' index at nVertices+1.
!>  Also note that vertexIndicesField%array is just the local indices 1:nVertices.
!-----------------------------------------------------------------------
   subroutine mpas_calculate_barycentric_weights_for_points(meshPool, xPoint, yPoint, zPoint, searchVertex, baryCellsOnPoints, baryWeightsOnPoints, ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(:), intent(in) :: xPoint, yPoint, zPoint  !< Input: coordinate of the points for which barycentric weights should be calculated
      integer, dimension(:), intent(in) :: searchVertex  !< Input: list of the vertex id to be used as the search origin for each point listed in xPoint, yPoint, XPoint
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: meshPool  !< Input: Mesh information
      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, dimension(:,:), intent(out) :: baryCellsOnPoints  !< Output: An array of dimension (3, nPoints) that contains the 3 neighboring cells to be used for the barycentric interpolation/extrapolation
      real(kind=RKIND), dimension(:,:), intent(out) :: baryWeightsOnPoints  !< Output: An array of dimension (3, nPoints) that contains the 3 bary weights for each point
      integer, intent(out) :: ierr !< Output: error flag
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer, pointer :: vertexDegree, nCells, nVertices
      logical, pointer :: on_a_sphere
      integer :: nPoints, iPoint, baseVertex
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnVertex, verticesOnCell
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
      real(kind=RKIND), dimension(3,3) :: triangleVertices
      integer :: cell1, cell2, cell3, c, t, iCell, iTriangle
      real(kind=RKIND), dimension(3) :: point, triangle_a, triangle_b, triangle_c
      logical :: in_a_triangle
      real(kind=RKIND) :: this_triangle_distance, nearest_triangle_distance
      integer :: nearest_triangle_index
      integer :: err_tmp

      ierr = 0
      err_tmp = 0

      nPoints = size(xPoint)

      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

      if (vertexDegree /= 3) then
         write (stderrUnit,*) 'Error: Barycentric weights can only be calculated if vertexDegree is 3'
         ierr = 1
         return
      endif

      do iPoint = 1, nPoints
         point = (/ xPoint(iPoint), yPoint(iPoint), zPoint(iPoint) /)
         baseVertex = searchVertex(iPoint)

         in_a_triangle = .false.
         nearest_triangle_distance = huge(nearest_triangle_distance)
         nearest_triangle_index = -1

         ! Step 1 - Identify the 'owning' triangle, if possible
         ! Start with the triangle defined by cellsOnVertex for the searchVertex.
         baryCellsOnPoints(:, iPoint) = cellsOnVertex(:, baseVertex)   ! start by assigning cellsOnVertex as the triangle
         cell1 = cellsOnVertex(1, baseVertex)
         cell2 = cellsOnVertex(2, baseVertex)
         cell3 = cellsOnVertex(3, baseVertex)
         if (cell1 <= nCells .and. cell2 <= nCells .and. cell3 <= nCells) then
            triangleVertices(1,:) = (/ xCell(cell1), yCell(cell1), zCell(cell1) /)
            triangleVertices(2,:) = (/ xCell(cell2), yCell(cell2), zCell(cell2) /)
            triangleVertices(3,:) = (/ xCell(cell3), yCell(cell3), zCell(cell3) /)
            in_a_triangle = mpas_point_in_polygon(point, triangleVertices, on_a_sphere)
            nearest_triangle_distance = mpas_distance(point, (/ xVertex(baseVertex), yVertex(baseVertex), zVertex(baseVertex) /), on_a_sphere)
            nearest_triangle_index = baseVertex
         else
            in_a_triangle = .false.
            nearest_triangle_distance = huge(this_triangle_distance)
            nearest_triangle_index = baseVertex
         endif

         if (.not. in_a_triangle) then
            ! If that was the owning triangle, don't do anything because we've already assigned cellsOnVertex as the triangle
            ! If that was not the owning triangle, check the neighboring triangles
            neighborCellsLoop: do c = 1, 3
               iCell = cellsOnVertex(c, baseVertex)
               if (iCell <= nCells) then
                  do t = 1, nEdgesOnCell(iCell)
                     iTriangle = verticesOnCell(t, iCell)
                     if (iTriangle <= nVertices) then
                        cell1 = cellsOnVertex(1, iTriangle)
                        cell2 = cellsOnVertex(2, iTriangle)
                        cell3 = cellsOnVertex(3, iTriangle)
                        if (cell1 <= nCells .and. cell2 <= nCells .and. cell3 <= nCells) then
                           triangleVertices(1,:) = (/ xCell(cell1), yCell(cell1), zCell(cell1) /)
                           triangleVertices(2,:) = (/ xCell(cell2), yCell(cell2), zCell(cell2) /)
                           triangleVertices(3,:) = (/ xCell(cell3), yCell(cell3), zCell(cell3) /)
                           if (mpas_point_in_polygon(point, triangleVertices, on_a_sphere)) then
                              in_a_triangle = .true.
                              baryCellsOnPoints(:, iPoint) = cellsOnVertex(:, iTriangle)
                              exit neighborCellsLoop ! No need to keep searching
                           endif
                           this_triangle_distance = mpas_distance(point, (/ xVertex(iTriangle), yVertex(iTriangle), zVertex(iTriangle) /), on_a_sphere)
                        else
                           this_triangle_distance = huge(this_triangle_distance)
                        endif ! check that all 3 cells are valid
                        if (this_triangle_distance < nearest_triangle_distance) then
                           nearest_triangle_distance = this_triangle_distance
                           nearest_triangle_index = iTriangle
                        endif
                     endif ! check if triangle is valid
                  end do
               endif
            end do neighborCellsLoop
         endif

         if (.not. in_a_triangle) then
            ! If none of the neighboring triangles are the owning triangle, then
            ! here we revert back to the closest valid triangle (which will be extrapolation)
            baryCellsOnPoints(:, iPoint) = cellsOnVertex(:, nearest_triangle_index)
         endif


         ! Step 2 - Determine the barycentric weights for the identified 'owning' triangle
         cell1 = baryCellsOnPoints(1, iPoint)
         cell2 = baryCellsOnPoints(2, iPoint)
         cell3 = baryCellsOnPoints(3, iPoint)

         triangle_a = (/ xCell(cell1), yCell(cell1), zCell(cell1) /)
         triangle_b = (/ xCell(cell2), yCell(cell2), zCell(cell2) /)
         triangle_c = (/ xCell(cell3), yCell(cell3), zCell(cell3) /)

         call mpas_calculate_barycentric_weights(point, triangle_a, triangle_b, triangle_c, meshPool,     &
             baryWeightsOnPoints(1, iPoint), baryWeightsOnPoints(2, iPoint), baryWeightsOnPoints(3, iPoint), err_tmp)
         ierr = ior(ierr, err_tmp)
      enddo

   end subroutine mpas_calculate_barycentric_weights_for_points !}}}


!***********************************************************************
!
!  routine mpas_triangle_signed_area
!
!> \brief   Calculates area of a triangle, whether on a sphere or a plane.
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle whether on a sphere or a plane.
!>  Note this does not handle triangles spanning planar periodic meshes because mpas_triangle_signed_area_plane does not!
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area(a, b, c, meshPool)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to get the area
      type (mpas_pool_type), intent(in) :: meshPool  !< Input: Mesh information (to find out if on a sphere)
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: on_a_sphere
      real(kind=RKIND), dimension(3) :: normalvec
      real(kind=RKIND), pointer :: radius

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

      if (on_a_sphere) then
         mpas_triangle_signed_area = mpas_triangle_signed_area_sphere(a, b, c, radius)
      else
         normalvec = (/ 0, 0, 1 /)
         mpas_triangle_signed_area = mpas_triangle_signed_area_plane(a, b, c, normalvec)
      endif
   end function mpas_triangle_signed_area !}}}


!***********************************************************************
!
!  routine mpas_triangle_signed_area_plane
!
!> \brief   Calculates signed area of a triangle in a plane
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle in a plane.
!>  Uses cross product.  Signed area will be positive if the vertices are oriented counterclockwise.
!>  Note this does not handle triangles spanning periodic meshes!
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area_plane(a, b, c, normalvec)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to calculate the area
      real(kind=RKIND), dimension(3), intent(in) :: normalvec  !< Input: 3d vector indicating the normal direction for the plane for assigning a sign to the area
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3) :: ab, ac, crossprod, triangleNormal

      ab = b - a
      ac = c - a
      call mpas_cross_product_in_r3(ab, ac, crossprod)
      if (mpas_vec_mag_in_r3(crossprod) == 0.0_RKIND) then
         mpas_triangle_signed_area_plane = 0.0_RKIND
      else
         triangleNormal = crossprod / mpas_vec_mag_in_r3(crossprod)
         mpas_triangle_signed_area_plane = 0.5_RKIND * (mpas_vec_mag_in_r3(crossprod)) *  &
              sum(triangleNormal * normalvec)
      endif
   end function mpas_triangle_signed_area_plane !}}}


!***********************************************************************
!
!  routine mpas_triangle_signed_area_sphere
!
!> \brief   Calculates area of a triangle on a sphere
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle on the surface of a sphere.
!>  Uses the spherical analog of Heron's formula.
!>  Copied from mesh generator.  A CCW winding angle is positive.
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area_sphere(a, b, c, radius)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
      real(kind=RKIND), intent(in) :: radius  !< sphere radius
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND) :: ab, bc, ca, semiperim, tanqe
      real(kind=RKIND), dimension(3) :: ablen, aclen, Dlen

      ab = mpas_arc_length(a(1), a(2), a(3), b(1), b(2), b(3))/radius
      bc = mpas_arc_length(b(1), b(2), b(3), c(1), c(2), c(3))/radius
      ca = mpas_arc_length(c(1), c(2), c(3), a(1), a(2), a(3))/radius
      semiperim = 0.5 * (ab + bc + ca)

      tanqe = sqrt(max(0.0_RKIND,tan(0.5_RKIND * semiperim) * tan(0.5_RKIND * (semiperim - ab)) &
                   * tan(0.5_RKIND * (semiperim - bc)) * tan(0.5_RKIND * (semiperim - ca))))

      mpas_triangle_signed_area_sphere = 4.0_RKIND * radius * radius * atan(tanqe)

      ! computing correct signs (in similar fashion to mpas_sphere_angle)
      ablen(1) = b(1) - a(1)
      ablen(2) = b(2) - a(2)
      ablen(3) = b(3) - a(3)

      aclen(1) = c(1) - a(1)
      aclen(2) = c(2) - a(2)
      aclen(3) = c(3) - a(3)

      dlen(1) =   (ablen(2) * aclen(3)) - (ablen(3) * aclen(2))
      dlen(2) = -((ablen(1) * aclen(3)) - (ablen(3) * aclen(1)))
      dlen(3) =   (ablen(1) * aclen(2)) - (ablen(2) * aclen(1))

      if ((Dlen(1)*a(1) + Dlen(2)*a(2) + Dlen(3)*a(3)) < 0.0) then
        mpas_triangle_signed_area_sphere = -mpas_triangle_signed_area_sphere
      end if

   end function mpas_triangle_signed_area_sphere !}}}


!***********************************************************************
!
!  routine mpas_point_in_polygon
!
!> \brief   Hit test to determine if a point is inside of a polygon
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine determines if a point is inside of a polygon.
!>  This is difficult because floating point arithmetic prevents a precise
!>  determination.  A tolerance is used to allow the point to be within the
!>  the polygon within some tolerance.  This means it is possible for a point
!>  to be identified to be within multiple polygons.  However, it avoids the
!>  situation where a point on a edge could be 'orphaned' - determined
!>  to not belong to *any* polygons.
!-----------------------------------------------------------------------
   logical function mpas_point_in_polygon(point, polygonVertices, on_a_sphere)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: point  !< Input: 3d (x,y,z) point
      real(kind=RKIND), dimension(:,:), intent(in) :: polygonVertices  !< Input: 3d (x,y,z) points forming the polygon to test, second dimension should be 3
      logical, intent(in) :: on_a_sphere  !< Input: If on a sphere
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3) :: normal_vector, crossprod, vec1, vec2
      integer :: polygonDegree, i
      integer, dimension(:), allocatable :: vertexNeighborFwd
      real(kind=RKIND), parameter :: eps = 1.0e-12_RKIND

      if (on_a_sphere) then
         normal_vector = point
      else
         normal_vector = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      endif

      polygonDegree = size(polygonVertices, 1)
      allocate(vertexNeighborFwd(polygonDegree))
      vertexNeighborFwd  = (/ (i+1, i = 1, polygonDegree) /)
      vertexNeighborFwd(polygonDegree) = 1

      mpas_point_in_polygon = .true.
      do i = 1, polygonDegree
         vec1 = polygonVertices(vertexNeighborFwd(i),:) - polygonVertices(i,:)
         vec2 = point - polygonVertices(i,:)
         call mpas_cross_product_in_r3(vec1, vec2, crossprod)
         if (sum(crossprod * normal_vector) < (0.0_RKIND - eps)) then
            mpas_point_in_polygon = .false.
            exit  ! If the point is ouside one of the edges, then we need not look further.
         endif
      enddo

      deallocate(vertexNeighborFwd)

   end function mpas_point_in_polygon !}}}


!***********************************************************************
!
!  subroutine mpas_cells_to_points_using_baryweights
!
!> \brief   Converts a single layer scalar field from cells to specified points using barycentric weights.
!> \author  Matt Hoffman
!> \date    14 Jan 2015
!> \details
!>  This routine converts a single layer scalar field from cells to a provided set
!>  of point locations using barycentric weights.
!>  The weights should be calculated on init using the
!>  mpas_calculate_barycentric_weights_for_points routine above.  That routine
!>  will calculate the baryCellsOnPoints, baryWeightsOnPoints fields.  These may
!>  be stored in meshPool, but that is not assumed.
!>
!>  To use this routine to interpolate from cells to vertices you can call this
!>  routine like:
!>    call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
!>       baryWeightsOnVertex(:, 1:nVertices), upperSurface, upperSurfaceVertex(1:nVertices), err_tmp)
!>  Note that you should only pass in the indices up to nVertices for the fields
!>  with a nVertices dimension in this case -
!>  this routine will not perform special treatment of the extra 'garbage' index at nVertices+1.
!-----------------------------------------------------------------------
   subroutine mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnPoints, baryWeightsOnPoints, fieldCells, fieldPoints, ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      real (kind=RKIND), dimension(:), intent(in) :: &
         fieldCells    !< Input: field on cells
      integer, dimension(:,:), intent(in) :: &
         baryCellsOnPoints    !< Input: The three cells that should be used for barycentric interpolation (or possibly extrapolation) for each point.  First dimension should be 3.
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         baryWeightsOnPoints    !< Input: the weights corresponding to each cell in baryCellsOnPoints.  First dimension should be 3.
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(out) :: &
         fieldPoints    !< Output: field on points
      integer, intent(out) :: ierr
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer, pointer :: vertexDegree
      integer :: nPoints, iPoint, cell1, cell2, cell3

      ierr = 0

      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      if (vertexDegree /= 3) then
         write (stderrUnit,*) 'Error: Barycentric weights can only be calculated if vertexDegree is 3'
         ierr = 1
         return
      endif

      nPoints = size(fieldPoints)

      do iPoint = 1, nPoints  ! Loop over vertices
        cell1 = baryCellsOnPoints(1, iPoint)
        cell2 = baryCellsOnPoints(2, iPoint)
        cell3 = baryCellsOnPoints(3, iPoint)
        fieldPoints(iPoint) = baryWeightsOnPoints(1, iPoint) * fieldCells(cell1) +    &
                              baryWeightsOnPoints(2, iPoint) * fieldCells(cell2) +    &
                              baryWeightsOnPoints(3, iPoint) * fieldCells(cell3)
      enddo

   end subroutine mpas_cells_to_points_using_baryweights !}}}

!***********************************************************************
!
!  subroutine mpas_unit_test_triangle_signed_area_sphere
!
!> \brief   Simple unit test for testing the mpas_triangle_signed_area_sphere
!> \author  Phillip J. Wolfram
!> \date    08/05/2015
!> \details
!>  This routine tests the mpas_triangle_signed_area_sphere routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_triangle_signed_area_sphere(ierr)!{{{
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
      real (kind=RKIND), dimension(3) :: a, b, c
      real (kind=RKIND) :: sphereRadius
      real (kind=RKIND) :: eps = 1.0e-12_RKIND

      ierr = 0

      a = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      b = (/ 0.0_RKIND, 1.0_RKIND, 0.0_RKIND /)
      c = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      sphereRadius = 1.0_RKIND

      ! 1/8 of spherical area = pi*r*r/8 = pi/2 for r=1.0
      if (abs((mpas_triangle_signed_area_sphere(a,b,c,sphereRadius) - pii/2.0_RKIND)/(pii/2.0_RKIND)) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_triangle_signed_area_sphere: unit sphere area wrong.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_triangle_signed_area_sphere: unit sphere area ok - SUCCESS'
      endif

      ! 1/8 of spherical area = pi*r*r/8 = pi/2 for r=1.0
      if (abs((mpas_triangle_signed_area_sphere(a,c,b,sphereRadius) + pii/2.0_RKIND)/(pii/2.0_RKIND)) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_triangle_signed_area_sphere: unit sphere area sign wrong.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_triangle_signed_area_sphere: unit sphere area sign ok - SUCCESS'
      endif

      ! 1/8 of spherical area = pi*r*r/8 = r*r*pi/2 for r=sphereRadius
      sphereRadius = 6371220._RKIND
      a = a*sphereRadius
      b = b*sphereRadius
      c = c*sphereRadius
      if (abs((mpas_triangle_signed_area_sphere(a,b,c,sphereRadius) - pii/2.0_RKIND*sphereRadius*sphereRadius)/(pii/2.0_RKIND*sphereRadius*sphereRadius)) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_triangle_signed_area_sphere: radius sphere area wrong.'
         write(stderrUnit,*) mpas_triangle_signed_area_sphere(a,b,c,sphereRadius) - pii/2.0_RKIND*sphereRadius*sphereRadius
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_triangle_signed_area_sphere: radius sphere area ok - SUCCESS'
      endif


   end subroutine mpas_unit_test_triangle_signed_area_sphere!}}}

!***********************************************************************
!
!  subroutine mpas_unit_test_in_triangle
!
!> \brief   Simple unit test for testing the mpas_point_in_polygon routine.
!> \author  Matt Hoffman
!> \date    4 Feb 2015
!> \details
!>  This routine tests the mpas_point_in_polygon routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_in_triangle(ierr)!{{{
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
      real (kind=RKIND), dimension(3,3) :: trianglePoints
      real (kind=RKIND), dimension(3) :: point1, point2

      ierr = 0

      trianglePoints(1,:) = (/ 0.0, 0.0, 0.0 /)
      trianglePoints(2,:) = (/ 1.0, 0.0, 0.0 /)
      trianglePoints(3,:) = (/ 0.0, 1.0, 0.0 /)

      point1 = (/ 0.333, 0.333, 0.333 /)
      point2 = (/ -0.01, 0.5, 0.0 /)

      if (.not. mpas_point_in_polygon(point1, trianglePoints, .false.)) then
         write(stderrUnit,*) 'Error in mpas_unit_test_in_triangle: point1 calculated to be outside of triangle.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_in_triangle: point1 test - SUCCESS'
      endif

      if (mpas_point_in_polygon(point2, trianglePoints, .false.)) then
         write(stderrUnit,*) 'Error in mpas_unit_test_in_triangle: point2 calculated to be inside of triangle.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_in_triangle: point2 test - SUCCESS'
      endif

   end subroutine mpas_unit_test_in_triangle !}}}



!***********************************************************************
!
!  subroutine mpas_unit_test_bary_weights
!
!> \brief   Simple unit test for testing the mpas_calculate_barycentric_weights routine.
!> \author  Matt Hoffman
!> \date    4 Feb 2015
!> \details
!>  This routine tests the mpas_calculate_barycentric_weights routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_bary_weights(meshPool, ierr)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool  !< Input: Mesh information
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
      real (kind=RKIND), dimension(3,3) :: trianglePoints
      real (kind=RKIND), dimension(3) :: point1, point2
      real (kind=RKIND), dimension(3) :: weights
      real (kind=RKIND) :: eps = 1.0e-12_RKIND

      ierr = 0

      trianglePoints(1,:) = (/ 0.0, 0.0, 0.0 /)
      trianglePoints(2,:) = (/ 1.0, 0.0, 0.0 /)
      trianglePoints(3,:) = (/ 0.0, 1.0, 0.0 /)

      point1 = (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /)
      point2 = (/ 0.0, 0.0, 0.0 /)


      call mpas_calculate_barycentric_weights(point1, trianglePoints(1,:), trianglePoints(2,:), trianglePoints(3,:), meshPool,     &
             weights(1), weights(2), weights(3), ierr)
      if (maxval(abs(weights - (/ 0.33333333333333333, 0.33333333333333333, 0.33333333333333333 /) )) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_bary_weights: point1 weights have error greater than tolerance.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_bary_weights: point1 test - SUCCESS'
      endif

      call mpas_calculate_barycentric_weights(point2, trianglePoints(1,:), trianglePoints(2,:), trianglePoints(3,:), meshPool,     &
             weights(1), weights(2), weights(3), ierr)
      if (maxval(abs(weights - (/ 1.0, 0.0, 0.0 /) )) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_bary_weights: point1 weights have error greater than tolerance.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_bary_weights: point2 test - SUCCESS'
      endif

   end subroutine mpas_unit_test_bary_weights !}}}


!***********************************************************************
!
!  subroutine mpas_get_nearby_cell_index
!
!> \brief   Determine nearest cell ID for a particular location
!> \author  Phillip Wolfram
!> \date    01/21/2015
!> \details
!>  This routine returns the cell index for a particular location based on
!>  distance to a cell center.  Two approaches are employed: 1) a brute
!>  force solution where the minimum distance over all cells is
!>  determined, and 2) a greedy, nearest neighbor algorithm which
!>  searches over all the neighbors of the lastCell index and finds the
!>  nearest point.
!-----------------------------------------------------------------------
   subroutine mpas_get_nearby_cell_index(nCells, xc,yc,zc , xp,yp,zp, meshPool, lastCell, cellsOnCell, nEdgesOnCell) !{{{
      implicit none

      integer, intent(in) :: nCells                            ! number of cells
      integer, dimension(:,:), intent(in) :: cellsOnCell       ! cell connectivity
      integer, dimension(:), intent(in) :: nEdgesOnCell        ! number of edges bordering cell
      real (kind=RKIND), dimension(:), intent(in) :: xc,yc,zc  ! center cell locations
      real (kind=RKIND), intent(in) :: xp,yp,zp                ! point to find cell locations
      type (mpas_pool_type), intent(in), pointer :: meshPool   ! meshPool pointer

      real (kind=RKIND) :: pointRadius   ! magnitude of a point radius
      integer :: iPoint, aPoint, cellGuess, cellID
      real (kind=RKIND), dimension(3) :: xPoint
      real (kind=RKIND), dimension(3, nCells) :: xCell
      real (kind=RKIND) :: dx,dy,dz      ! differences between two points
      real (kind=RKIND) :: r2Min         ! minimum squared distance found thus far
      real (kind=RKIND) :: r2            ! squared distance
      logical, pointer :: on_a_sphere    ! flag designating if we are on a sphere
      logical, pointer :: is_periodic    ! flag designating if periodicity is important
      real(kind=RKIND), pointer :: x_period, y_period

      integer, intent(inout) :: lastCell                          ! last known cell number

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
      call mpas_pool_get_config(meshPool, 'x_period', x_period)
      call mpas_pool_get_config(meshPool, 'y_period', y_period)

      ! could also be last than 1 for greater use in other cases (used in general since first index is 1)
      if (lastCell < 1) then
         ! brute force solution
         !{{{
         ! normalize locations to same spherical shell (unit) for direct comparison !{{{
         if (on_a_sphere) then
            pointRadius = sqrt(xp*xp + yp*yp + zp*zp)
            xPoint = (/ xp/pointRadius , yp/pointRadius , zp/pointRadius /)
            do aPoint = 1, nCells
               pointRadius = sqrt(xc(aPoint)*xc(aPoint) + yc(aPoint)*yc(aPoint) + zc(aPoint)*zc(aPoint))
               xCell(:,aPoint) = (/ xc(aPoint), yc(aPoint), zc(aPoint) /) / pointRadius
            end do
         else
            xPoint = (/ xp,yp,zp /)
            do aPoint = 1, nCells
               if (is_periodic) then
                  xCell(:, aPoint) = (/ mpas_fix_periodicity(xc(aPoint), xp, x_period), mpas_fix_periodicity(yc(aPoint), yp, y_period), zc(aPoint) /)
               else
                  xCell(:, aPoint) = (/ xc(aPoint), yc(aPoint), zc(aPoint) /)
               end if
            end do
         end if
         !}}}

         ! brute force solution
         cellID = 1
         dx = xPoint(1) - xCell(1,cellID)
         dy = xPoint(2) - xCell(2,cellID)
         dz = xPoint(3) - xCell(3,cellID)
         r2Min = dx*dx + dy*dy + dz*dz

         ! loop over all the other cells to get minimum
         do aPoint = 2, nCells
            ! compute squared distance
            dx = xPoint(1) - xCell(1,aPoint)
            dy = xPoint(2) - xCell(2,aPoint)
            dz = xPoint(3) - xCell(3,aPoint)
            r2 = dx*dx + dy*dy + dz*dz
            if ( r2 < r2Min) then
               ! we have a new closest point
               cellID = aPoint
               r2Min = r2
            end if
         end do
         !}}}
      else
         ! greedy, nearest neighbor algorithm
         !{{{

         ! search neighbors via greedy algorithm
         !write(stderrUnit,*) 'get_cell_id', lastCell
         cellGuess = lastCell
         cellID = -1
         ! loop over closest cells to make this a greedy method
         do while (cellID /= cellGuess)
            !write(stderrUnit,*) cellGuess, cellID

            ! we have a known cell
            cellID = cellGuess

            ! normalize locations to same spherical shell (unit) for direct comparison !{{{
            if (on_a_sphere) then
               pointRadius = sqrt(xp*xp + yp*yp + zp*zp)
               xPoint = (/ xp/pointRadius , yp/pointRadius , zp/pointRadius /)
               ! for point itself
               pointRadius = sqrt(xc(cellID)*xc(cellID) + yc(cellID)*yc(cellID) + zc(cellID)*zc(cellID))
               xCell(:,cellID) = (/ xc(cellID), yc(cellID), zc(cellID) /) / pointRadius
               ! for point neighbors
               do iPoint = 1, nEdgesOnCell(cellID)
                  aPoint = cellsOnCell(iPoint, cellID)
                  if (aPoint > nCells) cycle
                  pointRadius = sqrt(xc(aPoint)*xc(aPoint) + yc(aPoint)*yc(aPoint) + zc(aPoint)*zc(aPoint))
                  xCell(:,aPoint) = (/ xc(aPoint), yc(aPoint), zc(aPoint) /) / pointRadius
               end do
            else
               xPoint = (/ xp,yp,zp /)
               ! for point itself
               if (is_periodic) then
                  xCell(:, cellID) = (/ mpas_fix_periodicity(xc(cellID), xp, x_period), mpas_fix_periodicity(yc(cellID), yp, y_period), zc(cellID) /)
               else
                  xCell(:, cellID) = (/ xc(cellID), yc(cellID), zc(cellID) /)
               end if
               ! for point neighbors
               do iPoint = 1, nEdgesOnCell(cellID)
                  aPoint = cellsOnCell(iPoint, cellID)
                  if (aPoint > nCells) cycle
                  if (is_periodic) then
                     xCell(:, aPoint) = (/ mpas_fix_periodicity(xc(aPoint), xp, x_period), mpas_fix_periodicity(yc(aPoint), yp, y_period), zc(aPoint) /)
                  else
                     xCell(:, aPoint) = (/ xc(aPoint), yc(aPoint), zc(aPoint) /)
                  end if
               end do
            end if
            !}}}

            dx = xPoint(1) - xCell(1,cellID)
            dy = xPoint(2) - xCell(2,cellID)
            dz = xPoint(3) - xCell(3,cellID)
            r2Min = dx*dx + dy*dy + dz*dz
            !write(stderrUnit,*) 'r2Min_1 =', r2Min, ' cellID = ', cellID

            ! loop over all the other on the existing cell to see if we have found a closer cell
            do iPoint = 1, nEdgesOnCell(cellGuess)
               aPoint = cellsOnCell(iPoint, cellID)
               if (aPoint > nCells) cycle
               ! compute squared distance
               dx = xPoint(1) - xCell(1,aPoint)
               dy = xPoint(2) - xCell(2,aPoint)
               dz = xPoint(3) - xCell(3,aPoint)
               r2 = dx*dx + dy*dy + dz*dz
               if ( r2 < r2Min) then
                  ! we have a new closest point
                  cellGuess = aPoint
                  r2Min = r2
               end if
            end do
            !write(stderrUnit,*) 'r2Min_2 = ' ,r2Min, ' cellguess = ', cellguess
         end do
         !}}}
      end if

      ! return value
      lastCell = cellID

      ! could also perform a hit-test to make sure we got the right cell just to be sure
      ! this should be done outside this function for clarity

   end subroutine mpas_get_nearby_cell_index !}}}


!***********************************************************************
!
!  function mpas_get_vertical_id
!
!> \brief Determines vertical level index for a location
!> \author  Phillip Wolfram
!> \date    05/13/2014
!> \details
!>  This function determines the nearest vertical level to a particular
!>  location within a column.  If the location is below the bottom
!>  vertical level this returns -1, if it is above the top level, 0.
!>  Levels are referenced to vertical cell centers.
!-----------------------------------------------------------------------
   function mpas_get_vertical_id(nLevels, zLoc, zMid) !{{{

      implicit none
      ! in / out variables
      integer :: mpas_get_vertical_id ! index of vertical level
      real (kind=RKIND), intent(in) :: zLoc                 ! z-level location of point
      integer, intent(in) :: nLevels                           ! number of vertical levels
      real (kind=RKIND), dimension(:), intent(in) :: zMid      ! elevation of cell centers

      ! local variables
      integer :: aLevel                                        ! particular level returned

      ! cases except for bottom cell, if location is above the surface this yields 0
      do aLevel = 1, nLevels-1
         if((zMid(aLevel) <= zLoc .and. zLoc <= zMid(aLevel+1)) .or. (zMid(aLevel+1) <= zLoc .and. zLoc <= zMid(aLevel))) then
            ! the point is bounded by the levels (store the bottom level)
            mpas_get_vertical_id = aLevel
            return
         end if
      end do

      if(zLoc < minval(zMid(1:nLevels))) then
         ! case where location is smallest value
         mpas_get_vertical_id = -1
         return
      else if (zLoc > maxval(zMid(1:nLevels))) then
         ! case where location is largest value
         mpas_get_vertical_id = 0
         return
      end if

   end function mpas_get_vertical_id !}}}


!***********************************************************************
!
!  function mpas_wachspress_coordinates
!
!> \brief Compute the barycentric Wachspress coordinates for a polygon
!> \author  Phillip Wolfram
!> \date    01/26/2015
!> \details
!>  Computes the barycentric Wachspress coordinates for a polygon with nVertices
!>  points in R3, vertCoords for a particular pointInterp with normalized radius.
!>  Follows Gillette, A., Rand, A., Bajaj, C., 2011.
!>  Error estimates for generalized barycentric interpolation.
!>  Advances in computational mathematics 37 (3), 417–439.
!>  Optimized version of mpas_wachspress_coordinates uses optional cached B_i areas
!------------------------------------------------------------------------
   function mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, meshPool, areaBin) !{{{
      implicit none

      ! input points
      integer, intent(in) :: nVertices
      real (kind=RKIND), dimension(3, nVertices), intent(in) :: vertCoords
      real (kind=RKIND), dimension(3), intent(in) :: pointInterp
      real (kind=RKIND), dimension(nVertices), optional, intent(in) :: areaBin
      type (mpas_pool_type), pointer :: meshPool
      ! output
      real (kind=RKIND), dimension(nVertices) :: mpas_wachspress_coordinates
      ! computational intermediates
      real (kind=RKIND), dimension(nVertices) :: wach       ! The wachpress area-product
      real (kind=RKIND) :: wach_total                       ! The wachpress total weight
      integer :: i, j                                       ! Loop indices
      integer :: im1, i0, ip1                               ! im1 = (i-1), i0 = i, ip1 = (i+1)

      ! triangle areas to compute wachspress coordinate
      real (kind=RKIND), dimension(nVertices) :: areaA
      real (kind=RKIND), dimension(nVertices) :: areaB

      logical, pointer :: on_a_sphere
      real(kind=RKIND), pointer :: sphere_radius
      real(kind=RKIND) :: radiusLocal

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
      if ( on_a_sphere ) then
         radiusLocal = sphere_radius
      else
         radiusLocal = 1.0_RKIND
      end if

      if (.not. present(areaBin)) then
        ! compute areas
        do i = 1, nVertices
           ! compute first area B_i
           ! get vertex indices
           im1 = mod(nVertices + i - 2, nVertices) + 1
           i0  = mod(nVertices + i - 1, nVertices) + 1
           ip1 = mod(nVertices + i    , nVertices) + 1

           ! precompute B_i areas
           ! always the same because B_i independent of xp,yp,zp
           ! (COULD CACHE AND USE RESULT FROM ARRAY FOR FURTHER OPTIMIZATION)
           areaB(i) = mpas_triangle_signed_area(vertCoords(:, im1), vertCoords(:, i0), vertCoords(:, ip1), meshPool)
        end do
      else
        ! assign areas
        do i = 1, nVertices
          areaB(i) = areaBin(i)
        end do
      end if

      ! compute areas
      do i = 1, nVertices
         ! compute first area B_i
         ! get vertex indices
         im1 = mod(nVertices + i - 2, nVertices) + 1
         i0  = mod(nVertices + i - 1, nVertices) + 1
         ip1 = mod(nVertices + i    , nVertices) + 1

         ! compute A_ij areas
         ! must be computed each time
         areaA(i0) = mpas_triangle_signed_area(pointInterp, vertCoords(:, i0), vertCoords(:, ip1), meshPool)

         ! precomputed B_i areas, cached
      end do


      ! for each vertex compute wachpress coordinate
      do i = 1, nVertices
         wach(i) = areaB(i)
         do j = (i + 1), (i + nVertices - 2)
            i0  = mod(nVertices + j - 1, nVertices) + 1
            ! accumulate products for A_ij subareas
            wach(i) = wach(i) * areaA(i0)
         end do
      end do

      ! get summed weights for normalization
      wach_total = 0
      do i = 1, nVertices
         wach_total = wach_total + wach(i)
      end do

      ! compute lambda
      mpas_wachspress_coordinates= 0.0_RKIND
      do i = 1, nVertices
         mpas_wachspress_coordinates(i) = wach(i)/wach_total
      end do

   end function mpas_wachspress_coordinates!}}}


!***********************************************************************
!
!  routine mpas_wachspress_interpolate
!
!> \brief Interpolate using barycentric Wachspress coordinates
!> \author  Phillip Wolfram
!> \date    03/27/2015
!> \details
!>  Interpolate using the barycentric Wachspress coordinates for a polygon with nVertices
!>  having values phi.
!------------------------------------------------------------------------
   real (kind=RKIND) function mpas_wachspress_interpolate(lambda, phi) !{{{
      implicit none

      ! input points
      real (kind=RKIND), dimension(:), intent(in) :: lambda   !< Input: Wachspress coordinate / weight
      real (kind=RKIND), dimension(:), intent(in) :: phi      !< Input: values at lambda weights
      ! output for function
      !real (kind=RKIND), intent(out) :: mpas_wachspress_interpolate

      mpas_wachspress_interpolate = sum(phi * lambda)

   end function mpas_wachspress_interpolate! }}}


!***********************************************************************
!
!  subroutine mpas_unit_test_wachspress_hexagon
!
!> \brief   Simple unit test for testing wachspress interpolation on a hexagon.
!> \author  Phillip J. Wolfram
!> \date    07/21/2015
!> \details
!>  This routine tests the mpas_wachspress_coordinates and
!>  mpas_wachspress_interpolate routines on a hexagon.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_wachspress_hexagon(ierr)!{{{
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
      integer, parameter :: nVertices = 6
      real (kind=RKIND), dimension(3, nVertices) :: vertCoords
      real (kind=RKIND), dimension(nVertices) :: vertValues
      real (kind=RKIND), dimension(3) :: pointInterp
      type (mpas_pool_type), pointer :: meshPool
      integer :: j
      real (kind=RKIND), dimension(nVertices) :: lambda
      real (kind=RKIND), dimension(3) :: velocity
      real (kind=RKIND) :: eps = 1.0e-12_RKIND
      real (kind=RKIND) :: theta, interpValue

      ierr = 0
#if MPAS_DEBUG
      write(stderrUnit,*) 'starting mpas_unit_test_wachspress_hexagon'
#endif
      ! set on a plane
      call mpas_pool_create_pool(meshPool)
      call mpas_pool_add_config(meshPool,'on_a_sphere', .False.)
      call mpas_pool_add_config(meshPool, 'sphere_radius', 1.0_RKIND)
      ! hexagon geometry values
      do j = 1,nVertices
        theta = 2*pii/nVertices * (j-1)
        vertCoords(1,j) = cos(theta)
        vertCoords(2,j) = sin(theta)
        vertCoords(3,j) = 0.0_RKIND
        vertValues(j) = 3.0_RKIND + vertCoords(1,j) + 2.0_RKIND*vertCoords(2,j)
      end do

      ! interpolation at hexagon center
      pointInterp(:) = (/0.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      ! test without precached areas
      lambda = mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, meshPool)
      interpValue = mpas_wachspress_interpolate(lambda, vertValues)
      if (abs(interpValue - 3.0_RKIND) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_wachspress_hexagon: test1 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_wachspress_hexagon:  test1 - SUCCESS'
      end if

      ! interpolation at vertex
      pointInterp(:) = (/1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      ! test without precached areas
      lambda = mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, meshPool)
      interpValue = mpas_wachspress_interpolate(lambda, vertValues)
      if (abs(interpValue - 4.0_RKIND) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_wachspress_hexagon: test2 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_wachspress_hexagon:  test2 - SUCCESS'
      end if

      call mpas_pool_destroy_pool(meshPool)
#if MPAS_DEBUG
      write(stderrUnit,*) 'finished mpas_unit_test_wachspress_hexagon'
#endif

   end subroutine mpas_unit_test_wachspress_hexagon!}}}


!***********************************************************************
!
!  subroutine mpas_unit_test_wachspress_triangle
!
!> \brief   Simple unit test for testing wachspress interpolation on a triangle.
!> \author  Phillip J. Wolfram
!> \date    07/21/2015
!> \details
!>  This routine tests the mpas_wachspress_coordinates and
!>  mpas_wachspress_interpolate routines on a triangle.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_wachspress_triangle(ierr)!{{{
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
      integer, parameter :: nVertices = 3
      real (kind=RKIND), dimension(3, nVertices) :: vertCoords
      real (kind=RKIND), dimension(nVertices) :: vertValues
      real (kind=RKIND), dimension(3) :: pointInterp
      type (mpas_pool_type), pointer :: meshPool
      integer :: j
      real (kind=RKIND), dimension(nVertices) :: lambda
      real (kind=RKIND), dimension(3) :: velocity
      real (kind=RKIND) :: eps = 1.0e-12_RKIND
      real (kind=RKIND) :: interpValue

      ierr = 0
#if MPAS_DEBUG
      write(stderrUnit,*) 'starting mpas_unit_test_wachspress_triangle'
#endif
      ! set on a plane
      call mpas_pool_create_pool(meshPool)
      call mpas_pool_add_config(meshPool,'on_a_sphere', .False.)
      call mpas_pool_add_config(meshPool, 'sphere_radius', 1.0_RKIND)
      ! hexagon geometry values
      vertCoords(:,1) = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      vertCoords(:,2) = (/ 0.0_RKIND, 1.0_RKIND, 0.0_RKIND /)
      vertCoords(:,3) = (/ 0.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      vertValues(:) = 3.0_RKIND + vertCoords(1,:) + 2.0_RKIND*vertCoords(2,:)

      ! interpolation at triangle interior
      pointInterp(:) = (/0.25_RKIND, 0.25_RKIND, 0.0_RKIND /)
      ! test without precached areas
      lambda = mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, meshPool)
      interpValue = mpas_wachspress_interpolate(lambda, vertValues)
      if (abs(interpValue - 3.75_RKIND) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_wachspress_triangle: test1 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_wachspress_triangle: test1 - SUCCESS'
      end if

      ! interpolation at vertex
      pointInterp(:) = (/1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      ! test without precached areas
      lambda = mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, meshPool)
      interpValue = mpas_wachspress_interpolate(lambda, vertValues)
      if (abs(interpValue - 4.0_RKIND) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_wachspress_triangle: test2 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_wachspress_triangle: test2 - SUCCESS'
      end if

      call mpas_pool_destroy_pool(meshPool)
#if MPAS_DEBUG
      write(stderrUnit,*) 'finished mpas_unit_test_wachspress_triangle'
#endif

   end subroutine mpas_unit_test_wachspress_triangle!}}}


!***********************************************************************
!
!  routine mpas_convert_xyz_velocity_to_latlon
!
!> \brief Determine the zonal and meridional velocity field
!> \author  Phillip Wolfram
!> \date    09/13/2014
!> \details
!>  Returns the zonal and meridional velocity field (Uzon,Umer) from
!>  velocity field defined in R3: (x,y,z) point with velocity
!>  vectors (Ux,Uy,Uz).
!------------------------------------------------------------------------
   subroutine mpas_convert_xyz_velocity_to_latlon(Uzon, Umer, xyzPoint, xyzVel) !{{{
      implicit none

      ! input
      real(kind=RKIND), dimension(3), intent(in) :: xyzPoint, xyzVel
      ! output
      real(kind=RKIND), intent(out) :: Uzon, Umer
      ! local
      real(kind=RKIND) :: Rxy, Rxyz, slon,clon,slat,clat

      ! test for singularities at the poles
      if(xyzPoint(1) == 0.0_RKIND .and. xyzPoint(2) == 0.0_RKIND) then
#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'Point is located at a pole, cannot convert accurately to lat/lon velocity!'
#endif
        ! velocity is undefined but return 0s for testing purposes
        Uzon = 0.0_RKIND
        Umer = 0.0_RKIND
        return
      end if

      ! compute geometric coordinate transform coefficients
      Rxy  = sqrt(xyzPoint(1)**2.0_RKIND + xyzPoint(2)**2.0_RKIND)
      Rxyz = sqrt(xyzPoint(1)**2.0_RKIND + xyzPoint(2)**2.0_RKIND + xyzPoint(3)**2.0_RKIND)
      slon = xyzPoint(2)/Rxy
      clon = xyzPoint(1)/Rxy
      slat = xyzPoint(3)/Rxyz
      clat = Rxy/Rxyz


      ! compute the zonal and meridional velocity fields
      Uzon = -slon*xyzVel(1) + clon*xyzVel(2)
      Umer = -slat*(clon*xyzVel(1) + slon*xyzVel(2)) + clat*xyzVel(3)

   end subroutine mpas_convert_xyz_velocity_to_latlon !}}}


!***********************************************************************
!
!  subroutine mpas_unit_test_velocity_conversion
!
!> \brief   Simple unit test for testing the mpas_convert_xyz_velocity_to_latlon routine.
!> \author  Phillip J. Wolfram
!> \date    07/15/2015
!> \details
!>  This routine tests the mpas_convert_xyz_velocity_to_latlon routine.
!-----------------------------------------------------------------------
   subroutine mpas_unit_test_velocity_conversion(ierr)!{{{
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
      real (kind=RKIND), dimension(3) :: sphereLocation
      real (kind=RKIND), dimension(3) :: velocity
      real (kind=RKIND) :: zonalVel, meridVel
      real (kind=RKIND) :: eps = 1.0e-12_RKIND

      ierr = 0

      ! equator case (merid vel)
      sphereLocation = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      velocity = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      call mpas_convert_xyz_velocity_to_latlon(zonalVel, meridVel, sphereLocation, velocity)
      if (abs(zonalVel - 0.0_RKIND ) > eps .or. &
          abs(meridVel - 1.0_RKIND ) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_velocity_conversion: test1 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_velocity_conversion: test1 - SUCCESS'
      end if

      ! equator case (merid vel)
      sphereLocation = (/ 0.0_RKIND, 1.0_RKIND, 0.0_RKIND /)
      velocity = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      call mpas_convert_xyz_velocity_to_latlon(zonalVel, meridVel, sphereLocation, velocity)
      if (abs(zonalVel - 0.0_RKIND ) > eps .or. &
          abs(meridVel - 1.0_RKIND ) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_velocity_conversion: test2 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_velocity_conversion: test2 - SUCCESS'
      end if

      ! pole case
      sphereLocation = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      velocity = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      call mpas_convert_xyz_velocity_to_latlon(zonalVel, meridVel, sphereLocation, velocity)
      if (abs(zonalVel - 0.0_RKIND ) > eps .or. &
          abs(meridVel - 0.0_RKIND ) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_velocity_conversion: test3 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_velocity_conversion: test3 - SUCCESS'
      end if

      ! equator case (zonal vel)
      sphereLocation = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      velocity = (/ 0.0_RKIND, 1.0_RKIND, 0.0_RKIND /)
      call mpas_convert_xyz_velocity_to_latlon(zonalVel, meridVel, sphereLocation, velocity)
      if (abs(zonalVel - 1.0_RKIND ) > eps .or. &
          abs(meridVel - 0.0_RKIND ) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_velocity_conversion: test3 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_velocity_conversion: test3 - SUCCESS'
      end if

      ! equator case (zonal vel)
      sphereLocation = (/ 0.0_RKIND, 1.0_RKIND, 0.0_RKIND /)
      velocity = (/ -1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      call mpas_convert_xyz_velocity_to_latlon(zonalVel, meridVel, sphereLocation, velocity)
      if (abs(zonalVel - 1.0_RKIND ) > eps .or. &
          abs(meridVel - 0.0_RKIND ) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_velocity_conversion: test4 - FAILED.'
         ierr = 1
      else
         write(stderrUnit,*) 'mpas_unit_test_velocity_conversion: test4 - SUCCESS'
      end if

   end subroutine mpas_unit_test_velocity_conversion!}}}


!***********************************************************************
!
!  routine mpas_spherical_linear_interp
!
!> \brief Spherical linear interpolation of a point p0 to p1
!> \author  Phillip Wolfram
!> \date    04/21/2014
!> \details
!>  Returns the point following a spherical linear interpolation from p0 to p1
!>     with alpha = 0 => p0 and alpha = 1 => p1.  All interpolated points are
!>     on great circle route between p0 and p1, assuming that p0 and p1 are
!>     on the same spherical shell.
!------------------------------------------------------------------------
   subroutine mpas_spherical_linear_interp(pInterp, p0, p1, alpha) !{{{

     implicit none
     ! input variables
     real (kind=RKIND), dimension(3), intent(in) :: p0
     real (kind=RKIND), dimension(3), intent(in) :: p1
     real (kind=RKIND), intent(in) :: alpha

     ! function output variable
     real (kind=RKIND), dimension(3), intent(out) :: pInterp

     ! local variables
     real (kind=RKIND) :: p0mag, p1mag, omega, dotProd
     real (kind=RKIND), dimension(size(p0)) :: p0scaled
     real (kind=RKIND), dimension(size(p1)) :: p1scaled
     real (kind=RKIND), parameter :: eps = 1e-14_RKIND

     ! compute omega for angle of subtended arc
     p0mag = sqrt(sum(p0 * p0))
     p1mag = sqrt(sum(p1 * p1))
     ! scale for stability
     p0scaled = p0/p0mag
     p1scaled = p1/p1mag
     dotProd = min(1.0_RKIND,max(-1.0_RKIND,sum(p0scaled * p1scaled)))
     ! may need to potentially handle round-off errors with max?
     omega = acos(dotProd)
     ! handle case where p0 \approx p1
     if(abs(omega) < eps ) then
        pInterp = p0
        return
     end if
     pInterp = sin( (1.0_RKIND-alpha) * omega) / sin(omega) * p0 &
             + sin(alpha * omega) / sin(omega) * p1

   end subroutine mpas_spherical_linear_interp !}}}

end module mpas_geometry_utils
