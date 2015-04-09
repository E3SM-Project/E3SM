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
   use mpas_configure
   use mpas_constants
   use mpas_io_units
   use mpas_vector_operations

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
!        write(stdoutUnit,*) ' i, indx ',i,indx(i)
!     end do
!
!     write(stdoutUnit,*) ' '

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

      ierr = 0

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
             baryWeightsOnPoints(1, iPoint), baryWeightsOnPoints(2, iPoint), baryWeightsOnPoints(3, iPoint), ierr)

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

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if (on_a_sphere) then
         mpas_triangle_signed_area = mpas_triangle_signed_area_sphere(a, b, c)
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
!>  Copied from mesh generator.
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area_sphere(a, b, c)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND) :: ab, bc, ca, semiperim, tanqe

      ab = mpas_arc_length(a(1), a(2), a(3), b(1), b(2), b(3))
      bc = mpas_arc_length(b(1), b(2), b(3), c(1), c(2), c(3))
      ca = mpas_arc_length(c(1), c(2), c(3), a(1), a(2), a(3))
      semiperim = 0.5 * (ab + bc + ca)

      tanqe = sqrt(tan(0.5_RKIND * semiperim) * tan(0.5_RKIND * (semiperim - ab)) &
                   * tan(0.5_RKIND * (semiperim - bc)) * tan(0.5_RKIND * (semiperim - ca)) )

      mpas_triangle_signed_area_sphere = 4.0_RKIND * atan(tanqe)
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
         write(stdoutUnit,*) 'mpas_unit_test_in_triangle: point1 test - SUCCESS'
      endif

      if (mpas_point_in_polygon(point2, trianglePoints, .false.)) then
         write(stderrUnit,*) 'Error in mpas_unit_test_in_triangle: point2 calculated to be inside of triangle.'
         ierr = 1
      else
         write(stdoutUnit,*) 'mpas_unit_test_in_triangle: point2 test - SUCCESS'
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
         write(stdoutUnit,*) 'mpas_unit_test_bary_weights: point1 test - SUCCESS'
      endif

      call mpas_calculate_barycentric_weights(point2, trianglePoints(1,:), trianglePoints(2,:), trianglePoints(3,:), meshPool,     &
             weights(1), weights(2), weights(3), ierr)
      if (maxval(abs(weights - (/ 1.0, 0.0, 0.0 /) )) > eps) then
         write(stderrUnit,*) 'Error in mpas_unit_test_bary_weights: point1 weights have error greater than tolerance.'
         ierr = 1
      else
         write(stdoutUnit,*) 'mpas_unit_test_bary_weights: point2 test - SUCCESS'
      endif

   end subroutine mpas_unit_test_bary_weights !}}}

!-------------------------------------------------------------

end module mpas_geometry_utils
