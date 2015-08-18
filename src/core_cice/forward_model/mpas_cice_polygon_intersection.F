! Copyright (c) 2013, Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
! cice_polygon_intersection
!
!> \brief cice polygon intersection module
!> \author David Lee and William Lipscomb
!> \date April 2015
!> \details
!> This module provides a routine for determining the intersection region of
!> two polygons. To be used by the incremental remapping advection scheme in CICE.
!
!-----------------------------------------------------------------------

module cice_polygon_intersection

  use mpas_derived_types
  use mpas_io_units

  implicit none

  private
  save

  public :: cice_polygon_intersection_cpp, cice_polygon_intersection_fortran

  ! Interface to the C++ function for performing the intersection calculation
  interface
     subroutine polygon_intersection ( px1, py1, pn1, px2, py2, pn2, pxInt, pyInt, pnIntMax ) bind ( c )
       use iso_c_binding
       real*8 :: px1(*)
       real*8 :: py1(*)
       integer*4 :: pn1
       real*8 :: px2(*)
       real*8 :: py2(*)
       integer*4 :: pn2
       real*8 :: pxInt(*)
       real*8 :: pyInt(*)
       integer*4 :: pnIntMax
     end subroutine polygon_intersection
  end interface

   logical, parameter :: verbose_polyint = .false.
!!   logical, parameter :: verbose_polyint = .true.
  
contains

!***********************************************************************
!
! routine cice_polygon_intersection_cpp
!
!> \brief MPAS cice polygon intersection routine
!> \author David Lee and William Lipscomb
!> \date April 2015
!> \details
!> This routine calls C++ code to calculate the x and y coordinates of the vertices 
!> of an intersecting polygon of two input polygons.
!>
!> nPolyVerts1 - number of vertices of the first input polygon
!> polyVerts1 - coordinates of the vertices of the first input polygon
!> nPolyVerts2 - number of vertices of the second input polygon
!> polyVerts2 - coordinates of the vertices of the second input polygon
!>
!> Output:
!> nPolyVertsO - number of vertices of the intersecting polygon (0 if no such polygon exists)
!> polyVertsO - coordinates of the vertices of the intersecting polygon
!-----------------------------------------------------------------------
  
  subroutine cice_polygon_intersection_cpp( nPolyVerts1, polyVerts1, &
                                            nPolyVerts2, polyVerts2, &
                                            maxPolyVerts,            &
                                            nPolyVertsO, polyVertsO )!{{{

    integer, intent(in) :: nPolyVerts1   !< Input: number of vertices of the first polygon
    integer, intent(in) :: nPolyVerts2   !< Input: number of vertices of the second polygon
    integer, intent(in) :: maxPolyVerts  !< Input: maximum number of vertices of the intersecting polygon
    integer, intent(out) :: nPolyVertsO  !< Output: number of vertices of the intersecting polygon

    real(kind=RKIND), dimension(2,nPolyVerts1), intent(in) :: polyVerts1   !< Input: x/y coordinates of vertices of the first polygon
    real(kind=RKIND), dimension(2,nPolyVerts2), intent(in) :: polyVerts2   !< Input: x/y coordinates of vertices of the second polygon
    real(kind=RKIND), dimension(2,maxPolyVerts), intent(out) :: polyVertsO !< Output: x/y coordinates of vertices of the intersecting polygon

    ! Initialize the number of vertices in the intersection polygon to the max value
    ! polygon_intersection_c will return the actual number of vertices in the intersection polygon 
    nPolyVertsO = maxPolyVerts

    ! call the c++ function
    call polygon_intersection_c( polyVerts1(1,:), polyVerts1(2,:), nPolyVerts1, &
                                 polyVerts2(1,:), polyVerts2(2,:), nPolyVerts2, &
                                 polyVertsO(1,:), polyVertsO(2,:), nPolyVertsO )

    ! Remove any degenerate points from the output polygon

    call remove_degenerate_vertices(nPolyVertsO,  &
                                    polyVertsO)

  end subroutine cice_polygon_intersection_cpp!}}}

!***********************************************************************
!
! routine cice_polygon_intersection_fortran
!
!> \brief MPAS cice polygon intersection routine written in Fortran 90
!> \author William Lipscomb
!> \date May 2015
!> \details
!> This routine calculates the x and y coordinates of the vertices of an intersecting polygon
!> of two input polygons, based on algorithm of O'Rourke et al. (1982). Unlike the algorithm
!> called above (polygon_intersection_c), it is written in Fortran.
!>
!> nPolyVerts1 - number of vertices of the first input polygon
!> polyVerts1 - coordinates of the vertices of the first input polygon
!> nPolyVerts2 - number of vertices of the second input polygon
!> polyVerts2 - coordinates of the vertices of the second input polygon
!>
!> Output:
!> nPolyVertsO - number of vertices of the intersecting polygon (0 if no such polygon exists)
!> polyVertsO - coordinates of the vertices of the intersecting polygon
!>
!> Reference: 
!>   J. O'Rourke, C. Chien, T. Olson, D. Naddor, 1982: A new linear algorithm for
!>      intersecting convex polygons, Comp. Graph. Image Proc., 19, 384-391.
!>
!-----------------------------------------------------------------------
  
  subroutine cice_polygon_intersection_fortran(&
       nPolyVerts1In, polyVerts1In, &
       nPolyVerts2In, polyVerts2In, &
       maxPolyVerts,                &
       !!             nPolyVertsO, polyVertsO )
       nPolyVertsO, polyVertsO, iEdge, etest )

    integer, intent(in) :: nPolyVerts1In   !< Input: number of vertices of the first polygon
    integer, intent(in) :: nPolyVerts2In   !< Input: number of vertices of the second polygon
    integer, intent(in) :: maxPolyVerts  !< Input: maximum number of vertices of the intersecting polygon
    integer, intent(out) :: nPolyVertsO  !< Output: number of vertices of the intersecting polygon

    real(kind=RKIND), dimension(2,nPolyVerts1In), intent(in) :: polyVerts1In   !< Input: x/y coordinates of vertices of the first polygon
    real(kind=RKIND), dimension(2,nPolyVerts2In), intent(in) :: polyVerts2In   !< Input: x/y coordinates of vertices of the second polygon
    real(kind=RKIND), dimension(2,maxPolyVerts), intent(out) :: polyVertsO !< Output: x/y coordinates of vertices of the intersecting polygon

    !WHL - debug
    integer, intent(in) :: iEdge, etest  ! write diagnostic output if iEdge = etest

    ! Local variables
    !TODO - Choose a more descriptive name for pint

    integer :: nPolyVerts1, nPolyVerts2     ! number of vertices of the 2 polygons
                                            ! can be less than input value if the input polygon is degenerate
 
    real(kind=RKIND), dimension(2,nPolyVerts1In) :: polyVerts1   ! x/y coordinates of vertices of the first polygon
    real(kind=RKIND), dimension(2,nPolyVerts2In) :: polyVerts2   ! x/y coordinates of vertices of the second polygon

    real(kind=RKIND), dimension(2) ::  &
         edge1Start, edge1End,  &   ! x/y coordinates of the two vertices of an edge of polygon 1
         edge2Start, edge2End,  &   ! x/y coordinates of the two vertices of an edge of polygon 2
         edge1Vector,           &   ! 2D vector from start point to end point of edge 1
         edge2Vector,           &   ! 2D vector from start point to end point of edge 2
         pint,                  &   ! x/y coordinates of an intersection point
         pintFirst,             &   ! x/y coordinates of the first intersection point
         pintPrev                   ! x/y coordinates of the previous intersection point

    real(kind=RKIND) ::  &
         crossProduct        ! z component of the cross product of 2 vectors in the x-y plane

    logical :: &
         edgeIntersect,    & ! true if two polygon edges intersect
         startIntersect      ! true when looking for the first intersection point

    integer ::  &
         insidePolygon,    & ! = 1 if first polygon is inside, = 2 if second polygon is inside
         nLoopMax,         & ! max number of times through the loop
         vertexIndex1, vertexIndex2    ! vertex index for end point of current edge of polygon 1 and 2

    integer :: count

    !WHL - debug
    integer :: n, m, p

    ! initialize
    nPolyVerts1 = nPolyVerts1In
    polyVerts1(:,:) = polyVerts1In(:,:)

    nPolyVerts2 = nPolyVerts2In
    polyVerts2(:,:) = polyVerts2In(:,:)

    polyVertsO(:,:) = 0.0_RKIND
    nPolyVertsO = 0

    ! Remove any degenerate points from input polygon 1
    ! This is the departure region, which may have degenerate vertices where velocity = 0

    call remove_degenerate_vertices(nPolyVerts1,  &
                                    polyVerts1)

    ! This call commented out based on the assumption that each cell is nondegenerate
!    call remove_degenerate_vertices(nPolyVerts2,  &
!                                    polyVerts2)


    !WHL - Make sure both input polygons have at least 3 vertices

    if (nPolyVerts1 < 3 .or. nPolyVerts2 < 3) then  ! no work to do here
       return
    endif

    ! Initialize the polygon intersection algorithm
    pintFirst(:) = huge(0.0_RKIND)
    pintPrev(:) = huge(0.0_RKIND)
    insidePolygon = 0
    startIntersect = .true.
 
    nLoopMax = 2 * (nPolyVerts1 + nPolyVerts2)

    ! Choose a starting edge for each polygon
    ! Arbitrarily, put edgeStart at 1st vertex and edgeEnd at 2nd vertex for each polygon
    ! vertexIndex is the index of the end point of the edge

    vertexIndex1 = 2
    edge1Start(:)  = polyVerts1(:,vertexIndex1-1)
    edge1End(:)    = polyVerts1(:,vertexIndex1)
    edge1Vector(:) = edge1End(:) - edge1Start(:)

    vertexIndex2 = 2
    edge2Start(:)  = polyVerts2(:,vertexIndex2-1)
    edge2End(:)    = polyVerts2(:,vertexIndex2)
    edge2Vector(:) = edge2End(:) - edge2Start(:)

    !WHL - debug
    if (iEdge == etest) then
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'In cice_polygon_intersection, edge =', iEdge
       if (verbose_polyint) write(stderrUnit,*) 'nPolyVerts 1, 2 =', nPolyVerts1, nPolyVerts2
    endif

    ! loop around polygons (at most two times per polygon edge)

    do count = 1, nLoopMax

       ! Check whether the current edges intersect
       call find_line_intersection(&
            edge1Start,    edge1End,  &
            edge2Start,    edge2End,  &
            edgeIntersect,            &
            pint)

       !WHL - debug
       if (iEdge == etest .and. verbose_polyint) then
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'Edge 1 start:', Edge1Start(:)
          write(stderrUnit,*) 'Edge 1 end:', Edge1End(:)
          write(stderrUnit,*) 'Edge 2 start:', Edge2Start(:)
          write(stderrUnit,*) 'Edge 2 end:', Edge2End(:)
          write(stderrUnit,*) 'intersection:', pint(:)
          write(stderrUnit,*) 'edgeIntersect =', edgeIntersect
       endif
       
       if (edgeIntersect) then

          ! Check whether we should exit.  Do this if (1) this intersection point was not found on the previous iteration
          ! and (2) we have returned to the first intersection point.

          if (.not.identical_point(pint, pintPrev) .and. identical_point(pint, pintFirst)) then  

             !WHL - debug
             if (iEdge == etest) then
                if (verbose_polyint) then
                   write(stderrUnit,*) 'pint, pintFirst:', pint(:), pintFirst(:)
                   write(stderrUnit,*) 'Exiting, count =', count
                   write(stderrUnit,*) ' '
                endif
             endif
             
             exit    ! done computing vertices of the intersection polygon
             
          endif      ! back to first intersection point

          ! save this intersection point if it is the first one
          if (startIntersect) then
             pintFirst(:) = pint(:)
             startIntersect = .false.
          endif
          
          ! add this intersection point to the intersection polygon, if not already added
          if (.not. identical_point(pint, pintPrev)) then   ! not the same as the previous point

             nPolyVertsO = nPolyVertsO + 1
             polyVertsO(:,nPolyVertsO) = pint(:)

             ! Set pintPrev to this point
             ! This is to make sure we do not double-count this point (which could otherwise happen in the degenerate case
             !  of an intersection point lying on a polygon edge).
             pintPrev(:) = pint(:)

             !WHL - debug
             if (iEdge == etest .and. verbose_polyint) then
                write(stderrUnit,*) 'Add intersection point, nVertsO, pint =', nPolyVertsO, pint(:)
             endif

          endif  ! not identical to pintPrev
          
          ! check whether the current endpoint of polygon 1 is in the half plane of the current edge of polygon 2
          if (point_in_half_plane(edge1End, edge2Start, edge2End)) then
             insidePolygon = 1
          else
             insidePolygon = 2
          endif
          
          !WHL - debug
          if (iEdge == etest .and. verbose_polyint) then
             write(stderrUnit,*) 'insidePolygon =', insidePolygon
          endif

       endif        ! edgeIntersect

       ! compute cross product of the current edges
       ! (i.e., the z component of the cross product of two vectors in the x-y plane)

       call cross_product_2d(edge2Vector, edge1Vector, crossProduct)

       ! advance along polygon 1 or 2 based on the sign of the cross product
       !  and whether the endpoint of edge 1 is in the half plane of edge 2

       if (crossProduct >= 0.0_RKIND) then

          if (point_in_half_plane(edge1End, edge2Start, edge2End)) then

             ! advance along polygon 2
             call advance_polygon_edge(&
                  2,              &
                  insidePolygon,  &
                  vertexIndex2,   &
                  edge2Vector,    &
                  edge2Start,     &
                  edge2End,       &
                  polyVerts2,     &
                  nPolyVerts2,    &
                  polyVertsO,     &
                  maxPolyVerts,   &
                  !                 nPolyVertsO)
                  nPolyVertsO,    iEdge, etest)

             !WHL - debug
             if (iEdge == etest .and. verbose_polyint) then
                write(stderrUnit,*) 'cp >= 0, in half-plane, advance polygon 2'
             endif

          else  ! not in half plane

             !WHL - debug
             if (iEdge == etest .and. verbose_polyint) then
                write(stderrUnit,*) 'cp >= 0, not in half-plane, advance polygon 1'
             endif

             ! advance along polygon 1
             call advance_polygon_edge(&
                  1,              &
                  insidePolygon,  &
                  vertexIndex1,   &
                  edge1Vector,    &
                  edge1Start,     &
                  edge1End,       &
                  polyVerts1,     &
                  nPolyVerts1,    &
                  polyVertsO,     &
                  maxPolyVerts,   &
                  !                 nPolyVertsO)
                  nPolyVertsO,    iEdge, etest)

          endif   ! point in half plane

       else   ! crossProduct < 0

          if (point_in_half_plane(edge2End, edge1Start, edge1End)) then

             !WHL - debug
             if (iEdge == etest .and. verbose_polyint) then
                write(stderrUnit,*) 'cp < 0, in half-plane, advance polygon 1'
             endif

             ! advance along polygon 1
             call advance_polygon_edge(&
                  1,              &
                  insidePolygon,  &
                  vertexIndex1,   &
                  edge1Vector,    &
                  edge1Start,     &
                  edge1End,       &
                  polyVerts1,     &
                  nPolyVerts1,    &
                  polyVertsO,     &
                  maxPolyVerts,   &
                  !                 nPolyVertsO)
                  nPolyVertsO,    iEdge, etest)

          else   ! not in half plane

             !WHL - debug
             if (iEdge == etest .and. verbose_polyint) then
                write(stderrUnit,*) 'cp < 0, not in half-plane, advance polygon 2'
             endif

             ! advance along polygon 2
             call advance_polygon_edge(&
                  2,              &
                  insidePolygon,  &
                  vertexIndex2,   &
                  edge2Vector,    &
                  edge2Start,     &
                  edge2End,       &
                  polyVerts2,     &
                  nPolyVerts2,    &
                  polyVertsO,     &
                  maxPolyVerts,   &
                  !                 nPolyVertsO)
                  nPolyVertsO,    iEdge, etest)

          endif

       endif  ! crossProduct >= 0
       
    enddo     ! loop over polygon edges

    !WHL - debug
    if (iEdge == etest .and. verbose_polyint) then
       if (count > nLoopMax) then
          write(stderrUnit,*) 'No intersection'
       endif
    endif

    ! Remove any degenerate points from the output polygon

    call remove_degenerate_vertices(nPolyVertsO,  &
                                    polyVertsO)

  end subroutine cice_polygon_intersection_fortran

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine remove_degenerate_vertices
!
!> \brief  remove degenerate vertices from a polygon
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine checks to see if any two polygon vertices are identical,
!>  and if so removes the extra vertex.
!-----------------------------------------------------------------------

  subroutine remove_degenerate_vertices(nPolyVerts,  &
                                        polyVerts)

    integer, intent(inout)  :: nPolyVerts   !< Input/output: number of vertices of polygon
                                            !  May include degenerate vertices on input (but not output)

    real(kind=RKIND), dimension(:,:), intent(inout) :: polyVerts  !< Input/output: x/y coordinates of polygon vertices

    integer :: m, n, p

    if (nPolyVerts >= 2) then

       n = 1

       do while (n < nPolyVerts)
          m = n + 1

          do while (m <= nPolyVerts)

             if (identical_point(polyVerts(:,n), polyVerts(:,m))) then                
                ! remove this point from the polyVerts array and move the other points up
                do p = m, nPolyVerts-1
                   polyVerts(:,p) = polyVerts(:,p+1)
                enddo
                polyVerts(:,nPolyVerts) = 0.0_RKIND
                nPolyVerts = nPolyVerts - 1
             else
                m = m + 1
             endif
       
          enddo  ! m <= nPolyVerts
          
          n = n + 1

       enddo  ! n < nPolyVerts
       
    endif   ! nPolyVerts >= 2 

  end subroutine remove_degenerate_vertices

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine advance_polygon_edge
!
!> \brief  advance to the next edge in the polygon intersection algorithm
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine advances to the next edge of either polygon 1 or polygon 2
!>  in the polygon intersection algorithm
!-----------------------------------------------------------------------

  subroutine advance_polygon_edge(&
       polygonNumber,    &
       insidePolygon,    &
       vertexNumber,     &
       edgeVector,       &
       edgeStart,        &
       edgeEnd,          &
       polyVerts,        &
       nPolyVerts,       &
       polyVertsO,       &
       maxPolyVerts,     &
       !             nPolyVertsO)
       nPolyVertsO,    iEdge, etest)
       
    integer, intent(in) ::  &
         polygonNumber,    &  !< Input: input polygon (either 1 or 2) for edge advancement
         insidePolygon,    &  !< Input: polygon currently on the inside, either 1 or 2
         nPolyVerts,       &  !< Input: number of vertices of input polygon
         maxPolyVerts         !< Input: max number of vertices of intersection polygon

    integer, intent(inout) ::  &
         vertexNumber,     &  !< Input/output: vertex index of the end point of the edge
         nPolyVertsO          !< Input/output: cumulative number of vertices of the intersection polygon

    real(kind=RKIND), dimension(2), intent(inout) ::  &
         edgeStart, edgeEnd, &  !< Input/output:  x/y coordinates of start point and end point of edge
                                !                 old edge on input, new edge on output
         edgeVector             !< Input/output:  2D vector from start point to end point of edge
                                !                 old edge on input, new edge on output

    real(kind=RKIND), dimension(2,nPolyVerts), intent(in) ::  &
         polyVerts              !< Input: x/y coordinates of vertices input polygon

    real(kind=RKIND), dimension(2,maxPolyVerts), intent(inout) ::  &
         polyVertsO             !< Input/output: x/y coordinates of vertices of intersection polygon

    !WHL - debug
    integer, intent(in) :: iEdge, etest  ! write diagnostic output if iEdge = etest

    ! If the input polygon is currently the inside polygon, then add the end point of the current edge
    ! to the intersection polygon.
    ! Note: This procedure may add intersection points that have already been counted.
    !       These must be removed later.

    if (insidePolygon == polygonNumber) then

       nPolyVertsO = nPolyVertsO + 1
       PolyVertsO(:,nPolyVertsO) = edgeEnd(:)

       !WHL - debug
       if (iEdge == etest .and. verbose_polyint) then
          write(stderrUnit,*) 'Add intersection point, nVertsO, pint =', nPolyVertsO, edgeEnd(:)
       endif
       
    endif

    ! Advance from the current edge to the next edge of the input polygon
    vertexNumber = vertexNumber + 1
    if (vertexNumber > nPolyVerts) vertexNumber = 1
    edgeStart(:) = edgeEnd(:)
    edgeEnd(:) = polyVerts(:,vertexNumber)
    edgeVector(:) = edgeEnd(:) - edgeStart(:)

  end subroutine advance_polygon_edge

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine find_line_intersection
!
!> \brief  find the intersection point of two line segments in a plane
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given the Cartesian coordinates of the endpoints of 2 line segments in a plane, 
!>  find the point where the two line segments intersect (if they intersect).
!>  If the lines extending these segments intersect (even if the segment themselves
!>  do not), this routine returns the vertices of their intersection point.
!-----------------------------------------------------------------------

  !TODO - Change pint to a better name

  subroutine find_line_intersection(&
       p1, p2,        &
       p3, p4,        &
       lineIntersect, &
       pint,          &
       scalar)

    real(kind=RKIND), dimension(2), intent(in) :: &
         p1, p2,     & !< Input: x and y coordinates of the endpoints of the first line segment
         p3, p4        !< Input: x and y coordinates of the endpoints of the second line segment
        
    logical, intent(out) :: &
         lineIntersect !< Output: true if the lines intersect, else false
                       !    Note: If the lines share one or both endpoints, they do not intersect

    real(kind=RKIND), dimension(2), intent(out) :: &
         pint          !< Output:  x and y coordinates of the intersection point
                       !           returns (huge,huge) if the lines are parallel

    real(kind=RKIND), intent(out), optional :: &
         scalar        !< Output: a number describing where the intersection point lies on the second line
                       !          0 <= scalar <= 1 if the IP lies on the segment itself
                       !          scalar < 0 if the IP lies nearer p3 than p4
                       !          scalar > 1 if the IP lies nearer p4 than p3

    real(kind=RKIND) :: x1, x2, x3, x4, y1, y2, y3, y4
    real(kind=RKIND) :: rx, ry, sx, sy, t1, t2, rsCross, rsCrossMin

    ! initialize
    lineIntersect = .false.

    x1 = p1(1); y1 = p1(2)
    x2 = p2(1); y2 = p2(2)
    x3 = p3(1); y3 = p3(2)
    x4 = p4(1); y4 = p4(2)
    
    ! vector r = (x2,y2) - (x1,y1)
    rx = x2 - x1
    ry = y2 - y1

    ! vector s = (x4,y4) - (x3,y3)
    sx = x4 - x3
    sy = y4 - y3

    ! z component of cross product, r x s
    rsCross = rx*sy - ry*sx

    !TODO - Adjust this tolerance for single precision
    ! Min value for cross product (if smaller than this, the lines are deemed to be parallel)
    rsCrossMin = 1.e-10_RKIND * abs(rx*sy + ry*sx)

    ! Compute t1 and t2 such that (x1,y1) + t1*(rx,ry) = (x3,y3) + t2*(sx,sy)
    ! If both t1 and t2 are in the range [0,1], then the line segments intersect

    if (abs(rsCross) > rsCrossMin) then

       t1 = (sy*(x3 - x1) - sx*(y3-y1)) / rsCross
       t2 = (ry*(x3 - x1) - rx*(y3-y1)) / rsCross
 
       pint(1) = x1 + t1*rx
       pint(2) = y1 + t1*ry

       !WHL - Main IR module uses > and < rather than >= and <=, to avoid creating zero-area triangles.
       !      This algorithm, however, requires that if a vertex of one segment lies
       !       exactly on the other segment, the segments are deemed to intersect.
!!       if (t1 > 0.0_RKIND .and. t1 < 1.0_RKIND .and. t2 > 0.0_RKIND .and. t2 < 1.0_RKIND) then
       if (t1 >= 0.0_RKIND .and. t1 <= 1.0_RKIND .and. t2 >= 0.0_RKIND .and. t2 <= 1.0_RKIND) then
          lineIntersect = .true.
       endif

       if (present(scalar)) scalar = t2

    else   ! r x s = 0; the two lines are parallel

       pint(1) = huge(0.0_RKIND)
       pint(2) = huge(0.0_RKIND)

       if (present(scalar)) scalar = huge(0.0_RKIND)

    endif

  end subroutine find_line_intersection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cross_product_2d
!
!> \brief  compute the scalar value of the z component of the cross product
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given the components of two 2D vectors in the x-y plane, return the 
!>  scalar equal to the z component of cross product
!
!-----------------------------------------------------------------------

  subroutine cross_product_2d(&
       vector1, vector2,      &
       scalarOut)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         vector1, vector2        !< Input: x/y coordinates of two 2D Cartesian vectors

    real(kind=RKIND), intent(out) ::  &
         scalarOut               !< Output: z component of the cross product

    scalarOut = vector1(1)*vector2(2) - vector1(2)*vector2(1)

  end subroutine cross_product_2d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine point_in_half_plane
!
!> \brief  determine whether a given point is in the half-plane of a given line
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given the 2D Cartesian coordinates of two points on a line, along with a
!>  third point, determine whether the third point lies in the left half-plane
!>  of the line.
!
!-----------------------------------------------------------------------

  logical function point_in_half_plane(point, lineStart, lineEnd)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         point          ! x/y coordinates of an input point

    real(kind=RKIND), dimension(2), intent(in) ::  &
         lineStart, &  ! x/y coordinates of starting point on a line segment
         lineEnd     ! x/y coordinates of ending point on a line segment
    
    ! local variables

    real(kind=RKIND), dimension(2) :: &
         vector1, vector2   ! 2D work vectors in the x/y plane

    real(kind=RKIND) ::  &
         crossProduct       ! z component of the cross product

    ! vector from start point to end point of the input line
    vector1(:) = lineEnd(:) - lineStart(:)

    ! vector from start point of line to the input point
    vector2(:) = point(:) - lineStart(:)

    ! cross product
    call cross_product_2d(vector1, vector2, crossProduct)

    ! determine whether the point is in the half plane of the line
    if (crossProduct >= 0.0_RKIND) then
       point_in_half_plane = .true.
    else
       point_in_half_plane = .false.
    endif

  end function point_in_half_plane

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine identical_point
!
!> \brief  determine whether two points on a plane are identical
!> \author William Lipscomb
!> \date   May 2015
!> \details 
!>  This routine determines whether two points on a plane have identical
!>  x and y coordinates (within some tolerance).
!-----------------------------------------------------------------------

  logical function identical_point(p1, p2, toleranceIn)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         p1, p2      ! x/y coordinates of input points

    !TODO - Pass in tolerance based on grid spacing?
    real(kind=RKIND), optional :: toleranceIn   ! max difference between two numbers considered equal

    !TODO - Set appropriate tolerance for single precision
    real(kind=RKIND) :: tolerance

    if (present(toleranceIn)) then
       tolerance = toleranceIn
    else
       tolerance = 1.0e-08_RKIND   ! this value seems to work for standard test cases
    endif

    if ( abs(p1(1)-p2(1)) < tolerance .and. abs(p1(2)-p2(2)) < tolerance) then
       identical_point = .true.
    else
       identical_point = .false.
    endif
    
  end function identical_point


end module cice_polygon_intersection
