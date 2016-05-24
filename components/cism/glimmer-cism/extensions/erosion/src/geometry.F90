!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   geometry.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! module for handling polygon stuff
! - point in/outside polygon
! - intersection of two concave polygons
! - area of polygon

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module geometry
  use glimmer_global
  use glimmer_coordinates
  ! number of dimensions
  integer :: geom_unit = 6
  integer, parameter :: geom_dim = 2
  real(kind=dp), parameter, private :: tolerance = 1d-7

  private :: advance,add_intersect

  ! polygon type
  type geom_poly
     integer n
     type(coord_point), dimension(:), pointer :: poly
  end type geom_poly

  interface poly_add_vert
     module procedure poly_add_vert_d, poly_add_vert_s
  end interface
contains
  function poly_new(n)
    implicit none
    integer, intent(in) :: n
    type(geom_poly) :: poly_new

    poly_new%n = 0
    allocate(poly_new%poly(n))
  end function poly_new

  subroutine poly_delete(poly)
    implicit none
    type(geom_poly) :: poly
    deallocate(poly%poly)
  end subroutine poly_delete

  subroutine poly_add_vert_d(poly,x,y)
    implicit none
    type(geom_poly) :: poly
    real(kind=dp), intent(in) :: x,y

#ifdef DEBUG
    if (poly%n .eq. size(poly%poly)) then
       write(geom_unit,*) 'Error (',__FILE__,__LINE__,'): Cannot add extra vertices'
       write(geom_unit,*) '  current polygon:'
       call poly_print(poly,geom_unit)
    end if
#endif
    poly%n = poly%n + 1
    poly%poly(poly%n)%pt(1) = x
    poly%poly(poly%n)%pt(2) = y
  end subroutine poly_add_vert_d

  subroutine poly_add_vert_s(poly,x,y)
    implicit none
    type(geom_poly) :: poly
    real(kind=sp), intent(in) :: x,y

    call poly_add_vert_d(poly,real(x,kind=dp),real(y,kind=dp))
  end subroutine poly_add_vert_s

  function subvec(a,b)
    ! vector substraction
    implicit none
    type(coord_point), intent(in) :: a,b
    type(coord_point) subvec

    subvec%pt(:) = a%pt(:) - b%pt(:)
  end function subvec

  subroutine poly_print(poly, unit)
    implicit none
    type(geom_poly), intent(in) :: poly
    integer, intent(in) :: unit
    integer i
    
    do i=1, poly%n
       write(unit,*) i, poly%poly(i)%pt(1), poly%poly(i)%pt(2)
    end do
  end subroutine poly_print

  function triangle_area2(a, b, c)
    ! returns 2 the area of the triangle a, b, c
    ! positive if a, b, c are counterclockwise
    ! negative if a, b, c are clockwise
    implicit none
    type(coord_point), intent(in) :: a,b,c
    real(kind=dp) triangle_area2
    
    triangle_area2 = a%pt(1)*b%pt(2) - a%pt(2)*b%pt(1) + &
         a%pt(2)*c%pt(1) - a%pt(1)*c%pt(2) + &
         b%pt(1)*c%pt(2) - c%pt(1)*b%pt(2)
  end function triangle_area2

  function poly_area2(poly)
    ! returns 2 the area of the polygon poly
    implicit none
    type(geom_poly), intent(in) :: poly
    real(kind=dp) poly_area2
    integer i

    poly_area2 = 0.
    do i =1, poly%n - 1
       poly_area2 = poly_area2 + triangle_area2(poly%poly(1), poly%poly(i), poly%poly(i+1))
    end do
  end function poly_area2

  function left(a,b,c)
    ! determines if a,b,c are ordered ccw
    implicit none
    type(coord_point), intent(in) :: a,b,c
    logical left

    left = (triangle_area2(a,b,c) .gt. 0.d0)
  end function left
  
  function lefton(a,b,c)
    ! determines if a,b,c are ordered ccw
    implicit none
    type(coord_point), intent(in) :: a,b,c
    logical lefton

    lefton = (triangle_area2(a,b,c) .ge. 0.d0)
  end function lefton
   
  function collinear(a,b,c)
    ! are a,b,c on same line
    implicit none
    type(coord_point), intent(in) :: a,b,c
    logical collinear

    collinear = (abs(triangle_area2(a,b,c)) .le. tolerance)
  end function collinear
  
  function intersectProp(a,b,c,d)
    ! does [a,b] intersect [c,d]
    implicit none
    type(coord_point), intent(in) :: a,b,c,d
    logical intersectProp

    ! eliminate special cases
    if (collinear(a,b,c) .or. collinear(a,b,d) .or. collinear(c,d,a) .or. collinear(c,d,b)) then
       intersectProp = .false.
       return
    end if

    intersectProp = (left(a,b,c) .neqv. left(a,b,d)) .and. (left(c,d,a) .neqv. left(c,d,b))
  end function intersectProp

  function between(a,b,c)
    ! returns true iff (a,b,c) are colin and c lies on closed segment [a,b]
    implicit none
    type(coord_point), intent(in) :: a,b,c
    logical between
    
    if (.not.collinear(a,b,c)) then
       between = .false.
       return
    end if
    
    if (abs(a%pt(1)-b%pt(1)).gt.tolerance) then
       between = ((a%pt(1) .le. c%pt(1)) .and. (c%pt(1) .le. b%pt(1))) &
            .or. ((a%pt(1) .ge. c%pt(1)) .and. (c%pt(1) .ge. b%pt(1)))
    else
       between = ((a%pt(2) .le. c%pt(2)) .and. (c%pt(2) .le. b%pt(2))) &
            .or. ((a%pt(2) .ge. c%pt(2)) .and. (c%pt(2) .ge. b%pt(2)))
    end if
  end function between

  function intersect(a,b,c,d)
    ! returns true iff segments [a,b] and [c,d] intersect properly or improperly
    implicit none
    type(coord_point), intent(in) :: a,b,c,d
    logical intersect
    if (intersectProp(a,b,c,d)) then
       intersect = .true.
       return
    else if (between(a,b,c) .or. between(a,b,d) &
         .or. between(c,d,a) .or. between(c,d,b)) then
       intersect = .true.
       return
    else
       intersect = .false.
       return
    end if
  end function intersect

  function diagonalie(i,j,poly)
    ! returns true iff (vi,vj) is a proper internal or external diagonal of p
    implicit none
    integer, intent(in) :: i,j
    type(geom_poly), intent(in) :: poly
    logical diagonalie
    ! local vars
    integer k,k1

    ! loop over edges of poly
    do k=1,poly%n
       k1 = mod(k,poly%n)+1
       ! skip edges incident to i or j
       if (.not. ( (k.eq.i) .or. (k1.eq.i) .or. (k.eq.j) .or. (k1.eq.j))) then
          if (intersect(poly%poly(i),poly%poly(j),poly%poly(k),poly%poly(k1))) then
             diagonalie = .false.
          end if
       end if
    end do
    diagonalie = .true.
  end function diagonalie

  function incone(i,j,poly)
    ! returns true iff the diagonal (vi,vj) is strictly internal to the poly in the neighbourhood of the endpoint i
    implicit none
    integer, intent(in) :: i,j
    type(geom_poly), intent(in) :: poly
    logical incone
    ! local vars
    integer i1, in1

    i1 = mod(i,poly%n)+1
    in1 = mod(i+poly%n-2,poly%n) + 1
    
    ! if p(i) is convex vertex
    if (lefton(poly%poly(in1), poly%poly(i), poly%poly(i1))) then
       incone = left(poly%poly(i), poly%poly(j), poly%poly(in1)) .and. &
            left(poly%poly(j), poly%poly(i), poly%poly(i1))
       return
    else
       ! assume not colin, i.e. poly(i) is reflex
       incone = .not. (left(poly%poly(i), poly%poly(j), poly%poly(i1)) .and. &
            left(poly%poly(j), poly%poly(i), poly%poly(in1)))
       return
    end if
  end function incone

  function diagonal(i,j,poly)
    ! returns true iff the diagonal (vi,vj) is proper internal diag
    implicit none
    integer, intent(in) :: i,j
    type(geom_poly), intent(in) :: poly
    logical diagonal
    
    diagonal = incone(i,j,poly) .and. diagonalie(i,j,poly)
  end function diagonal

  function intersection(a,b,c,d,p)
    ! calculates intersection of segment [a,b] with [c,d], stores it in p and returns true if they intersect
    implicit none
    type(coord_point), intent(in) :: a,b,c,d
    logical intersection
    type(coord_point), intent(out) :: p

    real(kind=dp) s,t,denom

    denom = &
         a%pt(1) * (d%pt(2)-c%pt(2)) + &
         b%pt(1) * (c%pt(2)-d%pt(2)) + &
         d%pt(1) * (b%pt(2)-a%pt(2)) + &
         c%pt(1) * (a%pt(2)-b%pt(2))

    ! if denom is 0 then return false since segments are parallel (although they might overlap
    if (abs(denom) .le. tolerance) then
       intersection = .false.
       return
    end if

    s = ( &
         a%pt(1) * (d%pt(2)-c%pt(2)) + &
         c%pt(1) * (a%pt(2)-d%pt(2)) + &
         d%pt(1) * (c%pt(2)-a%pt(2)) &
         ) / denom

    t = - ( &
         a%pt(1) * (c%pt(2)-b%pt(2)) + &
         b%pt(1) * (a%pt(2)-c%pt(2)) + &
         c%pt(1) * (b%pt(2)-a%pt(2)) &
         ) / denom
         
    p%pt(:) = a%pt(:) + s*(b%pt(:) - a%pt(:))

    if (s.ge.0 .and. s.le.1. .and. t.ge.0. .and. t.le.1.) then
       intersection = .true.
       return
    end if

    intersection = .false.
  end function intersection
    
  function is_anti_clock(poly)
    ! check if verticies in poly are ordered anti clockwise
    type(geom_poly), intent(in) :: poly
    logical is_anti_clock

    integer i

    do i=1,poly%n-2
       is_anti_clock = left(poly%poly(i),poly%poly(i+1),poly%poly(i+2))
       if (.not.is_anti_clock) exit
    end do
  end function is_anti_clock

  function quad_is_convex(quad)
    ! is the quad concave?
    implicit none
    type(geom_poly), intent(in) :: quad
    logical quad_is_convex

#ifdef DEBUG
    if (quad%n .ne. 4) then
       write(*,*) 'Error (',__FILE__,__LINE__,'): quad has not four vertices, ',quad%n
    end if
#endif
    
    quad_is_convex = intersectProp(quad%poly(1),quad%poly(3),quad%poly(2),quad%poly(4))
  end function quad_is_convex

  function point_in_cpoly(pnt, poly)
    ! test if pnt is in concave polygon
    implicit none
    type(coord_point), intent(in) :: pnt
    type(geom_poly), intent(in) :: poly
    logical point_in_cpoly
    
    integer i

    point_in_cpoly = lefton(poly%poly(poly%n), poly%poly(1), pnt)
    do i=1,poly%n-1
       point_in_cpoly = point_in_cpoly .and. lefton(poly%poly(i), poly%poly(i+1), pnt)
    end do
  end function point_in_cpoly

  function ConvexIntersect(polya, polyb)
    ! intersection of two convex polygons
    implicit none
    type(geom_poly), intent(in) :: polya, polyb
    type(geom_poly) ConvexIntersect

    type(coord_point) :: p,q ! directed edges on polya and polyb
    type(coord_point) :: origin
    type(coord_point) :: pq ! intersection point
    integer a,a1,b,b1  ! indices on polya and polyb, and a-1, b-1
    integer i          ! loop counter
    real(kind=dp) bha,ahb    ! b in H(a); a in H(b)
    integer adva, advb ! number of advances since first inetersection
    real(kind=dp) cross
    integer inflag, count
    logical round 

    ! allocate a new polygon
    ConvexIntersect = poly_new(2*max(polya%n,polyb%n))
    round = .false.
    origin%pt = 0.
    a = 1
    b = 1
    i = 1
    inflag = -1
    adva = 0
    advb = 0
    count = 0
    do 
       count = count + 1
       a1 = mod(a+polya%n-2,polya%n)+1
       b1 = mod(b+polyb%n-2,polyb%n)+1
       p = subvec(polya%poly(a), polya%poly(a1))
       q = subvec(polyb%poly(b), polyb%poly(b1))
       cross = triangle_area2(origin, p,q)
       bha = triangle_area2(polya%poly(a1),polya%poly(a),polyb%poly(b))
       ahb = triangle_area2(polyb%poly(b1),polyb%poly(b),polya%poly(a))

#ifdef DEBUG_GEOM
       write(*,*) 'New loop:   ',count,a,b,cross,ahb,bha
#endif

       if (intersection(polya%poly(a1), polya%poly(a), polyb%poly(b1), polyb%poly(b), pq)) then
          if (inflag .eq. -1) then
             adva = 0
             advb = 0
          end if

          round = add_intersect(ConvexIntersect,pq)
#ifdef DEBUG_GEOM
          write(*,*) '   ', count,'intersection'
#endif
          if (ahb .gt. 0.) then
             inflag = 1
          end if
          if (bha .gt. 0) then
             inflag = 2
          end if
       end if

       ! advance rules
       ! special case a and b are collinear
       if (abs(cross).le.tolerance .and. abs(ahb).le.tolerance .and. abs(bha).le.tolerance) then
#ifdef DEBUG_GEOM  
          write(*,*) '   ', count,'spec'
#endif
          if (inflag .eq. 1) then
             call advance(b,advb,polyb%n,.false.,polyb%poly(b),ConvexIntersect,round)
          else
             call advance(a,adva,polya%n,.false.,polya%poly(a),ConvexIntersect,round)
          end if

          ! general rules
       else if (cross .ge. 0.) then
          if (bha .gt. 0) then
#ifdef DEBUG_GEOM  
             write(*,*) '   adv: ', count,'+ a', inflag
#endif
             call advance(a,adva,polya%n,(inflag.eq.1),polya%poly(a),ConvexIntersect,round)

          else
#ifdef DEBUG_GEOM  
             write(*,*) '   adv: ',count,'+ b', inflag
#endif
             call advance(b,advb,polyb%n,(inflag.eq.2),polyb%poly(b),ConvexIntersect,round)
          end if

       else
          if (ahb .gt. 0) then
#ifdef DEBUG_GEOM  
             write(*,*) '   adv: ',count,'- b', inflag
#endif
             call advance(b,advb,polyb%n,(inflag.eq.2),polyb%poly(b),ConvexIntersect,round)
          else
#ifdef DEBUG_GEOM  
             write(*,*) '   adv: ',count,'- a', inflag
#endif
             call advance(a,adva,polya%n,(inflag.eq.1),polya%poly(a),ConvexIntersect,round)
          end if
       end if

#ifdef DEBUG_GEOM  
       write(*,*) 'State :  ', count,adva,advb
#endif

       if ((adva.ge.polya%n .and. advb.ge.polyb%n) .or. adva.eq.2*polya%n .or. advb.eq.2*polyb%n .or. round) then
          exit
       end if
    end do
  end function ConvexIntersect

  subroutine advance(a, adv, n, flag, pta, poly, round)
    ! utility routine for ConvexIntersect
    implicit none
    integer, intent(inout) :: a,adv
    integer, intent(in) :: n
    logical, intent(in) :: flag
    type(coord_point), intent(in) :: pta
    type(geom_poly) :: poly
    logical round

    a = mod(a,n)+1
    adv = adv + 1
    if (flag) then
       round = add_intersect(poly,pta)
    end if
  end subroutine advance

  function add_intersect(poly,point)
    ! avoiding repeated vertices, return true iff at the starting point
    implicit none
    type(geom_poly) :: poly
    type(coord_point), intent(in) :: point
    logical add_intersect

    add_intersect = .false.
    if (poly%n .gt. 1) then
       if (any(abs(poly%poly(1)%pt-point%pt).gt.tolerance)) then
          call poly_add_vert(poly,point%pt(1),point%pt(2))
          return
       else
          add_intersect = .true.
          return
       end if
    end if
    if (poly%n .gt. 0) then ! more than one vertex
       ! check previous point
       if (any(abs(poly%poly(poly%n)%pt-point%pt).gt.tolerance)) then
          call poly_add_vert(poly,point%pt(1),point%pt(2))
          return
       end if
    end if
    call poly_add_vert(poly,point%pt(1),point%pt(2))
  end function add_intersect
end module geometry
