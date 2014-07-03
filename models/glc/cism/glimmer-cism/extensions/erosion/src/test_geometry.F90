!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_geometry.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! testing geometry module

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program testgeom
  use geometry
  implicit none
  type(geom_poly) :: polygon, poly2
  type(coord_point) :: test
  logical res

  integer i,j, i1, j1

  polygon = poly_new(10)
  poly2 = poly_new(4)
  
  call poly_add_vert(polygon,10.,1.)
  call poly_add_vert(polygon,15.,21.)
  call poly_add_vert(polygon,8.,25.)
  call poly_add_vert(polygon,5.,19.)

  call poly_add_vert(poly2,8.,-1.)
  call poly_add_vert(poly2,12.,-1.)
  call poly_add_vert(poly2,12.,3.)
  call poly_add_vert(poly2,8.,3.)

!  call poly_add_vert(poly2,7.,2.)
!  call poly_add_vert(poly2,12.,2.)
!  call poly_add_vert(poly2,12.,24.)
!  call poly_add_vert(poly2,7.,24.)

  call poly_print(polygon,6)

  write(*,*) 0.5*poly_area2(polygon)

  write(*,*) left(polygon%poly(1),polygon%poly(2),polygon%poly(3))

  res = quad_is_convex(polygon)
  write(*,*) res

  test%pt(1) = 9.
  test%pt(2) = 10.
  write(*,*) point_in_cpoly(test,polygon), test%pt(1), test%pt(2)

  test%pt(1) = 10.
  test%pt(2) = 1.
  write(*,*) point_in_cpoly(test,polygon), test%pt(1), test%pt(2)
  
  test%pt(1) = 9.
  test%pt(2) = 4.
  write(*,*) point_in_cpoly(test,polygon), test%pt(1), test%pt(2)

  test%pt(1) = 12.5
  test%pt(2) = 11.
  write(*,*) point_in_cpoly(test,polygon), test%pt(1), test%pt(2)

  do i=1,polygon%n
     i1 = mod(i+polygon%n-2,polygon%n)+1
     do j=1,poly2%n
        j1 = mod(j+poly2%n-2,poly2%n)+1
        res = intersection(polygon%poly(i1),polygon%poly(i),poly2%poly(j1),poly2%poly(j),test)
        write(*,*) i,i1,j,j1,res , test%pt(1), test%pt(2)
     end do
  end do
  
end program testgeom
