!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_map_proj4.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> Generates proj4 strings from projection data type.
!! Not used in GLIMMER at present.
module glimmer_map_proj4

  use glimmer_map_types

  implicit none

  private
  public :: glimmap_proj4

  integer, parameter :: proj4len=100

contains

  !> Returns a proj4 parameter string for a given set of projection parameters
  !! \return Pointer to array of projection parameter strings
  function glimmap_proj4(proj)
    
    use glimmer_log

    implicit none
    character(len=proj4len), dimension(:), pointer :: glimmap_proj4
    type(glimmap_proj) :: proj !*FD Projection of interest

    if (.not.glimmap_allocated(proj)) then
       call write_log('No known projection found!',GM_WARNING)
       return
    end if

    if (associated(proj%laea)) then
       glimmap_proj4 => glimmap_proj4_laea(proj%laea)
       return
    else if (associated(proj%aea)) then
       glimmap_proj4 => glimmap_proj4_aea(proj%aea)
       return
    else if (associated(proj%lcc)) then
       glimmap_proj4 => glimmap_proj4_lcc(proj%lcc)
       return
    else if (associated(proj%stere)) then
       glimmap_proj4 => glimmap_proj4_stere(proj%stere)
       return
    else
       call write_log('No known projection found!',GM_WARNING)
    end if
  end function glimmap_proj4

  !------------------------------------------------------------------
  ! private converters to proj4 strings
  !------------------------------------------------------------------

  !> Returns a proj4 parameter string for a stereographic projection
  function glimmap_proj4_stere(stere)
    implicit none
    character(len=proj4len), dimension(:), pointer :: glimmap_proj4_stere
    type(proj_stere) :: stere

    allocate(glimmap_proj4_stere(6))
    write(glimmap_proj4_stere(1),*) 'proj=stere'
    write(glimmap_proj4_stere(2),*) 'lon_0=',stere%longitude_of_central_meridian
    write(glimmap_proj4_stere(3),*) 'lat_0=',stere%latitude_of_projection_origin
    if (stere%pole/=0) then
       if (stere%standard_parallel /= 0) then
          write(glimmap_proj4_stere(4),*) 'lat_ts=',stere%standard_parallel
       else
          write(glimmap_proj4_stere(4),*) 'k_0=',stere%scale_factor_at_proj_origin
       end if
    else
       write(glimmap_proj4_stere(4),*) 'k_0=',stere%scale_factor_at_proj_origin
    end if
    write(glimmap_proj4_stere(5),*) 'x_0=',stere%false_easting
    write(glimmap_proj4_stere(6),*) 'y_0=',stere%false_northing
  end function glimmap_proj4_stere

  !> Returns a proj4 parameter string for a Lambert azimuthal equal area projection
  function glimmap_proj4_laea(laea)
    implicit none
    character(len=proj4len), dimension(:), pointer :: glimmap_proj4_laea
    type(proj_laea) :: laea

    allocate(glimmap_proj4_laea(5))
    write(glimmap_proj4_laea(1),*) 'proj=laea'
    write(glimmap_proj4_laea(2),*) 'lon_0=',laea%longitude_of_central_meridian
    write(glimmap_proj4_laea(3),*) 'lat_0=',laea%latitude_of_projection_origin
    write(glimmap_proj4_laea(4),*) 'x_0=',laea%false_easting
    write(glimmap_proj4_laea(5),*) 'y_0=',laea%false_northing
  end function glimmap_proj4_laea

  !> Returns a proj4 parameter string for a Lambert azimuthal equal area projection
  function glimmap_proj4_aea(aea)
    implicit none
    character(len=proj4len), dimension(:), pointer :: glimmap_proj4_aea
    type(proj_aea) :: aea

    allocate(glimmap_proj4_aea(7))
    write(glimmap_proj4_aea(1),*) 'proj=aea'
    write(glimmap_proj4_aea(2),*) 'lon_0=',aea%longitude_of_central_meridian
    write(glimmap_proj4_aea(3),*) 'lat_0=',aea%latitude_of_projection_origin
    write(glimmap_proj4_aea(4),*) 'lat_1=',aea%standard_parallel(1)
    write(glimmap_proj4_aea(5),*) 'lat_2=',aea%standard_parallel(2)
    write(glimmap_proj4_aea(6),*) 'x_0=',aea%false_easting
    write(glimmap_proj4_aea(7),*) 'y_0=',aea%false_northing
  end function glimmap_proj4_aea

  !> Returns a proj4 parameter string for a Lambert conformal conic projection
  function glimmap_proj4_lcc(lcc)
    implicit none
    character(len=proj4len), dimension(:), pointer :: glimmap_proj4_lcc
    type(proj_lcc) :: lcc

    allocate(glimmap_proj4_lcc(7))
    write(glimmap_proj4_lcc(1),*) 'proj=lcc'
    write(glimmap_proj4_lcc(2),*) 'lon_0=',lcc%longitude_of_central_meridian
    write(glimmap_proj4_lcc(3),*) 'lat_0=',lcc%latitude_of_projection_origin
    write(glimmap_proj4_lcc(4),*) 'lat_1=',lcc%standard_parallel(1)
    write(glimmap_proj4_lcc(5),*) 'lat_2=',lcc%standard_parallel(2)
    write(glimmap_proj4_lcc(6),*) 'x_0=',lcc%false_easting
    write(glimmap_proj4_lcc(7),*) 'y_0=',lcc%false_northing
  end function glimmap_proj4_lcc  

end module glimmer_map_proj4
