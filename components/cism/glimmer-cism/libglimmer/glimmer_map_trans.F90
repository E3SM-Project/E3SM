!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_map_trans.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!> convert between projections
module glimmer_map_trans

  use glimmer_map_types
  use glimmer_global, only: dp
  implicit none

  private
  public :: glimmap_ll_to_xy, glimmap_xy_to_ll, loncorrect


contains

  !> Convert lat-long coordinates to grid coordinates.
  !!
  !! The subroutine returns the x-y coordinates as real values,
  !! non-integer values indicating a position between grid-points.
  subroutine glimmap_ll_to_xy(lon,lat,x,y,proj,grid)

    use glimmer_log
    use glimmer_coordinates

    implicit none

    real(dp),intent(in)  :: lon !< The location of the point in lat-lon space (Longitude)
    real(dp),intent(in)  :: lat !< The location of the point in lat-lon space (Latitude)
    real(dp),intent(out) :: x   !< The location of the point in x-y space (x coordinate)
    real(dp),intent(out) :: y   !< The location of the point in x-y space (y coordinate)
    type(glimmap_proj),    intent(in) :: proj !< The projection being used
    type(coordsystem_type),intent(in) :: grid !< the grid definition

    real(dp) :: xx,yy ! These are real-space distances in meters

    if (associated(proj%laea)) then
       call glimmap_laea(lon,lat,xx,yy,proj%laea)
    else if (associated(proj%aea)) then
       call glimmap_aea(lon,lat,xx,yy,proj%aea)
    else if (associated(proj%lcc)) then
       call glimmap_lcc(lon,lat,xx,yy,proj%lcc)
    else if (associated(proj%stere)) then
       call glimmap_stere(lon,lat,xx,yy,proj%stere)
    else
       call write_log('No known projection found!',GM_WARNING)
    end if

    ! Now convert the real-space distances to grid-points using the grid type

    call space2grid(xx,yy,x,y,grid)

  end subroutine glimmap_ll_to_xy

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Convert grid coordinates to lat-lon coordinates. 
  !!
  !! The subroutine returns the lat-lon coordinates as real values,
  !! non-integer values indicating a position between grid-points.
  subroutine glimmap_xy_to_ll(lon,lat,x,y,proj,grid)


    use glimmer_log
    use glimmer_coordinates

    implicit none

    real(dp),intent(out) :: lon !< The location of the point in lat-lon space (Longitude)
    real(dp),intent(out) :: lat !< The location of the point in lat-lon space (Latitude)
    real(dp),intent(in)  :: x   !< The location of the point in x-y space (x coordinate)
    real(dp),intent(in)  :: y   !< The location of the point in x-y space (y coordinate)
    type(glimmap_proj),    intent(in) :: proj !< The projection being used
    type(coordsystem_type),intent(in) :: grid !< the grid definition

    real(dp) :: xx,yy ! These are real-space distances in meters

    ! First convert grid-point space to real space

    call grid2space(xx,yy,x,y,grid)

    if (associated(proj%laea)) then
       call glimmap_ilaea(lon,lat,xx,yy,proj%laea)
    else if (associated(proj%aea)) then
       call glimmap_iaea(lon,lat,xx,yy,proj%aea)
    else if (associated(proj%lcc)) then
       call glimmap_ilcc(lon,lat,xx,yy,proj%lcc)
    else if (associated(proj%stere)) then
       call glimmap_istere(lon,lat,xx,yy,proj%stere)
    else
       call write_log('No known projection found!',GM_WARNING)
    end if

    lon=loncorrect(lon,0.d0)

  end subroutine glimmap_xy_to_ll
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! PRIVATE subroutines follow
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Lambert azimuthal equal area projection
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Forward transformation: lat-lon -> x-y of Lambert azimuthal equal area projection
  subroutine glimmap_laea(lon,lat,x,y,params)


    use glimmer_log

    real(dp),intent(in)  :: lon !< longitude
    real(dp),intent(in)  :: lat !< latitude
    real(dp),intent(out) :: x   !< x
    real(dp),intent(out) :: y   !< y
    type(proj_laea),intent(in) :: params !< projection parameters

    real(dp) :: sin_lat,cos_lat,sin_lon,cos_lon,c,dlon,dlat,tmp,k
    character(80) :: errtxt

    dlon = lon-params%longitude_of_central_meridian

    ! Check domain of longitude

    dlon = loncorrect(dlon,-180.d0)

    ! Convert to radians and calculate sine and cos

    dlon = dlon*D2R ; dlat = lat*D2R

    call sincos(dlon,sin_lon,cos_lon);
    call sincos(dlat,sin_lat,cos_lat);
    c = cos_lat * cos_lon

    ! Mapping transformation

    tmp = 1.d0 + params%sinp * sin_lat + params%cosp * c

    if (tmp > 0.d0) then
       k = EQ_RAD * sqrt (2.d0 / tmp)
       x = k * cos_lat * sin_lon
       y = k * (params%cosp * sin_lat - params%sinp * c)
    else
       write(errtxt,*)'LAEA projection error:',lon,lat,params%latitude_of_projection_origin
       call write_log(trim(errtxt),GM_FATAL,__FILE__,__LINE__)
    endif

    ! Apply false eastings and northings

    x = x + params%false_easting
    y = y + params%false_northing

  end subroutine glimmap_laea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Inverse transformation: lat-lon -> x-y of Lambert azimuthal equal area projection
  subroutine glimmap_ilaea(lon,lat,x,y,params)

    use glimmer_log

    real(dp),intent(out) :: lon !< longitude
    real(dp),intent(out) :: lat !< latitude
    real(dp),intent(in)  :: x   !< x
    real(dp),intent(in)  :: y   !< y
    type(proj_laea),intent(in) :: params !< projection parameters

    real(dp) :: rho,c,sin_c,cos_c,xx,yy
    character(80) :: errtxt

    xx=x ; yy=y

    ! Account for false eastings and northings

    xx = xx - params%false_easting
    yy = yy - params%false_northing

    rho=hypot (xx,yy)

    if (abs(rho) < CONV_LIMIT) then
       ! If very near the centre of the map...
       lat = params%latitude_of_projection_origin
       lon = params%longitude_of_central_meridian
    else
       c = 2.d0 * asin(0.5d0 * rho * i_EQ_RAD)
       call sincos (c, sin_c, cos_c)
       lat = asin (cos_c * params%sinp + (yy * sin_c * params%cosp / rho)) * R2D
       select case(params%pole)
       case(1)
          lon = params%longitude_of_central_meridian + R2D * atan2 (xx, -yy)
       case(-1)
          lon = params%longitude_of_central_meridian + R2D * atan2 (xx, yy)
       case(0)
          lon = params%longitude_of_central_meridian + &
               R2D * atan2 (xx * sin_c, (rho * params%cosp * cos_c - yy * params%sinp * sin_c))
       case default
          write(errtxt,*)'Inverse LAEA projection error:',params%pole
          call write_log(trim(errtxt),GM_FATAL,__FILE__,__LINE__)
       end select
    endif

  end subroutine glimmap_ilaea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Albers equal area conic projection
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Forward transformation: lat-lon -> x-y of Albers equal area conic projection
  subroutine glimmap_aea(lon,lat,x,y,params)

    real(dp),intent(in)  :: lon !< longitude
    real(dp),intent(in)  :: lat !< latitude
    real(dp),intent(out) :: x   !< x
    real(dp),intent(out) :: y   !< y
    type(proj_aea),intent(in) :: params !< projection parameters

    real(dp) :: dlon,theta,sint,cost,rho

    dlon = lon-params%longitude_of_central_meridian

    ! Check domain of longitude

    dlon = loncorrect(dlon,-180.d0)
    theta = params%n * dlon * D2R
    call sincos(theta,sint,cost)

    rho = params%i_n*sqrt(params%c - 2.0*params%n*sin(lat*D2R))

    x = EQ_RAD * rho * sint
    y = EQ_RAD * (params%rho0_R - rho * cost)

    ! Apply false eastings and northings

    x = x + params%false_easting
    y = y + params%false_northing

  end subroutine glimmap_aea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Inverse transformation: lat-lon -> x-y of Albers equal area conic projection
  subroutine glimmap_iaea(lon,lat,x,y,params)

    real(dp),intent(out) :: lon !< longitude
    real(dp),intent(out) :: lat !< latitude
    real(dp),intent(in)  :: x   !< x
    real(dp),intent(in)  :: y   !< y
    type(proj_aea),intent(in) :: params !< projection parameters

    real(dp) :: xx,yy,rho,theta

    xx=x ; yy=y

    ! Account for false eastings and northings

    xx = xx - params%false_easting
    yy = yy - params%false_northing

    rho = sqrt(xx**2.d0 + (params%rho0 - yy)**2.d0)
    if (params%n > 0.d0) then
       theta = atan2(xx,(params%rho0-yy))
    else
       theta = atan2(-xx,(yy-params%rho0))
    end if

    lat = asin((params%c-(rho*params%n/EQ_RAD)**2.d0)*0.5d0*params%i_n)*R2D
    lon = params%longitude_of_central_meridian+R2D*theta*params%i_n

  end subroutine glimmap_iaea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Lambert conformal conic projection
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Forward transformation: lat-lon -> x-y of Lambert conformal conic projection
  subroutine glimmap_lcc(lon,lat,x,y,params)

    real(dp),intent(in)  :: lon !< longitude
    real(dp),intent(in)  :: lat !< latitude
    real(dp),intent(out) :: x   !< x
    real(dp),intent(out) :: y   !< y
    type(proj_lcc),intent(in) :: params !< projection parameters

    real(dp) :: dlon,rho,theta,sint,cost

    dlon = lon-params%longitude_of_central_meridian

    ! Check domain of longitude

    dlon = loncorrect(dlon,-180.d0)
    rho  = EQ_RAD * params%f/(tan(M_PI_4+lat*D2R/2.d0))**params%n
    theta = params%n*dlon*D2R
    call sincos(theta,sint,cost)

    x = rho * sint
    y = params%rho0 - rho * cost

    ! Apply false eastings and northings

    x = x + params%false_easting
    y = y + params%false_northing

  end subroutine glimmap_lcc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Inverse transformation: lat-lon -> x-y of Lambert conformal conic projection
  subroutine glimmap_ilcc(lon,lat,x,y,params)

    real(dp),intent(out) :: lon !< longitude
    real(dp),intent(out) :: lat !< latitude
    real(dp),intent(in)  :: x   !< x
    real(dp),intent(in)  :: y   !< y
    type(proj_lcc),intent(in) :: params !< projection parameters

    real(dp) :: xx,yy,rho,theta

    xx=x ; yy=y

    ! Account for false eastings and northings

    xx = xx - params%false_easting
    yy = yy - params%false_northing

    rho = sign(sqrt(xx**2.d0 + (params%rho0-yy)**2.d0),params%n)
    if (params%n > 0.d0) then
       theta = atan2(xx,(params%rho0-yy))
    else
       theta = atan2(-xx,(yy-params%rho0))
    end if

    if (abs(rho) < CONV_LIMIT) then
       lat = sign(real(90.d0,kind=dp),params%n)
    else
       lat = R2D * (2.d0 * atan((EQ_RAD*params%f/rho)**params%i_n) - M_PI_2)
    end if

    lon = params%longitude_of_central_meridian+R2D*theta*params%i_n

  end subroutine glimmap_ilcc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Stereographic projection
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Forward transformation: lat-lon -> x-y of Stereographic projection
  subroutine glimmap_stere(lon,lat,x,y,params)

    use glimmer_log

    real(dp),intent(in)  :: lon !< longitude
    real(dp),intent(in)  :: lat !< latitude
    real(dp),intent(out) :: x   !< x
    real(dp),intent(out) :: y   !< y
    type(proj_stere),intent(in) :: params !< projection parameters

    real(dp) :: dlon,k,dlat,slat,clat,slon,clon
    character(80) :: errtxt

    dlon = lon-params%longitude_of_central_meridian

    ! Check domain of longitude

    dlon = loncorrect(dlon,-180.d0)
    dlon = dlon * D2R
    dlat = lat  * D2R
    call sincos(dlon,slon,clon)

    select case(params%pole)
    case(1)  ! North pole
       x =  2.d0 * params%k0 * tan(M_PI_4 - dlat/2.d0)*slon
       y = -2.d0 * params%k0 * tan(M_PI_4 - dlat/2.d0)*clon
    case(-1)  ! South pole
       x = 2.d0 * params%k0 * tan(M_PI_4 + dlat/2.d0)*slon
       y = 2.d0 * params%k0 * tan(M_PI_4 + dlat/2.d0)*clon
    case(0)  ! Oblique
       call sincos(dlat,slat,clat)
       if (params%equatorial) then
          k = 2.d0 * params%k0 / (1.d0 + clat*clon)
          y = k * slat
       else
          k = 2.d0 * params%k0 / (1.d0 + params%sinp*slat + params%cosp*clat*clon)
          y = k * (params%cosp*slat - params%sinp*clat*clon)
       end if
       x = k * clat * slon
    case default 
       write(errtxt,*)'Stereographic projection error:',params%pole
       call write_log(trim(errtxt),GM_FATAL,__FILE__,__LINE__)
    end select

    ! Apply false eastings and northings

    x = x + params%false_easting
    y = y + params%false_northing

  end subroutine glimmap_stere

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Inverse transformation: lat-lon -> x-y of Stereographic projection
  subroutine glimmap_istere(lon,lat,x,y,params)

    real(dp),intent(out) :: lon !< longitude
    real(dp),intent(out) :: lat !< latitude
    real(dp),intent(in)  :: x   !< x
    real(dp),intent(in)  :: y   !< y
    type(proj_stere),intent(in) :: params !< projection parameters

    real(dp) :: xx,yy,rho,c,sinc,cosc

    xx=x ; yy=y

    ! Account for false eastings and northings

    xx = xx - params%false_easting
    yy = yy - params%false_northing

    rho = hypot(xx,yy)

    if (abs(rho)<CONV_LIMIT) then
       lon = params%longitude_of_central_meridian
       lat = params%latitude_of_projection_origin
    else
       c = 2.d0 * atan(rho*0.5d0*params%ik0)
       call sincos(c,sinc,cosc)
       select case(params%pole)
       case(0)
          lat = asin(cosc*params%sinp+(yy*sinc*params%cosp/rho))*R2D
          lon = params%longitude_of_central_meridian + &
               atan2(xx * sinc,(rho*params%cosp*cosc-yy*params%sinp*sinc))*R2D
       case(1)
          lat = asin(cosc)*R2D
          lon = params%longitude_of_central_meridian + &
               atan2(xx,-yy) * R2D
       case(-1)
          lat = asin(-cosc)*R2D
          lon = params%longitude_of_central_meridian + &
               atan2(xx,yy) * R2D
       end select
    end if

  end subroutine glimmap_istere

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Utility routines
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> transform from grid to space
  subroutine grid2space(x,y,gx,gy,coordsys)

    use glimmer_coordinates

    implicit none

    real(dp),intent(out) :: x  !< x-location in real space
    real(dp),intent(out) :: y  !< y-location in real space
    real(dp),intent(in)  :: gx !< x-location in grid space
    real(dp),intent(in)  :: gy !< y-location in grid space
    type(coordsystem_type), intent(in) :: coordsys  !< coordinate system

    x=coordsys%origin%pt(1) + real(gx - 1)*coordsys%delta%pt(1)
    y=coordsys%origin%pt(2) + real(gy - 1)*coordsys%delta%pt(2)

  end subroutine grid2space

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> convert from space to grid
  subroutine space2grid(x,y,gx,gy,coordsys)

    use glimmer_coordinates

    implicit none

    real(dp),intent(in)  :: x  !< x-location in real space
    real(dp),intent(in)  :: y  !< y-location in real space
    real(dp),intent(out) :: gx !< x-location in grid space
    real(dp),intent(out) :: gy !< y-location in grid space
    type(coordsystem_type), intent(in) :: coordsys  !< coordinate system

    gx = 1.d0 + (x - coordsys%origin%pt(1))/coordsys%delta%pt(1)
    gy = 1.d0 + (y - coordsys%origin%pt(2))/coordsys%delta%pt(2)

  end subroutine space2grid

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Calculates the sin and cos of an angle.
  subroutine sincos(a,s,c)

    implicit none

    real(dp),intent(in)  :: a !< Input value (radians).
    real(dp),intent(out) :: s !< sin(a)
    real(dp),intent(out) :: c !< cos(a)

    s = sin(a)
    c = cos(a)

  end subroutine sincos

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Normalises a value of longitude to the range starting at min degrees.
  !! \return The normalised value of longitude.
  real(dp) function loncorrect(lon,minimum)

    real(dp),intent(in) :: lon     !< The longitude under consideration (degrees east)
    real(dp),intent(in) :: minimum !< The lower end of the output range (degrees east)

    real(dp) :: maximum

    loncorrect = lon
    maximum = minimum + 360.d0

    do while (loncorrect >= maximum)
       loncorrect = loncorrect-360.d0
    enddo

    do while (loncorrect < minimum)
       loncorrect = loncorrect+360.d0
    enddo

  end function loncorrect

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> compute \f$\sqrt{x^2+y^2}\f$
  real(dp) function hypot(x,y)


    implicit none

    real(dp),intent(in) :: x !< One input value
    real(dp),intent(in) :: y !< Another input value

    hypot=sqrt(x*x+y*y)

  end function hypot

end module glimmer_map_trans
