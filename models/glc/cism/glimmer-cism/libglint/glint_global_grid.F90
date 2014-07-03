!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_global_grid.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_global_grid

  use glimmer_global
  use glimmer_physcon, only: pi

  implicit none

  real(rk),parameter :: pi2 = 2.0*pi

  ! ------------------------------------------------------------
  ! GLOBAL_GRID derived type
  ! ------------------------------------------------------------

  type global_grid

     !*FD Contains parameters specifying the global grid configuration.

     ! Dimensions of grid ---------------------------------------

     integer :: nx = 0  !*FD Number of points in the $x$-direction.
     integer :: ny = 0  !*FD Number of points in the $y$-direction.
     integer :: nec = 1 !*FD Number of elevation classes

     ! Locations of grid-points ---------------------------------

     real(rk),pointer,dimension(:) :: lats      => null() 
     !*FD Latitudinal locations of data-points in global fields (degrees)
     real(rk),pointer,dimension(:) :: lons      => null() 
     !*FD Longitudinal locations of data-points in global fields (degrees)

     ! Locations of grid-box boundaries -------------------------

     real(rk),pointer,dimension(:) :: lat_bound => null() 
     !*FD Latitudinal boundaries of data-points in global fields (degrees)
     real(rk),pointer,dimension(:) :: lon_bound => null() 
     !*FD Longitudinal boundaries of data-points in global fields (degrees)

     ! Areas of grid-boxes --------------------------------------

     real(rk),pointer,dimension(:,:) :: box_areas => null() 
     !*FD The areas of the grid-boxes (m$^2$). This is a two-dimensional array to take
     !*FD account of the possibility of a grid irregularly spaced in longitude

     integer,pointer,dimension(:,:) :: mask => null()
     !*FD This mask = 1 where the global data are valid, = 0 elsewhere
     !*FD (e.g., over ocean, where the GCM does not compute a surface mass balance)

  end type global_grid

  private pi  !TODO - Remove this line?

  interface min
     module procedure grid_min
  end interface

  interface operator(==)
     module procedure grid_equiv
  end interface

  interface operator(/=)
     module procedure grid_nequiv
  end interface

  interface operator(>)
     module procedure grid_greater_than
  end interface

  interface operator(<)
     module procedure grid_less_than
  end interface

  interface grid_alloc
     module procedure grid_alloc_2d,grid_alloc_3d
  end interface

contains

  subroutine new_global_grid(grid, lons, lats, lonb, latb, radius, correct, nec, mask)

    use glimmer_log

    !*FD Initialises a new global grid type

    type(global_grid),              intent(inout) :: grid !*FD The grid to be initialised
    real(rk),dimension(:),          intent(in)    :: lons !*FD Longitudinal positions of grid-points (degrees)
    real(rk),dimension(:),          intent(in)    :: lats !*FD Latitudinal positions of grid-points (degrees)
    real(rk),dimension(:), optional,intent(in)    :: lonb !*FD Longitudinal boundaries of grid-boxes (degrees)
    real(rk),dimension(:), optional,intent(in)    :: latb !*FD Latitudinal boundaries of grid-boxes (degrees)
    real(rk),              optional,intent(in)    :: radius  !*FD  The radius of the Earth (m)
    logical,               optional,intent(in)    :: correct !*FD  Set to correct for boundaries (default is .true.)
    integer,               optional,intent(in)    :: nec     !*FD  Number of elevation classes    
    integer,dimension(:,:),optional,intent(in)    :: mask    !*FD  Mask indicating where global data are valid 

    ! Internal variables

    real(rk) :: radea=1.0
    integer :: i,j
    logical :: cor

    ! Deal with optional non-correction

    if (present(correct)) then
       cor=correct
    else
       cor=.true.
    end if

    ! Check to see if things are allocated, and if so, deallocate them

    if (associated(grid%lats))      deallocate(grid%lats)
    if (associated(grid%lons))      deallocate(grid%lons)
    if (associated(grid%lat_bound)) deallocate(grid%lat_bound)
    if (associated(grid%lon_bound)) deallocate(grid%lon_bound)
    if (associated(grid%box_areas)) deallocate(grid%box_areas)
    if (associated(grid%mask))      deallocate(grid%mask)

    ! Find size of grid

    grid%nx=size(lons) ; grid%ny=size(lats)

    ! Allocate arrays

    allocate(grid%lons(grid%nx))
    allocate(grid%lats(grid%ny))
    allocate(grid%lon_bound(grid%nx+1))
    allocate(grid%lat_bound(grid%ny+1))
    allocate(grid%box_areas(grid%nx,grid%ny))
    allocate(grid%mask(grid%nx,grid%ny))

    ! Check dimensions of boundary arrays, if supplied

    if (present(lonb)) then
       if (.not.size(lonb)==grid%nx+1) then
          call write_log('Lonb mismatch in new_global_grid',GM_FATAL,__FILE__,__LINE__)
       endif
    endif

    if (present(latb)) then
       if (.not.size(latb)==grid%ny+1) then
          call write_log('Latb mismatch in new_global_grid',GM_FATAL,__FILE__,__LINE__)
       endif
    endif

    ! Copy lats and lons over

    grid%lats=lats
    grid%lons=lons

    ! In the following code, there are some things that only work correctly if 
    ! size(grid%lats) > 0 or size(grid%lons) > 0... we need to explicitly check those
    ! conditions to handle the multiple processor case, where tasks other than main have
    ! a size-0 global grid

    ! Check to see if we have polar points, and fudge if necessary
    ! By moving 1/20 of the way towards the next equatorward point.

    if (size(grid%lats) > 0) then
       if (abs(grid%lats(1)-90)<1e-8)       grid%lats(1)       =  90.0-(grid%lats(1)        -grid%lats(2))/20.0
       if (abs(grid%lats(grid%ny)+90)<1e-8) grid%lats(grid%ny) = -90.0+(grid%lats(grid%ny-1)-grid%lats(grid%ny))/20.0
    end if

    ! Calculate boundaries if necessary

    if (size(grid%lons) > 0) then
       if (present(lonb)) then
          grid%lon_bound=lonb
       else
          call calc_bounds_lon(lons,grid%lon_bound,cor)
       endif
    end if

    if (size(grid%lats) > 0) then
       if (present(latb)) then
          grid%lat_bound=latb
       else
          call calc_bounds_lat(lats,grid%lat_bound)
       endif
    end if

    ! Set radius of earth if necessary

    if (present(radius)) radea=radius

    ! Calculate areas of grid-boxes

    do i=1,grid%nx
       do j=1,grid%ny
          grid%box_areas(i,j)=delta_lon(grid%lon_bound(i),grid%lon_bound(i+1))*radea**2* &
               (sin_deg(grid%lat_bound(j))-sin_deg(grid%lat_bound(j+1)))
       enddo
    enddo

    ! Specify mask

    if (present(mask)) then
       grid%mask(:,:) = mask(:,:)
    else
       grid%mask(:,:) = 1   ! assume global data are valid everywhere 
    endif

    ! Set number of elevation classes

    if (present(nec)) then
       grid%nec = nec
    else
       grid%nec = 1
    endif 

  end subroutine new_global_grid

  !-----------------------------------------------------------------------------

  subroutine copy_global_grid(in,out)

    !*FD Copies a global grid type.

    type(global_grid),intent(in)  :: in  !*FD Input grid
    type(global_grid),intent(out) :: out !*FD Output grid

    ! Copy dimensions

    out%nx  = in%nx
    out%ny  = in%ny
    out%nec = in%nec

    ! Check to see if arrays are allocated, then deallocate and
    ! reallocate accordingly.

    if (associated(out%lats))      deallocate(out%lats)
    if (associated(out%lons))      deallocate(out%lons)
    if (associated(out%lat_bound)) deallocate(out%lat_bound)
    if (associated(out%lon_bound)) deallocate(out%lon_bound)
    if (associated(out%box_areas)) deallocate(out%box_areas)
    if (associated(out%mask))      deallocate(out%mask)


    allocate(out%lons(out%nx))
    allocate(out%lats(out%ny))
    allocate(out%lon_bound(out%nx+1))
    allocate(out%lat_bound(out%ny+1))
    allocate(out%box_areas(out%nx,out%ny))
    allocate(out%mask(out%nx,out%ny))

    ! Copy data

    out%lons     =in%lons
    out%lats     =in%lats
    out%lat_bound=in%lat_bound
    out%lon_bound=in%lon_bound
    out%box_areas=in%box_areas
    out%mask     =in%mask

  end subroutine copy_global_grid

  !-----------------------------------------------------------------------------

  subroutine get_grid_dims(grid,nx,ny,nec)

    type(global_grid),intent(in)  :: grid
    integer,intent(out) :: nx,ny
    integer,intent(out), optional :: nec

    nx=grid%nx
    ny=grid%ny

    if (present(nec)) nec=grid%nec
 
  end subroutine get_grid_dims

  !-----------------------------------------------------------------------------

  subroutine print_grid(grid,no_gp)

    type(global_grid),intent(in)  :: grid
    logical,optional,intent(in) :: no_gp

    logical :: ng

    if (present(no_gp)) then
       ng=no_gp
    else
       ng=.false.
    end if

    print*,'Grid parameters:'
    print*,'----------------'
    print*
    print*,'nx=', grid%nx
    print*,'ny=', grid%ny
    print*,'nec=',grid%nec
    if (.not.ng) then
       print*
       print*,'longitudes:'
       print*,grid%lons
       print*
       print*,'latitudes:'
       print*,grid%lats
       print*
       print*,'longitude boundaries:'
       print*,grid%lon_bound
       print*
       print*,'latitude boundaries:'
       print*,grid%lat_bound
    end if

  end subroutine print_grid

  !-----------------------------------------------------------------------------

  subroutine calc_bounds_lon(lons,lonb,correct)

    !*FD Calculates the longitudinal boundaries between
    !*FD global grid-boxes. Note that we assume that the boundaries lie 
    !*FD half-way between the points, although 
    !*FD this isn't strictly true for a Gaussian grid.

    implicit none

    real(rk),dimension(:),intent(in)  :: lons    !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: lonb    !*FD boundaries of grid-boxes (degrees)
    logical,              intent(in)  :: correct !*FD Set to correct for longitudinal grid boundary

    integer :: nxg,i

    nxg=size(lons)

    ! Longitudes

    do i=1,nxg-1
       lonb(i+1)=mid_lon(lons(i),lons(i+1),correct)
    enddo

    lonb(1)=mid_lon(lons(nxg),lons(1),correct)
    lonb(nxg+1)=lonb(1)

  end subroutine calc_bounds_lon

  !---------------------------------------------------------------------------------

  subroutine calc_bounds_lat(lat,latb)

    !*FD Calculates the boundaries between
    !*FD global grid-boxes. Note that we assume that the boundaries lie 
    !*FD half-way between the 
    !*FD points, both latitudinally and longitudinally, although 
    !*FD this isn't strictly true for a Gaussian grid.

    implicit none

    real(rk),dimension(:),intent(in)  :: lat  !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: latb !*FD boundaries of grid-boxes (degrees)

    integer :: nyg,j

    nyg=size(lat)

    ! Latitudes first - we assume the boundaries of the first and 
    ! last boxes coincide with the poles. Not sure how to
    ! handle it if they don't...

    latb(1)=90.0
    latb(nyg+1)=-90.0

    do j=2,nyg
       latb(j)=lat(j-1)-(lat(j-1)-lat(j))/2.0
    enddo

  end subroutine calc_bounds_lat

  !-------------------------------------------------------------

  real(rk) function mid_lon(a,b,correct)

    use glimmer_log

    !*FD Calculates the mid-point between two longitudes.
    !*FD \texttt{a} must be west of \texttt{b}.

    real(rk),intent(in) :: a,b
    logical :: correct

    real(rk) :: aa,bb,out

    aa=a ; bb=b

    if (aa>360.0.or.aa<0.0  .or. &
         bb>360.0.or.bb<0.0) then
       call write_log('Out of range in mid_lon',GM_FATAL,__FILE__,__LINE__)
    endif

    if (aa>bb) aa=aa-360.0

    out=aa+((bb-aa)/2.0)

    if (correct) then
       do
          if (out<=360.0) exit
          out=out-360.0
       end do

       do
          if (out>=0.0) exit
          out=out+360.0
       end do
    end if

    mid_lon=out

  end function mid_lon

  !-------------------------------------------------------------

  real(rk) function sin_deg(a)

    !*FD Calculate sin(a), where a is in degrees

    real(rk) :: a
    real(rk) :: aa

    aa=pi*a/180.0

    do
       if (aa<=pi2) exit
       aa=aa-pi2
    end do

    do
       if (aa>=0.0) exit
       aa=aa+pi2
    end do

    sin_deg=sin(aa)

  end function sin_deg

  !-------------------------------------------------------------

  real(rk) function delta_lon(a,b)

    real(rk) :: a,b
    real(rk) :: aa,bb,dl

    aa=a ; bb=b

    do
       dl=bb-aa
       if (dl>=0.0) exit
       aa=aa-360.0
    end do

    delta_lon=dl*pi/180.0

  end function delta_lon

  !-------------------------------------------------------------

  logical function grid_greater_than(a,b)

    type(global_grid),intent(in) :: a,b

    if (a%nx*a%ny > b%nx*b%ny) then
       grid_greater_than=.true.
    else
       grid_greater_than=.false.
    end if

  end function grid_greater_than

  !-------------------------------------------------------------

  logical function grid_less_than(a,b)

    type(global_grid),intent(in) :: a,b

    if (a%nx*a%ny < b%nx*b%ny) then
       grid_less_than=.true.
    else
       grid_less_than=.false.
    end if

  end function grid_less_than

  !-------------------------------------------------------------

  function grid_min(a,b)

    type(global_grid),intent(in) :: a,b
    type(global_grid) :: grid_min

    if (a>b) then
       grid_min=b
    else
       grid_min=a
    endif

  end function grid_min

  !-------------------------------------------------------------

  logical function grid_equiv(a,b)

    type(global_grid),intent(in)  :: a,b

    if (.not.check_associated(a).or. &
         .not.check_associated(b)) then
       grid_equiv=.false.
       return
    end if

    if (a%nx /= b%nx .or. a%ny /= b%ny .or. a%nec /= b%nec) then
       grid_equiv=.false.
       return
    end if

    ! N.B. Only checks grid-box centres and not boundaries

    if (any(a%lats /= b%lats).or. &
        any(a%lons /= b%lons).or. &
        any(a%mask /= b%mask).or. &
        any(a%box_areas /= b%box_areas)) then
       grid_equiv=.false.
       return
    end if

    grid_equiv=.true.

  end function grid_equiv

  !-------------------------------------------------------------

  logical function grid_nequiv(a,b)

    type(global_grid),intent(in)  :: a,b

    grid_nequiv=.not.grid_equiv(a,b)

  end function grid_nequiv

  !-------------------------------------------------------------

  logical function check_associated(a)

    type(global_grid),intent(in)  :: a

    if (associated(a%lats).and. &
         associated(a%lons).and. &
         associated(a%lat_bound).and. &
         associated(a%lon_bound).and. &
         associated(a%mask).and. &
         associated(a%box_areas)) then
       check_associated=.true.
    else
       check_associated=.false.
    end if

  end function check_associated

  !-------------------------------------------------------------

  subroutine grid_alloc_3d(array,grid,d3)

    real(rk),dimension(:,:,:),pointer :: array
    type(global_grid),intent(in)  :: grid
    integer,intent(in) :: d3

    if (associated(array)) deallocate(array)

    allocate(array(grid%nx,grid%ny,d3))

  end subroutine grid_alloc_3d

  !--------------------------------------------------------------

  subroutine grid_alloc_2d(array,grid)

    real(rk),dimension(:,:),pointer :: array
    type(global_grid),intent(in)  :: grid

    if (associated(array)) deallocate(array)

    allocate(array(grid%nx,grid%ny))

  end subroutine grid_alloc_2d

end module glint_global_grid
