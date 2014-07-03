!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_searchcircle.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!> improved algorithm for integrating a 2 dimensional array over large circles
!! this used for calculating continentality
!!
!! \author Magnus Hagdorn
module searchcircle
  
  implicit none

  type searchdata
     logical :: initialised = .false.
     integer :: radius                              !< search radius
     integer, pointer, dimension(:) :: ipos         !< positions on quater circle which will be moved along
     integer :: istart, jstart                      !< starting position of grid to be processed, will default to usally 1
     integer :: isize, jsize                        !< size of array to be processed
     real :: total_area
     real, pointer, dimension(:,:) :: sarray        !< array to be searched (expanded to include outside points
     real, pointer, dimension(:,:) :: weight        !< reciprocal weights
  end type searchdata

  !MAKE_RESTART

contains

  !> initialise search circle data structure
  !!
  !! \return initialised data type
  function sc_initdata(radius,istart,jstart,isize,jsize,searchgrid)
    implicit none
    integer, intent(in) :: radius                !< radius of search radius
    integer, intent(in) :: istart,jstart         !< starting position of grid to be processed
    integer, intent(in) :: isize,jsize           !< size of array to be processed
    real, dimension(:,:), optional :: searchgrid !< used for determining bounds of grid to be searched if not present, the bounds are assumed to be the same as the resultgrid
    
    type(searchdata) :: sc_initdata

    ! local variables
    real, allocatable, dimension(:) :: area
    real ::  area_temp
    integer i,j,intrad,ii,jj
    integer si_start,si_end,sj_start,sj_end,si_size,sj_size

    ! filling structure
    sc_initdata%radius = radius
    sc_initdata%istart = istart
    sc_initdata%jstart = jstart
    sc_initdata%isize  = isize
    sc_initdata%jsize  = jsize
    ! allocating data
    allocate(sc_initdata%sarray(1-radius:isize+radius,1-radius:jsize+radius))
    allocate(sc_initdata%weight(isize, jsize))
    allocate(sc_initdata%ipos(radius+1))

    if (present(searchgrid)) then
       si_start = lbound(searchgrid,1)
       sj_start = lbound(searchgrid,2)
       si_end = ubound(searchgrid,1)
       sj_end = ubound(searchgrid,2)
       si_size = si_end - si_start + 1
       sj_size = sj_end - sj_start + 1
    else
       si_start = 1
       sj_start = 1
       si_end = isize
       sj_end = jsize
       si_size = isize
       sj_size = jsize
    end if

    ! initialising data
    ! mask
    sc_initdata%sarray = 0.
    sc_initdata%ipos = 0
    ! weights
    ! calculate integral over quater circle
    allocate(area(0:radius))
    area(0) = radius    
    do i=1,radius
       sc_initdata%ipos(i) = int(sqrt(real(radius*radius-i*i)))
       area(i) = area(i-1)+real(sc_initdata%ipos(i))
    end do
    sc_initdata%total_area = 1.+4.*area(radius)

    ! complaining if search circle does not fit
    if (si_size < 2*radius+2 .and. sj_size < 2*radius+2) then
       ! internal sums
       sc_initdata%weight(1+radius:isize-radius, 1+radius:jsize-radius) = sc_initdata%total_area
       do j=jstart,jstart+jsize-1
          !left
          do i=istart,istart+radius-1
             area_temp = 0.
             do jj=max(sj_start,j-radius),min(sj_end,j+radius)
                intrad = int(sqrt(real(radius*radius-(jj-j)*(jj-j))))
                do ii=max(si_start,i-intrad),min(si_end,i+intrad)
                   area_temp = area_temp + 1.
                end do
             end do
             sc_initdata%weight(i-istart+1,j-jstart+1) = area_temp
          end do
          !right
          do i=istart+isize-1-radius,istart+isize-1
             area_temp = 0.
             do jj=max(sj_start,j-radius),min(sj_end,j+radius)
                intrad = int(sqrt(real(radius*radius-(jj-j)*(jj-j))))
                do ii=max(si_start,i-intrad),min(si_end,i+intrad)
                   area_temp = area_temp + 1.
                end do
             end do
             sc_initdata%weight(i-istart+1,j-jstart+1) = area_temp
          end do
       end do
       ! lower
       do j=jstart,jstart+radius-1
          do i=istart+radius,istart+isize-1-radius
             area_temp = 0.
             do jj=max(sj_start,j-radius),min(sj_end,j+radius)
                intrad = int(sqrt(real(radius*radius-(jj-j)*(jj-j))))
                do ii=max(si_start,i-intrad),min(si_end,i+intrad)
                   area_temp = area_temp + 1.
                end do
             end do
             sc_initdata%weight(i-istart+1,j-jstart+1) = area_temp
          end do
       end do
       ! upper
       do j=jstart+jsize-1-radius,jstart+jsize-1
          do i=istart+radius,istart+isize-1-radius
             area_temp = 0.
             do jj=max(sj_start,j-radius),min(sj_end,j+radius)
                intrad = int(sqrt(real(radius*radius-(jj-j)*(jj-j))))
                do ii=max(si_start,i-intrad),min(si_end,i+intrad)
                   area_temp = area_temp + 1.
                end do
             end do
             sc_initdata%weight(i-istart+1,j-jstart+1) = area_temp
          end do
       end do
    else
       do j=jstart,jstart+jsize-1
          do i=istart,istart+isize-1
             area_temp = 0.
             do jj=max(sj_start,j-radius),min(sj_end,j+radius)
                intrad = int(sqrt(real(radius*radius-(jj-j)*(jj-j))))
                do ii=max(si_start,i-intrad),min(si_end,i+intrad)
                   area_temp = area_temp + 1.
                end do
             end do
             sc_initdata%weight(i-istart+1,j-jstart+1) = area_temp
          end do
       end do
    end if

    sc_initdata%weight = sc_initdata%total_area/sc_initdata%weight

    sc_initdata%initialised = .true.
  end function sc_initdata


  !> do the search
  !!
  !! \bug cony does not match at boundary. no idea what is going on...
  subroutine sc_search(sdata,searchgrid,resultgrid)
    implicit none
    type(searchdata) :: sdata !< the search circle type
    real, dimension(:,:), intent(in) :: searchgrid !< the input mesh
    real, dimension(:,:), intent(out) :: resultgrid !< the result mesh

    ! local variables
    integer i,j,ii,jj,intrad
    integer :: istart,iend,jstart,jend

    if (.not.sdata%initialised) then
       write(*,*) 'Error (searchcircle), module is not initialised'
       stop
    end if

    ! checking grid sizes
    if (any(shape(resultgrid) /= (/sdata%isize,sdata%jsize/))) then
       write(*,*) 'Error (searchcircle), size of result grid does not match: ',shape(resultgrid),(/sdata%isize,sdata%jsize/)
       stop
    end if
    !filling search array
    sdata%sarray = 0.
    istart = max(1, sdata%istart-sdata%radius)
    iend   = min(size(searchgrid,1),sdata%istart+sdata%isize+sdata%radius-1)
    jstart = max(1, sdata%jstart-sdata%radius)
    jend   = min(size(searchgrid,2),sdata%jstart+sdata%jsize+sdata%radius-1)

    sdata%sarray(1+istart-sdata%istart:iend-sdata%istart+1, 1+jstart-sdata%jstart:jend-sdata%jstart+1) = &
         searchgrid(istart:iend, jstart:jend)
    resultgrid = 0.

    ! loop over grid
    do j=1,sdata%jsize
       ! do the full circle
       i=1  
       do jj=j-sdata%radius,j+sdata%radius
          intrad = int(sqrt(real(sdata%radius*sdata%radius-(jj-j)*(jj-j))))
          do ii=i-intrad,i+intrad
             resultgrid(i,j) = resultgrid(i,j) + sdata%sarray(ii,jj)
          end do
       end do
       
       ! loop over the remaing columns in the current row
       do i=2,sdata%isize
          resultgrid(i,j) = resultgrid(i-1,j) - sdata%sarray(i-sdata%radius,j) + sdata%sarray(i+sdata%radius,j)
          do jj=1,sdata%radius
             resultgrid(i,j) = resultgrid(i,j) - sdata%sarray(i-sdata%ipos(jj),j+jj) + sdata%sarray(i+sdata%ipos(jj),j+jj)
             resultgrid(i,j) = resultgrid(i,j) - sdata%sarray(i-sdata%ipos(jj),j-jj) + sdata%sarray(i+sdata%ipos(jj),j-jj)
          end do
       end do
    end do
    
    ! applying weights
    resultgrid = resultgrid * sdata%weight
  end subroutine sc_search

end module searchcircle
