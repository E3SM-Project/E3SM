!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_ts.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!> handling time series
!!
!! \author Magnus Hagdorn
!! \date 2006
!!
!! this module provides support for reading in tabulated ASCII data such as
!! time series. The data can then be accessed through functions which
!! interpolated the data.
module glimmer_ts

  use glimmer_global, only: dp
  implicit none

  !> time series derived type
  type glimmer_tseries
     integer :: numt=0                              !< number of times in time series
     integer :: numv=1                              !< number of values per time
     integer :: current=1                           !< current position in ts
     real(dp), dimension(:), pointer :: times=>NULL()   !< array holding times
     real(dp), dimension(:,:), pointer :: values=>NULL()!< array holding values
  end type glimmer_tseries

  interface glimmer_ts_step
     module procedure glimmer_ts_step_array, glimmer_ts_step_scalar
  end interface

  interface glimmer_ts_linear
     module procedure glimmer_ts_linear_array,glimmer_ts_linear_scalar
  end interface

  private :: get_i

contains

  !> read tabulated ASCII file
  subroutine glimmer_read_ts(ts,fname,numv)
    use glimmer_log
    use glimmer_global, only : msg_length
    implicit none
    type(glimmer_tseries) :: ts           !< time series data
    character(len=*), intent(in) :: fname !< read from this file
    integer, intent(in),optional :: numv  !< number of values per time
    
    ! local variables
    real(dp) :: d1,d2,fact=1.
    integer i,j,ios
    character(len=msg_length) :: message

    if (present(numv)) then
       ts%numv = numv
    else
       ts%numv = 1
    end if

    open(99,file=trim(fname),status='old',iostat=ios)
    
    if (ios.ne.0) then
       call write_log('Error opening file: '//trim(fname),type=GM_FATAL)
    end if

    ! find number of times and checking if ts is strictly monotonic
    ios = 0
    d1 = 1
    do
       d2 = d1
       read(99,*,iostat=ios) d1
       d1 = fact*d1
       if (ios.ne.0) then
          exit
       end if
       ts%numt = ts%numt + 1
       if (ts%numt.eq.1) then
          cycle
       else if (ts%numt.eq.2) then
          if (d1 > d2) then
             fact = 1.
          else if (d1 < d2) then
             fact = -1.
             d1 = -d1
          else
             write(message,*) 'Error, time series in file: '//trim(fname)//' is not monotonic line: ',ts%numt
             call write_log(message,type=GM_FATAL)
          end if
       else
          if (d1 <= d2) then
             write(message,*) 'Error, time series in file: '//trim(fname)//' is not monotonic line: ',ts%numt
             call write_log(message,type=GM_FATAL)
          end if
       end if
    end do
    rewind(99)

    allocate(ts%times(ts%numt))
    allocate(ts%values(ts%numv,ts%numt))
    ! read data
    do i=1,ts%numt
       read(99,*) ts%times(i),(ts%values(j,i),j=1,ts%numv)
    end do
    close(99)
  end subroutine glimmer_read_ts

  !> interpolate time series by stepping
  subroutine glimmer_ts_step_array(ts,time,value)
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !< time series data
    real(dp), intent(in)      :: time   !< time value to get
    real(dp), dimension(:)    :: value  !< interpolated value
       
    integer i

    if (size(value).ne.ts%numv) then
       call write_log('Error, wrong number of values',GM_FATAL,__FILE__,__LINE__)
    end if

    i = get_i(ts,time)
    if (i.eq.-1) then
       i = 1
    else if (i.eq.ts%numt+1) then
       i = ts%numt
    end if

    value = ts%values(:,i)
  end subroutine glimmer_ts_step_array

  !> interpolate time series by stepping
  subroutine glimmer_ts_step_scalar(ts,time,value)
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !< time series data
    real(dp), intent(in)      :: time   !< time value to get
    real(dp)                  :: value  !< interpolated value
       
    integer i

    i = get_i(ts,time)
    if (i.eq.-1) then
       i = 1
    else if (i.eq.ts%numt+1) then
       i = ts%numt
    end if

    value = ts%values(1,i)
  end subroutine glimmer_ts_step_scalar
 
  !> linear interpolate time series 
  subroutine glimmer_ts_linear_array(ts,time,value)
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !< time series data
    real(dp), intent(in)      :: time   !< time value to get
    real(dp), dimension(:)    :: value  !< interpolated value
       
    integer i
    real(dp),dimension(size(value)) :: slope

    if (size(value).ne.ts%numv) then
       call write_log('Error, wrong number of values',GM_FATAL,__FILE__,__LINE__)
    end if
    
    i = get_i(ts,time)
    if (i.eq.-1) then
       value(:) = ts%values(:,1)
    else if (i.eq.ts%numt+1) then
       value(:) = ts%values(:,ts%numt)
    else
       slope(:) = (ts%values(:,i+1)-ts%values(:,i))/(ts%times(i+1)-ts%times(i))
       value(:) = ts%values(:,i) + slope(:)*(time-ts%times(i))
    end if
  end subroutine glimmer_ts_linear_array

  !> linear interpolate time series
  subroutine glimmer_ts_linear_scalar(ts,time,value)
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !< time series data
    real(dp), intent(in)      :: time   !< time value to get
    real(dp)                  :: value  !< interpolated value
       
    integer i
    real(dp) :: slope

    i = get_i(ts,time)
    if (i.eq.-1) then
       value = ts%values(1,1)
    else if (i.eq.ts%numt+1) then
       value = ts%values(1,ts%numt)
    else
       slope = (ts%values(1,i+1)-ts%values(1,i))/(ts%times(i+1)-ts%times(i))
       value = ts%values(1,i) + slope*(time-ts%times(i))
    end if
  end subroutine glimmer_ts_linear_scalar
  
  !> get find the index
  function get_i(ts,time)
    implicit none
    type(glimmer_tseries) :: ts     !< time series data
    real(dp), intent(in)      :: time   !< time value to get
    integer get_i
    integer upper,lower

    ! BC
    if (time <= ts%times(1)) then
       get_i = -1
       return
    end if
    if (time >= ts%times(ts%numt)) then
       get_i = ts%numt + 1
       return
    end if
    ! first try if the interpolated value is around the last value
    ts%current=min(ts%current,ts%numt-1)
    if (time >= ts%times(ts%current) .and. time < ts%times(ts%current+1)) then
       get_i = ts%current
       return
    end if
    ! this didn't work, let's try the next interval
    ts%current=ts%current+1
    if (time >= ts%times(ts%current) .and. time < ts%times(ts%current+1)) then
       get_i = ts%current
       return
    end if
    ! nope, let's do a Newton search
    lower = 1
    upper = ts%numt
    do
       ts%current = lower+int((upper-lower)/2.)
       if (time >= ts%times(ts%current) .and. time < ts%times(ts%current+1)) then
          get_i = ts%current
          return
       end if
       if (time > ts%times(ts%current)) then
          lower = ts%current
       else
          upper = ts%current
       end if
    end do
  end function get_i

end module glimmer_ts
