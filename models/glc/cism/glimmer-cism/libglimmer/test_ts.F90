!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_ts.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Move this program to another directory so that this directory contains no programs?

program test_ts

  !*FD testing the time series module
  use glimmer_ts
  use glimmer_global, only: dp
  implicit none

  character(len=50) :: fname
  integer numv
  real(dp), dimension(:), allocatable :: val1,val2
  real(dp) :: time
  type(glimmer_tseries) :: ts

  write(*,*) 'Enter name containing time series'
  read(*,*) fname
  write(*,*) 'Enter number of values per time'
  read(*,*) numv
  
  allocate(val1(numv))
  allocate(val2(numv))
  
  call glimmer_read_ts(ts,fname,numv)
  
  do
     write(*,*) 'Enter new time'
     read(*,*) time
     call glimmer_ts_step(ts,time,val1)
     call glimmer_ts_linear(ts,time,val2)
     write(*,*) val1,val2
  end do
end program test_ts
