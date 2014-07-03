!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_config.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program testconfig
  !*FD testing config module
  !*FD written by Magnus Hagdorn, May 2004
  !*FD extended by Ian Rutt, August 2004
  use glimmer_config
  implicit none

  character(len=30) fname,sname,vname,vtype
  type(ConfigSection), pointer :: config,section

  character(len=100) :: charval
  integer :: intval
  real :: realval
  real, dimension(:), pointer :: realarray

  write(*,*) 'Enter name of configuration file'
  read(*,*) fname

  call ConfigRead(fname,config)

  call PrintConfig(config)

  do

     write(*,*)'Choose a section:'
     read(*,'(A)') sname
     if (trim(sname)=='*') exit
     call GetSection(config,section,sname)
     if (.not.associated(section)) then
        write(*,*) 'No such section'
        cycle
     end if

     do

        write(*,*)'Choose a value name:'
        read(*,*) vname
        if (trim(vname)=="*") exit
        write(*,*) 'Real(r), Integer(i), character(c) or real array (a)?'
        read(*,*) vtype
     
        select case(trim(vtype))
        case('i')
           call GetValue(section,trim(vname),intval)
           write(*,*)intval
        case('r')
           call GetValue(section,trim(vname),realval)
           write(*,*)realval
        case('c')
           call GetValue(section,trim(vname),charval)
           write(*,*)charval
        case('a')
           call GetValue(section,trim(vname),realarray)
           write(*,*) realarray
        end select
     end do
  end do

end program testconfig
