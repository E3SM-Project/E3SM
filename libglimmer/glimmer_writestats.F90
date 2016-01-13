!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_writestats.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glimmer_writestats
  !> F90 wrapper to gc_writestats
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009

  implicit none

contains

  subroutine glimmer_write_stats(resname, cfgname,wallTime)
    use glimmer_global, only : dp
    use parallel, only: main_task
    implicit none
    character(len=*), intent(in) :: resname !< name of the output result file
    character(len=*), intent(in) :: cfgname !< name of configuration file
    real(kind=dp), intent(in)    :: wallTime!< elapsed wall clock tine in seconds

    if (main_task) call gf_writestats(resname,cfgname,wallTime)
  end subroutine glimmer_write_stats

end module glimmer_writestats
