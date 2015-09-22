!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_writestats.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glimmer_writestats
  !> F90 wrapper to gc_writestats
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009

  implicit none

contains

  subroutine glimmer_write_stats(resname, cfgname,wallTime)
    use glimmer_global, only : dp
    implicit none
    character(len=*), intent(in) :: resname !< name of the output result file
    character(len=*), intent(in) :: cfgname !< name of configuration file
    real(kind=dp), intent(in)    :: wallTime!< elapsed wall clock tine in seconds

    call gf_writestats(resname,cfgname,wallTime)
  end subroutine glimmer_write_stats

end module glimmer_writestats
