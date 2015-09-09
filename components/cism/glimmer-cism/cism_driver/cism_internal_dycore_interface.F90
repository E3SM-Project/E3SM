!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_internal_dycore_interface.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module cism_internal_dycore_interface

contains

  subroutine cism_internal_dycore_interface(model)

 
  use parallel
  use glimmer_global, only:rk
  use glide
  use glissade
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_vector_ew, &
                                horiz_bcs_stag_vector_ns, horiz_bcs_stag_scalar

  use glide_diagnostics

  use glimmer_to_dycore

  implicit none

  type(glide_global_type), intent(inout) :: model

  real(kind=sp) cur_time, time_inc

  ! for external dycore:
  integer*4 dycore_model_index
  integer argc
  integer*4 p_index

  end subroutine cism_internal_dycore_interface

end module cism_internal_dycore_interface
