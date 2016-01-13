!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_external_dycore_interface.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_external_dycore_interface

contains

subroutine cism_init_external_dycore(external_dycore_type,model)
 
  use parallel
  use glimmer_global
  use glide
  use glissade
  use eismint_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init

  use glide_diagnostics

#if defined CISM_HAS_BISICLES || defined CISM_HAS_FELIX
#define CISM_HAS_EXTERNAL_DYCORE 1
#endif

#ifdef CISM_HAS_EXTERNAL_DYCORE
  use glimmer_to_dycore
#endif


  implicit none

  integer*4 :: external_dycore_type
  type(glide_global_type), intent(inout) :: model

  real(kind=dp) :: cur_time, time_inc

  ! for external dycore:
  integer*4 external_dycore_model_index
  ! integer argc
  integer*4 p_index


#ifdef CISM_HAS_EXTERNAL_DYCORE
  ! print *,"Initializing external dycore interface."
  call gtd_init_dycore_interface()

  call parallel_barrier()
  ! print *,"Initializing external dycore."
  call gtd_init_dycore(model,external_dycore_model_index)
  model%options%external_dycore_model_index = external_dycore_model_index
  call parallel_barrier()
#else
  print *,"ERROR: The program was not built with an external dynamic core."
#endif

end subroutine cism_init_external_dycore


subroutine cism_run_external_dycore(external_dycore_model_index,cur_time,time_inc)
  use parallel
  use glimmer_global
  use glide
  use glissade
  use eismint_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init

  use glide_diagnostics

#if defined CISM_HAS_BISICLES || defined CISM_HAS_FELIX
#define CISM_HAS_EXTERNAL_DYCORE 1
#endif

#ifdef CISM_HAS_EXTERNAL_DYCORE
  use glimmer_to_dycore
#endif

  integer*4 external_dycore_model_index
  real(kind=dp) :: cur_time, time_inc

#ifdef CISM_HAS_EXTERNAL_DYCORE
!  dycore_model_index = this_rank + 1
  external_dycore_model_index = 1

  call parallel_barrier()
  ! print *,"Running external dycore."
  call gtd_run_dycore(external_dycore_model_index,cur_time,time_inc)
  ! print *,"Completed Dycore Run."
  call parallel_barrier()
#else
  print *,"ERROR: The program was not built with an external dynamic core."
#endif

end subroutine cism_run_external_dycore

subroutine cism_finalize_external_dycore(external_dycore_type,model)
 
  use parallel
  use glimmer_global
  use glide
  use glissade
  use eismint_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init

  use glide_diagnostics

#if defined CISM_HAS_BISICLES || defined CISM_HAS_FELIX
#define CISM_HAS_EXTERNAL_DYCORE 1
#endif

#ifdef CISM_HAS_EXTERNAL_DYCORE
  use glimmer_to_dycore
#endif


  implicit none

  integer*4 :: external_dycore_type
  type(glide_global_type), intent(inout) :: model

  real(kind=dp) :: cur_time, time_inc

  ! for external dycore:
  integer*4 external_dycore_model_index
  ! integer argc
  integer*4 p_index


#ifdef CISM_HAS_EXTERNAL_DYCORE
   external_dycore_model_index = 1

  call parallel_barrier()
  ! print *,"Finalizing external dycore."
  call gtd_delete_dycore(external_dycore_model_index)
  model%options%external_dycore_model_index = external_dycore_model_index
  call parallel_barrier()
#else
  print *,"ERROR: The program was not built with an external dynamic core."
#endif

end subroutine cism_finalize_external_dycore


end module cism_external_dycore_interface
