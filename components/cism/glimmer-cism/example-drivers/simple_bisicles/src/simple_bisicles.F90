!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   simple_bisicles.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program simple_bisicles
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases

!  use parallel, only: nsub, ewub, nslb, ewlb, global_nsn, global_ewn, own_ewn, own_nsn, lhalo, uhalo, &
!                        rank => this_rank, ewtasks => ProcsEW, tasks, comm
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
  use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
                               horiz_bcs_stag_scalar

  use glide_diagnostics

  use glimmer_to_dycore

  use glide_nc_custom

  implicit none


  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  real(kind=sp) cur_time, time_inc

  integer :: tstep_count

  ! for external dycore:
  integer*4 dycore_model_index
  integer argc
  integer*4 p_index

  call parallel_initialise

  call glimmer_GetCommandline()

  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))

  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  call ConfigRead(commandline_configname,config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise GLIDE
  call glide_config(model,config)
!  call simple_initialise(climate,config)

  ! replace glide_initialise with glissade_initialise,
  ! for consistent halos in serial mode:
  call glissade_initialise(model)

  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart
  time_inc = model%numerics%tinc

!  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

  dycore_model_index = this_rank + 1
  
  print *,"Initializing external dycore interface."
  call gtd_init_dycore_interface()

!  print*, "acab before dycore:"
!  print*, model%climate%acab

  call parallel_barrier()
  print *,"Initializing external dycore."
  call gtd_init_dycore(model,dycore_model_index)
  call parallel_barrier()

!  print*, "acab after dycore:"
!  print*, model%climate%acab


  
! write nc files
  call glide_io_writeall(model, model, time=time)         

  print *,"Running external dycore."
  do while (cur_time < model%numerics%tend)
      print *, "CISM timestep: time = ", cur_time, ", dt = ", time_inc
      call gtd_run_dycore(dycore_model_index,cur_time,time_inc)
!     write nc files                                           
      call glide_io_writeall(model, model, time=time)          
  end do
  print *,"Completed Dycore Run."
  call parallel_barrier()

  ! finalise 
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

  call parallel_finalise
  
end program simple_bisicles
