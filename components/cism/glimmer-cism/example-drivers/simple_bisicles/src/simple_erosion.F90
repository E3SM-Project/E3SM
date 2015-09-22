!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   simple_erosion.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program simple_erosion
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use erosion
  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(erosion_type) :: er           ! erosion
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate

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
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call er_initialise(er,config,model)
  call CheckSections(config)

  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart
  do while(time.le.model%numerics%tend)
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     
     call glide_tstep_p1(model,time)
     call er_tstep(er,model)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now
     time = time + model%numerics%tinc
  end do

  ! finalise GLIDE
  call er_finalise(er)
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log
end program simple_erosion
