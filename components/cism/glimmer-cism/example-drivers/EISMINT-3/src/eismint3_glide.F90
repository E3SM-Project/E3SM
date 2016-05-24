!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eismint3_glide.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program eismint3_glide

  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT 3 Greenland test case. Adapted from simple_glide.f90

!WHLTSTEP
!  use glimmer_global, only:rk,fname_length
  use glimmer_global, only: dp, fname_length
  use glide
  use glimmer_log
  use glimmer_config
  use eismint3_forcing
  use eismint3_io
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glide_nc_custom, only : glide_nc_fillall
  implicit none

  type(glide_global_type) :: model        ! model instance
  type(eismint3_climate) :: climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(dp) time
  real(dp) t1,t2
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
  call glide_initialise(model)
  call eismint3_initialise(climate,config,model)
  call eismint3_io_createall(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  do
     call eismint3_clim(climate,model)
     call glide_tstep_p1(model,time)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     call eismint3_io_writeall(climate,model)
     ! override masking stuff for now
     time = time + model%numerics%tinc
     if (time.gt.model%numerics%tend) exit
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

end program eismint3_glide
