!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   simple_bisicles_old.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#define HAVE_MPI=1

program simple_bisicles
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats

  use glimmer_to_dycore

  implicit none


#ifdef GPTL
#include <gptl.inc>
#endif
#ifdef PAPI
#include <f90papi.h>
#endif

#ifdef HAVE_MPI
#include <mpif.h>

  integer npes,ierr,iam
#endif

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time

  ! make time vars single precision:
  real(kind=sp) cur_time, time_inc

  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

!  integer*4 dycore_type ! 0=BISICLES, 1=Octree
  integer*4 dycore_model_index
  integer argc

 ! print *,"Calling testGradient from Fortran..."
 ! call cmain(argc)

!  integer argv
!  print *,"Calling testLin from Fortran..."

!  call testLin(argc)

!  print *,"Calling testVelOpf from Fortran..."
!  call testvelopf(argc)

#ifdef HAVE_MPI
print *,"In HAVE_MPI, calling MPI_Init "
 call MPI_Init(ierr)

print *,"In HAVE_MPI, past MPI_Init"

 call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,iam,ierr)
print *,"In HAVE_MPI, mpes, iam = ",npes,iam

!  if (iam == 0) then

#endif

print *,"Running simple_bisicles, glimmer-cism-lanl version."


  ! start gptl
#ifdef GPTL
  ret = gptlsetoption (gptlprint_method,gptlfull_tree)
  ret = gptlsetoption (PAPI_FP_OPS, 1)
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('total')
#endif




  call glimmer_GetCommandline()

  print *,"config_fname: ",commandline_configname
  print *,"results_fname: ",commandline_resultsname

!  print *,"Calling linear solver from Fortran..."   
!  call test_linear_solver_call()  
!  print *,"Returned from linear solver call"   


!  print *,"Calling test_vel_op from Fortran..."   
!  call test_vel_op()  
!  print *,"Returned from test_vel_op call"   

  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! read configuration
!print *,"Calling ConfigRead..."
  call ConfigRead(commandline_configname,config)
print *,"Returned from ConfigRead call"  
  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise GLIDE
  call glide_config(model,config)

  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart
  time_inc = model%numerics%tinc

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

print *,"In simple_bisicles, calling init_dycore_interface..."  
  cur_time = time
  call gtd_init_dycore_interface()

print *,"In simple_bisicles, topg shape = ",size(shape(model%geometry%topg)), &
        shape(model%geometry%topg)

print *,"In simple_bisicles, btrc shape = ",size(shape(model%velocity%btrc)), &
        shape(model%velocity%btrc)  

!print *,model%velocity%btrc(1:4,1:4)

  call gtd_init_dycore(model,dycore_model_index)

  do while(time.le.model%numerics%tend)
     print *,"In simple_bisicles while loop, time = ",time

!     call glide_tstep_p1(model,time)
!     call glide_tstep_p2(model)
    
     ! run the BISICLES dynamic core model in place of tsteps p1 and p2:
     cur_time = time
     call gtd_run_dycore(dycore_model_index,cur_time,time_inc)

     call glide_tstep_p3(model)

     ! override masking stuff for now
     time = time + model%numerics%tinc
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)
  end do

  ! finalise GLIDE
  call glide_finalise(model)

  print *,"Calling gtd_delete_dycore..."
  call gtd_delete_dycore(dycore_model_index)
  print *,"Completed gtd_delete_dycore."

  print *,"Calling gtd_delete_dycore_interface..."
  call gtd_delete_dycore_interface();
  print *,"Completed gtd_delete_dycore_interface."

  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
!  call glimmer_writestats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log


  ! stop gptl
#ifdef GPTL
  ret = gptlstop ('total')
  ret = gptlpr (0)
  ret = gptlfinalize ()
#endif

#ifdef HAVE_MPI

!  endif   ! iam == 0 endif

!  call MPI_Barrier() 
print *,"In HAVE_MPI, calling MPI_Finalize "
   call MPI_Finalize(ierr)
#endif
 
end program simple_bisicles
