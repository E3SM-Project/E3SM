!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   gcm_cism_interface.F90 - part of the Community Ice Sheet Model (CISM)  
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

! from glide_types.F90:
!  integer, parameter :: DYCORE_GLIDE = 0     ! old shallow-ice dycore from Glimmer
!  integer, parameter :: DYCORE_GLAM = 1      ! Payne-Price finite-difference solver
!  integer, parameter :: DYCORE_GLISSADE = 2  ! prototype finite-element solver
!  integer, parameter :: DYCORE_ALBANYFELIX = 3  ! External Albany-Felix finite-element solver
!  integer, parameter :: DYCORE_BISICLES = 4     ! BISICLES external dycore

module gcm_cism_interface

  use parallel
  use glint_commandline
  use glide
  use cism_front_end

  use glint_example_clim
  use glint_main
  use gcm_to_cism_glint


  integer, parameter :: GCM_MINIMAL_MODEL = 0
  integer, parameter :: GCM_DATA_MODEL = 1
  integer, parameter :: GCM_CESM = 2

contains

subroutine gci_init_interface(which_gcm,g2c)
  use parallel
  use glint_commandline
  use glimmer_config
  use glide
  use glide_types
 
  use cism_front_end 

  integer, intent(in) :: which_gcm
  type(gcm_to_cism_type) :: g2c   ! holds everything

  integer :: whichdycore
  type(ConfigSection), pointer :: config  ! configuration stuff
  type(ConfigSection), pointer :: section  !< pointer to the section to be checked

  ! call parallel_initialise
  
  ! get the CISM dycore to be used:
  call glint_GetCommandline()
  call open_log(unit=50, fname=logname(commandline_configname))
  call ConfigRead(commandline_configname,config)
  call GetSection(config,section,'options')
  call GetValue(section,'dycore',whichdycore)
  if (main_task) print *,'CISM dycore type (0=Glide, 1=Glam, 2=Glissade, 3=AlbanyFelix, 4 = BISICLES) = ', whichdycore  

  ! check to see if running minimal GCM or data GCM.  Still need to add CESM GCM:
  call GetSection(config,section,'GLINT climate')

  if (associated(section)) then
    g2c%which_gcm = GCM_DATA_MODEL
  else 
    g2c%which_gcm = GCM_MINIMAL_MODEL
  end if
  if (main_task) print *,'g2c%which_gcm (1 = data, 2 = minimal) = ',g2c%which_gcm

  select case (g2c%which_gcm)
    case (GCM_MINIMAL_MODEL)
      if (main_task) print*, 'call cism_init_dycore'
      call cism_init_dycore(g2c%glide_model)
 
    case (GCM_DATA_MODEL)
      if (main_task) print*, 'call g2c_glint_init'
      call g2c_glint_init(g2c)

    case (GCM_CESM)
      ! call gcm_glint_GetCommandline_proxy()
      ! call g2c_glint_init(g2c) 

    case default
      if (main_task) print *,"Error -- unknown GCM type."
  end select

end subroutine gci_init_interface   

subroutine gci_run_model(g2c)

  type(gcm_to_cism_type) :: g2c 

  logical :: finished = .false.

  do while (.not. finished)
    select case (g2c%which_gcm)
      case (GCM_MINIMAL_MODEL)
        ! call gcm_update_model(gcm_model,cism_model)
!        if (main_task) print *,"In gci_run_model, calling cism_run_dycore"
        call cism_run_dycore(g2c%glide_model)

      case (GCM_DATA_MODEL,GCM_CESM)
!        if (main_task) print *,"In gci_run_model, calling g2c_glint_run"
        call g2c_glint_run(g2c)
        call g2c_glint_climate_time_step(g2c)
      case default
    end select
    finished = (gci_finished(g2c))
  end do
end subroutine gci_run_model


! gci_finished is used to test status of GCM
function gci_finished(g2c) result(finished)

  type(gcm_to_cism_type) :: g2c
  logical :: finished
 
  select case (g2c%which_gcm)
    case (GCM_MINIMAL_MODEL)
      finished = .true.

    case (GCM_DATA_MODEL,GCM_CESM)
      call g2c_glint_check_finished(g2c,finished)
    case default
  end select
  !if (main_task) print *,"In gci_finished, finished = ",finished  

end function gci_finished


subroutine gci_finalize_interface(g2c)

  type(gcm_to_cism_type) :: g2c

  select case (g2c%which_gcm)
    case (GCM_MINIMAL_MODEL)
      call cism_finalize_dycore(g2c%glide_model)
 
    case (GCM_DATA_MODEL)
      call g2c_glint_end(g2c)

    case (GCM_CESM)
      ! call g2c_glint_end(g2c)
    case default
  end select

end subroutine gci_finalize_interface


end module gcm_cism_interface
