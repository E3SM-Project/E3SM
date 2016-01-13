!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_driver.F90 - part of the Community Ice Sheet Model (CISM)  
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
program cism_driver

  use parallel
!  use glimmer_commandline
!  use glide
  use gcm_cism_interface
  use parallel

  integer :: which_gcm = GCM_DATA_MODEL
  type(gcm_to_cism_type) :: g2c

  if (command_argument_count() == 0) then
     print *,""
     print *,"Call cism_driver with either 1 or 2 arguments. Examples:"
     print *,"cism_driver ice_sheet.config"
     print *,"cism_driver ice_sheet.config climate.config"
     print *,""
     stop
  end if

  call parallel_initialise

  call gci_init_interface(which_gcm,g2c)
  call gci_run_model(g2c)
  call gci_finalize_interface(g2c)

  call parallel_finalise
end program cism_driver
