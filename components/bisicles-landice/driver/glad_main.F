!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_main.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glad_main

  ! This module provides an interface to GCMs in the case where fields have already been
  ! downscaled to the ice sheet grid (and the GCM does its own upscaling from the ice
  ! sheet grid to the land grid).
  !
  ! This only provides code for the SMB case, not for the PDD case.

  use glimmer_global, only: dp, fname_length
  use glad_type
  ! TJL TODO - This needs to be fixed
  ! use glad_constants
  ! use glimmer_config
  ! use glimmer_filenames, only : process_path
  ! use parallel, only: main_task
  ! use glad_input_averages, only : get_av_start_time, accumulate_averages, &
  !      calculate_averages, reset_glad_input_averages
  
  ! use glimmer_paramets, only: stdout, GLC_DEBUG

  implicit none
  private
  
  ! ------------------------------------------------------------
  ! glad_params derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type, public :: glad_params 

     !> Derived type containing parameters relevant to all instances of 
     !> the model - i.e. those parameters which pertain to the global model. 

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !> Number of ice model instances
     character(fname_length),pointer,dimension(:) :: config_fnames => null()    ! array of config filenames
     type(glad_instance),pointer,dimension(:) :: instances  => null() !> Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     integer  :: tstep_mbal = 1        !> Mass-balance timestep (hours)
     integer  :: start_time            !> Time of first call to glad (hours)
     integer  :: time_step             !> Calling timestep of global model (hours)

     ! Parameters that can be set by the GCM calling Glad

     logical  :: gcm_restart = .false. !> If true, restart the model from a GCM restart file
     character(fname_length) :: gcm_restart_file   !> Name of restart file
     integer  :: gcm_fileunit = 99     !> Fileunit specified by GCM for reading config files

  end type glad_params

end module glad_main

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
