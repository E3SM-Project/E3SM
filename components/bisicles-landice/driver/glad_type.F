!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_type.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#define NCO outfile%nc
#define NCI infile%nc

module glad_type

  !> contains type definitions for GLAD

  use glimmer_global, only: dp
  ! TJL TODO - This needs to be fixed
  ! use glide_types
  ! use glad_input_averages, only : glad_input_averages_type, initialize_glad_input_averages
  ! use glad_mbal_coupling, only : glad_mbc
  ! use glad_output_fluxes, only : glad_output_fluxes_type, initialize_glad_output_fluxes
  
  implicit none

  ! Constants that describe the options available

  ! basic Glad options

  integer, parameter :: EVOLVE_ICE_FALSE = 0   ! do not let the ice sheet evolve
                                               ! (hold the ice state fixed at initial condition)
  integer, parameter :: EVOLVE_ICE_TRUE  = 1   ! let the ice sheet evolve

  integer, parameter :: ZERO_GCM_FLUXES_FALSE = 0 ! send true fluxes to the GCM
  integer, parameter :: ZERO_GCM_FLUXES_TRUE  = 1 ! zero out all fluxes sent to the GCM

  !TODO - Add other Glad options here to avoid hardwiring of case numbers?

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glad_instance

     !> Derived type holding information about ice model instance. 
     !> Note that variables used for downscaling & upscaling are only valid on the main task,
     !> since all downscaling and upscaling is done there.

     ! TJL TODO - This needs to be fixed
     ! type(coordsystem_type)           :: lgrid              !> Local grid for interfacing with glide (grid on this task)
     !                                                        !> (WJS: Note that origin may be incorrect with multiple tasks;
     !                                                        !>  as far as I can tell, this isn't currently a problem)
     ! type(glad_input_averages_type)   :: glad_inputs        !> Time-averaged inputs from the climate model
     ! type(glide_global_type)          :: model              !> The instance and all its arrays.
     ! character(fname_length)          :: paramfile          !> The name of the configuration file.
     integer                          :: ice_tstep          !> Ice timestep in hours
     integer                          :: mbal_tstep         !> Mass-balance timestep in hours
     integer                          :: mbal_accum_time    !> Accumulation time for mass-balance (hours)
                                                            !> (defaults to ice time-step)
     integer                          :: ice_tstep_multiply=1 !> Ice time multiplier (non-dimensional)
     integer                          :: n_icetstep         !> Number of ice time-steps per mass-balance accumulation
     real(dp)                         :: glide_time         !> Time as seen by glide (years)
     integer                          :: next_time          !> The next time we expect to be called (hours)

     ! Climate inputs, on the local grid -------------------------

     real(dp),dimension(:,:),pointer :: artm => null() !> Annual mean air temperature
     real(dp),dimension(:,:),pointer :: acab => null() !> Annual mass balance (m/y water equiv)

     ! Arrays to accumulate mass-balance quantities --------------

     ! TJL TODO - This needs to be fixed
     ! type(glad_mbc) :: mbal_accum
     
     ! Climate options -------------------------------------------

     integer :: evolve_ice = 1

     !> Whether the ice sheet can evolve:
     !> \begin{description}
     !> \item[0] The ice sheet cannot evolve; hold fixed at initial state
     !> \item[1] The ice sheet can evolve

     logical :: test_coupling = .false.

     integer :: zero_gcm_fluxes = ZERO_GCM_FLUXES_FALSE
     
     !> Whether to zero out the fluxes (e.g., calving flux) sent to the GCM
     !> \begin{description}
     !> \item[0] send true fluxes to the GCM
     !> \item[1] zero out all fluxes sent to the GCM
     !> \end{description}

     ! Latitude & longitude of model grid points
     real(dp), dimension(:,:), pointer :: lat(:,:) => null()
     real(dp), dimension(:,:), pointer :: lon(:,:) => null()

     ! TJL TODO - This needs to be fixed
     ! Fields for averaging dycore output
     ! type(glad_output_fluxes_type) :: glad_output_fluxes
     real(dp), dimension(:,:), pointer :: hflx_tavg => null()  ! conductive heat flux at top surface (W m-2)
     real(dp), dimension(:,:), pointer :: rofi_tavg => null()  ! solid ice runoff (kg m-2 s-1)
     real(dp), dimension(:,:), pointer :: rofl_tavg => null()  ! liquid runoff from basal/interior melting (kg m-2 s-1)
     
     ! Pointers to file input and output

     ! TJL TODO - This needs to be fixed
     ! type(glimmer_nc_output),pointer :: out_first => null() !> first element of linked list defining netCDF outputs
     ! type(glimmer_nc_input), pointer :: in_first => null()  !> first element of linked list defining netCDF inputs

  end type glad_instance
end module glad_type
