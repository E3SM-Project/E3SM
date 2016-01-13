!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_output_fluxes.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glad_output_fluxes

  !> This module defines a type and related operations for working with fluxes output to
  !> the GCM. Its main purpose is to produce temporal averages of these outputs.

  ! (Most of the code here used to be in glint_upscale.F90)

  use glimmer_global, only: dp

  implicit none
  private

  type, public :: glad_output_fluxes_type
     private

     integer  :: av_count_output       ! step counter
     
     real(dp), dimension(:,:), pointer :: hflx_sum => null()  ! conductive heat flux at top surface (W m-2)
     real(dp), dimension(:,:), pointer :: rofi_sum => null()  ! solid ice runoff (kg m-2 s-1)
     real(dp), dimension(:,:), pointer :: rofl_sum => null()  ! liquid runoff from basal/interior melting (kg m-2 s-1)
  end type glad_output_fluxes_type

  public :: initialize_glad_output_fluxes   ! Initialize a glad_output_fluxes instance
  public :: accumulate_output_fluxes        ! Accumulate one time step's contribution to fluxes
  public :: calculate_average_output_fluxes ! Compute and return time-average fluxes
  public :: reset_output_fluxes             ! Reset output_fluxes state to start a new averaging period
  
contains

  subroutine initialize_glad_output_fluxes(output_fluxes, ewn, nsn)
    ! Initialize a glad_output_fluxes instance

    type(glad_output_fluxes_type), intent(inout) :: output_fluxes

    ! dimensions of local grid
    integer, intent(in) :: ewn
    integer, intent(in) :: nsn

    allocate(output_fluxes%rofi_sum(ewn,nsn))
    allocate(output_fluxes%rofl_sum(ewn,nsn))
    allocate(output_fluxes%hflx_sum(ewn,nsn))

    call reset_output_fluxes(output_fluxes)

  end subroutine initialize_glad_output_fluxes

  subroutine accumulate_output_fluxes(output_fluxes, model)
    ! Given the calving, basal melting, and conductive heat flux fields from the dycore,
    ! accumulate contributions to the rofi, rofl, and hflx fields to be sent to the coupler.

    use glimmer_paramets, only: thk0, tim0
    use glimmer_physcon, only : rhoi
    use glide_types, only : glide_global_type

    type(glad_output_fluxes_type), intent(inout) :: output_fluxes
    type(glide_global_type), intent(in) :: model

    output_fluxes%av_count_output = output_fluxes%av_count_output + 1

    !--------------------------------------------------------------------
    ! Accumulate solid runoff (calving)
    !--------------------------------------------------------------------
                       
    ! Note on units: model%calving%calving_thck has dimensionless ice thickness units
    !                Multiply by thk0 to convert to meters of ice
    !                Multiply by rhoi to convert to kg/m^2 water equiv.
    !                Divide by (dt*tim0) to convert to kg/m^2/s

    ! Convert to kg/m^2/s
    output_fluxes%rofi_sum(:,:) = output_fluxes%rofi_sum(:,:)  &
         + model%calving%calving_thck(:,:) * thk0 * rhoi / (model%numerics%dt * tim0)

    !--------------------------------------------------------------------
    ! Accumulate liquid runoff (basal melting)
    ! Note: This is basal melting for grounded ice only.
    !       Basal melting for floating ice will typically be an input from the coupler, not an output.
    !--------------------------------------------------------------------
                       
    ! Note on units: model%temper%bmlt has dimensionless units of ice thickness per unit time
    !                Multiply by thk0/tim0 to convert to meters ice per second
    !                Multiply by rhoi to convert to kg/m^2/s water equiv.

    ! Convert to kg/m^2/s
    output_fluxes%rofl_sum(:,:) = output_fluxes%rofl_sum(:,:)  &
         + model%temper%bmlt_ground(:,:) * thk0/tim0 * rhoi

    !--------------------------------------------------------------------
    ! Accumulate basal heat flux
    !--------------------------------------------------------------------

    ! Note on units: model%temper%ucondflx has units of W/m^2, positive down
    !                Flip the sign so that hflx is positive up.

    output_fluxes%hflx_sum(:,:) = output_fluxes%hflx_sum(:,:) &
         - model%temper%ucondflx(:,:)

  end subroutine accumulate_output_fluxes

  subroutine calculate_average_output_fluxes(output_fluxes, rofi_tavg, rofl_tavg, hflx_tavg)
    ! Compute and return time-average fluxes

    type(glad_output_fluxes_type), intent(in) :: output_fluxes
    real(dp), dimension(:,:), intent(out) :: rofi_tavg  ! average solid ice runoff (kg m-2 s-1)
    real(dp), dimension(:,:), intent(out) :: rofl_tavg  ! average liquid runoff from basal/interior melting (kg m-2 s-1)
    real(dp), dimension(:,:), intent(out) :: hflx_tavg  ! average conductive heat flux at top surface (W m-2)
    
    if (output_fluxes%av_count_output > 0) then
       rofi_tavg(:,:) = output_fluxes%rofi_sum(:,:) / real(output_fluxes%av_count_output,dp)
       rofl_tavg(:,:) = output_fluxes%rofl_sum(:,:) / real(output_fluxes%av_count_output,dp)
       hflx_tavg(:,:) = output_fluxes%hflx_sum(:,:) / real(output_fluxes%av_count_output,dp)
    else
       rofi_tavg(:,:) = 0.d0
       rofl_tavg(:,:) = 0.d0
       hflx_tavg(:,:) = 0.d0
    end if

  end subroutine calculate_average_output_fluxes

  subroutine reset_output_fluxes(output_fluxes)
    ! Reset output_fluxes state to start a new averaging period

    type(glad_output_fluxes_type), intent(inout) :: output_fluxes

    output_fluxes%av_count_output = 0
    output_fluxes%rofi_sum(:,:) = 0.d0
    output_fluxes%rofl_sum(:,:) = 0.d0
    output_fluxes%hflx_sum(:,:) = 0.d0
  end subroutine reset_output_fluxes
    
end module glad_output_fluxes
