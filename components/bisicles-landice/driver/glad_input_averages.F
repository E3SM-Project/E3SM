!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_input_averages.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glad_input_averages

  !> This module defines a type and related operations for working with inputs from the
  !> GCM. Its main purpose is to produce temporal averages of these inputs.

  ! Note that this module has some functionality in common with glad_mbal_coupling, but
  ! they are used at different stages in the time loops.

  ! NOTE(wjs, 2015-03-17) Most or all of the functionality here could be removed if we
  ! performed all of the necessary temporal averaging in the climate model, with coupling
  ! to the CISM code only happening once per mass balance time step. If we do that, we
  ! should probably add checks to ensure that the model is really just being called when
  ! it's time for a mass balance time step.

  use glimmer_global, only : dp
  use glimmer_paramets, only: GLC_DEBUG, stdout
  use glimmer_log
  use parallel, only : main_task
  
  implicit none
  private
  
  type, public :: glad_input_averages_type
     private

     integer :: av_start_time = 0  ! Value of time from the last occasion averaging was restarted (hours)
     integer :: av_steps      = 0  ! Number of times glimmer has been called in current round of averaging
     integer :: next_av_start = 0  ! Time when we expect next averaging to start
     logical :: new_av = .true.    ! Set to true if the next correct call starts a new averaging round
     
     real(dp),pointer,dimension(:,:) :: tot_qsmb => null()  ! running total surface mass balance (kg m-2 s-1)
     real(dp),pointer,dimension(:,:) :: tot_tsfc => null()  ! running total surface temperature (deg C)
  
  end type glad_input_averages_type

  public :: initialize_glad_input_averages
  public :: get_av_start_time
  public :: accumulate_averages
  public :: calculate_averages
  public :: reset_glad_input_averages
  
contains

  subroutine initialize_glad_input_averages(glad_inputs, ewn, nsn, next_av_start)
    ! Initialize a glad_inputs instance
    
    type(glad_input_averages_type), intent(inout) :: glad_inputs

    ! dimensions of local grid
    integer, intent(in) :: ewn
    integer, intent(in) :: nsn

    ! Starting time of next averaging period (hours)
    integer, intent(in) :: next_av_start

    allocate(glad_inputs%tot_qsmb(ewn,nsn));  glad_inputs%tot_qsmb = 0.d0
    allocate(glad_inputs%tot_tsfc(ewn,nsn));  glad_inputs%tot_tsfc = 0.d0

    glad_inputs%next_av_start = next_av_start
  end subroutine initialize_glad_input_averages

  integer function get_av_start_time(glad_inputs)
    ! Get value of time from the last occasion averaging was restarted (hours)
    type(glad_input_averages_type), intent(in) :: glad_inputs

    get_av_start_time = glad_inputs%av_start_time
  end function get_av_start_time
    
  subroutine accumulate_averages(glad_inputs, qsmb, tsfc, time)
    ! Accumulate averages based on one set of inputs.
    !
    ! Should be called every time we have new inputs from the climate model.
    
    type(glad_input_averages_type), intent(inout) :: glad_inputs
    real(dp),dimension(:,:),intent(in)  :: qsmb     ! flux of glacier ice (kg/m^2/s)
    real(dp),dimension(:,:),intent(in)  :: tsfc     ! surface ground temperature (C)
    integer, intent(in) :: time  ! Current model time
    
    if (glad_inputs%new_av) then
       call start_new_averaging_period(glad_inputs, time)
    end if

    glad_inputs%tot_qsmb(:,:) = glad_inputs%tot_qsmb(:,:) + qsmb(:,:)
    glad_inputs%tot_tsfc(:,:) = glad_inputs%tot_tsfc(:,:) + tsfc(:,:)
    
    glad_inputs%av_steps = glad_inputs%av_steps + 1
    
  end subroutine accumulate_averages

  subroutine calculate_averages(glad_inputs, qsmb, tsfc)
    ! Calculate averages over the averaging period
    type(glad_input_averages_type), intent(in) :: glad_inputs
    real(dp), dimension(:,:), intent(out) :: qsmb  ! average surface mass balance (kg m-2 s-1)
    real(dp), dimension(:,:), intent(out) :: tsfc  ! average surface temperature (deg C)
    
    qsmb(:,:) = glad_inputs%tot_qsmb(:,:) / real(glad_inputs%av_steps,dp)
    tsfc(:,:) = glad_inputs%tot_tsfc(:,:) / real(glad_inputs%av_steps,dp)
  end subroutine calculate_averages

  subroutine reset_glad_input_averages(glad_inputs, next_av_start)
    ! Resets this glad_inputs instance
    !
    ! Should be called at the end of an averaging period, in order to prepare for the
    ! next averaging period
    type(glad_input_averages_type), intent(inout) :: glad_inputs
    integer, intent(in) :: next_av_start  ! start time for next averaging period (hours)

    glad_inputs%tot_qsmb(:,:) = 0.d0
    glad_inputs%tot_tsfc(:,:) = 0.d0

    glad_inputs%av_steps      = 0
    glad_inputs%new_av        = .true.
    glad_inputs%next_av_start = next_av_start
  end subroutine reset_glad_input_averages
    
  subroutine start_new_averaging_period(glad_inputs, time)
    ! Should be called the first time accumulate_averages is called for a new averaging
    ! period. Sets some flags appropriately in this case.
    !
    ! Also performs some error checking to make sure we're not calling GLAD at an
    ! unexpected time.

    type(glad_input_averages_type), intent(inout) :: glad_inputs
    integer, intent(in) :: time  ! Current model time

    character(len=100) :: message
    
    if (GLC_DEBUG .and. main_task) then
       write (stdout,*) 'Accumulating averages, current time (hr) =', time
       write (stdout,*) 'av_start_time =', glad_inputs%av_start_time
       write (stdout,*) 'next_av_start =', glad_inputs%next_av_start
       write (stdout,*) 'new_av =', glad_inputs%new_av
    end if

    if (time == glad_inputs%next_av_start) then
       glad_inputs%av_start_time = time
       glad_inputs%new_av = .false.
    else
       write(message,*) 'Unexpected calling of GLAD at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if

  end subroutine start_new_averaging_period

end module glad_input_averages
