!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_mbal_coupling.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glad_mbal_coupling

  use glimmer_config
  use glimmer_global, only : dp

  implicit none

  ! Module to handle the accumulation of inputs.

  ! Note that this module has some functionality in common with glad_input_averages, but
  ! they are used at different stages in the time loops.

  type glad_mbc
     real(dp),dimension(:,:),pointer :: acab_save  => null() ! used to accumulate mass-balance
     real(dp),dimension(:,:),pointer :: artm_save  => null() ! used to average air-temperature
     real(dp),dimension(:,:),pointer :: acab       => null() ! Instantaneous mass-balance
     real(dp),dimension(:,:),pointer :: artm       => null() ! Instantaneous air temperature
     integer :: av_count  = 0 ! Counter for averaging inputs
     logical :: new_accum = .true.
     integer :: start_time    ! the time we started averaging (hours)
     integer :: tstep ! Timestep of mass-balance scheme in hours
  end type glad_mbc

contains

  subroutine glad_mbc_init(params,lgrid)

    ! Initialize the glad_mbc structure ('params').

    ! NOTE(wjs, 2015-03-19) In glint, when using SMB coupling, this was done in
    ! glint_downscale.F90: glint_init_input_gcm (rather than in glint_mbc_init). However,
    ! for simplicity and modularity, I am moving operations like this that act on glad_mbc
    ! into this glad_mbal_coupling module.

    use glimmer_coordinates
    use glad_constants, only : years2hours

    type(glad_mbc)  :: params
    type(coordsystem_type) :: lgrid

    ! Deallocate if necessary

    if (associated(params%acab_save))  deallocate(params%acab_save)
    if (associated(params%artm_save))  deallocate(params%artm_save)
    if (associated(params%acab))       deallocate(params%acab)
    if (associated(params%artm))       deallocate(params%artm)

    ! Allocate arrays and zero

    call coordsystem_allocate(lgrid,params%acab_save);  params%acab_save = 0.d0
    call coordsystem_allocate(lgrid,params%artm_save);  params%artm_save = 0.d0
    call coordsystem_allocate(lgrid,params%acab);       params%acab = 0.d0
    call coordsystem_allocate(lgrid,params%artm);       params%artm = 0.d0

    ! Set default mass balance time step
    !
    ! This is the default value that was being used in glint for the MASS_BALANCE_GCM
    ! scheme (some other schemes used different defaults)
    params%tstep = nint(years2hours)   ! mbal tstep = 1 year
    
  end subroutine glad_mbc_init

!++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_accumulate_input_gcm(params, time, acab, artm)

    ! In glint, this was done in glint_downscale.F90
    
    type(glad_mbc)  :: params
    integer :: time

    real(dp),dimension(:,:),intent(in) :: acab   ! Surface mass balance (m)
    real(dp),dimension(:,:),intent(in) :: artm   ! Mean air temperature (degC)

    ! Things to do the first time

    if (params%new_accum) then

       params%new_accum = .false.
       params%av_count  = 0

       ! Initialise

       params%acab_save = 0.d0
       params%artm_save = 0.d0
       params%start_time = time

    end if

    params%av_count = params%av_count + 1

    ! Accumulate

    params%acab_save = params%acab_save + acab
    params%artm_save = params%artm_save + artm

    ! Copy instantaneous fields

    params%acab = acab
    params%artm = artm

  end subroutine glad_accumulate_input_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_average_input_gcm(params, dt, acab, artm)

    ! In glint, this was done in glint_downscale.F90

    use glad_constants, only: hours2years

    type(glad_mbc)  :: params
    integer,                intent(in)    :: dt     !> mbal accumulation time (hours)
    real(dp),dimension(:,:),intent(out)   :: artm   !> Mean air temperature (degC)
    real(dp),dimension(:,:),intent(out)   :: acab   !> Mass-balance (m/yr)

    if (.not. params%new_accum) then
       params%artm_save = params%artm_save / real(params%av_count,dp)
    end if
    artm  = params%artm_save

    ! Note: acab_save has units of m, but acab has units of m/yr
    acab  = params%acab_save / real(dt*hours2years,dp)

    params%new_accum = .true.

  end subroutine glad_average_input_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++

  
end module glad_mbal_coupling

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
