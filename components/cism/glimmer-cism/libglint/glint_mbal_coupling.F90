!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_mbal_coupling.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_mbal_coupling

  use glint_mbal
  use glimmer_config

  implicit none

  !*FD Module to handle the accumulation of inputs and calculation of mass-balance

  type glint_mbc
     real(dp),dimension(:,:),pointer :: prcp_save  => null() !*FD used to accumulate precip
     real(dp),dimension(:,:),pointer :: ablt_save  => null() !*FD used to accumulate ablation
     real(dp),dimension(:,:),pointer :: acab_save  => null() !*FD used to accumulate mass-balance
     real(dp),dimension(:,:),pointer :: artm_save  => null() !*FD used to average air-temperature
     real(dp),dimension(:,:),pointer :: snowd      => null() !*FD Keeps track of snow depth
     real(dp),dimension(:,:),pointer :: siced      => null() !*FD Keeps track of superimposed ice depth 
     real(dp),dimension(:,:),pointer :: prcp       => null() !*FD Instantaneous precip
     real(dp),dimension(:,:),pointer :: ablt       => null() !*FD Instantaneous ablation
     real(dp),dimension(:,:),pointer :: acab       => null() !*FD Instantaneous mass-balance
     real(dp),dimension(:,:),pointer :: artm       => null() !*FD Instantaneous air temperature
     real(dp),dimension(:,:),pointer :: xwind      => null() !*FD Instantaneous x-wind
     real(dp),dimension(:,:),pointer :: ywind      => null() !*FD Instantaneous y-wind
     real(dp),dimension(:,:),pointer :: humidity   => null() !*FD Instantaneous humidity
     real(dp),dimension(:,:),pointer :: SWdown     => null() !*FD Instantaneous sw down
     real(dp),dimension(:,:),pointer :: LWdown     => null() !*FD Instantaneous lw down
     real(dp),dimension(:,:),pointer :: Psurf      => null() !*FD Instantaneous psurf
     real(dp),dimension(:,:),pointer :: snowd_save => null() !*FD Saves snow depth
     real(dp),dimension(:,:),pointer :: siced_save => null() !*FD Saves superimposed ice depth 
     integer :: av_count  = 0 !*FD Counter for averaging temperature input
     logical :: new_accum = .true.
     integer :: start_time    !*FD the time we started averaging (hours)
     type(glint_mbal_params) :: mbal
  end type glint_mbc

contains

  subroutine glint_mbc_init(params,lgrid,config,whichacab,snowd,siced,nx,ny,dx)

    use glimmer_coordinates

    type(glint_mbc)  :: params
    type(coordsystem_type) :: lgrid
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    integer          :: whichacab
    real(dp),dimension(:,:),intent(in) :: snowd !*FD Initial snow-depth field
    real(dp),dimension(:,:),intent(in) :: siced !*FD Initial superimposed ice field
    integer          :: nx,ny  !*FD grid dimensions (for SMB)
    real(dp)         :: dx     !* Grid length (for SMB)

    ! Deallocate if necessary

    if (associated(params%prcp_save))  deallocate(params%prcp_save)
    if (associated(params%ablt_save))  deallocate(params%ablt_save)
    if (associated(params%acab_save))  deallocate(params%acab_save)
    if (associated(params%artm_save))  deallocate(params%artm_save)
    if (associated(params%snowd))      deallocate(params%snowd)
    if (associated(params%siced))      deallocate(params%siced)
    if (associated(params%prcp))       deallocate(params%prcp)
    if (associated(params%ablt))       deallocate(params%ablt)
    if (associated(params%acab))       deallocate(params%acab)
    if (associated(params%artm))       deallocate(params%artm)
    if (associated(params%xwind))      deallocate(params%xwind)
    if (associated(params%ywind))      deallocate(params%ywind)
    if (associated(params%humidity))   deallocate(params%humidity)
    if (associated(params%SWdown))     deallocate(params%SWdown)
    if (associated(params%LWdown))     deallocate(params%LWdown)
    if (associated(params%Psurf))      deallocate(params%Psurf)
    if (associated(params%snowd_save)) deallocate(params%snowd_save)
    if (associated(params%siced_save)) deallocate(params%siced_save)

    ! Allocate arrays and zero

    call coordsystem_allocate(lgrid,params%prcp_save);  params%prcp_save = 0.d0
    call coordsystem_allocate(lgrid,params%ablt_save);  params%ablt_save = 0.d0
    call coordsystem_allocate(lgrid,params%acab_save);  params%acab_save = 0.d0
    call coordsystem_allocate(lgrid,params%artm_save);  params%artm_save = 0.d0
    call coordsystem_allocate(lgrid,params%snowd);      params%snowd = 0.d0
    call coordsystem_allocate(lgrid,params%siced);      params%siced = 0.d0
    call coordsystem_allocate(lgrid,params%prcp);       params%prcp = 0.d0
    call coordsystem_allocate(lgrid,params%ablt);       params%ablt = 0.d0
    call coordsystem_allocate(lgrid,params%acab);       params%acab = 0.d0
    call coordsystem_allocate(lgrid,params%artm);       params%artm = 0.d0
    call coordsystem_allocate(lgrid,params%xwind);      params%xwind = 0.d0
    call coordsystem_allocate(lgrid,params%ywind);      params%ywind = 0.d0
    call coordsystem_allocate(lgrid,params%humidity);   params%humidity = 0.d0
    call coordsystem_allocate(lgrid,params%SWdown);     params%SWdown = 0.d0
    call coordsystem_allocate(lgrid,params%LWdown);     params%LWdown = 0.d0
    call coordsystem_allocate(lgrid,params%Psurf);      params%Psurf = 0.d0
    call coordsystem_allocate(lgrid,params%snowd_save); params%snowd_save = 0.d0
    call coordsystem_allocate(lgrid,params%siced_save); params%siced_save = 0.d0

    ! Initialise the mass-balance scheme and other components

    call glint_mbal_init(params%mbal,config,whichacab,nx,ny,dx)

    ! Copy snow and ice depths if relevant

    if (mbal_has_snow_model(params%mbal)) then
       params%snowd = snowd
       params%siced = siced
    end if

  end subroutine glint_mbc_init

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate_mbal(params, time,     artm,    arng,    prcp,       &
                                   snowd,  siced,    xwind,   ywind,   local_orog, &
                                   thck,   humidity, SWdown,  LWdown,  Psurf)

    type(glint_mbc)  :: params
    integer :: time
    real(dp),dimension(:,:),intent(inout) :: artm      !*FD Mean air temperature (degC)
    real(dp),dimension(:,:),intent(in) :: arng         !*FD Air temperature half-range (degC)
    real(dp),dimension(:,:),intent(inout) :: prcp      !*FD Precipitation (m)
    real(dp),dimension(:,:),intent(in) :: snowd        !*FD Snow depth (m)
    real(dp),dimension(:,:),intent(in) :: siced        !*FD Superimposed ice depth (m)
    real(dp),dimension(:,:),intent(in) :: xwind        !*FD $x$-component of surface winds (m/s)
    real(dp),dimension(:,:),intent(in) :: ywind        !*FD $y$-component of surface winds (m/s)
    real(dp),dimension(:,:),intent(in) :: local_orog   !*FD Local orography (m)
    real(dp),dimension(:,:),intent(in) :: thck         !*FD Ice thickness (m)
    real(dp),dimension(:,:),intent(in) :: humidity     !*FD Relative humidity (%)
    real(dp),dimension(:,:),intent(in) :: SWdown       !*FD Downwelling shortwave (W/m^2)
    real(dp),dimension(:,:),intent(in) :: LWdown       !*FD Downwelling longwave (W/m^2)
    real(dp),dimension(:,:),intent(in) :: Psurf        !*FD Surface pressure (Pa)

    real(dp),dimension(size(artm,1),size(artm,2)) :: ablt,acab

    ! Things to do the first time

    if (params%new_accum) then

       params%new_accum = .false.
       params%av_count  = 0

       ! Initialise 

       params%snowd = snowd
       params%siced = siced

       params%prcp_save = 0.d0
       params%ablt_save = 0.d0
       params%acab_save = 0.d0
       params%artm_save = 0.d0

       params%start_time = time

    end if

    params%av_count = params%av_count+1

    call glint_mbal_calc(params%mbal,                      &
                         artm,  arng,  prcp,               &
                         (local_orog>0.d0 .or. thck>0.d0), &
                         params%snowd, params%siced,       &
                         ablt,  acab,  thck,               &
                         xwind, ywind, humidity,           &
                         SWdown,LWdown,Psurf) 

    ! Accumulate

    params%prcp_save = params%prcp_save + prcp
    params%ablt_save = params%ablt_save + ablt
    params%acab_save = params%acab_save + acab
    params%artm_save = params%artm_save + artm

    ! Copy instantaneous fields

    params%prcp = prcp
    params%ablt = ablt
    params%acab = acab
    params%artm = artm
    params%xwind = xwind
    params%ywind = ywind
    params%humidity = humidity
    params%SWdown = SWdown
    params%LWdown = LWdown
    params%Psurf = Psurf

  end subroutine glint_accumulate_mbal

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_average_mbal(params,artm,prcp,ablt,acab,snowd,siced,dt)

    use glint_constants, only: hours2years

    type(glint_mbc)  :: params
    real(dp),dimension(:,:),intent(out)   :: artm   !*FD Mean air temperature (degC)
    real(dp),dimension(:,:),intent(out)   :: prcp   !*FD Precipitation (m/yr)
    real(dp),dimension(:,:),intent(out)   :: ablt   !*FD Ablation (m/yr)
    real(dp),dimension(:,:),intent(out)   :: acab   !*FD Mass-balance (m/yr)
    real(dp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(dp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)
    integer,                intent(in)    :: dt     !*FD accumulation time in hours

    if (.not. params%new_accum) then
       params%artm_save = params%artm_save / real(params%av_count,dp)
    end if

    params%new_accum=.true.

    artm  = params%artm_save
    prcp  = params%prcp_save/real(dt*hours2years,dp)
    ablt  = params%ablt_save/real(dt*hours2years,dp)
    acab  = params%acab_save/real(dt*hours2years,dp)
    snowd = params%snowd
    siced = params%siced

    where (snowd < 0.d0) snowd = 0.d0
    where (siced < 0.d0) siced = 0.d0

  end subroutine glint_average_mbal

!++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_mbal_coupling

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
