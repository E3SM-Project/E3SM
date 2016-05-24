!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_daily_pdd.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_daily_pdd

  !*FD The daily PDD model.
  !*FD N.B. all quantities, inputs and outputs are in m water equivalent.
  !*FD PDD factors are no longer converted to ice equivalent.

  use glimmer_global, only: dp 
  use glimmer_physcon, only: pi,scyr,rhow,rhoi

  implicit none

  type glint_daily_pdd_params

     !*FD Holds parameters for daily positive-degree-day mass-balance
     !*FD calculation. 

     !TODO - Change to dp
     real(dp) :: wmax        = 0.6d0       !*FD Fraction of firn that must be ice before run-off occurs
     real(dp) :: pddfac_ice  = 0.008d0     !*FD PDD factor for ice (m water day$^{-1}$ $^{\circ}C$^{-1}$)
     real(dp) :: pddfac_snow = 0.003d0     !*FD PDD factor for snow (m water day$^{-1}$ $^{\circ}C$^{-1}$)
     real(dp) :: rain_threshold = 1.d0     !*FD Threshold for precip melting (degC)
     integer  :: whichrain = 1             !*FD method for determining whether precip is rain or snow
     real(dp) :: tau0 = 10.d0*scyr         !*FD Snow densification timescale (seconds)
     real(dp) :: constC = 0.0165d0         !*FD Snow density profile factor C (m$^{-1}$)
     real(dp) :: firnbound = 0.872d0       !*FD Ice-firn boundary as fraction of density of ice
     real(dp) :: snowdensity = 300.d0      !*FD Density of fresh snow ($\mathrm{kg}\,\mathrm{m}^{-3}$)
     real(dp) :: tstep = 24.d0*60.d0*60.d0 !*FD Scheme time-step (seconds)
     real(dp) :: a1,a2,a3                  !*FD Factors for relaxation of depth

  end type glint_daily_pdd_params

  real(dp),parameter :: one_over_pi=1.d0/pi

  private
  public :: glint_daily_pdd_params, glint_daily_pdd_init, glint_daily_pdd_mbal

contains

  subroutine glint_daily_pdd_init(params,config)

    use glimmer_physcon, only : rhow, rhoi
    use glimmer_config
  
    type(glint_daily_pdd_params) :: params
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file

    ! Read the config file and output to log

    call daily_pdd_readconfig(params,config)
    call daily_pdd_printconfig(params)

    params%a1 = params%tstep/params%tau0
    params%a2 = 1.d0 - params%a1/2.d0
    params%a3 = 1.d0 + params%a1/2.d0

  end subroutine glint_daily_pdd_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine daily_pdd_readconfig(params,config)

    !*FD Reads in configuration data for the daily PDD scheme.

    use glimmer_config

    type(glint_daily_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    real(dp) :: tau0

    tau0 = params%tau0/scyr

    call GetSection(config,section,'GLINT daily pdd')
    if (associated(section)) then
       call GetValue(section,'wmax',params%wmax)
       call GetValue(section,'pddfac_ice',params%pddfac_ice)
       call GetValue(section,'pddfac_snow',params%pddfac_snow)
       call GetValue(section,'rain_threshold',params%rain_threshold)
       call GetValue(section,'whichrain',params%whichrain)
       call GetValue(section,'tau0',tau0)
       call GetValue(section,'constC',params%constC)
       call GetValue(section,'firnbound',params%firnbound)
       call GetValue(section,'snowdensity',params%snowdensity)
    end if

    params%tau0 = tau0*scyr

  end subroutine daily_pdd_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine daily_pdd_printconfig(params)

    use glimmer_log

    type(glint_daily_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    character(len=100) :: message

    call write_log_div

    call write_log('GLINT daily PDD Scheme parameters:')
    call write_log('-----------------------------------')
    write(message,*) 'Snow refreezing fraction',params%wmax
    call write_log(message)
    write(message,*) 'PDD factor for ice',params%pddfac_ice
    call write_log(message)
    write(message,*) 'PDD factor for snow',params%pddfac_snow
    call write_log(message)
    write(message,*) 'Rain threshold temperature',params%rain_threshold,' degC'
    call write_log(message)
    write(message,*) 'Rain/snow partition method',params%whichrain
    call write_log(message)
    write(message,*) 'Snow densification time-scale',params%tau0/scyr,' years'
    call write_log(message)
    write(message,*) 'Snow density equilibrium profile factor',params%constC,' m^-1'
    call write_log(message)
    write(message,*) 'Ice-firn boundary as fraction of density of ice',params%firnbound
    call write_log(message)
    write(message,*) 'Density of fresh snow',params%snowdensity,'kg m^-3'
    call write_log(message)
    call write_log('')

  end subroutine daily_pdd_printconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_daily_pdd_mbal(params,artm,arng,prcp,snowd,siced,ablt,acab,landsea)

    type(glint_daily_pdd_params)          :: params !*FD Daily PDD scheme parameters
    real(dp),dimension(:,:),intent(in)    :: artm   !*FD Daily mean air-temperature ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)    :: arng   !*FD Daily temperature half-range ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)    :: prcp   !*FD Daily precipitation (m)
    real(dp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(dp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)
    real(dp),dimension(:,:),intent(out)   :: ablt   !*FD Daily ablation (m)
    real(dp),dimension(:,:),intent(out)   :: acab   !*FD Daily mass-balance (m)
    logical, dimension(:,:),intent(in)    :: landsea !*FD Land-sea mask (land is TRUE)

    real(dp),dimension(size(prcp,1),size(prcp,2)) :: rain ! Daily rain
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: degd ! Degree-day
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: giced ! Temporary array for glacial ice
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: old_snow,old_sice

    integer :: nx,ny,i,j

    nx=size(prcp,1) ; ny=size(prcp,2)

    old_snow = snowd
    old_sice = siced

    rain = rainorsnw(params%whichrain,   &
                     artm,  arng,  prcp, &
                     params%rain_threshold)
    degd = finddegdays(artm,arng)
    giced = 0.d0

    do i=1,nx
       do j=1,ny
          if (landsea(i,j)) then
             call degdaymodel(params,snowd(i,j),siced(i,j),giced(i,j),degd(i,j),rain(i,j),prcp(i,j)) 
             acab(i,j) = snowd(i,j)+siced(i,j)+giced(i,j)-old_snow(i,j)-old_sice(i,j)
             ablt(i,j) = max(prcp(i,j)-acab(i,j), 0.d0)
             call firn_densify(params,snowd(i,j),siced(i,j))
          else
             ablt(i,j)  = prcp(i,j)
             acab(i,j)  = 0.d0
             snowd(i,j) = 0.d0
             siced(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine glint_daily_pdd_mbal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(dp) function finddegdays(localartm,localarng)

    !*FD Finds the degree-day as a function of mean daily temperature and
    !*FD half-range. The calculation is made on the assumption that
    !*FD daily temperatures vary as $T = T_{0} - \Delta T * \cos(\theta)$.

    real(dp), intent(in) :: localartm !*FD Mean daily temperature (degC)
    real(dp), intent(in) :: localarng !*FD Daily temperture half-range (degC)

    real(dp) :: time

    if (localartm - localarng > 0.d0) then
       finddegdays = localartm
    else if (localartm + localarng < 0.d0) then  
       finddegdays = 0.d0
    else
       time = acos(localartm / localarng) 
       finddegdays = (localartm*(pi-time) + localarng*sin(time)) * one_over_pi
    end if

  end function finddegdays

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(dp) function rainorsnw(which,localartm,localarng,localprcp,threshold)

    !*FD Determines a value for snow precipitation, dependent on air temperature and half-range.
    !*FD Takes precipitation as input and returns amount of rain.

    integer, intent(in) :: which     !*FD Selects method of calculation
    real(dp),intent(in) :: localartm !*FD Air temperature (degC)
    real(dp),intent(in) :: localarng !*FD Temperature half-range (degC)
    real(dp),intent(in) :: localprcp !*FD Precipitation (arbitrary units)
    real(dp),intent(in) :: threshold !*FD Snow/rain threshold (degC)

    real(dp) :: acosarg

    select case(which)

    case(1)

       if (localarng == 0.d0) then
          ! There is no sinusoidal variation
          if (localartm > threshold) then
             rainorsnw = localprcp
          else
             rainorsnw = 0.d0
          end if
       else
          ! Assume sinusoidal variation
          acosarg = (localartm - threshold) / localarng

          if (acosarg <= -1.d0) then
             rainorsnw = 0.d0
          else if (acosarg >= 1.d0) then
             rainorsnw = localprcp
          else
             rainorsnw = localprcp * (1.d0 - one_over_pi * acos(acosarg))
          end if

       end if

    case default

       ! Just use mean temperature to determine if precip is snow or rain

       if (localartm > threshold) then
          rainorsnw = localprcp
       else
          rainorsnw = 0.d0
       end if

    end select

  end function rainorsnw

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine degdaymodel(params,snowdepth,sicedepth,gicedepth,degd,rain,prcp)

    !*FD Applies the positive degree-day model.
    !*FD The output of this subroutine is in the variables \texttt{snowdepth},
    !*FD \texttt{sicedepth}, and \texttt{gicedepth}, which give the new depths
    !*FD of snow, superimposed ice and glacial ice, respectively.

    type(glint_daily_pdd_params) :: params !*FD PDD parameters
    real(dp), intent(inout) :: snowdepth     !*FD Snow depth (m)
    real(dp), intent(inout) :: sicedepth     !*FD Superimposed ice depth (m)
    real(dp), intent(inout) :: gicedepth     !*FD Glacial ice depth (m)
    real(dp), intent(in)    :: degd          !*FD The degree-day (degC day)
    real(dp), intent(in)    :: prcp          !*FD Total precip (m)
    real(dp), intent(in)    :: rain          !*FD Rain (m)

    real(dp) :: potablt, wfrac, rfrez

    !-------------------------------------------------------------------------
    ! Assume snowfall goes into snow, and rainfall goes into superimposed ice

    snowdepth = snowdepth + prcp - rain
    sicedepth = sicedepth + rain

    ! Calculate amount of superimposed ice needed before runoff can occur:
    ! a fraction of the total depth of the firn layer.

    wfrac = params%wmax * (snowdepth + sicedepth)

    ! Initial potential ablation of snow

    potablt = degd * params%pddfac_snow
    
    ! Calculate possible amount of refreezing. We need to do this to 
    ! take into account of the fact that there may already be more superimposed
    ! ice than is required for runoff

    rfrez  = max(0.d0, wfrac-sicedepth)

    ! Start off trying to ablate snow, and add it to superimposed ice

    if (potablt < snowdepth) then 

       ! If we have enough snow to consume all potential ablation
       ! Melt that amount of snow
       snowdepth = snowdepth - potablt

       ! Refreeze up to the limit
       sicedepth = sicedepth + min(potablt,rfrez)

       ! We've used all the potential ablation, so set to zero
       potablt = 0.d0

    else

       ! If there isn't enough snow to use all the potential ablation
       ! Set potential ablation to what remains - we use this later.
       ! For this section, the ablation is the whole snow depth.
       potablt=potablt-snowdepth

       ! Refreeze up to the limit
       sicedepth = sicedepth+min(snowdepth,rfrez)

       ! We've ablated all the snow, so set to zero
       snowdepth = 0.d0

    endif

    ! If we have any potential ablation left, use it to melt ice, 
    ! first from the firn, and then glacial ice

    if (potablt > 0.d0) then 
       potablt = params%pddfac_ice * (potablt - snowdepth) / params%pddfac_snow       
       if (potablt < sicedepth) then
          sicedepth = sicedepth-potablt
        else
          potablt = potablt-sicedepth
          sicedepth = 0.d0
          gicedepth = gicedepth-potablt
       end if
    end if

  end subroutine degdaymodel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine firn_densify(params,snowdepth,sicedepth)

    type(glint_daily_pdd_params) :: params   !*FD PDD parameters
    real(dp), intent(inout) :: snowdepth     !*FD Snow depth (m)
    real(dp), intent(inout) :: sicedepth     !*FD Superimposed ice depth (m)

    real(dp) :: fracice,fracsnow,firndepth,equfdepth,newdepth

    firndepth = snowdepth+sicedepth
    fracice   = sicedepth/firndepth
    fracsnow  = snowdepth/firndepth

    if (fracsnow /= 0.d0) then
       equfdepth=(1.d0/(params%constC*rhow))*( &
            rhoi*log((fracsnow*(rhoi-params%snowdensity))/(rhoi*(1.d0-params%firnbound))) &
            +(fracice-params%firnbound)*rhoi+params%snowdensity*fracsnow)
    else
       equfdepth = -1.d0
    end if

    if (equfdepth >= 0.d0 .and. equfdepth < firndepth) then
       newdepth  = (params%a1*equfdepth+params%a2*firndepth)/params%a3
       snowdepth = fracsnow*newdepth
       sicedepth = fracice*newdepth
    end if

  end subroutine firn_densify

end module glint_daily_pdd
