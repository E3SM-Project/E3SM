! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Original version
! Apr 2015 - D. Swales - Modified for RTTOVv11.3
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,       ONLY: wp
  USE MOD_COSP_CONFIG,  ONLY: RTTOV_MAX_CHANNELS,rttovDir
  use mod_cosp_rttov,   only: platform,satellite,sensor,nChannels,iChannel,coef_rttov,   &
                              coef_scatt,opts,opts_scatt,construct_rttov_coeffilename,   &
                              construct_rttov_scatfilename
  IMPLICIT NONE
#include "rttov_read_coefs.interface"
#include "rttov_read_scattcoeffs.interface"

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type rttov_in
     integer,pointer :: &
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols,     & ! Number of subcolumns
          month           ! Month (needed for surface emissivity calculation)
     real(wp),pointer :: &
          zenang,       & ! Satellite zenith angle
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
     real(wp),dimension(:),pointer :: &
          surfem          ! Surface emissivities for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,       & ! Surface height
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t_skin,       & ! Surface skin temperature
          p_surf,       & ! Surface pressure
          t2m,          & ! 2 m Temperature
          q2m,          & ! 2 m Specific humidity
          lsmask,       & ! land-sea mask
          latitude,     & ! Latitude
          longitude,    & ! Longitude
          seaice          ! Sea-ice? 
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure @ model levels
          ph,           & ! Pressure @ model half levels
          t,            & ! Temperature 
          q,            & ! Specific humidity
          o3              ! Ozone
     
     ! These fields below are needed ONLY for the RTTOV all-sky brightness temperature
     real(wp),dimension(:,:),pointer :: &
          tca,          & ! Cloud fraction
          cldIce,       & ! Cloud ice
          cldLiq,       & ! Cloud liquid
          fl_rain,      & ! Precipitation flux (startiform+convective rain) (kg/m2/s)
          fl_snow         ! Precipitation flux (stratiform+convective snow)
  end type rttov_in
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(NchanIN,platformIN,satelliteIN,instrumentIN,channelsIN)
    integer,intent(in) :: & 
         NchanIN,     & ! Number of channels
         platformIN,  & ! Satellite platform
         satelliteIN, & ! Satellite
         instrumentIN   ! Instrument
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         channelsIN     ! RTTOV channels

    ! Local variables
    character(len=256) :: coef_file,scat_file
    integer :: errorstatus
    
    ! Initialize fields in module memory (cosp_rttovXX.F90)
    nChannels  = NchanIN
    platform   = platformIN 
    satellite  = satelliteIN 
    sensor     = instrumentIN 
    iChannel   = channelsIN

    ! Options common to RTTOV clear-sky Tb calculation
    opts%interpolation%addinterp  = .true.  ! allow interpolation of input profile
    opts%rt_all%use_q2m           = .true.
    opts%config%do_checkinput     = .false.
    opts%config%verbose           = .false.
    opts%rt_all%addrefrac         = .true.  ! include refraction in path calc
    opts%interpolation%reg_limit_extrap = .true.
    ! Options common to RTTOV clear-sky Tb calculation
    opts_scatt%config%do_checkinput    = .false.
    opts_scatt%config%verbose          = .false.
    opts_scatt%config%apply_reg_limits = .true.
    opts_scatt%interp_mode             = 1
    opts_scatt%reg_limit_extrap        = .true.
    opts_scatt%use_q2m                 = .true.
    opts%rt_mw%clw_data                = .true. 
        
    ! Read in scattering coefficient file.
    coef_file = trim(rttovDir)//"rtcoef_rttov11/rttov7pred54L/"// &
         trim(construct_rttov_coeffilename(platform,satellite,sensor))
    call rttov_read_coefs(errorstatus,coef_rttov, opts, file_coef=trim(coef_file))

    ! Read in scattering (clouds+aerosol) coefficient file. *ONLY NEEDED IF DOING RTTOV ALL-SKY.*
    !scat_file = trim(rttovDir)//"rtcoef_rttov11/cldaer/"//&
    !     trim(construct_rttov_scatfilename(platform,satellite,sensor))
    ! Can't pass filename to rttov_read_scattcoeffs!!!!!
    !call rttov_read_scattcoeffs (errorstatus, coef_rttov%coef, coef_scatt,)
 
  END SUBROUTINE COSP_RTTOV_INIT
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
