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
  USE COSP_KINDS, ONLY: wp
  IMPLICIT NONE
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
  SUBROUTINE COSP_RTTOV_INIT()

  END SUBROUTINE COSP_RTTOV_INIT
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
