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
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE MOD_COSP_RTTOV
  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : RTTOV_MAX_CHANNELS,N_HYDRO,rttovDir
  use cosp_phys_constants, only : mdry=>amd,mO3=>amO3,mco2=>amCO2,mCH4=>amCH4,           &
                                  mn2o=>amN2O,mco=>amCO
  IMPLICIT NONE

  ! Module parameters
  integer, parameter :: maxlim =  10000
  real(wp),parameter :: eps    =  0.622

  ! Initialization parameters
  integer :: &
       platform,   & ! RTTOV platform
       sensor,     & ! RTTOV instrument
       satellite,  & ! RTTOV satellite
       nChannels     ! Number of channels
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       iChannel      ! RTTOV channel numbers

CONTAINS
  subroutine rttov_column(nPoints,nLevels,nSubCols,q,p,t,o3,ph,h_surf,u_surf,v_surf,     &
                          p_surf,t_skin,t2m,q2m,lsmask,lon,lat,seaice,co2,ch4,n2o,co,    &
                          zenang,lCleanup,                                               &
                          ! Outputs
                          Tb,error,                                                      &
                          ! Optional arguments for surface emissivity calculation.
                          surfem,month,                                                  &
                          ! Optional arguments for all-sky calculation.
                          tca,ciw,clw,rain,snow)
    ! Inputs
    integer,intent(in) :: &
         nPoints, & ! Number of gridpoints
         nLevels, & ! Number of vertical levels
         nSubCols   ! Number of subcolumns
    real(wp),intent(in) :: &
         co2,     & ! CO2 mixing ratio (kg/kg)
         ch4,     & ! CH4 mixing ratio (kg/kg)
         n2o,     & ! N2O mixing ratio (kg/kg)
         co,      & ! CO mixing ratio (kg/kg)
         zenang     ! Satellite zenith angle
    real(wp),dimension(nPoints),intent(in) :: &
         h_surf,  & ! Surface height (m)
         u_surf,  & ! Surface u-wind (m/s)
         v_surf,  & ! Surface v-wind (m/s)
         p_surf,  & ! Surface pressure (Pa)
         t_skin,  & ! Skin temperature (K)
         t2m,     & ! 2-meter temperature (K)
         q2m,     & ! 2-meter specific humidity (kg/kg)
         lsmask,  & ! Land/sea mask
         lon,     & ! Longitude (deg)
         lat,     & ! Latitude (deg)
         seaice     ! Seaice fraction (0-1)
    real(wp),dimension(nPoints,nLevels),intent(in) :: &
         q,       & ! Specific humidity (kg/kg)
         p,       & ! Pressure(Pa)
         t,       & ! Temperature (K)
         o3         ! Ozone
    real(wp),dimension(nPoints,nLevels+1),intent(in) :: &
         ph         ! Pressure @ half-levels (Pa)
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types

    ! Optional inputs (Needed for surface emissivity calculation)
    integer,optional :: &
         month      ! Month (needed to determine table to load)
    real(wp),dimension(nChannels),optional :: &
         surfem     ! Surface emissivity for each RTTOV channel

    ! Optional inputs (Needed for all-sky calculation)
    real(wp),dimension(nPoints,nLevels),optional :: &
         tca       ! Total column cloud amount (0-1)
    real(wp),dimension(nPoints,nSubCols,nLevels),optional :: &
         ciw,    & ! Cloud ice
         clw,    & ! Cloud liquid
         rain,   & ! Precipitation flux (kg/m2/s)
         snow      ! Precipitation flux (kg/m2/s)

    ! Outputs
    real(wp),dimension(nPoints,nChannels) :: &
         Tb        ! RTTOV brightness temperature.
    character(len=128) :: &
         error     ! Error messages (only populated if error encountered)
    
  end subroutine rttov_column
  function construct_rttov_coeffilename(platform,satellite,instrument)
    ! Inputs
    integer,intent(in) :: platform,satellite,instrument
    ! Outputs
    character(len=256) :: construct_rttov_coeffilename
    ! Local variables
    character(len=256) :: coef_file
    integer :: error

    ! Initialize
    error = 0
    
    ! Platform
    if (platform .eq. 1)  coef_file = 'rtcoef_noaa_'
    if (platform .eq. 10) coef_file = 'rtcoef_metop_'
    if (platform .eq. 11) coef_file = 'rtcoef_envisat_'
    if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
       error=error+1
       write ( *,* ) 'Unsupported platform ID ',platform
       return
    endif

    ! Satellite
    if (satellite .lt. 10) then
       coef_file = trim(coef_file) // char(satellite+48)
    else if (satellite .lt. 100) then
       coef_file = trim(coef_file) // char(int(satellite/10)+48)
       coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
    else
       error=error+1
       write ( *,* ) 'Unsupported satellite number ',satellite
       return
    endif

    ! Sensor
    if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
    if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
    if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
    if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
       error=error+1
       write ( *,* ) 'Unsupported sensor number ', sensor
       return
    endif

    if (error .eq. 0) construct_rttov_coeffilename=coef_file
    
  end function construct_rttov_coeffilename
  function construct_rttov_scatfilename(platform,satellite,instrument)
    ! Inputs
    integer,intent(in) :: platform,satellite,instrument
    ! Outputs
    character(len=256) :: construct_rttov_scatfilename
    ! Local variables
    character(len=256) :: coef_file
    integer :: error

    ! Initialize
    error = 0
    
    ! Platform
    if (platform .eq. 1)  coef_file = 'sccldcoef_noaa_'
    if (platform .eq. 10) coef_file = 'sccldcoef_metop_'
    if (platform .eq. 11) coef_file = 'sccldcoef_envisat_'
    if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
       error=error+1
       write ( *,* ) 'Unsupported platform ID ',platform
       return
    endif

    ! Satellite
    if (satellite .lt. 10) then
       coef_file = trim(coef_file) // char(satellite+48)
    else if (satellite .lt. 100) then
       coef_file = trim(coef_file) // char(int(satellite/10)+48)
       coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
    else
       error=error+1
       write ( *,* ) 'Unsupported satellite number ',satellite
       return
    endif

    ! Sensor
    if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
    if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
    if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
    if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
       error=error+1
       write ( *,* ) 'Unsupported sensor number ', sensor
       return
    endif

    if (error .eq. 0) construct_rttov_scatfilename=coef_file
    
  end function construct_rttov_scatfilename
END MODULE MOD_COSP_RTTOV
