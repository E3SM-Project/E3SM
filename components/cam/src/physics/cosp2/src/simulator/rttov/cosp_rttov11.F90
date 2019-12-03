! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2016, Regents of the University of Colorado
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
! March 2016 - M. Johnston - Original version
! April 2016 - D. Swales   - Modified for use in COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_cosp_rttov
  use rttov_const,         only : errorstatus_success, errorstatus_fatal
  use rttov_types,         only : rttov_options,rttov_coefs,profile_type,                &
                                  transmission_type,radiance_type,rttov_chanprof,        &
                                  rttov_emissivity,profile_cloud_type,rttov_scatt_coef,  &
                                  rttov_options_scatt
  use rttov_const,         only : surftype_sea, surftype_land, surftype_seaice
  use rttov_unix_env,      only : rttov_exit
  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : RTTOV_MAX_CHANNELS,N_HYDRO,rttovDir
  use cosp_phys_constants, only : mdry=>amd,mO3=>amO3,mco2=>amCO2,mCH4=>amCH4,           &
                                  mn2o=>amN2O,mco=>amCO
  implicit none
#include "rttov_direct.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_read_coefs.interface"
#include "rttov_get_emis.interface"
#include "rttov_boundaryconditions.interface"

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

  ! Scattering coefficients (read in once during initialization)
  type(rttov_coefs) :: &
       coef_rttov
  type(rttov_scatt_coef) :: &
       coef_scatt
  ! RTTOV setup and options (set during initialization)
  type(rttov_options) :: &
       opts     ! defaults to everything optional switched off
  type(rttov_options_scatt) :: &
       opts_scatt
contains

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_column
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    ! Local variables
    integer :: &
         nloop,rmod,il,istart,istop,za,i,j,subcol,errorstatus,npts_it
    integer,dimension(60) :: &
         alloc_status
    real(wp),dimension(nPoints) :: &
         sh_surf
    real(wp),dimension(nPoints,nLevels) :: &
         sh,totalice
    real(wp),dimension(nSubCols,nPoints,nChannels) :: &
         Tbs ! Subcolumn brightness temperature
    logical :: &
         use_totalice, mmr_snowrain, cfrac
    logical :: &
         lallSky, & ! Control for type of brightness temperature calculation
                    ! (False(default) => clear-sky brightness temperature, True => All-sky)
         lsfcEmis   ! Control for surface emissivity calculation (true => compute surface emissivity,
                    ! provided that the field "month" is available)

#include "rttov_read_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_setup_emis_atlas.interface"
#include "rttov_deallocate_emis_atlas.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_boundaryconditions.interface"

    ! Initialize some things
    totalice   = 0._wp
    Tbs(:,:,:) = 0._wp
    Tb(:,:)    = 0._wp
    error      = ''

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Setup for call to RTTOV
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! First, check to see if we are doing an all-sky or clear-sky calculation brightness
    ! temperature
    lallSky = .false.
    if (present(tca) .and. present(clw) .and. present(ciw) .and. present(rain)           &
        .and. present(snow)) lallSky=.true.

    ! Check to see if we need to compute the surface emissivity (defualt is to compute
    ! surface emissivity using the atlas tables)
    lsfcEmis = .true.
    if (present(surfem)) lsfcEmis = .false.
    
    ! We also need the month for the emissivity atlas, so check...
    if (.not. present(month)) lsfcEmis = .false.

    if (lsfcEmis .eq. .false. .and. .not. present(surfem)) then
       error = 'ERROR (rttov_column): User did not provide surface emissivity and did not '//&
               'request the surface emissivity to be calculated!!!'
       return
    endif

    ! Convert specific humidity to ppmv
    sh       = ( q   / ( q   + eps * ( 1._wp - q ) ) ) * 1e6   
    sh_surf  = ( q2m / ( q2m + eps * ( 1._wp - q2m ) ) ) * 1e6

    ! Settings unique to all-sky call.
    use_totalice = .false.
    mmr_snowrain = .true.
    cfrac        = .true.
    opts_scatt%lusercfrac = cfrac

    ! RTTOV can handle only about 100 profiles at a time (fixme: check this with roger),
    ! so we are putting a loop of 100
    nloop  =  npoints / maxlim
    rmod   =  mod( npoints, maxlim )
    if( rmod .ne. 0 ) then
       nloop = nloop + 1
    endif

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialize emissivity atlas data for chosen sensor.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call rttov_setup_emis_atlas(errorstatus,opts,month,coef_rttov,path=trim(rttovDir)//"emis_data/")
    if (errorstatus /= errorstatus_success) then
       error = 'ERROR (rttov_column): Error reading emis atlas data!'
       return
    endif
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Some quality control prior to RTTOV call
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Ensure the options and coefficients are consistent
    if(opts_scatt%config%do_checkinput) then
       call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
       if (errorstatus /= errorstatus_success) then
          error =  'ERROR (rttov_column): Error when checking input data!'
          return
       endif
    endif
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call to RTTOV
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Looping over maxlim number of profiles
    do il = 1, nloop
       istart  =  (il - 1) * maxlim + 1
       istop   =  min(il * maxlim, npoints)     
       if( ( il .eq. nloop ) .and. ( rmod .ne. 0 ) ) then
          npts_it = rmod
       else
          npts_it = maxlim
       endif
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Clear-sky brightness temperature
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       if (.not. lallSky) then
          call rttov_multprof(nChannels,iChannel,surfem,npts_it,nLevels,platform,         &
                              satellite,sensor,opts,coef_rttov,zenang,                    &
                              p(istart:istop,:)/100._wp,t(istart:istop,:),                &
                              sh(istart:istop,:),(mdry/mo3)*o3(istart:istop,:)*1e6,       &
                              (mdry/mco2)*co2*1e6,(mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,&
                              (mdry/mco)*co*1e6,h_surf(istart:istop),u_surf(istart:istop),&
                              v_surf(istart:istop),t_skin(istart:istop),                  &
                              p_surf(istart:istop)/100.,t2m(istart:istop),                &
                              sh_surf(istart:istop),lsmask(istart:istop),                 &
                              seaice(istart:istop),lat(istart:istop),lon(istart:istop),   &
                              Tb(istart:istop,:))  
       endif
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! All-sky brightness temperature
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (lallSky) then
          ! Loop over all subcolumns
          do subcol = 1, nSubCols	
             ! Call RTTOV
             call cosp_rttov_mwscatt(nChannels,iChannel,surfem,nPoints,nlevels,platform,  &
                                     satellite,sensor,opts,opts_scatt,coef_rttov,         &
                                     coef_scatt,zenang,p(istart:istop,:)/100._wp,         &
                                     ph(istart:istop,:)/100._wp,t(istart:istop, :),       &
                                     sh(istart:istop, :),                                 &
                                     (mdry/mo3)*o3(istart:istop,:)*1e6,                   &
                                     clw(istart:istop,subcol,:),                          &
                                     ciw(istart:istop,subcol,:),tca(istart:istop, :),     &
                                     totalice(istart:istop,:),snow(istart:istop,subcol,:),& 
                                     rain(istart:istop,subcol,:),(mdry/mco2)*co2*1e6,     &
                                     (mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,             &
                                     (mdry/mco)*co*1e6,h_surf(istart:istop),              &
                                     u_surf(istart:istop),v_surf(istart:istop),           &
                                     t_skin(istart:istop), p_surf(istart:istop)/100.,     &
                                     t2m(istart:istop),sh_surf(istart:istop),             &
                                     lsmask(istart:istop),seaice(istart:istop),           &
                                     lat(istart:istop),lon(istart:istop), use_totalice,   &
                                     mmr_snowrain,cfrac,Tbs(subcol,istart:istop,:))
          enddo 
       endif 
    enddo

    ! For all-sky calculation we need to average together all of the cloudy subcolumns.
    if (lallSky) then
       do subcol = 1, nSubCols
          Tb = Tb + tbs(subcol,:,:)
       enddo
       Tb = Tb/nSubCols
    endif
    
    ! Free up space
    if (lCleanup) then
       call rttov_dealloc_coefs(errorstatus,coef_rttov)
       call rttov_deallocate_emis_atlas(coef_rttov)
       if (lallSky) call rttov_dealloc_scattcoeffs(coef_scatt)    
    endif
  end subroutine rttov_column

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_multprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine rttov_multprof( &
       nch_in,     & ! number of channels
       ichan_in,   & ! channel indices
       surfem_in,  & ! surface emissivity values
       prf_num_in, & ! number of profiles to simulate
       nlevels_in, & ! number of pressure levels
       plat_in,    & ! platform number
       sat_in,     & ! satellite number
       sens_in,    & ! instrument number
       opts,       &
       coef_rttov, &
       zenang_in,  & ! zenith angle
       p_in,       & ! pressure [hpa]
       t_in,       & ! temperature [ k ]
       q_in,       & ! specific humidity [ ppmv ]
       o3_in,      & ! ozone vmr [ ppmv ]
       co2_in,     & ! co2 vmr [ ppmv ] *this is a single value*
       ch4_in,     & ! ch4 vmr [ ppmv ] *this is a single value*
       n2o_in,     & ! n2o vmr [ ppmv ] *this is a single value*
       co_in,      & ! co vmr [ ppmv ]  *this is a single value*
       h_surf,     & ! surface height [ m ]
       u_surf,     & ! u wind at 10 m  [ m/s ]
       v_surf,     & ! v wind at 10 m [ m/s ]
       t_skin,     & ! skin temperatre [ k ]
       p_surf,     & ! surface pressure
       t_surf,     & ! 1.5 m temperature [ k ]
       q_surf,     & ! 1.5 m specific humidity [ ppmv ]
       lsmask,     & ! land sea mask
       seaice,     & ! seaice fraction
       latitude,   & ! latitude [ deg north ]
       longitude,  & ! longitude [ deg east ]
       tbs         & ! brightness temperature [ k ] (output)
       )

    !------ input arguments. no rttov kinds should be used here -----------------
    integer, intent(in)  :: nch_in             ! number of channels to be computed
    integer, intent(in)  :: ichan_in(nch_in)   ! indices of selected channels
    real(wp), intent(in)     :: surfem_in(nch_in)  ! surface emissivities for the channels
    integer, intent(in)  :: prf_num_in
    integer, intent(in)  :: nlevels_in
    integer, intent(in)  :: plat_in  ! satellite platform
    integer, intent(in)  :: sat_in   ! satellite number
    integer, intent(in)  :: sens_in  ! satellite sensor
    real(wp), intent(in)     :: zenang_in          ! satellite zenith angle

    type(rttov_options)  :: opts
    type(rttov_coefs)    :: coef_rttov

    real(wp), intent(in)     :: p_in(prf_num_in, nlevels_in)  ! pressure profiles
    real(wp), intent(in)     :: t_in(prf_num_in, nlevels_in)  ! temperature profiles
    real(wp), intent(in)     :: q_in(prf_num_in, nlevels_in)  ! humidity profiles
    real(wp), intent(in)     :: o3_in(prf_num_in, nlevels_in) ! ozone profiles

    ! the following trace gases contain constant values
    real(wp), intent(in) ::  co2_in ! carbon dioxide
    real(wp), intent(in) ::  ch4_in ! methane
    real(wp), intent(in) ::  n2o_in ! n2o
    real(wp), intent(in) ::  co_in  ! carbon monoxide
    real(wp), intent(in) ::  h_surf(prf_num_in)         ! surface height
    real(wp), intent(in) ::  u_surf(prf_num_in)         ! u component of surface wind
    real(wp), intent(in) ::  v_surf(prf_num_in)         ! v component of surface wind
    real(wp), intent(in) ::  t_skin(prf_num_in)         ! surface skin temperature
    real(wp), intent(in) ::  p_surf(prf_num_in)                  ! surface pressure
    real(wp), intent(in) ::  t_surf(prf_num_in)         ! 1.5 m temperature
    real(wp), intent(in) ::  q_surf(prf_num_in)         ! 1.5 m specific humidity
    real(wp), intent(in) ::  lsmask(prf_num_in)         ! land-sea mask
    real(wp), intent(in) ::  seaice(prf_num_in)         ! sea-ice fraction
    real(wp), intent(in) ::  latitude(prf_num_in)       ! latitude
    real(wp), intent(in) ::  longitude(prf_num_in)      ! longitude

    real(wp), intent(inout) :: tbs(prf_num_in, nch_in)  ! tbs (in the right format)

    !------ local variables. use only rttov kinds or derived types.
    !       logical variables are declared with the same kind
    !       as integers, as they are affected inthe same way by flags like -qintsize=8

    !     type(rttov_options)                  :: opts           ! options structure
    !     type(rttov_coefs),       allocatable :: coefs(:)       ! coefficients structure
    type(rttov_chanprof),    allocatable :: chanprof(:)    ! input channel/profile list
    type(profile_type),      allocatable :: profiles(:)    ! input profiles
    logical,      allocatable :: calcemis(:)    ! flag to indicate calculation of emissivity within rttov
    type(rttov_emissivity),  allocatable :: emissivity(:)  ! input/output surface emissivity
    type(transmission_type)              :: transmission   ! output transmittances
    type(radiance_type)                  :: radiance       ! output radiances

    integer, allocatable :: instrument(:,:) ! instrument id (3 x n_instruments)
    integer, allocatable :: nchan(:) ! number of channels per instrument
    integer, allocatable :: ichan(:,:)   ! channel list per instrument

    integer :: asw
    integer :: mxchn
    integer :: nrttovid     ! maximum number of instruments
    integer :: no_id        ! instrument loop index
    integer :: i, j, jch
    integer :: nprof  ! number of calls to rttov
    integer :: nch ! intermediate variable
    integer :: errorstatus
    integer :: ich, ich_temp, nchanprof, nchannels, chan
    integer :: alloc_status(60)

    real(wp),    allocatable :: input_emissivity (:)

    character (len=14)  :: nameofroutine = 'rttov_multprof'

    logical :: refrac, solrad, laerosl, lclouds, lsun, all_channels

    ! local variables for input arguments that need type casting to avoid type-mismatch with
    ! rttov kinds. this happens with some compiler flags (-qintsize=8).
    integer  :: prof_num
    integer  :: nlevels

    ! --------------------------------------------------------------------------
    ! 0. initialise cosp-specific things
    ! --------------------------------------------------------------------------

    ! type-casting of input arguments that need to be passed to rttov
    prof_num = prf_num_in
    nlevels  = nlevels_in
    nprof = prof_num

    ! currently we plan to calculate only 1 instrument per call
    nrttovid  =  1
    mxchn  =  nch_in

    errorstatus     = 0
    alloc_status(:) = 0

    !     allocate(coefs(nrttovid), stat = alloc_status(1))

    !     allocate(instrument(3, nrttovid), stat = alloc_status(4))

    !maximum number of channels allowed for one instrument is mxchn
    !    allocate(surfem(nch_in, nrttovid), stat = alloc_status(11))
    allocate(ichan(nch_in, nrttovid), stat = alloc_status(12))
    call rttov_error('ichan mem allocation error for profile array' , lalloc = .true.)


    do no_id = 1, nrttovid
       ichan(:, no_id)   = ichan_in
    enddo

    asw = 1 ! switch for allocation passed into rttov subroutines

    ! allocate input profile arrays
    allocate(profiles(nprof), stat = alloc_status(1))
    call rttov_error('Profile mem allocation error' , lalloc = .true.)

    call rttov_alloc_prof(     &
         errorstatus,             &
         nprof,                   &
         profiles,                &
         nlevels,                 &
         opts,                    &
         asw,                     &
         coefs = coef_rttov,      &
         init = .true.)
    call rttov_error('Profile 2 mem allocation error' , lalloc = .true.)
    ! --------------------------------------------------------------------------
    ! 5. store profile data in profile type
    ! --------------------------------------------------------------------------
    do i = 1, nprof
       profiles(i)%p(:) =  p_in(i, :)
       profiles(i)%t(:) =  t_in(i, :)
       profiles(i)%q(:) =  q_in(i, :)

       where(profiles(i)%q(:) < 1e-4)
          profiles(i)%q(:) = 1e-4
       end where

       profiles(i)%cfraction  =  0.
       profiles(i)%ctp        =  500.

       ! 2m parameters
       profiles(i)%s2m%p  =  p_surf(i)
       profiles(i)%s2m%t  =  t_surf(i)
       profiles(i)%s2m%q  =  q_surf(i)
       profiles(i)%s2m%u  =  u_surf(i) ! dar: hard-coded at 2ms-1?
       profiles(i)%s2m%v  =  v_surf(i) ! dar: hard-coded at 2ms-1?
       profiles(i)%s2m%wfetc  =  10000. ! dar: default?

       ! skin variables for emissivity calculations
       profiles(i)%skin%t          =  t_skin(i)

       ! fastem coefficients - for mw calculations
       profiles(i)%skin%fastem(1)  =  3.0
       profiles(i)%skin%fastem(2)  =  5.0
       profiles(i)%skin%fastem(3)  =  15.0
       profiles(i)%skin%fastem(4)  =  0.1
       profiles(i)%skin%fastem(5)  =  0.3

       profiles(i)%zenangle      = zenang_in ! pass in from cosp

       profiles(i)%azangle       = 0. ! hard-coded in rttov9 int

       profiles(i)%latitude      = latitude(i)
       profiles(i)%longitude     = longitude(i)
       profiles(i)%elevation     = h_surf(i)

       profiles(i)%sunzenangle   = 0. ! hard-coded in rttov9 int
       profiles(i)%sunazangle    = 0. ! hard-coded in rttov9 int

       ! surface type
       ! land-sea mask indicates proportion of land in grid
       if (lsmask(i) < 0.5) then
          profiles(i)%skin%surftype  = surftype_sea
       else
          profiles(i)%skin%surftype  = surftype_land
       endif
       ! sea-ice fraction
       if (seaice(i) >= 0.5) then
          profiles(i)%skin%surftype  = surftype_seaice
       endif

       ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
       profiles(i)%skin%watertype = 1
       profiles(i) %idg         = 0.
       profiles(i) %ish         = 0.
    enddo
    ! end of 5.

    ich_temp = 1
    nchannels = nch_in
    do no_id = 1, nrttovid

       ! --------------------------------------------------------------------------
       ! 3. build the list of profile/channel indices in chanprof
       ! --------------------------------------------------------------------------

       allocate(nchan(nprof))     ! number of channels per profile
       nchan(:) = size(ichan(:,no_id))  ! = nch_in

       ! size of chanprof array is total number of channels over all profiles
       ! square in this case - here same channels done for all profiles
       nchanprof = sum(nchan(:))

       ! pack channels and input emissivity arrays
       allocate(chanprof(nchanprof))
       !      allocate(emis(nchanprof))
       chanprof(:)%chan =0
       
       nch = 0
       do j = 1, nprof
          do jch = 1, nchan(j)
             nch = nch + 1
             chanprof(nch)%prof = j
             if(ichan(jch, no_id) < 1) then
                errorstatus = errorstatus_fatal
                call rttov_error('Sensor channel number must be 1 or greater' , lalloc = .true.)
             else
                chanprof(nch)%chan = ichan(jch, no_id)
             endif
          enddo
       enddo
       ! end of 3.

       ! allocate output radiance arrays
       call rttov_alloc_rad( &
            errorstatus,        &
            nchanprof,          &
            radiance,           &
            nlevels - 1,   & ! nlayers
            asw)
       call rttov_error('allocation error for radiance arrays' , lalloc = .true.)

       ! allocate transmittance structure
       call rttov_alloc_transmission( &
            errorstatus,             &
            transmission,            &
            nlevels - 1,        &
            nchanprof,               &
            asw,                     &
            init=.true.)
       call rttov_error('allocation error for transmission arrays' , lalloc = .true.)

       ! allocate arrays for surface emissivity
       allocate(calcemis(nchanprof), stat=alloc_status(1))
       allocate(emissivity(nchanprof), stat=alloc_status(2))
       call rttov_error('mem allocation error for emissivity arrays' , lalloc = .true.)

       call rttov_get_emis(       &
            & errorstatus, &
            & opts,        &
            & chanprof,    &
            & profiles,    &
            & coef_rttov,  &
                                !& resolution=resolution, &  ! *** MW atlas native
                                !                  resolution is 0.25 degree lat/lon; if you know better
                                !                  value for satellite footprint (larger than this) then
                                !                  you can specify it here
            & emissivity=emissivity(:)%emis_in)
       !            & emissivity(:)%emis_in)

       call rttov_error('Get emissivity error' , lalloc = .true.)
       calcemis(:) = .false.
       ! calculate emissivity for missing and ocean location (fastem)
       where (emissivity(:)%emis_in <= 0.0)
          calcemis(:) = .true.
       endwhere

       call rttov_direct(         &
            errorstatus,               &! out
            chanprof,                  &
            opts,                      &
            profiles,                  &! in
            coef_rttov,              &! in
            transmission,              &! out
            radiance,                  &
            calcemis = calcemis,       &! in
            emissivity = emissivity)    ! inout
       call rttov_error('rttov_direct error', lalloc = .true.)

       tbs(1:prof_num , ich_temp:ich_temp + size(ichan(:,no_id)) - 1) = &
            transpose(reshape(radiance%bt(1:nchanprof), (/ size(ichan(:,no_id)), prof_num/) ))

       ich_temp = ich_temp + size(ichan(:,no_id))

       ! --------------------------------------------------------------------------
       ! 8. deallocate all rttov arrays and structures
       ! --------------------------------------------------------------------------
       deallocate (nchan,       stat=alloc_status(3))
       deallocate (chanprof,    stat=alloc_status(4))
       deallocate (emissivity,  stat=alloc_status(5))
       deallocate (calcemis,    stat=alloc_status(6))
       call rttov_error('rttov array deallocation error', lalloc = .true.)

       asw = 0 ! switch for deallocation passed into rttov subroutines

       ! deallocate radiance arrays
       call rttov_alloc_rad(errorstatus, nchannels, radiance, nlevels - 1, asw)
       call rttov_error('radiance deallocation error', lalloc = .true.)

       ! deallocate transmission arrays
       call rttov_alloc_transmission(errorstatus, transmission, nlevels - 1, nchannels, asw)
       call rttov_error('transmission deallocation error', lalloc = .true.)

    enddo

    ! deallocate profile arrays
    call rttov_alloc_prof(errorstatus, nprof, profiles, nlevels, opts, asw)
    call rttov_error('profile deallocation error', lalloc = .true.)

    deallocate(profiles, stat=alloc_status(1))
    call rttov_error('mem deallocation error for profile array', lalloc= .true.)

  contains

    subroutine rttov_error(msg, lalloc)
      character(*) :: msg
      logical  :: lalloc

      if(lalloc) then
         if (any(alloc_status /= 0)) then
            write(*,*) msg
            errorstatus = 1
            call rttov_exit(errorstatus)
         endif
      else
         if (errorstatus /= errorstatus_success) then
            write(*,*) msg
            call rttov_exit(errorstatus)
         endif
      endif
    end subroutine rttov_error

  end subroutine rttov_multprof
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !----------------- subroutine cosp_rttov_mwscatt ---------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_mwscatt(&
       nch_in,     & ! number of channels
       ichan_in,   & ! channel indices
       surfem_in,  & ! surface emissivity values
       prf_num_in, & ! number of profiles to simulate
       nlevels_in, & ! number of pressure levels
       plat_in,    & ! platform number
       sat_in,     & ! satellite number
       sens_in,    & ! instrument number
       opts,       &
       opts_scatt, &
       coef_rttov, &
       coef_scatt, &
       zenang_in,  & ! zenith angle
       p_in,       & ! pressure [hpa]
       ph_in,      & ! pressure on half levels [hpa]
       t_in,       & ! temperature [ k ]
       q_in,       & ! specific humidity [ ppmv ]
       o3_in,      & ! ozone vmr [ ppmv ]
       clw_in,     & ! cloud water [0-1]
       ciw_in,     & ! cloud ice [0-1]
       cc_in,      & ! effective cloud fraction [0-1]
       totalice_in,& ! total ice, except snow [kg/kg] or [kg/m2/s]
       sp_in,      & ! solid precip with snow [kg/kg] or [kg/m2/s]
       rain_in,    & ! total liquid water [kg/kg] or [kg/m2/s]
       co2_in,     & ! co2 vmr [ ppmv ] *this is a single value*
       ch4_in,     & ! ch4 vmr [ ppmv ] *this is a single value*
       n2o_in,     & ! n2o vmr [ ppmv ] *this is a single value*
       co_in,      & ! co vmr [ ppmv ]  *this is a single value*
       h_surf,     & ! surface height [ m ]
       u_surf,     & ! u wind at 10 m  [ m/s ]
       v_surf,     & ! v wind at 10 m [ m/s ]
       t_skin,     & ! skin temperatre [ k ]
       p_surf,     & ! surface pressure
       t_surf,     & ! 1.5 m temperature [ k ]
       q_surf,     & ! 1.5 m specific humidity [ ppmv ]
       lsmask,     & ! land sea mask
       seaice,     & ! seaice fraction
       latitude,   & ! latitude [ deg north ]
       longitude,  & ! longitude [ deg east ]
       use_totalice,& ! separate ice and snow, or total ice hydrometeor types
       mmr_snowrain,& ! set units for snow and rain: if true units are kg/kg (the default)
       cfrac,       & ! opts_scatt%lusercfrac=true., supply the effective cloud fraction
       tbs          & ! brightness temperature [ k ] (output)
       )





    implicit none

    !------ input arguments. no rttov kinds should be used here -----------------
    integer, intent(in)  :: nch_in             ! number of channels to be computed
    integer, intent(in)  :: ichan_in(nch_in)   ! indices of selected channels
    real(wp), intent(in)     :: surfem_in(nch_in)  ! surface emissivities for the channels
    integer, intent(in)  :: prf_num_in
    integer, intent(in)  :: nlevels_in
    integer, intent(in)  :: plat_in  ! satellite platform
    integer, intent(in)  :: sat_in   ! satellite number
    integer, intent(in)  :: sens_in  ! satellite sensor
    real(wp), intent(in)     :: zenang_in          ! satellite zenith angle

    type(rttov_options)       :: opts
    type(rttov_options_scatt) :: opts_scatt
    type(rttov_coefs)         :: coef_rttov
    type(rttov_scatt_coef)    :: coef_scatt

    real(wp), intent(in)     :: p_in(prf_num_in, nlevels_in)  ! pressure profiles
    real(wp), intent(in)     :: t_in(prf_num_in, nlevels_in)  ! temperature profiles
    real(wp), intent(in)     :: q_in(prf_num_in, nlevels_in)  ! humidity profiles
    real(wp), intent(in)     :: o3_in(prf_num_in, nlevels_in) ! ozone profiles
    real(wp), intent(in)     :: clw_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: ciw_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: cc_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: totalice_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: sp_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: rain_in(prf_num_in, nlevels_in)
    real(wp), intent(in)     :: ph_in(prf_num_in, nlevels_in+1)

    ! the following trace gases contain constant values
    real(wp), intent(in) ::  co2_in ! carbon dioxide
    real(wp), intent(in) ::  ch4_in ! methane
    real(wp), intent(in) ::  n2o_in ! n2o
    real(wp), intent(in) ::  co_in  ! carbon monoxide
    real(wp), intent(in) ::  h_surf(prf_num_in)         ! surface height
    real(wp), intent(in) ::  u_surf(prf_num_in)         ! u component of surface wind
    real(wp), intent(in) ::  v_surf(prf_num_in)         ! v component of surface wind
    real(wp), intent(in) ::  t_skin(prf_num_in)         ! surface skin temperature
    real(wp), intent(in) ::  p_surf(prf_num_in)                  ! surface pressure
    real(wp), intent(in) ::  t_surf(prf_num_in)         ! 1.5 m temperature
    real(wp), intent(in) ::  q_surf(prf_num_in)         ! 1.5 m specific humidity
    real(wp), intent(in) ::  lsmask(prf_num_in)         ! land-sea mask
    real(wp), intent(in) ::  seaice(prf_num_in)         ! seaice fraction
    real(wp), intent(in) ::  latitude(prf_num_in)       ! latitude
    real(wp), intent(in) ::  longitude(prf_num_in)      ! longitude
    logical, intent(in) :: cfrac, use_totalice, mmr_snowrain

    real(wp), intent(inout) :: tbs(prf_num_in, nch_in)  ! tbs (in the right format)
    !****************** local variables **********************************************
    logical       , allocatable :: calcemis   (:)
    type(rttov_emissivity)   , allocatable :: emissivity  (:)
    integer       , allocatable :: frequencies (:)
    type(rttov_chanprof)     , allocatable :: chanprof    (:) ! channel and profile indices
    type(profile_type)       , allocatable :: profiles    (:)
    type(profile_cloud_type) , allocatable :: cld_profiles(:)

    integer, allocatable :: ichan(:,:)   ! channel list per instrument

    integer         :: errorstatus
    type (radiance_type)       :: radiance
    ! 	type (rttov_options)       :: opts     ! defaults to everything optional switched off
    ! 	type (rttov_options_scatt) :: opts_scatt
    ! 	type (rttov_coefs)         :: coef_rttov
    ! 	type (rttov_scatt_coef)    :: coef_scatt

    ! 	integer, allocatable :: instrument (:,:)
    integer :: j,k,asw
    integer :: nchanxnprof, ninstruments
    real(wp)     :: zenangle
    character (len=256) :: outstring
    integer :: alloc_status(60)

#include "rttov_init_rad.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_scatt.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_get_emis.interface"
#include "rttov_boundaryconditions.interface"

    errorstatus     = 0
    alloc_status(:) = 0
    ninstruments = 1 ! number of sensors or platforms

    allocate(ichan(nch_in, ninstruments), stat = alloc_status(3))

    do j = 1, ninstruments
       ichan(:, j) = ichan_in
    enddo

    nchanxnprof =  prf_num_in * nch_in            ! total channels to simulate * profiles

    allocate (chanprof(nchanxnprof))
    allocate (frequencies(nchanxnprof))
    allocate (emissivity(nchanxnprof))
    allocate (calcemis(nchanxnprof))
    allocate (profiles(prf_num_in))
    allocate (cld_profiles(prf_num_in))

    ! request rttov / fastem to calculate surface emissivity
    calcemis  = .true.
    emissivity % emis_in = 0.0

    ! setup indices
    call rttov_scatt_setupindex ( 	&
	 & prf_num_in,      			& ! in
	 & nch_in,          			& ! in
	 & coef_rttov%coef, 			& ! in
	 & nchanxnprof,       			& ! in
	 & chanprof,        			& ! out
	 & frequencies)       		      ! out

    ! allocate profiles (input) and radiance (output) structures
    asw = 1
    call rttov_alloc_prof( errorstatus,prf_num_in,profiles,nlevels_in,opts,asw, init = .true.)
    call rttov_alloc_scatt_prof(prf_num_in,cld_profiles, nlevels_in, .false., 1, init = .true.)
    call rttov_alloc_rad(errorstatus,nchanxnprof,radiance,nlevels_in-1,asw)

    ! fill the profile structures with data
    do j = 1, prf_num_in
       profiles(j)%latitude    = latitude(j)
       profiles(j)%longitude   = longitude(j)
       profiles(j)%elevation   = h_surf(j)
       profiles(j)%sunzenangle = 0.0 ! hard-coded in rttov9 int
       profiles(j)%sunazangle  = 0.0 ! hard-coded in rttov9 int
       profiles(j)%azangle     = 0.0
       profiles(j)%zenangle    = zenang_in
       profiles(j)%s2m%t       = t_surf(j)
       profiles(j)%s2m%q       = q_surf(j)
       profiles(j)%s2m%u       = u_surf(j)
       profiles(j)%s2m%v       = v_surf(j)
       profiles(j)%s2m%wfetc   = 10000.
       profiles(j)%skin%t      = t_skin(j)
       profiles(j)%skin%watertype  = 1 ! ocean water
       if (lsmask(j) < 0.5) then
          profiles(j)%skin%surftype  = surftype_sea
       else
          profiles(j)%skin%surftype  = surftype_land
       endif
       if (seaice(j) >= 0.5) then
          profiles(j)%skin%surftype  = surftype_seaice
       endif
       profiles(j)%skin%fastem(1)  =  3.0
       profiles(j)%skin%fastem(2)  =  5.0
       profiles(j)%skin%fastem(3)  =  15.0
       profiles(j)%skin%fastem(4)  =  0.1
       profiles(j)%skin%fastem(5)  =  0.3
       profiles(j)%cfraction       = 0.0
       profiles(j)%ctp 	           = 500.0 ! not used but still required by rttov
       profiles(j)%p(:)            = p_in(j,:)
       profiles(j)%t(:)            = t_in(j,:)
       profiles(j)%q(:)	           = q_in(j,:)
       profiles(j)%idg             = 0.
       profiles(j)%ish             = 0.
       where(profiles(j)%q(:) < 1e-4)
          profiles(j)%q(:) = 1e-4
       end where
       cld_profiles(j)%ph(:)   = ph_in(j,:)
       cld_profiles(j)%cc(:)   = cc_in(j,:)
       cld_profiles(j)%clw(:)  = clw_in(j,:)
       cld_profiles(j)%ciw(:)  = ciw_in(j,:)
       cld_profiles(j)%rain(:) = rain_in(j,:)
       cld_profiles(j)%sp(:)   = sp_in(j,:)
       profiles(j)%s2m%p       = cld_profiles(j)%ph(nlevels_in+1)
    enddo

    call rttov_get_emis(  &
         & errorstatus, &
         & opts,       &
         & chanprof,   &
         & profiles,   &
         & coef_rttov,      &
         !               & resolution=resolution, &  ! *** MW atlas native resolution is
         !                 0.25 degree lat/lon; if you know better value for satellite
         !                 footprint (larger than this) then you can specify it here
         & emissivity=emissivity(:)%emis_in)
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'In COSP_RTTOV11: Error RTTOV_GET_EMIS!'
       call rttov_exit(errorstatus)
    endif

    calcemis(:) = .false.
    where (emissivity(:)%emis_in <= 0.)
       calcemis(:) = .true.
    endwhere

    call rttov_scatt (&
         & errorstatus,         &! out
         & opts_scatt,          &! in
         & nlevels_in,          &! in
         & chanprof,            &! in
         & frequencies,         &! in
         & profiles,            &! in
         & cld_profiles,        &! in
         & coef_rttov,          &! in
         & coef_scatt,          &! in
         & calcemis,           &! in
         & emissivity,          &! in
         & radiance)             ! out

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'In COSP_RTTOV11: Error RTTOV_SCATT!'
       call rttov_exit(errorstatus)
    endif

    !write(*,*) 'Checking emissivities: ', maxval(emissivity(:)%emis_out), \
    !         minval(emissivity(:)%emis_out)
    tbs(1:prf_num_in,1:1+size(ichan(:,1))-1) = &
         transpose(reshape(radiance%bt(1:nchanxnprof),(/ size(ichan(:,1)),prf_num_in/) ))

    ! deallocate all storage
    asw = 0
    ! 	call rttov_dealloc_coefs(errorstatus,coef_rttov)
    ! 	call rttov_dealloc_scattcoeffs(coef_scatt)
    call rttov_alloc_prof(errorstatus,prf_num_in,profiles,nlevels_in,opts,asw)
    call rttov_alloc_scatt_prof(prf_num_in,cld_profiles,nlevels_in,.false.,asw)
    call rttov_alloc_rad(errorstatus,nchanxnprof,radiance,nlevels_in-1,asw)
    deallocate(ichan,chanprof,frequencies,emissivity,calcemis)  !instrument,
    !***************************************************************************
    !-------- end section --------
    !***************************************************************************
  end subroutine cosp_rttov_mwscatt
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
  
end module mod_cosp_rttov
