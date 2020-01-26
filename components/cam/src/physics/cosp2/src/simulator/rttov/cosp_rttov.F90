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
! History:
! Aug 2008 - V. John - Initial version
! Jun 2010 - A. Bodas-Salcedo - Conversion to module and tidy up
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: RTTOV_MAX_CHANNELS
  USE RTTOV_CONST,     only: errorstatus_fatal, errorstatus_warning, errorstatus_success
  USE RTTOV_TYPES,     only: rttov_coef, profile_type, transmission_type,                &
                             radiance_type, rttov_coef_scatt_ir, rttov_optpar_ir
  USE PARKIND1, Only : jpim, jprb
  IMPLICIT NONE

  ! Include subroutine interfaces
  include "rttov_errorreport.interface"
  include "rttov_setup.interface"
  include "rttov_errorhandling.interface"
  include "rttov_direct.interface"
  include "rttov_alloc_prof.interface"
  include "rttov_alloc_rad.interface"
  include "rttov_dealloc_coef.interface"
  
  ! Fields set during initialization
  integer :: &
     nch_in,  & ! Number of RTTOV channels
     plat_in, & ! RTTOV platform
     sat_in,  & ! RTTOV satellite
     sens_in    ! RTTOV instrument
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
     ichan_in   ! RTTOV channel indices
   
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                      SUBROUTINE RTTOV_PIXEL
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE RTTOV_subcolumn(surfem_in, prf_num_in, nlevels_in,         &
                         zenang_in, p_in,t_in, q_in, o3_in, co2_in, &
                         ch4_in, n2o_in, co_in, h_surf, u_surf,     &
                         v_surf, t_skin, p_surf, t_surf, q_surf,    &
                         lsmask, latitude, tbs)   
    ! INPUTS
    integer,intent(in) :: &
         prf_num_in,      & ! Number of profiles to simulate
         nlevels_in         ! Number of pressure levels
    real(wp),intent(in) :: &
         zenang_in,       & ! Satellite zenith angle
         co2_in,          & ! Carbon dioxide 
         ch4_in,          & ! Methane 
         n2o_in,          & ! n2o 
         co_in              ! Carbon monoxide
    real(wp),intent(in),dimension(nch_in) :: &
         surfem_in          ! Surface emissivities for the channels
    real(wp),intent(in),dimension(prf_num_in) :: &
         h_surf,          & ! Surface height
         u_surf,          & ! U component of surface wind
         v_surf,          & ! V component of surface wind
         t_skin,          & ! Surface skin temperature
         p_surf,          & ! Surface pressure
         t_surf,          & ! 1.5 m Temperature
         q_surf,          & ! 1.5 m Specific humidity
         lsmask,          & ! land-sea mask
         latitude           ! Latitude
    real(wp),intent(in),dimension(prf_num_in,nlevels_in) :: &
         p_in,            & ! Pressure profiles  
         t_in,            & ! Temperature profiles
         q_in,            & ! Humidity profiles
         o3_in              ! Ozone profiles
    
    ! OUTPUTS
    real(wp),intent(inout),dimension(prf_num_in,nch_in) :: &
         tbs                ! Tbs (in the right format)
    
    ! LOCAL VARIABLES
    type(transmission_type)                             :: transmission
    type(radiance_type)                                 :: radiance
    type(rttov_coef ),        allocatable, dimension(:) :: coef        
    type(profile_type),       allocatable, dimension(:) :: profiles
    type(rttov_coef_scatt_ir),allocatable, dimension(:) :: coef_scatt_ir
    type(rttov_optpar_ir),    allocatable, dimension(:) :: optp
    Integer(Kind=jpim), Allocatable, dimension(:,:) :: &
         instrument,         & ! Instrument id
         nchan,              & ! Number of channels per instrument and profile
         ichan                 ! Channel list per instrument
    Integer(Kind=jpim), Allocatable, dimension(:) :: &
         nchan1,             & ! Number of channels per instrument and profile
         nchannels,          & ! Number of channels per instrument 
         ifull,              & ! Full test (with TL,AD,K) per instrument
         nprof,              & ! Number of profiles per instrument
         nsurf,              & ! Surface id number per instrument
         nwater,             & ! Water id number per instrument
         channels,           & ! Channel list per instrument*profiles
         lprofiles,          & !
         rttov_errorstatus,  & ! rttov error return code
         setup_errorstatus     ! Setup return code
    Integer(Kind=jpim) :: &
         nprofiles,iref,isun,asw,mxchn,i,j,jch,errorstatus,io_status,ioout,interp, &
         Err_Unit,           & ! Logical error unit (<0 for default)
         verbosity_level,    & ! (<0 for default)
         nrttovid,           & ! Maximum number of instruments
         no_id,              & ! Instrument loop index
         Nprofs,             & ! Number of calls to RTTOV
         nch                   ! Intermediate variable
    Integer(Kind=jpim), Parameter :: &
         jpnav  =  31,       & ! Number of profile variables
         jpnsav =  6,        & ! Number of surface air variables
         jpnssv =  6,        & ! Number of skin variables
         jpncv  =  2,        & ! Number of cloud variables
         sscvar = jpnsav+jpnssv+jpncv ! Number of surface,skin,cloud vars
    Integer(Kind=jpim),dimension(60) :: &
         alloc_status    
    Real(Kind=jprb) :: &
         zenang, azang, sunzang, sunazang
    Real(kind=jprb), allocatable, dimension(:) :: &
         emissivity,fresnrefl,input_emissivity     
    Real(kind=jprb), allocatable, dimension(:,:) :: &
         surfem
    Real(Kind=jprb),dimension(nch_in*prf_num_in) :: &
         tbs_temp              ! A temporary variable to hold Tbs
    Character (len=3)  :: &
         cref, csun
    Character (len=80) :: &
         errMessage
    Character (len=14)  :: &
         NameOfRoutine = 'rttov_multprof'
    Logical :: addinterp,refrac,solrad,laerosl,lclouds,lsun,all_channels
    Logical,Allocatable,dimension(:) :: &
         calcemis
    integer(kind=jpim) :: prof_num,nlevels
    
    ! Type-casting of input arguments that need to be passed to RTTOV
    prof_num = prf_num_in
    nlevels  = nlevels_in
    
    ! Unit numbers for input/output
    IOOUT = 2
    
    ! Curretly we plan to calculate only 1 instrument per call
    nrttovid        =  1
    mxchn           =  nch_in
    errorstatus     = 0
    alloc_status(:) = 0
    all_channels    = .false.
    sunzang         = 0._jprb
    sunazang        = 0._jprb
    
    ! Initialise error management with default value for
    ! the error unit number and Fatal error message output
    Err_unit        = -1
    verbosity_level = 0
    
    ! All error message output
    call rttov_errorhandling( Err_unit, verbosity_level, print_checkinput_warnings=.false. )
    io_status  = 0
    errmessage = ''
    
    ! Assigning the zenith angle
    zenang = zenang_in
    
    ! Allocate
    allocate (coef(nrttovid),         stat = alloc_status(1))
    allocate (coef_scatt_ir(nrttovid),stat = alloc_status(2))
    allocate (optp(nrttovid),         stat = alloc_status(3))
    allocate (instrument(3,nrttovid), stat = alloc_status(4))
    allocate (ifull(nrttovid),        stat = alloc_status(5))
    allocate (nprof(nrttovid),        stat = alloc_status(6))
    allocate (nsurf(nrttovid),        stat = alloc_status(7))
    allocate (nwater(nrttovid),       stat = alloc_status(8))
    allocate (nchannels(nrttovid),    stat = alloc_status(9))
    allocate (nchan1(nrttovid),       stat = alloc_status(10))
    allocate (surfem(mxchn,nrttovid), stat = alloc_status(11))
    allocate (ichan (mxchn,nrttovid), stat = alloc_status(12))
    If( any(alloc_status /= 0) ) then
       errorstatus = errorstatus_fatal
       Write( errMessage, '( "mem allocation error")' )
       Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
       Stop
    End If
    surfem(:,:)  =  0.0_JPRB
    ichan(:,:)   =  0
    
    !!! FIXME: Shall we get rid of this loop? We use only one instrument
    DO NO_ID = 1, NRTTOVID    
       instrument(1,no_id) = plat_in
       instrument(2,no_id) = sat_in
       instrument(3,no_id) = sens_in
                
       !! Forward model only (0) or TL and AD (1) or K (2)?'
       !  This version supports only Forward model
       IFULL(no_id) = 0
       
       ! Number of profiles to test per call
       NPROF(no_id) = prof_num
       nprofiles    = NPROF(no_id)
       
       ! Total number of profiles to process
       NPROFS = prof_num
       NPROFS = NPROFS / NPROF(no_id)   ! Number of calls to RTTOV

       ! Check whether it is OK to use ocean all the time
       NWATER(no_id) = 1 ! Water type (0=fresh water, 1=ocean water)
       
       ! Set up channel numbers
       allocate (nchan(nprof(no_id),nrttovid),stat= alloc_status(3))
       nchan(1:nprof(no_id),no_id) = nch_in
       ichan(:, 1)   =  ichan_in
       surfem(:, 1)  =  surfem_in
       
       ! nchan(1,no_id) is now the real number of channels selected
       do j = 1 , nprof(no_id)
          nchan(j,no_id) = nch_in  
       enddo
       
       ! Compute channels*profiles
       nchannels(no_id) = 0
       Do j = 1 , nprof(no_id)
          nchannels(no_id) = nchannels(no_id) + nchan (j,no_id)
       End Do
       nchan1(no_id) = nchan(1,no_id)
    END DO
    
    ! Do you want clouds or aerosol?
    laerosl = .False.
    lclouds = .False.
    
    !#########################################################
    ! Beginning of rttov_setup test
    !#########################################################
    alloc_status = 0
    allocate ( setup_errorstatus(nrttovid),stat= alloc_status(1))
    If( any(alloc_status /= 0) ) then
       errorstatus = errorstatus_fatal
       Write( errMessage, '( "mem allocation error for errorsetup")' )
       Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
       Stop
    End If
    
    If (all_channels)Then
       Call rttov_setup (       &
            setup_errorstatus,  &! out
            Err_unit,           &! in
            verbosity_level,    &! in
            nrttovid,           &! in
            laerosl,            &! in
            lclouds,            &! in
            coef,               &! out
            coef_scatt_ir,      &! out
            optp,               &
            instrument)         ! in
    Else
       Call rttov_setup (         &
            setup_errorstatus,  &! out
            Err_unit,           &! in
            verbosity_level,    &! in
            nrttovid,           &! in
            laerosl,            &! in
            lclouds,            &! in
            coef,               &! out
            coef_scatt_ir,      &! out
            optp,               &
            instrument,         &! in
            ichan              ) ! in Optional 
    Endif
    
    if(any(setup_errorstatus(:) /= errorstatus_success ) ) then
       print*, 'rttov_setup fatal error'
       stop
    endif
    
    deallocate( setup_errorstatus ,stat=alloc_status(1))
    If( any(alloc_status /= 0) ) then
       errorstatus = errorstatus_fatal
       Write( errMessage, '( "mem deallocation error for setup_errorstatus")' )
       Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
       Stop
    End If
    
    DO no_id = 1, NRTTOVID
       if( any(coef(no_id)%ff_val_chn( : ) /= 1 )) then
          WRITE(*,*) ' some requested channels have bad validity parameter'
          do i = 1, nchan1(no_id)
             write(*,*) i, coef(no_id)%ff_val_chn(i)
          end do
       endif
    End Do

    DO no_id = 1, NRTTOVID

       !#########################################################
       ! Allocate memory for RTTOV_DIRECT
       !#########################################################
       allocate( rttov_errorstatus(nprof(no_id)),stat= alloc_status(3))

       ! Allocate profiles
       allocate( profiles(nprof(no_id)),stat= alloc_status(1))
       
       ! Allow profile interpolation
       interp  =  1
       
       if(interp == 0) addinterp = .false.
       if(interp == 1) addinterp = .true.
       asw = 1
       
       call rttov_alloc_prof (    &
            errorstatus,          &
            nprof(no_id),         &
            profiles,             &
            nlevels,              &
            coef_scatt_ir(no_id), &
            asw,                  &
            addclouds = lclouds,  &
            addaerosl = laerosl,  &
            init = .true.         )
        
       
       Do j = 1 , nprof(no_id)
          profiles(j) % nlevels =  nlevels
       Enddo
       
       alloc_status = 0_jpim
       ! number of channels per RTTOV call is only nchannels
       allocate(lprofiles       (nchannels(no_id)), stat = alloc_status(9))
       allocate(channels        (nchannels(no_id)), stat = alloc_status(10))
       allocate(emissivity      (nchannels(no_id)), stat = alloc_status(12))
       allocate(fresnrefl       (nchannels(no_id)), stat = alloc_status(13))
       allocate(input_emissivity(nchannels(no_id)), stat = alloc_status(14))
       allocate(calcemis        (nchannels(no_id)), stat = alloc_status(15))
       
       ! allocate transmittance arrays with number of channels
       allocate( transmission % tau_layers (profiles(1) % nlevels,nchannels(no_id) ), &
            stat= alloc_status(11))
       allocate( transmission % tau_total (nchannels(no_id) )                       , &
            stat= alloc_status(12))
       
       If( Any(alloc_status /= 0) ) Then
          errorstatus = errorstatus_fatal
          Write( errMessage, '( "allocation of transmission")' )
          Call Rttov_ErrorReport (errorstatus_fatal, errMessage, NameOfRoutine)
          !IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT',1,ZHOOK_HANDLE)
          Stop
       End If
       
       transmission % tau_layers = 0._jprb
       transmission % tau_total = 0._jprb
       
       ! Allocate radiance results arrays with number of channels
       asw = 1 ! allocate
       call rttov_alloc_rad (errorstatus,nchannels(no_id),radiance, &
            profiles(1)%nlevels,asw)
       
       AZANG  =  0
       ISUN   =  0
       IREF   =  1
       
       if(iref==0)then
          cref='NO'
          refrac=.False.
       else if(iref==1)then
          cref='YES'
          refrac=.True.
       endif
       
       if(sunzang<=87._JPRB)then
          solrad=.True.
       else
          solrad=.False.
       endif
       
       if(isun==1)then
          lsun=.true.
          if(sunzang<=87._JPRB)then
             csun='YES'
             solrad=.True.
          else
             csun='NO'
             solrad=.False.
          endif
       else
          csun='NO'
          solrad=.False.
       endif
       
       do i = 1, NPROF(no_id)
          profiles(i) %  p(:)  =  p_in(i, :)
          profiles(i) %  t(:)  =  t_in(i, :)
          profiles(i) %  q(:)  =  q_in(i, :)
          profiles(i) % o3(:)  =  o3_in(i, :)
          profiles(i) % co2(:) =  co2_in
          profiles(i) % ch4(:) =  ch4_in
          profiles(i) % n2o(:) =  n2o_in
          profiles(i) % co(:)  =  co_in
          profiles(i) % ozone_Data  =  .False.
          profiles(i) % co2_Data    =  .True.
          profiles(i) % n2o_data    =  .True.
          profiles(i) % ch4_Data    =  .True.
          profiles(i) % co_Data     =  .True.
          
          !FIXME: Make Cloud variables as passing ones if we go for all sky
          profiles(i) % cfraction  =  0.
          profiles(i) % ctp        =  500.
          profiles(i) % clw_Data   =  .False.
       
          ! 2m parameters 
          profiles(i) % s2m % p  =  p_surf(i)
          profiles(i) % s2m % t  =  t_in(i, 1)
          profiles(i) % s2m % q  =  q_in(i, 1)
          profiles(i) % s2m % u  =  2
          profiles(i) % s2m % v  =  2
          
          ! Skin variables for emissivity calculations
          profiles(i) % skin % t          =  t_skin(i)
          profiles(i) % skin % fastem(1)  =  3.0
          profiles(i) % skin % fastem(2)  =  5.0
          profiles(i) % skin % fastem(3)  =  15.0
          profiles(i) % skin % fastem(4)  =  0.1
          profiles(i) % skin % fastem(5)  =  0.3
          
          profiles(i) % zenangle      = zenang
          profiles(i) % azangle       = azang
          profiles(i) % latitude      = latitude(i)
          profiles(i) % elevation     = h_surf(i)
          profiles(i) % sunzenangle   = SUNZANG
          profiles(i) % sunazangle    = SUNAZANG
          profiles(i) % addsolar      = solrad
          profiles(i) % addrefrac     = refrac
          ! surface type
          profiles(i) % skin % surftype  = lsmask(i)
          !! FIXME: Check this one
          profiles(i) % skin % watertype = nwater(no_id)
          profiles(i) % aer_data   = laerosl
          profiles(i) % cld_data   = lclouds
          profiles(i) %idg         = 0._jprb
          profiles(i) %ish         = 0._jprb
          if( lclouds ) then
             profiles(i) %cloud(:,:)  = 0._jprb
             profiles(i) %cfrac(:,:)  = 0._jprb
          endif
       enddo
       
       ! Build the list of channels/profiles indices
       emissivity(:) = 0.0_JPRB
       channels(:) = 0_jpim
       lprofiles(:) = 0_jpim
       nch = 0_jpim
       Do j = 1 , nprof(no_id)
          DO  jch = 1,nchan1(no_id)
             nch = nch +1_jpim
             lprofiles ( nch ) = j
             if (all_channels)then
                channels( nch ) = ichan(jch,no_id)
             else 
                channels( nch ) = jch
             endif
             emissivity( nch ) = surfem(jch,no_id)
          End Do
       End Do
       
       input_emissivity(:) = emissivity(:)
       calcemis(:) = emissivity(:) < 0.01_JPRB
       
       ! FIXME: Check this one with Roger
       do j = 1 , NPROFS
          call rttov_direct(           &
               & rttov_errorstatus,    &! out
               & nprof(no_id),         &! in
               & nchannels(no_id),     &! in
               & channels,             &! in
               & lprofiles,            &! in
               & addinterp,            &! in
               & profiles,             &! in
               & coef(no_id),          &! in
               & coef_scatt_ir(no_id), &
               & optp(no_id)         , &
               & lsun,                 &! in
               & laerosl,              &! in
               & lclouds,              &! in
               & calcemis,             &! in
               & emissivity,           &! inout
               & transmission,         &! out
               & radiance  )            ! inout 
       enddo
       
       ! Initialising tbs array
       tbs(:, :)    =  0.0
       tbs_temp(:)  =  0.0
       tbs_temp     =  radiance%bt
       
       do i = 1, prof_num
          tbs(i, :)  =  tbs_temp((i-1)*nch_in+1:i*nch_in)
       enddo
       
       
       If ( any( rttov_errorstatus(:) == errorstatus_warning ) ) Then
          Do j = 1, nprof(no_id)
             If ( rttov_errorstatus(j) == errorstatus_warning ) Then
                write ( ioout, * ) 'rttov warning for profile', j
             End If
          End Do
       End If
       
       If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
          Do j = 1, nprof(no_id)
             If ( rttov_errorstatus(j) == errorstatus_fatal ) Then
                write ( ioout, * ) 'rttov error for profile',j
             End If
          End Do
          Stop
       End If
       
       ! Deallocate
       ! number of channels per RTTOV call is only nchannels
       deallocate( channels   ,stat=alloc_status(2))
       deallocate( lprofiles  ,stat=alloc_status(3))
       deallocate( emissivity ,stat=alloc_status(4))
       deallocate( fresnrefl  ,stat=alloc_status(5))
       deallocate( calcemis   ,stat=alloc_status(6))
       deallocate( input_emissivity ,stat= alloc_status(14))
       If( any(alloc_status /= 0) ) then
          errorstatus = errorstatus_fatal
          Write( errMessage, '( "mem deallocation error for channels etc")' )
          Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
          Stop
       End If
       
       asw = 0 ! deallocate radiance arrays
       call rttov_alloc_rad (errorstatus,nchan1(no_id),radiance,profiles(1) % nlevels,asw)
       If(errorstatus /= errorstatus_success) Then
          Write( errMessage, '( "deallocation error for radiances")' )
          Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
       Endif
       
       ! deallocate transmittances
       Deallocate( transmission % tau_total   ,stat= alloc_status(7))
       Deallocate( transmission % tau_layers  ,stat= alloc_status(8))
       If(errorstatus /= errorstatus_success) Then
          Write( errMessage, '( "deallocation error for transmittances")' )
          Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
          Stop
       Endif
       
       asw = 0 ! deallocate profile arrays
       call rttov_alloc_prof (errorstatus,nprof(no_id),profiles,profiles(1)%nlevels,coef_scatt_ir(no_id),asw,&
            & addclouds = lclouds, addaerosl = laerosl )
       deallocate( profiles,stat=alloc_status(1))
       If( any(alloc_status /= 0) ) then
          errorstatus = errorstatus_fatal
          Write( errMessage, '( "mem deallocation error for profiles")' )
          Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
          Stop
       End If
       
    EndDo
    
    DO no_id = 1, NRTTOVID 
       Call rttov_dealloc_coef (errorstatus, coef(no_id),coef_scatt_ir(no_id),optp(no_id))
       If(errorstatus /= errorstatus_success) Then
          Write( errMessage, '( "deallocation error for coeffs")' )
          Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
       Endif
    EndDo
    
  END SUBROUTINE RTTOV_SUBCOLUMN
END MODULE MOD_COSP_RTTOV
