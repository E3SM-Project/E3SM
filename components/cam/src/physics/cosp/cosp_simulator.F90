! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 07:08:38 -0700 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_simulator.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jan 2013 - G. Cesana - Add new variables linked to the lidar cloud phase 
!

#include "cosp_defs.h" 
MODULE MOD_COSP_SIMULATOR
  USE MOD_COSP_CONSTANTS, ONLY: I_RADAR, I_LIDAR, I_ISCCP, I_MISR, I_MODIS, &
                                I_RTTOV, I_STATS !, tsim
  USE MOD_COSP_TYPES
  USE MOD_COSP_RADAR
  USE MOD_COSP_LIDAR
  USE MOD_COSP_ISCCP_SIMULATOR
  USE MOD_COSP_MODIS_SIMULATOR
  USE MOD_COSP_MISR_SIMULATOR
#ifdef RTTOV
  USE MOD_COSP_RTTOV_SIMULATOR
#endif
  USE MOD_COSP_STATS
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef RTTOV
SUBROUTINE COSP_SIMULATOR(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
SUBROUTINE COSP_SIMULATOR(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif

  ! Arguments
  type(cosp_gridbox),intent(inout) :: gbx      ! Grid-box inputs
  type(cosp_subgrid),intent(in) :: sgx      ! Subgrid inputs
  type(cosp_sghydro),intent(in) :: sghydro  ! Subgrid info for hydrometeors
  type(cosp_config),intent(in)  :: cfg      ! Configuration options
  type(cosp_vgrid),intent(in)   :: vgrid    ! Information on vertical grid of stats
  type(cosp_sgradar),intent(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),intent(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccp),intent(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misr),intent(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modis),intent(inout)   :: modis   ! Output from MODIS simulator
#ifdef RTTOV
  type(cosp_rttov),intent(inout)    :: rttov    ! Output from RTTOV
#endif
  type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator
  ! Local variables
  ! Disable the use of timer, which causes runtime problem for pgi and pgi_acc
  ! it can be renabled if deemed necessary, by uncommenting all related to isim
  ! and tsim here and in cosp_constants.F90
  ! integer :: i,j,k,isim
  integer :: i,j,k
  logical :: inconsistent
  ! Timing variables
  integer :: t0,t1

  t0 = 0
  t1 = 0

  inconsistent=.false.
!   do k=1,gbx%Nhydro
!   do j=1,gbx%Nlevels
!   do i=1,gbx%Npoints
!     if ((gbx%mr_hydro(i,j,k)>0.0).and.(gbx%Reff(i,j,k)<=0.0)) inconsistent=.true.
!   enddo
!   enddo
!   enddo
!  if (inconsistent)  print *, '%%%% COSP_SIMULATOR: inconsistency in mr_hydro and Reff'


  !+++++++++ Radar model ++++++++++
! isim = I_RADAR
  if (cfg%Lradar_sim) then
!   call system_clock(t0)
    call cosp_radar(gbx,sgx,sghydro,sgradar)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++ Lidar model ++++++++++
! isim = I_LIDAR
  if (cfg%Llidar_sim) then
!   call system_clock(t0)
    call cosp_lidar(gbx,sgx,sghydro,sglidar)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++ ISCCP simulator ++++++++++
! isim = I_ISCCP
  if (cfg%Lisccp_sim) then
!   call system_clock(t0)
    call cosp_isccp_simulator(gbx,sgx,isccp)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++ MISR simulator ++++++++++
! isim = I_MISR
  if (cfg%Lmisr_sim) then
!   call system_clock(t0)
    call cosp_misr_simulator(gbx,sgx,misr)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++ MODIS simulator ++++++++++
! isim = I_MODIS
  if (cfg%Lmodis_sim) then
!   call system_clock(t0)
    call cosp_modis_simulator(gbx,sgx,sghydro,isccp, modis)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++ RTTOV ++++++++++ 
! isim = I_RTTOV
#ifdef RTTOV
  if (cfg%Lrttov_sim) then 
!   call system_clock(t0)
    call cosp_rttov_simulator(gbx,rttov)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif
#endif

  !+++++++++++ Summary statistics +++++++++++
! isim = I_STATS
  if (cfg%Lstats) then
!   call system_clock(t0)
    call cosp_stats(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
!   call system_clock(t1)
!   tsim(isim) = tsim(isim) + (t1 -t0)
  endif

  !+++++++++++ Change of units after computation of statistics +++++++++++
  ! This avoids using UDUNITS in CMOR

  ! Cloud fractions from 1 to %
  if (cfg%Lclcalipso) then
    where(stlidar%lidarcld /= R_UNDEF) stlidar%lidarcld = stlidar%lidarcld*100.0
  endif
  if (cfg%Lcltcalipso.OR.cfg%Lcllcalipso.OR.cfg%Lclmcalipso.OR.cfg%Lclhcalipso) then
    where(stlidar%cldlayer /= R_UNDEF) stlidar%cldlayer = stlidar%cldlayer*100.0
  endif
  if (cfg%Lclcalipso2) then
    where(stradar%lidar_only_freq_cloud /= R_UNDEF) stradar%lidar_only_freq_cloud = stradar%lidar_only_freq_cloud*100.0
  endif

  if (cfg%Lcltcalipsoliq.OR.cfg%Lcllcalipsoliq.OR.cfg%Lclmcalipsoliq.OR.cfg%Lclhcalipsoliq.OR. &
      cfg%Lcltcalipsoice.OR.cfg%Lcllcalipsoice.OR.cfg%Lclmcalipsoice.OR.cfg%Lclhcalipsoice.OR. &
      cfg%Lcltcalipsoun.OR.cfg%Lcllcalipsoun.OR.cfg%Lclmcalipsoun.OR.cfg%Lclhcalipsoun ) then
    where(stlidar%cldlayerphase /= R_UNDEF) stlidar%cldlayerphase = stlidar%cldlayerphase*100.0
  endif
  if (cfg%Lclcalipsoliq.OR.cfg%Lclcalipsoice.OR.cfg%Lclcalipsoun) then
    where(stlidar%lidarcldphase /= R_UNDEF) stlidar%lidarcldphase = stlidar%lidarcldphase*100.0
  endif
  if (cfg%Lclcalipsotmp.OR.cfg%Lclcalipsotmpliq.OR.cfg%Lclcalipsotmpice.OR.cfg%Lclcalipsotmpun) then
    where(stlidar%lidarcldtmp /= R_UNDEF) stlidar%lidarcldtmp = stlidar%lidarcldtmp*100.0
  endif

  if (cfg%Lcltisccp) then
     where(isccp%totalcldarea /= R_UNDEF) isccp%totalcldarea = isccp%totalcldarea*100.0
  endif  
  if (cfg%Lclisccp) then
    where(isccp%fq_isccp /= R_UNDEF) isccp%fq_isccp = isccp%fq_isccp*100.0
  endif

  if (cfg%LclMISR) then
    where(misr%fq_MISR /= R_UNDEF) misr%fq_MISR = misr%fq_MISR*100.0
  endif

  if (cfg%Lcltlidarradar) then
    where(stradar%radar_lidar_tcc /= R_UNDEF) stradar%radar_lidar_tcc = stradar%radar_lidar_tcc*100.0
  endif

  if (cfg%Lclmodis) then
    where(modis%Optical_Thickness_vs_Cloud_Top_Pressure /= R_UNDEF) modis%Optical_Thickness_vs_Cloud_Top_Pressure = &
                                                      modis%Optical_Thickness_vs_Cloud_Top_Pressure*100.0
  endif
  if (cfg%Lcltmodis) then
     where(modis%Cloud_Fraction_Total_Mean /= R_UNDEF) modis%Cloud_Fraction_Total_Mean = modis%Cloud_Fraction_Total_Mean*100.0
  endif
  if (cfg%Lclwmodis) then
     where(modis%Cloud_Fraction_Water_Mean /= R_UNDEF) modis%Cloud_Fraction_Water_Mean = modis%Cloud_Fraction_Water_Mean*100.0
  endif
  if (cfg%Lclimodis) then
     where(modis%Cloud_Fraction_Ice_Mean /= R_UNDEF) modis%Cloud_Fraction_Ice_Mean = modis%Cloud_Fraction_Ice_Mean*100.0
  endif

  if (cfg%Lclhmodis) then
     where(modis%Cloud_Fraction_High_Mean /= R_UNDEF) modis%Cloud_Fraction_High_Mean = modis%Cloud_Fraction_High_Mean*100.0
  endif
  if (cfg%Lclmmodis) then
     where(modis%Cloud_Fraction_Mid_Mean /= R_UNDEF) modis%Cloud_Fraction_Mid_Mean = modis%Cloud_Fraction_Mid_Mean*100.0
  endif
  if (cfg%Lcllmodis) then
     where(modis%Cloud_Fraction_Low_Mean /= R_UNDEF) modis%Cloud_Fraction_Low_Mean = modis%Cloud_Fraction_Low_Mean*100.0
  endif

  ! Change pressure from hPa to Pa.
  if (cfg%Lboxptopisccp) then
    where(isccp%boxptop /= R_UNDEF) isccp%boxptop = isccp%boxptop*100.0
  endif
  if (cfg%Lpctisccp) then
    where(isccp%meanptop /= R_UNDEF) isccp%meanptop = isccp%meanptop*100.0
  endif


END SUBROUTINE COSP_SIMULATOR

END MODULE MOD_COSP_SIMULATOR
