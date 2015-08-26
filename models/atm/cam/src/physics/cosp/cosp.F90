! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
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

#include "cosp_defs.h"
MODULE MOD_COSP
  USE MOD_COSP_TYPES
  USE MOD_COSP_SIMULATOR
  USE MOD_COSP_MODIS_SIMULATOR
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP ---------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef RTTOV
SUBROUTINE COSP(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
SUBROUTINE COSP(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif
  ! Arguments
  integer,intent(in) :: overlap !  overlap type in SCOPS: 1=max, 2=rand, 3=max/rand
  integer,intent(in) :: Ncolumns
  type(cosp_config),intent(in) :: cfg   ! Configuration options
  type(cosp_vgrid),intent(in) :: vgrid   ! Information on vertical grid of stats
  type(cosp_gridbox),intent(inout) :: gbx
  type(cosp_subgrid),intent(inout),target :: sgx   ! Subgrid info  !!++added JEK per jet 5-11-2010
  !!type(cosp_subgrid),intent(inout) :: sgx   ! Subgrid info
  type(cosp_sgradar),intent(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),intent(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccp),intent(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misr),intent(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modis),intent(inout)   :: modis   ! Output from MODIS simulator
#ifdef RTTOV
  type(cosp_rttov),intent(inout)   :: rttov   ! Output from RTTOV
#endif
  type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator

  ! Local variables 
  integer :: Npoints   ! Number of gridpoints
  integer :: Nlevels   ! Number of levels
  integer :: Nhydro    ! Number of hydrometeors
  integer :: Niter     ! Number of calls to cosp_simulator
  integer :: i_first,i_last ! First and last gridbox to be processed in each iteration
  integer :: i,Ni
  integer,dimension(2) :: ix,iy
  logical :: reff_zero
  real :: maxp,minp
  integer,dimension(:),allocatable :: & ! Dimensions nPoints
                  seed    !  It is recommended that the seed is set to a different value for each model
                          !  gridbox it is called on, as it is possible that the choice of the same 
                          !  seed value every time may introduce some statistical bias in the results, 
                          !  particularly for low values of NCOL.
  ! Types used in one iteration
  type(cosp_gridbox) :: gbx_it
  type(cosp_subgrid) :: sgx_it
  type(cosp_vgrid)   :: vgrid_it
  type(cosp_sgradar) :: sgradar_it
  type(cosp_sglidar) :: sglidar_it
  type(cosp_isccp)   :: isccp_it
  type(cosp_modis)   :: modis_it
  type(cosp_misr)    :: misr_it
#ifdef RTTOV
  type(cosp_rttov)   :: rttov_it
#endif
  type(cosp_radarstats) :: stradar_it
  type(cosp_lidarstats) :: stlidar_it
  
  !++++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

!++++++++++ Apply sanity checks to inputs ++++++++++
  call cosp_check_input('longitude',gbx%longitude,min_val=0.0,max_val=360.0)
  call cosp_check_input('latitude',gbx%latitude,min_val=-90.0,max_val=90.0)
  call cosp_check_input('dlev',gbx%dlev,min_val=0.0)
  call cosp_check_input('p',gbx%p,min_val=0.0)
  call cosp_check_input('ph',gbx%ph,min_val=0.0)
  call cosp_check_input('T',gbx%T,min_val=0.0)
  call cosp_check_input('q',gbx%q,min_val=0.0)
  call cosp_check_input('sh',gbx%sh,min_val=0.0)
  call cosp_check_input('dtau_s',gbx%dtau_s,min_val=0.0)
  call cosp_check_input('dtau_c',gbx%dtau_c,min_val=0.0)
  call cosp_check_input('dtau_s_snow',gbx%dtau_s_snow,min_val=0.0)
  call cosp_check_input('dem_s',gbx%dem_s,min_val=0.0,max_val=1.0)
  call cosp_check_input('dem_c',gbx%dem_c,min_val=0.0,max_val=1.0)
  call cosp_check_input('dem_s_snow',gbx%dem_s_snow,min_val=0.0) !!+jek, removed max_val=1.0 per steve klein
  ! Point information (Npoints)
  call cosp_check_input('land',gbx%land,min_val=0.0,max_val=1.0)
  call cosp_check_input('psfc',gbx%psfc,min_val=0.0)
  call cosp_check_input('sunlit',gbx%sunlit,min_val=0.0,max_val=1.0)
  call cosp_check_input('skt',gbx%skt,min_val=0.0)
  ! TOTAL and CONV cloud fraction for SCOPS
  call cosp_check_input('tca',gbx%tca,min_val=0.0,max_val=1.0)
  call cosp_check_input('cca',gbx%cca,min_val=0.0,max_val=1.0)
  ! Precipitation fluxes on model levels
  call cosp_check_input('rain_ls',gbx%rain_ls,min_val=0.0)
  call cosp_check_input('rain_cv',gbx%rain_cv,min_val=0.0)
  call cosp_check_input('snow_ls',gbx%snow_ls,min_val=0.0)
  call cosp_check_input('snow_cv',gbx%snow_cv,min_val=0.0)
  call cosp_check_input('grpl_ls',gbx%grpl_ls,min_val=0.0)
  ! Hydrometeors concentration and distribution parameters
  call cosp_check_input('mr_hydro',gbx%mr_hydro,min_val=0.0)
  ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
  call cosp_check_input('Reff',gbx%Reff,min_val=0.0)
  reff_zero=.true.
  if (any(gbx%Reff > 1.e-8)) then
     reff_zero=.false.
      ! reff_zero == .false.
      !     and gbx%use_reff == .true.   Reff use in radar and lidar
      !     and reff_zero    == .false.  Reff use in lidar and set to 0 for radar
  endif
  if ((.not. gbx%use_reff) .and. (reff_zero)) then ! No Reff in radar. Default in lidar
        gbx%Reff = DEFAULT_LIDAR_REFF
        print *, '---------- COSP WARNING ------------'
        print *, ''
        print *, 'Using default Reff in lidar simulations'
        print *, ''
        print *, '----------------------------------'
  endif
  
  ! Aerosols concentration and distribution parameters
  call cosp_check_input('conc_aero',gbx%conc_aero,min_val=0.0)
  ! Checks for CRM mode
  if (Ncolumns == 1) then
     if (gbx%use_precipitation_fluxes) then
        print *, '---------- COSP ERROR ------------'
        print *, ''
        print *, 'Use of precipitation fluxes not supported in CRM mode (Ncolumns=1)'
        print *, ''
        print *, '----------------------------------'
        stop
     endif
     if ((maxval(gbx%dtau_c) > 0.0).or.(maxval(gbx%dem_c) > 0.0)) then
        print *, '---------- COSP ERROR ------------'
        print *, ''
        print *, ' dtau_c > 0.0 or dem_c > 0.0. In CRM mode (Ncolumns=1), '
        print *, ' the optical depth (emmisivity) of all clouds must be '
        print *, ' passed through dtau_s (dem_s)'
        print *, ''
        print *, '----------------------------------'
        stop
     endif
  endif

   ! We base the seed in the decimal part of the surface pressure.
   allocate(seed(Npoints))
   seed = int(gbx%psfc) ! This is to avoid division by zero when Npoints = 1   
      ! Roj Oct/2008 ... Note: seed value of 0 caused me some problems + I want to 
      ! randomize for each call to COSP even when Npoints ==1
   minp = minval(gbx%psfc)
   maxp = maxval(gbx%psfc)
   !!if (Npoints .gt. 1) seed=int((gbx%psfc-minp)/(maxp-minp)*100000) + 1
   !! change per Yuying Zhang to make off-line and in-line COSP v1.3 results consistent.
   if (Npoints .gt. 1) seed=(gbx%psfc-int(gbx%psfc))*1000000
  
   if (gbx%Npoints_it >= gbx%Npoints) then ! One iteration gbx%Npoints
#ifdef RTTOV
        call cosp_iter(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
        call cosp_iter(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif
   else ! Several iterations to save memory
        Niter = gbx%Npoints/gbx%Npoints_it ! Integer division
        if (Niter*gbx%Npoints_it < gbx%Npoints) Niter = Niter + 1
        do i=1,Niter
            i_first = (i-1)*gbx%Npoints_it + 1
            i_last  = i_first + gbx%Npoints_it - 1
            i_last  = min(i_last,gbx%Npoints)
            Ni = i_last - i_first + 1
            if (i == 1) then
                ! Allocate types for all but last iteration
                call construct_cosp_gridbox(gbx%time,gbx%time_bnds,gbx%radar_freq,gbx%surface_radar,gbx%use_mie_tables, &
                                            gbx%use_gas_abs,gbx%do_ray,gbx%melt_lay,gbx%k2,Ni,Nlevels, &
                                            Ncolumns,N_HYDRO,gbx%Nprmts_max_hydro, &
                                            gbx%Naero,gbx%Nprmts_max_aero,Ni,gbx%lidar_ice_type,gbx%isccp_top_height, &
                                            gbx%isccp_top_height_direction,gbx%isccp_overlap,gbx%isccp_emsfc_lw, &
                                            gbx%use_precipitation_fluxes,gbx%use_reff, &
                                            gbx%plat,gbx%sat,gbx%inst,gbx%nchan,gbx%ZenAng, &
                                            gbx%Ichan(1:gbx%nchan),gbx%surfem(1:gbx%nchan), &
                                            gbx%co2,gbx%ch4,gbx%n2o,gbx%co, &
                                            gbx_it)
                call construct_cosp_vgrid(gbx_it,vgrid%Nlvgrid,vgrid%use_vgrid,vgrid%csat_vgrid,vgrid_it)
                call construct_cosp_subgrid(Ni, Ncolumns, Nlevels, sgx_it)
                call construct_cosp_sgradar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,sgradar_it)
                call construct_cosp_sglidar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar_it)
                call construct_cosp_isccp(cfg,Ni,Ncolumns,Nlevels,isccp_it)
                call construct_cosp_modis(cfg, Ni, modis_it)
                call construct_cosp_misr(cfg,Ni,misr_it)
#ifdef RTTOV
                call construct_cosp_rttov(Ni,gbx%nchan,rttov_it)
#endif
                call construct_cosp_radarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar_it)
                call construct_cosp_lidarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar_it)
            elseif (i == Niter) then ! last iteration
                call free_cosp_gridbox(gbx_it,.true.)
                call free_cosp_subgrid(sgx_it)
                call free_cosp_vgrid(vgrid_it)
                call free_cosp_sgradar(sgradar_it)
                call free_cosp_sglidar(sglidar_it)
                call free_cosp_isccp(isccp_it)
                call free_cosp_modis(modis_it)
                call free_cosp_misr(misr_it)
#ifdef RTTOV
                call free_cosp_rttov(rttov_it)
#endif
                call free_cosp_radarstats(stradar_it)
                call free_cosp_lidarstats(stlidar_it)
                ! Allocate types for iterations
                call construct_cosp_gridbox(gbx%time,gbx%time_bnds,gbx%radar_freq,gbx%surface_radar,gbx%use_mie_tables, &
                                            gbx%use_gas_abs,gbx%do_ray,gbx%melt_lay,gbx%k2,Ni,Nlevels, &
                                            Ncolumns,N_HYDRO,gbx%Nprmts_max_hydro, &
                                            gbx%Naero,gbx%Nprmts_max_aero,Ni,gbx%lidar_ice_type,gbx%isccp_top_height, &
                                            gbx%isccp_top_height_direction,gbx%isccp_overlap,gbx%isccp_emsfc_lw, &
                                            gbx%use_precipitation_fluxes,gbx%use_reff, &
                                            gbx%plat,gbx%sat,gbx%inst,gbx%nchan,gbx%ZenAng, &
                                            gbx%Ichan(1:gbx%nchan),gbx%surfem(1:gbx%nchan), &
                                            gbx%co2,gbx%ch4,gbx%n2o,gbx%co, &
                                            gbx_it)
                ! --- Copy arrays without Npoints as dimension ---
                gbx_it%dist_prmts_hydro = gbx%dist_prmts_hydro
                gbx_it%dist_type_aero   = gbx_it%dist_type_aero
                call construct_cosp_vgrid(gbx_it,vgrid%Nlvgrid,vgrid%use_vgrid,vgrid%csat_vgrid,vgrid_it)
                call construct_cosp_subgrid(Ni, Ncolumns, Nlevels, sgx_it)
                call construct_cosp_sgradar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,sgradar_it)
                call construct_cosp_sglidar(cfg,Ni,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar_it)
                call construct_cosp_isccp(cfg,Ni,Ncolumns,Nlevels,isccp_it)
                call construct_cosp_modis(cfg,Ni, modis_it)
                call construct_cosp_misr(cfg,Ni,misr_it)
#ifdef RTTOV 
                call construct_cosp_rttov(Ni,gbx%nchan,rttov_it) 
#endif 
                call construct_cosp_radarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar_it)
                call construct_cosp_lidarstats(cfg,Ni,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar_it)
            endif
            ! --- Copy sections of arrays with Npoints as dimension ---
            ix=(/i_first,i_last/)
            iy=(/1,Ni/)
            call cosp_gridbox_cpsection(ix,iy,gbx,gbx_it)
              ! These serve as initialisation of *_it types
            call cosp_subgrid_cpsection(ix,iy,sgx,sgx_it)
            if (cfg%Lradar_sim) call cosp_sgradar_cpsection(ix,iy,sgradar,sgradar_it)
            if (cfg%Llidar_sim) call cosp_sglidar_cpsection(ix,iy,sglidar,sglidar_it)
            if (cfg%Lisccp_sim) call cosp_isccp_cpsection(ix,iy,isccp,isccp_it)
            if (cfg%Lmodis_sim) call cosp_modis_cpsection(ix,iy,modis,modis_it)
            if (cfg%Lmisr_sim)  call cosp_misr_cpsection(ix,iy,misr,misr_it)
#ifdef RTTOV 
            if (cfg%Lrttov_sim) call cosp_rttov_cpsection(ix,iy,rttov,rttov_it) 
#endif
            if (cfg%Lradar_sim) call cosp_radarstats_cpsection(ix,iy,stradar,stradar_it)
            if (cfg%Llidar_sim) call cosp_lidarstats_cpsection(ix,iy,stlidar,stlidar_it)
            print *,'---------ix: ',ix
#ifdef RTTOV
            call cosp_iter(overlap,seed(ix(1):ix(2)),cfg,vgrid_it,gbx_it,sgx_it,sgradar_it, &
                           sglidar_it,isccp_it,misr_it,modis_it,rttov_it,stradar_it,stlidar_it)
#else
            call cosp_iter(overlap,seed(ix(1):ix(2)),cfg,vgrid_it,gbx_it,sgx_it,sgradar_it, &
                           sglidar_it,isccp_it,misr_it,modis_it,stradar_it,stlidar_it)
#endif
            ! --- Copy results to output structures ---
            ix=(/1,Ni/)
            iy=(/i_first,i_last/)
            call cosp_subgrid_cpsection(ix,iy,sgx_it,sgx)
            if (cfg%Lradar_sim) call cosp_sgradar_cpsection(ix,iy,sgradar_it,sgradar)
            if (cfg%Llidar_sim) call cosp_sglidar_cpsection(ix,iy,sglidar_it,sglidar)
            if (cfg%Lisccp_sim) call cosp_isccp_cpsection(ix,iy,isccp_it,isccp)
            if (cfg%Lmodis_sim) call cosp_modis_cpsection(ix,iy,modis_it,modis)
            if (cfg%Lmisr_sim)  call cosp_misr_cpsection(ix,iy,misr_it,misr)
#ifdef RTTOV 
            if (cfg%Lrttov_sim) call cosp_rttov_cpsection(ix,iy,rttov_it,rttov) 
#endif 
            if (cfg%Lradar_sim) call cosp_radarstats_cpsection(ix,iy,stradar_it,stradar)
            if (cfg%Llidar_sim) call cosp_lidarstats_cpsection(ix,iy,stlidar_it,stlidar)
        enddo
        ! Deallocate types
        call free_cosp_gridbox(gbx_it,.true.)
        call free_cosp_subgrid(sgx_it)
        call free_cosp_vgrid(vgrid_it)
        call free_cosp_sgradar(sgradar_it)
        call free_cosp_sglidar(sglidar_it)
        call free_cosp_isccp(isccp_it)
        call free_cosp_modis(modis_it)
        call free_cosp_misr(misr_it)
#ifdef RTTOV 
        call free_cosp_rttov(rttov_it) 
#endif
        call free_cosp_radarstats(stradar_it)
        call free_cosp_lidarstats(stlidar_it)
   endif
   deallocate(seed)

    
END SUBROUTINE COSP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_ITER ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef RTTOV
SUBROUTINE COSP_ITER(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
SUBROUTINE COSP_ITER(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif
  ! Arguments
  integer,intent(in) :: overlap !  overlap type in SCOPS: 1=max, 2=rand, 3=max/rand
  integer,dimension(:),intent(inout) :: seed
  type(cosp_config),intent(in) :: cfg   ! Configuration options
  type(cosp_vgrid),intent(in) :: vgrid   ! Information on vertical grid of stats
  type(cosp_gridbox),intent(inout) :: gbx
  type(cosp_subgrid),intent(inout) :: sgx   ! Subgrid info
  type(cosp_sgradar),intent(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),intent(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccp),intent(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misr),intent(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modis),intent(inout)   :: modis   ! Output from MODIS simulator
#ifdef RTTOV
  type(cosp_rttov),intent(inout)   :: rttov   ! Output from RTTOV
#endif
  type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator

  ! Local variables 
  integer :: Npoints   ! Number of gridpoints
  integer :: Ncolumns  ! Number of subcolumns
  integer :: Nlevels   ! Number of levels
  integer :: Nhydro    ! Number of hydrometeors
  integer :: i,j,k
  integer :: I_HYDRO 
  real,dimension(:,:),pointer :: column_frac_out ! Array with one column of frac_out
  real,dimension(:,:),pointer :: column_prec_out ! Array with one column of prec_frac !modified by ZYY
  integer :: scops_debug=0    !  set to non-zero value to print out inputs for debugging in SCOPS
  real,dimension(:, :),allocatable :: cca_scops,ls_p_rate,cv_p_rate, &
                     tca_scops ! Cloud cover in each model level (HORIZONTAL gridbox fraction) of total cloud.
                               ! Levels are from TOA to SURFACE. (nPoints, nLev)
  real,dimension(:,:),allocatable :: frac_ls,prec_ls,frac_cv,prec_cv ! Cloud/Precipitation fraction in each model level
                                                                     ! Levels are from SURFACE to TOA
  real,dimension(:,:),allocatable :: rho ! (Npoints, Nlevels). Atmospheric density
  type(cosp_sghydro) :: sghydro   ! Subgrid info for hydrometeors en each iteration

  
  !++++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Ncolumns = gbx%Ncolumns
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro
    
  !++++++++++ Climate/NWP mode ++++++++++  
  if (Ncolumns > 1) then
        !++++++++++ Subgrid sampling ++++++++++
        ! Allocate arrays before calling SCOPS
        allocate(frac_ls(Npoints,Nlevels),frac_cv(Npoints,Nlevels),prec_ls(Npoints,Nlevels),prec_cv(Npoints,Nlevels))
        allocate(tca_scops(Npoints,Nlevels),cca_scops(Npoints,Nlevels), &
                ls_p_rate(Npoints,Nlevels),cv_p_rate(Npoints,Nlevels))
        ! Initialize to zero
        frac_ls=0.0
        prec_ls=0.0
        frac_cv=0.0
        prec_cv=0.0
        ! Cloud fractions for SCOPS from TOA to SFC
        tca_scops = gbx%tca(:,Nlevels:1:-1)
        cca_scops = gbx%cca(:,Nlevels:1:-1)
        
        ! Call to SCOPS
        ! strat and conv arrays are passed with levels from TOA to SURFACE.
        call scops(Npoints,Nlevels,Ncolumns,seed,tca_scops,cca_scops,overlap,sgx%frac_out,scops_debug)
        
        ! temporarily use prec_ls/cv to transfer information about precipitation flux into prec_scops
        if(gbx%use_precipitation_fluxes) then
            ls_p_rate(:,Nlevels:1:-1)=gbx%rain_ls(:,1:Nlevels)+gbx%snow_ls(:,1:Nlevels)+gbx%grpl_ls(:,1:Nlevels)
            cv_p_rate(:,Nlevels:1:-1)=gbx%rain_cv(:,1:Nlevels)+gbx%snow_cv(:,1:Nlevels)
        else
            ls_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_LSRAIN)+ &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSSNOW)+ &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSGRPL)
            cv_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_CVRAIN)+ &
                                      gbx%mr_hydro(:,1:Nlevels,I_CVSNOW)
        endif
        
        call prec_scops(Npoints,Nlevels,Ncolumns,ls_p_rate,cv_p_rate,sgx%frac_out,sgx%prec_frac)
        
        ! Precipitation fraction
        do j=1,Npoints,1
        do k=1,Nlevels,1
            do i=1,Ncolumns,1
                if (sgx%frac_out (j,i,Nlevels+1-k) == I_LSC) frac_ls(j,k)=frac_ls(j,k)+1.
                if (sgx%frac_out (j,i,Nlevels+1-k) == I_CVC) frac_cv(j,k)=frac_cv(j,k)+1.
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 1) prec_ls(j,k)=prec_ls(j,k)+1.
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 2) prec_cv(j,k)=prec_cv(j,k)+1.
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 3) then
                    prec_cv(j,k)=prec_cv(j,k)+1.
                    prec_ls(j,k)=prec_ls(j,k)+1.
                endif
            enddo  !i
            frac_ls(j,k)=frac_ls(j,k)/Ncolumns
            frac_cv(j,k)=frac_cv(j,k)/Ncolumns
            prec_ls(j,k)=prec_ls(j,k)/Ncolumns
            prec_cv(j,k)=prec_cv(j,k)/Ncolumns
        enddo  !k
        enddo  !j
        
        !adjust grid-box mean snow properties to local properties
        !convert longwave optical depth to longwave emissivity
        !n.b.: both prec_ls and gbx% are ordered from bottom to top at this point (I believe!?!)
        do j=1,Npoints,1
        do k=1,Nlevels,1
            if (prec_ls(j,k) .ne. 0. .and. gbx%dtau_s_snow(j,k) .gt. 0.) then
                   gbx%dtau_s_snow(j,k) = gbx%dtau_s_snow(j,k)/prec_ls(j,k)
            end if
            if (prec_ls(j,k) .ne. 0. .and. gbx%dem_s_snow(j,k) .gt. 0.) then
                   gbx%dem_s_snow(j,k) = gbx%dem_s_snow(j,k)/prec_ls(j,k)
                   gbx%dem_s_snow(j,k) = 1. - exp ( -1. * gbx%dem_s_snow(j,k))
            end if
        enddo !k
        enddo !j

         ! Levels from SURFACE to TOA.
        if (Npoints*Ncolumns*Nlevels < 10000) then
            sgx%frac_out(1:Npoints,:,1:Nlevels)  = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
            sgx%prec_frac(1:Npoints,:,1:Nlevels) = sgx%prec_frac(1:Npoints,:,Nlevels:1:-1)
        else
            ! This is done within a loop (unvectorized) over nPoints to save memory
            do j=1,Npoints
                sgx%frac_out(j,:,1:Nlevels)  = sgx%frac_out(j,:,Nlevels:1:-1)
                sgx%prec_frac(j,:,1:Nlevels) = sgx%prec_frac(j,:,Nlevels:1:-1)
            enddo
        endif
       
       ! Deallocate arrays that will no longer be used
        deallocate(tca_scops,cca_scops,ls_p_rate,cv_p_rate)

        ! Populate the subgrid arrays
        call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)
        do k=1,Ncolumns
            !--------- Mixing ratios for clouds and Reff for Clouds and precip -------
            column_frac_out => sgx%frac_out(:,k,:)
            where (column_frac_out == I_LSC)     !+++++++++++ LS clouds ++++++++
                sghydro%mr_hydro(:,k,:,I_LSCLIQ) = gbx%mr_hydro(:,:,I_LSCLIQ)
                sghydro%mr_hydro(:,k,:,I_LSCICE) = gbx%mr_hydro(:,:,I_LSCICE)
                
                sghydro%Reff(:,k,:,I_LSCLIQ)     = gbx%Reff(:,:,I_LSCLIQ)
                sghydro%Reff(:,k,:,I_LSCICE)     = gbx%Reff(:,:,I_LSCICE)
              ! sghydro%Reff(:,k,:,I_LSRAIN)     = gbx%Reff(:,:,I_LSRAIN)
              ! sghydro%Reff(:,k,:,I_LSSNOW)     = gbx%Reff(:,:,I_LSSNOW)
              ! sghydro%Reff(:,k,:,I_LSGRPL)     = gbx%Reff(:,:,I_LSGRPL)
            elsewhere (column_frac_out == I_CVC) !+++++++++++ CONV clouds ++++++++
                sghydro%mr_hydro(:,k,:,I_CVCLIQ) = gbx%mr_hydro(:,:,I_CVCLIQ) 
                sghydro%mr_hydro(:,k,:,I_CVCICE) = gbx%mr_hydro(:,:,I_CVCICE) 
                
                sghydro%Reff(:,k,:,I_CVCLIQ)     = gbx%Reff(:,:,I_CVCLIQ) 
                sghydro%Reff(:,k,:,I_CVCICE)     = gbx%Reff(:,:,I_CVCICE) 
              ! sghydro%Reff(:,k,:,I_CVRAIN)     = gbx%Reff(:,:,I_CVRAIN) 
              ! sghydro%Reff(:,k,:,I_CVSNOW)     = gbx%Reff(:,:,I_CVSNOW) 
            end where
                        !modified by ZYY on Apr
            column_prec_out => sgx%prec_frac(:,k,:)

            where ((column_prec_out == 1) .or. (column_prec_out == 3) )  !++++ LS precip ++++
                sghydro%Reff(:,k,:,I_LSRAIN)     = gbx%Reff(:,:,I_LSRAIN)
                sghydro%Reff(:,k,:,I_LSSNOW)     = gbx%Reff(:,:,I_LSSNOW)
                sghydro%Reff(:,k,:,I_LSGRPL)     = gbx%Reff(:,:,I_LSGRPL)
            elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3)) !++++ CONV precip ++++
                sghydro%Reff(:,k,:,I_CVRAIN)     = gbx%Reff(:,:,I_CVRAIN) 
                sghydro%Reff(:,k,:,I_CVSNOW)     = gbx%Reff(:,:,I_CVSNOW) 
            end where
            !end of modification of ZYY on Apr

            !--------- Precip -------
            if (.not. gbx%use_precipitation_fluxes) then
                where (column_frac_out == I_LSC)  !+++++++++++ LS Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_LSRAIN) = gbx%mr_hydro(:,:,I_LSRAIN)
                    sghydro%mr_hydro(:,k,:,I_LSSNOW) = gbx%mr_hydro(:,:,I_LSSNOW)
                    sghydro%mr_hydro(:,k,:,I_LSGRPL) = gbx%mr_hydro(:,:,I_LSGRPL)
                elsewhere (column_frac_out == I_CVC) !+++++++++++ CONV Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_CVRAIN) = gbx%mr_hydro(:,:,I_CVRAIN) 
                    sghydro%mr_hydro(:,k,:,I_CVSNOW) = gbx%mr_hydro(:,:,I_CVSNOW) 
                end where 
            endif
        enddo
        ! convert the mixing ratio and precipitation flux from gridbox mean to the fraction-based values
        do k=1,Nlevels
            do j=1,Npoints
                !--------- Clouds -------
                if (frac_ls(j,k) .ne. 0.) then
                    sghydro%mr_hydro(j,:,k,I_LSCLIQ) = sghydro%mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                    sghydro%mr_hydro(j,:,k,I_LSCICE) = sghydro%mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
                endif
                if (frac_cv(j,k) .ne. 0.) then
                    sghydro%mr_hydro(j,:,k,I_CVCLIQ) = sghydro%mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                    sghydro%mr_hydro(j,:,k,I_CVCICE) = sghydro%mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
                endif
                !--------- Precip -------
                if (gbx%use_precipitation_fluxes) then
                    if (prec_ls(j,k) .ne. 0.) then
                        gbx%rain_ls(j,k) = gbx%rain_ls(j,k)/prec_ls(j,k)
                        gbx%snow_ls(j,k) = gbx%snow_ls(j,k)/prec_ls(j,k)
                        gbx%grpl_ls(j,k) = gbx%grpl_ls(j,k)/prec_ls(j,k)
                    endif
                    if (prec_cv(j,k) .ne. 0.) then
                        gbx%rain_cv(j,k) = gbx%rain_cv(j,k)/prec_cv(j,k)
                        gbx%snow_cv(j,k) = gbx%snow_cv(j,k)/prec_cv(j,k)
                    endif
                else
                    if (prec_ls(j,k) .ne. 0.) then
                        sghydro%mr_hydro(j,:,k,I_LSRAIN) = sghydro%mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSSNOW) = sghydro%mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSGRPL) = sghydro%mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                    endif
                    if (prec_cv(j,k) .ne. 0.) then
                        sghydro%mr_hydro(j,:,k,I_CVRAIN) = sghydro%mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                        sghydro%mr_hydro(j,:,k,I_CVSNOW) = sghydro%mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                    endif
                endif  
            enddo !k
        enddo !j
        deallocate(frac_ls,prec_ls,frac_cv,prec_cv)
        
        if (gbx%use_precipitation_fluxes) then
            ! Density
            allocate(rho(Npoints,Nlevels))
            I_HYDRO = I_LSRAIN
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1., &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                    gbx%rain_ls,sghydro%mr_hydro(:,:,:,I_HYDRO))
            I_HYDRO = I_LSSNOW
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1., &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                    gbx%snow_ls,sghydro%mr_hydro(:,:,:,I_HYDRO))
            I_HYDRO = I_CVRAIN
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,2., &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                    gbx%rain_cv,sghydro%mr_hydro(:,:,:,I_HYDRO))
            I_HYDRO = I_CVSNOW
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,2., &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                    gbx%snow_cv,sghydro%mr_hydro(:,:,:,I_HYDRO))
            I_HYDRO = I_LSGRPL
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1., &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                    gbx%grpl_ls,sghydro%mr_hydro(:,:,:,I_HYDRO))
            if(allocated(rho)) deallocate(rho)
        endif
   !++++++++++ CRM mode ++++++++++
   else
      call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)
      sghydro%mr_hydro(:,1,:,:) = gbx%mr_hydro
      sghydro%Reff(:,1,:,:) = gbx%Reff
      !--------- Clouds -------
      where ((gbx%dtau_s > 0.0))
             sgx%frac_out(:,1,:) = 1  ! Subgrid cloud array. Dimensions (Npoints,Ncolumns,Nlevels)
      endwhere
   endif ! Ncolumns > 1
  
   !++++++++++ Simulator ++++++++++
#ifdef RTTOV
    call cosp_simulator(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
    call cosp_simulator(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif

    ! Deallocate subgrid arrays
    call free_cosp_sghydro(sghydro)
END SUBROUTINE COSP_ITER

END MODULE MOD_COSP
