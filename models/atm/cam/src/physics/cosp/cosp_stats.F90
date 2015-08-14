! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 07:08:38 -0700 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_stats.F90 $
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
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
! Jun 2010 - T. Yokohata, T. Nishimura and K. Ogochi - Added NEC SXs optimisations
! Jan 2013 - G. Cesana        - Added betaperp and temperature arguments 
!                             - Added phase 3D/3Dtemperature/Map output variables in diag_lidar 
!
!
#include "cosp_defs.h" 
MODULE MOD_COSP_STATS
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  USE MOD_LLNL_STATS
  USE MOD_LMD_IPSL_STATS
  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_STATS ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_STATS(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)

   ! Input arguments
   type(cosp_gridbox),intent(in) :: gbx
   type(cosp_subgrid),intent(in) :: sgx
   type(cosp_config),intent(in)  :: cfg
   type(cosp_sgradar),intent(in) :: sgradar
   type(cosp_sglidar),intent(in) :: sglidar
   type(cosp_vgrid),intent(in)   :: vgrid
   ! Output arguments
   type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics for radar
   type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics for lidar

   ! Local variables
   integer :: Npoints  !# of grid points
   integer :: Nlevels  !# of levels
   integer :: Nhydro   !# of hydrometeors
   integer :: Ncolumns !# of columns
   integer :: Nlr
   logical :: ok_lidar_cfad = .false.
   real,dimension(:,:,:),allocatable :: Ze_out,betatot_out,betamol_in,betamol_out,ph_in,ph_out
   real,dimension(:,:),allocatable :: ph_c,betamol_c
   real,dimension(:,:,:),allocatable ::  betaperptot_out, temp_in, temp_out 
   real,dimension(:,:),allocatable :: temp_c

   Npoints  = gbx%Npoints
   Nlevels  = gbx%Nlevels
   Nhydro   = gbx%Nhydro
   Ncolumns = gbx%Ncolumns
   Nlr      = vgrid%Nlvgrid

   if (cfg%LcfadLidarsr532) ok_lidar_cfad=.true.

   if (vgrid%use_vgrid) then ! Statistics in a different vertical grid
        allocate(Ze_out(Npoints,Ncolumns,Nlr),betatot_out(Npoints,Ncolumns,Nlr), &
                 betamol_in(Npoints,1,Nlevels),betamol_out(Npoints,1,Nlr),betamol_c(Npoints,Nlr), &
                 ph_in(Npoints,1,Nlevels),ph_out(Npoints,1,Nlr),ph_c(Npoints,Nlr))
        Ze_out = 0.0
        betatot_out  = 0.0
        betamol_out= 0.0
        betamol_c  = 0.0
        ph_in(:,1,:)  = gbx%ph(:,:)
        ph_out  = 0.0
        ph_c    = 0.0
        allocate(betaperptot_out(Npoints,Ncolumns,Nlr),temp_in(Npoints,1,Nlevels),temp_out(Npoints,1,Nlr), &
                 temp_c(Npoints,Nlr))
        betaperptot_out = 0.0
        temp_in = 0.0
        temp_out = 0.0
        temp_c = 0.0

        !++++++++++++ Radar CFAD ++++++++++++++++
        if (cfg%Lradar_sim) then
            call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,gbx%zlev_half,sgradar%Ze_tot, &
                                           Nlr,vgrid%zl,vgrid%zu,Ze_out,log_units=.true.)
            stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,Ze_out, &
                                        DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH)
        endif
        !++++++++++++ Lidar CFAD ++++++++++++++++
        if (cfg%Llidar_sim) then
            betamol_in(:,1,:) = sglidar%beta_mol(:,:)
            call cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,gbx%zlev_half,betamol_in, &
                                           Nlr,vgrid%zl,vgrid%zu,betamol_out)
            call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,gbx%zlev_half,sglidar%beta_tot, &
                                           Nlr,vgrid%zl,vgrid%zu,betatot_out)

            temp_in(:,1,:) = gbx%T(:,:)
            call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,gbx%zlev_half,sglidar%betaperp_tot, &
                                           Nlr,vgrid%zl,vgrid%zu,betaperptot_out)
            call cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,gbx%zlev_half,temp_in, &
                                           Nlr,vgrid%zl,vgrid%zu,temp_out)
            temp_c(:,:) = temp_out(:,1,:)

            call cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,gbx%zlev_half,ph_in, &
                                           Nlr,vgrid%zl,vgrid%zu,ph_out)
            ph_c(:,:) = ph_out(:,1,:)
            betamol_c(:,:) = betamol_out(:,1,:)
            ! Stats from lidar_stat_summary
            call diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                            ,temp_c,betatot_out,betaperptot_out,betamol_c,sglidar%refl,gbx%land,ph_c &
                            ,LIDAR_UNDEF,ok_lidar_cfad &
                            ,stlidar%cfad_sr,stlidar%srbval &
                            ,LIDAR_NCAT,stlidar%lidarcld,stlidar%lidarcldphase &
                            ,stlidar%cldlayer,stlidar%cldlayerphase,stlidar%lidarcldtmp &
                            ,stlidar%parasolrefl)
        endif

        !++++++++++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++++++++++
        if (cfg%Lradar_sim.and.cfg%Llidar_sim) call cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr, &
                                    temp_c,betatot_out,betaperptot_out,betamol_c,Ze_out, &
                                    stradar%lidar_only_freq_cloud,stradar%radar_lidar_tcc, &
				    stradar%radar_tcc,stradar%radar_tcc_2) !+JEK
        deallocate(temp_out,temp_c,betaperptot_out)

        ! Deallocate arrays at coarse resolution
        deallocate(Ze_out,betatot_out,betamol_in,betamol_out,betamol_c,ph_in,ph_out,ph_c)
   else ! Statistics in model levels
        !++++++++++++ Radar CFAD ++++++++++++++++
        if (cfg%Lradar_sim) stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,sgradar%Ze_tot, &
                                        DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH)
        !++++++++++++ Lidar CFAD ++++++++++++++++
        ! Stats from lidar_stat_summary
        if (cfg%Llidar_sim) call diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                        ,sglidar%temp_tot,sglidar%beta_tot,sglidar%betaperp_tot,sglidar%beta_mol,sglidar%refl,gbx%land,gbx%ph &
                        ,LIDAR_UNDEF,ok_lidar_cfad &
                        ,stlidar%cfad_sr,stlidar%srbval &
                        ,LIDAR_NCAT,stlidar%lidarcld,stlidar%lidarcldphase,stlidar%cldlayer,stlidar%cldlayerphase &
                        ,stlidar%lidarcldtmp,stlidar%parasolrefl)
        !++++++++++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++++++++++
        if (cfg%Lradar_sim.and.cfg%Llidar_sim) call cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr, &
                                    sglidar%temp_tot,sglidar%beta_tot,sglidar%betaperp_tot,sglidar%beta_mol,sgradar%Ze_tot, &
                                    stradar%lidar_only_freq_cloud,stradar%radar_lidar_tcc, &
				    stradar%radar_tcc,stradar%radar_tcc_2)  !+JEK
   endif
   ! Replace undef
   where (stlidar%cfad_sr   == LIDAR_UNDEF) stlidar%cfad_sr   = R_UNDEF
   where (stlidar%lidarcld  == LIDAR_UNDEF) stlidar%lidarcld  = R_UNDEF
   where (stlidar%cldlayer  == LIDAR_UNDEF) stlidar%cldlayer  = R_UNDEF
   where (stlidar%parasolrefl == LIDAR_UNDEF) stlidar%parasolrefl = R_UNDEF
   where (stlidar%cldlayerphase  == LIDAR_UNDEF) stlidar%cldlayerphase  = R_UNDEF
   where (stlidar%lidarcldphase  == LIDAR_UNDEF) stlidar%lidarcldphase  = R_UNDEF
   where (stlidar%lidarcldtmp  == LIDAR_UNDEF) stlidar%lidarcldtmp  = R_UNDEF

END SUBROUTINE COSP_STATS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,Nglevels,newgrid_bot,newgrid_top,r,log_units)
   implicit none
   ! Input arguments
   integer,intent(in) :: Npoints  !# of grid points
   integer,intent(in) :: Nlevels  !# of levels
   integer,intent(in) :: Ncolumns !# of columns
   real,dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   real,dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
   integer,intent(in) :: Nglevels  !# levels in the new grid
   real,dimension(Nglevels),intent(in) :: newgrid_bot ! Lower boundary of new levels  [m]
   real,dimension(Nglevels),intent(in) :: newgrid_top ! Upper boundary of new levels  [m]
   logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
   ! Output
   real,dimension(Npoints,Ncolumns,Nglevels),intent(out) :: r ! Variable on new grid

   ! Local variables
   integer :: i,j,k
   logical :: lunits
   integer :: l
   real :: w ! Weight
   real :: dbb, dtb, dbt, dtt ! Distances between edges of both grids
   integer :: Nw  ! Number of weights
   real :: wt  ! Sum of weights
   real,dimension(Nlevels) :: oldgrid_bot,oldgrid_top ! Lower and upper boundaries of model grid
   real :: yp ! Local copy of y at a particular point.
              ! This allows for change of units.

   lunits=.false.
   if (present(log_units)) lunits=log_units

   r = 0.0

   do i=1,Npoints
     ! Calculate tops and bottoms of new and old grids
     oldgrid_bot = zhalf(i,:)
     oldgrid_top(1:Nlevels-1) = oldgrid_bot(2:Nlevels)
     oldgrid_top(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     l = 0 ! Index of level in the old grid
     ! Loop over levels in the new grid
     do k = 1,Nglevels
       Nw = 0 ! Number of weigths
       wt = 0.0 ! Sum of weights
       ! Loop over levels in the old grid and accumulate total for weighted average
       do
         l = l + 1
         w = 0.0 ! Initialise weight to 0
         ! Distances between edges of both grids
         dbb = oldgrid_bot(l) - newgrid_bot(k)
         dtb = oldgrid_top(l) - newgrid_bot(k)
         dbt = oldgrid_bot(l) - newgrid_top(k)
         dtt = oldgrid_top(l) - newgrid_top(k)
         if (dbt >= 0.0) exit ! Do next level in the new grid
         if (dtb > 0.0) then
           if (dbb <= 0.0) then
             if (dtt <= 0) then
               w = dtb
             else
               w = newgrid_top(k) - newgrid_bot(k)
             endif
           else
             if (dtt <= 0) then
               w = oldgrid_top(l) - oldgrid_bot(l)
             else
               w = -dbt
             endif
           endif
           ! If layers overlap (w/=0), then accumulate
           if (w /= 0.0) then
             Nw = Nw + 1
             wt = wt + w
             do j=1,Ncolumns
               if (lunits) then
                 if (y(i,j,l) /= R_UNDEF) then
                   yp = 10.0**(y(i,j,l)/10.0)
                 else
                   yp = 0.0
                 endif
               else
                 yp = y(i,j,l)
               endif
               r(i,j,k) = r(i,j,k) + w*yp
             enddo
           endif
         endif
       enddo
       l = l - 2
       if (l < 1) l = 0
       ! Calculate average in new grid
       if (Nw > 0) then
         do j=1,Ncolumns
           r(i,j,k) = r(i,j,k)/wt
         enddo
       endif
     enddo
   enddo

   ! Set points under surface to R_UNDEF, and change to dBZ if necessary
   do k=1,Nglevels
     do j=1,Ncolumns
       do i=1,Npoints
         if (newgrid_top(k) > zhalf(i,1)) then ! Level above model bottom level
           if (lunits) then
             if (r(i,j,k) <= 0.0) then
               r(i,j,k) = R_UNDEF
             else
               r(i,j,k) = 10.0*log10(r(i,j,k))
             endif
           endif
         else ! Level below surface
           r(i,j,k) = R_GROUND
         endif
       enddo
     enddo
   enddo

END SUBROUTINE COSP_CHANGE_VERTICAL_GRID

END MODULE MOD_COSP_STATS
