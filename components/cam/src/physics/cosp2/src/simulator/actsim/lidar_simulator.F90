! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Centre National de la Recherche Scientifique
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
! May 2007: ActSim code of M. Chiriaco and H. Chepfer rewritten by S. Bony
!
! May 2008, H. Chepfer:
! - Units of pressure inputs: Pa 
! - Non Spherical particles : LS Ice NS coefficients, CONV Ice NS coefficients
! - New input: ice_type (0=ice-spheres ; 1=ice-non-spherical)
!
! June 2008, A. Bodas-Salcedo:
! - Ported to Fortran 90 and optimisation changes
!
! August 2008, J-L Dufresne:
! - Optimisation changes (sum instructions suppressed)
!
! October 2008, S. Bony,  H. Chepfer and J-L. Dufresne :  
! - Interface with COSP v2.0:
!      cloud fraction removed from inputs
!      in-cloud condensed water now in input (instead of grid-averaged value)
!      depolarisation diagnostic removed
!      parasol (polder) reflectances (for 5 different solar zenith angles) added
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - Modification of the integration of the lidar equation.
! - change the cloud detection threshold
!
! April 2008, A. Bodas-Salcedo:
! - Bug fix in computation of pmol and pnorm of upper layer
!
! April 2008, J-L. Dufresne
! - Bug fix in computation of pmol and pnorm, thanks to Masaki Satoh: a factor 2 
! was missing. This affects the ATB values but not the cloud fraction. 
!
! January 2013, G. Cesana and H. Chepfer:
! - Add the perpendicular component of the backscattered signal (pnorm_perp_tot) in the arguments
! - Add the temperature for each levels (temp) in the arguments
! - Add the computation of the perpendicular component of the backscattered lidar signal 
! Reference: Cesana G. and H. Chepfer (2013): Evaluation of the cloud water phase
! in a climate model using CALIPSO-GOCCP, J. Geophys. Res., doi: 10.1002/jgrd.50376
!
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_lidar_simulator
  USE COSP_KINDS,         ONLY: wp
  USE MOD_COSP_CONFIG,    ONLY: SR_BINS,S_CLD,S_ATT,S_CLD_ATT,R_UNDEF,calipso_histBsct,  &
                                use_vgrid,vgrid_zl,vgrid_zu
  USE MOD_COSP_STATS,     ONLY: COSP_CHANGE_VERTICAL_GRID,hist1d
  implicit none
  
  ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
  ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
  ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
  ! (Cesana and Chepfer, JGR, 2013).
  !
  ! Relationship between ATBice and ATBperp,ice for ice particles:
  !                ATBperp,ice = Alpha*ATBice 
  ! Relationship between ATBice and ATBperp,ice for liquid particles:
  !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
  real(wp) :: &
       alpha,beta,gamma    

contains
  ! ######################################################################################
  ! SUBROUTINE lidar_subcolumn
  ! Inputs with a vertical dimensions (nlev) should ordered in along the vertical 
  ! dimension from TOA-2-SFC, for example: varIN(nlev) is varIN @ SFC. 
  ! ######################################################################################
  subroutine lidar_subcolumn(npoints,ncolumns,nlev,beta_mol,tau_mol,betatot,tautot,    &
                         betatot_ice,tautot_ice,betatot_liq,tautot_liq,                  &
                         pmol,pnorm,pnorm_perp_tot)

    ! INPUTS
    INTEGER,intent(in) :: & 
         npoints,      & ! Number of gridpoints
         ncolumns,     & ! Number of subcolumns
         nlev            ! Number of levels
    REAL(WP),intent(in),dimension(npoints,nlev) :: &
         beta_mol,     & ! Molecular backscatter coefficient
         tau_mol         ! Molecular optical depth

    REAL(WP),intent(in),dimension(npoints,ncolumns,nlev)       :: &
         betatot,      & ! 
         tautot,       & ! Optical thickess integrated from top
         betatot_ice,  & ! Backscatter coefficient for ice particles
         betatot_liq,  & ! Backscatter coefficient for liquid particles
         tautot_ice,   & ! Total optical thickness of ice
         tautot_liq      ! Total optical thickness of liq

    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,nlev) :: &
         pmol            ! Molecular attenuated backscatter lidar signal power(m^-1.sr^-1)
    REAL(WP),intent(out),dimension(npoints,ncolumns,nlev) :: &
         pnorm,        & ! Molecular backscatter signal power (m^-1.sr^-1)
         pnorm_perp_tot  ! Perpendicular lidar backscattered signal power

    ! LOCAL VARIABLES
    INTEGER :: k,icol
    REAL(WP),dimension(npoints) :: &
         tautot_lay        !
    REAL(WP),dimension(npoints,ncolumns,nlev) :: &
         pnorm_liq,      & ! Lidar backscattered signal power for liquid
         pnorm_ice,      & ! Lidar backscattered signal power for ice
         pnorm_perp_ice, & ! Perpendicular lidar backscattered signal power for ice
         pnorm_perp_liq, & ! Perpendicular lidar backscattered signal power for liq
         beta_perp_ice,  & ! Perpendicular backscatter coefficient for ice
         beta_perp_liq     ! Perpendicular backscatter coefficient for liquid    

    ! ####################################################################################
    ! *) Molecular signal
    ! ####################################################################################
    call cmp_backsignal(nlev,npoints,beta_mol(1:npoints,1:nlev),&
                        tau_mol(1:npoints,1:nlev),pmol(1:npoints,1:nlev))
                        
    ! ####################################################################################
    ! PLANE PARRALLEL FIELDS
    ! ####################################################################################
    do icol=1,ncolumns
       ! #################################################################################
       ! *) Total Backscatter signal
       ! #################################################################################
       call cmp_backsignal(nlev,npoints,betatot(1:npoints,icol,1:nlev),&
            tautot(1:npoints,icol,1:nlev),pnorm(1:npoints,icol,1:nlev))
       ! #################################################################################
       ! *) Ice/Liq Backscatter signal
       ! #################################################################################
       ! Computation of the ice and liquid lidar backscattered signal (ATBice and ATBliq)
       ! Ice only
       call cmp_backsignal(nlev,npoints,betatot_ice(1:npoints,icol,1:nlev),&
                      tautot_ice(1:npoints,icol,1:nlev),&
                      pnorm_ice(1:npoints,icol,1:nlev))
       ! Liquid only
       call cmp_backsignal(nlev,npoints,betatot_liq(1:npoints,icol,1:nlev),&
                      tautot_liq(1:npoints,icol,1:nlev),&
                      pnorm_liq(1:npoints,icol,1:nlev))
    enddo

    ! ####################################################################################
    ! PERDENDICULAR FIELDS
    ! ####################################################################################
    do icol=1,ncolumns

       ! #################################################################################
       ! *) Ice/Liq Perpendicular Backscatter signal
       ! #################################################################################
       ! Computation of ATBperp,ice/liq from ATBice/liq including the multiple scattering 
       ! contribution (Cesana and Chepfer 2013, JGR)
       do k=1,nlev
          ! Ice particles
          pnorm_perp_ice(1:npoints,icol,k) = Alpha * pnorm_ice(1:npoints,icol,k)

          ! Liquid particles
          pnorm_perp_liq(1:npoints,icol,k) = 1000._wp*Beta*pnorm_liq(1:npoints,icol,k)**2+&
               Gamma*pnorm_liq(1:npoints,icol,k) 
       enddo
  
       ! #################################################################################
       ! *) Computation of beta_perp_ice/liq using the lidar equation
       ! #################################################################################
       ! Ice only
       call cmp_beta(nlev,npoints,pnorm_perp_ice(1:npoints,icol,1:nlev),&
              tautot_ice(1:npoints,icol,1:nlev),beta_perp_ice(1:npoints,icol,1:nlev))        
 
       ! Liquid only
       call cmp_beta(nlev,npoints,pnorm_perp_liq(1:npoints,icol,1:nlev),&
            tautot_liq(1:npoints,icol,1:nlev),beta_perp_liq(1:npoints,icol,1:nlev))
          
       ! #################################################################################
       ! *) Perpendicular Backscatter signal
       ! #################################################################################
       ! Computation of the total perpendicular lidar signal (ATBperp for liq+ice)
       ! Upper layer
       WHERE(tautot(1:npoints,icol,1) .gt. 0)
          pnorm_perp_tot(1:npoints,icol,1) = (beta_perp_ice(1:npoints,icol,1)+           &
               beta_perp_liq(1:npoints,icol,1)-                                          &
               (beta_mol(1:npoints,1)/(1._wp+1._wp/0.0284_wp))) /                        &
               (2._wp*tautot(1:npoints,icol,1))*                                         &
               (1._wp-exp(-2._wp*tautot(1:npoints,icol,1)))
       ELSEWHERE
          pnorm_perp_tot(1:npoints,icol,1) = 0._wp
       ENDWHERE                                                  
             
       ! Other layers
       do k=2,nlev
          ! Optical thickness of layer k
          tautot_lay(1:npoints) = tautot(1:npoints,icol,k)-tautot(1:npoints,icol,k-1) 

          ! The perpendicular component of the molecular backscattered signal (Betaperp) 
          ! has been taken into account two times (once for liquid and once for ice). 
          ! We remove one contribution using 
          ! Betaperp=beta_mol(:,k)/(1+1/0.0284)) [bodhaine et al. 1999] in the following 
          ! equations:
          WHERE (pnorm(1:npoints,icol,k) .eq. 0)
             pnorm_perp_tot(1:npoints,icol,k)=0._wp
          ELSEWHERE
             where(tautot_lay(1:npoints) .gt. 0.)
                pnorm_perp_tot(1:npoints,icol,k) = (beta_perp_ice(1:npoints,icol,k)+     &
                   beta_perp_liq(1:npoints,icol,k)-(beta_mol(1:npoints,k)/(1._wp+1._wp/  &
                   0.0284_wp)))*EXP(-2._wp*tautot(1:npoints,icol,k-1))/                  &
                   (2._wp*tautot_lay(1:npoints))* (1._wp-EXP(-2._wp*tautot_lay(1:npoints)))
             elsewhere
                pnorm_perp_tot(1:npoints,icol,k) = (beta_perp_ice(1:npoints,icol,k)+     &
                   beta_perp_liq(1:npoints,icol,k)-(beta_mol(1:npoints,k)/(1._wp+1._wp/  &
                   0.0284_wp)))*EXP(-2._wp*tautot(1:npoints,icol,k-1))
             endwhere 
          ENDWHERE
       END DO
    enddo
  end subroutine lidar_subcolumn

  ! ######################################################################################
  ! SUBROUTINE lidar_column
  ! ######################################################################################
  subroutine lidar_column(npoints,ncol,nlevels,llm,max_bin,tmp, pnorm,                   &
                           pnorm_perp, pmol, pplay, ok_lidar_cfad, ncat, cfad2,    &
                           lidarcld, lidarcldphase, cldlayer, zlev, zlev_half,           &
                           cldlayerphase, lidarcldtmp)
    integer,parameter :: &
         nphase = 6 ! Number of cloud layer phase types

    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal grid points
         ncol,    & ! Number of subcolumns
         nlevels, & ! Number of vertical layers (OLD grid)
         llm,     & ! Number of vertical layers (NEW grid)
         max_bin, & ! Number of bins for SR CFADs
         ncat       ! Number of cloud layer types (low,mid,high,total)
    real(wp),intent(in),dimension(npoints,ncol,Nlevels) :: &
         pnorm,   & ! Lidar ATB
         pnorm_perp ! Lidar perpendicular ATB
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         pmol,    & ! Molecular ATB
         pplay,   & ! Pressure on model levels (Pa)
         tmp        ! Temperature at each levels
    logical,intent(in) :: &
         ok_lidar_cfad ! True if lidar CFAD diagnostics need to be computed
    real(wp),intent(in),dimension(npoints,nlevels) :: &
         zlev        ! Model full levels
    real(wp),intent(in),dimension(npoints,nlevels+1) :: &
         zlev_half   ! Model half levels
         
    ! Outputs
    real(wp),intent(inout),dimension(npoints,llm) :: &
         lidarcld      ! 3D "lidar" cloud fraction
    real(wp),intent(inout),dimension(npoints,ncat) :: &
         cldlayer      ! "lidar" cloud layer fraction (low, mid, high, total)
    real(wp),intent(inout),dimension(npoints,llm,nphase) :: &
         lidarcldphase ! 3D "lidar" phase cloud fraction
    real(wp),intent(inout),dimension(npoints,40,5) :: &
         lidarcldtmp   ! 3D "lidar" phase cloud fraction as a function of temp
    real(wp),intent(inout),dimension(npoints,ncat,nphase) :: &
         cldlayerphase ! "lidar" phase low mid high cloud fraction
    real(wp),intent(inout),dimension(npoints,max_bin,llm) :: &
         cfad2         ! CFADs of SR

    ! Local Variables
    integer :: ic,i,j
    real(wp),dimension(npoints,ncol,llm) :: &
         x3d
    real(wp),dimension(npoints,llm) :: &
         x3d_c,pnorm_c
    real(wp)  :: &
         xmax
    real(wp),dimension(npoints,1,Nlevels) :: t_in,ph_in,betamol_in
    real(wp),dimension(npoints,ncol,llm)  :: pnormFlip,pnorm_perpFlip
    real(wp),dimension(npoints,1,llm)     :: tmpFlip,pplayFlip,betamolFlip

    ! Vertically regrid input data
    if (use_vgrid) then 
       t_in(:,1,:)=tmp(:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            t_in,llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),tmpFlip(:,1,llm:1:-1))
       ph_in(:,1,:) = pplay(:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            ph_in,llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),pplayFlip(:,1,llm:1:-1))
       betamol_in(:,1,:) = pmol(:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            betamol_in,llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),betamolFlip(:,1,llm:1:-1))
       call cosp_change_vertical_grid(Npoints,Ncol,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            pnorm(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),pnormFlip(:,:,llm:1:-1))
       call cosp_change_vertical_grid(Npoints,Ncol,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            pnorm_perp(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),pnorm_perpFlip(:,:,llm:1:-1))
    endif

    ! Initialization (The histogram bins, are set up during initialization and the
    ! maximum value is used as the upper bounds.)
    xmax = maxval(calipso_histBsct)

    ! Compute LIDAR scattering ratio
    if (use_vgrid) then
       do ic = 1, ncol
          pnorm_c = pnormFlip(:,ic,:)
          where ((pnorm_c .lt. xmax) .and. (betamolFlip(:,1,:) .lt. xmax) .and.          &
                (betamolFlip(:,1,:) .gt. 0.0 ))
             x3d_c = pnorm_c/betamolFlip(:,1,:)
          elsewhere
             x3d_c = R_UNDEF
          end where
          x3d(:,ic,:) = x3d_c
       enddo
       ! Diagnose cloud fractions for subcolumn lidar scattering ratios
       CALL COSP_CLDFRAC(npoints,ncol,llm,ncat,nphase,tmpFlip,x3d,pnormFlip,             &
                         pnorm_perpFlip,pplayFlip,S_att,S_cld,S_cld_att,R_UNDEF,         &
                         lidarcld,cldlayer,lidarcldphase,cldlayerphase,lidarcldtmp)                           
    else
       do ic = 1, ncol
          pnorm_c = pnorm(:,ic,:)
          where ((pnorm_c.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 ))
             x3d_c = pnorm_c/pmol
          elsewhere
             x3d_c = R_UNDEF
          end where
          x3d(:,ic,:) = x3d_c
       enddo
       ! Diagnose cloud fractions for subcolumn lidar scattering ratios
       CALL COSP_CLDFRAC(npoints,ncol,nlevels,ncat,nphase,tmp,x3d,pnorm,pnorm_perp,pplay,&
                         S_att,S_cld,S_cld_att,R_UNDEF,lidarcld,cldlayer,lidarcldphase,  &
                         cldlayerphase,lidarcldtmp)
    endif

    ! CFADs
    if (ok_lidar_cfad) then
       ! CFADs of subgrid-scale lidar scattering ratios
       do i=1,Npoints
          do j=1,llm
             cfad2(i,:,j) = hist1D(ncol,x3d(i,:,j),SR_BINS,calipso_histBsct)
          enddo
       enddo
       where(cfad2 .ne. R_UNDEF) cfad2=cfad2/ncol

    endif 
    
    ! Unit conversions
    where(lidarcld /= R_UNDEF)      lidarcld      = lidarcld*100._wp
    where(cldlayer /= R_UNDEF)      cldlayer      = cldlayer*100._wp
    where(cldlayerphase /= R_UNDEF) cldlayerphase = cldlayerphase*100._wp
    where(lidarcldphase /= R_UNDEF) lidarcldphase = lidarcldphase*100._wp
    where(lidarcldtmp /= R_UNDEF)   lidarcldtmp   = lidarcldtmp*100._wp
   
  end subroutine lidar_column

  ! ######################################################################################
  ! The subroutines below compute the attenuated backscatter signal and the lidar 
  ! backscatter coefficients using eq (1) from doi:0094-8276/08/2008GL034207
  ! ######################################################################################
  subroutine cmp_backsignal(nlev,npoints,beta,tau,pnorm)
    ! INPUTS
    integer, intent(in) :: nlev,npoints
    real(wp),intent(in),dimension(npoints,nlev) :: beta,tau

    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev) :: pnorm

    ! Internal Variables
    real(wp), dimension(npoints) :: tautot_lay
    integer :: k

    ! Uppermost layer 
    pnorm(:,1) = beta(:,1) / (2._wp*tau(:,1)) * (1._wp-exp(-2._wp*tau(:,1)))

    ! Other layers
    do k=2,nlev
       tautot_lay(:) = tau(:,k)-tau(:,k-1) 
       WHERE (tautot_lay(:) .gt. 0.)
          pnorm(:,k) = beta(:,k)*EXP(-2._wp*tau(:,k-1)) /&
               (2._wp*tautot_lay(:))*(1._wp-EXP(-2._wp*tautot_lay(:)))
       ELSEWHERE
          ! This must never happen, but just in case, to avoid div. by 0
          pnorm(:,k) = beta(:,k) * EXP(-2._wp*tau(:,k-1))
       END WHERE

    END DO
  end subroutine cmp_backsignal

  subroutine cmp_beta(nlev,npoints,pnorm,tau,beta)
    ! INPUTS
    integer, intent(in) :: nlev,npoints
    real(wp),intent(in),dimension(npoints,nlev) :: pnorm,tau

    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev) :: beta

    ! Internal Variables
    real(wp), dimension(npoints) :: tautot_lay
    integer :: k
    real(wp) :: epsrealwp

    epsrealwp = epsilon(1._wp)
    beta(:,1) = pnorm(:,1) * (2._wp*tau(:,1))/(1._wp-exp(-2._wp*tau(:,1)))
    do k=2,nlev
       tautot_lay(:) = tau(:,k)-tau(:,k-1)       
       WHERE ( EXP(-2._wp*tau(:,k-1)) .gt. epsrealwp )
          WHERE (tautot_lay(:) .gt. 0.)
             beta(:,k) = pnorm(:,k)/ EXP(-2._wp*tau(:,k-1))* &
                  (2._wp*tautot_lay(:))/(1._wp-exp(-2._wp*tautot_lay(:)))
          ELSEWHERE
             beta(:,k)=pnorm(:,k)/EXP(-2._wp*tau(:,k-1))
          END WHERE
       ELSEWHERE
          beta(:,k)=pnorm(:,k)/epsrealwp
       END WHERE
    ENDDO

  end subroutine cmp_beta
    ! ####################################################################################
    ! SUBROUTINE cosp_cldfrac
    ! Conventions: Ncat must be equal to 4
    ! ####################################################################################
    SUBROUTINE COSP_CLDFRAC(Npoints,Ncolumns,Nlevels,Ncat,Nphase,tmp,x,ATB,ATBperp,      &
                               pplay,S_att,S_cld,S_cld_att,undef,lidarcld,cldlayer,      &
                               lidarcldphase,cldlayerphase,lidarcldtemp)
    ! Parameters
    integer,parameter :: Ntemp=40 ! indice of the temperature vector
    real(wp),parameter,dimension(Ntemp+1) :: &
       tempmod = [0.0,   183.15,186.15,189.15,192.15,195.15,198.15,201.15,204.15,207.15, &
                  210.15,213.15,216.15,219.15,222.15,225.15,228.15,231.15,234.15,237.15, &
                  240.15,243.15,246.15,249.15,252.15,255.15,258.15,261.15,264.15,267.15, &
                  270.15,273.15,276.15,279.15,282.15,285.15,288.15,291.15,294.15,297.15, &
                  473.15]
         
    ! Polynomial coefficient of the phase discrimination line used to separate liquid from ice
    ! (Cesana and Chepfer, JGR, 2013)
    ! ATBperp = ATB^5*alpha50 + ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50 + ATB*epsilon50 + zeta50
    real(wp),parameter :: &
       alpha50   = 9.0322e+15_wp,  & !
       beta50    = -2.1358e+12_wp, & !
       gamma50   = 173.3963e06_wp, & !
       delta50   = -3.9514e03_wp,  & !
       epsilon50 = 0.2559_wp,      & !
       zeta50    = -9.4776e-07_wp    ! 
       
	! Inputs
    integer,intent(in) :: &
       Npoints,  & ! Number of gridpoints
       Ncolumns, & ! Number of subcolumns
       Nlevels,  & ! Number of vertical levels
       Ncat,     & ! Number of cloud layer types
       Nphase      ! Number of cloud layer phase types
	               ! [ice,liquid,undefined,false ice,false liquid,Percent of ice]
    real(wp),intent(in) :: &
       S_att,    & !
       S_cld,    & !
       S_cld_att,& ! New threshold for undefine cloud phase detection
       undef       ! Undefined value
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: &
       x,        & ! 
       ATB,      & ! 3D attenuated backscatter
       ATBperp     ! 3D attenuated backscatter (perpendicular)
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
       tmp,      & ! Temperature   
       pplay       ! Pressure

	! Outputs
    real(wp),intent(out),dimension(Npoints,Ntemp,5) :: &
       lidarcldtemp  ! 3D Temperature 1=tot,2=ice,3=liq,4=undef,5=ice/ice+liq
    real(wp),intent(out),dimension(Npoints,Nlevels,Nphase) :: &
       lidarcldphase ! 3D cloud phase fraction
    real(wp),intent(out),dimension(Npoints,Nlevels) :: &
       lidarcld      ! 3D cloud fraction
    real(wp),intent(out),dimension(Npoints,Ncat) :: &
       cldlayer      ! Low, middle, high, total cloud fractions
    real(wp),intent(out),dimension(Npoints,Ncat,Nphase) :: &
       cldlayerphase ! Low, middle, high, total cloud fractions for ice liquid and undefine phase    
    
    ! Local variables
    integer  :: &
       ip, k, iz, ic, ncol, nlev, i, itemp, toplvlsat 
    real(wp) :: &
       p1,checktemp, ATBperp_tmp,checkcldlayerphase, checkcldlayerphase2
    real(wp),dimension(Npoints,Nlevels) :: &
       nsub,lidarcldphasetmp   
    real(wp),dimension(Npoints,Ntemp) :: &
       sumlidarcldtemp,lidarcldtempind
    real(wp),dimension(Npoints,Ncolumns,Ncat) :: &
       cldlay,nsublay   
    real(wp),dimension(Npoints,Ncat) :: &
       nsublayer,cldlayerphasetmp,cldlayerphasesum
    real(wp),dimension(Npoints,Ncolumns,Nlevels) :: &   
       tmpi, & ! Temperature of ice cld
       tmpl, & ! Temperature of liquid cld
       tmpu, & ! Temperature of undef cld
       cldy, & ! 
       srok    !
    real(wp),dimension(Npoints,Ncolumns,Ncat,Nphase) :: &
       cldlayphase ! subgrided low mid high phase cloud fraction
             
    ! ####################################################################################
	! 1) Initialize    
    ! ####################################################################################
    lidarcld              = 0._wp
    nsub                  = 0._wp
    cldlay                = 0._wp
    nsublay               = 0._wp
    ATBperp_tmp           = 0._wp
    lidarcldphase(:,:,:)  = 0._wp
    cldlayphase(:,:,:,:)  = 0._wp
    cldlayerphase(:,:,:)  = 0._wp
    tmpi(:,:,:)           = 0._wp
    tmpl(:,:,:)           = 0._wp
    tmpu(:,:,:)           = 0._wp
    cldlayerphasesum(:,:) = 0._wp
    lidarcldtemp(:,:,:)   = 0._wp
    lidarcldtempind(:,:)  = 0._wp
    sumlidarcldtemp(:,:)  = 0._wp
    lidarcldphasetmp(:,:) = 0._wp
    toplvlsat             = 0

    ! ####################################################################################
    ! 2) Cloud detection
    ! ####################################################################################
    do k=1,Nlevels
       ! Cloud detection at subgrid-scale:
       where ((x(:,:,k) .gt. S_cld) .and. (x(:,:,k) .ne. undef) )
          cldy(:,:,k)=1._wp
       elsewhere
          cldy(:,:,k)=0._wp
       endwhere
       
       ! Number of usefull sub-columns:
       where ((x(:,:,k) .gt. S_att) .and. (x(:,:,k) .ne. undef) )
          srok(:,:,k)=1._wp
       elsewhere
          srok(:,:,k)=0._wp
       endwhere
    enddo    
    
    ! ####################################################################################
    ! 3) Grid-box 3D cloud fraction and layered cloud fractions(ISCCP pressure categories)
    ! ####################################################################################
    lidarcld = 0._wp
    nsub     = 0._wp
    cldlay   = 0._wp
    nsublay  = 0._wp
    do k=1,Nlevels
       do ic = 1, Ncolumns
          do ip = 1, Npoints
          
             ! Computation of the cloud fraction as a function of the temperature instead
             ! of height, for ice,liquid and all clouds
             if(srok(ip,ic,k).gt.0.)then
                do itemp=1,Ntemp
                   if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
                      lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1._wp
                   endif
                enddo
             endif
             
             if(cldy(ip,ic,k).eq.1.)then
                do itemp=1,Ntemp 
                   if( (tmp(ip,k) .ge. tempmod(itemp)).and.(tmp(ip,k) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1._wp
                   endif
                enddo
             endif

             iz=1
             p1 = pplay(ip,k)
             if ( p1.gt.0. .and. p1.lt.(440._wp*100._wp)) then ! high clouds
                iz=3
             else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid clouds
                iz=2
             endif
             
             cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
             cldlay(ip,ic,4)  = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
             lidarcld(ip,k)   = lidarcld(ip,k) + cldy(ip,ic,k)
             
             nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
             nsublay(ip,ic,4)  = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
             nsub(ip,k)        = nsub(ip,k) + srok(ip,ic,k)
             
          enddo
       enddo
    enddo   
    
    ! Grid-box 3D cloud fraction
    where ( nsub(:,:).gt.0.0 )
       lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
    elsewhere
       lidarcld(:,:) = undef
    endwhere
    
    ! Layered cloud fractions
    cldlayer  = 0._wp
    nsublayer = 0._wp
    do iz = 1, Ncat
       do ic = 1, Ncolumns
          cldlayer(:,iz)  = cldlayer(:,iz)  + cldlay(:,ic,iz)
          nsublayer(:,iz) = nsublayer(:,iz) + nsublay(:,ic,iz)
       enddo
    enddo
    where (nsublayer(:,:) .gt. 0.0)
       cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
    elsewhere
       cldlayer(:,:) = undef
    endwhere
              
    ! ####################################################################################
    ! 4) Grid-box 3D cloud Phase
    ! ####################################################################################
    
    ! ####################################################################################
    ! 4.1) For Cloudy pixels with 8.16km < z < 19.2km
    ! ####################################################################################
    do ncol=1,Ncolumns
       do i=1,Npoints          
          do nlev=1,23 ! from 19.2km until 8.16km
               p1 = pplay(1,nlev)

             ! Avoid zero values
             if( (cldy(i,ncol,nlev).eq.1.) .and. (ATBperp(i,ncol,nlev).gt.0.) )then
                ! Computation of the ATBperp along the phase discrimination line
                ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                     (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                     ATB(i,ncol,nlev)*epsilon50 + zeta50     
                ! ########################################################################
                ! 4.1.a) Ice: ATBperp above the phase discrimination line
                ! ########################################################################
                if((ATBperp(i,ncol,nlev)-ATBperp_tmp) .ge. 0.)then ! Ice clouds

                   ! ICE with temperature above 273,15°K = Liquid (false ice)
                   if(tmp(i,nlev) .gt. 273.15) then ! Temperature above 273,15 K
                     ! Liquid: False ice corrected by the temperature to Liquid
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp ! False ice detection ==> added to Liquid
                                    
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,5) = lidarcldphase(i,nlev,5)+1._wp ! Keep the information "temperature criterium used"                      
                                                                              ! to classify the phase cloud
                      cldlayphase(i,ncol,4,2) = 1. ! tot cloud
                      if (p1 .gt. 0. .and. p1.lt.(440._wp*100._wp)) then ! high cloud
                         cldlayphase(i,ncol,3,2) = 1._wp
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then ! mid cloud
                         cldlayphase(i,ncol,2,2) = 1._wp
                      else ! low cloud
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                      cldlayphase(i,ncol,4,5) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,5) = 1._wp
                      ! Middle cloud
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,5) = 1._wp
                      ! Low cloud
                      else 
                         cldlayphase(i,ncol,1,5) = 1._wp
                      endif
                   else
                      ! ICE with temperature below 273,15°K
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud 
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then 
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                ! ########################################################################
                ! 4.1.b) Liquid: ATBperp below the phase discrimination line
                ! ########################################################################
                else
                   ! Liquid with temperature above 231,15°K
                   if(tmp(i,nlev) .gt. 231.15_wp) then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                   else
                      ! Liquid with temperature below 231,15°K = Ice (false liquid)
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp ! false liquid detection ==> added to ice
                      lidarcldphase(i,nlev,4) = lidarcldphase(i,nlev,4)+1._wp
                      cldlayphase(i,ncol,4,4) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,4) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,4) = 1._wp
                      ! Low cloud
                      else
                         cldlayphase(i,ncol,1,4) = 1._wp
                      endif
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                endif ! end of discrimination condition
             endif ! end of cloud condition
          enddo ! end of altitude loop

          ! ##############################################################################
          ! 4.2) For Cloudy pixels with 0km < z < 8.16km
          ! ##############################################################################
          toplvlsat = 0
          do nlev=24,Nlevels! from 8.16km until 0km
             p1 = pplay(i,nlev)

             if((cldy(i,ncol,nlev) .eq. 1.) .and. (ATBperp(i,ncol,nlev) .gt. 0.) )then
                ! Computation of the ATBperp of the phase discrimination line
                ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                     (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                     ATB(i,ncol,nlev)*epsilon50 + zeta50
                ! ########################################################################
                ! 4.2.a) Ice: ATBperp above the phase discrimination line
                ! ########################################################################
                ! ICE with temperature above 273,15°K = Liquid (false ice)
                if((ATBperp(i,ncol,nlev)-ATBperp_tmp) .ge. 0.)then ! Ice clouds
                   if(tmp(i,nlev) .gt. 273.15)then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp ! false ice ==> liq
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,5) = lidarcldphase(i,nlev,5)+1._wp
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then 
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud
                      else 
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                      
                      cldlayphase(i,ncol,4,5) = 1. ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,5) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,5) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,5) = 1._wp
                      endif
                   else
                      ! ICE with temperature below 273,15°K
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp
                     tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt.(680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                   
                ! ########################################################################
                ! 4.2.b) Liquid: ATBperp below the phase discrimination line
                ! ########################################################################
                else
                   ! Liquid with temperature above 231,15°K
                   if(tmp(i,nlev) .gt. 231.15)then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                   else
                      ! Liquid with temperature below 231,15°K = Ice (false liquid)
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp ! false liq ==> ice
                      lidarcldphase(i,nlev,4) = lidarcldphase(i,nlev,4)+1._wp ! false liq ==> ice
                      cldlayphase(i,ncol,4,4) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,4) = 1._wp
                      ! Middle   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,4) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,4) = 1._wp
                      endif
                      
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                endif ! end of discrimination condition
                
                toplvlsat=0
                
                ! Find the level of the highest cloud with SR>30
                if(x(i,ncol,nlev) .gt. S_cld_att) then ! SR > 30.
                    toplvlsat = nlev+1
                    goto 99
                endif
             endif ! end of cloud condition
          enddo ! end of altitude loop
99        continue
          
          ! ##############################################################################
          ! Undefined phase: For a cloud located below another cloud with SR>30
          ! see Cesana and Chepfer 2013 Sect.III.2
          ! ##############################################################################
          if(toplvlsat.ne.0) then
             do nlev = toplvlsat,Nlevels
                p1 = pplay(i,nlev)
                if(cldy(i,ncol,nlev).eq.1.)then
                   lidarcldphase(i,nlev,3) = lidarcldphase(i,nlev,3)+1._wp
                   tmpu(i,ncol,nlev)       = tmp(i,nlev)
                   cldlayphase(i,ncol,4,3) = 1._wp ! tot cloud
                   ! High cloud
                   if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                      cldlayphase(i,ncol,3,3) = 1._wp
                   ! Middle cloud   
                   else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                      cldlayphase(i,ncol,2,3) = 1._wp
                   ! Low cloud   
                   else
                      cldlayphase(i,ncol,1,3) = 1._wp
                   endif
                endif
             enddo
          endif
          toplvlsat=0
       enddo
    enddo
     
    ! ####################################################################################
    ! Computation of final cloud phase diagnosis
    ! ####################################################################################

    ! Compute the Ice percentage in cloud = ice/(ice+liq) as a function of the occurrences
    lidarcldphasetmp(:,:) = lidarcldphase(:,:,1)+lidarcldphase(:,:,2);
    WHERE (lidarcldphasetmp(:,:) .gt. 0.)
       lidarcldphase(:,:,6)=lidarcldphase(:,:,1)/lidarcldphasetmp(:,:)
    ELSEWHERE
       lidarcldphase(:,:,6) = undef
    ENDWHERE
    
    ! Compute Phase 3D Cloud Fraction
    !WHERE (nsub(:,Nlevels:1:-1) .gt. 0.0 )
    WHERE (nsub(:,:) .gt. 0.0 )  
       lidarcldphase(:,:,1)=lidarcldphase(:,:,1)/nsub(:,:)
       lidarcldphase(:,:,2)=lidarcldphase(:,:,2)/nsub(:,:)
       lidarcldphase(:,:,3)=lidarcldphase(:,:,3)/nsub(:,:)
       lidarcldphase(:,:,4)=lidarcldphase(:,:,4)/nsub(:,:)
       lidarcldphase(:,:,5)=lidarcldphase(:,:,5)/nsub(:,:)
    ELSEWHERE
       lidarcldphase(:,:,1) = undef
       lidarcldphase(:,:,2) = undef
       lidarcldphase(:,:,3) = undef
       lidarcldphase(:,:,4) = undef
       lidarcldphase(:,:,5) = undef
    ENDWHERE

    ! Compute Phase low mid high cloud fractions
    do iz = 1, Ncat
       do i=1,Nphase-3
          do ic = 1, Ncolumns
             cldlayerphase(:,iz,i)  = cldlayerphase(:,iz,i)  + cldlayphase(:,ic,iz,i)
             cldlayerphasesum(:,iz) = cldlayerphasesum(:,iz) + cldlayphase(:,ic,iz,i)
          enddo
       enddo
    enddo
    do iz = 1, Ncat
       do i=4,5
          do ic = 1, Ncolumns
             cldlayerphase(:,iz,i) = cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)
          enddo
       enddo
    enddo
    
    ! Compute the Ice percentage in cloud = ice/(ice+liq)
    cldlayerphasetmp(:,:)=cldlayerphase(:,:,1)+cldlayerphase(:,:,2)
    WHERE (cldlayerphasetmp(:,:).gt. 0.)
       cldlayerphase(:,:,6)=cldlayerphase(:,:,1)/cldlayerphasetmp(:,:)
    ELSEWHERE
       cldlayerphase(:,:,6) = undef
    ENDWHERE
    
    do i=1,Nphase-1
       WHERE ( cldlayerphasesum(:,:).gt.0.0 )
          cldlayerphase(:,:,i) = (cldlayerphase(:,:,i)/cldlayerphasesum(:,:)) * cldlayer(:,:)
       ENDWHERE
    enddo
    
    do i=1,Npoints
       do iz=1,Ncat
          checkcldlayerphase=0.
          checkcldlayerphase2=0.
          if (cldlayerphasesum(i,iz) .gt. 0.0 )then
             do ic=1,Nphase-3
                checkcldlayerphase = checkcldlayerphase+cldlayerphase(i,iz,ic)
             enddo
             checkcldlayerphase2 = cldlayer(i,iz)-checkcldlayerphase
             if((checkcldlayerphase2 .gt. 0.01) .or. (checkcldlayerphase2 .lt. -0.01) ) print *, checkcldlayerphase,cldlayer(i,iz)
          endif
       enddo
    enddo
    
    do i=1,Nphase-1
       WHERE (nsublayer(:,:) .eq. 0.0)
          cldlayerphase(:,:,i) = undef
       ENDWHERE
    enddo
 
    ! Compute Phase 3D as a function of temperature
    do nlev=1,Nlevels
       do ncol=1,Ncolumns
          do i=1,Npoints
             do itemp=1,Ntemp
                if(tmpi(i,ncol,nlev).gt.0.)then
                   if((tmpi(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpi(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,2)=lidarcldtemp(i,itemp,2)+1._wp
                   endif
                elseif(tmpl(i,ncol,nlev) .gt. 0.)then
                   if((tmpl(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpl(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,3)=lidarcldtemp(i,itemp,3)+1._wp
                   endif
                elseif(tmpu(i,ncol,nlev) .gt. 0.)then
                   if((tmpu(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpu(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,4)=lidarcldtemp(i,itemp,4)+1._wp
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
    
    ! Check temperature cloud fraction
    do i=1,Npoints
       do itemp=1,Ntemp
          checktemp=lidarcldtemp(i,itemp,2)+lidarcldtemp(i,itemp,3)+lidarcldtemp(i,itemp,4)
          !if(checktemp .NE. lidarcldtemp(i,itemp,1))then
          !   print *, i,itemp
          !   print *, lidarcldtemp(i,itemp,1:4)
          !endif
          
       enddo
    enddo
    
    ! Compute the Ice percentage in cloud = ice/(ice+liq)
    sumlidarcldtemp(:,:)=lidarcldtemp(:,:,2)+lidarcldtemp(:,:,3)    
    WHERE(sumlidarcldtemp(:,:) .gt. 0.)
       lidarcldtemp(:,:,5)=lidarcldtemp(:,:,2)/sumlidarcldtemp(:,:)
    ELSEWHERE
       lidarcldtemp(:,:,5)=undef
    ENDWHERE
    
    do i=1,4
       WHERE(lidarcldtempind(:,:) .gt. 0.)
          lidarcldtemp(:,:,i) = lidarcldtemp(:,:,i)/lidarcldtempind(:,:)
       ELSEWHERE
          lidarcldtemp(:,:,i) = undef
       ENDWHERE
    enddo
    
    RETURN
  END SUBROUTINE COSP_CLDFRAC

end module mod_lidar_simulator
