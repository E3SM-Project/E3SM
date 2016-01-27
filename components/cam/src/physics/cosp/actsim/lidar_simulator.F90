! Copyright (c) 2009, Centre National de la Recherche Scientifique
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 07:08:38 -0700 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/actsim/lidar_simulator.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the LMD/IPSL/CNRS/UPMC nor the names of its
!       contributors may be used to endorse or promote products derived from this software without 
!       specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      
      SUBROUTINE lidar_simulator(npoints,nlev,npart,nrefl &
                , undef &
                , pres, presf, temp &
                , q_lsliq, q_lsice, q_cvliq, q_cvice &
                , ls_radliq, ls_radice, cv_radliq, cv_radice &
                , ice_type, pmol, pnorm, pnorm_perp_tot,tautot, refl &
		, ls_radsnow, q_lssnow, prec_frac ) !+JEK+ZYY
!
!---------------------------------------------------------------------------------
! Purpose: To compute lidar signal from model-simulated profiles of cloud water
!          and cloud fraction in each sub-column of each model gridbox.
!
! References: 
! Chepfer H., S. Bony, D. Winker, M. Chiriaco, J.-L. Dufresne, G. Seze (2008),
! Use of CALIPSO lidar observations to evaluate the cloudiness simulated by a 
! climate model, Geophys. Res. Lett., 35, L15704, doi:10.1029/2008GL034207.
!
! Previous references:
! Chiriaco et al, MWR, 2006; Chepfer et al., MWR, 2007
!
! Contacts: Helene Chepfer (chepfer@lmd.polytechnique.fr), Sandrine Bony (bony@lmd.jussieu.fr)
!
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
!---------------------------------------------------------------------------------
!
! Inputs:
!  npoints  : number of horizontal points
!  nlev : number of vertical levels
!  npart: numberb of cloud meteors (stratiform_liq, stratiform_ice, conv_liq, conv_ice). 
!        (stratiform snow). !+JEK+ZYY
!        Currently npart must be 4 or 5 !+JEK+ZYY
!  nrefl: number of solar zenith angles for parasol reflectances
!  pres : pressure in the middle of atmospheric layers (full levels): Pa
!  presf: pressure in the interface of atmospheric layers (half levels): Pa
!     presf(..,1) : surface pressure ; presf(..,nlev+1)= TOA pressure
!  temp : temperature of atmospheric layers: K
!  q_lsliq: LS sub-column liquid water mixing ratio (kg/kg)
!  q_lsice: LS sub-column ice water mixing ratio (kg/kg)
!  q_cvliq: CONV sub-column liquid water mixing ratio (kg/kg)
!  q_cvice: CONV sub-column ice water mixing ratio (kg/kg)
!  q_lssnow: LS sub-column snow mixing ratio (kg/kg) !+JEK+ZYY
!  ls_radliq: effective radius of LS liquid particles (meters)
!  ls_radice: effective radius of LS ice particles (meters)
!  cv_radliq: effective radius of CONV liquid particles (meters)
!  cv_radice: effective radius of CONV ice particles (meters)
!  ls_radsnow:effective radius of LS snow particles (meters) !+JEK+ZYY
!  prec_frac: precip cover in each sub-column of the gridbox (output from prec_scops) !+JEK+ZYY
!  ice_type : ice particle shape hypothesis (ice_type=0 for spheres, ice_type=1 
!             for non spherical particles)
!
! Outputs:
!  pmol : molecular attenuated backscatter lidar signal power (m^-1.sr^-1)
!  pnorm: total attenuated backscatter lidar signal power (m^-1.sr^-1)
!  pnorm_perp_tot: perpendicular attenuated backscatter lidar signal power (m^-1.sr^-1)
!  tautot: optical thickess integrated from top to level z
!  refl : parasol(polder) reflectance
!
! Version 1.0 (June 2007)
! Version 1.1 (May 2008)
! Version 1.2 (June 2008)
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
!---------------------------------------------------------------------------------

      IMPLICIT NONE
      REAL :: SRsat
      PARAMETER (SRsat = 0.01) ! threshold full attenuation 

      LOGICAL ok_parasol
      PARAMETER (ok_parasol=.true.)  ! set to .true. if you want to activate parasol reflectances

      INTEGER i, k
      
      INTEGER INDX_LSLIQ,INDX_LSICE,INDX_CVLIQ,INDX_CVICE,INDX_LSSNOW  !+JEK+ZYY
      PARAMETER (INDX_LSLIQ=1,INDX_LSICE=2,INDX_CVLIQ=3,INDX_CVICE=4,INDX_LSSNOW=5) !+JEK+ZYY
! inputs:
      INTEGER npoints,nlev,npart,ice_type
      INTEGER nrefl
      real undef                 ! undefined value
      REAL pres(npoints,nlev)    ! pressure full levels
      REAL presf(npoints,nlev+1) ! pressure half levels
      REAL q_lsliq(npoints,nlev), q_lsice(npoints,nlev)
      REAL q_cvliq(npoints,nlev), q_cvice(npoints,nlev)
      REAL ls_radliq(npoints,nlev), ls_radice(npoints,nlev)
      REAL cv_radliq(npoints,nlev), cv_radice(npoints,nlev)
      REAL ls_radsnow(npoints,nlev)  !+JEK+ZYY
      REAL q_lssnow(npoints,nlev)    !+JEK+ZYY
      REAL prec_frac(npoints,nlev)   !+JEK+ZYY

! outputs (for each subcolumn):

      REAL pmol(npoints,nlev)  ! molecular backscatter signal power (m^-1.sr^-1)
      REAL pnorm(npoints,nlev) ! total lidar backscatter signal power (m^-1.sr^-1)
      REAL tautot(npoints,nlev)! optical thickess integrated from top
      REAL refl(npoints,nrefl)! parasol reflectance ! parasol

! actsim variables:

      REAL km, rdiffm, Qscat, Cmol
      PARAMETER (Cmol = 6.2446e-32) ! depends on wavelength
      PARAMETER (km = 1.38e-23)     ! Boltzmann constant (J/K)

      PARAMETER (rdiffm = 0.7)      ! multiple scattering correction parameter
      PARAMETER (Qscat = 2.0)       ! particle scattering efficiency at 532 nm

      REAL rholiq, rhoice
      PARAMETER (rholiq=1.0e+03)     ! liquid water (kg/m3)
      PARAMETER (rhoice=0.5e+03)     ! ice (kg/m3)

      REAL pi, rhopart(npart)
      REAL polpart(npart,5)  ! polynomial coefficients derived for spherical and non spherical
                             ! particules
      
      real, parameter :: TAUTOT_SMALL = 1.0e-8 ! smallest layer tau value in calculating beta_perp
                                               ! to avoid divide by zero error.
                                               ! This value works for double precision
      real, parameter :: TAUTOT_MAX = -0.5*log(tiny(TAUTOT_SMALL)) ! maximum allowed attenuation from
                                           ! one layer to TOA in beta_perp calculations to avoid 
                                           ! float underflow

!   grid-box variables:
      REAL rad_part(npoints,nlev,npart)
      REAL rhoair(npoints,nlev), zheight(npoints,nlev+1)
      REAL beta_mol(npoints,nlev), alpha_mol(npoints,nlev)
      REAL kp_part(npoints,nlev,npart)

!   sub-column variables:
      REAL qpart(npoints,nlev,npart) ! mixing ratio particles in each subcolumn
      REAL alpha_part(npoints,nlev,npart)
      REAL tau_mol_lay(npoints)  ! temporary variable, moL. opt. thickness of layer k
      REAL tau_mol(npoints,nlev) ! optical thickness between TOA and bottom of layer k
      REAL tau_part(npoints,nlev,npart)
      REAL betatot(npoints,nlev)
      REAL tautot_lay(npoints)   ! temporary variable, total opt. thickness of layer k
!     Optical thickness from TOA to surface for Parasol
     REAL tautot_S_liq(npoints),tautot_S_ice(npoints)     ! for liq and ice clouds


! Local variables
      REAL Alpha, Beta, Gamma  ! Polynomial coefficient for ATBperp computation
      REAL temp(npoints,nlev)                   ! temperature of layer k
      REAL betatot_ice(npoints,nlev)    ! backscatter coefficient for ice particles
      REAL beta_perp_ice(npoints,nlev)  ! perpendicular backscatter coefficient for ice
      REAL betatot_liq(npoints,nlev)    ! backscatter coefficient for liquid particles
      REAL beta_perp_liq(npoints,nlev)  ! perpendicular backscatter coefficient for liq
      REAL tautot_ice(npoints,nlev)     ! total optical thickness of ice
      REAL tautot_liq(npoints,nlev)     ! total optical thickness of liq
      REAL tautot_lay_ice(npoints)    ! total optical thickness of ice in the layer k
      REAL tautot_lay_liq(npoints)    ! total optical thickness of liq in the layer k
      REAL pnorm_liq(npoints,nlev)    ! lidar backscattered signal power for liquid
      REAL pnorm_ice(npoints,nlev)    ! lidar backscattered signal power for ice
      REAL pnorm_perp_ice(npoints,nlev) ! perpendicular lidar backscattered signal power for ice
      REAL pnorm_perp_liq(npoints,nlev) ! perpendicular lidar backscattered signal power for liq

! Output variable
      REAL pnorm_perp_tot (npoints,nlev) ! perpendicular lidar backscattered signal power

!------------------------------------------------------------
!---- 0. Initialisation :
!------------------------------------------------------------
betatot_ice(:,:)=0
betatot_liq(:,:)=0
beta_perp_ice(:,:)=0
beta_perp_liq(:,:)=0
tautot_ice(:,:)=0
tautot_liq(:,:)=0
tautot_lay_ice(:)=0;
tautot_lay_liq(:)=0;
pnorm_liq(:,:)=0
pnorm_ice(:,:)=0
pnorm_perp_ice(:,:)=0
pnorm_perp_liq(:,:)=0
pnorm_perp_tot(:,:)=0


! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the ATBperpendicular
! as a function of the ATB for ice or liquid cloud particles derived from CALIPSO-GOCCP
! observations at 120m vertical grid (Cesana and Chepfer, JGR, 2013).
!
! Relationship between ATBice and ATBperp,ice for ice particles
!  ATBperp,ice = Alpha*ATBice 
         Alpha = 0.2904

! Relationship between ATBice and ATBperp,ice for liquid particles
!  ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
         Beta = 0.4099
         Gamma = 0.009

!------------------------------------------------------------
!---- 1. Preliminary definitions and calculations :
!------------------------------------------------------------

!! commented out for snow by +JEK+ZYY
!      if ( npart .ne. 4 ) then
!        print *,'Error in lidar_simulator, npart should be 4, not',npart
!        stop
!      endif

      pi = dacos(-1.D0)

! Polynomial coefficients for spherical liq/ice particles derived from Mie theory.
! Polynomial coefficients for non spherical particles derived from a composite of
! Ray-tracing theory for large particles (e.g. Noel et al., Appl. Opt., 2001)
! and FDTD theory for very small particles (Yang et al., JQSRT, 2003).

! We repeat the same coefficients for LS and CONV cloud to make code more readable
!*     LS Liquid water coefficients:
         polpart(INDX_LSLIQ,1) =  2.6980e-8
         polpart(INDX_LSLIQ,2) = -3.7701e-6
         polpart(INDX_LSLIQ,3) =  1.6594e-4
         polpart(INDX_LSLIQ,4) = -0.0024
         polpart(INDX_LSLIQ,5) =  0.0626
!*     LS Ice coefficients: 
      if (ice_type.eq.0) then     
         polpart(INDX_LSICE,1) = -1.0176e-8
         polpart(INDX_LSICE,2) =  1.7615e-6
         polpart(INDX_LSICE,3) = -1.0480e-4
         polpart(INDX_LSICE,4) =  0.0019
         polpart(INDX_LSICE,5) =  0.0460
      endif
!*     LS Ice NS coefficients: 
      if (ice_type.eq.1) then 
         polpart(INDX_LSICE,1) = 1.3615e-8
         polpart(INDX_LSICE,2) = -2.04206e-6
         polpart(INDX_LSICE,3) = 7.51799e-5
         polpart(INDX_LSICE,4) = 0.00078213
         polpart(INDX_LSICE,5) = 0.0182131
      endif
!*     LS SNOW coefficients: !+JEK+ZYY
         polpart(INDX_LSSNOW,1) = 1.3615e-8   !+JEK+ZYY
         polpart(INDX_LSSNOW,2) = -2.04206e-6 !+JEK+ZYY
         polpart(INDX_LSSNOW,3) = 7.51799e-5  !+JEK+ZYY
         polpart(INDX_LSSNOW,4) = 0.00078213  !+JEK+ZYY
         polpart(INDX_LSSNOW,5) = 0.0182131   !+JEK+ZYY
!*     CONV Liquid water coefficients:
         polpart(INDX_CVLIQ,1) =  2.6980e-8
         polpart(INDX_CVLIQ,2) = -3.7701e-6
         polpart(INDX_CVLIQ,3) =  1.6594e-4
         polpart(INDX_CVLIQ,4) = -0.0024
         polpart(INDX_CVLIQ,5) =  0.0626
!*     CONV Ice coefficients: 
      if (ice_type.eq.0) then 
         polpart(INDX_CVICE,1) = -1.0176e-8
         polpart(INDX_CVICE,2) =  1.7615e-6
         polpart(INDX_CVICE,3) = -1.0480e-4
         polpart(INDX_CVICE,4) =  0.0019
         polpart(INDX_CVICE,5) =  0.0460
      endif
      if (ice_type.eq.1) then
         polpart(INDX_CVICE,1) = 1.3615e-8
         polpart(INDX_CVICE,2) = -2.04206e-6
         polpart(INDX_CVICE,3) = 7.51799e-5
         polpart(INDX_CVICE,4) = 0.00078213
         polpart(INDX_CVICE,5) = 0.0182131
      endif

! density:
!*    clear-sky air:
      rhoair = pres/(287.04*temp)

!*    liquid/ice particules:
      rhopart(INDX_LSLIQ) = rholiq
      rhopart(INDX_LSICE) = rhoice
      rhopart(INDX_CVLIQ) = rholiq
      rhopart(INDX_CVICE) = rhoice
      rhopart(INDX_LSSNOW) = rhoice/2. !+JEK+ZYY

! effective radius particles:
      rad_part(:,:,INDX_LSLIQ) = ls_radliq(:,:)
      rad_part(:,:,INDX_LSICE) = ls_radice(:,:)
      rad_part(:,:,INDX_CVLIQ) = cv_radliq(:,:)
      rad_part(:,:,INDX_CVICE) = cv_radice(:,:)
      rad_part(:,:,INDX_LSSNOW) = ls_radsnow(:,:) !+JEK+ZYY
      rad_part(:,:,:)=MAX(rad_part(:,:,:),0.)
      rad_part(:,:,:)=MIN(rad_part(:,:,:),70.0e-6)
      ls_radsnow(:,:)=MAX(ls_radsnow(:,:),0.) !+JEK+ZYY
      ls_radsnow(:,:)=MIN(ls_radsnow(:,:),1000.e-6) !+JEK+ZYY

! altitude at half pressure levels:
      zheight(:,1) = 0.0
      do k = 2, nlev+1
        zheight(:,k) = zheight(:,k-1) &
                  -(presf(:,k)-presf(:,k-1))/(rhoair(:,k-1)*9.81)
      enddo

!------------------------------------------------------------
!---- 2. Molecular alpha and beta:
!------------------------------------------------------------

      beta_mol = pres/km/temp * Cmol
      alpha_mol = 8.0*pi/3.0 * beta_mol

!------------------------------------------------------------
!---- 3. Particles alpha and beta:
!------------------------------------------------------------

! polynomes kp_lidar derived from Mie theory:
      do i = 1, npart
       where ( rad_part(:,:,i).gt.0.0)
         kp_part(:,:,i) = &
            polpart(i,1)*(rad_part(:,:,i)*1e6)**4 &
          + polpart(i,2)*(rad_part(:,:,i)*1e6)**3 &
          + polpart(i,3)*(rad_part(:,:,i)*1e6)**2 &
          + polpart(i,4)*(rad_part(:,:,i)*1e6) &
          + polpart(i,5)
        elsewhere
         kp_part(:,:,i) = 0.
        endwhere
      enddo
      
! mixing ratio particules in each subcolumn:
          qpart(:,:,INDX_LSLIQ) = q_lsliq(:,:) ! oct08
          qpart(:,:,INDX_LSICE) = q_lsice(:,:) ! oct08
          qpart(:,:,INDX_CVLIQ) = q_cvliq(:,:) ! oct08
          qpart(:,:,INDX_CVICE) = q_cvice(:,:) ! oct08
          qpart(:,:,INDX_LSSNOW) = q_lssnow(:,:) !+JEK+ZYY

! alpha of particles in each subcolumn:
      do i = 1, npart-1  !!+JEK+ZYY
        where ( rad_part(:,:,i).gt.0.0)
          alpha_part(:,:,i) = 3.0/4.0 * Qscat &
                 * rhoair(:,:) * qpart(:,:,i) &
                 / (rhopart(i) * rad_part(:,:,i) )
        elsewhere
          alpha_part(:,:,i) = 0.
        endwhere
      enddo

        where ( ls_radsnow(:,:).gt.0.0)      !+JEK+ZYY
          alpha_part(:,:,5) = 3.0/4.0 * Qscat & 
                 * rhoair(:,:) * qpart(:,:,5) &
                 / (rhopart(5) * ls_radsnow(:,:) ) !+JEK+ZYY
        elsewhere
          alpha_part(:,:,5) = 0. !+JEK+ZYY
        endwhere !+JEK+ZYY


!------------------------------------------------------------
!---- 4.1 Total Backscatter signal:
!------------------------------------------------------------

! optical thickness (molecular):
!     opt. thick of each layer
      tau_mol(:,1:nlev) = alpha_mol(:,1:nlev) &
         & *(zheight(:,2:nlev+1)-zheight(:,1:nlev))
!     opt. thick from TOA
      DO k = nlev-1, 1, -1
        tau_mol(:,k) = tau_mol(:,k) + tau_mol(:,k+1)
      ENDDO

! optical thickness (particles):

      tau_part = rdiffm * alpha_part
      DO i = 1, npart
!       opt. thick of each layer
        tau_part(:,:,i) = tau_part(:,:,i) &
           & * (zheight(:,2:nlev+1)-zheight(:,1:nlev) )
!       opt. thick from TOA
        DO k = nlev-1, 1, -1 
          tau_part(:,k,i) = tau_part(:,k,i) + tau_part(:,k+1,i)
        ENDDO
      ENDDO

! molecular signal:
!      Upper layer 
       pmol(:,nlev) = beta_mol(:,nlev) / (2.*tau_mol(:,nlev)) &
            & * (1.-exp(-2.0*tau_mol(:,nlev)))
!      Other layers
       DO k= nlev-1, 1, -1
        tau_mol_lay(:) = tau_mol(:,k)-tau_mol(:,k+1) ! opt. thick. of layer k
        WHERE (tau_mol_lay(:).GT.0.)
          pmol(:,k) = beta_mol(:,k) * EXP(-2.0*tau_mol(:,k+1)) / (2.*tau_mol_lay(:)) &
            & * (1.-exp(-2.0*tau_mol_lay(:)))
        ELSEWHERE
!         This must never happend, but just in case, to avoid div. by 0
          pmol(:,k) = beta_mol(:,k) * EXP(-2.0*tau_mol(:,k+1))
        END WHERE
      END DO

! Total signal (molecular + particules):
!
!
! For performance reason on vector computers, the 2 following lines should not be used
! and should be replace by the later one.
!      betatot(:,:) = beta_mol(:,:) + sum(kp_part*alpha_part,dim=3)
!      tautot(:,:)  = tau_mol(:,:)  + sum(tau_part,dim=3)
      betatot(:,:) = beta_mol(:,:)
      tautot(:,:)  = tau_mol(:,:)
      DO i = 1, npart
           betatot(:,:) = betatot(:,:) + kp_part(:,:,i)*alpha_part(:,:,i)
           tautot(:,:) = tautot(:,:)  + tau_part(:,:,i)
      ENDDO ! i
!
!     Upper layer 
      pnorm(:,nlev) = betatot(:,nlev) / (2.*tautot(:,nlev)) &
            & * (1.-exp(-2.0*tautot(:,nlev)))

!     Other layers
      DO k= nlev-1, 1, -1
          tautot_lay(:) = tautot(:,k)-tautot(:,k+1) ! optical thickness of layer k
        WHERE (tautot_lay(:).GT.0.)
          pnorm(:,k) = betatot(:,k) * EXP(-2.0*tautot(:,k+1)) / (2.*tautot_lay(:)) &
               & * (1.-EXP(-2.0*tautot_lay(:)))
        ELSEWHERE
!         This must never happend, but just in case, to avoid div. by 0
          pnorm(:,k) = betatot(:,k) * EXP(-2.0*tautot(:,k+1))
        END WHERE
      END DO

!------------------------------------------------------------
!---- 4.2 Ice/Liq Backscatter signal:
!------------------------------------------------------------

! Contribution of the molecular to beta
      betatot_ice(:,:) = beta_mol(:,:)
      betatot_liq(:,:) = beta_mol(:,:)

      tautot_ice(:,:) = tau_mol(:,:)
      tautot_liq(:,:) = tau_mol(:,:)

      DO i = 2, npart,2
           betatot_ice(:,:) = betatot_ice(:,:)+ kp_part(:,:,i)*alpha_part(:,:,i)
           tautot_ice(:,:) = tautot_ice(:,:)  + tau_part(:,:,i)
      ENDDO ! i
      DO i = 1, npart,2
           betatot_liq(:,:) = betatot_liq(:,:)+ kp_part(:,:,i)*alpha_part(:,:,i)
           tautot_liq(:,:) = tautot_liq(:,:)  + tau_part(:,:,i)
      ENDDO ! i


! Computation of the ice and liquid lidar backscattered signal (ATBice and ATBliq)
!     Ice only
!     Upper layer
      pnorm_ice(:,nlev) = betatot_ice(:,nlev) / (2.*tautot_ice(:,nlev)) &
            & * (1.-exp(-2.0*tautot_ice(:,nlev)))

      DO k= nlev-1, 1, -1
          tautot_lay_ice(:) = tautot_ice(:,k)-tautot_ice(:,k+1) 
        WHERE (tautot_lay_ice(:).GT.0.)
         pnorm_ice(:,k)=betatot_ice(:,k)*EXP(-2.0*tautot_ice(:,k+1))/(2.*tautot_lay_ice(:)) &
               & * (1.-EXP(-2.0*tautot_lay_ice(:)))
        ELSEWHERE
         pnorm_ice(:,k)=betatot_ice(:,k)*EXP(-2.0*tautot_ice(:,k+1))
        END WHERE
      ENDDO

!     Liquid only
!     Upper layer
      pnorm_liq(:,nlev) = betatot_liq(:,nlev) / (2.*tautot_liq(:,nlev)) &
            & * (1.-exp(-2.0*tautot_liq(:,nlev)))

      DO k= nlev-1, 1, -1
          tautot_lay_liq(:) = tautot_liq(:,k)-tautot_liq(:,k+1) 
        WHERE (tautot_lay_liq(:).GT.0.)
          pnorm_liq(:,k)=betatot_liq(:,k)*EXP(-2.0*tautot_liq(:,k+1))/(2.*tautot_lay_liq(:)) &
               & * (1.-EXP(-2.0*tautot_lay_liq(:)))
        ELSEWHERE
          pnorm_liq(:,k)=betatot_liq(:,k)*EXP(-2.0*tautot_liq(:,k+1))
        END WHERE
      ENDDO


! Computation of ATBperp,ice/liq from ATBice/liq including the multiple scattering 
! contribution (Cesana and Chepfer 2013, JGR)
!  ATBperp,ice = Alpha*ATBice 
!  ATBperp,liq = Beta*ATBliq^2 + Gamma*ATBliq

      DO k= nlev, 1, -1
	      pnorm_perp_ice(:,k) = Alpha * pnorm_ice(:,k) ! Ice particles
	      pnorm_perp_liq(:,k) = 1000*Beta * pnorm_liq(:,k)**2 + Gamma * pnorm_liq(:,k) ! Liquid particles
      ENDDO

! Computation of beta_perp_ice/liq using the lidar equation
!     Ice only
!     Upper layer 
      beta_perp_ice(:,nlev) = pnorm_perp_ice(:,nlev) * (2.*tautot_ice(:,nlev)) &
            & / (1.-exp(-2.0*tautot_ice(:,nlev)))

      DO k= nlev-1, 1, -1
         tautot_lay_ice(:) = tautot_ice(:,k)-tautot_ice(:,k+1)
         where (tautot_ice(:,k+1) .LT. TAUTOT_MAX)
            WHERE (tautot_lay_ice(:).GT.TAUTOT_SMALL)
               ! Need to compare against TAUTOT_SMALL because when
               ! tautot_lay_ice gets small enough, value of 
               ! [(1.-exp(-2.0*tautot_lay_ice(:))] approaches zero
               beta_perp_ice(:,k) = pnorm_perp_ice(:,k)/ EXP(-2.0*tautot_ice(:,k+1)) * (2.*tautot_lay_ice(:)) &
                                    & / (1.-exp(-2.0*tautot_lay_ice(:)))

            ELSEWHERE
               beta_perp_ice(:,k)=pnorm_perp_ice(:,k)/EXP(-2.0*tautot_ice(:,k+1))
            END WHERE
         else where
            ! If attenuation to TOA is too large, exp(-2.0*tautot_ice(:,k+1))
            ! goes to zero and the above code has a divide by zero error. A
            ! simple fix is to set the beta values of such layers to zero.
            beta_perp_ice(:,k) = 0
         end where
      ENDDO

!     Liquid only
!     Upper layer 
      beta_perp_liq(:,nlev) = pnorm_perp_liq(:,nlev) * (2.*tautot_liq(:,nlev)) &
            & / (1.-exp(-2.0*tautot_liq(:,nlev)))

      DO k= nlev-1, 1, -1
          tautot_lay_liq(:) = tautot_liq(:,k)-tautot_liq(:,k+1) 
          where (tautot_liq(:,k+1) .LT. TAUTOT_MAX) 
              WHERE (tautot_lay_liq(:).GT. TAUTOT_SMALL)
              ! Need to compare against TAUTOT_SMALL because when
              ! tautot_lay_liq gets small enough, value of 
              ! [(1.-exp(-2.0*tautot_lay_liq(:))] approaches zero
                 beta_perp_liq(:,k) = pnorm_perp_liq(:,k)/ EXP(-2.0*tautot_liq(:,k+1)) * (2.*tautot_lay_liq(:)) &
                                      & / (1.-exp(-2.0*tautot_lay_liq(:)))
   
              ELSEWHERE
                 beta_perp_liq(:,k)=pnorm_perp_liq(:,k)/EXP(-2.0*tautot_liq(:,k+1))
              END WHERE
           else where
              ! If attenuation to TOA is too large, exp(-2.0*tautot_liq(:,k+1))
              ! goes to zero and the above code has a divide by zero error. A
              ! simple fix is to set the beta values of such layers to zero.
              beta_perp_liq(:,k) = 0
           end where
      ENDDO



!------------------------------------------------------------
!---- 4.3 Perpendicular Backscatter signal:
!------------------------------------------------------------

! Computation of the total perpendicular lidar signal (ATBperp for liq+ice)
!     Upper layer 
    WHERE(tautot(:,nlev).GT.0)
          pnorm_perp_tot(:,nlev) = &
              (beta_perp_ice(:,nlev)+beta_perp_liq(:,nlev)-(beta_mol(:,nlev)/(1+1/0.0284))) / (2.*tautot(:,nlev)) &
              & * (1.-exp(-2.0*tautot(:,nlev)))
    ELSEWHERE
    pnorm_perp_tot(:,nlev) = 0.
    ENDWHERE

!     Other layers
      DO k= nlev-1, 1, -1
          tautot_lay(:) = tautot(:,k)-tautot(:,k+1) ! optical thickness of layer k

          ! The perpendicular component of the molecular backscattered signal (Betaperp) has been 
          ! taken into account two times (once for liquid and once for ice). 
          ! We remove one contribution using 
          ! Betaperp=beta_mol(:,k)/(1+1/0.0284)) [bodhaine et al. 1999] in the following equations:
            WHERE (pnorm(:,k).eq.0)
                  pnorm_perp_tot(:,k)=0.
                  ELSEWHERE
                    WHERE (tautot_lay(:).GT.0.)
                      pnorm_perp_tot(:,k) = &
                          (beta_perp_ice(:,k)+beta_perp_liq(:,k)-(beta_mol(:,k)/(1+1/0.0284))) * &
                          EXP(-2.0*tautot(:,k+1)) / (2.*tautot_lay(:)) &
                          & * (1.-EXP(-2.0*tautot_lay(:)))
                    ELSEWHERE
          !         This must never happen, but just in case, to avoid div. by 0
                      pnorm_perp_tot(:,k) = &
                           (beta_perp_ice(:,k)+beta_perp_liq(:,k)-(beta_mol(:,k)/(1+1/0.0284))) * &
                          EXP(-2.0*tautot(:,k+1))
                    END WHERE
            ENDWHERE

      END DO

!-------- End computation Lidar --------------------------

!---------------------------------------------------------
!  Parasol/Polder module
!
!  Purpose : Compute reflectance for one particular viewing direction
!  and 5 solar zenith angles (calculation valid only over ocean)
! ---------------------------------------------------------

! initialization:
    refl(:,:) = 0.0

! activate parasol calculations:
    if (ok_parasol) then

!     Optical thickness from TOA to surface
      tautot_S_liq = 0.
      tautot_S_ice = 0.
      tautot_S_liq(:) = tautot_S_liq(:) &
         + tau_part(:,1,1) + tau_part(:,1,3)
      tautot_S_ice(:) = tautot_S_ice(:) &
         + tau_part(:,1,2) + tau_part(:,1,4) &
         + tau_part(:,1,5) !+JEK+ZYY

      call parasol(npoints,nrefl,undef  &
                 ,tautot_S_liq,tautot_S_ice &
                 ,refl)

    endif ! ok_parasol

  END SUBROUTINE lidar_simulator
!
!---------------------------------------------------------------------------------
!
  SUBROUTINE parasol(npoints,nrefl,undef  &
                       ,tautot_S_liq,tautot_S_ice  &
                       ,refl)
!---------------------------------------------------------------------------------
! Purpose: To compute Parasol reflectance signal from model-simulated profiles 
!          of cloud water and cloud fraction in each sub-column of each model 
!          gridbox.
!
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - optimization for vectorization
!
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
!---------------------------------------------------------------------------------

    IMPLICIT NONE

! inputs
    INTEGER npoints              ! Number of horizontal gridpoints
    INTEGER nrefl                ! Number of angles for which the reflectance 
                                 ! is computed. Can not be greater then ntetas
    REAL undef                   ! Undefined value. Currently not used
    REAL tautot_S_liq(npoints)   ! liquid water cloud optical thickness, 
                                   ! integrated from TOA to surface
    REAL tautot_S_ice(npoints)   ! same for ice water clouds only
! outputs
    REAL refl(npoints,nrefl)     ! Parasol reflectances
!
! Local variables
    REAL tautot_S(npoints)       ! cloud optical thickness, from TOA to surface
    REAL frac_taucol_liq(npoints), frac_taucol_ice(npoints)

    REAL pi
!   look up table variables:
    INTEGER ny, it
    INTEGER ntetas, nbtau        ! number of angle and of optical thickness
                                   ! of the look-up table
    PARAMETER (ntetas=5, nbtau=7)
    REAL aa(ntetas,nbtau-1), ab(ntetas,nbtau-1)
    REAL ba(ntetas,nbtau-1), bb(ntetas,nbtau-1)  
    REAL tetas(ntetas),tau(nbtau)                        
    REAL r_norm(ntetas)
    REAL rlumA(ntetas,nbtau), rlumB(ntetas,nbtau)       
    REAL rlumA_mod(npoints,5), rlumB_mod(npoints,5) 

    DATA tau   /0., 1., 5., 10., 20., 50., 100./
    DATA tetas /0., 20., 40., 60., 80./
    
! Look-up table for spherical liquid particles:
    DATA (rlumA(1,ny),ny=1,nbtau) /0.03, 0.090886, 0.283965, &
     0.480587, 0.695235, 0.908229, 1.0 /
    DATA (rlumA(2,ny),ny=1,nbtau) /0.03, 0.072185, 0.252596, &
      0.436401,  0.631352, 0.823924, 0.909013 /
    DATA (rlumA(3,ny),ny=1,nbtau) /0.03, 0.058410, 0.224707, &
      0.367451,  0.509180, 0.648152, 0.709554 /
    DATA (rlumA(4,ny),ny=1,nbtau) /0.03, 0.052498, 0.175844, &
      0.252916,  0.326551, 0.398581, 0.430405 /
    DATA (rlumA(5,ny),ny=1,nbtau) /0.03, 0.034730, 0.064488, &
      0.081667,  0.098215, 0.114411, 0.121567 /

! Look-up table for ice particles:
    DATA (rlumB(1,ny),ny=1,nbtau) /0.03, 0.092170, 0.311941, &
       0.511298, 0.712079 , 0.898243 , 0.976646 /
    DATA (rlumB(2,ny),ny=1,nbtau) /0.03, 0.087082, 0.304293, &
       0.490879,  0.673565, 0.842026, 0.912966 /
    DATA (rlumB(3,ny),ny=1,nbtau) /0.03, 0.083325, 0.285193, &
      0.430266,  0.563747, 0.685773,  0.737154 /
    DATA (rlumB(4,ny),ny=1,nbtau) /0.03, 0.084935, 0.233450, &
      0.312280, 0.382376, 0.446371, 0.473317 /
    DATA (rlumB(5,ny),ny=1,nbtau) /0.03, 0.054157, 0.089911, &
      0.107854, 0.124127, 0.139004, 0.145269 /

!--------------------------------------------------------------------------------
! Lum_norm=f(tetaS,tau_cloud) derived from adding-doubling calculations
!        valid ONLY ABOVE OCEAN (albedo_sfce=5%)
!        valid only in one viewing direction (theta_v=30�, phi_s-phi_v=320�)
!        based on adding-doubling radiative transfer computation
!        for tau values (0 to 100) and for tetas values (0 to 80)
!        for 2 scattering phase functions: liquid spherical, ice non spherical

    IF ( nrefl.GT. ntetas ) THEN
        PRINT *,'Error in lidar_simulator, nrefl should be less then ',ntetas,' not',nrefl
        STOP
    ENDIF

    rlumA_mod=0
    rlumB_mod=0
!
    pi = ACOS(-1.0)
    r_norm(:)=1./ COS(pi/180.*tetas(:))
!
    tautot_S_liq(:)=MAX(tautot_S_liq(:),tau(1))
    tautot_S_ice(:)=MAX(tautot_S_ice(:),tau(1))
    tautot_S(:) = tautot_S_ice(:) + tautot_S_liq(:)
!
! relative fraction of the opt. thick due to liquid or ice clouds
    WHERE (tautot_S(:) .GT. 0.)
        frac_taucol_liq(:) = tautot_S_liq(:) / tautot_S(:)
        frac_taucol_ice(:) = tautot_S_ice(:) / tautot_S(:)
    ELSEWHERE
        frac_taucol_liq(:) = 1.
        frac_taucol_ice(:) = 0.
    END WHERE
    tautot_S(:)=MIN(tautot_S(:),tau(nbtau))
!
! Linear interpolation :

    DO ny=1,nbtau-1
! microphysics A (liquid clouds) 
      aA(:,ny) = (rlumA(:,ny+1)-rlumA(:,ny))/(tau(ny+1)-tau(ny))
      bA(:,ny) = rlumA(:,ny) - aA(:,ny)*tau(ny)
! microphysics B (ice clouds)
      aB(:,ny) = (rlumB(:,ny+1)-rlumB(:,ny))/(tau(ny+1)-tau(ny))
      bB(:,ny) = rlumB(:,ny) - aB(:,ny)*tau(ny)
    ENDDO
!
    DO it=1,ntetas
      DO ny=1,nbtau-1
        WHERE (tautot_S(:).GE.tau(ny).AND.tautot_S(:).LE.tau(ny+1))
            rlumA_mod(:,it) = aA(it,ny)*tautot_S(:) + bA(it,ny)
            rlumB_mod(:,it) = aB(it,ny)*tautot_S(:) + bB(it,ny)
        END WHERE
      END DO
    END DO
!
    DO it=1,ntetas
      refl(:,it) = frac_taucol_liq(:) * rlumA_mod(:,it) &
         + frac_taucol_ice(:) * rlumB_mod(:,it)
! normalized radiance -> reflectance: 
      refl(:,it) = refl(:,it) * r_norm(it)
    ENDDO

    RETURN
  END SUBROUTINE parasol
