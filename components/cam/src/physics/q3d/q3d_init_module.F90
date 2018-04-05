MODULE q3d_init_module
! Contains the programs related to initializing the channels

      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t

      USE parmsld, only: channel_l,channel_w,nk1,nk2,nk3,ntracer,nhalo,nhalo_adv
      USE mphysics_variables, only: rho_int,rhol,dzl,pl0,pil0 
      
      USE constld                    

!JUNG #ifdef FIXTHIS
! Need to be able to compile all of RRTMG
      USE trace_gases, only: trace_gas_input
!JUNG #endif
      
      USE rrtm_grid,   only: nrad_rad => nrad 
      USE rrtm_vars,   only: pres,presi

      USE bound_channel_module, only: bound_channel
      USE bound_extra,          only: bound_normal,bound_vert,bound_sync
      USE halo_q,               only: halo_correc_q
      USE halo_vort,            only: vort_comm_pre,vort_comm_post
      USE halo_z,               only: halo_correc_z     
      USE vort_3d_module,       only: zeta_diag
      USE wind_module,          only: wind_3d
      USE z_coord,              only: coords_2d
      
IMPLICIT NONE
PRIVATE

PUBLIC :: init_common,init_channel

CONTAINS

!===================================================================================
   SUBROUTINE INIT_common (isRestart) 
!===================================================================================
!  Initialize the common parameters and vertical profiles

      logical, intent(in) :: isRestart  ! .true. iff restart run
      
      ! Local
      INTEGER :: K
      
      REAL (kind=r8) :: PZERO = 100000._r8     ! [Pa]
      REAL (kind=r8) :: RBCP,CPBR,ALPHAW
      REAL (KIND=r8), DIMENSION(NK3) :: TV,ALPHA
      
      IF (.NOT. isRestart) THEN
      
!***********************************************************************************
!                             SET DEFAULT VALUES in CONSTLD
!***********************************************************************************

!------------------------
! Logical Parameters
!------------------------
      NOSFX      = .TRUE.               ! T, no surface flux
      NOTOPO     = .TRUE.               ! T, no topography
      NOMAP      = .FALSE.              ! T, no (horizontal) mapping
      ZCONST     = .TRUE.               ! T, use of constant-depth vertical grids  
      PHYSICS    = .FALSE.              ! T, include physical processes
      TURBULENCE = .FALSE.              ! T, include turbulence
      RADCODE    = .FALSE.              ! T, include radiation (RRTMG)  
      MICROCODE  = .FALSE.              ! T, include microphysics
      BUOY       = .TRUE.               ! T, include buoyancy
      NOTHERM    = .FALSE.              ! T, no thermodynamical processes

!------------------------      
! REAL Parameters
!------------------------
      DT      = 10.0_r8                 ! integration timestep of crm (s)
      DT_GCM  = 300.0_r8                ! dt of gcm (s) 
      
      DZ      = 1000.0_r8               ! grid size in vertical direction (m)
      DZSQ    = DZ * DZ

!      (Specified in Mapping_channel)
!      DX      =                        ! crm grid size in x-direction (m) 
!      DY      =                        ! crm grid size in y-direction (m)  
!      DXSQ    =
!      DYSQ    =
!      DXDY    =  
!      DX_GCM  =                        ! gcm grid size in x-direction (m)
!      DY_GCM  =                        ! gcm grid size in y-direction (m) 
      
      DZ1     = 500.0_r8                ! the lowest layer depth in z (m)  
      DOMAIN  = 30000.0_r8              ! vertical domain forstretched grids
      ZB      = 0.0_r8                  ! surface height (m)
      
      PI      = 4.0_r8 * DATAN(1.0_r8)  ! 3.14159... 
      GRAV    = 9.80616_r8              ! gravitational constant (m/s^2)
      HLF     = 2.500D+6                ! latent heat of vaporization (J/kg)
      HLM     = 3.336D+5                ! latent heat of fusion (J/kg) 
      CP      = 1004.5_r8               ! specific heat of dry air at constant pressure (J/kg/K)
      GAS_CONST = 287.0_r8              ! dry air gas constant
      DELTA   = 0.608_r8                ! 0.608
      VK      = 0.4_r8                  ! von Karman constant
      
      ZRSEA   = 2.0D-04                 ! surface roughness over sea (m)
      ZRLAND  = 1.0D-02                 ! surface roughness over land (m)
      
      SST     = 298.0_r8                ! sea surface temperature (K)

      PSFC    = 100000.0_r8             ! surface pressure (Pa)

      REARTH  = 6371220.0_r8            ! radius of the earth (m)
      OMEGA   = 7.292D-5                ! earth rotation rate (s^-1)  
      RAD2DEG = 180.0_r8/PI             ! deg. of 1 rad (i.e., 180/pi) 

      ALADV   = 1.0_r8                  ! alpha-value in advection

!      (Specified in Mapping_channel, after dx & dy are determined)      
!      WRXMU   =                        ! (mu/dt) in the relaxation method 
           
      DTRAD   = 3600.0_r8               ! time interval between radiation calls (s) 
           
      CRAD    = 1.0_r8/(3600._r8)       ! Rayleigh-type damping coefficient (s^-1)
      GWDHT   = 15000.0_r8              ! GWD above this height (m)
      
      DTHMAX  = 0.5_r8                  ! coefficient for random perturbation (K)
      Z1PERT  = 0.0_r8                  ! z_low for random perturbation  (m)
      Z2PERT  = 12500.0_r8              ! z_high for random perturbation (m)

!------------------------          
! Integer Parameters
!------------------------
      ITT_MAX = INT(dt_gcm/dt)          ! # of CRM time steps per one GCM time step 

! The below may not be needed.      
!!      NRESTART = 999999               ! restart writing interval
!!      NXS      = 2                    ! accumulation frequency for averaging
!!      NXSAVG   = 2                    ! averaging period for output 
!!      NFLPRT   = 999999               ! frequency of writing the full output
!!      NOUT     = 12                   ! # of data outputs

      ITINIT  = 0                       ! starting time for random perturbation 
      ITSTOP  = 0                       ! ending time for random perturbation
                       
      NSFLUX  = 6                       ! frequency of surface flux calculation
      NRAD    = INT (DTRAD/DT)          ! frequency of radiation calculation
      
      NITERW  = 10                      ! # of iterations for w
      NITERXY = 50                      ! # of iterations for psi and chi
      
!***********************************************************************************
!                        Define vertical grids: zt, zz, fnt, fnz 
!***********************************************************************************
      IF (ZCONST) THEN
!       Use of constant-depth vertical grids
        CZ2 = 0.0_r8
        CZ1 = 1.0_r8
      ELSE
!       Use of stretched vertical grids
        CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
        CZ1 = 1.0_r8 - CZ2 * DOMAIN
      ENDIF

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

!***********************************************************************************
! Define main vertical profiles: rho, rhoz, pbar, pibar, fn1, fn2, thbar, qvbar
! Need more work: how to prescribe thbar and pibar?
!***********************************************************************************
      DO K = 1, NK3
        QVBAR(K) = 0.0_r8
      ENDDO

      CPBR = CP / gas_const
      RBCP = 1.0_r8 / CPBR

      PISFC = ( PSFC / PZERO ) ** RBCP
      PBAR(1)  = PSFC   ! [Pa]

!     Prepare THBAR & PIBAR: temparary setting
      CALL PRESCRIBE_IC (THBAR,PIBAR)

      TV(1) = THBAR(1)*PIBAR(1)
      DO K = 2, NK3
       PBAR(K) = PZERO * PIBAR(K) ** CPBR
       TV(K)   = THBAR(K)*PIBAR(K)
      ENDDO
      
      PIBARZ(1) = PIBAR(1)
      PBARZ(1)  = PBAR(1)

      DO K = 2, nk2
        PIBARZ(K)= PIBAR(K) + (PIBAR(K+1)-PIBAR(K))/(ZT(K+1)-ZT(K)) * (ZZ(K)-ZT(K))
      ENDDO
      DO K = 2, nk2
        PBARZ(K) = PSFC * PIBARZ(K) ** ( cp / gas_const )
      ENDDO      

!     PROFILES OF RHO, RHOZ
      DO K = 2, NK3
        ALPHA(K) = gas_const * TV(K) / PBAR(K)
      ENDDO

      DO K = 2, NK2
       ALPHAW = (ALPHA(K) + ALPHA(K+1)) / 2.0_r8
       RHOZ(K) = 1.0_r8 / ALPHAW
      ENDDO

!     SURFACE DENSITY
      RHOZ(1) = PSFC / (gas_const * TV(1))

!     DENSITY FOR k=1/2
      RHO(1) = RHOZ(1)
      DO K = 2, NK3
        RHO(K) = 1.0_r8 / ALPHA(K)
      ENDDO

! TEMP----------------
!      DO K = 1, NK3
!        RHO(K) = 1.
!      ENDDO
!      DO K = 1, NK2
!        RHOZ(K) = 1.
!      ENDDO
! TEMP----------------

      DO K=2,NK2
       FN1(K)=RHO(K+1)*FNZ(K)/FNT(K+1)  ! map factor ratio, used in vort_3d_module
       FN2(K)=RHO(K)*FNZ(K)/FNT(K)      ! map factor ratio, used in vort_3d_module
      ENDDO 

!     Find the vertical index of the lowest layer where GWD is applied.
      IF (CRAD .NE. 0.0_r8) THEN
        K = 2
        DO WHILE (ZT(K) .LT. GWDHT)
        K = K + 1
        ENDDO
        KLOW_GWD = K            
      ENDIF

    ENDIF   ! .NOT. isRestart 

!***********************************************************************************
!                 Define parameters and vertical profiles used in PHYSICS
!***********************************************************************************
    
!--------------------------
    IF (physics) THEN
!--------------------------    

      ! Initialize Parameters for Microphysics    
      CALL initialize_microphysics  
      
      DO K=1,nk1
        RHOL(K) = RHO(K+1)          
        PL0(K)  = PBAR(K+1)
        PIL0(K) = PIBAR(K+1)
      ENDDO
      
      DO K=1,nk1
        DZL(K) = ZZ(K+1) - ZZ(K)
      ENDDO
      
      RHO_INT(0)   = RHOL(1) 
      RHO_INT(nk1) = RHOL(nk1)
      DO K=1,nk1-1   
        RHO_INT(K) = &
          (DZL(K)*RHOL(K) + DZL(K+1)*RHOL(K+1))/(2.0_r8*(ZT(K+2) - ZT(K+1)))
      ENDDO
       
      !-------------------- 
      IF (radcode) THEN
      !--------------------
      ! Initialize Parameters for Radiation   
      
      ! Get vertical grid parameters and define vertical pressure levels
      PRES(:)  = PBAR(2:nk2) * 1.0D-2    ! [mb]
      PRESI(:) = PBARZ(:) * 1.0D-2

      NRAD_rad = NRAD

      CALL TRACE_GAS_INPUT(channel_l, channel_w, nk1, PBAR(2:nk2), PBARZ)
      ! Horizontally uniform profiles are used.      

      !--------------------
      ENDIF  ! RADCODE
      !--------------------
      
!--------------------------      
    ENDIF  ! PHYSICS   
!--------------------------    

   CONTAINS
   
   SUBROUTINE PRESCRIBE_IC (THBAR_LO,PIBAR_LO)
   IMPLICIT NONE
!  Just temporary use.    

      REAL (KIND=r8),DIMENSION(nk3),INTENT(OUT) :: THBAR_LO,PIBAR_LO

!     Vertical domain size = 25 km      
      INTEGER, PARAMETER :: KM = 2500
      INTEGER, PARAMETER :: KM_TROP_THETA = 900,      &
                            KM1 = KM_TROP_THETA + 30, &
                            KM2 = KM1 + 30,           &
                            KM3 = KM2 + 750
      
      REAL (kind=r8), PARAMETER :: NT_SQUARE  = 0.00015_r8, &
                                   NT1_SQUARE = 0.00050_r8, &
                                   NT2_SQUARE = 0.00055_r8, &
                                   NT3_SQUARE = 0.00055_r8, &
                                   NS_SQUARE  = 0.00080_r8

      REAL (kind=r8), PARAMETER :: THETA00 = 300.0_r8    
      REAL (kind=r8), PARAMETER :: ZS = 0.0_r8, DZ_LO = 10.0_r8
      REAL (kind=r8), PARAMETER :: EXF_S = 1.0_r8     
!     Assume: p=1000mb at zs
         
      REAL (kind=r8), DIMENSION(0:KM) :: Z,THETA0,THETA_BAR,EXF,TEMP
      
      INTEGER :: LGAP_Z,LOC_K(NK2)
      INTEGER :: ITERATION,K,KLO
      
!--------------------------------------      
!     Set the local vertical grids
!--------------------------------------
      Z(0) = ZS
      DO K=1,KM
       Z(K) = Z(K-1) + DZ_LO
      ENDDO

!-----------------------------------------      
!     Set the vertical profile of theta0
!-----------------------------------------
      DO K=0,KM_TROP_THETA
         THETA0(K) = THETA00*DEXP(NT_SQUARE*(Z(K)-ZS)/GRAV)
      ENDDO

      DO K=KM_TROP_THETA+1,KM1
         THETA0(K) = THETA0(KM_TROP_THETA)*DEXP(NT1_SQUARE*(Z(K)-Z(KM_TROP_THETA))/GRAV)
      ENDDO

      DO K=KM1+1,KM2
         THETA0(K) = THETA0(KM1)*DEXP(NT2_SQUARE*(Z(K)-Z(KM1))/GRAV)
      ENDDO

      DO K=KM2+1,KM3
         THETA0(K) = THETA0(KM2)*DEXP(NT3_SQUARE*(Z(K)-Z(KM2))/GRAV)
      ENDDO

      DO K=KM3+1,KM
         THETA0(K) = THETA0(KM3)*DEXP(NS_SQUARE*(Z(K)-Z(KM3))/GRAV)
      ENDDO

!     Vertical Smoothing 
      DO ITERATION=1,1000
       TEMP = THETA0
       DO K=2,KM-2
        TEMP(K) = &
        (THETA0(K-2)+THETA0(K-1)+THETA0(K)+THETA0(K+1)+THETA0(K+2))/5.0_r8
       ENDDO
       THETA0 = TEMP
      ENDDO
      
!-----------------------------------------      
!     Set up theta
!-----------------------------------------
      DO K=0,KM
       THETA_BAR(K) = THETA0(K)
      ENDDO
      
      EXF(0) = EXF_S
      DO K=1,KM
       EXF(K) = THETA_BAR(K-1)*EXF(K-1)/THETA_BAR(K)
      ENDDO      
!---------------------------------------------------------      
!     PREPARE IC (following the resolution of the model)
!---------------------------------------------------------    
      LGAP_Z  = IDINT(DZ/DZ_LO)
       
      LOC_K(1) = 0
      LOC_K(2) = IDINT(0.5_r8*DZ/DZ_LO)
      DO K = 3,NK2
       LOC_K(K) = LOC_K(K-1) + LGAP_Z 
      ENDDO
      
!=====================      
      DO K = 2,NK2
!=====================       
       KLO  = LOC_K(K) 
       THBAR_LO(K) = THETA_BAR(KLO)
       PIBAR_LO(K) = EXF(KLO)
!=====================               
      ENDDO  ! K-loop 
!=====================       
      THBAR_LO(1) = THETA_BAR(0)
      PIBAR_LO(1) = EXF(0)     
       
!      THBAR_LO(1) = THBAR_LO(2)
!      PIBAR_LO(1) = PIBAR_LO(2)
            
      THBAR_LO(NK3) = THBAR_LO(NK2)
      PIBAR_LO(NK3) = PIBAR_LO(NK2)
             
   END SUBROUTINE prescribe_ic 
   
   END SUBROUTINE init_common

!===================================================================================
   SUBROUTINE INIT_channel (channel)
!===================================================================================
!  Initialize the variables of channel segments

      type(channel_t), intent(inout) :: channel  ! channel data

      ! Local
      REAL (KIND=r8) :: ZRSEA = 2.0D-04    ! temporary setting 
      
      INTEGER :: I,J,K,nt,num_seg,mi1,mj1,mim,mip,mjm,mjp,mip_c,mjp_c

! Calculate various coefficients    
!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 


      ! Calculate the coefficients used in turbulence and diffusion      
      DO J = 1, mj1
      DO I = 1, mi1  
      channel%seg(num_seg)%COEFX(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXDY*4.0_r8)
      
      channel%seg(num_seg)%COEFX_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXDY*4.0_r8)

      channel%seg(num_seg)%COEFX_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXDY*4.0_r8)

      channel%seg(num_seg)%COEFX_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXDY*4.0_r8)      
      ENDDO
      ENDDO

      ! Calculate the coefficients used in relaxation_3d     
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1 
          channel%seg(num_seg)%AGAU_CO(I,J,K-1) = -channel%seg(num_seg)%RG_T(I,J) &
                                                  *FNZ(K)*FNT(K)/(DZSQ*RHO(K))
                 
          channel%seg(num_seg)%BGAU_CO(I,J,K-1) = &
                     (WRXMU+(channel%seg(num_seg)%RGG_U(1,I,J)                    &
                            +channel%seg(num_seg)%RGG_U(1,I-1,J))/DXSQ            &
                           +(channel%seg(num_seg)%RGG_V(3,I,J)                    &
                            +channel%seg(num_seg)%RGG_V(3,I,J-1))/DYSQ)/RHOZ(K)   &
            +channel%seg(num_seg)%RG_T(I,J)*FNZ(K)*(FNT(K+1)/RHO(K+1)+FNT(K)/RHO(K))/DZSQ
           
          channel%seg(num_seg)%CGAU_CO(I,J,K-1) = -channel%seg(num_seg)%RG_T(I,J) &
                                                  *FNZ(K)*FNT(K+1)/(DZSQ*RHO(K+1))
        ENDDO
       ENDDO  
      ENDDO      

      ! Calculate the coefficients used in relaxation_2d        
      DO J = 1, mj1
       DO I = 1, mi1 
          channel%seg(num_seg)%RX2D_C1(I,J) = (channel%seg(num_seg)%RGG_V(1,I,J)         &
                                              +channel%seg(num_seg)%RGG_V(1,I+1,J))/DXSQ &
                                            + (channel%seg(num_seg)%RGG_U(3,I,J)         &
                                              +channel%seg(num_seg)%RGG_U(3,I,J+1))/DYSQ
                 
          channel%seg(num_seg)%RX2D_C2(I,J) = (channel%seg(num_seg)%RGG_U(1,I,J)         &
                                              +channel%seg(num_seg)%RGG_U(1,I-1,J))/DXSQ &
                                            + (channel%seg(num_seg)%RGG_V(3,I,J)         &
                                              +channel%seg(num_seg)%RGG_V(3,I,J-1))/DYSQ
       ENDDO  
      ENDDO   
!******************************  
      ENDDO   ! num_seg 
!******************************   

! Prepare characteristics of topography  
      CALL TOPO_PREPARE (channel)

! Initialize basic segment data 
!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

      mim = channel%seg(num_seg)%mim  
      mip = channel%seg(num_seg)%mip
      
      mjm = channel%seg(num_seg)%mjm
      mjp = channel%seg(num_seg)%mjp
      
      mip_c = channel%seg(num_seg)%mip_c
      mjp_c = channel%seg(num_seg)%mjp_c  

! Surface fluxes: initialization is not really necessary.
      IF (PHYSICS) THEN     
      DO J = 1,mj1
       DO I = 1,mi1
        channel%seg(num_seg)%SPREC(I,J) = 0.0_r8
        
        channel%seg(num_seg)%WTH(I,J) = 0.0_r8
        channel%seg(num_seg)%WQV(I,J) = 0.0_r8
       ENDDO
      ENDDO   

      DO J = 1,mjp_c
       DO I = 1,mip_c
        channel%seg(num_seg)%UW(I,J) = 0.0_r8
        channel%seg(num_seg)%WV(I,J) = 0.0_r8
        
        channel%seg(num_seg)%UW_CON(I,J) = 0.0_r8
        channel%seg(num_seg)%WV_CON(I,J) = 0.0_r8
       ENDDO
      ENDDO         
      ENDIF     

! Surface information     
!-------------------------------------------------------------------
!      LOCEAN is used in the calculation of surface fluxes. 
!      GWET is used in the calculation of surface fluxes. 
!      TG is used in radiation and surface fluxes.
!      ZROUGH is used in the calculation of turbulence coefficient.
!-------------------------------------------------------------------
      DO J = 1, mjp_c
       DO I = 1, mip_c 
        channel%seg(num_seg)%LOCEAN(I,J) = .TRUE.
       ENDDO
      ENDDO  

      DO J = mjm, mjp
       DO I = mim, mip 
        channel%seg(num_seg)%TG(I,J) = SST 
        channel%seg(num_seg)%ZROUGH(I,J) = ZRSEA
        channel%seg(num_seg)%GWET(I,J) = 1.0_r8
       ENDDO
      ENDDO  
!      CALL BOUND_NORMAL  (nhalo,channel,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)      
!      CALL BOUND_CHANNEL (nhalo_adv,channel,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)
!      IF (.not.nomap) &
!      CALL HALO_CORREC_Q (nhalo_adv,TG=.TRUE.,ZROUGH=.TRUE.,GWET=.TRUE.)

! Thermodynamic Fields (Temporary setting: zero deviation)
      
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D_bg(I,J,K)
         channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D_bg(I,J,K)
        ENDDO 
       ENDDO
      ENDDO

      IF (PHYSICS) THEN 
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D_bg(I,J,K)
         channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D_bg(I,J,K)
         channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D_bg(I,J,K)
         channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D_bg(I,J,K)
         channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D_bg(I,J,K)
        ENDDO 
       ENDDO
      ENDDO
      ENDIF

      DO NT = 1, ntracer 
      DO K = 2, NK2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%QT3D(I,J,K,nt) = channel%seg(num_seg)%QT3D_bg(I,J,K,nt)
        ENDDO 
       ENDDO
      ENDDO
      ENDDO
            
!------------------------
      IF (PHYSICS) THEN
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
                          
      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE., &
                       QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
!------------------------
      ELSE
!------------------------
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE.,QT3D=.TRUE.)
!------------------------
      ENDIF  ! PHYSICS
!------------------------            

! Vorticity Components (Temporary setting: zero deviation)

      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX_bg(I,J,K)
         channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)
      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
            
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ(I,J,nk2) = channel%seg(num_seg)%Z3DZ_bg(I,J,1)
       ENDDO
      ENDDO
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)   
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.)

      CALL ZETA_DIAG (channel)

! Diagnose wind components

      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%W3D(I,J,  1) = 0.0_r8
        channel%seg(num_seg)%W3D(I,J,nk2) = 0.0_r8
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%W3D(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      CALL BOUND_NORMAL  (nhalo,channel,W3D=.TRUE.)
      CALL BOUND_CHANNEL (nhalo,channel,W3D=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Q (nhalo,channel,W3D=.TRUE.) 

      DO K = 1, NK2
       DO J = mjm, mjp
        DO I = mim, mip
         channel%seg(num_seg)%W3Dnm1(I,J,K) = channel%seg(num_seg)%W3D(I,J,K)
        ENDDO
       ENDDO
      ENDDO 

      ! PSI & CHI are only for deviation parts      
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSI(I,J) = 0.0_r8
         channel%seg(num_seg)%CHI(I,J) = 0.0_r8
       ENDDO
      ENDDO  
      
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSInm1(I,J) = channel%seg(num_seg)%PSI(I,J)
         channel%seg(num_seg)%CHInm1(I,J) = channel%seg(num_seg)%CHI(I,J)
       ENDDO
      ENDDO                 
           
      CALL WIND_3D  (.TRUE., .TRUE., channel, NITER2D=5000)  
      
!******************************  
      ENDDO   ! num_seg 
!******************************   
      
      CONTAINS

!==========================================      
      SUBROUTINE TOPO_PREPARE (channel)
!==========================================           
!     INPUT : topographic data (height)
!     OUTPUT: associated parameters  
!
!     For a while, "NOTOPO=.T." setting is used because there is no input.
     
      type(channel_t), intent(inout) :: channel  ! channel data

      ! Local variables
      INTEGER :: mi1,mj1,mim,mim_a,mim_c,mip,mip_a,mjm,mjm_a,mjm_c,mjp,mjp_a
      INTEGER :: I,J,K

      INTEGER,DIMENSION(4)    :: KLOWQ_MAX
      INTEGER,DIMENSION(nk2,4):: NTOPOQ,NTOPOU,NTOPOV
      
      REAL (kind=r8) :: HEIGHT

!*************************************************************************
      IF (NOTOPO) then
!*************************************************************************
!===============================    
      DO num_seg = 1, 4
!=============================== 
      mim   = channel%seg(num_seg)%mim    
      mip   = channel%seg(num_seg)%mip  
      mim_a = channel%seg(num_seg)%mim_a   
      mip_a = channel%seg(num_seg)%mip_a  
    
      mjm   = channel%seg(num_seg)%mjm      
      mjp   = channel%seg(num_seg)%mjp    
      mjm_a = channel%seg(num_seg)%mjm_a  
      mjp_a = channel%seg(num_seg)%mjp_a         

      ! HX: level index of topography
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%HX(I,J) = 1
       ENDDO
      ENDDO

      ! KLOWQ: lowest air layer for Q
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%KLOWQ_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWU: lowest air level for U, ETA
      DO J = mjm_a, mjp_a
       DO I = mim, mip_a
        channel%seg(num_seg)%KLOWU_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWV: lowest air level for V, KSI
      DO J = mjm, mjp_a
       DO I = mim_a, mip_a
        channel%seg(num_seg)%KLOWV_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO

      channel%seg(num_seg)%KLOWQ_MAX_GLOB = 2

!-----------------------------------------
! For vorticity advection with mountain
!-----------------------------------------
      channel%seg(num_seg)%DM_KSI_XE = 0
      channel%seg(num_seg)%DM_KSI_XC = 0
      channel%seg(num_seg)%DM_KSI_YS = 0

      channel%seg(num_seg)%DM_ETA_XW = 0
      channel%seg(num_seg)%DM_ETA_YN = 0
      channel%seg(num_seg)%DM_ETA_YC = 0
!--------------------------------
! For turbulence with mountain
!--------------------------------
      channel%seg(num_seg)%DM_Q_TE = 0
      channel%seg(num_seg)%DM_Q_TW = 0
      channel%seg(num_seg)%DM_Q_TN = 0
      channel%seg(num_seg)%DM_Q_TS = 0

      channel%seg(num_seg)%DM_KSI_TE = 0
      channel%seg(num_seg)%DM_KSI_TW = 0
      channel%seg(num_seg)%DM_KSI_TN = 0
      channel%seg(num_seg)%DM_KSI_TS = 0

      channel%seg(num_seg)%DM_ETA_TE = 0
      channel%seg(num_seg)%DM_ETA_TW = 0
      channel%seg(num_seg)%DM_ETA_TN = 0
      channel%seg(num_seg)%DM_ETA_TS = 0

!===============================  
      ENDDO   ! num_seg 
!===============================   

!*************************************************************************
      ELSE   ! NOTOPO
!*************************************************************************
! Prepare information associated with the topographic data

!===============================   
      DO num_seg = 1, 4
!=============================== 
      mim   = channel%seg(num_seg)%mim    
      mip   = channel%seg(num_seg)%mip  
      mim_a = channel%seg(num_seg)%mim_a   
      mip_a = channel%seg(num_seg)%mip_a  
    
      mjm   = channel%seg(num_seg)%mjm      
      mjp   = channel%seg(num_seg)%mjp    
      mjm_a = channel%seg(num_seg)%mjm_a  
      mjp_a = channel%seg(num_seg)%mjp_a      
           
      ! Read Orographic Data (from INPUT): currently, specify 0.0 km       
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%TOPOZ(I,J) = 0.0_r8
       ENDDO 
      ENDDO 

      ! Halo filling for TOPO
!-----------------------------------------------
!      CALL BOUND_NORMAL  (nhalo,channel,TOPOZ=.TRUE.)
!      CALL BOUND_CHANNEL (nhalo,channel,TOPOZ=.TRUE.)
!
!      IF (.not.nomap) &
!      CALL HALO_CORREC_Q (nhalo,channel,TOPOZ=.TRUE.)
!-----------------------------------------------

      ! HX (i,j): level index of mountain top (at cell center)
      ! HX=1: no mountain   e.g., HX=4: mountain top is at k=4
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%HX(I,J) = 1
       ENDDO
      ENDDO

      DO J = mjm, mjp
       DO I = mim, mip
        HEIGHT = channel%seg(num_seg)%TOPOZ(I,J)

        IF (HEIGHT.GT.ZZ(1)) THEN
          K = 2
          DO WHILE (HEIGHT.GE.ZZ(K))
           K=K+1
          ENDDO
          channel%seg(num_seg)%HX(I,J) = &
                 (K-1)+NINT((HEIGHT-ZZ(K-1))/(ZZ(K)-ZZ(K-1)))
        ENDIF
       ENDDO
      ENDDO

      ! KLOWQ: lowest air layer for Q
      DO J = mjm, mjp
       DO I = mim, mip
        channel%seg(num_seg)%KLOWQ_IJ(I,J) = channel%seg(num_seg)%HX(I,J)+1
       ENDDO
      ENDDO
      ! KLOWU: lowest air level for U, ETA
      DO J = mjm_a, mjp_a
       DO I = mim, mip_a
        channel%seg(num_seg)%KLOWU_IJ(I,J) = &
        MAX0((channel%seg(num_seg)%HX(I,J)+1),(channel%seg(num_seg)%HX(I+1,J)+1))
       ENDDO
      ENDDO
      ! KLOWV: lowest air level for V, KSI
      DO J = mjm, mjp_a
       DO I = mim_a, mip_a
        channel%seg(num_seg)%KLOWV_IJ(I,J) = &
        MAX0((channel%seg(num_seg)%HX(I,J)+1),(channel%seg(num_seg)%HX(I,J+1)+1))
       ENDDO
      ENDDO

      KLOWQ_MAX(num_seg) = MAXVAL(channel%seg(num_seg)%KLOWQ_IJ(1:mi1,1:mj1))

!=============================== 
      ENDDO   ! num_seg 
!===============================   

      channel%seg(:)%KLOWQ_MAX_GLOB = MAXVAL(KLOWQ_MAX(1:4))
       
!===============================   
      DO num_seg = 1, 4
!===============================
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      mim_c = channel%seg(num_seg)%mim_c  
      mjm_c = channel%seg(num_seg)%mjm_c      
      
!---------------------------------------- 
! For vorticity advection with mountain
!---------------------------------------- 
      ! Variables are defined at the vorticity Flux locations
      channel%seg(num_seg)%DM_KSI_XE = 0
      channel%seg(num_seg)%DM_KSI_XC = 0
      channel%seg(num_seg)%DM_KSI_YS = 0

      channel%seg(num_seg)%DM_ETA_XW = 0
      channel%seg(num_seg)%DM_ETA_YN = 0
      channel%seg(num_seg)%DM_ETA_YC = 0

      ! Used in meridional advection for ksi & Zonal advection for eta
      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I,J)) channel%seg(num_seg)%DM_KSI_YS(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J)) channel%seg(num_seg)%DM_ETA_XW(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      ! Used in zonal advection for ksi & Meridional advection for eta
      DO J = mim_c, mj1
       DO I = mim_c, mi1
        DO K=2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I,J).OR.K.LE.channel%seg(num_seg)%HX(I,J+1))     &
                                      channel%seg(num_seg)%DM_KSI_XC(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I+2,J).OR.K.LE.channel%seg(num_seg)%HX(I+2,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_XE(I,J,K) = 1

         IF (K.LE.channel%seg(num_seg)%HX(I,J).OR.K.LE.channel%seg(num_seg)%HX(I+1,J))     &
                                      channel%seg(num_seg)%DM_ETA_YC(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+2).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+2)) &
                                      channel%seg(num_seg)%DM_ETA_YN(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

!---------------------------------------------------------------------
! For turbulence with mountain
!---------------------------------------------------------------------
      channel%seg(num_seg)%DM_Q_TE = 0
      channel%seg(num_seg)%DM_Q_TW = 0
      channel%seg(num_seg)%DM_Q_TN = 0
      channel%seg(num_seg)%DM_Q_TS = 0

      channel%seg(num_seg)%DM_KSI_TE = 0
      channel%seg(num_seg)%DM_KSI_TW = 0
      channel%seg(num_seg)%DM_KSI_TN = 0
      channel%seg(num_seg)%DM_KSI_TS = 0

      channel%seg(num_seg)%DM_ETA_TE = 0
      channel%seg(num_seg)%DM_ETA_TW = 0
      channel%seg(num_seg)%DM_ETA_TN = 0
      channel%seg(num_seg)%DM_ETA_TS = 0

      DO J = 1, mj1
       DO I = 1, mi1
        DO K=2,channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+1,J)) channel%seg(num_seg)%DM_Q_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J)) channel%seg(num_seg)%DM_Q_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+1)) channel%seg(num_seg)%DM_Q_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1)) channel%seg(num_seg)%DM_Q_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2, channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+1,J).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J).OR.K.LE.channel%seg(num_seg)%HX(I-1,J+1)) &
                                      channel%seg(num_seg)%DM_KSI_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+2)) channel%seg(num_seg)%DM_KSI_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1)) channel%seg(num_seg)%DM_KSI_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO J = mjm_c, mj1
       DO I = mim_c, mi1
        DO K = 2, channel%seg(num_seg)%KLOWQ_MAX_GLOB-1
         IF (K.LE.channel%seg(num_seg)%HX(I+2,J)) channel%seg(num_seg)%DM_ETA_TE(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I-1,J)) channel%seg(num_seg)%DM_ETA_TW(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J+1).OR.K.LE.channel%seg(num_seg)%HX(I+1,J+1)) &
                                      channel%seg(num_seg)%DM_ETA_TN(I,J,K) = 1
         IF (K.LE.channel%seg(num_seg)%HX(I,J-1).OR.K.LE.channel%seg(num_seg)%HX(I+1,J-1)) &
                                      channel%seg(num_seg)%DM_ETA_TS(I,J,K) = 1
        ENDDO
       ENDDO
      ENDDO

      DO K = 2, NK2
       NTOPOQ(K,num_seg) = 0
       NTOPOU(K,num_seg) = 0
       NTOPOV(K,num_seg) = 0
       DO J = 1, MJ1
        DO I = 1, MI1
         IF (K.GE.channel%seg(num_seg)%KLOWQ_IJ(I,J)) THEN
          NTOPOQ(K,num_seg) = NTOPOQ(K,num_seg)+1
         ENDIF
         IF (K.GE.channel%seg(num_seg)%KLOWU_IJ(I,J)) THEN
          NTOPOU(K,num_seg) = NTOPOU(K,num_seg)+1
         ENDIF
         IF (K.GE.channel%seg(num_seg)%KLOWV_IJ(I,J)) THEN
          NTOPOV(K,num_seg) = NTOPOV(K,num_seg)+1
         ENDIF
        ENDDO
       ENDDO
      ENDDO  
      
!=============================== 
      ENDDO   ! num_seg 
!===============================   
      
      DO K = 2, NK2
       channel%seg(1)%NTOPOQ_GLOB(K) = NTOPOQ(K,1)+NTOPOQ(K,2)+NTOPOQ(K,3)+NTOPOQ(K,4)
       channel%seg(1)%NTOPOU_GLOB(K) = NTOPOU(K,1)+NTOPOU(K,2)+NTOPOU(K,3)+NTOPOU(K,4)
       channel%seg(1)%NTOPOV_GLOB(K) = NTOPOV(K,1)+NTOPOV(K,2)+NTOPOV(K,3)+NTOPOV(K,4)
      ENDDO
             
      DO num_seg = 2, 4
       DO K = 2, NK2
        channel%seg(num_seg)%NTOPOQ_GLOB(K) = channel%seg(1)%NTOPOQ_GLOB(K)
        channel%seg(num_seg)%NTOPOU_GLOB(K) = channel%seg(1)%NTOPOU_GLOB(K)
        channel%seg(num_seg)%NTOPOV_GLOB(K) = channel%seg(1)%NTOPOV_GLOB(K)
       ENDDO
      ENDDO 

!*************************************************************************
      ENDIF   ! NO_TOPO
!*************************************************************************

      END SUBROUTINE topo_prepare 
      
   END SUBROUTINE init_channel

END MODULE q3d_init_module
