MODULE q3d_init_module
! Contains the programs related to initializing the x-arrays
! (Not finished yet)

      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
      USE vvm_data_types, only: channel_t

      USE parmsld, only: channel_l,channel_w,nk1,nk2,nk3
      USE mphysics_variables, only: rho_int,rhol,dzl,pl0,pil0 
      
      USE constld                    

#ifdef FIXTHIS
! Need to be able to compile all of RRTMG
USE trace_gases, only: trace_gas_input
USE rrtm_grid, only: nrad_rad => nrad, &
                     nrestart_rad => nrestart, &
                     masterproc,dtn,nstat    
USE rrtm_vars, only: pres,presi
#endif
      
IMPLICIT NONE
PRIVATE

PUBLIC :: init_common,init_channel

CONTAINS

!===================================================================================
   SUBROUTINE INIT_common 
!===================================================================================
!  Initialize the common parameters and vertical profiles

      INTEGER :: K
      REAL (KIND=dbl_kind),DIMENSION(nk2) :: PBARZ,PIBARZ  ! used for radcode initialization
      
      ! Local
      REAL (kind=dbl_kind) :: PZERO = 100000._dbl_kind     ! [Pa]
      REAL (kind=dbl_kind) :: RBCP,CPBR,ALPHAW
      REAL (KIND=dbl_kind), DIMENSION(NK3) :: TV,ALPHA
      
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
      DT      = 10.0_dbl_kind           ! integration timestep of crm (s)
      DT_GCM  = 300.0_dbl_kind          ! dt of gcm (s) 
      
      DZ      = 1000._dbl_kind          ! grid size in vertical direction (m)
      DZSQ    = DZ * DZ

!      Obtained in Mapping_channel
!      DX      =                        ! crm grid size in x-direction (m) 
!      DY      =                        ! crm grid size in y-direction (m)  
!      DXSQ    =
!      DYSQ    =
!      DXDY    =  
!      DX_GCM  =                        ! gcm grid size in x-direction (m)
!      DY_GCM  =                        ! gcm grid size in y-direction (m) 
      
      DZ1     = 500._dbl_kind           ! the lowest layer depth in z (m)  
      DOMAIN  = 30000._dbl_kind         ! vertical domain for z'-selection
      ZB      = 0.0_dbl_kind            ! surface height (m)
      
      PI      = 4._8 * atan(1._8)       ! 3.14159... 
      GRAV    = 9.80616_dbl_kind        ! gravitational constant (m/s^2)
      HLF     = 2.500E+6_dbl_kind       ! latent heat of vaporization (J/kg)
      HLM     = 3.336E+5_dbl_kind       ! latent heat of fusion (J/kg) 
      CP      = 1004.5_dbl_kind         ! specific heat of dry air at constant pressure (J/kg/K)
      GAS_CONST = 287.04_dbl_kind       ! dry air gas constant
      DELTA   = 0.608_dbl_kind          ! 0.608
      VK      = 0.4_dbl_kind            ! von Karman's constant
      
      ZRSEA   = 2.0E-04                 ! surface roughness over sea (m)
      ZRLAND  = 1.0E-02                 ! surface roughness over land (m)
      
!!      SST   = 298.0_dbl_kind          ! sea surface temperature (K)

      PSFC    = 100000._dbl_kind        ! surface pressure (Pa)

      REARTH  = 6371220_dbl_kind        ! radius of the earth (m)
      OMEGA   = 7.292E-5_dbl_kind       ! earth rotation rate (s^-1)  
      RAD2DEG = 180._dbl_kind/PI        ! deg. of 1 rad (i.e., 180/pi) 

      ALADV   = 1._dbl_kind             ! alpha-value in advection

!      Obtained in Mapping_channel (after dx & dy are determined)      
!      WRXMU   =                        ! (mu/dt) in the relaxation method 
           
      DTRAD   = 3600._dbl_kind          ! time interval between radiation calls (s) 
           
      CRAD    = 1._dbl_kind/(3600._dbl_kind)  ! Rayleigh-type damping coefficient (s^-1)
      GWDHT   = 15000.0_dbl_kind              ! GWD above this height (m)
      
      DTHMAX  = 0.5_dbl_kind            ! coefficient for random perturbation (K)
      Z1PERT  = 0.0_dbl_kind            ! z_low for random perturbation  (m)
      Z2PERT  = 12500._dbl_kind         ! z_high for random perturbation (m)

!------------------------          
! Integer Parameters
!------------------------
      ITT_MAX = INT(dt_gcm/dt)          ! # of CRM time steps per one GCM time step 
      
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
        CZ2 = D0_0
        CZ1 = D1_0
      ELSE
!       Use of stretched vertical grids
        CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
        CZ1 = D1_0 - CZ2 * DOMAIN
      ENDIF

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

!***********************************************************************************
! Define main vertical profiles: rho, rhoz, pbar, pibar, fn1, fn2, thbar, qvbar
! Need more work: prescribe thbar and pibar?
!***********************************************************************************
      DO K = 1, NK3
        QVBAR(K) = D0_0
      ENDDO

      CPBR = CP / gas_const
      RBCP = 1. / CPBR

      PISFC = ( PSFC / PZERO ) ** RBCP
      PBAR(1)  = PSFC   ! [Pa]

!     Prepare THBAR & PIBAR
!!!!      CALL PRESCRIBE_IC (THBAR,PIBAR,TH3D,PREU_U,PREU_V,PREV_U,PREV_V)

      TV(1) = THBAR(1)*PIBAR(1)
      DO K = 2, NK3
       PBAR(K) = PZERO * PIBAR(K) ** CPBR
       TV(K)   = THBAR(K)*PIBAR(K)
      ENDDO

!     PROFILES OF RHO, RHOZ
      DO K = 2, NK3
        ALPHA(K) = gas_const * TV(K) / PBAR(K)
      ENDDO

      DO K = 2, NK2
       ALPHAW = (ALPHA(K) + ALPHA(K+1)) / 2.
       RHOZ(K) = 1. / ALPHAW
      ENDDO

!     SURFACE DENSITY
      RHOZ(1) = PSFC / (gas_const * TV(1))

!     DENSITY FOR k=1/2
      RHO(1) = RHOZ(1)
      DO K = 2, NK3
        RHO(K) = 1. / ALPHA(K)
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

!***********************************************************************************
!                 Define parameters and vertical profiles used in PHYSICS
!***********************************************************************************

      IF (CRAD .NE. D0_0) THEN
        K = 2
        DO WHILE (ZT(K) .LT. GWDHT)
        K = K + 1
        ENDDO
        KLOW_GWD = K                    ! index of the lowest layer where GWD is applied
      ENDIF

!--------------------------
    IF (physics) THEN
!--------------------------    

      ! Initialize Parameters for Microphysics    
      CALL initialize_microphysics (dt) 
      
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
          (DZL(K)*RHOL(K) + DZL(K+1)*RHOL(K+1))/(D2_0*(ZT(K+2) - ZT(K+1)))
      ENDDO

#ifdef FIXTHIS
! Need to be able to compile all of RRTMG  
       
      !-------------------- 
      IF (radcode) THEN
      !--------------------
      ! Initialize Parameters for Radiation 

!!!      masterproc = my_task == 0     Q: how to determine masterproc?      
      
      ! Get vertical grid parameters and define vertical pressure levels
      PIBARZ(1) = PIBAR(1)
      PBARZ(1)  = PBAR(1)

      DO K = 2, nk2
        PIBARZ(K)= PIBAR(K) + (PIBAR(K+1)-PIBAR(K)) / (ZT(K+1)-ZT(K)) * (ZZ(K)-ZT(K))
      ENDDO
      DO K = 2, nk2
        PBARZ(K) = psfc * PIBARZ(K) ** ( cp / gas_const )
      ENDDO

      PRES(:)  = PBAR(2:nk2) * 1.D-2    ! [mb]
      PRESI(:) = PBARZ(:) * 1.D-2

      DTN  = DT
      NSTAT = 100000
      NRAD_rad = NRAD
      NRESTART_rad = 0
      ! Only call tracesini if nrestart==0 (Initialize trace gas vertical profiles)
      ! However, currently, the profiles are read by trace_gases.f90.
      ! "call tracesini" is commented out. (Seems, no use of "nrestart" for radiation)

      ! Read in trace gases (fort 91) by master process and broadcast them.
      CALL TRACE_GAS_INPUT(channel_l, channel_w, nk1, PBAR(2:nk2), PBARZ)
#endif

      !--------------------
      ENDIF  ! RADCODE
      !--------------------
!--------------------------      
    ENDIF  ! PHYSICS   
!--------------------------    

   END SUBROUTINE init_common

!===================================================================================
   SUBROUTINE INIT_channel (channel)
!===================================================================================
!  Initialize the variables of channel segments

      type(channel_t), intent(inout) :: channel  ! Channel data

      INTEGER :: I,J,K
      INTEGER :: num_seg,mi1,mj1 
    
!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 


!     Calculate the coefficients used in turbulence and diffusion      
      DO J = 1, mj1
      DO I = 1, mi1  
      channel%seg(num_seg)%COEFX(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA(I,J) = 1./(channel%seg(num_seg)%RG_T(I,J)*DXDY*4.)
      
      channel%seg(num_seg)%COEFX_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_K(I,J) = 1./(channel%seg(num_seg)%RG_V(I,J)*DXDY*4.)

      channel%seg(num_seg)%COEFX_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_E(I,J) = 1./(channel%seg(num_seg)%RG_U(I,J)*DXDY*4.)

      channel%seg(num_seg)%COEFX_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXSQ)
      channel%seg(num_seg)%COEFY_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DYSQ)
      channel%seg(num_seg)%COEFA_Z(I,J) = 1./(channel%seg(num_seg)%RG_Z(I,J)*DXDY*4.)      
      ENDDO
      ENDDO

!     Calculate the coefficients used in relaxation_3d     
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

!     Calculate the coefficients used in relaxation_2d        
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

      CALL TOPO_PREPARE (channel)
      
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
      
      REAL (kind=dbl_kind) :: HEIGHT

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
        channel%seg(num_seg)%TOPOZ(I,J) = D0_0
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
       channel%seg(:)%NTOPOQ_GLOB(K) = NTOPOQ(K,1)+NTOPOQ(K,2)+NTOPOQ(K,3)+NTOPOQ(K,4)
       channel%seg(:)%NTOPOU_GLOB(K) = NTOPOU(K,1)+NTOPOU(K,2)+NTOPOU(K,3)+NTOPOU(K,4)
       channel%seg(:)%NTOPOV_GLOB(K) = NTOPOV(K,1)+NTOPOV(K,2)+NTOPOV(K,3)+NTOPOV(K,4)
      ENDDO 

!*************************************************************************
      ENDIF   ! NO_TOPO
!*************************************************************************

      END SUBROUTINE topo_prepare 
      
   END SUBROUTINE init_channel

END MODULE q3d_init_module
