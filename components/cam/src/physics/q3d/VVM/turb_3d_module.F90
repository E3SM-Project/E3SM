MODULE turb_3d_module
! Calculate turbulence effect
! (SFLUX_2D and SFLUX_3D are removed. Need to know the sfc fluxes.)

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t
USE parmsld,        only: ntracer,nk1,nk2,nk3

USE constld, only: d0_0,d1_0,d2_0,d3_0,d4_0,d8_0,d0_5,d0_25, &
                   dt,dx,dy,dz,grav,vk,zt,fnz,fnt,rhoz,rho, &
                   nosfx,dxsq,dysq,dzsq,dxdy,physics           

IMPLICIT NONE
PRIVATE

PUBLIC :: turb_3d,turb_3d_therm,turb_3d_vort

CONTAINS

! Local Subroutines:
!---------------------------------------------------------------------
! SUBROUTINE turb_3d       : Calculates eddy viscosity and diffusivity coefficients
! SUBROUTINE turb_3d_therm : Calculates turbulent tendency of thermodynamic variables
! SUBROUTINE turb_3d_vort  : Calculates turbulent tendency of vorticity
!
!     FUNCTION TURB_H      : horizontal diffusion
!     FUNCTION TURB_V      : vertical diffusion
!---------------------------------------------------------------------
!     D0_0  = 0.0_dbl_kind
!     D1_0  = 1.0_dbl_kind
!     D2_0  = 2.0_dbl_kind
!     D3_0  = 3.0_dbl_kind
!     D4_0  = 4.0_dbl_kind
!     D8_0  = 8.0_dbl_kind
!     D0_5  = 0.5_dbl_kind
!     D0_25 = 0.25_dbl_kind

!=======================================================================
   SUBROUTINE TURB_3D (channel)
!=======================================================================
!  Calculates eddy viscosity and diffusivity coefficients

      type(channel_t), intent(inout) :: channel   ! channel data
      
!     Local variables
      REAL (KIND=dbl_kind) :: CRITMN,CRITMX,DELD,POWE,RAMD0S,RINUM
      REAL (KIND=dbl_kind), DIMENSION(nk2) :: DDX,DDY,TERM
      
      INTEGER :: I,J,K,ITYPE,KLOW,num_seg
      INTEGER :: mi1,mim,mim_c,mim_ce,mip,mip_c,mj1,mjm,mjm_c,mjm_ce,mjp,mjp_c    
      
      REAL (kind=dbl_kind), PARAMETER :: RNU = 100._dbl_kind
      REAL (kind=dbl_kind), PARAMETER :: RKD = 0.21_dbl_kind
      REAL (kind=dbl_kind), PARAMETER :: RKB = 0.7_dbl_kind

      POWE=D1_0/D3_0
      DELD=(DX*DY*DZ)**POWE

      CRITMN=10._dbl_kind
      CRITMX=0.8_dbl_kind*DELD*DELD/DT

      RAMD0S=(0.23_dbl_kind*DELD)**2

!**************************************************************************** 
      DO num_seg = 1, 4
!**************************************************************************** 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim    
      mim_c  = channel%seg(num_seg)%mim_c
      mim_ce = channel%seg(num_seg)%mim_ce 
      mip    = channel%seg(num_seg)%mip    
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm    = channel%seg(num_seg)%mjm    
      mjm_c  = channel%seg(num_seg)%mjm_c 
      mjm_ce = channel%seg(num_seg)%mjm_ce   
      mjp    = channel%seg(num_seg)%mjp    
      mjp_c  = channel%seg(num_seg)%mjp_c        
      
!====================================
! Specify eddy viscosity coefficient
!====================================
!     Strain tensor components (deformation & shear terms)
!     (Calculated with covariant components of wind -- Is this right?)

      DO K = 1, NK2

       DO J = mjm_ce,mjp_c
       DO I = mim_c,mip_c
        channel%seg(num_seg)%DEFYZ(I,J,K) = &
          (channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))         &
                  /(DY*channel%seg(num_seg)%RG_V(I,J))                                &
        + (channel%seg(num_seg)%U3DY_CO(I,J,K+1)-channel%seg(num_seg)%U3DY_CO(I,J,K)) &
                  *FNZ(K)/(DZ*channel%seg(num_seg)%RG_V(I,J))   
       ENDDO
       ENDDO

       DO J = mjm_c,mjp_c
       DO I = mim_ce,mip_c
        channel%seg(num_seg)%DEFXZ(I,J,K) = &
          (channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))         & 
                  /(DX*channel%seg(num_seg)%RG_U(I,J))                                &
        + (channel%seg(num_seg)%U3DX_CO(I,J,K+1)-channel%seg(num_seg)%U3DX_CO(I,J,K)) &
                  *FNZ(K)/(DZ*channel%seg(num_seg)%RG_U(I,J)) 
       ENDDO
       ENDDO

       DO J = mjm_ce,mjp_c
       DO I = mim_ce,mip_c
        channel%seg(num_seg)%DEFXY(I,J,K) = &
          (channel%seg(num_seg)%U3DY_CO(I+1,J,K)-channel%seg(num_seg)%U3DY_CO(I,J,K)) &
                  /(DX*channel%seg(num_seg)%RG_Z(I,J))                                &
        + (channel%seg(num_seg)%U3DX_CO(I,J+1,K)-channel%seg(num_seg)%U3DX_CO(I,J,K)) &
                  /(DY*channel%seg(num_seg)%RG_Z(I,J))
       ENDDO
       ENDDO

       DO J = mjm_c,mjp_c
       DO I = mim_c,mip_c
        channel%seg(num_seg)%DEFXX(I,J,K) = &
          (channel%seg(num_seg)%U3DX_CO(I,J,K)-channel%seg(num_seg)%U3DX_CO(I-1,J,K)) &
                  /(SQRT(channel%seg(num_seg)%GG_T(3,I,J))*DX)
        channel%seg(num_seg)%DEFYY(I,J,K) = &
          (channel%seg(num_seg)%U3DY_CO(I,J,K)-channel%seg(num_seg)%U3DY_CO(I,J-1,K)) &
                  /(SQRT(channel%seg(num_seg)%GG_T(1,I,J))*DY)
        channel%seg(num_seg)%DEFZZ(I,J,K) = &
          (FNT(K)*(channel%seg(num_seg)%W3D(I,J,K)-channel%seg(num_seg)%W3D(I,J,K-1))/DZ)
       ENDDO
       ENDDO

      ENDDO  ! K-loop

      DO J=mjm_c,mjp_c
      DO I=mim_c,mip_c
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

      DO K=KLOW,NK2
        TERM(K) = channel%seg(num_seg)%DEFXY(I-1,J-1,K  )**2   &
                + channel%seg(num_seg)%DEFXY(I-1,J,K  )**2     &
                + channel%seg(num_seg)%DEFXY(I  ,J-1,K  )**2   &
                + channel%seg(num_seg)%DEFXY(I  ,J,K  )**2     &
                + channel%seg(num_seg)%DEFXZ(I-1,J  ,K  )**2   &
                + channel%seg(num_seg)%DEFXZ(I  ,J,K  )**2     &
                + channel%seg(num_seg)%DEFXZ(I-1,J  ,K-1)**2   &
                + channel%seg(num_seg)%DEFXZ(I  ,J,K-1)**2     &
                + channel%seg(num_seg)%DEFYZ(I  ,J-1,K  )**2   &
                + channel%seg(num_seg)%DEFYZ(I  ,J,K  )**2     &
                + channel%seg(num_seg)%DEFYZ(I  ,J-1,K-1)**2   &
                + channel%seg(num_seg)%DEFYZ(I  ,J,K-1)**2
      ENDDO
      DO K=KLOW,NK2
        TERM(K) = D0_25*TERM(K) &
                + D2_0*(channel%seg(num_seg)%DEFXX(I,J,K)**2 &
                       +channel%seg(num_seg)%DEFYY(I,J,K)**2 &
                       +channel%seg(num_seg)%DEFZZ(I,J,K)**2) 
      ENDDO
      
      DO K=1,KLOW-1
        channel%seg(num_seg)%RKM(I,J,K) = D0_0
        channel%seg(num_seg)%RKH(I,J,K) = D0_0
      ENDDO

      DO K=KLOW,NK2
        DDY(K) = GRAV*FNZ(K)*(channel%seg(num_seg)%TH3D(I,J,K+1)  &
                               -channel%seg(num_seg)%TH3D(I,J,K)) &
                        /(DZ*(channel%seg(num_seg)%TH3D(I,J,K+1)  &
                              +channel%seg(num_seg)%TH3D(I,J,K))) &
               + GRAV*FNZ(K-1)*(channel%seg(num_seg)%TH3D(I,J,K)  &
                             -channel%seg(num_seg)%TH3D(I,J,K-1)) &
                        /(DZ*(channel%seg(num_seg)%TH3D(I,J,K)    &
                             +channel%seg(num_seg)%TH3D(I,J,K-1)))
      ENDDO
      DO K=KLOW,NK2
        DDX(K) = RAMD0S*(VK*VK*(ZT(K)+channel%seg(num_seg)%ZROUGH(I,J))**2) &
               /(RAMD0S+(VK*VK*(ZT(K)+channel%seg(num_seg)%ZROUGH(I,J))**2))
      ENDDO

      DO K=KLOW,NK2

      IF(TERM(K).EQ.D0_0) THEN
         channel%seg(num_seg)%RKM(I,J,K) = D0_0
         channel%seg(num_seg)%RKH(I,J,K) = D0_0
      ELSE
         RINUM = DDY(K)/TERM(K)

         IF(RINUM.LT.D0_0) THEN
            channel%seg(num_seg)%RKM(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*SQRT(D1_0 - 16._dbl_kind*RINUM)
            channel%seg(num_seg)%RKH(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*1.4_dbl_kind*SQRT(D1_0 - 40._dbl_kind*RINUM)
         ELSE IF(RINUM.LT.D0_25) THEN
            channel%seg(num_seg)%RKM(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*(D1_0 - D4_0*RINUM)**4
            channel%seg(num_seg)%RKH(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*1.4_dbl_kind*(D1_0 - 1.2_dbl_kind*RINUM) &
                                                 *(D1_0 - D4_0*RINUM)**4
         ELSE
            channel%seg(num_seg)%RKM(I,J,K) = D0_0
            channel%seg(num_seg)%RKH(I,J,K) = D0_0
         ENDIF
      ENDIF

      ENDDO  ! K-loop

      DO K=KLOW,NK2
        channel%seg(num_seg)%RKM(I,J,K) = max(channel%seg(num_seg)%RKM(I,J,K),CRITMN)
        channel%seg(num_seg)%RKH(I,J,K) = max(channel%seg(num_seg)%RKH(I,J,K),CRITMN)
      ENDDO
      DO K=KLOW,NK2
        channel%seg(num_seg)%RKM(I,J,K) = min(channel%seg(num_seg)%RKM(I,J,K),CRITMX)
        channel%seg(num_seg)%RKH(I,J,K) = min(channel%seg(num_seg)%RKH(I,J,K),CRITMX)
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!==================================================
! Calculate the surface fluxes
!===================================================

      IF(NOSFX) THEN
        DO J=1,mj1
         DO I=1,mi1
          channel%seg(num_seg)%UW(I,J)  = D0_0
          channel%seg(num_seg)%WV(I,J)  = D0_0
          channel%seg(num_seg)%WTH(I,J) = D0_0
          channel%seg(num_seg)%WQV(I,J) = D0_0
         ENDDO
        ENDDO 
      ELSE
!     How to get these data ?
      ENDIF
          
!*******************************
      ENDDO   ! num_seg 
!******************************* 

   END SUBROUTINE turb_3d

!=======================================================================
   SUBROUTINE TURB_3D_VORT (channel)
!=======================================================================
!  Calculate turbulent tendencies of vorticity components
!
!  PLAN: FU_TB,FV_TB should be separately diagnosed (TURB_3D_MOMENT)          
!-----------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local variables
      INTEGER :: I,J,K,KLOWQ_MAX,KLOWU,KLOWV,num_seg,mi1,mj1

      REAL (KIND=dbl_kind) :: COEFX_K,COEFY_K,COEFA_K
      REAL (KIND=dbl_kind) :: COEFX_E,COEFY_E,COEFA_E
      REAL (KIND=dbl_kind) :: COEFX_Z,COEFY_Z,COEFA_Z 

      REAL (KIND=dbl_kind) :: KVAL_P,KVAL_M
      REAL (KIND=dbl_kind) :: FXP,FXM,FYP,FYM
      REAL (KIND=dbl_kind) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=dbl_kind) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0
      
!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB
      
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

      DO J = 1, mj1
      DO I = 1, mi1

      COEFX_K = channel%seg(num_seg)%COEFX_K(I,J)
      COEFY_K = channel%seg(num_seg)%COEFY_K(I,J)
      COEFA_K = channel%seg(num_seg)%COEFA_K(I,J)

      COEFX_E = channel%seg(num_seg)%COEFX_E(I,J)
      COEFY_E = channel%seg(num_seg)%COEFY_E(I,J)
      COEFA_E = channel%seg(num_seg)%COEFA_E(I,J)

      COEFX_Z = channel%seg(num_seg)%COEFX_Z(I,J)
      COEFY_Z = channel%seg(num_seg)%COEFY_Z(I,J)
      COEFA_Z = channel%seg(num_seg)%COEFA_Z(I,J)
      
!-----------------------------------------------------------
!                         Z3DX (& V)
!-----------------------------------------------------------
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
       
      DO K = 1, KLOWV-1
       channel%seg(num_seg)%FZXTB(I,J,K) = D0_0
       channel%seg(num_seg)%FV_TB(I,J,K) = D0_0
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------

      DO K = KLOWV, KLOWQ_MAX-1
!     Influenced by topography

! For ksi -----------------------------------------------------------------
       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)       &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/D8_0 
                 
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)        &
                 +channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I-1,J,K+1)      &
                 +channel%seg(num_seg)%RKM(I,J,K+1)        &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K+1)    &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1))/D8_0 
                 
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1 .OR.   &
          channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I,J+1,K)  &
                      +channel%seg(num_seg)%RKM(I,J+1,K+1))

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,K) &
                      +channel%seg(num_seg)%RKM(I,J,K+1))

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1 .OR.    &
          channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)
       ENDIF

       channel%seg(num_seg)%FZXTB(I,J,K) = &
                    TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J  ,K),        & 
                           channel%seg(num_seg)%Z3DX(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J+1,K))                                         

! For v -----------------------------------------------------------------

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0        
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/D4_0 
                 
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P        
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)        &
                 +channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J+1,K))/D4_0 
                 
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1 .OR.   &
          channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = channel%seg(num_seg)%RKM(I,J+1,K)

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K)

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1 .OR.    &
          channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)
       ENDIF
                           
       channel%seg(num_seg)%FV_TB(I,J,K) = &
                    TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%U3DY(I-1,J-1,K),        &
                           channel%seg(num_seg)%U3DY(I,J-1,K),          &
                           channel%seg(num_seg)%U3DY(I+1,J-1,K),        &
                           channel%seg(num_seg)%U3DY(I-1,J  ,K),        & 
                           channel%seg(num_seg)%U3DY(I,J  ,K),          &
                           channel%seg(num_seg)%U3DY(I+1,J  ,K),        &
                           channel%seg(num_seg)%U3DY(I-1,J+1,K),        &
                           channel%seg(num_seg)%U3DY(I,J+1,K),          &
                           channel%seg(num_seg)%U3DY(I+1,J+1,K))                                                           

      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

! For ksi -----------------------------------------------------------------
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           & 
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)       &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/D8_0 
                 
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I-1,J,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K+1)     &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1))/D8_0                  

          FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
          FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
          FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P

          FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
          FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
          FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M

          FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
          FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)

        KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I,J+1,K)  &
                      +channel%seg(num_seg)%RKM(I,J+1,K+1))
        KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,K)  &
                      +channel%seg(num_seg)%RKM(I,J,K+1))

          FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
          FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
          FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P

          FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
          FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
          FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

          FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
          FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)

      channel%seg(num_seg)%FZXTB(I,J,K) = &
                   TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                          FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                          FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                          channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                          channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                          channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                          channel%seg(num_seg)%Z3DX(I-1,J,K),          &
                          channel%seg(num_seg)%Z3DX(I,J,K),            &
                          channel%seg(num_seg)%Z3DX(I+1,J,K),          &
                          channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                          channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                          channel%seg(num_seg)%Z3DX(I+1,J+1,K))                               

! For v -------------------------------------------------------------------          
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           & 
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/D4_0 
                 
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K))/D4_0                  

          FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
          FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
          FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P

          FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
          FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
          FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M

          FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
          FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)          
        
        KVAL_P = channel%seg(num_seg)%RKM(I,J+1,K)
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K)

          FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
          FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
          FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P

          FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
          FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
          FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

          FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
          FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)

      channel%seg(num_seg)%FV_TB(I,J,K) = &
                   TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                          FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                          FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                          channel%seg(num_seg)%U3DY(I-1,J-1,K),        &
                          channel%seg(num_seg)%U3DY(I,J-1,K),          &
                          channel%seg(num_seg)%U3DY(I+1,J-1,K),        &
                          channel%seg(num_seg)%U3DY(I-1,J,K),          &
                          channel%seg(num_seg)%U3DY(I,J,K),            &
                          channel%seg(num_seg)%U3DY(I+1,J,K),          &
                          channel%seg(num_seg)%U3DY(I-1,J+1,K),        &
                          channel%seg(num_seg)%U3DY(I,J+1,K),          &
                          channel%seg(num_seg)%U3DY(I+1,J+1,K))                       
        
      ENDDO  ! K-LOOP (high level: mountain-free region)
      
!------------------------
! VERTICAL DIFFUSION
!------------------------
      DO K = KLOWV, NK1

! For ksi --------------------------------      
       IF (K.EQ.KLOWV) THEN
          FXP = D0_5*FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I,J+1,K+1))*FNZ(K)/(DZSQ*RHOZ(K))
               
          FXM = D0_0
       ELSE
          FXP = D0_5*FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I,J+1,K+1))*FNZ(K)/(DZSQ*RHOZ(K))
               
          FXM = D0_5*FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)          &
               +channel%seg(num_seg)%RKM(I,J+1,K))*FNZ(K)/(DZSQ*RHOZ(K))
       ENDIF

       channel%seg(num_seg)%FZXTB(I,J,K) = channel%seg(num_seg)%FZXTB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%Z3DX(I,J,K-1), &
                                                  channel%seg(num_seg)%Z3DX(I,J,K),   &
                                                  channel%seg(num_seg)%Z3DX(I,J,K+1))
                                                  
! For v --------------------------------      
       IF (K.EQ.KLOWV) THEN
          FXP = D0_25*(FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
                                         +channel%seg(num_seg)%RKM(I,J+1,K+1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))  &                   
                     *FNT(K)/(DZSQ*RHO(K))
               
          FXM = D0_0
       ELSE
          FXP = D0_25*(FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
                                         +channel%seg(num_seg)%RKM(I,J+1,K+1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))  &                   
                     *FNT(K)/(DZSQ*RHO(K))
               
          FXM = D0_25*(FNT(K-1)*RHO(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)    &
                                         +channel%seg(num_seg)%RKM(I,J+1,K-1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))  &                   
                     *FNT(K)/(DZSQ*RHO(K))
       ENDIF

       channel%seg(num_seg)%FV_TB(I,J,K) = channel%seg(num_seg)%FV_TB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%U3DY(I,J,K-1), &
                                                  channel%seg(num_seg)%U3DY(I,J,K),   &
                                                  channel%seg(num_seg)%U3DY(I,J,K+1))
                                                                                                    
      ENDDO

!     Add the surface fluxes
      channel%seg(num_seg)%FZXTB(I,J,KLOWV)= channel%seg(num_seg)%FZXTB(I,J,KLOWV) &
            + FNZ(KLOWV)*FNT(KLOWV)*channel%seg(num_seg)%WV(I,J)/(RHO(KLOWV)*DZSQ)   

      channel%seg(num_seg)%FV_TB(I,J,KLOWV)= channel%seg(num_seg)%FV_TB(I,J,KLOWV) &
            + FNT(KLOWV)*channel%seg(num_seg)%WV(I,J)/(RHO(KLOWV)*DZ)   

!-----------------------------------------------------------
!                            Z3DY (& U)
!-----------------------------------------------------------
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = 1, KLOWU-1
       channel%seg(num_seg)%FZYTB(I,J,K) = D0_0
       channel%seg(num_seg)%FU_TB(I,J,K) = D0_0
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO K = KLOWU, KLOWQ_MAX-1
!     Influenced by topography

! For eta ------------------------------------------------------------------- 
       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I+1,J,K) &
                      +channel%seg(num_seg)%RKM(I+1,J,K+1))
                     
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,K) &
                      +channel%seg(num_seg)%RKM(I,J,K+1))
                     
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J,K+1)        &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/D8_0     
                     
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)         & 
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)       &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J-1,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K+1)     &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1))/D8_0      
                 
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       channel%seg(num_seg)%FZYTB(I,J,K) = &
                    TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))         

! For u -------------------------------------------------------------------

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        KVAL_P = channel%seg(num_seg)%RKM(I+1,J,K)
                     
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K) 
                     
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/D4_0     
                     
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)         & 
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)       &
                 +channel%seg(num_seg)%RKM(I+1,J,K))/D4_0      
                 
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       channel%seg(num_seg)%FU_TB(I,J,K) = &
                    TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%U3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%U3DX(I,J-1,K),          &
                           channel%seg(num_seg)%U3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%U3DX(I-1,J  ,K),        &
                           channel%seg(num_seg)%U3DX(I,J  ,K),          &
                           channel%seg(num_seg)%U3DX(I+1,J  ,K),        &
                           channel%seg(num_seg)%U3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%U3DX(I,J+1,K),          &
                           channel%seg(num_seg)%U3DX(I+1,J+1,K))         

      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

! For eta -------------------------------------------------------------------

        KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I+1,J,K) &
                      +channel%seg(num_seg)%RKM(I+1,J,K+1))
                     
        KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,K) &
                      +channel%seg(num_seg)%RKM(I,J,K+1))

        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)

        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I,J+1,K)          &
                 +channel%seg(num_seg)%RKM(I+1,J,K)          &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I,J,K+1)          &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)        &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/D8_0   
                       
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)          &
                 +channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J-1,K+1)        &
                 +channel%seg(num_seg)%RKM(I,J,K+1)          &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1))/D8_0  

        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        

        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)                                   

       channel%seg(num_seg)%FZYTB(I,J,K) = &
                    TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,    &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I-1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I,J,K),           &
                           channel%seg(num_seg)%Z3DY(I+1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))       

! For u -------------------------------------------------------------------

        KVAL_P = channel%seg(num_seg)%RKM(I+1,J,K)
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K)

        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)

        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I,J+1,K)          &
                 +channel%seg(num_seg)%RKM(I+1,J,K)          &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/D4_0   
                       
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)          &
                 +channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K))/D4_0  

        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        

        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)                                   

       channel%seg(num_seg)%FU_TB(I,J,K) = &
                    TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,    &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                           channel%seg(num_seg)%U3DX(I-1,J-1,K),       &
                           channel%seg(num_seg)%U3DX(I,J-1,K),         &
                           channel%seg(num_seg)%U3DX(I+1,J-1,K),       &
                           channel%seg(num_seg)%U3DX(I-1,J,K),         &
                           channel%seg(num_seg)%U3DX(I,J,K),           &
                           channel%seg(num_seg)%U3DX(I+1,J,K),         &
                           channel%seg(num_seg)%U3DX(I-1,J+1,K),       &
                           channel%seg(num_seg)%U3DX(I,J+1,K),         &
                           channel%seg(num_seg)%U3DX(I+1,J+1,K))       
      ENDDO  ! K-LOOP (high level: mountain-free region)

!------------------------
! VERTICAL DIFFUSION
!------------------------
      DO K = KLOWU, NK1

! For eta ------------------------------------
       IF (K.EQ.KLOWU) THEN
          FXP = D0_5*FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I+1,J,K+1))*FNZ(K)/(DZSQ*RHOZ(K))
          FXM = D0_0
       ELSE
          FXP = D0_5*FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I+1,J,K+1))*FNZ(K)/(DZSQ*RHOZ(K))
               
          FXM = D0_5*FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)          &
               +channel%seg(num_seg)%RKM(I+1,J,K))*FNZ(K)/(DZSQ*RHOZ(K))
       ENDIF

       channel%seg(num_seg)%FZYTB(I,J,K) = channel%seg(num_seg)%FZYTB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%Z3DY(I,J,K-1), & 
                                                  channel%seg(num_seg)%Z3DY(I,J,K),   &
                                                  channel%seg(num_seg)%Z3DY(I,J,K+1))

! For u ------------------------------------
       IF (K.EQ.KLOWU) THEN
          FXP = D0_25*(FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
                                         +channel%seg(num_seg)%RKM(I+1,J,K+1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))  &
                     *FNT(K)/(DZSQ*RHO(K))
          FXM = D0_0
       ELSE
          FXP = D0_25*(FNT(K+1)*RHO(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
                                         +channel%seg(num_seg)%RKM(I+1,J,K+1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))  &
                     *FNT(K)/(DZSQ*RHO(K))
               
          FXM = D0_25*(FNT(K-1)*RHO(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)    &
                                         +channel%seg(num_seg)%RKM(I+1,J,K-1)) &
                          +FNT(K)*RHO(K)*(channel%seg(num_seg)%RKM(I,J,K)      &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))  &
                     *FNT(K)/(DZSQ*RHO(K))
       ENDIF

       channel%seg(num_seg)%FU_TB(I,J,K) = channel%seg(num_seg)%FU_TB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%U3DX(I,J,K-1), & 
                                                  channel%seg(num_seg)%U3DX(I,J,K),   &
                                                  channel%seg(num_seg)%U3DX(I,J,K+1))                                                  
      ENDDO

!     Add the surface fluxes
      channel%seg(num_seg)%FZYTB(I,J,KLOWU) = channel%seg(num_seg)%FZYTB(I,J,KLOWU) &
             -FNZ(KLOWU)*FNT(KLOWU)*channel%seg(num_seg)%UW(I,J)/(RHO(KLOWU)*DZSQ)

      channel%seg(num_seg)%FU_TB(I,J,KLOWU) = channel%seg(num_seg)%FU_TB(I,J,KLOWU) &
             +FNT(KLOWU)*channel%seg(num_seg)%UW(I,J)/(RHO(KLOWU)*DZ)
 
!-----------------------------------------------------------
!                           Z3DZ (TOP LAYER)
!-----------------------------------------------------------
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I+1,J,NK2) &
                    +channel%seg(num_seg)%RKM(I+1,J+1,NK2))
      KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,NK2) &
                    +channel%seg(num_seg)%RKM(I,J+1,NK2))

        FXP    = channel%seg(num_seg)%RGG_V(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_U(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_U(2,I+1,J)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_V(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_U(2,I-1,J)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_U(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_U(2,I,J)*(KVAL_P-KVAL_M)

      KVAL_P = D0_5*(channel%seg(num_seg)%RKM(I,J+1,NK2) &
                    +channel%seg(num_seg)%RKM(I+1,J+1,NK2))
      KVAL_M = D0_5*(channel%seg(num_seg)%RKM(I,J,NK2)   &
                    +channel%seg(num_seg)%RKM(I+1,J,NK2))

        FYP    = channel%seg(num_seg)%RGG_U(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_V(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_V(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_U(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_V(2,I,J-1)*KVAL_M

        FXP_Y0 = channel%seg(num_seg)%RGG_V(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)

      channel%seg(num_seg)%FZTOPB(I,J) = &
                  TURB_H(COEFX_Z,COEFY_Z,COEFA_Z,FXP,FXM,FYP,FYM,    &
                         FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                         FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                         channel%seg(num_seg)%Z3DZ(I-1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J-1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I-1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I,J,NK2),         &
                         channel%seg(num_seg)%Z3DZ(I+1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I-1,J+1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J+1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J+1,NK2))   

!------------------------
! VERTICAL DIFFUSION
!------------------------
      FXP = D0_0
      FXM = (channel%seg(num_seg)%RKM(I,J,NK2)      &
            +channel%seg(num_seg)%RKM(I+1,J,NK2)    & 
            +channel%seg(num_seg)%RKM(I,J+1,NK2)    &
            +channel%seg(num_seg)%RKM(I+1,J+1,NK2)  &
            +channel%seg(num_seg)%RKM(I,J,NK1)      &
            +channel%seg(num_seg)%RKM(I+1,J,NK1)    &
            +channel%seg(num_seg)%RKM(I,J+1,NK1)    &
            +channel%seg(num_seg)%RKM(I+1,J+1,NK1)) &
           *FNZ(NK1)*RHOZ(NK1)*FNT(NK2)/(D8_0*DZSQ*RHO(NK2))

      channel%seg(num_seg)%FZTOPB(I,J) = channel%seg(num_seg)%FZTOPB(I,J)           &
                                       + TURB_V(FXP,FXM,                            &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK1), &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK2), &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK3))

      ENDDO   ! I-loop
      ENDDO   ! J-loop

!******************************  
      ENDDO   ! num_seg 
!****************************** 

   END SUBROUTINE turb_3d_vort

!=======================================================================
   SUBROUTINE TURB_3D_THERM (channel)
!=======================================================================   
!  Calculate turbulent tendencies of thermodynamic variables
 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER :: I,J,K,nt,klow,klowq_max,mi1,mj1,num_seg

      REAL (KIND=dbl_kind) :: COEFX,COEFY,COEFA

      REAL (KIND=dbl_kind) :: KVAL_P,KVAL_M
      REAL (KIND=dbl_kind) :: FXP,FXM,FYP,FYM
      REAL (KIND=dbl_kind) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=dbl_kind) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0

!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB
      
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      DO J = 1, mj1
      DO I = 1, mi1
      
      COEFX = channel%seg(num_seg)%COEFX(I,J)
      COEFY = channel%seg(num_seg)%COEFY(I,J)
      COEFA = channel%seg(num_seg)%COEFA(I,J)      
      
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

      DO K = KLOW, KLOWQ_MAX-1
!     Influenced by topography

        IF(channel%seg(num_seg)%DM_Q_TE(I,J,K).EQ.1) THEN
         FXP    = D0_0
         FYP_XP = D0_0
         FYM_XP = D0_0
        ELSE
         KVAL_P = D0_5*(channel%seg(num_seg)%RKH(I+1,J,K) &
                       +channel%seg(num_seg)%RKH(I,J,K))
         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)*KVAL_P
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)*KVAL_P
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_P
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TW(I,J,K).EQ.1) THEN
         FXM    = D0_0
         FYP_XM = D0_0
         FYM_XM = D0_0
        ELSE
         KVAL_M = D0_5*(channel%seg(num_seg)%RKH(I,J,K) &
                       +channel%seg(num_seg)%RKH(I-1,J,K))
         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)*KVAL_M
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)*KVAL_M
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)*KVAL_M
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TE(I,J,K).EQ.1 .OR.  &
           channel%seg(num_seg)%DM_Q_TW(I,J,K).EQ.1) THEN
         FYP_X0 = D0_0
         FYM_X0 = D0_0
        ELSE
         FYP_X0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)
         FYM_X0 = channel%seg(num_seg)%RGG_V(2,I,J-1)*(KVAL_P-KVAL_M)
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TN(I,J,K).EQ.1) THEN
         FYP    = D0_0
         FXP_YP = D0_0
         FXM_YP = D0_0
        ELSE
         KVAL_P = D0_5*(channel%seg(num_seg)%RKH(I,J+1,K) &
                       +channel%seg(num_seg)%RKH(I,J,K))
         FYP    = channel%seg(num_seg)%RGG_V(3,I,J)*KVAL_P
         FXP_YP = channel%seg(num_seg)%RGG_U(2,I,J+1)*KVAL_P
         FXM_YP = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_P
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TS(I,J,K).EQ.1) THEN
         FYM    = D0_0
         FXP_YM = D0_0
         FXM_YM = D0_0
        ELSE
         KVAL_M = D0_5*(channel%seg(num_seg)%RKH(I,J,K) &
                       +channel%seg(num_seg)%RKH(I,J-1,K))
         FYM    = channel%seg(num_seg)%RGG_V(3,I,J-1)*KVAL_M
         FXP_YM = channel%seg(num_seg)%RGG_U(2,I,J-1)*KVAL_M
         FXM_YM = channel%seg(num_seg)%RGG_U(2,I-1,J-1)*KVAL_M
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TN(I,J,K).EQ.1 .OR.  &
           channel%seg(num_seg)%DM_Q_TS(I,J,K).EQ.1) THEN
         FXP_Y0 = D0_0
         FXM_Y0 = D0_0
        ELSE
         FXP_Y0 = channel%seg(num_seg)%RGG_U(2,I,J)*(KVAL_P-KVAL_M)
         FXM_Y0 = channel%seg(num_seg)%RGG_U(2,I-1,J)*(KVAL_P-KVAL_M)
        ENDIF

      channel%seg(num_seg)%THAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                            channel%seg(num_seg)%TH3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%TH3D(I,J-1,K),          &
                            channel%seg(num_seg)%TH3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%TH3D(I-1,J,K),          &
                            channel%seg(num_seg)%TH3D(I,J,K),            &
                            channel%seg(num_seg)%TH3D(I+1,J,K),          &
                            channel%seg(num_seg)%TH3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%TH3D(I,J+1,K),          &
                            channel%seg(num_seg)%TH3D(I+1,J+1,K))

      channel%seg(num_seg)%QVAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                            channel%seg(num_seg)%QV3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QV3D(I,J-1,K),          &
                            channel%seg(num_seg)%QV3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QV3D(I-1,J,K),          &
                            channel%seg(num_seg)%QV3D(I,J,K),            &
                            channel%seg(num_seg)%QV3D(I+1,J,K),          &
                            channel%seg(num_seg)%QV3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QV3D(I,J+1,K),          &
                            channel%seg(num_seg)%QV3D(I+1,J+1,K))

!------------------------
      IF (PHYSICS) THEN
!------------------------
      channel%seg(num_seg)%QCAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &  
                            channel%seg(num_seg)%QC3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QC3D(I,J-1,K),          &
                            channel%seg(num_seg)%QC3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QC3D(I-1,J,K),          &
                            channel%seg(num_seg)%QC3D(I,J,K),            &
                            channel%seg(num_seg)%QC3D(I+1,J,K),          &
                            channel%seg(num_seg)%QC3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QC3D(I,J+1,K),          &
                            channel%seg(num_seg)%QC3D(I+1,J+1,K))

      channel%seg(num_seg)%QIAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &   
                            channel%seg(num_seg)%QI3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QI3D(I,J-1,K),          &
                            channel%seg(num_seg)%QI3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QI3D(I-1,J,K),          &
                            channel%seg(num_seg)%QI3D(I,J,K),            &
                            channel%seg(num_seg)%QI3D(I+1,J,K),          &
                            channel%seg(num_seg)%QI3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QI3D(I,J+1,K),          &
                            channel%seg(num_seg)%QI3D(I+1,J+1,K))
!------------------------
      ENDIF  ! PHYSICS
!------------------------

      DO nt = 1, ntracer
      channel%seg(num_seg)%QTAD3(I,J,K,nt) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &                 
                            channel%seg(num_seg)%QT3D(I-1,J-1,K,nt),     &
                            channel%seg(num_seg)%QT3D(I,J-1,K,nt),       &
                            channel%seg(num_seg)%QT3D(I+1,J-1,K,nt),     &
                            channel%seg(num_seg)%QT3D(I-1,J,K,nt),       &
                            channel%seg(num_seg)%QT3D(I,J,K,nt),         &
                            channel%seg(num_seg)%QT3D(I+1,J,K,nt),       &
                            channel%seg(num_seg)%QT3D(I-1,J+1,K,nt),     &
                            channel%seg(num_seg)%QT3D(I,J+1,K,nt),       &
                            channel%seg(num_seg)%QT3D(I+1,J+1,K,nt))
      ENDDO

      ENDDO   ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK2
!     Not influenced by topography (upper atmosphere)

         KVAL_P = D0_5*(channel%seg(num_seg)%RKH(I+1,J,K) &
                       +channel%seg(num_seg)%RKH(I,J,K))
         KVAL_M = D0_5*(channel%seg(num_seg)%RKH(I,J,K)   &
                       +channel%seg(num_seg)%RKH(I-1,J,K))

         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)*KVAL_P
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)*KVAL_P
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_P

         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)*KVAL_M
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)*KVAL_M
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)*KVAL_M

         FYP_X0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)
         FYM_X0 = channel%seg(num_seg)%RGG_V(2,I,J-1)*(KVAL_P-KVAL_M)

         KVAL_P = D0_5*(channel%seg(num_seg)%RKH(I,J+1,K) &
                       +channel%seg(num_seg)%RKH(I,J,K))
         KVAL_M = D0_5*(channel%seg(num_seg)%RKH(I,J,K)   &
                       +channel%seg(num_seg)%RKH(I,J-1,K))

         FYP    = channel%seg(num_seg)%RGG_V(3,I,J)*KVAL_P
         FXP_YP = channel%seg(num_seg)%RGG_U(2,I,J+1)*KVAL_P
         FXM_YP = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_P

         FYM    = channel%seg(num_seg)%RGG_V(3,I,J-1)*KVAL_M
         FXP_YM = channel%seg(num_seg)%RGG_U(2,I,J-1)*KVAL_M
         FXM_YM = channel%seg(num_seg)%RGG_U(2,I-1,J-1)*KVAL_M

         FXP_Y0 = channel%seg(num_seg)%RGG_U(2,I,J)*(KVAL_P-KVAL_M)
         FXM_Y0 = channel%seg(num_seg)%RGG_U(2,I-1,J)*(KVAL_P-KVAL_M)

      channel%seg(num_seg)%THAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                            channel%seg(num_seg)%TH3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%TH3D(I,J-1,K),          &
                            channel%seg(num_seg)%TH3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%TH3D(I-1,J,K),          &
                            channel%seg(num_seg)%TH3D(I,J,K),            &
                            channel%seg(num_seg)%TH3D(I+1,J,K),          &
                            channel%seg(num_seg)%TH3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%TH3D(I,J+1,K),          &
                            channel%seg(num_seg)%TH3D(I+1,J+1,K))

      channel%seg(num_seg)%QVAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &  
                            channel%seg(num_seg)%QV3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QV3D(I,J-1,K),          &
                            channel%seg(num_seg)%QV3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QV3D(I-1,J,K),          &
                            channel%seg(num_seg)%QV3D(I,J,K),            &
                            channel%seg(num_seg)%QV3D(I+1,J,K),          &
                            channel%seg(num_seg)%QV3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QV3D(I,J+1,K),          &
                            channel%seg(num_seg)%QV3D(I+1,J+1,K))

!------------------------
      IF (PHYSICS) THEN
!------------------------
      channel%seg(num_seg)%QCAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                            channel%seg(num_seg)%QC3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QC3D(I,J-1,K),          &
                            channel%seg(num_seg)%QC3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QC3D(I-1,J,K),          &
                            channel%seg(num_seg)%QC3D(I,J,K),            &
                            channel%seg(num_seg)%QC3D(I+1,J,K),          &
                            channel%seg(num_seg)%QC3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QC3D(I,J+1,K),          &
                            channel%seg(num_seg)%QC3D(I+1,J+1,K))

      channel%seg(num_seg)%QIAD3(I,J,K) = &
                     TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,           &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   & 
                            channel%seg(num_seg)%QI3D(I-1,J-1,K),        &
                            channel%seg(num_seg)%QI3D(I,J-1,K),          &
                            channel%seg(num_seg)%QI3D(I+1,J-1,K),        &
                            channel%seg(num_seg)%QI3D(I-1,J,K),          &
                            channel%seg(num_seg)%QI3D(I,J,K),            &        
                            channel%seg(num_seg)%QI3D(I+1,J,K),          &
                            channel%seg(num_seg)%QI3D(I-1,J+1,K),        &
                            channel%seg(num_seg)%QI3D(I,J+1,K),          &
                            channel%seg(num_seg)%QI3D(I+1,J+1,K))
!------------------------
      ENDIF  ! PHYSICS
!------------------------

      DO nt = 1, ntracer
      channel%seg(num_seg)%QTAD3(I,J,K,nt) &
                   = TURB_H(COEFX,COEFY,COEFA,FXP,FXM,FYP,FYM,            &
                            FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,    &
                            FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,    &               
                            channel%seg(num_seg)%QT3D(I-1,J-1,K,nt),      &
                            channel%seg(num_seg)%QT3D(I,J-1,K,nt),        &
                            channel%seg(num_seg)%QT3D(I+1,J-1,K,nt),      &
                            channel%seg(num_seg)%QT3D(I-1,J,K,nt),        &
                            channel%seg(num_seg)%QT3D(I,J,K,nt),          &
                            channel%seg(num_seg)%QT3D(I+1,J,K,nt),        &
                            channel%seg(num_seg)%QT3D(I-1,J+1,K,nt),      &
                            channel%seg(num_seg)%QT3D(I,J+1,K,nt),        &
                            channel%seg(num_seg)%QT3D(I+1,J+1,K,nt))
      ENDDO

      ENDDO  ! K-LOOP (high level: mountain-free region)                         

!-------------------------------------
! VERTICAL DIFFUSION
!-------------------------------------
      DO K = 1, KLOW-1
       channel%seg(num_seg)%THAD3(I,J,K) = D0_0
       channel%seg(num_seg)%QVAD3(I,J,K) = D0_0
       channel%seg(num_seg)%QCAD3(I,J,K) = D0_0
       channel%seg(num_seg)%QIAD3(I,J,K) = D0_0
      ENDDO

      DO K = KLOW, NK2

       IF (K.EQ.KLOW) THEN
          FXP = D0_5*FNT(K)*FNZ(K)*RHOZ(K)*(channel%seg(num_seg)%RKH(I,J,K+1)   &
               +channel%seg(num_seg)%RKH(I,J,K))/(RHO(K)*DZSQ)
          FXM = D0_0
       ELSE IF (K.EQ.NK2) THEN
          FXP = D0_0
          FXM = D0_5*FNT(K)*FNZ(K-1)*RHOZ(K-1)*(channel%seg(num_seg)%RKH(I,J,K) &
               +channel%seg(num_seg)%RKH(I,J,K-1))/(RHO(K)*DZSQ)
       ELSE
          FXP = D0_5*FNT(K)*FNZ(K)*RHOZ(K)*(channel%seg(num_seg)%RKH(I,J,K+1)   &
               +channel%seg(num_seg)%RKH(I,J,K))/(RHO(K)*DZSQ)
               
          FXM = D0_5*FNT(K)*FNZ(K-1)*RHOZ(K-1)*(channel%seg(num_seg)%RKH(I,J,K) &
               +channel%seg(num_seg)%RKH(I,J,K-1))/(RHO(K)*DZSQ)
       ENDIF

       channel%seg(num_seg)%THAD3(I,J,K) = channel%seg(num_seg)%THAD3(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%TH3D(I,J,K-1), &
                                                  channel%seg(num_seg)%TH3D(I,J,K),   &
                                                  channel%seg(num_seg)%TH3D(I,J,K+1))

       channel%seg(num_seg)%QVAD3(I,J,K) = channel%seg(num_seg)%QVAD3(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%QV3D(I,J,K-1), &
                                                  channel%seg(num_seg)%QV3D(I,J,K),   &
                                                  channel%seg(num_seg)%QV3D(I,J,K+1))                        
!------------------------
      IF (PHYSICS) THEN
!------------------------
       channel%seg(num_seg)%QCAD3(I,J,K) = channel%seg(num_seg)%QCAD3(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%QC3D(I,J,K-1), &
                                                  channel%seg(num_seg)%QC3D(I,J,K),   &
                                                  channel%seg(num_seg)%QC3D(I,J,K+1))
                                                  
       channel%seg(num_seg)%QIAD3(I,J,K) = channel%seg(num_seg)%QIAD3(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%QI3D(I,J,K-1), &
                                                  channel%seg(num_seg)%QI3D(I,J,K),   &
                                                  channel%seg(num_seg)%QI3D(I,J,K+1))  
!------------------------
      ENDIF  ! PHYSICS
!------------------------

       DO nt = 1, ntracer
       channel%seg(num_seg)%QTAD3(I,J,K,nt) = channel%seg(num_seg)%QTAD3(I,J,K,nt)    &
                                            + TURB_V(FXP,FXM,                         &
                                              channel%seg(num_seg)%QT3D(I,J,K-1,nt),  &
                                              channel%seg(num_seg)%QT3D(I,J,K,nt),    &
                                              channel%seg(num_seg)%QT3D(I,J,K+1,nt))
       ENDDO

      ENDDO  ! K-LOOP

!     Add the surface fluxes

      channel%seg(num_seg)%THAD3(I,J,KLOW) = channel%seg(num_seg)%THAD3(I,J,KLOW)  &
                + FNT(KLOW)*channel%seg(num_seg)%WTH(I,J)/(RHO(KLOW)*DZ)  
      channel%seg(num_seg)%QVAD3(I,J,KLOW) = channel%seg(num_seg)%QVAD3(I,J,KLOW)  &
                + FNT(KLOW)*channel%seg(num_seg)%WQV(I,J)/(RHO(KLOW)*DZ)     

      ENDDO   ! I-Loop
      ENDDO   ! J-Loop

!******************************  
      ENDDO   ! num_seg 
!****************************** 

   END SUBROUTINE turb_3d_therm

!=======================================================================       
      FUNCTION TURB_H(CX,CY,CA,KXP,KXM,KYP,KYM, &
                      KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0, &
                      KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0, &
                      A1,A2,A3,A4,A5,A6,A7,A8,A9) &
               RESULT(TURB_H_result)
!=======================================================================              
!     Calculate the turbulent effects (horizontal)
!     CX, CY, CA: coefficients
!     KXP and others: C (see the document)
!     A: target variable
!     A7(I-1,J+1,K),A8(I,J+1,K),A9(I+1,J+1,K)
!     A4(I-1,J  ,K),A5(I,J  ,K),A6(I+1,J  ,K)
!     A1(I-1,J-1,K),A2(I,J-1,K),A3(I+1,J-1,K)
!------------------------------------------------------------------
      USE shr_kind_mod, only: dbl_kind => shr_kind_r8
      IMPLICIT NONE

      REAL (KIND=dbl_kind), INTENT(IN) :: CX,CY,CA,KXP,KXM,KYP,KYM
      REAL (KIND=dbl_kind), INTENT(IN) :: KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0
      REAL (KIND=dbl_kind), INTENT(IN) :: KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0
      REAL (KIND=dbl_kind), INTENT(IN) :: A1,A2,A3,A4,A5,A6,A7,A8,A9
      REAL (KIND=dbl_kind) :: TURB_H_result

      TURB_H_result = CX*(KXP*(A6-A5)-KXM*(A5-A4))      &
                    + CY*(KYP*(A8-A5)-KYM*(A5-A2))      &
                    + CA*(KYP_XP*(A9-A6)+KYM_XP*(A6-A3) &
                         +KYP_X0*(A8-A5)+KYM_X0*(A5-A2) &
                         -KYP_XM*(A7-A4)-KYM_XM*(A4-A1) &
                         +KXP_YP*(A9-A8)+KXM_YP*(A8-A7) &
                         +KXP_Y0*(A6-A5)+KXM_Y0*(A5-A4) &
                         -KXP_YM*(A3-A2)-KXM_YM*(A2-A1))

      END FUNCTION turb_h

!=======================================================================      
      FUNCTION TURB_V(COP,COM,A1,A2,A3) RESULT(TURB_V_result)
!=======================================================================    
!     Calculate the turbulent effects (vertical)
!     COP, COM: coefficients
!     A1(I,J,K-1),A2(I,J,K), A3(I,J,K+1) : target variables
!-----------------------------------------------------------
      USE shr_kind_mod, only: dbl_kind => shr_kind_r8
      IMPLICIT NONE

      REAL (KIND=dbl_kind), INTENT(IN) :: COP,COM,A1,A2,A3
      REAL (KIND=dbl_kind) :: TURB_V_result

      TURB_V_result = COP*(A3-A2)-COM*(A2-A1)

      END FUNCTION turb_v

END MODULE turb_3d_module
