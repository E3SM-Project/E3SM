MODULE turb_3d_module
! Calculate turbulence effect
! Add LocalP option: use local density for physics

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t
USE parmsld,        only: ntracer,nk1,nk2,nk3

USE constld, only: nsflux,dt,dx,dy,dz,grav,vk,zt,fnz,fnt,rhoz,rho, &
                   pibar,pibarz,pbarz,zz, &
                   nosfx,dxsq,dysq,dzsq,dxdy,physics,localp 
                   
USE utils, only: rll2cov,rll2con,con2rll                             

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
!
! JUNG: Pay attention to the way of shear calculation in Richardson number.
!
! SUBROUTINE SAM_SURFACE   : Calculate surface momentum fluxes 
!        SUBROUTINE OCEFLX
!        SUBROUTINE LANDFLX
!        FUNCTION sam_qsatw
!        FUNCTION sam_esatw
!==========================================================================    


!=======================================================================
   SUBROUTINE TURB_3D (ITT, channel)
!=======================================================================
!  Calculates eddy viscosity and diffusivity coefficients

      INTEGER, INTENT(IN) :: ITT   ! time step count (used for calling surface flux)

      type(channel_t), intent(inout) :: channel   ! channel data
      
!     Local variables
      REAL (KIND=r8) :: CRITMN,CRITMX,DELD,RAMD0S,RINUM
      REAL (KIND=r8), DIMENSION(nk2) :: DDX,DDY,TERM
      
      INTEGER :: I,J,K,ITYPE,KLOW,num_seg
      INTEGER :: mi1,mim,mim_c,mim_ce,mip,mip_c,mj1,mjm,mjm_c,mjm_ce,mjp,mjp_c    

      DELD=(DX*DY*DZ)**(1.0/3.0)
      RAMD0S=(0.23_r8*DELD)**2
      
      CRITMN=10.0_r8
      CRITMX=0.8_r8*DELD*DELD/DT

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
!     (Calculated with covariant components of wind)

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

      ENDDO  ! K-loop
      
      DO K = 2, NK2

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
        TERM(K) = 0.25_r8*TERM(K) &
                + 2.0_r8*(channel%seg(num_seg)%DEFXX(I,J,K)**2 &
                         +channel%seg(num_seg)%DEFYY(I,J,K)**2 &
                         +channel%seg(num_seg)%DEFZZ(I,J,K)**2) 
      ENDDO
      
      DO K=1,KLOW-1
        channel%seg(num_seg)%RKM(I,J,K) = 0.0_r8
        channel%seg(num_seg)%RKH(I,J,K) = 0.0_r8
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

      IF(TERM(K).EQ.0.0_r8) THEN
         channel%seg(num_seg)%RKM(I,J,K) = 0.0_r8
         channel%seg(num_seg)%RKH(I,J,K) = 0.0_r8
      ELSE
         RINUM = DDY(K)/TERM(K)

         IF(RINUM.LT.0.0_r8) THEN
            channel%seg(num_seg)%RKM(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*SQRT(1.0_r8 - 16.0_r8*RINUM)
            channel%seg(num_seg)%RKH(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*1.4_r8*SQRT(1.0_r8 - 40.0_r8*RINUM)
         ELSE IF(RINUM.LT.0.25_r8) THEN
            channel%seg(num_seg)%RKM(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*(1.0_r8 - 4.0_r8*RINUM)**4
            channel%seg(num_seg)%RKH(I,J,K) = &
                SQRT(TERM(K))*DDX(K)*1.4_r8*(1.0_r8 - 1.2_r8*RINUM) &
                                           *(1.0_r8 - 4.0_r8*RINUM)**4
         ELSE
            channel%seg(num_seg)%RKM(I,J,K) = 0.0_r8
            channel%seg(num_seg)%RKH(I,J,K) = 0.0_r8
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
      
      IF (LocalP) THEN
       DO K=KLOW,NK2
        channel%seg(num_seg)%RKM(I,J,K) = channel%seg(num_seg)%RKM(I,J,K) &
                                         *channel%seg(num_seg)%RHO_bg(I,J,K)
        channel%seg(num_seg)%RKH(I,J,K) = channel%seg(num_seg)%RKH(I,J,K) &
                                         *channel%seg(num_seg)%RHO_bg(I,J,K)
       ENDDO            
      ELSE
       DO K=KLOW,NK2
        channel%seg(num_seg)%RKM(I,J,K) = channel%seg(num_seg)%RKM(I,J,K)*RHO(K)
        channel%seg(num_seg)%RKH(I,J,K) = channel%seg(num_seg)%RKH(I,J,K)*RHO(K)
       ENDDO      
      ENDIF

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!==================================================
! Calculate the surface fluxes
!===================================================

      IF(NOSFX) THEN
        DO J=1,mj1
         DO I=1,mi1
          channel%seg(num_seg)%WTH(I,J) = 0.0_r8
          channel%seg(num_seg)%WQV(I,J) = 0.0_r8
          channel%seg(num_seg)%UW(I,J)  = 0.0_r8
          channel%seg(num_seg)%WV(I,J)  = 0.0_r8
          
          channel%seg(num_seg)%UW_CON(I,J)  = 0.0_r8
          channel%seg(num_seg)%WV_CON(I,J)  = 0.0_r8
         ENDDO
        ENDDO 
      ENDIF
          
!*******************************
      ENDDO   ! num_seg 
!******************************* 

      IF(.NOT.NOSFX) THEN 
       IF (MOD((ITT-1),nsflux) .EQ. 0) CALL SAM_SURFACE(channel)
      ENDIF 

   END SUBROUTINE turb_3d

!=======================================================================
   SUBROUTINE TURB_3D_VORT (channel)
!=======================================================================
!  Calculate turbulent tendencies of vorticity components
!  Calculate turbulent tendencies of momentum (contravariant components)        
!-----------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local variables
      INTEGER :: I,J,K,num_seg,mi1,mj1
      INTEGER :: KLOWQ_MAX,KLOWU,KLOWV,KLOWQ,KLOWQ_P

      REAL (KIND=r8) :: COEFX_K,COEFY_K,COEFA_K
      REAL (KIND=r8) :: COEFX_E,COEFY_E,COEFA_E
      REAL (KIND=r8) :: COEFX_Z,COEFY_Z,COEFA_Z 

      REAL (KIND=r8) :: KVAL_P,KVAL_M
      REAL (KIND=r8) :: FXP,FXM,FYP,FYM
      REAL (KIND=r8) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=r8) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0
      
      REAL (KIND=r8) :: TEMP1,TEMP2
      REAL (KIND=r8) :: RHO_lay(nk2),RHO_int(nk2),FLUX_LO(nk2) 
      
!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB
      
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

!=============================
      DO J = 1, mj1
      DO I = 1, mi1
!=============================      

      COEFX_K = channel%seg(num_seg)%COEFX_K(I,J)
      COEFY_K = channel%seg(num_seg)%COEFY_K(I,J)
      COEFA_K = channel%seg(num_seg)%COEFA_K(I,J)

      COEFX_E = channel%seg(num_seg)%COEFX_E(I,J)
      COEFY_E = channel%seg(num_seg)%COEFY_E(I,J)
      COEFA_E = channel%seg(num_seg)%COEFA_E(I,J)

      COEFX_Z = channel%seg(num_seg)%COEFX_Z(I,J)
      COEFY_Z = channel%seg(num_seg)%COEFY_Z(I,J)
      COEFA_Z = channel%seg(num_seg)%COEFA_Z(I,J)
      
      KLOWQ = channel%seg(num_seg)%KLOWQ_IJ(I,J)
      
      IF (LocalP) THEN
        DO K = 2,nk2
         RHO_lay(K) = channel%seg(num_seg)%RHO_bg(I,J,K)     
        ENDDO 
        DO K = 1,nk2
         RHO_int(K) = channel%seg(num_seg)%RHOint_bg(I,J,K)    
        ENDDO          
      ELSE
        DO K = 2,nk2
         RHO_lay(K) = RHO(K)     
        ENDDO 
        DO K = 1,nk2
         RHO_int(K) = RHOZ(K)     
        ENDDO    
      ENDIF   
            
!-----------------------------------------------------------
!                         Z3DX (& V)
!-----------------------------------------------------------
      KLOWV   = channel%seg(num_seg)%KLOWV_IJ(I,J)
      KLOWQ_P = channel%seg(num_seg)%KLOWQ_IJ(I,J+1)
       
      DO K = 1, KLOWV-1
       channel%seg(num_seg)%FZXTB(I,J,K) = 0.0_r8
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------

      DO K = KLOWV, KLOWQ_MAX-1
!     Influenced by topography

! For ksi -----------------------------------------------------------------
       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = 0.0_r8
        FYP_XP = 0.0_r8
        FYM_XP = 0.0_r8
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)       &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/8.0_r8 
                 
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = 0.0_r8
        FYP_XM = 0.0_r8
        FYM_XM = 0.0_r8
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)        &
                 +channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I-1,J,K+1)      &
                 +channel%seg(num_seg)%RKM(I,J,K+1)        &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K+1)    &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1))/8.0_r8 
                 
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1 .OR.   &
          channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FYP_X0 = 0.0_r8
        FYM_X0 = 0.0_r8
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = 0.0_r8
        FXP_YP = 0.0_r8
        FXM_YP = 0.0_r8
       ELSE
        KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I,J+1,K)  &
                        +channel%seg(num_seg)%RKM(I,J+1,K+1))

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = 0.0_r8
        FXP_YM = 0.0_r8
        FXM_YM = 0.0_r8
       ELSE
        KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,K) &
                        +channel%seg(num_seg)%RKM(I,J,K+1))

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1 .OR.    &
          channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = 0.0_r8
        FXM_Y0 = 0.0_r8
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

! Diabatic_effect: turbulence (in terms of momentum)
! For v -----------------------------------------------------------------

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = 0.0_r8
        FYP_XP = 0.0_r8
        FYM_XP = 0.0_r8        
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/4.0_r8 
                 
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P        
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = 0.0_r8
        FYP_XM = 0.0_r8
        FYM_XM = 0.0_r8
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)        &
                 +channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J+1,K))/4.0_r8 
                 
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1 .OR.   &
          channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FYP_X0 = 0.0_r8
        FYM_X0 = 0.0_r8
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = 0.0_r8
        FXP_YP = 0.0_r8
        FXM_YP = 0.0_r8
       ELSE
        KVAL_P = channel%seg(num_seg)%RKM(I,J+1,K)

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = 0.0_r8
        FXP_YM = 0.0_r8
        FXM_YM = 0.0_r8
       ELSE
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K)

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1 .OR.    &
          channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = 0.0_r8
        FXM_Y0 = 0.0_r8
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)
       ENDIF
                           
       FLUX_LO(K) =  TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
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
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/8.0_r8 
                 
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I-1,J,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K+1)     &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1))/8.0_r8                  

          FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
          FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
          FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P

          FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
          FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
          FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M

          FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
          FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)

        KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I,J+1,K)  &
                        +channel%seg(num_seg)%RKM(I,J+1,K+1))
        KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,K)  &
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

      ENDDO  ! K-LOOP (high level: mountain-free region)
      
      DO K = KLOWQ_MAX, NK2  
! Diabatic_effect: turbulence (in terms of momentum)
! For v -------------------------------------------------------------------          
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)           & 
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J+1,K)         &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/4.0_r8 
                 
        KVAL_M = (channel%seg(num_seg)%RKM(I-1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I-1,J+1,K)       &
                 +channel%seg(num_seg)%RKM(I,J+1,K))/4.0_r8                  

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

       FLUX_LO(K) = TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
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
          FXP = 0.5_r8*FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I,J+1,K+1))*FNZ(K)/DZSQ
               
          FXM = 0.0_r8
       ELSE
          FXP = 0.5_r8*FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I,J+1,K+1))*FNZ(K)/DZSQ
               
          FXM = 0.5_r8*FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)        &
               +channel%seg(num_seg)%RKM(I,J+1,K))*FNZ(K)/DZSQ
       ENDIF

       channel%seg(num_seg)%FZXTB(I,J,K) = channel%seg(num_seg)%FZXTB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%Z3DX(I,J,K-1), &
                                                  channel%seg(num_seg)%Z3DX(I,J,K),   &
                                                  channel%seg(num_seg)%Z3DX(I,J,K+1))
      ENDDO  ! k-loop

      DO K = KLOWV, NK2      
! Diabatic_effect: turbulence (in terms of momentum)                                                  
! For v --------------------------------      
       IF (K.EQ.KLOWV) THEN
          FXP = 0.25_r8*(FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)             &
                                           +channel%seg(num_seg)%RKM(I,J+1,K+1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))    &                   
                       *FNT(K)/DZSQ
               
          FXM = 0.0_r8
       ELSE IF (K.EQ.NK2) THEN 
          FXP = 0.0_r8
               
          FXM = 0.25_r8*(FNT(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)             &
                                           +channel%seg(num_seg)%RKM(I,J+1,K-1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))    &                   
                       *FNT(K)/DZSQ       
       ELSE
          FXP = 0.25_r8*(FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)             &
                                           +channel%seg(num_seg)%RKM(I,J+1,K+1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))    &                   
                       *FNT(K)/DZSQ
               
          FXM = 0.25_r8*(FNT(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)             &
                                           +channel%seg(num_seg)%RKM(I,J+1,K-1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I,J+1,K)))    &                   
                       *FNT(K)/DZSQ
       ENDIF

       FLUX_LO(K) = FLUX_LO(K) + TURB_V(FXP,FXM,                            &
                                        channel%seg(num_seg)%U3DY(I,J,K-1), &
                                        channel%seg(num_seg)%U3DY(I,J,K),   &
                                        channel%seg(num_seg)%U3DY(I,J,K+1))                                                                                                    
      ENDDO  ! k-loop

      ! WV (covariant components) vs. WV_CON (contravariant components)   
      IF (KLOWQ_P .EQ. KLOWQ) THEN
        TEMP1 = channel%seg(num_seg)%WV(I,J)+channel%seg(num_seg)%WV(I,J+1)
        TEMP2 = channel%seg(num_seg)%WV_CON(I,J)+channel%seg(num_seg)%WV_CON(I,J+1)      
      ELSE IF (KLOWQ_P .GT. KLOWQ) THEN
        TEMP1 = channel%seg(num_seg)%WV(I,J+1)
        TEMP2 = channel%seg(num_seg)%WV_CON(I,J+1)
      ELSE
        TEMP1 = channel%seg(num_seg)%WV(I,J)
        TEMP2 = channel%seg(num_seg)%WV_CON(I,J)
      ENDIF
      
!     Add the effect of surface fluxes to ksi-equation
      channel%seg(num_seg)%FZXTB(I,J,KLOWV)= channel%seg(num_seg)%FZXTB(I,J,KLOWV) &
                               + RHO_int(KLOWV)*FNZ(KLOWV)*FNT(KLOWV)*0.5_r8*TEMP1 &
                               /(DZSQ*channel%seg(num_seg)%RG_V(I,J)*RHO_lay(KLOWV)) 

!     Add the effect of surface fluxes to v-moment 
      FLUX_LO(KLOWV) = FLUX_LO(KLOWV) + FNT(KLOWV)*0.5_r8*TEMP2/DZ                     

      DO K = KLOWV, NK1
        channel%seg(num_seg)%FZXTB(I,J,K)  = channel%seg(num_seg)%FZXTB(I,J,K) &
                                           / RHO_int(K)                         
      ENDDO  

      DO K = KLOWV, NK2
        channel%seg(num_seg)%FV_DIA(I,J,K) = channel%seg(num_seg)%FV_DIA(I,J,K) &
                                           + FLUX_LO(K)/RHO_lay(K)                                 
      ENDDO 
      
!-----------------------------------------------------------
!                            Z3DY (& U)
!-----------------------------------------------------------
      KLOWU   = channel%seg(num_seg)%KLOWU_IJ(I,J)
      KLOWQ_P = channel%seg(num_seg)%KLOWQ_IJ(I+1,J)

      DO K = 1, KLOWU-1
       channel%seg(num_seg)%FZYTB(I,J,K) = 0.0_r8
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO K = KLOWU, KLOWQ_MAX-1
!     Influenced by topography

! For eta ------------------------------------------------------------------- 
       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = 0.0_r8
        FYP_XP = 0.0_r8
        FYM_XP = 0.0_r8
       ELSE
        KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I+1,J,K) &
                        +channel%seg(num_seg)%RKM(I+1,J,K+1))
                     
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = 0.0_r8
        FYP_XM = 0.0_r8
        FYM_XM = 0.0_r8
       ELSE
        KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,K) &
                        +channel%seg(num_seg)%RKM(I,J,K+1))
                     
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FYP_X0 = 0.0_r8
        FYM_X0 = 0.0_r8
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = 0.0_r8
        FXP_YP = 0.0_r8
        FXM_YP = 0.0_r8
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K)      &
                 +channel%seg(num_seg)%RKM(I,J,K+1)        &
                 +channel%seg(num_seg)%RKM(I,J+1,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/8.0_r8     
                     
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = 0.0_r8
        FXP_YM = 0.0_r8
        FXM_YM = 0.0_r8
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)         & 
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)       &
                 +channel%seg(num_seg)%RKM(I+1,J,K)         &
                 +channel%seg(num_seg)%RKM(I,J-1,K+1)       &
                 +channel%seg(num_seg)%RKM(I,J,K+1)         &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K+1)     &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1))/8.0_r8      
                 
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = 0.0_r8
        FXM_Y0 = 0.0_r8
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

! Diabatic_effect: turbulence (in terms of momentum)
! For u -------------------------------------------------------------------

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = 0.0_r8
        FYP_XP = 0.0_r8
        FYM_XP = 0.0_r8
       ELSE
        KVAL_P = channel%seg(num_seg)%RKM(I+1,J,K)
                     
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = 0.0_r8
        FYP_XM = 0.0_r8
        FYM_XM = 0.0_r8
       ELSE
        KVAL_M = channel%seg(num_seg)%RKM(I,J,K) 
                     
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FYP_X0 = 0.0_r8
        FYM_X0 = 0.0_r8
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = 0.0_r8
        FXP_YP = 0.0_r8
        FXM_YP = 0.0_r8
       ELSE
        KVAL_P = (channel%seg(num_seg)%RKM(I,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J+1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/4.0_r8     
                     
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = 0.0_r8
        FXP_YM = 0.0_r8
        FXM_YM = 0.0_r8
       ELSE
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)         & 
                 +channel%seg(num_seg)%RKM(I,J,K)           &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)       &
                 +channel%seg(num_seg)%RKM(I+1,J,K))/4.0_r8      
                 
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = 0.0_r8
        FXM_Y0 = 0.0_r8
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       FLUX_LO(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,     &
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
!     mountain-free region

! For eta -------------------------------------------------------------------

        KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I+1,J,K) &
                        +channel%seg(num_seg)%RKM(I+1,J,K+1))
                     
        KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,K) &
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
                 +channel%seg(num_seg)%RKM(I+1,J+1,K+1))/8.0_r8   
                       
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)          &
                 +channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K)          &
                 +channel%seg(num_seg)%RKM(I,J-1,K+1)        &
                 +channel%seg(num_seg)%RKM(I,J,K+1)          &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K+1)      &
                 +channel%seg(num_seg)%RKM(I+1,J,K+1))/8.0_r8  

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
      
      ENDDO  ! K-LOOP
      
      DO K = KLOWQ_MAX, NK2
! Diabatic_effect: turbulence (in terms of momentum)
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
                 +channel%seg(num_seg)%RKM(I+1,J+1,K))/4.0_r8   
                       
        KVAL_M = (channel%seg(num_seg)%RKM(I,J-1,K)          &
                 +channel%seg(num_seg)%RKM(I,J,K)            &
                 +channel%seg(num_seg)%RKM(I+1,J-1,K)        &
                 +channel%seg(num_seg)%RKM(I+1,J,K))/4.0_r8  

        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        

        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)                                   

       FLUX_LO(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,    &
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
      ENDDO  ! k-loop

!------------------------
! VERTICAL DIFFUSION
!------------------------
      DO K = KLOWU, NK1

! For eta ------------------------------------
       IF (K.EQ.KLOWU) THEN
          FXP = 0.5_r8*FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I+1,J,K+1))*FNZ(K)/DZSQ
          FXM = 0.0_r8
       ELSE
          FXP = 0.5_r8*FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)    &
               +channel%seg(num_seg)%RKM(I+1,J,K+1))*FNZ(K)/DZSQ
               
          FXM = 0.5_r8*FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)        &
               +channel%seg(num_seg)%RKM(I+1,J,K))*FNZ(K)/DZSQ
       ENDIF

       channel%seg(num_seg)%FZYTB(I,J,K) = channel%seg(num_seg)%FZYTB(I,J,K)          &
                                         + TURB_V(FXP,FXM,                            &
                                                  channel%seg(num_seg)%Z3DY(I,J,K-1), & 
                                                  channel%seg(num_seg)%Z3DY(I,J,K),   &
                                                  channel%seg(num_seg)%Z3DY(I,J,K+1))
      ENDDO  ! k-loop
                                                  
      DO K = KLOWU, NK2
! Diabatic_effect: turbulence (in terms of momentum)
! For u ------------------------------------
       IF (K.EQ.KLOWU) THEN
          FXP = 0.25_r8*(FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)             &
                                           +channel%seg(num_seg)%RKM(I+1,J,K+1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))    &
                       *FNT(K)/DZSQ
          FXM = 0.0_r8
       ELSE IF (K.EQ.NK2) THEN
          FXP = 0.0_r8
               
          FXM = 0.25_r8*(FNT(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)             &
                                           +channel%seg(num_seg)%RKM(I+1,J,K-1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))    &
                       *FNT(K)/DZSQ       
       ELSE
          FXP = 0.25_r8*(FNT(K+1)*(channel%seg(num_seg)%RKM(I,J,K+1)             &
                                           +channel%seg(num_seg)%RKM(I+1,J,K+1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))    &
                       *FNT(K)/DZSQ
               
          FXM = 0.25_r8*(FNT(K-1)*(channel%seg(num_seg)%RKM(I,J,K-1)             &
                                           +channel%seg(num_seg)%RKM(I+1,J,K-1)) &
                          +FNT(K)*(channel%seg(num_seg)%RKM(I,J,K)               &
                                         +channel%seg(num_seg)%RKM(I+1,J,K)))    &
                       *FNT(K)/DZSQ
       ENDIF

       FLUX_LO(K) = FLUX_LO(K) + TURB_V(FXP,FXM,                            &
                                        channel%seg(num_seg)%U3DX(I,J,K-1), & 
                                        channel%seg(num_seg)%U3DX(I,J,K),   &
                                        channel%seg(num_seg)%U3DX(I,J,K+1))                                                  
      ENDDO  ! k-loop


     ! UW (covariant components)  vs.  UW_CON(contravariant components)  
      IF (KLOWQ_P .EQ. KLOWQ) THEN
        TEMP1 = channel%seg(num_seg)%UW(I,J)+channel%seg(num_seg)%UW(I+1,J)
        TEMP2 = channel%seg(num_seg)%UW_CON(I,J)+channel%seg(num_seg)%UW_CON(I+1,J)     
      ELSE IF (KLOWQ_P .GT. KLOWQ) THEN
        TEMP1 = channel%seg(num_seg)%UW(I+1,J)
        TEMP2 = channel%seg(num_seg)%UW_CON(I+1,J)
      ELSE
        TEMP1 = channel%seg(num_seg)%UW(I,J) 
        TEMP2 = channel%seg(num_seg)%UW_CON(I,J)                      
      ENDIF
           
!     Add the effect of surface fluxes to eta-equation           
      channel%seg(num_seg)%FZYTB(I,J,KLOWU)= channel%seg(num_seg)%FZYTB(I,J,KLOWU) &
                               - RHO_int(KLOWU)*FNZ(KLOWU)*FNT(KLOWU)*0.5_r8*TEMP1 &
                               /(DZSQ*channel%seg(num_seg)%RG_U(I,J)*RHO_lay(KLOWU))  

!     Add the effect of surface fluxes to u-moment
      FLUX_LO(KLOWU) = FLUX_LO(KLOWU) + FNT(KLOWU)*0.5_r8*TEMP2/DZ                                                                                 

      DO K = KLOWU, NK1
        channel%seg(num_seg)%FZYTB(I,J,K)  = channel%seg(num_seg)%FZYTB(I,J,K) &
                                           / RHO_int(K)                        
      ENDDO    
      DO K = KLOWU, NK2
        channel%seg(num_seg)%FU_DIA(I,J,K) = channel%seg(num_seg)%FU_DIA(I,J,K) &  
                                           + FLUX_LO(K)/RHO_lay(K)                                
      ENDDO                       
!-----------------------------------------------------------
!                           Z3DZ (TOP LAYER)
!-----------------------------------------------------------
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I+1,J,NK2) &
                      +channel%seg(num_seg)%RKM(I+1,J+1,NK2))
      KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,NK2) &
                      +channel%seg(num_seg)%RKM(I,J+1,NK2))

        FXP    = channel%seg(num_seg)%RGG_V(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_U(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_U(2,I+1,J)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_V(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_U(2,I-1,J)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_U(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_U(2,I,J)*(KVAL_P-KVAL_M)

      KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKM(I,J+1,NK2) &
                      +channel%seg(num_seg)%RKM(I+1,J+1,NK2))
      KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKM(I,J,NK2)   &
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
      FXP = 0.0_r8
      FXM = (channel%seg(num_seg)%RKM(I,J,NK2)      &
            +channel%seg(num_seg)%RKM(I+1,J,NK2)    & 
            +channel%seg(num_seg)%RKM(I,J+1,NK2)    &
            +channel%seg(num_seg)%RKM(I+1,J+1,NK2)  &
            +channel%seg(num_seg)%RKM(I,J,NK1)      &
            +channel%seg(num_seg)%RKM(I+1,J,NK1)    &
            +channel%seg(num_seg)%RKM(I,J+1,NK1)    &
            +channel%seg(num_seg)%RKM(I+1,J+1,NK1)) &
           *FNZ(NK1)*FNT(NK2)/(8.0_r8*DZSQ)

      channel%seg(num_seg)%FZTOPB(I,J) = channel%seg(num_seg)%FZTOPB(I,J)           &
                                       + TURB_V(FXP,FXM,                            &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK1), &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK2), &
                                                channel%seg(num_seg)%Z3DZ(I,J,NK3))
      
      channel%seg(num_seg)%FZTOPB(I,J) = channel%seg(num_seg)%FZTOPB(I,J) &
                                       / RHO_lay(NK2)                                          

!=============================
      ENDDO   ! I-loop
      ENDDO   ! J-loop
!=============================      

!******************************  
      ENDDO   ! num_seg 
!****************************** 

   END SUBROUTINE turb_3d_vort

!=======================================================================
   SUBROUTINE TURB_3D_THERM (channel)
!=======================================================================   
!  Calculate turbulent tendencies of thermodynamic variables
 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local
      LOGICAL :: CAL_QT = .FALSE.    ! Calculate the turbulent effect of QT?
      
      INTEGER :: I,J,K,nt,klow,klowq_max,mi1,mj1,num_seg

      REAL (KIND=r8) :: COEFX,COEFY,COEFA

      REAL (KIND=r8) :: KVAL_P,KVAL_M
      REAL (KIND=r8) :: FXP,FXM,FYP,FYM
      REAL (KIND=r8) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=r8) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0
      
      REAL (KIND=r8) :: RHO_lay(nk2)

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
      
      IF (LocalP) THEN
        DO K = 2,nk2
         RHO_lay(K) = channel%seg(num_seg)%RHO_bg(I,J,K)     
        ENDDO     
      ELSE
        DO K = 2,nk2
         RHO_lay(K) = RHO(K)     
        ENDDO 
      ENDIF   
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

      DO K = KLOW, KLOWQ_MAX-1
!     Influenced by topography

        IF(channel%seg(num_seg)%DM_Q_TE(I,J,K).EQ.1) THEN
         FXP    = 0.0_r8
         FYP_XP = 0.0_r8
         FYM_XP = 0.0_r8
        ELSE
         KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKH(I+1,J,K) &
                         +channel%seg(num_seg)%RKH(I,J,K))
         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)*KVAL_P
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)*KVAL_P
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_P
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TW(I,J,K).EQ.1) THEN
         FXM    = 0.0_r8
         FYP_XM = 0.0_r8
         FYM_XM = 0.0_r8
        ELSE
         KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKH(I,J,K) &
                         +channel%seg(num_seg)%RKH(I-1,J,K))
         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)*KVAL_M
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)*KVAL_M
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)*KVAL_M
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TE(I,J,K).EQ.1 .OR.  &
           channel%seg(num_seg)%DM_Q_TW(I,J,K).EQ.1) THEN
         FYP_X0 = 0.0_r8
         FYM_X0 = 0.0_r8
        ELSE
         FYP_X0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)
         FYM_X0 = channel%seg(num_seg)%RGG_V(2,I,J-1)*(KVAL_P-KVAL_M)
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TN(I,J,K).EQ.1) THEN
         FYP    = 0.0_r8
         FXP_YP = 0.0_r8
         FXM_YP = 0.0_r8
        ELSE
         KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKH(I,J+1,K) &
                         +channel%seg(num_seg)%RKH(I,J,K))
         FYP    = channel%seg(num_seg)%RGG_V(3,I,J)*KVAL_P
         FXP_YP = channel%seg(num_seg)%RGG_U(2,I,J+1)*KVAL_P
         FXM_YP = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_P
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TS(I,J,K).EQ.1) THEN
         FYM    = 0.0_r8
         FXP_YM = 0.0_r8
         FXM_YM = 0.0_r8
        ELSE
         KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKH(I,J,K) &
                         +channel%seg(num_seg)%RKH(I,J-1,K))
         FYM    = channel%seg(num_seg)%RGG_V(3,I,J-1)*KVAL_M
         FXP_YM = channel%seg(num_seg)%RGG_U(2,I,J-1)*KVAL_M
         FXM_YM = channel%seg(num_seg)%RGG_U(2,I-1,J-1)*KVAL_M
        ENDIF

        IF(channel%seg(num_seg)%DM_Q_TN(I,J,K).EQ.1 .OR.  &
           channel%seg(num_seg)%DM_Q_TS(I,J,K).EQ.1) THEN
         FXP_Y0 = 0.0_r8
         FXM_Y0 = 0.0_r8
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

      IF (CAL_QT) THEN
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
      ENDIF   ! CAL_QT

      ENDDO   ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK2
!     Not influenced by topography (upper atmosphere)

         KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKH(I+1,J,K) &
                         +channel%seg(num_seg)%RKH(I,J,K))
         KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKH(I,J,K)   &
                         +channel%seg(num_seg)%RKH(I-1,J,K))

         FXP    = channel%seg(num_seg)%RGG_U(1,I,J)*KVAL_P
         FYP_XP = channel%seg(num_seg)%RGG_V(2,I+1,J)*KVAL_P
         FYM_XP = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_P

         FXM    = channel%seg(num_seg)%RGG_U(1,I-1,J)*KVAL_M
         FYP_XM = channel%seg(num_seg)%RGG_V(2,I-1,J)*KVAL_M
         FYM_XM = channel%seg(num_seg)%RGG_V(2,I-1,J-1)*KVAL_M

         FYP_X0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)
         FYM_X0 = channel%seg(num_seg)%RGG_V(2,I,J-1)*(KVAL_P-KVAL_M)

         KVAL_P = 0.5_r8*(channel%seg(num_seg)%RKH(I,J+1,K) &
                         +channel%seg(num_seg)%RKH(I,J,K))
         KVAL_M = 0.5_r8*(channel%seg(num_seg)%RKH(I,J,K)   &
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

      IF (CAL_QT) THEN
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
      ENDIF  ! CAL_QT

      ENDDO  ! K-LOOP (high level: mountain-free region)                         

!-------------------------------------
! VERTICAL DIFFUSION
!-------------------------------------
      DO K = 1, KLOW-1
       channel%seg(num_seg)%THAD3(I,J,K) = 0.0_r8
       channel%seg(num_seg)%QVAD3(I,J,K) = 0.0_r8
       channel%seg(num_seg)%QCAD3(I,J,K) = 0.0_r8
       channel%seg(num_seg)%QIAD3(I,J,K) = 0.0_r8
      ENDDO

      DO K = KLOW, NK2

       IF (K.EQ.KLOW) THEN
          FXP = 0.5_r8*FNT(K)*FNZ(K)*(channel%seg(num_seg)%RKH(I,J,K+1)   &
               +channel%seg(num_seg)%RKH(I,J,K))/DZSQ
          FXM = 0.0_r8
       ELSE IF (K.EQ.NK2) THEN
          FXP = 0.0_r8
          FXM = 0.5_r8*FNT(K)*FNZ(K-1)*(channel%seg(num_seg)%RKH(I,J,K) &
               +channel%seg(num_seg)%RKH(I,J,K-1))/DZSQ
       ELSE
          FXP = 0.5_r8*FNT(K)*FNZ(K)*(channel%seg(num_seg)%RKH(I,J,K+1)   &
               +channel%seg(num_seg)%RKH(I,J,K))/DZSQ
               
          FXM = 0.5_r8*FNT(K)*FNZ(K-1)*(channel%seg(num_seg)%RKH(I,J,K) &
               +channel%seg(num_seg)%RKH(I,J,K-1))/DZSQ
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
                                                  
       channel%seg(num_seg)%THAD3(I,J,K) = channel%seg(num_seg)%THAD3(I,J,K) &
                                         / RHO_lay(K)
       channel%seg(num_seg)%QVAD3(I,J,K) = channel%seg(num_seg)%QVAD3(I,J,K) &
                                         / RHO_lay(K)                                         
                                                                                               
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
       
       channel%seg(num_seg)%QCAD3(I,J,K) = channel%seg(num_seg)%QCAD3(I,J,K) &
                                         / RHO_lay(K)
       channel%seg(num_seg)%QIAD3(I,J,K) = channel%seg(num_seg)%QIAD3(I,J,K) &
                                         / RHO_lay(K)                                          
!------------------------
      ENDIF  ! PHYSICS
!------------------------

       IF (CAL_QT) THEN
       DO nt = 1, ntracer
       channel%seg(num_seg)%QTAD3(I,J,K,nt) = channel%seg(num_seg)%QTAD3(I,J,K,nt)    &
                                            + TURB_V(FXP,FXM,                         &
                                              channel%seg(num_seg)%QT3D(I,J,K-1,nt),  &
                                              channel%seg(num_seg)%QT3D(I,J,K,nt),    &
                                              channel%seg(num_seg)%QT3D(I,J,K+1,nt))
       channel%seg(num_seg)%QTAD3(I,J,K,nt) = channel%seg(num_seg)%QTAD3(I,J,K,nt) &
                                            / RHO_lay(K)                                      
       ENDDO
       ENDIF  ! CAL_QT
       
      ENDDO  ! K-LOOP

!     Add the surface fluxes

      channel%seg(num_seg)%THAD3(I,J,KLOW) = channel%seg(num_seg)%THAD3(I,J,KLOW)  &
                + FNT(KLOW)*channel%seg(num_seg)%WTH(I,J)/(RHO_lay(KLOW)*DZ)  
      channel%seg(num_seg)%QVAD3(I,J,KLOW) = channel%seg(num_seg)%QVAD3(I,J,KLOW)  &
                + FNT(KLOW)*channel%seg(num_seg)%WQV(I,J)/(RHO_lay(KLOW)*DZ)     

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
      IMPLICIT NONE

      REAL (KIND=r8), INTENT(IN) :: CX,CY,CA,KXP,KXM,KYP,KYM
      REAL (KIND=r8), INTENT(IN) :: KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0
      REAL (KIND=r8), INTENT(IN) :: KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0
      REAL (KIND=r8), INTENT(IN) :: A1,A2,A3,A4,A5,A6,A7,A8,A9
      REAL (KIND=r8) :: TURB_H_result

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
      IMPLICIT NONE

      REAL (KIND=r8), INTENT(IN) :: COP,COM,A1,A2,A3
      REAL (KIND=r8) :: TURB_V_result

      TURB_V_result = COP*(A3-A2)-COM*(A2-A1)

      END FUNCTION turb_v

!========================================================================== 
   SUBROUTINE SAM_SURFACE (channel)
!==========================================================================    
!  1. Calculate surface momentum fluxes in RLL coordinates.   
!  2. Convert the output to covariant components (UW and WV).
!  3. Convert the output to contravariant components (UW_CON and WV_CON)
!
!  Add LocalP option:
!-------------------------------------------------------------------------- 
!  NOTE: Surface flux calculation from SAM.
!  (LANDFLX) A limit on fluxes to avoid too large fluxes over large surface roughness 
!            CRM time step "dt" is used for this.
!  (OCEFLX)  zbot should be > 10m; "grav" is used. 
!
!  Temporary use of channel%seg(num_seg)%TERM1 
!  SFX(I,J): sensible heat flux <WT> (Km/s)
!  LFX(I,J): latent heat flux   <WQ> (m/s)
!  MFX(I,J): momentum flux at t-point <WV> (m/s)
!  U_H_rll(I,J): U_H in RLL coordinates (t-point)
!  V_H_rll(I,J): V_H in RLL coordinates (t-point)
!-------------------------------------------------------------------------- 
   type(channel_t), intent(inout) :: channel   ! channel data
   
   ! local variables
   REAL(KIND=r8), PARAMETER :: SALT_FACTOR = 0.981_r8
   
   REAL(KIND=r8) :: HGT   ! height h, [m]
   REAL(KIND=r8) :: T_H   ! pot. temperature at height h, [K]
   REAL(KIND=r8) :: Q_H   ! water vapor mixing ratio at height h, [g/g]
   REAL(KIND=r8) :: U_H   ! zonal wind at height h, [m/s]
   REAL(KIND=r8) :: V_H   ! merid wind at height h, [m/s]
   
   REAL(KIND=r8) :: TA_H  ! temperature at height h, [K]
   
   REAL(KIND=r8) :: T_S   ! pot. Temperature at surface, [K]
   REAL(KIND=r8) :: Q_S   ! saturated water vapor mixing ratio at surface, [g/g]
   
   REAL(KIND=r8) :: TA_S  ! temperature at surface, [K]
   REAL(kind=r8) :: Pmb_S ! surface pressure in [mb]
   
   REAL(KIND=r8) :: Z0    ! surface roughness length, [m]
   REAL(KIND=r8) :: GWET  ! surface (ground) wetness
     
   INTEGER num_seg,mi1,mj1,mip_c,mjp_c,I,J,KLOW   
          
   REAL(KIND=r8) :: VMAP(2),VMAP_COV(2),VMAP_CON(2)                                                 

   REAL(KIND=r8), DIMENSION(:,:), allocatable :: SFX,LFX,MFX,U_H_rll,V_H_rll
   REAL(KIND=r8), DIMENSION(:,:), allocatable :: PIBAR_KLOW,RHO_KLOW,PIBARZ_KLOW_1, &
                                                 PBARZ_KLOW_1,RHOZ_KLOW_1    ! JUNG_LocalP

!**************************************************************************** 
   DO num_seg = 1, 4
!**************************************************************************** 
   mi1   = channel%seg(num_seg)%mi1  ! x-size of channel segment
   mj1   = channel%seg(num_seg)%mj1  ! y-size of channel segment 
   
   mip_c = channel%seg(num_seg)%mip_c 
   mjp_c = channel%seg(num_seg)%mjp_c 

   allocate(SFX(mip_c,mjp_c))
   allocate(LFX(mip_c,mjp_c))
   allocate(MFX(mip_c,mjp_c))
   allocate(U_H_rll(mip_c,mjp_c))
   allocate(V_H_rll(mip_c,mjp_c))
   
   allocate(PIBAR_KLOW(mip_c,mjp_c))
   allocate(RHO_KLOW(mip_c,mjp_c))
   allocate(PIBARZ_KLOW_1(mip_c,mjp_c))
   allocate(PBARZ_KLOW_1(mip_c,mjp_c))
   allocate(RHOZ_KLOW_1(mip_c,mjp_c))
 
   DO J = 1, mjp_c
    DO I = 1, mip_c
     KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
     IF (LocalP) THEN
       PIBAR_KLOW(I,J)    = channel%seg(num_seg)%PImid_bg(I,J,KLOW)
       PIBARZ_KLOW_1(I,J) = channel%seg(num_seg)%PIint_bg(I,J,KLOW-1) 
       PBARZ_KLOW_1(I,J)  = channel%seg(num_seg)%Pint_bg(I,J,KLOW-1) 
       RHO_KLOW(I,J)      = channel%seg(num_seg)%RHO_bg(I,J,KLOW) 
       RHOZ_KLOW_1(I,J)   = RHO_KLOW(I,J)
     ELSE
       PIBAR_KLOW(I,J)    = PIBAR(KLOW)
       PIBARZ_KLOW_1(I,J) = PIBARZ(KLOW-1) 
       PBARZ_KLOW_1(I,J)  = PBARZ(KLOW-1)
       RHO_KLOW(I,J)      = RHO(KLOW)
       RHOZ_KLOW_1(I,J)   = RHO_KLOW(I,J)
     ENDIF
    ENDDO
   ENDDO  
         
   DO J = 1, mjp_c
    DO I = 1, mip_c
     KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
     
     HGT  = ZT(KLOW) - ZZ(KLOW-1)
     T_H  = channel%seg(num_seg)%TH3D(I,J,KLOW)
     Q_H  = channel%seg(num_seg)%QV3D(I,J,KLOW)
     
     U_H  = 0.5_r8*(channel%seg(num_seg)%U3DX(I-1,J,KLOW) &
                   +channel%seg(num_seg)%U3DX(I,J,KLOW)) 
     V_H  = 0.5_r8*(channel%seg(num_seg)%U3DY(I,J-1,KLOW) &
                   +channel%seg(num_seg)%U3DY(I,J,KLOW)) 
     
     ! Convert a contravariant wind vector to RLL wind vector. 
     VMAP = CON2RLL(channel%seg(num_seg)%AM_T(1,I,J),channel%seg(num_seg)%AM_T(2,I,J), &
                    channel%seg(num_seg)%AM_T(3,I,J),channel%seg(num_seg)%AM_T(4,I,J), &
                    U_H,V_H)      
                  
     U_H_rll(I,J) = VMAP(1)                      
     V_H_rll(I,J) = VMAP(2)  
     
     TA_H = channel%seg(num_seg)%TH3D(I,J,KLOW)*PIBAR_KLOW(I,J)
     
     T_S   = channel%seg(num_seg)%TG(I,J)/PIBARZ_KLOW_1(I,J)
     TA_S  = channel%seg(num_seg)%TG(I,J)
     Pmb_S = PBARZ_KLOW_1(I,J)*0.01_r8
     
     GWET  = channel%seg(num_seg)%GWET(I,J)

     IF(channel%seg(num_seg)%LOCEAN(I,J)) THEN

       Q_S = GWET*SALT_FACTOR*SAM_QSATW(TA_S,Pmb_S)
       CALL OCEFLX(RHO_KLOW(I,J),VMAP(1),VMAP(2),TA_H,Q_H,T_H,HGT,T_S,Q_S, &
                   SFX(I,J),LFX(I,J),MFX(I,J))
          
     ELSE
     
       Z0  = channel%seg(num_seg)%ZROUGH(I,J)
       Q_S = GWET*SAM_QSATW(TA_S,Pmb_S)
       CALL LANDFLX(T_H,T_S,Q_H,Q_S,VMAP(1),VMAP(2),HGT,Z0,    &
                    SFX(I,J),LFX(I,J),MFX(I,J))

     ENDIF  
     
    ENDDO  ! i-loop
   ENDDO   ! j-loop 

   DO J = 1, mj1
    DO I = 1, mi1
       channel%seg(num_seg)%WTH(I,J) = SFX(I,J)*RHOZ_KLOW_1(I,J)/PIBARZ_KLOW_1(I,J)
       channel%seg(num_seg)%WQV(I,J) = LFX(I,J)*RHOZ_KLOW_1(I,J)
    ENDDO
   ENDDO   

!#ifdef JUNG_TEST

   ! UW and WV are defined at t-points (Values in RLL coordinates)
   DO J = 1, mjp_c
    DO I = 1, mip_c
       channel%seg(num_seg)%UW(I,J) = U_H_rll(I,J)*MFX(I,J)*RHOZ_KLOW_1(I,J)  
       channel%seg(num_seg)%WV(I,J) = V_H_rll(I,J)*MFX(I,J)*RHOZ_KLOW_1(I,J)  
    ENDDO
   ENDDO 

   DO J = 1, mjp_c
    DO I = 1, mip_c
       ! RLL to Covariant components  
       VMAP_COV = RLL2COV(channel%seg(num_seg)%AM_T(1,I,J),channel%seg(num_seg)%AM_T(2,I,J), &
                          channel%seg(num_seg)%AM_T(3,I,J),channel%seg(num_seg)%AM_T(4,I,J), &
                          channel%seg(num_seg)%UW(I,J),channel%seg(num_seg)%WV(I,J))
       ! RLL to Contravariant components                     
       VMAP_CON = RLL2CON(channel%seg(num_seg)%AMI_T(1,I,J),channel%seg(num_seg)%AMI_T(2,I,J), &
                          channel%seg(num_seg)%AMI_T(3,I,J),channel%seg(num_seg)%AMI_T(4,I,J), &
                          channel%seg(num_seg)%UW(I,J),channel%seg(num_seg)%WV(I,J))                   
                       
       channel%seg(num_seg)%UW(I,J) = VMAP_COV(1)
       channel%seg(num_seg)%WV(I,J) = VMAP_COV(2) 
       
       channel%seg(num_seg)%UW_CON(I,J) = VMAP_CON(1)
       channel%seg(num_seg)%WV_CON(I,J) = VMAP_CON(2) 
    ENDDO
   ENDDO   

!#else
#ifdef JUNG_TEST

   DO J = 1, mjp_c
    DO I = 1, mip_c
       channel%seg(num_seg)%UW(I,J) = 0.0_r8
       channel%seg(num_seg)%WV(I,J) = 0.0_r8 
       
       channel%seg(num_seg)%UW_CON(I,J) = 0.0_r8
       channel%seg(num_seg)%WV_CON(I,J) = 0.0_r8 
    ENDDO
   ENDDO   

#endif     

   deallocate(SFX)
   deallocate(LFX)
   deallocate(MFX)
   deallocate(U_H_rll)
   deallocate(V_H_rll)
   
   deallocate(PIBAR_KLOW)
   deallocate(RHO_KLOW)
   deallocate(PIBARZ_KLOW_1)
   deallocate(PBARZ_KLOW_1)
   deallocate(RHOZ_KLOW_1)
      
!*******************************
      ENDDO   ! num_seg 
!******************************* 

   CONTAINS

   SUBROUTINE OCEFLX(rbot,ubot,vbot,tbot,qbot,thbot,zbot,ts,qs,shf,lhf,taus)

! Compute ocean to atmosphere surface fluxes of sensible, latent heat and stress components:
!
! Assume:
!   1) Neutral 10m drag coeff: 
!         cdn = .0027/U10N + .000142 + .0000764 U10N
!   2) Neutral 10m stanton number: 
!         ctn = .0327 sqrt(cdn), unstable
!         ctn = .0180 sqrt(cdn), stable
!   3) Neutral 10m dalton number:  
!         cen = .0346 sqrt(cdn)
!
! Note:
!   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   2) wind speeds should all be above a minimum speed (umin)
!--------------------------Code History---------------------------------
! Original:      Bill Large/M.Vertenstein, Sep. 1995 for CCM3.5
! Standardized:  L. Buja,     Feb 1996
! Reviewed:      B. Briegleb, March 1996
! Adopted for LES by Marat Khairoutdinov, July 1998
!------------------------------Arguments--------------------------------
! Input arguments                        
      real(KIND=r8), INTENT(IN) :: rbot    ! Density at bottom model level
      real(KIND=r8), INTENT(IN) :: ubot    ! Bottom level u wind
      real(KIND=r8), INTENT(IN) :: vbot    ! Bottom level v wind
      real(KIND=r8), INTENT(IN) :: tbot    ! Bottom level temperature
      real(KIND=r8), INTENT(IN) :: qbot    ! Bottom level specific humidity
      real(KIND=r8), INTENT(IN) :: thbot   ! Bottom level potential temperature ---> SEE NOTE BELOW.
      real(KIND=r8), INTENT(IN) :: zbot    ! Bottom level height above surface
      real(KIND=r8), INTENT(IN) :: ts      ! Surface temperature (Jung: potential temp.)
      real(KIND=r8), INTENT(IN) :: qs      ! Surface saturation specific humidity

! Output arguments                                 
      real(KIND=r8), INTENT(OUT) :: shf     ! Initial sensible heat flux (Km/s)
      real(KIND=r8), INTENT(OUT) :: lhf     ! Initial latent heat flux (m/s)
      real(KIND=r8), INTENT(OUT) :: taus    ! surface stress factor (m/s)

!---------------------------Local variables-----------------------------
      real(KIND=r8) ustar          ! ustar
      real(KIND=r8) tstar          ! tstar
      real(KIND=r8) qstar          ! qstar
      real(KIND=r8) u10n           ! neutral 10 m wind speed over ocean
      real(KIND=r8) vmag           ! Surface wind magnitude
      real(KIND=r8) thvbot         ! Bottom lev virtual potential temp
      real(KIND=r8) delt           ! potential T difference (K)
      real(KIND=r8) delq           ! specific humidity difference (kg/kg)
      real(KIND=r8) rdn            ! sqrt of neutral exchange coeff (momentum)
      real(KIND=r8) rhn            ! sqrt of neutral exchange coeff (heat)
      real(KIND=r8) ren            ! sqrt of neutral exchange coeff (tracers)          
      real(KIND=r8) rd             ! sqrt of exchange coeff (momentum)
      real(KIND=r8) rh             ! sqrt of exchange coeff (heat)
      real(KIND=r8) re             ! sqrt of exchange coeff (tracers)
      real(KIND=r8) hol            ! Ref hgt (10m) / monin-obukhov length
      real(KIND=r8) xsq            ! Temporary variable
      real(KIND=r8) xqq            ! Temporary variable
      real(KIND=r8) alz            ! ln(zbot/z10)
      
      real(KIND=r8) tau            ! Reference height stress
      real(KIND=r8) psimh          ! Stability funct at ref lev (momentum)
      real(KIND=r8) psixh          ! Stability funct at ref lev (heat & tracers) 
      real(KIND=r8) stable         ! Stability factor
      real(KIND=r8) bn             ! exchange coef funct for interpolation
      real(KIND=r8) bh             ! exchange coef funct for interpolation
      real(KIND=r8) fac            ! interpolation factor
      real(KIND=r8) ln0            ! log factor for interpolation
      real(KIND=r8) ln3            ! log factor for interpolation
      
      real(KIND=r8) gravit	     ! =9.81 m/s2
      real(KIND=r8) xkar           ! Von Karman constant (=0.4)
      
      real(KIND=r8), parameter :: zref = 10.0_r8  ! 10m reference height
      real(KIND=r8), parameter :: umin = 1.0_r8   ! Minimum wind speed at bottom level

!--------------------------Statement functions--------------------------

      real(KIND=r8) psimhu         ! Unstable part of psimh
      real(KIND=r8) psixhu         ! Unstable part of psixh
      real(KIND=r8) cdn            ! Neutral drag coeff at bottom model level
      real(KIND=r8) xd             ! Dummy argument
      real(KIND=r8) Umps           ! Wind velocity (m/sec)

      cdn(Umps)  = 0.0027_r8/Umps + 0.000142_r8 + 0.0000764_r8 * Umps
      psimhu(xd) = log((1.0_r8+xd*(2.0_r8+xd))*(1.0_r8+xd*xd)/8.0_r8) &
                 - 2.0_r8*atan(xd) + 1.571_r8
      psixhu(xd) = 2.0_r8 * log((1.0_r8 + xd*xd)/2.0_r8)

!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------  
!	     gravit = 9.81_r8
	     gravit = grav
	     xkar = vk

         if(zbot.le.zref) then
          print*,'oceflx: zbot should be > 10m'
         end if

         vmag   = max(umin, sqrt(ubot**2 + vbot**2))
         thvbot = thbot  * (1.0_r8 + 0.606_r8*qbot )
         delt   = thbot  - ts
         delq   = qbot  - qs  
         alz    = log(zbot /zref) 
        
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
! Initial guess for roots of neutral exchange coefficients, 
! assume z/L=0. and u10n is approximated by vmag.
! Stable if (thbot > ts ).

         stable = 0.5_r8 + sign(0.5_r8 , delt)
         rdn  = sqrt(cdn(vmag))
         rhn  = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8 
         ren  = 0.0346_r8 

! Initial guess of ustar, tstar and qstar

         ustar = rdn*vmag
         tstar = rhn*delt
         qstar = ren*delq

! Compute stability and evaluate all stability functions
! Stable if (thbot > ts or hol > 0 )

         hol = xkar*gravit*zbot*(tstar/thvbot+qstar/(1.0_r8/0.606_r8+qbot))/ustar**2
         hol = sign( min(abs(hol),10.0_r8), hol )
         stable = 0.5_r8 + sign(0.5_r8 , hol)
         xsq   = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
         xqq   = sqrt(xsq)
         psimh = -5.0_r8 * hol * stable + (1.0_r8-stable)*psimhu(xqq)
         psixh = -5.0_r8 * hol * stable + (1.0_r8-stable)*psixhu(xqq)

! Shift 10m neutral wind speed using old rdn coefficient

         rd   = rdn / (1.0_r8+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn

! Update the neutral transfer coefficients at 10m and neutral stability

         rdn = sqrt(cdn(u10n))
         ren = 0.0346_r8
         rhn = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8 

! Shift all coeffs to measurement height and stability

         rd = rdn / (1.0_r8+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.0_r8+rhn/xkar*(alz-psixh)) 
         re = ren / (1.0_r8+ren/xkar*(alz-psixh))

! Update ustar, tstar, qstar using updated, shifted coeffs 

         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 

!---------------------------------------------------------------
! Second iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------

! Recompute stability & evaluate all stability functions  
! Stable if (thbot > ts or hol > 0 )
 
         hol = xkar*gravit*zbot*(tstar/thvbot+qstar/(1.0_r8/0.606_r8+qbot))/ustar**2
         hol = sign( min(abs(hol),10.0_r8), hol )
         stable = 0.5_r8 + sign(0.5_r8 , hol)
         xsq   = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
         xqq   = sqrt(xsq)
         psimh = -5.0_r8 * hol * stable + (1.0_r8-stable)*psimhu(xqq)
         psixh = -5.0_r8 * hol * stable + (1.0_r8-stable)*psixhu(xqq)

! Shift 10m neutral wind speed using old rdn coefficient

         rd   = rdn / (1.0_r8+rdn/xkar*(alz-psimh))
         u10n = vmag * rd / rdn

! Update the neutral transfer coefficients at 10m and neutral stability

         rdn = sqrt(cdn(u10n))
         ren = 0.0346_r8
         rhn = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8 

! Shift all coeffs to measurement height and stability

         rd = rdn / (1.0_r8+rdn/xkar*(alz-psimh)) 
         rh = rhn / (1.0_r8+rhn/xkar*(alz-psixh)) 
         re = ren / (1.0_r8+ren/xkar*(alz-psixh))

!---------------------------------------------------------------
! Compute the fluxes
!---------------------------------------------------------------

! Update ustar, tstar, qstar using updated, shifted coeffs 

         ustar = rd * vmag 
         tstar = rh * delt 
         qstar = re * delq 

! Compute surface stress components

         tau   =  rbot  * ustar * ustar 
         
         taus  = - tau / vmag / rbot     ! this will be interpolated to wind points

! Compute heat flux components at current surface temperature
! (Define positive latent and sensible heat as upwards into the atm)

         shf  = -tau * tstar / ustar / rbot
         lhf  = -tau * qstar / ustar / rbot

END SUBROUTINE oceflx

!------------------------------------------------------------------------
SUBROUTINE LANDFLX(th, ts, qh, qs, uh, vh, h, z0, shf, lhf, taus)
! Monin-Obukhov Similarity for land
! Coded using description document  LSM4 in CESM (Marat Khairoutdinov, Oct 2016)
!------------------------------------------------------------------------
! Input:
real(KIND=r8), INTENT(IN) :: th   ! pot. temperature at height h, K
real(KIND=r8), INTENT(IN) :: ts   ! pot. Temperature at z0, K
real(KIND=r8), INTENT(IN) :: qh   ! vapor at height h, g/g
real(KIND=r8), INTENT(IN) :: qs   ! saturated vapor at z0, g/g
real(KIND=r8), INTENT(IN) :: uh   ! zonal wind at height h, m/s
real(KIND=r8), INTENT(IN) :: vh   ! merid wind at height h, m/s
real(KIND=r8), INTENT(IN) :: h    ! height h, m
real(KIND=r8), INTENT(IN) :: z0   ! surfase roughness, m

! Output:
real(KIND=r8), INTENT(OUT) :: shf   ! sensible heat flux (K m/s)
real(KIND=r8), INTENT(OUT) :: lhf   ! latent heat flux (m/s)
real(KIND=r8), INTENT(OUT) :: taus  ! zonal surface stress factor(m/s)

! Local
real(KIND=r8) r     ! bulk Richardson number
real(KIND=r8) pii, zodym, zodyh, vel
real(KIND=r8) ustar, tstar, qstar
real(KIND=r8) x, x0, xm0, xh0, xsi, xsim0, xsih0, xsi1, fm, fh, error
real(KIND=r8) psim1, psim2, psim3, psim4, xx, yy, fm0
real(KIND=r8) psih1, psih2, psih3, psih4
real(KIND=r8) z0h, zTh, zt0
real(KIND=r8), parameter :: xsim = -1.574_r8
real(KIND=r8), parameter :: xsih = -0.465_r8
real(KIND=r8), parameter :: xm = sqrt(sqrt((1.0_r8-16.0_r8*xsim)))
real(KIND=r8), parameter :: xh = sqrt(sqrt((1.0_r8-16.0_r8*xsih)))
real(KIND=r8), parameter :: errormax = 0.01_r8
integer, parameter :: nitermax = 10
integer niter

xx(yy) = sqrt(sqrt((1.0_r8-16.0_r8*yy)))

! unstable: -1.574<xsi<0
!---------------------------------
psim1(x,x0)= 2.0_r8*log((1.0_r8+x)/(1.0_r8+x0)) &
           + log((1.0_r8+x*x)/(1.0_r8+x0*x0))-2.0_r8*(atan(x)-atan(x0))
           
psih1(x,x0)=2.0_r8*log((1.0_r8+x*x)/(1.0_r8+x0*x0))
!---------------------------------

! very unstable:  xsi < -1.574
!---------------------------------
psim2(xsi,xsim0,xm0) = log(xsim/xsim0)-psim1(xm,xm0)+1.14_r8*((-xsi)**0.3333_r8-(-xsim)**0.3333_r8)

psih2(xsi,xsih0,xh0) = log(xsih/xsih0)-psih1(xh,xh0)+0.8_r8*((-xsi)**0.3333_r8-(-xsih)**0.3333_r8)
!---------------------------------

! stable: 0 < xsi < 1
!---------------------------------
psim3(xsi,xsim0) = -5.0_r8*(xsi-xsim0)

psih3(xsi,xsih0) = -5.0_r8*(xsi-xsih0)
!---------------------------------

! very stable: 0 < xsi < 1
!---------------------------------
psim4(xsi,xsim0) = log(xsi**5/xsim0)+5.0_r8*(1.-xsim0)+xsi-1.0_r8

psih4(xsi,xsih0) = log(xsi**5/xsih0)+5.0_r8*(1.-xsih0)+xsi-1.0_r8
!---------------------------------

vel = sqrt(max(0.5_r8,uh**2+vh**2))
r=max(-10.0_r8,min(0.19_r8,9.81_r8*((th-ts)/ts*(1.0_r8+0.61_r8*qh)+0.61_r8*(qh-qs))*h/vel**2))

! initial guess:

z0h=z0/h
zt0 = 0.135_r8*z0  ! roughness length for scalars
zTh=zt0/h       ! (assume ln(z0/z0t) = 2.)
zodym = log(1.0_r8/z0h)
zodyh = log(1.0_r8/zTh)

if(r.gt.0.) then
 xsi = r*zodym/(1.0_r8-5.0_r8*r)
else
 xsi = r*zodym
end if

niter = 0
error = 1000
do while (error.gt.errormax.and.niter.lt.nitermax)

xsi1 = xsi
niter = niter + 1
xsim0 = z0h*xsi
xsih0 = zTh*xsi 

if(xsi.lt.-0.01_r8) then
  if(xsi.ge.xsim) then
    x = xx(xsi)
    x0 = xx(xsim0)
    fm = zodym-psim1(x,x0)
  else
    xm0 = xx(xsim0)
    fm = psim2(xsi,xsim0,xm0)
  end if
  if(xsi.ge.xsih) then
    x = xx(xsi)
    x0 = xx(xsih0)
    fh = zodyh-psih1(x,x0)
  else
    xh0 = xx(xsih0)
    fh = psih2(xsi,xsih0,xh0)
  end if
elseif(xsi.gt.0.01_r8) then
  if(xsi.le.1.0_r8) then
    fm = zodym-psim3(xsi,xsim0)
    fh = zodyh-psih3(xsi,xsih0)
  else
    fm = psim4(xsi,xsim0)
    fh = psih4(xsi,xsih0)
  end if
else
  fm = zodym
  fh = zodyh
end if
xsi = r*fm*fm/fh
error = abs(xsi-xsi1)

end do

! limit fh and fh to avoid too large fluxes especially over large surface
! roughness. Basically, make the maximum slowdown of the velocity not bigger
! than ! 50% in one timestep.

fm0 = sqrt(0.4_r8**2*vel*dt/0.5_r8/h)
fm = max(fm0,fm)
fh = max(fh,fh/fm*fm0)

shf = 0.4_r8**2*vel/(fm*fh)*(ts-th)
lhf = 0.4_r8**2*vel/(fm*fh)*(qs-qh)
taus=-0.4_r8**2*vel/(fm*fm)

END SUBROUTINE landflx

!------------------------------------------------               
REAL(kind=r8) FUNCTION sam_qsatw(t,p)
! Saturation vapor mixing ratio. 
! Based on Flatau et.al, (JAM, 1992:1507)
!------------------------------------------------
REAL(kind=r8) :: t	! temperature (K)
REAL(kind=r8) :: p	! pressure    (mb)
REAL(kind=r8) :: esat

esat = sam_esatw(t)
sam_qsatw = 0.622_r8 * esat/max(esat,p-esat)

END

!------------------------------------------------
REAL(kind=r8) FUNCTION sam_esatw(t)
! Saturation vapor pressure.
! Based on Flatau et.al, (JAM, 1992:1507)
!------------------------------------------------
REAL(kind=r8) :: t	! temperature (K)
REAL(kind=r8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
       6.11239921_r8,    0.443987641_r8,   0.142986287D-1,  &
       0.264847430D-3,   0.302950461D-5,   0.206739458D-7,  &
       0.640689451D-10, -0.952447341D-13, -0.976195544D-15/
REAL rt

rt = max(-80.0_r8,t-273.16_r8)
sam_esatw = a0 + rt*(a1+rt*(a2+rt*(a3+rt*(a4+rt*(a5+rt*(a6+rt*(a7+a8*rt))))))) 

END

   END SUBROUTINE sam_surface      

END MODULE turb_3d_module
