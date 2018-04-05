MODULE advec_3d_module

      USE shr_kind_mod, only: r8 => shr_kind_r8
      USE parmsld,      only: nk1,nk2,nk3
      USE constld,      only: dx,dy,dz,dt,aladv,rho,rhoz,fnt,fnz,zz

IMPLICIT NONE
PRIVATE

PUBLIC :: advec_3d

CONTAINS

! Local Subroutines:
!------------------------------------------------------------------------
! SUBROUTINE advec_3d
!------------------------------------------------------------------------

!========================================================================
      SUBROUTINE ADVEC_3D (mi1,mim,mim_a,mim_c,mip,mip_c,                 &
                           mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,                 &
                           Q,TERMA,u3dx,u3dy,w3d,rg_t,rg_u,rg_v,klowq_ij, &
                           POSITIVE,TPHGT,DIV_CAL,TERM_ADD)
!========================================================================
!     3-D Advection of thermodynamic variables
!
!     INPUT : Q (quantity to be advected), wind components, and parameters
!             POSITIVE=.T. (use positive definite scheme)
!             DIV_CAL =.T. (calculate the mass convergence rate)
!             TPHGT    : prescribe tropopause height for vertical advection 
!             TERM_ADD : used for converting flx conv form to advec form
!
!     OUTPUT: TERMA (flux convergence/advective tendency)
!-------------------------------------------------------------------------  
      INTEGER, INTENT(IN) :: mi1,mim,mim_a,mim_c,mip,mip_c, &
                             mj1,mjm,mjm_a,mjm_c,mjp,mjp_c
                               
      REAL (KIND=r8), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(IN)  :: Q 
      REAL (KIND=r8), DIMENSION(mi1,mj1,nk2),         INTENT(OUT) :: TERMA 
      
      ! Wind field
      REAL (KIND=r8), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(IN)  :: U3DX,U3DY 
      REAL (KIND=r8), DIMENSION(mim:mip,mjm:mjp,NK2), INTENT(IN)  :: W3D 
      ! Mapping information
      REAL (KIND=r8), DIMENSION(mim:mip,mjm:mjp),     INTENT(IN)  :: RG_T,RG_U,RG_V 
      ! Topography information
      INTEGER, DIMENSION(mim:mip,mjm:mjp), INTENT(IN)  :: KLOWQ_IJ
      
      LOGICAL, INTENT(IN) :: POSITIVE

      REAL (KIND=r8), OPTIONAL, INTENT(IN) :: TPHGT     ! tropopause height
      LOGICAL, OPTIONAL, INTENT(IN) :: DIV_CAL
      REAL (KIND=r8), DIMENSION(mi1,mj1,nk2),OPTIONAL,INTENT(IN) :: TERM_ADD

!     Local variables ---------
      REAL (KIND=r8) :: TEMPX(mim_a:mip_c,mj1,NK2),TEMPY(mi1,mjm_a:mjp_c,NK2)
      REAL (KIND=r8) :: UPI(mim_a:mip_c,mjm_a:mjp_c,NK2) 
      REAL (KIND=r8) :: UMI(mim_a:mip_c,mjm_a:mjp_c,NK2)   
      REAL (KIND=r8) :: UPSRI(mim_a:mip_c,mjm_a:mjp_c,NK2)
      REAL (KIND=r8) :: UMSRI(mim_a:mip_c,mjm_a:mjp_c,NK2)
      REAL (KIND=r8) :: FLXI(0:mi1,0:mj1,NK2)

!     For positive definite advection
      REAL (KIND=r8) :: EPSIL = 1.0d-30
      REAL (KIND=r8) :: P_X(mim_c:mip_c,mj1,nk2)
      REAL (KIND=r8) :: P_Y(mi1,mjm_c:mjp_c,nk2)
      REAL (KIND=r8) :: P_Z(mi1,mj1,nk1)
      REAL (KIND=r8) :: A0,GAMMA_P,GAMMA_HAT_P,GAMMA_M,GAMMA_HAT_M, &
                              BETA_P,BETA_HAT_P,BETA_M,BETA_HAT_M

      LOGICAL :: USE_MU = .FALSE.                        
      REAL (KIND=r8) :: COEF_P,COEF_M,MU_P,MU_M,ALPHA_P,ALPHA_M            

      INTEGER :: I, J, K, KLOW

!****************************************
      IF (.not.PRESENT(DIV_CAL)) THEN
!****************************************

!----------------------------
      IF (POSITIVE) THEN
!----------------------------
      DO K=2,NK2
       DO J=1,MJ1
        DO I=mim_c,mip_c
         P_X(I,J,K)=(Q(I-1,J,K)-2.0_r8*Q(I,J,K)+Q(I+1,J,K))**2 + EPSIL
        ENDDO
       ENDDO
       DO J=mjm_c,mjp_c
        DO I=1,MI1
         P_Y(I,J,K)=(Q(I,J-1,K)-2.0_r8*Q(I,J,K)+Q(I,J+1,K))**2 + EPSIL
        ENDDO
       ENDDO
      ENDDO

      DO K=klow+1,nk1
       DO J=1,MJ1
        DO I=1,MI1
         P_Z(I,J,K)=(Q(I,J,K-1)-2.0_r8*Q(I,J,K)+Q(I,J,K+1))**2 + EPSIL
        ENDDO
       ENDDO
      ENDDO
!----------------------------
      ENDIF   ! POSITIVE
!----------------------------

!****************************************
      ENDIF
!****************************************

!=====================================================================
!                            Zonal advection
!=====================================================================
      DO K = 2, NK2
       DO J = 1, MJ1
        DO I = mim_a,mip_c
         TEMPX(I,J,K)=U3DX(I,J,K)*rho(K)*RG_U(I,J)
        ENDDO
       ENDDO
      ENDDO

!*********************************
      IF (PRESENT(DIV_CAL)) THEN
!*********************************

      DO K=2,NK2
       DO J=1,MJ1
        DO I=0,MI1
         FLXI(I,J,K) = 0.5_r8*TEMPX(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))
        ENDDO
       ENDDO
      ENDDO

!*********************************
      ELSE
!*********************************

      DO K=2,NK2
       DO J=1,MJ1
        DO I=mim_a,mip_c
         UPI(I,J,K)=0.5_r8*(TEMPX(I,J,K)+ABS(TEMPX(I,J,K)))
         UMI(I,J,K)=0.5_r8*(TEMPX(I,J,K)-ABS(TEMPX(I,J,K)))
        ENDDO
       ENDDO
      ENDDO

      DO K=2,NK2
       DO J=1,mj1
        DO I=mim_a,mip_c
         UPSRI(I,J,K)=SQRT(UPI(I,J,K))
         UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
        ENDDO
       ENDDO
      ENDDO

      DO K=2,NK2
       DO J=1,MJ1
        DO I=0,MI1

!       -------------------
        IF (POSITIVE) THEN
!       -------------------
         IF (USE_MU) THEN
           MU_P    = ABS(U3DX(I,J,K))*dt/dx
           ALPHA_P = (1.0_r8 + MU_P)/6.0_r8
           COEF_P  = (1.0_r8 - 2.0_r8*ALPHA_P)/(2.0_r8*ALPHA_P)

           MU_M    = ABS(U3DX(I+1,J,K))*dt/dx
           ALPHA_M = (1.0_r8 + MU_M)/6.0_r8
           COEF_M  = (1.0_r8 - 2.0_r8*ALPHA_M)/(2.0_r8*ALPHA_M)
         ELSE
           COEF_P  =  2.0_r8
           COEF_M  =  2.0_r8
         ENDIF

         a0 = MAX(Q(I,J,K),0.0_r8)*MAX(Q(I+1,J,K),0.0_r8)

         GAMMA_P     = (P_X(I,J,K)**2)/(P_X(I,J,K)**2 + a0)
         GAMMA_HAT_P = GAMMA_P
         BETA_P      = 1.0_r8 + COEF_P*GAMMA_P
         BETA_HAT_P  = 1.0_r8 - GAMMA_HAT_P

         GAMMA_M     = (P_X(I+1,J,K)**2)/(P_X(I+1,J,K)**2 + a0)
         GAMMA_HAT_M = GAMMA_M
         BETA_M      = 1.0_r8 + COEF_M*GAMMA_M
         BETA_HAT_M  = 1.0_r8 - GAMMA_HAT_M
!       -------------------
        ELSE
!       -------------------
         BETA_P      = 1.0_r8
         BETA_HAT_P  = 1.0_r8
         BETA_M      = 1.0_r8
         BETA_HAT_M  = 1.0_r8
!       -------------------
        ENDIF
!       -------------------

         FLXI(I,J,K) = &
                 0.5_r8*TEMPX(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))                      &
         - aladv*(UPI(I,J,K)*BETA_P*(Q(I+1,J,K)-Q(I,J,K))                       &
               -UPSRI(I,J,K)*UPSRI(I-1,J,K)*BETA_HAT_P*(Q(I,J,K)-Q(I-1,J,K))    &
               +UMI(I,J,K)*BETA_M*(Q(I,J,K)-Q(I+1,J,K))                         &
               +UMSRI(I,J,K)*UMSRI(I+1,J,K)*BETA_HAT_M*(Q(I+1,J,K)-Q(I+2,J,K))) &
               /6.0_r8
        ENDDO
       ENDDO
      ENDDO

!*********************************
      ENDIF  ! DIV_CAL
!*********************************

      DO K=2,NK2
       DO J=1,MJ1
        DO I=1,MI1
         terma(I,J,K)=-(FLXI(I,J,K)-FLXI(I-1,J,K))/(RG_T(I,J)*dx)
        ENDDO
       ENDDO
      ENDDO

!=====================================================================
!                          Meridional advection
!=====================================================================
      DO K = 2, NK2
       DO J = mjm_a, mjp_c
        DO I = 1, MI1
        TEMPY(I,J,K)=U3DY(I,J,K)*rho(K)*RG_V(I,J)
        ENDDO
       ENDDO
      ENDDO

!*********************************
      IF (PRESENT(DIV_CAL)) THEN
!*********************************
      DO K=2,NK2
       DO J=0,MJ1
        DO I=1,MI1
         FLXI(I,J,K) = 0.5_r8*TEMPY(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))
        ENDDO
       ENDDO
      ENDDO
!*********************************
      ELSE
!*********************************

      DO K=2,NK2
       DO J=mjm_a,mjp_c
        DO I=1,mi1
         UPI(I,J,K)=0.5_r8*(TEMPY(I,J,K)+ABS(TEMPY(I,J,K)))
         UMI(I,J,K)=0.5_r8*(TEMPY(I,J,K)-ABS(TEMPY(I,J,K)))
        ENDDO
       ENDDO
      ENDDO

      DO K=2,NK2
       DO J=mjm_a,mjp_c
        DO i=1,mi1
         UPSRi(i,J,K)=SQRT(UPi(i,J,K))
         UMSRi(i,J,K)=SQRT(ABS(UMi(i,J,K)))
        ENDDO
       ENDDO
      ENDDO

      DO K=2,NK2
       DO J=0,MJ1
        DO I=1,MI1

!       -------------------
        IF (POSITIVE) THEN
!       -------------------
         IF (USE_MU) THEN
           MU_P    = ABS(U3DY(I,J,K))*dt/dy
           ALPHA_P = (1.0_r8 + MU_P)/6.0_r8
           COEF_P  = (1.0_r8 - 2.0_r8*ALPHA_P)/(2.0_r8*ALPHA_P)

           MU_M    = ABS(U3DY(I,J+1,K))*dt/dy
           ALPHA_M = (1.0_r8 + MU_M)/6.0_r8
           COEF_M  = (1.0_r8 - 2.0_r8*ALPHA_M)/(2.0_r8*ALPHA_M)
         ELSE
           COEF_P  =  2.0_r8
           COEF_M  =  2.0_r8
         ENDIF

         a0 = MAX(Q(I,J,K),0.0_r8)*MAX(Q(I,J+1,K),0.0_r8)

         GAMMA_P     = (P_Y(I,J,K)**2)/(P_Y(I,J,K)**2 + a0)
         GAMMA_HAT_P = GAMMA_P
         BETA_P      = 1.0_r8 + COEF_P*GAMMA_P
         BETA_HAT_P  = 1.0_r8 - GAMMA_HAT_P

         GAMMA_M     = (P_Y(I,J+1,K)**2)/(P_Y(I,J+1,K)**2 + a0)
         GAMMA_HAT_M = GAMMA_M
         BETA_M      = 1.0_r8 + COEF_M*GAMMA_M
         BETA_HAT_M  = 1.0_r8 - GAMMA_HAT_M
!       -------------------
        ELSE
!       -------------------
         BETA_P      = 1.0_r8
         BETA_HAT_P  = 1.0_r8
         BETA_M      = 1.0_r8
         BETA_HAT_M  = 1.0_r8
!       -------------------
        ENDIF
!       -------------------

         FLXI(I,J,K) = &
              0.5_r8*TEMPY(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))                         &
         - aladv*(UPI(I,J,K)*BETA_P*(Q(I,J+1,K)-Q(I,J,K))                       &
               -UPSRI(I,J,K)*UPSRI(I,J-1,K)*BETA_HAT_P*(Q(I,J,K)-Q(I,J-1,K))    &
               +UMI(I,J,K)*BETA_M*(Q(I,J,K)-Q(I,J+1,K))                         &
               +UMSRI(I,J,K)*UMSRI(I,J+1,K)*BETA_HAT_M*(Q(I,J+1,K)-Q(I,J+2,K))) &
               /6.0_r8
        ENDDO
       ENDDO
      ENDDO

!*********************************
      ENDIF  ! DIV_CAL
!*********************************

      DO K=2,NK2
       DO J=1,MJ1
        DO I=1,MI1
         terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(RG_T(I,J)*dy)
        ENDDO
       ENDDO
      ENDDO

!=====================================================================
!                             Vertical advection
!=====================================================================
      DO J=1,MJ1
      DO I=1,MI1

      KLOW = KLOWQ_IJ(I,J)

      DO K=KLOW,NK1
       TEMPX(I,J,K)=W3D(I,J,K)*rhoz(K)
      ENDDO

!*********************************
      IF (PRESENT(DIV_CAL)) THEN
!*********************************

      DO K=KLOW+1,NK1-1
       FLXI(I,J,K) = 0.5_r8*TEMPX(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))
      ENDDO

!     LEVEL: K=NK1
      IF(TEMPX(I,J,NK1).GE.0.0_r8) THEN
       FLXI(I,J,NK1) = 0.5_r8*TEMPX(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
      ELSE
       FLXi(I,J,NK1)=0.5_r8*TEMPX(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
       IF (POSITIVE) THEN
        FLXI(I,J,NK1)=MAX(FLXI(I,J,NK1),-dz*rho(NK2)*Q(I,J,NK2) &
                     /(fnt(NK2)*dt))
       ENDIF
      ENDIF

!     LEVEL: K=KLOW
      IF(TEMPX(I,J,KLOW).GE.0.0_r8) THEN
       FLXI(I,J,KLOW)=0.5_r8*TEMPX(I,J,KLOW)*(Q(I,J,KLOW+1)+Q(I,J,KLOW))
       IF (POSITIVE) THEN
        FLXI(I,J,KLOW)=MIN(FLXI(I,J,KLOW),dz*rho(KLOW)*Q(I,J,KLOW) &
                      /(fnt(KLOW)*dt))
       ENDIF
      ELSE
       FLXI(I,J,KLOW) = 0.5_r8*TEMPX(I,J,KLOW)*(Q(I,J,KLOW+1)+Q(I,J,KLOW))
      ENDIF

!*********************************
      ELSE
!*********************************

      DO K=KLOW,NK1
       UPI(I,J,K)=0.5_r8*(TEMPX(I,J,K)+ABS(TEMPX(I,J,K)))
       UMI(I,J,K)=0.5_r8*(TEMPX(I,J,K)-ABS(TEMPX(I,J,K)))
      ENDDO
      DO K=KLOW,NK1
       UPSRI(I,J,K)=SQRT(UPI(I,J,K))
       UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
      ENDDO

      DO K=KLOW+1,NK1-1

!       -------------------
        IF (POSITIVE) THEN
!       -------------------

         IF (USE_MU) THEN
           MU_P    = fnz(K)*ABS(W3D(I,J,K))*dt/dz
           ALPHA_P = (1.0_r8 + MU_P)/6.0_r8
           COEF_P  = (1.0_r8 - 2.0_r8*ALPHA_P)/(2.0_r8*ALPHA_P)

           MU_M    = fnz(K+1)*ABS(W3D(I,J,K+1))*dt/dz
           ALPHA_M = (1.0_r8 + MU_M)/6.0_r8
           COEF_M  = (1.0_r8 - 2.0_r8*ALPHA_M)/(2.0_r8*ALPHA_M)
         ELSE
           COEF_P  =  2.0_r8
           COEF_M  =  2.0_r8
         ENDIF

         a0 = MAX(Q(I,J,K),0.0_r8)*MAX(Q(I,J,K+1),0.0_r8)

         GAMMA_P     = (P_Z(I,J,K)**2)/(P_Z(I,J,K)**2 + a0)
         GAMMA_HAT_P = GAMMA_P
         BETA_P      = 1.0_r8 + COEF_P*GAMMA_P
         BETA_HAT_P  = 1.0_r8 - GAMMA_HAT_P

         GAMMA_M     = (P_Z(I,J,K+1)**2)/(P_Z(I,J,K+1)**2 + a0)
         GAMMA_HAT_M = GAMMA_M
         BETA_M      = 1.0_r8 + COEF_M*GAMMA_M
         BETA_HAT_M  = 1.0_r8 - GAMMA_HAT_M

!       -------------------
        ELSE
!       -------------------

         BETA_P      = 1.0_r8
         BETA_HAT_P  = 1.0_r8
         BETA_M      = 1.0_r8
         BETA_HAT_M  = 1.0_r8

!       -------------------
        ENDIF
!       -------------------
       
       IF (PRESENT(TPHGT).AND.ZZ(K).GE.TPHGT) THEN
         FLXI(I,J,K) = 0.5_r8*TEMPX(I,J,K)*(Q(I,J,K+1)+Q(I,J,K)) 
       ELSE 
       FLXI(I,J,K) = &
            0.5_r8*TEMPX(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))                         &
       - aladv*(UPI(I,J,K)*BETA_P*(Q(I,J,K+1)-Q(I,J,K))                       &
             -UPSRI(I,J,K)*UPSRI(I,J,K-1)*BETA_HAT_P*(Q(I,J,K)-Q(I,J,K-1))    &
             +UMI(I,J,K)*BETA_M*(Q(I,J,K)-Q(I,J,K+1))                         &
             +UMSRI(I,J,K)*UMSRI(I,J,K+1)*BETA_HAT_M*(Q(I,J,K+1)-Q(I,J,K+2))) &
             /6.0_r8
       ENDIF      
      ENDDO

!     LEVEL: K=NK1
!=====================================
      IF(TEMPX(I,J,NK1).GE.0.0_r8) THEN
!=====================================
!       -------------------
        IF (POSITIVE) THEN
!       -------------------
         IF (USE_MU) THEN
           MU_P    = fnz(NK1)*ABS(W3D(I,J,NK1))*dt/dz
           ALPHA_P = (1.0_r8 + MU_P)/6.0_r8
           COEF_P  = (1.0_r8 - 2.0_r8*ALPHA_P)/(2.0_r8*ALPHA_P)
         ELSE
           COEF_P  =  2.0_r8
         ENDIF

         a0 = MAX(Q(I,J,NK1),0.0_r8)*MAX(Q(I,J,NK2),0.0_r8)

         GAMMA_P     = P_Z(I,J,NK1)**2/(P_Z(I,J,NK1)**2+a0)
         GAMMA_HAT_P = GAMMA_P
         BETA_P      = 1.0_r8 + COEF_P*GAMMA_P
         BETA_HAT_P  = 1.0_r8 - GAMMA_HAT_P
!       -------------------
        ELSE
!       -------------------

         BETA_P      = 1.0_r8
         BETA_HAT_P  = 1.0_r8
!       -------------------
        ENDIF
!       -------------------
         
         IF (PRESENT(TPHGT).AND.ZZ(NK1).GE.TPHGT) THEN
           FLXI(I,J,NK1) = 0.5_r8*TEMPX(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
         ELSE 
         FLXI(I,J,NK1) = &
            0.5_r8*TEMPX(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))                       &
          - aladv*(UPI(I,J,NK1)*BETA_P*(Q(I,J,NK2)-Q(I,J,NK1))                        &
               -UPSRI(I,J,NK1)*UPSRI(I,J,NK1-1)*BETA_HAT_P*(Q(I,J,NK1)-Q(I,J,NK1-1))) &
               /6.0_r8
         ENDIF      
!=====================================
      ELSE
!=====================================

         FLXi(I,J,NK1)=0.5_r8*TEMPX(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
         IF (POSITIVE) THEN
          FLXI(I,J,NK1)=MAX(FLXI(I,J,NK1),-dz*rho(NK2)*Q(I,J,NK2) &
                       /(fnt(NK2)*dt))
         ENDIF
!=====================================
      ENDIF
!=====================================

!     LEVEL: K=KLOW
!=====================================
      IF(TEMPX(I,J,KLOW).GE.0.0_r8) THEN
!=====================================
         FLXI(I,J,KLOW)=0.5_r8*TEMPX(I,J,KLOW)*(Q(I,J,KLOW+1)+Q(I,J,KLOW))
         IF (POSITIVE) THEN
          FLXI(I,J,KLOW)=MIN(FLXI(I,J,KLOW),dz*rho(KLOW)*Q(I,J,KLOW) &
                        /(fnt(KLOW)*dt))
         ENDIF
!=====================================
      ELSE
!=====================================
!       -------------------
        IF (POSITIVE) THEN
!       -------------------
         IF (USE_MU) THEN
           MU_M    = fnz(KLOW+1)*ABS(W3D(I,J,KLOW+1))*dt/dz
           ALPHA_M = (1.0_r8 + MU_M)/6.0_r8
           COEF_M  = (1.0_r8 - 2.0_r8*ALPHA_M)/(2.0_r8*ALPHA_M)
         ELSE
           COEF_M  =  2.0_r8
         ENDIF

         a0 = MAX(Q(I,J,KLOW),0.0_r8)*MAX(Q(I,J,KLOW+1),0.0_r8)

         GAMMA_M     = P_Z(I,J,KLOW+1)**2/(P_Z(I,J,KLOW+1)**2+a0)
         GAMMA_HAT_M = GAMMA_M
         BETA_M      = 1.0_r8 + COEF_M*GAMMA_M
         BETA_HAT_M  = 1.0_r8 - GAMMA_HAT_M
!       -------------------
        ELSE
!       -------------------
         BETA_M      = 1.0_r8
         BETA_HAT_M  = 1.0_r8
!       -------------------
        ENDIF  ! POSITIVE
!       -------------------

         FLXI(I,J,KLOW) = &
           0.5_r8*TEMPX(I,J,KLOW)*(Q(I,J,KLOW+1)+Q(I,J,KLOW))         &
         - aladv*(UMI(I,J,KLOW)*BETA_M*(Q(I,J,KLOW)-Q(I,J,KLOW+1))    &
               +UMSRI(I,J,KLOW)*UMSRI(I,J,KLOW+1)*BETA_HAT_M*(Q(I,J,KLOW+1)-Q(I,J,KLOW+2))) &
               /6.0_r8
!=====================================
      ENDIF
!=====================================

!*********************************
      ENDIF  ! DIV_CAL
!*********************************

      FLXI(I,J,KLOW-1) = 0.0_r8
      FLXI(I,J,NK2) = 0.0_r8

      DO K=KLOW,NK2
       TERMA(I,J,K)=TERMA(I,J,K)-(FLXI(I,J,K)-FLXI(I,J,K-1))*fnt(K)/dz
      ENDDO

      DO K=KLOW,NK2
       TERMA(I,J,K)=TERMA(I,J,K)/rho(K)
      ENDDO

      IF (PRESENT(TERM_ADD)) THEN

      DO K=KLOW,NK2
       TERMA(I,J,K)=TERMA(I,J,K)-Q(I,J,K)*TERM_ADD(I,J,K)
      ENDDO

      ENDIF

      DO K=1,KLOW-1
       TERMA(I,J,K) = 0.0_r8
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

   END SUBROUTINE advec_3d

END MODULE advec_3d_module
