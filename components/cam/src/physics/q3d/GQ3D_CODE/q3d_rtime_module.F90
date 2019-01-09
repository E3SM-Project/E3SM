MODULE q3d_rtime_module
! Contains the program that calculates the relaxation time scale (CRM fields toward GCM)

      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t

      USE parmsld,      only: nVGCM_seg,netsz,nk1,nk2,nk3
      USE constld,      only: dx_gcm,dy_gcm

IMPLICIT NONE
PRIVATE

PUBLIC :: cal_rtime

CONTAINS

! Local Subroutines:
!--------------------------------------------------------------------------------
! Subroutine cal_rtime : Calculate tau_rx (time-scale of relaxation) 
!
!      Subroutine csect_avg  : Average over each channel-segment 
!      Subroutine bound_g1   : Fill halos (face edge) for group-1 channels
!      Subroutine bound_g2   : Fill halos (face edge) for group-2 channels
!      Subroutine bound_g3   : Fill halos (face edge) for group-3 channels
!      Subroutine csect_dist : Distribute data on vGCM grid to each CRM grid 
!--------------------------------------------------------------------------------

!================================================================================
   SUBROUTINE CAL_RTIME (channel)
!================================================================================
!  Calculate tau_rx (time-scale of relaxation)
!  channel width is fixed as 1 (i.e., # of prognostic grid = 1: chn = 1).
!
!  JUNG: In the 1-D interpolation, mapping is not considered at the cube edge.
!
!  TOPOGRAPHY (coupling strength) needs to be recovered.
!-------------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel  ! Channel data

      ! Local
      LOGICAL :: INTERPOL = .FALSE.
      INTEGER mi1,mim,mip,mj1,mjm,mjp
      INTEGER k,num_seg,num_gcm

      REAL (KIND=r8), DIMENSION(nk2,0:nVGCM_seg+1,4) :: MN_SPEED
      REAL (KIND=r8), DIMENSION(nk2,0:nVGCM_seg+1,4) :: MN_SPEED_Z
      REAL (KIND=r8), DIMENSION(1,0:nVGCM_seg+1,4)   :: MN_SPEED_ZZ

!================================
      DO num_seg = 1, 4
!================================
      
      mi1 = channel%seg(num_seg)%mi1     ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1     ! y-size of channel segment (temp setting)

      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip

      mjm = channel%seg(num_seg)%mjm
      mjp = channel%seg(num_seg)%mjp

      ! Calculate mean speed at vGCM point
      CALL CSECT_AVG (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,      &
                      mi1,mim,mip,mj1,mjm,mjp,channel%seg(num_seg)%U3DX,        &
                      channel%seg(num_seg)%U3DY,channel%seg(num_seg)%AM_T,      &
                      MN_SPEED(:,1:nVGCM_seg,num_seg))

      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1
        MN_SPEED_Z(K,num_gcm,num_seg) = &
                0.5_r8*(MN_SPEED(K,num_gcm,num_seg)+MN_SPEED(K+1,num_gcm,num_seg))
       ENDDO
        MN_SPEED_Z(1,num_gcm,num_seg)   = 0.0_r8
        MN_SPEED_Z(nk2,num_gcm,num_seg) = 0.0_r8
      ENDDO

      DO num_gcm = 1,nVGCM_seg
        MN_SPEED_ZZ(1,num_gcm,num_seg) = MN_SPEED(nk2,num_gcm,num_seg)
      ENDDO

!================================
      ENDDO   ! num_seg
!================================

      IF (INTERPOL) THEN
      ! Filling the halos along the channel: no mapping

      SELECT CASE (channel%num_chg)
      CASE(1)
        CALL BOUND_G1 (nk1,MN_SPEED(2:nk2,:,:))
        CALL BOUND_G1 (nk1-1,MN_SPEED_Z(2:nk1,:,:))
        CALL BOUND_G1 (1,MN_SPEED_ZZ)

      CASE(2)
        CALL BOUND_G2 (nk1,MN_SPEED(2:nk2,:,:))
        CALL BOUND_G2 (nk1-1,MN_SPEED_Z(2:nk1,:,:))
        CALL BOUND_G2 (1,MN_SPEED_ZZ)

      CASE(3)
        CALL BOUND_G3 (nk1,MN_SPEED(2:nk2,:,:))
        CALL BOUND_G3 (nk1-1,MN_SPEED_Z(2:nk1,:,:))
        CALL BOUND_G3 (1,MN_SPEED_ZZ)
      END SELECT

      ENDIF

!================================
      DO num_seg = 1, 4
!================================
      ! TAU_RX_TQ
      CALL CSECT_DIST (INTERPOL,0.0_r8,0.0_r8,channel%seg(num_seg)%lcen,      &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,   &
                       mi1,mj1,nk1,MN_SPEED(2:nk2,0:nVGCM_seg+1,num_seg),     &
                       channel%seg(num_seg)%TAU_RX_TQ(:,:,2:nk2))

      ! TAU_RX_ZX
      CALL CSECT_DIST (INTERPOL,0.0_r8,0.5_r8,channel%seg(num_seg)%lcen,      &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,   &
                       mi1,mj1,nk1-1,MN_SPEED_Z(2:nk1,0:nVGCM_seg+1,num_seg), &
                       channel%seg(num_seg)%TAU_RX_ZX(:,:,2:nk1))

      ! TAU_RX_ZY
      CALL CSECT_DIST (INTERPOL,0.5_r8,0.0_r8,channel%seg(num_seg)%lcen,      &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,   &
                       mi1,mj1,nk1-1,MN_SPEED_Z(2:nk1,0:nVGCM_seg+1,num_seg), &
                       channel%seg(num_seg)%TAU_RX_ZY(:,:,2:nk1))

      ! TAU_RX_ZZ
      CALL CSECT_DIST (INTERPOL,0.5_r8,0.5_r8,channel%seg(num_seg)%lcen,    &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                       mi1,mj1,1,MN_SPEED_ZZ(1,0:nVGCM_seg+1,num_seg),      &
                       channel%seg(num_seg)%TAU_RX_ZZ)
!================================
      ENDDO   ! num_seg
!================================                     

      CONTAINS

       SUBROUTINE CSECT_AVG (lsta,lend,mi1,mim,mip,mj1,mjm,mjp,U,V,AM_T,AVG)
!      input  (U & V) : U & V
!      output (AVG)   : Mean speed average over a vGCM-cell at q-layers (k=2~nk2)
       
       USE utils, only: con2cov
       
       IMPLICIT NONE

       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment

       integer, intent(in) :: mim,mip,mjm,mjp

       real (kind=r8), INTENT(IN)  :: U(mim:mip,mjm:mjp,nk3)     ! contravariant wind component
       real (kind=r8), INTENT(IN)  :: V(mim:mip,mjm:mjp,nk3)     ! contravariant wind component
       real (kind=r8), INTENT(IN)  :: AM_T(4,mim:mip,mjm:mjp)    ! transformation matrix 
       
       real (kind=r8), INTENT(OUT) :: AVG(nk2,nVGCM_seg)

       ! Local
       integer n_gcm,k,npo,chl,chn
       real (kind=r8) :: VAL(netsz),TEMP,TEMP_U,TEMP_V,TEMP_U_cov,TEMP_V_cov 
       REAL (KIND=r8) :: VMAP(2)

       chn = 1
!---------------------------------
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------

       DO n_gcm = 1,nVGCM_seg
        DO k = 2,nk2
         DO npo = 1, netsz
          chl = lsta(n_gcm) - 1 + npo
          TEMP_U = 0.5_r8*(U(chl,chn,K)+U(chl-1,chn,K)) 
          TEMP_V = 0.5_r8*(V(chl,chn,K)+V(chl,chn-1,K))
          
          VMAP = CON2COV(AM_T(1,chl,chn),AM_T(2,chl,chn),AM_T(3,chl,chn),AM_T(4,chl,chn),  &
                         TEMP_U,TEMP_V)
                         
          TEMP_U_cov = VMAP(1)
          TEMP_V_cov = VMAP(2)              
          
          TEMP = TEMP_U*TEMP_U_cov + TEMP_V*TEMP_V_cov

          VAL(npo) = SQRT(TEMP)
         ENDDO
         AVG(k,n_gcm) = SUM(VAL)/FLOAT(netsz)
        ENDDO
       ENDDO
              
!---------------------------------
       else                  ! y-array
!---------------------------------

       DO n_gcm = 1,nVGCM_seg
        DO k = 2,nk2
         DO npo = 1, netsz
          chl = lsta(n_gcm) - 1 + npo
          TEMP_U = 0.5_r8*(U(chn,chl,K)+U(chn-1,chl,K)) 
          TEMP_V = 0.5_r8*(V(chn,chl,K)+V(chn,chl-1,K))

          VMAP = CON2COV(AM_T(1,chn,chl),AM_T(2,chn,chl),AM_T(3,chn,chl),AM_T(4,chn,chl),  &
                         TEMP_U,TEMP_V)
                         
          TEMP_U_cov = VMAP(1)
          TEMP_V_cov = VMAP(2)              
          
          TEMP = TEMP_U*TEMP_U_cov + TEMP_V*TEMP_V_cov
          
          VAL(npo) = SQRT(TEMP)
         ENDDO
         AVG(k,n_gcm) = SUM(VAL)/FLOAT(netsz)
        ENDDO
       ENDDO
       
!---------------------------------
       endif
!---------------------------------

       END SUBROUTINE csect_avg

      SUBROUTINE BOUND_G1 (kdim,A)
!     Filling the halos along the channel (Group1)

      INTEGER,INTENT(IN) :: kdim
      REAL (KIND=r8),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, kdim
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,nVGCM_seg,4)

          A(K,chl_p,2) = A(K,1,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,1,4)
          A(K,chl_m,3) = A(K,nVGCM_seg,2)

          A(K,chl_p,4) = A(K,1,1)
          A(K,chl_m,4) = A(K,nVGCM_seg,3)
        ENDDO

      END SUBROUTINE bound_g1

      SUBROUTINE BOUND_G2 (kdim,A)
!     Filling the halos along the channel (Group2)

      INTEGER,INTENT(IN) :: kdim
      REAL (KIND=r8),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, kdim
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,1,4)

          A(K,chl_p,2) = A(K,1,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,nVGCM_seg,4)
          A(K,chl_m,3) = A(K,nVGCM_seg,2)

          A(K,chl_p,4) = A(K,nVGCM_seg,3)
          A(K,chl_m,4) = A(K,1,1)
        ENDDO

      END SUBROUTINE bound_g2

      SUBROUTINE BOUND_G3 (kdim,A)
!     Filling the halos along the channel (Group3)

      INTEGER, INTENT(IN) :: kdim    ! vertical size of A
      REAL (KIND=r8),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A

      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0

        DO K = 1, kdim
          A(K,chl_p,1) = A(K,1,2)
          A(K,chl_m,1) = A(K,1,4)

          A(K,chl_p,2) = A(K,nVGCM_seg,3)
          A(K,chl_m,2) = A(K,nVGCM_seg,1)

          A(K,chl_p,3) = A(K,nVGCM_seg,2)
          A(K,chl_m,3) = A(K,nVGCM_seg,4)

          A(K,chl_p,4) = A(K,1,3)
          A(K,chl_m,4) = A(K,1,1)
        ENDDO

      END SUBROUTINE bound_g3

       SUBROUTINE CSECT_DIST (INTERPOL,x_cor,y_cor,lcen,lsta,lend,mi1,mj1,kdim,A,AOUT)
!      input  (A)    : MN_SPEED
!      output (AOUT) : Relaxation time scale

       IMPLICIT NONE

       LOGICAL, INTENT(IN) :: INTERPOL
       REAL (kind=r8), INTENT(IN) :: x_cor,y_cor  ! depending on the input variable location
       REAL (kind=r8), INTENT(IN) :: lcen(nVGCM_seg)  ! center CRM point of a vGCM cell
       INTEGER, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       INTEGER, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell

       INTEGER, INTENT(IN) :: mi1               ! x-size of a channel-segment
       INTEGER, INTENT(IN) :: mj1               ! y-size of a channel-segment
       INTEGER, INTENT(IN) :: kdim              ! vertical size of a channel-segment

       REAL (kind=r8), INTENT(IN)  :: A(kdim,0:nVGCM_seg+1)
       REAL (kind=r8), INTENT(OUT) :: AOUT(mi1,mj1,kdim)

       ! Local
       INTEGER ng,k,chl,chn,npo
       REAL (KIND=r8) :: dist,dist1,dist2

       chn = 1
!---------------------------------
       if (mi1.gt.mj1) then
!---------------------------------  x-array
       DO ng = 1, nVGCM_seg

        DO chl = lsta(ng), lend(ng)

        IF (INTERPOL) THEN

         DIST1 = lcen(ng) - (float(chl) + x_cor)
         IF (DIST1.GE.0.0_r8) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF

         DO k = 1, kdim
          AOUT(chl,chn,K) = (DIST*A(K,NPO)+DIST2*A(K,ng))/FLOAT(netsz)
         ENDDO

         DO k = 1, kdim
          IF (AOUT(chl,chn,K).NE.0.0_r8) THEN
            AOUT(chl,chn,K) = DX_GCM/AOUT(chl,chn,K)
          ENDIF
         ENDDO

        ELSE   ! INTERPOL

         DO k = 1, kdim
           IF (A(K,ng).EQ.0.0_r8) THEN
             AOUT(chl,chn,K) = 0.0_r8
           ELSE
             AOUT(chl,chn,K) = DX_GCM/A(K,ng)
           ENDIF
         ENDDO

        ENDIF  ! INTERPOL

        ENDDO  ! chl-loop
       ENDDO   ! ng-loop
!---------------------------------
       else                  ! y-array
!---------------------------------
       DO ng = 1, nVGCM_seg

        DO chl = lsta(ng), lend(ng)

        IF (INTERPOL) THEN

         DIST1 = lcen(ng) - (float(chl) + y_cor)
         IF (DIST1.GE.0.0_r8) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF

         DO k = 1, kdim
          AOUT(chn,chl,K) = (DIST*A(K,NPO)+DIST2*A(K,ng))/FLOAT(netsz)
         ENDDO

         DO k = 1, kdim
          IF (AOUT(chn,chl,K).NE.0.0_r8) THEN
            AOUT(chn,chl,K) = DY_GCM/AOUT(chn,chl,K)
          ENDIF
         ENDDO

        ELSE   ! INTERPOL

         DO k = 1, kdim
           IF (A(K,ng).EQ.0.0_r8) THEN
             AOUT(chn,chl,K) = 0.0_r8
           ELSE
             AOUT(chn,chl,K) = DY_GCM/A(K,ng)
           ENDIF
         ENDDO

        ENDIF  ! INTERPOL

        ENDDO  ! chl-loop
       ENDDO   ! ng-loop

!---------------------------------
       endif
!---------------------------------

       END SUBROUTINE csect_dist

   END SUBROUTINE cal_rtime

END MODULE q3d_rtime_module
