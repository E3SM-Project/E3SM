MODULE q3d_rtime_module
! Contains the program that calculates the relaxation time scale (CRM fields toward GCM)

      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
      USE vvm_data_types, only: channel_t
      
      USE parmsld,      only: nVGCM_seg,netsz,nk1,nk2,nk3
      USE constld,      only: d0_0,d0_5,dx_gcm 
      
IMPLICIT NONE
PRIVATE

PUBLIC :: cal_rtime

CONTAINS

!================================================================================
   SUBROUTINE CAL_RTIME (channel)
!================================================================================
!  Prepare tau_rx (time-scale of relaxation)
!  channel width is fixed as 1 (i.e., # of prognostic grid = 1: chn = 1).
!
!  PLAN to do:
!  Should it be the wind in RLL coordinates (vector handling)?
!  TOPOGRAPHY (coupling strength) needs to be recovered.
!-------------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel  ! Channel data

#ifdef FIXTHIS

      ! Local 
      LOGICAL :: INTERPOL = .FALSE.
      INTEGER mi1,mim,mip,mj1,mjm,mjp
      INTEGER k,num_seg,num_gcm

      REAL (KIND=dbl_kind), DIMENSION(nk2,0:nVGCM_seg+1,4) :: MN_SPEED
      REAL (KIND=dbl_kind), DIMENSION(nk2,0:nVGCM_seg+1,4) :: MN_SPEED_Z
      REAL (KIND=dbl_kind), DIMENSION(1,0:nVGCM_seg+1,4)   :: MN_SPEED_ZZ
        
!=======================================================================   
      DO num_seg = 1, 4
!=======================================================================
      
      mi1 = channel%seg(num_seg)%mi1     ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1     ! y-size of channel segment (temp setting)

      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip 
      
      mjm = channel%seg(num_seg)%mjm 
      mjp = channel%seg(num_seg)%mjp   
      
      ! Calculate mean speed at vGCM point      
      CALL CSECT_AVG (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                      mi1,mim,mip,mj1,mjm,mjp,channel%seg(num_seg)%U3DX,   &
                      channel%seg(num_seg)%U3DY,MN_SPEED(:,1:nVGCM_seg,num_seg))

      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1     
        MN_SPEED_Z(K,num_gcm,num_seg) = &
                D0_5*(MN_SPEED(K,num_gcm,num_seg)+MN_SPEED(K+1,num_gcm,num_seg))
       ENDDO 
        MN_SPEED_Z(1,num_gcm,num_seg)   = D0_0
        MN_SPEED_Z(nk2,num_gcm,num_seg) = D0_0
      ENDDO
      
      DO num_gcm = 1,nVGCM_seg    
        MN_SPEED_ZZ(1,num_gcm,num_seg) = MN_SPEED(nk2,num_gcm,num_seg)
      ENDDO      
              
!=======================================================================   
      ENDDO   ! num_seg 
!=======================================================================

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
      
      ! TAU_RX_TQ
      CALL CSECT_DIST (INTERPOL,D0_0,D0_0,channel%seg(num_seg)%lcen,        &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                       mi1,mj1,nk1,MN_SPEED(2:nk2,:,num_seg),               &
                       channel%seg(num_seg)%TAU_RX_TQ)
       
      ! TAU_RX_ZX
      CALL CSECT_DIST (INTERPOL,D0_0,D0_5,channel%seg(num_seg)%lcen,        &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                       mi1,mj1,nk1-1,MN_SPEED_Z(2:nk1,:,num_seg),           &
                       channel%seg(num_seg)%TAU_RX_ZX)
       
      ! TAU_RX_ZY
      CALL CSECT_DIST (INTERPOL,D0_5,D0_0,channel%seg(num_seg)%lcen,        &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                       mi1,mj1,nk1-1,MN_SPEED_Z(2:nk1,:,num_seg),           &
                       channel%seg(num_seg)%TAU_RX_ZY)
                           
      ! TAU_RX_ZZ
      CALL CSECT_DIST (INTERPOL,D0_5,D0_5,channel%seg(num_seg)%lcen,        &
                       channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend, &
                       mi1,mj1,1,MN_SPEED_ZZ(:,:,num_seg),                  &
                       channel%seg(num_seg)%TAU_RX_ZZ)                     

      CONTAINS
      
       SUBROUTINE CSECT_AVG (lsta,lend,mi1,mim,mip,mj1,mjm,mjp,U,V,AVG)
!      input  (U & V) : U & V      
!      output (AVG)   : Mean speed average over a vGCM-cell at q-layers (k=2~nk2)  

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       integer, intent(in) :: mim,mip,mjm,mjp
       
       real (kind=dbl_kind), INTENT(IN)  :: U(mim:mip,mjm:mjp,nk3)
       real (kind=dbl_kind), INTENT(IN)  :: V(mim:mip,mjm:mjp,nk3)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k,chp
       real (kind=dbl_kind) :: SUM,TEMP
       
!---------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk2
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          TEMP = D0_5*(U(CHP,1,K)**2 + U(CHP-1,1,K)**2)   &
               + D0_5*(V(CHP,1,K)**2 + V(CHP,0,K)**2)
               
          SUM = SUM + SQRT(TEMP) 
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO   
!---------------------------------       
       else                  ! y-array
!---------------------------------   
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk2
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          TEMP = D0_5*(U(1,CHP,K)**2 + U(0,CHP,K)**2)   &
               + D0_5*(V(1,CHP,K)**2 + V(1,CHP-1,K)**2)
                
          SUM = SUM + SQRT(TEMP) 
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO  
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_avg 
      
      SUBROUTINE BOUND_G1 (kdim,A)
!     Filling the halos along the channel (Group1) 

      INTEGER,INTENT(IN) :: kdim     
      REAL (KIND=dbl_kind),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
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
      REAL (KIND=dbl_kind),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
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
      REAL (KIND=dbl_kind),DIMENSION(kdim,0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
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
       REAL (kind=dbl_kind), INTENT(IN) :: x_cor,y_cor  ! depending on the input variable location
       REAL (kind=dbl_kind), INTENT(IN) :: lcen(nVGCM_seg)  ! center CRM point of a vGCM cell 
       INTEGER, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       INTEGER, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       
       INTEGER, INTENT(IN) :: mi1               ! x-size of a channel-segment
       INTEGER, INTENT(IN) :: mj1               ! y-size of a channel-segment
       INTEGER, INTENT(IN) :: kdim              ! vertical size of a channel-segment
       
       REAL (kind=dbl_kind), INTENT(IN)  :: A(kdim,0:nVGCM_seg+1)
       REAL (kind=dbl_kind), INTENT(OUT) :: AOUT(mi1,mj1,kdim)
       
       ! Local
       INTEGER ng,k,chl,chn,npo
       REAL (KIND=dbl_kind) :: dist,dist1,dist2
       
       chn = 1
!---------------------------------       
       if (mi1.gt.mj1) then  
!---------------------------------  x-array     
       DO ng = 1, nVGCM_seg 

        DO chl = lsta(ng), lend(ng)

        IF (INTERPOL) THEN
         
         DIST1 = lcen(ng) - (float(chl) + x_cor)
         IF (DIST1.GE.D0_0) THEN
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
          IF (AOUT(chl,chn,K).NE.D0_0) THEN
            AOUT(chl,chn,K) = DX_GCM/AOUT(chl,chn,K)
          ENDIF
         ENDDO 
         
        ELSE   ! INTERPOL

         DO k = 1, kdim
           IF A(K,ng).EQ.D0_0) THEN
             AOUT(chl,chn,K) = D0_0
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

         DIST1 = lcen - (float(chl) + y_cor)
         IF (DIST1.GE.D0_0) THEN
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
          IF (AOUT(chn,chl,K).NE.D0_0) THEN
            AOUT(chn,chl,K) = DX_GCM/AOUT(chn,chl,K)
          ENDIF
         ENDDO 
         
        ELSE   ! INTERPOL        

         DO k = 1, kdim
           IF A(K,ng).EQ.D0_0) THEN
             AOUT(chn,chl,K) = D0_0
           ELSE
             AOUT(chn,chl,K) = DX_GCM/A(K,ng)
           ENDIF
         ENDDO 
        
        ENDIF  ! INTERPOL 
        
        ENDDO  ! chl-loop  
       ENDDO   ! ng-loop   
        
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_dist     

#endif
   END SUBROUTINE cal_rtime

END MODULE q3d_rtime_module
