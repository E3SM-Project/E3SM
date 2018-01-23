MODULE q3d_effect_module
! Contains the programs for calculating CRM effects that will be provided to the GCM

      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
      USE vvm_data_types, only: channel_t
      
      USE parmsld, only: nVGCM_seg,netsz,ntracer,nk1,nk2,nk3
      USE constld, only: d0_0,d0_5,d0_25,dz,physics,fnt,rho,rhoz
      
IMPLICIT NONE
PRIVATE

PUBLIC :: ini_tendency,    & ! initialize tendencies at each CRM timestep
          ini_crm_eft,     & ! initialize for cal_eddy_trans & cal_tendency
          cal_eddy_trans,  & ! calculate mean eddy transport effects
          cal_tendency,    & ! calculate mean diabatic tendencies
          ini_sfc_eft,     & ! initialize for cal_sfc_eft
          cal_sfc_eft        ! calculate mean sfc fluxes

CONTAINS

!================================================================================
   SUBROUTINE INI_TENDENCY (channel)
!================================================================================
!  Initialize total diabatic tendencies at each CRM time step 
!------------------------------------------------------------------------
!  IN/OUT: FTH3D_DIA, FQV3D_DIA, FQT3D_DIA, FU_DIA, FV_DIA,
!          FQC3D_DIA, FQI3D_DIA, FQR3D_DIA, FQS3D_DIA, FQG3D_DIA (channel)
!------------------------------------------------------------------------      
      type(channel_t), intent(inout) :: channel   ! Channel data      

#ifndef CAM
      
      INTEGER num_seg

!===============================   
      DO num_seg = 1, 4
!=============================== 

      ! total diabatic effects 
      channel%seg(num_seg)%FTH3D_DIA = D0_0
      channel%seg(num_seg)%FQV3D_DIA = D0_0
      channel%seg(num_seg)%FQC3D_DIA = D0_0
      channel%seg(num_seg)%FQI3D_DIA = D0_0
      channel%seg(num_seg)%FQR3D_DIA = D0_0
      channel%seg(num_seg)%FQS3D_DIA = D0_0
      channel%seg(num_seg)%FQG3D_DIA = D0_0     
      channel%seg(num_seg)%FQT3D_DIA = D0_0
      
      channel%seg(num_seg)%FU_DIA = D0_0
      channel%seg(num_seg)%FV_DIA = D0_0
      
      ! Microphysical effects      
      channel%seg(num_seg)%THAD_MICRO = D0_0
      channel%seg(num_seg)%QVAD_MICRO = D0_0
      channel%seg(num_seg)%QCAD_MICRO = D0_0
      channel%seg(num_seg)%QIAD_MICRO = D0_0
      channel%seg(num_seg)%QSAD_MICRO = D0_0
      channel%seg(num_seg)%QGAD_MICRO = D0_0
      channel%seg(num_seg)%QRAD_MICRO = D0_0  
      
      !Turbulent effects      
      channel%seg(num_seg)%THAD3 = D0_0
      channel%seg(num_seg)%QVAD3 = D0_0
      channel%seg(num_seg)%QCAD3 = D0_0
      channel%seg(num_seg)%QIAD3 = D0_0
      channel%seg(num_seg)%QTAD3 = D0_0
      
      ! Radiative effects 
      channel%seg(num_seg)%FTHRAD = D0_0    
      
!=======================================================================   
      ENDDO   ! num_seg 
!=======================================================================
#endif

      END SUBROUTINE ini_tendency
      
!========================================================================
   SUBROUTINE INI_CRM_EFT(channel)
!========================================================================
! Called before the CRM time marching (initialization for time average)
! E1 - vertical eddy fluxes (finally, flux convergence)
! E2 - diabatic tendencies (physics)
!------------------------------------------------------------------------
! OUT: TH_E1, QV_E1, QT_E1, U_E1, V_E1, QC_E1, QI_E1, QR_E1, QS_E1, QG_E1,
!      TH_E2, QV_E2, QT_E2, U_E2, V_E2, QC_E2, QI_E2, QR_E2, QS_E2, QG_E2 
!      (channel)
!-------------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel   ! Channel data
      
      INTEGER num_seg,num_gcm,nt,k

#ifndef CAM

!===============================   
      DO num_seg = 1, 4
!=============================== 
      DO num_gcm = 1,nVGCM_seg   ! # of vGCM points in a segment
      
      DO k = 1,nk2
       channel%seg(num_seg)%TH_E1(k,num_gcm) = D0_0
       channel%seg(num_seg)%QV_E1(k,num_gcm) = D0_0

       DO nt = 1,ntracer
        channel%seg(num_seg)%QT_E1(k,num_gcm,nt) = D0_0
       ENDDO

       channel%seg(num_seg)%U_E1(k,num_gcm)  = D0_0
       channel%seg(num_seg)%V_E1(k,num_gcm)  = D0_0
      ENDDO
      
      IF (physics) THEN 
       DO k = 1,nk2
        channel%seg(num_seg)%QC_E1(k,num_gcm) = D0_0
        channel%seg(num_seg)%QI_E1(k,num_gcm) = D0_0
        channel%seg(num_seg)%QR_E1(k,num_gcm) = D0_0
        channel%seg(num_seg)%QS_E1(k,num_gcm) = D0_0
        channel%seg(num_seg)%QG_E1(k,num_gcm) = D0_0
       ENDDO     
      ENDIF      

      DO k = 1,nk2
       channel%seg(num_seg)%TH_E2(k,num_gcm) = D0_0
       channel%seg(num_seg)%QV_E2(k,num_gcm) = D0_0

       DO nt = 1,ntracer
        channel%seg(num_seg)%QT_E2(k,num_gcm,nt) = D0_0
       ENDDO

       channel%seg(num_seg)%U_E2(k,num_gcm)  = D0_0
       channel%seg(num_seg)%V_E2(k,num_gcm)  = D0_0
      ENDDO

      IF (physics) THEN 
       DO k = 1,nk2
        channel%seg(num_seg)%QC_E2(k,num_gcm) = D0_0
        channel%seg(num_seg)%QI_E2(k,num_gcm) = D0_0
        channel%seg(num_seg)%QR_E2(k,num_gcm) = D0_0
        channel%seg(num_seg)%QS_E2(k,num_gcm) = D0_0
        channel%seg(num_seg)%QG_E2(k,num_gcm) = D0_0
       ENDDO     
      ENDIF 
      
      ENDDO   ! num_gcm
      
!===============================
      ENDDO   ! num_seg 
!===============================         

#endif

   END SUBROUTINE ini_crm_eft

!================================================================================
   SUBROUTINE CAL_EDDY_TRANS (INUM,MAKE_AVG,channel)
!================================================================================
!  Calculate the vertical eddy flux & finally the vertical flux convergence
!  Average the fluxes over space (channel) and time (a GCM timestep)
!------------------------------------------------------------------------
!  IN:  TH3D_ED, QV3D_ED, QT3D_ED, U3DX_ED, U3DY_ED, W3D_ED, 
!       QC3D_ED, QI3D_ED, QR3D_ED, QS3D_ED, QG3D_ED (channel)
!
!  IN/OUT: TH_E1, QV_E1, QT_E1, U_E1, V_E1,
!          QC_E1, QI_E1, QR_E1, QS_E1, QG_E1 (channel)
!
!  channel width is fixed as 1 (i.e., # of prognostic grid =1).
!
!  PLAN to do:
!  Mapping should be included for vector components (u & v).
!  The modification for "TOPOGRAPHY" is removed and should be added.
!------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: INUM      ! # of data for time averaging
      LOGICAL, INTENT(IN) :: MAKE_AVG  ! if true, make an average
      
      type(channel_t), intent(inout) :: channel   ! Channel data

#ifndef CAM

      integer mi1,mj1
      integer k,nt,num_seg,num_gcm

      REAL (KIND=dbl_kind) :: AVG(nk2,nVGCM_seg)
      
!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1 = channel%seg(num_seg)%mi1     ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1     ! y-size of channel segment (temp setting)

!----------------------------------------------------------------------
!     1. Calculate vertical eddy fluxes averaged over a vGCM cell-size   
!---------------------------------------------------------------------- 
     
      ! EDDY FLUX: w'(TH)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%TH3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%TH_E1(K,num_gcm) = &
               channel%seg(num_seg)%TH_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO  
      
      ! EDDY FLUX: w'(QV)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QV3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QV_E1(K,num_gcm) = &
               channel%seg(num_seg)%QV_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
             
      IF (physics) THEN 

      ! EDDY FLUX: w'(QC)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QC3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QC_E1(K,num_gcm) = &
               channel%seg(num_seg)%QC_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ! EDDY FLUX: w'(QI)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QI3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QI_E1(K,num_gcm) = &
               channel%seg(num_seg)%QI_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ! EDDY FLUX: w'(QR)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QR3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QR_E1(K,num_gcm) = &
               channel%seg(num_seg)%QR_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ! EDDY FLUX: w'(QS)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QS3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QS_E1(K,num_gcm) = &
               channel%seg(num_seg)%QS_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ! EDDY FLUX: w'(QG)'
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QG3D_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QG_E1(K,num_gcm) = &
               channel%seg(num_seg)%QG_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
                              
      ENDIF  ! PHYSICS
      
      ! EDDY FLUX: w'(QT)'
      DO nt = 1,ntracer
      
      CALL CSECT_AVG_Q (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%QT3D_ED(:,:,:,nt), &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%QT_E1(K,num_gcm,nt) = &
               channel%seg(num_seg)%QT_E1(K,num_gcm,nt) + AVG(K,num_gcm)
       ENDDO
      ENDDO  
      
      ENDDO    

      ! EDDY FLUX: w'(u)'
      CALL CSECT_AVG_U (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%U3DX_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%U_E1(K,num_gcm) = &
               channel%seg(num_seg)%U_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ! EDDY FLUX: w'(v)'
      CALL CSECT_AVG_V (channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%U3DY_ED,  &
                        channel%seg(num_seg)%W3D_ED,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk1                
         channel%seg(num_seg)%V_E1(K,num_gcm) = &
               channel%seg(num_seg)%V_E1(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO              
      
!******************************************************
      IF (MAKE_AVG) THEN
!******************************************************

!-----------------------------------------------------------------
!     2. Calculate mean eddy fluxes (average over a vGCM timestep)
!-----------------------------------------------------------------

      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk1
        channel%seg(num_seg)%TH_E1(K,num_gcm)= &  
             channel%seg(num_seg)%TH_E1(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QV_E1(K,num_gcm)= &
             channel%seg(num_seg)%QV_E1(K,num_gcm)/FLOAT(INUM)
       ENDDO
      ENDDO 

      IF (physics) THEN 
      
      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk1
        channel%seg(num_seg)%QC_E1(K,num_gcm)= &
             channel%seg(num_seg)%QC_E1(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QI_E1(K,num_gcm)= &
             channel%seg(num_seg)%QI_E1(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QR_E1(K,num_gcm)= &
             channel%seg(num_seg)%QR_E1(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QS_E1(K,num_gcm)= &
             channel%seg(num_seg)%QS_E1(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QG_E1(K,num_gcm)= &
             channel%seg(num_seg)%QV_E1(K,num_gcm)/FLOAT(INUM)                                                    
       ENDDO
      ENDDO 
      
      ENDIF
      
      DO NT=1,ntracer
       DO num_gcm = 1,nVGCM_seg
        DO K = 2,nk1
         channel%seg(num_seg)%QT_E1(K,num_gcm,nt)= &
             channel%seg(num_seg)%QT_E1(K,num_gcm,nt)/FLOAT(INUM)
        ENDDO
       ENDDO         
      ENDDO
         
      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk1
        channel%seg(num_seg)%U_E1(K,num_gcm)= &
             channel%seg(num_seg)%U_E1(K,num_gcm)/FLOAT(INUM)
         
        channel%seg(num_seg)%V_E1(K,num_gcm)= &
             channel%seg(num_seg)%V_E1(K,num_gcm)/FLOAT(INUM)
       ENDDO
      ENDDO 
      
!-----------------------------------------------------------------
!     3. Calculate the mean vertical flux convergence
!-----------------------------------------------------------------

      ! TH
      CALL VERT_CONV(channel%seg(num_seg)%TH_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%TH_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      ! QV
      CALL VERT_CONV(channel%seg(num_seg)%QV_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QV_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      IF (physics) THEN 

      ! QC
      CALL VERT_CONV(channel%seg(num_seg)%QC_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QC_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      ! QI
      CALL VERT_CONV(channel%seg(num_seg)%QI_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QI_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      ! QR
      CALL VERT_CONV(channel%seg(num_seg)%QR_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QR_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      ! QS
      CALL VERT_CONV(channel%seg(num_seg)%QS_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QS_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  

      ! QG
      CALL VERT_CONV(channel%seg(num_seg)%QG_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QG_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  
      
      ENDIF

      ! QT      
      DO nt = 1, ntracer

      CALL VERT_CONV(channel%seg(num_seg)%QT_E1(:,:,nt),AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%QT_E1(K,num_gcm,nt) = AOUT(K,num_gcm)
       ENDDO
      ENDDO        
      
      ENDDO

      ! U
      CALL VERT_CONV(channel%seg(num_seg)%U_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%U_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO  
      
      ! V
      CALL VERT_CONV(channel%seg(num_seg)%V_E1,AOUT)
                     
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2                
         channel%seg(num_seg)%V_E1(K,num_gcm) = AOUT(K,num_gcm)
       ENDDO
      ENDDO        
            
!********************************************
      ENDIF      ! MAKE_AVG
!********************************************

!===============================
      ENDDO   ! num_seg 
!=============================== 

      CONTAINS

       SUBROUTINE CSECT_AVG_Q(lsta,lend,mi1,mj1,A,W,AVG)
!      input  (A)   : eddy field (q) transported by w'
!             (W)   : w' 
!      output (AVG) : vertical flux average over a vGCM-cell at w-level (k=2~nk1)  

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       real (kind=dbl_kind), INTENT(IN)  :: A(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(IN)  :: W(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k,chp
       real (kind=dbl_kind) :: SUM
       
!---------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_5*(A(CHP,1,K)+A(CHP,1,K+1))*W(CHP,1,K)
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO      
!---------------------------------       
       else                 ! y-array
!---------------------------------   
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_5*(A(1,CHP,K)+A(1,CHP,K+1))*W(1,CHP,K)
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO  
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_avg_q  
       
       SUBROUTINE CSECT_AVG_U(lsta,lend,mi1,mj1,A,W,AVG)
!      input  (A)   : eddy field (U) transported by w'
!             (W)   : w' 
!      output (AVG) : vertical flux average over a vGCM-cell at w-level (k=2~nk1) 

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       real (kind=dbl_kind), INTENT(IN)  :: A(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(IN)  :: W(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k,chp
       real (kind=dbl_kind) :: SUM
       
!---------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_25*(A(CHP,1,K)+A(CHP,1,K+1) &
                            +A(CHP-1,1,K)+A(CHP-1,1,K+1))*W(CHP,1,K) 
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO                       
!---------------------------------       
       else                  ! y-array
!---------------------------------   
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_25*(A(1,CHP,K)+A(1,CHP,K+1) &
                            +A(0,CHP,K)+A(0,CHP,K+1))*W(1,CHP,K)   
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO  
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_avg_u  
       
       SUBROUTINE CSECT_AVG_V(lsta,lend,mi1,mj1,A,W,AVG)
!      input  (A)   : eddy field (V) transported by w'
!             (W)   : w' 
!      output (AVG) : vertical flux average over a vGCM-cell at w-level (k=2~nk1) 

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       real (kind=dbl_kind), INTENT(IN)  :: A(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(IN)  :: W(0:mi1,0:mj1,nk2)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k,chp
       real (kind=dbl_kind) :: SUM
       
!---------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_25*(A(CHP,1,K)+A(CHP,1,K+1) &
                            +A(CHP,0,K)+A(CHP,0,K+1))*W(CHP,1,K)
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO                 
!---------------------------------       
       else                 ! y-array
!---------------------------------   
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk1
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + D0_25*(A(1,CHP,K)+A(1,CHP,K+1) &
                            +A(1,CHP-1,K)+A(1,CHP-1,K+1))*W(1,CHP,K)  
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO  
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_avg_v   
       
       SUBROUTINE VERT_CONV(A,AOUT)
       
!      input  (A)   : eddy fluxes at level position
!      output (AOUT): vertical flux convergence at layer position    
!      DZ: grid size in z-(vertical) direction (m)

       IMPLICIT NONE            
      
       real (kind=dbl_kind), INTENT(IN)  :: A(nk2,nVGCM_seg)
       real (kind=dbl_kind), INTENT(OUT) :: AOUT(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k
       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk2
         AOUT(k,n_gcm) = &
          -FNT(K)*(RHOZ(K)*A(k,n_gcm)-RHOZ(K-1)*A(k-1,n_gcm))/(RHO(K)*DZ)
        ENDDO
       ENDDO  
             
       END SUBROUTINE vert_conv                     
                 
#endif

   END SUBROUTINE cal_eddy_trans

!================================================================================
   SUBROUTINE CAL_TENDENCY (INUM,MAKE_AVG,channel)
!================================================================================
!  calculate the mean diabatic effects due to physical processes
!  averaging over space (a vGCM cell size) and time (a vGCM timestep)
!------------------------------------------------------------------------
!  IN:  FTH3D_DIA, FQV3D_DIA, FQT3D_DIA, FU_DIA, FV_DIA, 
!       FQC3D_DIA, FQI3D_DIA, FQR3D_DIA, FQS3D_DIA, FQG3D_DIA (channel)        
!
!  IN/OUT: TH_E2, QV_E2, QT_E2, U_E2, V_E2,
!          QC_E2, QI_E2, QR_E2, QS_E2, QG_E2 (channel)
!
!  channel width is fixed as 1 (i.e., # of prognostic grid =1).
!
!  PLAN to do:
!  Mapping should be included for vector components (u & v).
!  The modification for "TOPOGRAPHY" is removed and should be added.
!------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: INUM       ! # of data for time averaging
      LOGICAL, INTENT(IN) :: MAKE_AVG   ! if true, make an average
      
      type(channel_t), intent(inout) :: channel   ! Channel data      

#ifndef CAM
      integer mi1,mj1
      integer k,nt,num_seg,num_gcm

      REAL (KIND=dbl_kind) :: AVG(nk2,nVGCM_seg)

!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1 = channel%seg(num_seg)%mi1     ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1     ! y-size of channel segment (temp setting)
      
      ! Diabatic effect: TH
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FTH3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%TH_E2(K,num_gcm) = &
               channel%seg(num_seg)%TH_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO 
      
      ! Diabatic effect: QV
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQV3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QV_E2(K,num_gcm) = &
               channel%seg(num_seg)%QV_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
     
      IF (physics) THEN

      ! Diabatic effect: QC
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQC3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QC_E2(K,num_gcm) = &
               channel%seg(num_seg)%QC_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

      ! Diabatic effect: QI
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQI3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QI_E2(K,num_gcm) = &
               channel%seg(num_seg)%QI_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

      ! Diabatic effect: QR
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQR3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QR_E2(K,num_gcm) = &
               channel%seg(num_seg)%QR_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

      ! Diabatic effect: QS
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQS3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QS_E2(K,num_gcm) = &
               channel%seg(num_seg)%QS_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

      ! Diabatic effect: QG
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQG3D_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QG_E2(K,num_gcm) = &
               channel%seg(num_seg)%QG_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        
      
      ENDIF       

      DO nt=1,ntracer
      
      ! Diabatic effect: QT
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FQT3D_DIA(:,:,:,nt),AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%QT_E2(K,num_gcm,nt) = &
               channel%seg(num_seg)%QT_E2(K,num_gcm,nt) + AVG(K,num_gcm)
       ENDDO
      ENDDO  
      
      ENDDO

      ! Diabatic effect: U (U3DX)
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FU_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%U_E2(K,num_gcm) = &
               channel%seg(num_seg)%U_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

      ! Diabatic effect: V (U3DY)
      CALL CSECT_AVG_3D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%FV_DIA,AVG)
       
      DO num_gcm = 1,nVGCM_seg
       DO K = 2,nk2               
         channel%seg(num_seg)%V_E2(K,num_gcm) = &
               channel%seg(num_seg)%V_E2(K,num_gcm) + AVG(K,num_gcm)
       ENDDO
      ENDDO        

!******************************************************
       IF (MAKE_AVG) THEN
!      Make time averages over a timestep of GCM
!******************************************************
      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk2
        channel%seg(num_seg)%TH_E2(K,num_gcm)= &  
             channel%seg(num_seg)%TH_E2(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QV_E2(K,num_gcm)= &
             channel%seg(num_seg)%QV_E2(K,num_gcm)/FLOAT(INUM)
       ENDDO
      ENDDO 

      IF (physics) THEN 
      
      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk2
        channel%seg(num_seg)%QC_E2(K,num_gcm)= &
             channel%seg(num_seg)%QC_E2(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QI_E2(K,num_gcm)= &
             channel%seg(num_seg)%QI_E2(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QR_E2(K,num_gcm)= &
             channel%seg(num_seg)%QR_E2(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QS_E2(K,num_gcm)= &
             channel%seg(num_seg)%QS_E2(K,num_gcm)/FLOAT(INUM)
        channel%seg(num_seg)%QG_E2(K,num_gcm)= &
             channel%seg(num_seg)%QV_E2(K,num_gcm)/FLOAT(INUM)                                                    
       ENDDO
      ENDDO 
      
      ENDIF
      
      DO NT=1,ntracer
       DO num_gcm = 1,nVGCM_seg
        DO K = 2,nk2
         channel%seg(num_seg)%QT_E2(K,num_gcm,nt)= &
             channel%seg(num_seg)%QT_E2(K,num_gcm,nt)/FLOAT(INUM)
        ENDDO
       ENDDO         
      ENDDO
         
      DO num_gcm = 1,nVGCM_seg 
       DO K = 2,nk2
        channel%seg(num_seg)%U_E2(K,num_gcm)= &
             channel%seg(num_seg)%U_E2(K,num_gcm)/FLOAT(INUM)
         
        channel%seg(num_seg)%V_E2(K,num_gcm)= &
             channel%seg(num_seg)%V_E2(K,num_gcm)/FLOAT(INUM)
       ENDDO
      ENDDO 
!********************************************
      ENDIF      ! MAKE_AVG
!********************************************

!=======================================================================   
      ENDDO   ! num_seg 
!=======================================================================
      
      CONTAINS
      
       SUBROUTINE CSECT_AVG_3D(lsta,lend,mi1,mj1,A,AVG)
!      input  (A)   : 3D target variable        
!      output (AVG) : Average over a vGCM-cell at q, u, v - layers (k=2~nk2)  

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       real (kind=dbl_kind), INTENT(IN)  :: A(mi1,mj1,nk2)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nk2,nVGCM_seg)
       
       ! Local
       integer n_gcm,k,chp
       real (kind=dbl_kind) :: SUM
       
!---------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!---------------------------------       
       DO n_gcm = 1,nVGCM_seg 
        DO k = 2,nk2
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + A(CHP,1,K)
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
          SUM = SUM + A(1,CHP,K)
         ENDDO       
         AVG(k,n_gcm) = SUM/FLOAT(netsz)
        ENDDO
       ENDDO  
!---------------------------------       
       endif
!---------------------------------                
       
       END SUBROUTINE csect_avg_3d 
#endif
   END SUBROUTINE cal_tendency

!========================================================================
   SUBROUTINE INI_SFC_EFT(channel)
!========================================================================
!  Called before the CRM time marching (initialization for time average)
!  surface precipitation and surface fluxes
!------------------------------------------------------------------------
!  OUT: SPREC_E0, WTH_E0, WQV_E0, UW_E0, WV_E0 (channel)
!
!       Add more later ?
!-------------------------------------------------------------------------
      type(channel_t), intent(inout) :: channel   ! Channel data
      
#ifndef CAM

      ! Local
      integer num_seg,num_gcm

      IF (physics) THEN

!===============================   
      DO num_seg = 1, 4
!=============================== 
      DO num_gcm = 1,nVGCM_seg   ! # of vGCM points in a segment 
       channel%seg(num_seg)%SPREC_E0(num_gcm) = D0_0

       channel%seg(num_seg)%WTH_E0(num_gcm) = D0_0
       channel%seg(num_seg)%WQV_E0(num_gcm) = D0_0
       channel%seg(num_seg)%UW_E0(num_gcm)  = D0_0
       channel%seg(num_seg)%WV_E0(num_gcm)  = D0_0
      ENDDO 
!===============================   
      ENDDO
!=============================== 

      ENDIF  ! PHYSICS

#endif

   END SUBROUTINE ini_sfc_eft

!================================================================================
   SUBROUTINE CAL_SFC_EFT (INUM,MAKE_AVG,channel)
!================================================================================
!  calculate the surface precipitation & surface fluxes
!  Average (over a channel segment & a GCM timestep) is provided to a GCM grid
!------------------------------------------------------------------------
!  IN: SPREC, WTH, WQV, UW, WV (channel)
!
!  IN/OUT: SPREC_E0, WTH_E0, WQV_E0, UW_E0, WV_E0 (channel)
!
!  channel width is fixed as 1 (i.e., # of prognostic grid =1).
!
!  PLAN to do:
!  Mapping should be included for vector components (u & v).
!------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: INUM      ! # of data for time averaging
      LOGICAL, INTENT(IN) :: MAKE_AVG  ! if true, make an average

      type(channel_t), intent(inout) :: channel   ! Channel data
      
#ifndef CAM
      
      ! Local
      integer mi1,mj1
      integer num_seg,num_gcm
      REAL (KIND=dbl_kind) :: AVG(nVGCM_seg)

      IF (physics) THEN
      
!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment (temp setting)
         
!     Make space average over a channel segment

      ! SPREC
      CALL CSECT_AVG_2D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%SPREC,AVG)
                        
      DO num_gcm = 1,nVGCM_seg
       channel%seg(num_seg)%SPREC_E0(num_gcm) = &
         channel%seg(num_seg)%SPREC_E0(num_gcm) + AVG(num_gcm)
      ENDDO 
      
      ! WTH
      CALL CSECT_AVG_2D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%WTH,AVG)
                        
      DO num_gcm = 1,nVGCM_seg
       channel%seg(num_seg)%WTH_E0(num_gcm) = &
         channel%seg(num_seg)%WTH_E0(num_gcm) + AVG(num_gcm)
      ENDDO
      
      ! WQV
      CALL CSECT_AVG_2D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%WQV,AVG)
                        
      DO num_gcm = 1,nVGCM_seg
       channel%seg(num_seg)%WQV_E0(num_gcm) = &
         channel%seg(num_seg)%WQV_E0(num_gcm) + AVG(num_gcm)
      ENDDO
      
      ! UW
      CALL CSECT_AVG_2D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%UW(1:mi1,1:mj1),AVG)
                        
      DO num_gcm = 1,nVGCM_seg
       channel%seg(num_seg)%UW_E0(num_gcm) = &
         channel%seg(num_seg)%UW_E0(num_gcm) + AVG(num_gcm)
      ENDDO
      
      ! WV
      CALL CSECT_AVG_2D(channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,mi1,mj1, &
                        channel%seg(num_seg)%WV(1:mi1,1:mj1),AVG)
                        
      DO num_gcm = 1,nVGCM_seg
       channel%seg(num_seg)%WV_E0(num_gcm) = &
         channel%seg(num_seg)%WV_E0(num_gcm) + AVG(num_gcm)
      ENDDO                         

!******************************************************
      IF (MAKE_AVG) THEN
!     Make time averages over a timestep of vGCM
!******************************************************
      DO num_gcm = 1,nVGCM_seg
      
      channel%seg(num_seg)%SPREC_E0(num_gcm) = &
        channel%seg(num_seg)%SPREC_E0(num_gcm)/FLOAT(INUM)

      channel%seg(num_seg)%WTH_E0(num_gcm) =  &
        channel%seg(num_seg)%WTH_E0(num_gcm)/FLOAT(INUM)
        
      channel%seg(num_seg)%WQV_E0(num_gcm) =  &
        channel%seg(num_seg)%WQV_E0(num_gcm)/FLOAT(INUM)
        
      channel%seg(num_seg)%UW_E0(num_gcm) =   &
        channel%seg(num_seg)%UW_E0(num_gcm)/FLOAT(INUM)
        
      channel%seg(num_seg)%WV_E0(num_gcm) =   &
        channel%seg(num_seg)%WV_E0(num_gcm)/FLOAT(INUM)
      
      ENDDO
!*****************
      ENDIF
!*****************
!===============================   
      ENDDO
!=============================== 

      ENDIF  ! PHYSICS

      CONTAINS
      
       SUBROUTINE CSECT_AVG_2D(lsta,lend,mi1,mj1,A,AVG)
!      input  (A)   : 2D target variable        
!      output (AVG) : Average over a vGCM-cell  

       IMPLICIT NONE 
                    
       integer, INTENT(IN) :: lsta(nVGCM_seg)   ! starting CRM point of a vGCM cell
       integer, INTENT(IN) :: lend(nVGCM_seg)   ! ending CRM point of a vGCM cell 
       integer, INTENT(IN) :: mi1               ! x-size of a channel-segment
       integer, INTENT(IN) :: mj1               ! y-size of a channel-segment
       
       real (kind=dbl_kind), INTENT(IN)  :: A(mi1,mj1)
       real (kind=dbl_kind), INTENT(OUT) :: AVG(nVGCM_seg)
       
       ! Local
       integer n_gcm,chp
       real (kind=dbl_kind) :: SUM
       
!--------------------------------------       
       if (mi1.gt.mj1) then  ! x-array
!--------------------------------------       
       DO n_gcm = 1,nVGCM_seg 
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + A(CHP,1)
         ENDDO       
         AVG(n_gcm) = SUM/FLOAT(netsz)
       ENDDO     
!--------------------------------------       
       else                  ! y-array
!--------------------------------------   
       DO n_gcm = 1,nVGCM_seg 
         SUM = D0_0 
         DO CHP = ista(n_gcm),iend(n_gcm)
          SUM = SUM + A(1,CHP)
         ENDDO       
         AVG(n_gcm) = SUM/FLOAT(netsz)
       ENDDO  
!--------------------------------------         
       endif
!--------------------------------------               
       END SUBROUTINE csect_avg_2d 
       
#endif

   END SUBROUTINE cal_sfc_eft

END MODULE q3d_effect_module
