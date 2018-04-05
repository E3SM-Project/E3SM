MODULE damping
! Rayleigh-type Gravity Wave Damping in the upper layers 

      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t 
      
      USE parmsld, only: nk1,nk2
      USE constld, only: crad,klow_gwd,dt,zz,zt,thbar,qvbar,physics

IMPLICIT NONE
PRIVATE
   
PUBLIC :: damping_vort,damping_therm

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE damping_vort   : damping on vorticity
! SUBROUTINE damping_therm  : damping on theta, qv, qc, and qi
!-----------------------------------------------------------------------

!=======================================================================
   SUBROUTINE DAMPING_VORT (channel)
!=======================================================================   
!  Gravity Wave Damping (Rayleigh-type) on vorticity components 

      type(channel_t), INTENT(INOUT) :: channel   ! Channel data

!     Local variables 
      INTEGER :: I, J, K, mi1, mj1, num_seg
      REAL (KIND=r8) :: CGR

!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
!     Gravity Wave Damping: if z > GWDHT     

!     horizontal vorticity change due to Gravity Wave Damping    
      DO K = klow_gwd, nk1      
        CGR = CRAD*(ZZ(K)-ZT(klow_gwd))/(ZT(NK2)-ZT(klow_gwd))       
        DO J = 1, mj1
         DO I = 1, mi1
           channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K) &
                                            - DT*CGR*channel%seg(num_seg)%Z3DX(I,J,K)
           channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K) &
                                            - DT*CGR*channel%seg(num_seg)%Z3DY(I,J,K)
         ENDDO
        ENDDO
      ENDDO
      
!     vertical vorticity change due to Gravity Wave Damping (at the uppermost layer)
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                                           - DT*CRAD*channel%seg(num_seg)%Z3DZ(I,J,NK2)
       ENDDO
      ENDDO

!     Horizontal Momentum tendency due to Gravity Wave Damping  
      DO K = klow_gwd, nk2
        CGR = CRAD*(ZT(K)-ZT(klow_gwd))/(ZT(NK2)-ZT(klow_gwd))
        DO J = 1, mj1
          DO I = 1, mi1
           channel%seg(num_seg)%FU_DIA(I,J,K) = channel%seg(num_seg)%FU_DIA(I,J,K) &
                                              - CGR*channel%seg(num_seg)%U3DX(I,J,K)
           channel%seg(num_seg)%FV_DIA(I,J,K) = channel%seg(num_seg)%FV_DIA(I,J,K) &
                                              - CGR*channel%seg(num_seg)%U3DY(I,J,K)
          ENDDO
        ENDDO
      ENDDO

!*************************   
      ENDDO  ! num_seg 
!*************************
      
   END SUBROUTINE damping_vort

!=======================================================================
   SUBROUTINE DAMPING_THERM (channel)
!=======================================================================   
!  (Rayleigh-type) Gravity Wave Damping on (th, qv, qc, qi)
!-----------------------------------------------------------------------  
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
!     Local variables   
      INTEGER :: num_seg,mi1,mj1,I,J,K
      REAL (KIND=r8) :: CGR,TEMP_TH,TEMP_QV,TEMP_QC,TEMP_QI

!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
                 
!     Gravity Wave Damping: if z > GWDHT   
      DO K = klow_gwd, nk2     
         CGR =  CRAD*(ZT(K)-ZT(klow_gwd))/(ZT(NK2)-ZT(klow_gwd))
         
         DO J = 1, mj1
          DO I = 1, mi1
           TEMP_TH = - CGR*(channel%seg(num_seg)%TH3D(I,J,K)-THBAR(K))
           TEMP_QV = - CGR*(channel%seg(num_seg)%QV3D(I,J,K)-QVBAR(K))

!          Update the fields
           channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                                            + DT*TEMP_TH
           channel%seg(num_seg)%QV3D(I,J,K) = channel%seg(num_seg)%QV3D(I,J,K) &
                                            + DT*TEMP_QV
                                            
!          Update the diabatic tendency: damping effects                                            
           channel%seg(num_seg)%FTH3D_DIA(I,J,K) = channel%seg(num_seg)%FTH3D_DIA(I,J,K) &
                                                 + TEMP_TH
           channel%seg(num_seg)%FQV3D_DIA(I,J,K) = channel%seg(num_seg)%FQV3D_DIA(I,J,K) &
                                                 + TEMP_QV                                              
          ENDDO
         ENDDO 
          
         IF (PHYSICS) THEN    
         DO J = 1, mj1
          DO I = 1, mi1    
           TEMP_QC = - CGR*channel%seg(num_seg)%QC3D(I,J,K)
           TEMP_QI = - CGR*channel%seg(num_seg)%QI3D(I,J,K)

!          Update the fields
           channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K) &
                                            + DT*TEMP_QC
           channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K) &
                                            + DT*TEMP_QI
                                            
!          Update the diabatic tendency: damping effects                                                
           channel%seg(num_seg)%FQC3D_DIA(I,J,K) = channel%seg(num_seg)%FQC3D_DIA(I,J,K) &
                                                 + TEMP_QC
           channel%seg(num_seg)%FQI3D_DIA(I,J,K) = channel%seg(num_seg)%FQI3D_DIA(I,J,K) &
                                                 + TEMP_QI                                           
          ENDDO
         ENDDO
         ENDIF     
      
      ENDDO  ! k-loop

!*************************   
      ENDDO  ! num_seg 
!*************************

   END SUBROUTINE damping_therm

END MODULE damping
