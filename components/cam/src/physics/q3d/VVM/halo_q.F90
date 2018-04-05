MODULE halo_q
! Map the halo values from neighboring cube surface (data at the q-point). 

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t 

USE parmsld,        only: ntracer,nk1,nk2,nhalo
  
! Functions being called
USE utils, only: fintrp,cubicint

IMPLICIT NONE
PRIVATE

PUBLIC :: halo_correc_q    ! correct the halo values from other cube faces 
          
CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------
! SUBROUTINE halo_correc_q
!
! SUBROUTINE HALO_COR_ALPHA : Map the halo-values from other faces 
!                      through 1-D interpolation (ALONG THE ALPHA AXIS)
! SUBROUTINE HALO_COR_BETA  : Map the halo-values from other faces 
!                      through 1-D interpolation (ALONG THE BETA AXIS)
!----------------------------------------------------------------------

!=======================================================================      
      SUBROUTINE HALO_CORREC_Q (num_halo,channel, &
                                TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D,QT3D, &
                                W3D,CHI,TOPOZ,ZROUGH,GWET,TG)
!=======================================================================
!     NUM_HALO can be different among the variables,
!     but, the horizontal array size of the variables is fixed, (mim:mip,mjm:mjp)

      INTEGER, INTENT(IN) :: num_halo             ! Number of halo points to be mapped 
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D,QT3D, &
                                       W3D,CHI,TOPOZ,ZROUGH,GWET,TG
              
      ! Local        
      INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp,nt,nlenI,nlenJ 
      INTEGER :: MODE1= 1, MODE2= 1
      ! 1 (Linear interpolation) 2 (Cubic Spline interpolation)
      ! MODE2 is used for positive definite variables
      
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip   
      
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp    

      nlenJ = mj1 + 2*nhalo
      nlenI = mi1 + 2*nhalo      
!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE1,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%th3d(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE1,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%th3d(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QV3D
!---------------------------------
      IF (PRESENT(QV3D)) THEN
        if (QV3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QV3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QV3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QC3D
!---------------------------------
      IF (PRESENT(QC3D)) THEN
        if (QC3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QC3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QC3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QI3D
!---------------------------------
      IF (PRESENT(QI3D)) THEN
        if (QI3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QI3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QI3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QR3D
!---------------------------------
      IF (PRESENT(QR3D)) THEN
        if (QR3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QR3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QR3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QS3D
!---------------------------------
      IF (PRESENT(QS3D)) THEN
        if (QS3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QS3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QS3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QG3D
!---------------------------------
      IF (PRESENT(QG3D)) THEN
        if (QG3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QG3D(:,:,2:nk2))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QG3D(:,:,2:nk2))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! QT3D
!---------------------------------
      IF (PRESENT(QT3D)) THEN
        if (QT3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array     
          DO nt = 1,ntracer   
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%QT3D(:,:,2:nk2,nt))   
          ENDDO                           
        ELSE
          ! y-array
          DO nt = 1,ntracer
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1,num_halo,MODE2,    &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%QT3D(:,:,2:nk2,nt))
          ENDDO 
        ENDIF                      
        endif                        
      ENDIF                                          
!---------------------------------
! W3D
!---------------------------------
      IF (PRESENT(W3D)) THEN
        if (W3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlenJ,nk1-1,num_halo,MODE1,  &
                              channel%seg(num_seg)%jhalo_loc_t,                &
                              channel%seg(num_seg)%cbeta_t,                    &
                              channel%seg(num_seg)%hbeta_t,                    &
                              channel%seg(num_seg)%W3D(:,:,2:nk1))         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlenI,nk1-1,num_halo,MODE1,  &
                              channel%seg(num_seg)%ihalo_loc_t,                &
                              channel%seg(num_seg)%calpha_t,                   &
                              channel%seg(num_seg)%halpha_t,                   &
                              channel%seg(num_seg)%W3D(:,:,2:nk1))
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! CHI
!---------------------------------
      IF (PRESENT(CHI)) THEN
        if (CHI) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlenJ,num_halo,MODE1,      &
                                channel%seg(num_seg)%jhalo_loc_t,              &
                                channel%seg(num_seg)%cbeta_t,                  &
                                channel%seg(num_seg)%hbeta_t,                  &
                                channel%seg(num_seg)%CHI)         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlenI,num_halo,MODE1,      &
                                channel%seg(num_seg)%ihalo_loc_t,              &
                                channel%seg(num_seg)%calpha_t,                 &
                                channel%seg(num_seg)%halpha_t,                 &
                                channel%seg(num_seg)%CHI)
        ENDIF                      
        endif                        
      ENDIF 
!---------------------------------
! TG
!---------------------------------
      IF (PRESENT(TG)) THEN
        if (TG) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlenJ,num_halo,MODE1,      &
                                channel%seg(num_seg)%jhalo_loc_t,              &
                                channel%seg(num_seg)%cbeta_t,                  &
                                channel%seg(num_seg)%hbeta_t,                  &
                                channel%seg(num_seg)%TG)         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlenI,num_halo,MODE1,      &
                                channel%seg(num_seg)%ihalo_loc_t,              &
                                channel%seg(num_seg)%calpha_t,                 &
                                channel%seg(num_seg)%halpha_t,                 &
                                channel%seg(num_seg)%TG)
        ENDIF                      
        endif                        
      ENDIF 
      
! Probably, the below will not be used.      
!---------------------------------
! TOPOZ
!---------------------------------
      IF (PRESENT(TOPOZ)) THEN
        if (TOPOZ) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlenJ,num_halo,MODE2,      &
                                channel%seg(num_seg)%jhalo_loc_t,              &
                                channel%seg(num_seg)%cbeta_t,                  &
                                channel%seg(num_seg)%hbeta_t,                  &
                                channel%seg(num_seg)%TOPOZ)         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlenI,num_halo,MODE2,      &
                                channel%seg(num_seg)%ihalo_loc_t,              &
                                channel%seg(num_seg)%calpha_t,                 &
                                channel%seg(num_seg)%halpha_t,                 &
                                channel%seg(num_seg)%TOPOZ)
        ENDIF                      
        endif                        
      ENDIF 
!---------------------------------
! ZROUGH
!---------------------------------
      IF (PRESENT(ZROUGH)) THEN
        if (ZROUGH) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlenJ,num_halo,MODE2,      &
                                channel%seg(num_seg)%jhalo_loc_t,              &
                                channel%seg(num_seg)%cbeta_t,                  &
                                channel%seg(num_seg)%hbeta_t,                  &
                                channel%seg(num_seg)%ZROUGH)         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlenI,num_halo,MODE2,      &
                                channel%seg(num_seg)%ihalo_loc_t,              &
                                channel%seg(num_seg)%calpha_t,                 &
                                channel%seg(num_seg)%halpha_t,                 &
                                channel%seg(num_seg)%ZROUGH)
        ENDIF                      
        endif                        
      ENDIF 
!---------------------------------
! GWET
!---------------------------------
      IF (PRESENT(GWET)) THEN
        if (GWET) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlenJ,num_halo,MODE2,      &
                                channel%seg(num_seg)%jhalo_loc_t,              &
                                channel%seg(num_seg)%cbeta_t,                  &
                                channel%seg(num_seg)%hbeta_t,                  &
                                channel%seg(num_seg)%GWET)         
        ELSE
          ! y-array
          CALL HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlenI,num_halo,MODE2,      &
                                channel%seg(num_seg)%ihalo_loc_t,              &
                                channel%seg(num_seg)%calpha_t,                 &
                                channel%seg(num_seg)%halpha_t,                 &
                                channel%seg(num_seg)%GWET)
        ENDIF                      
        endif                        
      ENDIF                         
           
!********************************   
       ENDDO  ! num_seg 
!******************************** 

      END SUBROUTINE halo_correc_q

!=======================================================================  
      SUBROUTINE HALO_COR_ALPH (mj1,mjm,mjp,mim,mip,nlen,KDIMN,mhalo,MODE, &
                                IHALO_LOC,CALPHA,HALPHA,A)
!======================================================================= 
!     (ALONG THE ALPHA AXIS) Mapping halo-values through 1-D interpolation. 

      INTEGER, INTENT(IN) :: mj1,mjm,mjp,mim,mip,nlen  ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN        ! vertical array size  
      INTEGER, INTENT(IN) :: mhalo        ! number of halo points for interpolation (.le. nhalo) 
      INTEGER, INTENT(IN) :: MODE         ! interpolation type (1: linear  2: cubic spline)        
      
      ! i-index of lower points of given values for linear interpolation      
      INTEGER,DIMENSION(mim:mip,nhalo),INTENT(IN) :: IHALO_LOC 
      
      ! alpha values of the regular q-grid points   
      REAL (KIND=r8),DIMENSION(mim:mip),INTENT(IN) :: CALPHA
        
      ! alpha values of the target points 
      REAL (KIND=r8),DIMENSION(mim:mip,nhalo),INTENT(IN) :: HALPHA
      
      ! target variable
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(NLEN) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=r8),DIMENSION(mim:mip,KDIMN) :: HMV,HPV  
      INTEGER :: I,J,K,I1VAL,I2VAL
      
      IF (MODE.EQ.2)  THEN
        DO I = mim,mip  
         XVAL(I+nhalo) = CALPHA(I)
        ENDDO 
      ENDIF
      
!=======================================       
      SELECT CASE (MODE)
!=======================================  
!================================= 
!     Linear interpolation
      CASE(1)
!=================================    
      DO J=1,mhalo   
        
         DO I=mim,mip
          I1VAL = IHALO_LOC(I,J)
          I2VAL = I1VAL + 1
       
          DO K=1,KDIMN 
           HMV(I,K)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,1-J,K)   &
                                          ,CALPHA(I2VAL),A(I2VAL,1-J,K))
           HPV(I,K)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,MJ1+J,K) &
                                          ,CALPHA(I2VAL),A(I2VAL,MJ1+J,K)) 
          ENDDO                        
         ENDDO                                                                     

        DO K=1,KDIMN 
         DO I=mim,mip
          A(I,1-J,K)   = HMV(I,K)
          A(I,mj1+J,K) = HPV(I,K)
         ENDDO
        ENDDO  
       
      ENDDO  ! J-loop
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!================================= 
      DO J=1,mhalo
      
         DO I=mim,mip
          XNEW(I+nhalo) = HALPHA(I,J)
         ENDDO   
         
         DO K=1,KDIMN
         
          DO I=mim,mip
           YVAL(I+nhalo) = A(I,1-J,K)
          ENDDO                                                                           
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO I=mim,mip
           HMV(I,K) = YNEW(I+nhalo) 
          ENDDO 
          DO I=mim,mip
           YVAL(I+nhalo) = A(I,MJ1+J,K)
          ENDDO                                                                           
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO I=mim,mip
           HPV(I,K) = YNEW(I+nhalo) 
          ENDDO 
                    
         ENDDO  ! k-loop

         DO K=1,KDIMN 
          DO I=mim,mip
           A(I,1-J,K)   = HMV(I,K)
           A(I,mj1+J,K) = HPV(I,K)
          ENDDO
         ENDDO  
          
      ENDDO   ! J-loop 
!================================= 
      END SELECT
!=================================

      END SUBROUTINE halo_cor_alph     

!=======================================================================  
      SUBROUTINE HALO_COR_ALPH_2 (mj1,mjm,mjp,mim,mip,nlen,mhalo,MODE, &
                                  IHALO_LOC,CALPHA,HALPHA,A)
!======================================================================= 
!     Same as HALO_COR_ALPH, but used for a 2D variable.  

      INTEGER, INTENT(IN) :: mj1,mjm,mjp,mim,mip,nlen  ! horizontal array sizes 
      INTEGER, INTENT(IN) :: mhalo        ! number of halo points for interpolation (.le. nhalo) 
      INTEGER, INTENT(IN) :: MODE         ! interpolation type (1: linear  2: cubic spline)        
      
      ! i-index of lower points of given values for linear interpolation      
      INTEGER,DIMENSION(mim:mip,nhalo),INTENT(IN) :: IHALO_LOC 
      
      ! alpha values of the regular q-grid points   
      REAL (KIND=r8),DIMENSION(mim:mip),INTENT(IN) :: CALPHA
        
      ! alpha values of the target points 
      REAL (KIND=r8),DIMENSION(mim:mip,nhalo),INTENT(IN) :: HALPHA
      
      ! target variable
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(NLEN) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=r8),DIMENSION(mim:mip) :: HMV,HPV  
      INTEGER :: I,J,I1VAL,I2VAL
      
      IF (MODE.EQ.2)  THEN
        DO I = mim,mip  
         XVAL(I+nhalo) = CALPHA(I)
        ENDDO 
      ENDIF
      
!=======================================       
      SELECT CASE (MODE)
!=======================================  
!================================= 
!     Linear interpolation
      CASE(1)
!=================================    
      DO J=1,mhalo   
        
         DO I=mim,mip
          I1VAL = IHALO_LOC(I,J)
          I2VAL = I1VAL + 1
       
           HMV(I)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,1-J)   &
                                        ,CALPHA(I2VAL),A(I2VAL,1-J))
           HPV(I)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,MJ1+J) &
                                        ,CALPHA(I2VAL),A(I2VAL,MJ1+J))                        
         ENDDO                                                                     

         DO I=mim,mip
          A(I,1-J)   = HMV(I)
          A(I,mj1+J) = HPV(I)
         ENDDO
       
      ENDDO  ! J-loop
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!================================= 
      DO J=1,mhalo
      
         DO I=mim,mip
          XNEW(I+nhalo) = HALPHA(I,J)
         ENDDO   
         
          DO I=mim,mip
           YVAL(I+nhalo) = A(I,1-J)
          ENDDO                                                                           
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO I=mim,mip
           HMV(I) = YNEW(I+nhalo) 
          ENDDO 
          DO I=mim,mip
           YVAL(I+nhalo) = A(I,MJ1+J)
          ENDDO                                                                           
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO I=mim,mip
           HPV(I) = YNEW(I+nhalo) 
          ENDDO 

          DO I=mim,mip
           A(I,1-J)   = HMV(I)
           A(I,mj1+J) = HPV(I)
          ENDDO 
          
      ENDDO   ! J-loop 
!================================= 
      END SELECT
!=================================
      
      END SUBROUTINE halo_cor_alph_2           

!=======================================================================          
      SUBROUTINE HALO_COR_BETA (mi1,mim,mip,mjm,mjp,nlen,KDIMN,mhalo,MODE, &
                                JHALO_LOC,CBETA,HBETA,A)
!=======================================================================  
!     (ALONG THE BETA AXIS) Mapping halo-values through 1-D interpolation. 
    
      INTEGER, INTENT(IN) :: mi1,mim,mip,mjm,mjp,nlen  ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN        ! vertical array size  
      INTEGER, INTENT(IN) :: mhalo        ! number of halo points for interpolation (.le. nhalo) 
      INTEGER, INTENT(IN) :: MODE         ! interpolation type (1: linear  2: cubic spline)        
      
      ! j-index of lower points of given values for linear interpolation      
      INTEGER,DIMENSION(nhalo,mjm:mjp),INTENT(IN) :: JHALO_LOC 
      
      ! beta values of the regular q-grid points       
      REAL (KIND=r8),DIMENSION(mjm:mjp),INTENT(IN) :: CBETA
      
      ! beta values of the target points 
      REAL (KIND=r8),DIMENSION(nhalo,mjm:mjp),INTENT(IN) :: HBETA
      
      ! target variable
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(NLEN) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=r8),DIMENSION(mjm:mjp,KDIMN) :: HMV,HPV    
      INTEGER :: I,J,K,J1VAL,J2VAL
      
      IF (MODE.EQ.2)  THEN
        DO J = mjm,mjp  
         XVAL(J+nhalo) = CBETA(J)
        ENDDO 
      ENDIF

!=======================================       
      SELECT CASE (MODE)
!======================================= 
       
!================================= 
!     Linear interpolation
      CASE(1)
!================================= 
      DO I=1,mhalo
 
         DO J=mjm,mjp
          J1VAL = JHALO_LOC(I,J)
          J2VAL = J1VAL+1
       
          DO K=1,KDIMN 
           HMV(J,K)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(1-I,J1VAL,K)   &
                                         ,CBETA(J2VAL),A(1-I,J2VAL,K))
           HPV(J,K)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(mi1+I,J1VAL,K) &
                                         ,CBETA(J2VAL),A(mi1+I,J2VAL,K)) 
          ENDDO                           
         ENDDO      

         DO K=1,KDIMN 
          DO J=mjm,mjp
           A(1-I,J,K)   = HMV(J,K)
           A(mi1+I,J,K) = HPV(J,K)
          ENDDO
         ENDDO  
       
       ENDDO  ! I-loop
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================   
      DO I=1,mhalo

         DO J = mjm,mjp  
          XNEW(J+nhalo) = HBETA(I,J)
         ENDDO        
       
         DO K=1,KDIMN
          DO J=mjm,mjp
           YVAL(J+nhalo) = A(1-I,J,K)
          ENDDO                                                                  
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO J=mjm,mjp
           HMV(J,K) = YNEW(J+nhalo) 
          ENDDO       
          DO J=mjm,mjp
           YVAL(J+nhalo) = A(MI1+I,J,K)
          ENDDO                                                                  
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO J=mjm,mjp
           HPV(J,K) = YNEW(J+nhalo) 
          ENDDO                       
         ENDDO  ! k-loop
 
         DO K=1,KDIMN 
          DO J=mjm,mjp
           A(1-I,J,K)   = HMV(J,K)
           A(mi1+I,J,K) = HPV(J,K)
          ENDDO
         ENDDO  
        
      ENDDO  ! I-loop
!================================= 
      END SELECT
!=================================  
      
      END SUBROUTINE halo_cor_beta
      
!=======================================================================          
      SUBROUTINE HALO_COR_BETA_2 (mi1,mim,mip,mjm,mjp,nlen,mhalo,MODE, &
                                  JHALO_LOC,CBETA,HBETA,A)
!=======================================================================  
!     Same as HALO_COR_BETA, but used for a 2D variable.
    
      INTEGER, INTENT(IN) :: mi1,mim,mip,mjm,mjp,nlen  ! horizontal array sizes
      INTEGER, INTENT(IN) :: mhalo        ! number of halo points for interpolation (.le. nhalo) 
      INTEGER, INTENT(IN) :: MODE         ! interpolation type (1: linear  2: cubic spline)        
      
      ! j-index of lower points of given values for linear interpolation      
      INTEGER,DIMENSION(nhalo,mjm:mjp),INTENT(IN) :: JHALO_LOC 
      
      ! beta values of the regular q-grid points       
      REAL (KIND=r8),DIMENSION(mjm:mjp),INTENT(IN) :: CBETA
      
      ! beta values of the target points 
      REAL (KIND=r8),DIMENSION(nhalo,mjm:mjp),INTENT(IN) :: HBETA
      
      ! target variable
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(NLEN) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=r8),DIMENSION(mjm:mjp) :: HMV,HPV    
      INTEGER :: I,J,J1VAL,J2VAL
      
      IF (MODE.EQ.2)  THEN
        DO J = mjm,mjp  
         XVAL(J+nhalo) = CBETA(J)
        ENDDO 
      ENDIF

!=======================================       
      SELECT CASE (MODE)
!======================================= 
       
!================================= 
!     Linear interpolation
      CASE(1)
!================================= 
      DO I=1,mhalo
 
         DO J=mjm,mjp
          J1VAL = JHALO_LOC(I,J)
          J2VAL = J1VAL+1
       
           HMV(J)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(1-I,J1VAL)   &
                                       ,CBETA(J2VAL),A(1-I,J2VAL))
           HPV(J)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(mi1+I,J1VAL) &
                                       ,CBETA(J2VAL),A(mi1+I,J2VAL))                           
         ENDDO      

         DO J=mjm,mjp
           A(1-I,J)   = HMV(J)
           A(mi1+I,J) = HPV(J)
         ENDDO  
       
       ENDDO  ! I-loop
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================   
      DO I=1,mhalo

         DO J = mjm,mjp  
          XNEW(J+nhalo) = HBETA(I,J)
         ENDDO        
       
          DO J=mjm,mjp
           YVAL(J+nhalo) = A(1-I,J)
          ENDDO                                                                  
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO J=mjm,mjp
           HMV(J) = YNEW(J+nhalo) 
          ENDDO       
          DO J=mjm,mjp
           YVAL(J+nhalo) = A(MI1+I,J)
          ENDDO                                                                  
          CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
          DO J=mjm,mjp
           HPV(J) = YNEW(J+nhalo) 
          ENDDO                       
         
          DO J=mjm,mjp
           A(1-I,J)   = HMV(J)
           A(mi1+I,J) = HPV(J)
          ENDDO  
        
      ENDDO  ! I-loop
!================================= 
      END SELECT
!=================================  
      
      END SUBROUTINE halo_cor_beta_2      
      
END MODULE halo_q