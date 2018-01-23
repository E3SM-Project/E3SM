MODULE q_chk_module
! 1. Filling holes of qc, qi, qr, qs, qg (vertical borrowing only)
! 2. Call saturation adjustment

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk2,nk3,nhalo_adv
USE constld, only: d0_0,dt,rho,fnt,pbar,pibar,nomap

! Subroutines being called
USE bound_channel_module, only: bound_channel
USE bound_extra, only: bound_normal,bound_vert
USE halo_q, only: halo_correc_q

IMPLICIT NONE
PRIVATE

PUBLIC :: q_chk_3d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE q_chk_3d  : call FILL_VERT and call SATURATION
! SUBROUTINE FILL_VERT : filling the negative values (holes)
!-----------------------------------------------------------------------

!=======================================================================
      SUBROUTINE Q_CHK_3D (channel)
!======================================================================= 
      type(channel_t), intent(inout) :: channel   ! channel data

      ! Local variables
      REAL (KIND=dbl_kind) :: PBAR8(nk2),PIBAR8(nk2),INVFNT, &
                              TSTAR,QVSTAR,QCSTAR,QISTAR,QFREEZE
      REAL (KIND=dbl_kind), DIMENSION(nk2) :: TEMP_TH,TEMP_QV,TEMP_QC,TEMP_QI
      
      LOGICAL :: FILL_CHK = .FALSE.
      INTEGER, DIMENSION(4) :: QCNC1,QINC1,QRNC1,QSNC1,QGNC1   ! check purpose
      INTEGER, DIMENSION(4) :: QCNC2,QINC2,QRNC2,QSNC2,QGNC2   ! check purpose
      INTEGER :: QCcount1,QIcount1,QRcount1,QScount1,QGcount1  ! check purpose
      INTEGER :: QCcount2,QIcount2,QRcount2,QScount2,QGcount2  ! check purpose      
      
      INTEGER :: I,J,K,KLOW
      INTEGER :: num_seg,mi1,mj1,mim,mip,mjm,mjp

      ! Used for saturation adjustment
      DO K = 2, NK2
       PBAR8(K)  = PBAR(K) * 0.01_dbl_kind
       PIBAR8(K) = PIBAR(K)
      ENDDO
            
!*******************************************************************   
      DO num_seg = 1, 4
!*******************************************************************  
      mi1  = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim  = channel%seg(num_seg)%mim  
      mip  = channel%seg(num_seg)%mip   

      mj1  = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm  = channel%seg(num_seg)%mjm   
      mjp  = channel%seg(num_seg)%mjp 
                        
!---------------------------------
!     Filling Holes
!---------------------------------
      DO K = 2, nk2
       invfnt = rho(K)/fnt(k)
       DO J = mjm, mjp
        DO I = mim, mip
         channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K)*invfnt
         channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K)*invfnt
         channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D(I,J,K)*invfnt
         channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D(I,J,K)*invfnt
         channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D(I,J,K)*invfnt
        ENDDO
       ENDDO
      ENDDO

      CALL FILL_VERT(channel%seg(num_seg)%QC3D,mi1,mj1,mim,mip,mjm,mjp,nk3, &
                     QCNC1(num_seg),QCNC2(num_seg))
      CALL FILL_VERT(channel%seg(num_seg)%QI3D,mi1,mj1,mim,mip,mjm,mjp,nk3, &
                     QINC1(num_seg),QINC2(num_seg))
      CALL FILL_VERT(channel%seg(num_seg)%QR3D,mi1,mj1,mim,mip,mjm,mjp,nk3, &
                     QRNC1(num_seg),QRNC2(num_seg))
      CALL FILL_VERT(channel%seg(num_seg)%QS3D,mi1,mj1,mim,mip,mjm,mjp,nk3, &
                     QSNC1(num_seg),QSNC2(num_seg))
      CALL FILL_VERT(channel%seg(num_seg)%QG3D,mi1,mj1,mim,mip,mjm,mjp,nk3, &
                     QGNC1(num_seg),QGNC2(num_seg))

      DO K = 2, NK2
       invfnt = fnt(k)/rho(K)
       DO J = mjm, mjp
        DO I = mim ,mip
         channel%seg(num_seg)%QC3D(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K)*invfnt
         channel%seg(num_seg)%QI3D(I,J,K) = channel%seg(num_seg)%QI3D(I,J,K)*invfnt
         channel%seg(num_seg)%QS3D(I,J,K) = channel%seg(num_seg)%QS3D(I,J,K)*invfnt
         channel%seg(num_seg)%QG3D(I,J,K) = channel%seg(num_seg)%QG3D(I,J,K)*invfnt
         channel%seg(num_seg)%QR3D(I,J,K) = channel%seg(num_seg)%QR3D(I,J,K)*invfnt
        ENDDO
       ENDDO
      ENDDO

!---------------------------------------
!     Saturation adjustment
!---------------------------------------
      DO J = 1, mj1
      DO I = 1, mi1
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
      
      DO K = KLOW, nk2
        TEMP_TH(K) = channel%seg(num_seg)%TH3D(I,J,K)
        TEMP_QV(K) = channel%seg(num_seg)%QV3D(I,J,K)
        TEMP_QC(K) = channel%seg(num_seg)%QC3D(I,J,K)
        TEMP_QI(K) = channel%seg(num_seg)%QI3D(I,J,K)
      ENDDO
       
      DO K = KLOW, nk2
        TSTAR  = TEMP_TH(K) * PIBAR8(K)
        QVSTAR = TEMP_QV(K)
        QCSTAR = TEMP_QC(K)
        QISTAR = TEMP_QI(K)
        
        CALL SATURATION(I,J,K,TSTAR,PBAR8(K),QVSTAR,QCSTAR,QISTAR,QFREEZE)
        
        channel%seg(num_seg)%TH3D(I,J,K) = TSTAR / PIBAR8(K)
        channel%seg(num_seg)%QV3D(I,J,K) = QVSTAR
        channel%seg(num_seg)%QC3D(I,J,K) = QCSTAR
        channel%seg(num_seg)%QI3D(I,J,K) = QISTAR
      ENDDO
      DO K=2,KLOW-1
       channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,KLOW)
      ENDDO

!     Diabatic Tendency: Saturation Adjustment Effect   
      DO K = KLOW, NK2
        channel%seg(num_seg)%FTH3D_DIA(I,J,K) = &
                          channel%seg(num_seg)%FTH3D_DIA(I,J,K) &
                       + (channel%seg(num_seg)%TH3D(I,J,K)-TEMP_TH(K))/DT
        channel%seg(num_seg)%FQV3D_DIA(I,J,K) = &
                          channel%seg(num_seg)%FQV3D_DIA(I,J,K) &
                       + (channel%seg(num_seg)%QV3D(I,J,K)-TEMP_QV(K))/DT       
        channel%seg(num_seg)%FQC3D_DIA(I,J,K) = &
                          channel%seg(num_seg)%FQC3D_DIA(I,J,K) &
                       + (channel%seg(num_seg)%QC3D(I,J,K)-TEMP_QC(K))/DT       
        channel%seg(num_seg)%FQI3D_DIA(I,J,K) = &
                          channel%seg(num_seg)%FQI3D_DIA(I,J,K) &
                       + (channel%seg(num_seg)%QI3D(I,J,K)-TEMP_QI(K))/DT   
      ENDDO

      ENDDO  ! i-loop
      ENDDO  ! j-loop
!**************************   
      ENDDO  ! num_seg 
!**************************  

      IF (FILL_CHK) THEN
        QCcount1 = SUM(QCNC1)
        QIcount1 = SUM(QINC1)
        QRcount1 = SUM(QRNC1)
        QScount1 = SUM(QSNC1)
        QGcount1 = SUM(QGNC1)
        
        QCcount2 = SUM(QCNC2)
        QIcount2 = SUM(QINC2)
        QRcount2 = SUM(QRNC2)
        QScount2 = SUM(QSNC2)
        QGcount2 = SUM(QGNC2)
        WRITE(6,*) 'FILL_CHK  QC: N1 & N2 =',QCcount1,QCcount2
        WRITE(6,*) 'FILL_CHK  QI: N1 & N2 =',QIcount1,QIcount2
        WRITE(6,*) 'FILL_CHK  QR: N1 & N2 =',QRcount1,QRcount2
        WRITE(6,*) 'FILL_CHK  QS: N1 & N2 =',QScount1,QScount2
        WRITE(6,*) 'FILL_CHK  QG: N1 & N2 =',QGcount1,QGcount2
      ENDIF
      
      CALL BOUND_NORMAL  (nhalo,channel,TH3D=.TRUE.,QV3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
                          
      CALL BOUND_CHANNEL (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)
      IF (.not.nomap) &
      CALL HALO_CORREC_Q (nhalo_adv,channel,TH3D=.TRUE.,QV3D=.TRUE., &
                          QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

      CALL BOUND_VERT (channel,TH3D=.TRUE.,QV3D=.TRUE., &
                       QC3D=.TRUE.,QI3D=.TRUE.,QR3D=.TRUE.,QS3D=.TRUE.,QG3D=.TRUE.)

   END SUBROUTINE q_chk_3d

!=======================================================================
   SUBROUTINE FILL_VERT(A,mi1,mj1,mim,mip,mjm,mjp,ksize,NCOUNT1,NCOUNT2)
!=======================================================================
      INTEGER,  INTENT(IN) :: mi1,mj1,mim,mip,mjm,mjp,ksize  
      REAL(KIND=dbl_kind), INTENT(INOUT) :: A(mim:mip,mjm:mjp,ksize)
      INTEGER, INTENT(OUT) :: NCOUNT1,NCOUNT2

      ! Local variables
      INTEGER :: I,J,K
      REAL(KIND=dbl_kind) :: DDX(mi1,mj1,ksize)
      REAL(KIND=dbl_kind) :: FRAC,PL5,PU5,SUMA,TEM

!     INITIALIZE BORROW (DDX) TO ZERO.
      DDX = D0_0

!-------------------------------------------
!     CALCULATE APPORTIONMENT OF BORROWING
!-------------------------------------------
      NCOUNT1 = 0

      DO K = 2, ksize-1
       DO J = 1, mj1
        DO I = 1, mi1

        IF ( A(I,J,K) .GE. D0_0 ) CYCLE

!       BORROW LOCALLY ONLY FROM POSITIVE POINTS.
        PL5 = MAX( D0_0, A(I,J,K-1) )
        PU5 = MAX( D0_0, A(I,J,K+1) )

        IF ( K == 2 ) PL5 = D0_0
        IF ( K == ksize-1 ) PU5 = D0_0

        SUMA = PL5 + PU5

!       IS THERE ENOUGH LOCALLY TO FILL HOLE?
        TEM = SUMA + A(I,J,K)

        IF ( TEM .le. D0_0) THEN
!         NOT ENOUGH, SO TAKE EVERYTHING AVAILABLE LOCALLY
          FRAC = 1.
        ELSE
!         ENOUGH EXISTS LOCALLY TO FILL.
          FRAC = - A(I,J,K) / SUMA
        ENDIF

!       CALCULATE LOCAL BORROWING
        DDX(I,J,K-1) = DDX(I,J,K-1) + FRAC * PL5
        DDX(I,J,K+1) = DDX(I,J,K+1) + FRAC * PU5

!       FILL HOLE ( RESET TO ZERO ).
        A(I,J,K) = D0_0
        NCOUNT1 = NCOUNT1 + 1

        ENDDO
       ENDDO
      ENDDO

!     DO LOCAL BORROWING
      DO K = 2, ksize-1
       DO J = 1, mj1
        DO I = 1, mi1
         A(I,J,K) = A(I,J,K) - DDX(I,J,K)
        ENDDO
       ENDDO
      ENDDO

!     Set the rest zero
      NCOUNT2 = 0
      DO K = 2, ksize-1
       DO J = 1, mj1
        DO I = 1, mi1

         IF ( A(I,J,K) .GE. D0_0 ) CYCLE
         A(I,J,K) = D0_0
         NCOUNT2 = NCOUNT2 + 1

        ENDDO
       ENDDO
      ENDDO
 
   END SUBROUTINE fill_vert

END MODULE q_chk_module
