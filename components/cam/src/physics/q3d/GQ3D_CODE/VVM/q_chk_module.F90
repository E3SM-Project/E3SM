MODULE q_chk_module
! 1. Filling holes of qc, qi, qr, qs, qg (vertical borrowing only)
! 2. Call saturation adjustment
!
! Add "LocalP" option: using local pressure to calculate Temp. 

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk2,nk3,nhalo,nhalo_adv
USE constld, only: dt,rho,fnt,pbar,pibar,nomap,localp, &
                   rad2deg  ! DEBUG

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
      REAL (KIND=r8) :: PBAR8(nk2),PIBAR8(nk2),INVFNT, &
                        TSTAR,QVSTAR,QCSTAR,QISTAR,QFREEZE
      REAL (KIND=r8), DIMENSION(nk2) :: TEMP_TH,TEMP_QV,TEMP_QC,TEMP_QI
      
      LOGICAL :: FILL_CHK = .FALSE.   ! DEBUG
      
      INTEGER, DIMENSION(4) :: QCNC1,QINC1,QRNC1,QSNC1,QGNC1   ! check purpose
      INTEGER, DIMENSION(4) :: QCNC2,QINC2,QRNC2,QSNC2,QGNC2   ! check purpose
      INTEGER :: QCcount1,QIcount1,QRcount1,QScount1,QGcount1  ! check purpose
      INTEGER :: QCcount2,QIcount2,QRcount2,QScount2,QGcount2  ! check purpose      
      
      INTEGER :: I,J,K,KLOW
      INTEGER :: num_seg,mi1,mj1,mim,mip,mjm,mjp
      
      CHARACTER (LEN=50) :: title,qctxt,qitxt,qrtxt,qstxt,qgtxt
                                 
      ! Used for saturation adjustment
      IF(.NOT.LocalP) THEN
       DO K = 2, NK2
        PBAR8(K)  = PBAR(K) * 0.01_r8
        PIBAR8(K) = PIBAR(K)
       ENDDO
      ENDIF
            
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
      
      IF (FILL_CHK) THEN
        write(unit=title,fmt='(a8,i2)') 'num_seg=',num_seg
        
        write(unit=qctxt,fmt='(a13,2i5)') 'qc: N1 & N2 =',qcnc1(num_seg),qcnc2(num_seg)
        write(unit=qitxt,fmt='(a13,2i5)') 'qi: N1 & N2 =',qinc1(num_seg),qinc2(num_seg)
        write(unit=qrtxt,fmt='(a13,2i5)') 'qr: N1 & N2 =',qrnc1(num_seg),qrnc2(num_seg)
        write(unit=qstxt,fmt='(a13,2i5)') 'qs: N1 & N2 =',qsnc1(num_seg),qsnc2(num_seg)
        write(unit=qgtxt,fmt='(a13,2i5)') 'qg: N1 & N2 =',qgnc1(num_seg),qgnc2(num_seg)  
        
        CALL CHK_WRITE (title,channel)  
        CALL CHK_WRITE (qctxt,channel)
        CALL CHK_WRITE (qitxt,channel)
        CALL CHK_WRITE (qrtxt,channel)
        CALL CHK_WRITE (qstxt,channel)
        CALL CHK_WRITE (qgtxt,channel)    
      ENDIF      

!---------------------------------------
!     Saturation adjustment
!---------------------------------------
      
      DO J = 1, mj1
      DO I = 1, mi1
      KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
      
      IF (LocalP) THEN
       DO K = 2, NK2
        PBAR8(K)  = channel%seg(num_seg)%Pmid_bg(I,J,K) * 0.01_r8
        PIBAR8(K) = channel%seg(num_seg)%PImid_bg(I,J,K)
       ENDDO
      ENDIF

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
      REAL(KIND=r8), INTENT(INOUT) :: A(mim:mip,mjm:mjp,ksize)
      INTEGER, INTENT(OUT) :: NCOUNT1,NCOUNT2

      ! Local variables
      INTEGER :: I,J,K
      REAL(KIND=r8) :: DDX(mi1,mj1,ksize)
      REAL(KIND=r8) :: FRAC,PL5,PU5,SUMA,TEM

!     INITIALIZE BORROW (DDX) TO ZERO.
      DDX = 0.0_r8

!-------------------------------------------
!     CALCULATE APPORTIONMENT OF BORROWING
!-------------------------------------------
      NCOUNT1 = 0

      DO K = 2, ksize-1
       DO J = 1, mj1
        DO I = 1, mi1

        IF ( A(I,J,K) .GE. 0.0_r8 ) CYCLE

!       BORROW LOCALLY ONLY FROM POSITIVE POINTS.
        PL5 = MAX( 0.0_r8, A(I,J,K-1) )
        PU5 = MAX( 0.0_r8, A(I,J,K+1) )

        IF ( K == 2 ) PL5 = 0.0_r8
        IF ( K == ksize-1 ) PU5 = 0.0_r8

        SUMA = PL5 + PU5

!       IS THERE ENOUGH LOCALLY TO FILL HOLE?
        TEM = SUMA + A(I,J,K)

        IF ( TEM .le. 0.0_r8) THEN
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
        A(I,J,K) = 0.0_r8
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

         IF ( A(I,J,K) .GE. 0.0_r8 ) CYCLE
         A(I,J,K) = 0.0_r8
         NCOUNT2 = NCOUNT2 + 1

        ENDDO
       ENDDO
      ENDDO
 
   END SUBROUTINE fill_vert
   
!=======================================================================      
      SUBROUTINE CHK_WRITE (title,channel, &
                            TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D)
!=======================================================================
!     NUM_HALO can be different among the variables,
!     but, the horizontal array size of the variables is fixed, (mim:mip,mjm:mjp)
     
      CHARACTER (LEN=*), INTENT(IN)  :: title
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D
              
      ! Local        
      INTEGER :: num_seg,mi1,mj1,ifortw 
      INTEGER :: ilo1,ilo2,ilo3,klo
      
      ifortw = channel%num_chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65

      write(ifortw,*) TRIM(title)
       
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      klo  = 2
      ilo1 = 1 
      ilo2 = 187
      ilo3 = 375

!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'TH', &
                             channel%seg(num_seg)%TH3D(ilo1,1,klo), &
                             channel%seg(num_seg)%TH3D(ilo2,1,klo), &
                             channel%seg(num_seg)%TH3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'TH', &
                             channel%seg(num_seg)%TH3D(1,ilo1,klo), &
                             channel%seg(num_seg)%TH3D(1,ilo2,klo), &
                             channel%seg(num_seg)%TH3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QV', &
                             channel%seg(num_seg)%QV3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QV3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QV3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QV', &
                             channel%seg(num_seg)%QV3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QV3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QV3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QC', &
                             channel%seg(num_seg)%QC3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QC3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QC3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QC', &
                             channel%seg(num_seg)%QC3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QC3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QC3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QI', &
                             channel%seg(num_seg)%QI3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QI3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QI3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QI', &
                             channel%seg(num_seg)%QI3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QI3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QI3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QR', &
                             channel%seg(num_seg)%QR3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QR3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QR3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QR', &
                             channel%seg(num_seg)%QR3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QR3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QR3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QS', &
                             channel%seg(num_seg)%QS3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QS3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QS3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QS', &
                             channel%seg(num_seg)%QS3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QS3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QS3D(1,ilo3,klo)
          
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
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QG', &
                             channel%seg(num_seg)%QG3D(ilo1,1,klo), &
                             channel%seg(num_seg)%QG3D(ilo2,1,klo), &
                             channel%seg(num_seg)%QG3D(ilo3,1,klo)
        ELSE
          ! y-array
          write(ifortw,'(i2,2x,a2,2x,3E15.3)') num_seg,'QG', &
                             channel%seg(num_seg)%QG3D(1,ilo1,klo), &
                             channel%seg(num_seg)%QG3D(1,ilo2,klo), &
                             channel%seg(num_seg)%QG3D(1,ilo3,klo)
          
        ENDIF                      
        endif                        
      ENDIF                              
      


!**********************************************************************   
      ENDDO   
!**********************************************************************

      END SUBROUTINE chk_write    
      
!=======================================================================      
      SUBROUTINE CHK_title (title,chn)
!=======================================================================
      CHARACTER (LEN=*), INTENT(IN)  :: title
      INTEGER, INTENT(IN) :: chn
              
      ! Local        
      INTEGER :: ifortw 
      
      ifortw = chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65
      
      write(ifortw,*) TRIM(title),',I,J,K,num_seg,lon_deg,lat_deg,var'
        
      END SUBROUTINE chk_title         
      
!=======================================================================      
      SUBROUTINE CHK_VALUE (chn,ival,jval,kval,nval,lon,lat,var,var1,var2)
!=======================================================================
      INTEGER, INTENT(IN) :: chn,ival,jval,kval,nval
      REAL (KIND=r8), INTENT(IN) :: lon,lat,var
      REAL (KIND=r8), OPTIONAL, INTENT(IN) :: var1,var2
              
      ! Local        
      INTEGER :: ifortw 
      REAL (KIND=r8) :: lon_deg,lat_deg
      
      ifortw = chn
      if (ifortw == 5) ifortw=55
      if (ifortw == 6) ifortw=65
      
      lon_deg = lon*rad2deg
      lat_deg = lat*rad2deg

      IF (present(var1)) THEN
        if (present(var2)) then
        write(ifortw,'(4i5,2f8.1,2x,3E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var,var1,var2
        else
        write(ifortw,'(4i5,2f8.1,2x,2E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var,var1
        endif
      ELSE 
        write(ifortw,'(4i5,2f8.1,2x,E15.5)') ival,jval,kval,nval,lon_deg,lat_deg,var
      ENDIF
        
      END SUBROUTINE chk_value                  
   
END MODULE q_chk_module
