MODULE halo_z
! Routines that correct the halopoint values from neighboring face.
! Target: 2-D data at the z-point (on a given height) 

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nhalo,nk2

! Function
USE utils, only: fintrp,cubicint

IMPLICIT NONE
PRIVATE

PUBLIC :: halo_correc_z

CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------
! SUBROUTINE halo_correc_z
!
! SUBROUTINE HALO_X : data arrangement for x-channel segments
! SUBROUTINE HALO_Y : data arrangement for y-channel segments
!
! SUBROUTINE HALO_COR_ALPHA : Correct the halo-values from other faces 
!                      through 1-D interpolation (along the alpha axis)
! SUBROUTINE HALO_COR_BETA  : Correct the halo-values from other faces 
!                      through 1-D interpolation (along the beta axis)
!----------------------------------------------------------------------

!=======================================================================  
      SUBROUTINE HALO_CORREC_Z (num_halo,channel,Z3DZ,PSI)
!=======================================================================         
!     Correcting the data on z-point, which are from other cube faces.
!     attention: num_halo must be >=2 in order to use an extrapolation in halo_x and halo_y.

      INTEGER,INTENT(IN) :: num_halo  ! # of halo points to be corrected
      type(channel_t), INTENT(INOUT) :: channel   ! channel data
      LOGICAL, INTENT(IN), OPTIONAL  :: Z3DZ,PSI
      
      ! Local
      INTEGER :: MODE = 1
!     1 (Linear interpolation) 2 (Cubic Spline interpolation)
      
      INTEGER :: num_seg,mi1,mj1,mim,mip,mjm,mjp,nlen

!************************************  
      DO num_seg = 1, 4
!************************************
      mi1 = channel%seg(num_seg)%mi1  
      mim = channel%seg(num_seg)%mim     
      mip = channel%seg(num_seg)%mip  
    
      mj1 = channel%seg(num_seg)%mj1  
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp 

!--------------------------------------- X-channel
      IF (mi1.GT.mj1) THEN  
!--------------------------------------- X-channel
      nlen = mj1+2*nhalo 
      
      IF (PRESENT(Z3DZ)) THEN
       if (Z3DZ) then
          CALL HALO_X (num_halo,mi1,mim,mip,mjm,mjp,channel%seg(num_seg)%nface,  &
                       channel%seg(num_seg)%Z3DZ(:,:,nk2))
                       
          CALL HALO_COR_BETA (MODE,num_halo,nlen,mi1,mim,mip,mjm,mjp, &
                              channel%seg(num_seg)%CBETA_Z,       &
                              channel%seg(num_seg)%HBETA_Z,       &
                              channel%seg(num_seg)%JHALO_LOC_Z,   &
                              channel%seg(num_seg)%Z3DZ(:,:,nk2))             
       endif
      ENDIF

      IF (PRESENT(PSI)) THEN
       if (PSI) then
          CALL HALO_X (num_halo,mi1,mim,mip,mjm,mjp,channel%seg(num_seg)%nface,  &
                       channel%seg(num_seg)%PSI)
                       
          CALL HALO_COR_BETA (MODE,num_halo,nlen,mi1,mim,mip,mjm,mjp, &
                              channel%seg(num_seg)%CBETA_Z,       &
                              channel%seg(num_seg)%HBETA_Z,       &
                              channel%seg(num_seg)%JHALO_LOC_Z,   &
                              channel%seg(num_seg)%PSI) 
       endif
      ENDIF      


!--------------------------------------- Y-channel    
      ELSE   
!--------------------------------------- Y-channel
      nlen = mi1+2*nhalo 
      
      IF (PRESENT(Z3DZ)) THEN
       if (Z3DZ) then
          CALL HALO_Y (num_halo,mj1,mim,mip,mjm,mjp,channel%seg(num_seg)%nface,  &
                       channel%seg(num_seg)%Z3DZ(:,:,nk2))
                       
          CALL HALO_COR_ALPHA (MODE,num_halo,nlen,mj1,mim,mip,mjm,mjp,  &
                               channel%seg(num_seg)%CALPHA_Z,       &
                               channel%seg(num_seg)%HALPHA_Z,       &
                               channel%seg(num_seg)%IHALO_LOC_Z,    &
                               channel%seg(num_seg)%Z3DZ(:,:,nk2))             
       endif
      ENDIF

      IF (PRESENT(PSI)) THEN
       if (PSI) then
          CALL HALO_Y (num_halo,mj1,mim,mip,mjm,mjp,channel%seg(num_seg)%nface,  &
                       channel%seg(num_seg)%PSI)
                       
          CALL HALO_COR_ALPHA (MODE,num_halo,nlen,mj1,mim,mip,mjm,mjp,  &
                               channel%seg(num_seg)%CALPHA_Z,       &
                               channel%seg(num_seg)%HALPHA_Z,       &
                               channel%seg(num_seg)%IHALO_LOC_Z,    &
                               channel%seg(num_seg)%PSI) 
       endif
      ENDIF      
!-----------------
      ENDIF       
!-----------------    

!************************************  
      ENDDO  ! num_seg
!************************************

      END SUBROUTINE halo_correc_z

      SUBROUTINE HALO_X (num_halo,mi1,mim,mip,mjm,mjp,nface,A)
!=======================================================================       
!     Arrange data before calling HALO_COR_BETA (for variables at zeta-point)
!     (Needed because of the staggered distribution of variables.)
!------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: num_halo         ! number of halo points for interpolation
      INTEGER, INTENT(IN) :: mi1,mim,mip,mjm,mjp  ! channel segment size
      INTEGER, INTENT(IN) :: nface            ! number of face where the segment is placed
      
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp) :: TVAR
      INTEGER :: J,NPO
      
      TVAR(:,:) = A(:,:)

! Halo data arrangement      
!=========================================== 
      SELECT CASE (NFACE)
!===========================================  
!=================   
      CASE(5)
!=================         
!        data shift on the beta-direction (correct halos at Edge_M)                        
         DO NPO=1,num_halo
          DO J=mjm,mjp-1
           A(1-NPO,J) = TVAR(1-NPO,J+1)
          ENDDO
           A(1-NPO,mjp) = 1.5_r8*TVAR(1-NPO,mjp)-0.5_r8*TVAR(1-NPO,mjp-1)
         ENDDO               
!=================   
      CASE(6)
!=================  
!        data shift on the beta-direction (correct halos at Edge_P)         
         DO NPO=1,num_halo
          DO J=mjm,mjp-1
           A(mi1+NPO,J) = TVAR(mi1+NPO,J+1)
          ENDDO
           A(mi1+NPO,mjp) = 1.5_r8*TVAR(mi1+NPO,mjp)-0.5_r8*TVAR(mi1+NPO,mjp-1)
         ENDDO 

!        Data shift on the alpha-direction (correct halos at Edge_M)  
         NPO=1
         DO J=mjm,mjp
           A(1-NPO,J) = 0.5_r8*(1.5_r8*TVAR(1-NPO,J)-0.5_r8*TVAR(1-NPO-1,J) &
                               +1.5_r8*TVAR(1,J)-0.5_r8*TVAR(2,J))
                               
!          To avoid an extrapolation: 
!           A(1-NPO,J) = 0.5_r8*(TVAR(1-NPO,J)+TVAR(1,J))
         ENDDO            
         DO NPO=2,num_halo
          DO J=mjm,mjp 
           A(1-NPO,J) = TVAR(1-NPO+1,J)
          ENDDO 
         ENDDO
!===========================================      
      END SELECT   
!=========================================== 

      END SUBROUTINE halo_x

      SUBROUTINE HALO_Y (num_halo,mj1,mim,mip,mjm,mjp,nface,A)
!=======================================================================       
!     Arrange data before calling HALO_COR_ALPHA (for variables at zeta-point)
!     (Needed because of the staggered distribution of variables.)
!------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: num_halo         ! number of halo points for interpolation
      INTEGER, INTENT(IN) :: mj1,mim,mip,mjm,mjp  ! channel segment size
      INTEGER, INTENT(IN) :: nface            ! number of face where the segment is placed
      
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp) :: TVAR 
      INTEGER :: I,NPO
      
      TVAR(:,:) = A(:,:)
      
! Halo data arrangement       
!=========================================== 
      SELECT CASE (nface)
!=========================================== 
!=================     
      CASE(2)
!=================  
!        data shift on the alpha-direction (correct halos at Edge_M)                        
         DO NPO=1,num_halo
          DO I=mim,mip-1
           A(I,1-NPO) = TVAR(I+1,1-NPO)
          ENDDO
           A(mip,1-NPO) = 1.5_r8*TVAR(mip,1-NPO)-0.5_r8*TVAR(mip-1,1-NPO)
         ENDDO
!================= 
      CASE(3)
!================= 
!        data shift on the alpha-direction (correct halos at Edge_P)          
         DO NPO=1,nhalo
          DO I=mim,mip-1
           A(I,mj1+NPO) = TVAR(I+1,mj1+NPO)
          ENDDO
           A(mip,mj1+NPO) = 1.5_r8*TVAR(mip,mj1+NPO)-0.5_r8*TVAR(mip-1,mj1+NPO)
         ENDDO 
            
!        data shift on the alpha-direction (correct halos at Edge_M)                       
         DO NPO=1,num_halo
          DO I=mim,mip-1
           A(I,1-NPO) = TVAR(I+1,1-NPO)
          ENDDO
           A(mip,1-NPO) = 1.5_r8*TVAR(mip,1-NPO)-0.5_r8*TVAR(mip-1,1-NPO)
           
          DO I=mim,mip
           TVAR(I,1-NPO) = A(I,1-NPO)
          ENDDO 
         ENDDO
         
!        Data shift on the beta-direction (correct halos at Edge_M)  
         NPO=1
         DO I=mim,mip
           A(I,1-NPO) = 0.5_r8*(1.5_r8*TVAR(I,1-NPO)-0.5_r8*TVAR(I,1-NPO-1) &
                               +1.5_r8*TVAR(I,1)-0.5_r8*TVAR(I,2))

!          To avoid an extrapolation:
!           A(I,1-NPO) = 0.5_r8*(TVAR(I,1-NPO)+TVAR(I,1))
         ENDDO            
         DO NPO=2,num_halo
          DO I=mim,mip 
           A(I,1-NPO) = TVAR(I,1-NPO+1)
          ENDDO 
         ENDDO                
!=================  
      CASE(4)
!=================  
!        data shift on the alpha-direction (correct halos at Edge_P)         
         DO NPO=1,num_halo
          DO I=mim,mip-1
           A(I,mj1+NPO) = TVAR(I+1,mj1+NPO)
          ENDDO
           A(mip,mj1+NPO) = 1.5_r8*TVAR(mip,mj1+NPO)-0.5_r8*TVAR(mip-1,mj1+NPO)
         ENDDO  

!        Data shift on the beta-direction (correct halos at Edge_M)  
         NPO=1
         DO I=mim,mip
           A(I,1-NPO) = 0.5_r8*(1.5_r8*TVAR(I,1-NPO)-0.5_r8*TVAR(I,1-NPO-1) &
                               +1.5_r8*TVAR(I,1)-0.5_r8*TVAR(I,2))
                               
!          To avoid an extrapolation:
!           A(I,1-NPO) = 0.5_r8*(TVAR(I,1-NPO)+TVAR(I,1))
         ENDDO            
         DO NPO=2,num_halo
          DO I=mim,mip 
           A(I,1-NPO) = TVAR(I,1-NPO+1)
          ENDDO 
         ENDDO                
!=================   
      CASE(5)
!================= 
!        data shift on the alpha-direction (correct halos at Edge_P)          
         DO NPO=1,nhalo
          DO I=mim,mip-1
           A(I,mj1+NPO) = TVAR(I+1,mj1+NPO)
          ENDDO
           A(mip,mj1+NPO) = 1.5_r8*TVAR(mip,mj1+NPO)-0.5_r8*TVAR(mip-1,mj1+NPO)
         ENDDO 
!=================    
      CASE(6)
!=================   
!        data shift on the alpha-direction (correct halos at Edge_M)                        
         DO NPO=1,num_halo
          DO I=mim,mip-1
           A(I,1-NPO) = TVAR(I+1,1-NPO)
          ENDDO
           A(mip,1-NPO) = 1.5_r8*TVAR(mip,1-NPO)-0.5_r8*TVAR(mip-1,1-NPO)
           
          DO I=mim,mip
           TVAR(I,1-NPO) = A(I,1-NPO)
          ENDDO 
         ENDDO

!        Data shift on the beta-direction (correct halos at Edge_M)  
         NPO=1
         DO I=mim,mip
           A(I,1-NPO) = 0.5_r8*(1.5_r8*TVAR(I,1-NPO)-0.5_r8*TVAR(I,1-NPO-1) &
                               +1.5_r8*TVAR(I,1)-0.5_r8*TVAR(I,2)) 
                                     
!          To avoid an extrapolation in the direction crossing the edge:
!           A(I,1-NPO) = 0.5_r8*(TVAR(I,1-NPO)+TVAR(I,1))          
         ENDDO            
         DO NPO=2,num_halo
          DO I=mim,mip 
           A(I,1-NPO) = TVAR(I,1-NPO+1)
          ENDDO 
         ENDDO
!===========================================   
      END SELECT   
!===========================================  

      END SUBROUTINE halo_y
            
      SUBROUTINE HALO_COR_BETA (MODE,mhalo,nlen,mi1,mim,mip,mjm,mjp, &
                                CBETA,HBETA,JHALO_LOC,A)
!=======================================================================       
!     Correct the halo-values from other faces through 1-D interpolation (along the beta axis)
!------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: MODE             ! interpolation type (1: linear  2: cubic spline)
      INTEGER, INTENT(IN) :: mhalo            ! number of halo points for interpolation
      INTEGER, INTENT(IN) :: nlen,mi1,mim,mip,mjm,mjp  ! channel segment size

      REAL (KIND=r8),DIMENSION(mjm:mjp),INTENT(IN) :: &
           CBETA       ! beta values (coordinates) of the regular zeta-grid points 
      REAL (KIND=r8),DIMENSION(nhalo,mjm:mjp),INTENT(IN) :: &
           HBETA       ! beta values of the target halo points
      INTEGER,DIMENSION(nhalo,mjm:mjp),INTENT(IN) ::  &
           JHALO_LOC   ! j-index of lower points of given values for linear interpolation
      
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      INTEGER :: I,J,IPO,J1VAL
      REAL (KIND=r8),DIMENSION(nlen) :: XVAL,YVAL,XNEW,YNEW
      
      IF (MODE.EQ.2)  THEN
        DO J = mjm,mjp  
         XVAL(J+nhalo) = CBETA(J)
        ENDDO 
      ENDIF
      
!---------------------------------------------------------------------------
! EDGE_M CALCULATION
!---------------------------------------------------------------------------  
!     No correction is needed for the 1st mhalo array in EDGE_M direction 
!     because they are on the common edge between two faces.    
!=======================================       
      SELECT CASE (MODE)
!=======================================
!================================= 
!     Linear interpolation
      CASE(1)
!=================================
      DO I=2,mhalo
       IPO = I-1
       
         DO J=mjm,mjp
          J1VAL = JHALO_LOC(IPO,J)
          YNEW(J+nhalo) = FINTRP (1,HBETA(IPO,J),CBETA(J1VAL)  ,A(1-I,J1VAL) &
                                                ,CBETA(J1VAL+1),A(1-I,J1VAL+1))                                             
         ENDDO 

         DO J=mjm,mjp
          A(1-I,J) = YNEW(J+nhalo)
         ENDDO            
      ENDDO
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!================================= 
      DO I=2,mhalo
       IPO = I-1
           
         DO J=mjm,mjp
          YVAL(J+nhalo) = A(1-I,J)
          XNEW(J+nhalo) = HBETA(IPO,J) 
         ENDDO                                                                  
         CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)

         DO J=mjm,mjp
          A(1-I,J) = YNEW(J+nhalo)
         ENDDO         
      ENDDO         
!================================= 
      END SELECT
!=================================
     
!---------------------------------------------------------------------------
! EDGE_P CALCULATION
!---------------------------------------------------------------------------    

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
          YNEW(J+nhalo) = FINTRP (1,HBETA(I,J),CBETA(J1VAL)  ,A(MI1+I,J1VAL) &
                                              ,CBETA(J1VAL+1),A(MI1+I,J1VAL+1))  
         ENDDO 
        DO J=mjm,mjp
         A(MI1+I,J) = YNEW(J+nhalo)
        ENDDO         
                                                                                                    
      ENDDO
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================    
      DO I=1,mhalo

         DO J=mjm,mjp
          YVAL(J+nhalo) = A(MI1+I,J)
          XNEW(J+nhalo) = HBETA(I,J) 
         ENDDO                                                                  
         CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)  

        DO J=mjm,mjp
         A(MI1+I,J) = YNEW(J+nhalo)
        ENDDO         
            
      ENDDO               
!================================= 
      END SELECT
!=================================        
         
      END SUBROUTINE halo_cor_beta
                    
!=======================================================================       
      SUBROUTINE HALO_COR_ALPHA (MODE,mhalo,nlen,mj1,mim,mip,mjm,mjp, &
                                 CALPHA,HALPHA,IHALO_LOC,A)
!=======================================================================         
!     Correct the halo-values from other faces through 1-D interpolation (along the alpha axis)
!
!     MODE: interpolation type (1: linear  2: cubic spline)
!     EDGE_M=.T. (boundary patch - direction) vs EDGE_P=.T. (boundary patch + direction)
!     mhalo: number of halo points for interpolation (.le. nhalo) 
!     CALPHA: alpha values of the regular zeta-grid points  
!     HALPHA: alpha values of the target points
!     IHALO_LOC: i-index of lower points of given values for linear interpolation
!     A : target variable
!------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: MODE             ! interpolation type (1: linear  2: cubic spline)
      INTEGER, INTENT(IN) :: mhalo            ! number of halo points for interpolation
      INTEGER, INTENT(IN) :: nlen,mj1,mim,mip,mjm,mjp  ! channel segment size
      
      REAL (KIND=r8),DIMENSION(mim:mip),INTENT(IN) :: &
           CALPHA      ! alpha values (coordinates) of the regular zeta-grid points 
      REAL (KIND=r8),DIMENSION(mim:mip,nhalo),INTENT(IN) :: &
           HALPHA      ! alpha values of the target halo points
      INTEGER,DIMENSION(mim:mip,nhalo),INTENT(IN) ::  &
           IHALO_LOC   ! i-index of lower points of given values for linear interpolation
      
      REAL (KIND=r8),DIMENSION(mim:mip,mjm:mjp),INTENT(INOUT) :: A
      
      ! Local
      INTEGER :: I,J,JPO,I1VAL
      REAL (KIND=r8),DIMENSION(nlen) :: XVAL,YVAL,XNEW,YNEW
      
      IF (MODE.EQ.2)  THEN
        DO I = mim,mip  
         XVAL(I+nhalo) = CALPHA(I)
        ENDDO 
      ENDIF
      
!---------------------------------------------------------------------------
! EDGE_M CALCULATION
!---------------------------------------------------------------------------   
!     No correction is needed for the 1st mhalo array in EDGE_M direction 
!     because they are on the common edge between two faces.          
!=======================================       
      SELECT CASE (MODE)
!=======================================       
!================================= 
!     Linear interpolation
      CASE(1)
!=================================  
      DO J=2,mhalo
       JPO = J-1  
         DO I=mim,mip
          I1VAL = IHALO_LOC(I,JPO)
          YNEW(I+nhalo) = FINTRP (1,HALPHA(I,JPO),CALPHA(I1VAL)  ,A(I1VAL,1-J)   &
                                                 ,CALPHA(I1VAL+1),A(I1VAL+1,1-J))                             
         ENDDO  
         DO I=mim,mip
          A(I,1-J) = YNEW(I+nhalo)
         ENDDO                                                                           
      ENDDO
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!================================= 
      DO J=2,mhalo
       JPO = J-1 
         DO I=mim,mip
          YVAL(I+nhalo) = A(I,1-J)
          XNEW(I+nhalo) = HALPHA(I,JPO) 
         ENDDO                                                                  
         CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
         DO I=mim,mip
          A(I,1-J) = YNEW(I+nhalo)
         ENDDO         
      ENDDO
!================================= 
      END SELECT
!=================================          
            
!---------------------------------------------------------------------------
! EDGE_P CALCULATION
!---------------------------------------------------------------------------      
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
          YNEW(I+nhalo) = FINTRP (1,HALPHA(I,J),CALPHA(I1VAL)  ,A(I1VAL,MJ1+J) &
                                               ,CALPHA(I1VAL+1),A(I1VAL+1,MJ1+J))                                
         ENDDO
         DO I=mim,mip
          A(I,MJ1+J) = YNEW(I+nhalo) 
         ENDDO
      ENDDO
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================  
      DO J=1,mhalo
         DO I=mim,mip
          YVAL(I+nhalo) = A(I,MJ1+J)
          XNEW(I+nhalo) = HALPHA(I,J)
         ENDDO                                                                  
         CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)    
         DO I=mim,mip
          A(I,MJ1+J) = YNEW(I+nhalo) 
         ENDDO
      ENDDO       
!================================= 
      END SELECT
!=================================  
      
      END SUBROUTINE halo_cor_alpha     
      
END MODULE halo_z           