MODULE force_3d_module
! Implements a random perturbation to the potential temperature

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk2,nk3
USE constld, only: d0_0,zt,rad2deg

! Subroutines being called
USE utils,   only: indexr,ran2

IMPLICIT NONE
PRIVATE

PUBLIC :: add_rperturb

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE add_rperturb  : Add random perturbation to theta
!-----------------------------------------------------------------------

!=======================================================================
   SUBROUTINE ADD_RPERTURB (A,Z1,Z2,channel)
!=======================================================================   
!  Prepare a random perturbation 
!  (temporary use of channel%seg(num_seg)%term1)
!
!  JH: Currently, TOPOGRAPHY is not considered.
!-----------------------------------------------------------------------
      REAL (KIND=dbl_kind), INTENT(IN)  :: A,Z1,Z2
      
      type(channel_t), intent(inout) :: channel   ! channel data
      
!     Local variables
      LOGICAL :: LF
      REAL (KIND=dbl_kind) :: SUMV(NK2,4),NDATA(NK2,4),AMEAN(NK2)
      REAL (KIND=dbl_kind) :: AMAX(NK2),AMIN(NK2),AMAX_all(NK2),AMIN_all(NK2)
      REAL (KIND=dbl_kind) :: AMAX0,AMIN0,ADD1,ADD2

      INTEGER :: NSUM,IDUM,I,J,K,K1,K2,num_seg,mi1,mj1,chl,chn
      
      chn = 1    ! Assume the channel width is 1 (i.e., mi1 or mj1 = 1).

      ! Set the lateral boundary for applying the perturbation.
      ADD1 = -10._dbl_kind/RAD2DEG
      ADD2 =  10._dbl_kind/RAD2DEG

! Get vertical indices
      K1 = INDEXR(Z1,NK3,ZT,LF) + 1
      K2 = INDEXR(Z2,NK3,ZT,LF)

!******************************  
      DO num_seg = 1, 4
!******************************
      mi1 = channel%seg(num_seg)%mi1
      mj1 = channel%seg(num_seg)%mj1   

! Select random deviates between -1 and 1
      IDUM = 5
      DO K = K1, K2
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%term1(I,J,K) = RAN2(IDUM)
        ENDDO
       ENDDO
      ENDDO

! Calculate the sum of perturbation at each channel segment and level
      IF (mi1.gt.mj1) THEN
        ! x-channel segment
        DO K = K1, K2
         SUMV(k,num_seg) = D0_0
         DO chl = 1, MI1
          SUMV(k,num_seg) = SUMV(k,num_seg) + channel%seg(num_seg)%term1(chl,chn,K)
         ENDDO
          NDATA(k,num_seg) = float(mi1)
        ENDDO  
      ELSE
        ! y-channel segment
        DO K = K1, K2
         SUMV(k,num_seg) = D0_0
         DO chl = 1, MJ1
          SUMV(k,num_seg) = SUMV(k,num_seg) + channel%seg(num_seg)%term1(chn,chl,K)
         ENDDO
          NDATA(k,num_seg) = float(mj1)
        ENDDO        
      ENDIF
!***********************  
      ENDDO  ! num_seg 
!***********************

! Calculate the perturbation mean at each level
      
      DO K = K1, K2 
       AMEAN(K) = (SUMV(k,1)+SUMV(k,2)+SUMV(k,3)+SUMV(k,4)) &
                 /(NDATA(k,1)+NDATA(k,2)+NDATA(k,3)+NDATA(k,4))
      ENDDO 

      AMAX_all(:) = D0_0
      AMIN_all(:) = D0_0
       
!******************************  
      DO num_seg = 1, 4
!******************************
      mi1 = channel%seg(num_seg)%mi1
      mj1 = channel%seg(num_seg)%mj1   
      
      DO K = K1, K2
       DO J = 1, MJ1
        DO I = 1, MI1
         channel%seg(num_seg)%term1(I,J,K) = channel%seg(num_seg)%term1(I,J,K) - AMEAN(K)
        ENDDO
       ENDDO
      ENDDO
      
! Find the max & min values at each level

      DO K = K1, K2
       AMAX(k) = MAXVAL(channel%seg(num_seg)%term1(1:mi1,1:mj1,k))
       AMIN(k) = MINVAL(channel%seg(num_seg)%term1(1:mi1,1:mj1,k))
      ENDDO 
      
      DO K = K1, K2
       AMAX_all(K) = MAX(AMAX_all(K),AMAX(k))
       AMIN_all(K) = MIN(AMIN_all(K),AMIN(k))
      ENDDO        
!******************************  
      ENDDO  ! num_seg 
!******************************

      AMAX0 = maxval(AMAX_all(k1:k2))
      AMIN0 = minval(AMIN_all(k1:k2))

! Normalize to perturbation magnitude
!******************************  
      DO num_seg = 1, 4
!******************************
      mi1 = channel%seg(num_seg)%mi1
      mj1 = channel%seg(num_seg)%mj1   

      DO K = K1, K2
       DO J = 1, MJ1
        DO I = 1, MI1
         IF (channel%seg(num_seg)%term1(I,J,K).LT.d0_0) &
             channel%seg(num_seg)%term1(I,J,K) = -channel%seg(num_seg)%term1(I,J,K)*A/AMIN0
         IF (channel%seg(num_seg)%term1(I,J,K).GT.d0_0) &
             channel%seg(num_seg)%term1(I,J,K) =  channel%seg(num_seg)%term1(I,J,K)*A/AMAX0
        ENDDO
       ENDDO
      ENDDO

      DO K = K1, K2
       DO J = 1, MJ1
        DO I = 1, MI1
         IF (RLAT_T(I,J).GE.ADD1 .AND. RLAT_T(I,J).LE.ADD2) THEN
          channel%seg(num_seg)%TH3D(I,J,K) = channel%seg(num_seg)%TH3D(I,J,K) &
                                           + channel%seg(num_seg)%term1(I,J,K)
         ENDIF
        ENDDO
       ENDDO
      ENDDO

!***********************  
      ENDDO  ! num_seg 
!***********************

   END SUBROUTINE add_rperturb

END MODULE force_3d_module
