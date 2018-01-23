MODULE halo_vort
! Routines that correct the halo-point vector components from neighboring face.

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk1,nk2,nhalo
USE constld, only: d2_0,d3_0,d0_5

! Functions used
USE utils,   only: fintrp,cubicint

IMPLICIT NONE
PRIVATE

PUBLIC :: vort_comm_pre,  &  ! procedure before data communication
          vort_comm_post     ! procedure after data communication

CONTAINS

! Local Subroutines:
!---------------------------------------------------------------------------
! SUBROUTINE vort_comm_pre  : replace face edge values 
!                             with the ones in RLL before sending data 
! SUBROUTINE vort_comm_post : recover the data after receiving data
!                             from RLL to CUBE 
!
! Correct the halo-values from other faces
! SUBROUTINE HALO_COR_VORT_ALPHA : along the alpha axis
! SUBROUTINE HALO_COR_VORT_BETA  : along the beta axis
!
! SUBROUTINE HALO_ALPHA : called from HALO_COR_VORT_ALPHA   
! SUBROUTINE HALO_BETA  : called from HALO_COR_VORT_BETA
!---------------------------------------------------------------------------
!     D2_0  = 2.0_dbl_kind
!     D3_0  = 3.0_dbl_kind
!     D0_5  = 0.5_dbl_kind

!======================================================================= 
      SUBROUTINE VORT_COMM_PRE (channel)
!=======================================================================  
!     Procedure before data exchange between faces 
!     
!     HALO_COR_VORT_ALPHA & HALO_COR_VORT_BETA:
!     1. obtain inner halo values (for sending) at q-points
!     2. find the corresponding values at the halo points of neighboring face
!     3. map the values of 2 to obtain the values on RLL grids
!
!     Replace Z3DX and Z3DY only at the edge of cube faces 
!--------------------------------------------------------------------------------      
      type(channel_t), intent(inout) :: channel   ! channel data
                
      INTEGER :: MODE = 1
!     1 (Linear interpolation) 2 (Cubic Spline interpolation)
    
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: Z3DX_ALPHA,Z3DY_ALPHA
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: Z3DX_BETA,Z3DY_BETA 
      
      INTEGER :: I,J,K,NPO
      INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp

!****************************************************************************   
      DO num_seg = 1, 4
!**************************************************************************** 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim    
      mip    = channel%seg(num_seg)%mip    
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm    = channel%seg(num_seg)%mjm    
      mjp    = channel%seg(num_seg)%mjp          
      
      ALLOCATE(Z3DX_BETA(2*nhalo,mjm:mjp,NK1-1))
      ALLOCATE(Z3DY_BETA(2*nhalo,mjm:mjp,NK1-1))
      
      ALLOCATE(Z3DX_ALPHA(mim:mip,2*nhalo,NK1-1))
      ALLOCATE(Z3DY_ALPHA(mim:mip,2*nhalo,NK1-1))
   
!     Save the original data      
      DO K = 2, nk1 
       DO J = mjm, mjp
        DO I = mim, mip
         channel%seg(num_seg)%Z3DX0(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K)
         channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K)
        ENDDO
       ENDDO
      ENDDO
         
!=======================================
      IF (mi1.GT.mj1) THEN  
!=======================================  x-array     
 
      CALL HALO_COR_VORT_BETA (mi1,mim,mip,mjm,mjp,nk1-1,MODE,       &
                               channel%seg(num_seg)%AM_VORT_BETA_T,  &
                               channel%seg(num_seg)%CBETA_T,         &
                               channel%seg(num_seg)%HBETA_T,         &
                               channel%seg(num_seg)%JHALO_LOC_T,     &
                               channel%seg(num_seg)%Z3DX(:,:,2:NK1), &
                               channel%seg(num_seg)%Z3DY(:,:,2:NK1), &
                               Z3DX_BETA,Z3DY_BETA)
                              
        DO I = 1, nhalo
         DO K = 2, NK1 
          DO J = mjm, mjp
           channel%seg(num_seg)%Z3DX(I,J,K) = Z3DX_BETA(I,J,K-1)
           channel%seg(num_seg)%Z3DY(I,J,K) = Z3DY_BETA(I,J,K-1)
          ENDDO
         ENDDO
        ENDDO  
      
        DO I = mi1-nhalo+1, mi1
         NPO = I - mi1 + 2*nhalo
         DO K = 2, NK1
          DO J = mjm, mjp
           channel%seg(num_seg)%Z3DX(I,J,K) = Z3DX_BETA(NPO,J,K-1)
           channel%seg(num_seg)%Z3DY(I,J,K) = Z3DY_BETA(NPO,J,K-1)
          ENDDO
         ENDDO
        ENDDO     
             
!======================================= 
      ELSE    
!=======================================  y-array

      CALL HALO_COR_VORT_ALPHA (mj1,mim,mip,mjm,mjp,nk1-1,MODE,       &
                                channel%seg(num_seg)%AM_VORT_ALPHA_T, &
                                channel%seg(num_seg)%CALPHA_T,        &
                                channel%seg(num_seg)%HALPHA_T,        &
                                channel%seg(num_seg)%IHALO_LOC_T,     &  
                                channel%seg(num_seg)%Z3DX(:,:,2:NK1), &
                                channel%seg(num_seg)%Z3DY(:,:,2:NK1), &
                                Z3DX_ALPHA,Z3DY_ALPHA)
      
        DO J = 1, nhalo
         DO K = 2, NK1 
          DO I = mim, mip
           channel%seg(num_seg)%Z3DX(I,J,K) = Z3DX_ALPHA(I,J,K-1)
           channel%seg(num_seg)%Z3DY(I,J,K) = Z3DY_ALPHA(I,J,K-1)
          ENDDO
         ENDDO
        ENDDO    
      
        DO J = mj1-nhalo+1, mj1
         NPO = J - mj1 + 2*nhalo
         DO K = 2, NK1
          DO I = mim, mip
           channel%seg(num_seg)%Z3DX(I,J,K) = Z3DX_ALPHA(I,NPO,K-1)
           channel%seg(num_seg)%Z3DY(I,J,K) = Z3DY_ALPHA(I,NPO,K-1)
          ENDDO
         ENDDO
        ENDDO                
!=======================================       
      ENDIF       
!=======================================   

      DEALLOCATE(Z3DX_BETA)
      DEALLOCATE(Z3DY_BETA)
      
      DEALLOCATE(Z3DX_ALPHA)
      DEALLOCATE(Z3DY_ALPHA)

!***************************************   
      ENDDO  ! num_seg 
!*************************************** 
      
      END SUBROUTINE vort_comm_pre   

!=======================================================================            
      SUBROUTINE VORT_COMM_POST (channel)
!=======================================================================  
!     Procedure after receiving vorticity components from other processes 
!     
!     1. Recover own internal values
!     2. Remap the halo values at q-points (RLL --> CUBE) 
!     3. Obtain the halo values at ksi- and eta-points
!-------------------------------------------------------------------------------- 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: Z3DX_ALPHA,Z3DY_ALPHA
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: Z3DX_BETA,Z3DY_BETA 
      
      INTEGER, PARAMETER :: INTERPOL_ITER = 10 
      INTEGER :: I,J,K,IPOM,IPOP,JPOM,JPOP,ITER
      INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp
      
!****************************************************************************   
      DO num_seg = 1, 4
!**************************************************************************** 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim    
      mip    = channel%seg(num_seg)%mip    
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm    = channel%seg(num_seg)%mjm    
      mjp    = channel%seg(num_seg)%mjp          
      
      ALLOCATE(Z3DX_BETA(2*nhalo,mjm:mjp,NK1-1))
      ALLOCATE(Z3DY_BETA(2*nhalo,mjm:mjp,NK1-1))
      
      ALLOCATE(Z3DX_ALPHA(mim:mip,2*nhalo,NK1-1))
      ALLOCATE(Z3DY_ALPHA(mim:mip,2*nhalo,NK1-1))
            
!=========================================  X-channel
      IF (mi1.GT.mj1) THEN  
!=========================================  X-channel

!--------------------------------------
!     1. Recover own internal values    
!--------------------------------------      
      DO K = 2, NK1
       DO J = mjm, mjp
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX0(I,J,K)
         channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY0(I,J,K)
        ENDDO
       ENDDO
      ENDDO
!--------------------------------------------------------------        
!     2. Remap the vorticity components (RLL --> Cube grids) 
!--------------------------------------------------------------
      DO K = 2, NK1
       DO J = mjm, mjp
        DO I = 1, nhalo
         IPOM =   1 - I
         IPOP = MI1 + I
         
         Z3DX_BETA(I,J,K-1) = &
           channel%seg(num_seg)%AMI_T(1,IPOM,J)*channel%seg(num_seg)%Z3DX(IPOM,J,K) &
         + channel%seg(num_seg)%AMI_T(2,IPOM,J)*channel%seg(num_seg)%Z3DY(IPOM,J,K)
         
         Z3DY_BETA(I,J,K-1) = &
           channel%seg(num_seg)%AMI_T(3,IPOM,J)*channel%seg(num_seg)%Z3DX(IPOM,J,K) &
         + channel%seg(num_seg)%AMI_T(4,IPOM,J)*channel%seg(num_seg)%Z3DY(IPOM,J,K)
         
         Z3DX_BETA(I+nhalo,J,K-1) = &
           channel%seg(num_seg)%AMI_T(1,IPOP,J)*channel%seg(num_seg)%Z3DX(IPOP,J,K) &
         + channel%seg(num_seg)%AMI_T(2,IPOP,J)*channel%seg(num_seg)%Z3DY(IPOP,J,K)
         
         Z3DY_BETA(I+nhalo,J,K-1) = &
           channel%seg(num_seg)%AMI_T(3,IPOP,J)*channel%seg(num_seg)%Z3DX(IPOP,J,K) &
         + channel%seg(num_seg)%AMI_T(4,IPOP,J)*channel%seg(num_seg)%Z3DY(IPOP,J,K)
        ENDDO
       ENDDO  
      ENDDO
!--------------------------------------------------------------      
!     3. Obtain the vorticity components from q-point values 
!--------------------------------------------------------------
!     nhalo <-- 1 |    CUBE FACE     | nhalo+1 --> 2*nhalo 

!     Process including Western Edge          
       DO I = 1, nhalo
        IPOM = 1 - I
        DO K = 2, NK1
         DO J = mjm, mjp-1
          channel%seg(num_seg)%Z3DX(IPOM,J,K) = D0_5*(Z3DX_BETA(I,J,K-1)+Z3DX_BETA(I,J+1,K-1))
         ENDDO
          channel%seg(num_seg)%Z3DX(IPOM,mjp,K) = D2_0*Z3DX_BETA(I,mjp,K-1) &
                                                - channel%seg(num_seg)%Z3DX(IPOM,mjp-1,K)
        ENDDO  
       ENDDO
         
       DO I = 2, nhalo
        IPOM = 1 - I
        DO K = 2, NK1 
         DO J = mjm, mjp 
          channel%seg(num_seg)%Z3DY(IPOM,J,K) = D0_5*(Z3DY_BETA(I,J,K-1)+Z3DY_BETA(I-1,J,K-1))
         ENDDO
        ENDDO  
       ENDDO         
       DO K = 2, NK1 
        DO J = mjm, mjp        
         channel%seg(num_seg)%Z3DY(0,J,K) = (D2_0*Z3DY_BETA(1,J,K-1) &
                                          + channel%seg(num_seg)%Z3DY(1,J,K))/D3_0
        ENDDO  
       ENDDO
        
!     Process including Estern Edge   
       DO I = nhalo+1, 2*nhalo
        IPOP = MI1 + I - nhalo
        DO K = 2, NK1
         DO J = mjm, mjp-1
          channel%seg(num_seg)%Z3DX(IPOP,J,K) =   &
                 D0_5*(Z3DX_BETA(I,J,K-1)+Z3DX_BETA(I,J+1,K-1))
         ENDDO
          channel%seg(num_seg)%Z3DX(IPOP,mjp,K) = &
                 D2_0*Z3DX_BETA(I,mjp,K-1)-channel%seg(num_seg)%Z3DX(IPOP,mjp-1,K)
        ENDDO  
       ENDDO
         
       DO I = nhalo+1, 2*nhalo-1
        IPOP = MI1 + I - nhalo
        DO K = 2, NK1 
         DO J = mjm, mjp 
          channel%seg(num_seg)%Z3DY(IPOP,J,K) = &
                 D0_5*(Z3DY_BETA(I,J,K-1)+Z3DY_BETA(I+1,J,K-1))
         ENDDO
        ENDDO  
       ENDDO
       
       DO ITER = 1, INTERPOL_ITER   
       DO K = 2, NK1 
        DO J = mjm, mjp 
         Z3DY_BETA(nhalo+1,J,K-1) = D0_5*(channel%seg(num_seg)%Z3DY(MI1,J,K) &
                                  + channel%seg(num_seg)%Z3DY(MI1+1,J,K))
        ENDDO
       ENDDO  
       DO K = 2, NK1 
        DO J = mjm, mjp 
         channel%seg(num_seg)%Z3DY(MI1+1,J,K) = &
                D0_5*(Z3DY_BETA(nhalo+1,J,K-1)+Z3DY_BETA(nhalo+2,J,K-1))
        ENDDO
       ENDDO 
       ENDDO  ! ITER              
       
       I    = nhalo + nhalo
       IPOP = MI1   + nhalo        
       DO K = 2, NK1 
        DO J = mjm, mjp 
         channel%seg(num_seg)%Z3DY(IPOP,J,K) = &
                D2_0*Z3DY_BETA(I,J,K-1)-channel%seg(num_seg)%Z3DY(IPOP-1,J,K)
        ENDDO  
       ENDDO
             
!=========================================  Y-channel
      ELSE    
!=========================================  Y-channel

!--------------------------------------
!     1. Recover own internal values     
!--------------------------------------     
      DO K = 2, NK1
       DO I = mim, mip
        DO J = 1, mj1
         channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX0(I,J,K-1)
         channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY0(I,J,K-1)
        ENDDO
       ENDDO
      ENDDO
!-----------------------------------------------------                 
!     2. Remap the vorticity (RLL --> Cube grids) 
!-----------------------------------------------------
      DO  K= 2, NK1
       DO I = mim, mip
        DO J = 1, nhalo  
         JPOM =   1 - J
         JPOP = MJ1 + J
         
         Z3DX_ALPHA(I,J,K-1) = &
           channel%seg(num_seg)%AMI_T(1,I,JPOM)*channel%seg(num_seg)%Z3DX(I,JPOM,K) &
         + channel%seg(num_seg)%AMI_T(2,I,JPOM)*channel%seg(num_seg)%Z3DY(I,JPOM,K)
         
         Z3DY_ALPHA(I,J,K-1) = &
           channel%seg(num_seg)%AMI_T(3,I,JPOM)*channel%seg(num_seg)%Z3DX(I,JPOM,K) &
         + channel%seg(num_seg)%AMI_T(4,I,JPOM)*channel%seg(num_seg)%Z3DY(I,JPOM,K)
         
         Z3DX_ALPHA(I,J+nhalo,K-1) = &
           channel%seg(num_seg)%AMI_T(1,I,JPOP)*channel%seg(num_seg)%Z3DX(I,JPOP,K) &
         + channel%seg(num_seg)%AMI_T(2,I,JPOP)*channel%seg(num_seg)%Z3DY(I,JPOP,K)
         
         Z3DY_ALPHA(I,J+nhalo,K-1) = &
           channel%seg(num_seg)%AMI_T(3,I,JPOP)*channel%seg(num_seg)%Z3DX(I,JPOP,K) &
         + channel%seg(num_seg)%AMI_T(4,I,JPOP)*channel%seg(num_seg)%Z3DY(I,JPOP,K)
        ENDDO
       ENDDO  
      ENDDO

!--------------------------------------------------------------      
!     3. Obtain the vorticity components from q-point values 
!--------------------------------------------------------------

!     Process including Southern Edge          
       DO J = 1, nhalo
        JPOM = 1 - J
        DO K = 2, NK1
         DO I = mim, mip-1
          channel%seg(num_seg)%Z3DY(I,JPOM,K) = &
                 D0_5*(Z3DY_ALPHA(I,J,K-1)+Z3DY_ALPHA(I+1,J,K-1))
         ENDDO
          channel%seg(num_seg)%Z3DY(mip,JPOM,K) = D2_0*Z3DY_ALPHA(mip,J,K-1) &
                                                - channel%seg(num_seg)%Z3DY(mip-1,JPOM,K)
        ENDDO  
       ENDDO
         
       DO J = 2, nhalo
        JPOM = 1 - J
        DO K = 2, NK1 
         DO I = mim, mip 
          channel%seg(num_seg)%Z3DX(I,JPOM,K) = &
                 D0_5*(Z3DX_ALPHA(I,J,K-1)+Z3DX_ALPHA(I,J-1,K-1))
         ENDDO
        ENDDO  
       ENDDO         
       DO K = 2, NK1 
        DO I = mim, mip 
         channel%seg(num_seg)%Z3DX(I,0,K) = (D2_0*Z3DX_ALPHA(I,1,K-1) &
                                          + channel%seg(num_seg)%Z3DX(I,1,K))/D3_0
        ENDDO  
       ENDDO

!     Process including Northern Edge   
       DO J = nhalo+1, 2*nhalo
        JPOP = MJ1 + J - nhalo
        DO K = 2, NK1
         DO I = mim, mip-1
          channel%seg(num_seg)%Z3DY(I,JPOP,K) = &
                 D0_5*(Z3DY_ALPHA(I,J,K-1)+Z3DY_ALPHA(I+1,J,K-1))
         ENDDO
          channel%seg(num_seg)%Z3DY(mip,JPOP,K) = D2_0*Z3DY_ALPHA(mip,J,K-1) &
                                                - channel%seg(num_seg)%Z3DY(mip-1,JPOP,K)
        ENDDO  
       ENDDO
         
       DO J = nhalo+1, 2*nhalo-1
        JPOP = MJ1 + J - nhalo
        DO K = 2, NK1 
         DO I = mim, mip 
          channel%seg(num_seg)%Z3DX(I,JPOP,K) = &
                 D0_5*(Z3DX_ALPHA(I,J,K-1)+Z3DX_ALPHA(I,J+1,K-1))
         ENDDO
        ENDDO  
       ENDDO  

       DO ITER = 1, INTERPOL_ITER   
        DO K = 2, NK1 
         DO I = mim, mip 
          Z3DX_ALPHA(I,nhalo+1,K-1) = D0_5*(channel%seg(num_seg)%Z3DX(I,mj1,K) &
                                    + channel%seg(num_seg)%Z3DX(I,mj1+1,K))
         ENDDO
        ENDDO  
        DO K = 2, NK1 
         DO I = mim, mip 
          channel%seg(num_seg)%Z3DX(I,mj1+1,K) = &
                 D0_5*(Z3DX_ALPHA(I,nhalo+1,K-1)+Z3DX_ALPHA(I,nhalo+2,K-1))
         ENDDO
        ENDDO           
       ENDDO  ! ITER    
                     
       J    = nhalo + nhalo
       JPOP = MJ1   + nhalo        
       DO K = 2, NK1 
        DO I = mim, mip 
         channel%seg(num_seg)%Z3DX(I,JPOP,K) = &
                D2_0*Z3DX_ALPHA(I,J,K-1)-channel%seg(num_seg)%Z3DX(I,JPOP-1,K)
        ENDDO  
       ENDDO
!=========================================  Y-channel     
      ENDIF       
!=========================================  Y-channel

      DEALLOCATE(Z3DX_BETA)
      DEALLOCATE(Z3DY_BETA)
      
      DEALLOCATE(Z3DX_ALPHA)
      DEALLOCATE(Z3DY_ALPHA)

!***************************************   
      ENDDO  ! num_seg 
!*************************************** 
      
      END SUBROUTINE vort_comm_post  
                        
!=======================================================================  
      SUBROUTINE HALO_COR_VORT_ALPHA (mj1,mim,mip,mjm,mjp,KDIMN,MODE,   &
                                      AM_VORT_ALPHA,CALPHA,HALPHA,IHALO_LOC, &
                                      AKSI,AETA,AKSI0,AETA0) 
!=======================================================================        
      INTEGER, INTENT(IN) :: mj1,mim,mip,mjm,mjp,nlen  ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN   ! vertical array size
      INTEGER, INTENT(IN) :: MODE    ! interpolation type (1: linear  2: cubic spline)

      REAL(KIND=dbl_kind),DIMENSION(4,mim:mip,nhalo*2),INTENT(IN) :: &
         AM_VORT_ALPHA  ! transformation matrix for a vector in q-point 
      REAL(KIND=dbl_kind),DIMENSION(mim:mip),INTENT(IN) ::           &
         CALPHA         ! curvilinear coordinates, alpha [rad]
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,nhalo),INTENT(IN) ::     &
         HALPHA         ! location of halo points: alpha value [rad] 
      INTEGER,DIMENSION(mim:mip,nhalo),INTENT(IN) ::                 &
         IHALO_LOC     ! location of halo points: nearby i-index   


      REAL(KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(IN) ::  &
         AKSI,AETA     ! original field
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2*nhalo,kdimn),INTENT(OUT) :: &
         AKSI0,AETA0   ! corrected field for communication  

      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2*nhalo,kdimn) :: AKSI_ALPHA,AETA_ALPHA
      INTEGER :: I,J,K,NPO

!     1. Obtain the vorticity components at q-point 
   
      DO K=1,KDIMN

       DO J=2,nhalo
        DO I=mim,mip
         AKSI_ALPHA(I,J,K) = D0_5*(AKSI(I,J,K)+AKSI(I,J-1,K))
        ENDDO 
       ENDDO  
       DO I=mim,mip
        AKSI_ALPHA(I,1,K)  = D2_0*AKSI(I,1,K)-AKSI_ALPHA(I,2,K)
       ENDDO
              
       DO J=mj1-nhalo+1,mj1
        NPO = J - mj1 + 2*nhalo      
        DO I=mim,mip
         AKSI_ALPHA(I,NPO,K) = D0_5*(AKSI(I,J,K)+AKSI(I,J-1,K))
        ENDDO 
       ENDDO          
       
       DO J=1,nhalo
        DO I=mim+1,mip         
         AETA_ALPHA(I,J,K)  = D0_5*(AETA(I,J,K)+AETA(I-1,J,K))
        ENDDO
        AETA_ALPHA(mim,J,K) = D2_0*AETA(mim,J,K)-AETA_ALPHA(mim+1,J,K) 
       ENDDO 
       
       DO J=mj1-nhalo+1,mj1
        NPO = J - mj1 + 2*nhalo
        DO I=mim+1,mip         
         AETA_ALPHA(I,NPO,K)  = D0_5*(AETA(I,J,K)+AETA(I-1,J,K))
        ENDDO 
        AETA_ALPHA(mim,NPO,K) = D2_0*AETA(mim,J,K)-AETA_ALPHA(mim+1,NPO,K)
       ENDDO 
                             
      ENDDO  ! k-loop

!     2. Obtain the vorticity components needed for the neighboring face through 1-D interpolation
      
      nlen = mi1+2*nhalo
      CALL HALO_ALPHA (mim,mip,nlen,KDIMN,MODE,CALPHA,HALPHA,IHALO_LOC,AKSI_ALPHA)
      CALL HALO_ALPHA (mim,mip,nlen,KDIMN,MODE,CALPHA,HALPHA,IHALO_LOC,AETA_ALPHA)  
        
!     3. Map the vorticity to obtain the values on RLL grids

      DO K=1,KDIMN 
       DO I=mim,mip
       
        DO J=1,nhalo
         AKSI0(I,J,K) = AM_VORT_ALPHA(1,I,J)*AKSI_ALPHA(I,J,K) &
                      + AM_VORT_ALPHA(2,I,J)*AETA_ALPHA(I,J,K)
         AETA0(I,J,K) = AM_VORT_ALPHA(3,I,J)*AKSI_ALPHA(I,J,K) &
                      + AM_VORT_ALPHA(4,I,J)*AETA_ALPHA(I,J,K)
        ENDDO
        DO J=nhalo+1,2*nhalo
         AKSI0(I,J,K) = AM_VORT_ALPHA(1,I,J)*AKSI_ALPHA(I,J,K) &
                      + AM_VORT_ALPHA(2,I,J)*AETA_ALPHA(I,J,K)
         AETA0(I,J,K) = AM_VORT_ALPHA(3,I,J)*AKSI_ALPHA(I,J,K) &
                      + AM_VORT_ALPHA(4,I,J)*AETA_ALPHA(I,J,K)
        ENDDO

       ENDDO  
      ENDDO 
                  
      END SUBROUTINE halo_cor_vort_alpha 

!=======================================================================      
      SUBROUTINE HALO_COR_VORT_BETA (mi1,mim,mip,mjm,mjp,KDIMN,MODE,     &
                                     AM_VORT_BETA,CBETA,HBETA,JHALO_LOC, &
                                     AKSI,AETA,AKSI0,AETA0) 
!=======================================================================        
      INTEGER, INTENT(IN) :: mi1,mim,mip,mjm,mjp  ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN   ! vertical array size
      INTEGER, INTENT(IN) :: MODE    ! interpolation type (1: linear  2: cubic spline)

      REAL(KIND=dbl_kind),DIMENSION(4,nhalo*2,mjm:mjp),INTENT(IN) :: &
         AM_VORT_BETA  ! transformation matrix for a vector in q-point 
      REAL(KIND=dbl_kind),DIMENSION(mjm:mjp),INTENT(IN) ::           &
         CBETA         ! curvilinear coordinates, beta [rad]
      REAL(KIND=dbl_kind),DIMENSION(nhalo,mjm:mjp),INTENT(IN) ::     &
         HBETA         ! location of halo points: beta value [rad] 
      INTEGER,DIMENSION(nhalo,mjm:mjp),INTENT(IN) ::                 &
         JHALO_LOC     ! location of halo points: nearby j-index   
                  
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(IN) ::  &
         AKSI,AETA     ! original field
      REAL(KIND=dbl_kind),DIMENSION(2*nhalo,mjm:mjp,kdimn),INTENT(OUT) :: &
         AKSI0,AETA0   ! corrected field for communication  

      REAL(KIND=dbl_kind),DIMENSION(2*nhalo,mjm:mjp,kdimn) :: AKSI_BETA,AETA_BETA 
      INTEGER :: nlen,I,J,K,NPO
             
!     1. Obtain the vorticity components at q-point (every processes) 
   
      DO K=1,KDIMN

       DO I=1,nhalo
        DO J=mjm+1,mjp         
         AKSI_BETA(I,J,K)  = D0_5*(AKSI(I,J,K)+AKSI(I,J-1,K))
        ENDDO
        AKSI_BETA(I,mjm,K) = D2_0*AKSI(I,mjm,K)-AKSI_BETA(I,mjm+1,K) 
       ENDDO 
       
       DO I=mi1-nhalo+1,mi1
        NPO = I - mi1 + 2*nhalo
        DO J=mjm+1,mjp         
         AKSI_BETA(NPO,J,K)  = D0_5*(AKSI(I,J,K)+AKSI(I,J-1,K))
        ENDDO
        AKSI_BETA(NPO,mjm,K) = D2_0*AKSI(I,mjm,K)-AKSI_BETA(NPO,mjm+1,K)  
       ENDDO 
                     
       DO I=2,nhalo
        DO J=mjm,mjp
         AETA_BETA(I,J,K) = D0_5*(AETA(I,J,K)+AETA(I-1,J,K))
        ENDDO 
       ENDDO  
       DO J=mjm,mjp
        AETA_BETA(1,J,K) = D2_0*AETA(1,J,K)-AETA_BETA(2,J,K)
       ENDDO
              
       DO I=mi1-nhalo+1,mi1
        NPO = I - mi1 + 2*nhalo      
        DO J=mjm,mjp
         AETA_BETA(NPO,J,K) = D0_5*(AETA(I,J,K)+AETA(I-1,J,K))
        ENDDO 
       ENDDO          
                
      ENDDO  ! k-loop

!     2. Obtain the vorticity components needed for the neighboring face through 1-D interpolation
      
      nlen = mj1+2*nhalo
      CALL HALO_BETA (mjm,mjp,nlen,KDIMN,MODE,CBETA,HBETA,JHALO_LOC,AKSI_BETA)
      CALL HALO_BETA (mjm,mjp,nlen,KDIMN,MODE,CBETA,HBETA,JHALO_LOC,AETA_BETA)
      
!     3. map the vorticity components to obtain the values on RLL grids

      DO K=1,KDIMN 
       DO J=mjm,mjp
       
        DO I=1,nhalo
         AKSI0(I,J,K) = AM_VORT_BETA(1,I,J)*AKSI_BETA(I,J,K) &
                      + AM_VORT_BETA(2,I,J)*AETA_BETA(I,J,K)
         AETA0(I,J,K) = AM_VORT_BETA(3,I,J)*AKSI_BETA(I,J,K) &
                      + AM_VORT_BETA(4,I,J)*AETA_BETA(I,J,K)
        ENDDO
        DO I=nhalo+1,2*nhalo
         AKSI0(I,J,K) = AM_VORT_BETA(1,I,J)*AKSI_BETA(I,J,K) &
                      + AM_VORT_BETA(2,I,J)*AETA_BETA(I,J,K)
         AETA0(I,J,K) = AM_VORT_BETA(3,I,J)*AKSI_BETA(I,J,K) &
                      + AM_VORT_BETA(4,I,J)*AETA_BETA(I,J,K)
        ENDDO

       ENDDO  
      ENDDO

      END SUBROUTINE halo_cor_vort_beta      

!=======================================================================  
      SUBROUTINE HALO_ALPHA (mim,mip,nlen,KDIMN,MODE,CALPHA,HALPHA,IHALO_LOC,A)   
!=======================================================================           
!     (ALONG ALPHA AXIS) Determine halo-values through 1-D interpolation.

      INTEGER, INTENT(IN) :: mim,mip,nlen ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN        ! vertical array size  
      INTEGER, INTENT(IN) :: MODE         ! interpolation type (1: linear  2: cubic spline) 
      
      ! alpha values of the regular q-grid points
      REAL (KIND=dbl_kind),DIMENSION(mim:mip),INTENT(IN) ::        &
        CALPHA      ! alpha coordinates of the regular q-grid points
      REAL (KIND=dbl_kind),DIMENSION(mim:mip,nhalo),INTENT(IN) ::  &
        HALPHA      ! alpha values of the target halo points
      INTEGER,DIMENSION(mim:mip,nhalo),INTENT(IN) ::               &
        IHALO_LOC   ! i-index of CALPHA that is just less than (or equal to) HALPHA  
      
      ! target variable
      REAL (KIND=dbl_kind),DIMENSION(mim:mip,2*nhalo,KDIMN),INTENT(INOUT) :: A
      
      ! Local variables 
      REAL (KIND=dbl_kind),DIMENSION(nlen) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=dbl_kind),DIMENSION(mim:mip,KDIMN) :: HMV,HPV  
      INTEGER :: I,J,K,NPO,I1VAL,I2VAL
            
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
      DO J=1,nhalo     
       NPO = (2*nhalo+1) - J     
              
       DO I=mim,mip
        I1VAL = IHALO_LOC(I,J)
        I2VAL = I1VAL + 1
        DO K=1,KDIMN 
         HMV(I,K)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,J,K)   &
                                        ,CALPHA(I2VAL),A(I2VAL,J,K))
         HPV(I,K)= FINTRP (1,HALPHA(I,J),CALPHA(I1VAL),A(I1VAL,NPO,K) &
                                        ,CALPHA(I2VAL),A(I2VAL,NPO,K)) 
        ENDDO                              
       ENDDO  
       DO K=1,KDIMN 
        DO I=mim,mip
         A(I,J,K)   = HMV(I,K)
         A(I,NPO,K) = HPV(I,K)
        ENDDO
       ENDDO  
              
      ENDDO  ! J-loop     

!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================               
      DO J=1,nhalo     
       NPO = (2*nhalo+1) - J
         
       DO I=mim,mip
        XNEW(I+nhalo) = HALPHA(I,J)
       ENDDO   
         
       DO K=1,KDIMN
        DO I=mim,mip
         YVAL(I+nhalo) = A(I,J,K)
        ENDDO                                                                           
        CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
        DO I=mim,mip
         HMV(I,K) = YNEW(I+nhalo) 
        ENDDO 
        DO I=mim,mip
         YVAL(I+nhalo) = A(I,NPO,K)
        ENDDO                                                                           
        CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
        DO I=mim,mip
         HPV(I,K) = YNEW(I+nhalo) 
        ENDDO             
       ENDDO    
       DO K=1,KDIMN 
        DO I=mim,mip
         A(I,J,K)   = HMV(I,K)
         A(I,NPO,K) = HPV(I,K)
        ENDDO
       ENDDO  
            
      ENDDO  ! J-loop   
      
!================================= 
      END SELECT
!=======================================      
      
      END SUBROUTINE halo_alpha     

!=======================================================================      
      SUBROUTINE HALO_BETA (mjm,mjp,nlen,KDIMN,MODE,CBETA,HBETA,JHALO_LOC,A)
!=======================================================================        
!     (ALONG BETA AXIS) Determine halo-values through 1-D interpolation.

      INTEGER, INTENT(IN) :: mjm,mjp,nlen  ! horizontal array sizes
      INTEGER, INTENT(IN) :: KDIMN         ! vertical array size  
      INTEGER, INTENT(IN) :: MODE          ! interpolation type (1: linear  2: cubic spline)        
   
      REAL (KIND=dbl_kind),DIMENSION(mjm:mjp),INTENT(IN) ::        &
        CBETA       ! beta coordinates of the regular q-grid points
      REAL (KIND=dbl_kind),DIMENSION(nhalo,mjm:mjp),INTENT(IN) ::  &
        HBETA       ! beta values of the target halo points 
      INTEGER,DIMENSION(nhalo,mjm:mjp),INTENT(IN) ::               &
        JHALO_LOC   ! j-index of CBETA that is just less than (or equal to) HBETA     
      
      ! target variable
      REAL (KIND=dbl_kind),DIMENSION(2*nhalo,mjm:mjp,KDIMN),INTENT(INOUT) :: A
      
      ! Local variables 
      REAL (KIND=dbl_kind),DIMENSION(nlen) :: XVAL,YVAL,XNEW,YNEW
      REAL (KIND=dbl_kind),DIMENSION(mjm:mjp,KDIMN) :: HMV,HPV
      INTEGER :: I,J,K,NPO,J1VAL,J2VAL
     
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
      DO I=1,nhalo
       NPO = (2*nhalo+1) - I
        
       DO J=mjm,mjp
        J1VAL = JHALO_LOC(I,J)
        J2VAL = J1VAL + 1
        DO K=1,KDIMN 
         HMV(J,K)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(I,J1VAL,K)   &
                                       ,CBETA(J2VAL),A(I,J2VAL,K))
         HPV(J,K)= FINTRP (1,HBETA(I,J),CBETA(J1VAL),A(NPO,J1VAL,K) &
                                       ,CBETA(J2VAL),A(NPO,J2VAL,K)) 
        ENDDO                               
       ENDDO   
       DO K=1,KDIMN 
        DO J=mjm,mjp
          A(I,J,K)   = HMV(J,K)
          A(NPO,J,K) = HPV(J,K)
        ENDDO
       ENDDO  
         
      ENDDO  ! I-loop    
     
!=================================      
!     Cubic spline interpolation
      CASE(2) 
!=================================                                                                                 
      DO I=1,nhalo
       NPO = (2*nhalo+1) - I
   
       DO J = mjm,mjp  
        XNEW(J+nhalo) = HBETA(I,J)
       ENDDO        
       DO K=1,KDIMN
        DO J=mjm,mjp
         YVAL(J+nhalo) = A(I,J,K)
        ENDDO                                                                  
        CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
        DO J=mjm,mjp
         HMV(J,K) = YNEW(J+nhalo) 
        ENDDO       
        DO J=mjm,mjp
         YVAL(J+nhalo) = A(NPO,J,K)
        ENDDO                                                                  
        CALL CUBICINT (NLEN,NLEN,XVAL,YVAL,XNEW,YNEW)
        DO J=mjm,mjp
         HPV(J,K) = YNEW(J+nhalo) 
        ENDDO                        
       ENDDO 
       DO K=1,KDIMN 
        DO J=mjm,mjp
          A(I,J,K)   = HMV(J,K)
          A(NPO,J,K) = HPV(J,K)
        ENDDO
       ENDDO  
         
      ENDDO  ! I-loop    

!================================= 
      END SELECT
!=================================      
      
      END SUBROUTINE halo_beta
      
END MODULE halo_vort