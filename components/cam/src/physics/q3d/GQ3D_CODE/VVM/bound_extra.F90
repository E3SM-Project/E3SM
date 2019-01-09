module bound_extra

! this module contains routines that specify the boundary values in the normal and  
! vertical directions of the channel domain.

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t 

USE parmsld,        only: ntracer,nk1,nk2,nk3 

IMPLICIT NONE
PRIVATE

PUBLIC ::  bound_normal,bound_vert,bound_sync    
           
CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------------------------
! SUBROUTINE bound_normal : specify the boundary values in the normal direction (Q3D algorithm)
! SUBROUTINE bound_vert   : specify the boundary values in the vertical direction 
! SUBROUTINE HALO_NS      : called from bound_normal (for X-channel) 
! SUBROUTINE HALO_NS_2    : called from bound_normal (for X-channel) 
! SUBROUTINE HALO_EW      : called from bound_normal (for Y-channel)
! SUBROUTINE HALO_EW_2    : called from bound_normal (for Y-channel)
! SUBROUTINE HALO_TB1     : called from bound_vert
! SUBROUTINE HALO_TB2     : called from bound_vert
!-----------------------------------------------------------------------------------------

!=========================================================================================
      SUBROUTINE BOUND_NORMAL (num_halo,channel,TH3D,QV3D,QT3D,  &              
                               QC3D,QI3D,QR3D,QS3D,QG3D,         &
                               W3D,Z3DX,Z3DY,Z3DZ,               &
                               PSI,CHI,TOPOZ,ZROUGH,GWET,TG)              
!=========================================================================================                               
!     Specify the ghost values in the normal direction to the channel.
!     Number of halo points to be specified can be different (NUM_HALO),
!     but, the horizontal array size of the variables is fixed, i.e., (mim:mip,mjm:mjp)
   
      INTEGER, INTENT(IN) :: num_halo
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D,QT3D, &
                                       W3D,Z3DX,Z3DY,Z3DZ, &
                                       PSI,CHI,TOPOZ,ZROUGH,GWET,TG
              
      INTEGER :: num_seg,mi1,mim,mip,mj1,mjm,mjp,nt 

!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip   
      
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp    
      
!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%TH3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%TH3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%TH3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%TH3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QV3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QV3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QV3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QV3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QC3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QC3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QC3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QC3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QI3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QI3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QI3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QI3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QR3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QR3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QR3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QR3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QS3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QS3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QS3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QS3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QG3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QG3D_bg(:,:,2:NK2))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QG3D(:,:,2:NK2),       &
                        channel%seg(num_seg)%QG3D_bg(:,:,2:NK2))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QT3D(:,:,2:NK2,nt),    &
                        channel%seg(num_seg)%QT3D_bg(:,:,2:NK2,nt))
          ENDDO              
        ELSE
          ! y-array
          DO nt = 1,ntracer
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1,       &
                        channel%seg(num_seg)%QT3D(:,:,2:NK2,nt),    &
                        channel%seg(num_seg)%QT3D_bg(:,:,2:NK2,nt))
          ENDDO              
        ENDIF                      
        endif                        
      ENDIF
!---------------------------------
! Z3DX
!---------------------------------
      IF (PRESENT(Z3DX)) THEN
        if (Z3DX) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,     &
                        channel%seg(num_seg)%Z3DX(:,:,2:NK1),       &
                        channel%seg(num_seg)%Z3DX_bg(:,:,2:NK1))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,     &
                        channel%seg(num_seg)%Z3DX(:,:,2:NK1),       &
                        channel%seg(num_seg)%Z3DX_bg(:,:,2:NK1))
        ENDIF                      
        endif                        
      ENDIF 
!---------------------------------
! Z3DY
!---------------------------------
      IF (PRESENT(Z3DY)) THEN
        if (Z3DY) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,     &
                        channel%seg(num_seg)%Z3DY(:,:,2:NK1),       &
                        channel%seg(num_seg)%Z3DY_bg(:,:,2:NK1))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,     &
                        channel%seg(num_seg)%Z3DY(:,:,2:NK1),       &
                        channel%seg(num_seg)%Z3DY_bg(:,:,2:NK1))
        ENDIF                      
        endif                        
      ENDIF       
!---------------------------------
! Z3DZ
!---------------------------------
      IF (PRESENT(Z3DZ)) THEN
        if (Z3DZ) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,1,       &
                        channel%seg(num_seg)%Z3DZ(:,:,NK2),       &
                        channel%seg(num_seg)%Z3DZ_bg(:,:,1))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,1,       &
                        channel%seg(num_seg)%Z3DZ(:,:,NK2),       &
                        channel%seg(num_seg)%Z3DZ_bg(:,:,1))
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
          CALL HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,    &
                        channel%seg(num_seg)%W3D(:,:,2:NK1),       &
                        channel%seg(num_seg)%W3D_bg(:,:,2:NK1))
        ELSE
          ! y-array
          CALL HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,NK1-1,    &
                        channel%seg(num_seg)%W3D(:,:,2:NK1),       &
                        channel%seg(num_seg)%W3D_bg(:,:,2:NK1))
        ENDIF                      
        endif                        
      ENDIF

! 2-D fields with no background      
!---------------------------------
! PSI
!---------------------------------
      IF (PRESENT(PSI)) THEN
        if (PSI) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%PSI)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%PSI)
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
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%CHI)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%CHI)
        ENDIF                      
        endif                        
      ENDIF  
!---------------------------------
! TOPOZ
!---------------------------------
      IF (PRESENT(TOPOZ)) THEN
        if (TOPOZ) then
        IF (mi1.GT.mj1) THEN
          ! x-array        
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%TOPOZ)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%TOPOZ)
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
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%ZROUGH)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%ZROUGH)
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
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%GWET)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%GWET)
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
          CALL HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%TG)
        ELSE
          ! y-array
          CALL HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,channel%seg(num_seg)%TG)
        ENDIF                      
        endif                        
      ENDIF  
      
!********************************   
      ENDDO  ! num_seg 
!******************************** 

      end subroutine bound_normal

!=========================================================================================      
      SUBROUTINE BOUND_VERT (channel,TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D,QT3D, &
                             W3D,Z3DX,Z3DY)
!=========================================================================================                             
!     Specify the boundary values at the model top and bottom 
     
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      
      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,QC3D,QI3D,QR3D,QS3D,QG3D,QT3D, &
                                       W3D,Z3DX,Z3DY
              
      INTEGER :: num_seg,mim,mip,mjm,mjp,nt 

!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mim = channel%seg(num_seg)%mim
      mip = channel%seg(num_seg)%mip   
      
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp 
      
!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%TH3D)     
        endif                        
      ENDIF
!---------------------------------
! QV3D
!---------------------------------
      IF (PRESENT(QV3D)) THEN
        if (QV3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QV3D)     
        endif                        
      ENDIF
!---------------------------------
! QC3D
!---------------------------------
      IF (PRESENT(QC3D)) THEN
        if (QC3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QC3D)     
        endif                        
      ENDIF
!---------------------------------
! QI3D
!---------------------------------
      IF (PRESENT(QI3D)) THEN
        if (QI3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QI3D)     
        endif                        
      ENDIF
!---------------------------------
! QR3D
!---------------------------------
      IF (PRESENT(QR3D)) THEN
        if (QR3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QR3D)     
        endif                        
      ENDIF
!---------------------------------
! QS3D
!---------------------------------
      IF (PRESENT(QS3D)) THEN
        if (QS3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QS3D)     
        endif                        
      ENDIF
!---------------------------------
! QG3D
!---------------------------------
      IF (PRESENT(QG3D)) THEN
        if (QG3D) then
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QG3D)     
        endif                        
      ENDIF
!---------------------------------
! QT3D
!---------------------------------
      IF (PRESENT(QT3D)) THEN
        if (QT3D) then
          DO nt = 1,ntracer
          CALL HALO_TB1 (mim,mip,mjm,mjp,nk3,channel%seg(num_seg)%QT3D(:,:,:,nt))    
          ENDDO 
        endif                        
      ENDIF    
                                            
!---------------------------------
! Z3DX
!---------------------------------
      IF (PRESENT(Z3DX)) THEN
        if (Z3DX) then
          CALL HALO_TB2 (mim,mip,mjm,mjp,nk2,channel%seg(num_seg)%Z3DX)     
        endif                        
      ENDIF
!---------------------------------
! Z3DY
!---------------------------------
      IF (PRESENT(Z3DY)) THEN
        if (Z3DY) then
          CALL HALO_TB2 (mim,mip,mjm,mjp,nk2,channel%seg(num_seg)%Z3DY)     
        endif                        
      ENDIF      
!---------------------------------
! W3D
!---------------------------------
      IF (PRESENT(W3D)) THEN
        if (W3D) then
          CALL HALO_TB2 (mim,mip,mjm,mjp,nk2,channel%seg(num_seg)%W3D)     
        endif                        
      ENDIF

!********************************   
      ENDDO  ! num_seg 
!******************************** 
       
      end subroutine bound_vert      

!=========================================================================================
      SUBROUTINE BOUND_SYNC (channel,Z3DZ,PSI)              
!=========================================================================================                               
!     Specify the ghost values in the normal direction to the channel.
!     Number of halo points to be specified can be different (NUM_HALO),
!     but, the horizontal array size of the variables is fixed, i.e., (mim:mip,mjm:mjp)
   
      type(channel_t), INTENT(INOUT) :: channel   ! Channel data
      LOGICAL, INTENT(IN), OPTIONAL  :: Z3DZ,PSI

      ! Local
      INTEGER :: mi1_2,mim_2,mip_2,mjm_2,mjp_2
      INTEGER :: mj1_3,mim_3,mip_3,mjm_3,mjp_3,mi1_3
      INTEGER :: mj1_4,mim_4,mip_4,mjm_4,mjp_4

      ! Size of segment 2
      mi1_2 = channel%seg(2)%mi1 
      
      mim_2 = channel%seg(2)%mim
      mip_2 = channel%seg(2)%mip   
      mjm_2 = channel%seg(2)%mjm    
      mjp_2 = channel%seg(2)%mjp    

      ! Size of segment 3 
      mi1_3 = channel%seg(3)%mi1  
      mj1_3 = channel%seg(3)%mj1 
       
      mim_3 = channel%seg(3)%mim
      mip_3 = channel%seg(3)%mip   
      mjm_3 = channel%seg(3)%mjm    
      mjp_3 = channel%seg(3)%mjp    

      ! Size of segment 4
      mj1_4 = channel%seg(4)%mj1 
        
      mim_4 = channel%seg(4)%mim
      mip_4 = channel%seg(4)%mip   
      mjm_4 = channel%seg(4)%mjm    
      mjp_4 = channel%seg(4)%mjp    

!---------------------------------
! Z3DZ
!---------------------------------
      IF (PRESENT(Z3DZ)) THEN
        if (Z3DZ) then
        SELECT CASE (channel%num_chg)
        CASE(2)
          CALL SYNC_GROUP2 (mj1_3,mim_3,mip_3,mjm_3,mjp_3,         &
                            mj1_4,mim_4,mip_4,mjm_4,mjp_4,mi1_3,1, & 
                            channel%seg(3)%Z3DZ(:,:,NK2),channel%seg(4)%Z3DZ(:,:,NK2))

        CASE(3)
          CALL SYNC_GROUP3 (mi1_2,mim_2,mip_2,mjm_2,mjp_2,    &
                            mj1_3,mim_3,mip_3,mjm_3,mjp_3,1,  &
                            channel%seg(2)%Z3DZ(:,:,NK2),channel%seg(3)%Z3DZ(:,:,NK2))
        END SELECT
        endif                        
      ENDIF 

! 2-D fields with no background      
!---------------------------------
! PSI
!---------------------------------
      IF (PRESENT(PSI)) THEN
        if (PSI) then
        SELECT CASE (channel%num_chg)
        CASE(2)
          CALL SYNC_GROUP2_2 (mj1_3,mim_3,mip_3,mjm_3,mjp_3,       &
                              mj1_4,mim_4,mip_4,mjm_4,mjp_4,mi1_3, &
                              channel%seg(3)%PSI,channel%seg(4)%PSI)

        CASE(3)
          CALL SYNC_GROUP3_2 (mi1_2,mim_2,mip_2,mjm_2,mjp_2,    &
                              mj1_3,mim_3,mip_3,mjm_3,mjp_3,    &
                              channel%seg(2)%PSI,channel%seg(3)%PSI)
        END SELECT
        endif                        
      ENDIF  

      CONTAINS
      
      SUBROUTINE SYNC_GROUP2 (mj1_3,mim_3,mip_3,mjm_3,mjp_3, &
                              mj1_4,mim_4,mip_4,mjm_4,mjp_4, &
                              mi1_3,ksize,a3,a4)
!=========================================================================================      
!     Synchronize the edge values of F5 (Seg3) and F3 (Seg4) in GROUP 2 channels 
      
      INTEGER, INTENT(IN) :: mj1_3,mim_3,mip_3,mjm_3,mjp_3, &
                             mj1_4,mim_4,mip_4,mjm_4,mjp_4,mi1_3,ksize
      REAL (kind=r8), DIMENSION(mim_3:mip_3,mjm_3:mjp_3,ksize), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_4:mip_4,mjm_4:mjp_4,ksize), INTENT(INOUT) :: A4
      
      ! Local
      REAL (kind=r8) :: TEMP
      INTEGER :: chn,chn_rev,k               

!---------------------------- 
      ! Seg 3 & Seg 4
!----------------------------
        DO K = 1, ksize 
         DO chn = mim_3, mip_3-1
          chn_rev = mi1_3 - chn 
          TEMP = (A3(chn,mj1_3,K)+A4(chn_rev,mj1_4,K))/2.0_r8
          A3(chn,mj1_3,K)     = TEMP
          A4(chn_rev,mj1_4,K) = TEMP 
         ENDDO
        ENDDO   

      end subroutine sync_group2 
      
      SUBROUTINE SYNC_GROUP2_2 (mj1_3,mim_3,mip_3,mjm_3,mjp_3, &
                                mj1_4,mim_4,mip_4,mjm_4,mjp_4, &
                                mi1_3,a3,a4)
!=========================================================================================      
!     Synchronize the edge values of F5 (Seg3) and F3 (Seg4) in GROUP 2 channels 
      
      INTEGER, INTENT(IN) :: mj1_3,mim_3,mip_3,mjm_3,mjp_3, &
                             mj1_4,mim_4,mip_4,mjm_4,mjp_4,mi1_3
      REAL (kind=r8), DIMENSION(mim_3:mip_3,mjm_3:mjp_3), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_4:mip_4,mjm_4:mjp_4), INTENT(INOUT) :: A4
      
      ! Local
      REAL (kind=r8) :: TEMP
      INTEGER :: chn,chn_rev               

!---------------------------- 
      ! Seg 3 & Seg 4
!----------------------------
         DO chn = mim_3, mip_3-1
          chn_rev = mi1_3 - chn 
          TEMP = (A3(chn,mj1_3)+A4(chn_rev,mj1_4))/2.0_r8
          A3(chn,mj1_3)     = TEMP
          A4(chn_rev,mj1_4) = TEMP 
         ENDDO  

      end subroutine sync_group2_2      
      
      SUBROUTINE SYNC_GROUP3 (mi1_2,mim_2,mip_2,mjm_2,mjp_2, &
                              mj1_3,mim_3,mip_3,mjm_3,mjp_3,ksize,a2,a3)
!=========================================================================================      
!     Synchronize the edge values of F5 (Seg2) and F2 (Seg3) in GROUP 3 channels 
      
      INTEGER, INTENT(IN) :: mi1_2,mim_2,mip_2,mjm_2,mjp_2, &
                             mj1_3,mim_3,mip_3,mjm_3,mjp_3,ksize       
      REAL (kind=r8), DIMENSION(mim_2:mip_2,mjm_2:mjp_2,ksize), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_3:mip_3,mjm_3:mjp_3,ksize), INTENT(INOUT) :: A3
      
      ! Local
      REAL (kind=r8) :: TEMP
      INTEGER :: chn,k               

!---------------------------- 
      ! Seg 2 & Seg 3
!----------------------------
        DO K = 1, ksize 
         DO chn = mjm_2, mjp_2
          TEMP = (A2(mi1_2,chn,K)+A3(chn,mj1_3,K))/2.0_r8
          
          A2(mi1_2,chn,K) = TEMP
          A3(chn,mj1_3,K) = TEMP
         ENDDO
        ENDDO   
             
      end subroutine sync_group3  
      
      SUBROUTINE SYNC_GROUP3_2 (mi1_2,mim_2,mip_2,mjm_2,mjp_2, &
                                mj1_3,mim_3,mip_3,mjm_3,mjp_3,a2,a3)
!=========================================================================================      
!     Synchronize the edge values of F5 (Seg2) and F2 (Seg3) in GROUP 3 channels 
      
      INTEGER, INTENT(IN) :: mi1_2,mim_2,mip_2,mjm_2,mjp_2, &
                             mj1_3,mim_3,mip_3,mjm_3,mjp_3       
      REAL (kind=r8), DIMENSION(mim_2:mip_2,mjm_2:mjp_2), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_3:mip_3,mjm_3:mjp_3), INTENT(INOUT) :: A3
      
      ! Local
      REAL (kind=r8) :: TEMP
      INTEGER :: chn               

!---------------------------- 
      ! Seg 2 & Seg 3
!----------------------------
         DO chn = mjm_2, mjp_2
          TEMP = (A2(mi1_2,chn)+A3(chn,mj1_3))/2.0_r8
          
          A2(mi1_2,chn) = TEMP
          A3(chn,mj1_3) = TEMP
         ENDDO
             
      end subroutine sync_group3_2      
      
      end subroutine bound_sync        
            
!=========================================================================================      
      SUBROUTINE HALO_NS (mi1,mim,mip,mj1,mjm,mjp,num_halo,ksize,a,abg)
!=========================================================================================      
!     specify the halo values at the northern and southern sides of the x-channel 
!     (deviation is uniform across the channel)
      
      INTEGER, INTENT(IN) :: mi1,mim,mip,mj1,mjm,mjp,num_halo,ksize       
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(IN) :: Abg
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(INOUT) :: A
      
      INTEGER :: I,J,K,JM,JP       

         DO J = 1,num_halo
          JM = 1 - J
          JP = mj1 + J
          DO K = 1, ksize 
           DO I = 1, mi1
            A(I,JM,K) = Abg(I,JM,K) + (A(I,1,K) - Abg(I,1,K))
            A(I,JP,K) = Abg(I,JP,K) + (A(I,mj1,K) - Abg(I,mj1,K))
           ENDDO
          ENDDO  
         ENDDO 
                    
      end subroutine halo_ns
      
!=========================================================================================      
      SUBROUTINE HALO_NS_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,a)
!=========================================================================================      
!     Same as HALO_NS, but for a 2D variable (whole value is uniform across the channel)
      
      INTEGER, INTENT(IN) :: mi1,mim,mip,mj1,mjm,mjp,num_halo       
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp), INTENT(INOUT) :: A
      
      INTEGER :: I,J,JM,JP       

         DO J = 1,num_halo
          JM = 1 - J
          JP = mj1 + J
           DO I = 1, mi1
            A(I,JM) = A(I,1)
            A(I,JP) = A(I,mj1)
           ENDDO 
         ENDDO 
                    
      end subroutine halo_ns_2      

!=========================================================================================      
      SUBROUTINE HALO_EW (mi1,mim,mip,mj1,mjm,mjp,num_halo,ksize,a,abg)
!=========================================================================================      
!     specify the halo values at the eastern and western sides of the y-channel 
!     (deviation is uniform across the channel)
      
      INTEGER, INTENT(IN) :: mi1,mim,mip,mj1,mjm,mjp,num_halo,ksize     
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(IN) :: Abg
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(INOUT) :: A
      
      INTEGER :: I,J,K,IM,IP    

         DO I = 1,NUM_HALO
          IM = 1 - I
          IP = mi1 + I
          DO K = 1, ksize 
           DO J = 1, mj1
            A(IM,J,K) = Abg(IM,J,K) + (A(1,J,K) - Abg(1,J,K)) 
            A(IP,J,K) = Abg(IP,J,K) + (A(mi1,J,K) - Abg(mi1,J,K))
           ENDDO
          ENDDO
         ENDDO  
            
      end subroutine halo_ew
      
!=========================================================================================      
      SUBROUTINE HALO_EW_2 (mi1,mim,mip,mj1,mjm,mjp,num_halo,a)
!=========================================================================================      
!     Same as HALO_EW, but for a 2D variable (whole value is uniform across the channel)
      
      INTEGER, INTENT(IN) :: mi1,mim,mip,mj1,mjm,mjp,num_halo     
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp), INTENT(INOUT) :: A
      
      INTEGER :: I,J,IM,IP    

         DO I = 1,NUM_HALO
          IM = 1 - I
          IP = mi1 + I
           DO J = 1, mj1
            A(IM,J) = A(1,J)  
            A(IP,J) = A(mi1,J) 
           ENDDO
         ENDDO  
            
      end subroutine halo_ew_2      

!=========================================================================================
      SUBROUTINE HALO_TB1 (mim,mip,mjm,mjp,ksize,A)
!=========================================================================================      
!     specify the boundary values at the model top and bottom 
      
      INTEGER, INTENT(IN) :: mim,mip,mjm,mjp,ksize       
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(INOUT) :: A
      
      INTEGER :: I,J    
      
      DO J = mjm,mjp
       DO I = mim,mip
        A(I,J,1)     = A(I,J,2)
        A(I,J,KSIZE) = A(I,J,KSIZE-1)
       ENDDO
      ENDDO
      
      end subroutine halo_tb1

!=========================================================================================
      SUBROUTINE HALO_TB2 (mim,mip,mjm,mjp,ksize,A)
!=========================================================================================      
!     specify the boundary values at the model top and bottom 
      
      INTEGER, INTENT(IN) :: mim,mip,mjm,mjp,ksize       
      REAL (kind=r8), DIMENSION(mim:mip,mjm:mjp,ksize), INTENT(INOUT) :: A
      
      INTEGER :: I,J    
      
      DO J = mjm,mjp
       DO I = mim,mip
        A(I,J,1)     = 0.0_r8
        A(I,J,KSIZE) = 0.0_r8
       ENDDO
      ENDDO
      
      end subroutine halo_tb2      
            
end module bound_extra
            