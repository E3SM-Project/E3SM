MODULE q3d_eddy_module

      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
      USE vvm_data_types, only: channel_t
      
      USE parmsld, only: ntracer,nk2
      USE constld, only: physics
      
IMPLICIT NONE
PRIVATE

PUBLIC :: eddy_prepare

CONTAINS

!=======================================================================
   SUBROUTINE EDDY_PREPARE(channel)
!=======================================================================
!  Prepare eddy fields to calculate the eddy transport effects
!  U3DX_ED & U3DY_ED (0:mi1,0:mj1)   Others (1:mi1,1:mj1)

      type(channel_t), intent(inout) :: channel   ! channel data
   
      integer mi1,mim,mip,mj1,mjm,mjp
      integer num_seg,nt

#ifdef FIXTHIS
!=======================================================================   
      DO num_seg = 1, 4
!=======================================================================
      mi1 = channel%seg(num_seg)%mi1     ! x-size of channel segment (temp setting)
      mj1 = channel%seg(num_seg)%mj1     ! y-size of channel segment (temp setting)

      mim = channel%seg(num_seg)%mim 
      mip = channel%seg(num_seg)%mip 
      mjm = channel%seg(num_seg)%mjm
      mjp = channel%seg(num_seg)%mjp

      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%TH3D(:,:,1:nk2), &
                     channel%seg(num_seg)%TH3D_BG,         &
                     channel%seg(num_seg)%TH3D_ED)
                                                   
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QV3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QV3D_BG,         &
                     channel%seg(num_seg)%QV3D_ED)

      IF (PHYSICS) THEN 
        
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QC3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QC3D_BG,         &
                     channel%seg(num_seg)%QC3D_ED)
                                                   
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QI3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QI3D_BG,         &
                     channel%seg(num_seg)%QI3D_ED)

      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QR3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QR3D_BG,         &
                     channel%seg(num_seg)%QR3D_ED)
                                                   
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QS3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QS3D_BG,         &
                     channel%seg(num_seg)%QS3D_ED)

      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%QG3D(:,:,1:nk2), &
                     channel%seg(num_seg)%QG3D_BG,         &
                     channel%seg(num_seg)%QG3D_ED)
                                                   
      ENDIF ! PHYSICS                                             
             
      DO nt = 1,ntracer                                             
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,               &
                     channel%seg(num_seg)%QT3D(:,:,1:nk2,nt), &
                     channel%seg(num_seg)%QT3D_BG(:,:,:,nt),  &
                     channel%seg(num_seg)%QT3D_ED(:,:,:,nt))
      ENDDO               
                                                   
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,0,            &
                     channel%seg(num_seg)%U3DX(:,:,1:nk2), &
                     channel%seg(num_seg)%U3DX_BG,         &
                     channel%seg(num_seg)%U3DX_ED)           
                     
      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,0,            &
                     channel%seg(num_seg)%U3DY(:,:,1:nk2), &
                     channel%seg(num_seg)%U3DY_BG,         &
                     channel%seg(num_seg)%U3DY_ED)                                                                        

      CALL GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,1,            &
                     channel%seg(num_seg)%W3D,             &
                     channel%seg(num_seg)%W3D_BG,          &
                     channel%seg(num_seg)%W3D_ED)

!=======================================================================   
      ENDDO   ! num_seg 
!=======================================================================
      
    CONTAINS

    SUBROUTINE GET_EDDY (mi1,mim,mip,mj1,mjm,mjp,nbeg,A,A_BG,A_ED)
    
    INTEGER, INTENT(IN) :: mi1,mim,mip,mj1,mjm,mjp,nbeg
    REAL (KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,nk2),INTENT(IN) :: A
    REAL (KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,nk2),INTENT(IN) :: A_BG
    REAL (KIND=dbl_kind),DIMENSION(0:mi1,0:mj1,nk2),INTENT(OUT)    :: A_ED
    
    INTEGER I,J,K
    
    DO K = 1, nk2
     DO J = nbeg, mj1
      DO I = nbeg, mi1
       A_ED(I,J,K) = A(I,J,K) - A_BG(I,J,K)
      ENDDO
     ENDDO
    ENDDO   
    
    END SUBROUTINE get_eddy

#endif

   END SUBROUTINE eddy_prepare
  
END MODULE q3d_eddy_module