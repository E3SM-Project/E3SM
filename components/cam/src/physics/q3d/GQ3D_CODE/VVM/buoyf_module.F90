MODULE buoyf_module
! Calculates buoyancy terms in vorticity equations

      USE shr_kind_mod,   only: r8 => shr_kind_r8
      USE vvm_data_types, only: channel_t

      USE parmsld, only: nk1,nk2
      USE constld, only: thbar,dx,dy,grav,delta,physics

IMPLICIT NONE
PRIVATE

PUBLIC :: buoyf_3d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE buoyf_3d: Calculates buoyancy terms in vorticity equations
!-----------------------------------------------------------------------

!=======================================================================
   SUBROUTINE BUOYF_3D ( channel )
!=======================================================================
!  Temporary use of channel%seg(num_seg)%TERM1

      type(channel_t), intent(inout) :: channel   ! channel data
     
!     Local variables 
      INTEGER :: I, J, K, KLOWU, KLOWV
      INTEGER :: num_seg,mi1,mj1,mim_c,mip_c,mjm_c,mjp_c

!*******************************************************************   
      DO num_seg = 1, 4
!*******************************************************************  
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment 
      mim_c  = channel%seg(num_seg)%mim_c     
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment  
      mjm_c  = channel%seg(num_seg)%mjm_c    
      mjp_c  = channel%seg(num_seg)%mjp_c         

      IF (PHYSICS) THEN
        ! temporary use of term1
        DO K = 2, nk2
         DO J = mjm_c, mjp_c
          DO I = mim_c, mip_c
           channel%seg(num_seg)%TERM1(I,J,K) = channel%seg(num_seg)%QC3D(I,J,K) &
                                             + channel%seg(num_seg)%QR3D(I,J,K) &
                                             + channel%seg(num_seg)%QI3D(I,J,K) &
                                             + channel%seg(num_seg)%QS3D(I,J,K) &
                                             + channel%seg(num_seg)%QG3D(I,J,K)
          ENDDO
         ENDDO
        ENDDO
      ELSE
        channel%seg(num_seg)%TERM1(:,:,:) = 0.0_r8
      ENDIF

      DO J = 1, mj1
       DO I = 1, mi1
        klowu = channel%seg(num_seg)%KLOWU_IJ(I,J)
        klowv = channel%seg(num_seg)%KLOWV_IJ(I,J)

        DO K = klowv, nk1
         channel%seg(num_seg)%FZXBU(I,J,K) = &
                             (channel%seg(num_seg)%TH3D(I,J+1,K)              &
                             -channel%seg(num_seg)%TH3D(I,J,K))/THBAR(K)      &
                            +(channel%seg(num_seg)%TH3D(I,J+1,K+1)            &
                             -channel%seg(num_seg)%TH3D(I,J,K+1))/THBAR(K+1)  &
                      +DELTA*(channel%seg(num_seg)%QV3D(I,J+1,K)              &
                             -channel%seg(num_seg)%QV3D(I,J,K)                &
                             +channel%seg(num_seg)%QV3D(I,J+1,K+1)            &
                             -channel%seg(num_seg)%QV3D(I,J,K+1))             &
            -(channel%seg(num_seg)%TERM1(I,J+1,K)-channel%seg(num_seg)%TERM1(I,J,K) &
             +channel%seg(num_seg)%TERM1(I,J+1,K+1)-channel%seg(num_seg)%TERM1(I,J,K+1))
        ENDDO
        DO K = klowv, nk1
         channel%seg(num_seg)%FZXBU(I,J,K) = &
         channel%seg(num_seg)%FZXBU(I,J,K)*GRAV/(2.0_r8*DY*channel%seg(num_seg)%RG_V(I,J))
        ENDDO

        DO K = klowu, nk1
         channel%seg(num_seg)%FZYBU(I,J,K) = &
                             (channel%seg(num_seg)%TH3D(I+1,J,K)              &
                             -channel%seg(num_seg)%TH3D(I,J,K))/THBAR(K)      &
                            +(channel%seg(num_seg)%TH3D(I+1,J,K+1)            &
                             -channel%seg(num_seg)%TH3D(I,J,K+1))/THBAR(K+1)  &
                      +DELTA*(channel%seg(num_seg)%QV3D(I+1,J,K)              &
                             -channel%seg(num_seg)%QV3D(I,J,K)                &
                             +channel%seg(num_seg)%QV3D(I+1,J,K+1)            &
                             -channel%seg(num_seg)%QV3D(I,J,K+1))             &
            -(channel%seg(num_seg)%TERM1(I+1,J,K)-channel%seg(num_seg)%TERM1(I,J,K) &
             +channel%seg(num_seg)%TERM1(I+1,J,K+1)-channel%seg(num_seg)%TERM1(I,J,K+1))
        ENDDO

        DO K = klowu, nk1
         channel%seg(num_seg)%FZYBU(I,J,K) = &
        -channel%seg(num_seg)%FZYBU(I,J,K)*GRAV/(2.0_r8*DX*channel%seg(num_seg)%RG_U(I,J))
        ENDDO

       ENDDO  ! I-loop
      ENDDO  ! J-loop

!****************************  
      ENDDO  ! num_seg
!****************************

   END SUBROUTINE buoyf_3d

END MODULE buoyf_module
