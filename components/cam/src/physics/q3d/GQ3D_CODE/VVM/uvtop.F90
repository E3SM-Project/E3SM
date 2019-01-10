MODULE uvtop
! Update horizontal wind at a certain layer 
! (currently, at the lowest layer) 

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk2
USE constld, only: dx,dy

! Subroutine being called
USE elliptic, only: relax_2d  

IMPLICIT NONE
PRIVATE

PUBLIC :: uvcal_3d  
   
CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------------
! SUBROUTINE uvcal_3d : Update horizontal wind at the lowest model layer
!----------------------------------------------------------------------------

!=======================================================================
   SUBROUTINE UVCAL_3D (channel,NITER, ittv)
!=======================================================================
      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER, INTENT(IN), OPTIONAL :: NITER      ! number of iteration for 2D elliptic equation, 
                                                  ! which is different from the basic setting 
                                                  ! "niterxy" (test purpose)    
                                                  
      INTEGER, INTENT(IN), OPTIONAL :: ittv    ! JUNG_DEBUG                                                

!     Local variables
      REAL (KIND=r8), DIMENSION(:,:), ALLOCATABLE :: UROT
      REAL (KIND=r8), DIMENSION(:,:), ALLOCATABLE :: VROT
      
      INTEGER :: I,J,num_seg,mim,mip,mim_a,mip_a,mjm,mjp,mjm_a,mjp_a
      
!     Solve the 2-D elliptic equations: only deviation fields

      ! Type 1: PSI
      IF (.not.PRESENT(NITER)) THEN
      
       if (present(ittv)) then
         write(90,*) 'UVTOP b4 psi'
         CALL RELAX_2D (1,channel,ittv=ittv)
         write(90,*) 'UVTOP a4 psi'
       else
         CALL RELAX_2D (1,channel)
       endif
      ELSE
       CALL RELAX_2D (1,channel,NITER2D=NITER)
      ENDIF
      
      ! Type 2: CHI
      IF (.not.PRESENT(NITER)) THEN
      
       if (present(ittv)) then
       write(91,*) 'UVTOP b4 chi'
       CALL RELAX_2D (2,channel,ittv=ittv)
       write(91,*) 'UVTOP a4 chi'
       else
       CALL RELAX_2D (2,channel)
       endif
      ELSE
       CALL RELAX_2D (2,channel,NITER2D=NITER)
      ENDIF

!=======================================================================   
      DO num_seg = 1, 4
!======================================================================= 
      mim   = channel%seg(num_seg)%mim    
      mip   = channel%seg(num_seg)%mip  
      mim_a = channel%seg(num_seg)%mim_a    
      mip_a = channel%seg(num_seg)%mip_a  
      
      mjm   = channel%seg(num_seg)%mjm   
      mjp   = channel%seg(num_seg)%mjp           
      mjm_a = channel%seg(num_seg)%mjm_a   
      mjp_a = channel%seg(num_seg)%mjp_a           
      
      ALLOCATE(UROT(mim:mip,mjm_a:mjp))
      ALLOCATE(VROT(mim_a:mip,mjm:mjp))
      
!     rotational part (deviation) 
      DO J = mjm_a, mjp  
       DO I = mim, mip
        UROT(I,J) = -(channel%seg(num_seg)%PSI(I,J)-channel%seg(num_seg)%PSI(I,J-1)) &
                    /(channel%seg(num_seg)%RG_U(I,J)*DY) 
       ENDDO
      ENDDO
      DO J = mjm, mjp  
       DO I = mim_a, mip
        VROT(I,J) =  (channel%seg(num_seg)%PSI(I,J)-channel%seg(num_seg)%PSI(I-1,J)) &
                    /(channel%seg(num_seg)%RG_V(I,J)*DX)      
       ENDDO
      ENDDO      

!-----------------------------------
!     1. Contravariant components
!-----------------------------------
      DO J = mjm_a, mjp_a  
       DO I = mim, mip_a
        channel%seg(num_seg)%U3DX(I,J,2) = channel%seg(num_seg)%U3DX_BG(I,J,2)             &
        +UROT(I,J)+channel%seg(num_seg)%GCONT_U(1,I,J)*(channel%seg(num_seg)%CHI(I+1,J)    &
                                                       -channel%seg(num_seg)%CHI(I,J))/DX  &
        +(channel%seg(num_seg)%GCONT_V(2,I,J)*(channel%seg(num_seg)%CHI(I,J+1)             &
                                              -channel%seg(num_seg)%CHI(I,J))              &
         +channel%seg(num_seg)%GCONT_V(2,I,J-1)*(channel%seg(num_seg)%CHI(I,J)             &
                                                -channel%seg(num_seg)%CHI(I,J-1))          &
         +channel%seg(num_seg)%GCONT_V(2,I+1,J)*(channel%seg(num_seg)%CHI(I+1,J+1)         &
                                                -channel%seg(num_seg)%CHI(I+1,J))          &
         +channel%seg(num_seg)%GCONT_V(2,I+1,J-1)*(channel%seg(num_seg)%CHI(I+1,J)         &
                                                  -channel%seg(num_seg)%CHI(I+1,J-1)))     &
         /(4.0_r8*DY)           
       ENDDO
      ENDDO
      
      DO J = mjm, mjp_a  
       DO I = mim_a, mip_a
        channel%seg(num_seg)%U3DY(I,J,2) = channel%seg(num_seg)%U3DY_BG(I,J,2)            &
        +VROT(I,J)+channel%seg(num_seg)%GCONT_V(3,I,J)*(channel%seg(num_seg)%CHI(I,J+1)   &
                                                       -channel%seg(num_seg)%CHI(I,J))/DY &
        +(channel%seg(num_seg)%GCONT_U(2,I,J)*(channel%seg(num_seg)%CHI(I+1,J)            &
                                              -channel%seg(num_seg)%CHI(I,J))             &
         +channel%seg(num_seg)%GCONT_U(2,I-1,J)*(channel%seg(num_seg)%CHI(I,J)            &
                                                -channel%seg(num_seg)%CHI(I-1,J))         &
         +channel%seg(num_seg)%GCONT_U(2,I,J+1)*(channel%seg(num_seg)%CHI(I+1,J+1)        &
                                                -channel%seg(num_seg)%CHI(I,J+1))         &
         +channel%seg(num_seg)%GCONT_U(2,I-1,J+1)*(channel%seg(num_seg)%CHI(I,J+1)        &
                                                -channel%seg(num_seg)%CHI(I-1,J+1)))      &
         /(4.0_r8*DX)            
       ENDDO
      ENDDO      

!-------------------------------------
!     2. Covariant components
!-------------------------------------     
      DO J = mjm_a, mjp
       DO I = mim_a, mip_a
        channel%seg(num_seg)%U3DX_CO(I,J,2) = channel%seg(num_seg)%U3DX_CO_BG(I,J,2)  &
                  +channel%seg(num_seg)%GG_U(3,I,J)*UROT(I,J)                         &
                  -(channel%seg(num_seg)%GG_V(2,I+1,J)*VROT(I+1,J)                    &
                   +channel%seg(num_seg)%GG_V(2,I,J)*VROT(I,J)                        &
                   +channel%seg(num_seg)%GG_V(2,I+1,J-1)*VROT(I+1,J-1)                &
                   +channel%seg(num_seg)%GG_V(2,I,J-1)*VROT(I,J-1))/4.0_r8            &
                  +(channel%seg(num_seg)%CHI(I+1,J)-channel%seg(num_seg)%CHI(I,J))/DX        
       ENDDO
      ENDDO          
  
      DO J = mjm_a, mjp_a  
       DO I = mim_a, mip
        channel%seg(num_seg)%U3DY_CO(I,J,2) = channel%seg(num_seg)%U3DY_CO_BG(I,J,2)  &
                  +channel%seg(num_seg)%GG_V(1,I,J)*VROT(I,J)                         &
                  -(channel%seg(num_seg)%GG_U(2,I,J+1)*UROT(I,J+1)                    &
                   +channel%seg(num_seg)%GG_U(2,I,J)*UROT(I,J)                        &
                   +channel%seg(num_seg)%GG_U(2,I-1,J+1)*UROT(I-1,J+1)                &
                   +channel%seg(num_seg)%GG_U(2,I-1,J)*UROT(I-1,J))/4.0_r8            &
                  +(channel%seg(num_seg)%CHI(I,J+1)-channel%seg(num_seg)%CHI(I,J))/DY                 
       ENDDO
      ENDDO      
      
      DEALLOCATE(UROT)
      DEALLOCATE(VROT)
      
!===============================  
      ENDDO  ! num_seg 
!=============================== 
      
   END SUBROUTINE uvcal_3d   

END MODULE uvtop
