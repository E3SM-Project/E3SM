MODULE wind_module
! diagnose wind components

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: GQ3D_DATA_T

USE parmsld, only: nk1,nk2,nk3
USE constld, only: d0_0,d4_0,dx,dy,dz,fnz

! Subroutines being called
USE elliptic, only: relax_3d
USE uvtop,    only: uvtop_3d

IMPLICIT NONE
PRIVATE

PUBLIC :: wind_3d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE wind_3d : diagnose wind components
!-----------------------------------------------------------------------
!    D0_0  = 0.0_dbl_kind
!    D4_0  = 4.0_dbl_kind

!=======================================================================
   SUBROUTINE WIND_3D (N1,N2,UCHANGE,WCHANGE,channel,NITER2D)
!=======================================================================   
!  Diagnose three components of wind
!
!        w3d, z3dx, z3dy are good for nhalo
!        u3dx and u3dy are good for nhalo_adv
!
!        nhalo_cal = 1       (fixed)             mim_c:mip_c
!        nhalo_adv = 2       (might be changed)  mim_a:mip_a
!        nhalo = nhalo_adv+1 (fixed)             mim:mip
!        (mim_ce = mim_c-1; mip_ce = mip_c+1)
!
!  u and v are calculated for different sizes because they are 
!  diagnostically calculated from a limited dataset.
!
!  PLAN: modify the wind calculation (Q3D: close the dynamics)
!------------------------------------------------------------------------ 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER, INTENT(IN) :: N1, N2
      LOGICAL, INTENT(IN) :: UCHANGE, WCHANGE

      INTEGER, INTENT(IN), OPTIONAL :: NITER2D   
      ! number of iteration for 2D elliptic equation, 
      ! which is different from the basic setting "niterxy" (test purpose) 

!     Local Variables
      REAL (KIND=dbl_kind), DIMENSION(nk1) :: ADD_U,ADD_V

      INTEGER :: I,J,K,KLOWU,KLOWV,num_seg
      INTEGER :: mim,mim_a,mim_ce,mip_a,mip_c,mip_ce  
      INTEGER :: mjm,mjm_a,mjm_ce,mjp_a,mjp_c,mjp_ce 
!------------------------------------------------------------------------      

      IF (WCHANGE) CALL RELAX_3D (channel)

      IF (UCHANGE) THEN

      ! UPDATE wind at the model top
      IF (.not.PRESENT(NITER2D)) THEN
        CALL UVTOP_3D (channel)
      ELSE
        CALL UVTOP_3D (channel,NITER=NITER2D)
      ENDIF

      ! UPDATE wind below the model top
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mim    = channel%seg(num_seg)%mim
      mim_a  = channel%seg(num_seg)%mim_a 
      mim_ce = channel%seg(num_seg)%mim_ce
      
      mip_a  = channel%seg(num_seg)%mip_a
      mip_c  = channel%seg(num_seg)%mip_c
      mip_ce = channel%seg(num_seg)%mip_ce    
      
      mjm    = channel%seg(num_seg)%mjm 
      mjm_a  = channel%seg(num_seg)%mjm_a   
      mjm_ce = channel%seg(num_seg)%mjm_ce   
         
      mjp_a  = channel%seg(num_seg)%mjp_a 
      mjp_c  = channel%seg(num_seg)%mjp_c
      mjp_ce = channel%seg(num_seg)%mjp_ce     
      
!----------------------------------------------------------
!     1. Contravariant components of horizontal wind
!     (calculated with contravariant vorticity components)
!----------------------------------------------------------
!---------------
! U-component
!---------------

      DO J = mjm_a, mjp_a
      DO I = mim, mip_a
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = KLOWU, nk1
       ADD_U(K) = &
       (channel%seg(num_seg)%GCONT_V(2,I+1,J)*(channel%seg(num_seg)%W3D(I+1,J+1,K)    &
                                              -channel%seg(num_seg)%W3D(I+1,J,K))     &
       +channel%seg(num_seg)%GCONT_V(2,I+1,J-1)*(channel%seg(num_seg)%W3D(I+1,J,K)    &
                                                -channel%seg(num_seg)%W3D(I+1,J-1,K)) &
       +channel%seg(num_seg)%GCONT_V(2,I,J)*(channel%seg(num_seg)%W3D(I,J+1,K)        &
                                            -channel%seg(num_seg)%W3D(I,J,K))         &
       +channel%seg(num_seg)%GCONT_V(2,I,J-1)*(channel%seg(num_seg)%W3D(I,J,K)        &
                                              -channel%seg(num_seg)%W3D(I,J-1,K)))    &
       /(4.dbl_kind*DY)                                                               &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(D4_0*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = nk1, KLOWU, -1
       channel%seg(num_seg)%U3DX(I,J,K) = channel%seg(num_seg)%U3DX(I,J,K+1)          &
       -(channel%seg(num_seg)%GCONT_U(1,I,J)                                          &
           *((channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))/DX   &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K))        &
         +ADD_U(K))*DZ/FNZ(K)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = D0_0
      ENDDO

      channel%seg(num_seg)%U3DX(I,J,1)   = channel%seg(num_seg)%U3DX(I,J,2)
      channel%seg(num_seg)%U3DX(I,J,nk3) = channel%seg(num_seg)%U3DX(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!---------------
! V-component
!---------------

      DO J = mjm, mjp_a
      DO I = mim_a, mip_a

      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K=KLOWV,nk1
       ADD_V(K) = &
       (channel%seg(num_seg)%GCONT_U(2,I,J+1)*(channel%seg(num_seg)%W3D(I+1,J+1,K)    &
                                              -channel%seg(num_seg)%W3D(I,J+1,K))     &
       +channel%seg(num_seg)%GCONT_U(2,I-1,J+1)*(channel%seg(num_seg)%W3D(I,J+1,K)    &
                                                -channel%seg(num_seg)%W3D(I-1,J+1,K)) &
       +channel%seg(num_seg)%GCONT_U(2,I,J)*(channel%seg(num_seg)%W3D(I+1,J,K)        &
                                            -channel%seg(num_seg)%W3D(I,J,K))         &
       +channel%seg(num_seg)%GCONT_U(2,I-1,J)*(channel%seg(num_seg)%W3D(I,J,K)        &
                                              -channel%seg(num_seg)%W3D(I-1,J,K)))    &
       /(4.dbl_kind*DX)                                                               &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(D4_0*channel%seg(num_seg)%RG_V(I,J))

      ENDDO
      DO K = nk1, KLOWV, -1
       channel%seg(num_seg)%U3DY(I,J,K) = channel%seg(num_seg)%U3DY(I,J,K+1)          &
       -(channel%seg(num_seg)%GCONT_V(3,I,J)                                          &
           *((channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))/DY   &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K))        &
         +ADD_V(K))*DZ/FNZ(K)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = D0_0
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!---------------------------------------------------------------------
!     2. Covariant components of horizontal wind

!     (Used in a limited way: calculation of deformation and vorticities)
!---------------------------------------------------------------------
!---------------
! U-component
!---------------

      DO J = mjm_ce, mjp_ce
       DO I = mim_ce, mip_c
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

        DO K = nk1, KLOWU, -1
         channel%seg(num_seg)%U3DX_CO(I,J,K) = channel%seg(num_seg)%U3DX_CO(I,J,K+1)       &
                 -((channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))/DX  &
                   +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K))*DZ/FNZ(K)
        ENDDO

        DO K = 2, KLOWU-1
         channel%seg(num_seg)%U3DX_CO(I,J,K) =  &                                          
                channel%seg(num_seg)%GG_U(3,I,J)*channel%seg(num_seg)%U3DX(I,J,K)           &
              - (channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%U3DY(I,J,K)          &
                +channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)      &
                +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)      &
                +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K)) &
                /D4_0
        ENDDO

        channel%seg(num_seg)%U3DX_CO(I,J,1)   = channel%seg(num_seg)%U3DX_CO(I,J,2)
        channel%seg(num_seg)%U3DX_CO(I,J,nk3) = channel%seg(num_seg)%U3DX_CO(I,J,nk2)
       ENDDO
      ENDDO

!---------------
! V-component
!---------------

      DO J = mjm_ce, mjp_c
       DO I = mim_ce, mip_ce
        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

        DO K = nk1, KLOWV, -1
         channel%seg(num_seg)%U3DY_CO(I,J,K) = channel%seg(num_seg)%U3DY_CO(I,J,K+1)       &
                 -((channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))/DY  &
                   -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K))*DZ/FNZ(K)
        ENDDO

        DO K = 2, KLOWV -1
         channel%seg(num_seg)%U3DY_CO(I,J,K) = &
                channel%seg(num_seg)%GG_V(1,I,J)*channel%seg(num_seg)%U3DY(I,J,K)           &
              - (channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%U3DX(I,J,K)          &
                +channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K)) &
                /D4_0
        ENDDO

        channel%seg(num_seg)%U3DY_CO(I,J,1)   = channel%seg(num_seg)%U3DY_CO(I,J,2)
        channel%seg(num_seg)%U3DY_CO(I,J,nk3) = channel%seg(num_seg)%U3DY_CO(I,J,nk2)
       ENDDO
      ENDDO

!**********************************************************************   
      ENDDO  ! num_seg
!**********************************************************************

      ENDIF   ! UCHANGE

      END SUBROUTINE wind_3d

END MODULE wind_module
