MODULE wind_module
! diagnose wind components

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk1,nk2,nk3
USE constld, only: dx,dy,dz,fnz,fnt,rho,rhoz

! Subroutines being called
USE elliptic, only: relax_3d
USE uvtop,    only: uvcal_3d       

IMPLICIT NONE
PRIVATE

PUBLIC :: wind_3d_v2,wind_3d_v1,wind_3d_v0

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE wind_3d : diagnose wind components
!-----------------------------------------------------------------------

!=======================================================================
   SUBROUTINE WIND_3D_V2 (UCHANGE,WCHANGE,channel,NITER2D, ITTV)
!=======================================================================  
!  Version 2
!  JUNG: Explicitly apply the definition of zeta (across the channel)
!        Explicitly apply the continuity equation (across the channel)
!        Calculate the covariant components based on the analytic equation. 
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
!  JUNG: The way of calculating wind at ghost points is different from LQ3D.
!------------------------------------------------------------------------ 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      LOGICAL, INTENT(IN) :: UCHANGE, WCHANGE

      INTEGER, INTENT(IN), OPTIONAL :: NITER2D   
      ! number of iteration for 2D elliptic equation, 
      ! which is different from the basic setting "niterxy" (test purpose) 

      INTEGER, INTENT(IN), OPTIONAL :: ITTV   ! JUNG_DEBUG
       
!     Local Variables
      REAL (KIND=r8), DIMENSION(nk1) :: ADD_U,ADD_V
      
      integer :: ifortw   ! remove

      INTEGER :: I,J,K,KLOWU,KLOWV,num_seg
      INTEGER :: mi1,mim,mip,mim_c,mim_a,mip_a
      INTEGER :: mj1,mjm,mjp,mjm_c,mjm_a,mjp_a
!------------------------------------------------------------------------    
      if (present(ittv)) then  
        ifortw = channel%num_chn 
        if (ifortw == 5) ifortw=55
        if (ifortw == 6) ifortw=65
      endif

      if (present(ittv)) write(ifortw,*) 'In WIND_3D b4 RELAX'

      IF (WCHANGE) CALL RELAX_3D (channel)

      if (present(ittv)) write(ifortw,*) 'In WIND_3D a4 RELAX'
      
      IF (UCHANGE) THEN

      ! UPDATE wind at the model top
      IF (.not.PRESENT(NITER2D)) THEN
        IF (present(ittv)) then
          CALL UVCAL_3D (channel,ittv=ittv)
        else
          CALL UVCAL_3D (channel)
        endif
      ELSE
        CALL UVCAL_3D (channel,NITER=NITER2D)
      ENDIF

      if (present(ittv)) write(ifortw,*) 'In WIND_3D  a4 UVLOW'
      
      ! UPDATE wind above the lowest model layer
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1    = channel%seg(num_seg)%mi1
      mj1    = channel%seg(num_seg)%mj1
      
      mim    = channel%seg(num_seg)%mim
      mjm    = channel%seg(num_seg)%mjm 
      
      mip    = channel%seg(num_seg)%mip
      mjp    = channel%seg(num_seg)%mjp  
      
      mim_c  = channel%seg(num_seg)%mim_c
      mjm_c  = channel%seg(num_seg)%mjm_c 
      
      mim_a  = channel%seg(num_seg)%mim_a 
      mjm_a  = channel%seg(num_seg)%mjm_a
      
      mip_a  = channel%seg(num_seg)%mip_a
      mjp_a  = channel%seg(num_seg)%mjp_a 

!========================================================== ! x-array 
      IF (mi1 .gt. mj1) THEN
!==========================================================

! 1. Covariant components of horizontal wind
!    (Explicitly satisfy the relationship between the wind and Zeta)
!---------------
! V-component (Covariant components)
!---------------
      DO J = mjm_a, mjp_a
       DO I = mim_a, mip
        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
        DO K = KLOWV+1, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = channel%seg(num_seg)%U3DY_CO(I,J,K-1)           &
                 +((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY  &
                   -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K-1))*DZ/FNZ(K-1)          
        ENDDO
       ENDDO
      ENDDO
!---------------
! U-component (Covariant components)
!---------------
      DO J = 1, mj1
       DO I = mim_a, mip_a
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
        DO K = KLOWU+1, nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) = channel%seg(num_seg)%U3DX_CO(I,J,K-1)           &
                 +((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX  &
                   +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K-1))*DZ/FNZ(K-1)          
        ENDDO
       ENDDO
      ENDDO
      
      ! across the channel (north +): Apply Zeta-definition
      DO J = mj1+1, mjp
       DO I = mim_a, mip_a
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

        DO K = KLOWU+1, nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) = channel%seg(num_seg)%U3DX_CO(I,J-1,K)         &
        +DY*((channel%seg(num_seg)%U3DY_CO(I+1,J-1,K)-channel%seg(num_seg)%U3DY_CO(I,J-1,K))/DX  &
            -channel%seg(num_seg)%RG_Z(I,J-1)*channel%seg(num_seg)%Z3DZ(I,J-1,K))            
        ENDDO
       ENDDO
      ENDDO
      ! across the channel (south -): Apply Zeta-definition
      DO J = 0, mjm_a, -1
       DO I = mim_a, mip_a
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

        DO K = KLOWU+1, nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) = channel%seg(num_seg)%U3DX_CO(I,J+1,K)         &
        -DY*((channel%seg(num_seg)%U3DY_CO(I+1,J,K)-channel%seg(num_seg)%U3DY_CO(I,J,K))/DX  &
            -channel%seg(num_seg)%RG_Z(I,J)*channel%seg(num_seg)%Z3DZ(I,J,K))
        ENDDO
       ENDDO
      ENDDO                  

! 2. Contravariant components of horizontal wind
!    (Explicitly satisfy the continuity equation)
!----------------------------------------------
!  U-component (Contravariant components)
!----------------------------------------------
      ! Calculate new eta-value (Z3DY0) 
      DO I = mim_a, mip_a
      
      DO J = mjm_a, 0
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
        DO K = KLOWU, nk1
        channel%seg(num_seg)%Z3DY0(I,J,K) = &
         - (channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))/DX      &
         + (channel%seg(num_seg)%U3DX_CO(I,J,K+1)-channel%seg(num_seg)%U3DX_CO(I,J,K)) &
          *FNZ(K)/DZ
        ENDDO
        DO K = KLOWU, nk1
        channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%Z3DY0(I,J,K)  &
                                           /channel%seg(num_seg)%RG_U(I,J)
        ENDDO
      ENDDO      
         
      DO J = 1, mj1
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
        DO K = KLOWU, nk1
        channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K)
        ENDDO
      ENDDO
      
      DO J = mj1+1, mjp_a
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
        DO K = KLOWU, nk1
        channel%seg(num_seg)%Z3DY0(I,J,K) = &
         - (channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))/DX      &
         + (channel%seg(num_seg)%U3DX_CO(I,J,K+1)-channel%seg(num_seg)%U3DX_CO(I,J,K)) &
          *FNZ(K)/DZ
        ENDDO
        DO K = KLOWU, nk1
        channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%Z3DY0(I,J,K)  &
                                           /channel%seg(num_seg)%RG_U(I,J)
        ENDDO
      ENDDO
      
      ENDDO   ! I-loop      

      ! original method, but using Z3DY0    
      DO J = mjm_a, mjp_a
      DO I = mim_a, mip_a
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
       /(4.0_r8*DY)                                                                   &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = KLOWU+1,nk2
       channel%seg(num_seg)%U3DX(I,J,K) = channel%seg(num_seg)%U3DX(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_U(1,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY0(I,J,K-1))       &
         +ADD_U(K-1))*DZ/FNZ(K-1)   
      ENDDO

!     kinematic boundary condition
      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DX(I,J,1)   = channel%seg(num_seg)%U3DX(I,J,2)
      channel%seg(num_seg)%U3DX(I,J,nk3) = channel%seg(num_seg)%U3DX(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!----------------------------------------------
!  V-component (Contravariant components)
!----------------------------------------------
      ! original method
      DO J = 1, mj1
      DO I = mim_c, mip_a
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
       /(4.0_r8*DX)                                                                   &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_V(I,J))

      ENDDO
           
      DO K =KLOWV+1, nk2
       channel%seg(num_seg)%U3DY(I,J,K) = channel%seg(num_seg)%U3DY(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_V(3,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K-1))        &
         +ADD_V(K-1))*DZ/FNZ(K-1)  
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (north +): Apply the continuity eq. 
      DO J = mj1+1, mjp_a
      DO I = mim_c, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  &     
          channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)  &           
         -DY*FNT(K)*channel%seg(num_seg)%RG_T(I,J)                             &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J,K)                   &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J,K-1))/(DZ*RHO(K))  &
         -DY*(channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%U3DX(I,J,K)  & 
             -channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K))/DX       
      ENDDO
      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  &     
              channel%seg(num_seg)%U3DY(I,J,K)/channel%seg(num_seg)%RG_V(I,J)
      ENDDO      
                 
      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
      ! across the channel (south -): Apply the continuity eq. 
      DO J = 0, mjm, -1
      DO I = mim_c, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  & 
          channel%seg(num_seg)%RG_V(I,J+1)*channel%seg(num_seg)%U3DY(I,J+1,K)      &                
         +DY*FNT(K)*channel%seg(num_seg)%RG_T(I,J+1)                               &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J+1,K)                     &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J+1,K-1))/(DZ*RHO(K))    &
         +DY*(channel%seg(num_seg)%RG_U(I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)  & 
             -channel%seg(num_seg)%RG_U(I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K))/DX            
      ENDDO
      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  & 
              channel%seg(num_seg)%U3DY(I,J,K)/channel%seg(num_seg)%RG_V(I,J)  
      ENDDO
                  
      ENDDO  ! I-loop
      ENDDO  ! J-loop      

!     Boundary Values
      DO J = mjm, mjp_a
      DO I = mim_c, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
      
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!  Boundary Values of Covariant components
!  JH: If topography is included, correct the below.
!---------------
! U-component (Covariant components)
!---------------      
!      DO J = mjm_a, mjp_a
!       DO I = mim_c, mip_c
!        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
!        DO K = 2, KLOWU-1
!         channel%seg(num_seg)%U3DX_CO(I,J,K) =  &                                          
!                channel%seg(num_seg)%GG_U(3,I,J)*channel%seg(num_seg)%U3DX(I,J,K)           &
!              - (channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%U3DY(I,J,K)          &
!                +channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)      &
!                +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)      &
!                +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K)) &
!                /4.0_r8
!        ENDDO
!       ENDDO
!      ENDDO      
!---------------
! V-component (Covariant components)
!---------------
!      DO J = mjm_a, mjp_c
!       DO I = mim_c, mip_a
!        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
!
!        DO K = 2, KLOWV-1
!         channel%seg(num_seg)%U3DY_CO(I,J,K) = &
!                channel%seg(num_seg)%GG_V(1,I,J)*channel%seg(num_seg)%U3DY(I,J,K)           &
!              - (channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%U3DX(I,J,K)          &
!                +channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)      &
!                +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)      &
!                +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K)) &
!                /4.0_r8
!        ENDDO    
!       ENDDO
!      ENDDO

!========================================================== ! y-array 
      ELSE
!==========================================================

! 1. Covariant components of horizontal wind
!    (Explicitly satisfy the relationship between the wind and Zeta)
!---------------
! U-component (Covariant components)
!---------------
      DO J = mjm_a, mjp
       DO I = mim_a, mip_a
        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
        DO K = KLOWU+1,nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) = channel%seg(num_seg)%U3DX_CO(I,J,K-1)           &
                 +((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX  &
                   +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K-1))*DZ/FNZ(K-1)
        ENDDO
       ENDDO
      ENDDO

!---------------
! v-component (Covariant components)
!---------------
      DO J = mjm_a, mjp_a
       DO I = 1, mi1
        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
        DO K = KLOWV+1, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = channel%seg(num_seg)%U3DY_CO(I,J,K-1)           &
                 +((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY  &
                   -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K-1))*DZ/FNZ(K-1)
        ENDDO
       ENDDO
      ENDDO

      ! across the channel (east +): Apply Zeta-definition
      DO J = mjm_a, mjp_a
       DO I = mi1+1, mip
        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

        DO K = KLOWV+1, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = channel%seg(num_seg)%U3DY_CO(I-1,J,K)             &
        +DX*((channel%seg(num_seg)%U3DX_CO(I-1,J+1,K)-channel%seg(num_seg)%U3DX_CO(I-1,J,K))/DY  &
            +channel%seg(num_seg)%RG_Z(I-1,J)*channel%seg(num_seg)%Z3DZ(I-1,J,K))
        ENDDO
       ENDDO
      ENDDO
      ! across the channel (south -): Apply Zeta-definition
      DO J = mjm_a, mjp_a
       DO I = 0, mim_a, -1
        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

        DO K = KLOWV+1, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = channel%seg(num_seg)%U3DY_CO(I+1,J,K)         &
        -DX*((channel%seg(num_seg)%U3DX_CO(I,J+1,K)-channel%seg(num_seg)%U3DX_CO(I,J,K))/DY  &
            +channel%seg(num_seg)%RG_Z(I,J)*channel%seg(num_seg)%Z3DZ(I,J,K))
        ENDDO
       ENDDO
      ENDDO            

! 2. Contravariant components of horizontal wind
!    (Explicitly satisfy the continuity equation)
!----------------------------------------------
!  V-component (Contravariant components)
!----------------------------------------------
      ! Calculate new ksi-value (Z3DX0) 
      DO J = mjm_a, mjp_a
      
        DO I = mim_a, 0
         KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
         DO K=KLOWV,nk1
           channel%seg(num_seg)%Z3DX0(I,J,K) = &
              (channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))/DY      &
            - (channel%seg(num_seg)%U3DY_CO(I,J,K+1)-channel%seg(num_seg)%U3DY_CO(I,J,K)) &
             *FNZ(K)/DZ
         ENDDO
         DO K=KLOWV,nk1
           channel%seg(num_seg)%Z3DX0(I,J,K) = channel%seg(num_seg)%Z3DX0(I,J,K) &
                                              /channel%seg(num_seg)%RG_V(I,J)
         ENDDO
        ENDDO 
        
        DO I = 1, mi1
         KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
         DO K=KLOWV,nk1
           channel%seg(num_seg)%Z3DX0(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K)
         ENDDO
        ENDDO  
        
        DO I = mi1+1, mip_a   
         KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
         DO K=KLOWV,nk1
           channel%seg(num_seg)%Z3DX0(I,J,K) = &
              (channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))/DY &
            - (channel%seg(num_seg)%U3DY_CO(I,J,K+1)-channel%seg(num_seg)%U3DY_CO(I,J,K)) &
             *FNZ(K)/DZ
         ENDDO
         DO K=KLOWV,nk1
           channel%seg(num_seg)%Z3DX0(I,J,K) = channel%seg(num_seg)%Z3DX0(I,J,K) &
                                              /channel%seg(num_seg)%RG_V(I,J)
         ENDDO
        ENDDO 
        
      ENDDO  ! J-loop      

      ! original method, but using Z3DX0    
      DO J = mjm_a, mjp_a
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
       /(4.0_r8*DX)                                                                   &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_V(I,J))
      ENDDO
           
      DO K =KLOWV+1, nk2
       channel%seg(num_seg)%U3DY(I,J,K) = channel%seg(num_seg)%U3DY(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_V(3,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX0(I,J,K-1))       &
         +ADD_V(K-1))*DZ/FNZ(K-1)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!----------------------------------------------
!  U-component (Contravariant components)
!----------------------------------------------
      ! original method    
      DO J = mjm_c, mjp_a
      DO I = 1, mi1
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
       /(4.0_r8*DY)                                                                   &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = KLOWU+1,nk2
       channel%seg(num_seg)%U3DX(I,J,K) = channel%seg(num_seg)%U3DX(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_U(1,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K-1))        &
         +ADD_U(K-1))*DZ/FNZ(K-1)
      ENDDO
      
      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (East +): Apply the continuity eq.
      DO J = mjm_c, mjp_a
      DO I = mi1+1, mip_a
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
      
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = & 
          channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)   &
         -DX*FNT(K)*channel%seg(num_seg)%RG_T(I,J)                              &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J,K)                    &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J,K-1))/(DZ*RHO(K))   &
         -DX*(channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%U3DY(I,J,K)   &
             -channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K))/DY  
      ENDDO
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = &
              channel%seg(num_seg)%U3DX(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
      ENDDO
      
      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (West -): Apply the continuity eq.
      DO J = mjm_c, mjp_a
      DO I = 0, mim, -1
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
      
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) =  &                 
          channel%seg(num_seg)%RG_U(I+1,J)*channel%seg(num_seg)%U3DX(I+1,J,K)      & 
         +DX*FNT(K)*channel%seg(num_seg)%RG_T(I+1,J)                               &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I+1,J,K)                     &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I+1,J,K-1))/(DZ*RHO(K))    &
         +DX*(channel%seg(num_seg)%RG_V(I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)  &
             -channel%seg(num_seg)%RG_V(I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K))/DY  
    
      ENDDO
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = &                 
              channel%seg(num_seg)%U3DX(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!     Boundary Values
      DO J = mjm_c, mjp_a
      DO I = mim, mip_a
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DX(I,J,1)   = channel%seg(num_seg)%U3DX(I,J,2)
      channel%seg(num_seg)%U3DX(I,J,nk3) = channel%seg(num_seg)%U3DX(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!  Boundary Values of Covariant components
!  JH: If topography is included, correct the below.
!---------------
! U-component (Covariant components)
!---------------      
!      DO J = mjm_a, mjp
!       DO I = mim_a, mip_a
!        KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
!
!        DO K = 2, KLOWU-1
!         channel%seg(num_seg)%U3DX_CO(I,J,K) =  &                                          
!                channel%seg(num_seg)%GG_U(3,I,J)*channel%seg(num_seg)%U3DX(I,J,K)           &
!              - (channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%U3DY(I,J,K)          &
!                +channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)      &
!                +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)      &
!                +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K)) &
!                /4.0_r8
!        ENDDO
!       ENDDO
!      ENDDO      
!---------------
! V-component (Covariant components)
!---------------
!      DO J = mjm_a, mjp_a
!       DO I = mim_a, mip
!        KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
!
!        DO K = 2, KLOWV-1
!         channel%seg(num_seg)%U3DY_CO(I,J,K) = &
!                channel%seg(num_seg)%GG_V(1,I,J)*channel%seg(num_seg)%U3DY(I,J,K)           &
!              - (channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%U3DX(I,J,K)          &
!                +channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)      &
!                +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)      &
!                +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K)) &
!                /4.0_r8
!        ENDDO     
!       ENDDO
!      ENDDO

!==========================================================
      ENDIF 
!========================================================== 
   
!  Boundary Values of Covariant components
!---------------
! U-component (Covariant components)
!---------------      
      DO J = mjm_a, mjp
       DO I = mim_a, mip_a
        channel%seg(num_seg)%U3DX_CO(I,J,1)   = channel%seg(num_seg)%U3DX_CO(I,J,2)
        channel%seg(num_seg)%U3DX_CO(I,J,nk3) = channel%seg(num_seg)%U3DX_CO(I,J,nk2)
       ENDDO
      ENDDO
!---------------
! V-component (Covariant components)
!---------------
      DO J = mjm_a, mjp_a
       DO I = mim_a, mip
        channel%seg(num_seg)%U3DY_CO(I,J,1)   = channel%seg(num_seg)%U3DY_CO(I,J,2)
        channel%seg(num_seg)%U3DY_CO(I,J,nk3) = channel%seg(num_seg)%U3DY_CO(I,J,nk2)        
       ENDDO
      ENDDO
                        
!**********************************************************************   
      ENDDO  ! num_seg
!**********************************************************************

      ENDIF   ! UCHANGE

      END SUBROUTINE wind_3d_v2

!=======================================================================
   SUBROUTINE WIND_3D_V1 (UCHANGE,WCHANGE,channel,NITER2D, ITTV)
!=======================================================================   
!  version_1
!  JUNG: Explicitly apply the continuity equation (across the channel)
!        Calculate the covariant components based on the contravariant components.
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
!  JUNG: The way of calculating wind at ghost points is different from LQ3D.
!------------------------------------------------------------------------ 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      LOGICAL, INTENT(IN) :: UCHANGE, WCHANGE

      INTEGER, INTENT(IN), OPTIONAL :: NITER2D   
      ! number of iteration for 2D elliptic equation, 
      ! which is different from the basic setting "niterxy" (test purpose) 

      INTEGER, INTENT(IN), OPTIONAL :: ITTV   ! JUNG_DEBUG
       
!     Local Variables
      REAL (KIND=r8), DIMENSION(nk1) :: ADD_U,ADD_V
      
      integer :: ifortw   ! remove

      INTEGER :: I,J,K,KLOWU,KLOWV,num_seg
      INTEGER :: mi1,mim,mip,mim_c,mip_c,mim_a,mip_a
      INTEGER :: mj1,mjm,mjp,mjm_c,mjp_c,mjm_a,mjp_a
!------------------------------------------------------------------------    
      if (present(ittv)) then  
        ifortw = channel%num_chn 
        if (ifortw == 5) ifortw=55
        if (ifortw == 6) ifortw=65
      endif

      if (present(ittv)) write(ifortw,*) 'In WIND_3D b4 RELAX'

      IF (WCHANGE) CALL RELAX_3D (channel)

      if (present(ittv)) write(ifortw,*) 'In WIND_3D a4 RELAX'
      
      IF (UCHANGE) THEN

      ! UPDATE wind at the model top
      IF (.not.PRESENT(NITER2D)) THEN
        IF (present(ittv)) then
          CALL UVCAL_3D (channel,ittv=ittv)
        else
          CALL UVCAL_3D (channel)
        endif
      ELSE
        CALL UVCAL_3D (channel,NITER=NITER2D)
      ENDIF

      if (present(ittv)) write(ifortw,*) 'In WIND_3D  a4 UVLOW'
      
      ! UPDATE wind above the lowest model layer
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1    = channel%seg(num_seg)%mi1
      mj1    = channel%seg(num_seg)%mj1
      
      mim    = channel%seg(num_seg)%mim
      mjm    = channel%seg(num_seg)%mjm 
      
      mip    = channel%seg(num_seg)%mip
      mjp    = channel%seg(num_seg)%mjp  
      
      mim_c  = channel%seg(num_seg)%mim_c
      mjm_c  = channel%seg(num_seg)%mjm_c 

      mip_c  = channel%seg(num_seg)%mip_c
      mjp_c  = channel%seg(num_seg)%mjp_c 
            
      mim_a  = channel%seg(num_seg)%mim_a 
      mjm_a  = channel%seg(num_seg)%mjm_a
      
      mip_a  = channel%seg(num_seg)%mip_a
      mjp_a  = channel%seg(num_seg)%mjp_a 

!========================================================== ! x-array 
      IF (mi1 .gt. mj1) THEN
!==========================================================
! 1. Contravariant components of horizontal wind
!    (Explicitly satisfy the continuity equation)
!----------------------------------------------
!  U-component (Contravariant components)
!----------------------------------------------
      ! original method    
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
       /(4.0_r8*DY)                                                                   &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = KLOWU+1,nk2
       channel%seg(num_seg)%U3DX(I,J,K) = channel%seg(num_seg)%U3DX(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_U(1,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K-1))        &
         +ADD_U(K-1))*DZ/FNZ(K-1)   
      ENDDO

!     kinematic boundary condition
      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DX(I,J,1)   = channel%seg(num_seg)%U3DX(I,J,2)
      channel%seg(num_seg)%U3DX(I,J,nk3) = channel%seg(num_seg)%U3DX(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!----------------------------------------------
!  V-component (Contravariant components)
!----------------------------------------------
      ! original method
      DO J = 1, mj1
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
       /(4.0_r8*DX)                                                                   &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_V(I,J))

      ENDDO
           
      DO K =KLOWV+1, nk2
       channel%seg(num_seg)%U3DY(I,J,K) = channel%seg(num_seg)%U3DY(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_V(3,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K-1))        &
         +ADD_V(K-1))*DZ/FNZ(K-1)  
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (north +): Apply the continuity eq. 
      DO J = mj1+1, mjp_a
      DO I = mim_a, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  &     
          channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)  &           
         -DY*FNT(K)*channel%seg(num_seg)%RG_T(I,J)                             &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J,K)                   &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J,K-1))/(DZ*RHO(K))  &
         -DY*(channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%U3DX(I,J,K)  & 
             -channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K))/DX       
      ENDDO
      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  &     
              channel%seg(num_seg)%U3DY(I,J,K)/channel%seg(num_seg)%RG_V(I,J)
      ENDDO      
                 
      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
      ! across the channel (south -): Apply the continuity eq. 
      DO J = 0, mjm, -1
      DO I = mim_a, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  & 
          channel%seg(num_seg)%RG_V(I,J+1)*channel%seg(num_seg)%U3DY(I,J+1,K)      &                
         +DY*FNT(K)*channel%seg(num_seg)%RG_T(I,J+1)                               &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J+1,K)                     &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J+1,K-1))/(DZ*RHO(K))    &
         +DY*(channel%seg(num_seg)%RG_U(I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)  & 
             -channel%seg(num_seg)%RG_U(I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K))/DX            
      ENDDO
      DO K=KLOWV+1,nk2
      channel%seg(num_seg)%U3DY(I,J,K) =  & 
              channel%seg(num_seg)%U3DY(I,J,K)/channel%seg(num_seg)%RG_V(I,J)  
      ENDDO
                  
      ENDDO  ! I-loop
      ENDDO  ! J-loop      

!     Boundary Values
      DO J = mjm, mjp_a
      DO I = mim_a, mip_a
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
      
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!========================================================== ! y-array 
      ELSE
!==========================================================

! 1. Contravariant components of horizontal wind
!    (Explicitly satisfy the continuity equation)
!----------------------------------------------
!  V-component (Contravariant components)
!----------------------------------------------

      ! original method    
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
       /(4.0_r8*DX)                                                                   &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_V(I,J))
      ENDDO
           
      DO K =KLOWV+1, nk2
       channel%seg(num_seg)%U3DY(I,J,K) = channel%seg(num_seg)%U3DY(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_V(3,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I,J+1,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DY &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K-1))       &
         +ADD_V(K-1))*DZ/FNZ(K-1)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop
      
!----------------------------------------------
!  U-component (Contravariant components)
!----------------------------------------------
      ! original method    
      DO J = mjm_a, mjp_a
      DO I = 1, mi1
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
       /(4.0_r8*DY)                                                                   &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = KLOWU+1,nk2
       channel%seg(num_seg)%U3DX(I,J,K) = channel%seg(num_seg)%U3DX(I,J,K-1)            &
       +(channel%seg(num_seg)%GCONT_U(1,I,J)                                            &
           *((channel%seg(num_seg)%W3D(I+1,J,K-1)-channel%seg(num_seg)%W3D(I,J,K-1))/DX &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K-1))        &
         +ADD_U(K-1))*DZ/FNZ(K-1)
      ENDDO
      
      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (East +): Apply the continuity eq.
      DO J = mjm_a, mjp_a
      DO I = mi1+1, mip_a
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
      
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = & 
          channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)   &
         -DX*FNT(K)*channel%seg(num_seg)%RG_T(I,J)                              &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I,J,K)                    &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I,J,K-1))/(DZ*RHO(K))   &
         -DX*(channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%U3DY(I,J,K)   &
             -channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K))/DY  
      ENDDO
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = &
              channel%seg(num_seg)%U3DX(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
      ENDDO
      
      ENDDO  ! I-loop
      ENDDO  ! J-loop

      ! across the channel (West -): Apply the continuity eq.
      DO J = mjm_a, mjp_a
      DO I = 0, mim, -1
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
      
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) =  &                 
          channel%seg(num_seg)%RG_U(I+1,J)*channel%seg(num_seg)%U3DX(I+1,J,K)      & 
         +DX*FNT(K)*channel%seg(num_seg)%RG_T(I+1,J)                               &
                   *(RHOZ(K)*channel%seg(num_seg)%W3D(I+1,J,K)                     &
                    -RHOZ(K-1)*channel%seg(num_seg)%W3D(I+1,J,K-1))/(DZ*RHO(K))    &
         +DX*(channel%seg(num_seg)%RG_V(I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)  &
             -channel%seg(num_seg)%RG_V(I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K))/DY  
    
      ENDDO
      DO K = KLOWU+1, nk2
      channel%seg(num_seg)%U3DX(I,J,K) = &                 
              channel%seg(num_seg)%U3DX(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!     Boundary Values
      DO J = mjm_a, mjp_a
      DO I = mim, mip_a
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DX(I,J,1)   = channel%seg(num_seg)%U3DX(I,J,2)
      channel%seg(num_seg)%U3DX(I,J,nk3) = channel%seg(num_seg)%U3DX(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!==========================================================
      ENDIF 
!========================================================== 
   
!---------------------------------------------------------------------
!     2. Covariant components: Obtained from contravariant components 
!        Used: calculation of deformation (turbulence) and vorticities (relax3d)
!---------------------------------------------------------------------
!     Subroutine UVCAL_3D calculates the covariant components at k = 2.
 
!---------------
! U-component
!---------------
      DO J = mjm_a, mjp_a
       DO I = mim_a, mip_c
        DO K = 3, nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) =  &                                          
                channel%seg(num_seg)%GG_U(3,I,J)*channel%seg(num_seg)%U3DX(I,J,K)           &
              - (channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%U3DY(I,J,K)          &
                +channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)      &
                +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)      &
                +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K)) &
                /4.0_r8
        ENDDO  
        channel%seg(num_seg)%U3DX_CO(I,J,1)   = channel%seg(num_seg)%U3DX_CO(I,J,2)
        channel%seg(num_seg)%U3DX_CO(I,J,nk3) = channel%seg(num_seg)%U3DX_CO(I,J,nk2)                          
       ENDDO
      ENDDO          

!---------------
! V-component
!---------------  
      DO J = mjm_a, mjp_c  
       DO I = mim_a, mip_a
        DO K = 3, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = &
                channel%seg(num_seg)%GG_V(1,I,J)*channel%seg(num_seg)%U3DY(I,J,K)           &
              - (channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%U3DX(I,J,K)          &
                +channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K)) &
                /4.0_r8
        ENDDO          
        channel%seg(num_seg)%U3DY_CO(I,J,1)   = channel%seg(num_seg)%U3DY_CO(I,J,2)
        channel%seg(num_seg)%U3DY_CO(I,J,nk3) = channel%seg(num_seg)%U3DY_CO(I,J,nk2)                          
       ENDDO
      ENDDO      
      
!**********************************************************************   
      ENDDO  ! num_seg
!**********************************************************************

      ENDIF   ! UCHANGE

      END SUBROUTINE wind_3d_v1

!=======================================================================
   SUBROUTINE WIND_3D_V0 (UCHANGE,WCHANGE,channel,NITER2D, ITTV)
!=======================================================================   
!  version_0
!  JUNG: No application of the continuity equation & zeta definition (across the channel)
!        Calculate the covariant components based on the contravariant components.
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
!  JUNG: The way of calculating wind at ghost points is different from LQ3D.
!------------------------------------------------------------------------ 
      type(channel_t), intent(inout) :: channel   ! channel data
      
      LOGICAL, INTENT(IN) :: UCHANGE, WCHANGE

      INTEGER, INTENT(IN), OPTIONAL :: NITER2D   
      ! number of iteration for 2D elliptic equation, 
      ! which is different from the basic setting "niterxy" (test purpose) 

      INTEGER, INTENT(IN), OPTIONAL :: ITTV   ! JUNG_DEBUG
       
!     Local Variables
      REAL (KIND=r8), DIMENSION(nk1) :: ADD_U,ADD_V
      
      integer :: ifortw   ! remove

      INTEGER :: I,J,K,KLOWU,KLOWV,num_seg
      INTEGER :: mim,mim_a,mim_ce,mip_a,mip_c,mip_ce  
      INTEGER :: mjm,mjm_a,mjm_ce,mjp_a,mjp_c,mjp_ce 
!------------------------------------------------------------------------      
      if (present(ittv)) then  
        ifortw = channel%num_chn 
        if (ifortw == 5) ifortw=55
        if (ifortw == 6) ifortw=65
      endif

      if (present(ittv)) write(ifortw,*) 'In WIND_3D b4 RELAX'

      IF (WCHANGE) CALL RELAX_3D (channel)

      if (present(ittv)) write(ifortw,*) 'In WIND_3D a4 RELAX'
      
      IF (UCHANGE) THEN

      ! UPDATE wind at the model top
      IF (.not.PRESENT(NITER2D)) THEN
        IF (present(ittv)) then
          CALL UVCAL_3D (channel,ittv=ittv)
        else
          CALL UVCAL_3D (channel)
        endif
      ELSE
        CALL UVCAL_3D (channel,NITER=NITER2D)
      ENDIF

      if (present(ittv)) write(ifortw,*) 'In WIND_3D  a4 UVLOW'
      
      ! UPDATE wind above the lowest model layer
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
       /(4.0_r8*DY)                                                                   &
     - (channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K)         &
       +channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%Z3DX(I,J,K)             &
       +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX(I+1,J-1,K)     &
       +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%Z3DX(I,J-1,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_U(I,J))
      ENDDO

      DO K = KLOWU,nk1
       channel%seg(num_seg)%U3DX(I,J,K+1) = channel%seg(num_seg)%U3DX(I,J,K)          &
       +(channel%seg(num_seg)%GCONT_U(1,I,J)                                          &
           *((channel%seg(num_seg)%W3D(I+1,J,K)-channel%seg(num_seg)%W3D(I,J,K))/DX   &
             +channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K))        &
         +ADD_U(K))*DZ/FNZ(K)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWU-1
       channel%seg(num_seg)%U3DX(I,J,K) = 0.0_r8
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
       /(4.0_r8*DX)                                                                   &
     + (channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K)         &
       +channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%Z3DY(I,J,K)             &
       +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY(I-1,J+1,K)     &
       +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%Z3DY(I-1,J,K))        &
       /(4.0_r8*channel%seg(num_seg)%RG_V(I,J))

      ENDDO
           
      DO K =KLOWV, nk1
       channel%seg(num_seg)%U3DY(I,J,K+1) = channel%seg(num_seg)%U3DY(I,J,K)          &
       +(channel%seg(num_seg)%GCONT_V(3,I,J)                                          &
           *((channel%seg(num_seg)%W3D(I,J+1,K)-channel%seg(num_seg)%W3D(I,J,K))/DY   &
             -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K))        &
         +ADD_V(K))*DZ/FNZ(K)
      ENDDO
      
!     kinematic boundary condition
      DO K = 2, KLOWV-1
       channel%seg(num_seg)%U3DY(I,J,K) = 0.0_r8
      ENDDO

      channel%seg(num_seg)%U3DY(I,J,1)   = channel%seg(num_seg)%U3DY(I,J,2)
      channel%seg(num_seg)%U3DY(I,J,nk3) = channel%seg(num_seg)%U3DY(I,J,nk2)

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!---------------------------------------------------------------------
!     2. Covariant components: Obtained from contravariant components 
!        Used: calculation of deformation (turbulence) and vorticities (relax3d)
!---------------------------------------------------------------------
!     Subroutine UVCAL_3D calculates the covariant components at k = 2.
 
!---------------
! U-component
!---------------
      DO J = mjm_a, mjp_a
       DO I = mim_a, mip_c
        DO K = 3, nk2
         channel%seg(num_seg)%U3DX_CO(I,J,K) =  &                                          
                channel%seg(num_seg)%GG_U(3,I,J)*channel%seg(num_seg)%U3DX(I,J,K)           &
              - (channel%seg(num_seg)%GG_V(2,I,J)*channel%seg(num_seg)%U3DY(I,J,K)          &
                +channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%U3DY(I+1,J,K)      &
                +channel%seg(num_seg)%GG_V(2,I,J-1)*channel%seg(num_seg)%U3DY(I,J-1,K)      &
                +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%U3DY(I+1,J-1,K)) &
                /4.0_r8
        ENDDO  
        channel%seg(num_seg)%U3DX_CO(I,J,1)   = channel%seg(num_seg)%U3DX_CO(I,J,2)
        channel%seg(num_seg)%U3DX_CO(I,J,nk3) = channel%seg(num_seg)%U3DX_CO(I,J,nk2)                          
       ENDDO
      ENDDO          

!---------------
! V-component
!---------------  
      DO J = mjm_a, mjp_c  
       DO I = mim_a, mip_a
        DO K = 3, nk2
         channel%seg(num_seg)%U3DY_CO(I,J,K) = &
                channel%seg(num_seg)%GG_V(1,I,J)*channel%seg(num_seg)%U3DY(I,J,K)           &
              - (channel%seg(num_seg)%GG_U(2,I,J)*channel%seg(num_seg)%U3DX(I,J,K)          &
                +channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%U3DX(I,J+1,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J)*channel%seg(num_seg)%U3DX(I-1,J,K)      &
                +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%U3DX(I-1,J+1,K)) &
                /4.0_r8
        ENDDO          
        channel%seg(num_seg)%U3DY_CO(I,J,1)   = channel%seg(num_seg)%U3DY_CO(I,J,2)
        channel%seg(num_seg)%U3DY_CO(I,J,nk3) = channel%seg(num_seg)%U3DY_CO(I,J,nk2)                          
       ENDDO
      ENDDO      
                  
!**********************************************************************   
      ENDDO  ! num_seg
!**********************************************************************

      ENDIF   ! UCHANGE

      END SUBROUTINE wind_3d_v0
      
END MODULE wind_module
