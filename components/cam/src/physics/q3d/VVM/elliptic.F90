module elliptic
! Solve 2D and 3D elliptic equations by iteration method

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk1,nk2,nhalo,nhalo_cal
USE constld, only: wrxmu,niterw,niterxy,dx,dy,dz, &
                   rho,rhoz,fnt,fnz,nomap,dxsq,dysq,dxdy

! Subroutines being called
USE bound_channel_module, only: bound_channel
USE bound_extra, only: bound_normal
USE halo_q,      only: halo_correc_q 
USE halo_z,      only: halo_correc_z

IMPLICIT NONE
PRIVATE

PUBLIC :: relax_3d, relax_2d

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE relax_3d  : Solve 3D elliptic equations by iteration method
! SUBROUTINE relax_2d  : Solve 2D elliptic equations by iteration method
! SUBROUTINE GAUSS_LD  : Gauss elimination
!-----------------------------------------------------------------------

!=======================================================================
      SUBROUTINE RELAX_3D (channel)
!=======================================================================    
! Solve the 3D elliptic equation for w by iteration method
!
! Temporary use of channel%seg(num_seg)%Z3DX0, channel%seg(num_seg)%Z3DY0, 
! and channel%seg(num_seg)%TERM1  
          
      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local variables
      REAL (kind=r8), DIMENSION(NK2) :: &
         AGAU,        & ! coefficient a in gaussian elimination
         BGAU,        & ! coefficient b in gaussian elimination
         CGAU           ! coefficient c in gaussian elimination         

      REAL (kind=r8) :: RHSV(NK2),ROUT(NK2)
      
      INTEGER :: I,J,K,ITER,KLOW,KLOWU,KLOWV
      INTEGER :: num_seg,mi1,mj1,mim_c,mip_c,mjm_c,mjp_c

!     Pre-process of iteration
!-------------------------------   
      DO num_seg = 1, 4
!------------------------------- 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim_c  = channel%seg(num_seg)%mim_c    
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm_c  = channel%seg(num_seg)%mjm_c   
      mjp_c  = channel%seg(num_seg)%mjp_c           
      
      ! Prepare vorticity components with topographic boundary conditions     
      DO J = mjm_c,mjp_c
      DO I = mim_c,mip_c

       KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
       DO K = KLOWV, NK2
        channel%seg(num_seg)%Z3DX0(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K)
       ENDDO
       DO K = 2, KLOWV-1
        channel%seg(num_seg)%Z3DX0(I,J,K) = &
          (channel%seg(num_seg)%W3D(I,J+1,K)                                     &
          -channel%seg(num_seg)%W3D(I,J,K))/(DY*channel%seg(num_seg)%RG_V(I,J))  &
        - (channel%seg(num_seg)%U3DY_CO(I,J,K+1)                                 &
          -channel%seg(num_seg)%U3DY_CO(I,J,K))*FNZ(K)/(DZ*channel%seg(num_seg)%RG_V(I,J))
       ENDDO

       KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)
       DO K = KLOWU,NK2
        channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K)
       ENDDO
       DO K = 2, KLOWU-1
        channel%seg(num_seg)%Z3DY0(I,J,K) = &
          (channel%seg(num_seg)%U3DX_CO(I,J,K+1)                                            &
          -channel%seg(num_seg)%U3DX_CO(I,J,K))*FNZ(K)/(DZ*channel%seg(num_seg)%RG_U(I,J))  &
        - (channel%seg(num_seg)%W3D(I+1,J,K)                                                &
          -channel%seg(num_seg)%W3D(I,J,K))/(DX*channel%seg(num_seg)%RG_U(I,J))
       ENDDO

      ENDDO
      ENDDO
      
      ! Calculate the forcing term
      DO K = 2, NK1
       DO J = 1, MJ1
        DO I = 1, MI1
         channel%seg(num_seg)%TERM1(I,J,K) = &
           (channel%seg(num_seg)%GG_U(1,I,J)*channel%seg(num_seg)%Z3DY0(I,J,K)                  &
           -channel%seg(num_seg)%GG_U(1,I-1,J)*channel%seg(num_seg)%Z3DY0(I-1,J,K))/DX          &        
          -(channel%seg(num_seg)%GG_V(3,I,J)*channel%seg(num_seg)%Z3DX0(I,J,K)                  &
           -channel%seg(num_seg)%GG_V(3,I,J-1)*channel%seg(num_seg)%Z3DX0(I,J-1,K))/DY          &
          -(channel%seg(num_seg)%GG_V(2,I+1,J)*channel%seg(num_seg)%Z3DX0(I+1,J,K)              &
           +channel%seg(num_seg)%GG_V(2,I+1,J-1)*channel%seg(num_seg)%Z3DX0(I+1,J-1,K)          &
           -channel%seg(num_seg)%GG_V(2,I-1,J)*channel%seg(num_seg)%Z3DX0(I-1,J,K)              &
           -channel%seg(num_seg)%GG_V(2,I-1,J-1)*channel%seg(num_seg)%Z3DX0(I-1,J-1,K))         &
          /(4.0_r8*DX)                                                                     &
          +(channel%seg(num_seg)%GG_U(2,I,J+1)*channel%seg(num_seg)%Z3DY0(I,J+1,K)              &
           +channel%seg(num_seg)%GG_U(2,I-1,J+1)*channel%seg(num_seg)%Z3DY0(I-1,J+1,K)          &
           -channel%seg(num_seg)%GG_U(2,I,J-1)*channel%seg(num_seg)%Z3DY0(I,J-1,K)              &
           -channel%seg(num_seg)%GG_U(2,I-1,J-1)*channel%seg(num_seg)%Z3DY0(I-1,J-1,K))         &
          /(4.0_r8*DY) 
        ENDDO
       ENDDO
      ENDDO 

!     linear extrapolation of initial guess
      DO K = 2, NK1
       DO J = mjm_c, mjp_c
        DO I = mim_c, mip_c
         channel%seg(num_seg)%Z3DX0(I,J,K) = &
           2.0_r8*channel%seg(num_seg)%W3D(I,J,K)-channel%seg(num_seg)%W3DNM1(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      DO K = 1, NK2
       DO J = mjm_c, mjp_c
        DO I = mim_c, mip_c
         channel%seg(num_seg)%W3DNM1(I,J,K) = channel%seg(num_seg)%W3D(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      DO K = 2, NK1
       DO J = mjm_c, mjp_c
        DO I = mim_c, mip_c
         channel%seg(num_seg)%W3D(I,J,K) = channel%seg(num_seg)%Z3DX0(I,J,K)
        ENDDO
       ENDDO
      ENDDO      

!-------------------------------   
      ENDDO  ! num_seg
!------------------------------- 

!******************************************
      DO ITER = 1, NITERW
!******************************************

!-------------------------------  
      DO num_seg = 1, 4
!------------------------------- 
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim_c  = channel%seg(num_seg)%mim_c    
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm_c  = channel%seg(num_seg)%mjm_c   
      mjp_c  = channel%seg(num_seg)%mjp_c         
          
      DO K = 2, NK1
       DO J =  mjm_c,mjp_c
        DO I = mim_c,mip_c
         channel%seg(num_seg)%Z3DY0(I,J,K) = channel%seg(num_seg)%W3D(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      
      DO J = 1, mj1
       DO I = 1, mi1
       KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)
       
        DO K = 2, NK1
         AGAU(K) = channel%seg(num_seg)%AGAU_CO(I,J,K-1)
         BGAU(K) = channel%seg(num_seg)%BGAU_CO(I,J,K-1)
         CGAU(K) = channel%seg(num_seg)%CGAU_CO(I,J,K-1)
        ENDDO
        
        ! Z3DY0: temporary value of w3d, term1: forcing
        DO K = KLOW, NK1
         RHSV(K) = WRXMU*channel%seg(num_seg)%Z3DY0(I,J,K) + channel%seg(num_seg)%TERM1(I,J,K) &
       + (channel%seg(num_seg)%RGG_U(1,I,J)*channel%seg(num_seg)%Z3DY0(I+1,J,K)           &
         +channel%seg(num_seg)%RGG_U(1,I-1,J)*channel%seg(num_seg)%Z3DY0(I-1,J,K))/DXSQ   &
       + (channel%seg(num_seg)%RGG_V(3,I,J)*channel%seg(num_seg)%Z3DY0(I,J+1,K)           &
         +channel%seg(num_seg)%RGG_V(3,I,J-1)*channel%seg(num_seg)%Z3DY0(I,J-1,K))/DYSQ   &
       + (channel%seg(num_seg)%RGG_V(2,I+1,J)                                             &
           *(channel%seg(num_seg)%Z3DY0(I+1,J+1,K)-channel%seg(num_seg)%Z3DY0(I+1,J,K))   &
         +channel%seg(num_seg)%RGG_V(2,I+1,J-1)                                           &
           *(channel%seg(num_seg)%Z3DY0(I+1,J,K)-channel%seg(num_seg)%Z3DY0(I+1,J-1,K))   &
         -channel%seg(num_seg)%RGG_V(2,I-1,J)                                             &
           *(channel%seg(num_seg)%Z3DY0(I-1,J+1,K)-channel%seg(num_seg)%Z3DY0(I-1,J,K))   &
         -channel%seg(num_seg)%RGG_V(2,I-1,J-1)                                           &
           *(channel%seg(num_seg)%Z3DY0(I-1,J,K)-channel%seg(num_seg)%Z3DY0(I-1,J-1,K))   &  
         +channel%seg(num_seg)%RGG_U(2,I,J+1)                                             &
           *(channel%seg(num_seg)%Z3DY0(I+1,J+1,K)-channel%seg(num_seg)%Z3DY0(I,J+1,K))   &
         +channel%seg(num_seg)%RGG_U(2,I-1,J+1)                                           &
           *(channel%seg(num_seg)%Z3DY0(I,J+1,K)-channel%seg(num_seg)%Z3DY0(I-1,J+1,K))   &
         -channel%seg(num_seg)%RGG_U(2,I,J-1)                                             &
           *(channel%seg(num_seg)%Z3DY0(I+1,J-1,K)-channel%seg(num_seg)%Z3DY0(I,J-1,K))   &
         -channel%seg(num_seg)%RGG_U(2,I-1,J-1)                                           &
           *(channel%seg(num_seg)%Z3DY0(I,J-1,K)-channel%seg(num_seg)%Z3DY0(I-1,J-1,K)))  &
         /(4.0_r8*DXDY)   
        ENDDO 
                
        CALL GAUSS_LD(KLOW,NK2,AGAU,BGAU,CGAU,RHSV,ROUT)
        
        DO K = 2, NK1
         channel%seg(num_seg)%W3D(I,J,K) = ROUT(K)/RHOZ(K)
        ENDDO
              
       ENDDO    !  i-loop
      ENDDO    !  j-loop  

!-------------------------------   
      ENDDO  ! num_seg 
!------------------------------- 
      CALL BOUND_NORMAL  (nhalo,channel,W3D=.TRUE.)
      
      IF (ITER.LT.NITERW) THEN
        CALL BOUND_CHANNEL (nhalo_cal,channel,W3D=.TRUE.)
        IF (.not.nomap) CALL HALO_CORREC_Q (nhalo_cal,channel,W3D=.TRUE.) 
      ELSE
        CALL BOUND_CHANNEL (nhalo,channel,W3D=.TRUE.)
        IF (.not.nomap) CALL HALO_CORREC_Q (nhalo,channel,W3D=.TRUE.)
      ENDIF 
      
!******************************************
      ENDDO    ! iteration-loop   
!******************************************

      CONTAINS

!=======================================================================      
      SUBROUTINE GAUSS_LD(KMIN,KP2,AN,BN,CN,ZM,PM)
!=======================================================================      
!
! bn(2)*pm(2)+cn(2)*pm(3)                            = zm(2)
! an(3)*pm(2)+bn(3)*pm(3)+cn(3)*pm(4)                = zm(3)
!           an(k)*pm(k-1)+bn(k)*pm(k)+cn(k)*pn(k+1)  = zm(k)
!
!            an(kp2-1)*pm(kp2-2)+bn(kp2-1)*pm(kp2-1) = zm(kp2-1)
!
!     tridiagonal system [an,bn,cn].pm = zm
!     loop 20: a=l.u, l=(modif-bn in diag, an in subdiag)
!                     u=(ones in diag, modif cn in supdiag)
!     loop 30: u.pm = g = zm \ l, stored it on pm
!     loop 40: pm = g \ u, stored it on pm

      integer, intent(in) :: kp2,kmin
      real (kind=r8), intent(in) :: zm(kp2),an(kp2),bn(kp2),cn(kp2)
      real (kind=r8), intent(out):: pm(kp2)

      real (kind=r8) :: bn_new(kp2),cn_new(kp2),pm_temp(kp2)
      integer :: k

      cn_new(kmin) = cn(kmin) / bn(kmin)
      do k = kmin+1, kp2-1
      bn_new(k) = bn(k) - an(k) * cn_new(k-1)
      cn_new(k) = cn(k) / bn_new(k)
      enddo

      pm_temp(kmin) = zm(kmin) / bn(kmin)

      do k = kmin+1, kp2-1
      pm_temp(k) = ( zm(k) - an(k) * pm_temp(k-1) ) / bn_new(k)
      enddo

      pm(kp2-1) = pm_temp(kp2-1)
      do k = kp2-2, kmin, -1
      pm(k) = pm_temp(k) - cn_new(k) * pm(k+1)
      enddo
      do k = 1, kmin-1
      pm(k) = 0.0_r8
      enddo

      end subroutine gauss_ld

      end subroutine relax_3d

!=======================================================================
      SUBROUTINE RELAX_2D (NTYPE,channel,NITER2D)
!=======================================================================   
!     Solve 2D elliptic equations by iteration method
!
!     Temporary use of channel%seg(num_seg)%Z3DX0, channel%seg(num_seg)%Z3DY0, 
!     and channel%seg(num_seg)%TERM1  
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NTYPE
!     ntype=1: PSI (z-point value)    ntype=2: CHI (q-point value)

      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER, INTENT(IN), OPTIONAL  :: NITER2D   ! number of iteration for 2D elliptic equation, 
                                                  ! which is different from the basic setting 
                                                  ! "niterxy" (test purpose)   
     
      ! Local variables
      LOGICAL :: USE_MU = .TRUE.
      REAL (kind=r8) :: COEF0  
      
      INTEGER :: I,J,ITER,NITER
      INTEGER :: num_seg,mi1,mj1,mim,mip,mjm,mjp,mim_c,mip_c,mjm_c,mjp_c
      
      IF (PRESENT(NITER2D)) THEN
         NITER = NITER2D
      ELSE
         NITER = NITERXY
      ENDIF

      IF (USE_MU) THEN
        COEF0 = WRXMU
      ELSE
        COEF0 = 0.0_r8
      ENDIF   

!===============================================================================         
      SELECT CASE (NTYPE)
!=============================================      

!=============================================      
      CASE(1)
!     PSI calculation (zeta-point)     
!============================================= 

!     Pre-process of iteration
!-------------------------------
      DO num_seg = 1, 4
!------------------------------- 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim    
      mip = channel%seg(num_seg)%mip  
    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm = channel%seg(num_seg)%mjm   
      mjp = channel%seg(num_seg)%mjp           

      ! Calculate the forcing term
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%TERM1(I,J,1) = channel%seg(num_seg)%RG_Z(I,J) &
                *(channel%seg(num_seg)%Z3DZ(I,J,NK2)-channel%seg(num_seg)%Z3DZ_BG(I,J,1))
       ENDDO
      ENDDO     

      ! linear extrapolation of initial guess (temporary use of z3dx0) 
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%Z3DX0(I,J,1) = 2.0_r8*channel%seg(num_seg)%PSI(I,J) &
                                           - channel%seg(num_seg)%PSINM1(I,J)
       ENDDO
      enddo

      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSINM1(I,J) = channel%seg(num_seg)%PSI(I,J) 
       ENDDO
      enddo 
      
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%PSI(I,J) = channel%seg(num_seg)%Z3DX0(I,J,1)
       ENDDO
      enddo            
!-------------------------------   
      ENDDO  ! num_seg
!------------------------------- 

!******************************************
      DO ITER = 1,NITER
!******************************************
!-------------------------------   
      DO num_seg = 1, 4
!------------------------------- 
      mi1   = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim_c = channel%seg(num_seg)%mim_c    
      mip_c = channel%seg(num_seg)%mip_c  
    
      mj1   = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm_c = channel%seg(num_seg)%mjm_c   
      mjp_c = channel%seg(num_seg)%mjp_c          
      
      DO J =  mjm_c,mjp_c
       DO I = mim_c,mip_c
         channel%seg(num_seg)%Z3DY0(I,J,1) = channel%seg(num_seg)%PSI(I,J)
       ENDDO
      ENDDO

      ! Z3DY0(I,J,1): temporary psi value, TERM1(I,J,1): forcing
      DO J = 1, mj1
       DO I = 1, mi1
         channel%seg(num_seg)%PSI(I,J) = COEF0*channel%seg(num_seg)%Z3DY0(I,J,1)         &
                                       - channel%seg(num_seg)%TERM1(I,J,1)               &
         + (channel%seg(num_seg)%RGG_V(1,I+1,J)*channel%seg(num_seg)%Z3DY0(I+1,J,1)      &
           +channel%seg(num_seg)%RGG_V(1,I,J)*channel%seg(num_seg)%Z3DY0(I-1,J,1))/DXSQ  &
         + (channel%seg(num_seg)%RGG_U(3,I,J+1)*channel%seg(num_seg)%Z3DY0(I,J+1,1)      &
           +channel%seg(num_seg)%RGG_U(3,I,J)*channel%seg(num_seg)%Z3DY0(I,J-1,1))/DYSQ  &
         + (channel%seg(num_seg)%RGG_U(2,I+1,J+1)*(channel%seg(num_seg)%Z3DY0(I+1,J+1,1) &
                                                  -channel%seg(num_seg)%Z3DY0(I+1,J,1))  &
           +channel%seg(num_seg)%RGG_U(2,I+1,J)*(channel%seg(num_seg)%Z3DY0(I+1,J,1)     &
                                                -channel%seg(num_seg)%Z3DY0(I+1,J-1,1))  &
           -channel%seg(num_seg)%RGG_U(2,I-1,J+1)*(channel%seg(num_seg)%Z3DY0(I-1,J+1,1) &
                                                  -channel%seg(num_seg)%Z3DY0(I-1,J,1))  &
           -channel%seg(num_seg)%RGG_U(2,I-1,J)*(channel%seg(num_seg)%Z3DY0(I-1,J,1)     &
                                                -channel%seg(num_seg)%Z3DY0(I-1,J-1,1))  & 
           +channel%seg(num_seg)%RGG_V(2,I+1,J+1)*(channel%seg(num_seg)%Z3DY0(I+1,J+1,1) &
                                                  -channel%seg(num_seg)%Z3DY0(I,J+1,1))  &
           +channel%seg(num_seg)%RGG_V(2,I,J+1)*(channel%seg(num_seg)%Z3DY0(I,J+1,1)     &
                                                -channel%seg(num_seg)%Z3DY0(I-1,J+1,1))  &
           -channel%seg(num_seg)%RGG_V(2,I+1,J-1)*(channel%seg(num_seg)%Z3DY0(I+1,J-1,1) &
                                                  -channel%seg(num_seg)%Z3DY0(I,J-1,1))  &
           -channel%seg(num_seg)%RGG_V(2,I,J-1)*(channel%seg(num_seg)%Z3DY0(I,J-1,1)     &
                                                -channel%seg(num_seg)%Z3DY0(I-1,J-1,1))) &
           /(4.0_r8*DXDY) 
       ENDDO         
      ENDDO
      
     DO J = 1, mj1
       DO I = 1, mi1
         channel%seg(num_seg)%PSI(I,J) = channel%seg(num_seg)%PSI(I,J) &
                                       /(COEF0+channel%seg(num_seg)%RX2D_C1(I,J))
       ENDDO         
      ENDDO
!-------------------------------   
      ENDDO  ! num_seg 
!------------------------------- 
      
      CALL BOUND_NORMAL  (nhalo,channel,PSI=.TRUE.)                      
      if (iter.lt.niter) then
        CALL BOUND_CHANNEL (nhalo_cal,channel,PSI=.TRUE.)
        if (.not.nomap) CALL HALO_CORREC_Z (nhalo_cal,channel,PSI=.TRUE.)
      else
        
        CALL BOUND_CHANNEL (nhalo,channel,PSI=.TRUE.)
        if (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,PSI=.TRUE.) 
      endif 
      
!******************************************
      ENDDO    ! iteration-loop   
!****************************************** 
         
!=============================================           
      CASE(2)
!     CHI calculation (q-point) 
!=============================================

!     Pre-process of iteration
!-------------------------------   
      DO num_seg = 1, 4
!------------------------------- 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim = channel%seg(num_seg)%mim    
      mip = channel%seg(num_seg)%mip  
    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm = channel%seg(num_seg)%mjm   
      mjp = channel%seg(num_seg)%mjp           

      ! Calculate the forcing term 
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%TERM1(I,J,1) = channel%seg(num_seg)%RG_T(I,J)*FNT(NK2)*RHOZ(NK1) &
                                           *(channel%seg(num_seg)%W3D(I,J,NK1) &
                                            -channel%seg(num_seg)%W3D_BG(I,J,NK1))/(RHO(NK2)*DZ)
       ENDDO
      ENDDO   
      
      ! linear extrapolation of initial guess           
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%Z3DX0(I,J,1) = 2.0_r8*channel%seg(num_seg)%CHI(I,J) &
                                           - channel%seg(num_seg)%CHINM1(I,J)
       ENDDO
      enddo
      
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%CHINM1(I,J) = channel%seg(num_seg)%CHI(I,J) 
       ENDDO
      enddo 
      
      DO J = mjm, mjp
       DO I = mim, mip
         channel%seg(num_seg)%CHI(I,J) = channel%seg(num_seg)%Z3DX0(I,J,1)
       ENDDO
      enddo            
!-------------------------------   
      ENDDO  ! num_seg
!------------------------------- 

!******************************************
      DO ITER = 1,NITER
!******************************************
!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1   = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim_c = channel%seg(num_seg)%mim_c    
      mip_c = channel%seg(num_seg)%mip_c  
    
      mj1   = channel%seg(num_seg)%mj1  ! y-size of channel segment    
      mjm_c = channel%seg(num_seg)%mjm_c   
      mjp_c = channel%seg(num_seg)%mjp_c          
      
      DO J =  mjm_c,mjp_c
       DO I = mim_c,mip_c
         channel%seg(num_seg)%Z3DY0(I,J,1) = channel%seg(num_seg)%CHI(I,J)
       ENDDO
      ENDDO
      
      ! Z3DY0(I,J,1): temporary chi value, TERM1(I,J,1): forcing 
      DO J = 1, mj1
       DO I = 1, mi1
         channel%seg(num_seg)%CHI(I,J) = COEF0*channel%seg(num_seg)%Z3DY0(I,J,1)           &
                                       - channel%seg(num_seg)%TERM1(I,J,1)                 &
         + (channel%seg(num_seg)%RGG_U(1,I,J)*channel%seg(num_seg)%Z3DY0(I+1,J,1)          &
           +channel%seg(num_seg)%RGG_U(1,I-1,J)*channel%seg(num_seg)%Z3DY0(I-1,J,1))/DXSQ  &
         + (channel%seg(num_seg)%RGG_V(3,I,J)*channel%seg(num_seg)%Z3DY0(I,J+1,1)          &
           +channel%seg(num_seg)%RGG_V(3,I,J-1)*channel%seg(num_seg)%Z3DY0(I,J-1,1))/DYSQ  &
         + (channel%seg(num_seg)%RGG_V(2,I+1,J)*(channel%seg(num_seg)%Z3DY0(I+1,J+1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I+1,J,1))            &
           +channel%seg(num_seg)%RGG_V(2,I+1,J-1)*(channel%seg(num_seg)%Z3DY0(I+1,J,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I+1,J-1,1))          &
           -channel%seg(num_seg)%RGG_V(2,I-1,J)*(channel%seg(num_seg)%Z3DY0(I-1,J+1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I-1,J,1))            &
           -channel%seg(num_seg)%RGG_V(2,I-1,J-1)*(channel%seg(num_seg)%Z3DY0(I-1,J,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I-1,J-1,1))          & 
           +channel%seg(num_seg)%RGG_U(2,I,J+1)*(channel%seg(num_seg)%Z3DY0(I+1,J+1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I,J+1,1))            &
           +channel%seg(num_seg)%RGG_U(2,I-1,J+1)*(channel%seg(num_seg)%Z3DY0(I,J+1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I-1,J+1,1))          &
           -channel%seg(num_seg)%RGG_U(2,I,J-1)*(channel%seg(num_seg)%Z3DY0(I+1,J-1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I,J-1,1))            &
           -channel%seg(num_seg)%RGG_U(2,I-1,J-1)*(channel%seg(num_seg)%Z3DY0(I,J-1,1)     &
                                          -channel%seg(num_seg)%Z3DY0(I-1,J-1,1)))         &
           /(4.0_r8*DXDY) 
       ENDDO         
      ENDDO
      
      DO J = 1,mj1
       DO I = 1,mi1
         channel%seg(num_seg)%CHI(I,J) = channel%seg(num_seg)%CHI(I,J) &
                                       /(COEF0+channel%seg(num_seg)%RX2D_C2(I,J))
       ENDDO         
      ENDDO  
       
!-------------------------------   
      ENDDO  ! num_seg 
!-------------------------------
      
      CALL BOUND_NORMAL  (nhalo,channel,CHI=.TRUE.)
      if (iter.lt.niter) then
        CALL BOUND_CHANNEL (nhalo_cal,channel,CHI=.TRUE.)
        if (.not.nomap) CALL HALO_CORREC_Q (nhalo_cal,channel,CHI=.TRUE.)
      else
        CALL BOUND_CHANNEL (nhalo,channel,CHI=.TRUE.)
        if (.not.nomap) CALL HALO_CORREC_Q (nhalo,channel,CHI=.TRUE.) 
      endif 
            
!******************************************
      ENDDO    ! iteration-loop   
!******************************************                  

!=============================================
      END SELECT
!===============================================================================      
                    
      end subroutine relax_2d
      
end module elliptic
