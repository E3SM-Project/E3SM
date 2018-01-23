MODULE vort_3d_module
! vorticity prediction

USE shr_kind_mod,   only: dbl_kind => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld, only: nk1,nk2,nk3,nhalo,nVGCM_seg,netsz
USE constld, only: d0_0,d1_0,d6_0,d8_0,d0_5,d0_25, &
                   a,b,dx,dy,dz,dt,aladv,rho,rhoz,fnt,fnz,fn1,fn2, &
                   dxsq,dysq,dxdy,crad,nomap,turbulence

! Subroutines being called
USE bound_channel_module, only: bound_channel  ! HALO DATA COMMUNICATION
USE bound_extra,    only: bound_normal,bound_vert
USE halo_vort,      only: vort_comm_pre,vort_comm_post
USE halo_z,         only: halo_correc_z
USE damping,        only: damping_vort

IMPLICIT NONE
PRIVATE

PUBLIC :: vort_3d,vort_3d_corec,zeta_diag

CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------------
! SUBROUTINE vort_3d       : predict vorticity
! SUBROUTINE vort_3d_corec : revise vorticity
! SUBROUTINE zeta_diag     : zeta (z3dz) diagnosis  at k < nk2
!
! SUBROUTINE RKSI_3D       : ksi  (z3dx) prediction
! SUBROUTINE RETA_3D       : eta  (z3dy) prediction
! SUBROUTINE ZETA_3D       : zeta (z3dz) prediction at k = nk2
! SUBROUTINE VADVEC_H1     : horizontal advection 1
! SUBROUTINE VADVEC_H2     : horizontal advection 2
! SUBROUTINE VADVEC_V      : vertical advection
! SUBROUTINE ABM_3D        : update z3dx and z3dy from dynamical processes
! SUBROUTINE ABM_3D_TURB   : update z3dx and z3dy from turbulence
! SUBROUTINE RELAXATION    : relax z3dx, z3dy, z3dz toward background (Q3D algorithm)
!            SUBROUTINE BOUND_G1_V, SUBROUTINE BOUND_G2_V, SUBROUTINE BOUND_G3_V
!            SUBROUTINE BOUND_G1_S, SUBROUTINE BOUND_G2_S, SUBROUTINE BOUND_G3_S               
! SUBROUTINE HDIFFUSION_SMAGORINSKY : nonlinear numerical diffusion (will not be used)
! SUBROUTINE LDIFFUSION2            : linear numerical diffusion
! FUNCTION   TURB_H        : called from LDIFFUSION2 & HDIFFUSION_SMAGORINSKY   
!----------------------------------------------------------------------------
!    D0_0  = 0.0_dbl_kind
!    D1_0  = 1.0_dbl_kind
!    D6_0  = 6.0_dbl_kind
!    D8_0  = 8.0_dbl_kind
!    D0_5  = 0.5_dbl_kind
!    D0_25 = 0.25_dbl_kind

!=======================================================================
   SUBROUTINE VORT_3D ( N1, N2, channel )
!=======================================================================
!  ALL CALCULATIONS ASSOCIATED WITH VORTICITY.

      INTEGER, INTENT(IN) :: N1,N2                ! previous and current time index (AB scheme)       
      type(channel_t), intent(inout) :: channel   ! channel data 
            
      ! Local variables
      INTEGER :: I, J, K, num_seg, mim, mip, mjm, mjp

!-------------------------------------------
!     CALCULATING THE VORTICITY TENDENCY
!-------------------------------------------

!************************************  
      DO num_seg = 1, 4
!************************************
      mim = channel%seg(num_seg)%mim     
      mip = channel%seg(num_seg)%mip  
    
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp 
      
      do k = 1, nk2
        do j = mjm, mjp
          do i = mim, mip
            channel%seg(num_seg)%Z3DX(I,J,K)=channel%seg(num_seg)%Z3DX(I,J,K)/RHOZ(K)
            channel%seg(num_seg)%Z3DY(I,J,K)=channel%seg(num_seg)%Z3DY(I,J,K)/RHOZ(K)
          enddo
        enddo
      enddo
      do k = 1, nk3
        do j = mjm, mjp
          do i = mim, mip
            channel%seg(num_seg)%Z3DZ(I,J,K) = &
           (channel%seg(num_seg)%Z3DZ(I,J,K)+channel%seg(num_seg)%FVAL(I,J))/RHO(K)
          enddo
        enddo
      enddo
!******************   
      ENDDO  
!******************     

      CALL RKSI_3D ( N1, N2, channel )
      CALL RETA_3D ( N1, N2, channel )

      CALL ZETA_3D ( N1, N2, channel )

!************************************  
      DO num_seg = 1, 4
!************************************
      mim = channel%seg(num_seg)%mim     
      mip = channel%seg(num_seg)%mip  
    
      mjm = channel%seg(num_seg)%mjm    
      mjp = channel%seg(num_seg)%mjp 
      
      do k = 1, nk2
        do j = mjm, mjp
          do i = mim, mip
            channel%seg(num_seg)%Z3DX(I,J,K)=channel%seg(num_seg)%Z3DX(I,J,K)*RHOZ(K)
            channel%seg(num_seg)%Z3DY(I,J,K)=channel%seg(num_seg)%Z3DY(I,J,K)*RHOZ(K)
          enddo
        enddo
      enddo
      do k = 1, nk3
        do j = mjm, mjp
          do i = mim, mip
            channel%seg(num_seg)%Z3DZ(I,J,K) = &
            channel%seg(num_seg)%Z3DZ(I,J,K)*RHO(K)-channel%seg(num_seg)%FVAL(I,J)
          enddo
        enddo
      enddo
!******************   
      ENDDO  
!****************** 

!--------------------------------------------------------
!     UPDATE: Z3DX & Z3DY
!--------------------------------------------------------
      CALL ABM_3D ( N1, N2, channel )

!-------------------
!     HALO UPDATE
!-------------------
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)

      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

!--------------------------------------------------------
!     UPDATE: Z3DZ at k = nk2
!--------------------------------------------------------
!************************************  
      DO num_seg = 1, 4
!************************************
      mi1 = channel%seg(num_seg)%mi1     
      mj1 = channel%seg(num_seg)%mj1    
        
      DO J = 1, mj1
      DO I = 1, mi1
       channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2)   &
                                          + A*channel%seg(num_seg)%FZTOP(I,J,N2) &
                                          + B*channel%seg(num_seg)%FZTOP(I,J,N1)
      ENDDO
      ENDDO
!******************   
      ENDDO  
!******************       

!----------------------------------
!     HALO UPDATE: Z3DZ at k = NK2
!     Good for halo = nhalo (>=2)
!----------------------------------
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)   
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.)

!--------------------------------------
!     UPDATE: Z3DZ at k < NK2
!--------------------------------------

      CALL ZETA_DIAG (channel)

   END SUBROUTINE vort_3d

!=======================================================================
   SUBROUTINE VORT_3D_COREC (channel)
!======================================================================= 
!  Update vorticity components (turbulence, upper layer damping, 
!  relaxation (coupling), and numerical diffusion (if necessary)
   
      type(channel_t), INTENT(INOUT) :: channel   ! channel data
      
      ! Local variables
      LOGICAL :: NDIFFU = .FALSE.
      LOGICAL :: LDIFFU = .TRUE.

      INTEGER :: I,J,K,num_seg,mi1,mj1
      
!=================================================================================
      IF (TURBULENCE) THEN
!=================================================================================
!     UPDATE: Z3DX & Z3DY
      CALL ABM_3D_TURB ( channel )

!**************************  
      DO num_seg = 1, 4
!**************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
!     UPDATE: Z3DZ at k=NK2
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                                           + DT*channel%seg(num_seg)%FZTOPB(I,J)
       ENDDO
      ENDDO

!     Diabatic_effect: turbulence (in terms of momentum)
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = 1, mi1
          channel%seg(num_seg)%FU_DIA(I,J,K) = channel%seg(num_seg)%FU_DIA(I,J,K) &
                                             + channel%seg(num_seg)%FU_TB(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      DO K=2,NK1
       DO J = 1, mj1
        DO I = 1, mi1
          channel%seg(num_seg)%FV_DIA(I,J,K) = channel%seg(num_seg)%FV_DIA(I,J,K) &
                                             + channel%seg(num_seg)%FV_TB(I,J,K)
        ENDDO
       ENDDO
      ENDDO
!************************  
      ENDDO  ! num_seg 
!************************
!========================
      ENDIF  
!========================

!=================================================================================
      IF (CRAD.NE.0.) CALL DAMPING_VORT (channel)
!=================================================================================

!=================================================================================
!     RELAXATION EFFECT (z3dx, z3dy, z3dz)
!=================================================================================
      CALL RELAXATION (channel)

!     -------------------
!     Z3DX & Z3DY
!     -------------------
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)

      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)      

!     -------------------
!     Z3DZ at k=NK2
!     -------------------
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)   
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.)      

!     -------------------
!     Z3DZ at k < NK2
!     -------------------
      CALL ZETA_DIAG (channel)

!=================================================================================
      IF (LDIFFU) THEN
!=================================================================================
      CALL LDIFFUSION2 (channel)

!     -------------------
!     Z3DX & Z3DY
!     -------------------
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)

      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)  

!     -------------------
!     Z3DZ at k=NK2
!     -------------------
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)   
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.)    
      
!     -------------------
!     Z3DZ at k < NK2
!     -------------------
      CALL ZETA_DIAG (channel) 

!==============
      ENDIF
!==============

!=================================================================================
      IF (NDIFFU) THEN  ! Smagorinsky type diffusion
!=================================================================================
      CALL HDIFFUSION_SMAGORINSKY (channel) 

!     -------------------
!     Z3DX & Z3DY
!     -------------------
      CALL BOUND_NORMAL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)

      IF (.NOT.nomap) CALL VORT_COMM_PRE (channel)
      CALL BOUND_CHANNEL (nhalo,channel,Z3DX=.TRUE.,Z3DY=.TRUE.)
      IF (.NOT.nomap) CALL VORT_COMM_POST (channel)

      CALL BOUND_VERT (channel,Z3DX=.TRUE.,Z3DY=.TRUE.)  

!     -------------------
!     Z3DZ at k=NK2
!     -------------------
      CALL BOUND_NORMAL  (nhalo,channel,Z3DZ=.TRUE.)
      CALL BOUND_SYNC (channel,Z3DZ=.TRUE.)   
      CALL BOUND_CHANNEL (nhalo,channel,Z3DZ=.TRUE.)
      IF (.not.nomap) CALL HALO_CORREC_Z (nhalo,channel,Z3DZ=.TRUE.) 

!     -------------------
!     Z3DZ at k < NK2
!     -------------------
      CALL ZETA_DIAG (channel) 

!==============
      ENDIF
!==============

   END SUBROUTINE vort_3d_corec

!=======================================================================
   SUBROUTINE RKSI_3D ( N1, N2, channel )
!=======================================================================
!  Predict ksi (z3dx), x-component of vorticity
!     mass-weight factors: fn1 & fn2
!     FN1(K)=RHO(K+1)*FNZ(K)/FNT(K+1)
!     FN2(K)=RHO(K)*FNZ(K)/FNT(K)
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: N1,N2                ! previous and current time index (AB scheme)       
      type(channel_t), intent(inout) :: channel   ! channel data 

      ! Local variables
      LOGICAL :: MODFLUX
      REAL (KIND=dbl_kind), DIMENSION(NK2) :: WWND,QVER,TERM1,TERM2,TERM3
      REAL (KIND=dbl_kind), DIMENSION(NK1) :: FINALZ
      
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: WND  
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: FINAL

      INTEGER :: I,J,K,L,klowv
      INTEGER :: num_seg,klowq_max
      INTEGER :: mi1,mim,mip,mim_a,mip_a,mim_c,mip_c,mim_ce
      INTEGER :: mj1,mjm,mjp,mjm_a,mjp_a,mjm_c,mjp_c,mjm_ce

      L = N2

!*******************************************************************   
      DO num_seg = 1, 4
!*******************************************************************  
      klowq_max = channel%seg(num_seg)%klowq_max_glob
      
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim  
      mip    = channel%seg(num_seg)%mip
      mim_a  = channel%seg(num_seg)%mim_a 
      mip_a  = channel%seg(num_seg)%mip_a   
      mim_c  = channel%seg(num_seg)%mim_c 
      mip_c  = channel%seg(num_seg)%mip_c 
      mim_ce = channel%seg(num_seg)%mim_ce 
       
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm    = channel%seg(num_seg)%mjm   
      mjp    = channel%seg(num_seg)%mjp 
      mjm_a  = channel%seg(num_seg)%mjm_a 
      mjp_a  = channel%seg(num_seg)%mjp_a  
      mjm_c  = channel%seg(num_seg)%mjm_c       
      mjp_c  = channel%seg(num_seg)%mjp_c 
      mjm_ce = channel%seg(num_seg)%mjm_ce         
        
      ALLOCATE(FINAL(mi1,mj1))
      ALLOCATE(WND(mim_a:mip_c,mjm_a:mjp_c))
      
!-----------------------
!     ZONAL ADVECTION
!-----------------------
      MODFLUX = .TRUE.
      ! k-loop for the lower atmosphere with topography
      DO K = 2, klowq_max-1

      DO J = 1, mj1
      DO I = mim_a, mip_c
      WND(I,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DX(I,J+1,K+1)   &
                               +channel%seg(num_seg)%U3DX(I,J,K+1))    &
                       +FN2(K)*(channel%seg(num_seg)%U3DX(I,J+1,K)     &
                               +channel%seg(num_seg)%U3DX(I,J,K)))     &
                      *channel%seg(num_seg)%RG_Z(I,J)
      ENDDO
      ENDDO

      CALL VADVEC_H1 (mi1,mim,mim_a,mim_c,mip,mip_c,                 &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,                 &
                      channel%seg(num_seg)%DM_KSI_XE(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_TE(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_XC(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_TW(:,:,K),         &
                      channel%seg(num_seg)%Z3DX(:,:,K),              &
                      WND,DX,DIRECTION='X',MODFLUX,FINAL)

      DO J = 1, mj1
      DO  I= 1, mi1
       channel%seg(num_seg)%FZX(I,J,K,L) = FINAL(I,J)
      ENDDO
      ENDDO

      ENDDO  

      MODFLUX = .FALSE.
      ! k-loop for free atmosphere
      DO K = klowq_max, NK1

      DO J = 1, mj1
      DO I = mim_a, mip_c
      WND(I,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DX(I,J+1,K+1)   &
                               +channel%seg(num_seg)%U3DX(I,J,K+1))    &
                       +FN2(K)*(channel%seg(num_seg)%U3DX(I,J+1,K)     &
                               +channel%seg(num_seg)%U3DX(I,J,K)))     &
                      *channel%seg(num_seg)%RG_Z(I,J)
      ENDDO
      ENDDO

      CALL VADVEC_H1 (mi1,mim,mim_a,mim_c,mip,mip_c,                 &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,                 &
                      channel%seg(num_seg)%DM_KSI_XE(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_TE(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_XC(:,:,K),         &
                      channel%seg(num_seg)%DM_KSI_TW(:,:,K),         &
                      channel%seg(num_seg)%Z3DX(:,:,K),              &
                      WND,DX,DIRECTION='X',MODFLUX,FINAL)

      DO J = 1, mj1
      DO I = 1, mi1
      channel%seg(num_seg)%FZX(I,J,K,L) = FINAL(I,J)
      ENDDO
      ENDDO

      ENDDO  

!--------------------------
!     MERIDIONAL ADVECTION
!--------------------------
      MODFLUX = .TRUE.
      ! k-loop for the lower atmosphere with topography
      DO K= 2, klowq_max-1

      DO J = mjm_a, mjp_c
      DO i = 1, mi1
      WND(i,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DY(I,J+1,K+1)   &
                               +channel%seg(num_seg)%U3DY(I,J,K+1))    &
                       +FN2(K)*(channel%seg(num_seg)%U3DY(I,J+1,K)     &
                               +channel%seg(num_seg)%U3DY(I,J,K)))     &
                      *channel%seg(num_seg)%RG_T(I,J+1)
      ENDDO
      ENDDO

      CALL VADVEC_H2 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_KSI_TN(:,:,K),  &
                      channel%seg(num_seg)%DM_KSI_YS(:,:,K),  &
                      channel%seg(num_seg)%Z3DX(:,:,K),       &
                      WND,DY,DIRECTION='Y',MODFLUX,FINAL)

      DO J = 1, mj1
      DO I = 1, mi1
      channel%seg(num_seg)%FZX(I,J,K,L) = channel%seg(num_seg)%FZX(I,J,K,L) &
                                        + FINAL(i,J)
      ENDDO
      ENDDO

      ENDDO  

      MODFLUX = .FALSE.
      ! k-loop for free atmosphere 
      DO K = klowq_max, NK1

      DO J = mjm_a, mjp_c
      DO i = 1, mi1
      WND(i,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DY(I,J+1,K+1)   &
                               +channel%seg(num_seg)%U3DY(I,J,K+1))    &
                       +FN2(K)*(channel%seg(num_seg)%U3DY(I,J+1,K)     &
                               +channel%seg(num_seg)%U3DY(I,J,K)))     &
                      *channel%seg(num_seg)%RG_T(I,J+1)
      ENDDO
      ENDDO

      CALL VADVEC_H2 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_KSI_TN(:,:,K),  &
                      channel%seg(num_seg)%DM_KSI_YS(:,:,K),  &
                      channel%seg(num_seg)%Z3DX(:,:,K),       &
                      WND,DY,DIRECTION='Y',MODFLUX,FINAL)

      DO J = 1, mj1
      DO I = 1, mi1
      channel%seg(num_seg)%FZX(I,J,K,L) = channel%seg(num_seg)%FZX(I,J,K,L) &
                                        + FINAL(i,J)
      ENDDO
      ENDDO

      ENDDO  ! k-loop for the upper atmosphere

      DO K= 2, nk1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%FZX(I,J,K,L) = &
            channel%seg(num_seg)%FZX(I,J,K,L)/channel%seg(num_seg)%RG_V(I,J)
        ENDDO
       ENDDO
      ENDDO

!--------------------------
!     VERTICAL ADVECTION
!--------------------------
      DO J = 1, mj1
      DO I = 1, mi1
      klowv = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K = klowv-1, nk1
        WWND(K) = D0_25*(RHOZ(K)*(channel%seg(num_seg)%W3D(I,J,K)       &
                                 +channel%seg(num_seg)%W3D(I,J+1,K))    &
                        +RHOZ(K+1)*(channel%seg(num_seg)%W3D(I,J,K+1)   &
                                   +channel%seg(num_seg)%W3D(I,J+1,K+1)))
      ENDDO
      DO K = klowv-1, nk2
        QVER(K) = channel%seg(num_seg)%Z3DX(I,J,K)
      ENDDO

      CALL VADVEC_V (WWND,QVER,FNZ,nk1,nk2,nk3,klowv-1,DZ,FINALZ)

      DO K = klowv, nk1
        channel%seg(num_seg)%FZX(I,J,K,L) = channel%seg(num_seg)%FZX(I,J,K,L) &
                                          + FINALZ(K)
      ENDDO
      DO K = 1, klowv-1
        channel%seg(num_seg)%FZX(I,J,K,L) = D0_0
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!=========================================
      DO J = 1, mj1
      DO I = 1, mi1
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)
!=========================================

!-------------------------------------
!     STRETCHING TERM
!     (3.25) in Jung & Arakawa (2005)
!-------------------------------------
      DO K = klowv, nk1
      TERM1(K) = (FN1(K)*(channel%seg(num_seg)%U3DX(I,J+1,K+1)     &
                         -channel%seg(num_seg)%U3DX(I-1,J+1,K+1))  &
                 +FN2(K)*(channel%seg(num_seg)%U3DX(I,J+1,K)       &
                         -channel%seg(num_seg)%U3DX(I-1,J+1,K)))   &
                 *(channel%seg(num_seg)%Z3DX(I,J,K)                &
                         +channel%seg(num_seg)%Z3DX(I,J+1,K))      &
               + (FN1(K)*(channel%seg(num_seg)%U3DX(I,J,K+1)       &
                         -channel%seg(num_seg)%U3DX(I-1,J,K+1))    &
                 +FN2(K)*(channel%seg(num_seg)%U3DX(I,J,K)         &
                         -channel%seg(num_seg)%U3DX(I-1,J,K)))     &
                 *(channel%seg(num_seg)%Z3DX(I,J,K)                &
                         +channel%seg(num_seg)%Z3DX(I,J-1,K))                   
      ENDDO
      DO K = klowv, nk1
       TERM1(K) = TERM1(K)/(D8_0*DX)
      ENDDO

!-------------------------------------------
!     TWISTING TERMS & f-TERMS
!     (Deformations are not used anymore)
!-------------------------------------------
      DO K = klowv, nk1
      TERM2(K) = (channel%seg(num_seg)%Z3DY(I,J+1,K)                        &
                                      +channel%seg(num_seg)%Z3DY(I,J,K))    &
                *(FN2(K)*(channel%seg(num_seg)%U3DX(I,J+1,K)                &
                         -channel%seg(num_seg)%U3DX(I,J,K))                 &
                 +FN1(K)*(channel%seg(num_seg)%U3DX(I,J+1,K+1)              &
                         -channel%seg(num_seg)%U3DX(I,J,K+1)))              &
               + (channel%seg(num_seg)%Z3DY(I-1,J+1,K)                      &
                                      +channel%seg(num_seg)%Z3DY(I-1,J,K))  &
                *(FN2(K)*(channel%seg(num_seg)%U3DX(I-1,J+1,K)              &
                         -channel%seg(num_seg)%U3DX(I-1,J,K))               &
                 +FN1(K)*(channel%seg(num_seg)%U3DX(I-1,J+1,K+1)            &
                         -channel%seg(num_seg)%U3DX(I-1,J,K+1)))  
      ENDDO
      DO K = klowv, nk1
       TERM2(K) = TERM2(K)/(D8_0*DY)
      ENDDO

      DO K = klowv, nk1
      TERM3(K) = (channel%seg(num_seg)%U3DX(I-1,J,K+1)   &
                 -channel%seg(num_seg)%U3DX(I-1,J,K)     &
                 +channel%seg(num_seg)%U3DX(I-1,J+1,K+1) &
                 -channel%seg(num_seg)%U3DX(I-1,J+1,K))  &
                *(channel%seg(num_seg)%Z3DZ(I-1,J,K)     &
                                      +channel%seg(num_seg)%Z3DZ(I-1,J,K+1))  &
               + (channel%seg(num_seg)%U3DX(I,J,K+1)     &
                 -channel%seg(num_seg)%U3DX(I,J,K)       &
                 +channel%seg(num_seg)%U3DX(I,J+1,K+1)   &
                 -channel%seg(num_seg)%U3DX(I,J+1,K))    &
                *(channel%seg(num_seg)%Z3DZ(I,J,K)       &
                                      +channel%seg(num_seg)%Z3DZ(I,J,K+1))
      ENDDO
      DO K = klowv, nk1
       TERM3(K) = RHOZ(K)*FNZ(K)*TERM3(K)/(D8_0*DZ)
      ENDDO

!-------------------------------------------
!     TOTAL SUM
!-------------------------------------------
      DO K = klowv, nk1
       channel%seg(num_seg)%FZX(I,J,K,L) = channel%seg(num_seg)%FZX(I,J,K,L) &
                                         + TERM1(K)+TERM2(K)+TERM3(K)
      ENDDO

!=========================================
      ENDDO  ! I-loop
      ENDDO  ! J-loop
!=========================================

      DEALLOCATE(FINAL)
      DEALLOCATE(WND)
!*******************************   
      ENDDO  ! num_seg 
!*******************************  
         
   END SUBROUTINE rksi_3d

!=======================================================================
   SUBROUTINE RETA_3D ( N1, N2, channel )
!=======================================================================
!  Predict eta (z3dy)      
!     mass-weight factors: fn1 & fn2
!     FN1(K)=RHO(K+1)*FNZ(K)/FNT(K+1)
!     FN2(K)=RHO(K)*FNZ(K)/FNT(K)
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: N1,N2                ! previous and current time index (AB scheme)       
      type(channel_t), intent(inout) :: channel   ! channel data 
      
! Local variables
      LOGICAL :: MODFLUX
      REAL (KIND=dbl_kind), DIMENSION(NK2) :: WWND,QVER,TERM1,TERM2,TERM3
      REAL (KIND=dbl_kind), DIMENSION(NK1) :: FINALZ
      
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: WND  
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: FINAL
      
      INTEGER :: I,J,K,L,klowu
      INTEGER :: num_seg,klowq_max
      INTEGER :: mi1,mim,mip,mim_a,mip_a,mim_c,mip_c,mim_ce
      INTEGER :: mj1,mjm,mjp,mjm_a,mjp_a,mjm_c,mjp_c,mjm_ce

      L = N2

!*******************************************************************   
      DO num_seg = 1, 4
!*******************************************************************  
      klowq_max = channel%seg(num_seg)%klowq_max_glob
      
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mim    = channel%seg(num_seg)%mim  
      mip    = channel%seg(num_seg)%mip
      mim_a  = channel%seg(num_seg)%mim_a 
      mip_a  = channel%seg(num_seg)%mip_a   
      mim_c  = channel%seg(num_seg)%mim_c 
      mip_c  = channel%seg(num_seg)%mip_c 
      mim_ce = channel%seg(num_seg)%mim_ce 
       
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      mjm    = channel%seg(num_seg)%mjm   
      mjp    = channel%seg(num_seg)%mjp 
      mjm_a  = channel%seg(num_seg)%mjm_a 
      mjp_a  = channel%seg(num_seg)%mjp_a  
      mjm_c  = channel%seg(num_seg)%mjm_c       
      mjp_c  = channel%seg(num_seg)%mjp_c 
      mjm_ce = channel%seg(num_seg)%mjm_ce         
        
      ALLOCATE(FINAL(mi1,mj1))
      ALLOCATE(WND(mim_a:mip_c,mjm_a:mjp_c))
      
!-----------------------
!     ZONAL ADVECTION
!-----------------------
      MODFLUX = .TRUE.
      ! k-loop for the lower atmosphere with topography
      DO K = 2, klowq_max-1

      DO j = 1, mj1
      DO I = mim_a, mip_c
      WND(I,j) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DX(I+1,J,K+1)  &
                               +channel%seg(num_seg)%U3DX(I,J,K+1))   &
                       +FN2(K)*(channel%seg(num_seg)%U3DX(I+1,J,K)    &
                               +channel%seg(num_seg)%U3DX(I,J,K)))    &
                      *channel%seg(num_seg)%RG_T(I+1,J)
      ENDDO
      ENDDO

      CALL VADVEC_H2 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_ETA_TE(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_XW(:,:,K),  &
                      channel%seg(num_seg)%Z3DY(:,:,K),       &
                      WND,DX,DIRECTION='X',MODFLUX,FINAL)

      DO J = 1, mj1
      DO I = 1, mi1
       channel%seg(num_seg)%FZY(I,J,K,L) = FINAL(I,j)
      ENDDO
      ENDDO

      ENDDO  

      MODFLUX = .FALSE.
      ! k-loop for free atmosphere
      DO K = klowq_max, nk1

      DO j = 1, mj1
      DO I = mim_a, mip_c
      WND(I,j) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DX(I+1,J,K+1)  &
                               +channel%seg(num_seg)%U3DX(I,J,K+1))   &
                       +FN2(K)*(channel%seg(num_seg)%U3DX(I+1,J,K)    &
                               +channel%seg(num_seg)%U3DX(I,J,K)))    &
                      *channel%seg(num_seg)%RG_T(I+1,J)
      ENDDO
      ENDDO

      CALL VADVEC_H2 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_ETA_TE(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_XW(:,:,K),  &
                      channel%seg(num_seg)%Z3DY(:,:,K),       &
                      WND,DX,DIRECTION='X',MODFLUX,FINAL)

      DO J = 1, mj1
      DO I = 1, mi1
       channel%seg(num_seg)%FZY(I,J,K,L) = FINAL(I,j)
      ENDDO
      ENDDO

      ENDDO  ! k-loop for the upper atmosphere

!--------------------------
!     MERIDIONAL ADVECTION
!--------------------------
      MODFLUX=.TRUE.
      ! k-loop for the lower atmosphere with topography
      DO K= 2, klowq_max-1

      DO J = mjm_a, mjp_c
      DO i = 1, mi1
      WND(i,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DY(I+1,J,K+1)  &
                               +channel%seg(num_seg)%U3DY(I,J,K+1))   &
                       +FN2(K)*(channel%seg(num_seg)%U3DY(I+1,J,K)    &
                               +channel%seg(num_seg)%U3DY(I,J,K)))    &
                      *channel%seg(num_seg)%RG_Z(I,J)
      ENDDO
      ENDDO

      CALL VADVEC_H1 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_ETA_YN(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_TN(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_YC(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_TS(:,:,K),  &
                      channel%seg(num_seg)%Z3DY(:,:,K),       &
                      WND,DY,DIRECTION='Y',MODFLUX,FINAL)

      DO J = 1, mj1
      DO i = 1, mi1
       channel%seg(num_seg)%FZY(I,J,K,L) = channel%seg(num_seg)%FZY(I,J,K,L) &
                                         + FINAL(i,J)
      ENDDO
      ENDDO

      ENDDO  

      MODFLUX = .FALSE.
      ! k-loop for free atmosphere 
      DO K= klowq_max, nk1

      DO J = mjm_a, mjp_c
      DO i = 1, mi1
      WND(i,J) = D0_25*(FN1(K)*(channel%seg(num_seg)%U3DY(I+1,J,K+1)  &
                               +channel%seg(num_seg)%U3DY(I,J,K+1))   &
                       +FN2(K)*(channel%seg(num_seg)%U3DY(I+1,J,K)    &
                               +channel%seg(num_seg)%U3DY(I,J,K)))    &
                      *channel%seg(num_seg)%RG_Z(I,J)
      ENDDO
      ENDDO

      CALL VADVEC_H1 (mi1,mim,mim_a,mim_c,mip,mip_c,          &
                      mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,          &
                      channel%seg(num_seg)%DM_ETA_YN(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_TN(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_YC(:,:,K),  &
                      channel%seg(num_seg)%DM_ETA_TS(:,:,K),  &
                      channel%seg(num_seg)%Z3DY(:,:,K),       &
                      WND,DY,DIRECTION='Y',MODFLUX,FINAL)

      DO J = 1, mj1
      DO i = 1, mi1
       channel%seg(num_seg)%FZY(I,J,K,L) = channel%seg(num_seg)%FZY(I,J,K,L) &
                                         + FINAL(i,J)
      ENDDO
      ENDDO

      ENDDO 

      DO K= 2, nk1
       DO J = 1, mj1
        DO I = 1, mi1
         channel%seg(num_seg)%FZY(I,J,K,L) = &
                channel%seg(num_seg)%FZY(I,J,K,L)/channel%seg(num_seg)%RG_U(I,J)
        ENDDO
       ENDDO
      ENDDO

!--------------------------
!     VERTICAL ADVECTION
!--------------------------
      DO J = 1, mj1
      DO I = 1, mi1
      klowu = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = klowu-1, nk1
        WWND(K) = D0_25*(RHOZ(K)*(channel%seg(num_seg)%W3D(I+1,J,K)     &
                                 +channel%seg(num_seg)%W3D(I,J,K))      &
                        +RHOZ(K+1)*(channel%seg(num_seg)%W3D(I+1,J,K+1) &
                                   +channel%seg(num_seg)%W3D(I,J,K+1)))              
      ENDDO
      DO K = klowu-1, nk2
        QVER(K) = channel%seg(num_seg)%Z3DY(I,J,K)
      ENDDO

      CALL VADVEC_V (WWND,QVER,FNZ,NK1,NK2,NK3,KLOWU-1,DZ,FINALZ)

      DO K = klowu, nk1
        channel%seg(num_seg)%FZY(I,J,K,L) = channel%seg(num_seg)%FZY(I,J,K,L) &
                                          + FINALZ(K)
      ENDDO
      DO K = 1, klowu-1
        channel%seg(num_seg)%FZY(I,J,K,L) = D0_0
      ENDDO

      ENDDO  ! I-loop
      ENDDO  ! J-loop

!=========================================
      DO J = 1, mj1
      DO I = 1, mi1
      klowu = channel%seg(num_seg)%KLOWU_IJ(I,J)
!=========================================

!-------------------------
!     STRETCHING TERM
!-------------------------
      DO K = klowu, nk1
      TERM1(K) = (FN1(K)*(channel%seg(num_seg)%U3DY(I,J,K+1)             &
                         -channel%seg(num_seg)%U3DY(I,J-1,K+1))          &
                 +FN2(K)*(channel%seg(num_seg)%U3DY(I,J,K)               &
                         -channel%seg(num_seg)%U3DY(I,J-1,K)))           &
                *(channel%seg(num_seg)%Z3DY(I-1,J,K)                     &
                                      +channel%seg(num_seg)%Z3DY(I,J,K)) &
                 +(FN1(K)*(channel%seg(num_seg)%U3DY(I+1,J,K+1)          &
                          -channel%seg(num_seg)%U3DY(I+1,J-1,K+1))       &
                  +FN2(K)*(channel%seg(num_seg)%U3DY(I+1,J,K)            &
                          -channel%seg(num_seg)%U3DY(I+1,J-1,K)))        &
                *(channel%seg(num_seg)%Z3DY(I,J,K)                       &
                                      +channel%seg(num_seg)%Z3DY(I+1,J,K))
      ENDDO
      DO K = klowu, nk1
       TERM1(K) = TERM1(K)/(D8_0*DY)
      ENDDO

!-------------------------------------------
!     TWISTING TERMS & f-TERMS
!     (Deformations are not used anymore)
!-------------------------------------------
      DO K = klowu, nk1
       TERM2(K) = (channel%seg(num_seg)%Z3DX(I+1,J,K)                      &
                                       +channel%seg(num_seg)%Z3DX(I,J,K))  &
                 *(FN2(K)*(channel%seg(num_seg)%U3DY(I+1,J,K)              &
                          -channel%seg(num_seg)%U3DY(I,J,K))               &
                  +FN1(K)*(channel%seg(num_seg)%U3DY(I+1,J,K+1)            &
                          -channel%seg(num_seg)%U3DY(I,J,K+1)))            &
                + (channel%seg(num_seg)%Z3DX(I+1,J-1,K)                    &
                                      +channel%seg(num_seg)%Z3DX(I,J-1,K)) &
                 *(FN2(K)*(channel%seg(num_seg)%U3DY(I+1,J-1,K)            &
                          -channel%seg(num_seg)%U3DY(I,J-1,K))             &
                  +FN1(K)*(channel%seg(num_seg)%U3DY(I+1,J-1,K+1)          &
                          -channel%seg(num_seg)%U3DY(I,J-1,K+1)))  
      ENDDO
      DO K = KLOWU,NK1
       TERM2(K) = TERM2(K)/(D8_0*DX)
      ENDDO

      DO K = klowu, nk1
       TERM3(K) = (channel%seg(num_seg)%U3DY(I,J-1,K+1)    &
                  -channel%seg(num_seg)%U3DY(I,J-1,K)      &
                  +channel%seg(num_seg)%U3DY(I+1,J-1,K+1)  &
                  -channel%seg(num_seg)%U3DY(I+1,J-1,K))   &
                 *(channel%seg(num_seg)%Z3DZ(I,J-1,K)      &
                                       +channel%seg(num_seg)%Z3DZ(I,J-1,K+1))  &
                + (channel%seg(num_seg)%U3DY(I,J,K+1)      &
                  -channel%seg(num_seg)%U3DY(I,J,K)        &
                  +channel%seg(num_seg)%U3DY(I+1,J,K+1)    &
                  -channel%seg(num_seg)%U3DY(I+1,J,K))     &
                 *(channel%seg(num_seg)%Z3DZ(I,J,K)        &
                                       +channel%seg(num_seg)%Z3DZ(I,J,K+1))
      ENDDO
      DO K = klowu, nk1
       TERM3(K) = RHOZ(K)*FNZ(K)*TERM3(K)/(8.*DZ)
      ENDDO

!-------------------------------------------
!     TOTAL EFFECT
!-------------------------------------------
      DO K = klowu, nk1
       channel%seg(num_seg)%FZY(I,J,K,L) = channel%seg(num_seg)%FZY(I,J,K,L) &
                                         + TERM1(K)+TERM2(K)+TERM3(K)
      ENDDO

!=========================================
      ENDDO  ! I-loop
      ENDDO  ! J-loop
!=========================================

      DEALLOCATE(FINAL)
      DEALLOCATE(WND)
      
!*******************************   
      ENDDO  ! num_seg 
!*******************************  

   END SUBROUTINE reta_3d

!=======================================================================
   SUBROUTINE ZETA_3D (N1, N2, channel)
!=======================================================================
!  Predict zeta (z3dz)
      
      INTEGER, INTENT(IN) :: N1,N2                ! previous and current time index (AB scheme)       
      type(channel_t), intent(inout) :: channel   ! channel data 
      
      ! Local variables
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: FLX
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: UPi,UMi     
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: UPSRi,UMSRi
      REAL (KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: TEMP

      REAL (KIND=dbl_kind) :: UPi_V,UPSRi_V,UPSR2i_V
      REAL (KIND=dbl_kind) :: FACTOR1,FACTOR2,TERM1,TERM2,TERM3
      
      INTEGER :: I,J,K,L
      INTEGER :: num_seg,mi1,mim_a,mip_c,mj1,mjm_a,mjp_c

      L = N2

!*******************************************************************   
      DO num_seg = 1, 4
!*******************************************************************  
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment  
      mim_a  = channel%seg(num_seg)%mim_a    
      mip_c  = channel%seg(num_seg)%mip_c 
      
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment   
      mjm_a  = channel%seg(num_seg)%mjm_a        
      mjp_c  = channel%seg(num_seg)%mjp_c      
      
      ALLOCATE(FLX(0:mi1,0:mj1))
      ALLOCATE(UPi(mim_a:mip_c,mjm_a:mjp_c))
      ALLOCATE(UMi(mim_a:mip_c,mjm_a:mjp_c))
      ALLOCATE(UPSRi(mim_a:mip_c,mjm_a:mjp_c))
      ALLOCATE(UMSRi(mim_a:mip_c,mjm_a:mjp_c))
      ALLOCATE(TEMP(mim_a:mip_c,mjm_a:mjp_c))
      
!===================
!     ADVECTION
!===================
! Zonal advection

      DO J = 1, mj1
       DO i = mim_a, mip_c
        TEMP(I,J) = D0_25*RHO(NK2)*(channel%seg(num_seg)%U3DX(I,J+1,NK2)   &
                                   +channel%seg(num_seg)%U3DX(I,J,NK2)     &
                                   +channel%seg(num_seg)%U3DX(I+1,J+1,NK2) &
                                   +channel%seg(num_seg)%U3DX(I+1,J,NK2))  &
                                  *channel%seg(num_seg)%RG_V(I+1,J)
       ENDDO 
      ENDDO

      DO J = 1, mj1
       DO i = mim_a, mip_c
        UPi(I,j) = D0_5*(TEMP(I,J)+ABS(TEMP(I,J)))
        UMi(I,j) = D0_5*(TEMP(I,J)-ABS(TEMP(I,J)))
       ENDDO
      ENDDO      

      DO J = 1, mj1
       DO i = mim_a, mip_c
        UPSRi(I,j) = SQRT(UPi(I,j)) 
        UMSRi(I,j) = SQRT(ABS(UMi(I,j)))
       ENDDO
      ENDDO

      DO J = 1, mj1
       DO i = 0, mi1
        FLX(I,j) = D0_5*TEMP(I,j)*(channel%seg(num_seg)%Z3DZ(I+1,J,NK2)   &
                                  +channel%seg(num_seg)%Z3DZ(I,J,NK2))    &
      - ALADV*(UPi(I,j)*(channel%seg(num_seg)%Z3DZ(I+1,J,NK2)                     &
                        -channel%seg(num_seg)%Z3DZ(I,J,NK2))                      &
              -UPSRi(I,j)*UPSRi(I-1,j)*(channel%seg(num_seg)%Z3DZ(I,J,NK2)        &
                                       -channel%seg(num_seg)%Z3DZ(I-1,J,NK2))     &
              +UMi(I,j)*(channel%seg(num_seg)%Z3DZ(I,J,NK2)                       &
                        -channel%seg(num_seg)%Z3DZ(I+1,J,NK2))                    &
              +UMSRi(I,j)*UMSRi(I+1,j)*(channel%seg(num_seg)%Z3DZ(I+1,J,NK2)      &
                                -channel%seg(num_seg)%Z3DZ(I+2,J,NK2)))/D6_0
       ENDDO
      ENDDO

      DO J = 1, mj1
       DO i = 1, mi1
        channel%seg(num_seg)%FZTOP(I,J,L) = &
                 -(FLX(I,j)-FLX(I-1,j))/(channel%seg(num_seg)%RG_Z(I,J)*DX)
       ENDDO
      ENDDO

! Meridional advection

      DO j = mjm_a, mjp_c
       DO I = 1, mi1
        TEMP(i,J) = D0_25*RHO(NK2)*(channel%seg(num_seg)%U3DY(I+1,J,NK2)   &
                                   +channel%seg(num_seg)%U3DY(I,J,NK2)     &
                                   +channel%seg(num_seg)%U3DY(I+1,J+1,NK2) &
                                   +channel%seg(num_seg)%U3DY(I,J+1,NK2))  &
                                  *channel%seg(num_seg)%RG_U(I,J+1)
       ENDDO
      ENDDO

      DO j = mjm_a, mjp_c
       DO I = 1, mi1
        UPi(i,J) = D0_5*(TEMP(i,J)+ABS(TEMP(i,J)))
        UMi(i,J) = D0_5*(TEMP(i,J)-ABS(TEMP(i,J)))
       ENDDO
      ENDDO

      DO j = mjm_a, mjp_c
       DO I = 1, mi1
        UPSRi(i,J) = SQRT(UPi(i,J))
        UMSRi(i,J) = SQRT(ABS(UMi(i,J)))
       ENDDO
      ENDDO

      DO j = 0, mj1
       DO I = 1, mi1
        FLX(i,J) = D0_5*TEMP(i,J)*(channel%seg(num_seg)%Z3DZ(I,J+1,NK2)  &
                                  +channel%seg(num_seg)%Z3DZ(I,J,NK2))   &
    -ALADV*(UPi(i,J)*(channel%seg(num_seg)%Z3DZ(I,J+1,NK2)                       &
                     -channel%seg(num_seg)%Z3DZ(I,J,NK2))                        &
           -UPSRi(i,J)*UPSRi(i,J-1)*(channel%seg(num_seg)%Z3DZ(I,J,NK2)          &
                                    -channel%seg(num_seg)%Z3DZ(I,J-1,NK2))       &
           +UMi(i,J)*(channel%seg(num_seg)%Z3DZ(I,J,NK2)                         &
                     -channel%seg(num_seg)%Z3DZ(I,J+1,NK2))                      &
           +UMSRi(i,J)*UMSRi(i,J+1)*(channel%seg(num_seg)%Z3DZ(I,J+1,NK2)        &
                                    -channel%seg(num_seg)%Z3DZ(I,J+2,NK2)))/D6_0
       ENDDO
      ENDDO

      DO j = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%FZTOP(I,J,L) = channel%seg(num_seg)%FZTOP(I,J,L) &
                    -(FLX(i,J)-FLX(i,J-1))/(channel%seg(num_seg)%RG_Z(I,J)*DY)
       ENDDO
      ENDDO

! Vertical advection

      DO J = 1, mj1
       DO I = 1, mi1
        TERM1 = D0_25*RHOZ(NK1)*(channel%seg(num_seg)%W3D(I,J,NK1)   &
                                +channel%seg(num_seg)%W3D(I+1,J,NK1) &
                                +channel%seg(num_seg)%W3D(I,J+1,NK1) &
                                +channel%seg(num_seg)%W3D(I+1,J+1,NK1))
        TERM2 = D0_25*RHOZ(NK2-2)*(channel%seg(num_seg)%W3D(I,J,NK2-2)    &
                                  +channel%seg(num_seg)%W3D(I+1,J,NK2-2)  &
                                  +channel%seg(num_seg)%W3D(I,J+1,NK2-2)  &
                                  +channel%seg(num_seg)%W3D(I+1,J+1,NK2-2))

        UPi_V    = D0_5*(TERM1+ABS(TERM1))
        UPSRi_V  = SQRT(UPi_V)
        UPSR2i_V = SQRT(D0_5*(TERM2+ABS(TERM2)))
      
        IF(TERM1.GE.D0_0) THEN
         FLX(i,J) = D0_5*TERM1*(channel%seg(num_seg)%Z3DZ(I,J,NK2)   &
                               +channel%seg(num_seg)%Z3DZ(I,J,NK1))  &
                  -ALADV*(UPi_V*(channel%seg(num_seg)%Z3DZ(I,J,NK2)          &
                                   -channel%seg(num_seg)%Z3DZ(I,J,NK1))      &
           -UPSRi_V*UPSR2i_V*(channel%seg(num_seg)%Z3DZ(I,J,NK1)             &
                                   -channel%seg(num_seg)%Z3DZ(I,J,NK1-1)))/D6_0
        ELSE
         FLX(i,J) = D0_5*TERM1*(channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                               +channel%seg(num_seg)%Z3DZ(I,J,NK1))
        END IF
       ENDDO
      ENDDO

      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%FZTOP(I,J,L) = channel%seg(num_seg)%FZTOP(I,J,L) &
                                          + FLX(i,J)*FNT(NK2)/DZ
!       Using W(NK2)=0
       ENDDO
      ENDDO

!======================
!     STRETCHING TERM
!     (f-term is included)
!======================
      DO J = 1, mj1
       DO I = 1, mi1
        TERM1 = &
        (channel%seg(num_seg)%W3D(I,J,NK1)+channel%seg(num_seg)%W3D(I,J+1,NK1))     &
       *(channel%seg(num_seg)%Z3DZ(I-1,J,NK2)+channel%seg(num_seg)%Z3DZ(I,J,NK2))   &
      + (channel%seg(num_seg)%W3D(I+1,J,NK1)+channel%seg(num_seg)%W3D(I+1,J+1,NK1)) &
       *(channel%seg(num_seg)%Z3DZ(I,J,NK2)+channel%seg(num_seg)%Z3DZ(I+1,J,NK2))
        
        channel%seg(num_seg)%FZTOP(I,J,L) = channel%seg(num_seg)%FZTOP(I,J,L) &
                                          - RHO(NK2)*TERM1*FNT(NK2)/(D8_0*DZ) 
       ENDDO
      ENDDO

!======================
!     TWISTING TERMS
!======================
      FACTOR1 = RHOZ(NK2)*FNT(NK2)/FNZ(NK2)
      FACTOR2 = RHOZ(NK1)*FNT(NK2)/FNZ(NK1)

      DO J = 1, mj1
       DO I = 1, mi1
        TERM2 = FACTOR1*(channel%seg(num_seg)%Z3DX(I,J,NK2)      &
                        +channel%seg(num_seg)%Z3DX(I+1,J,NK2))   &
                       *(channel%seg(num_seg)%W3D(I+1,J+1,NK2)   &
                        -channel%seg(num_seg)%W3D(I,J+1,NK2)     &
                        +channel%seg(num_seg)%W3D(I+1,J,NK2)     &
                        -channel%seg(num_seg)%W3D(I,J,NK2))      &
              + FACTOR2*(channel%seg(num_seg)%Z3DX(I,J,NK1)      &
                        +channel%seg(num_seg)%Z3DX(I+1,J,NK1))   &
                       *(channel%seg(num_seg)%W3D(I+1,J+1,NK1)   &
                        -channel%seg(num_seg)%W3D(I,J+1,NK1)     &
                        +channel%seg(num_seg)%W3D(I+1,J,NK1)     &
                        -channel%seg(num_seg)%W3D(I,J,NK1))
       
        channel%seg(num_seg)%FZTOP(I,J,L) = channel%seg(num_seg)%FZTOP(I,J,L) &
                                          + TERM2/(D8_0*DX)
       ENDDO
      ENDDO

      DO J = 1, mj1
       DO I = 1, mi1
         TERM3 = FACTOR1*(channel%seg(num_seg)%Z3DY(I,J,NK2)      &
                         +channel%seg(num_seg)%Z3DY(I,J+1,NK2))   &
                        *(channel%seg(num_seg)%W3D(I+1,J+1,NK2)   &
                         -channel%seg(num_seg)%W3D(I+1,J,NK2)     &
                         +channel%seg(num_seg)%W3D(I,J+1,NK2)     &
                         -channel%seg(num_seg)%W3D(I,J,NK2))      &
               + FACTOR2*(channel%seg(num_seg)%Z3DY(I,J,NK1)      &
                         +channel%seg(num_seg)%Z3DY(I,J+1,NK1))   &
                        *(channel%seg(num_seg)%W3D(I+1,J+1,NK1)   &
                         -channel%seg(num_seg)%W3D(I+1,J,NK1)     &
                         +channel%seg(num_seg)%W3D(I,J+1,NK1)     &
                         -channel%seg(num_seg)%W3D(I,J,NK1))
      
        channel%seg(num_seg)%FZTOP(I,J,L) = channel%seg(num_seg)%FZTOP(I,J,L) &
                                          + TERM3/(D8_0*DY)
       ENDDO
      ENDDO

      DEALLOCATE(FLX)
      DEALLOCATE(UPi)
      DEALLOCATE(UMi)
      DEALLOCATE(UPSRi)
      DEALLOCATE(UMSRi)
      DEALLOCATE(TEMP)

!*******************************************************************   
      ENDDO  ! num_seg 
!*******************************************************************        
      
   END SUBROUTINE zeta_3d

!=======================================================================
   SUBROUTINE ZETA_DIAG (channel)
!=======================================================================
!  Diagnose zeta (z3dz)  k < nk2
!  Data needed only for (mim_c:mip_c, mjm_c:mjp_c, k)

      type(channel_t), intent(inout) :: channel   ! channel data
      
      INTEGER :: I,J,K,num_seg,mim_c,mip_c,mjm_c,mjp_c

!===============================   
      DO num_seg = 1, 4
!=============================== 
      mim_c  = channel%seg(num_seg)%mim_c    
      mip_c  = channel%seg(num_seg)%mip_c 
        
      mjm_c  = channel%seg(num_seg)%mjm_c        
      mjp_c  = channel%seg(num_seg)%mjp_c  
      
      DO K = nk1, 1, -1
      DO J = mjm_c, mjp_c
      DO I = mim_c, mip_c
       channel%seg(num_seg)%Z3DZ(I,J,K) = channel%seg(num_seg)%Z3DZ(I,J,K+1)    &
          +(channel%seg(num_seg)%RG_V(I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,K) &
           -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,K))    &
          *DZ/(channel%seg(num_seg)%RG_Z(I,J)*DX*FNZ(K))                        &
          +(channel%seg(num_seg)%RG_U(I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,K) &
           -channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,K))    &
          *DZ/(channel%seg(num_seg)%RG_Z(I,J)*DY*FNZ(K))
      ENDDO
      ENDDO
      ENDDO

      DO J = mjm_c, mjp_c
      DO I = mim_c, mip_c
       channel%seg(num_seg)%Z3DZ(I,J,NK3) = channel%seg(num_seg)%Z3DZ(I,J,NK2)    &
          -(channel%seg(num_seg)%RG_V(I+1,J)*channel%seg(num_seg)%Z3DX(I+1,J,NK2) &
           -channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%Z3DX(I,J,NK2))    &
          *DZ/(channel%seg(num_seg)%RG_Z(I,J)*DX*FNZ(NK2))                        &
          -(channel%seg(num_seg)%RG_U(I,J+1)*channel%seg(num_seg)%Z3DY(I,J+1,NK2) &
           -channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%Z3DY(I,J,NK2))    &
          *DZ/(channel%seg(num_seg)%RG_Z(I,J)*DY*FNZ(NK2))
      ENDDO
      ENDDO
!===============================   
      ENDDO  ! num_seg 
!===============================

   END SUBROUTINE zeta_diag

!=======================================================================      
   SUBROUTINE VADVEC_H1 (mi1,mim,mim_a,mim_c,mip,mip_c,    &
                         mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,    &
                         DM_P2,DM_P1,DM_M1,DM_M2,QVAL,     &
                         WINDV,DGRID,direction,MODFLUX,FINAL)
!=======================================================================                         
      INTEGER, INTENT(IN) :: mi1,mim,mim_a,mim_c,mip,mip_c
      INTEGER, INTENT(IN) :: mj1,mjm,mjm_a,mjm_c,mjp,mjp_c
      INTEGER, INTENT(IN), DIMENSION(mim_c:mi1,mjm_c:mj1) :: DM_P2,DM_P1
      INTEGER, INTENT(IN), DIMENSION(mim_c:mi1,mjm_c:mj1) :: DM_M2,DM_M1
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,mjm:mjp) :: QVAL   
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim_a:mip_c,mjm_a:mjp_c) :: WINDV
      REAL (KIND=dbl_kind), INTENT(IN) :: DGRID
      CHARACTER (LEN=*), INTENT(IN) :: direction
      LOGICAL, INTENT(IN) :: MODFLUX
      
      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(mi1,mj1) :: FINAL   
      
      ! Local variables
      REAL (KIND=dbl_kind) :: TEMPX(mim_a:mip_c,mj1),TEMPY(mi1,mjm_a:mjp_c)
      REAL (KIND=dbl_kind), DIMENSION(mim_a:mip_c,mjm_a:mjp_c) :: UP,UM,UPSR,UMSR
      REAL (KIND=dbl_kind), DIMENSION(mim_c:mi1,mjm_c:mj1) :: FLX
      INTEGER ::  I, J

!*******************************************
    IF ( direction == 'X') then
!*******************************************
      DO J = 1, mj1
      DO I = mim_a, mip_c
       TEMPX(I,j) = WINDV(I,j)
      ENDDO
      ENDDO

      DO J = 1, mj1
      DO I = mim_a, mip_c
        UP(I,j) = D0_5*(TEMPX(I,j)+ABS(TEMPX(I,j)))
        UM(I,j) = D0_5*(TEMPX(I,j)-ABS(TEMPX(I,j)))
      ENDDO
      ENDDO

      DO j = 1, mj1
      DO I = mim_a, mip_c
        UPSR(I,j) = SQRT(UP(I,j))
        UMSR(I,j) = SQRT(ABS(UM(I,j)))
      ENDDO
      ENDDO

!----------------------------------
      IF (MODFLUX) THEN
!----------------------------------
      DO J = 1, mJ1
      DO I = mim_c, mi1

      IF (DM_P1(I,J).EQ.1.OR.DM_M1(I,J).EQ.1) THEN
       FLX(I,j) = D0_0
      ELSE

       IF(WINDV(I,j).LE.D0_0) THEN
         IF (DM_P2(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                &
                           +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
         ENDIF
       ELSE
         IF (DM_M2(I,J).EQ.1) THEN
          FLX(I,j)=D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                &
                           +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
         ENDIF
       ENDIF

      ENDIF

      ENDDO
      ENDDO
!----------------------------------
      ELSE
!----------------------------------
      DO J = 1, mJ1
      DO I = mim_c, mi1
        FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))               &
                 - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                &
                         -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
                         +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                &
                         +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
      ENDDO
      ENDDO
!----------------------------------
      ENDIF
!----------------------------------

      DO j = 1, mj1
      DO I = 1, mi1
       FINAL(I,j) = -(FLX(I,j)-FLX(I-1,j))/DGRID
      ENDDO
      ENDDO

!*******************************************
    ELSE
!   direction == 'Y'
!*******************************************
      DO J = mjm_a, mjp_c
      DO I = 1, mi1
       TEMPY(I,j) = WINDV(I,j)
      ENDDO
      ENDDO

      DO J = mjm_a, mjp_c
      DO I = 1, mi1
        UP(I,j) = D0_5*(TEMPY(I,j)+ABS(TEMPY(I,j)))
        UM(I,j) = D0_5*(TEMPY(I,j)-ABS(TEMPY(I,j)))
      ENDDO
      ENDDO

      DO j = mjm_a, mjp_c
      DO I = 1, mi1
        UPSR(I,j) = SQRT(UP(I,j))
        UMSR(I,j) = SQRT(ABS(UM(I,j)))
      ENDDO
      ENDDO
!----------------------------------
      IF (MODFLUX) THEN
!----------------------------------
      DO J = mjm_c, mJ1
      DO I = 1, mi1

      IF (DM_P1(I,J).EQ.1.OR.DM_M1(I,J).EQ.1) THEN
       FLX(I,j) = D0_0
      ELSE

       IF(WINDV(i,J).LE.D0_0) THEN
         IF (DM_P2(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                &
                           +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
         ENDIF
       ELSE
         IF (DM_M2(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                &
                           +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
         ENDIF
       ENDIF

      ENDIF

      ENDDO
      ENDDO
!----------------------------------
      ELSE
!----------------------------------
      DO J = mjm_c, mJ1
      DO I = 1, mi1
        FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))               &
                 - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                &
                         -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
                         +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                &
                         +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
      ENDDO
      ENDDO
!----------------------------------
      ENDIF
!----------------------------------

      DO j = 1, mj1
      DO I = 1, mi1
       FINAL(I,j) = -(FLX(I,j)-FLX(I,j-1))/DGRID
      ENDDO
      ENDDO

!*******************************************
    ENDIF
!*******************************************

   END SUBROUTINE vadvec_h1

!=======================================================================
   SUBROUTINE VADVEC_H2 (mi1,mim,mim_a,mim_c,mip,mip_c,    &
                         mj1,mjm,mjm_a,mjm_c,mjp,mjp_c,    &
                         DM_P1,DM_M1,QVAL, &
                         WINDV,DGRID,direction,MODFLUX,FINAL)
!=======================================================================    
      INTEGER, INTENT(IN) :: mi1,mim,mim_a,mim_c,mip,mip_c
      INTEGER, INTENT(IN) :: mj1,mjm,mjm_a,mjm_c,mjp,mjp_c                     
      INTEGER, INTENT(IN), DIMENSION(mim_c:mi1,mjm_c:mj1) :: DM_P1
      INTEGER, INTENT(IN), DIMENSION(mim_c:mi1,mjm_c:mj1) :: DM_M1
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,mjm:mjp) :: QVAL 
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim_a:mip_c,mjm_a:mjp_c) :: WINDV 
      REAL (KIND=dbl_kind), INTENT(IN) :: DGRID
      CHARACTER (LEN=*), INTENT(IN) :: direction
      LOGICAL, INTENT(IN) :: MODFLUX
      
      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(mi1,mj1) :: FINAL     

      ! Local variables
      REAL (KIND=dbl_kind) :: TEMPX(mim_a:mip_c,mj1),TEMPY(mi1,mjm_a:mjp_c)
      REAL (KIND=dbl_kind), DIMENSION(mim_a:mip_c,mjm_a:mjp_c) :: UP,UM,UPSR,UMSR
      REAL (KIND=dbl_kind), DIMENSION(mim_c:mi1,mjm_c:mj1) :: FLX
      INTEGER ::  I, J

!*******************************************
    IF ( direction == 'X') then
!*******************************************
      DO J = 1, mj1
      DO I = mim_a, mip_c
       TEMPX(I,j) = WINDV(I,j)
      ENDDO
      ENDDO

      DO J = 1, mj1
      DO I = mim_a, mip_c
        UP(I,j) = D0_5*(TEMPX(I,j)+ABS(TEMPX(I,j)))
        UM(I,j) = D0_5*(TEMPX(I,j)-ABS(TEMPX(I,j)))
      ENDDO
      ENDDO

      DO j = 1, mj1
      DO I = mim_a, mip_c
        UPSR(I,j) = SQRT(UP(I,j))
        UMSR(I,j) = SQRT(ABS(UM(I,j)))
      ENDDO
      ENDDO

!----------------------------------
      IF (MODFLUX) THEN
!----------------------------------
      DO J = 1, mJ1
      DO I = mim_c, mi1

      IF(WINDV(i,J).LE.D0_0) THEN
         IF (DM_P1(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                &
                           +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
         ENDIF
      ELSE
         IF (DM_M1(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                &
                           +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
         ENDIF
      ENDIF

      ENDDO
      ENDDO
!----------------------------------
      ELSE
!----------------------------------
      DO J = 1, mJ1
      DO I = mim_c, mi1
        FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))                &
                 - ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))                 &
                         -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))   &
                         +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))                 &
                         +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/D6_0
      ENDDO
      ENDDO
!----------------------------------
      ENDIF
!----------------------------------

      DO j = 1, mj1
      DO I = 1, mi1
       FINAL(I,j) = -(FLX(I,j)-FLX(I-1,j))/DGRID
      ENDDO
      ENDDO

!*******************************************
    ELSE
!   direction == 'Y'
!*******************************************
      DO J = mjm_a, mjp_c
      DO I = 1, mi1
       TEMPY(I,j) = WINDV(I,j)
      ENDDO
      ENDDO

      DO J = mjm_a, mjp_c
      DO I = 1, mi1
        UP(I,j) = D0_5*(TEMPY(I,j)+ABS(TEMPY(I,j)))
        UM(I,j) = D0_5*(TEMPY(I,j)-ABS(TEMPY(I,j)))
      ENDDO
      ENDDO

      DO j = mjm_a, mjp_c
      DO I = 1, mi1
        UPSR(I,j) = SQRT(UP(I,j))
        UMSR(I,j) = SQRT(ABS(UM(I,j)))
      ENDDO
      ENDDO
!----------------------------------
      IF (MODFLUX) THEN
!----------------------------------
      DO J = mjm_c, mJ1
      DO I = 1, mi1

       IF(WINDV(i,J).LE.D0_0) THEN
         IF (DM_P1(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                &
                           +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
         ENDIF
       ELSE
         IF (DM_M1(I,J).EQ.1) THEN
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))
         ELSE
          FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))               &
                   - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                &
                           -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
                           +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                &
                           +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
         ENDIF
       ENDIF

      ENDDO
      ENDDO
!----------------------------------
      ELSE
!----------------------------------
      DO J = mjm_c, mJ1
      DO I = 1, mi1 
        FLX(I,j) = D0_5*WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))                &
                 - ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))                 &
                         -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))   &
                         +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))                 &
                         +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/D6_0
      ENDDO
      ENDDO
!----------------------------------
      ENDIF
!----------------------------------

      DO j = 1, mj1
      DO I = 1, mi1
       FINAL(I,j) = -(FLX(I,j)-FLX(I,j-1))/DGRID
      ENDDO
      ENDDO

!*******************************************
    ENDIF
!*******************************************

   END SUBROUTINE vadvec_h2

!=======================================================================
   SUBROUTINE VADVEC_V(WINDV,QVAL,FNZ,KZ,KZP1,KZP2,KMIN,DZ,FINAL)
!=======================================================================
      INTEGER, INTENT(IN) :: KZ,KZP1,KZP2,KMIN
      REAL (KIND=DBL_KIND), INTENT(IN), DIMENSION(KZP1) :: WINDV,QVAL
      REAL (KIND=DBL_KIND), INTENT(IN), DIMENSION(KZP2) :: FNZ
      REAL (KIND=DBL_KIND), INTENT(IN) :: DZ
      REAL (KIND=DBL_KIND), INTENT(OUT), DIMENSION(KZ) :: FINAL

      ! Local variables
      REAL (KIND=dbl_kind), DIMENSION(kz) :: UP,UM,UPSR,UMSR,FLX
      INTEGER :: K

      DO K = KMIN, KZ
       UP(K) = D0_5*(WINDV(K)+ABS(WINDV(K)))
       UM(K) = D0_5*(WINDV(K)-ABS(WINDV(K)))
      ENDDO
      DO K = KMIN, KZ
       UPSR(K) = SQRT(UP(K))
       UMSR(K) = SQRT(ABS(UM(K)))
      ENDDO

      DO K = KMIN+1, KZ-1
       FLX(K) = D0_5*WINDV(K)*(QVAL(K+1)+QVAL(K))              &
              - ALADV*(UP(K)*(QVAL(K+1)-QVAL(K))               &
                      -UPSR(K)*UPSR(K-1)*(QVAL(K)-QVAL(K-1))   &
                      +UM(K)*(QVAL(K)-QVAL(K+1))               &
                      +UMSR(K)*UMSR(K+1)*(QVAL(K+1)-QVAL(K+2)))/D6_0
      ENDDO
      IF(WINDV(KZ).GE.D0_0) THEN
        FLX(KZ) = D0_5*WINDV(KZ)*(QVAL(KZ+1)+QVAL(KZ))         &
                - ALADV*(UP(KZ)*(QVAL(KZ+1)-QVAL(KZ))          &
                        -UPSR(KZ)*UPSR(KZ-1)*(QVAL(KZ)-QVAL(KZ-1)))/D6_0
      ELSE
        FLX(KZ) = D0_5*WINDV(KZ)*(QVAL(KZ+1)+QVAL(KZ))
      ENDIF

      IF(WINDV(KMIN).GE.D0_0) THEN
        FLX(KMIN) = D0_5*WINDV(KMIN)*(QVAL(KMIN+1)+QVAL(KMIN))
      ELSE
        FLX(KMIN) = D0_5*WINDV(KMIN)*(QVAL(KMIN+1)+QVAL(KMIN))              &
                  - ALADV*(UM(KMIN)*(QVAL(KMIN)-QVAL(KMIN+1))               &
                          +UMSR(KMIN)*UMSR(KMIN+1)*(QVAL(KMIN+1)-QVAL(KMIN+2)))/D6_0
      ENDIF

      DO K = KMIN+1, KZ
       FINAL(K) = -(FLX(K)-FLX(K-1))*FNZ(K)/DZ
      ENDDO

   END SUBROUTINE vadvec_v

!=======================================================================
   SUBROUTINE ABM_3D ( N1, N2, channel )
!=======================================================================   
!  Updating the vorticity components 

      INTEGER, INTENT(IN) :: N1,N2                ! previous and current time index (AB scheme)       
      type(channel_t), intent(inout) :: channel   ! channel data 

      ! Local
      INTEGER :: I,J,K,klowu,klowv,num_seg,mi1,mj1

!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
! JH  Applying free-slip condition at the surface.

      DO J = 1, mj1
      DO I = 1, mi1
      klowu = channel%seg(num_seg)%KLOWU_IJ(I,J)
      klowv = channel%seg(num_seg)%KLOWV_IJ(I,J)

! Predict vorticity
       DO K = klowv, nk1
        channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K)     &
                                         + DT*channel%seg(num_seg)%FZXBU(I,J,K) &
                                         + A*channel%seg(num_seg)%FZX(I,J,K,N2) &
                                         + B*channel%seg(num_seg)%FZX(I,J,K,N1)          
       ENDDO
       DO K = 2, klowv-1
        channel%seg(num_seg)%Z3DX(I,J,K) = D0_0
       ENDDO

       DO K = klowu,nk1
        channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K)     &
                                         + DT*channel%seg(num_seg)%FZYBU(I,J,K) &
                                         + A*channel%seg(num_seg)%FZY(I,J,K,N2) &
                                         +B*channel%seg(num_seg)%FZY(I,J,K,N1)            
       ENDDO
       DO K = 2, klowu-1
        channel%seg(num_seg)%Z3DY(I,J,K) = D0_0
       ENDDO

      ENDDO
      ENDDO
!===============================   
      ENDDO  ! num_seg 
!===============================

   END SUBROUTINE abm_3d

!=======================================================================
   SUBROUTINE ABM_3D_TURB ( channel )
!======================================================================= 
!  Updating the vorticity components
    
      type(channel_t), intent(inout) :: channel   ! channel data 

      ! Local
      INTEGER :: I,J,K,klowu,klowv,num_seg,mi1,mj1


!===============================   
      DO num_seg = 1, 4
!=============================== 
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
      
      DO J = 1, mj1
      DO I = 1, mi1
      klowu = channel%seg(num_seg)%KLOWU_IJ(I,J)
      klowv = channel%seg(num_seg)%KLOWV_IJ(I,J)

       DO K = klowv, nk1
        channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K)  &
                                         + DT*channel%seg(num_seg)%FZXTB(I,J,K)
       ENDDO

       DO K = klowu, nk1
        channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K)  &
                                         + DT*channel%seg(num_seg)%FZYTB(I,J,K)
       ENDDO

      ENDDO
      ENDDO
!===============================   
      ENDDO  ! num_seg 
!===============================

  END SUBROUTINE abm_3d_turb

!=======================================================================
      SUBROUTINE RELAXATION (channel)
!=======================================================================   
!     Relax the vGCM-cell mean of deviation toward zero
!     (prognostic) channel width is assumed to be 1 (chn=1)
!
!     PLAN: In the averaging, TOPOGRAPHY effect should be included.
!           (Any mapping problem?)
!------------------------------------------------------------------------ 
      type(channel_t), INTENT(INOUT) :: channel   ! channel data
       
!     Local variables
      LOGICAL :: INTERPOL = .FALSE.
      
      REAL (KIND=dbl_kind), DIMENSION(4,2,4) :: AM
      REAL (KIND=dbl_kind), DIMENSION(4,2,4) :: AMI
      
      REAL (KIND=dbl_kind),DIMENSION(nk1-1,0:nVGCM_seg+1,4) :: MN_ZX
      REAL (KIND=dbl_kind),DIMENSION(nk1-1,0:nVGCM_seg+1,4) :: MN_ZY
      REAL (KIND=dbl_kind),DIMENSION(0:nVGCM_seg+1,4) :: MN_ZZ
      
      REAL (KIND=dbl_kind) :: XMN(netsz)
      
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_ZX
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: TEMP_ZY
      REAL (KIND=dbl_kind), DIMENSION(:,:),   ALLOCATABLE :: TEMP_ZZ
      
      INTEGER :: mi1,mj1,k
      INTEGER :: num_seg,ng,lsta,lend,lcen_m,lcen_p,chl,chn,npo

      chn = 1    !  Assume (prognostic) channel width is 1.
      
! Calculate the vGCM-cell average of deviation.     
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment  
      
      IF (INTERPOL) THEN
      AMI(:,:,num_seg) = channel%seg(num_seg)%AMI_vGCM_halo(:,:)
      ! Inverse of the transformation matrix specially prepared 
      ! at the 1st halos of the vGCM in each channel segment.
      ! 1st dimension: matrix index (1) 11, (2) 12, (3) 21, (4) 22
      ! 2nd dimension: (1) minus direction (2) plus direction of halo points
      ! 3rd dimension: segment index
      ENDIF
      
!     Calculate the segment average of deviation
!---------------------------- 
      IF (mi1.GT.mj1) THEN
!---------------------------- x-channel seg.
        IF (INTERPOL) THEN
        lcen_m = INT(channel%seg(num_seg)%lcen(1))
        lcen_p = INT(channel%seg(num_seg)%lcen(nVGCM_seg))
        
        DO npo = 1, 4
          AM(npo,1,num_seg) = channel%seg(num_seg)%AM_T(npo,lcen_m,chn)
          AM(npo,2,num_seg) = channel%seg(num_seg)%AM_T(npo,lcen_p,chn)
        ENDDO 
        ENDIF
        
      DO ng = 1, nVGCM_seg
        
        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng)
        
        DO K = 2, nk1
        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DX(chl,chn,K) &
                  - channel%seg(num_seg)%Z3DX_BG(chl,chn,K) 
        ENDDO
        MN_ZX(K-1,ng,num_seg) = SUM(XMN)/FLOAT(netsz)

        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DY(chl,chn,K) &
                  - channel%seg(num_seg)%Z3DY_BG(chl,chn,K) 
        ENDDO
        MN_ZY(K-1,ng,num_seg) = SUM(XMN)/FLOAT(netsz)
        ENDDO

        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DZ(chl,chn,nk2) &
                  - channel%seg(num_seg)%Z3DZ_BG(chl,chn,1) 
        ENDDO
        MN_ZZ(ng,num_seg) = SUM(XMN)/FLOAT(netsz)

      ENDDO  ! ng-loop
!---------------------------- 
      ELSE
!---------------------------- y-channel seg.
        IF (INTERPOL) THEN 
        lcen_m = INT(channel%seg(num_seg)%lcen(1))
        lcen_p = INT(channel%seg(num_seg)%lcen(nVGCM_seg))
        
        DO npo = 1, 4
          AM(npo,1,num_seg) = channel%seg(num_seg)%AM_T(npo,chn,lcen_m)
          AM(npo,2,num_seg) = channel%seg(num_seg)%AM_T(npo,chn,lcen_p)
        ENDDO 
        ENDIF
        
      DO ng = 1, nVGCM_seg
        
        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng) 
        
        DO K = 2, nk1
        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DX(chn,chl,K) &
                  - channel%seg(num_seg)%Z3DX_BG(chn,chl,K) 
        ENDDO
        MN_ZX(K-1,ng,num_seg) = SUM(XMN)/FLOAT(netsz)

        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DY(chn,chl,K) &
                  - channel%seg(num_seg)%Z3DY_BG(chn,chl,K) 
        ENDDO
        MN_ZY(K-1,ng,num_seg) = SUM(XMN)/FLOAT(netsz)
        ENDDO

        DO chl = lsta, lend
         XMN(chl) = channel%seg(num_seg)%Z3DZ(chn,chl,nk2) &
                  - channel%seg(num_seg)%Z3DZ_BG(chn,chl,1) 
        ENDDO
        MN_ZZ(ng,num_seg) = SUM(XMN)/FLOAT(netsz)

      ENDDO  ! ng-loop
!----------------------------
      ENDIF
!----------------------------   
      
!****************************   
      ENDDO  ! num_seg 
!****************************

      IF (INTERPOL) THEN
      ! Filling the halos along the channel: vector mapping 
      
      SELECT CASE (channel%num_chg)
      CASE(1)
        CALL BOUND_G1_V (AM,AMI,MN_ZX,MN_ZY)
        CALL BOUND_G1_S (MN_ZZ)
        
      CASE(2)
        CALL BOUND_G2_V (AM,AMI,MN_ZX,MN_ZY)
        CALL BOUND_G2_S (MN_ZZ)
        
      CASE(3)
        CALL BOUND_G3_V (AM,AMI,MN_ZX,MN_ZY)
        CALL BOUND_G3_S (MN_ZZ)
      END SELECT
      
      ENDIF
      
!**********************************************************************   
      DO num_seg = 1, 4
!**********************************************************************
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment    
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

      ALLOCATE(TEMP_ZX(mi1,mj1,nk1-1))
      ALLOCATE(TEMP_ZY(mi1,mj1,nk1-1))
      ALLOCATE(TEMP_ZZ(mi1,mj1))
      
! Interpolation of vGCM values to CRM grid points      
!----------------------------
      IF (mi1.GT.mj1) THEN
!---------------------------- x-channel seg.
      DO ng = 1, nVGCM_seg
        
        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng)
        lcen = channel%seg(num_seg)%lcen(ng)
        
        DO chl = lsta, lend 

        IF (INTERPOL) THEN
         
         DIST1 = lcen - float(chl)
         IF (DIST1.GE.D0_0) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF 
          
         DO K = 1, nk1-1       
          TEMP_ZX(chl,chn,K) = &
            (DIST*MN_ZX(K,NPO,num_seg)+DIST2*MN_ZX(K,ng,num_seg))/FLOAT(netsz)
          TEMP_ZY(chl,chn,K) = &
            (DIST*MN_ZY(K,NPO,num_seg)+DIST2*MN_ZY(K,ng,num_seg))/FLOAT(netsz)
         ENDDO 
                
          TEMP_ZZ(chl,chn) = &
            (DIST*MN_ZZ(NPO,num_seg)+DIST2*MN_ZZ(ng,num_seg))/FLOAT(netsz)
            
        ELSE  ! INTERPOL

         DO K = 1, nk1-1       
          TEMP_ZX(chl,chn,K) = MN_ZX(K,ng,num_seg)
          TEMP_ZY(chl,chn,K) = MN_ZY(K,ng,num_seg)
         ENDDO 
          TEMP_ZZ(chl,chn)   = MN_ZZ(ng,num_seg)
                    
        ENDIF  ! INTERPOL    
         
        ENDDO  ! chl-loop 

      ENDDO  ! ng-loop
!----------------------------
      ELSE
!---------------------------- y-channel seg.
      DO ng = 1, nVGCM_seg
        
        lsta = channel%seg(num_seg)%lsta(ng)
        lend = channel%seg(num_seg)%lend(ng)
        lcen = channel%seg(num_seg)%lcen(ng)
        
        DO chl = lsta, lend 

        IF (INTERPOL) THEN
        
         DIST1 = lcen - float(chl)
         IF (DIST1.GE.D0_0) THEN
           DIST2 = FLOAT(netsz) - DIST1
           DIST  = DIST1
           NPO   = ng - 1
         ELSE
           DIST2 = FLOAT(netsz) + DIST1
           DIST  = ABS(DIST1)
           NPO   = ng + 1
         ENDIF 
          
         DO K = 1, nk1       
          TEMP_ZX(chn,chl,K) = &
            (DIST*MN_ZX(K,NPO,num_seg)+DIST2*MN_ZX(K,ng,num_seg))/FLOAT(netsz)
          TEMP_ZY(chn,chl,K) = &
            (DIST*MN_ZY(K,NPO,num_seg)+DIST2*MN_ZY(K,ng,num_seg))/FLOAT(netsz)
         ENDDO        
          TEMP_ZZ(chn,chl) = &
            (DIST*MN_ZZ(NPO,num_seg)+DIST2*MN_ZZ(ng,num_seg))/FLOAT(netsz)

        ELSE  ! INTERPOL

         DO K = 1, nk1       
          TEMP_ZX(chn,chl,K) = MN_ZX(K,ng,num_seg)
          TEMP_ZY(chn,chl,K) = MN_ZY(K,ng,num_seg)
         ENDDO        
          TEMP_ZZ(chn,chl) = MN_ZZ(ng,num_seg)
                    
        ENDIF ! INTERPOL
        
        ENDDO  ! chl-loop 

      ENDDO  ! ng-loop
!----------------------------
      ENDIF
!---------------------------- 

!     Apply the relaxation (topography is not considered)
      DO K = 2, NK1
       DO J = 1, MJ1
        DO I = 1, MI1

          IF (channel%seg(num_seg)%TAU_RX_ZX(I,J,K).NE.D0_0) THEN
           channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K) &
                    -DT*TEMP_ZX(I,J,K-1)/channel%seg(num_seg)%TAU_RX_ZX(I,J,K)
          ENDIF

          IF (channel%seg(num_seg)%TAU_RX_ZY(I,J,K).NE.D0_0) THEN
           channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K) &
                    -DT*TEMP_ZY(I,J,K-1)/channel%seg(num_seg)%TAU_RX_ZY(I,J,K)
          ENDIF

        ENDDO
       ENDDO
      ENDDO

       DO J = 1, MJ1
        DO I = 1, MI1
          IF (channel%seg(num_seg)%TAU_RX_ZZ(I,J,1).NE.D0_0) THEN
           channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                    -DT*TEMP_ZZ(I,J)/channel%seg(num_seg)%TAU_RX_ZZ(I,J,1)
          ENDIF
        ENDDO
       ENDDO

      DEALLOCATE(TEMP_ZX)
      DEALLOCATE(TEMP_ZY)
      DEALLOCATE(TEMP_ZZ)

!**********************************************************************   
      ENDDO  ! num_seg 
!**********************************************************************

      CONTAINS

      SUBROUTINE BOUND_G1_V (AM,AMI,U,V)
!     Vector: Filling the halos along the channel (Group1) 

      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AM    ! transformation matrix at the center of vGCM cell
                  ! 1st dimension (matrix index); 2nd (1st & last vGCM cells); 3rd (seg. index)
      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AMI   ! inverse of transformation matrix at the center of vGCM cell
                  ! 1st (matrix index); 2nd (-/+ halo vGCM cells); 3rd (seg. index)
                   
      REAL (KIND=dbl_kind),DIMENSION(nk1-1,0:nVGCM_seg+1,4),INTENT(INOUT) :: U,V
      REAL (KIND=dbl_kind) :: U_rll,V_rll
      INTEGER :: chl_p,chl_m,K

        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        DO K = 1, nk1-1
          U_rll = (AM(1,1,2)*U(K,1,2) + AM(2,1,2)*V(K,1,2))
          V_rll = (AM(3,1,2)*U(K,1,2) + AM(4,1,2)*V(K,1,2))
          U(K,chl_p,1) = (AMI(1,2,1)*U_rll + AMI(2,2,1)*V_rll)
          V(K,chl_p,1) = (AMI(3,2,1)*U_rll + AMI(4,2,1)*V_rll)

          U_rll = (AM(1,2,4)*U(K,nVGCM_seg,4) + AM(2,2,4)*V(K,nVGCM_seg,4))
          V_rll = (AM(3,2,4)*U(K,nVGCM_seg,4) + AM(4,2,4)*V(K,nVGCM_seg,4))
          U(K,chl_m,1) = (AMI(1,1,1)*U_rll + AMI(2,1,1)*V_rll)
          V(K,chl_m,1) = (AMI(3,1,1)*U_rll + AMI(4,1,1)*V_rll)          
          
          U_rll = (AM(1,1,3)*U(K,1,3) + AM(2,1,3)*V(K,1,3))
          V_rll = (AM(3,1,3)*U(K,1,3) + AM(4,1,3)*V(K,1,3))
          U(K,chl_p,2) = (AMI(1,2,2)*U_rll + AMI(2,2,2)*V_rll)
          V(K,chl_p,2) = (AMI(3,2,2)*U_rll + AMI(4,2,2)*V_rll)           

          U_rll = (AM(1,2,1)*U(K,nVGCM_seg,1) + AM(2,2,1)*V(K,nVGCM_seg,1))
          V_rll = (AM(3,2,1)*U(K,nVGCM_seg,1) + AM(4,2,1)*V(K,nVGCM_seg,1))
          U(K,chl_m,2) = (AMI(1,1,2)*U_rll + AMI(2,1,2)*V_rll)
          V(K,chl_m,2) = (AMI(3,1,2)*U_rll + AMI(4,1,2)*V_rll)               
        
          U_rll = (AM(1,1,4)*U(K,1,4) + AM(2,1,4)*V(K,1,4))
          V_rll = (AM(3,1,4)*U(K,1,4) + AM(4,1,4)*V(K,1,4))
          U(K,chl_p,3) = (AMI(1,2,3)*U_rll + AMI(2,2,3)*V_rll)
          V(K,chl_p,3) = (AMI(3,2,3)*U_rll + AMI(4,2,3)*V_rll)      
        
          U_rll = (AM(1,2,2)*U(K,nVGCM_seg,2) + AM(2,2,2)*V(K,nVGCM_seg,2))
          V_rll = (AM(3,2,2)*U(K,nVGCM_seg,2) + AM(4,2,2)*V(K,nVGCM_seg,2))
          U(K,chl_m,3) = (AMI(1,1,3)*U_rll + AMI(2,1,3)*V_rll)
          V(K,chl_m,3) = (AMI(3,1,3)*U_rll + AMI(4,1,3)*V_rll)      
          
          U_rll = (AM(1,1,1)*U(K,1,1) + AM(2,1,1)*V(K,1,1))
          V_rll = (AM(3,1,1)*U(K,1,1) + AM(4,1,1)*V(K,1,1))
          U(K,chl_p,4) = (AMI(1,2,4)*U_rll + AMI(2,2,4)*V_rll)
          V(K,chl_p,4) = (AMI(3,2,4)*U_rll + AMI(4,2,4)*V_rll)      

          U_rll = (AM(1,2,3)*U(K,nVGCM_seg,3) + AM(2,2,3)*V(K,nVGCM_seg,3))
          V_rll = (AM(3,2,3)*U(K,nVGCM_seg,3) + AM(4,2,3)*V(K,nVGCM_seg,3))
          U(K,chl_m,4) = (AMI(1,1,4)*U_rll + AMI(2,1,4)*V_rll)
          V(K,chl_m,4) = (AMI(3,1,4)*U_rll + AMI(4,1,4)*V_rll)      
        ENDDO  
            
      END SUBROUTINE bound_g1_v
      
      SUBROUTINE BOUND_G2_V (AM,AMI,U,V)
!     Vector: Filling the halos along the channel (Group2) 

      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AM    ! transformation matrix at the center of vGCM cell
                  ! 1st dimension (matrix index); 2nd (1st & last vGCM cells); 3rd (seg. index)
      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AMI   ! inverse of transformation matrix at the center of vGCM cell
                  ! 1st (matrix index); 2nd (-/+ halo vGCM cells); 3rd (seg. index)
                   
      REAL (KIND=dbl_kind),DIMENSION(nk1-1,0:nVGCM_seg+1,4),INTENT(INOUT) :: U,V
      REAL (KIND=dbl_kind) :: U_rll,V_rll
      INTEGER :: chl_p,chl_m,K
      
        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        DO K = 1, nk1-1
          U_rll = (AM(1,1,2)*U(K,1,2) + AM(2,1,2)*V(K,1,2))
          V_rll = (AM(3,1,2)*U(K,1,2) + AM(4,1,2)*V(K,1,2))
          U(K,chl_p,1) = (AMI(1,2,1)*U_rll + AMI(2,2,1)*V_rll)
          V(K,chl_p,1) = (AMI(3,2,1)*U_rll + AMI(4,2,1)*V_rll)

          U_rll = (AM(1,1,4)*U(K,1,4) + AM(2,1,4)*V(K,1,4))
          V_rll = (AM(3,1,4)*U(K,1,4) + AM(4,1,4)*V(K,1,4))
          U(K,chl_m,1) = (AMI(1,1,1)*U_rll + AMI(2,1,1)*V_rll)
          V(K,chl_m,1) = (AMI(3,1,1)*U_rll + AMI(4,1,1)*V_rll)

          U_rll = (AM(1,1,3)*U(K,1,3) + AM(2,1,3)*V(K,1,3))
          V_rll = (AM(3,1,3)*U(K,1,3) + AM(4,1,3)*V(K,1,3))
          U(K,chl_p,2) = (AMI(1,2,2)*U_rll + AMI(2,2,2)*V_rll)
          V(K,chl_p,2) = (AMI(3,2,2)*U_rll + AMI(4,2,2)*V_rll)                  

          U_rll = (AM(1,2,1)*U(K,nVGCM_seg,1) + AM(2,2,1)*V(K,nVGCM_seg,1))
          V_rll = (AM(3,2,1)*U(K,nVGCM_seg,1) + AM(4,2,1)*V(K,nVGCM_seg,1))
          U(K,chl_m,2) = (AMI(1,1,2)*U_rll + AMI(2,1,2)*V_rll)
          V(K,chl_m,2) = (AMI(3,1,2)*U_rll + AMI(4,1,2)*V_rll) 
          
          U_rll = (AM(1,2,4)*U(K,nVGCM_seg,4) + AM(2,2,4)*V(K,nVGCM_seg,4))
          V_rll = (AM(3,2,4)*U(K,nVGCM_seg,4) + AM(4,2,4)*V(K,nVGCM_seg,4))
          U(K,chl_p,3) = (AMI(1,2,3)*U_rll + AMI(2,2,3)*V_rll)
          V(K,chl_p,3) = (AMI(3,2,3)*U_rll + AMI(4,2,3)*V_rll)                       
        
          U_rll = (AM(1,2,2)*U(K,nVGCM_seg,2) + AM(2,2,2)*V(K,nVGCM_seg,2))
          V_rll = (AM(3,2,2)*U(K,nVGCM_seg,2) + AM(4,2,2)*V(K,nVGCM_seg,2))
          U(K,chl_m,3) = (AMI(1,1,3)*U_rll + AMI(2,1,3)*V_rll)
          V(K,chl_m,3) = (AMI(3,1,3)*U_rll + AMI(4,1,3)*V_rll) 

          U_rll = (AM(1,2,3)*U(K,nVGCM_seg,3) + AM(2,2,3)*V(K,nVGCM_seg,3))
          V_rll = (AM(3,2,3)*U(K,nVGCM_seg,3) + AM(4,2,3)*V(K,nVGCM_seg,3))
          U(K,chl_p,4) = (AMI(1,2,4)*U_rll + AMI(2,2,4)*V_rll)
          V(K,chl_p,4) = (AMI(3,2,4)*U_rll + AMI(4,2,4)*V_rll)                       
        
          U_rll = (AM(1,1,1)*U(K,1,1) + AM(2,1,1)*V(K,1,1))
          V_rll = (AM(3,1,1)*U(K,1,1) + AM(4,1,1)*V(K,1,1))
          U(K,chl_m,4) = (AMI(1,1,4)*U_rll + AMI(2,1,4)*V_rll)
          V(K,chl_m,4) = (AMI(3,1,4)*U_rll + AMI(4,1,4)*V_rll)     
        ENDDO  
            
      END SUBROUTINE bound_g2_v

      SUBROUTINE BOUND_G3_V (AM,AMI,U,V)
!     Vector: Filling the halos along the channel (Group3)   

      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AM    ! transformation matrix at the center of vGCM cell
                  ! 1st dimension (matrix index); 2nd (1st & last vGCM cells); 3rd (seg. index)
      REAL (KIND=dbl_kind),DIMENSION(4,2,4), INTENT(IN) :: &
            AMI   ! inverse of transformation matrix at the center of vGCM cell
                  ! 1st (matrix index); 2nd (-/+ halo vGCM cells); 3rd (seg. index)
                     
      REAL (KIND=dbl_kind),DIMENSION(nk1-1,0:nVGCM_seg+1,4),INTENT(INOUT) :: U,V
      REAL (KIND=dbl_kind) :: U_rll,V_rll
      INTEGER :: chl_p,chl_m,K
      
        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        DO K = 1, nk1-1
          U_rll = (AM(1,1,2)*U(K,1,2) + AM(2,1,2)*V(K,1,2))
          V_rll = (AM(3,1,2)*U(K,1,2) + AM(4,1,2)*V(K,1,2))
          U(K,chl_p,1) = (AMI(1,2,1)*U_rll + AMI(2,2,1)*V_rll)
          V(K,chl_p,1) = (AMI(3,2,1)*U_rll + AMI(4,2,1)*V_rll)

          U_rll = (AM(1,1,4)*U(K,1,4) + AM(2,1,4)*V(K,1,4))
          V_rll = (AM(3,1,4)*U(K,1,4) + AM(4,1,4)*V(K,1,4))
          U(K,chl_m,1) = (AMI(1,1,1)*U_rll + AMI(2,1,1)*V_rll)
          V(K,chl_m,1) = (AMI(3,1,1)*U_rll + AMI(4,1,1)*V_rll)

          U_rll = (AM(1,2,3)*U(K,nVGCM_seg,3) + AM(2,2,3)*V(K,nVGCM_seg,3))
          V_rll = (AM(3,2,3)*U(K,nVGCM_seg,3) + AM(4,2,3)*V(K,nVGCM_seg,3))
          U(K,chl_p,2) = (AMI(1,2,2)*U_rll + AMI(2,2,2)*V_rll)
          V(K,chl_p,2) = (AMI(3,2,2)*U_rll + AMI(4,2,2)*V_rll)

          U_rll = (AM(1,2,1)*U(K,nVGCM_seg,1) + AM(2,2,1)*V(K,nVGCM_seg,1))
          V_rll = (AM(3,2,1)*U(K,nVGCM_seg,1) + AM(4,2,1)*V(K,nVGCM_seg,1))
          U(K,chl_m,2) = (AMI(1,1,2)*U_rll + AMI(2,1,2)*V_rll)
          V(K,chl_m,2) = (AMI(3,1,2)*U_rll + AMI(4,1,2)*V_rll)

          U_rll = (AM(1,2,2)*U(K,nVGCM_seg,2) + AM(2,2,2)*V(K,nVGCM_seg,2))
          V_rll = (AM(3,2,2)*U(K,nVGCM_seg,2) + AM(4,2,2)*V(K,nVGCM_seg,2))
          U(K,chl_p,3) = (AMI(1,2,3)*U_rll + AMI(2,2,3)*V_rll)
          V(K,chl_p,3) = (AMI(3,2,3)*U_rll + AMI(4,2,3)*V_rll)

          U_rll = (AM(1,2,4)*U(K,nVGCM_seg,4) + AM(2,2,4)*V(K,nVGCM_seg,4))
          V_rll = (AM(3,2,4)*U(K,nVGCM_seg,4) + AM(4,2,4)*V(K,nVGCM_seg,4))
          U(K,chl_m,3) = (AMI(1,1,3)*U_rll + AMI(2,1,3)*V_rll)
          V(K,chl_m,3) = (AMI(3,1,3)*U_rll + AMI(4,1,3)*V_rll)

          U_rll = (AM(1,1,3)*U(K,1,3) + AM(2,1,3)*V(K,1,3))
          V_rll = (AM(3,1,3)*U(K,1,3) + AM(4,1,3)*V(K,1,3))
          U(K,chl_p,4) = (AMI(1,2,4)*U_rll + AMI(2,2,4)*V_rll)
          V(K,chl_p,4) = (AMI(3,2,4)*U_rll + AMI(4,2,4)*V_rll)

          U_rll = (AM(1,1,1)*U(K,1,1) + AM(2,1,1)*V(K,1,1))
          V_rll = (AM(3,1,1)*U(K,1,1) + AM(4,1,1)*V(K,1,1))
          U(K,chl_m,4) = (AMI(1,1,4)*U_rll + AMI(2,1,4)*V_rll)
          V(K,chl_m,4) = (AMI(3,1,4)*U_rll + AMI(4,1,4)*V_rll)
        ENDDO  
            
      END SUBROUTINE bound_g3_v
            
      SUBROUTINE BOUND_G1_S (A)
!     Scalar: Filling the halos along the channel (Group1)      
      REAL (KIND=dbl_kind),DIMENSION(0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
      INTEGER :: chl_p,chl_m

        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        A(chl_p,1) = A(1,2)
        A(chl_m,1) = A(nVGCM_seg,4)
          
        A(chl_p,2) = A(1,3)  
        A(chl_m,2) = A(nVGCM_seg,1)
          
        A(chl_p,3) = A(1,4)  
        A(chl_m,3) = A(nVGCM_seg,2)
          
        A(chl_p,4) = A(1,1)  
        A(chl_m,4) = A(nVGCM_seg,3)  
            
      END SUBROUTINE bound_g1_s

      SUBROUTINE BOUND_G2_S (A)
!     Scalar: Filling the halos along the channel (Group2)      
      REAL (KIND=dbl_kind),DIMENSION(0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
      INTEGER :: chl_p,chl_m

        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        A(chl_p,1) = A(1,2)
        A(chl_m,1) = A(1,4)
          
        A(chl_p,2) = A(1,3)  
        A(chl_m,2) = A(nVGCM_seg,1)
          
        A(chl_p,3) = A(nVGCM_seg,4)  
        A(chl_m,3) = A(nVGCM_seg,2)
          
        A(chl_p,4) = A(nVGCM_seg,3)  
        A(chl_m,4) = A(1,1)
            
      END SUBROUTINE bound_g2_s

      SUBROUTINE BOUND_G3_S (A)
!     Scalar: Filling the halos along the channel (Group3)      
      REAL (KIND=dbl_kind),DIMENSION(0:nVGCM_seg+1,4),INTENT(INOUT) :: A
      
      INTEGER :: chl_p,chl_m
  
        chl_p = nVGCM_seg + 1
        chl_m = 0
        
        A(chl_p,1) = A(1,2)
        A(chl_m,1) = A(1,4)
          
        A(chl_p,2) = A(nVGCM_seg,3)  
        A(chl_m,2) = A(nVGCM_seg,1)
          
        A(chl_p,3) = A(nVGCM_seg,2)  
        A(chl_m,3) = A(nVGCM_seg,4)
          
        A(chl_p,4) = A(1,3)  
        A(chl_m,4) = A(1,1) 
            
      END SUBROUTINE bound_g3_s
            
      END SUBROUTINE relaxation

!=======================================================================
      SUBROUTINE HDIFFUSION_SMAGORINSKY (channel)
!=======================================================================      
!     Smagorinsky-type (nonlinear) numerical diffusion (horizontal)

      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local variables
      REAL (KIND=dbl_kind) :: FAC_NLD = D1_0
      REAL (KIND=dbl_kind) :: CRITV   = D0_5
      
      REAL (KIND=dbl_kind), DIMENSION(nk2) :: RKV,TEMP_ZX,TEMP_ZY
      REAL (KIND=dbl_kind) :: TEMP_ZZ
      
      REAL (KIND=dbl_kind) :: COEFX_K,COEFY_K,COEFA_K
      REAL (KIND=dbl_kind) :: COEFX_E,COEFY_E,COEFA_E
      REAL (KIND=dbl_kind) :: COEFX_Z,COEFY_Z,COEFA_Z 

      REAL (KIND=dbl_kind) :: KVAL_P,KVAL_M
      REAL (KIND=dbl_kind) :: FXP,FXM,FYP,FYM
      REAL (KIND=dbl_kind) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=dbl_kind) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0
      
      REAL (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE :: KVAL

      INTEGER :: I,J,K,KLOWU,KLOWV
      INTEGER :: KLOWQ_MAX,num_seg,mi1,mj1,mim_c,mip_c,mjm_c,mjp_c 

!**************************************************************************** 
      DO num_seg = 1, 4
!**************************************************************************** 
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB
      
      mi1    = channel%seg(num_seg)%mi1  ! x-size of channel segment  
      mim_c  = channel%seg(num_seg)%mim_c
      mip_c  = channel%seg(num_seg)%mip_c  
    
      mj1    = channel%seg(num_seg)%mj1  ! y-size of channel segment   
      mjm_c  = channel%seg(num_seg)%mjm_c     
      mjp_c  = channel%seg(num_seg)%mjp_c    

      ALLOCATE(KVAL(mim_c:mip_c,mjm_c:mjp_c,nk2)) 

!     Calculate the diffusion coefficient      
      DO J = mjm_c, mjp_c
      DO I = mim_c, mip_c
       KLOW = channel%seg(num_seg)%KLOWQ_IJ(I,J)

       DO K = KLOW, NK2
         RKV(K) = D0_25*(channel%seg(num_seg)%DEFXY(I-1,J-1,K)**2 &
                        +channel%seg(num_seg)%DEFXY(I-1,J,K)**2   &
                        +channel%seg(num_seg)%DEFXY(I,J-1,K)**2   &
                        +channel%seg(num_seg)%DEFXY(I,J,K)**2)    &
        +(channel%seg(num_seg)%DEFXX(I,J,K)-channel%seg(num_seg)%DEFYY(I,J,K))**2
       ENDDO
       DO K = 1, KLOW-1
         RKV(K) = D0_0
       ENDDO
       
       DO K = 1, NK2
         KVAL(I,J,K) = FAC_NLD*channel%seg(num_seg)%RG_T(I,J)*DXSQ*SQRT(RKV(K))
       ENDDO

       DO K = 1, NK2
         KVAL(I,J,K) = MIN(KVAL(I,J,K),channel%seg(num_seg)%RG_T(I,J)*DXSQ*CRITV/DT)
       ENDDO

      ENDDO  
      ENDDO      

      DO J = 1, mj1
      DO I = 1, mi1
      
      COEFX_K = channel%seg(num_seg)%COEFX_K(I,J)
      COEFY_K = channel%seg(num_seg)%COEFY_K(I,J)
      COEFA_K = channel%seg(num_seg)%COEFA_K(I,J)

      COEFX_E = channel%seg(num_seg)%COEFX_E(I,J)
      COEFY_E = channel%seg(num_seg)%COEFY_E(I,J)
      COEFA_E = channel%seg(num_seg)%COEFA_E(I,J)

      COEFX_Z = channel%seg(num_seg)%COEFX_Z(I,J)
      COEFY_Z = channel%seg(num_seg)%COEFY_Z(I,J)
      COEFA_Z = channel%seg(num_seg)%COEFA_Z(I,J)
      
!-----------------------------------------------------------
!                                 Z3DX
!-----------------------------------------------------------
       KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

       DO K = 1, KLOWV-1
        TEMP_ZX(K) = D0_0
       ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
     DO K = KLOWV, KLOWQ_MAX-1
!     Influenced by topography

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        KVAL_P = (KVAL(I,J,K)+KVAL(I+1,J,K)+KVAL(I,J+1,K)       &
                 +KVAL(I+1,J+1,K)+KVAL(I,J,K+1)+KVAL(I+1,J,K+1) &
                 +KVAL(I,J+1,K+1)+KVAL(I+1,J+1,K+1))/D8_0 
                 
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = (KVAL(I-1,J,K)+KVAL(I,J,K)+KVAL(I-1,J+1,K)     &
                 +KVAL(I,J+1,K)+KVAL(I-1,J,K+1)+KVAL(I,J,K+1)   &
                 +KVAL(I-1,J+1,K+1)+KVAL(I,J+1,K+1))/D8_0 
                 
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1 .OR.   &
          channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = D0_5*(KVAL(I,J+1,K)+KVAL(I,J+1,K+1))

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = D0_5*(KVAL(I,J,K)+KVAL(I,J,K+1))

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1 .OR.    &
          channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)
       ENDIF

       TEMP_ZX(K) = TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J  ,K),        & 
                           channel%seg(num_seg)%Z3DX(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J+1,K))                                                  

      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

        KVAL_P = (KVAL(I,J,K)+KVAL(I+1,J,K)+KVAL(I,J+1,K)         &
                 +KVAL(I+1,J+1,K)+KVAL(I,J,K+1)+KVAL(I+1,J,K+1)   &
                 +KVAL(I,J+1,K+1)+KVAL(I+1,J+1,K+1))/D8_0 
                 
        KVAL_M = (KVAL(I-1,J,K)+KVAL(I,J,K)+KVAL(I-1,J+1,K)       &
                 +KVAL(I,J+1,K)+KVAL(I-1,J,K+1)+KVAL(I,J,K+1)     &
                 +KVAL(I-1,J+1,K+1)+KVAL(I,J+1,K+1))/D8_0 

          FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)*KVAL_P
          FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
          FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)*KVAL_P

          FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)*KVAL_M
          FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)*KVAL_M
          FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)*KVAL_M

          FYP_X0 = channel%seg(num_seg)%RGG_T(2,I,J+1)*(KVAL_P-KVAL_M)
          FYM_X0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)

        KVAL_P = D0_5*(KVAL(I,J+1,K)+KVAL(I,J+1,K+1))
        KVAL_M = D0_5*(KVAL(I,J,K)+KVAL(I,J,K+1))

          FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)*KVAL_P
          FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)*KVAL_P
          FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)*KVAL_P

          FYM    = channel%seg(num_seg)%RGG_T(3,I,J)*KVAL_M
          FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)*KVAL_M
          FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

          FXP_Y0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
          FXM_Y0 = channel%seg(num_seg)%RGG_Z(2,I-1,J)*(KVAL_P-KVAL_M)

       TEMP_ZX(K) = TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J,K),          &
                           channel%seg(num_seg)%Z3DX(I,J,K),            &
                           channel%seg(num_seg)%Z3DX(I+1,J,K),          &
                           channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J+1,K))                                  

      ENDDO  ! K-LOOP (high level: mountain-free region)

      DO K = 2, NK1
        channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K) &
                                         + DT*TEMP_ZX(K)
      ENDDO
!-----------------------------------------------------------
!                                Z3DY
!-----------------------------------------------------------
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = 1, KLOWU-1
        TEMP_ZY(K) = D0_0
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO K = KLOWU, KLOWQ_MAX-1
!     Influenced by topography

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        KVAL_P = D0_5*(KVAL(I+1,J,K)+KVAL(I+1,J,K+1))
                     
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        KVAL_M = D0_5*(KVAL(I,J,K)+KVAL(I,J,K+1))
                     
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FYP_X0 = D0_0
        FYM_X0 = D0_0
       ELSE
        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        KVAL_P = (KVAL(I,J,K)+KVAL(I,J+1,K)+KVAL(I+1,J,K)        &
                 +KVAL(I+1,J+1,K)+KVAL(I,J,K+1)+KVAL(I,J+1,K+1)  &
                 +KVAL(I+1,J,K+1)+KVAL(I+1,J+1,K+1))/D8_0     
                     
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        KVAL_M = (KVAL(I,J-1,K)+KVAL(I,J,K)+KVAL(I+1,J-1,K)       &
                 +KVAL(I+1,J,K)+KVAL(I,J-1,K+1)+KVAL(I,J,K+1)     &
                 +KVAL(I+1,J-1,K+1)+KVAL(I+1,J,K+1))/D8_0      
                 
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        
       ENDIF

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1 .OR. &
          channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FXP_Y0 = D0_0
        FXM_Y0 = D0_0
       ELSE
        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)
       ENDIF

       TEMP_ZY(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))         
                           
      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

        KVAL_P = D0_5*(KVAL(I+1,J,K)+KVAL(I+1,J,K+1))
                     
        KVAL_M = D0_5*(KVAL(I,J,K)+KVAL(I,J,K+1))

        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_Z(2,I,J)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_Z(2,I,J-1)*(KVAL_P-KVAL_M)

        KVAL_P = (KVAL(I,J,K)+KVAL(I,J+1,K)+KVAL(I+1,J,K)          &
                 +KVAL(I+1,J+1,K)+KVAL(I,J,K+1)+KVAL(I,J+1,K+1)    &
                 +KVAL(I+1,J,K+1)+KVAL(I+1,J+1,K+1))/D8_0   
                       
        KVAL_M = (KVAL(I,J-1,K)+KVAL(I,J,K)+KVAL(I+1,J-1,K)        &
                 +KVAL(I+1,J,K)+KVAL(I,J-1,K+1)+KVAL(I,J,K+1)      &
                 +KVAL(I+1,J-1,K+1)+KVAL(I+1,J,K+1))/D8_0  

        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)*KVAL_M                                        

        FXP_Y0 = channel%seg(num_seg)%RGG_T(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_T(2,I,J)*(KVAL_P-KVAL_M)                                   


       TEMP_ZY(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,    &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I-1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I,J,K),           &
                           channel%seg(num_seg)%Z3DY(I+1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))         

      ENDDO  ! K-LOOP (high level: mountain-free region)

      DO K = 2, NK1
        channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K) &
                                         + DT*TEMP_ZY(K)
      ENDDO

!-----------------------------------------------------------
!                           Z3DZ (TOP LAYER)
!-----------------------------------------------------------
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      KVAL_P = D0_5*(KVAL(I+1,J,NK2)+KVAL(I+1,J+1,NK2))
      KVAL_M = D0_5*(KVAL(I,J,NK2)+KVAL(I,J+1,NK2))

        FXP    = channel%seg(num_seg)%RGG_V(1,I+1,J)*KVAL_P
        FYP_XP = channel%seg(num_seg)%RGG_U(2,I+1,J+1)*KVAL_P
        FYM_XP = channel%seg(num_seg)%RGG_U(2,I+1,J)*KVAL_P

        FXM    = channel%seg(num_seg)%RGG_V(1,I,J)*KVAL_M
        FYP_XM = channel%seg(num_seg)%RGG_U(2,I-1,J+1)*KVAL_M
        FYM_XM = channel%seg(num_seg)%RGG_U(2,I-1,J)*KVAL_M

        FYP_X0 = channel%seg(num_seg)%RGG_U(2,I,J+1)*(KVAL_P-KVAL_M)
        FYM_X0 = channel%seg(num_seg)%RGG_U(2,I,J)*(KVAL_P-KVAL_M)

      KVAL_P = D0_5*(KVAL(I,J+1,NK2)+KVAL(I+1,J+1,NK2))
      KVAL_M = D0_5*(KVAL(I,J,NK2)+KVAL(I+1,J,NK2))

        FYP    = channel%seg(num_seg)%RGG_U(3,I,J+1)*KVAL_P
        FXP_YP = channel%seg(num_seg)%RGG_V(2,I+1,J+1)*KVAL_P
        FXM_YP = channel%seg(num_seg)%RGG_V(2,I,J+1)*KVAL_P

        FYM    = channel%seg(num_seg)%RGG_U(3,I,J)*KVAL_M
        FXP_YM = channel%seg(num_seg)%RGG_V(2,I+1,J-1)*KVAL_M
        FXM_YM = channel%seg(num_seg)%RGG_V(2,I,J-1)*KVAL_M

        FXP_Y0 = channel%seg(num_seg)%RGG_V(2,I+1,J)*(KVAL_P-KVAL_M)
        FXM_Y0 = channel%seg(num_seg)%RGG_V(2,I,J)*(KVAL_P-KVAL_M)

        TEMP_ZZ = TURB_H(COEFX_Z,COEFY_Z,COEFA_Z,FXP,FXM,FYP,FYM,    &
                         FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                         FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                         channel%seg(num_seg)%Z3DZ(I-1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J-1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I-1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I,J,NK2),         &
                         channel%seg(num_seg)%Z3DZ(I+1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I-1,J+1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J+1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J+1,NK2))  

        channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                                           + DT*TEMP_ZZ

      ENDDO   ! I-loop
      ENDDO   ! J-loop

      DEALLOCATE(KVAL) 
      
!******************************  
      ENDDO   ! num_seg 
!****************************** 

   END SUBROUTINE hdiffusion_smagorinsky

!=======================================================================   
      SUBROUTINE LDIFFUSION2 (channel)
!=======================================================================      
!     Numerical horizontal diffusion (linear 2nd-order)

      type(channel_t), intent(inout) :: channel   ! channel data
      
      ! Local variables
      REAL (KIND=dbl_kind) :: COEF_FAC = 400._dbl_kind

      REAL (KIND=dbl_kind),DIMENSION(nk2) :: TEMP_ZX,TEMP_ZY
      REAL (KIND=dbl_kind) :: TEMP_ZZ

      REAL (KIND=dbl_kind) :: COEFX_K,COEFY_K,COEFA_K
      REAL (KIND=dbl_kind) :: COEFX_E,COEFY_E,COEFA_E
      REAL (KIND=dbl_kind) :: COEFX_Z,COEFY_Z,COEFA_Z

      REAL (KIND=dbl_kind) :: KVAL,KVAL0
      REAL (KIND=dbl_kind) :: FXP,FXM,FYP,FYM
      REAL (KIND=dbl_kind) :: FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0
      REAL (KIND=dbl_kind) :: FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0
      
      INTEGER :: I,J,K,KLOWQ_MAX,KLOWU,KLOWV,num_seg,mi1,mj1

      KVAL0 = DXSQ/(4.*DT)     ! 2dx-wave is removed over one timestep  
      KVAL  = KVAL0/COEF_FAC   ! 2dx-wave is removed over COEF_FAC timesteps 

!*************************************************************************   
      DO num_seg = 1, 4
!************************************************************************* 
      KLOWQ_MAX = channel%seg(num_seg)%KLOWQ_MAX_GLOB
      
      mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
      mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 

      DO J = 1, mj1
      DO I = 1, mi1

      COEFX_K = channel%seg(num_seg)%COEFX_K(I,J)
      COEFY_K = channel%seg(num_seg)%COEFY_K(I,J)
      COEFA_K = channel%seg(num_seg)%COEFA_K(I,J)

      COEFX_E = channel%seg(num_seg)%COEFX_E(I,J)
      COEFY_E = channel%seg(num_seg)%COEFY_E(I,J)
      COEFA_E = channel%seg(num_seg)%COEFA_E(I,J)

      COEFX_Z = channel%seg(num_seg)%COEFX_Z(I,J)
      COEFY_Z = channel%seg(num_seg)%COEFY_Z(I,J)
      COEFA_Z = channel%seg(num_seg)%COEFA_Z(I,J)
!-----------------------------------------------------------
!                                 Z3DX
!-----------------------------------------------------------
      KLOWV = channel%seg(num_seg)%KLOWV_IJ(I,J)

      DO K = 1, KLOWV-1
       TEMP_ZX(K) = D0_0
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO K = KLOWV, KLOWQ_MAX-1
!     Influenced by topography

       IF(channel%seg(num_seg)%DM_KSI_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)
       ENDIF

        FYP_X0 = D0_0
        FYM_X0 = D0_0
       

       IF(channel%seg(num_seg)%DM_KSI_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)
       ENDIF

       IF(channel%seg(num_seg)%DM_KSI_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)
       ENDIF

        FXP_Y0 = D0_0
        FXM_Y0 = D0_0

       TEMP_ZX(K) = TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J  ,K),        & 
                           channel%seg(num_seg)%Z3DX(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J+1,K))                                       

      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

        FXP    = channel%seg(num_seg)%RGG_Z(1,I,J)
        FYP_XP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)
        FYM_XP = channel%seg(num_seg)%RGG_T(2,I+1,J)

        FXM    = channel%seg(num_seg)%RGG_Z(1,I-1,J)
        FYP_XM = channel%seg(num_seg)%RGG_T(2,I-1,J+1)
        FYM_XM = channel%seg(num_seg)%RGG_T(2,I-1,J)

        FYP_X0 = D0_0
        FYM_X0 = D0_0

        FYP    = channel%seg(num_seg)%RGG_T(3,I,J+1)
        FXP_YP = channel%seg(num_seg)%RGG_Z(2,I,J+1)
        FXM_YP = channel%seg(num_seg)%RGG_Z(2,I-1,J+1)

        FYM    = channel%seg(num_seg)%RGG_T(3,I,J)
        FXP_YM = channel%seg(num_seg)%RGG_Z(2,I,J-1)
        FXM_YM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)

        FXP_Y0 = D0_0
        FXM_Y0 = D0_0

       TEMP_ZX(K) = TURB_H(COEFX_K,COEFY_K,COEFA_K,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DX(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DX(I-1,J,K),          &
                           channel%seg(num_seg)%Z3DX(I,J,K),            &
                           channel%seg(num_seg)%Z3DX(I+1,J,K),          &
                           channel%seg(num_seg)%Z3DX(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DX(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DX(I+1,J+1,K))                                     

      ENDDO  ! K-LOOP (high level: mountain-free region)

      DO K = 2, NK1
        channel%seg(num_seg)%Z3DX(I,J,K) = channel%seg(num_seg)%Z3DX(I,J,K) &
                                         + DT*KVAL*TEMP_ZX(K)
      ENDDO
!-----------------------------------------------------------
!                                Z3DY
!-----------------------------------------------------------
      KLOWU = channel%seg(num_seg)%KLOWU_IJ(I,J)

      DO K = 1, KLOWU-1
       TEMP_ZY(K) = D0_0
      ENDDO

!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
      DO K = KLOWU, KLOWQ_MAX-1
!     Influenced by topography

       IF(channel%seg(num_seg)%DM_ETA_TE(I,J,K).EQ.1) THEN
        FXP    = D0_0
        FYP_XP = D0_0
        FYM_XP = D0_0
       ELSE
        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TW(I,J,K).EQ.1) THEN
        FXM    = D0_0
        FYP_XM = D0_0
        FYM_XM = D0_0
       ELSE
        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)
       ENDIF

        FYP_X0 = D0_0
        FYM_X0 = D0_0

       IF(channel%seg(num_seg)%DM_ETA_TN(I,J,K).EQ.1) THEN
        FYP    = D0_0
        FXP_YP = D0_0
        FXM_YP = D0_0
       ELSE
        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)
       ENDIF
       IF(channel%seg(num_seg)%DM_ETA_TS(I,J,K).EQ.1) THEN
        FYM    = D0_0
        FXP_YM = D0_0
        FXM_YM = D0_0
       ELSE
        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)                                       
       ENDIF

        FXP_Y0 = D0_0
        FXM_Y0 = D0_0

       TEMP_ZY(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,     &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,   &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,   &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I,J  ,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J  ,K),        &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),        &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),          &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))        
                           
      ENDDO  ! K-LOOP (Low level: mountainous region)

      DO K = KLOWQ_MAX, NK1
!     Not influenced by topography

        FXP    = channel%seg(num_seg)%RGG_T(1,I+1,J)
        FYP_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J)
        FYM_XP = channel%seg(num_seg)%RGG_Z(2,I+1,J-1)

        FXM    = channel%seg(num_seg)%RGG_T(1,I,J)
        FYP_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J)
        FYM_XM = channel%seg(num_seg)%RGG_Z(2,I-1,J-1)

        FYP_X0 = D0_0
        FYM_X0 = D0_0   

        FYP    = channel%seg(num_seg)%RGG_Z(3,I,J)
        FXP_YP = channel%seg(num_seg)%RGG_T(2,I+1,J+1)
        FXM_YP = channel%seg(num_seg)%RGG_T(2,I,J+1)

        FYM    = channel%seg(num_seg)%RGG_Z(3,I,J-1)
        FXP_YM = channel%seg(num_seg)%RGG_T(2,I+1,J-1)
        FXM_YM = channel%seg(num_seg)%RGG_T(2,I,J-1)                                      

        FXP_Y0 = D0_0
        FXM_Y0 = D0_0                               

       TEMP_ZY(K) = TURB_H(COEFX_E,COEFY_E,COEFA_E,FXP,FXM,FYP,FYM,    &
                           FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                           FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                           channel%seg(num_seg)%Z3DY(I-1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J-1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J-1,K),       &
                           channel%seg(num_seg)%Z3DY(I-1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I,J,K),           &
                           channel%seg(num_seg)%Z3DY(I+1,J,K),         &
                           channel%seg(num_seg)%Z3DY(I-1,J+1,K),       &
                           channel%seg(num_seg)%Z3DY(I,J+1,K),         &
                           channel%seg(num_seg)%Z3DY(I+1,J+1,K))              

      ENDDO  ! K-LOOP (high level: mountain-free region)

      DO K = 2, NK1
        channel%seg(num_seg)%Z3DY(I,J,K) = channel%seg(num_seg)%Z3DY(I,J,K) &
                                         + DT*KVAL*TEMP_ZY(K)
      ENDDO

!-----------------------------------------------------------
!                           Z3DZ (TOP LAYER)
!-----------------------------------------------------------
!-------------------------
! HORIZONTAL DIFFUSION
!-------------------------
        FXP    = channel%seg(num_seg)%RGG_V(1,I+1,J)
        FYP_XP = channel%seg(num_seg)%RGG_U(2,I+1,J+1)
        FYM_XP = channel%seg(num_seg)%RGG_U(2,I+1,J)

        FXM    = channel%seg(num_seg)%RGG_V(1,I,J)
        FYP_XM = channel%seg(num_seg)%RGG_U(2,I-1,J+1)
        FYM_XM = channel%seg(num_seg)%RGG_U(2,I-1,J)

        FYP_X0 = D0_0
        FYM_X0 = D0_0

        FYP    = channel%seg(num_seg)%RGG_U(3,I,J+1)
        FXP_YP = channel%seg(num_seg)%RGG_V(2,I+1,J+1)
        FXM_YP = channel%seg(num_seg)%RGG_V(2,I,J+1)

        FYM    = channel%seg(num_seg)%RGG_U(3,I,J)
        FXP_YM = channel%seg(num_seg)%RGG_V(2,I+1,J-1)
        FXM_YM = channel%seg(num_seg)%RGG_V(2,I,J-1)

        FXP_Y0 = D0_0
        FXM_Y0 = D0_0

        TEMP_ZZ = TURB_H(COEFX_Z,COEFY_Z,COEFA_Z,FXP,FXM,FYP,FYM,    &
                         FYP_XP,FYM_XP,FYP_XM,FYM_XM,FYP_X0,FYM_X0,  &
                         FXP_YP,FXM_YP,FXP_YM,FXM_YM,FXP_Y0,FXM_Y0,  &
                         channel%seg(num_seg)%Z3DZ(I-1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J-1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J-1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I-1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I,J,NK2),         &
                         channel%seg(num_seg)%Z3DZ(I+1,J,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I-1,J+1,NK2),     &
                         channel%seg(num_seg)%Z3DZ(I,J+1,NK2),       &
                         channel%seg(num_seg)%Z3DZ(I+1,J+1,NK2))   

        channel%seg(num_seg)%Z3DZ(I,J,NK2) = channel%seg(num_seg)%Z3DZ(I,J,NK2) &
                                           + DT*KVAL*TEMP_ZZ
      
      ENDDO   ! I-loop
      ENDDO   ! J-loop

!******************************  
      ENDDO   ! num_seg 
!****************************** 

   END SUBROUTINE ldiffusion2

!=======================================================================   
      FUNCTION TURB_H(CX,CY,CA,KXP,KXM,KYP,KYM, &
                      KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0, &
                      KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0, &
                      A1,A2,A3,A4,A5,A6,A7,A8,A9) &
               RESULT(TURB_H_result)
!=======================================================================                                  
!     Calculate the turbulent effects (horizontal)
!     CX, CY, CA: coefficients
!     KXP and others: C (see the document)
!     A: target variable
!     A7(I-1,J+1,K),A8(I,J+1,K),A9(I+1,J+1,K)
!     A4(I-1,J  ,K),A5(I,J  ,K),A6(I+1,J  ,K)
!     A1(I-1,J-1,K),A2(I,J-1,K),A3(I+1,J-1,K)
!------------------------------------------------------------------

      USE shr_kind_mod, only: dbl_kind => shr_kind_r8
      IMPLICIT NONE

      REAL (KIND=dbl_kind), INTENT(IN) :: CX,CY,CA,KXP,KXM,KYP,KYM
      REAL (KIND=dbl_kind), INTENT(IN) :: KYP_XP,KYM_XP,KYP_XM,KYM_XM,KYP_X0,KYM_X0
      REAL (KIND=dbl_kind), INTENT(IN) :: KXP_YP,KXM_YP,KXP_YM,KXM_YM,KXP_Y0,KXM_Y0
      REAL (KIND=dbl_kind), INTENT(IN) :: A1,A2,A3,A4,A5,A6,A7,A8,A9
      REAL (KIND=dbl_kind) :: TURB_H_result

      TURB_H_result = CX*(KXP*(A6-A5)-KXM*(A5-A4))      &
                    + CY*(KYP*(A8-A5)-KYM*(A5-A2))      &
                    + CA*(KYP_XP*(A9-A6)+KYM_XP*(A6-A3) &
                         +KYP_X0*(A8-A5)+KYM_X0*(A5-A2) &
                         -KYP_XM*(A7-A4)-KYM_XM*(A4-A1) &
                         +KXP_YP*(A9-A8)+KXM_YP*(A8-A7) &
                         +KXP_Y0*(A6-A5)+KXM_Y0*(A5-A4) &
                         -KXP_YM*(A3-A2)-KXM_YM*(A2-A1))

      END FUNCTION TURB_H

END MODULE vort_3d_module
