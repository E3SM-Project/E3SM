MODULE q3d_bg_module
! Contains the programs related to calculating the background fields

      USE shr_kind_mod,    only: r8 => shr_kind_r8
      USE vGCM_data_types, only: vGCM_state_t
      USE vvm_data_types,  only: channel_t

      USE parmsld, only: nhalo_vgcm,nvgcm_seg,nlevel, &
                         ntracer,nk1,nk2,nhalo,netsz,netsz_sq
      USE constld, only: dx,dy,dz,fnz,physics
      
IMPLICIT NONE
PRIVATE

PUBLIC :: cal_bg

CONTAINS

!===================================================================================
   SUBROUTINE CAL_BG (vGCM_st, channel)
!===================================================================================
!  calculate background fields (pay attention to the vertical size of Z3DZ_bg)

!  Check out the consistency between Ucon (Vcon) & U3DX (U3DY) and cov. components
!  Assume that the vertical levels are common between vGCM_st (1:nk1) and channel (2:nk2)
!  Vertical interpolation should be considered.
!  Change the vertical range of w3d for interpolation.

   type(vGCM_state_t), intent(in)    :: vGCM_st(4)  ! vGCM state
   type(channel_t),    intent(inout) :: channel     ! channel data
   
   ! Local
   integer num_seg,mi1,mj1,mim,mip,mim_a,mip_a,mjm,mjp,mjm_a,mjp_a
   integer nt,I,J,K  

!JUNG #ifdef FIXTHIS

!********************************************************************  
   DO num_seg = 1, 4
!******************************************************************** 
   mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
   mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment 
   
   mim   = channel%seg(num_seg)%mim    
   mip   = channel%seg(num_seg)%mip 
   mim_a = channel%seg(num_seg)%mim_a 
   mip_a = channel%seg(num_seg)%mip_a    
   mjm   = channel%seg(num_seg)%mjm    
   mjp   = channel%seg(num_seg)%mjp 
   mjm_a = channel%seg(num_seg)%mjm_a
   mjp_a = channel%seg(num_seg)%mjp_a    

!====================================================================
   IF (mi1.GT.mj1) THEN
!================================================== x-channel segment   

!-----------------------------------------------------------------
!  1. Prepare background fields on CRM grids (interpolation)
!-----------------------------------------------------------------
   ! 1st parameter: .T. positive definite variable 
   ! 2nd parameter: variable location (1: q-point, 2: u-point, 3: v-point)
   ! 3rd parameter: # of halos to be prepared along the channel segment
   !                (# of halos to be prepared across the channel is fixed as nhalo)
 
   ! Thermodynamic variables (No need to get halos in x-direction)   
   CALL GCM2Q3D_X (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%TH(:,:,1:nk1),                              &
                   channel%seg(num_seg)%TH3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_X (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QV(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QV3D_bg(:,:,2:nk2))
    
   if (physics) then                
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QC(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QC3D_bg(:,:,2:nk2))
                   
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QI(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QI3D_bg(:,:,2:nk2))
                   
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QR(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QR3D_bg(:,:,2:nk2))
                   
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QS(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QS3D_bg(:,:,2:nk2))
                   
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QG(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QG3D_bg(:,:,2:nk2))
   endif      
    
   DO NT = 1,ntracer          
   CALL GCM2Q3D_X (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QT(:,:,1:nk1,nt),                           &
                   channel%seg(num_seg)%QT3D_bg(:,:,2:nk2,nt))
   ENDDO                

   ! Contravariant components of horizontal wind (k= 2:nk1): used for eddy calculation          
   CALL GCM2Q3D_X (.FALSE.,2,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucon(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DX_bg(:,:,2:nk1))
                   
   CALL GCM2Q3D_X (.FALSE.,3,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcon(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DY_bg(:,:,2:nk1))

   ! Contravariant components of horizontal wind (k= nk2): used for top-wind calculation
   CALL GCM2Q3D_X (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucon(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DX_bg(:,:,nk2:nk2))
                   
   CALL GCM2Q3D_X (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcon(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DY_bg(:,:,nk2:nk2))
                   
   ! Covariant components of horizontal wind (k= 2:nk1): used for vorticity calculation                 
   CALL GCM2Q3D_X (.FALSE.,2,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucov(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DX_co_bg(:,:,2:nk1))
                   
   CALL GCM2Q3D_X (.FALSE.,3,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcov(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DY_co_bg(:,:,2:nk1))
                   
   ! Covariant components of horizontal wind (k= nk2): used for top-wind calculation                  
   CALL GCM2Q3D_X (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucov(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DX_co_bg(:,:,nk2:nk2))
                   
   CALL GCM2Q3D_X (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcov(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DY_co_bg(:,:,nk2:nk2))
   
   ! Vertical velocity: used for vorticity calculation                 
   CALL GCM2Q3D_X (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%W(:,:,1:nk1),                               &
                   channel%seg(num_seg)%W3D_bg(:,:,2:nk2))

!--------------------------------------------------------------------------------
!  2. Calculate the vorticity components from the covariant components of wind
!--------------------------------------------------------------------------------

!---------------------------- 
!     Z3DX_bg calculation
!----------------------------
      DO K = 2, NK1
       DO J = mjm, mjp_a
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DX_bg(I,J,K) = &
                    (channel%seg(num_seg)%W3D_bg(I,J+1,K)      &
                    -channel%seg(num_seg)%W3D_bg(I,J,K))/DY    &
                  - (channel%seg(num_seg)%U3DY_CO_bg(I,J,K+1)  &
                    -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))*FNZ(K)/DZ         
        ENDDO
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = mjm, mjp_a
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DX_bg(I,J,K) = &
           channel%seg(num_seg)%Z3DX_bg(I,J,K)/channel%seg(num_seg)%RG_V(I,J)
        ENDDO
       ENDDO
      ENDDO

       ! computational boundary condition
       DO J = mjm, mjp_a
        DO I = 1, mi1 
         channel%seg(num_seg)%Z3DX_bg(I,J,  1) = 0.0_r8
         channel%seg(num_seg)%Z3DX_bg(I,J,nk2) = 0.0_r8       
       ENDDO
      ENDDO
 
      ! filling the halo data
      DO K = 1, NK2
       DO I = 1, mi1
         channel%seg(num_seg)%Z3DX_bg(I,mjp,K) = channel%seg(num_seg)%Z3DX_bg(I,mjp_a,K)
       ENDDO
      ENDDO   
!----------------------------                     
!     Z3DY_bg calculation
!----------------------------
      DO K = 2, NK1
       DO J = mjm, mjp
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DY_bg(I,J,K) = &
                  - (channel%seg(num_seg)%W3D_bg(I+1,J,K)      &
                    -channel%seg(num_seg)%W3D_bg(I,J,K))/DX    &
                  + (channel%seg(num_seg)%U3DX_CO_bg(I,J,K+1)  &
                    -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))*FNZ(K)/DZ                 
        ENDDO              
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = mjm, mjp
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DY_bg(I,J,K) = &
           channel%seg(num_seg)%Z3DY_bg(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
        ENDDO
       ENDDO
      ENDDO        

       ! computational boundary condition
       DO J = mjm, mjp
        DO I = 1, mi1
         channel%seg(num_seg)%Z3DY_bg(I,J,  1) = 0.0_r8
         channel%seg(num_seg)%Z3DY_bg(I,J,nk2) = 0.0_r8       
        ENDDO
       ENDDO
!----------------------------------------      
!     Z3DZ_bg calculation (only at nk2)   
!----------------------------------------
      DO J = mjm_a, mjp_a
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,NK2)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,NK2))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,NK2)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,NK2))/DY
       ENDDO
      ENDDO   
      DO J = mjm_a, mjp_a
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,1)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
      ENDDO   

      ! filling the halo data
      DO I = 1, mi1
       channel%seg(num_seg)%Z3DZ_bg(I,mjm,1) = channel%seg(num_seg)%Z3DZ_bg(I,mjm_a,1)
       channel%seg(num_seg)%Z3DZ_bg(I,mjp,1) = channel%seg(num_seg)%Z3DZ_bg(I,mjp_a,1)
      ENDDO

!====================================================================   
   ELSE
!================================================== y-channel segment   

!-----------------------------------------------------------------
!  1. Prepare background fields on CRM grids (interpolation)
!----------------------------------------------------------------- 
   ! 1st parameter: .T. positive definite variable 
   ! 2nd parameter: variable location (1: q-point, 2: u-point, 3: v-point)
   ! 3rd parameter: # of halos to be prepared along the channel segment
   !                (# of halos to be prepared across the channel is fixed as nhalo)
  
   ! Thermodynamic variables (No need to get halos in the channel direction)    
   CALL GCM2Q3D_Y (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%TH(:,:,1:nk1),                              &
                   channel%seg(num_seg)%TH3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QV(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QV3D_bg(:,:,2:nk2))

   if (physics) then 
   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QC(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QC3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QI(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QI3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QR(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QR3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QS(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QS3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QG(:,:,1:nk1),                              &
                   channel%seg(num_seg)%QG3D_bg(:,:,2:nk2))                     
   endif
 
   DO NT = 1,ntracer    
   CALL GCM2Q3D_Y (.FALSE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%QT(:,:,1:nk1,nt),                           &
                   channel%seg(num_seg)%QT3D_bg(:,:,2:nk2,nt))
   ENDDO
   
   ! Contravariant components of horizontal wind (k= 2:nk1): used for eddy calculation     
   CALL GCM2Q3D_Y (.FALSE.,2,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucon(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DX_bg(:,:,2:nk1))

   CALL GCM2Q3D_Y (.FALSE.,3,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcon(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DY_bg(:,:,2:nk1))

   ! Contravariant components of horizontal wind (k= nk2): used for top-wind calculation
   CALL GCM2Q3D_Y (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucon(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DX_bg(:,:,nk2:nk2))

   CALL GCM2Q3D_Y (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcon(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DY_bg(:,:,nk2:nk2))

   ! Covariant components of horizontal wind (k= 2:nk1): used for vorticity calculation  
   CALL GCM2Q3D_Y (.FALSE.,2,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucov(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DX_co_bg(:,:,2:nk1))

   CALL GCM2Q3D_Y (.FALSE.,3,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcov(:,:,1:nk1-1),                            &
                   channel%seg(num_seg)%U3DY_co_bg(:,:,2:nk1))

   ! Covariant components of horizontal wind (k= nk2): used for top-wind calculation          
   CALL GCM2Q3D_Y (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Ucov(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DX_co_bg(:,:,nk2:nk2))

   CALL GCM2Q3D_Y (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st(num_seg)%Vcov(:,:,nk1:nk1),                            &
                   channel%seg(num_seg)%U3DY_co_bg(:,:,nk2:nk2))

   ! Vertical velocity: used for vorticity calculation
   CALL GCM2Q3D_Y (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st(num_seg)%W(:,:,1:nk1),                               &
                   channel%seg(num_seg)%W3D_bg(:,:,2:nk2))

!---------------------------------------------------------------------------------
!  2. Calculate the vorticity components from the covariant components of wind
!---------------------------------------------------------------------------------    
!----------------------------
!     Z3DX_bg calculation
!----------------------------
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = mim, mip
         channel%seg(num_seg)%Z3DX_bg(I,J,K) = &
                    (channel%seg(num_seg)%W3D_bg(I,J+1,K)      &
                    -channel%seg(num_seg)%W3D_bg(I,J,K))/DY    &
                  - (channel%seg(num_seg)%U3DY_CO_bg(I,J,K+1)  &
                    -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))*FNZ(K)/DZ         
        ENDDO
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = mim, mip
         channel%seg(num_seg)%Z3DX_bg(I,J,K) = &
           channel%seg(num_seg)%Z3DX_bg(I,J,K)/channel%seg(num_seg)%RG_V(I,J)
        ENDDO
       ENDDO
      ENDDO

       ! computational boundary condition
       DO J = 1, mj1
        DO I = mim, mip 
         channel%seg(num_seg)%Z3DX_bg(I,J,  1) = 0.0_r8
         channel%seg(num_seg)%Z3DX_bg(I,J,nk2) = 0.0_r8       
       ENDDO
      ENDDO
 
!----------------------------                     
!     Z3DY_bg calculation
!----------------------------
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = mim, mip_a
         channel%seg(num_seg)%Z3DY_bg(I,J,K) = &
                  - (channel%seg(num_seg)%W3D_bg(I+1,J,K)      &
                    -channel%seg(num_seg)%W3D_bg(I,J,K))/DX    &
                  + (channel%seg(num_seg)%U3DX_CO_bg(I,J,K+1)  &
                    -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))*FNZ(K)/DZ                 
        ENDDO              
       ENDDO
      ENDDO
      DO K = 2, NK1
       DO J = 1, mj1
        DO I = mim, mip_a
         channel%seg(num_seg)%Z3DY_bg(I,J,K) = &
           channel%seg(num_seg)%Z3DY_bg(I,J,K)/channel%seg(num_seg)%RG_U(I,J)
        ENDDO
       ENDDO
      ENDDO        

       ! computational boundary condition
       DO J = 1, mj1
        DO I = mim, mip_a
         channel%seg(num_seg)%Z3DY_bg(I,J,  1) = 0.0_r8
         channel%seg(num_seg)%Z3DY_bg(I,J,nk2) = 0.0_r8       
       ENDDO
      ENDDO
      
      ! filling the halo data
      DO K = 1, NK2
       DO J = 1, mj1
         channel%seg(num_seg)%Z3DY_bg(mip,J,K) = channel%seg(num_seg)%Z3DY_bg(mip_a,J,K)
       ENDDO
      ENDDO   
!---------------------------------------             
!     Z3DZ_bg calculation (only at nk2)   
!---------------------------------------
      DO J = 1, mj1
       DO I = mim_a, mip_a
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,NK2)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,NK2))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,NK2)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,NK2))/DY
       ENDDO
      ENDDO   
      DO J = 1, mj1
       DO I = mim_a, mip_a
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,1)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
      ENDDO   

      ! filling the halo data
      DO J = 1, mj1
       channel%seg(num_seg)%Z3DZ_bg(mim,J,1) = channel%seg(num_seg)%Z3DZ_bg(mim_a,J,1)
       channel%seg(num_seg)%Z3DZ_bg(mip,J,1) = channel%seg(num_seg)%Z3DZ_bg(mip_a,J,1)
      ENDDO
!====================================================================  
   ENDIF
!==================================================================== 

!********************************************************************   
   ENDDO   ! num_seg 
!********************************************************************

!JUNG #endif

   END SUBROUTINE cal_bg

!JUNG #ifdef FIXTHIS
!===================================================================================
   SUBROUTINE GCM2Q3D_X (POSITIVE,LOC_V,mhalo_l,NFACE,mim,mip,mjm,mjp,kdim, &
                         lsta,lend,lcen,A_in,A_bg)
!===================================================================================
!  Obtain background fields of the CRM through an interpolation of vGCM values 
!  (for an x-channel segment) 
!
!  INPUT : vGCM state variable (A_in) with a fixed data shape
!  OUTPUT: corresponding background field (A_bg)

   LOGICAL, INTENT(IN) :: POSITIVE           ! if .T., a positive definite variable
   INTEGER, INTENT(IN) :: LOC_V              ! 1 (q-point)  2 (u-point)  3 (v-point)
   INTEGER, INTENT(IN) :: mhalo_l            ! # of halo along the channel needed
   INTEGER, INTENT(IN) :: NFACE              ! # of cube-face, to which the segment belongs
   INTEGER, INTENT(IN) :: mim,mip,mjm,mjp    ! Horizontal sizes of a channel segment 
   INTEGER, INTENT(IN) :: kdim               ! Vertical size of input data
   INTEGER, INTENT(IN) :: lsta(nVGCM_seg)    ! Starting CRM-grid index in each vGCM cell
   INTEGER, INTENT(IN) :: lend(nVGCM_seg)    ! Ending CRM-grid index in each vGCM cell
   REAL(kind=r8), INTENT(IN) :: lcen(nVGCM_seg)  ! Center CRM-grid index in each vGCM cell
   REAL(kind=r8), INTENT(IN) ::  & 
       A_in(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)
       
   REAL(kind=r8), INTENT(OUT) :: A_bg(mim:mip,mjm:mjp,kdim)    
       
   ! Local variables    
   REAL(kind=r8) :: A(1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,-nhalo_vGCM:nhalo_vGCM,kdim)
   REAL(kind=r8) :: GRADI(0:nVGCM_seg+1,-1:1,kdim)
   REAL(kind=r8) :: GRADJ(0:nVGCM_seg+1,-1:1,kdim)
   REAL(KIND=r8) :: VAL0,VAL1,VAL2,VAL3
   REAL(KIND=r8) :: V01,V02,V03,V04,V05,V06,V07,V08
   REAL(KIND=r8) :: XVAL,YVAL,X0,X_COR,Y_COR
   
   INTEGER :: I,J,K,IPO,JPO,chl,chn,lsta0,lend0,lcen_int
   
   SELECT CASE (LOC_V)
   CASE(1)  
     ! q-point variable
     X_COR = 0.0_r8
     Y_COR = 0.0_r8
   CASE(2)
     ! u-point variable
     X_COR = 0.5_r8
     Y_COR = 0.0_r8
   CASE(3)
     ! v-point variable
     X_COR = 0.0_r8
     Y_COR = 0.5_r8
   END SELECT
   
!  Prepare data on the local Cartesian coordinates (for x-channels)
   CALL PREP_X (NFACE, kdim, A_in, A)

!  Prepare the 1st-order horizontal gradients
   DO K = 1, kdim
    DO J = -1, 1
     DO I = 0, nVGCM_seg+1
      GRADI(I,J,K) = 0.5_r8*(A(I+1,J,K)-A(I-1,J,K))
      GRADJ(I,J,K) = 0.5_r8*(A(I,J+1,K)-A(I,J-1,K))
     ENDDO
    ENDDO
   ENDDO

!=================================================
   DO I = 1, nVGCM_seg
!=================================================   
    lsta0 = lsta(I)
    lend0 = lend(I)
    
    IF (I.EQ.1) lsta0 = lsta0 - mhalo_l
    IF (I.EQ.nVGCM_seg) lend0 = lend0 + mhalo_l
    
    lcen_int = INT(lcen(I))
    
    X0 = FLOAT(netsz) - lcen(I)

!   --------------------
!   Section 1 (Left/UP)
!   --------------------
    IPO = I-1
    JPO = 0
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)

    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
        
    DO chl = lsta0, lcen_int
     XVAL = X0 + FLOAT(chl) + X_COR
     DO chn = 1, mjp
       YVAL = FLOAT(CHN-1) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))                         
     ENDDO
    ENDDO
    ENDDO 
    
!   --------------------
!   Section 2 (Left/DOWN)
!   --------------------
    IPO = I-1
    JPO = -1
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)
          
    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
    
    DO chl = lsta0, lcen_int
     XVAL = X0 + FLOAT(chl) + X_COR
     DO chn = mjm, 0
       YVAL = FLOAT(CHN-1) + FLOAT(netsz) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))                         
     ENDDO
    ENDDO
    ENDDO     

!   --------------------
!   Section 3 (RIGHT/UP)
!   --------------------
    IPO = I
    JPO = 0
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)

    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
        
    DO chl = lcen_int+1,lend0
     XVAL = FLOAT(chl) - lcen(I) + X_COR
     DO chn = 1, mjp
       YVAL = FLOAT(CHN-1) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))                         
     ENDDO
    ENDDO
    ENDDO 

!   --------------------
!   Section 4 (RIGHT/DOWN)
!   --------------------
    IPO = I
    JPO = -1
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)
          
    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))

    DO chl = lcen_int+1,lend0
     XVAL = FLOAT(chl) - lcen(I) + X_COR
     DO chn = mjm, 0
       YVAL = FLOAT(CHN-1) + FLOAT(netsz) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))                         
     ENDDO
    ENDDO
    ENDDO     

!=================================================
   ENDDO   ! I-loop  
!=================================================   

   CONTAINS
   
   SUBROUTINE PREP_X (NFACE, kdim, A_in, A_out)
!-----------------------------------------------   
!  Prepare data arrange on the local Cartesian coordinates (for x-channels)
!  Assume: nLevel of vGCM data = nk1 of CRM (common vertical coordinates) 
!  Need to consider vertical interpolation later
!-----------------------------------------------   
   INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
   INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data 
   
   REAL(kind=r8), INTENT(IN)  :: &     
       A_in(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)
   REAL(kind=r8), INTENT(OUT) :: &
       A_out(1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,-nhalo_vGCM:nhalo_vGCM,kdim)

   ! Local
   INTEGER :: I,J,K
    
   IF (NFACE .NE. 6) THEN
      DO K = 1, kdim
       DO J = -nhalo_vGCM,nhalo_vGCM
        DO I = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
          A_out(I,J,K) = A_in(-J,I,K)
        ENDDO
       ENDDO 
      ENDDO  
   ELSE
      DO K = 1, kdim
       DO J = -nhalo_vGCM,nhalo_vGCM
        DO I = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
          A_out(I,J,K) = A_in(J,-I,K)
        ENDDO
       ENDDO 
      ENDDO       
   ENDIF
        
   END SUBROUTINE prep_x
   
   END SUBROUTINE gcm2q3d_x

!===================================================================================
   SUBROUTINE GCM2Q3D_Y (POSITIVE,LOC_V,mhalo_l,NFACE,mim,mip,mjm,mjp,kdim, &
                         lsta,lend,lcen,A_in,A_bg)
!===================================================================================
!  Obtain background fields of the CRM through an interpolation of vGCM values 
!  (for an y-channel segment) 
!
!  INPUT : vGCM state variable (A_in) with a fixed data shape
!  OUTPUT: corresponding background field (A_bg)

   LOGICAL, INTENT(IN) :: POSITIVE           ! if .T., a positive definite variable
   INTEGER, INTENT(IN) :: LOC_V              ! 1 (q-point)  2 (u-point)  3 (v-point)
   INTEGER, INTENT(IN) :: mhalo_l            ! # of halo along the channel needed
   INTEGER, INTENT(IN) :: NFACE              ! # of cube-face, to which the segment belongs
   INTEGER, INTENT(IN) :: mim,mip,mjm,mjp    ! Horizontal sizes of a channel segment 
   INTEGER, INTENT(IN) :: kdim               ! Vertical size of input data
   INTEGER, INTENT(IN) :: lsta(nVGCM_seg)    ! Starting CRM-grid index in each vGCM cell
   INTEGER, INTENT(IN) :: lend(nVGCM_seg)    ! Ending CRM-grid index in each vGCM cell
   REAL(kind=r8), INTENT(IN) :: lcen(nVGCM_seg)  ! Center CRM-grid index in each vGCM cell
   REAL(kind=r8), INTENT(IN) ::  & 
       A_in(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)
       
   REAL(kind=r8), INTENT(OUT) :: A_bg(mim:mip,mjm:mjp,kdim)    
       
   ! Local variables    
   REAL(kind=r8) :: A(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)
   REAL(kind=r8) :: GRADI(-1:1,0:nVGCM_seg+1,kdim)
   REAL(kind=r8) :: GRADJ(-1:1,0:nVGCM_seg+1,kdim)
   REAL(KIND=r8) :: VAL0,VAL1,VAL2,VAL3
   REAL(KIND=r8) :: V01,V02,V03,V04,V05,V06,V07,V08
   REAL(KIND=r8) :: XVAL,YVAL,X0,X_COR,Y_COR,Y0
   
   INTEGER :: I,J,K,IPO,JPO,chl,chn,lsta0,lend0,lcen_int
   
   SELECT CASE (LOC_V)
   CASE(1)  
     ! q-point variable
     X_COR = 0.0_r8
     Y_COR = 0.0_r8
   CASE(2)
     ! u-point variable
     X_COR = 0.5_r8
     Y_COR = 0.0_r8
   CASE(3)
     ! v-point variable
     X_COR = 0.0_r8
     Y_COR = 0.5_r8
   END SELECT
   
!  Prepare data on the local Cartesian coordinates (for y-channels)
   CALL PREP_Y (NFACE, kdim, A_in, A)

!  Prepare the 1st-order horizontal gradients
   DO K = 1, kdim
    DO J = 0, nVGCM_seg+1
     DO I = -1, 1
      GRADI(I,J,K) = 0.5_r8*(A(I+1,J,K)-A(I-1,J,K))
      GRADJ(I,J,K) = 0.5_r8*(A(I,J+1,K)-A(I,J-1,K))
     ENDDO
    ENDDO
   ENDDO

!=================================================
   DO J = 1, nVGCM_seg
!=================================================   
    lsta0 = lsta(J)
    lend0 = lend(J)
    
    IF (J.EQ.1) lsta0 = lsta0 - mhalo_l
    IF (J.EQ.nVGCM_seg) lend0 = lend0 + mhalo_l 
    
    lcen_int = INT(lcen(J))
    
    Y0 = FLOAT(netsz) - lcen(J)

!   --------------------
!   Section 1 (Left/DOWN)
!   --------------------
    IPO = -1
    JPO = J-1
     
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)

    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
        
    DO chl = lsta0, lcen_int
     YVAL = Y0 + FLOAT(chl) + Y_COR
     DO chn = mim, 0
       XVAL = FLOAT(CHN-1) + FLOAT(netsz) + X_COR
       A_bg(chn,chl,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chn,chl,K) = MAX(0.0_r8,A_bg(chn,chl,K))                         
     ENDDO
    ENDDO
    ENDDO 
    
!   --------------------
!   Section 2 (RIGHT/DOWN)
!   --------------------
    IPO = 0
    JPO = J-1
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)
          
    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
    
    DO chl = lsta0, lcen_int
     YVAL = Y0 + FLOAT(chl) + Y_COR
     DO chn = 1, mip
       XVAL = FLOAT(CHN-1) + X_COR 
       A_bg(chn,chl,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chn,chl,K) = MAX(0.0_r8,A_bg(chn,chl,K))                         
     ENDDO
    ENDDO
    ENDDO     

!   --------------------
!   Section 3 (LEFT/UP)
!   --------------------
    IPO = -1
    JPO = J
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)

    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
        
    DO chl = lcen_int+1,lend0
     YVAL = FLOAT(chl) - lcen(J) + Y_COR
     DO chn = mim, 0
       XVAL = FLOAT(CHN-1) + FLOAT(netsz) + X_COR
       A_bg(chn,chl,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chn,chl,K) = MAX(0.0_r8,A_bg(chn,chl,K))                         
     ENDDO
    ENDDO
    ENDDO 
    
!   --------------------
!   Section 4 (RIGHT/UP)
!   --------------------
    IPO = 0
    JPO = J
    
    DO K = 1, kdim
    VAL0 = A(IPO,   JPO,   K)
    VAL1 = A(IPO+1, JPO,   K)
    VAL2 = A(IPO+1, JPO+1, K)
    VAL3 = A(IPO,   JPO+1, K)
          
    V01 = (VAL1 - VAL0 - GRADI(IPO+1, JPO,   K))
    V02 = (VAL3 - VAL0 - GRADJ(IPO,   JPO+1, K))
    V03 = (VAL1 - VAL0 - GRADI(IPO,   JPO,   K))
    V04 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO+1, K))
    V05 = (VAL2 - VAL3 - GRADI(IPO,   JPO+1, K))
    V06 = (VAL2 - VAL1 - GRADJ(IPO+1, JPO,   K))
    V07 = (VAL2 - VAL3 - GRADI(IPO+1, JPO+1, K))
    V08 = (VAL3 - VAL0 - GRADJ(IPO,   JPO,   K))
    
    DO chl = lcen_int+1,lend0
     YVAL = FLOAT(chl) - lcen(J) + Y_COR
     DO chn = 1, mip
       XVAL = FLOAT(CHN-1) + X_COR
       A_bg(chn,chl,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)
                              
       IF (POSITIVE) A_bg(chn,chl,K) = MAX(0.0_r8,A_bg(chn,chl,K))                         
     ENDDO
    ENDDO
    ENDDO     

!=================================================
   ENDDO   ! J-loop
!=================================================   

   CONTAINS

!-----------------------------------------------
   SUBROUTINE PREP_Y (NFACE, kdim, A_in, A_out)
!  Prepare data on the local Cartesian coordinates (for y-channels)
!  Assume: nLevel of vGCM data = nk1 of CRM (common vertical coordinates) 
!  Need to consider vertical interpolation later
!-----------------------------------------------   
   INTEGER, INTENT(IN) :: NFACE     ! # of cube-face, to which the data belongs
   INTEGER, INTENT(IN) :: kdim      ! Vertical size of input data 
   
   real(kind=r8), INTENT(IN)  :: &     
       A_in(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)
   real(kind=r8), INTENT(OUT) :: &
       A_out(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,kdim)

   ! Local
   INTEGER :: I,J,K
    
   IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO K = 1, kdim
       DO J = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
        DO I = -nhalo_vGCM,nhalo_vGCM
          A_out(I,J,K) = A_in(-I,-J,K)
        ENDDO
       ENDDO 
      ENDDO  
   ELSE
      DO K = 1, kdim
       DO J = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
        DO I = -nhalo_vGCM,nhalo_vGCM
          A_out(I,J,K) = A_in(I,J,K)
        ENDDO
       ENDDO 
      ENDDO       
   ENDIF
        
   END SUBROUTINE prep_y   

   END SUBROUTINE gcm2q3d_y
   
   FUNCTION QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                    VAL0,VAL1,VAL2,VAL3,N,NSQ) RESULT (QUADINT_result)

   REAL (KIND=r8) :: QUADINT_result
   REAL (KIND=r8), INTENT(IN) :: XVAL,YVAL  ! x and y locations relative to the GCM grid (i,j)
   REAL (KIND=r8), INTENT(IN) :: V01,V02,V03,V04,V05,V06,V07,V08 ! terms in the interpolation ftn
   REAL (KIND=r8), INTENT(IN) :: VAL0,VAL1,VAL2,VAL3 ! terms in the interpolation ftn
   INTEGER, INTENT(IN) :: N,NSQ

   QUADINT_result = &
       (XVAL-N)*(YVAL-N)*(V01*XVAL*XVAL+V02*YVAL*YVAL+VAL0*NSQ) &
      +XVAL*(YVAL-N)*(V03*(XVAL-N)*(XVAL-N)-V04*YVAL*YVAL-VAL1*NSQ) &
      -XVAL*YVAL*(V05*(XVAL-N)*(XVAL-N)+V06*(YVAL-N)*(YVAL-N)-VAL2*NSQ) &
      -(XVAL-N)*YVAL*(V07*XVAL*XVAL-V08*(YVAL-N)*(YVAL-N)+VAL3*NSQ)

   QUADINT_result=QUADINT_result/FLOAT(NSQ*NSQ)

   END FUNCTION quadint
!JUNG #endif

END MODULE q3d_bg_module
