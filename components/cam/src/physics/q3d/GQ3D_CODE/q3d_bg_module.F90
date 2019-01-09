MODULE q3d_bg_module
! Contains the programs related to calculating the background fields

      USE shr_kind_mod,    only: r8 => shr_kind_r8
      USE vGCM_data_types, only: vGCM_state_t,vGCM_map_t
      USE vvm_data_types,  only: channel_t

      USE parmsld, only: nhalo_vgcm,nvgcm_seg,nlevel, &
                         ntracer,nk1,nk2,nk1_ext,nk2_ext,nhalo,netsz,netsz_sq
      USE constld, only: dx,dy,dz,fnz,fnt,pibar,pibarz,zz,zt,rho,rhoz,grav, &
                         cappa,gas_const,pzero,physics,val_missing

      USE utils,   only: indexr,fintrp,spline,splint,con2cov,rll2con,con2rll
      
IMPLICIT NONE
PRIVATE

  type, private :: vGCM_vgrid
    real(kind=r8), allocatable :: Pint(:,:,:) ! pressure at interfaces (used in the extended layers) 
    real(kind=r8), allocatable :: TH(:,:,:)   ! potential temp. (used in the extended layers)  
             
    real(kind=r8), allocatable :: QV(:,:,:)   ! water vapor mixing ratio          
    real(kind=r8), allocatable :: QC(:,:,:)   ! cloud liquid water mixing ratio      
    real(kind=r8), allocatable :: QI(:,:,:)   ! cloud ice water mixing ratio    
    real(kind=r8), allocatable :: QR(:,:,:)   ! rain mixing ratio
    real(kind=r8), allocatable :: QS(:,:,:)   ! snow mixing ratio       
    real(kind=r8), allocatable :: QG(:,:,:)   ! graupel mixing ratio    
    real(kind=r8), allocatable :: QT(:,:,:,:) ! tracer mixing ratio  
    real(kind=r8), allocatable :: Ucon(:,:,:) ! contravariant zonal velocity component  
    real(kind=r8), allocatable :: Vcon(:,:,:) ! contravariant meridional velocity component   
    real(kind=r8), allocatable :: W(:,:,:)    ! vertical velocity
    
    real(kind=r8), allocatable :: PI_G(:,:)   ! pi-value at mid-layer of GCM (at vGCM point)    
    real(kind=r8), allocatable :: ZM_G(:,:)   ! zm-value at mid-layer of GCM (at vGCM point) 
                                       
  end type vGCM_vgrid
    
  type(vGCM_vgrid) :: vGCM_st_new(4)          ! vGCM state on VVM vertical grids
    
PUBLIC :: cal_bg

CONTAINS

! Local Subroutines:
!----------------------------------------------------------------------------------
! Subroutine cal_bg : Calculate background fields
!
!   Subroutine vgcm_prepare : Vertical interpolation of vGCM state variables
!      Subroutine vindex    : Find the vertical grid index for a linear interpolation
!      Subroutine vgcm_layers_to_vvm_layers 
!      Subroutine vgcm_layers_to_vvm_interfaces 
!      Subroutine vgcm_interfaces_to_vvm_interfaces
!   Subroutine gcm2q3d_x    : Horizontal interpolation for x-channel
!      Subroutine prep_x    : data flip
!   Subroutine gcm2q3d_y    : Horizontal interpolation for y-channel
!      Subroutine prep_y    : data flip  
!   Function quadint
!----------------------------------------------------------------------------------

!===================================================================================
   SUBROUTINE CAL_BG (vGCM_st, vGCM_map, channel)
!===================================================================================
!  Calculate background fields from the CAM-SE states. 
!  Pay attention to 
!  the vertical size of Z3DZ_bg = 1
!  the vertical size of TH3D_bg = nk2_ext
!  the vertical size of others  = nk2

   type(vGCM_state_t), intent(in)    :: vGCM_st(4)   ! vGCM state
   type(vGCM_map_t),   intent(in)    :: vGCM_map(4)  ! vGCM map
   type(channel_t),    intent(inout) :: channel      ! channel data
   
   ! Local
   REAL(KIND=r8) :: VMAP(2),TMP
   LOGICAL :: W3D_bg_INTERPOL = .FALSE.  
   
   integer mi1,mim,mip,mim_a,mip_a,mim_c,mip_c
   integer mj1,mjm,mjp,mjm_a,mjp_a,mjm_c,mjp_c
   integer num_seg,nt,I,J,K,nn,IPOM,IPOP,JPOM,JPOP

   DO num_seg = 1, 4
     allocate(vGCM_st_new(num_seg)%Pint(-nhalo_vGCM:nhalo_vGCM, &
                                        1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2_ext))                        
     allocate(vGCM_st_new(num_seg)%TH(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2_ext))
                                      
     allocate(vGCM_st_new(num_seg)%QV(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QC(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QI(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QR(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QS(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QG(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%QT(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2,ntracer))    
     allocate(vGCM_st_new(num_seg)%Ucon(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%Vcon(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2))    
     allocate(vGCM_st_new(num_seg)%W(-nhalo_vGCM:nhalo_vGCM, &
                                      1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nk2)) 
                                                         
     allocate(vGCM_st_new(num_seg)%PI_G(nLevel,nVGCM_seg))   
     allocate(vGCM_st_new(num_seg)%ZM_G(nLevel,nVGCM_seg))     
          
     IF (.NOT.PHYSICS) THEN
       vGCM_st_new(num_seg)%QC(:,:,:) = 0.0_r8
       vGCM_st_new(num_seg)%QI(:,:,:) = 0.0_r8
       vGCM_st_new(num_seg)%QR(:,:,:) = 0.0_r8
       vGCM_st_new(num_seg)%QS(:,:,:) = 0.0_r8
       vGCM_st_new(num_seg)%QG(:,:,:) = 0.0_r8                         
     ENDIF                                      
     
     vGCM_st_new(num_seg)%QT(:,:,:,:) = 0.0_r8                                                                                                                                                                                                                                                     
   ENDDO

   ! Prepare vGCM state variables on vertical grids of VVM (vertical interpolation)
   CALL VGCM_PREPARE (vGCM_st,vGCM_map,vGCM_st_new)
                          
!********************************************************************
   DO num_seg = 1, 4
!********************************************************************
   
   mi1 = channel%seg(num_seg)%mi1  ! x-size of channel segment
   mj1 = channel%seg(num_seg)%mj1  ! y-size of channel segment

   mim   = channel%seg(num_seg)%mim
   mip   = channel%seg(num_seg)%mip
   mim_a = channel%seg(num_seg)%mim_a
   mip_a = channel%seg(num_seg)%mip_a
   mim_c = channel%seg(num_seg)%mim_c
   mip_c = channel%seg(num_seg)%mip_c
      
   mjm   = channel%seg(num_seg)%mjm
   mjp   = channel%seg(num_seg)%mjp
   mjm_a = channel%seg(num_seg)%mjm_a
   mjp_a = channel%seg(num_seg)%mjp_a
   mjm_c = channel%seg(num_seg)%mjm_c
   mjp_c = channel%seg(num_seg)%mjp_c

   channel%seg(num_seg)%PI_G(:,:) = vGCM_st_new(num_seg)%PI_G(:,:) 
   channel%seg(num_seg)%ZM_G(:,:) = vGCM_st_new(num_seg)%ZM_G(:,:) 

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

   CALL GCM2Q3D_X (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim_c,mip_c,mjm_c,mjp_c,nk2_ext, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,              &
                   channel%seg(num_seg)%lcen,                                        &
                   vGCM_st_new(num_seg)%Pint(:,:,1:nk2_ext),                         &
                   channel%seg(num_seg)%Pint_bg(:,:,1:nk2_ext))
                                      
   CALL GCM2Q3D_X (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk2_ext-1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,              &
                   channel%seg(num_seg)%lcen,                                        &
                   vGCM_st_new(num_seg)%TH(:,:,2:nk2_ext),                           &
                   channel%seg(num_seg)%TH3D_bg(:,:,2:nk2_ext))
   
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,       &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QV(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QV3D_bg(:,:,2:nk2))               

   if (physics) then
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,       &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QC(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QC3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,       &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QI(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QI3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QR(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QR3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QS(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QS3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QG(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QG3D_bg(:,:,2:nk2))
   endif

   DO nt = 1,ntracer
   CALL GCM2Q3D_X (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QT(:,:,2:nk2,nt),                       &
                   channel%seg(num_seg)%QT3D_bg(:,:,2:nk2,nt))
   ENDDO

!  U3DX_bg contravariant component ---------------------------------  
   !  at u-point  (k=2:nk2): 
   CALL GCM2Q3D_X (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Ucon(:,:,2:nk2),                         &
                   channel%seg(num_seg)%U3DX_bg(:,:,2:nk2))
   
   ! at v-point   (k=2:2): (temporary use of Z3DX0) 
   CALL GCM2Q3D_X (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Ucon(:,:,2:2),                           &
                   channel%seg(num_seg)%Z3DX0(:,:,2:2))     
                   
   ! at v-point   (k=3:nk2): (temporary use of Z3DX0) 
   CALL GCM2Q3D_X (.FALSE.,3,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Ucon(:,:,3:nk2),                         &
                   channel%seg(num_seg)%Z3DX0(:,:,3:nk2))                        
!-------------------------------------------------------------------
   
!  U3DY_bg contravariant component ---------------------------------
   ! at v-point   (k=2:nk2): used for eddy calculation  
   CALL GCM2Q3D_X (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Vcon(:,:,2:nk2),                         &
                   channel%seg(num_seg)%U3DY_bg(:,:,2:nk2))
     
   ! at u-point   (k=2:2): (temporary use of Z3DY0)    
   CALL GCM2Q3D_X (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Vcon(:,:,2:2),                           &
                   channel%seg(num_seg)%Z3DY0(:,:,2:2))      

   ! at u-point   (k=3:nk2): (temporary use of Z3DY0)    
   CALL GCM2Q3D_X (.FALSE.,2,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%Vcon(:,:,3:nk2),                         &
                   channel%seg(num_seg)%Z3DY0(:,:,3:nk2))      

!---------------------------------

  IF (W3D_bg_INTERPOL) THEN
   ! Vertical velocity: used for vorticity calculation
   CALL GCM2Q3D_X (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,          &
                   channel%seg(num_seg)%lcen,                                    &
                   vGCM_st_new(num_seg)%W(:,:,2:nk1),                            &
                   channel%seg(num_seg)%W3D_bg(:,:,2:nk1))
  ENDIF                

!-----------------------------------------------------------------------
!     1. Prepare background fields of Pressure, Exner ftn, Density
!-----------------------------------------------------------------------
      DO J = mjm_c, mjp_c
       DO I = mim_c, mip_c
        
        channel%seg(num_seg)%Pmid_bg(I,J,1) = channel%seg(num_seg)%Pint_bg(I,J,1) 
        DO K = 2, NK2_ext
         channel%seg(num_seg)%Pmid_bg(I,J,K) = &
           0.5*(channel%seg(num_seg)%Pint_bg(I,J,K)+channel%seg(num_seg)%Pint_bg(I,J,K-1))
        ENDDO
        DO K = 1, NK2_ext 
         channel%seg(num_seg)%PIint_bg(I,J,K) = &
                         (channel%seg(num_seg)%Pint_bg(I,J,K)/pzero)**cappa
        
         channel%seg(num_seg)%PImid_bg(I,J,K) = &
                         (channel%seg(num_seg)%Pmid_bg(I,J,K)/pzero)**cappa
        ENDDO        
        
        DO K = 2, NK2_ext
         channel%seg(num_seg)%RHO_bg(I,J,K) = channel%seg(num_seg)%Pmid_bg(I,J,K) &
        /(gas_const*channel%seg(num_seg)%TH3D_bg(I,J,K)*channel%seg(num_seg)%PImid_bg(I,J,K))
        ENDDO

        channel%seg(num_seg)%RHOint_bg(I,J,1) = channel%seg(num_seg)%RHO_bg(I,J,2)
        channel%seg(num_seg)%RHOint_bg(I,J,NK2_ext) = channel%seg(num_seg)%RHO_bg(I,J,NK2_ext)
        DO K = 2, NK1_ext
         channel%seg(num_seg)%RHOint_bg(I,J,K) = &
                   ((ZZ(K)-ZZ(K-1))*channel%seg(num_seg)%RHO_bg(I,J,K)    &
                   +(ZZ(K+1)-ZZ(K))*channel%seg(num_seg)%RHO_bg(I,J,K+1)) &
                   /(2.0_r8*(ZT(K+1)-ZT(K)))
        ENDDO        
                        
       ENDDO
      ENDDO 
      
!--------------------------------------------------------------------------------
!     2. Calculate the BG covariant components of wind 
!--------------------------------------------------------------------------------
      ! For UV_CAL at k=2
      DO K = 2, 2
       DO J = mjm, mjp
        DO I = mim, mip
          VMAP = CON2COV(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_CO_bg(I,J,k) = VMAP(1) 
          
          VMAP = CON2COV(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_CO_bg(I,J,k) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO       
      
      DO K = 3, NK2
       DO J = mjm, mjp
        DO I = 0, mi1+1
          VMAP = CON2COV(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_CO_bg(I,J,k) = VMAP(1) 
          
          VMAP = CON2COV(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_CO_bg(I,J,k) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO                  
      
      IF (.NOT.W3D_bg_INTERPOL) THEN
!--------------------------------------------------------------------------------
!     3. Calculate the background vertical velocity
!--------------------------------------------------------------------------------

      channel%seg(num_seg)%W3D_bg(:,:,1)   = 0.0_r8
      channel%seg(num_seg)%W3D_bg(:,:,nk2) = 0.0_r8

!     Obtain (rho*W3D_bg) based on mass continuity       
      DO J = mjm_a, mjp
       DO I = mim_c, mip_c
        DO K = 2, nk1
         channel%seg(num_seg)%W3D_bg(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K-1)             &
     - DZ*RHO(K)*((channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%U3DX_bg(I,J,K)          & 
                  -channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX_bg(I-1,J,K))/DX  &
                 +(channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%U3DY_bg(I,J,K)          &
                  -channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY_bg(I,J-1,K))/DY) &
                /(channel%seg(num_seg)%RG_T(I,J)*FNT(K))  
        ENDDO
       ENDDO
      ENDDO
       
!     Obtain (W3D_bg) 
      DO J = mjm_a, mjp
       DO I = mim_c, mip_c
        DO K = 2, NK1
         channel%seg(num_seg)%W3D_bg(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K)/RHOZ(K) 
        ENDDO
       ENDDO
      ENDDO

!     Filling data
      DO I = mim_c, mip_c
       DO K = 2, NK1
        channel%seg(num_seg)%W3D_bg(I,mjm,K) = channel%seg(num_seg)%W3D_bg(I,mjm_a,K)
       ENDDO
      ENDDO    

      ENDIF  ! .NOT.W3D_bg_INTERPOL
      
!--------------------------------------------------------------------------------
!     4.0 Calculate the rll components of wind  (For debugging only) 
!--------------------------------------------------------------------------------
      channel%seg(num_seg)%U3DX_ll_bg(:,:,1) = val_missing
      channel%seg(num_seg)%U3DY_ll_bg(:,:,1) = val_missing
      DO K = 2, NK2
       DO J = mjm, mjp 
        DO nn = 2,nhalo
         IPOM = 1 - nn
         IPOP = mi1 + nn
         channel%seg(num_seg)%U3DX_ll_bg(IPOM,J,K) = val_missing
         channel%seg(num_seg)%U3DX_ll_bg(IPOP,J,K) = val_missing 
         channel%seg(num_seg)%U3DY_ll_bg(IPOM,J,K) = val_missing
         channel%seg(num_seg)%U3DY_ll_bg(IPOP,J,K) = val_missing
        ENDDO        
       ENDDO
      ENDDO          
      
      DO K = 2, NK2
       DO J = mjm, mjp
        DO I = 0, mi1+1
          VMAP = CON2RLL(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_ll_bg(I,J,K) = VMAP(1) 
          
          VMAP = CON2RLL(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_ll_bg(I,J,K) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO              

!--------------------------------------------------------------------------------
!    4. Calculate the vorticity components from the covariant components of wind
!--------------------------------------------------------------------------------

!----------------------------
!     Z3DX_bg calculation
!----------------------------

      DO K = 1, NK2
       DO J = mjm, mjp 
        DO nn = 1,nhalo
         IPOM = 1 - nn
         IPOP = mi1 + nn
         channel%seg(num_seg)%Z3DX_bg(IPOM,J,K) = val_missing
         channel%seg(num_seg)%Z3DX_bg(IPOP,J,K) = val_missing 
        ENDDO        
       ENDDO
      ENDDO          
      
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
      DO K = 1, NK2
       DO J = mjm, mjp 
        DO nn = 1,nhalo
         IPOM = 1 - nn
         IPOP = mi1 + nn
         channel%seg(num_seg)%Z3DY_bg(IPOM,J,K) = val_missing
         channel%seg(num_seg)%Z3DY_bg(IPOP,J,K) = val_missing 
        ENDDO        
       ENDDO
      ENDDO          

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
!     Z3DZ_bg calculation 
!----------------------------------------
! (used for z3dz prediction)
      K = NK2
      DO J = mjm, mjp 
       DO nn = 1,nhalo
         IPOM = 1 - nn
         IPOP = mi1 + nn
         channel%seg(num_seg)%Z3DZ_bg(IPOM,J,1) = val_missing
         channel%seg(num_seg)%Z3DZ_bg(IPOP,J,1) = val_missing        
       ENDDO
      ENDDO  
      
      DO J = mjm, mjp_a
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,K)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,K)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))/DY
       ENDDO
      ENDDO
      DO J = mjm, mjp_a
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,1)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
      ENDDO

      ! filling the halo data
      DO I = 1, mi1
       channel%seg(num_seg)%Z3DZ_bg(I,mjp,1) = channel%seg(num_seg)%Z3DZ_bg(I,mjp_a,1)
      ENDDO

! (used for RELAX_2D)   
      K = 2       
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,2) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,K)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,K)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))/DY
       ENDDO
      ENDDO
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,2) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,2)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
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

   CALL GCM2Q3D_Y (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim_c,mip_c,mjm_c,mjp_c,nk2_ext, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,              &
                   channel%seg(num_seg)%lcen,                                        &
                   vGCM_st_new(num_seg)%Pint(:,:,1:nk2_ext),                         &
                   channel%seg(num_seg)%Pint_bg(:,:,1:nk2_ext))
                   
   CALL GCM2Q3D_Y (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk2_ext-1, &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,              &
                   channel%seg(num_seg)%lcen,                                        &
                   vGCM_st_new(num_seg)%TH(:,:,2:nk2_ext),                           &
                   channel%seg(num_seg)%TH3D_bg(:,:,2:nk2_ext))
                   
   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,      &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QV(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QV3D_bg(:,:,2:nk2))
   
   if (physics) then
   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,       &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QC(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QC3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,       &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,             &
                   channel%seg(num_seg)%lcen,                                       &
                   vGCM_st_new(num_seg)%QI(:,:,2:nk2),                              &
                   channel%seg(num_seg)%QI3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QR(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QR3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QS(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QS3D_bg(:,:,2:nk2))

   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QG(:,:,2:nk2),                          &
                   channel%seg(num_seg)%QG3D_bg(:,:,2:nk2))
   endif

   DO NT = 1,ntracer
   CALL GCM2Q3D_Y (.TRUE.,1,0,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,   &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,         &
                   channel%seg(num_seg)%lcen,                                   &
                   vGCM_st_new(num_seg)%QT(:,:,2:nk2,nt),                       &
                   channel%seg(num_seg)%QT3D_bg(:,:,2:nk2,nt))
   ENDDO

! U3DX_bg contravariant component ---------------------------------  
   !  at u-point  (k=2:nk2): 
   CALL GCM2Q3D_Y (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Ucon(:,:,2:nk2),                          &
                   channel%seg(num_seg)%U3DX_bg(:,:,2:nk2))   
                   
   !  at v-point  (k= 2:2): (temporary use of Z3DX0)   
   CALL GCM2Q3D_Y (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Ucon(:,:,2:2),                            &
                   channel%seg(num_seg)%Z3DX0(:,:,2:2))

   !  at v-point  (k= 3:nk2): (temporary use of Z3DX0)   
   CALL GCM2Q3D_Y (.FALSE.,3,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Ucon(:,:,3:nk2),                          &
                   channel%seg(num_seg)%Z3DX0(:,:,3:nk2))

!-------------------------------------------------------------------

!  U3DY_bg contravariant component ---------------------------------
   ! at v-point   (k= 2:nk2):                   
   CALL GCM2Q3D_Y (.FALSE.,3,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Vcon(:,:,2:nk2),                          &
                   channel%seg(num_seg)%U3DY_bg(:,:,2:nk2))       
                      
   ! at u-point   (k= 2:2): (temporary use of Z3DY0)                     
   CALL GCM2Q3D_Y (.FALSE.,2,nhalo,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Vcon(:,:,2:2),                            &
                   channel%seg(num_seg)%Z3DY0(:,:,2:2))
                   
   ! at u-point   (k= 3:nk2): (temporary use of Z3DY0)                     
   CALL GCM2Q3D_Y (.FALSE.,2,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%Vcon(:,:,3:nk2),                          &
                   channel%seg(num_seg)%Z3DY0(:,:,3:nk2))                   
!-------------------------------------------------------------------                   
 
  IF (W3D_bg_INTERPOL) THEN 
   ! Vertical velocity: used for vorticity calculation
   CALL GCM2Q3D_Y (.FALSE.,1,1,channel%seg(num_seg)%NFACE,mim,mip,mjm,mjp,nk1-1,  &
                   channel%seg(num_seg)%lsta,channel%seg(num_seg)%lend,           &
                   channel%seg(num_seg)%lcen,                                     &
                   vGCM_st_new(num_seg)%W(:,:,2:nk1),                             &
                   channel%seg(num_seg)%W3D_bg(:,:,2:nk1))
  ENDIF
  
!--------------------------------------------------------------------
!     1. Prepare background fields of Pressure, Exner ftn, Density
!--------------------------------------------------------------------
      DO J = mjm_c, mjp_c
       DO I = mim_c, mip_c
       
        channel%seg(num_seg)%Pmid_bg(I,J,1) = channel%seg(num_seg)%Pint_bg(I,J,1) 
        DO K = 2, NK2_ext
         channel%seg(num_seg)%Pmid_bg(I,J,K) = &
           0.5*(channel%seg(num_seg)%Pint_bg(I,J,K)+channel%seg(num_seg)%Pint_bg(I,J,K-1))
        ENDDO
        DO K = 1, NK2_ext 
         channel%seg(num_seg)%PIint_bg(I,J,K) = &
                         (channel%seg(num_seg)%Pint_bg(I,J,K)/pzero)**cappa
        
         channel%seg(num_seg)%PImid_bg(I,J,K) = &
                         (channel%seg(num_seg)%Pmid_bg(I,J,K)/pzero)**cappa
        ENDDO        

        DO K = 2, NK2_ext
         channel%seg(num_seg)%RHO_bg(I,J,K) = channel%seg(num_seg)%Pmid_bg(I,J,K) &
        /(gas_const*channel%seg(num_seg)%TH3D_bg(I,J,K)*channel%seg(num_seg)%PImid_bg(I,J,K))
        ENDDO
        
        channel%seg(num_seg)%RHOint_bg(I,J,1) = channel%seg(num_seg)%RHO_bg(I,J,2)
        channel%seg(num_seg)%RHOint_bg(I,J,NK2_ext) = channel%seg(num_seg)%RHO_bg(I,J,NK2_ext)
        DO K = 2, NK1_ext
         channel%seg(num_seg)%RHOint_bg(I,J,K) = &
                   ((ZZ(K)-ZZ(K-1))*channel%seg(num_seg)%RHO_bg(I,J,K)    &
                   +(ZZ(K+1)-ZZ(K))*channel%seg(num_seg)%RHO_bg(I,J,K+1)) &
                   /(2.0_r8*(ZT(K+1)-ZT(K)))
        ENDDO        
                        
       ENDDO
      ENDDO     
 
!--------------------------------------------------------------------------------
!     2. Calculate the BG covariant components of wind
!--------------------------------------------------------------------------------
      ! For UV_CAL at k=2
      DO K = 2, 2
       DO J = mjm, mjp
        DO I = mim, mip
         
          VMAP = CON2COV(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_CO_bg(I,J,K) = VMAP(1) 
          
          VMAP = CON2COV(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_CO_bg(I,J,K) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO        
      
      DO K = 3, NK2
       DO J = 0, mj1+1
        DO I = mim, mip
         
          VMAP = CON2COV(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_CO_bg(I,J,K) = VMAP(1) 
          
          VMAP = CON2COV(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_CO_bg(I,J,K) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO                    
       
      IF (.NOT.W3D_bg_INTERPOL) THEN 
!--------------------------------------------------------------------------------
!     3. Calculate the background vertical velocity
!--------------------------------------------------------------------------------

      channel%seg(num_seg)%W3D_bg(:,:,1)   = 0.0_r8
      channel%seg(num_seg)%W3D_bg(:,:,nk2) = 0.0_r8

!     Obtain (rho*W3D_bg) based on mass continuity       
      DO J = mjm_c, mjp_c
       DO I = mim_a, mip
        DO K = 2, NK1
         channel%seg(num_seg)%W3D_bg(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K-1)             &
     - DZ*RHO(K)*((channel%seg(num_seg)%RG_U(I,J)*channel%seg(num_seg)%U3DX_bg(I,J,K)          & 
                  -channel%seg(num_seg)%RG_U(I-1,J)*channel%seg(num_seg)%U3DX_bg(I-1,J,K))/DX  &
                 +(channel%seg(num_seg)%RG_V(I,J)*channel%seg(num_seg)%U3DY_bg(I,J,K)          &
                  -channel%seg(num_seg)%RG_V(I,J-1)*channel%seg(num_seg)%U3DY_bg(I,J-1,K))/DY) &
                /(channel%seg(num_seg)%RG_T(I,J)*FNT(K))  
        ENDDO
       ENDDO
      ENDDO
       
!     Obtain (W3D_bg)
      DO J = mjm_c, mjp_c
       DO I = mim_a, mip
        DO K = 2, NK1
         channel%seg(num_seg)%W3D_bg(I,J,K) = channel%seg(num_seg)%W3D_bg(I,J,K)/RHOZ(K) 
        ENDDO
       ENDDO
      ENDDO

      DO J = mjm_c, mjp_c
       DO K = 2, NK1
        channel%seg(num_seg)%W3D_bg(mim,J,K) = channel%seg(num_seg)%W3D_bg(mim_a,J,K)
       ENDDO
      ENDDO

      ENDIF  ! .NOT.W3D_bg_INTERPOL
!--------------------------------------------------------------------------------
!     4.0 Calculate the rll components of wind (For debugging only)
!--------------------------------------------------------------------------------
      channel%seg(num_seg)%U3DX_ll_bg(:,:,1) = val_missing
      channel%seg(num_seg)%U3DY_ll_bg(:,:,1) = val_missing
      
      DO K = 2, NK2
       DO I = mim, mip 
        DO nn = 2,nhalo
         JPOM = 1 - nn
         JPOP = mj1 + nn
         channel%seg(num_seg)%U3DX_ll_bg(I,JPOM,K) = val_missing
         channel%seg(num_seg)%U3DX_ll_bg(I,JPOP,K) = val_missing 
         channel%seg(num_seg)%U3DY_ll_bg(I,JPOM,K) = val_missing
         channel%seg(num_seg)%U3DY_ll_bg(I,JPOP,K) = val_missing
        ENDDO        
       ENDDO
      ENDDO     
      
      DO K = 2, NK2
       DO J = 0, mj1+1
        DO I = mim, mip
         
          VMAP = CON2RLL(channel%seg(num_seg)%AM_U(1,i,j),    &
                         channel%seg(num_seg)%AM_U(2,i,j),    &
                         channel%seg(num_seg)%AM_U(3,i,j),    &
                         channel%seg(num_seg)%AM_U(4,i,j),    &
                         channel%seg(num_seg)%U3DX_bg(i,j,k), &
                         channel%seg(num_seg)%Z3DY0(i,j,k))
          channel%seg(num_seg)%U3DX_ll_bg(I,J,K) = VMAP(1) 
          
          VMAP = CON2RLL(channel%seg(num_seg)%AM_V(1,i,j),    &
                         channel%seg(num_seg)%AM_V(2,i,j),    &
                         channel%seg(num_seg)%AM_V(3,i,j),    &
                         channel%seg(num_seg)%AM_V(4,i,j),    &
                         channel%seg(num_seg)%Z3DX0(i,j,k),   &
                         channel%seg(num_seg)%U3DY_bg(i,j,k))
          channel%seg(num_seg)%U3DY_ll_bg(I,J,K) = VMAP(2)
        ENDDO
       ENDDO
      ENDDO              
                                                      
!---------------------------------------------------------------------------------
!     4. Calculate the vorticity components from the covariant components of wind
!---------------------------------------------------------------------------------
!----------------------------
!     Z3DX_bg calculation
!----------------------------
    
      DO K = 1, NK2
       DO I = mim, mip 
        DO nn = 1,nhalo
         JPOM = 1 - nn
         JPOP = mj1 + nn
         channel%seg(num_seg)%Z3DX_bg(I,JPOM,K) = val_missing
         channel%seg(num_seg)%Z3DX_bg(I,JPOP,K) = val_missing 
        ENDDO        
       ENDDO
      ENDDO     
      
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
      DO K = 1, NK2
       DO I = mim, mip 
        DO nn = 1,nhalo
         JPOM = 1 - nn
         JPOP = mj1 + nn
         channel%seg(num_seg)%Z3DY_bg(I,JPOM,K) = val_missing
         channel%seg(num_seg)%Z3DY_bg(I,JPOP,K) = val_missing 
        ENDDO        
       ENDDO
      ENDDO     

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
!     Z3DZ_bg calculation 
!---------------------------------------
! (used for prediction of z3dz)
      K = NK2
      DO I = mim, mip 
       DO nn = 1,nhalo
         JPOM = 1 - nn
         JPOP = mj1 + nn
         channel%seg(num_seg)%Z3DZ_bg(I,JPOM,1) = val_missing
         channel%seg(num_seg)%Z3DZ_bg(I,JPOP,1) = val_missing       
       ENDDO
      ENDDO     
      
      DO J = 1, mj1
       DO I = mim, mip_a
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,K)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,K)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))/DY
       ENDDO
      ENDDO
      DO J = 1, mj1
       DO I = mim, mip_a
        channel%seg(num_seg)%Z3DZ_bg(I,J,1) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,1)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
      ENDDO

      ! filling the halo data
      DO J = 1, mj1
       channel%seg(num_seg)%Z3DZ_bg(mip,J,1) = channel%seg(num_seg)%Z3DZ_bg(mip_a,J,1)
      ENDDO
      
! (Used for RELAX_2D)    
      K = 2
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,2) = &
                   (channel%seg(num_seg)%U3DY_CO_bg(I+1,J,K)     &
                   -channel%seg(num_seg)%U3DY_CO_bg(I,J,K))/DX   &
                 - (channel%seg(num_seg)%U3DX_CO_bg(I,J+1,K)     &
                   -channel%seg(num_seg)%U3DX_CO_bg(I,J,K))/DY
       ENDDO
      ENDDO
      DO J = 1, mj1
       DO I = 1, mi1
        channel%seg(num_seg)%Z3DZ_bg(I,J,2) = &
          channel%seg(num_seg)%Z3DZ_bg(I,J,2)/channel%seg(num_seg)%RG_Z(I,J)
       ENDDO
      ENDDO
  
!====================================================================
   ENDIF
!====================================================================

!********************************************************************
   ENDDO   ! num_seg
!********************************************************************

   DO num_seg = 1, 4
     deallocate(vGCM_st_new(num_seg)%Pint) 
     deallocate(vGCM_st_new(num_seg)%TH)
     deallocate(vGCM_st_new(num_seg)%QV)   
     deallocate(vGCM_st_new(num_seg)%QC)    
     deallocate(vGCM_st_new(num_seg)%QI)    
     deallocate(vGCM_st_new(num_seg)%QR)    
     deallocate(vGCM_st_new(num_seg)%QS)    
     deallocate(vGCM_st_new(num_seg)%QG)    
     deallocate(vGCM_st_new(num_seg)%QT)    
     deallocate(vGCM_st_new(num_seg)%Ucon)    
     deallocate(vGCM_st_new(num_seg)%Vcon)    
     deallocate(vGCM_st_new(num_seg)%W)
                   
     deallocate(vGCM_st_new(num_seg)%PI_G) 
     deallocate(vGCM_st_new(num_seg)%ZM_G)                                                                                                                                                                                                                                                                                 
   ENDDO
   
   END SUBROUTINE cal_bg

!===================================================================================
   SUBROUTINE VGCM_PREPARE (vGCM_st,vGCM_map,vGCM_stnew)   
!===================================================================================
!  Prepare vGCM state variables on vertical grids of VVM through interpolation
!  (Convert Temp. to theta at nLevel and interpolate the theta to CRM layers)  
!  Pay attention to 
!  the vertical size of Z3DZ_bg = 1 (=2  JUNG_UVLOW)
!  the vertical size of TH3D_bg = nk2_ext
!  the vertical size of others  = nk2

!  check unit of Pint [Pa] 
     
    type(vGCM_state_t), intent(in)  :: vGCM_st(4)      ! vGCM state
    type(vGCM_map_t),   intent(in)  :: vGCM_map(4)     ! vGCM map
    type(vGCM_vgrid), intent(inout) :: vGCM_stnew(4)   ! vGCM state
    
    ! Local
    REAL(kind=r8) :: ZM_mid(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)

    REAL(kind=r8) :: P_mid(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: Rho_mid(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: Pi_mid(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: Pi_int(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel+1)
    
    REAL(kind=r8) :: TH(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: W(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)

    REAL(kind=r8) :: Ucon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: Vcon(-nhalo_vGCM:nhalo_vGCM,1-nhalo_vGCM:nVGCM_seg+nhalo_vGCM,nLevel)
    REAL(kind=r8) :: VEC_CON(2)  
    
    INTEGER :: KVAL(nk2),num_seg,nt,I,J,K,Ibeg,Iend
        
!********************************************************************
    DO num_seg = 1, 4
!********************************************************************    
      DO J = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
      
      Ibeg = -nhalo_vGCM
      Iend = nhalo_vGCM      
      IF (J .EQ. (1-nhalo_vGCM) .OR. J .EQ. (nVGCM_seg+nhalo_vGCM)) THEN
        Ibeg = -nhalo_vGCM + 1
        Iend =  nhalo_vGCM - 1
      ENDIF  

      ! convert wind on RLL to contravariant components (on GCM grids) 
      DO I = Ibeg,Iend 
       DO K = 1, nLevel      
        VEC_CON = RLL2CON (vGCM_map(num_seg)%AMI(1,I,J),vGCM_map(num_seg)%AMI(2,I,J), &
                           vGCM_map(num_seg)%AMI(3,I,J),vGCM_map(num_seg)%AMI(4,I,J), &
                           vGCM_st(num_seg)%U(I,J,K),vGCM_st(num_seg)%V(I,J,K))
                           
        Ucon(I,J,K) = VEC_CON(1)
        Vcon(I,J,K) = VEC_CON(2) 
       ENDDO
      ENDDO      

      ! calculate Exner ftn & lnp at interfcae (on GCM grids) 
      DO I = Ibeg,Iend       
       DO K = 1, nLevel+1
         Pi_int(I,J,K) = (vGCM_st(num_seg)%Pint(I,J,K)/PZERO)**cappa  
       ENDDO
      ENDDO

      ! calculate p at mid-layer (on GCM grids)
      DO I = Ibeg,Iend       
       DO K = 1, nLevel
         P_mid(I,J,K) = 0.5_r8*(vGCM_st(num_seg)%Pint(I,J,K)+vGCM_st(num_seg)%Pint(I,J,K+1))
       ENDDO
      ENDDO  

      ! calculate Exner ftn & density at mid-layer (on GCM grids)                   
      DO I = Ibeg,Iend       
       DO K = 1, nLevel                          
         Pi_mid(I,J,K) = (P_mid(I,J,K)/PZERO)**cappa
         Rho_mid(I,J,K) = P_mid(I,J,K)/(gas_const*vGCM_st(num_seg)%T(I,J,K))
       ENDDO
      ENDDO
      
      ! calculate geopotential height at mid-layer (on GCM grids)                
      DO I = Ibeg,Iend       
       DO K = 1, nLevel                          
         ZM_mid(I,J,K) = 0.5_r8*(vGCM_st(num_seg)%zm_int(I,J,K)+vGCM_st(num_seg)%zm_int(I,J,K+1))
       ENDDO
      ENDDO  

      ! calculate potential (temp. & w) from (T & omega)
      DO I = Ibeg,Iend       
       DO K = 1, nLevel  
         TH(I,J,K) =  vGCM_st(num_seg)%T(I,J,K)/Pi_mid(I,J,K)    
         W(I,J,K)  = - vGCM_st(num_seg)%omega(I,J,K)/(grav*Rho_mid(I,J,K))
       ENDDO
      ENDDO 

      ENDDO  ! j-loop
      
!     Save Pi_G for the use in converting TH to/from T in other places
      DO K = 1, nLevel
       DO J = 1,nVGCM_seg
        vGCM_stnew(num_seg)%PI_G(K,J) = Pi_mid(0,J,K)
        vGCM_stnew(num_seg)%ZM_G(K,J) = ZM_mid(0,J,K)  
       ENDDO
      ENDDO    
      
!---------------------------------------------------------       
      ! vertical interpolation (column by column)      
!---------------------------------------------------------      
      DO J = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM

      Ibeg = -nhalo_vGCM
      Iend = nhalo_vGCM      
      IF (J .EQ. (1-nhalo_vGCM) .OR. J .EQ. (nVGCM_seg+nhalo_vGCM)) THEN
        Ibeg = -nhalo_vGCM + 1
        Iend =  nhalo_vGCM - 1
      ENDIF  
      
      DO I = Ibeg,Iend
 
      CALL vGCM_i_to_VVM_i (nk2_ext,nLevel,Pi_int(I,J,:),                               &
                            vGCM_st(num_seg)%zm_int(I,J,:),zz,vGCM_stnew(num_seg)%Pint(I,J,:))
                                                                              
      CALL vGCM_l_to_VVM_l (nk2_ext,nLevel,TH(I,J,:),                                    &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2_ext),  &
                            vGCM_stnew(num_seg)%TH(I,J,:))
                 
      ! Use a cubic-spline interpolation & constrain the positive values.                                
      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QV(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%QV(I,J,:),positive=.true.)

      ! Use a linear interpolation for QC,QI,QR,QS,QG,QT to avoid a generation of negative value.
      CALL vIndex_z (nk2,nLevel,ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),Kval)   
      
      IF (PHYSICS) THEN
      
      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QC(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &                        
                            vGCM_stnew(num_seg)%QC(I,J,:),K1val=Kval)

      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QI(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%QI(I,J,:),K1val=Kval)  
                                      
      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QR(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%QR(I,J,:),K1val=Kval)  
      
      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QS(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%QS(I,J,:),K1val=Kval)  
                                      
      CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QG(I,J,:),                       &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%QG(I,J,:),K1val=Kval) 
                            
      ENDIF   ! PHYSICS                             
            
      DO nt = 1,ntracer
       CALL vGCM_l_to_VVM_l (nk2,nLevel,vGCM_st(num_seg)%QT(I,J,:,nt),                   &
                             ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),     &
                             vGCM_stnew(num_seg)%QT(I,J,:,nt),K1val=Kval)  
      ENDDO                                                                                                
                                                                                    
      CALL vGCM_l_to_VVM_l (nk2,nLevel,Ucon(I,J,:),                                      &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%Ucon(I,J,:))                                   
                               
      CALL vGCM_l_to_VVM_l (nk2,nLevel,Vcon(I,J,:),                                      &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zt(1:nk2),      &
                            vGCM_stnew(num_seg)%Vcon(I,J,:)) 
                                     
      CALL vGCM_l_to_VVM_i (nk2,nLevel,W(I,J,:),                                         &
                            ZM_mid(I,J,:),vGCM_st(num_seg)%zm_int(I,J,:),zz(1:nk2),      &
                            vGCM_stnew(num_seg)%W(I,J,:))                                       
           
      vGCM_stnew(num_seg)%W(I,J,1)   = 0.0_r8
      vGCM_stnew(num_seg)%W(I,J,nk2) = 0.0_r8   
      
      ! Convert: Exner ftn pi to p
      DO K = 1,nk2_ext
       vGCM_stnew(num_seg)%Pint(I,J,K) = PZERO*(vGCM_stnew(num_seg)%Pint(I,J,K)**(1.0_r8/cappa))                 
      ENDDO
             
      ENDDO  ! i-loop

      ENDDO  ! j-loop  

!********************************************************************
    ENDDO  !num_seg 
!********************************************************************
      
    CONTAINS

!===========================================================================================
      Subroutine vIndex (kdim,kmG,exnerl_G,exner_G,exnerl_C,Kval)
!===========================================================================================      
!     Find the vertical grid index for a linear interpolation.
 
!     exnerl_G, exner_G, exnerl: downward indexing (top to bottom) 
!     exnerl_C: upward indexing (bottom to top)

      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: exnerl_G   ! exner ftn of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: exner_G    ! exner ftn of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: exnerl_C   ! exner ftn of VVM (at mid-layer)
      
      INTEGER, DIMENSION(kdim), INTENT(OUT)  :: kval     

      ! LOCAL
      LOGICAL LF
      REAL(kind=r8), DIMENSION(0:kmG+1) :: exnerl
      INTEGER k
      
      do k=1,kmG
       exnerl(k) = exnerl_G(k)
      enddo
      exnerl(0)     = exner_G(1)
      exnerl(kmG+1) = exner_G(kmG+1)
      
      do k=2,kdim 
       kval(k) = INDEXR (exnerl_C(k),kmG+2,exnerl,LF) - 1 
      enddo      
      
      END SUBROUTINE vindex
      
!===========================================================================================
      Subroutine vIndex_z (kdim,kmG,hgtl_G,hgt_G,hgtl_C,Kval)
!===========================================================================================      
!     Find the vertical grid index for a linear interpolation.
 
!     hgtl_G, hgt_G: downward indexing (top to bottom) 
!     hgtl_C: upward indexing (bottom to top)

      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: hgtl_G   ! height of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: hgt_G    ! height of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgtl_C   ! height of VVM (at mid-layer)
      
      INTEGER, DIMENSION(kdim), INTENT(OUT)  :: kval     

      ! LOCAL
      LOGICAL LF
      REAL(kind=r8), DIMENSION(0:kmG+1) :: hgtl
      INTEGER k
      
      do k=1,kmG
       hgtl(k) = hgtl_G(kmG-k+1)
      enddo
      
      hgtl(0)     = hgt_G(kmG+1)
      if (hgt_G(1) .ge. hgtl_C(kdim)) then
        hgtl(kmG+1) = hgt_G(1) 
      else
        hgtl(kmG+1) = 1.001_r8*hgtl_C(kdim)
      endif         
      
      do k=2,kdim 
       kval(k) = INDEXR (hgtl_C(k),kmG+2,hgtl,LF) - 1 
      enddo      
      
      END SUBROUTINE vindex_z      
         
!===========================================================================================
      Subroutine vGCM_layers_to_VVM_layers (kdim,kmG,val_G,exnerl_G,exner_G, &
                                            exnerl_C,val_C,k1val)
!===========================================================================================      
!     Vertical interpolation of GCM values to CRM grids (layer to layer)
!     Val_G, exnerl, val: downward indexing (top to bottom) 
!     Val_C: upward indexing (bottom to top)
      
      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: val_G      ! target variable of GCM (at mid-layer)
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: exnerl_G   ! exner ftn of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: exner_G    ! exner ftn of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: exnerl_C   ! exner ftn of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C      ! target variable of VVM (at mid-layer)
      
      INTEGER, DIMENSION(kdim), OPTIONAL, INTENT(IN)  :: k1val
      
      ! LOCAL
      REAL(kind=r8), DIMENSION(0:kmG+1) :: exnerl,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k,k1,k2
   
      do k=1,kmG
       exnerl(k) = exnerl_G(k)
       val(k) = val_G(k)
      enddo
       
      ! Prepare boundary values 
      exnerl(0)     = exner_G(1)
      exnerl(kmG+1) = exner_G(kmG+1)
      
      val(0)     = val(1)   + 0.6_r8*(exnerl(0)-exnerl(1))  &
                                    *(val(1)-val(2))/(exnerl(1)-exnerl(2))                    
                             
      val(kmG+1) = val(kmG) + 0.6_r8*(exnerl(kmG+1)-exnerl(kmG))  &
                                    *(val(kmG)-val(kmG-1))/(exnerl(kmG)-exnerl(kmG-1))

      IF (PRESENT(k1val)) THEN
      ! linear interpolation
      
        DO k=2,kdim 
         K1 = k1val(k)
         K2 = K1 + 1
         val_C(k) = FINTRP(1,exnerl_C(k),exnerl(K1),val(K1),exnerl(K2),val(K2))
        ENDDO     

      ELSE                              
      ! cubic-spline interpolation
        yp1  =  0.5_r8*((val(1)-val(0))/(exnerl(1)-exnerl(0))  &
                       +(val(2)-val(1))/(exnerl(2)-exnerl(1)))

        ypkm =  0.5_r8*((val(kmG+1)-val(kmG))/(exnerl(kmG+1)-exnerl(kmG))   & 
                       +(val(kmG)-val(kmG-1))/(exnerl(kmG)-exnerl(kmG-1)))
      
        call spline(exnerl,val,kmG+2,yp1,ypkm,y2)

        do k=2,kdim 
         call splint(exnerl,val,y2,kmG+2,exnerl_C(k),val_C(k))
        enddo
        
      ENDIF

      END SUBROUTINE vGCM_layers_to_VVM_layers

!===========================================================================================
      Subroutine vGCM_layers_to_VVM_interfaces (kdim,kmG,val_G,exnerl_G,exner_G,exner_C,val_C)
!===========================================================================================  
!     This routine is meant to be used for w only.    
!     Vertical interpolation of GCM values to CRM grids (layer to interface)
!     Val_G, exnerl, val: downward indexing (top to bottom) 
!     Val_C: upward indexing (bottom to top)
      
      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: val_G      ! target variable of GCM (at mid-layer)
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: exnerl_G   ! exner ftn of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: exner_G    ! exner ftn of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: exner_C    ! exner ftn of VVM (at interface)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C      ! target variable of VVM (at interface)

      ! LOCAL
      REAL(kind=r8), DIMENSION(0:kmG+1) :: exnerl,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k
             
      do k=1,kmG
       exnerl(k) = exnerl_G(k)
       val(k) = val_G(k)
      enddo
       
      ! Prepare boundary values 
      exnerl(0)     = exner_G(1)
      exnerl(kmG+1) = exner_G(kmG+1)
      
      val(0)     = val(1)   + 0.6_r8*(exnerl(0)-exnerl(1))        &
                                    *(val(1)-val(2))/(exnerl(1)-exnerl(2))
                             
      val(kmG+1) = val(kmG) + 0.6_r8*(exnerl(kmG+1)-exnerl(kmG))  &
                                    *(val(kmG)-val(kmG-1))/(exnerl(kmG)-exnerl(kmG-1))

      yp1  =  0.5_r8*((val(1)-val(0))/(exnerl(1)-exnerl(0))        &
                     +(val(2)-val(1))/(exnerl(2)-exnerl(1)))

      ypkm =  0.5_r8*((val(kmG+1)-val(kmG))/(exnerl(kmG+1)-exnerl(kmG))   & 
                     +(val(kmG)-val(kmG-1))/(exnerl(kmG)-exnerl(kmG-1)))
      
      call spline(exnerl,val,kmG+2,yp1,ypkm,y2)
       
      do k=2,kdim-1 
       call splint(exnerl,val,y2,kmG+2,exner_C(k),val_C(k))
      enddo
       
      END SUBROUTINE vGCM_layers_to_VVM_interfaces    
             
!===========================================================================================
      SUBROUTINE vGCM_interfaces_to_VVM_interfaces (kdim,kmG,val_G,exner_G,exner_C,val_C)
!===========================================================================================
!     Vertical interpolation of GCM values to CRM grids (interface to interface)
!     Val_G, exner, val: downward indexing (top to bottom)  
!     Val_C: upward indexing (bottom to top)

      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: val_G     ! target variable of GCM (at interface)
      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: exner_G   ! exner ftn of GCM (at interface)

      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: exner_C   ! exner ftn of VVM (at interface)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C     ! target variable of VVM (at interface) 

      ! Local
      REAL(kind=r8), DIMENSION(0:kmG+2) :: exner,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k
      
      do k=1,kmG+1
       exner(k) = exner_G(k)
       val(k) = val_G(k)
      enddo

      exner(0)     = 0.9_r8*exner_G(1)
      exner(kmG+2) = 1.001_r8*exner_G(kmG+1)
      
      val(0)     = val(1)
      val(kmG+2) = val(kmG+1)

      yp1  = 0.5_r8*((val(1)-val(0))/(exner(1)-exner(0))    &
                    +(val(2)-val(1))/(exner(2)-exner(1)))

      ypkm = 0.5_r8*((val(kmG+2)-val(kmG+1))/(exner(kmG+2)-exner(kmG+1))   & 
                    +(val(kmG+1)-val(kmG))/(exner(kmG+1)-exner(kmG)))

      call spline(exner,val,kmG+3,yp1,ypkm,y2)

      do k=1,kdim 
       call splint(exner,val,y2,kmG+3,exner_C(k),val_C(k))
      enddo

      END SUBROUTINE vGCM_interfaces_to_VVM_interfaces

!===========================================================================================
      Subroutine vGCM_l_to_VVM_l (kdim,kmG,val_G,hgtl_G,hgt_G, &
                                  hgtl_C,val_C,k1val,positive)
!===========================================================================================      
!     Vertical interpolation of GCM values to CRM grids (layer to layer)
!     Val_G, hgtl_G, val_G: downward indexing (top to bottom) 
!     Val_C: upward indexing (bottom to top)
      
      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: val_G    ! target variable of GCM (at mid-layer)
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: hgtl_G   ! height of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: hgt_G    ! height of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgtl_C   ! height of VVM (at mid-layer)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C    ! target variable of VVM (at mid-layer)
      
      INTEGER, DIMENSION(kdim), OPTIONAL, INTENT(IN)  :: k1val
      LOGICAL, OPTIONAL, INTENT(IN)  :: positive
      
      ! LOCAL
      REAL(kind=r8), parameter :: fac = 1.0_r8
      REAL(kind=r8), DIMENSION(0:kmG+1) :: hgtl,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k,k1,k2
   
      ! Reverse the array
      do k=1,kmG
       hgtl(k) = hgtl_G(kmG-k+1)
       val(k)  = val_G(kmG-k+1)
      enddo

!-------------------------------------------------------------------------------
!     Prepare boundary values
!-------------------------------------------------------------------------------      
      hgtl(0) = hgt_G(kmG+1)
      val(0)  = val(1) + (hgtl(0)-hgtl(1))*fac*(val(1)-val(2))/(hgtl(1)-hgtl(2))
      
      if (hgt_G(1) .ge. hgtl_C(kdim)) then
        hgtl(kmG+1) = hgt_G(1) 
      else
        hgtl(kmG+1) = 1.001_r8*hgtl_C(kdim)
      endif
      
      ypkm =  0.5_r8*((val(kmG)-val(kmG-1))/(hgtl(kmG)-hgtl(kmG-1))   & 
                     +(val(kmG-1)-val(kmG-2))/(hgtl(kmG-1)-hgtl(kmG-2)))
                            
      val(kmG+1) = val(kmG) + (hgtl(kmG+1)-hgtl(kmG))*fac*ypkm
!-------------------------------------------------------------------------------  
                                           
      IF (PRESENT(k1val)) THEN
      ! linear interpolation
            
        DO k=2,kdim         
         K1 = k1val(k)
         K2 = K1 + 1
         val_C(k) = FINTRP(1,hgtl_C(k),hgtl(K1),val(K1),hgtl(K2),val(K2))
        ENDDO     

      ELSE                              
      ! cubic-spline interpolation
        yp1  =  0.5_r8*((val(1)-val(0))/(hgtl(1)-hgtl(0))  &
                       +(val(2)-val(1))/(hgtl(2)-hgtl(1)))

        ypkm =  0.5_r8*((val(kmG+1)-val(kmG))/(hgtl(kmG+1)-hgtl(kmG))   & 
                       +(val(kmG)-val(kmG-1))/(hgtl(kmG)-hgtl(kmG-1)))
      
        call spline(hgtl,val,kmG+2,yp1,ypkm,y2)

        do k=2,kdim 
         call splint(hgtl,val,y2,kmG+2,hgtl_C(k),val_C(k))
        enddo
        
      ENDIF
            
      if (present(positive)) then
        do k=2,kdim 
         val_C(k) = MAX(0.0_r8,val_C(k))
        enddo
      endif

      END SUBROUTINE vGCM_l_to_VVM_l

!===========================================================================================
      Subroutine vGCM_l_to_VVM_i (kdim,kmG,val_G,hgtl_G,hgt_G,hgt_C,val_C)
!===========================================================================================  
!     This routine is meant to be used for w only.    
!     Vertical interpolation of GCM values to CRM grids (layer to interface)
!     Val_G, hgtl_G, val_G: downward indexing (top to bottom) 
!     Val_C: upward indexing (bottom to top)
      
      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 
      
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: val_G    ! target variable of GCM (at mid-layer)
      REAL(kind=r8), DIMENSION(kmG),   INTENT(IN) :: hgtl_G   ! height of GCM (at mid-layer) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: hgt_G    ! height of GCM (at interface)
  
      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgt_C    ! height of VVM (at interface)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C    ! target variable of VVM (at interface)

      ! LOCAL
      REAL(kind=r8), parameter :: fac = 1.0_r8
      REAL(kind=r8), DIMENSION(0:kmG+1) :: hgtl,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k
         
      ! Reverse the array       
      do k=1,kmG
       hgtl(k) = hgtl_G(kmG-k+1)
       val(k) = val_G(kmG-k+1)
      enddo
       
!-------------------------------------------------------------------------------
!     Prepare boundary values
!-------------------------------------------------------------------------------      
      hgtl(0) = hgt_G(kmG+1)
      val(0)  = val(1) + (hgtl(0)-hgtl(1))*fac*(val(1)-val(2))/(hgtl(1)-hgtl(2))
      
      if (hgt_G(1) .ge. hgt_C(kdim-1)) then
        hgtl(kmG+1) = hgt_G(1) 
      else
        hgtl(kmG+1) = 1.001_r8*hgt_C(kdim-1)
      endif
      
      ypkm =  0.5_r8*((val(kmG)-val(kmG-1))/(hgtl(kmG)-hgtl(kmG-1))   & 
                     +(val(kmG-1)-val(kmG-2))/(hgtl(kmG-1)-hgtl(kmG-2)))
                            
      val(kmG+1) = val(kmG) + (hgtl(kmG+1)-hgtl(kmG))*fac*ypkm
!-------------------------------------------------------------------------------      

      yp1  =  0.5_r8*((val(1)-val(0))/(hgtl(1)-hgtl(0))               &
                     +(val(2)-val(1))/(hgtl(2)-hgtl(1)))

      ypkm =  0.5_r8*((val(kmG+1)-val(kmG))/(hgtl(kmG+1)-hgtl(kmG))   & 
                     +(val(kmG)-val(kmG-1))/(hgtl(kmG)-hgtl(kmG-1)))
      
      call spline(hgtl,val,kmG+2,yp1,ypkm,y2)
       
      do k=2,kdim-1 
       call splint(hgtl,val,y2,kmG+2,hgt_C(k),val_C(k))
      enddo
       
      END SUBROUTINE vGCM_l_to_VVM_i  
      
!===========================================================================================
      SUBROUTINE vGCM_i_to_VVM_i (kdim,kmG,val_G,hgt_G,hgt_C,val_C)
!===========================================================================================
!     Vertical interpolation of GCM values to CRM grids (interface to interface)
!     Val_G, exner_G, val_G: downward indexing (top to bottom)  
!     Val_C: upward indexing (bottom to top)

      INTEGER, INTENT(IN) :: kdim   ! vertical dimension of VAL_C (nk2 or nk2_total)
      INTEGER, INTENT(IN) :: kmG    ! vertical dimension of VAL_G (nLevel) 

      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: val_G   ! target variable of GCM (at interface)
      REAL(kind=r8), DIMENSION(kmG+1), INTENT(IN) :: hgt_G   ! height of GCM (at interface)

      REAL(kind=r8), DIMENSION(kdim), INTENT(IN)  :: hgt_C   ! height of VVM (at interface)
      REAL(kind=r8), DIMENSION(kdim), INTENT(OUT) :: val_C   ! target variable of VVM (at interface) 
      
      ! Local
      REAL(kind=r8), parameter :: fac = 1.0_r8   ! This routine is used only for pi interpolation
      
      REAL(kind=r8), DIMENSION(kmG+2) :: hgt,val,y2
      REAL(kind=r8) :: yp1,ypkm 
      INTEGER k
            
      do k=1,kmG+1
       hgt(k) = hgt_G(kmG-k+2)
       val(k) = val_G(kmG-k+2)
      enddo

!--------------------------------------------------------------------
!     Prepare upper boundary values
!--------------------------------------------------------------------       
      if (hgt(kmG+1) .ge. hgt_C(kdim)) then
        hgt(kmG+2) = 1.001_r8*hgt(kmG+1)
      else
        hgt(kmG+2) = 1.001_r8*hgt_C(kdim)
      endif
      
      ypkm = 0.5_r8*((val(kmG+1)-val(kmG))/(hgt(kmG+1)-hgt(kmG))       & 
                    +(val(kmG)-val(kmG-1))/(hgt(kmG)-hgt(kmG-1)))            
      
      val(kmG+2) = val(kmG+1) + (hgt(kmG+2)-hgt(kmG+1))*fac*ypkm 
!--------------------------------------------------------------------      
      
      yp1  = 0.5_r8*((val(2)-val(1))/(hgt(2)-hgt(1)))

      ypkm = 0.5_r8*((val(kmG+2)-val(kmG+1))/(hgt(kmG+2)-hgt(kmG+1))   & 
                    +(val(kmG+1)-val(kmG))/(hgt(kmG+1)-hgt(kmG)))

      call spline(hgt,val,kmG+2,yp1,ypkm,y2)

      do k=1,kdim 
       call splint(hgt,val,y2,kmG+2,hgt_C(k),val_C(k))
      enddo      
      
      END SUBROUTINE vGCM_i_to_VVM_i       
             
   END SUBROUTINE vgcm_prepare

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
   REAL(KIND=r8) :: XVAL,YVAL,X0,Y0,X_COR,Y_COR

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

!-------------------------------------------------------------    
!   Fill up the halos, for which bg values are not calculated. 
!-------------------------------------------------------------
    lsta0 = lsta(1) - mhalo_l
    lend0 = lend(nVGCM_seg) + mhalo_l
    
    DO K = 1, kdim
     DO chn = mjm, mjp
      DO chl = mim,lsta0-1
       A_bg(chl,chn,K) = val_missing
      ENDDO
      DO chl = lend0+1,mip
       A_bg(chl,chn,K) = val_missing
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

    X0 = FLOAT(netsz) - lcen(I) + X_COR
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
     XVAL = FLOAT(chl) + X0  
     DO chn = 1, mjp
       YVAL = FLOAT(chn-1) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)

       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))
     ENDDO
    ENDDO
    ENDDO

    Y0 = FLOAT(netsz) + Y_COR
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
     XVAL = FLOAT(chl) + X0 
     DO chn = mjm, 0
       YVAL = FLOAT(chn-1) + Y0
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)

       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))
     ENDDO
    ENDDO
    ENDDO

    X0 = - lcen(I) + X_COR
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
     XVAL = FLOAT(chl) + X0 
     DO chn = 1, mjp
       YVAL = FLOAT(chn-1) + Y_COR
       A_bg(chl,chn,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)

       IF (POSITIVE) A_bg(chl,chn,K) = MAX(0.0_r8,A_bg(chl,chn,K))
     ENDDO
    ENDDO
    ENDDO

    Y0 = FLOAT(netsz) + Y_COR
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
     XVAL = FLOAT(chl) + X0
     DO chn = mjm, 0
       YVAL = FLOAT(chn-1) + Y0
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
   INTEGER :: I,J,K,Irev

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
          Irev = nVGCM_seg - I + 1
                     
          A_out(I,J,K) = A_in(J,Irev,K)

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
   REAL(KIND=r8) :: XVAL,YVAL,X0,Y0,X_COR,Y_COR

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

!-------------------------------------------------------------    
!   Fill up the halos, for which bg values are not calculated. 
!-------------------------------------------------------------
    lsta0 = lsta(1) - mhalo_l
    lend0 = lend(nVGCM_seg) + mhalo_l
    
    DO K = 1, kdim
     DO chn = mim, mip
      DO chl = mjm,lsta0-1
       A_bg(chn,chl,K) = val_missing
      ENDDO
      DO chl = lend0+1,mjp
       A_bg(chn,chl,K) = val_missing
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

    Y0 = FLOAT(netsz) - lcen(J) + Y_COR
    X0 = FLOAT(netsz) + X_COR
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
     YVAL = FLOAT(chl) + Y0  
     DO chn = mim, 0
       XVAL = FLOAT(chn-1) + X0 
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
     YVAL = FLOAT(chl) + Y0 
     DO chn = 1, mip
       XVAL = FLOAT(chn-1) + X_COR
       A_bg(chn,chl,K) = QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                                 VAL0,VAL1,VAL2,VAL3,netsz,netsz_sq)

       IF (POSITIVE) A_bg(chn,chl,K) = MAX(0.0_r8,A_bg(chn,chl,K))
     ENDDO
    ENDDO
    ENDDO

    Y0 = - lcen(J) + Y_COR 
    X0 = FLOAT(netsz) + X_COR
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
     YVAL = FLOAT(chl) + Y0
     DO chn = mim, 0
       XVAL = FLOAT(chn-1) + X0
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
     YVAL = FLOAT(chl) + Y0
     DO chn = 1, mip
       XVAL = FLOAT(chn-1) + X_COR
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
   INTEGER :: I,J,K,Jrev

   IF (NFACE.EQ.2 .OR. NFACE.EQ.3) THEN
      DO K = 1, kdim
       DO J = 1-nhalo_vGCM,nVGCM_seg+nhalo_vGCM
        Jrev = nVGCM_seg - J + 1
        DO I = -nhalo_vGCM,nhalo_vGCM
        
          A_out(I,J,K) = A_in(-I,Jrev,K)

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
      
!===========================================================================================
   FUNCTION QUADINT(XVAL,YVAL,V01,V02,V03,V04,V05,V06,V07,V08, &
                    VAL0,VAL1,VAL2,VAL3,N,NSQ) RESULT (QUADINT_result)
!===========================================================================================                    

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

END MODULE q3d_bg_module
