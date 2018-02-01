module vvm_data_types

  USE shr_kind_mod, only: dbl_kind => shr_kind_r8
  USE parmsld,      only: nVGCM,nVGCM_seg,channel_seg_l, &
                          ntracer,nk1,nk2,nk3,nhalo 
  
  implicit none
  private

!***********************************************************************************************
! Define a type (1): Variables declared for a channel segment (covering 1/4 of a channel)
!***********************************************************************************************
  
  type, public :: channel_seg_t
    integer :: nface                  ! Number of face where the segment belongs
    integer :: nx_size,ny_size        ! Sizes of x- and y-domain of a channel segment
    integer :: mi1,mim,mim_c,mim_ce,mim_a,mip,mip_c,mip_ce,mip_a  
    integer :: mj1,mjm,mjm_c,mjm_ce,mjm_a,mjp,mjp_c,mjp_ce,mjp_a      
    
    integer :: lsta(nVGCM_seg)        ! starting CRM-point index of a vGCM cell (along a channel)
    integer :: lend(nVGCM_seg)        ! ending CRM-point index of a vGCM cell
    real (kind=dbl_kind) :: lcen(nVGCM_seg)  ! center CRM-point index of a vGCM cell (vGCM point)
!----------------- 
!   TIMEINFO
!-----------------     
    integer :: iyr                    ! Year of current time step
    integer :: imonth                 ! Month of current time step
    integer :: iday                   ! Calendar day of current time step
    
    real (kind=dbl_kind) :: rjday0    ! Fractional Julian day at the start of the model simulation
    real (kind=dbl_kind) :: rjday     ! Fractional Julian day of the current model time step
    real (kind=dbl_kind) :: utc_time  ! UTC time of the current model time step    

!-----------------  
!   BACKGROUND FIELDS
!-----------------     
    ! Counterparts of Prognostic 3-D thermodynamical variables     
    real (kind=dbl_kind), pointer :: TH3D_bg(:,:,:)   => NULL() ! potential temperature (K)
    real (kind=dbl_kind), pointer :: QV3D_bg(:,:,:)   => NULL() ! water vapor mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QC3D_bg(:,:,:)   => NULL() ! cloud water mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QI3D_bg(:,:,:)   => NULL() ! cloud ice mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QR3D_bg(:,:,:)   => NULL() ! rain mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QS3D_bg(:,:,:)   => NULL() ! snow mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QG3D_bg(:,:,:)   => NULL() ! graupel mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QT3D_bg(:,:,:,:) => NULL() ! tracer mixing ratio (kg/kg)  
       
    ! Counterparts of Prognostic 3-D dynamical variables
    real (kind=dbl_kind), pointer :: Z3DX_bg(:,:,:) => NULL() ! x-component of vorticity, dw/dy-dv/dz (1/s)
    real (kind=dbl_kind), pointer :: Z3DY_bg(:,:,:) => NULL() ! y-component of vorticity, du/dz-dw/dx (1/s)
    real (kind=dbl_kind), pointer :: Z3DZ_bg(:,:,:) => NULL() ! z-component of vorticity, dv/dx-du/dy (1/s)
    
    ! Counterparts of Diagnostic 3-D dynamical variables
    real (kind=dbl_kind), pointer :: U3DX_bg(:,:,:)    => NULL() ! (contravariant) zonal velocity, u (m/s)
    real (kind=dbl_kind), pointer :: U3DY_bg(:,:,:)    => NULL() ! (contravariant) meridional velocity, v (m/s)
    real (kind=dbl_kind), pointer :: U3DX_co_bg(:,:,:) => NULL() ! (covariant) zonal velocity, u (m/s)
    real (kind=dbl_kind), pointer :: U3DY_co_bg(:,:,:) => NULL() ! (covariant) meridional velocity, v (m/s)
    real (kind=dbl_kind), pointer :: W3D_bg(:,:,:)     => NULL() ! vertical velocity, w (m/s)
    
!-----------------  
!   CONST3D
!-----------------     
    ! Prognostic 3-D thermodynamical variables     
    real (kind=dbl_kind), pointer :: TH3D(:,:,:)   => NULL() ! potential temperature (K)
    real (kind=dbl_kind), pointer :: QV3D(:,:,:)   => NULL() ! water vapor mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QC3D(:,:,:)   => NULL() ! cloud water mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QI3D(:,:,:)   => NULL() ! cloud ice mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QR3D(:,:,:)   => NULL() ! rain mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QS3D(:,:,:)   => NULL() ! snow mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QG3D(:,:,:)   => NULL() ! graupel mixing ratio (kg/kg)
    real (kind=dbl_kind), pointer :: QT3D(:,:,:,:) => NULL() ! tracer mixing ratio (kg/kg)  
       
    ! Prognostic 3-D dynamical variables
    real (kind=dbl_kind), pointer :: Z3DX(:,:,:) => NULL() ! x-component of vorticity, dw/dy-dv/dz (1/s)
    real (kind=dbl_kind), pointer :: Z3DY(:,:,:) => NULL() ! y-component of vorticity, du/dz-dw/dx (1/s)
    real (kind=dbl_kind), pointer :: Z3DZ(:,:,:) => NULL() ! z-component of vorticity, dv/dx-du/dy (1/s)
    
    ! Diagnostic 3-D dynamical variables
    real (kind=dbl_kind), pointer :: U3DX(:,:,:)    => NULL() ! (contravariant) zonal velocity, u (m/s)
    real (kind=dbl_kind), pointer :: U3DY(:,:,:)    => NULL() ! (contravariant) meridional velocity, v (m/s)
    real (kind=dbl_kind), pointer :: U3DX_CO(:,:,:) => NULL() ! (covariant) zonal velocity, u (m/s)
    real (kind=dbl_kind), pointer :: U3DY_CO(:,:,:) => NULL() ! (covariant) meridional velocity, v (m/s)
    real (kind=dbl_kind), pointer :: W3D(:,:,:)     => NULL() ! vertical velocity, w (m/s)
    real (kind=dbl_kind), pointer :: W3DNM1(:,:,:)  => NULL() ! vertical velocity, w (m/s), in the previous timestep
    
    ! Diagnostic 2-D dynamical variables (uppermost layer)
    real (kind=dbl_kind), pointer :: PSI(:,:)    => NULL() ! stream function (m**2/s)
    real (kind=dbl_kind), pointer :: PSINM1(:,:) => NULL() ! stream function (m**2/s), in the previous timestep
    real (kind=dbl_kind), pointer :: CHI(:,:)    => NULL() ! velocity potential (m**2/s)
    real (kind=dbl_kind), pointer :: CHINM1(:,:) => NULL() ! velocity potential (m**2/s), in the previous timestep
    
    ! Diagnostic 2-D variables (top or surface) : Add more diagnostics (or remove some)
    real (kind=dbl_kind), pointer :: TG(:,:)     => NULL() ! ground temperature (K)
    real (kind=dbl_kind), pointer :: ZROUGH(:,:) => NULL() ! roughness length (m)
    real (kind=dbl_kind), pointer :: GWET(:,:)   => NULL() ! ground wetness
    
    real (kind=dbl_kind), pointer :: UW(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2
    real (kind=dbl_kind), pointer :: WV(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2
    
    real (kind=dbl_kind), pointer :: WTH(:,:) => NULL() ! surface heat (potential temp.) flux (kg/m**3)(K)(m/s)
    real (kind=dbl_kind), pointer :: WQV(:,:) => NULL() ! surface moisture flux (kg/m**3)(kg/kg)(m/s)
    
    real (kind=dbl_kind), pointer :: SPREC(:,:) => NULL() ! surface precipitation rate (kg/m**2/s): sprec*3600. (mm/hr)

!-----------------  
!   CONST3D_Extra
!-----------------  
    ! alternative space (used in halo_vort & relax_3d): multi-use
    real (kind=dbl_kind), pointer :: Z3DX0(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: Z3DY0(:,:,:) => NULL()  
    real (kind=dbl_kind), pointer :: TERM1(:,:,:) => NULL()     
    
    ! deformation & shear components (used in turbulence and nonlinear diffusion)  
    real (kind=dbl_kind), pointer :: DEFYZ(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: DEFXZ(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: DEFXY(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: DEFXX(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: DEFYY(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: DEFZZ(:,:,:) => NULL() 

    ! coefficients (used in turbulence and diffusion: initially calculated)  
    real (kind=dbl_kind), pointer :: COEFX(:,:) => NULL()    
    real (kind=dbl_kind), pointer :: COEFY(:,:) => NULL()    
    real (kind=dbl_kind), pointer :: COEFA(:,:) => NULL()       
    real (kind=dbl_kind), pointer :: COEFX_K(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFY_K(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFA_K(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFX_E(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFY_E(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFA_E(:,:) => NULL() 
    real (kind=dbl_kind), pointer :: COEFX_Z(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFY_Z(:,:) => NULL()  
    real (kind=dbl_kind), pointer :: COEFA_Z(:,:) => NULL()     

    ! coefficients (used in relaxation_3d: initially calculated)
    real (kind=dbl_kind), pointer :: AGAU_CO(:,:,:) => NULL()
    real (kind=dbl_kind), pointer :: BGAU_CO(:,:,:) => NULL()
    real (kind=dbl_kind), pointer :: CGAU_CO(:,:,:) => NULL()
    
    ! coefficients (used in relaxation_2d: initially calculated)
    real (kind=dbl_kind), pointer :: RX2D_C1(:,:) => NULL()
    real (kind=dbl_kind), pointer :: RX2D_C2(:,:) => NULL()
    
!-----------------  
!   CONST3D_TD    
!----------------- 
    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (advection).    
    real (kind=dbl_kind), pointer :: FTH3D(:,:,:,:)   => NULL() ! (K/s) & random perturbation 
    real (kind=dbl_kind), pointer :: FQV3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQC3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQI3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQR3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=dbl_kind), pointer :: FQS3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=dbl_kind), pointer :: FQG3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=dbl_kind), pointer :: FQT3D(:,:,:,:,:) => NULL()
    
    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (1/s/s) 
    ! (due to advection, stretching, twisting, and Coriolis effect)   
    real (kind=dbl_kind), pointer :: FZX(:,:,:,:) => NULL() ! x-component of vorticity  
    real (kind=dbl_kind), pointer :: FZY(:,:,:,:) => NULL() ! y-component of vorticity 
    real (kind=dbl_kind), pointer :: FZTOP(:,:,:) => NULL() ! z-component of vorticity only at the top layer
    
    ! Tendencies due to turbulence.
    real (kind=dbl_kind), pointer :: THAD3(:,:,:)   => NULL() ! (K/s)
    real (kind=dbl_kind), pointer :: QVAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QCAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QIAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=dbl_kind), pointer :: QRAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=dbl_kind), pointer :: QSAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=dbl_kind), pointer :: QGAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QTAD3(:,:,:,:) => NULL() ! (kg/kg/s)
    
    real (kind=dbl_kind), pointer :: FZXTB(:,:,:) => NULL() ! x-component of vorticity (1/s/s) 
    real (kind=dbl_kind), pointer :: FZYTB(:,:,:) => NULL() ! y-component of vorticity (1/s/s) 
    real (kind=dbl_kind), pointer :: FZTOPB(:,:)  => NULL() ! z-component of vorticity (1/s/s) 
    
    real (kind=dbl_kind), pointer :: FU_TB(:,:,:) => NULL() ! x-component of momentum (m/s/s)
    real (kind=dbl_kind), pointer :: FV_TB(:,:,:) => NULL() ! y-component of momentum (m/s/s)
    
    ! Tendencies due to buoyancy.
    real (kind=dbl_kind), pointer :: FZXBU(:,:,:) => NULL() ! x-component of vorticity (1/s/s)
    real (kind=dbl_kind), pointer :: FZYBU(:,:,:) => NULL() ! y-component of vorticity (1/s/s) 
    
    ! Tendencies due to microphysics.
    real (kind=dbl_kind), pointer :: THAD_MICRO(:,:,:) => NULL() ! (K/s)
    real (kind=dbl_kind), pointer :: QVAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QCAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QIAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QRAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QSAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: QGAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)   
!-----------------    
!   Q3D_COUPLE (coupling) 
!-----------------  
    ! Relaxation time scale, sec 
    real (kind=dbl_kind), pointer :: TAU_RX_TQ(:,:,:) => NULL() ! for T (Q) 
    real (kind=dbl_kind), pointer :: TAU_RX_ZX(:,:,:) => NULL() ! for z3dx
    real (kind=dbl_kind), pointer :: TAU_RX_ZY(:,:,:) => NULL() ! for z3dy
    real (kind=dbl_kind), pointer :: TAU_RX_ZZ(:,:,:) => NULL() ! for z3dz (top layer only)
    
    ! Tendency due to diabatic effect (coupling) 
    real (kind=dbl_kind), pointer :: FTH3D_DIA(:,:,:) => NULL() ! (K/s)
    real (kind=dbl_kind), pointer :: FQV3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQC3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQI3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQR3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQS3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=dbl_kind), pointer :: FQG3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    
    real (kind=dbl_kind), pointer :: FQT3D_DIA(:,:,:,:) => NULL()
    
    real (kind=dbl_kind), pointer :: FU_DIA(:,:,:) => NULL() ! x-component of momentum (m/s/s)
    real (kind=dbl_kind), pointer :: FV_DIA(:,:,:) => NULL() ! y-component of momentum (m/s/s)
!----------------- 
!   Q3D_EDDY (coupling) 
!-----------------  
    ! Eddy Components
    real (kind=dbl_kind), pointer :: TH3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QV3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QC3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QI3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QR3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QS3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QG3D_ED(:,:,:)   => NULL() 
    real (kind=dbl_kind), pointer :: QT3D_ED(:,:,:,:) => NULL() 
    
    real (kind=dbl_kind), pointer :: U3DX_ED(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: U3DY_ED(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: W3D_ED(:,:,:)  => NULL() 
    
!----------------- 
!   Q3D_EDDY_EFFECT (coupling) 
!----------------- 
    real (kind=dbl_kind), pointer :: TH_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QV_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QC_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QI_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QR_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QS_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QG_E1(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QT_E1(:,:,:) => NULL()
    
    real (kind=dbl_kind), pointer :: U_E1(:,:) => NULL()
    real (kind=dbl_kind), pointer :: V_E1(:,:) => NULL()
    
    real (kind=dbl_kind), pointer :: TH_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QV_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QC_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QI_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QR_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QS_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QG_E2(:,:)   => NULL()
    real (kind=dbl_kind), pointer :: QT_E2(:,:,:) => NULL()
    
    real (kind=dbl_kind), pointer :: U_E2(:,:) => NULL()
    real (kind=dbl_kind), pointer :: V_E2(:,:) => NULL()    

!----------------- 
!   Q3D_EDDY_EFFECT (coupling)  2D Variables: e.g., sfx fluxes
!----------------- 
    real (kind=dbl_kind), pointer :: SPREC_E0(:) => NULL()
    real (kind=dbl_kind), pointer :: WTH_E0(:)   => NULL()
    real (kind=dbl_kind), pointer :: WQV_E0(:)   => NULL()
    real (kind=dbl_kind), pointer :: UW_E0(:)    => NULL()
    real (kind=dbl_kind), pointer :: WV_E0(:)    => NULL()
    
!-----------------    
!   MAPPING        
!-----------------  
    integer :: NFACE,NFACE_M,NFACE_P  ! index of cube face, indices of neighbors (-, +)
    
    ! Base of curvilinear coordinates, alpha [rad]
    real (kind=dbl_kind), pointer :: CALPHA_T(:) => NULL() ! at T-point   
    real (kind=dbl_kind), pointer :: CALPHA_U(:) => NULL() ! at U-point 
    real (kind=dbl_kind), pointer :: CALPHA_V(:) => NULL() ! at V-point 
    real (kind=dbl_kind), pointer :: CALPHA_Z(:) => NULL() ! at Z-point     
    
    ! Base of curvilinear coordinates, beta [rad]
    real (kind=dbl_kind), pointer :: CBETA_T(:)  => NULL() ! at T-point 
    real (kind=dbl_kind), pointer :: CBETA_U(:)  => NULL() ! at U-point  
    real (kind=dbl_kind), pointer :: CBETA_V(:)  => NULL() ! at V-point  
    real (kind=dbl_kind), pointer :: CBETA_Z(:)  => NULL() ! at Z-point 
    
    ! Longitude [rad] 
    real (kind=dbl_kind), pointer :: RLON_T(:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: RLON_U(:,:) => NULL() ! at U-point
    real (kind=dbl_kind), pointer :: RLON_V(:,:) => NULL() ! at V-point
    real (kind=dbl_kind), pointer :: RLON_Z(:,:) => NULL() ! at Z-point 
    
    ! Latitude [rad] 
    real (kind=dbl_kind), pointer :: RLAT_T(:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: RLAT_U(:,:) => NULL() ! at U-point
    real (kind=dbl_kind), pointer :: RLAT_V(:,:) => NULL() ! at V-point
    real (kind=dbl_kind), pointer :: RLAT_Z(:,:) => NULL() ! at Z-point 
    
    ! Coefficient of Coriolis Force defined at Z-point  
    real (kind=dbl_kind), pointer :: FVAL(:,:) => NULL() 
    
    ! Jacobian of the transformation
    real (kind=dbl_kind), pointer :: RG_T (:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: RG_U (:,:) => NULL() ! at U-point
    real (kind=dbl_kind), pointer :: RG_V (:,:) => NULL() ! at V-point
    real (kind=dbl_kind), pointer :: RG_Z (:,:) => NULL() ! at Z-point
    
    ! Coefficients of the contravariant metric tensor
    real (kind=dbl_kind), pointer :: GCONT_T(:,:,:) => NULL() ! at T-point 
    real (kind=dbl_kind), pointer :: GCONT_U(:,:,:) => NULL() ! at U-point 
    real (kind=dbl_kind), pointer :: GCONT_V(:,:,:) => NULL() ! at V-point 
    real (kind=dbl_kind), pointer :: GCONT_Z(:,:,:) => NULL() ! at Z-point 
    
    ! Secondary coefficients: used for calculating deformation & vorticity
    real (kind=dbl_kind), pointer :: RGG_T(:,:,:) => NULL() ! at T-point 
    real (kind=dbl_kind), pointer :: RGG_U(:,:,:) => NULL() ! at U-point 
    real (kind=dbl_kind), pointer :: RGG_V(:,:,:) => NULL() ! at V-point 
    real (kind=dbl_kind), pointer :: RGG_Z(:,:,:) => NULL() ! at Z-point 
    
    real (kind=dbl_kind), pointer :: GG_T(:,:,:) => NULL() ! at T-point 
    real (kind=dbl_kind), pointer :: GG_U(:,:,:) => NULL() ! at U-point 
    real (kind=dbl_kind), pointer :: GG_V(:,:,:) => NULL() ! at V-point 
    real (kind=dbl_kind), pointer :: GG_Z(:,:,:) => NULL() ! at Z-point 
    
    ! Coefficients of Transformation matrix A
    real (kind=dbl_kind), pointer :: AM_T(:,:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: AM_U(:,:,:) => NULL() ! at U-point
    real (kind=dbl_kind), pointer :: AM_V(:,:,:) => NULL() ! at V-point
    real (kind=dbl_kind), pointer :: AM_Z(:,:,:) => NULL() ! at Z-point
    
    ! Coefficients of the inverse of A
    real (kind=dbl_kind), pointer :: AMI_T(:,:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: AMI_U(:,:,:) => NULL() ! at U-point
    real (kind=dbl_kind), pointer :: AMI_V(:,:,:) => NULL() ! at V-point
    real (kind=dbl_kind), pointer :: AMI_Z(:,:,:) => NULL() ! at Z-point
    
    ! Coefficients of the inverse of A at the 1st halos of vGCM at the Edges 
    real (kind=dbl_kind) :: AMI_vGCM_halo(4,2)
    
    ! Halo points: location & alpha [rad]
    integer, pointer :: IHALO_LOC_T(:,:) => NULL() ! at T-point
    integer, pointer :: IHALO_LOC_Z(:,:) => NULL() ! at Z-point
    
    real (kind=dbl_kind), pointer :: HALPHA_T(:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: HALPHA_Z(:,:) => NULL() ! at Z-point
    
    ! Halo points: location & beta [rad]
    integer, pointer :: JHALO_LOC_T(:,:) => NULL() ! at T-point
    integer, pointer :: JHALO_LOC_Z(:,:) => NULL() ! at Z-point
    
    real (kind=dbl_kind), pointer :: HBETA_T(:,:) => NULL() ! at T-point
    real (kind=dbl_kind), pointer :: HBETA_Z(:,:) => NULL() ! at Z-point
    
    ! Halo points: mapping rules for vector in T-point
    real (kind=dbl_kind), pointer :: AM_VORT_ALPHA(:,:,:) => NULL() 
    real (kind=dbl_kind), pointer :: AM_VORT_BETA(:,:,:)  => NULL() 
!-----------------    
!   TOPO
!----------------- 
    ! Mountain heights
    real (kind=dbl_kind), pointer :: TOPOZ(:,:) => NULL()

    ! Vertical level index of mountain top: if hx=1, no mountain
    integer, pointer :: HX(:,:) => NULL()
    
    integer KLOWQ_MAX_GLOB  ! Lowest mountain-free layer in a channel
    
    integer, pointer :: KLOWQ_IJ(:,:) => NULL() ! Lowest mountain-free layer for T (Q) 
    integer, pointer :: KLOWU_IJ(:,:) => NULL() ! Lowest mountain-free layer for U
    integer, pointer :: KLOWV_IJ(:,:) => NULL() ! Lowest mountain-free layer for V

    integer, pointer :: NTOPOQ_GLOB(:) => NULL() ! # of topo-free cell in a channel: may not be used 
    integer, pointer :: NTOPOU_GLOB(:) => NULL() ! # of topo-free cell in a channel: may not be used 
    integer, pointer :: NTOPOV_GLOB(:) => NULL() ! # of topo-free cell in a channel: may not be used 
    
    integer, pointer :: DM_Q_TE(:,:,:) => NULL() ! used for q turbulence
    integer, pointer :: DM_Q_TW(:,:,:) => NULL() ! used for q turbulence
    integer, pointer :: DM_Q_TN(:,:,:) => NULL() ! used for q turbulence 
    integer, pointer :: DM_Q_TS(:,:,:) => NULL() ! used for q turbulence
    
    integer, pointer :: DM_KSI_XE(:,:,:) => NULL() ! used for z3dx advection
    integer, pointer :: DM_KSI_XC(:,:,:) => NULL() ! used for z3dx advection
    integer, pointer :: DM_KSI_YS(:,:,:) => NULL() ! used for z3dx advection
    
    integer, pointer :: DM_KSI_TE(:,:,:) => NULL() ! used for z3dx turbulence
    integer, pointer :: DM_KSI_TW(:,:,:) => NULL() ! used for z3dx turbulence
    integer, pointer :: DM_KSI_TN(:,:,:) => NULL() ! used for z3dx turbulence
    integer, pointer :: DM_KSI_TS(:,:,:) => NULL() ! used for z3dx turbulence
    
    integer, pointer :: DM_ETA_YN(:,:,:) => NULL() ! used for z3dy advection
    integer, pointer :: DM_ETA_YC(:,:,:) => NULL() ! used for z3dy advection
    integer, pointer :: DM_ETA_XW(:,:,:) => NULL() ! used for z3dy advection
    
    integer, pointer :: DM_ETA_TE(:,:,:) => NULL() ! used for z3dy turbulence
    integer, pointer :: DM_ETA_TW(:,:,:) => NULL() ! used for z3dy turbulence
    integer, pointer :: DM_ETA_TN(:,:,:) => NULL() ! used for z3dy turbulence
    integer, pointer :: DM_ETA_TS(:,:,:) => NULL() ! used for z3dy turbulence
!----------------- 
!   RADOUTLD 
!----------------- 
    ! Diagnostic output from the radiation parameterization.
    real (kind=dbl_kind), pointer :: FTHRAD(:,:,:) => NULL() ! tendency of potential temp. due to radiation (K/s)
    
    real (kind=dbl_kind), pointer :: FULWO(:,:,:)   => NULL() ! upward longwave flux (W/m**2)
    real (kind=dbl_kind), pointer :: FDLWO(:,:,:)   => NULL() ! downward longwave flux (W/m**2)
    real (kind=dbl_kind), pointer :: FUSWO(:,:,:)   => NULL() ! upward shortwave flux (W/m**2)
    real (kind=dbl_kind), pointer :: FDSWO(:,:,:)   => NULL() ! downward shortwave flux (W/m**2)
    real (kind=dbl_kind), pointer :: DTRADLW(:,:,:) => NULL() ! longwave heating rate (K/s)
    real (kind=dbl_kind), pointer :: DTRADSW(:,:,:) => NULL() ! shortwave heating rate (K/s)
      
    real (kind=dbl_kind), pointer :: WPLIQ(:,:,:) => NULL() ! liquid water path (g/m**2)
    real (kind=dbl_kind), pointer :: WPICE(:,:,:) => NULL() ! ice water path (g/m**2)
    real (kind=dbl_kind), pointer :: RELIQ(:,:,:) => NULL() ! effective radius of water (microns)
    real (kind=dbl_kind), pointer :: REICE(:,:,:) => NULL() ! effective radius of ice (microns)
    
    real (kind=dbl_kind), pointer :: FULWTOA(:,:) => NULL() ! upward longwave flux at the top of the atmosphere (W/m**2)
    real (kind=dbl_kind), pointer :: FUSWTOA(:,:) => NULL() ! upward shortwave flux at the top of the atmosphere (W/m**2)
    real (kind=dbl_kind), pointer :: FDSWTOA(:,:) => NULL() ! downward shortwave flux at the top of the atmosphere (W/m**2)
    
    real (kind=dbl_kind), pointer :: OLR(:,:) => NULL() ! outgoing long wave radiation (W/m**2)
    
    real (kind=dbl_kind), pointer :: dtlwavg(:) => NULL() ! area mean of longwave heating rate (K/s)
    real (kind=dbl_kind), pointer :: dtswavg(:) => NULL() ! area mean of shortwave heating rate (K/s)
!----------------------------------------------------------------------- 
    
  end type channel_seg_t

!***********************************************************************************************
! Define a type (2): Combination of 4 channel segments (a channel)
!***********************************************************************************************
  type, public :: channel_t
        integer :: num_chg       ! the number of channel group-id (1, 2, or 3)
        type(channel_seg_t), public :: seg(4)
  end type  channel_t

!***********************************************************************************************
! Define variables with derived type
!***********************************************************************************************

  type(channel_t), public, pointer :: channel_data(:) => NULL()

!***********************************************************************************************
  
  public :: allocate_channel_data

CONTAINS

  subroutine allocate_channel_data(channel, nch, isize, jsize)

    ! Dummy argument
    integer,         intent(in)    :: nch           ! Number of channels to be allocated  
    type(channel_t), intent(inout) :: channel(nch)  ! Channel variable to be allocated
    integer,         intent(in)    :: isize(4,nch)  ! Sizes of x-domain for 4 segments
    integer,         intent(in)    :: jsize(4,nch)  ! Sizes of y-domain for 4 segments

    ! Local data
    integer :: i,j

    do j = 1, nch
    do i = 1, 4
      call allocate_channel_section(channel(j)%seg(i), isize(i,j), jsize(i,j))
    enddo
    enddo

  end subroutine allocate_channel_data

  subroutine allocate_channel_section (segment,isize,jsize)

    ! Dummy argument
    type(channel_seg_t), intent(inout) :: segment    ! Segment variable to be allocated
    integer,             intent(in)    :: isize      ! Size of x-domain of a segment
    integer,             intent(in)    :: jsize      ! Size of x-domain of a segment

    segment%nx_size = isize
    segment%ny_size = jsize

!----------------------------------------------------------------------- 
!   BACKGROUND FIELDS
!-----------------------------------------------------------------------    
    ! Counterparts of Prognostic 3-D thermodynamical variables     
    allocate(segment%TH3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QV3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QC3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QI3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QR3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QS3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QG3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%QT3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2,ntracer))
       
    ! Counterparts of Prognostic 3-D dynamical variables
    allocate(segment%Z3DX_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%Z3DY_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%Z3DZ_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,1))
    
    ! Counterparts of Diagnostic 3-D dynamical variables
    allocate(segment%U3DX_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DY_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DX_CO_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DY_CO_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%W3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
        
!-----------------------------------------------------------------------    
!   CONST3D
!-----------------------------------------------------------------------    
    ! Prognostic 3-D thermodynamical variables     
    allocate(segment%TH3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QV3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QC3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QI3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QR3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QS3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QG3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%QT3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3,ntracer))
       
    ! Prognostic 3-D dynamical variables
    allocate(segment%Z3DX(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%Z3DY(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%Z3DZ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    
    ! Diagnostic 3-D dynamical variables
    allocate(segment%U3DX(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%U3DY(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%U3DX_CO(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%U3DY_CO(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
    allocate(segment%W3D(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%W3DNM1(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
        
    ! Diagnostic 2-D dynamical variables (uppermost layer)
    allocate(segment%PSI(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%PSINM1(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%CHI(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%CHINM1(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Diagnostic 2-D variables (top or surface) : Add more diagnostics (or remove some)
    allocate(segment%TG(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%ZROUGH(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GWET(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    allocate(segment%UW(isize,jsize))
    allocate(segment%WV(isize,jsize))
    allocate(segment%WTH(isize,jsize))
    allocate(segment%WQV(isize,jsize))
    allocate(segment%SPREC(isize,jsize))
    
!-----------------  
!   CONST3D_Extra
!----------------- 
    ! alternative space (used in halo_vort & relax_3d): multi-use
    allocate(segment%Z3DX0(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%Z3DY0(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%TERM1(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk3))
 
    ! deformation & shear components (used in turbulence and nonlinear diffusion)  
    allocate(segment%DEFYZ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%DEFXZ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%DEFXY(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%DEFXX(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%DEFYY(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%DEFZZ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))

    ! coefficients (used in turbulence and diffusion: initially calculated)  
    allocate(segment%COEFX(isize,jsize))   
    allocate(segment%COEFY(isize,jsize))   
    allocate(segment%COEFA(isize,jsize))       
    allocate(segment%COEFX_K(isize,jsize))
    allocate(segment%COEFY_K(isize,jsize)) 
    allocate(segment%COEFA_K(isize,jsize)) 
    allocate(segment%COEFX_E(isize,jsize)) 
    allocate(segment%COEFY_E(isize,jsize)) 
    allocate(segment%COEFA_E(isize,jsize))
    allocate(segment%COEFX_Z(isize,jsize)) 
    allocate(segment%COEFY_Z(isize,jsize)) 
    allocate(segment%COEFA_Z(isize,jsize))    

    ! coefficients (used in relaxation_3d: initially calculated)
    allocate(segment%AGAU_CO(isize,jsize,nk1-1))
    allocate(segment%BGAU_CO(isize,jsize,nk1-1))
    allocate(segment%CGAU_CO(isize,jsize,nk1-1))   
    
    ! coefficients (used in relaxation_2d: initially calculated)
    allocate(segment%RX2D_C1(isize,jsize))
    allocate(segment%RX2D_C2(isize,jsize))
      
!-----------------------------------------------------------------------  
!   CONST3D_TD    
!----------------------------------------------------------------------- 
    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (advection).    
    allocate(segment%FTH3D(isize,jsize,nk2,2))
    allocate(segment%FQV3D(isize,jsize,nk2,2))
    allocate(segment%FQC3D(isize,jsize,nk2,2))
    allocate(segment%FQI3D(isize,jsize,nk2,2))
    allocate(segment%FQR3D(isize,jsize,nk2,2))
    allocate(segment%FQS3D(isize,jsize,nk2,2))
    allocate(segment%FQG3D(isize,jsize,nk2,2))
    allocate(segment%FQT3D(isize,jsize,nk2,2,ntracer))
    
    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (1/s/s) 
    ! (due to advection, stretching, twisting, and Coriolis effect)   
    allocate(segment%FZX(isize,jsize,nk2,2))
    allocate(segment%FZY(isize,jsize,nk2,2)) 
    allocate(segment%FZTOP(isize,jsize,2))
    
    ! Tendencies due to turbulence.
    allocate(segment%THAD3(isize,jsize,nk2))
    allocate(segment%QVAD3(isize,jsize,nk2))
    allocate(segment%QCAD3(isize,jsize,nk2))
    allocate(segment%QIAD3(isize,jsize,nk2))
    allocate(segment%QTAD3(isize,jsize,nk2,ntracer))
    
    allocate(segment%FZXTB(isize,jsize,nk2))
    allocate(segment%FZYTB(isize,jsize,nk2))
    allocate(segment%FZTOPB(isize,jsize))
    
    allocate(segment%FU_TB(isize,jsize,nk2))
    allocate(segment%FV_TB(isize,jsize,nk2))
    
    ! Tendencies due to buoyancy.
    allocate(segment%FZXBU(isize,jsize,nk2))
    allocate(segment%FZYBU(isize,jsize,nk2))
    
    ! Tendencies due to microphysics.
    allocate(segment%THAD_MICRO(isize,jsize,nk2))
    allocate(segment%QVAD_MICRO(isize,jsize,nk2))
    allocate(segment%QCAD_MICRO(isize,jsize,nk2))
    allocate(segment%QIAD_MICRO(isize,jsize,nk2))
    allocate(segment%QRAD_MICRO(isize,jsize,nk2))
    allocate(segment%QSAD_MICRO(isize,jsize,nk2))
    allocate(segment%QGAD_MICRO(isize,jsize,nk2))   
!-----------------------------------------------------------------------   
!   Q3D_COUPLE (coupling) 
!----------------------------------------------------------------------- 
    ! Relaxation time scale, sec 
    allocate(segment%TAU_RX_TQ(isize,jsize,nk2))
    allocate(segment%TAU_RX_ZX(isize,jsize,nk2))
    allocate(segment%TAU_RX_ZY(isize,jsize,nk2))
    allocate(segment%TAU_RX_ZZ(isize,jsize,1))
    
    ! Tendency due to diabatic effect (coupling) 
    allocate(segment%FTH3D_DIA(isize,jsize,nk2))
    allocate(segment%FQV3D_DIA(isize,jsize,nk2))
    allocate(segment%FQC3D_DIA(isize,jsize,nk2))
    allocate(segment%FQI3D_DIA(isize,jsize,nk2))
    allocate(segment%FQR3D_DIA(isize,jsize,nk2))
    allocate(segment%FQS3D_DIA(isize,jsize,nk2))
    allocate(segment%FQG3D_DIA(isize,jsize,nk2))
    allocate(segment%FQT3D_DIA(isize,jsize,nk2,ntracer))
    
    allocate(segment%FU_DIA(isize,jsize,nk2))
    allocate(segment%FV_DIA(isize,jsize,nk2))
!-----------------------------------------------------------------------
!   Q3D_EDDY (coupling) 
!----------------------------------------------------------------------- 
    ! Eddy Components
    allocate(segment%TH3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QV3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QC3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QI3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QR3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QS3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QG3D_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%QT3D_ED(0:isize,0:jsize,nk2,ntracer)) 
    
    allocate(segment%U3DX_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%U3DY_ED(0:isize,0:jsize,nk2)) 
    allocate(segment%W3D_ED(0:isize,0:jsize,nk2)) 
!-----------------------------------------------------------------------
!   Q3D_EDDY_EFFECT (coupling) 
!----------------------------------------------------------------------- 
    allocate(segment%TH_E1(nk2,nVGCM_seg))
    allocate(segment%QV_E1(nk2,nVGCM_seg))
    allocate(segment%QC_E1(nk2,nVGCM_seg))
    allocate(segment%QI_E1(nk2,nVGCM_seg))
    allocate(segment%QR_E1(nk2,nVGCM_seg))
    allocate(segment%QS_E1(nk2,nVGCM_seg))
    allocate(segment%QG_E1(nk2,nVGCM_seg))
    allocate(segment%QT_E1(nk2,nVGCM_seg,ntracer))
    
    allocate(segment%U_E1(nk2,nVGCM_seg))
    allocate(segment%V_E1(nk2,nVGCM_seg))
    
    allocate(segment%TH_E2(nk2,nVGCM_seg))
    allocate(segment%QV_E2(nk2,nVGCM_seg))
    allocate(segment%QC_E2(nk2,nVGCM_seg))
    allocate(segment%QI_E2(nk2,nVGCM_seg))
    allocate(segment%QR_E2(nk2,nVGCM_seg))
    allocate(segment%QS_E2(nk2,nVGCM_seg))
    allocate(segment%QG_E2(nk2,nVGCM_seg))
    allocate(segment%QT_E2(nk2,nVGCM_seg,ntracer))
    
    allocate(segment%U_E2(nk2,nVGCM_seg))
    allocate(segment%V_E2(nk2,nVGCM_seg)) 
    
!-----------------------------------------------------------------------
!   Q3D_EDDY_EFFECT (coupling)  2D Variables: e.g., sfx fluxes
!----------------------------------------------------------------------- 
    allocate(segment%SPREC_E0(nVGCM_seg))
    
    allocate(segment%WTH_E0(nVGCM_seg))
    allocate(segment%WQV_E0(nVGCM_seg))
    allocate(segment%UW_E0(nVGCM_seg))
    allocate(segment%WV_E0(nVGCM_seg))
        
!-----------------------------------------------------------------------   
!   MAPPING        
!----------------------------------------------------------------------- 
    ! Base of curvilinear coordinates, alpha [rad]
    allocate(segment%CALPHA_T(1-nhalo:isize+nhalo))  
    allocate(segment%CALPHA_U(1-nhalo:isize+nhalo)) 
    allocate(segment%CALPHA_V(1-nhalo:isize+nhalo))
    allocate(segment%CALPHA_Z(1-nhalo:isize+nhalo))     
    
    ! Base of curvilinear coordinates, beta [rad]
    allocate(segment%CBETA_T(1-nhalo:jsize+nhalo))
    allocate(segment%CBETA_U(1-nhalo:jsize+nhalo))  
    allocate(segment%CBETA_V(1-nhalo:jsize+nhalo))
    allocate(segment%CBETA_Z(1-nhalo:jsize+nhalo))
    
    ! Longitude [rad] 
    allocate(segment%RLON_T(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLON_U(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLON_V(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLON_Z(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Latitude [rad] 
    allocate(segment%RLAT_T(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLAT_U(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLAT_V(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RLAT_Z(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Coefficient of Coriolis Force defined at Z-point  
    allocate(segment%FVAL(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Jacobian of the transformation
    allocate(segment%RG_T (1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RG_U (1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RG_V (1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RG_Z (1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Coefficients of the contravariant metric tensor
    allocate(segment%GCONT_T(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GCONT_U(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GCONT_V(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GCONT_Z(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo)) 
    
    ! Secondary coefficients: used for calculating deformation & vorticity
    allocate(segment%RGG_T(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RGG_U(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RGG_V(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%RGG_Z(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    allocate(segment%GG_T(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GG_U(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo)) 
    allocate(segment%GG_V(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo)) 
    allocate(segment%GG_Z(3,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Coefficients of Transformation matrix A
    allocate(segment%AM_T(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AM_U(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AM_V(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AM_Z(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    !  Coefficients of the inverse of A
    allocate(segment%AMI_T(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AMI_U(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AMI_V(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%AMI_Z(4,1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    ! Halo points: location & alpha [rad]
    allocate(segment%IHALO_LOC_T(1-nhalo:isize+nhalo,nhalo))
    allocate(segment%IHALO_LOC_Z(1-nhalo:isize+nhalo,nhalo))
    
    allocate(segment%HALPHA_T(1-nhalo:isize+nhalo,nhalo))
    allocate(segment%HALPHA_Z(1-nhalo:isize+nhalo,nhalo))
    
    ! Halo points: location & beta [rad]
    allocate(segment%JHALO_LOC_T(nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%JHALO_LOC_Z(nhalo,1-nhalo:jsize+nhalo))
    
    allocate(segment%HBETA_T(nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%HBETA_Z(nhalo,1-nhalo:jsize+nhalo))
    
    ! Halo points: mapping rules for vector in T-point
    allocate(segment%AM_VORT_ALPHA(4,1-nhalo:isize+nhalo,nhalo*2))
    allocate(segment%AM_VORT_BETA(4,nhalo*2,1-nhalo:jsize+nhalo))
!-----------------------------------------------------------------------    
!   TOPO
!----------------------------------------------------------------------- 
    ! Mountain heights
    allocate(segment%TOPOZ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))

    ! Vertical level index of mountain top: if hx=1, no mountain
    allocate(segment%HX(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    
    allocate(segment%KLOWQ_IJ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%KLOWU_IJ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%KLOWV_IJ(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))

    allocate(segment%NTOPOQ_GLOB(nk2))
    allocate(segment%NTOPOU_GLOB(nk2)) 
    allocate(segment%NTOPOV_GLOB(nk2)) 
    
    allocate(segment%DM_Q_TE(isize,jsize,nk1))
    allocate(segment%DM_Q_TW(isize,jsize,nk1))
    allocate(segment%DM_Q_TN(isize,jsize,nk1)) 
    allocate(segment%DM_Q_TS(isize,jsize,nk1))
    
    allocate(segment%DM_KSI_XE(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_KSI_XC(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_KSI_YS(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    
    allocate(segment%DM_KSI_TE(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_KSI_TW(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_KSI_TN(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_KSI_TS(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    
    allocate(segment%DM_ETA_YN(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_ETA_YC(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_ETA_XW(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    
    allocate(segment%DM_ETA_TE(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_ETA_TW(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_ETA_TN(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
    allocate(segment%DM_ETA_TS(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk1))
!----------------------------------------------------------------------- 
!   RADOUTLD 
!----------------------------------------------------------------------- 
    ! Diagnostic output from the radiation parameterization.
    allocate(segment%FTHRAD(isize,jsize,nk2))
    
    allocate(segment%FULWO(isize,jsize,nk2))
    allocate(segment%FDLWO(isize,jsize,nk2))
    allocate(segment%FUSWO(isize,jsize,nk2))
    allocate(segment%FDSWO(isize,jsize,nk2))
    allocate(segment%DTRADLW(isize,jsize,nk2))
    allocate(segment%DTRADSW(isize,jsize,nk2))
      
    allocate(segment%WPLIQ(isize,jsize,nk1))
    allocate(segment%WPICE(isize,jsize,nk1))
    allocate(segment%RELIQ(isize,jsize,nk1))
    allocate(segment%REICE(isize,jsize,nk1))
    
    allocate(segment%FULWTOA(isize,jsize))
    allocate(segment%FUSWTOA(isize,jsize))
    allocate(segment%FDSWTOA(isize,jsize))
    
    allocate(segment%OLR(isize,jsize))
    
    allocate(segment%dtlwavg(nk2))
    allocate(segment%dtswavg(nk2))
!-----------------------------------------------------------------------     
 
    ! Allocate additional fields here

  end subroutine allocate_channel_section

end module vvm_data_types
