module vvm_data_types

! The vertical size of (TH3D_BG, TH_E1, TH_E2, FTH3D_DIA, and RADOUTLD variables) 
! is "nk2_ext", not "nk2".  

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE parmsld,      only: nVGCM,nVGCM_seg,channel_seg_l,    &
                          nlevel,ntracer,nk1,nk2,nk3,nhalo, &
                          nk1_ext,nk2_ext

  implicit none
  private

!***********************************************************************************************
! Define a type (1): Variables declared for a channel segment (covering 1/4 of a channel)
!***********************************************************************************************

  type, public :: channel_seg_t
  
!------------------------------------------
!   BasicInfo (RESTART: need to recalculate)
!------------------------------------------ 
    integer :: nface                  ! Number of face where the segment belongs
    integer :: nx_size,ny_size        ! Sizes of x- and y-domain of a channel segment
    integer :: mi1,mim,mim_c,mim_ce,mim_a,mip,mip_c,mip_ce,mip_a
    integer :: mj1,mjm,mjm_c,mjm_ce,mjm_a,mjp,mjp_c,mjp_ce,mjp_a

    integer, allocatable :: lsta(:)       ! starting CRM-point index of a vGCM cell (along a channel)
    integer, allocatable :: lend(:)       ! ending CRM-point index of a vGCM cell
    real(kind=r8), allocatable :: lcen(:) ! center CRM-point index of a vGCM cell (vGCM point)
    
!------------------------------------------
!   BACKGROUND FIELDS
!------------------------------------------
    ! Counterparts of Prognostic 3-D thermodynamical variables
    real (kind=r8), pointer :: TH3D_bg(:,:,:)   => NULL() ! potential temperature (K)
    real (kind=r8), pointer :: QV3D_bg(:,:,:)   => NULL() ! water vapor mixing ratio (kg/kg)
    real (kind=r8), pointer :: QC3D_bg(:,:,:)   => NULL() ! cloud water mixing ratio (kg/kg)
    real (kind=r8), pointer :: QI3D_bg(:,:,:)   => NULL() ! cloud ice mixing ratio (kg/kg)
    real (kind=r8), pointer :: QR3D_bg(:,:,:)   => NULL() ! rain mixing ratio (kg/kg)
    real (kind=r8), pointer :: QS3D_bg(:,:,:)   => NULL() ! snow mixing ratio (kg/kg)
    real (kind=r8), pointer :: QG3D_bg(:,:,:)   => NULL() ! graupel mixing ratio (kg/kg)
    real (kind=r8), pointer :: QT3D_bg(:,:,:,:) => NULL() ! tracer mixing ratio (kg/kg)

    ! Counterparts of Prognostic 3-D dynamical variables
    real (kind=r8), pointer :: Z3DX_bg(:,:,:) => NULL() ! x-component of vorticity, dw/dy-dv/dz (1/s)
    real (kind=r8), pointer :: Z3DY_bg(:,:,:) => NULL() ! y-component of vorticity, du/dz-dw/dx (1/s)
    real (kind=r8), pointer :: Z3DZ_bg(:,:,:) => NULL() ! z-component of vorticity, dv/dx-du/dy (1/s)

    ! Counterparts of Diagnostic 3-D dynamical variables
    real (kind=r8), pointer :: U3DX_bg(:,:,:)    => NULL() ! (contravariant) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY_bg(:,:,:)    => NULL() ! (contravariant) meridional velocity, v (m/s)
    real (kind=r8), pointer :: U3DX_co_bg(:,:,:) => NULL() ! (covariant) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY_co_bg(:,:,:) => NULL() ! (covariant) meridional velocity, v (m/s)
    real (kind=r8), pointer :: W3D_bg(:,:,:)     => NULL() ! vertical velocity, w (m/s)
    
!   Debugging only (Remove later)    
    real (kind=r8), pointer :: U3DX_ll_bg(:,:,:) => NULL() ! (lon/lat) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY_ll_bg(:,:,:) => NULL() ! (lon/lat) meridional velocity, v (m/s)
    
    real (kind=r8), pointer :: U3DX_ll(:,:,:) => NULL() ! (lon/lat) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY_ll(:,:,:) => NULL() ! (lon/lat) meridional velocity, v (m/s)

!----------------------------------------------------------
!   pi-value of GCM at vGCM points (Q3D coupling Info)
!----------------------------------------------------------
    real (kind=r8), pointer :: Pint_bg(:,:,:)   => NULL() ! pressure at interfaces (Pa)
    real (kind=r8), pointer :: Pmid_bg(:,:,:)   => NULL() ! pressure at mid-layers (Pa)
    real (kind=r8), pointer :: Piint_bg(:,:,:)  => NULL() ! Exner ftn at interfaces (Pa)
    real (kind=r8), pointer :: Pimid_bg(:,:,:)  => NULL() ! Exner ftn at mid-layers (Pa) 
    real (kind=r8), pointer :: RHO_bg(:,:,:)    => NULL() ! density at mid-layers (kg/m**3)   
    real (kind=r8), pointer :: RHOint_bg(:,:,:) => NULL() ! density at interfaces (kg/m**3)    
    
    real (kind=r8), pointer :: PI_G(:,:)   => NULL()
    real (kind=r8), pointer :: ZM_G(:,:)   => NULL()    
    
!------------------------------------------
!   CONST3D (RESTART: read data, except for LOCEAN)
!------------------------------------------
    ! Prognostic 3-D thermodynamical variables
    real (kind=r8), pointer :: TH3D(:,:,:)   => NULL() ! potential temperature (K)
    real (kind=r8), pointer :: QV3D(:,:,:)   => NULL() ! water vapor mixing ratio (kg/kg)
    real (kind=r8), pointer :: QC3D(:,:,:)   => NULL() ! cloud water mixing ratio (kg/kg)
    real (kind=r8), pointer :: QI3D(:,:,:)   => NULL() ! cloud ice mixing ratio (kg/kg)
    real (kind=r8), pointer :: QR3D(:,:,:)   => NULL() ! rain mixing ratio (kg/kg)
    real (kind=r8), pointer :: QS3D(:,:,:)   => NULL() ! snow mixing ratio (kg/kg)
    real (kind=r8), pointer :: QG3D(:,:,:)   => NULL() ! graupel mixing ratio (kg/kg)
    real (kind=r8), pointer :: QT3D(:,:,:,:) => NULL() ! tracer mixing ratio (kg/kg)

    ! Prognostic 3-D dynamical variables
    real (kind=r8), pointer :: Z3DX(:,:,:) => NULL() ! x-component of vorticity, dw/dy-dv/dz (1/s)
    real (kind=r8), pointer :: Z3DY(:,:,:) => NULL() ! y-component of vorticity, du/dz-dw/dx (1/s)
    real (kind=r8), pointer :: Z3DZ(:,:,:) => NULL() ! z-component of vorticity, dv/dx-du/dy (1/s)

    ! Diagnostic 3-D dynamical variables
    real (kind=r8), pointer :: U3DX(:,:,:)    => NULL() ! (contravariant) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY(:,:,:)    => NULL() ! (contravariant) meridional velocity, v (m/s)
    real (kind=r8), pointer :: U3DX_CO(:,:,:) => NULL() ! (covariant) zonal velocity, u (m/s)
    real (kind=r8), pointer :: U3DY_CO(:,:,:) => NULL() ! (covariant) meridional velocity, v (m/s)
    real (kind=r8), pointer :: W3D(:,:,:)     => NULL() ! vertical velocity, w (m/s)
    real (kind=r8), pointer :: W3DNM1(:,:,:)  => NULL() ! vertical velocity, w (m/s), in the previous timestep

    ! Diagnostic 2-D dynamical variables (uppermost layer)
    real (kind=r8), pointer :: PSI(:,:)    => NULL() ! stream function (m**2/s)
    real (kind=r8), pointer :: PSINM1(:,:) => NULL() ! stream function (m**2/s), in the previous timestep
    real (kind=r8), pointer :: CHI(:,:)    => NULL() ! velocity potential (m**2/s)
    real (kind=r8), pointer :: CHINM1(:,:) => NULL() ! velocity potential (m**2/s), in the previous timestep

    ! Diagnostic 2-D variables (top or surface) : Add more diagnostics (or remove some)
    logical, pointer :: LOCEAN(:,:) => NULL()  ! logical variable for being ocean

    real (kind=r8), pointer :: TG(:,:)     => NULL() ! ground temperature (K)
    real (kind=r8), pointer :: ZROUGH(:,:) => NULL() ! roughness length (m)
    real (kind=r8), pointer :: GWET(:,:)   => NULL() ! ground wetness

    real (kind=r8), pointer :: UW(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2
    real (kind=r8), pointer :: WV(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2
    real (kind=r8), pointer :: UW_CON(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2
    real (kind=r8), pointer :: WV_CON(:,:) => NULL() ! surface momentum flux (kg/m**3)(m/s)**2

    real (kind=r8), pointer :: WTH(:,:) => NULL() ! surface heat (potential temp.) flux (kg/m**3)(K)(m/s)
    real (kind=r8), pointer :: WQV(:,:) => NULL() ! surface moisture flux (kg/m**3)(kg/kg)(m/s)

    real (kind=r8), pointer :: SPREC(:,:) => NULL() ! surface precipitation rate (kg/m**2/s): sprec*3600. (mm/hr)

!--------------------------------------------------------
!   CONST3D_Extra_Initialized (RESTART: recalculated)
!--------------------------------------------------------
    ! coefficients (used in turbulence and diffusion: initially calculated)
    real (kind=r8), pointer :: COEFX(:,:) => NULL()
    real (kind=r8), pointer :: COEFY(:,:) => NULL()
    real (kind=r8), pointer :: COEFA(:,:) => NULL()
    real (kind=r8), pointer :: COEFX_K(:,:) => NULL()
    real (kind=r8), pointer :: COEFY_K(:,:) => NULL()
    real (kind=r8), pointer :: COEFA_K(:,:) => NULL()
    real (kind=r8), pointer :: COEFX_E(:,:) => NULL()
    real (kind=r8), pointer :: COEFY_E(:,:) => NULL()
    real (kind=r8), pointer :: COEFA_E(:,:) => NULL()
    real (kind=r8), pointer :: COEFX_Z(:,:) => NULL()
    real (kind=r8), pointer :: COEFY_Z(:,:) => NULL()
    real (kind=r8), pointer :: COEFA_Z(:,:) => NULL()

    ! coefficients (used in relaxation_3d: initially calculated)
    real (kind=r8), pointer :: AGAU_CO(:,:,:) => NULL()
    real (kind=r8), pointer :: BGAU_CO(:,:,:) => NULL()
    real (kind=r8), pointer :: CGAU_CO(:,:,:) => NULL()

    ! coefficients (used in relaxation_2d: initially calculated)
    real (kind=r8), pointer :: RX2D_C1(:,:) => NULL()
    real (kind=r8), pointer :: RX2D_C2(:,:) => NULL()
    
!------------------------------------------
!   CONST3D_Extra
!------------------------------------------
    ! alternative space (used in halo_vort & relax_3d): multi-use
    real (kind=r8), pointer :: Z3DX0(:,:,:) => NULL()
    real (kind=r8), pointer :: Z3DY0(:,:,:) => NULL()
    real (kind=r8), pointer :: TERM1(:,:,:) => NULL()

    ! deformation & shear components (used in turbulence and nonlinear diffusion)
    real (kind=r8), pointer :: DEFYZ(:,:,:) => NULL()
    real (kind=r8), pointer :: DEFXZ(:,:,:) => NULL()
    real (kind=r8), pointer :: DEFXY(:,:,:) => NULL()
    real (kind=r8), pointer :: DEFXX(:,:,:) => NULL()
    real (kind=r8), pointer :: DEFYY(:,:,:) => NULL()
    real (kind=r8), pointer :: DEFZZ(:,:,:) => NULL()

    ! turbulence coefficients
    real (kind=r8), pointer :: RKH(:,:,:) => NULL()
    real (kind=r8), pointer :: RKM(:,:,:) => NULL()

!------------------------------------------
!   CONST3D_TD
!------------------------------------------
    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (advection).
    real (kind=r8), pointer :: FTH3D(:,:,:,:)   => NULL() ! (K/s) & random perturbation
    real (kind=r8), pointer :: FQV3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQC3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQI3D(:,:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQR3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=r8), pointer :: FQS3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=r8), pointer :: FQG3D(:,:,:,:)   => NULL() ! (kg/kg/s) & falling with terminal velocity
    real (kind=r8), pointer :: FQT3D(:,:,:,:,:) => NULL()

    ! Tendencies used in Adams-Bashforth 2nd-order time scheme (1/s/s)
    ! (due to advection, stretching, twisting, and Coriolis effect)
    real (kind=r8), pointer :: FZX(:,:,:,:) => NULL() ! x-component of vorticity
    real (kind=r8), pointer :: FZY(:,:,:,:) => NULL() ! y-component of vorticity
    real (kind=r8), pointer :: FZTOP(:,:,:) => NULL() ! z-component of vorticity only at the top layer

    ! Tendencies due to turbulence.
    real (kind=r8), pointer :: THAD3(:,:,:)   => NULL() ! (K/s)
    real (kind=r8), pointer :: QVAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QCAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QIAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=r8), pointer :: QRAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=r8), pointer :: QSAD3(:,:,:)   => NULL() ! (kg/kg/s)
!     real (kind=r8), pointer :: QGAD3(:,:,:)   => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QTAD3(:,:,:,:) => NULL() ! (kg/kg/s)

    real (kind=r8), pointer :: FZXTB(:,:,:) => NULL() ! x-component of vorticity (1/s/s)
    real (kind=r8), pointer :: FZYTB(:,:,:) => NULL() ! y-component of vorticity (1/s/s)
    real (kind=r8), pointer :: FZTOPB(:,:)  => NULL() ! z-component of vorticity (1/s/s)

    ! Tendencies due to buoyancy.
    real (kind=r8), pointer :: FZXBU(:,:,:) => NULL() ! x-component of vorticity (1/s/s)
    real (kind=r8), pointer :: FZYBU(:,:,:) => NULL() ! y-component of vorticity (1/s/s)

    ! Tendencies due to microphysics.
    real (kind=r8), pointer :: THAD_MICRO(:,:,:) => NULL() ! (K/s)
    real (kind=r8), pointer :: QVAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QCAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QIAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QRAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QSAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: QGAD_MICRO(:,:,:) => NULL() ! (kg/kg/s)

!------------------------------------------
!   Q3D_COUPLE (coupling)
!------------------------------------------
    ! Relaxation time scale, sec
    real (kind=r8), pointer :: TAU_RX_TQ(:,:,:) => NULL() ! for T (Q)
    real (kind=r8), pointer :: TAU_RX_ZX(:,:,:) => NULL() ! for z3dx
    real (kind=r8), pointer :: TAU_RX_ZY(:,:,:) => NULL() ! for z3dy
    real (kind=r8), pointer :: TAU_RX_ZZ(:,:,:) => NULL() ! for z3dz (top layer only)

    ! Tendency due to diabatic effect (coupling)
    real (kind=r8), pointer :: FTH3D_DIA(:,:,:) => NULL() ! (K/s)
    real (kind=r8), pointer :: FQV3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQC3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQI3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQR3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQS3D_DIA(:,:,:) => NULL() ! (kg/kg/s)
    real (kind=r8), pointer :: FQG3D_DIA(:,:,:) => NULL() ! (kg/kg/s)

    real (kind=r8), pointer :: FQT3D_DIA(:,:,:,:) => NULL()

    real (kind=r8), pointer :: FU_DIA(:,:,:) => NULL() ! x-component of momentum (m/s/s)
    real (kind=r8), pointer :: FV_DIA(:,:,:) => NULL() ! y-component of momentum (m/s/s)

!------------------------------------------
!   Q3D_EDDY (coupling)
!------------------------------------------
    ! Eddy Components
    real (kind=r8), pointer :: TH3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QV3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QC3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QI3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QR3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QS3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QG3D_ED(:,:,:)   => NULL()
    real (kind=r8), pointer :: QT3D_ED(:,:,:,:) => NULL()

    real (kind=r8), pointer :: U3DX_ED(:,:,:) => NULL()
    real (kind=r8), pointer :: U3DY_ED(:,:,:) => NULL()
    real (kind=r8), pointer :: W3D_ED(:,:,:)  => NULL()
    
!------------------------------------------
!   Q3D_EDDY_EFFECT (coupling)
!------------------------------------------
    real (kind=r8), pointer :: TH_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QV_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QC_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QI_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QR_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QS_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QG_E1(:,:)   => NULL()
    real (kind=r8), pointer :: QT_E1(:,:,:) => NULL()

    real (kind=r8), pointer :: U_E1(:,:) => NULL()
    real (kind=r8), pointer :: V_E1(:,:) => NULL()

    real (kind=r8), pointer :: TH_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QV_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QC_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QI_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QR_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QS_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QG_E2(:,:)   => NULL()
    real (kind=r8), pointer :: QT_E2(:,:,:) => NULL()

    real (kind=r8), pointer :: U_E2(:,:) => NULL()
    real (kind=r8), pointer :: V_E2(:,:) => NULL()

!------------------------------------------
!   Q3D_OUT (coupling)  
!------------------------------------------
!   1. 2D Variables (Time-averaged values) 
    real (kind=r8), pointer :: SPREC_E0(:) => NULL()
    real (kind=r8), pointer :: WTH_E0(:)   => NULL()
    real (kind=r8), pointer :: WQV_E0(:)   => NULL()
    real (kind=r8), pointer :: UW_E0(:)    => NULL()
    real (kind=r8), pointer :: WV_E0(:)    => NULL()
    
!   2. 3D Variables (State values)
    real (kind=r8), pointer :: TH3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QV3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QC3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QI3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QR3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QS3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QG3D_sa(:,:)   => NULL()
    real (kind=r8), pointer :: QT3D_sa(:,:,:) => NULL()

    real (kind=r8), pointer :: U3DX_sa(:,:)   => NULL()
    real (kind=r8), pointer :: U3DY_sa(:,:)   => NULL()
    real (kind=r8), pointer :: W3D_sa(:,:)    => NULL()    

!------------------------------------------
!   MAPPING (RESTART: need to recalculate)
!------------------------------------------
    ! Base of curvilinear coordinates, alpha [rad]
    real (kind=r8), pointer :: CALPHA_T(:) => NULL() ! at T-point
    real (kind=r8), pointer :: CALPHA_U(:) => NULL() ! at U-point
    real (kind=r8), pointer :: CALPHA_V(:) => NULL() ! at V-point
    real (kind=r8), pointer :: CALPHA_Z(:) => NULL() ! at Z-point

    ! Base of curvilinear coordinates, beta [rad]
    real (kind=r8), pointer :: CBETA_T(:)  => NULL() ! at T-point
    real (kind=r8), pointer :: CBETA_U(:)  => NULL() ! at U-point
    real (kind=r8), pointer :: CBETA_V(:)  => NULL() ! at V-point
    real (kind=r8), pointer :: CBETA_Z(:)  => NULL() ! at Z-point

    ! Longitude [rad]
    real (kind=r8), pointer :: RLON_T(:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: RLON_U(:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: RLON_V(:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: RLON_Z(:,:) => NULL() ! at Z-point

    ! Latitude [rad]
    real (kind=r8), pointer :: RLAT_T(:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: RLAT_U(:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: RLAT_V(:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: RLAT_Z(:,:) => NULL() ! at Z-point

    ! Coefficient of Coriolis Force defined at Z-point
    real (kind=r8), pointer :: FVAL(:,:) => NULL()

    ! Jacobian of the transformation
    real (kind=r8), pointer :: RG_T (:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: RG_U (:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: RG_V (:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: RG_Z (:,:) => NULL() ! at Z-point

    ! Coefficients of the contravariant metric tensor
    real (kind=r8), pointer :: GCONT_T(:,:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: GCONT_U(:,:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: GCONT_V(:,:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: GCONT_Z(:,:,:) => NULL() ! at Z-point

    ! Secondary coefficients: used for calculating deformation & vorticity
    real (kind=r8), pointer :: RGG_T(:,:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: RGG_U(:,:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: RGG_V(:,:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: RGG_Z(:,:,:) => NULL() ! at Z-point

    real (kind=r8), pointer :: GG_T(:,:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: GG_U(:,:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: GG_V(:,:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: GG_Z(:,:,:) => NULL() ! at Z-point

    ! Coefficients of Transformation matrix A
    real (kind=r8), pointer :: AM_T(:,:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: AM_U(:,:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: AM_V(:,:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: AM_Z(:,:,:) => NULL() ! at Z-point

    ! Coefficients of the inverse of A
    real (kind=r8), pointer :: AMI_T(:,:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: AMI_U(:,:,:) => NULL() ! at U-point
    real (kind=r8), pointer :: AMI_V(:,:,:) => NULL() ! at V-point
    real (kind=r8), pointer :: AMI_Z(:,:,:) => NULL() ! at Z-point

    ! Coefficients of the inverse of A at the 1st halos of vGCM at the Edges
    real (kind=r8) :: AMI_vGCM_halo(4,2)

    ! Halo points: location & alpha [rad]
    integer, pointer :: IHALO_LOC_T(:,:) => NULL() ! at T-point
    integer, pointer :: IHALO_LOC_Z(:,:) => NULL() ! at Z-point

    real (kind=r8), pointer :: HALPHA_T(:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: HALPHA_Z(:,:) => NULL() ! at Z-point

    ! Halo points: location & beta [rad]
    integer, pointer :: JHALO_LOC_T(:,:) => NULL() ! at T-point
    integer, pointer :: JHALO_LOC_Z(:,:) => NULL() ! at Z-point

    real (kind=r8), pointer :: HBETA_T(:,:) => NULL() ! at T-point
    real (kind=r8), pointer :: HBETA_Z(:,:) => NULL() ! at Z-point

    ! Halo points: mapping rules for vector in T-point
    real (kind=r8), pointer :: AM_VORT_ALPHA(:,:,:) => NULL()
    real (kind=r8), pointer :: AM_VORT_BETA(:,:,:)  => NULL()
    
!------------------------------------------
!   TOPO (RESTART: need to recalculate)
!------------------------------------------
    ! Mountain heights
    real (kind=r8), pointer :: TOPOZ(:,:) => NULL()

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
!------------------------------------------
!   RADOUTLD
!------------------------------------------
    ! Diagnostic output from the radiation parameterization.
    real (kind=r8), pointer :: FTHRAD(:,:,:) => NULL() ! tendency of potential temp. due to radiation (K/s)

    real (kind=r8), pointer :: FULWO(:,:,:)   => NULL() ! upward longwave flux (W/m**2)
    real (kind=r8), pointer :: FDLWO(:,:,:)   => NULL() ! downward longwave flux (W/m**2)
    real (kind=r8), pointer :: FUSWO(:,:,:)   => NULL() ! upward shortwave flux (W/m**2)
    real (kind=r8), pointer :: FDSWO(:,:,:)   => NULL() ! downward shortwave flux (W/m**2)
    real (kind=r8), pointer :: DTRADLW(:,:,:) => NULL() ! longwave heating rate (K/s)
    real (kind=r8), pointer :: DTRADSW(:,:,:) => NULL() ! shortwave heating rate (K/s)

    real (kind=r8), pointer :: WPLIQ(:,:,:) => NULL() ! liquid water path (g/m**2)
    real (kind=r8), pointer :: WPICE(:,:,:) => NULL() ! ice water path (g/m**2)
    real (kind=r8), pointer :: RELIQ(:,:,:) => NULL() ! effective radius of water (microns)
    real (kind=r8), pointer :: REICE(:,:,:) => NULL() ! effective radius of ice (microns)

    real (kind=r8), pointer :: FULWTOA(:,:) => NULL() ! upward longwave flux at the top of the atmosphere (W/m**2)
    real (kind=r8), pointer :: FUSWTOA(:,:) => NULL() ! upward shortwave flux at the top of the atmosphere (W/m**2)
    real (kind=r8), pointer :: FDSWTOA(:,:) => NULL() ! downward shortwave flux at the top of the atmosphere (W/m**2)

    real (kind=r8), pointer :: OLR(:,:) => NULL() ! outgoing long wave radiation (W/m**2)

    real (kind=r8), pointer :: dtlwavg(:) => NULL() ! area mean of longwave heating rate (K/s)
    real (kind=r8), pointer :: dtswavg(:) => NULL() ! area mean of shortwave heating rate (K/s)
!-----------------------------------------------------------------------

  end type channel_seg_t

!***********************************************************************************************
! Define a type (2): Combination of 4 channel segments (a channel)
!***********************************************************************************************
  type, public :: channel_t
        logical :: isFirst_psi, isFirst_chi 
        integer :: num_chg     ! the number of channel group-id (1, 2, or 3) (needed for RESTART)
        integer :: num_chn     ! the number of channel
        type(channel_seg_t), public :: seg(4)
  end type  channel_t

!***********************************************************************************************
! Define variables with derived type
!***********************************************************************************************

  type(channel_t), public, pointer :: channel_data(:) => NULL()

!***********************************************************************************************

  public :: allocate_channel_data

CONTAINS

  subroutine allocate_channel_data(begchan, endchan, channel)

    ! Dummy argument
    integer,         intent(in)    :: begchan       ! First channel number to be allocated
    integer,         intent(in)    :: endchan       ! Last channel number to be allocated
    type(channel_t), intent(inout) :: channel(begchan:endchan)  ! Channel variable to be allocated

    ! Local data
    integer :: i,j

    if (endchan > 0) then
       ! There may not be any channels in this task.
       do j = begchan, endchan
          do i = 1, 4
             call allocate_channel_section(channel(j)%seg(i), channel(j)%seg(i)%nx_size, &
                                                              channel(j)%seg(i)%ny_size)
          end do
       end do
    end if

  end subroutine allocate_channel_data

  subroutine allocate_channel_section (segment,isize,jsize)

    ! Dummy argument
    type(channel_seg_t), intent(inout) :: segment    ! Segment variable to be allocated
    integer,             intent(in)    :: isize      ! Size of x-domain of a segment
    integer,             intent(in)    :: jsize      ! Size of x-domain of a segment

    segment%nx_size = isize
    segment%ny_size = jsize

    allocate(segment%lsta(nVGCM_seg))    
    allocate(segment%lend(nVGCM_seg))
    allocate(segment%lcen(nVGCM_seg))

!-----------------------------------------------------------------------
!   BACKGROUND FIELDS:
!   # of halo values needed across the channel = nhalo
!   
!   # of halo values needed along the channel
!   Thermodynamic variables     : none
!   U3DX_bg & U3DY_bg           :  1 (2:nk1)   nhalo (nk2)
!   U3DX_CO_bg & U3DY_CO_bg     :  1 (2:nk1)   nhalo (nk2)
!   W3D_bg                      :  1 
!   Z3DX_bg & Z3DY_BG & Z3DZ_BG : none
!-----------------------------------------------------------------------
    ! Counterparts of Prognostic 3-D thermodynamical variables
    allocate(segment%TH3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2_ext))
    
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
    
    allocate(segment%Z3DZ_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,2))      
    ! K = 1 : Used for predicting Z3DZ at the uppermost layer
    ! K = 2 : Used for solving the 2D elliptic equation at the lowest layer 

    ! Counterparts of Diagnostic 3-D dynamical variables
    allocate(segment%U3DX_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DY_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DX_CO_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DY_CO_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%W3D_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))

!   Debugging only (remove later)    
    allocate(segment%U3DX_ll_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2))
    allocate(segment%U3DY_ll_bg(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo,nk2)) 
    
    allocate(segment%U3DX_ll(isize,jsize,nk2))
    allocate(segment%U3DY_ll(isize,jsize,nk2))  

!-----------------------------------------------------------------------
!   pi-value of GCM at vGCM points (Q3D coupling Info)
!-----------------------------------------------------------------------
    allocate(segment%Pint_bg(0:isize+1,0:jsize+1,nk2_ext))
    allocate(segment%Pmid_bg(0:isize+1,0:jsize+1,nk2_ext))
    allocate(segment%PIint_bg(0:isize+1,0:jsize+1,nk2_ext))
    allocate(segment%PImid_bg(0:isize+1,0:jsize+1,nk2_ext))
    allocate(segment%RHO_bg(0:isize+1,0:jsize+1,nk2_ext))
    allocate(segment%RHOint_bg(0:isize+1,0:jsize+1,nk2_ext))
     
    allocate(segment%PI_G(nLevel,nVGCM_seg)) 
    allocate(segment%ZM_G(nLevel,nVGCM_seg))  
    
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
    allocate(segment%LOCEAN(1:isize+1,1:jsize+1))

    allocate(segment%TG(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%ZROUGH(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))
    allocate(segment%GWET(1-nhalo:isize+nhalo,1-nhalo:jsize+nhalo))

    allocate(segment%UW(isize+1,jsize+1))
    allocate(segment%WV(isize+1,jsize+1))
    allocate(segment%UW_CON(isize+1,jsize+1))
    allocate(segment%WV_CON(isize+1,jsize+1))

    allocate(segment%WTH(isize,jsize))
    allocate(segment%WQV(isize,jsize))
    allocate(segment%SPREC(isize,jsize))

!-----------------
!   CONST3D_Extra_Initialized
!-----------------
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

    ! turbulence coefficients
    allocate(segment%RKH(0:isize+1,0:jsize+1,nk2))
    allocate(segment%RKM(0:isize+1,0:jsize+1,nk2))

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
    allocate(segment%FTH3D_DIA(isize,jsize,nk2_ext))
    
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
    allocate(segment%TH_E1(nk2_ext,nVGCM_seg))
    
    allocate(segment%QV_E1(nk2,nVGCM_seg))
    allocate(segment%QC_E1(nk2,nVGCM_seg))
    allocate(segment%QI_E1(nk2,nVGCM_seg))
    allocate(segment%QR_E1(nk2,nVGCM_seg))
    allocate(segment%QS_E1(nk2,nVGCM_seg))
    allocate(segment%QG_E1(nk2,nVGCM_seg))
    allocate(segment%QT_E1(nk2,nVGCM_seg,ntracer))
    allocate(segment%U_E1(nk2,nVGCM_seg))
    allocate(segment%V_E1(nk2,nVGCM_seg))

    allocate(segment%TH_E2(nk2_ext,nVGCM_seg))
    
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
!   Q3D_OUT (coupling)  
!-----------------------------------------------------------------------
!   1. 2D Variables (Time-averaged values) 
    allocate(segment%SPREC_E0(nVGCM_seg))
    allocate(segment%WTH_E0(nVGCM_seg))

    allocate(segment%WQV_E0(nVGCM_seg))
    allocate(segment%UW_E0(nVGCM_seg))
    allocate(segment%WV_E0(nVGCM_seg))
    
!   2. 3D Variables (State values)
    allocate(segment%TH3D_sa(nk2,nVGCM_seg))
    allocate(segment%QV3D_sa(nk2,nVGCM_seg))
    allocate(segment%QC3D_sa(nk2,nVGCM_seg))
    allocate(segment%QI3D_sa(nk2,nVGCM_seg))
    allocate(segment%QR3D_sa(nk2,nVGCM_seg))
    allocate(segment%QS3D_sa(nk2,nVGCM_seg))
    allocate(segment%QG3D_sa(nk2,nVGCM_seg))
    allocate(segment%QT3D_sa(nk2,nVGCM_seg,ntracer))

    allocate(segment%U3DX_sa(nk2,nVGCM_seg))
    allocate(segment%U3DY_sa(nk2,nVGCM_seg))
    allocate(segment%W3D_sa(nk2,nVGCM_seg))   

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
!   RADOUTLD: Use extended vertical layers
!-----------------------------------------------------------------------
    ! Diagnostic output from the radiation parameterization.
    allocate(segment%FTHRAD(isize,jsize,nk2_ext))

    allocate(segment%FULWO(isize,jsize,nk2_ext))
    allocate(segment%FDLWO(isize,jsize,nk2_ext))
    allocate(segment%FUSWO(isize,jsize,nk2_ext))
    allocate(segment%FDSWO(isize,jsize,nk2_ext))
    allocate(segment%DTRADLW(isize,jsize,nk2_ext))
    allocate(segment%DTRADSW(isize,jsize,nk2_ext))

    allocate(segment%WPLIQ(isize,jsize,nk1_ext))
    allocate(segment%WPICE(isize,jsize,nk1_ext))
    allocate(segment%RELIQ(isize,jsize,nk1_ext))
    allocate(segment%REICE(isize,jsize,nk1_ext))

    allocate(segment%FULWTOA(isize,jsize))
    allocate(segment%FUSWTOA(isize,jsize))
    allocate(segment%FDSWTOA(isize,jsize))

    allocate(segment%OLR(isize,jsize))

    allocate(segment%dtlwavg(nk2_ext))
    allocate(segment%dtswavg(nk2_ext))
!-----------------------------------------------------------------------

    ! Allocate additional fields here

  end subroutine allocate_channel_section

end module vvm_data_types
