module radae
!------------------------------------------------------------------------------
!
! Description:
!
! Data and subroutines to calculate absorptivities and emissivity needed
! for the LW radiation calculation.
!
! Public interfaces are: 
!
! radae_init ------------ Initialization
! initialize_radbuffer -- Initialize the 3D abs/emis arrays.
! radabs ---------------- Compute absorptivities.
! radems ---------------- Compute emissivity.
! radtpl ---------------- Compute Temperatures and path lengths.
! radoz2 ---------------- Compute ozone path lengths.
! trcpth ---------------- Compute ghg path lengths.
!
! Author:  B. Collins
!
! $Id$
!------------------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use ppgrid,           only: pcols, pverp, begchunk, endchunk, pver
  use infnan,           only: posinf, assignment(=)
  use pmgrid,           only: plev, plevp
  use radconstants,     only: nlwbands, idx_LW_0650_0800, idx_LW_0500_0650, &
        idx_LW_1000_1200, idx_LW_0800_1000,  idx_LW_1200_2000
  use abortutils,       only: endrun
  use cam_logfile,      only: iulog
  use wv_saturation,    only: qsat_water


  implicit none

  save
!-----------------------------------------------------------------------------
! PUBLIC:: By default data and interfaces are private
!-----------------------------------------------------------------------------
  private
  public radabs, radems, radtpl, radae_init, initialize_radbuffer, radoz2, trcpth

  integer, public, parameter :: nbands = 2          ! Number of spectral bands
!
! Following data needed for restarts and in radclwmx
!
  real(r8), public, allocatable, target :: abstot_3d(:,:,:,:) ! Non-adjacent layer absorptivites
  real(r8), public, allocatable, target :: absnxt_3d(:,:,:,:) ! Nearest layer absorptivities
  real(r8), public, allocatable, target :: emstot_3d(:,:,:)   ! Total emissivity
  integer,  public :: ntoplw    ! top level to solve for longwave cooling

!-----------------------------------------------------------------------------
! PRIVATE:: The rest of the data is private to this module.
!-----------------------------------------------------------------------------
  real(r8) :: p0    ! Standard pressure (dynes/cm**2)
  real(r8) :: amd   ! Molecular weight of dry air (g/mol)
  real(r8) :: amco2 ! Molecular weight of co2   (g/mol)
  real(r8) :: mwo3  ! Molecular weight of O3 (g/mol)

  real(r8) :: gravit     ! acceleration due to gravity (m/s**2)
  real(r8) :: gravit_cgs ! acceleration due to gravity (cm/s**2)
  real(r8) :: rga        ! 1./gravit_cgs
  real(r8) :: epsilo     ! Ratio of mol. wght of H2O to dry air
  real(r8) :: omeps      ! 1._r8 - epsilo
  real(r8) :: sslp       ! Standard sea-level pressure (dynes/cm**2)
  real(r8) :: stebol_cgs ! Stefan-Boltzmann's constant (CGS)
  real(r8) :: rgsslp     ! 0.5/(gravit_cgs*sslp)
  real(r8) :: dpfo3      ! Voigt correction factor for O3
  real(r8) :: dpfco2     ! Voigt correction factor for CO2

  integer, parameter :: n_u = 25   ! Number of U in abs/emis tables
  integer, parameter :: n_p = 10   ! Number of P in abs/emis tables
  integer, parameter :: n_tp = 10  ! Number of T_p in abs/emis tables
  integer, parameter :: n_te = 21  ! Number of T_e in abs/emis tables
  integer, parameter :: n_rh = 7   ! Number of RH in abs/emis tables

  real(r8):: ah2onw(n_p, n_tp, n_u, n_te, n_rh)   ! absorptivity (non-window)
  real(r8):: eh2onw(n_p, n_tp, n_u, n_te, n_rh)   ! emissivity   (non-window)
  real(r8):: ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! absorptivity (window, for adjacent layers)
  real(r8):: cn_ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for absorptivity (window)
  real(r8):: cn_eh2ow(n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for emissivity   (window)
  real(r8):: ln_ah2ow(n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for absorptivity (window)
  real(r8):: ln_eh2ow(n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for emissivity   (window)
!
! Constant coefficients for water vapor overlap with trace gases.
! Reference: Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
  real(r8):: coefh(2,4) = reshape(  &
         (/ (/5.46557e+01_r8,-7.30387e-02_r8/), &
            (/1.09311e+02_r8,-1.46077e-01_r8/), &
            (/5.11479e+01_r8,-6.82615e-02_r8/), &
            (/1.02296e+02_r8,-1.36523e-01_r8/) /), (/2,4/) )
!
  real(r8):: coefj(3,2) = reshape( &
            (/ (/2.82096e-02_r8,2.47836e-04_r8,1.16904e-06_r8/), &
               (/9.27379e-02_r8,8.04454e-04_r8,6.88844e-06_r8/) /), (/3,2/) )
!
  real(r8):: coefk(3,2) = reshape( &
            (/ (/2.48852e-01_r8,2.09667e-03_r8,2.60377e-06_r8/) , &
               (/1.03594e+00_r8,6.58620e-03_r8,4.04456e-06_r8/) /), (/3,2/) )
  real(r8):: c16,c17,c26,c27,c28,c29,c30,c31
!
! Farwing correction constants for narrow-band emissivity model,
! introduced to account for the deficiencies in narrow-band model
! used to derive the emissivity; tuned with Arkings line-by-line
! calculations.  Just used for water vapor overlap with trace gases.
!
  real(r8):: fwcoef      ! Farwing correction constant
  real(r8):: fwc1,fwc2   ! Farwing correction constants 
  real(r8):: fc1         ! Farwing correction constant 
!
! Collins/Hackney/Edwards (C/H/E) & Collins/Lee-Taylor/Edwards (C/LT/E) 
!       H2O parameterization
!
! Notation:
! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
! absorptivity/emissivity in window are fit using an expression:
!
!      a/e = f_a/e * {1.0 - ln_a/e * cn_a/e} 
!
! absorptivity/emissivity in non-window are fit using:
! 
!      a/e = f_a/e * a/e_norm
!
! where
!      a/e = absorptivity/emissivity
! a/e_norm = absorptivity/emissivity normalized to 1
!    f_a/e = value of a/e as U->infinity = f(T_e) only
!   cn_a/e = continuum transmission
!   ln_a/e = line transmission
!
! spectral interval:
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1 (rotation and rotation-vibration)
!   2 = 800-1200 cm^-1                  (window)
!
  real(r8), parameter:: min_tp_h2o = 160.0_r8        ! min T_p for pre-calculated abs/emis 
  real(r8), parameter:: max_tp_h2o = 349.999999_r8   ! max T_p for pre-calculated abs/emis
  integer, parameter :: ntemp = 192 ! Number of temperatures in H2O sat. table for Tp
  integer, parameter :: o_fa = 6   ! Degree+1 of poly of T_e for absorptivity as U->inf.
  integer, parameter :: o_fe = 6   ! Degree+1 of poly of T_e for emissivity as U->inf.
!-----------------------------------------------------------------------------
! Data for f in C/H/E fit -- value of A and E as U->infinity
! New C/LT/E fit (Hitran 2K, CKD 2.4) -- no change
!     These values are determined by integrals of Planck functions or
!     derivatives of Planck functions only.
!-----------------------------------------------------------------------------
!
! fa/fe coefficients for 2 bands (0-800 & 1200-2200, 800-1200 cm^-1)
!
! Coefficients of polynomial for f_a in T_e
!
  real(r8), parameter:: fat(o_fa,nbands) = reshape( (/ &
       (/-1.06665373E-01_r8,  2.90617375E-02_r8, -2.70642049E-04_r8,   &   ! 0-800&1200-2200 cm^-1
          1.07595511E-06_r8, -1.97419681E-09_r8,  1.37763374E-12_r8/), &   !   0-800&1200-2200 cm^-1
       (/ 1.10666537E+00_r8, -2.90617375E-02_r8,  2.70642049E-04_r8,   &   ! 800-1200 cm^-1
         -1.07595511E-06_r8,  1.97419681E-09_r8, -1.37763374E-12_r8/) /) & !   800-1200 cm^-1
       , (/o_fa,nbands/) )
!
! Coefficients of polynomial for f_e in T_e
!
  real(r8), parameter:: fet(o_fe,nbands) = reshape( (/ & 
      (/3.46148163E-01_r8,  1.51240299E-02_r8, -1.21846479E-04_r8,   &   ! 0-800&1200-2200 cm^-1
        4.04970123E-07_r8, -6.15368936E-10_r8,  3.52415071E-13_r8/), &   !   0-800&1200-2200 cm^-1
      (/6.53851837E-01_r8, -1.51240299E-02_r8,  1.21846479E-04_r8,   &   ! 800-1200 cm^-1
       -4.04970123E-07_r8,  6.15368936E-10_r8, -3.52415071E-13_r8/) /) & !   800-1200 cm^-1
      , (/o_fa,nbands/) )
!
! Note: max values should be slightly underestimated to avoid index bound violations
!
  real(r8), parameter:: min_lp_h2o = -3.0_r8         ! min log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: min_p_h2o = 1.0e-3_r8        ! min log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: max_lp_h2o = -0.0000001_r8   ! max log_10(P) for pre-calculated abs/emis 
  real(r8), parameter:: dlp_h2o = 0.3333333333333_r8 ! difference in adjacent elements of lp_h2o
 
  real(r8), parameter:: dtp_h2o = 21.111111111111_r8 ! difference in adjacent elements of tp_h2o

  real(r8), parameter:: min_rh_h2o = 0.0_r8          ! min RH for pre-calculated abs/emis 
  real(r8), parameter:: max_rh_h2o = 1.19999999_r8   ! max RH for pre-calculated abs/emis 
  real(r8), parameter:: drh_h2o = 0.2_r8             ! difference in adjacent elements of RH

  real(r8), parameter:: min_te_h2o = -120.0_r8       ! min T_e-T_p for pre-calculated abs/emis 
  real(r8), parameter:: max_te_h2o = 79.999999_r8    ! max T_e-T_p for pre-calculated abs/emis 
  real(r8), parameter:: dte_h2o  = 10.0_r8           ! difference in adjacent elements of te_h2o

  real(r8), parameter:: min_lu_h2o = -8.0_r8         ! min log_10(U) for pre-calculated abs/emis 
  real(r8), parameter:: min_u_h2o  = 1.0e-8_r8       ! min pressure-weighted path-length
  real(r8), parameter:: max_lu_h2o =  3.9999999_r8   ! max log_10(U) for pre-calculated abs/emis 
  real(r8), parameter:: dlu_h2o  = 0.5_r8            ! difference in adjacent elements of lu_h2o

  real(r8), parameter:: g1(6)=(/0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/)
  real(r8), parameter :: g2(6)=(/14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/)
  real(r8), parameter :: g3(6)=(/26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/)
  real(r8), parameter :: g4(6)=(/0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/)
  real(r8), parameter :: ab(6)=(/3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/)
  real(r8), parameter :: bb(6)=(/-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/)
  real(r8), parameter :: abp(6)=(/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/)
  real(r8), parameter :: bbp(6)=(/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/)


! Public Interfaces
!====================================================================================
CONTAINS
!====================================================================================

subroutine radabs(lchnk   ,ncol    ,             &
   pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  , &
   s2c    ,tcg     ,w        ,h2otr   ,plco2   , &
   plh2o  ,co2t    ,tint     ,tlayr   ,plol    , &
   plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  , &
   un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  , &
   uco213 ,uco221  ,uco222   ,uco223  ,uptype  , &
   bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 , &
   abstot ,absnxt  ,plh2ob  ,wb       , &
   odap_aer ,aer_trn_ttl, co2mmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
! 
! Method: 
! h2o  ....  Uses nonisothermal emissivity method for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
!            Implementation updated by Collins, Hackney, and Edwards (2001)
!               using line-by-line calculations based upon Hitran 1996 and
!               CKD 2.1 for absorptivity and emissivity
!
!            Implementation updated by Collins, Lee-Taylor, and Edwards (2003)
!               using line-by-line calculations based upon Hitran 2000 and
!               CKD 2.4 for absorptivity and emissivity
!
! co2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019.
!            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
!            are also included.
!
! o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R.Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of CFC11
!
! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of CFC12
!
!
! Computes individual absorptivities for non-adjacent layers, accounting
! for band overlap, and sums to obtain the total; then, computes the
! nearest layer contribution.
! 
! Author: W. Collins (H2O absorptivity) and J. Kiehl
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns

   real(r8), intent(in) :: pbr(pcols,pver)            ! Prssr at mid-levels (dynes/cm2)
   real(r8), intent(in) :: pnm(pcols,pverp)           ! Prssr at interfaces (dynes/cm2)
   real(r8), intent(in) :: co2em(pcols,pverp)         ! Co2 emissivity function
   real(r8), intent(in) :: co2eml(pcols,pver)         ! Co2 emissivity function
   real(r8), intent(in) :: tplnka(pcols,pverp)        ! Planck fnctn level temperature
   real(r8), intent(in) :: s2c(pcols,pverp)           ! H2o continuum path length
   real(r8), intent(in) :: tcg(pcols,pverp)           ! H2o-mass-wgted temp. (Curtis-Godson approx.)
   real(r8), intent(in) :: w(pcols,pverp)             ! H2o prs wghted path
   real(r8), intent(in) :: h2otr(pcols,pverp)         ! H2o trnsmssn fnct for o3 overlap
   real(r8), intent(in) :: plco2(pcols,pverp)         ! Co2 prs wghted path length
   real(r8), intent(in) :: plh2o(pcols,pverp)         ! H2o prs wfhted path length
   real(r8), intent(in) :: co2t(pcols,pverp)          ! Tmp and prs wghted path length
   real(r8), intent(in) :: tint(pcols,pverp)          ! Interface temperatures
   real(r8), intent(in) :: tlayr(pcols,pverp)         ! K-1 level temperatures
   real(r8), intent(in) :: plol(pcols,pverp)          ! Ozone prs wghted path length
   real(r8), intent(in) :: plos(pcols,pverp)          ! Ozone path length
   real(r8), intent(in) :: pmln(pcols,pver)           ! Ln(pmidm1)
   real(r8), intent(in) :: piln(pcols,pverp)          ! Ln(pintm1)
   real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! Pressure weighted h2o path with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 
   real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! H2o path length with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 

! [fraction] absorbtion optical depth, cumulative from top
   real(r8), intent(in) :: odap_aer(pcols,pver,nlwbands)

! [fraction] Total transmission between interfaces k1 and k2
   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) 

!
! Trace gas variables
!
   real(r8), intent(in) :: co2mmr(pcols)              ! co2 column mean mass mixing ratio
   real(r8), intent(in) :: ucfc11(pcols,pverp)        ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)        ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)         ! N2O path length
   real(r8), intent(in) :: un2o1(pcols,pverp)         ! N2O path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)          ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)        ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)        ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)        ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)        ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)        ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)        ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uptype(pcols,pverp)        ! continuum path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)         ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)         ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)          ! pressure factor for ch4
   real(r8), intent(in) :: abplnk1(14,pcols,pverp)    ! non-nearest layer Planck factor
   real(r8), intent(in) :: abplnk2(14,pcols,pverp)    ! nearest layer factor
!
! Output arguments
!
   real(r8), intent(out) :: abstot(pcols,ntoplw:pverp,ntoplw:pverp) ! Total absorptivity
   real(r8), intent(out) :: absnxt(pcols,pver,4)      ! Total nearest layer absorptivity
!
!---------------------------Local variables-----------------------------
!
   integer i                   ! Longitude index
   integer k                   ! Level index
   integer k1                  ! Level index
   integer k2                  ! Level index
   integer kn                  ! Nearest level index
   integer wvl                 ! Wavelength index

   real(r8) abstrc(pcols)              ! total trace gas absorptivity
   real(r8) bplnk(14,pcols,4)          ! Planck functions for sub-divided layers
   real(r8) pnew(pcols)        ! Effective pressure for H2O vapor linewidth
   real(r8) pnewb(nbands)      ! Effective pressure for h2o linewidth w/
                               !    Hulst-Curtis-Godson correction for
                               !    each band
   real(r8) u(pcols)           ! Pressure weighted H2O path length
   real(r8) ub(nbands)         ! Pressure weighted H2O path length with
                               !    Hulst-Curtis-Godson correction for
                               !    each band
   real(r8) tbar(pcols,4)      ! Mean layer temperature
   real(r8) emm(pcols,4)       ! Mean co2 emissivity
   real(r8) o3emm(pcols,4)     ! Mean o3 emissivity
   real(r8) o3bndi             ! Ozone band parameter
   real(r8) temh2o(pcols,4)    ! Mean layer temperature equivalent to tbar
   real(r8) k21                ! Exponential coefficient used to calculate
!                              !  rotation band transmissvty in the 650-800
!                              !  cm-1 region (tr1)
   real(r8) k22                ! Exponential coefficient used to calculate
!                              !  rotation band transmissvty in the 500-650
!                              !  cm-1 region (tr2)
   real(r8) uc1(pcols)         ! H2o continuum pathlength in 500-800 cm-1
   real(r8) to3h2o(pcols)      ! H2o trnsmsn for overlap with o3
   real(r8) pi                 ! For co2 absorptivity computation
   real(r8) sqti(pcols)        ! Used to store sqrt of mean temperature
   real(r8) et                 ! Co2 hot band factor
   real(r8) et2                ! Co2 hot band factor squared
   real(r8) et4                ! Co2 hot band factor to fourth power
   real(r8) omet               ! Co2 stimulated emission term
   real(r8) f1co2              ! Co2 central band factor
   real(r8) f2co2(pcols)       ! Co2 weak band factor
   real(r8) f3co2(pcols)       ! Co2 weak band factor
   real(r8) t1co2(pcols)       ! Overlap factr weak bands on strong band
   real(r8) sqwp               ! Sqrt of co2 pathlength
   real(r8) f1sqwp(pcols)      ! Main co2 band factor
   real(r8) oneme              ! Co2 stimulated emission term
   real(r8) alphat             ! Part of the co2 stimulated emission term
   real(r8) co2vmr(pcols)      ! CO2 column mean vmr
   real(r8) rmw                ! ratio of molecular weights (air/co2)
   real(r8) wco2               ! Constants used to define co2 pathlength
   real(r8) posqt              ! Effective pressure for co2 line width
   real(r8) u7(pcols)          ! Co2 hot band path length
   real(r8) u8                 ! Co2 hot band path length
   real(r8) u9                 ! Co2 hot band path length
   real(r8) u13                ! Co2 hot band path length
   real(r8) rbeta7(pcols)      ! Inverse of co2 hot band line width par
   real(r8) rbeta8             ! Inverse of co2 hot band line width par
   real(r8) rbeta9             ! Inverse of co2 hot band line width par
   real(r8) rbeta13            ! Inverse of co2 hot band line width par
   real(r8) tpatha             ! For absorptivity computation
   real(r8) abso(pcols,4)      ! Absorptivity for various gases/bands
   real(r8) dtx(pcols)         ! Planck temperature minus 250 K
   real(r8) dty(pcols)         ! Path temperature minus 250 K
   real(r8) term7(pcols,2)     ! Kl_inf(i) in eq(r8) of table A3a of R&D
   real(r8) term8(pcols,2)     ! Delta kl_inf(i) in eq(r8)
   real(r8) tr1                ! Eqn(6) in table A2 of R&D for 650-800
   real(r8) tr10(pcols)        ! Eqn (6) times eq(4) in table A2
!                              !  of R&D for 500-650 cm-1 region
   real(r8) tr2                ! Eqn(6) in table A2 of R&D for 500-650
   real(r8) tr5                ! Eqn(4) in table A2 of R&D for 650-800
   real(r8) tr6                ! Eqn(4) in table A2 of R&D for 500-650
   real(r8) tr9(pcols)         ! Equation (6) times eq(4) in table A2
!                              !  of R&D for 650-800 cm-1 region
   real(r8) sqrtu(pcols)       ! Sqrt of pressure weighted h20 pathlength
   real(r8) fwk(pcols)         ! Equation(33) in R&D far wing correction
   real(r8) fwku(pcols)        ! GU term in eqs(1) and (6) in table A2
   real(r8) to3co2(pcols)      ! P weighted temp in ozone band model
   real(r8) dpnm(pcols)        ! Pressure difference between two levels
   real(r8) pnmsq(pcols,pverp) ! Pressure squared
   real(r8) dw(pcols)          ! Amount of h2o between two levels
   real(r8) uinpl(pcols,4)     ! Nearest layer subdivision factor
   real(r8) winpl(pcols,4)     ! Nearest layer subdivision factor
   real(r8) zinpl(pcols,4)     ! Nearest layer subdivision factor
   real(r8) pinpl(pcols,4)     ! Nearest layer subdivision factor
   real(r8) dplh2o(pcols)      ! Difference in press weighted h2o amount
   real(r8) r293               ! 1/293
   real(r8) r250               ! 1/250
   real(r8) r3205              ! Line width factor for o3 (see R&Di)
   real(r8) r300               ! 1/300
   real(r8) rsslp              ! Reciprocal of sea level pressure
   real(r8) r2sslp             ! 1/2 of rsslp
   real(r8) ds2c               ! Y in eq(7) in table A2 of R&D
   real(r8)  dplos             ! Ozone pathlength eq(A2) in R&Di
   real(r8) dplol              ! Presure weighted ozone pathlength
   real(r8) tlocal             ! Local interface temperature
   real(r8) beta               ! Ozone mean line parameter eq(A3) in R&Di
!                               (includes Voigt line correction factor)
   real(r8) rphat              ! Effective pressure for ozone beta
   real(r8) tcrfac             ! Ozone temperature factor table 1 R&Di
   real(r8) tmp1               ! Ozone band factor see eq(A1) in R&Di
   real(r8) u1                 ! Effective ozone pathlength eq(A2) in R&Di
   real(r8) realnu             ! 1/beta factor in ozone band model eq(A1)
   real(r8) tmp2               ! Ozone band factor see eq(A1) in R&Di
   real(r8) u2                 ! Effective ozone pathlength eq(A2) in R&Di
   real(r8) rsqti              ! Reciprocal of sqrt of path temperature
   real(r8) tpath              ! Path temperature used in co2 band model
   real(r8) tmp3               ! Weak band factor see K&B
   real(r8) rdpnmsq            ! Reciprocal of difference in press^2
   real(r8) rdpnm              ! Reciprocal of difference in press
   real(r8) p1                 ! Mean pressure factor
   real(r8) p2                 ! Mean pressure factor
   real(r8) dtym10             ! T - 260 used in eq(9) and (10) table A3a
   real(r8) dplco2             ! Co2 path length
   real(r8) te                 ! A_0 T factor in ozone model table 1 of R&Di
   real(r8) denom              ! Denominator in eq(r8) of table A3a of R&D
   real(r8) th2o(pcols)        ! transmission due to H2O
   real(r8) tco2(pcols)        ! transmission due to CO2
   real(r8) to3(pcols)         ! transmission due to O3
!
! Transmission terms for various spectral intervals:
!
   real(r8) trab2(pcols)       ! H2o   500 -  800 cm-1
   real(r8) absbnd             ! Proportional to co2 band absorptance
   real(r8) dbvtit(pcols,pverp)! Intrfc drvtv plnck fnctn for o3
   real(r8) dbvtly(pcols,pver) ! Level drvtv plnck fnctn for o3
!
! Variables for Collins/Hackney/Edwards (C/H/E) & 
!       Collins/Lee-Taylor/Edwards (C/LT/E) H2O parameterization

!
! Notation:
! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
   real(r8) fa               ! asymptotic value of abs. as U->infinity
   real(r8) a_star           ! normalized absorptivity for non-window
   real(r8) l_star           ! interpolated line transmission
   real(r8) c_star           ! interpolated continuum transmission

   real(r8) te1              ! emission temperature
   real(r8) te2              ! te^2
   real(r8) te3              ! te^3
   real(r8) te4              ! te^4
   real(r8) te5              ! te^5

   real(r8) log_u            ! log base 10 of U 
   real(r8) log_uc           ! log base 10 of H2O continuum path
   real(r8) log_p            ! log base 10 of P
   real(r8) t_p              ! T_p
   real(r8) t_e              ! T_e (offset by T_p)

   integer iu                ! index for log10(U)
   integer iu1               ! iu + 1
   integer iuc               ! index for log10(H2O continuum path)
   integer iuc1              ! iuc + 1
   integer ip                ! index for log10(P)
   integer ip1               ! ip + 1
   integer itp               ! index for T_p
   integer itp1              ! itp + 1
   integer ite               ! index for T_e
   integer ite1              ! ite + 1
   integer irh               ! index for RH
   integer irh1              ! irh + 1

   real(r8) dvar             ! normalized variation in T_p/T_e/P/U
   real(r8) uvar             ! U * diffusivity factor
   real(r8) uscl             ! factor for lineary scaling as U->0

   real(r8) wu               ! weight for U
   real(r8) wu1              ! 1 - wu
   real(r8) wuc              ! weight for H2O continuum path
   real(r8) wuc1             ! 1 - wuc
   real(r8) wp               ! weight for P
   real(r8) wp1              ! 1 - wp
   real(r8) wtp              ! weight for T_p
   real(r8) wtp1             ! 1 - wtp
   real(r8) wte              ! weight for T_e
   real(r8) wte1             ! 1 - wte
   real(r8) wrh              ! weight for RH
   real(r8) wrh1             ! 1 - wrh

   real(r8) w_0_0_           ! weight for Tp/Te combination
   real(r8) w_0_1_           ! weight for Tp/Te combination
   real(r8) w_1_0_           ! weight for Tp/Te combination
   real(r8) w_1_1_           ! weight for Tp/Te combination

   real(r8) w_0_00           ! weight for Tp/Te/RH combination
   real(r8) w_0_01           ! weight for Tp/Te/RH combination
   real(r8) w_0_10           ! weight for Tp/Te/RH combination
   real(r8) w_0_11           ! weight for Tp/Te/RH combination
   real(r8) w_1_00           ! weight for Tp/Te/RH combination
   real(r8) w_1_01           ! weight for Tp/Te/RH combination
   real(r8) w_1_10           ! weight for Tp/Te/RH combination
   real(r8) w_1_11           ! weight for Tp/Te/RH combination

   real(r8) w00_00           ! weight for P/Tp/Te/RH combination
   real(r8) w00_01           ! weight for P/Tp/Te/RH combination
   real(r8) w00_10           ! weight for P/Tp/Te/RH combination
   real(r8) w00_11           ! weight for P/Tp/Te/RH combination
   real(r8) w01_00           ! weight for P/Tp/Te/RH combination
   real(r8) w01_01           ! weight for P/Tp/Te/RH combination
   real(r8) w01_10           ! weight for P/Tp/Te/RH combination
   real(r8) w01_11           ! weight for P/Tp/Te/RH combination
   real(r8) w10_00           ! weight for P/Tp/Te/RH combination
   real(r8) w10_01           ! weight for P/Tp/Te/RH combination
   real(r8) w10_10           ! weight for P/Tp/Te/RH combination
   real(r8) w10_11           ! weight for P/Tp/Te/RH combination
   real(r8) w11_00           ! weight for P/Tp/Te/RH combination
   real(r8) w11_01           ! weight for P/Tp/Te/RH combination
   real(r8) w11_10           ! weight for P/Tp/Te/RH combination
   real(r8) w11_11           ! weight for P/Tp/Te/RH combination

   integer ib                ! spectral interval:
                             !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
                             !   2 = 800-1200 cm^-1


   real(r8) pch2o            ! H2O continuum path
   real(r8) fch2o            ! temp. factor for continuum
   real(r8) uch2o            ! U corresponding to H2O cont. path (window)

   real(r8) fdif             ! secant(zenith angle) for diffusivity approx.

   real(r8) sslp_mks         ! Sea-level pressure in MKS units
   real(r8) esx              ! saturation vapor pressure returned by qsat
   real(r8) qsx              ! saturation mixing ratio returned by qsat
   real(r8) pnew_mks         ! pnew in MKS units
   real(r8) q_path           ! effective specific humidity along path
   real(r8) rh_path          ! effective relative humidity along path

      integer bnd_idx        ! LW band index
      real(r8) aer_pth_dlt   ! [kg m-2] STRAER path between interface levels k1 and k2
      real(r8) aer_pth_ngh(pcols)
                             ! [kg m-2] STRAER path between neighboring layers
      real(r8) odap_aer_ttl  ! [fraction] Total path absorption optical depth
      real(r8) aer_trn_ngh(pcols,nlwbands) 
                             ! [fraction] Total transmission between 
                             !            nearest neighbor sub-levels
!
!--------------------------Statement function---------------------------
!
   real(r8) dbvt,t             ! Planck fnctn tmp derivative for o3
!
   dbvt(t)=(-2.8911366682e-4_r8+(2.3771251896e-6_r8+1.1305188929e-10_r8*t)*t)/ &
      (1.0_r8+(-6.1364820707e-3_r8+1.5550319767e-5_r8*t)*t)
!
!
!-----------------------------------------------------------------------
!
! Initialize
!
   do k2=1,4
      do k1=1,ntoplw-1
         absnxt(:,k1,k2) = posinf    ! set unused portions for lf95 restart write
      end do
   end do

   do k=ntoplw,pverp
      abstot(:,k,k) = posinf         ! set unused portions for lf95 restart write
   end do

   do k=ntoplw,pver
      do i=1,ncol
         dbvtly(i,k) = dbvt(tlayr(i,k+1))
         dbvtit(i,k) = dbvt(tint(i,k))
      end do
   end do
   rmw = amd/amco2
   do i=1,ncol
      dbvtit(i,pverp) = dbvt(tint(i,pverp))
      co2vmr(i) = co2mmr(i) * rmw
   end do
!
   r293    = 1._r8/293._r8
   r250    = 1._r8/250._r8
   r3205   = 1._r8/.3205_r8
   r300    = 1._r8/300._r8
   rsslp   = 1._r8/sslp
   r2sslp  = 1._r8/(2._r8*sslp)
!
!Constants for computing U corresponding to H2O cont. path
!
   fdif       = 1.66_r8
   sslp_mks   = sslp / 10.0_r8
!
! Non-adjacent layer absorptivity:
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O absorptivity
!
! 500-800 cm^-1 H2o continuum/line overlap already included
!                in abso(i,1).  This used to be in abso(i,4)
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! abso(i,4)   co2 15  micrometer band system
!
   do k=ntoplw,pverp
      do i=1,ncol
         pnmsq(i,k) = pnm(i,k)**2
         dtx(i) = tplnka(i,k) - 250._r8
      end do
   end do
!
! Non-nearest layer level loops
!
   do k1=pverp,ntoplw,-1
      do k2=pverp,ntoplw,-1
         if (k1 == k2) cycle
         do i=1,ncol
            dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
            u(i)      = abs(dplh2o(i))
            sqrtu(i)  = sqrt(u(i))
            ds2c      = abs(s2c(i,k1) - s2c(i,k2))
            dw(i)     = abs(w(i,k1) - w(i,k2))
            uc1(i)    = (ds2c + 1.7e-3_r8*u(i))*(1._r8 +  2._r8*ds2c)/(1._r8 + 15._r8*ds2c)
            pch2o     = ds2c
            pnew(i)   = u(i)/dw(i)
            pnew_mks  = pnew(i) * sslp_mks
!
! Changed effective path temperature to std. Curtis-Godson form
!
            tpatha = abs(tcg(i,k1) - tcg(i,k2))/dw(i)
            t_p = min(max(tpatha, min_tp_h2o), max_tp_h2o)

            call qsat_water(t_p, pnew_mks, esx, qsx)
!
! Compute effective RH along path
!
            q_path = dw(i) / abs(pnm(i,k1) - pnm(i,k2)) / rga
!
! Calculate effective u, pnew for each band using
!        Hulst-Curtis-Godson approximation:
! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
!           2nd edition, Oxford University Press, 1989.
! Effective H2O path (w)
!      eq. 6.24, p. 228
! Effective H2O path pressure (pnew = u/w):
!      eq. 6.29, p. 228
!
            ub(1) = abs(plh2ob(1,i,k1) - plh2ob(1,i,k2)) / psi(t_p,1)
            ub(2) = abs(plh2ob(2,i,k1) - plh2ob(2,i,k2)) / psi(t_p,2)
            
            pnewb(1) = ub(1) / abs(wb(1,i,k1) - wb(1,i,k2)) * phi(t_p,1)
            pnewb(2) = ub(2) / abs(wb(2,i,k1) - wb(2,i,k2)) * phi(t_p,2)

            dtx(i)      = tplnka(i,k2) - 250._r8
            dty(i)      = tpatha       - 250._r8

            fwk(i)  = fwcoef + fwc1/(1._r8 + fwc2*u(i))
            fwku(i) = fwk(i)*u(i)
!
! Define variables for C/H/E (now C/LT/E) fit
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O absorptivity
!
! Notation:
! U   = integral (P/P_0 dW)  
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
!
! Terms for asymptotic value of emissivity
!
            te1  = tplnka(i,k2)
            te2  = te1 * te1
            te3  = te2 * te1
            te4  = te3 * te1
            te5  = te4 * te1

!
!  Band-independent indices for lines and continuum tables
!
            dvar = (t_p - min_tp_h2o) / dtp_h2o
            itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
            itp1 = itp + 1
            wtp = dvar - floor(dvar)
            wtp1 = 1.0_r8 - wtp
            
            t_e = min(max(tplnka(i,k2)-t_p, min_te_h2o), max_te_h2o)
            dvar = (t_e - min_te_h2o) / dte_h2o
            ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
            ite1 = ite + 1
            wte = dvar - floor(dvar)
            wte1 = 1.0_r8 - wte
            
            rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
            dvar = (rh_path - min_rh_h2o) / drh_h2o
            irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
            irh1 = irh + 1
            wrh = dvar - floor(dvar)
            wrh1 = 1.0_r8 - wrh

            w_0_0_ = wtp  * wte
            w_0_1_ = wtp  * wte1
            w_1_0_ = wtp1 * wte 
            w_1_1_ = wtp1 * wte1
            
            w_0_00 = w_0_0_ * wrh
            w_0_01 = w_0_0_ * wrh1
            w_0_10 = w_0_1_ * wrh
            w_0_11 = w_0_1_ * wrh1
            w_1_00 = w_1_0_ * wrh
            w_1_01 = w_1_0_ * wrh1
            w_1_10 = w_1_1_ * wrh
            w_1_11 = w_1_1_ * wrh1

!
! H2O Continuum path for 0-800 and 1200-2200 cm^-1
!
!    Assume foreign continuum dominates total H2O continuum in these bands
!    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
!    Then the effective H2O path is just 
!         U_c = integral[ f(P) dW ]
!    where 
!           W = water-vapor mass and 
!        f(P) = dependence of foreign continuum on pressure 
!             = P / sslp
!    Then 
!         U_c = U (the same effective H2O path as for lines)
!
!
! Continuum terms for 800-1200 cm^-1
!
!    Assume self continuum dominates total H2O continuum for this band
!    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
!    Then the effective H2O self-continuum path is 
!         U_c = integral[ h(e,T) dW ]                        (*eq. 1*)
!    where 
!           W = water-vapor mass and 
!           e = partial pressure of H2O along path
!           T = temperature along path
!      h(e,T) = dependence of foreign continuum on e,T
!             = e / sslp * f(T)
!
!    Replacing
!           e =~ q * P / epsilo
!           q = mixing ratio of H2O
!     epsilo = 0.622
!
!    and using the definition
!           U = integral [ (P / sslp) dW ]
!             = (P / sslp) W                                 (homogeneous path)
!
!    the effective path length for the self continuum is
!         U_c = (q / epsilo) f(T) U                         (*eq. 2*)
!
!    Once values of T, U, and q have been calculated for the inhomogeneous
!        path, this sets U_c for the corresponding
!        homogeneous atmosphere.  However, this need not equal the
!        value of U_c' defined by eq. 1 for the actual inhomogeneous atmosphere
!        under consideration.
!
!    Solution: hold T and q constant, solve for U' that gives U_c' by
!        inverting eq. (2):
!
!        U' = (U_c * epsilo) / (q * f(T))
!
            fch2o = fh2oself(t_p) 
            uch2o = (pch2o * epsilo) / (q_path * fch2o)

!
! Band-dependent indices for non-window
!
            ib = 1

            uvar = ub(ib) * fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0_r8 - wu
            
            log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0_r8 - wp
         
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 
!
! Asymptotic value of absorptivity as U->infinity
!
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5

            a_star = &
                 ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
            abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
                                 aer_trn_ttl(i,k1,k2,ib)), &
                             0.0_r8), 1.0_r8)
!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
                         
!
! Band-dependent indices for window
!
            ib = 2

            uvar = ub(ib) * fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0_r8 - wu
            
            log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0_r8 - wp
         
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 

            log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
            dvar = (log_uc - min_lu_h2o) / dlu_h2o
            iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iuc1 = iuc + 1
            wuc = dvar - floor(dvar)
            wuc1 = 1.0_r8 - wuc
!
! Asymptotic value of absorptivity as U->infinity
!
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5

            l_star = &
                 ln_ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ln_ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ln_ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ln_ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ln_ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ln_ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ln_ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 

            c_star = &
                 cn_ah2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
            abso(i,ib) = min(max(fa * (1.0_r8 - l_star * c_star * &
                                 aer_trn_ttl(i,k1,k2,ib)), &
                             0.0_r8), 1.0_r8) 
!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif

         end do
!
! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
         do i=1,ncol
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8 + c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8 + c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8 + c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8 + c27*dty(i))
         end do
!
! 500 -  800 cm-1   h2o rotation band overlap with co2
!
         do i=1,ncol
            k21    = term7(i,1) + term8(i,1)/ &
               (1._r8 + (c30 + c31*(dty(i)-10._r8)*(dty(i)-10._r8))*sqrtu(i))
            k22    = term7(i,2) + term8(i,2)/ &
               (1._r8 + (c28 + c29*(dty(i)-10._r8))*sqrtu(i))
            tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
            tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
            tr1=tr1*aer_trn_ttl(i,k1,k2,idx_LW_0650_0800) 
!                                          ! H2O line+STRAER trn 650--800 cm-1
            tr2=tr2*aer_trn_ttl(i,k1,k2,idx_LW_0500_0650)
!                                          ! H2O line+STRAER trn 500--650 cm-1
            tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
            tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i)   = tr1*tr5
            tr10(i)  = tr2*tr6
            th2o(i) = tr10(i)
            trab2(i) = 0.65_r8*tr9(i) + 0.35_r8*tr10(i)
         end do
         if (k2 < k1) then
            do i=1,ncol
               to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
            end do
         else
            do i=1,ncol
               to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
            end do
         end if
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
         do i=1,ncol
            dpnm(i)  = pnm(i,k1) - pnm(i,k2)
            to3co2(i) = (pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/dpnm(i)
            te       = (to3co2(i)*r293)**.7_r8
            dplos    = plos(i,k1) - plos(i,k2)
            if (dplos == 0._r8) then
              abso(i,3) = 0._r8
              to3(i) = 1._r8
              write(iulog,*) 'radiation ozone error  ',k1,k2,plos(i,k1)
            else 
              dplol    = plol(i,k1) - plol(i,k2)
              u1       = 18.29_r8*abs(dplos)/te
              u2       = .5649_r8*abs(dplos)/te
              rphat    = dplol/dplos
              tlocal   = tint(i,k2)
              tcrfac   = sqrt(tlocal*r250)*te
              beta     = r3205*(rphat + dpfo3*tcrfac)
              realnu   = te/beta
              tmp1     = u1/sqrt(4._r8 + u1*(1._r8 + realnu))
              tmp2     = u2/sqrt(4._r8 + u2*(1._r8 + realnu))
              o3bndi    = 74._r8*te*log(1._r8 + tmp1 + tmp2)
              abso(i,3) = o3bndi*to3h2o(i)*dbvtit(i,k2)
              to3(i)   = 1.0_r8/(1._r8 + 0.1_r8*tmp1 + 0.1_r8*tmp2)
            endif
         end do
!
! abso(i,4)      co2 15  micrometer band system
!
         do i=1,ncol
            sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
            et        = exp(-480._r8/to3co2(i))
            sqti(i)   = sqrt(to3co2(i))
            rsqti     = 1._r8/sqti(i)
            et2       = et*et
            et4       = et2*et2
            omet      = 1._r8 - 1.5_r8*et2
            f1co2     = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
            f1sqwp(i) = f1co2*sqwp
            t1co2(i)  = 1._r8/(1._r8 + (245.18_r8*omet*sqwp*rsqti))
            oneme     = 1._r8 - et2
            alphat    = oneme**3*rsqti
            pi        = abs(dpnm(i))
            wco2      =  2.5221_r8*co2vmr(i)*pi*rga
            u7(i)     =  4.9411e4_r8*alphat*et2*wco2
            u8        =  3.9744e4_r8*alphat*et4*wco2
            u9        =  1.0447e5_r8*alphat*et4*et2*wco2
            u13       = 2.8388e3_r8*alphat*et4*wco2
            tpath     = to3co2(i)
            tlocal    = tint(i,k2)
            tcrfac    = sqrt(tlocal*r250*tpath*r300)
            posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp + dpfco2*tcrfac)*rsqti
            rbeta7(i) = 1._r8/(5.3228_r8*posqt)
            rbeta8    = 1._r8/(10.6576_r8*posqt)
            rbeta9    = rbeta7(i)
            rbeta13   = rbeta9
            f2co2(i)  = (u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i)))) + &
               (u8   /sqrt(4._r8 + u8*(1._r8 + rbeta8))) + &
               (u9   /sqrt(4._r8 + u9*(1._r8 + rbeta9)))
            f3co2(i)  = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
         end do
         if (k2 >= k1) then
            do i=1,ncol
               sqti(i) = sqrt(tlayr(i,k2))
            end do
         end if
!
         do i=1,ncol
            tmp1      = log(1._r8 + f1sqwp(i))
            tmp2      = log(1._r8 + f2co2(i))
            tmp3      = log(1._r8 + f3co2(i))
            absbnd    = (tmp1 + 2._r8*t1co2(i)*tmp2 + 2._r8*tmp3)*sqti(i)
            abso(i,4) = trab2(i)*co2em(i,k2)*absbnd
            tco2(i)   = 1._r8/(1.0_r8+10.0_r8*(u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i)))))
         end do
!
! Calculate absorptivity due to trace gases, abstrc
!
         call trcab(ncol     ,                            &
            k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
            un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
            uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
            bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
            s2c     ,uptype  ,u       ,abplnk1 ,tco2    , &
            th2o    ,to3     ,abstrc  , &
            aer_trn_ttl)
!
! Sum total absorptivity
!
         do i=1,ncol
            abstot(i,k1,k2) = abso(i,1) + abso(i,2) + &
               abso(i,3) + abso(i,4) + abstrc(i)
         end do
      end do ! do k2 = 
   end do ! do k1 = 
!
! Adjacent layer absorptivity:
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O absorptivity
!
! 500-800 cm^-1 H2o continuum/line overlap already included
!                in abso(i,1).  This used to be in abso(i,4)
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! abso(i,4)   co2 15  micrometer band system
!
! Nearest layer level loop
!
   do k2=pver,ntoplw,-1
      do i=1,ncol
         tbar(i,1)   = 0.5_r8*(tint(i,k2+1) + tlayr(i,k2+1))
         emm(i,1)    = 0.5_r8*(co2em(i,k2+1) + co2eml(i,k2))
         tbar(i,2)   = 0.5_r8*(tlayr(i,k2+1) + tint(i,k2))
         emm(i,2)    = 0.5_r8*(co2em(i,k2) + co2eml(i,k2))
         tbar(i,3)   = 0.5_r8*(tbar(i,2) + tbar(i,1))
         emm(i,3)    = emm(i,1)
         tbar(i,4)   = tbar(i,3)
         emm(i,4)    = emm(i,2)
         o3emm(i,1)  = 0.5_r8*(dbvtit(i,k2+1) + dbvtly(i,k2))
         o3emm(i,2)  = 0.5_r8*(dbvtit(i,k2) + dbvtly(i,k2))
         o3emm(i,3)  = o3emm(i,1)
         o3emm(i,4)  = o3emm(i,2)
         temh2o(i,1) = tbar(i,1)
         temh2o(i,2) = tbar(i,2)
         temh2o(i,3) = tbar(i,1)
         temh2o(i,4) = tbar(i,2)
         dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
      end do
!
!  Weighted Planck functions for trace gases
!
      do wvl = 1,14
         do i = 1,ncol
            bplnk(wvl,i,1) = 0.5_r8*(abplnk1(wvl,i,k2+1) + abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5_r8*(abplnk1(wvl,i,k2) + abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
         end do
      end do
      
      do i=1,ncol
         rdpnmsq    = 1._r8/(pnmsq(i,k2+1) - pnmsq(i,k2))
         rdpnm      = 1._r8/dpnm(i)
         p1         = .5_r8*(pbr(i,k2) + pnm(i,k2+1))
         p2         = .5_r8*(pbr(i,k2) + pnm(i,k2  ))
         uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
         uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
         uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
         uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
         winpl(i,1) = (.5_r8*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
         winpl(i,2) = (.5_r8*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
         winpl(i,3) = (.5_r8*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*rdpnm
         winpl(i,4) = (.5_r8*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*rdpnm
         tmp1       = 1._r8/(piln(i,k2+1) - piln(i,k2))
         tmp2       = piln(i,k2+1) - pmln(i,k2)
         tmp3       = piln(i,k2  ) - pmln(i,k2)
         zinpl(i,1) = (.5_r8*tmp2          )*tmp1
         zinpl(i,2) = (        - .5_r8*tmp3)*tmp1
         zinpl(i,3) = (.5_r8*tmp2 -    tmp3)*tmp1
         zinpl(i,4) = (   tmp2 - .5_r8*tmp3)*tmp1
         pinpl(i,1) = 0.5_r8*(p1 + pnm(i,k2+1))
         pinpl(i,2) = 0.5_r8*(p2 + pnm(i,k2  ))
         pinpl(i,3) = 0.5_r8*(p1 + pnm(i,k2  ))
         pinpl(i,4) = 0.5_r8*(p2 + pnm(i,k2+1))
      end do
      do kn=1,4
         do i=1,ncol
            u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
            sqrtu(i) = sqrt(u(i))
            dw(i)    = abs(w(i,k2) - w(i,k2+1))
            pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
            pnew_mks  = pnew(i) * sslp_mks
            t_p = min(max(tbar(i,kn), min_tp_h2o), max_tp_h2o)

            call qsat_water(t_p, pnew_mks, esx, qsx)

            q_path = dw(i) / ABS(dpnm(i)) / rga
            
            ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
            uc1(i)   = uinpl(i,kn)*ds2c
            pch2o    = uc1(i)
            uc1(i)   = (uc1(i) + 1.7e-3_r8*u(i))*(1._r8 +  2._r8*uc1(i))/(1._r8 + 15._r8*uc1(i))
            dtx(i)      = temh2o(i,kn) - 250._r8
            dty(i)      = tbar(i,kn) - 250._r8
            
            fwk(i)    = fwcoef + fwc1/(1._r8 + fwc2*u(i))
            fwku(i)   = fwk(i)*u(i)

            aer_trn_ngh(i, 1:nlwbands)= &
              exp(-fdif * uinpl(i,kn) * odap_aer(i, k2, 1:nlwbands ) )

!
! Define variables for C/H/E (now C/LT/E) fit
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O absorptivity
!
! Notation:
! U   = integral (P/P_0 dW)  
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
!
! Terms for asymptotic value of emissivity
!
            te1  = temh2o(i,kn)
            te2  = te1 * te1
            te3  = te2 * te1
            te4  = te3 * te1
            te5  = te4 * te1

!
! Indices for lines and continuum tables 
! Note: because we are dealing with the nearest layer,
!       the Hulst-Curtis-Godson corrections
!       for inhomogeneous paths are not applied.
!
            uvar = u(i)*fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0_r8 - wu
            
            log_p  = min(log10(max(pnew(i), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0_r8 - wp
            
            dvar = (t_p - min_tp_h2o) / dtp_h2o
            itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
            itp1 = itp + 1
            wtp = dvar - floor(dvar)
            wtp1 = 1.0_r8 - wtp
            
            t_e = min(max(temh2o(i,kn)-t_p,min_te_h2o),max_te_h2o)
            dvar = (t_e - min_te_h2o) / dte_h2o
            ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
            ite1 = ite + 1
            wte = dvar - floor(dvar)
            wte1 = 1.0_r8 - wte
            
            rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
            dvar = (rh_path - min_rh_h2o) / drh_h2o
            irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
            irh1 = irh + 1
            wrh = dvar - floor(dvar)
            wrh1 = 1.0_r8 - wrh
            
            w_0_0_ = wtp  * wte
            w_0_1_ = wtp  * wte1
            w_1_0_ = wtp1 * wte 
            w_1_1_ = wtp1 * wte1
            
            w_0_00 = w_0_0_ * wrh
            w_0_01 = w_0_0_ * wrh1
            w_0_10 = w_0_1_ * wrh
            w_0_11 = w_0_1_ * wrh1
            w_1_00 = w_1_0_ * wrh
            w_1_01 = w_1_0_ * wrh1
            w_1_10 = w_1_1_ * wrh
            w_1_11 = w_1_1_ * wrh1
            
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 

!
! Non-window absorptivity
!
            ib = 1
            
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5
            
            a_star = &
                 ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
            
            abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
                                 aer_trn_ngh(i,ib)), &
                             0.0_r8), 1.0_r8)

!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
            
!
! Window absorptivity
!
            ib = 2
            
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5
            
            a_star = &
                 ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
            
            abso(i,ib) = min(max(fa * (1.0_r8 - (1.0_r8 - a_star) * &
                                 aer_trn_ngh(i,ib)), &
                             0.0_r8), 1.0_r8)

!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
            
         end do
!
! Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
         do i=1,ncol
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8 + c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8 + c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8 + c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8 + c27*dty(i))
         end do
!
! 500 -  800 cm-1   h2o rotation band overlap with co2
!
         do i=1,ncol
            dtym10     = dty(i) - 10._r8
            denom      = 1._r8 + (c30 + c31*dtym10*dtym10)*sqrtu(i)
            k21        = term7(i,1) + term8(i,1)/denom
            denom      = 1._r8 + (c28 + c29*dtym10       )*sqrtu(i)
            k22        = term7(i,2) + term8(i,2)/denom
            tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
            tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
            tr1=tr1*aer_trn_ngh(i,idx_LW_0650_0800) 
!                                         ! H2O line+STRAER trn 650--800 cm-1
            tr2=tr2*aer_trn_ngh(i,idx_LW_0500_0650) 
!                                         ! H2O line+STRAER trn 500--650 cm-1
            tr5     = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
            tr6     = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i)  = tr1*tr5
            tr10(i) = tr2*tr6
            trab2(i)= 0.65_r8*tr9(i) + 0.35_r8*tr10(i)
            th2o(i) = tr10(i)
         end do
!
! abso(i,3)  o3  9.6 micrometer (nu3 and nu1 bands)
!
         do i=1,ncol
            te        = (tbar(i,kn)*r293)**.7_r8
            dplos     = abs(plos(i,k2+1) - plos(i,k2))
            u1        = zinpl(i,kn)*18.29_r8*dplos/te
            u2        = zinpl(i,kn)*.5649_r8*dplos/te
            tlocal    = tbar(i,kn)
            tcrfac    = sqrt(tlocal*r250)*te
            beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
            realnu    = te/beta
            tmp1      = u1/sqrt(4._r8 + u1*(1._r8 + realnu))
            tmp2      = u2/sqrt(4._r8 + u2*(1._r8 + realnu))
            o3bndi    = 74._r8*te*log(1._r8 + tmp1 + tmp2)
            abso(i,3) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
            to3(i)    = 1.0_r8/(1._r8 + 0.1_r8*tmp1 + 0.1_r8*tmp2)
         end do
!
! abso(i,4)   co2 15  micrometer band system
!
         do i=1,ncol
            dplco2   = plco2(i,k2+1) - plco2(i,k2)
            sqwp     = sqrt(uinpl(i,kn)*dplco2)
            et       = exp(-480._r8/tbar(i,kn))
            sqti(i)  = sqrt(tbar(i,kn))
            rsqti    = 1._r8/sqti(i)
            et2      = et*et
            et4      = et2*et2
            omet     = (1._r8 - 1.5_r8*et2)
            f1co2    = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
            f1sqwp(i)= f1co2*sqwp
            t1co2(i) = 1._r8/(1._r8 + (245.18_r8*omet*sqwp*rsqti))
            oneme    = 1._r8 - et2
            alphat   = oneme**3*rsqti
            pi       = abs(dpnm(i))*winpl(i,kn)
            wco2     = 2.5221_r8*co2vmr(i)*pi*rga
            u7(i)    = 4.9411e4_r8*alphat*et2*wco2
            u8       = 3.9744e4_r8*alphat*et4*wco2
            u9       = 1.0447e5_r8*alphat*et4*et2*wco2
            u13      = 2.8388e3_r8*alphat*et4*wco2
            tpath    = tbar(i,kn)
            tlocal   = tbar(i,kn)
            tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
            posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
            rbeta7(i)= 1._r8/(5.3228_r8*posqt)
            rbeta8   = 1._r8/(10.6576_r8*posqt)
            rbeta9   = rbeta7(i)
            rbeta13  = rbeta9
            f2co2(i) = u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i))) + &
                 u8   /sqrt(4._r8 + u8*(1._r8 + rbeta8)) + &
                 u9   /sqrt(4._r8 + u9*(1._r8 + rbeta9))
            f3co2(i) = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
            tmp1     = log(1._r8 + f1sqwp(i))
            tmp2     = log(1._r8 + f2co2(i))
            tmp3     = log(1._r8 + f3co2(i))
            absbnd   = (tmp1 + 2._r8*t1co2(i)*tmp2 + 2._r8*tmp3)*sqti(i)
            abso(i,4)= trab2(i)*emm(i,kn)*absbnd
            tco2(i)  = 1.0_r8/(1.0_r8+ 10.0_r8*u7(i)/sqrt(4._r8 + u7(i)*(1._r8 + rbeta7(i))))
         end do ! do i =
!
! Calculate trace gas absorptivity for nearest layer, abstrc
!
         call trcabn(ncol      ,                            &
              k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
              un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
              uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
              winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
              uptype  ,dw      ,s2c     ,u       ,pnew    , &
              abstrc  ,uinpl   , &
              aer_trn_ngh)
!
! Total next layer absorptivity:
!
         do i=1,ncol
            absnxt(i,k2,kn) = abso(i,1) + abso(i,2) + &
                 abso(i,3) + abso(i,4) + abstrc(i)
         end do
      end do ! do kn =
   end do ! do k2 =

   return
end subroutine radabs

!====================================================================================

subroutine radems(lchnk   ,ncol    ,                            &
                  s2c     ,tcg     ,w       ,tplnke  ,plh2o   , &
                  pnm     ,plco2   ,tint    ,tint4   ,tlayr   , &
                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211 ,uco212   , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  , &
                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  , &
                  plh2ob  ,wb      , &
                  aer_trn_ttl, co2mmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute emissivity for H2O, CO2, O3, CH4, N2O, CFC11 and CFC12
! 
! Method: 
! H2O  ....  Uses nonisothermal emissivity method for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
!
!            Implementation updated by Collins,Hackney, and Edwards 2001
!               using line-by-line calculations based upon Hitran 1996 and
!               CKD 2.1 for absorptivity and emissivity
!
!            Implementation updated by Collins, Lee-Taylor, and Edwards (2003)
!               using line-by-line calculations based upon Hitran 2000 and
!               CKD 2.4 for absorptivity and emissivity
!
! CO2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019. Also includes the effects
!            of the 9.4 and 10.4 micron bands of CO2.
!
! O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R. Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of CFC11
!
! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of CFC12
!
!
! Computes individual emissivities, accounting for band overlap, and
! sums to obtain the total.
!
! Author: W. Collins (H2O emissivity) and J. Kiehl
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns

   real(r8), intent(in) :: s2c(pcols,pverp)        ! H2o continuum path length
   real(r8), intent(in) :: tcg(pcols,pverp)        ! H2o-mass-wgted temp. (Curtis-Godson approx.)
   real(r8), intent(in) :: w(pcols,pverp)          ! H2o path length
   real(r8), intent(in) :: tplnke(pcols)           ! Layer planck temperature
   real(r8), intent(in) :: plh2o(pcols,pverp)      ! H2o prs wghted path length
   real(r8), intent(in) :: pnm(pcols,pverp)        ! Model interface pressure
   real(r8), intent(in) :: plco2(pcols,pverp)      ! Prs wghted path of co2
   real(r8), intent(in) :: tint(pcols,pverp)       ! Model interface temperatures
   real(r8), intent(in) :: tint4(pcols,pverp)      ! Tint to the 4th power
   real(r8), intent(in) :: tlayr(pcols,pverp)      ! K-1 model layer temperature
   real(r8), intent(in) :: tlayr4(pcols,pverp)     ! Tlayr to the 4th power
   real(r8), intent(in) :: plol(pcols,pverp)       ! Pressure wghtd ozone path
   real(r8), intent(in) :: plos(pcols,pverp)       ! Ozone path
   real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! Pressure weighted h2o path with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 
   real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! H2o path length with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) 
!                               ! [fraction] Total strat. aerosol
!                               ! transmission between interfaces k1 and k2  

!
! Trace gas variables
!
   real(r8), intent(in) :: co2mmr(pcols)           ! co2 column mean mass mixing ratio
   real(r8), intent(in) :: ucfc11(pcols,pverp)     ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)     ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)      ! N2O path length
   real(r8), intent(in) :: un2o1(pcols,pverp)      ! N2O path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)       ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
   real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type continuum path length
!
! Output arguments
!
   real(r8), intent(out) :: emstot(pcols,pverp)     ! Total emissivity
   real(r8), intent(out) :: co2em(pcols,pverp)      ! Layer co2 normalzd plnck funct drvtv
   real(r8), intent(out) :: co2eml(pcols,pver)      ! Intrfc co2 normalzd plnck func drvtv
   real(r8), intent(out) :: co2t(pcols,pverp)       ! Tmp and prs weighted path length
   real(r8), intent(out) :: h2otr(pcols,pverp)      ! H2o transmission over o3 band
   real(r8), intent(out) :: abplnk1(14,pcols,pverp) ! non-nearest layer Plack factor
   real(r8), intent(out) :: abplnk2(14,pcols,pverp) ! nearest layer factor

!
!---------------------------Local variables-----------------------------
!
   integer i                    ! Longitude index
   integer k                    ! Level index]
   integer k1                   ! Level index
!
! Local variables for H2O:
!
   real(r8) h2oems(pcols,pverp)     ! H2o emissivity
   real(r8) tpathe                  ! Used to compute h2o emissivity
   real(r8) dtx(pcols)              ! Planck temperature minus 250 K
   real(r8) dty(pcols)              ! Path temperature minus 250 K
!
! The 500-800 cm^-1 emission in emis(i,4) has been combined
!              into the 0-800 cm^-1 emission in emis(i,1)
!
   real(r8) emis(pcols,2)           ! H2O emissivity 
!
!
!
   real(r8) term7(pcols,2)          ! Kl_inf(i) in eq(r8) of table A3a of R&D
   real(r8) term8(pcols,2)          ! Delta kl_inf(i) in eq(r8)
   real(r8) tr1(pcols)              ! Equation(6) in table A2 for 650-800
   real(r8) tr2(pcols)              ! Equation(6) in table A2 for 500-650
   real(r8) tr3(pcols)              ! Equation(4) in table A2 for 650-800
   real(r8) tr4(pcols)              ! Equation(4),table A2 of R&D for 500-650
   real(r8) tr7(pcols)              ! Equation (6) times eq(4) in table A2
!                                      of R&D for 650-800 cm-1 region
   real(r8) tr8(pcols)              ! Equation (6) times eq(4) in table A2
!                                      of R&D for 500-650 cm-1 region
   real(r8) k21(pcols)              ! Exponential coefficient used to calc
!                                     rot band transmissivity in the 650-800
!                                     cm-1 region (tr1)
   real(r8) k22(pcols)              ! Exponential coefficient used to calc
!                                     rot band transmissivity in the 500-650
!                                     cm-1 region (tr2)
   real(r8) u(pcols)                ! Pressure weighted H2O path length
   real(r8) ub(nbands)              ! Pressure weighted H2O path length with
                                    !  Hulst-Curtis-Godson correction for
                                    !  each band
   real(r8) pnew                    ! Effective pressure for h2o linewidth
   real(r8) pnewb(nbands)           ! Effective pressure for h2o linewidth w/
                                    !  Hulst-Curtis-Godson correction for
                                    !  each band
   real(r8) uc1(pcols)              ! H2o continuum pathlength 500-800 cm-1
   real(r8) fwk                     ! Equation(33) in R&D far wing correction
   real(r8) troco2(pcols,pverp)     ! H2o overlap factor for co2 absorption
   real(r8) emplnk(14,pcols)        ! emissivity Planck factor
   real(r8) emstrc(pcols,pverp)     ! total trace gas emissivity
!
! Local variables for CO2:
!
   real(r8) co2vmr(pcols)            ! CO2 column mean vmr
   real(r8) rmw                      ! ratio of molecular weights (air/co2)
   real(r8) co2ems(pcols,pverp)      ! Co2 emissivity
   real(r8) co2plk(pcols)            ! Used to compute co2 emissivity
   real(r8) sum(pcols)               ! Used to calculate path temperature
   real(r8) t1i                      ! Co2 hot band temperature factor
   real(r8) sqti                     ! Sqrt of temperature
   real(r8) pi                       ! Pressure used in co2 mean line width
   real(r8) et                       ! Co2 hot band factor
   real(r8) et2                      ! Co2 hot band factor
   real(r8) et4                      ! Co2 hot band factor
   real(r8) omet                     ! Co2 stimulated emission term
   real(r8) ex                       ! Part of co2 planck function
   real(r8) f1co2                    ! Co2 weak band factor
   real(r8) f2co2                    ! Co2 weak band factor
   real(r8) f3co2                    ! Co2 weak band factor
   real(r8) t1co2                    ! Overlap factor weak bands strong band
   real(r8) sqwp                     ! Sqrt of co2 pathlength
   real(r8) f1sqwp                   ! Main co2 band factor
   real(r8) oneme                    ! Co2 stimulated emission term
   real(r8) alphat                   ! Part of the co2 stimulated emiss term
   real(r8) wco2                     ! Consts used to define co2 pathlength
   real(r8) posqt                    ! Effective pressure for co2 line width
   real(r8) rbeta7                   ! Inverse of co2 hot band line width par
   real(r8) rbeta8                   ! Inverse of co2 hot band line width par
   real(r8) rbeta9                   ! Inverse of co2 hot band line width par
   real(r8) rbeta13                  ! Inverse of co2 hot band line width par
   real(r8) tpath                    ! Path temp used in co2 band model
   real(r8) tmp1                     ! Co2 band factor
   real(r8) tmp2                     ! Co2 band factor
   real(r8) tmp3                     ! Co2 band factor
   real(r8) tlayr5                   ! Temperature factor in co2 Planck func
   real(r8) rsqti                    ! Reciprocal of sqrt of temperature
   real(r8) exm1sq                   ! Part of co2 Planck function
   real(r8) u7                       ! Absorber amt for various co2 band systems
   real(r8) u8                       ! Absorber amt for various co2 band systems
   real(r8) u9                       ! Absorber amt for various co2 band systems
   real(r8) u13                      ! Absorber amt for various co2 band systems
   real(r8) r250                     ! Inverse 250K
   real(r8) r300                     ! Inverse 300K
   real(r8) rsslp                    ! Inverse standard sea-level pressure
!
! Local variables for O3:
!
   real(r8) o3ems(pcols,pverp)       ! Ozone emissivity
   real(r8) dbvtt(pcols)             ! Tmp drvtv of planck fctn for tplnke
   real(r8) dbvt,fo3,t,ux,vx
   real(r8) te                       ! Temperature factor
   real(r8) u1                       ! Path length factor
   real(r8) u2                       ! Path length factor
   real(r8) phat                     ! Effecitive path length pressure
   real(r8) tlocal                   ! Local planck function temperature
   real(r8) tcrfac                   ! Scaled temperature factor
   real(r8) beta                     ! Absorption funct factor voigt effect
   real(r8) realnu                   ! Absorption function factor
   real(r8) o3bndi                   ! Band absorption factor
!
! Transmission terms for various spectral intervals:
!
   real(r8) absbnd                   ! Proportional to co2 band absorptance
   real(r8) tco2(pcols)              ! co2 overlap factor
   real(r8) th2o(pcols)              ! h2o overlap factor
   real(r8) to3(pcols)               ! o3 overlap factor
!
! Variables for new H2O parameterization
!
! Notation:
! U   = integral (P/P_0 dW)  eq. 15 in Ramanathan/Downey 1986
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
   real(r8) fe               ! asymptotic value of emis. as U->infinity
   real(r8) e_star           ! normalized non-window emissivity
   real(r8) l_star           ! interpolated line transmission
   real(r8) c_star           ! interpolated continuum transmission

   real(r8) te1              ! emission temperature
   real(r8) te2              ! te^2
   real(r8) te3              ! te^3
   real(r8) te4              ! te^4
   real(r8) te5              ! te^5

   real(r8) log_u            ! log base 10 of U 
   real(r8) log_uc           ! log base 10 of H2O continuum path
   real(r8) log_p            ! log base 10 of P
   real(r8) t_p              ! T_p
   real(r8) t_e              ! T_e (offset by T_p)

   integer iu                ! index for log10(U)
   integer iu1               ! iu + 1
   integer iuc               ! index for log10(H2O continuum path)
   integer iuc1              ! iuc + 1
   integer ip                ! index for log10(P)
   integer ip1               ! ip + 1
   integer itp               ! index for T_p
   integer itp1              ! itp + 1
   integer ite               ! index for T_e
   integer ite1              ! ite + 1
   integer irh               ! index for RH
   integer irh1              ! irh + 1

   real(r8) dvar             ! normalized variation in T_p/T_e/P/U
   real(r8) uvar             ! U * diffusivity factor
   real(r8) uscl             ! factor for lineary scaling as U->0

   real(r8) wu               ! weight for U
   real(r8) wu1              ! 1 - wu
   real(r8) wuc              ! weight for H2O continuum path
   real(r8) wuc1             ! 1 - wuc
   real(r8) wp               ! weight for P
   real(r8) wp1              ! 1 - wp
   real(r8) wtp              ! weight for T_p
   real(r8) wtp1             ! 1 - wtp
   real(r8) wte              ! weight for T_e
   real(r8) wte1             ! 1 - wte
   real(r8) wrh              ! weight for RH
   real(r8) wrh1             ! 1 - wrh

   real(r8) w_0_0_           ! weight for Tp/Te combination
   real(r8) w_0_1_           ! weight for Tp/Te combination
   real(r8) w_1_0_           ! weight for Tp/Te combination
   real(r8) w_1_1_           ! weight for Tp/Te combination

   real(r8) w_0_00           ! weight for Tp/Te/RH combination
   real(r8) w_0_01           ! weight for Tp/Te/RH combination
   real(r8) w_0_10           ! weight for Tp/Te/RH combination
   real(r8) w_0_11           ! weight for Tp/Te/RH combination
   real(r8) w_1_00           ! weight for Tp/Te/RH combination
   real(r8) w_1_01           ! weight for Tp/Te/RH combination
   real(r8) w_1_10           ! weight for Tp/Te/RH combination
   real(r8) w_1_11           ! weight for Tp/Te/RH combination

   real(r8) w00_00           ! weight for P/Tp/Te/RH combination
   real(r8) w00_01           ! weight for P/Tp/Te/RH combination
   real(r8) w00_10           ! weight for P/Tp/Te/RH combination
   real(r8) w00_11           ! weight for P/Tp/Te/RH combination
   real(r8) w01_00           ! weight for P/Tp/Te/RH combination
   real(r8) w01_01           ! weight for P/Tp/Te/RH combination
   real(r8) w01_10           ! weight for P/Tp/Te/RH combination
   real(r8) w01_11           ! weight for P/Tp/Te/RH combination
   real(r8) w10_00           ! weight for P/Tp/Te/RH combination
   real(r8) w10_01           ! weight for P/Tp/Te/RH combination
   real(r8) w10_10           ! weight for P/Tp/Te/RH combination
   real(r8) w10_11           ! weight for P/Tp/Te/RH combination
   real(r8) w11_00           ! weight for P/Tp/Te/RH combination
   real(r8) w11_01           ! weight for P/Tp/Te/RH combination
   real(r8) w11_10           ! weight for P/Tp/Te/RH combination
   real(r8) w11_11           ! weight for P/Tp/Te/RH combination

   integer ib                ! spectral interval:
                             !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
                             !   2 = 800-1200 cm^-1

   real(r8) pch2o            ! H2O continuum path
   real(r8) fch2o            ! temp. factor for continuum
   real(r8) uch2o            ! U corresponding to H2O cont. path (window)

   real(r8) fdif             ! secant(zenith angle) for diffusivity approx.

   real(r8) sslp_mks         ! Sea-level pressure in MKS units
   real(r8) esx              ! saturation vapor pressure returned by qsat
   real(r8) qsx              ! saturation mixing ratio returned by qsat
   real(r8) pnew_mks         ! pnew in MKS units
   real(r8) q_path           ! effective specific humidity along path
   real(r8) rh_path          ! effective relative humidity along path

!
!---------------------------Statement functions-------------------------
!
! Derivative of planck function at 9.6 micro-meter wavelength, and
! an absorption function factor:
!
!
   dbvt(t)=(-2.8911366682e-4_r8+(2.3771251896e-6_r8+1.1305188929e-10_r8*t)*t)/ &
           (1.0_r8+(-6.1364820707e-3_r8+1.5550319767e-5_r8*t)*t)
!
   fo3(ux,vx)=ux/sqrt(4._r8+ux*(1._r8+vx))
!
!
!
!-----------------------------------------------------------------------
!
! Initialize
!
   r250  = 1._r8/250._r8
   r300  = 1._r8/300._r8
   rsslp = 1._r8/sslp
   rmw   = amd/amco2
   do i=1,ncol
      co2vmr(i) = co2mmr(i) * rmw
   end do
!
! Constants for computing U corresponding to H2O cont. path
!
   fdif       = 1.66_r8
   sslp_mks   = sslp / 10.0_r8
!
! Planck function for co2
!
   do i=1,ncol
      ex             = exp(960._r8/tplnke(i))
      co2plk(i)      = 5.e8_r8/((tplnke(i)**4)*(ex - 1._r8))
      co2t(i,ntoplw) = tplnke(i)
      sum(i)         = co2t(i,ntoplw)*pnm(i,ntoplw)
   end do
   k = ntoplw
   do k1=pverp,ntoplw+1,-1
      k = k + 1
      do i=1,ncol
         sum(i)         = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
         ex             = exp(960._r8/tlayr(i,k1))
         tlayr5         = tlayr(i,k1)*tlayr4(i,k1)
         co2eml(i,k1-1) = 1.2e11_r8*ex/(tlayr5*(ex - 1._r8)**2)
         co2t(i,k)      = sum(i)/pnm(i,k)
      end do
   end do
!
! Initialize planck function derivative for O3
!
   do i=1,ncol
      dbvtt(i) = dbvt(tplnke(i))
   end do
!
! Calculate trace gas Planck functions
!
   call trcplk(ncol    ,                                     &
               tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
               abplnk2 )
         
   if ( ntoplw > 1 )then
      emstot(:ncol,:ntoplw-1) = 0._r8
   end if
          
!
! Interface loop
!
   do k1=ntoplw,pverp
!
! H2O emissivity
!
! emis(i,1)     0 -  800 cm-1   h2o rotation band
! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! emis(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O emissivity
!
!      emis(i,3)   = 0.0
!
! For the p type continuum
!
      do i=1,ncol
         u(i)        = plh2o(i,k1)
         pnew        = u(i)/w(i,k1)
         pnew_mks    = pnew * sslp_mks
!
! Apply scaling factor for 500-800 continuum
!
         uc1(i)      = (s2c(i,k1) + 1.7e-3_r8*plh2o(i,k1))*(1._r8 + 2._r8*s2c(i,k1))/ &
                       (1._r8 + 15._r8*s2c(i,k1))
         pch2o       = s2c(i,k1)
!
! Changed effective path temperature to std. Curtis-Godson form
!
         tpathe   = tcg(i,k1)/w(i,k1)
         t_p = min(max(tpathe, min_tp_h2o), max_tp_h2o)

         call qsat_water(t_p, pnew_mks, esx, qsx)

!
! Compute effective RH along path
!
         q_path = w(i,k1) / pnm(i,k1) / rga
!
! Calculate effective u, pnew for each band using
!        Hulst-Curtis-Godson approximation:
! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
!           2nd edition, Oxford University Press, 1989.
! Effective H2O path (w)
!      eq. 6.24, p. 228
! Effective H2O path pressure (pnew = u/w):
!      eq. 6.29, p. 228
!
         ub(1) = plh2ob(1,i,k1) / psi(t_p,1)
         ub(2) = plh2ob(2,i,k1) / psi(t_p,2)

         pnewb(1) = ub(1) / wb(1,i,k1) * phi(t_p,1)
         pnewb(2) = ub(2) / wb(2,i,k1) * phi(t_p,2)
!
!
!
         dtx(i) = tplnke(i) - 250._r8
         dty(i) = tpathe - 250._r8
!
! Define variables for C/H/E (now C/LT/E) fit
!
! emis(i,1)     0 -  800 cm-1   h2o rotation band
! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! emis(i,2)   800 - 1200 cm-1   h2o window
!
! Separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for H2O emissivity
!
! emis(i,3)   = 0.0
!
! Notation:
! U   = integral (P/P_0 dW)  
! P   = atmospheric pressure
! P_0 = reference atmospheric pressure
! W   = precipitable water path
! T_e = emission temperature
! T_p = path temperature
! RH  = path relative humidity
!
! Terms for asymptotic value of emissivity
!
         te1  = tplnke(i)
         te2  = te1 * te1
         te3  = te2 * te1
         te4  = te3 * te1
         te5  = te4 * te1
!
! Band-independent indices for lines and continuum tables
!
         dvar = (t_p - min_tp_h2o) / dtp_h2o
         itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
         itp1 = itp + 1
         wtp = dvar - floor(dvar)
         wtp1 = 1.0_r8 - wtp

         t_e = min(max(tplnke(i) - t_p, min_te_h2o), max_te_h2o)
         dvar = (t_e - min_te_h2o) / dte_h2o
         ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
         ite1 = ite + 1
         wte = dvar - floor(dvar)
         wte1 = 1.0_r8 - wte

         rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
         dvar = (rh_path - min_rh_h2o) / drh_h2o
         irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
         irh1 = irh + 1
         wrh = dvar - floor(dvar)
         wrh1 = 1.0_r8 - wrh

         w_0_0_ = wtp  * wte
         w_0_1_ = wtp  * wte1
         w_1_0_ = wtp1 * wte 
         w_1_1_ = wtp1 * wte1

         w_0_00 = w_0_0_ * wrh
         w_0_01 = w_0_0_ * wrh1
         w_0_10 = w_0_1_ * wrh
         w_0_11 = w_0_1_ * wrh1
         w_1_00 = w_1_0_ * wrh
         w_1_01 = w_1_0_ * wrh1
         w_1_10 = w_1_1_ * wrh
         w_1_11 = w_1_1_ * wrh1
!
! H2O Continuum path for 0-800 and 1200-2200 cm^-1
!
!    Assume foreign continuum dominates total H2O continuum in these bands
!    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
!    Then the effective H2O path is just 
!         U_c = integral[ f(P) dW ]
!    where 
!           W = water-vapor mass and 
!        f(P) = dependence of foreign continuum on pressure 
!             = P / sslp
!    Then 
!         U_c = U (the same effective H2O path as for lines)
!
!
! Continuum terms for 800-1200 cm^-1
!
!    Assume self continuum dominates total H2O continuum for this band
!    per Clough et al, JGR, v. 97, no. D14 (Oct 20, 1992), p. 15776
!    Then the effective H2O self-continuum path is 
!         U_c = integral[ h(e,T) dW ]                        (*eq. 1*)
!    where 
!           W = water-vapor mass and 
!           e = partial pressure of H2O along path
!           T = temperature along path
!      h(e,T) = dependence of foreign continuum on e,T
!             = e / sslp * f(T)
!
!    Replacing
!           e =~ q * P / epsilo
!           q = mixing ratio of H2O
!     epsilo = 0.622
!
!    and using the definition
!           U = integral [ (P / sslp) dW ]
!             = (P / sslp) W                                 (homogeneous path)
!
!    the effective path length for the self continuum is
!         U_c = (q / epsilo) f(T) U                         (*eq. 2*)
!
!    Once values of T, U, and q have been calculated for the inhomogeneous
!        path, this sets U_c for the corresponding
!        homogeneous atmosphere.  However, this need not equal the
!        value of U_c' defined by eq. 1 for the actual inhomogeneous atmosphere
!        under consideration.
!
!    Solution: hold T and q constant, solve for U' that gives U_c' by
!        inverting eq. (2):
!
!        U' = (U_c * epsilo) / (q * f(T))
!
         fch2o = fh2oself(t_p)
         uch2o = (pch2o * epsilo) / (q_path * fch2o)

!
! Band-dependent indices for non-window
!
         ib = 1

         uvar = ub(ib) * fdif
         log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
         dvar = (log_u - min_lu_h2o) / dlu_h2o
         iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iu1 = iu + 1
         wu = dvar - floor(dvar)
         wu1 = 1.0_r8 - wu
         
         log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
         dvar = (log_p - min_lp_h2o) / dlp_h2o
         ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
         ip1 = ip + 1
         wp = dvar - floor(dvar)
         wp1 = 1.0_r8 - wp

         w00_00 = wp  * w_0_00 
         w00_01 = wp  * w_0_01 
         w00_10 = wp  * w_0_10 
         w00_11 = wp  * w_0_11 
         w01_00 = wp  * w_1_00 
         w01_01 = wp  * w_1_01 
         w01_10 = wp  * w_1_10 
         w01_11 = wp  * w_1_11 
         w10_00 = wp1 * w_0_00 
         w10_01 = wp1 * w_0_01 
         w10_10 = wp1 * w_0_10 
         w10_11 = wp1 * w_0_11 
         w11_00 = wp1 * w_1_00 
         w11_01 = wp1 * w_1_01 
         w11_10 = wp1 * w_1_10 
         w11_11 = wp1 * w_1_11 

!
! Asymptotic value of emissivity as U->infinity
!
         fe = fet(1,ib) + &
              fet(2,ib) * te1 + &
              fet(3,ib) * te2 + &
              fet(4,ib) * te3 + &
              fet(5,ib) * te4 + &
              fet(6,ib) * te5

         e_star = &
              eh2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
              eh2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
              eh2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
              eh2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
              eh2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
              eh2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
              eh2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
              eh2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
              eh2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
              eh2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
              eh2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
              eh2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
              eh2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
              eh2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
              eh2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
              eh2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
              eh2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
              eh2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
              eh2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
              eh2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
              eh2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
              eh2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
              eh2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
              eh2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
              eh2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
              eh2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
              eh2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
              eh2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
              eh2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
              eh2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
              eh2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
              eh2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
         emis(i,ib) = min(max(fe * (1.0_r8 - (1.0_r8 - e_star) * &
                              aer_trn_ttl(i,k1,1,ib)), &
                          0.0_r8), 1.0_r8)
!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
         if (uvar < min_u_h2o) then
            uscl = uvar / min_u_h2o
            emis(i,ib) = emis(i,ib) * uscl
         endif

                      

!
! Band-dependent indices for window
!
         ib = 2

         uvar = ub(ib) * fdif
         log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
         dvar = (log_u - min_lu_h2o) / dlu_h2o
         iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iu1 = iu + 1
         wu = dvar - floor(dvar)
         wu1 = 1.0_r8 - wu
         
         log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
         dvar = (log_p - min_lp_h2o) / dlp_h2o
         ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
         ip1 = ip + 1
         wp = dvar - floor(dvar)
         wp1 = 1.0_r8 - wp

         w00_00 = wp  * w_0_00 
         w00_01 = wp  * w_0_01 
         w00_10 = wp  * w_0_10 
         w00_11 = wp  * w_0_11 
         w01_00 = wp  * w_1_00 
         w01_01 = wp  * w_1_01 
         w01_10 = wp  * w_1_10 
         w01_11 = wp  * w_1_11 
         w10_00 = wp1 * w_0_00 
         w10_01 = wp1 * w_0_01 
         w10_10 = wp1 * w_0_10 
         w10_11 = wp1 * w_0_11 
         w11_00 = wp1 * w_1_00 
         w11_01 = wp1 * w_1_01 
         w11_10 = wp1 * w_1_10 
         w11_11 = wp1 * w_1_11 

         log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
         dvar = (log_uc - min_lu_h2o) / dlu_h2o
         iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iuc1 = iuc + 1
         wuc = dvar - floor(dvar)
         wuc1 = 1.0_r8 - wuc
!
! Asymptotic value of emissivity as U->infinity
!
         fe = fet(1,ib) + &
              fet(2,ib) * te1 + &
              fet(3,ib) * te2 + &
              fet(4,ib) * te3 + &
              fet(5,ib) * te4 + &
              fet(6,ib) * te5

         l_star = &
              ln_eh2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
              ln_eh2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
              ln_eh2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
              ln_eh2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
              ln_eh2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
              ln_eh2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
              ln_eh2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
              ln_eh2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 

         c_star = &
              cn_eh2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
              cn_eh2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
              cn_eh2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
              cn_eh2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
         emis(i,ib) = min(max(fe * (1.0_r8 - l_star * c_star * &
                              aer_trn_ttl(i,k1,1,ib)), &
                          0.0_r8), 1.0_r8) 
!
! Invoke linear limit for scaling wrt u below min_u_h2o
!
         if (uvar < min_u_h2o) then
            uscl = uvar / min_u_h2o
            emis(i,ib) = emis(i,ib) * uscl
         endif

                      
!
! Compute total emissivity for H2O
!
         h2oems(i,k1) = emis(i,1)+emis(i,2)

      end do
!
!
!

      do i=1,ncol
         term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1._r8+c16*dty(i))
         term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1._r8+c17*dty(i))
         term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1._r8+c26*dty(i))
         term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1._r8+c27*dty(i))
      end do
      do i=1,ncol
!
! 500 -  800 cm-1   rotation band overlap with co2
!
         k21(i) = term7(i,1) + term8(i,1)/ &
                 (1._r8 + (c30 + c31*(dty(i)-10._r8)*(dty(i)-10._r8))*sqrt(u(i)))
         k22(i) = term7(i,2) + term8(i,2)/ &
                 (1._r8 + (c28 + c29*(dty(i)-10._r8))*sqrt(u(i)))
         fwk    = fwcoef + fwc1/(1._r8+fwc2*u(i))
         tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
         tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
         tr1(i)=tr1(i)*aer_trn_ttl(i,k1,1,idx_LW_0650_0800) 
!                                            ! H2O line+aer trn 650--800 cm-1
         tr2(i)=tr2(i)*aer_trn_ttl(i,k1,1,idx_LW_0500_0650) 
!                                            ! H2O line+aer trn 500--650 cm-1
         tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
         tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
         tr7(i) = tr1(i)*tr3(i)
         tr8(i) = tr2(i)*tr4(i)
         troco2(i,k1) = 0.65_r8*tr7(i) + 0.35_r8*tr8(i)
         th2o(i) = tr8(i)
      end do
!
! CO2 emissivity for 15 micron band system
!
      do i=1,ncol
         t1i    = exp(-480._r8/co2t(i,k1))
         sqti   = sqrt(co2t(i,k1))
         rsqti  = 1._r8/sqti
         et     = t1i
         et2    = et*et
         et4    = et2*et2
         omet   = 1._r8 - 1.5_r8*et2
         f1co2  = 899.70_r8*omet*(1._r8 + 1.94774_r8*et + 4.73486_r8*et2)*rsqti
         sqwp   = sqrt(plco2(i,k1))
         f1sqwp = f1co2*sqwp
         t1co2  = 1._r8/(1._r8 + 245.18_r8*omet*sqwp*rsqti)
         oneme  = 1._r8 - et2
         alphat = oneme**3*rsqti
         wco2   = 2.5221_r8*co2vmr(i)*pnm(i,k1)*rga
         u7     = 4.9411e4_r8*alphat*et2*wco2
         u8     = 3.9744e4_r8*alphat*et4*wco2
         u9     = 1.0447e5_r8*alphat*et4*et2*wco2
         u13    = 2.8388e3_r8*alphat*et4*wco2
!
         tpath  = co2t(i,k1)
         tlocal = tplnke(i)
         tcrfac = sqrt((tlocal*r250)*(tpath*r300))
         pi     = pnm(i,k1)*rsslp + 2._r8*dpfco2*tcrfac
         posqt  = pi/(2._r8*sqti)
         rbeta7 =  1._r8/( 5.3288_r8*posqt)
         rbeta8 = 1._r8/ (10.6576_r8*posqt)
         rbeta9 = rbeta7
         rbeta13= rbeta9
         f2co2  = (u7/sqrt(4._r8 + u7*(1._r8 + rbeta7))) + &
                  (u8/sqrt(4._r8 + u8*(1._r8 + rbeta8))) + &
                  (u9/sqrt(4._r8 + u9*(1._r8 + rbeta9)))
         f3co2  = u13/sqrt(4._r8 + u13*(1._r8 + rbeta13))
         tmp1   = log(1._r8 + f1sqwp)
         tmp2   = log(1._r8 +  f2co2)
         tmp3   = log(1._r8 +  f3co2)
         absbnd = (tmp1 + 2._r8*t1co2*tmp2 + 2._r8*tmp3)*sqti
         tco2(i)=1.0_r8/(1.0_r8+10.0_r8*(u7/sqrt(4._r8 + u7*(1._r8 + rbeta7))))
         co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
         ex     = exp(960._r8/tint(i,k1))
         exm1sq = (ex - 1._r8)**2
         co2em(i,k1) = 1.2e11_r8*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
      end do
!
! O3 emissivity
!
      do i=1,ncol
         h2otr(i,k1) = exp(-12._r8*s2c(i,k1))
          h2otr(i,k1)=h2otr(i,k1)*aer_trn_ttl(i,k1,1,idx_LW_1000_1200)
         te          = (co2t(i,k1)/293._r8)**.7_r8
         u1          = 18.29_r8*plos(i,k1)/te
         u2          = .5649_r8*plos(i,k1)/te
         phat        = plos(i,k1)/plol(i,k1)
         tlocal      = tplnke(i)
         tcrfac      = sqrt(tlocal*r250)*te
         beta        = (1._r8/.3205_r8)*((1._r8/phat) + (dpfo3*tcrfac))
         realnu      = (1._r8/beta)*te
         o3bndi      = 74._r8*te*(tplnke(i)/375._r8)*log(1._r8 + fo3(u1,realnu) + fo3(u2,realnu))
         o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
         to3(i)=1.0_r8/(1._r8 + 0.1_r8*fo3(u1,realnu) + 0.1_r8*fo3(u2,realnu))
      end do
!
!   Calculate trace gas emissivities
!
      call trcems(ncol    ,                                     &
                  k1      ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
                  bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
                  uco222  ,uco223  ,uptype  ,w       ,s2c     , &
                  u       ,emplnk  ,th2o    ,tco2    ,to3     , &
                  emstrc  , &
                  aer_trn_ttl)
!
! Total emissivity:
!
      do i=1,ncol
         emstot(i,k1) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1)  &
                        + emstrc(i,k1)
      end do
   end do ! End of interface loop

end subroutine radems

!====================================================================================

subroutine radtpl(ncol    ,                                     &
                  tnm     ,lwupcgs ,qnm     ,pnm     ,plco2   ,plh2o   , &
                  tplnka  ,s2c     ,tcg     ,w       ,tplnke  , &
                  tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    , &
                  piln    ,plh2ob  ,wb      ,co2mmr)
!--------------------------------------------------------------------
!
! Purpose:
! Compute temperatures and path lengths for longwave radiation
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: CCM1
!------------------------------Arguments-----------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: tnm(pcols,pver)      ! Model level temperatures
   real(r8), intent(in) :: lwupcgs(pcols)       ! Surface longwave up flux
   real(r8), intent(in) :: qnm(pcols,pver)      ! Model level specific humidity
   real(r8), intent(in) :: pnm(pcols,pverp)     ! Pressure at model interfaces (dynes/cm2)
   real(r8), intent(in) :: pmln(pcols,pver)     ! Ln(pmidm1)
   real(r8), intent(in) :: piln(pcols,pverp)    ! Ln(pintm1)
   real(r8), intent(in) :: co2mmr(pcols)        ! co2 column mean mass mixing ratio
!
! Output arguments
!
   real(r8), intent(out) :: plco2(pcols,pverp)   ! Pressure weighted co2 path
   real(r8), intent(out) :: plh2o(pcols,pverp)   ! Pressure weighted h2o path
   real(r8), intent(out) :: tplnka(pcols,pverp)  ! Level temperature from interface temperatures
   real(r8), intent(out) :: s2c(pcols,pverp)     ! H2o continuum path length
   real(r8), intent(out) :: tcg(pcols,pverp)     ! H2o-mass-wgted temp. (Curtis-Godson approx.)
   real(r8), intent(out) :: w(pcols,pverp)       ! H2o path length
   real(r8), intent(out) :: tplnke(pcols)        ! Equal to tplnka
   real(r8), intent(out) :: tint(pcols,pverp)    ! Layer interface temperature
   real(r8), intent(out) :: tint4(pcols,pverp)   ! Tint to the 4th power
   real(r8), intent(out) :: tlayr(pcols,pverp)   ! K-1 level temperature
   real(r8), intent(out) :: tlayr4(pcols,pverp)  ! Tlayr to the 4th power
   real(r8), intent(out) :: plh2ob(nbands,pcols,pverp)! Pressure weighted h2o path with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 
   real(r8), intent(out) :: wb(nbands,pcols,pverp)    ! H2o path length with 
                                                      !    Hulst-Curtis-Godson temp. factor 
                                                      !    for H2O bands 

!
!---------------------------Local variables--------------------------
!
   integer i                 ! Longitude index
   integer k                 ! Level index
   integer kp1               ! Level index + 1

   real(r8) repsil               ! Inver ratio mol weight h2o to dry air
   real(r8) dy                   ! Thickness of layer for tmp interp
   real(r8) dpnm                 ! Pressure thickness of layer
   real(r8) dpnmsq               ! Prs squared difference across layer
   real(r8) dw                   ! Increment in H2O path length
   real(r8) dplh2o               ! Increment in plh2o
   real(r8) cpwpl                ! Const in co2 mix ratio to path length conversn

!--------------------------------------------------------------------
!
   repsil = 1._r8/epsilo
!
! Compute co2 and h2o paths
!
   cpwpl = 0.5_r8/(gravit_cgs*p0)
   do i=1,ncol
      plh2o(i,ntoplw)  = rgsslp*qnm(i,ntoplw)*pnm(i,ntoplw)*pnm(i,ntoplw)
      plco2(i,ntoplw)  = co2mmr(i)*cpwpl*pnm(i,ntoplw)*pnm(i,ntoplw)
   end do
   do k=ntoplw,pver
      do i=1,ncol
         plh2o(i,k+1)  = plh2o(i,k) + rgsslp* &
                         (pnm(i,k+1)**2 - pnm(i,k)**2)*qnm(i,k)
         plco2(i,k+1)  = co2mmr(i)*cpwpl*pnm(i,k+1)**2
      end do
   end do
!
! Set the top and bottom intermediate level temperatures,
! top level planck temperature and top layer temp**4.
!
! Tint is lower interface temperature
! (not available for bottom layer, so use ground temperature)
!
   do i=1,ncol
      tint4(i,pverp)   = lwupcgs(i)/stebol_cgs
      tint(i,pverp)    = sqrt(sqrt(tint4(i,pverp)))
      tplnka(i,ntoplw) = tnm(i,ntoplw)
      tint(i,ntoplw)   = tplnka(i,ntoplw)
      tlayr4(i,ntoplw) = tplnka(i,ntoplw)**4
      tint4(i,ntoplw)  = tlayr4(i,ntoplw)
   end do
!
! Intermediate level temperatures are computed using temperature
! at the full level below less dy*delta t,between the full level
!
   do k=ntoplw+1,pver
      do i=1,ncol
         dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
         tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
         tint4(i,k) = tint(i,k)**4
      end do
   end do
!
! Now set the layer temp=full level temperatures and establish a
! planck temperature for absorption (tplnka) which is the average
! the intermediate level temperatures.  Note that tplnka is not
! equal to the full level temperatures.
!
   do k=ntoplw+1,pverp
      do i=1,ncol
         tlayr(i,k)  = tnm(i,k-1)
         tlayr4(i,k) = tlayr(i,k)**4
         tplnka(i,k) = .5_r8*(tint(i,k) + tint(i,k-1))
      end do
   end do
!
! Calculate tplank for emissivity calculation.
! Assume isothermal tplnke i.e. all levels=ttop.
!
   do i=1,ncol
      tplnke(i)       = tplnka(i,ntoplw)
      tlayr(i,ntoplw) = tint(i,ntoplw)
   end do
!
! Now compute h2o path fields:
!
   do i=1,ncol
!
! Changed effective path temperature to std. Curtis-Godson form
!
      tcg(i,ntoplw) = rga*qnm(i,ntoplw)*pnm(i,ntoplw)*tnm(i,ntoplw)
      w(i,ntoplw)   = sslp * (plh2o(i,ntoplw)*2._r8) / pnm(i,ntoplw)
!
! Hulst-Curtis-Godson scaling for H2O path
!
      wb(1,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),1)
      wb(2,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),2)
!
! Hulst-Curtis-Godson scaling for effective pressure along H2O path
!
      plh2ob(1,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),1)
      plh2ob(2,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),2)

      s2c(i,ntoplw) = plh2o(i,ntoplw)*fh2oself(tnm(i,ntoplw))*qnm(i,ntoplw)*repsil
   end do

   do k=ntoplw,pver
      do i=1,ncol
         dpnm       = pnm(i,k+1) - pnm(i,k)
         dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
         dw         = rga*qnm(i,k)*dpnm
         kp1        = k+1
         w(i,kp1)   = w(i,k) + dw
!
! Hulst-Curtis-Godson scaling for H2O path
!
         wb(1,i,kp1) = wb(1,i,k) + dw * phi(tnm(i,k),1)
         wb(2,i,kp1) = wb(2,i,k) + dw * phi(tnm(i,k),2)
!
! Hulst-Curtis-Godson scaling for effective pressure along H2O path
!
         dplh2o = plh2o(i,kp1) - plh2o(i,k)

         plh2ob(1,i,kp1) = plh2ob(1,i,k) + dplh2o * psi(tnm(i,k),1)
         plh2ob(2,i,kp1) = plh2ob(2,i,k) + dplh2o * psi(tnm(i,k),2)
!
! Changed effective path temperature to std. Curtis-Godson form
!
         tcg(i,kp1) = tcg(i,k) + dw*tnm(i,k)
         s2c(i,kp1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)* &
                      fh2oself(tnm(i,k))*qnm(i,k)*repsil
      end do
   end do

end subroutine radtpl

!====================================================================================

subroutine radae_init(gravx, epsilox, stebol, pstdx, mwdryx, mwco2x, mwo3x)
!
! Initialize radae module data
!
   use pio,          only: file_desc_t, var_desc_t, pio_inq_dimid, pio_inquire_dimension, &
                           pio_inquire_variable, pio_inq_varid, pio_get_var, pio_nowrite, &
                           pio_max_var_dims, pio_max_name, pio_closefile
   use cam_pio_utils,only: cam_pio_openfile
   use ioFileMod,    only: getfil
   use filenames,    only: absems_data

!
! Input variables
!
   real(r8), intent(in) :: gravx   ! Acceleration due to gravity (m/s**2)
   real(r8), intent(in) :: epsilox ! Ratio of mol. wght of H2O to dry air
   real(r8), intent(in) :: stebol  ! Stefan-Boltzmann's constant (MKS)
   real(r8), intent(in) :: pstdx   ! Standard pressure (pascals)
   real(r8), intent(in) :: mwdryx  ! Molecular weight of dry air 
   real(r8), intent(in) :: mwco2x  ! Molecular weight of carbon dioxide
   real(r8), intent(in) :: mwo3x   ! Molecular weight of ozone
!
!      Variables for loading absorptivity/emissivity
!
   type(file_desc_T) :: ncid_ae                ! NetCDF file id for abs/ems file

   integer pdimid                 ! pressure dimension id
   integer psize                  ! pressure dimension size

   integer tpdimid                ! path temperature dimension id
   integer tpsize                 ! path temperature size

   integer tedimid                ! emission temperature dimension id
   integer tesize                 ! emission temperature size

   integer udimid                 ! u (H2O path) dimension id
   integer usize                  ! u (H2O path) dimension size

   integer rhdimid                ! relative humidity dimension id
   integer rhsize                 ! relative humidity dimension size

   type(var_desc_t) ::    ah2onwid            ! var. id for non-wndw abs.
   type(var_desc_t) ::    eh2onwid            ! var. id for non-wndw ems.
   type(var_desc_t) ::    ah2owid             ! var. id for wndw abs. (adjacent layers)
   type(var_desc_t) :: cn_ah2owid             ! var. id for continuum trans. for wndw abs.
   type(var_desc_t) :: cn_eh2owid             ! var. id for continuum trans. for wndw ems.
   type(var_desc_t) :: ln_ah2owid             ! var. id for line trans. for wndw abs.
   type(var_desc_t) :: ln_eh2owid             ! var. id for line trans. for wndw ems.
   
   character*(PIO_MAX_NAME) tmpname! dummy variable for var/dim names
   character(len=256) locfn       ! local filename
   integer tmptype                ! dummy variable for variable type
   integer ndims                  ! number of dimensions
   integer dims(PIO_MAX_VAR_DIMS)  ! vector of dimension ids
   integer natt                   ! number of attributes

   integer ierr                  ! ierr flag returned from pio (pio handles errors internally so it is not checked)
!
! Constants to set
!
   gravit     = gravx
   gravit_cgs = 100._r8*gravx
   rga        = 1._r8/gravit_cgs
   epsilo     = epsilox
   omeps      = 1._r8 - epsilo
   sslp       = 1.013250e6_r8
   stebol_cgs = 1.e3_r8*stebol
   rgsslp     = 0.5_r8/(gravit_cgs*sslp)
   dpfo3      = 2.5e-3_r8
   dpfco2     = 5.0e-3_r8

   p0         = pstdx*10.0_r8
   amd        = mwdryx
   amco2      = mwco2x
   mwo3       = mwo3x
!
! Coefficients for h2o emissivity and absorptivity for overlap of H2O 
!    and trace gases.
!
   c16  = coefj(3,1)/coefj(2,1)
   c17  = coefk(3,1)/coefk(2,1)
   c26  = coefj(3,2)/coefj(2,2)
   c27  = coefk(3,2)/coefk(2,2)
   c28  = .5_r8
   c29  = .002053_r8
   c30  = .1_r8
   c31  = 3.0e-5_r8
!
! Initialize further longwave constants referring to far wing
! correction for overlap of H2O and trace gases; R&D refers to:
!
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
   fwcoef = .1_r8           ! See eq(33) R&D
   fwc1   = .30_r8          ! See eq(33) R&D
   fwc2   = 4.5_r8          ! See eq(33) and eq(34) in R&D
   fc1    = 2.6_r8          ! See eq(34) R&D

   call getfil(absems_data, locfn)
   call cam_pio_openfile(ncid_ae, locfn, PIO_NOWRITE)
   
   ierr =  pio_inq_dimid(ncid_ae, 'p', pdimid)
   ierr =  pio_inquire_dimension(ncid_ae, pdimid, len=psize)
   
   ierr =  pio_inq_dimid(ncid_ae, 'tp', tpdimid)
   ierr =  pio_inquire_dimension(ncid_ae, tpdimid, len=tpsize)
   
   ierr =  pio_inq_dimid(ncid_ae, 'te', tedimid)
   ierr =  pio_inquire_dimension(ncid_ae, tedimid, len=tesize)
   
   ierr =  pio_inq_dimid(ncid_ae, 'u', udimid)
   ierr =  pio_inquire_dimension(ncid_ae, udimid, len=usize)
   
   ierr =  pio_inq_dimid(ncid_ae, 'rh', rhdimid)
   ierr =  pio_inquire_dimension(ncid_ae, rhdimid, len=rhsize)
      
   if (psize    /= n_p  .or. &
        tpsize   /= n_tp .or. &
        usize    /= n_u  .or. &
        tesize   /= n_te .or. &
        rhsize   /= n_rh) then
      call endrun ('RADAEINI: dimensions for abs/ems do not match internal def.')
   endif
   
   ierr =  pio_inq_varid(ncid_ae, 'ah2onw',   ah2onwid)
   ierr =  pio_inq_varid(ncid_ae, 'eh2onw',   eh2onwid)
   ierr =  pio_inq_varid(ncid_ae, 'ah2ow',    ah2owid)
   ierr =  pio_inq_varid(ncid_ae, 'cn_ah2ow', cn_ah2owid)
   ierr =  pio_inq_varid(ncid_ae, 'cn_eh2ow', cn_eh2owid)
   ierr =  pio_inq_varid(ncid_ae, 'ln_ah2ow', ln_ah2owid)
   ierr =  pio_inq_varid(ncid_ae, 'ln_eh2ow', ln_eh2owid)
      
   ierr =  pio_inquire_variable(ncid_ae, ah2onwid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: non-wndw abs. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, eh2onwid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: non-wndw ems. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, ah2owid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: window abs. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, cn_ah2owid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: cont. trans for abs. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, cn_eh2owid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: cont. trans. for ems. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, ln_ah2owid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: line trans for abs. in file /= internal def.')
   endif
   ierr =  pio_inquire_variable(ncid_ae, ln_eh2owid, tmpname, tmptype, ndims, dims, natt)
   if (ndims /= 5 .or. &
        dims(1) /= pdimid    .or. &
        dims(2) /= tpdimid   .or. &
        dims(3) /= udimid    .or. &
        dims(4) /= tedimid   .or. &
        dims(5) /= rhdimid) then
      call endrun ('RADAEINI: line trans. for ems. in file /= internal def.')
   endif
   
   ierr =  pio_get_var (ncid_ae, ah2onwid,   ah2onw)
   ierr =  pio_get_var (ncid_ae, eh2onwid,   eh2onw)
   ierr =  pio_get_var (ncid_ae, ah2owid,    ah2ow)
   ierr =  pio_get_var (ncid_ae, cn_ah2owid, cn_ah2ow)
   ierr =  pio_get_var (ncid_ae, cn_eh2owid, cn_eh2ow)
   ierr =  pio_get_var (ncid_ae, ln_ah2owid, ln_ah2ow)
   ierr =  pio_get_var (ncid_ae, ln_eh2owid, ln_eh2ow)
      
   call pio_closefile(ncid_ae)

end subroutine radae_init

!====================================================================================

subroutine initialize_radbuffer
!
! Initialize radiation buffer data
!

  use ref_pres, only : pref_mid
  use phys_control, only : phys_getopts

  character(len=16) :: radiation_scheme
  integer :: k

! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
! longwave cooling to about 80 km (1 Pa)
   if (pref_mid(1) .lt. 0.1_r8) then
      do k = 1, plev
         if (pref_mid(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
   if (masterproc) then
      write(iulog,*) 'INITIALIZE_RADBUFFER: ntoplw =',ntoplw, ' pressure:',pref_mid(ntoplw)
   endif
 
   call phys_getopts(radiation_scheme_out=radiation_scheme)
   
  if(radiation_scheme.eq.'camrt') then
     allocate (abstot_3d(pcols,ntoplw:pverp,ntoplw:pverp,begchunk:endchunk))
     allocate (absnxt_3d(pcols,pver,4,begchunk:endchunk))
     allocate (emstot_3d(pcols,pverp,begchunk:endchunk))
     abstot_3d(:,:,:,:) = posinf
     absnxt_3d(:,:,:,:) = posinf
     emstot_3d(:,:,:) = posinf
  end if
  return
end subroutine initialize_radbuffer

!====================================================================================

subroutine radoz2(ncol, o3, pint, plol, plos)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes the path length integrals to the model interfaces given the
! ozone volume mixing ratio
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact J. Kiehl
! 
!------------------------------Input arguments--------------------------
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: o3(pcols,pver)      ! ozone mass mixing ratio
   real(r8), intent(in) :: pint(pcols,pverp)   ! Model interface pressures
!
!----------------------------Output arguments---------------------------
!
   real(r8), intent(out) :: plol(pcols,pverp)   ! Ozone prs weighted path length (cm)
   real(r8), intent(out) :: plos(pcols,pverp)   ! Ozone path length (cm)
!
!---------------------------Local workspace-----------------------------
!
   integer i                ! longitude index
   integer k                ! level index

   real(r8) :: v0    ! Volume of a gas at stp (m**3/kmol)
   real(r8) :: p0    ! Standard pressure (pascals)
   real(r8) :: cplos ! constant for ozone path length integral
   real(r8) :: cplol ! constant for ozone path length integral

!
!-----------------------------------------------------------------------
!*******************************************************************
! These hardwired constants need to be replaced with common values.
! They are here for testing infrastructure changes that should not
! change answers.
! Constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
   v0    = 22.4136_r8         ! Volume of a gas at stp (m**3/kmol)
   p0    = 0.1_r8*sslp        ! Standard pressure (pascals)
   cplos = v0/(mwo3*gravit)       *100.0_r8
   cplol = v0/(mwo3*gravit*p0)*0.5_r8*100.0_r8
!*******************************************************************
!
! Evaluate the ozone path length integrals to interfaces;
! factors of .1 and .01 to convert pressures from cgs to mks:
!
   do i=1,ncol
      plos(i,ntoplw) = 0.1_r8 *cplos*o3(i,ntoplw)*pint(i,ntoplw)
      plol(i,ntoplw) = 0.01_r8*cplol*o3(i,ntoplw)*pint(i,ntoplw)*pint(i,ntoplw)
   end do
   do k=ntoplw+1,pverp
      do i=1,ncol
         plos(i,k) = plos(i,k-1) + 0.1_r8*cplos*o3(i,k-1)*(pint(i,k) - pint(i,k-1))
         plol(i,k) = plol(i,k-1) + 0.01_r8*cplol*o3(i,k-1)* &
                    (pint(i,k)*pint(i,k) - pint(i,k-1)*pint(i,k-1))
      end do
   end do

end subroutine radoz2

!====================================================================================

subroutine trcpth(ncol                                        , &
                  tnm     ,pnm     ,cfc11   ,cfc12   ,n2o     , &
                  ch4     ,qnm     ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                  bch4    ,uptype  ,co2mmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate path lengths and pressure factors for CH4, N2O, CFC11
! and CFC12.
! 
! Method: 
! See CCM3 description for details
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: tnm(pcols,pver)      ! Model level temperatures
   real(r8), intent(in) :: pnm(pcols,pverp)     ! Pres. at model interfaces (dynes/cm2)
   real(r8), intent(in) :: qnm(pcols,pver)      ! h2o specific humidity
   real(r8), intent(in) :: cfc11(pcols,pver)    ! CFC11 mass mixing ratio
!
   real(r8), intent(in) :: cfc12(pcols,pver)    ! CFC12 mass mixing ratio
   real(r8), intent(in) :: n2o(pcols,pver)      ! N2O mass mixing ratio
   real(r8), intent(in) :: ch4(pcols,pver)      ! CH4 mass mixing ratio
   real(r8), intent(in) :: co2mmr(pcols)        ! co2 column mean mass mixing ratio

!
! Output arguments
!
   real(r8), intent(out) :: ucfc11(pcols,pverp)  ! CFC11 path length
   real(r8), intent(out) :: ucfc12(pcols,pverp)  ! CFC12 path length
   real(r8), intent(out) :: un2o0(pcols,pverp)   ! N2O path length
   real(r8), intent(out) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
   real(r8), intent(out) :: uch4(pcols,pverp)    ! CH4 path length
!
   real(r8), intent(out) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(out) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(out) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(out) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(out) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
!
   real(r8), intent(out) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(out) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(out) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(out) :: bch4(pcols,pverp)    ! pressure factor for ch4
   real(r8), intent(out) :: uptype(pcols,pverp)  ! p-type continuum path length

!
!---------------------------Local variables-----------------------------
!
   integer   i               ! Longitude index
   integer   k               ! Level index
!
   real(r8) co2fac(pcols,1)      ! co2 factor
   real(r8) alpha1(pcols)        ! stimulated emission term
   real(r8) alpha2(pcols)        ! stimulated emission term
   real(r8) rt(pcols)            ! reciprocal of local temperature
   real(r8) rsqrt(pcols)         ! reciprocal of sqrt of temp
!
   real(r8) pbar(pcols)          ! mean pressure
   real(r8) dpnm(pcols)          ! difference in pressure
   real(r8) diff                 ! diffusivity factor
!
!--------------------------Data Statements------------------------------
!
   data diff /1.66_r8/
!
!-----------------------------------------------------------------------
!
!  Calculate path lengths for the trace gases at model top
!

   do i = 1,ncol
      ucfc11(i,ntoplw) = 1.8_r8 * cfc11(i,ntoplw) * pnm(i,ntoplw) * rga
      ucfc12(i,ntoplw) = 1.8_r8 * cfc12(i,ntoplw) * pnm(i,ntoplw) * rga
      un2o0(i,ntoplw) = diff * 1.02346e5_r8 * n2o(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
      un2o1(i,ntoplw) = diff * 2.01909_r8 * un2o0(i,ntoplw) * exp(-847.36_r8/tnm(i,ntoplw))
      uch4(i,ntoplw)  = diff * 8.60957e4_r8 * ch4(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
      co2fac(i,1)     = diff * co2mmr(i) * pnm(i,ntoplw) * rga
      alpha1(i) = (1.0_r8 - exp(-1540.0_r8/tnm(i,ntoplw)))**3.0_r8/sqrt(tnm(i,ntoplw))
      alpha2(i) = (1.0_r8 - exp(-1360.0_r8/tnm(i,ntoplw)))**3.0_r8/sqrt(tnm(i,ntoplw))
      uco211(i,ntoplw) = 3.42217e3_r8 * co2fac(i,1) * alpha1(i) * exp(-1849.7_r8/tnm(i,ntoplw))
      uco212(i,ntoplw) = 6.02454e3_r8 * co2fac(i,1) * alpha1(i) * exp(-2782.1_r8/tnm(i,ntoplw))
      uco213(i,ntoplw) = 5.53143e3_r8 * co2fac(i,1) * alpha1(i) * exp(-3723.2_r8/tnm(i,ntoplw))
      uco221(i,ntoplw) = 3.88984e3_r8 * co2fac(i,1) * alpha2(i) * exp(-1997.6_r8/tnm(i,ntoplw))
      uco222(i,ntoplw) = 3.67108e3_r8 * co2fac(i,1) * alpha2(i) * exp(-3843.8_r8/tnm(i,ntoplw))
      uco223(i,ntoplw) = 6.50642e3_r8 * co2fac(i,1) * alpha2(i) * exp(-2989.7_r8/tnm(i,ntoplw))
      bn2o0(i,ntoplw) = diff * 19.399_r8 * pnm(i,ntoplw)**2.0_r8 * n2o(i,ntoplw) * &
                   1.02346e5_r8 * rga / (sslp*tnm(i,ntoplw))
      bn2o1(i,ntoplw) = bn2o0(i,ntoplw) * exp(-847.36_r8/tnm(i,ntoplw)) * 2.06646e5_r8
      bch4(i,ntoplw) = diff * 2.94449_r8 * ch4(i,ntoplw) * pnm(i,ntoplw)**2.0_r8 * rga * &
                  8.60957e4_r8 / (sslp*tnm(i,ntoplw))
      uptype(i,ntoplw) = diff * qnm(i,ntoplw) * pnm(i,ntoplw)**2.0_r8 *  &
                    exp(1800.0_r8*(1.0_r8/tnm(i,ntoplw) - 1.0_r8/296.0_r8)) * rga / sslp
   end do
!
! Calculate trace gas path lengths through model atmosphere
!
   do k = ntoplw,pver
      do i = 1,ncol
         rt(i) = 1._r8/tnm(i,k)
         rsqrt(i) = sqrt(rt(i))
         pbar(i) = 0.5_r8 * (pnm(i,k+1) + pnm(i,k)) / sslp
         dpnm(i) = (pnm(i,k+1) - pnm(i,k)) * rga
         alpha1(i) = diff * rsqrt(i) * (1.0_r8 - exp(-1540.0_r8/tnm(i,k)))**3.0_r8
         alpha2(i) = diff * rsqrt(i) * (1.0_r8 - exp(-1360.0_r8/tnm(i,k)))**3.0_r8
         ucfc11(i,k+1) = ucfc11(i,k) +  1.8_r8 * cfc11(i,k) * dpnm(i)
         ucfc12(i,k+1) = ucfc12(i,k) +  1.8_r8 * cfc12(i,k) * dpnm(i)
         un2o0(i,k+1) = un2o0(i,k) + diff * 1.02346e5_r8 * n2o(i,k) * rsqrt(i) * dpnm(i)
         un2o1(i,k+1) = un2o1(i,k) + diff * 2.06646e5_r8 * n2o(i,k) * &
                        rsqrt(i) * exp(-847.36_r8/tnm(i,k)) * dpnm(i)
         uch4(i,k+1) = uch4(i,k) + diff * 8.60957e4_r8 * ch4(i,k) * rsqrt(i) * dpnm(i)
         uco211(i,k+1) = uco211(i,k) + 1.15_r8*3.42217e3_r8 * alpha1(i) * &
                         co2mmr(i) * exp(-1849.7_r8/tnm(i,k)) * dpnm(i)
         uco212(i,k+1) = uco212(i,k) + 1.15_r8*6.02454e3_r8 * alpha1(i) * &
                         co2mmr(i) * exp(-2782.1_r8/tnm(i,k)) * dpnm(i)
         uco213(i,k+1) = uco213(i,k) + 1.15_r8*5.53143e3_r8 * alpha1(i) * &
                         co2mmr(i) * exp(-3723.2_r8/tnm(i,k)) * dpnm(i)
         uco221(i,k+1) = uco221(i,k) + 1.15_r8*3.88984e3_r8 * alpha2(i) * &
                         co2mmr(i) * exp(-1997.6_r8/tnm(i,k)) * dpnm(i)
         uco222(i,k+1) = uco222(i,k) + 1.15_r8*3.67108e3_r8 * alpha2(i) * &
                         co2mmr(i) * exp(-3843.8_r8/tnm(i,k)) * dpnm(i)
         uco223(i,k+1) = uco223(i,k) + 1.15_r8*6.50642e3_r8 * alpha2(i) * &
                         co2mmr(i) * exp(-2989.7_r8/tnm(i,k)) * dpnm(i)
         bn2o0(i,k+1) = bn2o0(i,k) + diff * 19.399_r8 * pbar(i) * rt(i) &
                        * 1.02346e5_r8 * n2o(i,k) * dpnm(i)
         bn2o1(i,k+1) = bn2o1(i,k) + diff * 19.399_r8 * pbar(i) * rt(i) &
                        * 2.06646e5_r8 * exp(-847.36_r8/tnm(i,k)) * n2o(i,k)*dpnm(i)
         bch4(i,k+1) = bch4(i,k) + diff * 2.94449_r8 * rt(i) * pbar(i) &
                       * 8.60957e4_r8 * ch4(i,k) * dpnm(i)
         uptype(i,k+1) = uptype(i,k) + diff *qnm(i,k) * &
                         exp(1800.0_r8*(1.0_r8/tnm(i,k) - 1.0_r8/296.0_r8)) * pbar(i) * dpnm(i)
      end do
   end do
!
   return
end subroutine trcpth



!====================================================================================
! Private Interfaces
!====================================================================================

function fh2oself( temp )
!
! Short function for H2O self-continuum temperature factor in
!   calculation of effective H2O self-continuum path length
!
! H2O Continuum: CKD 2.4
! Code for continuum: GENLN3
! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
!                     Transmittance and Radiance Model, Version 3.0 Description
!                     and Users Guide, NCAR/TN-367+STR, 147 pp.
!
! In GENLN, the temperature scaling of the self-continuum is handled
!    by exponential interpolation/extrapolation from observations at
!    260K and 296K by:
!
!         TFAC =  (T(IPATH) - 296.0)/(260.0 - 296.0) 
!         CSFFT = CSFF296*(CSFF260/CSFF296)**TFAC
!
! For 800-1200 cm^-1, (CSFF260/CSFF296) ranges from ~2.1 to ~1.9
!     with increasing wavenumber.  The ratio <CSFF260>/<CSFF296>,
!     where <> indicates average over wavenumber, is ~2.07
! 
! fh2oself is (<CSFF260>/<CSFF296>)**TFAC
!
   real(r8),intent(in) :: temp     ! path temperature
   real(r8) fh2oself               ! mean ratio of self-continuum at temp and 296K

   fh2oself = 2.0727484_r8**((296.0_r8 - temp) / 36.0_r8)
end function fh2oself

!====================================================================================

function phi(tpx,iband)
!
! History: First version for Hitran 1996 (C/H/E)
!          Current version for Hitran 2000 (C/LT/E)
! Short function for Hulst-Curtis-Godson temperature factors for
!   computing effective H2O path
! Line data for H2O: Hitran 2000, plus H2O patches v11.0 for 1341 missing
!                    lines between 500 and 2820 cm^-1.
!                    See cfa-www.harvard.edu/HITRAN
! Isotopes of H2O: all
! Line widths: air-broadened only (self set to 0)
! Code for line strengths and widths: GENLN3
! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
!                     Transmittance and Radiance Model, Version 3.0 Description
!                     and Users Guide, NCAR/TN-367+STR, 147 pp.
!
! Note: functions have been normalized by dividing by their values at
!       a path temperature of 160K
!
! spectral intervals:
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1
!   2 = 800-1200 cm^-1      
!
! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
!           2nd edition, Oxford University Press, 1989.
! Phi: function for H2O path
!      eq. 6.25, p. 228
!
   real(r8),intent(in):: tpx      ! path temperature
   integer, intent(in):: iband    ! band to process
   real(r8) phi                   ! phi for given band
   real(r8),parameter ::  phi_r0(nbands) = (/ 9.60917711E-01_r8, -2.21031342E+01_r8/)
   real(r8),parameter ::  phi_r1(nbands) = (/ 4.86076751E-04_r8,  4.24062610E-01_r8/)
   real(r8),parameter ::  phi_r2(nbands) = (/-1.84806265E-06_r8, -2.95543415E-03_r8/)
   real(r8),parameter ::  phi_r3(nbands) = (/ 2.11239959E-09_r8,  7.52470896E-06_r8/)

   phi = (((phi_r3(iband) * tpx) + phi_r2(iband)) * tpx + phi_r1(iband)) &
          * tpx + phi_r0(iband)
end function phi

!====================================================================================

function psi(tpx,iband)
!
! History: First version for Hitran 1996 (C/H/E)
!          Current version for Hitran 2000 (C/LT/E)
! Short function for Hulst-Curtis-Godson temperature factors for
!   computing effective H2O path
! Line data for H2O: Hitran 2000, plus H2O patches v11.0 for 1341 missing
!                    lines between 500 and 2820 cm^-1. 
!                    See cfa-www.harvard.edu/HITRAN
! Isotopes of H2O: all
! Line widths: air-broadened only (self set to 0)
! Code for line strengths and widths: GENLN3
! Reference: Edwards, D.P., 1992: GENLN2, A General Line-by-Line Atmospheric
!                     Transmittance and Radiance Model, Version 3.0 Description
!                     and Users Guide, NCAR/TN-367+STR, 147 pp.
!
! Note: functions have been normalized by dividing by their values at
!       a path temperature of 160K
!
! spectral intervals:
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1
!   2 = 800-1200 cm^-1      
!
! Formulae: Goody and Yung, Atmospheric Radiation: Theoretical Basis, 
!           2nd edition, Oxford University Press, 1989.
! Psi: function for pressure along path
!      eq. 6.30, p. 228
!
   real(r8),intent(in):: tpx      ! path temperature
   integer, intent(in):: iband    ! band to process
   real(r8) psi                   ! psi for given band
   real(r8),parameter ::  psi_r0(nbands) = (/ 5.65308452E-01_r8, -7.30087891E+01_r8/)
   real(r8),parameter ::  psi_r1(nbands) = (/ 4.07519005E-03_r8,  1.22199547E+00_r8/)
   real(r8),parameter ::  psi_r2(nbands) = (/-1.04347237E-05_r8, -7.12256227E-03_r8/)
   real(r8),parameter ::  psi_r3(nbands) = (/ 1.23765354E-08_r8,  1.47852825E-05_r8/)

   psi = (((psi_r3(iband) * tpx) + psi_r2(iband)) * tpx + psi_r1(iband)) * tpx + psi_r0(iband)

end function psi

!====================================================================================

subroutine trcab(ncol    ,                                     &
                 k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
                 un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                 uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                 bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
                 s2c     ,uptype  ,dplh2o  ,abplnk1 ,tco2    , &
                 th2o    ,to3     ,abstrc  , &
                 aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
! CFC12.
! 
! Method: 
! See CCM3 description for equations.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                     ! number of atmospheric columns
   integer, intent(in) :: k1,k2                    ! level indices
!
   real(r8), intent(in) :: to3co2(pcols)           ! pressure weighted temperature
   real(r8), intent(in) :: pnm(pcols,pverp)        ! interface pressures
   real(r8), intent(in) :: ucfc11(pcols,pverp)     ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)     ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)      ! N2O path length
!
   real(r8), intent(in) :: un2o1(pcols,pverp)      ! N2O path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)       ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)     ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)     ! CO2 9.4 micron band path length
!
   real(r8), intent(in) :: uco221(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)     ! CO2 10.4 micron band path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
!
   real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
   real(r8), intent(in) :: dw(pcols)               ! h2o path length
   real(r8), intent(in) :: pnew(pcols)             ! pressure
   real(r8), intent(in) :: s2c(pcols,pverp)        ! continuum path length
   real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type h2o path length
!
   real(r8), intent(in) :: dplh2o(pcols)           ! p squared h2o path length
   real(r8), intent(in) :: abplnk1(14,pcols,pverp) ! Planck factor
   real(r8), intent(in) :: tco2(pcols)             ! co2 transmission factor
   real(r8), intent(in) :: th2o(pcols)             ! h2o transmission factor
   real(r8), intent(in) :: to3(pcols)              ! o3 transmission factor

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) ! aer trn.

!
!  Output Arguments
!
   real(r8), intent(out) :: abstrc(pcols)           ! total trace gas absorptivity
!
!--------------------------Local Variables------------------------------
!
   integer  i,l                     ! loop counters

   real(r8) sqti(pcols)             ! square root of mean temp
   real(r8) du1                     ! cfc11 path length
   real(r8) du2                     ! cfc12 path length
   real(r8) acfc1                   ! cfc11 absorptivity 798 cm-1
   real(r8) acfc2                   ! cfc11 absorptivity 846 cm-1
!
   real(r8) acfc3                   ! cfc11 absorptivity 933 cm-1
   real(r8) acfc4                   ! cfc11 absorptivity 1085 cm-1
   real(r8) acfc5                   ! cfc12 absorptivity 889 cm-1
   real(r8) acfc6                   ! cfc12 absorptivity 923 cm-1
   real(r8) acfc7                   ! cfc12 absorptivity 1102 cm-1
!
   real(r8) acfc8                   ! cfc12 absorptivity 1161 cm-1
   real(r8) du01                    ! n2o path length
   real(r8) dbeta01                 ! n2o pressure factor
   real(r8) dbeta11                 !         "
   real(r8) an2o1                   ! absorptivity of 1285 cm-1 n2o band
!
   real(r8) du02                    ! n2o path length
   real(r8) dbeta02                 ! n2o pressure factor
   real(r8) an2o2                   ! absorptivity of 589 cm-1 n2o band
   real(r8) du03                    ! n2o path length
   real(r8) dbeta03                 ! n2o pressure factor
!
   real(r8) an2o3                   ! absorptivity of 1168 cm-1 n2o band
   real(r8) duch4                   ! ch4 path length
   real(r8) dbetac                  ! ch4 pressure factor
   real(r8) ach4                    ! absorptivity of 1306 cm-1 ch4 band
   real(r8) du11                    ! co2 path length
!
   real(r8) du12                    !       "
   real(r8) du13                    !       "
   real(r8) dbetc1                  ! co2 pressure factor
   real(r8) dbetc2                  ! co2 pressure factor
   real(r8) aco21                   ! absorptivity of 1064 cm-1 band
!
   real(r8) du21                    ! co2 path length
   real(r8) du22                    !       "
   real(r8) du23                    !       "
   real(r8) aco22                   ! absorptivity of 961 cm-1 band
   real(r8) tt(pcols)               ! temp. factor for h2o overlap factor
!
   real(r8) psi1                    !                 "
   real(r8) phi1                    !                 "
   real(r8) p1                      ! h2o overlap factor
   real(r8) w1                      !        "
   real(r8) ds2c(pcols)             ! continuum path length
!
   real(r8) duptyp(pcols)           ! p-type path length
   real(r8) tw(pcols,6)             ! h2o transmission factor
!   real(r8) g1(6)                   !         "
!   real(r8) g2(6)                   !         "
!   real(r8) g3(6)                   !         "
!
!   real(r8) g4(6)                   !         "
!   real(r8) ab(6)                   ! h2o temp. factor
!   real(r8) bb(6)                   !         "
!   real(r8) abp(6)                  !         "
!   real(r8) bbp(6)                  !         "
!
   real(r8) tcfc3                   ! transmission for cfc11 band
   real(r8) tcfc4                   ! transmission for cfc11 band
   real(r8) tcfc6                   ! transmission for cfc12 band
   real(r8) tcfc7                   ! transmission for cfc12 band
   real(r8) tcfc8                   ! transmission for cfc12 band
!
   real(r8) tlw                     ! h2o transmission
   real(r8) tch4                    ! ch4 transmission
!
!--------------------------Data Statements------------------------------
!
!   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
!   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
!   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
!   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
!   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
!   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
!   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
!   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
!
!--------------------------Statement Functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
!
!------------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(to3co2(i))
!
! h2o transmission
!
      tt(i) = abs(to3co2(i) - 250.0_r8)
      ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
      duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i) + bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i) + bb(l)*tt(i)*tt(i))
         p1 = pnew(i)*(psi1/phi1)/sslp
         w1 = dw(i)*phi1
         tw(i,l) = exp(-g1(l)*p1*(sqrt(1.0_r8 + g2(l)*(w1/p1)) - 1.0_r8) - &
                   g3(l)*ds2c(i)-g4(l)*duptyp(i))
      end do
   end do
!
   do i=1,ncol
      tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ttl(i,k1,k2,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
                       0.3_r8*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000)) 
      tw(i,2)=tw(i,2)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
      tw(i,3)=tw(i,3)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
      tw(i,4)=tw(i,4)*aer_trn_ttl(i,k1,k2,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
      tw(i,5)=tw(i,5)*aer_trn_ttl(i,k1,k2,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
      tw(i,6)=tw(i,6)*aer_trn_ttl(i,k1,k2,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
   end do                    ! end loop over lon
   do i = 1,ncol
      du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
      du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
!
! cfc transmissions
!
      tcfc3 = exp(-175.005_r8*du1)
      tcfc4 = exp(-1202.18_r8*du1)
      tcfc6 = exp(-5786.73_r8*du2)
      tcfc7 = exp(-2873.51_r8*du2)
      tcfc8 = exp(-2085.59_r8*du2)
!
! Absorptivity for CFC11 bands
!
      acfc1 =  50.0_r8*(1.0_r8 - exp(-54.09_r8*du1))*tw(i,1)*abplnk1(7,i,k2)
      acfc2 =  60.0_r8*(1.0_r8 - exp(-5130.03_r8*du1))*tw(i,2)*abplnk1(8,i,k2)
      acfc3 =  60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6*abplnk1(9,i,k2)
      acfc4 = 100.0_r8*(1.0_r8 - tcfc4)*tw(i,5)*abplnk1(10,i,k2)
!
! Absorptivity for CFC12 bands
!
      acfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*du2))*tw(i,3)*abplnk1(11,i,k2)
      acfc6 = 50.0_r8*(1.0_r8 - tcfc6)* tw(i,4) * abplnk1(12,i,k2)
      acfc7 = 80.0_r8*(1.0_r8 - tcfc7)* tw(i,5) * tcfc4*abplnk1(13,i,k2)
      acfc8 = 70.0_r8*(1.0_r8 - tcfc8)* tw(i,6) * abplnk1(14,i,k2)
!
! Emissivity for CH4 band 1306 cm-1
!
      tlw = exp(-1.0_r8*sqrt(dplh2o(i)))
      tlw=tlw*aer_trn_ttl(i,k1,k2,idx_LW_1200_2000)
      duch4 = abs(uch4(i,k1) - uch4(i,k2))
      dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
      ach4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(duch4,dbetac))*tlw*abplnk1(3,i,k2)
      tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(duch4,dbetac))
!
! Absorptivity for N2O bands
!
      du01 = abs(un2o0(i,k1) - un2o0(i,k2))
      du11 = abs(un2o1(i,k1) - un2o1(i,k2))
      dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
      dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
!
! 1285 cm-1 band
!
      an2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(du01,dbeta01) &
              + func(du11,dbeta11))*tlw*tch4*abplnk1(4,i,k2)
      du02 = 0.100090_r8*du01
      du12 = 0.0992746_r8*du11
      dbeta02 = 0.964282_r8*dbeta01
!
! 589 cm-1 band
!
      an2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(du02,dbeta02) + &
              func(du12,dbeta02))*th2o(i)*tco2(i)*abplnk1(5,i,k2)
      du03 = 0.0333767_r8*du01
      dbeta03 = 0.982143_r8*dbeta01
!
! 1168 cm-1 band
!
      an2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(du03,dbeta03))* &
              tw(i,6)*tcfc8*abplnk1(6,i,k2)
!
! Emissivity for 1064 cm-1 band of CO2
!
      du11 = abs(uco211(i,k1) - uco211(i,k2))
      du12 = abs(uco212(i,k1) - uco212(i,k2))
      du13 = abs(uco213(i,k1) - uco213(i,k2))
      dbetc1 = 2.97558_r8*abs(pnm(i,k1) + pnm(i,k2))/(2.0_r8*sslp*sqti(i))
      dbetc2 = 2.0_r8*dbetc1
      aco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(du11,dbetc1) &
              + func(du12,dbetc2) + func(du13,dbetc2)) &
              *to3(i)*tw(i,5)*tcfc4*tcfc7*abplnk1(2,i,k2)
!
! Emissivity for 961 cm-1 band
!
      du21 = abs(uco221(i,k1) - uco221(i,k2))
      du22 = abs(uco222(i,k1) - uco222(i,k2))
      du23 = abs(uco223(i,k1) - uco223(i,k2))
      aco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(du21,dbetc1) &
              + func(du22,dbetc1) + func(du23,dbetc2)) &
              *tw(i,4)*tcfc3*tcfc6*abplnk1(1,i,k2)
!
! total trace gas absorptivity
!
      abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
                  aco21 + aco22
   end do

end subroutine trcab

!====================================================================================

subroutine trcabn(ncol    ,                                     &
                  k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
                  winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
                  uptype  ,dw      ,s2c     ,up2     ,pnew    , &
                  abstrc  ,uinpl   , &
                  aer_trn_ngh)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate nearest layer absorptivity due to CH4, N2O, CFC11 and CFC12
! 
! Method: 
! Equations in CCM3 description
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: k2                    ! level index
   integer, intent(in) :: kn                    ! level index
!
   real(r8), intent(in) :: tbar(pcols,4)        ! pressure weighted temperature
   real(r8), intent(in) :: ucfc11(pcols,pverp)  ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)  ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)   ! N2O path length
   real(r8), intent(in) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
!
   real(r8), intent(in) :: uch4(pcols,pverp)    ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
!
   real(r8), intent(in) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(in) :: bplnk(14,pcols,4)    ! weighted Planck fnc. for absorptivity
   real(r8), intent(in) :: winpl(pcols,4)       ! fractional path length
   real(r8), intent(in) :: pinpl(pcols,4)       ! pressure factor for subdivided layer
!
   real(r8), intent(in) :: tco2(pcols)          ! co2 transmission
   real(r8), intent(in) :: th2o(pcols)          ! h2o transmission
   real(r8), intent(in) :: to3(pcols)           ! o3 transmission
   real(r8), intent(in) :: dw(pcols)            ! h2o path length
   real(r8), intent(in) :: pnew(pcols)          ! pressure factor
!
   real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum factor
   real(r8), intent(in) :: uptype(pcols,pverp)  ! p-type path length
   real(r8), intent(in) :: up2(pcols)           ! p squared path length
   real(r8), intent(in) :: uinpl(pcols,4)       ! Nearest layer subdivision factor
   real(r8), intent(in) :: aer_trn_ngh(pcols,nlwbands) 
                             ! [fraction] Total transmission between 
                             !            nearest neighbor sub-levels
!
!  Output Arguments
!
   real(r8), intent(out) :: abstrc(pcols)        ! total trace gas absorptivity

!
!--------------------------Local Variables------------------------------
!
   integer i,l                   ! loop counters
!
   real(r8) sqti(pcols)          ! square root of mean temp
   real(r8) rsqti(pcols)         ! reciprocal of sqti
   real(r8) du1                  ! cfc11 path length
   real(r8) du2                  ! cfc12 path length
   real(r8) acfc1                ! absorptivity of cfc11 798 cm-1 band
!
   real(r8) acfc2                ! absorptivity of cfc11 846 cm-1 band
   real(r8) acfc3                ! absorptivity of cfc11 933 cm-1 band
   real(r8) acfc4                ! absorptivity of cfc11 1085 cm-1 band
   real(r8) acfc5                ! absorptivity of cfc11 889 cm-1 band
   real(r8) acfc6                ! absorptivity of cfc11 923 cm-1 band
!
   real(r8) acfc7                ! absorptivity of cfc11 1102 cm-1 band
   real(r8) acfc8                ! absorptivity of cfc11 1161 cm-1 band
   real(r8) du01                 ! n2o path length
   real(r8) dbeta01              ! n2o pressure factors
   real(r8) dbeta11              !        "
!
   real(r8)  an2o1               ! absorptivity of the 1285 cm-1 n2o band
   real(r8) du02                 ! n2o path length
   real(r8) dbeta02              ! n2o pressure factor
   real(r8) an2o2                ! absorptivity of the 589 cm-1 n2o band
   real(r8) du03                 ! n2o path length
!
   real(r8) dbeta03              ! n2o pressure factor
   real(r8) an2o3                ! absorptivity of the 1168 cm-1 n2o band
   real(r8) duch4                ! ch4 path length
   real(r8) dbetac               ! ch4 pressure factor
   real(r8) ach4                 ! absorptivity of the 1306 cm-1 ch4 band
!
   real(r8) du11                 ! co2 path length
   real(r8) du12                 !       "
   real(r8) du13                 !       "
   real(r8) dbetc1               ! co2 pressure factor
   real(r8) dbetc2               ! co2 pressure factor
!
   real(r8) aco21                ! absorptivity of the 1064 cm-1 co2 band
   real(r8) du21                 ! co2 path length
   real(r8) du22                 !       "
   real(r8) du23                 !       "
   real(r8) aco22                ! absorptivity of the 961 cm-1 co2 band
!
   real(r8) tt(pcols)            ! temp. factor for h2o overlap
   real(r8) psi1                 !          "
   real(r8) phi1                 !          "
   real(r8) p1                   ! factor for h2o overlap
   real(r8) w1                   !          "
!
   real(r8) ds2c(pcols)          ! continuum path length
   real(r8) duptyp(pcols)        ! p-type path length
   real(r8) tw(pcols,6)          ! h2o transmission overlap
!   real(r8) g1(6)                ! h2o overlap factor
!   real(r8) g2(6)                !         "
!
!   real(r8) g3(6)                !         "
!   real(r8) g4(6)                !         "
!   real(r8) ab(6)                ! h2o temp. factor
!   real(r8) bb(6)                !         "
!   real(r8) abp(6)               !         "
!
!   real(r8) bbp(6)               !         "
   real(r8) tcfc3                ! transmission of cfc11 band
   real(r8) tcfc4                ! transmission of cfc11 band
   real(r8) tcfc6                ! transmission of cfc12 band
   real(r8) tcfc7                !         "
!
   real(r8) tcfc8                !         "
   real(r8) tlw                  ! h2o transmission
   real(r8) tch4                 ! ch4 transmission
!
!--------------------------Data Statements------------------------------
!
!   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
!   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
!   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
!   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
!   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
!   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
!   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
!   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
!
!--------------------------Statement Functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
!
!------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(tbar(i,kn))
      rsqti(i) = 1._r8 / sqti(i)
!
! h2o transmission
!
      tt(i) = abs(tbar(i,kn) - 250.0_r8)
      ds2c(i) = abs(s2c(i,k2+1) - s2c(i,k2))*uinpl(i,kn)
      duptyp(i) = abs(uptype(i,k2+1) - uptype(i,k2))*uinpl(i,kn)
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
         p1 = pnew(i) * (psi1/phi1) / sslp
         w1 = dw(i) * winpl(i,kn) * phi1
         tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0_r8+g2(l)*(w1/p1))-1.0_r8) &
                   - g3(l)*ds2c(i)-g4(l)*duptyp(i))
      end do
   end do
!
   do i=1,ncol
      tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ngh(i,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
                       0.3_r8*aer_trn_ngh(i,idx_LW_0800_1000))
      tw(i,2)=tw(i,2)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
      tw(i,3)=tw(i,3)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
      tw(i,4)=tw(i,4)*aer_trn_ngh(i,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
      tw(i,5)=tw(i,5)*aer_trn_ngh(i,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
      tw(i,6)=tw(i,6)*aer_trn_ngh(i,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
   end do                    ! end loop over lon

   do i = 1,ncol
!
      du1 = abs(ucfc11(i,k2+1) - ucfc11(i,k2)) * winpl(i,kn)
      du2 = abs(ucfc12(i,k2+1) - ucfc12(i,k2)) * winpl(i,kn)
!
! cfc transmissions
!
      tcfc3 = exp(-175.005_r8*du1)
      tcfc4 = exp(-1202.18_r8*du1)
      tcfc6 = exp(-5786.73_r8*du2)
      tcfc7 = exp(-2873.51_r8*du2)
      tcfc8 = exp(-2085.59_r8*du2)
!
! Absorptivity for CFC11 bands
!
      acfc1 = 50.0_r8*(1.0_r8 - exp(-54.09_r8*du1)) * tw(i,1)*bplnk(7,i,kn)
      acfc2 = 60.0_r8*(1.0_r8 - exp(-5130.03_r8*du1))*tw(i,2)*bplnk(8,i,kn)
      acfc3 = 60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6 * bplnk(9,i,kn)
      acfc4 = 100.0_r8*(1.0_r8 - tcfc4)* tw(i,5) * bplnk(10,i,kn)
!
! Absorptivity for CFC12 bands
!
      acfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*du2))*tw(i,3)*bplnk(11,i,kn)
      acfc6 = 50.0_r8*(1.0_r8 - tcfc6)*tw(i,4)*bplnk(12,i,kn)
      acfc7 = 80.0_r8*(1.0_r8 - tcfc7)* tw(i,5)*tcfc4 *bplnk(13,i,kn)
      acfc8 = 70.0_r8*(1.0_r8 - tcfc8)*tw(i,6)*bplnk(14,i,kn)
!
! Absorptivity for CH4 band 1306 cm-1
!
      tlw = exp(-1.0_r8*sqrt(up2(i)))
      tlw=tlw*aer_trn_ngh(i,idx_LW_1200_2000)
      duch4 = abs(uch4(i,k2+1) - uch4(i,k2)) * winpl(i,kn)
      dbetac = 2.94449_r8 * pinpl(i,kn) * rsqti(i) / sslp
      ach4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(duch4,dbetac)) * tlw * bplnk(3,i,kn)
      tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(duch4,dbetac))
!
! Absorptivity for N2O bands
!
      du01 = abs(un2o0(i,k2+1) - un2o0(i,k2)) * winpl(i,kn)
      du11 = abs(un2o1(i,k2+1) - un2o1(i,k2)) * winpl(i,kn)
      dbeta01 = 19.399_r8 *  pinpl(i,kn) * rsqti(i) / sslp
      dbeta11 = dbeta01
!
! 1285 cm-1 band
!
      an2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(du01,dbeta01) &
              + func(du11,dbeta11)) * tlw * tch4 * bplnk(4,i,kn)
      du02 = 0.100090_r8*du01
      du12 = 0.0992746_r8*du11
      dbeta02 = 0.964282_r8*dbeta01
!
! 589 cm-1 band
!
      an2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(du02,dbeta02) &
              +  func(du12,dbeta02)) * tco2(i) * th2o(i) * bplnk(5,i,kn)
      du03 = 0.0333767_r8*du01
      dbeta03 = 0.982143_r8*dbeta01
!
! 1168 cm-1 band
!
      an2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(du03,dbeta03)) * &
              tw(i,6) * tcfc8 * bplnk(6,i,kn)
!
! Absorptivity for 1064 cm-1 band of CO2
!
      du11 = abs(uco211(i,k2+1) - uco211(i,k2)) * winpl(i,kn)
      du12 = abs(uco212(i,k2+1) - uco212(i,k2)) * winpl(i,kn)
      du13 = abs(uco213(i,k2+1) - uco213(i,k2)) * winpl(i,kn)
      dbetc1 = 2.97558_r8 * pinpl(i,kn) * rsqti(i) / sslp
      dbetc2 = 2.0_r8 * dbetc1
      aco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(du11,dbetc1) &
              + func(du12,dbetc2) + func(du13,dbetc2)) &
              * to3(i) * tw(i,5) * tcfc4 * tcfc7 * bplnk(2,i,kn)
!
! Absorptivity for 961 cm-1 band of co2
!
      du21 = abs(uco221(i,k2+1) - uco221(i,k2)) * winpl(i,kn)
      du22 = abs(uco222(i,k2+1) - uco222(i,k2)) * winpl(i,kn)
      du23 = abs(uco223(i,k2+1) - uco223(i,k2)) * winpl(i,kn)
      aco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(du21,dbetc1) &
              + func(du22,dbetc1) + func(du23,dbetc2)) &
              * tw(i,4) * tcfc3 * tcfc6 * bplnk(1,i,kn)
!
! total trace gas absorptivity
!
      abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
                  aco21 + aco22
   end do

end subroutine trcabn

!====================================================================================

subroutine trcems(ncol    ,                                     &
                  k       ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
                  bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
                  uco222  ,uco223  ,uptype  ,w       ,s2c     , &
                  up2     ,emplnk  ,th2o    ,tco2    ,to3     , &
                  emstrc  , &
                  aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Calculate emissivity for CH4, N2O, CFC11 and CFC12 bands.
! 
! Method: 
!  See CCM3 Description for equations.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: co2t(pcols,pverp)    ! pressure weighted temperature
   real(r8), intent(in) :: pnm(pcols,pverp)     ! interface pressure
   real(r8), intent(in) :: ucfc11(pcols,pverp)  ! CFC11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)  ! CFC12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)   ! N2O path length
!
   real(r8), intent(in) :: un2o1(pcols,pverp)   ! N2O path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)    ! CH4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)  ! CO2 9.4 micron band path length
!
   real(r8), intent(in) :: uco221(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r8), intent(in) :: uptype(pcols,pverp)  ! continuum path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
!
   real(r8), intent(in) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)    ! pressure factor for ch4
   real(r8), intent(in) :: emplnk(14,pcols)     ! emissivity Planck factor
   real(r8), intent(in) :: th2o(pcols)          ! water vapor overlap factor
   real(r8), intent(in) :: tco2(pcols)          ! co2 overlap factor
!
   real(r8), intent(in) :: to3(pcols)           ! o3 overlap factor
   real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum path length
   real(r8), intent(in) :: w(pcols,pverp)       ! h2o path length
   real(r8), intent(in) :: up2(pcols)           ! pressure squared h2o path length
!
   integer, intent(in) :: k                 ! level index

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,nlwbands) ! aer trn.

!
!  Output Arguments
!
   real(r8), intent(out) :: emstrc(pcols,pverp)  ! total trace gas emissivity

!
!--------------------------Local Variables------------------------------
!
   integer i,l               ! loop counters
!
   real(r8) sqti(pcols)          ! square root of mean temp
   real(r8) ecfc1                ! emissivity of cfc11 798 cm-1 band
   real(r8) ecfc2                !     "      "    "   846 cm-1 band
   real(r8) ecfc3                !     "      "    "   933 cm-1 band
   real(r8) ecfc4                !     "      "    "   1085 cm-1 band
!
   real(r8) ecfc5                !     "      "  cfc12 889 cm-1 band
   real(r8) ecfc6                !     "      "    "   923 cm-1 band
   real(r8) ecfc7                !     "      "    "   1102 cm-1 band
   real(r8) ecfc8                !     "      "    "   1161 cm-1 band
   real(r8) u01                  ! n2o path length
!
   real(r8) u11                  ! n2o path length
   real(r8) beta01               ! n2o pressure factor
   real(r8) beta11               ! n2o pressure factor
   real(r8) en2o1                ! emissivity of the 1285 cm-1 N2O band
   real(r8) u02                  ! n2o path length
!
   real(r8) u12                  ! n2o path length
   real(r8) beta02               ! n2o pressure factor
   real(r8) en2o2                ! emissivity of the 589 cm-1 N2O band
   real(r8) u03                  ! n2o path length
   real(r8) beta03               ! n2o pressure factor
!
   real(r8) en2o3                ! emissivity of the 1168 cm-1 N2O band
   real(r8) betac                ! ch4 pressure factor
   real(r8) ech4                 ! emissivity of 1306 cm-1 CH4 band
   real(r8) betac1               ! co2 pressure factor
   real(r8) betac2               ! co2 pressure factor
!
   real(r8) eco21                ! emissivity of 1064 cm-1 CO2 band
   real(r8) eco22                ! emissivity of 961 cm-1 CO2 band
   real(r8) tt(pcols)            ! temp. factor for h2o overlap factor
   real(r8) psi1                 ! narrow band h2o temp. factor
   real(r8) phi1                 !             "
!
   real(r8) p1                   ! h2o line overlap factor
   real(r8) w1                   !          "
   real(r8) tw(pcols,6)          ! h2o transmission overlap
!   real(r8) g1(6)                ! h2o overlap factor
!   real(r8) g2(6)                !          "
!
!   real(r8) g3(6)                !          "
!   real(r8) g4(6)                !          "
!   real(r8) ab(6)                !          "
!   real(r8) bb(6)                !          "
!   real(r8) abp(6)               !          "
!
!   real(r8) bbp(6)               !          "
   real(r8) tcfc3                ! transmission for cfc11 band
   real(r8) tcfc4                !          "
   real(r8) tcfc6                ! transmission for cfc12 band
   real(r8) tcfc7                !          "
!
   real(r8) tcfc8                !          "
   real(r8) tlw                  ! h2o overlap factor
   real(r8) tch4                 ! ch4 overlap factor
!
!--------------------------Data Statements------------------------------
!
!   data g1 /0.0468556_r8,0.0397454_r8,0.0407664_r8,0.0304380_r8,0.0540398_r8,0.0321962_r8/
!   data g2 /14.4832_r8,4.30242_r8,5.23523_r8,3.25342_r8,0.698935_r8,16.5599_r8/
!   data g3 /26.1898_r8,18.4476_r8,15.3633_r8,12.1927_r8,9.14992_r8,8.07092_r8/
!   data g4 /0.0261782_r8,0.0369516_r8,0.0307266_r8,0.0243854_r8,0.0182932_r8,0.0161418_r8/
!   data ab /3.0857e-2_r8,2.3524e-2_r8,1.7310e-2_r8,2.6661e-2_r8,2.8074e-2_r8,2.2915e-2_r8/
!   data bb /-1.3512e-4_r8,-6.8320e-5_r8,-3.2609e-5_r8,-1.0228e-5_r8,-9.5743e-5_r8,-1.0304e-4_r8/
!   data abp/2.9129e-2_r8,2.4101e-2_r8,1.9821e-2_r8,2.6904e-2_r8,2.9458e-2_r8,1.9892e-2_r8/
!   data bbp/-1.3139e-4_r8,-5.5688e-5_r8,-4.6380e-5_r8,-8.0362e-5_r8,-1.0115e-4_r8,-8.8061e-5_r8/
!
!--------------------------Statement Functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0_r8 + u*(1.0_r8 + 1.0_r8 / b))
!
!-----------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(co2t(i,k))
!
! Transmission for h2o
!
      tt(i) = abs(co2t(i,k) - 250.0_r8)
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
         p1 = pnm(i,k) * (psi1/phi1) / sslp
         w1 = w(i,k) * phi1
         tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0_r8+g2(l)*(w1/p1))-1.0_r8) &
                   - g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
      end do
   end do

!     Overlap H2O tranmission with STRAER continuum in 6 trace gas 
!                 subbands

      do i=1,ncol
         tw(i,1)=tw(i,1)*(0.7_r8*aer_trn_ttl(i,k,1,idx_LW_0650_0800)+&! l=1: 0750--0820 cm-1
                          0.3_r8*aer_trn_ttl(i,k,1,idx_LW_0800_1000))
         tw(i,2)=tw(i,2)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=2: 0820--0880 cm-1
         tw(i,3)=tw(i,3)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=3: 0880--0900 cm-1
         tw(i,4)=tw(i,4)*aer_trn_ttl(i,k,1,idx_LW_0800_1000) ! l=4: 0900--1000 cm-1
         tw(i,5)=tw(i,5)*aer_trn_ttl(i,k,1,idx_LW_1000_1200) ! l=5: 1000--1120 cm-1
         tw(i,6)=tw(i,6)*aer_trn_ttl(i,k,1,idx_LW_1000_1200) ! l=6: 1120--1170 cm-1
      end do                    ! end loop over lon
!
   do i = 1,ncol
!
! transmission due to cfc bands
!
      tcfc3 = exp(-175.005_r8*ucfc11(i,k))
      tcfc4 = exp(-1202.18_r8*ucfc11(i,k))
      tcfc6 = exp(-5786.73_r8*ucfc12(i,k))
      tcfc7 = exp(-2873.51_r8*ucfc12(i,k))
      tcfc8 = exp(-2085.59_r8*ucfc12(i,k))
!
! Emissivity for CFC11 bands
!
      ecfc1 = 50.0_r8*(1.0_r8 - exp(-54.09_r8*ucfc11(i,k))) * tw(i,1) * emplnk(7,i)
      ecfc2 = 60.0_r8*(1.0_r8 - exp(-5130.03_r8*ucfc11(i,k)))* tw(i,2) * emplnk(8,i)
      ecfc3 = 60.0_r8*(1.0_r8 - tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
      ecfc4 = 100.0_r8*(1.0_r8 - tcfc4)*tw(i,5)*emplnk(10,i)
!
! Emissivity for CFC12 bands
!
      ecfc5 = 45.0_r8*(1.0_r8 - exp(-1272.35_r8*ucfc12(i,k)))*tw(i,3)*emplnk(11,i)
      ecfc6 = 50.0_r8*(1.0_r8 - tcfc6)*tw(i,4)*emplnk(12,i)
      ecfc7 = 80.0_r8*(1.0_r8 - tcfc7)*tw(i,5)* tcfc4 * emplnk(13,i)
      ecfc8 = 70.0_r8*(1.0_r8 - tcfc8)*tw(i,6) * emplnk(14,i)
!
! Emissivity for CH4 band 1306 cm-1
!
      tlw = exp(-1.0_r8*sqrt(up2(i)))

!     Overlap H2O vibration rotation band with STRAER continuum 
!             for CH4 1306 cm-1 and N2O 1285 cm-1 bands

            tlw=tlw*aer_trn_ttl(i,k,1,idx_LW_1200_2000)
      betac = bch4(i,k)/uch4(i,k)
      ech4 = 6.00444_r8*sqti(i)*log(1.0_r8 + func(uch4(i,k),betac)) *tlw * emplnk(3,i)
      tch4 = 1.0_r8/(1.0_r8 + 0.02_r8*func(uch4(i,k),betac))
!
! Emissivity for N2O bands
!
      u01 = un2o0(i,k)
      u11 = un2o1(i,k)
      beta01 = bn2o0(i,k)/un2o0(i,k)
      beta11 = bn2o1(i,k)/un2o1(i,k)
!
! 1285 cm-1 band
!
      en2o1 = 2.35558_r8*sqti(i)*log(1.0_r8 + func(u01,beta01) + &
              func(u11,beta11))*tlw*tch4*emplnk(4,i)
      u02 = 0.100090_r8*u01
      u12 = 0.0992746_r8*u11
      beta02 = 0.964282_r8*beta01
!
! 589 cm-1 band
!
      en2o2 = 2.65581_r8*sqti(i)*log(1.0_r8 + func(u02,beta02) + &
              func(u12,beta02)) * tco2(i) * th2o(i) * emplnk(5,i)
      u03 = 0.0333767_r8*u01
      beta03 = 0.982143_r8*beta01
!
! 1168 cm-1 band
!
      en2o3 = 2.54034_r8*sqti(i)*log(1.0_r8 + func(u03,beta03)) * &
              tw(i,6) * tcfc8 * emplnk(6,i)
!
! Emissivity for 1064 cm-1 band of CO2
!
      betac1 = 2.97558_r8*pnm(i,k) / (sslp*sqti(i))
      betac2 = 2.0_r8 * betac1
      eco21 = 3.7571_r8*sqti(i)*log(1.0_r8 + func(uco211(i,k),betac1) &
              + func(uco212(i,k),betac2) + func(uco213(i,k),betac2)) &
              * to3(i) * tw(i,5) * tcfc4 * tcfc7 * emplnk(2,i)
!
! Emissivity for 961 cm-1 band
!
      eco22 = 3.8443_r8*sqti(i)*log(1.0_r8 + func(uco221(i,k),betac1) &
              + func(uco222(i,k),betac1) + func(uco223(i,k),betac2)) &
              * tw(i,4) * tcfc3 * tcfc6 * emplnk(1,i)
!
! total trace gas emissivity
!
      emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 +ecfc6 + &
                    ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 + &
                    eco21 + eco22
   end do

end subroutine trcems

!====================================================================================

subroutine trcplk(ncol    ,                                     &
                  tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
                  abplnk2 )
!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Calculate Planck factors for absorptivity and emissivity of
!   CH4, N2O, CFC11 and CFC12
! 
! Method: 
!   Planck function and derivative evaluated at the band center.
! 
! Author: J. Kiehl
! 
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: tint(pcols,pverp)   ! interface temperatures
   real(r8), intent(in) :: tlayr(pcols,pverp)  ! k-1 level temperatures
   real(r8), intent(in) :: tplnke(pcols)       ! Top Layer temperature
!
! output arguments
!
   real(r8), intent(out) :: emplnk(14,pcols)         ! emissivity Planck factor
   real(r8), intent(out) :: abplnk1(14,pcols,pverp)  ! non-nearest layer Plack factor
   real(r8), intent(out) :: abplnk2(14,pcols,pverp)  ! nearest layer factor

!
!--------------------------Local Variables------------------------------
!
   integer wvl                   ! wavelength index
   integer i,k                   ! loop counters
!
   real(r8) f1(14)                   ! Planck function factor
   real(r8) f2(14)                   !        "
   real(r8) f3(14)                   !        "
!
!--------------------------Data Statements------------------------------
!
   data f1 /5.85713e8_r8,7.94950e8_r8,1.47009e9_r8,1.40031e9_r8,1.34853e8_r8, &
            1.05158e9_r8,3.35370e8_r8,3.99601e8_r8,5.35994e8_r8,8.42955e8_r8, &
            4.63682e8_r8,5.18944e8_r8,8.83202e8_r8,1.03279e9_r8/
   data f2 /2.02493e11_r8,3.04286e11_r8,6.90698e11_r8,6.47333e11_r8, &
            2.85744e10_r8,4.41862e11_r8,9.62780e10_r8,1.21618e11_r8, &
            1.79905e11_r8,3.29029e11_r8,1.48294e11_r8,1.72315e11_r8, &
            3.50140e11_r8,4.31364e11_r8/
   data f3 /1383.0_r8,1531.0_r8,1879.0_r8,1849.0_r8,848.0_r8,1681.0_r8, &
            1148.0_r8,1217.0_r8,1343.0_r8,1561.0_r8,1279.0_r8,1328.0_r8, &
            1586.0_r8,1671.0_r8/
!
!-----------------------------------------------------------------------
!
! Calculate emissivity Planck factor
!
   do wvl = 1,14
      do i = 1,ncol
         emplnk(wvl,i) = f1(wvl)/(tplnke(i)**4.0_r8*(exp(f3(wvl)/tplnke(i))-1.0_r8))
      end do
   end do
!
! Calculate absorptivity Planck factor for tint and tlayr temperatures
!
   do wvl = 1,14
      do k = ntoplw, pverp
         do i = 1, ncol
!
! non-nearlest layer function
!
            abplnk1(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tint(i,k)))  &
                               /(tint(i,k)**5.0_r8*(exp(f3(wvl)/tint(i,k))-1.0_r8)**2.0_r8)
!
! nearest layer function
!
            abplnk2(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tlayr(i,k))) &
                               /(tlayr(i,k)**5.0_r8*(exp(f3(wvl)/tlayr(i,k))-1.0_r8)**2.0_r8)
         end do
      end do
   end do

end subroutine trcplk


!====================================================================================

end module radae
