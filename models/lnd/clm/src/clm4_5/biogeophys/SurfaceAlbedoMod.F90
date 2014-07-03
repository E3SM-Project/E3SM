module SurfaceAlbedoMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceAlbedoMod
!
! !DESCRIPTION:
! Performs surface albedo calculations
!
! !PUBLIC TYPES:
  use clm_varcon  , only : istsoil
  use clm_varpar  , only : numrad, nlevcan
  use clm_varcon  , only : istcrop
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : nlevsno
  use SNICARMod   , only : sno_nbr_aer, SNICAR_RT, DO_SNO_AER, DO_SNO_OC

  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
!
! !PUBLIC DATA MEMBERS:
! The CLM default albice values are too high.
! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  real(r8), public  :: albice(numrad) = &       ! albedo land ice by waveband (1=vis, 2=nir)
                       (/ 0.80_r8, 0.55_r8 /)
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilAlbedo    ! Determine ground surface albedo
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceAlbedo
!
! !INTERFACE:
  subroutine SurfaceAlbedo(lbg, ubg, lbc, ubc, lbp, ubp, &
                           num_nourbanc, filter_nourbanc, &
                           num_nourbanp, filter_nourbanp, &
                           nextsw_cday, declinp1)
!
! !DESCRIPTION:
! Surface albedo and two-stream fluxes
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! Calculate sunlit and shaded fluxes as described by
! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
! a multi-layer canopy to calculate APAR profile 
! The calling sequence is:
! -> SurfaceAlbedo:   albedos for next time step
!    -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!    -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
!    -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
!    -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir,vis dif, nir dir, nir dif)
!

! !USES:
    use clmtype
    use shr_orb_mod
    use clm_time_manager, only : get_nstep
    use abortutils      , only : endrun
    use clm_varctl      , only : iulog, subgridflag

!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbg, ubg                   ! gridcell bounds
    integer , intent(in) :: lbc, ubc                   ! column bounds
    integer , intent(in) :: lbp, ubp                   ! pft bounds
    integer , intent(in) :: num_nourbanc               ! number of columns in non-urban filter
    integer , intent(in) :: filter_nourbanc(ubc-lbc+1) ! column filter for non-urban points
    integer , intent(in) :: num_nourbanp               ! number of pfts in non-urban filter
    integer , intent(in) :: filter_nourbanp(ubp-lbp+1) ! pft filter for non-urban points
    real(r8), intent(in) :: nextsw_cday                   ! calendar day at Greenwich (1.00, ..., days/year)
    real(r8), intent(in) :: declinp1                   ! declination angle (radians) for next time step
!
! !CALLED FROM:
! subroutine clm_driver1
! subroutine iniTimeVar
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrate to new data structures
! 8/20/03, Mariana Vertenstein: Vectorized routine
! 11/3/03, Peter Thornton: added decl(c) output for use in CN code.
! 03/28/08, Mark Flanner: added SNICAR, which required reversing the
!  order of calls to SNICAR_RT and SoilAlbedo and the location where
!  ground albedo is calculated
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    logical , pointer :: pactive(:)   ! true=>do computations on this pft (see reweightMod for details)
    integer , pointer :: pgridcell(:) ! gridcell of corresponding pft
    integer , pointer :: plandunit(:) ! index into landunit level quantities
    integer , pointer :: itypelun(:)  ! landunit type
    integer , pointer :: pcolumn(:)   ! column of corresponding pft
    integer , pointer :: cgridcell(:) ! gridcell of corresponding column
    real(r8), pointer :: lat(:)       ! gridcell latitude (radians)
    real(r8), pointer :: lon(:)       ! gridcell longitude (radians)
    real(r8), pointer :: tlai(:)      ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)      ! one-sided stem area index, no burying by snow
    real(r8), pointer :: elai(:)      ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)      ! one-sided stem area index with burying by snow
    real(r8), pointer :: h2osno(:)    ! snow water (mm H2O)
    real(r8), pointer :: rhol(:,:)    ! leaf reflectance: 1=vis, 2=nir
    real(r8), pointer :: rhos(:,:)    ! stem reflectance: 1=vis, 2=nir
    real(r8), pointer :: taul(:,:)    ! leaf transmittance: 1=vis, 2=nir
    real(r8), pointer :: taus(:,:)    ! stem transmittance: 1=vis, 2=nir
    integer , pointer :: ivt(:)       ! pft vegetation type
!
! local pointers toimplicit out arguments
!
    real(r8), pointer :: coszen(:)	    ! cosine of solar zenith angle
    real(r8), pointer :: albgrd(:,:)        ! ground albedo (direct)
    real(r8), pointer :: albgri(:,:)        ! ground albedo (diffuse)
    real(r8), pointer :: albd(:,:)          ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)          ! surface albedo (diffuse)
    real(r8), pointer :: fabd(:,:)          ! flux absorbed by canopy per unit direct flux
    real(r8), pointer :: fabd_sun(:,:)      ! flux absorbed by sunlit canopy per unit direct flux
    real(r8), pointer :: fabd_sha(:,:)      ! flux absorbed by shaded canopy per unit direct flux
    real(r8), pointer :: fabi(:,:)          ! flux absorbed by canopy per unit diffuse flux
    real(r8), pointer :: fabi_sun(:,:)      ! flux absorbed by sunlit canopy per unit diffuse flux
    real(r8), pointer :: fabi_sha(:,:)      ! flux absorbed by shaded canopy per unit diffuse flux
    real(r8), pointer :: ftdd(:,:)          ! down direct flux below canopy per unit direct flux
    real(r8), pointer :: ftid(:,:)          ! down diffuse flux below canopy per unit direct flux
    real(r8), pointer :: ftii(:,:)          ! down diffuse flux below canopy per unit diffuse flux
    real(r8), pointer :: vcmaxcintsun(:)    ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(r8), pointer :: vcmaxcintsha(:)    ! leaf to canopy scaling coefficient, shaded leaf vcmax
    integer , pointer :: ncan(:)            ! number of canopy layers
    integer , pointer :: nrad(:)            ! number of canopy layers, above snow for radiative transfer
    real(r8), pointer :: fabd_sun_z(:,:)    ! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabd_sha_z(:,:)    ! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabi_sun_z(:,:)    ! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabi_sha_z(:,:)    ! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fsun_z(:,:)        ! sunlit fraction of canopy layer
    real(r8), pointer :: tlai_z(:,:)        ! tlai increment for canopy layer
    real(r8), pointer :: tsai_z(:,:)        ! tsai increment for canopy layer
    real(r8), pointer :: decl(:)            ! solar declination angle (radians)
    real(r8), pointer :: frac_sno(:)        ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid water content (col,lyr) [kg/m2]
    real(r8), pointer :: h2osoi_ice(:,:)    ! ice lens content (col,lyr) [kg/m2]
    real(r8), pointer :: mss_cnc_bcphi(:,:) ! mass concentration of hydrophilic BC (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_bcpho(:,:) ! mass concentration of hydrophobic BC (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocphi(:,:) ! mass concentration of hydrophilic OC (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocpho(:,:) ! mass concentration of hydrophobic OC (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst1(:,:)  ! mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst2(:,:)  ! mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst3(:,:)  ! mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst4(:,:)  ! mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]
    real(r8), pointer :: albsod(:,:)        ! direct-beam soil albedo (col,bnd) [frc]
    real(r8), pointer :: albsoi(:,:)        ! diffuse soil albedo (col,bnd) [frc]
    real(r8), pointer :: flx_absdv(:,:)     ! direct flux absorption factor (col,lyr): VIS [frc]
    real(r8), pointer :: flx_absdn(:,:)     ! direct flux absorption factor (col,lyr): NIR [frc]
    real(r8), pointer :: flx_absiv(:,:)     ! diffuse flux absorption factor (col,lyr): VIS [frc]
    real(r8), pointer :: flx_absin(:,:)     ! diffuse flux absorption factor (col,lyr): NIR [frc]
    real(r8), pointer :: snw_rds(:,:)       ! snow grain radius (col,lyr) [microns]
    real(r8), pointer :: albgrd_pur(:,:)    ! pure snow ground albedo (direct)
    real(r8), pointer :: albgri_pur(:,:)    ! pure snow ground albedo (diffuse)
    real(r8), pointer :: albgrd_bc(:,:)     ! ground albedo without BC (direct)
    real(r8), pointer :: albgri_bc(:,:)     ! ground albedo without BC (diffuse)
    real(r8), pointer :: albgrd_oc(:,:)     ! ground albedo without OC (direct)
    real(r8), pointer :: albgri_oc(:,:)     ! ground albedo without OC (diffuse)
    real(r8), pointer :: albgrd_dst(:,:)    ! ground albedo without dust (direct)
    real(r8), pointer :: albgri_dst(:,:)    ! ground albedo without dust (diffuse)
    real(r8), pointer :: albsnd_hst(:,:)    ! snow albedo, direct, for history files (col,bnd) [frc]
    real(r8), pointer :: albsni_hst(:,:)    ! snow ground albedo, diffuse, for history files (col,bnd) [frc]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    real(r8), parameter :: mpe = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8) :: extkn                      ! nitrogen allocation coefficient
    integer  :: fp,fc,g,c,p,iv             ! indices
    integer  :: ib                         ! band index
    integer  :: ic                         ! 0=unit incoming direct; 1=unit incoming diffuse
    real(r8) :: wl(lbp:ubp)                ! fraction of LAI+SAI that is LAI
    real(r8) :: ws(lbp:ubp)                ! fraction of LAI+SAI that is SAI
    real(r8) :: dinc                       ! lai+sai increment for canopy layer
    real(r8) :: dincmax                    ! maximum lai+sai increment for canopy layer
    real(r8) :: dincmax_sum                ! cumulative sum of maximum lai+sai increment for canopy layer
    real(r8) :: laisum                     ! sum of canopy layer lai for error check
    real(r8) :: saisum                     ! sum of canopy layer sai for error check
    real(r8) :: blai(lbp:ubp)              ! lai buried by snow: tlai - elai
    real(r8) :: bsai(lbp:ubp)              ! sai buried by snow: tsai - esai
    real(r8) :: rho(lbp:ubp,numrad)        ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8) :: tau(lbp:ubp,numrad)        ! leaf/stem tran weighted by fraction LAI and SAI
    real(r8) :: albsnd(lbc:ubc,numrad)     ! snow albedo (direct)
    real(r8) :: albsni(lbc:ubc,numrad)     ! snow albedo (diffuse)
    real(r8) :: coszen_gcell(lbg:ubg)      ! cosine solar zenith angle for next time step (gridcell level)
    real(r8) :: coszen_col(lbc:ubc)        ! cosine solar zenith angle for next time step (pft level)
    real(r8) :: coszen_pft(lbp:ubp)        ! cosine solar zenith angle for next time step (pft level)
    integer  :: num_vegsol                 ! number of vegetated pfts where coszen>0
    integer  :: filter_vegsol(ubp-lbp+1)   ! pft filter where vegetated and coszen>0
    integer  :: num_novegsol               ! number of vegetated pfts where coszen>0
    integer  :: filter_novegsol(ubp-lbp+1) ! pft filter where vegetated and coszen>0
    integer, parameter :: nband =numrad    ! number of solar radiation waveband classes
    integer  :: flg_slr                    ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer  :: flg_snw_ice                ! flag for SNICAR (=1 when called from CLM, =2 when called from sea-ice)
    real(r8) :: albsnd_pur(lbc:ubc,numrad) ! direct pure snow albedo (radiative forcing)
    real(r8) :: albsni_pur(lbc:ubc,numrad) ! diffuse pure snow albedo (radiative forcing)
    real(r8) :: albsnd_bc(lbc:ubc,numrad)  ! direct snow albedo without BC (radiative forcing)
    real(r8) :: albsni_bc(lbc:ubc,numrad)  ! diffuse snow albedo without BC (radiative forcing)
    real(r8) :: albsnd_oc(lbc:ubc,numrad)  ! direct snow albedo without OC (radiative forcing)
    real(r8) :: albsni_oc(lbc:ubc,numrad)  ! diffuse snow albedo without OC (radiative forcing)
    real(r8) :: albsnd_dst(lbc:ubc,numrad) ! direct snow albedo without dust (radiative forcing)
    real(r8) :: albsni_dst(lbc:ubc,numrad) ! diffuse snow albedo without dust (radiative forcing)
    integer  :: i                          ! index for layers [idx]
    real(r8) :: flx_absd_snw(lbc:ubc,-nlevsno+1:1,numrad)   ! flux absorption factor for just snow (direct) [frc]
    real(r8) :: flx_absi_snw(lbc:ubc,-nlevsno+1:1,numrad)   ! flux absorption factor for just snow (diffuse) [frc]
    real(r8) :: foo_snw(lbc:ubc,-nlevsno+1:1,numrad)        ! dummy array for forcing calls
    real(r8) :: albsfc(lbc:ubc,numrad)                      ! albedo of surface underneath snow (col,bnd) 
    real(r8) :: h2osno_liq(lbc:ubc,-nlevsno+1:0)            ! liquid snow content (col,lyr) [kg m-2]
    real(r8) :: h2osno_ice(lbc:ubc,-nlevsno+1:0)            ! ice content in snow (col,lyr) [kg m-2]
    integer  :: snw_rds_in(lbc:ubc,-nlevsno+1:0)            ! snow grain size sent to SNICAR (col,lyr) [microns]
    real(r8) :: mss_cnc_aer_in_frc_pur(lbc:ubc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_bc(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)  ! mass concentration of aerosol species for BC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_oc(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)  ! mass concentration of aerosol species for OC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_dst(lbc:ubc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for dust forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_fdb(lbc:ubc,-nlevsno+1:0,sno_nbr_aer)     ! mass concentration of all aerosol species for feedback calculation (col,lyr,aer) [kg kg-1]
  !-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    lat => clm3%g%lat
    lon => clm3%g%lon

    ! Assign local pointers to derived subtypes components (landunit level)

    itypelun       => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    cgridcell      => clm3%g%l%c%gridcell
    h2osno         => clm3%g%l%c%cws%h2osno
    albgrd         => clm3%g%l%c%cps%albgrd
    albgri         => clm3%g%l%c%cps%albgri
    decl           => clm3%g%l%c%cps%decl 
    coszen         => clm3%g%l%c%cps%coszen 
    albsod         => clm3%g%l%c%cps%albsod
    albsoi         => clm3%g%l%c%cps%albsoi
    frac_sno       => clm3%g%l%c%cps%frac_sno
    flx_absdv      => clm3%g%l%c%cps%flx_absdv
    flx_absdn      => clm3%g%l%c%cps%flx_absdn
    flx_absiv      => clm3%g%l%c%cps%flx_absiv
    flx_absin      => clm3%g%l%c%cps%flx_absin
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    snw_rds        => clm3%g%l%c%cps%snw_rds
    albgrd_pur     => clm3%g%l%c%cps%albgrd_pur
    albgri_pur     => clm3%g%l%c%cps%albgri_pur
    albgrd_bc      => clm3%g%l%c%cps%albgrd_bc
    albgri_bc      => clm3%g%l%c%cps%albgri_bc
    albgrd_oc      => clm3%g%l%c%cps%albgrd_oc
    albgri_oc      => clm3%g%l%c%cps%albgri_oc
    albgrd_dst     => clm3%g%l%c%cps%albgrd_dst
    albgri_dst     => clm3%g%l%c%cps%albgri_dst
    mss_cnc_bcphi  => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho  => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_cnc_ocphi  => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho  => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_cnc_dst1   => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2   => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3   => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4   => clm3%g%l%c%cps%mss_cnc_dst4
    albsnd_hst     => clm3%g%l%c%cps%albsnd_hst
    albsni_hst     => clm3%g%l%c%cps%albsni_hst

    ! Assign local pointers to derived subtypes components (pft-level)

    pactive   => clm3%g%l%c%p%active
    plandunit => clm3%g%l%c%p%landunit
    pgridcell => clm3%g%l%c%p%gridcell
    pcolumn   => clm3%g%l%c%p%column
    albd      => clm3%g%l%c%p%pps%albd
    albi      => clm3%g%l%c%p%pps%albi
    fabd      => clm3%g%l%c%p%pps%fabd
    fabd_sun  => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha  => clm3%g%l%c%p%pps%fabd_sha
    fabi      => clm3%g%l%c%p%pps%fabi
    fabi_sun  => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha  => clm3%g%l%c%p%pps%fabi_sha
    ftdd      => clm3%g%l%c%p%pps%ftdd
    ftid      => clm3%g%l%c%p%pps%ftid
    ftii      => clm3%g%l%c%p%pps%ftii
    vcmaxcintsun => clm3%g%l%c%p%pps%vcmaxcintsun
    vcmaxcintsha => clm3%g%l%c%p%pps%vcmaxcintsha
    ncan      => clm3%g%l%c%p%pps%ncan
    nrad      => clm3%g%l%c%p%pps%nrad
    fabd_sun_z  => clm3%g%l%c%p%pps%fabd_sun_z
    fabd_sha_z  => clm3%g%l%c%p%pps%fabd_sha_z
    fabi_sun_z  => clm3%g%l%c%p%pps%fabi_sun_z
    fabi_sha_z  => clm3%g%l%c%p%pps%fabi_sha_z
    fsun_z      => clm3%g%l%c%p%pps%fsun_z
    tlai_z    => clm3%g%l%c%p%pps%tlai_z
    tsai_z    => clm3%g%l%c%p%pps%tsai_z
    tlai      => clm3%g%l%c%p%pps%tlai
    tsai      => clm3%g%l%c%p%pps%tsai
    elai      => clm3%g%l%c%p%pps%elai
    esai      => clm3%g%l%c%p%pps%esai
    ivt       => clm3%g%l%c%p%itype
    rhol      => pftcon%rhol
    rhos      => pftcon%rhos
    taul      => pftcon%taul
    taus      => pftcon%taus
    
    ! Cosine solar zenith angle for next time step

    do g = lbg, ubg
       coszen_gcell(g) = shr_orb_cosz (nextsw_cday, lat(g), lon(g), declinp1)
    end do

    ! Save coszen and declination values to  clm3 data structures for
    ! use in other places in the CN and urban code

    do c = lbc,ubc
       g = cgridcell(c)
       coszen_col(c) = coszen_gcell(g)
       coszen(c) = coszen_col(c)
       decl(c) = declinp1
    end do

    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
!      if (pwtgcell(p)>0._r8) then ! "if" added due to chg in filter definition
       g = pgridcell(p)
       coszen_pft(p) = coszen_gcell(g)
!      end if ! then removed for CNDV (and dyn. landuse?) cases to work
    end do

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          albsod(c,ib)     = 0._r8
          albsoi(c,ib)     = 0._r8
          albgrd(c,ib)     = 0._r8
          albgri(c,ib)     = 0._r8
          albgrd_pur(c,ib) = 0._r8
          albgri_pur(c,ib) = 0._r8
          albgrd_bc(c,ib)  = 0._r8
          albgri_bc(c,ib)  = 0._r8
          albgrd_oc(c,ib)  = 0._r8
          albgri_oc(c,ib)  = 0._r8
          albgrd_dst(c,ib) = 0._r8
          albgri_dst(c,ib) = 0._r8
          do i=-nlevsno+1,1,1
             flx_absdv(c,i) = 0._r8
             flx_absdn(c,i) = 0._r8
             flx_absiv(c,i) = 0._r8
             flx_absin(c,i) = 0._r8
          enddo
       end do
       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
!         if (pwtgcell(p)>0._r8) then ! "if" added due to chg in filter definition
          albd(p,ib) = 1._r8
          albi(p,ib) = 1._r8
          fabd(p,ib) = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib) = 0._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 0._r8
!         if (ib==1) then
!            ncan(p) = 0
!            nrad(p) = 0
!            do iv = 1, nlevcan
!               fabd_sun_z(p,iv) = 0._r8
!               fabd_sha_z(p,iv) = 0._r8
!               fabi_sun_z(p,iv) = 0._r8
!               fabi_sha_z(p,iv) = 0._r8
!               fsun_z(p,iv) = 0._r8
!               tlai_z(p,iv) = 0._r8
!               tsai_z(p,iv) = 0._r8
!            end do
!         end if
!         end if ! then removed for CNDV (and dyn. landuse?) cases to work
       end do
    end do

    ! SoilAlbedo called before SNICAR_RT
    ! so that reflectance of soil beneath snow column is known 
    ! ahead of time for snow RT calculation.

    ! Snow albedos
    ! Note that snow albedo routine will only compute nonzero snow albedos
    ! where h2osno> 0 and coszen > 0
    
    ! Ground surface albedos
    ! Note that ground albedo routine will only compute nonzero snow albedos
    ! where coszen > 0

    call SoilAlbedo(lbc, ubc, num_nourbanc, filter_nourbanc, &
                    coszen_col, albsnd, albsni) 

    ! set variables to pass to SNICAR.
    
    flg_snw_ice = 1   ! calling from CLM, not CSIM
    do c=lbc,ubc
       albsfc(c,:)     = albsoi(c,:)
       h2osno_liq(c,:) = h2osoi_liq(c,-nlevsno+1:0)
       h2osno_ice(c,:) = h2osoi_ice(c,-nlevsno+1:0)
       snw_rds_in(c,:) = nint(snw_rds(c,:))

       ! zero aerosol input arrays
       mss_cnc_aer_in_frc_pur(c,:,:) = 0._r8
       mss_cnc_aer_in_frc_bc(c,:,:)  = 0._r8
       mss_cnc_aer_in_frc_oc(c,:,:)  = 0._r8
       mss_cnc_aer_in_frc_dst(c,:,:) = 0._r8
       mss_cnc_aer_in_fdb(c,:,:)     = 0._r8
    end do

    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if (DO_SNO_AER) then
       mss_cnc_aer_in_fdb(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
       mss_cnc_aer_in_fdb(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)
       
       ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
       !  1) Knowledge of their optical properties is primitive
       !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow, 
       !     it has a negligible darkening effect.
       if (DO_SNO_OC) then
          mss_cnc_aer_in_fdb(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
          mss_cnc_aer_in_fdb(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
       endif
       
       mss_cnc_aer_in_fdb(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
       mss_cnc_aer_in_fdb(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
       mss_cnc_aer_in_fdb(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
       mss_cnc_aer_in_fdb(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)
    endif

! If radiative forcing is being calculated, first estimate clean-snow albedo
#if (defined SNICAR_FRC)

    ! 1. BC input array:
    !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
    mss_cnc_aer_in_frc_bc(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)
    if (DO_SNO_OC) then
       mss_cnc_aer_in_frc_bc(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
       mss_cnc_aer_in_frc_bc(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
    endif

    ! BC FORCING CALCULATIONS
    flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_bc, albsfc, albsnd_bc, foo_snw)
    
    flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_bc, albsfc, albsni_bc, foo_snw)

    ! 2. OC input array:
    !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
    if (DO_SNO_OC) then
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,5) = mss_cnc_dst1(lbc:ubc,:)
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,6) = mss_cnc_dst2(lbc:ubc,:)
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,7) = mss_cnc_dst3(lbc:ubc,:)
       mss_cnc_aer_in_frc_oc(lbc:ubc,:,8) = mss_cnc_dst4(lbc:ubc,:)
    
       ! OC FORCING CALCULATIONS
       flg_slr = 1; ! direct-beam
       call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                      coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                      mss_cnc_aer_in_frc_oc, albsfc, albsnd_oc, foo_snw)
    
       flg_slr = 2; ! diffuse
       call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                      coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                      mss_cnc_aer_in_frc_oc, albsfc, albsni_oc, foo_snw)
    endif
    
    ! 3. DUST input array:
    ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
    mss_cnc_aer_in_frc_dst(lbc:ubc,:,1) = mss_cnc_bcphi(lbc:ubc,:)
    mss_cnc_aer_in_frc_dst(lbc:ubc,:,2) = mss_cnc_bcpho(lbc:ubc,:)
    if (DO_SNO_OC) then
       mss_cnc_aer_in_frc_dst(lbc:ubc,:,3) = mss_cnc_ocphi(lbc:ubc,:)
       mss_cnc_aer_in_frc_dst(lbc:ubc,:,4) = mss_cnc_ocpho(lbc:ubc,:)
    endif
    
    ! DUST FORCING CALCULATIONS
    flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_dst, albsfc, albsnd_dst, foo_snw)
    
    flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_dst, albsfc, albsni_dst, foo_snw)

    ! 4. ALL AEROSOL FORCING CALCULATION
    ! (pure snow albedo)
    flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_pur, albsfc, albsnd_pur, foo_snw)
    
    flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_frc_pur, albsfc, albsni_pur, foo_snw)

#endif

    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
    flg_slr = 1; ! direct-beam
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_fdb, albsfc, albsnd, flx_absd_snw)

    flg_slr = 2; ! diffuse
    call SNICAR_RT(flg_snw_ice, lbc, ubc, num_nourbanc, filter_nourbanc,    &
                   coszen_col, flg_slr, h2osno_liq, h2osno_ice, snw_rds_in, &
                   mss_cnc_aer_in_fdb, albsfc, albsni, flx_absi_snw)

    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          if (coszen(c) > 0._r8) then
             ! ground albedo was originally computed in SoilAlbedo, but is now computed here
             ! because the order of SoilAlbedo and SNICAR_RT was switched for SNICAR.
             albgrd(c,ib) = albsod(c,ib)*(1._r8-frac_sno(c)) + albsnd(c,ib)*frac_sno(c)
             albgri(c,ib) = albsoi(c,ib)*(1._r8-frac_sno(c)) + albsni(c,ib)*frac_sno(c)

             ! albedos for radiative forcing calculations:
#if (defined SNICAR_FRC)
             ! BC forcing albedo
             albgrd_bc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_bc(c,ib)*frac_sno(c)
             albgri_bc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_bc(c,ib)*frac_sno(c)
             
             if (DO_SNO_OC) then
                ! OC forcing albedo
                albgrd_oc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_oc(c,ib)*frac_sno(c)
                albgri_oc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_oc(c,ib)*frac_sno(c)
             endif

             ! dust forcing albedo
             albgrd_dst(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_dst(c,ib)*frac_sno(c)
             albgri_dst(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_dst(c,ib)*frac_sno(c)

             ! pure snow albedo for all-aerosol radiative forcing
             albgrd_pur(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_pur(c,ib)*frac_sno(c)
             albgri_pur(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_pur(c,ib)*frac_sno(c)
#endif

             ! also in this loop (but optionally in a different loop for vectorized code)
             !  weight snow layer radiative absorption factors based on snow fraction and soil albedo
             !  (NEEDED FOR ENERGY CONSERVATION)
             do i = -nlevsno+1,1,1
              if (subgridflag == 0) then 
                if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                endif
             else
                if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                endif
             endif
             enddo
          endif
       enddo
    enddo

    ! for diagnostics, set snow albedo to spval over non-snow points 
    ! so that it is not averaged in history buffer
    ! (OPTIONAL)
    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          if ((coszen(c) > 0._r8) .and. (h2osno(c) > 0._r8)) then
             albsnd_hst(c,ib) = albsnd(c,ib)
             albsni_hst(c,ib) = albsni(c,ib)
          else
             albsnd_hst(c,ib) = 0._r8
             albsni_hst(c,ib) = 0._r8
          endif
       enddo
    enddo

    ! Create solar-vegetated filter for the following calculations

    num_vegsol = 0
    num_novegsol = 0
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
          if (coszen_pft(p) > 0._r8) then
             if ((itypelun(plandunit(p)) == istsoil .or.  &
                  itypelun(plandunit(p)) == istcrop     ) &
                 .and. (elai(p) + esai(p)) > 0._r8        &
                 .and. pactive(p)) then
                num_vegsol = num_vegsol + 1
                filter_vegsol(num_vegsol) = p
             else
                num_novegsol = num_novegsol + 1
                filter_novegsol(num_novegsol) = p
             end if
          end if
    end do

    ! Weight reflectance/transmittance by lai and sai
    ! Only perform on vegetated pfts where coszen > 0

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       wl(p) = elai(p) / max( elai(p)+esai(p), mpe )
       ws(p) = esai(p) / max( elai(p)+esai(p), mpe )
    end do

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          rho(p,ib) = max( rhol(ivt(p),ib)*wl(p) + rhos(ivt(p),ib)*ws(p), mpe )
          tau(p,ib) = max( taul(ivt(p),ib)*wl(p) + taus(ivt(p),ib)*ws(p), mpe )
       end do
    end do

    ! Diagnose number of canopy layers for radiative transfer, in increments of dincmax.
    ! Add to number of layers so long as cumulative leaf+stem area does not exceed total
    ! leaf+stem area. Then add any remaining leaf+stem area to next layer and exit the loop.
    ! Do this first for elai and esai (not buried by snow) and then for the part of the
    ! canopy that is buried by snow. 
    ! ------------------
    ! tlai_z = leaf area increment for a layer
    ! tsai_z = stem area increment for a layer
    ! nrad   = number of canopy layers above snow
    ! ncan   = total number of canopy layers
    !
    ! tlai_z summed from 1 to nrad = elai
    ! tlai_z summed from 1 to ncan = tlai

    ! tsai_z summed from 1 to nrad = esai
    ! tsai_z summed from 1 to ncan = tsai
    ! ------------------
    !
    ! Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol" 
    ! because layering is needed for all time steps regardless of radiation
    !
    ! Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered by
    ! nlevcan = 1

    dincmax = 0.25_r8
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)

       if (nlevcan == 1) then
          nrad(p) = 1
          ncan(p) = 1
          tlai_z(p,1) = elai(p)
          tsai_z(p,1) = esai(p)
       else if (nlevcan > 1) then
          if (elai(p)+esai(p) == 0._r8) then
             nrad(p) = 0
          else
             dincmax_sum = 0._r8
             do iv = 1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((elai(p)+esai(p))-dincmax_sum) > 1.e-06_r8) then
                   nrad(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                else
                   nrad(p) = iv
                   dinc = dincmax - (dincmax_sum - (elai(p)+esai(p)))
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                   exit
                end if
             end do

             ! Mimumum of 4 canopy layers

             if (nrad(p) < 4) then
                nrad(p) = 4
                do iv = 1, nrad(p)
                   tlai_z(p,iv) = elai(p) / nrad(p)
                   tsai_z(p,iv) = esai(p) / nrad(p)
                end do
             end if
          end if
       end if

       ! Error check: make sure cumulative of increments does not exceed total

       laisum = 0._r8
       saisum = 0._r8
       do iv = 1, nrad(p)
          laisum = laisum + tlai_z(p,iv)
          saisum = saisum + tsai_z(p,iv)
       end do
       if (abs(laisum-elai(p)) > 1.e-06_r8 .or. abs(saisum-esai(p)) > 1.e-06_r8) then
          write (iulog,*) 'multi-layer canopy error 01 in SurfaceAlbedo: ',nrad(p),elai(p),laisum,esai(p),saisum
          call endrun()
       end if

       ! Repeat to find canopy layers buried by snow

       if (nlevcan > 1) then
          blai(p) = tlai(p) - elai(p)
          bsai(p) = tsai(p) - esai(p)
          if (blai(p)+bsai(p) == 0._r8) then
             ncan(p) = nrad(p)
          else
             dincmax_sum = 0._r8
             do iv = nrad(p)+1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((blai(p)+bsai(p))-dincmax_sum) > 1.e-06_r8) then
                   ncan(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                else
                   ncan(p) = iv
                   dinc = dincmax - (dincmax_sum - (blai(p)+bsai(p)))
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                   exit
                end if
             end do
          end if

          ! Error check: make sure cumulative of increments does not exceed total

          laisum = 0._r8
          saisum = 0._r8
          do iv = 1, ncan(p)
             laisum = laisum + tlai_z(p,iv)
             saisum = saisum + tsai_z(p,iv)
          end do
          if (abs(laisum-tlai(p)) > 1.e-06_r8 .or. abs(saisum-tsai(p)) > 1.e-06_r8) then
             write (iulog,*) 'multi-layer canopy error 02 in SurfaceAlbedo: ',nrad(p),ncan(p)
             write (iulog,*) tlai(p),elai(p),blai(p),laisum,tsai(p),esai(p),bsai(p),saisum
             call endrun()
          end if
       end if

    end do

    ! Zero fluxes for active canopy layers

    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       do iv = 1, nrad(p)
          fabd_sun_z(p,iv) = 0._r8
          fabd_sha_z(p,iv) = 0._r8
          fabi_sun_z(p,iv) = 0._r8
          fabi_sha_z(p,iv) = 0._r8
          fsun_z(p,iv) = 0._r8
       end do
    end do

    ! Default leaf to canopy scaling coefficients, used when coszen <= 0.
    ! This is the leaf nitrogen profile integrated over the full canopy.
    ! Integrate exp(-kn*x) over x=0 to x=elai and assign to shaded canopy,
    ! because sunlit fraction is 0. Canopy scaling coefficients are set in
    ! TwoStream for coszen > 0. So kn must be set here and in TwoStream.

    extkn = 0.30_r8
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       if (nlevcan == 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn
          if (elai(p) > 0._r8) then
             vcmaxcintsha(p) = vcmaxcintsha(p) / elai(p)
          else
             vcmaxcintsha(p) = 0._r8
          end if
       else if (nlevcan > 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = 0._r8
       end if
    end do

    ! Calculate surface albedos and fluxes 
    ! Only perform on vegetated pfts where coszen > 0

    call TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, num_vegsol, coszen_pft, rho, tau)
       
    ! Determine values for non-vegetated pfts where coszen > 0

    do ib = 1,numrad
       do fp = 1,num_novegsol
          p = filter_novegsol(fp)
          c = pcolumn(p)
          fabd(p,ib) = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib) = 1._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 1._r8
          albd(p,ib) = albgrd(c,ib)
          albi(p,ib) = albgri(c,ib)
!         if (ib == 1) then
!            do iv = 1, nrad(p)
!               fabd_sun_z(p,iv) = 0._r8
!               fabd_sha_z(p,iv) = 0._r8
!               fabi_sun_z(p,iv) = 0._r8
!               fabi_sha_z(p,iv) = 0._r8
!               fsun_z(p,iv) = 0._r8
!            end do
!         end if
       end do
    end do

  end subroutine SurfaceAlbedo


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilAlbedo
!
! !INTERFACE:
  subroutine SoilAlbedo (lbc, ubc, num_nourbanc, filter_nourbanc, coszen, albsnd, albsni)
!
! !DESCRIPTION:
! Determine ground surface albedo, accounting for snow
!
! !USES:
    use clmtype
    use clm_varpar, only : numrad
    use clm_varcon, only : albsat, albdry, tfrz, istice, istice_mec
    use clm_varcon, only : istdlak
    use SLakeCon  , only : alblak
    use SLakeCon  , only : alblakwi, calb, lakepuddling
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                   ! column bounds
    integer , intent(in) :: num_nourbanc               ! number of columns in non-urban points in column filter
    integer , intent(in) :: filter_nourbanc(ubc-lbc+1) ! column filter for non-urban points
    real(r8), intent(in) :: coszen(lbc:ubc)            ! cos solar zenith angle next time step (column-level)
    real(r8), intent(in) :: albsnd(lbc:ubc,numrad)     ! snow albedo (direct)
    real(r8), intent(in) :: albsni(lbc:ubc,numrad)     ! snow albedo (diffuse)
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/5/02, Peter Thornton: Migrated to new data structures.
! 8/20/03, Mariana Vertenstein: Vectorized routine
! 03/28/08, Mark Flanner: changes for SNICAR
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: clandunit(:)    ! landunit of corresponding column
    integer , pointer :: ltype(:)        ! landunit type
    integer , pointer :: isoicol(:)      ! soil color class
    real(r8), pointer :: t_grnd(:)       ! ground temperature (Kelvin)
    real(r8), pointer :: h2osoi_vol(:,:) ! volumetric soil water [m3/m3]
    real(r8), pointer :: lake_icefrac(:,:)  ! mass fraction of lake layer that is frozen
    integer , pointer :: snl(:)             ! number of snow layers
!
! local pointers to original implicit out arguments
!
    real(r8), pointer:: albgrd(:,:)      ! ground albedo (direct)
    real(r8), pointer:: albgri(:,:)      ! ground albedo (diffuse)
    ! albsod and albsoi are now clm_type variables so they can be used by SNICAR.
    real(r8), pointer :: albsod(:,:)        ! soil albedo (direct)
    real(r8), pointer :: albsoi(:,:)        ! soil albedo (diffuse)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer, parameter :: nband =numrad ! number of solar radiation waveband classes
    integer  :: fc            ! non-urban filter column index
    integer  :: c,l           ! indices
    integer  :: ib            ! waveband number (1=vis, 2=nir)
    real(r8) :: inc           ! soil water correction factor for soil albedo
    ! albsod and albsoi are now clm_type variables so they can be used by SNICAR.
    !real(r8) :: albsod        ! soil albedo (direct)
    !real(r8) :: albsoi        ! soil albedo (diffuse)
    integer  :: soilcol       ! soilcolor
    real(r8) :: sicefr        ! Lake surface ice fraction (based on D. Mironov 2010)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit  => clm3%g%l%c%landunit
    isoicol    => clm3%g%l%c%cps%isoicol
    t_grnd     => clm3%g%l%c%ces%t_grnd
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
    albgrd     => clm3%g%l%c%cps%albgrd
    albgri     => clm3%g%l%c%cps%albgri
    albsod     => clm3%g%l%c%cps%albsod
    albsoi     => clm3%g%l%c%cps%albsoi
    snl        => clm3%g%l%c%cps%snl
    lake_icefrac => clm3%g%l%c%cws%lake_icefrac

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype

    ! Compute soil albedos

    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          if (coszen(c) > 0._r8) then
             l = clandunit(c)

             if (ltype(l) == istsoil .or. ltype(l) == istcrop)  then ! soil
                inc    = max(0.11_r8-0.40_r8*h2osoi_vol(c,1), 0._r8)
                soilcol = isoicol(c)
                ! changed from local variable to clm_type:
                !albsod = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                !albsoi = albsod
                albsod(c,ib) = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                albsoi(c,ib) = albsod(c,ib)
             else if (ltype(l) == istice .or. ltype(l) == istice_mec)  then  ! land ice
                ! changed from local variable to clm_type:
                !albsod = albice(ib)
                !albsoi = albsod
                albsod(c,ib) = albice(ib)
                albsoi(c,ib) = albsod(c,ib)
             ! unfrozen lake, wetland
             else if (t_grnd(c) > tfrz .or. (lakepuddling .and. ltype(l) == istdlak .and. t_grnd(c) == tfrz .and. &
                      lake_icefrac(c,1) < 1._r8 .and. lake_icefrac(c,2) > 0._r8) ) then

                albsod(c,ib) = 0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8)
                ! This expression is apparently from BATS according to Yongjiu Dai.

                ! The diffuse albedo should be an average over the whole sky of an angular-dependent direct expression.
                ! The expression above may have been derived to encompass both (e.g. Henderson-Sellers 1986),
                ! but I'll assume it applies more appropriately to the direct form for now.

                ! ZMS: Attn EK, currently restoring this for wetlands even though it is wrong in order to try to get
                ! bfb baseline comparison when no lakes are present. I'm assuming wetlands will be phased out anyway.
                if (ltype(l) == istdlak) then
                   albsoi(c,ib) = 0.10_r8
                else
                   albsoi(c,ib) = albsod(c,ib)
                end if

             else                                     ! frozen lake, wetland
                ! Introduce crude surface frozen fraction according to D. Mironov (2010)
                ! Attn EK: This formulation is probably just as good for "wetlands" if they are not phased out.
                ! Tenatively I'm restricting this to lakes because I haven't tested it for wetlands. But if anything
                ! the albedo should be lower when melting over frozen ground than a solid frozen lake.
                ! 
                if (ltype(l) == istdlak .and. .not. lakepuddling .and. snl(c) == 0) then
                    ! Need to reference snow layers here because t_grnd could be over snow or ice
                                      ! but we really want the ice surface temperature with no snow
                   sicefr = 1._r8 - exp(-calb * (tfrz - t_grnd(c))/tfrz)
                   albsod(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), &
                                  0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8))
                   albsoi(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), 0.10_r8)
                   ! Make sure this is no less than the open water albedo above.
                   ! Setting lake_melt_icealb(:) = alblak(:) in namelist reverts the melting albedo to the cold
                   ! snow-free value.
                else
                   albsod(c,ib) = alblak(ib)
                   albsoi(c,ib) = albsod(c,ib)
                end if
             end if

             ! Weighting is done in SurfaceAlbedo, after the call to SNICAR_RT
             ! This had to be done, because SoilAlbedo is called before SNICAR_RT, so at
             ! this point, snow albedo is not yet known.
          end if
       end do
    end do

  end subroutine SoilAlbedo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: TwoStream
!
! !INTERFACE:
  subroutine TwoStream (lbc, ubc, lbp, ubp, filter_vegsol, num_vegsol, coszen, rho, tau)
!
! !DESCRIPTION:
! Two-stream fluxes for canopy radiative transfer
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
! Calculate sunlit and shaded fluxes as described by
! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
! a multi-layer canopy to calculate APAR profile 
!
! !USES:
    use clmtype
    use clm_varpar, only : numrad, nlevcan
    use clm_varcon, only : omegas, tfrz, betads, betais
    use clm_varctl, only : iulog
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                 ! column bounds
    integer , intent(in)  :: lbp, ubp                 ! pft bounds
    integer , intent(in)  :: filter_vegsol(ubp-lbp+1) ! filter for vegetated pfts with coszen>0
    integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0
    real(r8), intent(in)  :: coszen(lbp:ubp)          ! cosine solar zenith angle for next time step
    real(r8), intent(in)  :: rho(lbp:ubp,numrad)      ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8), intent(in)  :: tau(lbp:ubp,numrad)      ! leaf/stem tran weighted by fraction LAI and SAI
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! Modified for speedup: Mariana Vertenstein, 8/26/02
! Vectorized routine: Mariana Vertenstein:  8/20/03
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    real(r8), pointer :: albgrd(:,:)   ! ground albedo (direct) (column-level)
    real(r8), pointer :: albgri(:,:)   ! ground albedo (diffuse)(column-level)
    real(r8), pointer :: t_veg(:)      ! vegetation temperature (Kelvin)
    real(r8), pointer :: fwet(:)       ! fraction of canopy that is wet (0 to 1)
    integer , pointer :: ivt(:)        ! pft vegetation type
    real(r8), pointer :: xl(:)         ! ecophys const - leaf/stem orientation index
    real(r8), pointer :: elai(:)       ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)       ! one-sided stem area index with burying by snow
    integer , pointer :: nrad(:)       ! number of canopy layers, above snow for radiative transfer
    real(r8), pointer :: tlai_z(:,:)   ! tlai increment for canopy layer
    real(r8), pointer :: tsai_z(:,:)   ! tsai increment for canopy layer
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: albd(:,:)     ! surface albedo (direct)
    real(r8), pointer :: albi(:,:)     ! surface albedo (diffuse)
    real(r8), pointer :: fabd(:,:)     ! flux absorbed by canopy per unit direct flux
    real(r8), pointer :: fabd_sun(:,:) ! flux absorbed by sunlit canopy per unit direct flux
    real(r8), pointer :: fabd_sha(:,:) ! flux absorbed by shaded canopy per unit direct flux
    real(r8), pointer :: fabi(:,:)     ! flux absorbed by canopy per unit diffuse flux
    real(r8), pointer :: fabi_sun(:,:) ! flux absorbed by sunlit canopy per unit diffuse flux
    real(r8), pointer :: fabi_sha(:,:) ! flux absorbed by shaded canopy per unit diffuse flux
    real(r8), pointer :: ftdd(:,:)     ! down direct flux below canopy per unit direct flx
    real(r8), pointer :: ftid(:,:)     ! down diffuse flux below canopy per unit direct flx
    real(r8), pointer :: ftii(:,:)     ! down diffuse flux below canopy per unit diffuse flx
    real(r8), pointer :: fabd_sun_z(:,:)! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabd_sha_z(:,:)! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabi_sun_z(:,:)! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fabi_sha_z(:,:)! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
    real(r8), pointer :: fsun_z(:,:)    ! sunlit fraction of canopy layer
    real(r8), pointer :: vcmaxcintsun(:) ! leaf to canopy scaling coefficient, sunlit leaf vcmax
    real(r8), pointer :: vcmaxcintsha(:) ! leaf to canopy scaling coefficient, shaded leaf vcmax
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: fp,p,c,iv        ! array indices
    integer  :: ib               ! waveband number
    real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
    real(r8) :: asu              ! single scattering albedo
    real(r8) :: chil(lbp:ubp)    ! -0.4 <= xl <= 0.6
    real(r8) :: gdir(lbp:ubp)    ! leaf projection in solar direction (0 to 1)
    real(r8) :: twostext(lbp:ubp)! optical depth of direct beam per unit leaf area
    real(r8) :: avmu(lbp:ubp)    ! average diffuse optical depth
    real(r8) :: omega(lbp:ubp,numrad)   ! fraction of intercepted radiation that is scattered (0 to 1)
    real(r8) :: omegal           ! omega for leaves
    real(r8) :: betai            ! upscatter parameter for diffuse radiation
    real(r8) :: betail           ! betai for leaves
    real(r8) :: betad            ! upscatter parameter for direct beam radiation
    real(r8) :: betadl           ! betad for leaves
    real(r8) :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
    real(r8) :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
    real(r8) :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
    real(r8) :: phi1,phi2,sigma                                   ! temporary
    real(r8) :: temp0(lbp:ubp),temp1,temp2(lbp:ubp)               ! temporary
    real(r8) :: t1                                                ! temporary
    real(r8) :: a1,a2                                             ! parameter for sunlit/shaded leaf radiation absorption
    real(r8) :: v,dv,u,du                                         ! temporary for flux derivatives
    real(r8) :: dh2,dh3,dh5,dh6,dh7,dh8,dh9,dh10                  ! temporary for flux derivatives
    real(r8) :: da1,da2                                           ! temporary for flux derivatives
    real(r8) :: d_ftid,d_ftii                                     ! ftid, ftii derivative with respect to lai+sai
    real(r8) :: d_fabd,d_fabi                                     ! fabd, fabi derivative with respect to lai+sai
    real(r8) :: d_fabd_sun,d_fabd_sha                             ! fabd_sun, fabd_sha derivative with respect to lai+sai
    real(r8) :: d_fabi_sun,d_fabi_sha                             ! fabi_sun, fabi_sha derivative with respect to lai+sai
    real(r8) :: laisum                                            ! cumulative lai+sai for canopy layer (at middle of layer)
    real(r8) :: extkb                                             ! direct beam extinction coefficient
    real(r8) :: extkn                                             ! nitrogen allocation coefficient
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    albgrd  => clm3%g%l%c%cps%albgrd
    albgri  => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn  => clm3%g%l%c%p%column
    fwet     => clm3%g%l%c%p%pps%fwet
    t_veg    => clm3%g%l%c%p%pes%t_veg
    ivt      => clm3%g%l%c%p%itype
    elai     => clm3%g%l%c%p%pps%elai
    esai     => clm3%g%l%c%p%pps%esai
    albd     => clm3%g%l%c%p%pps%albd
    albi     => clm3%g%l%c%p%pps%albi
    fabd     => clm3%g%l%c%p%pps%fabd
    fabd_sun => clm3%g%l%c%p%pps%fabd_sun
    fabd_sha => clm3%g%l%c%p%pps%fabd_sha
    fabi     => clm3%g%l%c%p%pps%fabi
    fabi_sun => clm3%g%l%c%p%pps%fabi_sun
    fabi_sha => clm3%g%l%c%p%pps%fabi_sha
    ftdd     => clm3%g%l%c%p%pps%ftdd
    ftid     => clm3%g%l%c%p%pps%ftid
    ftii     => clm3%g%l%c%p%pps%ftii
    xl       => pftcon%xl
    nrad     => clm3%g%l%c%p%pps%nrad
    fabd_sun_z  => clm3%g%l%c%p%pps%fabd_sun_z
    fabd_sha_z  => clm3%g%l%c%p%pps%fabd_sha_z
    fabi_sun_z  => clm3%g%l%c%p%pps%fabi_sun_z
    fabi_sha_z  => clm3%g%l%c%p%pps%fabi_sha_z
    fsun_z      => clm3%g%l%c%p%pps%fsun_z
    tlai_z   => clm3%g%l%c%p%pps%tlai_z
    tsai_z   => clm3%g%l%c%p%pps%tsai_z
    vcmaxcintsun => clm3%g%l%c%p%pps%vcmaxcintsun
    vcmaxcintsha => clm3%g%l%c%p%pps%vcmaxcintsha

    ! Calculate two-stream parameters that are independent of waveband:
    ! chil, gdir, twostext, avmu, and temp0 and temp2 (used for asu)

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       
       ! note that the following limit only acts on cosz values > 0 and less than 
       ! 0.001, not on values cosz = 0, since these zero have already been filtered
       ! out in filter_vegsol
       cosz = max(0.001_r8, coszen(p))

       chil(p) = min( max(xl(ivt(p)), -0.4_r8), 0.6_r8 )
       if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8
       phi1 = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
       phi2 = 0.877_r8 * (1._r8-2._r8*phi1)
       gdir(p) = phi1 + phi2*cosz
       twostext(p) = gdir(p)/cosz
       avmu(p) = ( 1._r8 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
       temp0(p) = gdir(p) + phi2*cosz
       temp1 = phi1*cosz
       temp2(p) = ( 1._r8 - temp1/temp0(p) * log((temp1+temp0(p))/temp1) )
    end do

   ! Loop over all wavebands to calculate for the full canopy the scattered fluxes 
   ! reflected upward and transmitted downward by the canopy and the flux absorbed by the 
   ! canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given 
   ! an underlying surface of known albedo.
   !
   ! Output:
   ! ------------------
   ! Direct beam fluxes
   ! ------------------
   ! albd       - Upward scattered flux above canopy (per unit direct beam flux)
   ! ftid       - Downward scattered flux below canopy (per unit direct beam flux)
   ! ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
   ! fabd       - Flux absorbed by canopy (per unit direct beam flux)
   ! fabd_sun   - Sunlit portion of fabd
   ! fabd_sha   - Shaded portion of fabd
   ! fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
   ! ------------------
   ! Diffuse fluxes
   ! ------------------
   ! albi       - Upward scattered flux above canopy (per unit diffuse flux)
   ! ftii       - Downward scattered flux below canopy (per unit diffuse flux)
   ! fabi       - Flux absorbed by canopy (per unit diffuse flux)
   ! fabi_sun   - Sunlit portion of fabi
   ! fabi_sha   - Shaded portion of fabi
   ! fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          c = pcolumn(p)

          ! Calculate two-stream parameters omega, betad, and betai.
          ! Omega, betad, betai are adjusted for snow. Values for omega*betad
          ! and omega*betai are calculated and then divided by the new omega
          ! because the product omega*betai, omega*betad is used in solution.
          ! Also, the transmittances and reflectances (tau, rho) are linear
          ! weights of leaf and stem values.

          omegal = rho(p,ib) + tau(p,ib)
          asu = 0.5_r8*omegal*gdir(p)/temp0(p) *temp2(p)
          betadl = (1._r8+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu
          betail = 0.5_r8 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) &
                 * ((1._r8+chil(p))/2._r8)**2) / omegal

          ! Adjust omega, betad, and betai for intercepted snow

          if (t_veg(p) > tfrz) then                             !no snow
             tmp0 = omegal
             tmp1 = betadl
             tmp2 = betail
          else
             tmp0 =   (1._r8-fwet(p))*omegal        + fwet(p)*omegas(ib)
             tmp1 = ( (1._r8-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
             tmp2 = ( (1._r8-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
          end if
          omega(p,ib) = tmp0           
          betad = tmp1 
          betai = tmp2  

          ! Common terms

          b = 1._r8 - omega(p,ib) + omega(p,ib)*betai
          c1 = omega(p,ib)*betai
          tmp0 = avmu(p)*twostext(p)
          d = tmp0 * omega(p,ib)*betad
          f = tmp0 * omega(p,ib)*(1._r8-betad)
          tmp1 = b*b - c1*c1
          h = sqrt(tmp1) / avmu(p)
          sigma = tmp0*tmp0 - tmp1
          p1 = b + avmu(p)*h
          p2 = b - avmu(p)*h
          p3 = b + tmp0
          p4 = b - tmp0

          ! Absorbed, reflected, transmitted fluxes per unit incoming radiation
          ! for full canopy

          t1 = min(h*(elai(p)+esai(p)), 40._r8)
          s1 = exp(-t1)
          t1 = min(twostext(p)*(elai(p)+esai(p)), 40._r8)
          s2 = exp(-t1)
          
          ! Direct beam 

          u1 = b - c1/albgrd(c,ib)
          u2 = b - c1*albgrd(c,ib)
          u3 = f + c1*albgrd(c,ib)
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

          albd(p,ib) = h1/sigma + h2 + h3
          ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1
          ftdd(p,ib) = s2
          fabd(p,ib) = 1._r8 - albd(p,ib) - (1._r8-albgrd(c,ib))*ftdd(p,ib) - (1._r8-albgri(c,ib))*ftid(p,ib)

          a1 = h1 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h2         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h3         * (1._r8 - s2/s1) / (twostext(p) - h)

          a2 = h4 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h5         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h6         * (1._r8 - s2/s1) / (twostext(p) - h)

          fabd_sun(p,ib) = (1._r8 - omega(p,ib)) * ( 1._r8 - s2 + 1._r8 / avmu(p) * (a1 + a2) )
          fabd_sha(p,ib) = fabd(p,ib) - fabd_sun(p,ib)

          ! Diffuse

          u1 = b - c1/albgri(c,ib)
          u2 = b - c1*albgri(c,ib)
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

          albi(p,ib) = h7 + h8
          ftii(p,ib) = h9*s1 + h10/s1
          fabi(p,ib) = 1._r8 - albi(p,ib) - (1._r8-albgri(c,ib))*ftii(p,ib)

          a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
          a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)

          fabi_sun(p,ib) = (1._r8 - omega(p,ib)) / avmu(p) * (a1 + a2)
          fabi_sha(p,ib) = fabi(p,ib) - fabi_sun(p,ib)

          ! Repeat two-stream calculations for each canopy layer to calculate derivatives. 
          ! tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are 
          ! calculated at the center of the layer. Derivatives are needed only for the 
          ! visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
          ! Derivatives are calculated first per unit lai+sai and then normalized for sunlit 
          ! or shaded fraction of canopy layer.

          ! Sun/shade big leaf code uses only one layer, with canopy integrated values from above
          ! and also canopy-integrated scaling coefficients

          if (ib == 1) then
             if (nlevcan == 1) then

                ! sunlit fraction of canopy
                fsun_z(p,1) = (1._r8 - s2) / t1

                ! absorbed PAR (per unit sun/shade lai+sai)
                laisum = elai(p)+esai(p)
                fabd_sun_z(p,1) = fabd_sun(p,ib) / (fsun_z(p,1)*laisum)
                fabi_sun_z(p,1) = fabi_sun(p,ib) / (fsun_z(p,1)*laisum)
                fabd_sha_z(p,1) = fabd_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)
                fabi_sha_z(p,1) = fabi_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)

                ! leaf to canopy scaling coefficients
                extkn = 0.30_r8
                extkb = twostext(p)
                vcmaxcintsun(p) = (1._r8 - exp(-(extkn+extkb)*elai(p))) / (extkn + extkb)
                vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn - vcmaxcintsun(p)
                if (elai(p) .gt. 0._r8) then
                  vcmaxcintsun(p) = vcmaxcintsun(p) / (fsun_z(p,1)*elai(p))
                  vcmaxcintsha(p) = vcmaxcintsha(p) / ((1._r8 - fsun_z(p,1))*elai(p))
                else
                  vcmaxcintsun(p) = 0._r8
                  vcmaxcintsha(p) = 0._r8
                end if

             else if (nlevcan > 1) then
                do iv = 1, nrad(p)

                ! Cumulative lai+sai at center of layer

                if (iv == 1) then
                   laisum = 0.5_r8 * (tlai_z(p,iv)+tsai_z(p,iv))
                else
                   laisum = laisum + 0.5_r8 * ((tlai_z(p,iv-1)+tsai_z(p,iv-1))+(tlai_z(p,iv)+tsai_z(p,iv)))
                end if

                ! Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction

                t1 = min(h*laisum, 40._r8)
                s1 = exp(-t1)
                t1 = min(twostext(p)*laisum, 40._r8)
                s2 = exp(-t1)
                fsun_z(p,iv) = s2

                ! ===============
                ! Direct beam
                ! ===============

                ! Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai

                u1 = b - c1/albgrd(c,ib)
                u2 = b - c1*albgrd(c,ib)
                u3 = f + c1*albgrd(c,ib)
                tmp2 = u1 - avmu(p)*h
                tmp3 = u1 + avmu(p)*h
                d1 = p1*tmp2/s1 - p2*tmp3*s1
                tmp4 = u2 + avmu(p)*h
                tmp5 = u2 - avmu(p)*h
                d2 = tmp4/s1 - tmp5*s1
                h1 = -d*p4 - c1*f
                tmp6 = d - h1*p3/sigma
                tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
                h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
                h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
                h4 = -f*p3 - c1*d
                tmp8 = h4/sigma
                tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
                h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
                h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

                a1 = h1 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
                   + h2         * (1._r8 - s2*s1) / (twostext(p) + h) &
                   + h3         * (1._r8 - s2/s1) / (twostext(p) - h)

                a2 = h4 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
                   + h5         * (1._r8 - s2*s1) / (twostext(p) + h) &
                   + h6         * (1._r8 - s2/s1) / (twostext(p) - h)

                ! Derivatives for h2, h3, h5, h6 and a1, a2

                v = d1
                dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

                u = tmp6 * tmp2 / s1 - p2 * tmp7
                du = h * tmp6 * tmp2 / s1 + twostext(p) * p2 * tmp7
                dh2 = (v * du - u * dv) / (v * v)

                u = -tmp6 * tmp3 * s1 + p1 * tmp7
                du = h * tmp6 * tmp3 * s1 - twostext(p) * p1 * tmp7
                dh3 = (v * du - u * dv) / (v * v)

                v = d2
                dv = h * tmp4 / s1 + h * tmp5 * s1

                u = -h4/sigma * tmp4 / s1 - tmp9
                du = -h * h4/sigma * tmp4 / s1 + twostext(p) * tmp9
                dh5 = (v * du - u * dv) / (v * v)

                u = h4/sigma * tmp5 * s1 + tmp9
                du = -h * h4/sigma * tmp5 * s1 - twostext(p) * tmp9
                dh6 = (v * du - u * dv) / (v * v)

                da1 = h1/sigma * s2*s2 + h2 * s2*s1 + h3 * s2/s1 &
                    + (1._r8 - s2*s1) / (twostext(p) + h) * dh2 &
                    + (1._r8 - s2/s1) / (twostext(p) - h) * dh3
                da2 = h4/sigma * s2*s2 + h5 * s2*s1 + h6 * s2/s1 &
                    + (1._r8 - s2*s1) / (twostext(p) + h) * dh5 &
                    + (1._r8 - s2/s1) / (twostext(p) - h) * dh6

                ! Flux derivatives

                d_ftid = -twostext(p)*h4/sigma*s2 - h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1
                d_fabd = -(dh2+dh3) + (1._r8-albgrd(c,ib))*twostext(p)*s2 - (1._r8-albgri(c,ib))*d_ftid
                d_fabd_sun = (1._r8 - omega(p,ib)) * (twostext(p)*s2 + 1._r8 / avmu(p) * (da1 + da2))
                d_fabd_sha = d_fabd - d_fabd_sun

                fabd_sun_z(p,iv) = max(d_fabd_sun, 0._r8)
                fabd_sha_z(p,iv) = max(d_fabd_sha, 0._r8)

                ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need 
                ! to normalize derivatives by sunlit or shaded fraction to get
                ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha

                fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / fsun_z(p,iv)
                fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))

                ! ===============
                ! Diffuse
                ! ===============

                ! Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai

                u1 = b - c1/albgri(c,ib)
                u2 = b - c1*albgri(c,ib)
                tmp2 = u1 - avmu(p)*h
                tmp3 = u1 + avmu(p)*h
                d1 = p1*tmp2/s1 - p2*tmp3*s1
                tmp4 = u2 + avmu(p)*h
                tmp5 = u2 - avmu(p)*h
                d2 = tmp4/s1 - tmp5*s1
                h7 = (c1*tmp2) / (d1*s1)
                h8 = (-c1*tmp3*s1) / d1
                h9 = tmp4 / (d2*s1)
                h10 = (-tmp5*s1) / d2

                a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
                a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)

                ! Derivatives for h7, h8, h9, h10 and a1, a2

                v = d1
                dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

                u = c1 * tmp2 / s1
                du = h * c1 * tmp2 / s1
                dh7 = (v * du - u * dv) / (v * v)

                u = -c1 * tmp3 * s1
                du = h * c1 * tmp3 * s1
                dh8 = (v * du - u * dv) / (v * v)

                v = d2
                dv = h * tmp4 / s1 + h * tmp5 * s1

                u = tmp4 / s1
                du = h * tmp4 / s1
                dh9 = (v * du - u * dv) / (v * v)

                u = -tmp5 * s1
                du = h * tmp5 * s1
                dh10 = (v * du - u * dv) / (v * v)

                da1 = h7*s2*s1 +  h8*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh7 + (1._r8-s2/s1)/(twostext(p)-h)*dh8
                da2 = h9*s2*s1 + h10*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh9 + (1._r8-s2/s1)/(twostext(p)-h)*dh10

                ! Flux derivatives

                d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1
                d_fabi = -(dh7+dh8) - (1._r8-albgri(c,ib))*d_ftii
                d_fabi_sun = (1._r8 - omega(p,ib)) / avmu(p) * (da1 + da2)
                d_fabi_sha = d_fabi - d_fabi_sun

                fabi_sun_z(p,iv) = max(d_fabi_sun, 0._r8)
                fabi_sha_z(p,iv) = max(d_fabi_sha, 0._r8)

                ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need 
                ! to normalize derivatives by sunlit or shaded fraction to get
                ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha

                fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / fsun_z(p,iv)
                fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))

                end do   ! end of canopy layer loop
             end if
          end if

       end do   ! end of pft loop
    end do   ! end of radiation band loop

  end subroutine TwoStream

end module SurfaceAlbedoMod
