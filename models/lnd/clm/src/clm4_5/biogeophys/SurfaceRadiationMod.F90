module SurfaceRadiationMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceRadiationMod
!
! !DESCRIPTION:
! Calculate solar fluxes absorbed by vegetation and ground surface
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varctl  , only: iulog

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 11/26/03, Peter Thornton: Added new routine for improved treatment of
!    sunlit/shaded canopy radiation.
! 4/26/05, Peter Thornton: Adopted the sun/shade algorithm as the default,
!    removed the old SurfaceRadiation(), and renamed SurfaceRadiationSunShade()
!    as SurfaceRadiation().
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRadiation
!
! !INTERFACE:
   subroutine SurfaceRadiation(lbp, ubp, num_nourbanp, filter_nourbanp)
!
! !DESCRIPTION: 
! Solar fluxes absorbed by vegetation and ground surface
! Note possible problem when land is on different grid than atmosphere.
! Land may have sun above the horizon (coszen > 0) but atmosphere may
! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.
! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
! land may have sun below horizon. This is okay because fabd, fabi,
! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
! the radiation is reflected. NDVI should equal zero in this case.
! However, the way the code is currently implemented this is only true
! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
!
! !USES:
     use clmtype
     use clm_atmlnd      , only : clm_a2l
     use clm_varpar      , only : numrad
     use clm_varcon      , only : spval, istsoil, degpsec, isecspday
     use clm_varcon      , only : istcrop
     use clm_varctl      , only : subgridflag
     use clm_time_manager, only : get_curr_date, get_step_size
     use clm_varpar      , only : nlevsno
     use SNICARMod       , only : DO_SNO_OC
     use abortutils      , only : endrun
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: lbp, ubp                   ! pft upper and lower bounds
     integer, intent(in) :: num_nourbanp               ! number of pfts in non-urban points in pft filter
     integer, intent(in) :: filter_nourbanp(ubp-lbp+1) ! pft filter for non-urban points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/18/02, Peter Thornton: Migrated to new data structures. Added a pft loop.
! 6/05/03, Peter Thornton: Modified sunlit/shaded canopy treatment. Original code
! had all radiation being absorbed in the sunlit canopy, and now the sunlit and shaded
! canopies are each given the appropriate fluxes.  There was also an inconsistency in
! the original code, where parsun was not being scaled by leaf area, and so represented
! the entire canopy flux.  This goes into Stomata (in CanopyFluxes) where it is assumed
! to be a flux per unit leaf area. In addition, the fpsn flux coming out of Stomata was
! being scaled back up to the canopy by multiplying by lai, but the input radiation flux was
! for the entire canopy to begin with.  Corrected this inconsistency in this version, so that
! the parsun and parsha fluxes going into canopy fluxes are per unit lai in the sunlit and
! shaded canopies.
! 6/9/03, Peter Thornton: Moved coszen from g%gps to c%cps to avoid problem
! with OpenMP threading over columns, where different columns hit the radiation
! time step at different times during execution.
! 6/10/03, Peter Thornton: Added constraint on negative tot_aid, instead of
! exiting with error. Appears to be happening only at roundoff level.
! 6/11/03, Peter Thornton: Moved calculation of ext inside if (coszen),
! and added check on laisun = 0 and laisha = 0 in calculation of sun_aperlai
! and sha_aperlai.
! 11/26/03, Peter Thornton: During migration to new vector code, created 
!   this as a new routine to handle sunlit/shaded canopy calculations.
! 03/28/08, Mark Flanner: Incorporated SNICAR, including absorbed solar radiation
!   in each snow layer and top soil layer, and optional radiative forcing calculation
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
     real(r8), pointer :: albsod(:,:)      ! direct-beam soil albedo (col,bnd) [frc]
     real(r8), pointer :: albsoi(:,:)      ! diffuse soil albedo (col,bnd) [frc]
     real(r8), pointer :: sabg_soil(:)     ! solar radiation absorbed by soil (W/m**2)
     real(r8), pointer :: sabg_snow(:)     ! solar radiation absorbed by snow (W/m**2)
     logical , pointer :: pactive(:)       ! true=>do computations on this pft (see reweightMod for details)
     integer , pointer :: ivt(:)           ! pft vegetation type
     integer , pointer :: pcolumn(:)       ! pft's column index
     integer , pointer :: pgridcell(:)     ! pft's gridcell index
     real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
     real(r8), pointer :: esai(:)          ! one-sided stem area index with burying by snow
     real(r8), pointer :: londeg(:)        ! longitude (degrees)
     real(r8), pointer :: latdeg(:)        ! latitude (degrees)
     real(r8), pointer :: coszen(:)	   ! cosine of solar zenith angle
     real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
     real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
     real(r8), pointer :: fabd(:,:)        ! flux absorbed by canopy per unit direct flux
     real(r8), pointer :: fabd_sun(:,:)    ! flux absorbed by sunlit canopy per unit direct flux
     real(r8), pointer :: fabd_sha(:,:)    ! flux absorbed by shaded canopy per unit direct flux
     real(r8), pointer :: fabi(:,:)        ! flux absorbed by canopy per unit diffuse flux
     real(r8), pointer :: fabi_sun(:,:)    ! flux absorbed by sunlit canopy per unit diffuse flux
     real(r8), pointer :: fabi_sha(:,:)    ! flux absorbed by shaded canopy per unit diffuse flux
     real(r8), pointer :: ftdd(:,:)        ! down direct flux below canopy per unit direct flux
     real(r8), pointer :: ftid(:,:)        ! down diffuse flux below canopy per unit direct flux
     real(r8), pointer :: ftii(:,:)        ! down diffuse flux below canopy per unit diffuse flux
     integer , pointer :: nrad(:)          ! number of canopy layers, above snow for radiative transfer
     real(r8), pointer :: fabd_sun_z(:,:)  ! absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabd_sha_z(:,:)  ! absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sun_z(:,:)  ! absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fabi_sha_z(:,:)  ! absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
     real(r8), pointer :: fsun_z(:,:)      ! sunlit fraction of canopy layer
     real(r8), pointer :: tlai_z(:,:)      ! tlai increment for canopy layer
     real(r8), pointer :: tsai_z(:,:)      ! tsai increment for canopy layer
     real(r8), pointer :: albgrd(:,:)      ! ground albedo (direct)
     real(r8), pointer :: albgri(:,:)      ! ground albedo (diffuse)
     real(r8), pointer :: albd(:,:)        ! surface albedo (direct)
     real(r8), pointer :: albi(:,:)        ! surface albedo (diffuse)
!
! local pointers to original implicit out arguments
!
     real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
     real(r8), pointer :: laisun(:)        ! sunlit leaf area
     real(r8), pointer :: laisha(:)        ! shaded leaf area
     real(r8), pointer :: laisun_z(:,:)    ! sunlit leaf area for canopy layer
     real(r8), pointer :: laisha_z(:,:)    ! shaded leaf area for canopy layer
     real(r8), pointer :: parsun_z(:,:)    ! absorbed PAR for sunlit leaves in canopy layer
     real(r8), pointer :: parsha_z(:,:)    ! absorbed PAR for shaded leaves in canopy layer
     real(r8), pointer :: sabg(:)          ! solar radiation absorbed by ground (W/m**2)
     real(r8), pointer :: sabv(:)          ! solar radiation absorbed by vegetation (W/m**2)
     real(r8), pointer :: fsa(:)           ! solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: fsa_r(:)         ! rural solar radiation absorbed (total) (W/m**2)
     integer , pointer :: ityplun(:)       ! landunit type
     integer , pointer :: plandunit(:)     ! index into landunit level quantities
     real(r8), pointer :: fsr(:)           ! solar radiation reflected (W/m**2)
     real(r8), pointer :: fsds_vis_d(:)    ! incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_d(:)    ! incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_i(:)    ! incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i(:)    ! incident diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_d(:)     ! reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_d(:)     ! reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_i(:)     ! reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_i(:)     ! reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_d_ln(:) ! incident direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsds_nir_d_ln(:) ! incident direct beam nir solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_vis_d_ln(:)  ! reflected direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_ln(:)  ! reflected direct beam nir solar rad at local noon (W/m**2)
     real(r8), pointer :: fsds_vis_i_ln(:) ! incident diffuse beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: parveg_ln(:)     ! absorbed par by vegetation at local noon (W/m**2)
     real(r8), pointer :: flx_absdv(:,:)   ! direct flux absorption factor (col,lyr): VIS [frc]
     real(r8), pointer :: flx_absdn(:,:)   ! direct flux absorption factor (col,lyr): NIR [frc]
     real(r8), pointer :: flx_absiv(:,:)   ! diffuse flux absorption factor (col,lyr): VIS [frc]
     real(r8), pointer :: flx_absin(:,:)   ! diffuse flux absorption factor (col,lyr): NIR [frc]
     integer , pointer :: snl(:)           ! negative number of snow layers [nbr]
     real(r8), pointer :: albgrd_pur(:,:)    ! pure snow ground albedo (direct)
     real(r8), pointer :: albgri_pur(:,:)    ! pure snow ground albedo (diffuse)
     real(r8), pointer :: albgrd_bc(:,:)     ! ground albedo without BC (direct) (col,bnd)
     real(r8), pointer :: albgri_bc(:,:)     ! ground albedo without BC (diffuse) (col,bnd)
     real(r8), pointer :: albgrd_oc(:,:)     ! ground albedo without OC (direct) (col,bnd)
     real(r8), pointer :: albgri_oc(:,:)     ! ground albedo without OC (diffuse) (col,bnd)
     real(r8), pointer :: albgrd_dst(:,:)    ! ground albedo without dust (direct) (col,bnd)
     real(r8), pointer :: albgri_dst(:,:)    ! ground albedo without dust (diffuse) (col,bnd)
     real(r8), pointer :: albsnd_hst(:,:)    ! snow albedo, direct, for history files (col,bnd) [frc]
     real(r8), pointer :: albsni_hst(:,:)    ! snow ground albedo, diffuse, for history files (col,bnd
     real(r8), pointer :: sabg_lyr(:,:)      ! absorbed radiative flux (pft,lyr) [W/m2]
     real(r8), pointer :: sabg_pen(:)        ! (rural) shortwave radiation penetrating top soisno layer [W/m2]
     real(r8), pointer :: sfc_frc_aer(:)     ! surface forcing of snow with all aerosols (pft) [W/m2]
     real(r8), pointer :: sfc_frc_bc(:)      ! surface forcing of snow with BC (pft) [W/m2]
     real(r8), pointer :: sfc_frc_oc(:)      ! surface forcing of snow with OC (pft) [W/m2]
     real(r8), pointer :: sfc_frc_dst(:)     ! surface forcing of snow with dust (pft) [W/m2]
     real(r8), pointer :: sfc_frc_aer_sno(:) ! surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer :: sfc_frc_bc_sno(:)  ! surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer :: sfc_frc_oc_sno(:)  ! surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer :: sfc_frc_dst_sno(:) ! surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer :: frac_sno(:)      ! fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: fsr_sno_vd(:)    ! reflected visible, direct radiation from snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsr_sno_nd(:)    ! reflected near-IR, direct radiation from snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsr_sno_vi(:)    ! reflected visible, diffuse radiation from snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsr_sno_ni(:)    ! reflected near-IR, diffuse radiation from snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsds_sno_vd(:)   ! incident visible, direct radiation on snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsds_sno_nd(:)   ! incident near-IR, direct radiation on snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsds_sno_vi(:)   ! incident visible, diffuse radiation on snow (for history files) (pft) [W/m2]
     real(r8), pointer :: fsds_sno_ni(:)   ! incident near-IR, diffuse radiation on snow (for history files) (pft) [W/m2]
     real(r8), pointer :: snow_depth(:)        ! snow height (m)

!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
     integer , parameter :: nband = numrad    ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06_r8   ! prevents overflow for division by zero
     integer  :: fp                  ! non-urban filter pft index
     integer  :: p                   ! pft index
     integer  :: c                   ! column index
     integer  :: l                   ! landunit index
     integer  :: g                   ! grid cell index
     integer  :: ib                  ! waveband number (1=vis, 2=nir)
     integer  :: iv                  ! canopy layer
     real(r8) :: absrad              ! absorbed solar radiation (W/m**2)
     real(r8) :: rnir                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: trd(lbp:ubp,numrad) ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri(lbp:ubp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(lbp:ubp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(lbp:ubp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     integer  :: local_secp1         ! seconds into current date in local time
     real(r8) :: dtime               ! land model time step (sec)
     integer  :: year,month,day,secs !  calendar info for current time step
     integer  :: i                   ! layer index [idx]
     real(r8) :: sabg_snl_sum        ! temporary, absorbed energy in all active snow layers [W/m2]
     real(r8) :: absrad_pur          ! temp: absorbed solar radiation by pure snow [W/m2]
     real(r8) :: absrad_bc           ! temp: absorbed solar radiation without BC [W/m2]
     real(r8) :: absrad_oc           ! temp: absorbed solar radiation without OC [W/m2]
     real(r8) :: absrad_dst          ! temp: absorbed solar radiation without dust [W/m2]
     real(r8) :: sabg_pur(lbp:ubp)   ! solar radiation absorbed by ground with pure snow [W/m2]
     real(r8) :: sabg_bc(lbp:ubp)    ! solar radiation absorbed by ground without BC [W/m2]
     real(r8) :: sabg_oc(lbp:ubp)    ! solar radiation absorbed by ground without OC [W/m2]
     real(r8) :: sabg_dst(lbp:ubp)   ! solar radiation absorbed by ground without dust [W/m2]
     real(r8) :: parveg(lbp:ubp)     ! absorbed par by vegetation (W/m**2)
     integer  :: nstep               ! time step number
!------------------------------------------------------------------------------

     ! Assign local pointers to multi-level derived type members (gridcell level)

     londeg        => clm3%g%londeg
     latdeg        => clm3%g%latdeg
     forc_solad    => clm_a2l%forc_solad
     forc_solai    => clm_a2l%forc_solai

     ! Assign local pointers to multi-level derived type members (landunit level)

     ityplun       => clm3%g%l%itype

     ! Assign local pointers to multi-level derived type members (column level)

     albgrd        => clm3%g%l%c%cps%albgrd
     albgri        => clm3%g%l%c%cps%albgri
     coszen        => clm3%g%l%c%cps%coszen

     ! Assign local pointers to derived type members (pft-level)

     albsod        => clm3%g%l%c%cps%albsod
     albsoi        => clm3%g%l%c%cps%albsoi
     sabg_soil     => clm3%g%l%c%p%pef%sabg_soil
     sabg_snow     => clm3%g%l%c%p%pef%sabg_snow
     pactive       => clm3%g%l%c%p%active
     plandunit     => clm3%g%l%c%p%landunit
     ivt           => clm3%g%l%c%p%itype
     pcolumn       => clm3%g%l%c%p%column
     pgridcell     => clm3%g%l%c%p%gridcell
     elai          => clm3%g%l%c%p%pps%elai
     esai          => clm3%g%l%c%p%pps%esai
     laisun        => clm3%g%l%c%p%pps%laisun
     laisha        => clm3%g%l%c%p%pps%laisha
     laisun_z      => clm3%g%l%c%p%pps%laisun_z
     laisha_z      => clm3%g%l%c%p%pps%laisha_z
     albd          => clm3%g%l%c%p%pps%albd
     albi          => clm3%g%l%c%p%pps%albi
     fabd          => clm3%g%l%c%p%pps%fabd
     fabd_sun      => clm3%g%l%c%p%pps%fabd_sun
     fabd_sha      => clm3%g%l%c%p%pps%fabd_sha
     fabi          => clm3%g%l%c%p%pps%fabi
     fabi_sun      => clm3%g%l%c%p%pps%fabi_sun
     fabi_sha      => clm3%g%l%c%p%pps%fabi_sha
     ftdd          => clm3%g%l%c%p%pps%ftdd
     ftid          => clm3%g%l%c%p%pps%ftid
     ftii          => clm3%g%l%c%p%pps%ftii
     nrad          => clm3%g%l%c%p%pps%nrad
     fabd_sun_z    => clm3%g%l%c%p%pps%fabd_sun_z
     fabd_sha_z    => clm3%g%l%c%p%pps%fabd_sha_z
     fabi_sun_z    => clm3%g%l%c%p%pps%fabi_sun_z
     fabi_sha_z    => clm3%g%l%c%p%pps%fabi_sha_z
     fsun_z        => clm3%g%l%c%p%pps%fsun_z
     tlai_z        => clm3%g%l%c%p%pps%tlai_z
     tsai_z        => clm3%g%l%c%p%pps%tsai_z
     fsun          => clm3%g%l%c%p%pps%fsun
     sabg          => clm3%g%l%c%p%pef%sabg
     sabv          => clm3%g%l%c%p%pef%sabv
     snow_depth        => clm3%g%l%c%cps%snow_depth
     fsa           => clm3%g%l%c%p%pef%fsa
     fsa_r         => clm3%g%l%c%p%pef%fsa_r
     fsr           => clm3%g%l%c%p%pef%fsr
     parsun_z      => clm3%g%l%c%p%pef%parsun_z
     parsha_z      => clm3%g%l%c%p%pef%parsha_z
     fsds_vis_d    => clm3%g%l%c%p%pef%fsds_vis_d
     fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
     fsds_vis_i    => clm3%g%l%c%p%pef%fsds_vis_i
     fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
     fsr_vis_d     => clm3%g%l%c%p%pef%fsr_vis_d
     fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
     fsr_vis_i     => clm3%g%l%c%p%pef%fsr_vis_i
     fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
     fsds_vis_d_ln => clm3%g%l%c%p%pef%fsds_vis_d_ln
     fsds_nir_d_ln => clm3%g%l%c%p%pef%fsds_nir_d_ln
     parveg_ln     => clm3%g%l%c%p%pef%parveg_ln
     fsds_vis_i_ln => clm3%g%l%c%p%pef%fsds_vis_i_ln
     fsr_vis_d_ln  => clm3%g%l%c%p%pef%fsr_vis_d_ln
     fsr_nir_d_ln  => clm3%g%l%c%p%pef%fsr_nir_d_ln
     
     ! Assign local pointers to derived type members (ecophysiological)

     frac_sno         => clm3%g%l%c%cps%frac_sno
     flx_absdv        => clm3%g%l%c%cps%flx_absdv
     flx_absdn        => clm3%g%l%c%cps%flx_absdn
     flx_absiv        => clm3%g%l%c%cps%flx_absiv
     flx_absin        => clm3%g%l%c%cps%flx_absin
     sabg_lyr         => clm3%g%l%c%p%pef%sabg_lyr
     sabg_pen         => clm3%g%l%c%p%pef%sabg_pen
     snl              => clm3%g%l%c%cps%snl
     sfc_frc_aer      => clm3%g%l%c%p%pef%sfc_frc_aer
     sfc_frc_aer_sno  => clm3%g%l%c%p%pef%sfc_frc_aer_sno
     albgrd_pur       => clm3%g%l%c%cps%albgrd_pur
     albgri_pur       => clm3%g%l%c%cps%albgri_pur
     sfc_frc_bc       => clm3%g%l%c%p%pef%sfc_frc_bc
     sfc_frc_bc_sno   => clm3%g%l%c%p%pef%sfc_frc_bc_sno
     albgrd_bc        => clm3%g%l%c%cps%albgrd_bc
     albgri_bc        => clm3%g%l%c%cps%albgri_bc
     sfc_frc_oc       => clm3%g%l%c%p%pef%sfc_frc_oc
     sfc_frc_oc_sno   => clm3%g%l%c%p%pef%sfc_frc_oc_sno
     albgrd_oc        => clm3%g%l%c%cps%albgrd_oc
     albgri_oc        => clm3%g%l%c%cps%albgri_oc
     sfc_frc_dst      => clm3%g%l%c%p%pef%sfc_frc_dst
     sfc_frc_dst_sno  => clm3%g%l%c%p%pef%sfc_frc_dst_sno
     albgrd_dst       => clm3%g%l%c%cps%albgrd_dst
     albgri_dst       => clm3%g%l%c%cps%albgri_dst
     albsnd_hst       => clm3%g%l%c%cps%albsnd_hst
     albsni_hst       => clm3%g%l%c%cps%albsni_hst
     fsr_sno_vd       => clm3%g%l%c%p%pef%fsr_sno_vd
     fsr_sno_nd       => clm3%g%l%c%p%pef%fsr_sno_nd
     fsr_sno_vi       => clm3%g%l%c%p%pef%fsr_sno_vi
     fsr_sno_ni       => clm3%g%l%c%p%pef%fsr_sno_ni
     fsds_sno_vd      => clm3%g%l%c%p%pef%fsds_sno_vd
     fsds_sno_nd      => clm3%g%l%c%p%pef%fsds_sno_nd
     fsds_sno_vi      => clm3%g%l%c%p%pef%fsds_sno_vi
     fsds_sno_ni      => clm3%g%l%c%p%pef%fsds_sno_ni

     ! Determine seconds off current time step
     
     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Initialize fluxes

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        ! was redundant b/c filter already included wt>0; 
        ! not redundant anymore with chg in filter definition
        if (pactive(p)) then
           sabg_soil(p)  = 0._r8
           sabg_snow(p)  = 0._r8
           sabg(p)       = 0._r8
           sabv(p)       = 0._r8
           fsa(p)        = 0._r8
           l = plandunit(p)
           if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
             fsa_r(p)      = 0._r8
           end if
           sabg_lyr(p,:) = 0._r8
           sabg_pur(p)   = 0._r8
           sabg_bc(p)    = 0._r8
           sabg_oc(p)    = 0._r8
           sabg_dst(p)   = 0._r8
           do iv = 1, nrad(p)
              parsun_z(p,iv) = 0._r8
              parsha_z(p,iv) = 0._r8
              laisun_z(p,iv) = 0._r8
              laisha_z(p,iv) = 0._r8
           end do
        end if
     end do 

     ! Loop over pfts to calculate laisun_z and laisha_z for each layer.
     ! Derive canopy laisun, laisha, and fsun from layer sums.
     ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
     ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        if (pactive(p)) then

           laisun(p) = 0._r8
           laisha(p) = 0._r8
           do iv = 1, nrad(p)
              laisun_z(p,iv) = tlai_z(p,iv) * fsun_z(p,iv)
              laisha_z(p,iv) = tlai_z(p,iv) * (1._r8 - fsun_z(p,iv))
              laisun(p) = laisun(p) + laisun_z(p,iv) 
              laisha(p) = laisha(p) + laisha_z(p,iv) 
           end do
           if (elai(p) > 0._r8) then
              fsun(p) = laisun(p) / elai(p)
           else
              fsun(p) = 0._r8
           end if
        
        end if
     end do
        
     ! Loop over nband wavebands
     do ib = 1, nband
        do fp = 1,num_nourbanp
           p = filter_nourbanp(fp)
           if (pactive(p)) then
              c = pcolumn(p)
              l = plandunit(p)
              g = pgridcell(p)

              ! Absorbed by canopy
              
              cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
              cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
              sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
              fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
              if (ib == 1) then
                parveg(p) = cad(p,ib) + cai(p,ib)
              end if
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
              end if

              ! Absorbed PAR profile through canopy
              ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
              ! are canopy integrated so that layer values equal big leaf values.
              
              if (ib == 1) then
                 do iv = 1, nrad(p)
                    parsun_z(p,iv) = forc_solad(g,ib)*fabd_sun_z(p,iv) + forc_solai(g,ib)*fabi_sun_z(p,iv)
                    parsha_z(p,iv) = forc_solad(g,ib)*fabd_sha_z(p,iv) + forc_solai(g,ib)*fabi_sha_z(p,iv)
                 end do
              end if

              ! Transmitted = solar fluxes incident on ground

              trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
              tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)

              ! Solar radiation absorbed by ground surface
              
              ! calculate absorbed solar by soil/snow separately
              absrad  = trd(p,ib)*(1._r8-albsod(c,ib)) + tri(p,ib)*(1._r8-albsoi(c,ib))
              sabg_soil(p) = sabg_soil(p) + absrad
              absrad  = trd(p,ib)*(1._r8-albsnd_hst(c,ib)) + tri(p,ib)*(1._r8-albsni_hst(c,ib))
              sabg_snow(p) = sabg_snow(p) + absrad
              absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
              sabg(p) = sabg(p) + absrad
              fsa(p)  = fsa(p)  + absrad
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + absrad
              end if
              if (snl(c) == 0) then
                 sabg_snow(p) = sabg(p)
                 sabg_soil(p) = sabg(p)
              endif
              ! if no subgrid fluxes, make sure to set both components equal to weighted average
              if (subgridflag == 0) then 
                 sabg_snow(p) = sabg(p)
                 sabg_soil(p) = sabg(p)
              endif

#if (defined SNICAR_FRC)
              ! Solar radiation absorbed by ground surface without BC
              absrad_bc = trd(p,ib)*(1._r8-albgrd_bc(c,ib)) + tri(p,ib)*(1._r8-albgri_bc(c,ib))
              sabg_bc(p) = sabg_bc(p) + absrad_bc

              ! Solar radiation absorbed by ground surface without OC
              absrad_oc = trd(p,ib)*(1._r8-albgrd_oc(c,ib)) + tri(p,ib)*(1._r8-albgri_oc(c,ib))
              sabg_oc(p) = sabg_oc(p) + absrad_oc

              ! Solar radiation absorbed by ground surface without dust
              absrad_dst = trd(p,ib)*(1._r8-albgrd_dst(c,ib)) + tri(p,ib)*(1._r8-albgri_dst(c,ib))
              sabg_dst(p) = sabg_dst(p) + absrad_dst

              ! Solar radiation absorbed by ground surface without any aerosols
              absrad_pur = trd(p,ib)*(1._r8-albgrd_pur(c,ib)) + tri(p,ib)*(1._r8-albgri_pur(c,ib))
              sabg_pur(p) = sabg_pur(p) + absrad_pur
#endif

           end if
        end do ! end of pft loop
     end do ! end nbands loop   

     !   compute absorbed flux in each snow layer and top soil layer,
     !   based on flux factors computed in the radiative transfer portion of SNICAR.

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
           if (pactive(p)) then
           c = pcolumn(p)
           l = plandunit(p)
           sabg_snl_sum = 0._r8

           ! CASE1: No snow layers: all energy is absorbed in top soil layer
           if (snl(c) == 0) then
              sabg_lyr(p,:) = 0._r8
              sabg_lyr(p,1) = sabg(p)
              sabg_snl_sum  = sabg_lyr(p,1)
   
           ! CASE 2: Snow layers present: absorbed radiation is scaled according to 
           ! flux factors computed by SNICAR
           else
              do i = -nlevsno+1,1,1
                 sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + flx_absdn(c,i)*trd(p,2) + &
                                 flx_absiv(c,i)*tri(p,1) + flx_absin(c,i)*tri(p,2)
                 ! summed radiation in active snow layers:
                 if (i >= snl(c)+1) then
                    sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
                 endif
              enddo
   
              ! Error handling: The situation below can occur when solar radiation is 
              ! NOT computed every timestep.
              ! When the number of snow layers has changed in between computations of the 
              ! absorbed solar energy in each layer, we must redistribute the absorbed energy
              ! to avoid physically unrealistic conditions. The assumptions made below are 
              ! somewhat arbitrary, but this situation does not arise very frequently. 
              ! This error handling is implemented to accomodate any value of the
              ! radiation frequency.
              ! change condition to match sabg_snow isntead of sabg
              if (abs(sabg_snl_sum-sabg_snow(p)) > 0.00001_r8) then
                 if (snl(c) == 0) then
                    sabg_lyr(p,-4:0) = 0._r8
                    sabg_lyr(p,1) = sabg(p)
                 elseif (snl(c) == -1) then
                    sabg_lyr(p,-4:-1) = 0._r8
                    sabg_lyr(p,0) = sabg_snow(p)*0.6_r8
                    sabg_lyr(p,1) = sabg_snow(p)*0.4_r8
                 else
                    sabg_lyr(p,:) = 0._r8
                    sabg_lyr(p,snl(c)+1) = sabg_snow(p)*0.75_r8
                    sabg_lyr(p,snl(c)+2) = sabg_snow(p)*0.25_r8
                 endif
              endif

              ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
              ! to prevent unrealistic timestep soil warming 
              if (subgridflag == 0) then 
                 if (snow_depth(c) < 0.10_r8) then
                    if (snl(c) == 0) then
                       sabg_lyr(p,-4:0) = 0._r8
                       sabg_lyr(p,1) = sabg(p)
                    elseif (snl(c) == -1) then
                       sabg_lyr(p,-4:-1) = 0._r8
                       sabg_lyr(p,0) = sabg(p)
                       sabg_lyr(p,1) = 0._r8
                    else
                       sabg_lyr(p,:) = 0._r8
                       sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                       sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                    endif
                 endif
              endif
           endif

           ! This situation should not happen:
           if (abs(sum(sabg_lyr(p,:))-sabg_snow(p)) > 0.00001_r8) then
              write(iulog,*) "SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation. pft = ",   &
                             p," Col= ", c, " Diff= ",sum(sabg_lyr(p,:))-sabg_snow(p), " sabg_snow(p)= ", sabg_snow(p), " sabg_sum(p)= ", &
                             sum(sabg_lyr(p,:)), " snl(c)= ", snl(c)
              write(iulog,*) "flx_absdv1= ", trd(p,1)*(1.-albgrd(c,1)), "flx_absdv2= ", sum(flx_absdv(c,:))*trd(p,1)
              write(iulog,*) "flx_absiv1= ", tri(p,1)*(1.-albgri(c,1))," flx_absiv2= ", sum(flx_absiv(c,:))*tri(p,1)
              write(iulog,*) "flx_absdn1= ", trd(p,2)*(1.-albgrd(c,2))," flx_absdn2= ", sum(flx_absdn(c,:))*trd(p,2)
              write(iulog,*) "flx_absin1= ", tri(p,2)*(1.-albgri(c,2))," flx_absin2= ", sum(flx_absin(c,:))*tri(p,2)
   
              write(iulog,*) "albgrd_nir= ", albgrd(c,2)
              write(iulog,*) "coszen= ", coszen(c)
              call endrun()
           endif

           ! Diagnostic: shortwave penetrating ground (e.g. top layer)
           if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
              sabg_pen(p) = sabg(p) - sabg_lyr(p, snl(c)+1)
           end if

#if (defined SNICAR_FRC)

           ! BC aerosol forcing (pft-level):
           sfc_frc_bc(p) = sabg(p) - sabg_bc(p)
   
           ! OC aerosol forcing (pft-level):
           if (DO_SNO_OC) then
              sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
           else
              sfc_frc_oc(p) = 0._r8
           endif
   
           ! dust aerosol forcing (pft-level):
           sfc_frc_dst(p) = sabg(p) - sabg_dst(p)
   
           ! all-aerosol forcing (pft-level):
           sfc_frc_aer(p) = sabg(p) - sabg_pur(p)        
           
           ! forcings averaged only over snow:
           if (frac_sno(c) > 0._r8) then
              sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
              sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
              sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
              sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
           else
              sfc_frc_bc_sno(p)  = spval
              sfc_frc_oc_sno(p)  = spval
              sfc_frc_dst_sno(p) = spval
              sfc_frc_aer_sno(p) = spval
           endif

#endif
        endif
     enddo

     ! Radiation diagnostics

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        if (pactive(p)) then
           g = pgridcell(p)
        
           ! NDVI and reflected solar radiation
           
           rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
           rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
           fsr(p) = rvis + rnir
           
           fsds_vis_d(p) = forc_solad(g,1)
           fsds_nir_d(p) = forc_solad(g,2)
           fsds_vis_i(p) = forc_solai(g,1)
           fsds_nir_i(p) = forc_solai(g,2)
           fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
           fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
           fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
           fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)
           
           local_secp1 = secs + nint((londeg(g)/degpsec)/dtime)*dtime
           local_secp1 = mod(local_secp1,isecspday)
           if (local_secp1 == isecspday/2) then
              fsds_vis_d_ln(p) = forc_solad(g,1)
              fsds_nir_d_ln(p) = forc_solad(g,2)
              fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
              fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
              fsds_vis_i_ln(p) = forc_solai(g,1)
              parveg_ln(p)     = parveg(p)
           else
              fsds_vis_d_ln(p) = spval
              fsds_nir_d_ln(p) = spval
              fsr_vis_d_ln(p) = spval
              fsr_nir_d_ln(p) = spval
              fsds_vis_i_ln(p) = spval
              parveg_ln(p)     = spval
           end if

           ! diagnostic variables (downwelling and absorbed radiation partitioning) for history files
           ! (OPTIONAL)
           c = pcolumn(p)
           if (snl(c) < 0) then
              fsds_sno_vd(p) = forc_solad(g,1)
              fsds_sno_nd(p) = forc_solad(g,2)
              fsds_sno_vi(p) = forc_solai(g,1)
              fsds_sno_ni(p) = forc_solai(g,2)

              fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
              fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
              fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
              fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
           else
              fsds_sno_vd(p) = spval
              fsds_sno_nd(p) = spval
              fsds_sno_vi(p) = spval
              fsds_sno_ni(p) = spval

              fsr_sno_vd(p) = spval
              fsr_sno_nd(p) = spval
              fsr_sno_vi(p) = spval
              fsr_sno_ni(p) = spval
           endif
        end if
     end do 

   end subroutine SurfaceRadiation

end module SurfaceRadiationMod
