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
     use clm_varcon      , only : istice_mec
     use clm_varcon      , only : istcrop
     use clm_time_manager, only : get_curr_date, get_step_size
     use clm_varpar      , only : nlevsno
     use SNICARMod       , only : DO_SNO_OC
     use abortutils      , only : endrun
     use clm_varctl      , only : use_snicar_frc
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
! 6/9/03, Peter Thornton: Moved coszen from gps to c%cps to avoid problem
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
     integer , pointer :: ivt(:)           ! pft vegetation type
     integer , pointer :: pcolumn(:)       ! pft's column index
     integer , pointer :: pgridcell(:)     ! pft's gridcell index
     real(r8), pointer :: pwtgcell(:)      ! pft's weight relative to corresponding gridcell
     real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
     real(r8), pointer :: esai(:)          ! one-sided stem area index with burying by snow
     real(r8), pointer :: londeg(:)        ! longitude (degrees)
     real(r8), pointer :: latdeg(:)        ! latitude (degrees)
     real(r8), pointer :: slasun(:)        ! specific leaf area for sunlit canopy, projected area basis (m^2/gC)
     real(r8), pointer :: slasha(:)        ! specific leaf area for shaded canopy, projected area basis (m^2/gC)
     real(r8), pointer :: gdir(:)	   ! leaf projection in solar direction (0 to 1)
     real(r8), pointer :: omega(:,:)       ! fraction of intercepted radiation that is scattered (0 to 1)
     real(r8), pointer :: coszen(:)	   ! cosine of solar zenith angle
     real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
     real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
     real(r8), pointer :: fabd(:,:)        ! flux absorbed by veg per unit direct flux
     real(r8), pointer :: fabi(:,:)        ! flux absorbed by veg per unit diffuse flux
     real(r8), pointer :: ftdd(:,:)        ! down direct flux below veg per unit dir flx
     real(r8), pointer :: ftid(:,:)        ! down diffuse flux below veg per unit dir flx
     real(r8), pointer :: ftii(:,:)        ! down diffuse flux below veg per unit dif flx
     real(r8), pointer :: albgrd(:,:)      ! ground albedo (direct)
     real(r8), pointer :: albgri(:,:)      ! ground albedo (diffuse)
     real(r8), pointer :: albd(:,:)        ! surface albedo (direct)
     real(r8), pointer :: albi(:,:)        ! surface albedo (diffuse)
     real(r8), pointer :: slatop(:)        ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real(r8), pointer :: dsladlai(:)      ! dSLA/dLAI, projected area basis [m^2/gC]
!
! local pointers to original implicit out arguments
!
     real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
     real(r8), pointer :: laisun(:)        ! sunlit leaf area
     real(r8), pointer :: laisha(:)        ! shaded leaf area
     real(r8), pointer :: sabg(:)          ! solar radiation absorbed by ground (W/m**2)
     real(r8), pointer :: sabv(:)          ! solar radiation absorbed by vegetation (W/m**2)
     real(r8), pointer :: fsa(:)           ! solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: fsa_r(:)         ! rural solar radiation absorbed (total) (W/m**2)
     integer , pointer :: ityplun(:)       ! landunit type
     integer , pointer :: plandunit(:)     ! index into landunit level quantities
     real(r8), pointer :: parsun(:)        ! average absorbed PAR for sunlit leaves (W/m**2)
     real(r8), pointer :: parsha(:)        ! average absorbed PAR for shaded leaves (W/m**2)
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
     real(r8), pointer :: eff_kid(:,:)     ! effective extinction coefficient for indirect from direct
     real(r8), pointer :: eff_kii(:,:)     ! effective extinction coefficient for indirect from indirect
     real(r8), pointer :: sun_faid(:,:)    ! fraction sun canopy absorbed indirect from direct
     real(r8), pointer :: sun_faii(:,:)    ! fraction sun canopy absorbed indirect from indirect
     real(r8), pointer :: sha_faid(:,:)    ! fraction shade canopy absorbed indirect from direct
     real(r8), pointer :: sha_faii(:,:)    ! fraction shade canopy absorbed indirect from indirect
     real(r8), pointer :: sun_add(:,:)     ! sun canopy absorbed direct from direct (W/m**2)
     real(r8), pointer :: tot_aid(:,:)     ! total canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sun_aid(:,:)     ! sun canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sun_aii(:,:)     ! sun canopy absorbed indirect from indirect (W/m**2)
     real(r8), pointer :: sha_aid(:,:)     ! shade canopy absorbed indirect from direct (W/m**2)
     real(r8), pointer :: sha_aii(:,:)     ! shade canopy absorbed indirect from indirect (W/m**2)
     real(r8), pointer :: sun_atot(:,:)    ! sun canopy total absorbed (W/m**2)
     real(r8), pointer :: sha_atot(:,:)    ! shade canopy total absorbed (W/m**2)
     real(r8), pointer :: sun_alf(:,:)     ! sun canopy total absorbed by leaves (W/m**2)
     real(r8), pointer :: sha_alf(:,:)     ! shade canopy total absored by leaves (W/m**2)
     real(r8), pointer :: sun_aperlai(:,:) ! sun canopy total absorbed per unit LAI (W/m**2)
     real(r8), pointer :: sha_aperlai(:,:) ! shade canopy total absorbed per unit LAI (W/m**2)
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
     real(r8), pointer :: snowdp(:)        ! snow height (m)

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
     real(r8) :: absrad              ! absorbed solar radiation (W/m**2)
     real(r8) :: rnir                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: laifra              ! leaf area fraction of canopy
     real(r8) :: trd(lbp:ubp,numrad) ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri(lbp:ubp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(lbp:ubp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(lbp:ubp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     real(r8) :: vai(lbp:ubp)        ! total leaf area index + stem area index, one sided
     real(r8) :: ext                 ! optical depth direct beam per unit LAI+SAI
     real(r8) :: t1, t2              ! temporary variables
     real(r8) :: cosz
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
!------------------------------------------------------------------------------

     ! Assign local pointers to multi-level derived type members (gridcell level)

     londeg        => grc%londeg
     latdeg        => grc%latdeg
     forc_solad    => clm_a2l%forc_solad
     forc_solai    => clm_a2l%forc_solai

     ! Assign local pointers to multi-level derived type members (landunit level)

     ityplun       => lun%itype

     ! Assign local pointers to multi-level derived type members (column level)

     albgrd        => cps%albgrd
     albgri        => cps%albgri
     coszen        => cps%coszen

     ! Assign local pointers to derived type members (pft-level)

     plandunit     => pft%landunit
     ivt           => pft%itype
     pcolumn       => pft%column
     pgridcell     => pft%gridcell
     pwtgcell      => pft%wtgcell
     elai          => pps%elai
     esai          => pps%esai
     slasun        => pps%slasun
     slasha        => pps%slasha
     gdir          => pps%gdir
     omega         => pps%omega
     laisun        => pps%laisun
     laisha        => pps%laisha
     fabd          => pps%fabd
     fabi          => pps%fabi
     ftdd          => pps%ftdd
     ftid          => pps%ftid
     ftii          => pps%ftii
     albd          => pps%albd
     albi          => pps%albi
     fsun          => pps%fsun
     sabg          => pef%sabg
     sabv          => pef%sabv
     snowdp        => cps%snowdp
     fsa           => pef%fsa
     fsa_r         => pef%fsa_r
     fsr           => pef%fsr
     parsun        => pef%parsun
     parsha        => pef%parsha
     fsds_vis_d    => pef%fsds_vis_d
     fsds_nir_d    => pef%fsds_nir_d
     fsds_vis_i    => pef%fsds_vis_i
     fsds_nir_i    => pef%fsds_nir_i
     fsr_vis_d     => pef%fsr_vis_d
     fsr_nir_d     => pef%fsr_nir_d
     fsr_vis_i     => pef%fsr_vis_i
     fsr_nir_i     => pef%fsr_nir_i
     fsds_vis_d_ln => pef%fsds_vis_d_ln
     fsds_nir_d_ln => pef%fsds_nir_d_ln
     fsr_vis_d_ln  => pef%fsr_vis_d_ln
     fsr_nir_d_ln  => pef%fsr_nir_d_ln
     eff_kid       => pps%eff_kid
     eff_kii       => pps%eff_kii
     sun_faid      => pps%sun_faid
     sun_faii      => pps%sun_faii
     sha_faid      => pps%sha_faid
     sha_faii      => pps%sha_faii
     sun_add       => pef%sun_add
     tot_aid       => pef%tot_aid
     sun_aid       => pef%sun_aid
     sun_aii       => pef%sun_aii
     sha_aid       => pef%sha_aid
     sha_aii       => pef%sha_aii
     sun_atot      => pef%sun_atot
     sha_atot      => pef%sha_atot
     sun_alf       => pef%sun_alf
     sha_alf       => pef%sha_alf
     sun_aperlai   => pef%sun_aperlai
     sha_aperlai   => pef%sha_aperlai
     
     ! Assign local pointers to derived type members (ecophysiological)

     slatop           => pftcon%slatop
     dsladlai         => pftcon%dsladlai
     frac_sno         => cps%frac_sno
     flx_absdv        => cps%flx_absdv
     flx_absdn        => cps%flx_absdn
     flx_absiv        => cps%flx_absiv
     flx_absin        => cps%flx_absin
     sabg_lyr         => pef%sabg_lyr
     snl              => cps%snl
     sfc_frc_aer      => pef%sfc_frc_aer
     sfc_frc_aer_sno  => pef%sfc_frc_aer_sno
     albgrd_pur       => cps%albgrd_pur
     albgri_pur       => cps%albgri_pur
     sfc_frc_bc       => pef%sfc_frc_bc
     sfc_frc_bc_sno   => pef%sfc_frc_bc_sno
     albgrd_bc        => cps%albgrd_bc
     albgri_bc        => cps%albgri_bc
     sfc_frc_oc       => pef%sfc_frc_oc
     sfc_frc_oc_sno   => pef%sfc_frc_oc_sno
     albgrd_oc        => cps%albgrd_oc
     albgri_oc        => cps%albgri_oc
     sfc_frc_dst      => pef%sfc_frc_dst
     sfc_frc_dst_sno  => pef%sfc_frc_dst_sno
     albgrd_dst       => cps%albgrd_dst
     albgri_dst       => cps%albgri_dst
     albsnd_hst       => cps%albsnd_hst
     albsni_hst       => cps%albsni_hst
     fsr_sno_vd       => pef%fsr_sno_vd
     fsr_sno_nd       => pef%fsr_sno_nd
     fsr_sno_vi       => pef%fsr_sno_vi
     fsr_sno_ni       => pef%fsr_sno_ni
     fsds_sno_vd      => pef%fsds_sno_vd
     fsds_sno_nd      => pef%fsds_sno_nd
     fsds_sno_vi      => pef%fsds_sno_vi
     fsds_sno_ni      => pef%fsds_sno_ni

     ! Determine seconds off current time step
     
     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Determine fluxes

     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        ! was redundant b/c filter already included wt>0; 
        ! not redundant anymore with chg in filter definition
        l = plandunit(p)
        !Note: Some glacier_mec pfts may have zero weight
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
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
        end if
     end do 

     ! Loop over pfts to calculate fsun, etc
     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        l = plandunit(p)
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           c = pcolumn(p)
           g = pgridcell(p)
        
           vai(p) = elai(p) + esai(p)
           if (coszen(c) > 0._r8 .and. elai(p) > 0._r8 .and. gdir(p) > 0._r8) then
              cosz = max(0.001_r8, coszen(c))
              ext = gdir(p)/cosz
              t1 = min(ext*elai(p), 40.0_r8)
              t2 = exp(-t1)
              fsun(p) = (1._r8-t2)/t1
              
              ! new control on low lai, to avoid numerical problems in
              ! calculation of slasun, slasha
              ! PET: 2/29/04
              
              if (elai(p) > 0.01_r8) then
                 laisun(p) = elai(p)*fsun(p)
                 laisha(p) = elai(p)*(1._r8-fsun(p))
                 
                 ! calculate the average specific leaf area for sunlit and shaded
                 ! canopies, when effective LAI > 0
                 slasun(p) = (t2*dsladlai(ivt(p))*ext*elai(p) + &
                              t2*dsladlai(ivt(p)) + &
                              t2*slatop(ivt(p))*ext - &
                              dsladlai(ivt(p)) - &
                              slatop(ivt(p))*ext) / &
                              (ext*(t2-1._r8))
                 slasha(p) = ((slatop(ivt(p)) + &
                             (dsladlai(ivt(p)) * elai(p)/2.0_r8)) * elai(p) - &
                             laisun(p)*slasun(p)) / laisha(p)
              else
                 ! special case for low elai
                 fsun(p) = 1._r8
                 laisun(p) = elai(p)
                 laisha(p) = 0._r8
                 slasun(p) = slatop(ivt(p))
                 slasha(p) = 0._r8
              end if
           else
              fsun(p)   = 0._r8
              laisun(p) = 0._r8
              laisha(p) = elai(p)
              slasun(p) = 0._r8
              slasha(p) = 0._r8
           end if
        end if
     end do
        
     ! Loop over nband wavebands
     do ib = 1, nband
        do fp = 1,num_nourbanp
           p = filter_nourbanp(fp)
           l = plandunit(p)
           if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
              c = pcolumn(p)
              g = pgridcell(p)
              
              ! Absorbed by canopy
              
              cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
              cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
              sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
              fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
              l = plandunit(p)
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
              end if
              
              ! Transmitted = solar fluxes incident on ground
              
              trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
              tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)
     
              ! Solar radiation absorbed by ground surface
              
              absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
              sabg(p) = sabg(p) + absrad
              fsa(p)  = fsa(p)  + absrad
              if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + absrad
              end if

              if (use_snicar_frc) then
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
              end if


              ! New sunlit.shaded canopy algorithm
              
              if (coszen(c) > 0._r8 .and. elai(p) > 0._r8 .and. gdir(p) > 0._r8 ) then
                 
                 ! 1. calculate flux of direct beam radiation absorbed in the 
                 ! sunlit canopy as direct (sun_add), and the flux of direct
                 ! beam radiation absorbed in the total canopy as indirect
                 
                 sun_add(p,ib) = forc_solad(g,ib) * (1._r8-ftdd(p,ib)) * (1._r8-omega(p,ib))
                 tot_aid(p,ib) = (forc_solad(g,ib) * fabd(p,ib)) - sun_add(p,ib)
                 
                 ! the following constraint set to catch round-off level errors
                 ! that can cause negative tot_aid
                 
                 tot_aid(p,ib) = max(tot_aid(p,ib), 0._r8)
                 
                 ! 2. calculate the effective extinction coefficients for indirect
                 ! transmission originating from direct and indirect streams,
                 ! using ftid and ftii
                 
                 !eff_kid(p,ib) = -(log(ftid(p,ib)))/vai(p)
                 !eff_kii(p,ib) = -(log(ftii(p,ib)))/vai(p)
                 
                 ! 3. calculate the fraction of indirect radiation being absorbed 
                 ! in the sunlit and shaded canopy fraction. Some of this indirect originates in
                 ! the direct beam and some originates in the indirect beam.

                 !sun_faid(p,ib) = 1.-exp(-eff_kid(p,ib) * vaisun(p))
                 !sun_faii(p,ib) = 1.-exp(-eff_kii(p,ib) * vaisun(p))
                 sun_faid(p,ib) = fsun(p)
                 sun_faii(p,ib) = fsun(p)
                 sha_faid(p,ib) = 1._r8-sun_faid(p,ib)
                 sha_faii(p,ib) = 1._r8-sun_faii(p,ib)

                 ! 4. calculate the total indirect flux absorbed by the sunlit
                 ! and shaded canopy based on these fractions and the fabd and
                 ! fabi from surface albedo calculations

                 sun_aid(p,ib) = tot_aid(p,ib) * sun_faid(p,ib)
                 sun_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sun_faii(p,ib)
                 sha_aid(p,ib) = tot_aid(p,ib) * sha_faid(p,ib)
                 sha_aii(p,ib) = forc_solai(g,ib)*fabi(p,ib)*sha_faii(p,ib)
                 
                 ! 5. calculate the total flux absorbed in the sunlit and shaded
                 ! canopy as the sum of these terms
                 
                 sun_atot(p,ib) = sun_add(p,ib) + sun_aid(p,ib) + sun_aii(p,ib)
                 sha_atot(p,ib) = sha_aid(p,ib) + sha_aii(p,ib)
                 
                 ! 6. calculate the total flux absorbed by leaves in the sunlit
                 ! and shaded canopies
                 
                 laifra = elai(p)/vai(p)
                 sun_alf(p,ib) = sun_atot(p,ib) * laifra
                 sha_alf(p,ib) = sha_atot(p,ib) * laifra
                 
                 ! 7. calculate the fluxes per unit lai in the sunlit and shaded
                 ! canopies
                 
                 if (laisun(p) > 0._r8) then
                    sun_aperlai(p,ib) = sun_alf(p,ib)/laisun(p)
                 else
                    sun_aperlai(p,ib) = 0._r8
                 endif
                 if (laisha(p) > 0._r8) then
                    sha_aperlai(p,ib) = sha_alf(p,ib)/laisha(p)
                 else
                    sha_aperlai(p,ib) = 0._r8
                 endif
             
              else   ! coszen = 0 or elai = 0
                 
                 sun_add(p,ib)     = 0._r8
                 tot_aid(p,ib)     = 0._r8
                 eff_kid(p,ib)     = 0._r8
                 eff_kii(p,ib)     = 0._r8
                 sun_faid(p,ib)    = 0._r8
                 sun_faii(p,ib)    = 0._r8
                 sha_faid(p,ib)    = 0._r8
                 sha_faii(p,ib)    = 0._r8
                 sun_aid(p,ib)     = 0._r8
                 sun_aii(p,ib)     = 0._r8
                 sha_aid(p,ib)     = 0._r8
                 sha_aii(p,ib)     = 0._r8
                 sun_atot(p,ib)    = 0._r8
                 sha_atot(p,ib)    = 0._r8
                 sun_alf(p,ib)     = 0._r8
                 sha_alf(p,ib)     = 0._r8
                 sun_aperlai(p,ib) = 0._r8
                 sha_aperlai(p,ib) = 0._r8
                 
              end if
           end if
        end do ! end of pft loop
     end do ! end nbands loop   

     
     !   compute absorbed flux in each snow layer and top soil layer,
     !   based on flux factors computed in the radiative transfer portion of SNICAR.
     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
           l = plandunit(p)
           if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           c = pcolumn(p)
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
              if (abs(sabg_snl_sum-sabg(p)) > 0.00001_r8) then
                 if (snl(c) == 0) then
                    sabg_lyr(p,-4:0) = 0._r8
                    sabg_lyr(p,1) = sabg(p)
                 elseif (snl(c) == -1) then
                    sabg_lyr(p,-4:-1) = 0._r8
                    sabg_lyr(p,0) = sabg(p)*0.6_r8
                    sabg_lyr(p,1) = sabg(p)*0.4_r8
                 else
                    sabg_lyr(p,:) = 0._r8
                    sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                    sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                 endif
              endif

              ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
              ! to prevent unrealistic timestep soil warming 
              if (snowdp(c) < 0.10_r8) then
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

           ! This situation should not happen:
           if (abs(sum(sabg_lyr(p,:))-sabg(p)) > 0.00001_r8) then
              write(iulog,*) "SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation. pft = ",   &
                             p," Col= ", c, " Diff= ",sum(sabg_lyr(p,:))-sabg(p), " sabg(p)= ", sabg(p), " sabg_sum(p)= ", &
                             sum(sabg_lyr(p,:)), " snl(c)= ", snl(c)
              write(iulog,*) "flx_absdv1= ", trd(p,1)*(1.-albgrd(c,1)), "flx_absdv2= ", sum(flx_absdv(c,:))*trd(p,1)
              write(iulog,*) "flx_absiv1= ", tri(p,1)*(1.-albgri(c,1))," flx_absiv2= ", sum(flx_absiv(c,:))*tri(p,1)
              write(iulog,*) "flx_absdn1= ", trd(p,2)*(1.-albgrd(c,2))," flx_absdn2= ", sum(flx_absdn(c,:))*trd(p,2)
              write(iulog,*) "flx_absin1= ", tri(p,2)*(1.-albgri(c,2))," flx_absin2= ", sum(flx_absin(c,:))*tri(p,2)
   
              write(iulog,*) "albgrd_nir= ", albgrd(c,2)
              write(iulog,*) "coszen= ", coszen(c)
              call endrun()
           endif

 
           if (use_snicar_frc) then

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

           end if
        endif
     enddo


     do fp = 1,num_nourbanp
        p = filter_nourbanp(fp)
        l = plandunit(p)
        if (pwtgcell(p)>0._r8 .or. ityplun(l)==istice_mec) then
           g = pgridcell(p)
        
           ! Final step of new sunlit/shaded canopy algorithm
           ! 8. calculate the total and per-unit-lai fluxes for PAR in the
           ! sunlit and shaded canopy leaf fractions
           
           parsun(p) = sun_aperlai(p,1)
           parsha(p) = sha_aperlai(p,1)
           
           ! The following code is duplicated from SurfaceRadiation
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
           else
              fsds_vis_d_ln(p) = spval
              fsds_nir_d_ln(p) = spval
              fsr_vis_d_ln(p) = spval
              fsr_nir_d_ln(p) = spval
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
