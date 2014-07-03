module VOCEmissionMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: VOCEmissionMod
!
! !DESCRIPTION:
! Volatile organic compound emission
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl,   only : iulog
  use abortutils,   only : endrun
  use clm_varpar,   only : numpft
  use pftvarcon ,   only : ndllf_evr_tmp_tree,  ndllf_evr_brl_tree,    &
                           ndllf_dcd_brl_tree,  nbrdlf_evr_trp_tree,   &
                           nbrdlf_evr_tmp_tree, nbrdlf_dcd_brl_shrub,  &
                           nbrdlf_dcd_trp_tree, nbrdlf_dcd_tmp_tree,   &
                           nbrdlf_dcd_brl_tree, nbrdlf_evr_shrub,      &
                           nc3_arctic_grass,    nc3crop,               &
                           nc4_grass,           noveg

  use shr_megan_mod,  only : shr_megan_megcomps_n, shr_megan_megcomp_t, shr_megan_linkedlist
  use shr_megan_mod,  only : shr_megan_mechcomps_n, shr_megan_mechcomps, shr_megan_mapped_emisfctrs
  use MEGANFactorsMod,only : Agro, Amat, Anew, Aold, betaT, ct1, ct2, LDF, Ceo

!
! !PUBLIC TYPES:
  implicit none
  save

  logical, parameter :: debug = .false.

!
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
  public :: VOCEmission_init
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
! !IROUTINE: VOCEmission
!
! !INTERFACE:
  subroutine VOCEmission (lbp, ubp, num_soilp, filter_soilp )
!
! ! NEW DESCRIPTION
! Volatile organic compound emission
! This code simulates volatile organic compound emissions following
! MEGAN (Model of Emissions of Gases and Aerosols from Nature) v2.1 
! for 20 compound classes. The original description of this
! algorithm (for isoprene only) can be found in Guenther et al., 2006
! (we follow equations 2-9, 16-17, 20 for explicit canopy).
! The model scheme came be described as:
!    E= epsilon * gamma * rho
! VOC flux (E) [ug m-2 h-1] is calculated from baseline emission
! factors (epsilon) [ug m-2 h-1] which are specified for each of the 16
! CLM PFTs (in input file) OR in the case of isoprene, from
! mapped EFs for each PFT which reflect species divergence of emissions,
! particularly in North America. 
! The emission activity factor (gamma) [unitless] for includes 
! dependence on PPFT, temperature, LAI, leaf age and soil moisture.
! For isoprene only we also include the effect of CO2 inhibition as
! described by Heald et al., 2009. 
! The canopy environment constant was calculated offline for CLM+CAM at 
! standard conditions.
! We assume that the escape efficiency (rho) here is unity following
! Guenther et al., 2006.
! A manuscript describing MEGAN 2.1 and the implementation in CLM is
! in preparation: Guenther, Heald et al., 2012
! Subroutine written to operate at the patch level.
!
! Input: <filename> to be read in with EFs and some parameters.  
!        Currently these are set in procedure init_EF_params
! Output: vocflx(shr_megan_mechcomps_n) !VOC flux [moles/m2/sec]
!
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_atmlnd   , only : clm_a2l
    use clmtype
    use domainMod,     only : ldomain 
    use clm_varcon   , only : spval
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_soilp                   ! number of columns in soil pft filter
    integer, intent(in) :: filter_soilp(num_soilp)     ! pft filter for soil
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
! 2/1/02: Peter Thornton: migration to new data structure
! 4/15/06: Colette L. Heald: modify for updated MEGAN model (Guenther et al., 2006)
! 4/29/11: Colette L. Heald: expand MEGAN to 20 compound classes
! 7 Feb 2012: Francis Vitt: Implemented capability to specify MEGAN emissions in namelist
!                           and read in MEGAN factors from file.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pgridcell(:)     ! gridcell index of corresponding pft
    integer , pointer :: pcolumn(:)       ! column index of corresponding pft
    integer , pointer :: ivt(:)           ! pft vegetation type for current
    real(r8), pointer :: t_veg(:)         ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
    real(r8), pointer :: clayfrac(:)      ! fraction of soil that is clay
    real(r8), pointer :: sandfrac(:)      ! fraction of soil that is sand
    real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation     (visible only)
    real(r8), pointer :: elai_p(:)        ! one-sided leaf area index from previous timestep
    real(r8), pointer :: t_veg24(:)       ! avg pft vegetation temperature for last 24 hrs
    real(r8), pointer :: t_veg240(:)      ! avg pft vegetation temperature for last 240 hrs
    real(r8), pointer :: fsun24(:)        ! sunlit fraction of canopy last 24 hrs
    real(r8), pointer :: fsun240(:)       ! sunlit fraction of canopy last 240 hrs
    real(r8), pointer :: forc_solad24(:)  ! direct beam radiation last 24hrs  (visible only)
    real(r8), pointer :: forc_solai24(:)  ! diffuse radiation  last 24hrs     (visible only)
    real(r8), pointer :: forc_solad240(:) ! direct beam radiation last 240hrs (visible only)
    real(r8), pointer :: forc_solai240(:) ! diffuse radiation  last 240hrs    (visible only)
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water (m3/m3)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice soil content (kg/m3)
    real(r8), pointer :: dz(:,:)          ! depth of layer (m)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: cisun(:)         ! sunlit intracellular CO2 (Pa)
    real(r8), pointer :: cisha(:)         ! shaded intracellular CO2 (Pa)
    real(r8), pointer :: forc_pbot(:)     ! atmospheric pressure (Pa)
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: vocflx(:,:)      ! VOC flux [moles/m2/sec]
    real(r8), pointer :: vocflx_tot(:)    ! VOC flux [moles/m2/sec]

    type(megan_out_type), pointer :: meg_out(:) ! fluxes for CLM history 

    real(r8), pointer :: gamma_out(:)
    real(r8), pointer :: gammaT_out(:)
    real(r8), pointer :: gammaP_out(:)
    real(r8), pointer :: gammaL_out(:)
    real(r8), pointer :: gammaA_out(:)
    real(r8), pointer :: gammaS_out(:)
    real(r8), pointer :: gammaC_out(:)

    real(r8), pointer :: Eopt_out(:)     
    real(r8), pointer :: topt_out(:)
    real(r8), pointer :: alpha_out(:)
    real(r8), pointer :: cp_out(:)
    real(r8), pointer :: paru_out(:)
    real(r8), pointer :: par24u_out(:)
    real(r8), pointer :: par240u_out(:)
    real(r8), pointer :: para_out(:)
    real(r8), pointer :: par24a_out(:)
    real(r8), pointer :: par240a_out(:)

!
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,g,c                ! indices
    real(r8) :: epsilon                 ! emission factor [ug m-2 h-1]
    real(r8) :: par_sun                 ! temporary
    real(r8) :: par24_sun               ! temporary
    real(r8) :: par240_sun              ! temporary
    real(r8) :: par_sha                 ! temporary
    real(r8) :: par24_sha               ! temporary
    real(r8) :: par240_sha              ! temporary
    real(r8) :: gamma                   ! activity factor (accounting for light, T, age, LAI conditions)
    real(r8) :: gamma_p                 ! activity factor for PPFD
    real(r8) :: gamma_l                 ! activity factor for PPFD & LAI
    real(r8) :: gamma_t                 ! activity factor for temperature
    real(r8) :: gamma_a                 ! activity factor for leaf age
    real(r8) :: gamma_sm                ! activity factor for soil moisture
    real(r8) :: gamma_c                 ! activity factor for CO2 (only isoprene)
    
    integer :: class_num, n_meg_comps, imech, imeg, ii
    character(len=16) :: mech_name

    real(r8) :: vocflx_meg(shr_megan_megcomps_n)
    type(shr_megan_megcomp_t), pointer :: meg_cmp

    real(r8) :: cp, alpha,  Eopt, topt  ! for history output

    ! factor used convert MEGAN units [micro-grams/m2/hr] to CAM srf emis units [g/m2/sec]
    real(r8), parameter :: megemis_units_factor = 1._r8/3600._r8/1.e6_r8

!    real(r8) :: root_depth(0:numpft)    ! Root depth [m]
!
!!-----------------------------------------------------------------------
!
!    ! root depth (m) (defined based on Zeng et al., 2001, cf Guenther 2006)
!    root_depth(noveg)                                     = 0._r8   ! bare-soil
!    root_depth(ndllf_evr_tmp_tree:ndllf_evr_brl_tree)     = 1.8_r8  ! evergreen tree
!    root_depth(ndllf_dcd_brl_tree)                        = 2.0_r8  ! needleleaf deciduous boreal tree
!    root_depth(nbrdlf_evr_trp_tree:nbrdlf_evr_tmp_tree)   = 3.0_r8  ! broadleaf evergreen tree
!    root_depth(nbrdlf_dcd_trp_tree:nbrdlf_dcd_brl_tree)   = 2.0_r8  ! broadleaf deciduous tree
!    root_depth(nbrdlf_evr_shrub:nbrdlf_dcd_brl_shrub)     = 2.5_r8  ! shrub
!    root_depth(nc3_arctic_grass:numpft)                   = 1.5_r8  ! grass/crop
!
!-----------------------------------------------------------------------
    if ( shr_megan_mechcomps_n < 1) return

    ! Assign local pointers to derived type members (gridcell-level)
    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai
    forc_pbot  => clm_a2l%forc_pbot

    ! Assign local pointers to derived subtypes components (column-level)
    h2osoi_vol       => cws%h2osoi_vol
    h2osoi_ice       => cws%h2osoi_ice
    dz               => cps%dz
    bsw              => cps%bsw
    watsat           => cps%watsat
    sucsat           => cps%sucsat

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell        => pft%gridcell
    pcolumn          => pft%column
    ivt              => pft%itype
    t_veg            => pes%t_veg
    fsun             => pps%fsun
    elai             => pps%elai
    clayfrac         => pps%clayfrac
    sandfrac         => pps%sandfrac

    cisun            => pps%cisun
    cisha            => pps%cisha

    vocflx           => pvf%vocflx
    vocflx_tot       => pvf%vocflx_tot
    meg_out          => pvf%meg

    gammaL_out       => pvf%gammaL_out
    gammaT_out       => pvf%gammaT_out
    gammaP_out       => pvf%gammaP_out
    gammaA_out       => pvf%gammaA_out
    gammaS_out       => pvf%gammaS_out
    gammaC_out       => pvf%gammaC_out
    gamma_out        => pvf%gamma_out

    Eopt_out         => pvf%Eopt_out
    topt_out         => pvf%topt_out
    alpha_out        => pvf%alpha_out
    cp_out           => pvf%cp_out
    paru_out         => pvf%paru_out
    par24u_out       => pvf%par24u_out
    par240u_out      => pvf%par240u_out
    para_out         => pvf%para_out
    par24a_out       => pvf%par24a_out
    par240a_out      => pvf%par240a_out

    t_veg24          => pvs%t_veg24
    t_veg240         => pvs%t_veg240
    forc_solad24     => pvs%fsd24
    forc_solad240    => pvs%fsd240
    forc_solai24     => pvs%fsi24
    forc_solai240    => pvs%fsi240
    fsun24           => pvs%fsun24
    fsun240          => pvs%fsun240
    elai_p           => pvs%elai_p

    ! initialize variables which get passed to the atmosphere
    vocflx(lbp:ubp,:) = 0._r8
    vocflx_tot(lbp:ubp) = 0._r8

    do imeg=1,shr_megan_megcomps_n
      meg_out(imeg)%flux_out(lbp:ubp) = 0._r8
    enddo
    
    gamma_out(lbp:ubp) = spval
    gammaP_out(lbp:ubp) = spval
    gammaT_out(lbp:ubp) = spval
    gammaA_out(lbp:ubp) = spval
    gammaS_out(lbp:ubp) = spval
    gammaL_out(lbp:ubp) = spval
    gammaC_out(lbp:ubp) = spval

    paru_out(lbp:ubp) = spval
    par24u_out(lbp:ubp) = spval
    par240u_out(lbp:ubp) = spval

    para_out(lbp:ubp) = spval
    par24a_out(lbp:ubp) = spval
    par240a_out(lbp:ubp) = spval

    alpha_out(lbp:ubp) = spval
    cp_out(lbp:ubp) = spval

    topt_out(lbp:ubp) = spval
    Eopt_out(lbp:ubp) = spval

    ! initalize to zero since this might not alway get set
    vocflx_meg(:) = 0._r8

    ! Begin loop over points
    !_______________________________________________________________________________
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       g = pgridcell(p)
       c = pcolumn(p)

       ! initialize EF
       epsilon=0._r8

       ! calculate VOC emissions for non-bare ground PFTs
       if (ivt(p) > 0) then 
          gamma=0._r8

          ! Calculate PAR: multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)
          !------------------------
          ! SUN:
          par_sun = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
          par24_sun = (forc_solad24(p) + fsun24(p) * forc_solai24(p)) * 4.6_r8
          par240_sun = (forc_solad240(p) + fsun240(p) * forc_solai240(p)) * 4.6_r8
          ! SHADE:
          par_sha = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
          par24_sha = ((1._r8 - fsun24(p)) * forc_solai24(p)) * 4.6_r8
          par240_sha = ((1._r8 - fsun240(p)) * forc_solai240(p)) * 4.6_r8

          ! Activity factor for LAI (Guenther et al., 2006): all species
          gamma_l = get_gamma_L(fsun240(p), elai(p))

          ! Activity factor for soil moisture: all species (commented out for now)
!          gamma_sm = get_gamma_SM(clayfrac(p), sandfrac(p), h2osoi_vol(c,:), h2osoi_ice(c,:), &
!               dz(c,:), bsw(c,:), watsat(c,:), sucsat(c,:), root_depth(ivt(p)))
          gamma_sm = 1.0_r8

          ! Loop through VOCs for light, temperature and leaf age activity factor & apply
          ! all final activity factors to baseline emission factors
          !_______________________________________________________________________________

          ! loop over megan compounds
          meg_cmp => shr_megan_linkedlist
          meg_cmp_loop: do while(associated(meg_cmp))
             imeg = meg_cmp%index

             ! set emis factor
             ! if specified, set EF for isoprene with mapped values
             if ( trim(meg_cmp%name) == 'isoprene' .and. shr_megan_mapped_emisfctrs) then
                epsilon = get_map_EF(ivt(p),g)
             else
                epsilon = meg_cmp%emis_factors(ivt(p))
             end if

             class_num = meg_cmp%class_number

             ! Activity factor for PPFD
             gamma_p = get_gamma_P(par_sun, par24_sun, par240_sun, par_sha, par24_sha, par240_sha, &
                  fsun(p), fsun240(p), forc_solad240(p),forc_solai240(p), LDF(class_num), cp, alpha)

             ! Activity factor for T
             gamma_t = get_gamma_T(t_veg240(p), t_veg24(p),t_veg(p), ct1(class_num), ct2(class_num),&
                                   betaT(class_num),LDF(class_num), Ceo(class_num), Eopt, topt)

             ! Activity factor for Leaf Age
             gamma_a = get_gamma_A(ivt(p), elai_p(p),elai(p),class_num)

             ! Activity factor for CO2 (only for isoprene)
             if (trim(meg_cmp%name) == 'isoprene') then 
                gamma_c = get_gamma_C(cisun(p),cisha(p),forc_pbot(g),fsun(p))
             else
                gamma_c = 1._r8
             end if

             ! Calculate total scaling factor
             gamma = gamma_l * gamma_sm * gamma_a * gamma_p * gamma_T * gamma_c

             if ( (gamma >=0.0_r8) .and. (gamma< 100._r8) ) then

                vocflx_meg(imeg) = epsilon * gamma * megemis_units_factor / meg_cmp%molec_weight ! moles/m2/sec

                ! assign to arrays for history file output (not weighted by landfrac)
                meg_out(imeg)%flux_out(p) = meg_out(imeg)%flux_out(p) &
                                          + epsilon * gamma * megemis_units_factor*1.e-3_r8 ! Kg/m2/sec

                if (imeg==1) then 
                   ! 
                   gamma_out(p)=gamma
                   gammaP_out(p)=gamma_p
                   gammaT_out(p)=gamma_t
                   gammaA_out(p)=gamma_a
                   gammaS_out(p)=gamma_sm
                   gammaL_out(p)=gamma_l
                   gammaC_out(p)=gamma_c

                   paru_out(p)=par_sun
                   par24u_out(p)=par24_sun
                   par240u_out(p)=par240_sun

                   para_out(p)=par_sha
                   par24a_out(p)=par24_sha
                   par240a_out(p)=par240_sha

                   alpha_out(p)=alpha
                   cp_out(p)=cp

                   topt_out(p)=topt
                   Eopt_out(p)=Eopt

                end if
             endif

             if (debug .and. gamma > 0.0_r8) then
                write(iulog,*) 'MEGAN: n, megan name, epsilon, gamma, vocflx: ', &
                     imeg, meg_cmp%name, epsilon, gamma, vocflx_meg(imeg), gamma_p,gamma_t,gamma_a,gamma_sm,gamma_l
             endif

             meg_cmp => meg_cmp%next_megcomp
          enddo meg_cmp_loop

          ! sum up the megan compound fluxes for the fluxes of chem mechanism compounds 
          do imech = 1,shr_megan_mechcomps_n
             n_meg_comps = shr_megan_mechcomps(imech)%n_megan_comps
             do imeg = 1,n_meg_comps ! loop over number of megan compounds that make up the nth mechanism compoud
                ii = shr_megan_mechcomps(imech)%megan_comps(imeg)%ptr%index
                vocflx(p,imech) = vocflx(p,imech) + vocflx_meg(ii)
             enddo
             vocflx_tot(p) = vocflx_tot(p) + vocflx(p,imech) ! moles/m2/sec
          enddo

       end if ! ivt(1:15 only)

    enddo ! fp 

  end subroutine VOCEmission
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! !IROUTINE: init_EF_params
!
! !INTERFACE:
  subroutine VOCEmission_init(  )
    ! Interface to set all input parameters for 20 VOC compound classes.
    ! including EFs for 16(+1 bare ground) PFTs.
    ! For now set all specified values, in future to be replaced with values read in from file.
    ! (heald, 04/27/11)

    use shr_megan_mod,   only : shr_megan_factors_file
    use MEGANFactorsMod, only : megan_factors_init, megan_factors_get
    
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!
! !USES
!
! !ARGUMENTS:
    implicit none

!    character(len=*),intent(in) :: filename
!
! !LOCAL VARIABLES:
!-----------------------------------------------------------------------

    integer :: nmech, nmeg
    type(shr_megan_megcomp_t), pointer :: meg_cmp

    integer  :: class_num
    real(r8) :: factors(numpft)
    real(r8) :: molec_wght

    if ( shr_megan_mechcomps_n < 1) return

    call megan_factors_init( shr_megan_factors_file )
    
    meg_cmp => shr_megan_linkedlist
    do while(associated(meg_cmp))
       allocate(meg_cmp%emis_factors(numpft))
       call megan_factors_get( trim(meg_cmp%name), factors, class_num, molec_wght )
       meg_cmp%emis_factors = factors
       meg_cmp%class_number = class_num
       meg_cmp%molec_weight = molec_wght
       meg_cmp => meg_cmp%next_megcomp
    enddo

  end subroutine VOCEmission_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_map_EF
!
! !INTERFACE:
function get_map_EF(ivt_in,g_in)
!
! Get mapped EF for isoprene
! Use gridded values for 6 PFTs specified by MEGAN following
! Guenther et al. (2006).  Map the numpft CLM PFTs to these 6.
! Units: [ug m-2 h-1] 
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!

! !USES:
  use clmtype
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:

    ! varibles in
    integer, intent(in) :: ivt_in
    integer, intent(in) :: g_in
    real(r8)            :: get_map_EF

    real(r8), pointer :: efisop(:,:)      ! emission factors for isoprene for each pft [ug m-2 h-1]

    ! assign local pointer
    efisop     => gve%efisop
!-----------------------------------------------------------------------
    get_map_EF = 0._r8
    
    if (     ivt_in == ndllf_evr_tmp_tree  &
         .or.     ivt_in == ndllf_evr_brl_tree) then     !fineleaf evergreen
       get_map_EF = efisop(2,g_in)
    else if (ivt_in == ndllf_dcd_brl_tree) then     !fineleaf deciduous
       get_map_EF = efisop(3,g_in)
    else if (ivt_in >= nbrdlf_evr_trp_tree &
         .and.    ivt_in <= nbrdlf_dcd_brl_tree) then    !broadleaf trees
       get_map_EF = efisop(1,g_in)
    else if (ivt_in >= nbrdlf_evr_shrub &
         .and.    ivt_in <= nbrdlf_dcd_brl_shrub) then   !shrubs
       get_map_EF = efisop(4,g_in)
    else if (ivt_in >= nc3_arctic_grass &
         .and.    ivt_in <= nc4_grass) then              !grass
       get_map_EF = efisop(5,g_in)
    else if (ivt_in >= nc3crop) then                !crops
       get_map_EF =efisop(6,g_in)
    end if

end function get_map_EF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_gamma_P
!
! !INTERFACE:
  function get_gamma_P(par_sun_in, par24_sun_in, par240_sun_in, par_sha_in, par24_sha_in, par240_sha_in, &
       fsun_in, fsun240_in, forc_solad240_in,forc_solai240_in, LDF_in, cp, alpha) 
  
! Activity factor for PPFD (Guenther et al., 2006): all light dependent species
!-------------------------
! With distinction between sunlit and shaded leafs, weight scalings by
! fsun and fshade 
! Scale total incident par by fraction of sunlit leaves (added on 1/2002)
    
! fvitt -- forc_solad240, forc_solai240 can be zero when CLM finidat is specified
!          which will cause par240 to be zero and produce NaNs via log(par240)
! dml   -- fsun240 can be equal to or greater than one before 10 day averages are
!           set on startup or if a new pft comes online during land cover change.
!           Avoid this problem by only doing calculations with fsun240 when fsun240 is
!           between 0 and 1
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: par_sun_in
    real(r8),intent(in) :: par24_sun_in
    real(r8),intent(in) :: par240_sun_in
    real(r8),intent(in) :: par_sha_in
    real(r8),intent(in) :: par24_sha_in
    real(r8),intent(in) :: par240_sha_in
    real(r8),intent(in) :: fsun_in
    real(r8),intent(in) :: fsun240_in
    real(r8),intent(in) :: forc_solad240_in
    real(r8),intent(in) :: forc_solai240_in
    real(r8),intent(in) :: LDF_in

    real(r8),intent(out) :: cp                      ! temporary
    real(r8),intent(out) :: alpha                   ! temporary
    real(r8) :: gamma_p_LDF             ! activity factor for PPFD
    real(r8) :: get_gamma_P             ! return value

    real(r8), parameter :: ca1 = 0.004_r8                  ! empirical coefficent for alpha
    real(r8), parameter :: ca2 = 0.0005_r8                 ! empirical coefficent for alpha
    real(r8), parameter :: ca3 = 0.0468_r8                 ! empirical coefficent for cp
    real(r8), parameter :: par0_sun = 200._r8              ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: par0_shade = 50._r8             ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: alpha_fix = 0.001_r8            ! empirical coefficient
    real(r8), parameter :: cp_fix = 1.21_r8                ! empirical coefficient
!
! local pointers to implicit in arguments
!
!-----------------------------------------------------------------------

  if ( (fsun240_in > 0._r8) .and. (fsun240_in < 1._r8) .and.  (forc_solad240_in > 0._r8) &
       .and. (forc_solai240_in > 0._r8)) then
     ! With alpha and cp calculated based on eq 6 and 7:
     ! Note indexing for accumulated variables is all at pft level
     ! SUN:
     alpha = ca1 - ca2 * log(par240_sun_in)
     cp = ca3 * exp(ca2 * (par24_sun_in-par0_sun))*par240_sun_in**(0.6_r8)
     gamma_p_LDF = fsun_in * ( cp * alpha * par_sun_in * (1._r8 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5_r8) )
     ! SHADE:
     alpha = ca1 - ca2 * log(par240_sha_in)
     cp = ca3 * exp(ca2 * (par_sha_in-par0_shade))*par240_sha_in**(0.6_r8)
     gamma_p_LDF = gamma_p_LDF + (1._r8-fsun_in) * (cp*alpha*par_sha_in*(1._r8 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5_r8))
  else
     ! With fixed alpha and cp (from MEGAN User's Guide):
     ! SUN: direct + diffuse  
     alpha = alpha_fix
     cp = cp_fix
     gamma_p_LDF = fsun_in * ( cp * alpha*par_sun_in * (1._r8 + alpha*alpha*par_sun_in*par_sun_in)**(-0.5_r8) )
     ! SHADE: diffuse 
     gamma_p_LDF = gamma_p_LDF + (1._r8-fsun_in) * (cp*alpha*par_sha_in*(1._r8 + alpha*alpha*par_sha_in*par_sha_in)**(-0.5_r8))
  end if
  
  ! Calculate total activity factor for PPFD accounting for light-dependent fraction
  get_gamma_P = (1._r8 - LDF_in) + LDF_in * gamma_p_LDF

end function get_gamma_P
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_gamma_L
!
! !INTERFACE:
function get_gamma_L(fsun240_in,elai_in)
!
! Activity factor for LAI (Guenther et al., 2006): all species
! Guenther et al., 2006 eq 3
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!

! !USES:
    use clm_varcon   , only : denice
    use clm_varpar   , only : nlevsoi
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: fsun240_in
    real(r8),intent(in) :: elai_in
    real(r8)            :: get_gamma_L             ! return value
 

    ! parameters
    real(r8), parameter :: cce = 0.30_r8                   ! factor to set emissions to unity @ std
    real(r8), parameter :: cce1 = 0.24_r8                  ! same as Cce but for non-accumulated vars
!-----------------------------------------------------------------------
    if ( (fsun240_in > 0.0_r8) .and. (fsun240_in < 1.e30_r8) ) then 
       get_gamma_L = cce * elai_in
    else
       get_gamma_L = cce1 * elai_in
    end if

end function get_gamma_L
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_gamma_SM
!
! !INTERFACE:
function get_gamma_SM(clayfrac_in, sandfrac_in, h2osoi_vol_in, h2osoi_ice_in, dz_in, &
     bsw_in, watsat_in, sucsat_in, root_depth_in)
!
! Activity factor for soil moisture (Guenther et al., 2006): all species
!----------------------------------
! Calculate the mean scaling factor throughout the root depth.
! wilting point potential is in units of matric potential (mm) 
! (1 J/Kg = 0.001 MPa, approx = 0.1 m)
! convert to volumetric soil water using equation 7.118 of the CLM4 Technical Note
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!

! !USES:
    use clm_varcon   , only : denice
    use clm_varpar   , only : nlevsoi
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: clayfrac_in
    real(r8),intent(in) :: sandfrac_in
    real(r8),intent(in) :: h2osoi_vol_in(nlevsoi)
    real(r8),intent(in) :: h2osoi_ice_in(nlevsoi)
    real(r8),intent(in) :: dz_in(nlevsoi)
    real(r8),intent(in) :: bsw_in(nlevsoi)
    real(r8),intent(in) :: watsat_in(nlevsoi)
    real(r8),intent(in) :: sucsat_in(nlevsoi)
    real(r8),intent(in) :: root_depth_in
 
    real(r8)            :: get_gamma_SM

    ! local variables
    integer  :: j
    real(r8) :: nl                      ! temporary number of soil levels
    real(r8) :: theta_ice               ! water content in ice in m3/m3
    real(r8) :: wilt                    ! wilting point in m3/m3
    real(r8) :: theta1                  ! temporary

    ! parameters
    real(r8), parameter :: deltheta1=0.06_r8               ! empirical coefficient
    real(r8), parameter :: smpmax = 2.57e5_r8              ! maximum soil matrix potential

    if ((clayfrac_in > 0) .and. (sandfrac_in > 0)) then 
       get_gamma_SM = 0._r8
       nl=0._r8
       
       do j = 1,nlevsoi
          if  (sum(dz_in(1:j)) < root_depth_in)  then
             theta_ice = h2osoi_ice_in(j)/(dz_in(j)*denice)
             wilt = ((smpmax/sucsat_in(j))**(-1._r8/bsw_in(j))) * (watsat_in(j) - theta_ice)
             theta1 = wilt + deltheta1
             if (h2osoi_vol_in(j) >= theta1) then 
                get_gamma_SM = get_gamma_SM + 1._r8
             else if ( (h2osoi_vol_in(j) > wilt) .and. (h2osoi_vol_in(j) < theta1) ) then
                get_gamma_SM = get_gamma_SM + ( h2osoi_vol_in(j) - wilt ) / deltheta1
             else
                get_gamma_SM = get_gamma_SM + 0._r8
             end if
             nl=nl+1._r8
          end if
       end do
       
       if (nl > 0._r8) then
          get_gamma_SM = get_gamma_SM/nl
       endif

       if (get_gamma_SM > 1.0_r8) then 
          write(iulog,*) 'healdSM > 1: gamma_SM, nl', get_gamma_SM, nl
          get_gamma_SM=1.0_r8
       endif

    else
       get_gamma_SM = 1.0_r8
    end if


end function get_gamma_SM
!-----------------------------------------------------------------------
  
!-----------------------------------------------------------------------
! FUNCTION: get_gamma_T
!
! !INTERFACE:
function get_gamma_T(t_veg240_in, t_veg24_in,t_veg_in, ct1_in, ct2_in, betaT_in, LDF_in, Ceo_in, Eopt, topt)

! Activity factor for temperature 
!--------------------------------
! Calculate both a light-dependent fraction as in Guenther et al., 2006 for isoprene
! of a max saturation type form. Also caculate a light-independent fraction of the
! form of an exponential. Final activity factor depends on light dependent fraction
! of compound type.
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: t_veg240_in
    real(r8),intent(in) :: t_veg24_in
    real(r8),intent(in) :: t_veg_in
    real(r8),intent(in) :: ct1_in
    real(r8),intent(in) :: ct2_in
    real(r8),intent(in) :: betaT_in
    real(r8),intent(in) :: LDF_in
    real(r8),intent(in) :: Ceo_in
    real(r8),intent(out) :: Eopt                    ! temporary 
    real(r8),intent(out) :: topt                    ! temporary 

    ! local variables
    real(r8) :: get_gamma_T
    real(r8) :: gamma_t_LDF             ! activity factor for temperature
    real(r8) :: gamma_t_LIF             ! activity factor for temperature
    real(r8) :: x                       ! temporary 

    ! parameters
    real(r8), parameter :: co1 = 313._r8                   ! empirical coefficient
    real(r8), parameter :: co2 = 0.6_r8                    ! empirical coefficient
    real(r8), parameter :: co4 = 0.05_r8                   ! empirical coefficient
    real(r8), parameter :: tstd0 = 297_r8                  ! std temperature [K]
    real(r8), parameter :: topt_fix = 317._r8              ! std temperature [K]
    real(r8), parameter :: Eopt_fix = 2.26_r8              ! empirical coefficient
    real(r8), parameter :: ct3 = 0.00831_r8                ! empirical coefficient (0.0083 in User's Guide)
    real(r8), parameter :: tstd = 303.15_r8                ! std temperature [K]
    real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
!-----------------------------------------------------------------------

    ! Light dependent fraction (Guenther et al., 2006)
    if ( (t_veg240_in > 0.0_r8) .and. (t_veg240_in < 1.e30_r8) ) then 
       ! topt and Eopt from eq 8 and 9:
       topt = co1 + (co2 * (t_veg240_in-tstd0))
       Eopt = Ceo_in * exp (co4 * (t_veg24_in-tstd0)) * exp(co4 * (t_veg240_in -tstd0))
    else
       topt = topt_fix
       Eopt = Eopt_fix
    endif
    x = ( (1._r8/topt) - (1._r8/(t_veg_in)) ) / ct3
    gamma_t_LDF = Eopt * ( ct2_in * exp(ct1_in * x)/(ct2_in - ct1_in * (1._r8 - exp(ct2_in * x))) )
    
    
    ! Light independent fraction (of exp(beta T) form)
    gamma_t_LIF = exp(betaT_in * (t_veg_in - tstd))
    
    ! Calculate total activity factor for light as a function of light-dependent fraction
    !--------------------------------
    get_gamma_T = (1-LDF_in)*gamma_T_LIF + LDF_in*gamma_T_LDF 

end function get_gamma_T
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_gamma_A
!
! !INTERFACE:
function get_gamma_A(ivt_in, elai_p_in,elai_in,nclass_in)

! Activity factor for leaf age (Guenther et al., 2006)
!-----------------------------
! If not CNDV elai is constant therefore gamma_a=1.0
! gamma_a set to unity for evergreens (PFTs 1, 2, 4, 5)
! Note that we assume here that the time step is shorter than the number of 
!days after budbreak required to induce isoprene emissions (ti=12 days) and 
! the number of days after budbreak to reach peak emission (tm=28 days)
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (4/27/11)
!
! !ARGUMENTS:
    implicit none
! !LOCAL VARIABLES:

    ! varibles in
    integer,intent(in)  :: ivt_in
    integer,intent(in)  :: nclass_in
    real(r8),intent(in) :: elai_p_in
    real(r8),intent(in) :: elai_in

    real(r8)            :: get_gamma_A

    ! local variables
    real(r8) :: elai_prev               ! lai for previous timestep
    real(r8) :: fnew, fgro, fmat, fold  ! fractions of leaves at different phenological stages
    !-----------------------------------------------------------------------
    if ( (ivt_in == ndllf_dcd_brl_tree) .or. (ivt_in >= nbrdlf_dcd_trp_tree) ) then  ! non-evergreen
       
       if ( (elai_p_in > 0.0_r8) .and. (elai_p_in < 1.e30_r8) )then 
          elai_prev = 2._r8*elai_p_in-elai_in  ! have accumulated average lai over last timestep
          if (elai_prev == elai_in) then
             fnew = 0.0_r8
             fgro = 0.0_r8
             fmat = 1.0_r8
             fold = 0.0_r8
          else if (elai_prev > elai_in) then
             fnew = 0.0_r8
             fgro = 0.0_r8
             fmat = 1.0_r8 - (elai_prev - elai_in)/elai_prev
             fold = (elai_prev - elai_in)/elai_prev
          else if (elai_prev < elai_in) then
             fnew = 1 - (elai_prev / elai_in)
             fgro = 0.0_r8
             fmat = (elai_prev / elai_in)
             fold = 0.0_r8
          end if
          
          get_gamma_A = fnew*Anew(nclass_in) + fgro*Agro(nclass_in) + fmat*Amat(nclass_in) + fold*Aold(nclass_in)

       else
          get_gamma_A = 1.0_r8
       end if
       
    else
       get_gamma_A = 1.0_r8
    end if
    

  end function get_gamma_A
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FUNCTION: get_gamma_C
!
! !INTERFACE:
function get_gamma_C(cisun_in,cisha_in,forc_pbot_in,fsun_in)

! Activity factor for instantaneous CO2 changes (Heald et al., 2009)
!-------------------------
! With distinction between sunlit and shaded leafs, weight scalings by
! fsun and fshade 
!
! !CALLED FROM: VOCEmission
!
! !REVISION HISTORY:
! Author: Colette L. Heald (11/30/11)
!
! !USES:
    use clm_varctl,    only : co2_ppmv      ! corresponds to CCSM_CO2_PPMV set in env_conf.xml
!
! !ARGUMENTS:
    implicit none
! !LOCAL VARIABLES:

    ! varibles in
    real(r8),intent(in) :: cisun_in
    real(r8),intent(in) :: cisha_in
    real(r8),intent(in) :: forc_pbot_in
    real(r8),intent(in) :: fsun_in

    real(r8)            :: get_gamma_C

    ! local variables
    real(r8)            :: IEmin            ! empirical coeff for CO2 
    real(r8)            :: IEmax            ! empirical coeff for CO2 
    real(r8)            :: ECi50            ! empirical coeff for CO2 
    real(r8)            :: Cislope          ! empirical coeff for CO2 
    real(r8)            :: fint             ! interpolation fraction for CO2
    real(r8)            :: ci               ! temporary sunlight/shade weighted cisun & cisha (umolCO2/mol)

    !-----------------------------------------------------------------------

    ! Determine long-term CO2 growth environment (ie. ambient CO2) and interpolate
    ! parameters
    if ( co2_ppmv < 400._r8 ) then
       IEmin   = 0.7301_r8
       IEmax   = 1.0542_r8
       ECi50   = 457._r8
       Cislope = 3.1665_r8
    else if ( (co2_ppmv > 400._r8) .and. (co2_ppmv < 500._r8) ) then
       fint    = (co2_ppmv - 400._r8)/100._r8
       IEmin   = fint*0.7301_r8 + (1._r8 - fint)*0.7034_r8
       IEmax   = fint*1.0542_r8 + (1._r8 - fint)*0.9897_r8
       ECi50   = fint*457._r8 + (1._r8 - fint)*472._r8
       Cislope = fint*3.1665_r8 + (1._r8 - fint)*3.0652_r8
    else if ( (co2_ppmv > 500._r8) .and. (co2_ppmv < 600._r8) ) then
       fint = (co2_ppmv - 500._r8)/100._r8
       IEmin   = fint*0.7034_r8 + (1._r8 - fint)*0.6768_r8
       IEmax   = fint*0.9897_r8 + (1._r8 - fint)*0.9253_r8
       ECi50   = fint*472._r8 + (1._r8 - fint)*488._r8
       Cislope = fint*3.0652_r8 + (1._r8 - fint)*2.9321_r8
    else if ( (co2_ppmv > 600._r8) .and. (co2_ppmv < 700._r8) ) then
       fint = (co2_ppmv - 600._r8)/100._r8
       IEmin   = fint*0.6768_r8 + (1._r8 - fint)*0.6500_r8
       IEmax   = fint*0.9253_r8 + (1._r8 - fint)*0.8611_r8
       ECi50   = fint*488._r8 + (1._r8 - fint)*508._r8
       Cislope = fint*2.9321_r8 + (1._r8 - fint)*2.7497_r8
    else if ( (co2_ppmv > 700._r8) .and. (co2_ppmv < 800._r8) ) then
       fint = (co2_ppmv - 700._r8)/100._r8
       IEmin   = fint*0.6500_r8 + (1._r8 - fint)*0.6063_r8
       IEmax   = fint*0.8611_r8 + (1._r8 - fint)*0.7976_r8
       ECi50   = fint*508._r8 + (1._r8 - fint)*575._r8
       Cislope = fint*2.7497_r8 + (1._r8 - fint)*2.3643_r8
    else if ( co2_ppmv > 800._r8 ) then
       IEmin   = 0.6063_r8
       IEmax   = 0.7976_r8
       ECi50   = 575._r8
       Cislope = 2.3643_r8
    end if

    ! Intracellular CO2 concentrations (ci) given in Pa, divide by atmos
    ! pressure to get mixing ratio (umolCO2/mol)
    if ( (cisun_in > 0._r8) .and. (cisha_in > 0._r8) .and. (forc_pbot_in > 0._r8) .and. (fsun_in > 0._r8) ) then
       ci = ( fsun_in*cisun_in + (1._r8-fsun_in)*cisha_in )/forc_pbot_in * 1.e6_r8
       get_gamma_C = IEmin + ( (IEmax-IEmin)/(1._r8+(ci/ECi50)**Cislope) )
    else if ( (cisha_in > 0._r8) .and. (forc_pbot_in > 0._r8) ) then
       ci = cisha_in/forc_pbot_in * 1.e6_r8
       get_gamma_C = IEmin + ( (IEmax-IEmin)/(1._r8+(ci/ECi50)**Cislope) )
    else if ( (cisun_in > 0._r8) .and. (forc_pbot_in > 0._r8) ) then
       ci = cisun_in/forc_pbot_in * 1.e6_r8
       get_gamma_C = IEmin + ( (IEmax-IEmin)/(1._r8+(ci/ECi50)**Cislope) )
    else
       get_gamma_C = 1._r8
    end if

  end function get_gamma_C
!-----------------------------------------------------------------------

end module VOCEmissionMod


