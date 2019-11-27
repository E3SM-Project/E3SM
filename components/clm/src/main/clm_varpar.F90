module clm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing CLM parameters
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varctl   , only: use_extralakelayers, use_vertsoilc, use_crop, use_betr
  use clm_varctl   , only: use_century_decomp, use_c13, use_c14, use_fates
  use clm_varctl   , only: iulog, create_crop_landunit, irrigate
  use clm_varctl   , only: use_vichydro
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  logical, public :: more_vertlayers = .false. ! true => run with more vertical soil layers

  ! Note - model resolution is read in from the surface dataset

  integer, parameter :: nlev_equalspace   = 15
  integer, parameter :: toplev_equalspace =  6
  integer            :: nlevsoi               ! number of hydrologically active soil layers
  integer            :: nlevsoifl             ! number of soil layers on input file
  integer            :: nlevgrnd              ! number of ground layers 
                                              ! (includes lower layers that are hydrologically inactive)
  integer            :: nlevurb               ! number of urban layers
  integer            :: nlevlak               ! number of lake layers
  integer            :: nlevdecomp            ! number of biogeochemically active soil layers
  integer            :: nlevdecomp_full       ! number of biogeochemical layers 
                                              ! (includes lower layers that are biogeochemically inactive)
  integer            :: nlevtrc_soil
  integer            :: nlevtrc_full
  
  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer, parameter :: ngases      =   3     ! CH4, O2, & CO2
  integer, parameter :: nlevcan     =   1     ! number of leaf layers in canopy layer
  integer, parameter :: nvegwcs     =   4     ! number of vegetation water conductance segments
  integer, parameter :: numwat      =   5     ! number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: nmonth      =   12    ! number of months in year for crop planting
  integer, parameter :: ivis        =   1     ! index for visible band
  integer, parameter :: inir        =   2     ! index for near-infrared band
  integer, parameter :: numsolar    =   2     ! number of solar type bands: direct, diffuse
  integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3     ! number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200     ! number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: mxpft       =  24     ! maximum number of PFT's for any mode;
  ! FIX(RF,032414) might we set some of these automatically from reading pft-physiology?
  integer, parameter :: numveg      =  16     ! number of veg types (without specific crop)
  integer, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang
  integer            :: nlayert               ! number of VIC soil layer + 3 lower thermal layers

  integer :: numpft      = mxpft   ! actual # of patches (without bare)
  integer :: numcft      =  10     ! actual # of crops
  logical :: crop_prog   = .true.  ! If prognostic crops is turned on
  integer :: maxpatch_urb= 5       ! max number of urban patches (columns) in urban landunit

  integer :: maxpatch_pft        ! max number of plant functional types in naturally vegetated landunit (namelist setting)

  integer, parameter :: nsoilorder  =  15     ! number of soil orders

  integer, parameter :: nlevslp = 11          ! number of slope percentile levels

  ! constants for decomposition cascade

  integer :: i_met_lit 
  integer :: i_cel_lit 
  integer :: i_lig_lit 
  integer :: i_cwd 

  integer :: ndecomp_pools
  integer :: ndecomp_cascade_transitions

  ! Indices used in surface file read and set in clm_varpar_init

  integer :: natpft_lb          ! In PFT arrays, lower bound of Patches on the natural veg landunit (i.e., bare ground index)
  integer :: natpft_ub          ! In PFT arrays, upper bound of Patches on the natural veg landunit
  integer :: natpft_size        ! Number of Patches on natural veg landunit (including bare ground)
  integer :: cft_lb             ! In PFT arrays, lower bound of Patches on the crop landunit
  integer :: cft_ub             ! In PFT arrays, upper bound of Patches on the crop landunit
  integer :: cft_size           ! Number of Patches on crop landunit

  integer :: maxpatch_glcmec    ! max number of elevation classes
  integer :: max_patch_per_col

  real(r8) :: mach_eps            ! machine epsilon
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters
  !
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine clm_varpar_init()
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    character(len=32) :: subname = 'clm_varpar_init'  ! subroutine name
    !------------------------------------------------------------------------------

    ! Crop settings and consistency checks

    if (use_crop) then
       numpft      = mxpft   ! actual # of patches (without bare)
       numcft      =  10     ! actual # of crops
       crop_prog   = .true.  ! If prognostic crops is turned on
    else
       numpft      = numveg  ! actual # of patches (without bare)
       numcft      =   2     ! actual # of crops
       crop_prog   = .false. ! If prognostic crops is turned on
    end if

    ! For arrays containing all Patches (natural veg & crop), determine lower and upper bounds
    ! for (1) Patches on the natural vegetation landunit (includes bare ground, and includes
    ! crops if create_crop_landunit=false), and (2) CFTs on the crop landunit (no elements
    ! if create_crop_landunit=false)

    if (create_crop_landunit) then
       natpft_size = (numpft + 1) - numcft    ! note that numpft doesn't include bare ground -- thus we add 1
       cft_size    = numcft
    else
       natpft_size = numpft + 1               ! note that numpft doesn't include bare ground -- thus we add 1
       cft_size    = 0
    end if

    natpft_lb = 0
    natpft_ub = natpft_lb + natpft_size - 1
    cft_lb = natpft_ub + 1
    cft_ub = cft_lb + cft_size - 1

    max_patch_per_col= max(numpft+1, numcft, maxpatch_urb)
    mach_eps       = epsilon(1.0_r8)

    nlevsoifl   =  10
    nlevurb     =  5
    if ( .not. more_vertlayers )then
       nlevsoi     =  nlevsoifl
       nlevgrnd    =  15
    else
       nlevsoi     =  8  + nlev_equalspace
       nlevgrnd    =  15 + nlev_equalspace
    end if

    if (use_vichydro) then
       nlayert     =  nlayer + (nlevgrnd -nlevsoi)
    endif

    ! here is a switch to set the number of soil levels for the biogeochemistry calculations.
    ! currently it works on either a single level or on nlevsoi and nlevgrnd levels
    if (use_vertsoilc) then
       nlevdecomp      = nlevsoi
       nlevdecomp_full = nlevgrnd
    else
       nlevdecomp      = 1
       nlevdecomp_full = 1
    end if
    
    nlevtrc_full   = nlevsoi
    if(use_betr) then
      nlevtrc_soil = nlevsoi
    endif

    if (.not. use_extralakelayers) then
       nlevlak     =  10     ! number of lake layers
    else
       nlevlak     =  25     ! number of lake layers (Yields better results for site simulations)
    end if

    if ( use_fates ) then
       i_cwd = 0
       i_met_lit = 1
       i_cel_lit = 2
       i_lig_lit = 3
       if (use_century_decomp) then
          ndecomp_pools = 6
          ndecomp_cascade_transitions = 8
       else
          ndecomp_pools = 7
          ndecomp_cascade_transitions = 7
       end if
       
    else

       if (use_century_decomp) then
          ndecomp_pools = 7
          ndecomp_cascade_transitions = 10
          i_met_lit = 1
          i_cel_lit = 2
          i_lig_lit = 3
          i_cwd = 4
       else
          ndecomp_pools = 8
          ndecomp_cascade_transitions = 9
          i_met_lit = 1
          i_cel_lit = 2
          i_lig_lit = 3
          i_cwd = 4
       end if

    endif
    


  end subroutine clm_varpar_init

end module clm_varpar
