module clm_varpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! -------------------------------------------------------
! Module Parameters
! -------------------------------------------------------

  logical, public :: more_vertlayers = .false. ! true => run with more vertical soil layers


! Note - model resolution is read in from the surface dataset

  integer, parameter :: nlev_equalspace   = 15
  integer, parameter :: toplev_equalspace =  6
  integer            :: nlevsoi               ! number of hydrologically active soil layers
  integer            :: nlevsoifl             ! number of soil layers on input file
  integer            :: nlevgrnd              ! number of ground layers 
                                              ! (includes lower layers that are hydrologically inactive)
  integer            :: nlevurb               ! number of urban layers

#ifndef EXTRALAKELAYERS
  integer, parameter :: nlevlak     =  10     ! number of lake layers
#else
  ! Yields better results for site simulations.
  integer, parameter :: nlevlak     =  25     ! number of lake layers
#endif
  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  ! For CH4 code
  integer, parameter :: ngases = 3 ! CH4, O2, & CO2
  integer, parameter :: nlevcan     =   1     ! number of leaf layers in canopy
  integer, parameter :: numwat      =   5     ! number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: ivis        =   1     ! index for visible band
  integer, parameter :: inir        =   2     ! index for near-infrared band
  integer, parameter :: numsolar    =   2     ! number of solar type bands: direct, diffuse
  integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3     ! number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200     ! number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: mxpft       =  24     ! maximum number of PFT's for any mode; might we set some of these automatically from reading pft-physiology?
  integer, parameter :: numveg      =  16     ! number of veg types (without specific crop)
#if (defined VICHYDRO)
  integer, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang
  integer, parameter :: nlayert     =   8     ! number of VIC soil layer + 3 lower thermal layers
#endif
#if (defined CROP)
  integer, parameter :: numpft      = mxpft   ! actual # of pfts (without bare)
  integer, parameter :: numcft      =  10     ! actual # of crops
#else
  integer, parameter :: numpft      = numveg  ! actual # of pfts (without bare)
  integer, parameter :: numcft      =   2     ! actual # of crops
#endif
  integer, parameter :: maxpatch_pft= MAXPATCH_PFT ! max number of plant functional types in naturally vegetated landunit
  integer, parameter :: numurbl     = 3       ! number of urban landunits

  integer            :: nlevdecomp                    ! number of biogeochemically active soil layers
  integer            :: nlevdecomp_full               ! number of biogeochemical layers (includes lower layers that are biogeochemically inactive)

! -------------------------------------------------------
! Module Varaibles (initialized in clm_varpar_init)
! -------------------------------------------------------

#ifndef CENTURY_DECOMP
  ! parameters for decomposition cascade
  integer, parameter :: ndecomp_pools = 8
  integer, parameter :: ndecomp_cascade_transitions = 9
  integer, parameter :: i_met_lit = 1
  integer, parameter :: i_cel_lit = 2
  integer, parameter :: i_lig_lit = 3
  integer, parameter :: i_cwd = 4
  integer, parameter :: nsompools = 4
#else
  ! parameters for decomposition cascade
  integer, parameter :: ndecomp_pools = 7
  integer, parameter :: ndecomp_cascade_transitions = 10
  integer, parameter :: i_met_lit = 1
  integer, parameter :: i_cel_lit = 2
  integer, parameter :: i_lig_lit = 3
  integer, parameter :: i_cwd = 4
  integer, parameter :: nsompools = 3
#endif

! Indices used in surface file read and set in clm_varpar_init

  integer :: maxpatch           ! max number of patches
  integer :: maxpatch_glcmec    ! max number of elevation classes
  integer :: maxpatch_urb       ! max number of urban pfts (columns) in urban landunit
  integer :: npatch_urban       ! number of urban pfts (columns) in urban landunit
  integer :: npatch_lake        ! number of lake pfts (columns) in lake landunit
  integer :: npatch_wet         ! number of wetland pfts (columns) in wetland landunit
  integer :: npatch_glacier     ! number of glacier pfts (columns) in glacier landunit
  integer :: npatch_glacier_mec ! number of glacier_mec pfts (columns) in glacier_mec landunit
  integer :: max_pft_per_gcell 
  integer :: max_pft_per_lu 
  integer :: max_pft_per_col
  integer :: npatch_urban_tbd
  integer :: npatch_urban_hd
  integer :: npatch_urban_md

  real(r8) :: mach_eps            ! machine epsilon

! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!-----------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varpar_init
!
! !INTERFACE:
  subroutine clm_varpar_init()
!
! !DESCRIPTION:
! This subroutine initializes parameters in clm_varpar
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by T Craig
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  maxpatch_urb   = 5
  npatch_urban_tbd = maxpatch_pft + 1
  npatch_urban_hd  = npatch_urban_tbd + maxpatch_urb
  npatch_urban_md  = npatch_urban_hd + maxpatch_urb
  npatch_lake      = npatch_urban_md + maxpatch_urb
  npatch_wet     = npatch_lake  + 1
  npatch_glacier = npatch_wet   + 1
  npatch_glacier_mec = npatch_glacier + maxpatch_glcmec
  maxpatch       = npatch_glacier_mec
  mach_eps       = epsilon(1.0_r8)

  max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb*numurbl + maxpatch_glcmec
#if (defined CROP)
  max_pft_per_gcell = max_pft_per_gcell +  numcft  
#endif
  max_pft_per_lu    = max(numpft+1, numcft, maxpatch_urb)
  max_pft_per_col   = max(numpft+1, numcft, maxpatch_urb)

  nlevsoifl   =  10
  nlevurb     =  5
  if ( .not. more_vertlayers )then
     nlevsoi     =  nlevsoifl
     nlevgrnd    =  15
  else
     nlevsoi     =  8  + nlev_equalspace
     nlevgrnd    =  15 + nlev_equalspace
  end if

  ! here is a switch to set the number of soil levels for the biogeochemistry calculations.
  ! currently it works on either a single level or on nlevsoi and nlevgrnd levels
#ifdef VERTSOILC
  nlevdecomp      = nlevsoi
  nlevdecomp_full = nlevgrnd
#else
  nlevdecomp      = 1
  nlevdecomp_full = 1
#endif

  end subroutine clm_varpar_init

!------------------------------------------------------------------------------
end module clm_varpar
