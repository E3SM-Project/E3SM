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

! Note - model resolution is read in from the surface dataset

  integer, parameter :: nlevsoi     =  10     ! number of hydrologically active soil layers
  integer, parameter :: nlevgrnd    =  15     ! number of ground layers (includes lower layers that are hydrologically inactive)
  integer, parameter :: nlevurb     = nlevgrnd! number of urban layers (must equal nlevgrnd right now)
  integer, parameter :: nlevlak     =  10     ! number of lake layers
  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer, parameter :: numwat      =   5     ! number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: ivis        =   1     ! index for visible band
  integer, parameter :: inir        =   2     ! index for near-infrared band
  integer, parameter :: numsolar    =   2     ! number of solar type bands: direct, diffuse
  integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3     ! number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200     ! number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: mxpft       =  20     ! maximum number of PFT's for any mode
  integer, parameter :: numveg      =  16     ! number of veg types (without specific crop)
#if (defined CROP)
  integer, parameter :: numpft      = mxpft   ! actual # of pfts (without bare)
  integer, parameter :: numcft      =   6     ! actual # of crops
#else
  integer, parameter :: numpft      = numveg  ! actual # of pfts (without bare)
  integer, parameter :: numcft      =   2     ! actual # of crops
#endif
  integer, parameter :: maxpatch_pft= MAXPATCH_PFT ! max number of plant functional types in naturally vegetated landunit

! -------------------------------------------------------
! Module Varaibles (initialized in clm_varpar_init)
! -------------------------------------------------------

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
  npatch_urban   = maxpatch_pft + 1
  npatch_lake    = npatch_urban + maxpatch_urb
  npatch_wet     = npatch_lake  + 1
  npatch_glacier = npatch_wet   + 1
  npatch_glacier_mec = npatch_glacier + maxpatch_glcmec
  maxpatch       = npatch_glacier_mec

  max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb + maxpatch_glcmec
#if (defined CROP)
  max_pft_per_gcell = max_pft_per_gcell +  numcft  
#endif
  max_pft_per_lu    = max(numpft+1, numcft, maxpatch_urb)
  max_pft_per_col   = max(numpft+1, numcft, maxpatch_urb)

  end subroutine clm_varpar_init

!------------------------------------------------------------------------------
end module clm_varpar
