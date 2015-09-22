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
! Define land surface 2-d grid. The model resolution is read in from the surface dataset
!
  integer             :: lsmlon              ! maximum number of longitude points on lsm grid
  integer             :: lsmlat              ! number of latitude points on lsm grid

! Define number of levels

  integer, parameter :: nlevsoi     =  10     ! number of hydrologically active soil layers
  integer, parameter :: nlevgrnd    =  15     ! number of ground layers (includes lower layers that are hydrologically inactive)
  integer, parameter :: nlevurb     = nlevgrnd! number of urban layers (must equal nlevgrnd right now)
  integer, parameter :: nlevlak     =  10     ! number of lake layers
  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers

! Define miscellaneous parameters

  integer, parameter :: numwat      =   5   ! number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numrad      =   2   ! number of solar radiation bands: vis, nir
  integer, parameter :: ivis        =   1   ! index for visible band
  integer, parameter :: inir        =   2   ! index for near-infrared band
  integer, parameter :: numsolar    =   2   ! number of solar type bands: direct, diffuse
  integer, parameter :: ndst        =   4   ! number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3   ! number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200   ! number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: nvoc        =   5   ! number of voc categories

! Define parameters for RTM river routing model

  integer            :: rtmlon        !number of rtm longitudes
  integer            :: rtmlat        !number of rtm latitudes

! Define indices used in surface file read
! maxpatch_pft     = max number of plant functional types in naturally vegetated landunit
! maxpatch_urb     = max number of urban pfts (columns) in urban landunit
! maxpatch_wet     = max number of wetland pfts (columns) in wetland landunit
! maxpatch_lake    = max number of lake pfts (columns) in lake landunit
! maxpatch_glacier = max number of glacier pfts (columns) in glacier landunit
! maxpatch_glcmec  = max number of glacier_mec pfts (columns) in glacier_mec landunit

  integer, parameter :: mxpft          = 20     ! maximum number of PFT's for any mode
  integer, parameter :: numveg         = 16     ! number of veg types (without specific crop)
#if (defined CROP)
  integer, parameter :: numpft         = mxpft  ! actual # of pfts (without bare)
  integer, parameter :: numcft         =  6     ! actual # of crops
#else
  integer, parameter :: numpft         = numveg ! actual # of pfts (without bare)
  integer, parameter :: numcft         =  2     ! actual # of crops
#endif
  integer, parameter :: maxpatch_urb   = 5
  integer            :: maxpatch_pft
#if defined(GLC_NEC_10)
  integer, parameter :: maxpatch_glcmec = 10
#elif defined(GLC_NEC_5)
  integer, parameter :: maxpatch_glcmec = 5
#elif defined(GLC_NEC_3)
  integer, parameter :: maxpatch_glcmec = 3
#elif defined(GLC_NEC_1)
  integer, parameter :: maxpatch_glcmec = 1
#else
  integer, parameter :: maxpatch_glcmec = 0
#endif
  integer            :: npatch_urban
  integer            :: npatch_lake 
  integer            :: npatch_wet  
  integer            :: npatch_glacier
  integer            :: npatch_glacier_mec
  integer            :: maxpatch    

! clm_varpar_init seems to do something similar; less prone to error to move
! these three lines there? (slevis)
#if (defined CROP)
  integer, parameter :: max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb + maxpatch_glcmec
#else
  integer, parameter :: max_pft_per_gcell = numpft+1 + 3 + maxpatch_urb + numcft + maxpatch_glcmec
#endif
  integer, parameter :: max_pft_per_lu    = max(numpft+1, numcft, maxpatch_urb)
  integer, parameter :: max_pft_per_col   = max(numpft+1, numcft, maxpatch_urb)

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
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  lsmlon         = 1
  lsmlat         = 1
  rtmlon         = 1
  rtmlat         = 1
  maxpatch_pft   = MAXPATCH_PFT
  npatch_urban   = maxpatch_pft + 1
  npatch_lake    = npatch_urban + maxpatch_urb
  npatch_wet     = npatch_lake  + 1
  npatch_glacier = npatch_wet   + 1
  npatch_glacier_mec = npatch_glacier + maxpatch_glcmec
  maxpatch       = npatch_glacier_mec

  end subroutine clm_varpar_init

!------------------------------------------------------------------------------
end module clm_varpar
