module clm_varpar
  ! !PUBLIC TYPES:
  implicit none
  save

  integer, parameter :: nlevsno     =   5     ! maximum number of snow layers
  integer, parameter :: nlevtrc_soil = 10
  integer  :: nlevsoi = 10
  integer  :: nlevgrnd = 15
  integer  :: maxpatch_glcmec= -1    !some dumb number for stand alone betr
  integer            :: nlevlak               ! number of lake layers
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: ngases      =   3     ! CH4, O2, & CO2, this is for centurybgc
  integer, parameter :: mxpft       =  24     ! maximum number of PFT's for any mode;
  integer, parameter :: nsoilorder  =  15     ! number of soil orders

  integer :: natpft_size        ! Number of Patches on natural veg landunit (including bare ground)
  integer :: nlevdecomp_full = nlevtrc_soil
  integer :: cft_size           ! Number of Patches on crop landunit
  integer :: numpft
  logical :: crop_prog   = .false.  ! If prognostic crops is turned on
  integer :: nlevdecomp  = nlevtrc_soil
  integer :: ndecomp_pools
  integer :: ndecomp_cascade_transitions
  integer :: i_met_lit = 1
  integer :: i_cel_lit = 2
  integer :: i_lig_lit = 3
  integer :: i_cwd = 4

  integer, parameter :: ivis        =   1     ! index for visible band
  integer, parameter :: inir        =   2     ! index for near-infrared band
  integer :: maxpatch_pft        ! max number of plant functional types in naturally vegetated landunit (namelist setting)
end module clm_varpar
