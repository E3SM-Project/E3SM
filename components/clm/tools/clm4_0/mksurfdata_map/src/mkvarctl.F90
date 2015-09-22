module mkvarctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkvarctl
!
! !DESCRIPTION:
! Module containing control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
  real(r8), public, parameter :: spval = 1.e36    ! special value

  logical, public    :: outnc_large_files     ! output files in 64-bit format for large files
  logical, public    :: outnc_double          ! output ALL data in files as 64-bit
  integer, public    :: outnc_dims = 2        ! only applicable to lat/lon grids
  logical, public    :: outnc_1d              ! true => output file is 1d  

  character(len= 32), public :: mksrf_gridnm   = ' '  ! name of grid to use on output file
  character(len=256), public :: mksrf_fgrid    = ' '  ! land grid file name to use 
  character(len=256), public :: mksrf_gridtype = ' '  ! land gridtype, global or reg
  character(len=256), public :: mksrf_fvegtyp  = ' '  ! vegetation data file name
  character(len=256), public :: mksrf_fsoitex  = ' '  ! soil texture data file name
  character(len=256), public :: mksrf_forganic = ' '  ! organic matter data file name
  character(len=256), public :: mksrf_fsoicol  = ' '  ! soil color data file name
  character(len=256), public :: mksrf_flakwat  = ' '  ! inland lake data file name
  character(len=256), public :: mksrf_fwetlnd  = ' '  ! inland wetlands data file name
  character(len=256), public :: mksrf_furban   = ' '  ! urban data file name
  character(len=256), public :: mksrf_firrig   = ' '  ! irrigated area data file name
  character(len=256), public :: mksrf_fglacier = ' '  ! glacier data file name
  character(len=256), public :: mksrf_furbtopo = ' '  ! urban topography data file name
  character(len=256), public :: mksrf_flndtopo = ' '  ! land topography data file name
  character(len=256), public :: mksrf_fmax     = ' '  ! fmax data file name
  character(len=256), public :: mksrf_flai     = ' '  ! lai data filename
  character(len=256), public :: mksrf_fdynuse  = ' '  ! ascii file containing names of dynamic land use files
  character(len=256), public :: mksrf_fvocef   = ' '  ! VOC Emission Factor data file name

  integer           , public :: numpft         = 16   ! number of plant types

  character(len=256), public :: map_fpft        = ' ' ! Mapping file for PFT
  character(len=256), public :: map_flakwat     = ' ' ! Mapping file for lake water
  character(len=256), public :: map_fwetlnd     = ' ' ! Mapping file for wetland water
  character(len=256), public :: map_fglacier    = ' ' ! Mapping file for glacier
  character(len=256), public :: map_fsoitex     = ' ' ! Mapping file for soil texture
  character(len=256), public :: map_fsoicol     = ' ' ! Mapping file for soil color
  character(len=256), public :: map_furban      = ' ' ! Mapping file for urban
  character(len=256), public :: map_furbtopo    = ' ' ! Mapping file for urban topography
  character(len=256), public :: map_flndtopo    = ' ' ! Mapping file for land topography
  character(len=256), public :: map_fmax        = ' ' ! Mapping file for soil frac max
  character(len=256), public :: map_forganic    = ' ' ! Mapping file for organic soil
  character(len=256), public :: map_fvocef      = ' ' ! Mapping file for VOC emission factors
  character(len=256), public :: map_flai        = ' ' ! Mapping file for LAI
  character(len=256), public :: map_fharvest    = ' ' ! Mapping file for harvesting
  character(len=256), public :: map_firrig      = ' ' ! Mapping file for irrigation
!
! Variables to override data read in with
! (This is mostly for single-point mode, but could be used for sensitivity studies)
!
  logical,  public   :: all_urban              ! output ALL data as 100% covered in urban
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
