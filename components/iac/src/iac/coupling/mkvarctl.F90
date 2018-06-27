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
  use mkvarpar    , only: numpft
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

  character(len= 32), public :: mksrf_gridnm   = ' '  ! name of grid to use on output file
  character(len=256), public :: mksrf_fgrid    = ' '  ! land grid file name to use 
  character(len=256), public :: mksrf_gridtype = ' '  ! land gridtype, global or reg
  character(len=256), public :: mksrf_fvegtyp  = ' '  ! vegetation data file name
  character(len=256), public :: mksrf_fsoitex  = ' '  ! soil texture data file name
  character(len=256), public :: mksrf_forganic = ' '  ! organic matter data file name
  character(len=256), public :: mksrf_fsoicol  = ' '  ! soil color data file name
  character(len=256), public :: mksrf_flanwat  = ' '  ! inland water data file name
  character(len=256), public :: mksrf_furban   = ' '  ! urban data file name
  character(len=256), public :: mksrf_firrig   = ' '  ! irrigated area data file name
  character(len=256), public :: mksrf_fglacier = ' '  ! glacier data file name
  character(len=256), public :: mksrf_fglctopo = ' '  ! glc topography data file name
  character(len=256), public :: mksrf_flndtopo = ' '  ! land topography data file name
  character(len=256), public :: mksrf_ffrac    = ' '  ! grid fraction data file name
  character(len=256), public :: mksrf_fmax     = ' '  ! fmax data file name
  character(len=256), public :: mksrf_flai     = ' '  ! lai data filename
  character(len=256), public :: mksrf_fdynuse  = ' '  ! ascii file containing names of dynamic land use files
  character(len=256), public :: mksrf_fvocef   = ' '  ! VOC Emission Factor data file name

  character(len=256), public :: map_fpft        = ' '     
  character(len=256), public :: map_flanwat     = ' '  
  character(len=256), public :: map_fglacier    = ' ' 
  character(len=256), public :: map_fsoitex     = ' '  
  character(len=256), public :: map_fsoicol     = ' '  
  character(len=256), public :: map_furban      = ' '   
  character(len=256), public :: map_fglctopo    = ' '    
  character(len=256), public :: map_flndtopo    = ' '    
  character(len=256), public :: map_fmax        = ' '    
  character(len=256), public :: map_forganic    = ' ' 
  character(len=256), public :: map_fvocef      = ' '   
  character(len=256), public :: map_fglcmec_t2g = ' ' 
  character(len=256), public :: map_fglcmec_g2g = ' ' 
  character(len=256), public :: map_flai        = ' ' 
  character(len=256), public :: map_fharvest    = ' ' 
  character(len=256), public :: map_firrig      = ' ' 
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
