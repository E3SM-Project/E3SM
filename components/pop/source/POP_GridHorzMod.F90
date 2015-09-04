!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_GridHorzMod

!BOP
! !MODULE: POP_GridHorzMod
!
! !DESCRIPTION:
!  This module contains the basic data types and common parameters for 
!  POP horizontal grids.  It defines supported horizontal grids and
!  basic functions for the grid type.
!
! !REVISION HISTORY:
!  SVN:$Id$
!  2007-02-09: Phil Jones
!              Initial implementation with some basic types

! !USES:

   use POP_KindsMod
   use POP_ErrorMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

! !PUBLIC DATA TYPES:

! !DEFINED PARAMETERS:

   !*** identifiers for commonly-used POP grids

   character ( 7), parameter, public :: POP_gridHorzKindUnknown  = 'Unknown'
   character ( 6), parameter, public :: POP_gridHorzKindLatLon   = 'LatLon'
   character (13), parameter, public :: POP_gridHorzKindDisplacedPole = 'DisplacedPole'
   character ( 7), parameter, public :: POP_gridHorzKindTripole  = 'Tripole'
   character ( 8), parameter, public :: POP_gridHorzKindGeodesic = 'Geodesic'

   !*** identifiers for types of staggering

   character ( 7), parameter, public :: POP_gridHorzStaggerUnknown = 'Unknown'
   character ( 1), parameter, public :: POP_gridHorzStaggerA       = 'A'
   character ( 4), parameter, public :: POP_gridHorzStaggerB_NE    = 'B_NE'
   character ( 4), parameter, public :: POP_gridHorzStaggerC_NE    = 'C_NE'
   character ( 1), parameter, public :: POP_gridHorzStaggerZ       = 'Z'

   !*** identifiers for common staggering locations

   character ( 7), parameter, public :: POP_gridHorzLocUnknown  = 'Unknown'
   character ( 6), parameter, public :: POP_gridHorzLocCenter   = 'Center'
   character ( 5), parameter, public :: POP_gridHorzLocNface    = 'NFace'
   character ( 5), parameter, public :: POP_gridHorzLocSface    = 'SFace'
   character ( 5), parameter, public :: POP_gridHorzLocEface    = 'EFace'
   character ( 5), parameter, public :: POP_gridHorzLocWface    = 'WFace'
   character ( 8), parameter, public :: POP_gridHorzLocNEcorner = 'NECorner'
   character ( 8), parameter, public :: POP_gridHorzLocNWcorner = 'NWCorner'
   character ( 8), parameter, public :: POP_gridHorzLocSEcorner = 'SECorner'
   character ( 8), parameter, public :: POP_gridHorzLocSWcorner = 'SWCorner'
   character ( 8), parameter, public :: POP_gridHorzLocNoUpdate = 'NoUpdate'

!EOP
!BOC
!EOC
!***********************************************************************

 !contains

!***********************************************************************
!BOP
! !IROUTINE: 
! !INTERFACE:


! !DESCRIPTION:
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!EOC

!***********************************************************************

 end module POP_GridHorzMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
