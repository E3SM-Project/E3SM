!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_GridVertMod

!BOP
! !MODULE: POP_GridVertMod
!
! !DESCRIPTION:
!  This module contains vertical grid data types and common vertical
!  grid parameters.
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

   character ( 7), parameter, public :: POP_gridVertKindUnknown = 'Unknown'
   character ( 5), parameter, public :: POP_gridVertKindZ       = 'Depth'
   !character ( 3), parameter, public :: POP_gridVertKindALE     = 'ALE'
   !character ( 6), parameter, public :: POP_gridVertKindIsopyc  = 'Isopyc'

   !*** identifiers for grid locations

   character ( 7), parameter, public :: POP_gridVertLocUnknown  = 'Unknown'
   character ( 4), parameter, public :: POP_gridVertLocZtop     = 'Ztop'
   character ( 4), parameter, public :: POP_gridVertLocZmid     = 'Zmid'
   character ( 4), parameter, public :: POP_gridVertLocZbot     = 'Zbot'
   !character ( 8), parameter, public :: POP_gridVertLocLayerTop = 'LayerTop'
   !character ( 8), parameter, public :: POP_gridVertLocLayerMid = 'LayerMid'
   !character ( 8), parameter, public :: POP_gridVertLocLayerBot = 'LayerBot'

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

 end module POP_GridVertMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
