#include "pilgrim.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: scatter --- wrapper for PILGRIM utility
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine scatter( input, decomp, output, incomm )
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
      use decompmodule, only : decomptype
#if defined( SPMD )
      use parutilitiesmodule, only : parscatter 
#endif
      implicit none

! !INPUT PARAMETERS:
      real(CPP_REAL8) input(*)              ! Input array (global)
      type (decomptype) :: decomp          ! Decomposition
      integer incomm                       ! Communicator

! !OUTPUT PARAMETERS:
      real(CPP_REAL8) output(*)             ! Output array (local)

! !DESCRIPTION:
!     Scatter the global input array (on PE 0) according to 
!     decomposition to the local output array.  This intermediate
!     routine is a way to trick the compiler into passing 1-D 
!     arrays to the parutilitiesmodule method parscatter.
!
! !REVISION HISTORY:
!   WS  00.11.28:  Creation
!   AAM 01.05.08:  Added communicator as input argument
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined( SPMD )
      call parscatter( incomm, 0, input, decomp, output )
#endif
      return
!EOC
      end subroutine scatter
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: gather --- wrapper for PILGRIM utility
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine gather( input, decomp, output )
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
      use decompmodule, only : decomptype
#if defined( SPMD )
      use parutilitiesmodule, only : commglobal, pargather
#endif
      implicit none

! !INPUT PARAMETERS:
      real(CPP_REAL8) input(*)              ! Input array (global)
      type (decomptype) :: decomp          ! Decomposition

! !OUTPUT PARAMETERS:
      real(CPP_REAL8) output(*)             ! Output array (local)

! !DESCRIPTION:
!     Gather the local input array according to decomposition
!     to the global output array (on PE 0).  This intermediate
!     routine is a way to trick the compiler into passing 1-D 
!     arrays to the parutilitiesmodule method pargather.
!
! !REVISION HISTORY:
!   WS  00.11.28:  Creation
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined( SPMD )
      call pargather( commglobal, 0, input, decomp, output )
#endif
      return
!EOC
      end subroutine gather
!-----------------------------------------------------------------------

#if defined( SPMD )
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: begintransfer --- wrapper for PILGRIM utility
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine begintransfer( pattern, input, output )
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
#if defined( SPMD )
      use parutilitiesmodule, only : parpatterntype, parbegintransfer, commglobal
#endif
      implicit none

! !INPUT PARAMETERS:
      type (parpatterntype) :: pattern     ! Decomposition
      real(CPP_REAL8) input(*)              ! Input array
! !OUTPUT PARAMETERS:
      real(CPP_REAL8) output(*)             ! Output array

! !DESCRIPTION:
!     Initiate an asynchronous collective transfer of the input
!     array to the output array as defined by the communication 
!     pattern.
!
! !REVISION HISTORY:
!   WS  01.03.11:  Creation
!EOP
!-----------------------------------------------------------------------
!BOC

#if defined( SPMD )
      call parbegintransfer( commglobal, pattern, input, output )
#endif
      return
!EOC
      end subroutine begintransfer
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: endtransfer --- wrapper for PILGRIM utility
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine endtransfer( pattern, input, output )
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
#if defined( SPMD )
      use parutilitiesmodule, only : parpatterntype, parendtransfer, commglobal
#endif
      implicit none

! !INPUT PARAMETERS:
      type (parpatterntype) :: pattern     ! Decomposition
      real(CPP_REAL8) input(*)              ! Input array

! !OUTPUT PARAMETERS:
      real(CPP_REAL8) output(*)             ! Output array

! !DESCRIPTION:
!     Complete an asynchronous collective transfer of the input
!     array to the output array as defined by the communication 
!     pattern.
!
! !REVISION HISTORY:
!   WS  01.03.11:  Creation
!EOP
!---------------------------------------------------------------------
!BOC

#if defined( SPMD )
      call parendtransfer( commglobal, pattern, input, output )
#endif
      return
!EOC
      end subroutine endtransfer
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: exchangevector --- wrapper for PILGRIM utility
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine exchangevector( incomm, inlen, input, outlen, output )
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
#if defined( SPMD )
      use parutilitiesmodule, only : parexchangevector
#endif
      implicit none

! !INPUT PARAMETERS:
      integer incomm                       ! Communicator
      integer  inlen(*)                    ! Input lengths per PE
      real(CPP_REAL8) input(*)              ! Input array (global)

! !OUTPUT PARAMETERS:
      integer  outlen(*)                   ! Output lengths per PE
      real(CPP_REAL8) output(*)             ! Output array (local)

! !DESCRIPTION:
!     Perform a synchronous collective transfer of the input vector,
!     blocked in segments to be sent to each PE in ascending order, and 
!     the lengths of the blocks given by inlen.  The routine returns the 
!     output block lengths (those received on the local PE) and
!     the vector output which is blocked by ascending PE order.
!
! !REVISION HISTORY:
!   WS  01.03.11:  Creation
!EOP
!---------------------------------------------------------------------
!BOC

#if defined( SPMD )
      call parexchangevector( incomm, inlen, input, outlen, output )
#endif
      return
!EOC
      end subroutine exchangevector
!---------------------------------------------------------------------
#endif
