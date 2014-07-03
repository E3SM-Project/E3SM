!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_BroadcastMod

!BOP
! !MODULE: POP_BroadcastMod
! !DESCRIPTION:
!  This module contains a set of broadcast routines for scalars
!  and arrays.  It supports broadcasting a value from one task
!  to all other tasks.  This particular module contains serial
!  implementations of these routines and typically perform no
!  operations since there is no need to broadcast what is already
!  known.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_BroadcastMod.F90 9 2006-07-14 23:02:40Z  $
!  2006-07-14: Phil Jones
!              New broadcast module following new name conventions
!                Added support for broadcasting arrays of all dimensions
!
! !USES:

   use POP_KindsMod
   use POP_CommMod
   use POP_ErrorMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: POP_Broadcast

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface POP_Broadcast
     module procedure POP_BroadcastScalarR8,    &
                      POP_BroadcastScalarR4,    &
                      POP_BroadcastScalarI4,    &
                      POP_BroadcastScalarLog,   &
                      POP_BroadcastScalarChar,  &
                      POP_BroadcastArrayR8D1,   &
                      POP_BroadcastArrayR4D1,   &
                      POP_BroadcastArrayI4D1,   &
                      POP_BroadcastArrayLogD1,  &
                      POP_BroadcastArrayCharD1, &
                      POP_BroadcastArrayR8D2,   &
                      POP_BroadcastArrayR4D2,   &
                      POP_BroadcastArrayI4D2,   &
                      POP_BroadcastArrayLogD2,  &
                      POP_BroadcastArrayR8D3,   &
                      POP_BroadcastArrayR4D3,   &
                      POP_BroadcastArrayI4D3,   &
                      POP_BroadcastArrayLogD3,  &
                      POP_BroadcastArrayR8D4,   &
                      POP_BroadcastArrayR4D4,   &
                      POP_BroadcastArrayI4D4,   &
                      POP_BroadcastArrayLogD4,  &
                      POP_BroadcastArrayR8D5,   &
                      POP_BroadcastArrayR4D5,   &
                      POP_BroadcastArrayI4D5,   &
                      POP_BroadcastArrayLogD5,  &
                      POP_BroadcastArrayR8D6,   &
                      POP_BroadcastArrayR4D6,   &
                      POP_BroadcastArrayI4D6,   &
                      POP_BroadcastArrayLogD6
   end interface

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastScalarR8
! !INTERFACE:

 subroutine POP_BroadcastScalarR8(scalar, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a scalar r8 variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), intent(inout) :: &
      scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!EOC

end subroutine POP_BroadcastScalarR8

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastScalarR4
! !INTERFACE:

subroutine POP_BroadcastScalarR4(scalar, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a scalar real variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), intent(inout) :: &
      scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastScalarR4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastScalarI4
! !INTERFACE:

subroutine POP_BroadcastScalarI4(scalar, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a scalar integer variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      scalar                ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastScalarI4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastScalarLog
! !INTERFACE:

subroutine POP_BroadcastScalarLog(scalar, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a scalar logical variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), intent(inout) :: &
     scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastScalarLog

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastScalarChar
! !INTERFACE:

subroutine POP_BroadcastScalarChar(scalar, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a scalar character variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   character (*), intent(inout) :: &
     scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastScalarChar

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D1
! !INTERFACE:

subroutine POP_BroadcastArrayR8D1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a vector r8 variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D1
! !INTERFACE:

subroutine POP_BroadcastArrayR4D1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real vector from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D1
! !INTERFACE:

subroutine POP_BroadcastArrayI4D1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts an integer vector from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD1
! !INTERFACE:

subroutine POP_BroadcastArrayLogD1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical vector from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayCharD1
! !INTERFACE:

subroutine POP_BroadcastArrayCharD1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a character vector from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   character (POP_charLength), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayCharD1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D2
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D2(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a r8 2d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D2

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D2
! !INTERFACE:

 subroutine POP_BroadcastArrayR4D2(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real 2d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D2

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D2
! !INTERFACE:

 subroutine POP_BroadcastArrayI4D2(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a 2d integer array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D2

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD2
! !INTERFACE:

 subroutine POP_BroadcastArrayLogD2(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical 2d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD2

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D3
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D3(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a double 3d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D3

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D3
! !INTERFACE:

 subroutine POP_BroadcastArrayR4D3(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real 3d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D3

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D3
! !INTERFACE:

 subroutine POP_BroadcastArrayI4D3(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts an integer 3d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D3

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD3
! !INTERFACE:

 subroutine POP_BroadcastArrayLogD3(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical 3d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD3

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D4
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D4(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a double 4d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D4
! !INTERFACE:

 subroutine POP_BroadcastArrayR4D4(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real 4d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D4
! !INTERFACE:

 subroutine POP_BroadcastArrayI4D4(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts an integer 4d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD4
! !INTERFACE:

 subroutine POP_BroadcastArrayLogD4(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical 4d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD4

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D5
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D5(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a double 5d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D5

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D5
! !INTERFACE:

 subroutine POP_BroadcastArrayR4D5(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real 5d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D5

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D5
! !INTERFACE:

 subroutine POP_BroadcastArrayI4D5(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts an integer 5d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D5

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD5
! !INTERFACE:

 subroutine POP_BroadcastArrayLogD5(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical 5d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD5

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D6
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D6(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a double 6d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR8D6

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR4D6
! !INTERFACE:

 subroutine POP_BroadcastArrayR4D6(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real 6d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayR4D6

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayI4D6
! !INTERFACE:

 subroutine POP_BroadcastArrayI4D6(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts an integer 6d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayI4D6

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayLogD6
! !INTERFACE:

 subroutine POP_BroadcastArrayLogD6(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a logical 6d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  for serial codes, nothing is required
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD6

!***********************************************************************

 end module POP_BroadcastMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
