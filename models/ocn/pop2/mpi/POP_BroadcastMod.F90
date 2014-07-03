!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_BroadcastMod

!BOP
! !MODULE: POP_BroadcastMod
! !DESCRIPTION:
!  This module contains routines for broadcasting scalar and array
!  data from one processor or task to the rest of the tasks. This
!  particular version contains MPI implementations of these routines.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_BroadcastMod.F90 66 2007-03-08 22:03:07Z pwjones $
!  2006-08-07: Phil Jones
!              Added new broadcast module following new name
!              conventions
!              Added vector character broadcast needed by NCAR
!
! !USES:

   use POP_KindsMod
   use POP_CommMod
   use POP_ErrorMod

   implicit none
   private
   save

   include 'mpif.h'

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
!  Broadcasts a scalar real8 variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), intent(inout) :: &
      scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   call MPI_BCAST(scalar, 1, POP_mpiR8, srcTask, &
                             POP_communicator, ierr)
   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarR8: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)
   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarR8: error in MPI barrier')
      return
   endif

!-----------------------------------------------------------------------
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
! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), intent(inout) :: &
      scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   call MPI_BCAST(scalar, 1, POP_mpiR4, srcTask, &
                             POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarR$: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarR4: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      scalar                ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   call MPI_BCAST(scalar, 1, MPI_INTEGER, srcTask, &
                             POP_communicator,ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarI4: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarI4: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), intent(inout) :: &
     scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     itmp,               &! local temporary
     ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   if (scalar) then
     itmp = 1
   else
     itmp = 0
   endif

   call MPI_BCAST(itmp, 1, MPI_INTEGER, srcTask, &
                           POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarLog: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastScalarLog: error in MPI Barrier')
      return
   endif

   if (itmp == 1) then
     scalar = .true.
   else
     scalar = .false.
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   character (*), intent(inout) :: &
     scalar               ! scalar to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     clength,            &! length of character
     ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   clength = len(scalar)

   call MPI_BCAST(scalar, clength, MPI_CHARACTER, srcTask, &
                                   POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                       'POP_BroadcastScalarChar: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                     'POP_BroadcastScalarChar: error in MPI Barrier')
      return
   endif

!--------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastScalarChar

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D1
! !INTERFACE:

subroutine POP_BroadcastArrayR8D1(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a vector real8 variable from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:), intent(inout) :: &
      array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nelements,       &! size of array
      ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D1: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D1: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nelements,       &! size of array to be broadcast
      ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D1: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D1: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D1: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D1: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:), allocatable :: &
      tmpArray             ! temporary array for MPI bcast

   integer (POP_i4) :: &
      nelements,          &! size of array to be broadcast
      ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(nelements))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD1: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                     'POP_BroadcastArrayLogD1: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   character (POP_charLength), dimension(:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements*POP_charLength, MPI_CHARACTER, &
                         srcTask, POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                       'POP_BroadcastArrayCharD1: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                     'POP_BroadcastArrayCharD1: error in MPI Barrier')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayCharD1

!***********************************************************************
!BOP
! !IROUTINE: POP_BroadcastArrayR8D2
! !INTERFACE:

 subroutine POP_BroadcastArrayR8D2(array, srcTask, errorCode)

! !DESCRIPTION:
!  Broadcasts a real8 2d array from one processor (srcTask)
!  to all other processors. This is a specific instance of the generic
!  POP\_Broadcast interface.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nelements,         &! size of array
      ierr                ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D2: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR8D2: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D2: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR4D2: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:), intent(inout) :: &
       array              ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D2: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)
   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayI4D2: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:), allocatable :: &
     tmpArray             ! temporary array for MPI bcast

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(size(array,dim=1),size(array,dim=2)))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD2: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayLogD2: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D3: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR8D3: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D3: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR4D3: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D3: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayI4D3: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:,:), allocatable :: &
     tmpArray            ! temporary array for MPI bcast

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(size(array,dim=1), &
                     size(array,dim=2), &
                     size(array,dim=3)))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD3: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD3: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D4: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR8D4: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D4: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR4D4: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D4: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayI4D4: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:,:,:), allocatable :: &
     tmpArray            ! temporary array for MPI bcast

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(size(array,dim=1), &
                     size(array,dim=2), &
                     size(array,dim=3), &
                     size(array,dim=4)))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD4: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD4: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D5: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR8D5: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D5: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR4D5: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D5: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayI4D5: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:,:,:,:), allocatable :: &
     tmpArray            ! temporary array for MPI bcast

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(size(array,dim=1), &
                     size(array,dim=2), &
                     size(array,dim=3), &
                     size(array,dim=4), &
                     size(array,dim=5)))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD5: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD5: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:,:), intent(inout) :: &
     array             ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR8, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR8D6: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR8D6: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
     srcTask              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, POP_mpiR4, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayR4D6: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayR4D6: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, srcTask, &
                                    POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayI4D6: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                      'POP_BroadcastArrayI4D6: error in MPI Barrier')
      return
   endif

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

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask             ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:,:,:,:), intent(inout) :: &
      array               ! array to be broadcast

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:,:,:,:,:), allocatable :: &
     tmpArray            ! temporary array for MPI bcast

   integer (POP_i4) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!-----------------------------------------------------------------------

   errorCode = POP_success

   nelements = size(array)
   allocate(tmpArray(size(array,dim=1), &
                     size(array,dim=2), &
                     size(array,dim=3), &
                     size(array,dim=4), &
                     size(array,dim=5), &
                     size(array,dim=6)))

   where (array)
     tmpArray = 1
   elsewhere
     tmpArray = 0
   end where

   call MPI_BCAST(tmpArray, nelements, MPI_INTEGER, srcTask, &
                                       POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD6: error in MPI Bcast')
      return
   endif

   call MPI_BARRIER(POP_communicator, ierr)

   if (ierr /= MPI_SUCCESS) then
      call POP_ErrorSet(errorCode, &
                        'POP_BroadcastArrayLogD6: error in MPI Barrier')
      return
   endif

   where (tmpArray == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(tmpArray)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_BroadcastArrayLogD6

!***********************************************************************

 end module POP_BroadcastMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
