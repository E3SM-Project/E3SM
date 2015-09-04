!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_GridDimMod

!BOP
! !MODULE: POP_GridDimMod
!
! !DESCRIPTION:
!  This module contains the a grid dimension data type and basic 
!  operations on that data type.  
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

   public :: POP_GridDimCreate, &
             POP_GridDimDestroy, &
             POP_GridDimGet

! !PUBLIC DATA TYPES:

   type, public :: POP_GridDim
      private

      integer (POP_i4) :: indx
   end type

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  actual grid dimension type and array of defined dimensions
!
!-----------------------------------------------------------------------

   type :: POP_GridDimType
      character(POP_charLength) :: name
      integer (POP_i4) :: length  ! 1 to n, but 0 means unlimited
      integer (POP_i4) :: start, stop, stride  ! For slicing and dicing
   end type
 
   integer (POP_i4), parameter :: POP_maxGridDims = 20
   
   type (POP_GridDimType), dimension(POP_maxGridDims) :: &
      POP_allGridDims              ! storage for all defined dims

   logical (POP_Logical), dimension(POP_maxGridDims) :: &
      POP_definedGridDims = .false.

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_GridDimCreate
! !INTERFACE:

   function POP_GridDimCreate (name, length, errorCode, &
                               start, stop, stride)     &
                               result(gridDim)

! !DESCRIPTION:
!  Creates a grid dimension type with basic dimension data. Coordinates
!  are not included in this type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      name                      ! name for dimension

   integer (POP_i4), intent(in) :: &
      length                    ! length of this dimension

   integer (POP_i4), intent(in), optional :: &
      start,                   &! optional starting index (default 1)
      stop,                    &! optional stopping index (default 1)
      stride                    ! optional stride (1 is default)

! !OUTPUT PARAMETERS:

   type (POP_GridDim) :: gridDim  ! created GridDim

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      indx              ! first free location in GridDim array

!-----------------------------------------------------------------------
!
!  initialize error code and return value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   gridDim%indx = 0

!-----------------------------------------------------------------------
!
!  find first free GridDim
!
!-----------------------------------------------------------------------

   indx = 0
   searchLoop: do n=1,POP_maxGridDims
      if (.not. POP_definedGridDims(n)) then
         indx = n
         exit searchLoop
      endif
   end do searchLoop

   if (indx == 0) then
      call POP_ErrorSet(errorCode, &
         'POP_GridDimCreate: defined GridDims exceeded maxGridDims')
      return
   endif

!-----------------------------------------------------------------------
!
!  set values
!
!-----------------------------------------------------------------------

   POP_allGridDims(indx)%name = ' '
   POP_allGridDims(indx)%name = name

   POP_allGridDims(indx)%length = length

   if (present(start)) then
      POP_allGridDims(indx)%start = start
   else
      POP_allGridDims(indx)%start = 1
   endif

   if (present(stop)) then
      POP_allGridDims(indx)%stop  = stop
   else
      POP_allGridDims(indx)%stop  = 1
   endif

   if (present(stride)) then
      POP_allGridDims(indx)%stride = stride
   else
      POP_allGridDims(indx)%stride = 1
   endif

!-----------------------------------------------------------------------
!
!  return index in GridDim type
!
!-----------------------------------------------------------------------

   gridDim%indx = indx

!-----------------------------------------------------------------------
!EOC

   end function POP_GridDimCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_GridDimDestroy
! !INTERFACE:

   subroutine POP_GridDimDestroy(GridDim, errorCode)

! !DESCRIPTION:
!  Creates a grid dimension type with basic dimension data. Coordinates
!  are not included in this type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_GridDim), intent(inout) :: GridDim ! GridDim to destroy

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      indx              ! index of GridDim in GridDim array

!-----------------------------------------------------------------------
!
!  initialize error code and return value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   indx = GridDim%indx
   GridDim%indx = 0

!-----------------------------------------------------------------------
!
!  clear values
!
!-----------------------------------------------------------------------

   POP_allGridDims(indx)%name = ' '
   POP_allGridDims(indx)%length = 0

   POP_allGridDims(indx)%start  = 0
   POP_allGridDims(indx)%stop   = 0
   POP_allGridDims(indx)%stride = 0

   POP_definedGridDims(indx) = .false.

!-----------------------------------------------------------------------
!EOC

   end subroutine POP_GridDimDestroy

!***********************************************************************
!BOP
! !IROUTINE: POP_GridDimGet
! !INTERFACE:

   subroutine POP_GridDimGet(GridDim, errorCode, name, &
                             length, start, stop, stride)

! !DESCRIPTION:
!  Retrieves dimension information from a grid dimension type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_GridDim), intent(in) :: GridDim ! GridDim containing info

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   integer (POP_i4), intent(out), optional :: &
      length, &! dimension length
      start,  &! start index of dimension
      stop,   &! stop  index of dimension
      stride   ! stride along dimension

   character (*), intent(out), optional :: &
      name     ! dimension name

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      indx              ! index of GridDim in GridDim array

!-----------------------------------------------------------------------
!
!  initialize error code and return value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   indx = GridDim%indx

   if (indx < 0 .or. indx > POP_maxGridDims) then
      call POP_ErrorSet(errorCode, &
           'POP_GridDimGet: attempt to retrieve info from bad Grid Dim')
      return
   endif

   if (.not. POP_definedGridDims(indx)) then
      call POP_ErrorSet(errorCode, &
        'POP_GridDimGet: attempt to retrieve from undefined Grid Dim')
      return
   endif

!-----------------------------------------------------------------------
!
!  set values if requested
!
!-----------------------------------------------------------------------

   if (present(name  )) name   = POP_allGridDims(indx)%name
   if (present(length)) length = POP_allGridDims(indx)%length
   if (present(start )) start  = POP_allGridDims(indx)%start
   if (present(stop  )) stop   = POP_allGridDims(indx)%stop
   if (present(stride)) stride = POP_allGridDims(indx)%stride

!-----------------------------------------------------------------------
!EOC

   end subroutine POP_GridDimGet

!***********************************************************************

 end module POP_GridDimMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
