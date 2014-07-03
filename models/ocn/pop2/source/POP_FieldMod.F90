!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_FieldMod

!BOP
! !MODULE: POP_FieldMod
!
! !DESCRIPTION:
!  This module contains a basic field class with a field structure
!  containing data and meta-data in the form of attributes.  Attributes
!  are stored in attribute arrays, but the interfaces recognize the
!  standard attributes shortName, longName, units, numDims, fieldLoc,
!  fieldKind, fieldDims, missingValue, validRangeMin, validRangeMax.
!
! !REVISION HISTORY:
!  SVN:$Id$
!  2007-02-08: Phil Jones
!              Initial implementation

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_GridDimMod
   use POP_GridHorzMod
   use POP_GridVertMod

   implicit none
   private
   save

! !PUBLIC TYPES:

   ! Full field data type

   type, public :: POP_Field
      private

      !--- special standard attributes (all other standard attributes
      !--- are stored in the attribute arrays, but are recognized by
      !--- the field interfaces

      type (POP_GridDim), dimension(:), pointer   :: fieldDims

      !--- attributes of each data type

      integer  (POP_i4) :: nAtts  ! total number of attributes
      
      integer  (POP_i4) :: nAttsChar ! num of defined char attributes
      character(POP_CharLength), dimension(:), pointer  :: attribNameChar
      character(POP_CharLength), dimension(:), pointer  :: attribValChar

      integer  (POP_i4) :: nAttsLog  ! num of defined logical attributes
      character(POP_CharLength), dimension(:), pointer  :: attribNameLog
      logical  (POP_Logical),    dimension(:), pointer  :: attribValLog

      integer  (POP_i4) :: nAttsI4   ! num of defined integer attributes
      character(POP_CharLength), dimension(:), pointer  :: attribNameI4
      integer  (POP_i4),         dimension(:), pointer  :: attribValI4

      integer  (POP_i4) :: nAttsR4   ! num of defined real attributes
      character(POP_CharLength), dimension(:), pointer  :: attribNameR4
      real     (POP_r4),         dimension(:), pointer  :: attribValR4

      integer  (POP_i4) :: nAttsR8   ! num of defined double attributes
      character(POP_CharLength), dimension(:), pointer  :: attribNameR8
      real     (POP_r8),         dimension(:), pointer  :: attribValR8

      !   Only one of these next nine pointers can be associated.
      !   The others must be nullified.  For convenience in
      !   initialization, these declarations are the last listed
      !   in this type.

      logical(POP_Logical), dimension(:,:,:), pointer     :: data2DLog
      logical(POP_Logical), dimension(:,:,:,:), pointer   :: data3DLog
      logical(POP_Logical), dimension(:,:,:,:,:), pointer :: data4DLog
      integer(POP_i4), dimension(:,:,:), pointer          :: data2DI4
      integer(POP_i4), dimension(:,:,:,:), pointer        :: data3DI4
      integer(POP_i4), dimension(:,:,:,:,:), pointer      :: data4DI4
      real(POP_r4), dimension(:,:,:), pointer             :: data2DR4
      real(POP_r4), dimension(:,:,:,:), pointer           :: data3DR4
      real(POP_r4), dimension(:,:,:,:,:), pointer         :: data4DR4
      real(POP_r8), dimension(:,:,:), pointer             :: data2DR8
      real(POP_r8), dimension(:,:,:,:), pointer           :: data3DR8
      real(POP_r8), dimension(:,:,:,:,:), pointer         :: data4DR8

   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_FieldCreate,           &
             POP_FieldDestroy,          &
             POP_FieldGetNumAttributes, &
             POP_FieldGetDimSize,       &
             POP_FieldAttributeSet,     &
             POP_FieldAttributeGet,     &
             POP_FieldAttachData,       &
             POP_FieldDetachData,       &
             POP_FieldGetData

! !DEFINED PARAMETERS:

   !*** identifiers for commonly-used POP field types

   character (7), parameter, public :: POP_fieldKindUnknown  = 'unknown'
   character (6), parameter, public :: POP_fieldKindScalar   = 'scalar'
   character (6), parameter, public :: POP_fieldKindVector   = 'vector'
   character (5), parameter, public :: POP_fieldKindAngle    = 'angle'
   character (8), parameter, public :: POP_fieldKindNoUpdate = 'noUpdate'

   !*** identifiers for commonly-used POP data types

   character (7), parameter, public :: POP_fieldDataTypeUnknown = 'unknown'
   character (7), parameter, public :: POP_fieldDataTypeLogical = 'logical'
   character (2), parameter, public :: POP_fieldDataTypeI4      = 'i4'
   character (2), parameter, public :: POP_fieldDataTypeR4      = 'r4'
   character (2), parameter, public :: POP_fieldDataTypeR8      = 'r8'

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  module data types
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  generic interfaces
!
!-----------------------------------------------------------------------

   interface POP_FieldAttributeSet
      module procedure POP_FieldAttributeSetChar, &
                       POP_FieldAttributeSetLog,  &
                       POP_FieldAttributeSetI4,   &
                       POP_FieldAttributeSetR4,   &
                       POP_FieldAttributeSetR8
   end interface 

   interface POP_FieldAttributeGet
      module procedure POP_FieldAttributeGetChar, &
                       POP_FieldAttributeGetLog,  &
                       POP_FieldAttributeGetI4,   &
                       POP_FieldAttributeGetR4,   &
                       POP_FieldAttributeGetR8,   &
                       POP_FieldAttributeGetDims
   end interface 

   interface POP_FieldAttachData
      module procedure POP_FieldAttachData2DLog, &
                       POP_FieldAttachData3DLog, &
                       POP_FieldAttachData4DLog, &
                       POP_FieldAttachData2DI4 , &
                       POP_FieldAttachData3DI4 , &
                       POP_FieldAttachData4DI4 , &
                       POP_FieldAttachData2DR4 , &
                       POP_FieldAttachData3DR4 , &
                       POP_FieldAttachData4DR4 , &
                       POP_FieldAttachData2DR8 , &
                       POP_FieldAttachData3DR8 , &
                       POP_FieldAttachData4DR8
   end interface 

   interface POP_FieldGetData
      module procedure POP_FieldGetData2DLog, &
                       POP_FieldGetData3DLog, &
                       POP_FieldGetData4DLog, &
                       POP_FieldGetData2DI4 , &
                       POP_FieldGetData3DI4 , &
                       POP_FieldGetData4DI4 , &
                       POP_FieldGetData2DR4 , &
                       POP_FieldGetData3DR4 , &
                       POP_FieldGetData4DR4 , &
                       POP_FieldGetData2DR8 , &
                       POP_FieldGetData3DR8 , &
                       POP_FieldGetData4DR8
   end interface 

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldCreate
! !INTERFACE:

   function POP_FieldCreate (errorCode, shortName, longName, units,   &
                             horzLoc, vertLoc, fieldKind, dataType,   &
                             numDims, fieldDims)                      &
            result(field)

! !DESCRIPTION:
!  Creates a field type with field metadata, if supplied.  The actual 
!  field data is attached or detached with a separate function call.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in), optional :: &
      shortName,               &! short (one word) name for field
      longName,                &! longer descriptive name for field
      units,                   &! units for field
      horzLoc,                 &! horizontal staggering location
      vertLoc,                 &! vertical   staggering location
      dataType,                &! field data type (log, i4, r4, r8)
      fieldKind                 ! field type (scalar,vector,angle)

   integer (POP_i4), intent(in), optional :: &
      numDims                   ! num of spatial dimensions for field

   type (POP_GridDim), dimension(:), intent(in), optional :: &
      fieldDims                 ! grid dimension descriptor for each dim

   !*** missingValue, validRange not included as they are dependent
   !*** on the data type - use the FieldAttributeSet call for those

! !OUTPUT PARAMETERS:

   type (POP_Field) :: field  ! created field type

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      istat,           &! error flag for allocate
      nDims             ! number of dimensions

!-----------------------------------------------------------------------
!
!  Initialize all components of data type to defaults
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   nullify(field%fieldDims)

   field%nAtts = 0

   field%nAttsChar = 0
   nullify(field%attribNameChar)
   nullify(field%attribValChar)

   field%nAttsLog  = 0
   nullify(field%attribNameLog)
   nullify(field%attribValLog)

   field%nAttsI4   = 0
   nullify(field%attribNameI4)
   nullify(field%attribValI4)

   field%nAttsR4   = 0
   nullify(field%attribNameR4)
   nullify(field%attribValR4)

   field%nAttsR8   = 0
   nullify(field%attribNameR8)
   nullify(field%attribValR8)

   nullify(field%data2DLog)
   nullify(field%data3DLog)
   nullify(field%data4DLog)
   nullify(field%data2DI4)
   nullify(field%data3DI4)
   nullify(field%data4DI4)
   nullify(field%data2DR4)
   nullify(field%data3DR4)
   nullify(field%data4DR4)
   nullify(field%data2DR8)
   nullify(field%data3DR8)
   nullify(field%data4DR8)

!-----------------------------------------------------------------------
!
!  now fill metadata with input values
!
!-----------------------------------------------------------------------

   if (present(shortName)) then
      call POP_FieldAttributeSet(field, 'shortName', trim(shortName), &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'shortName', 'unknown', &
                                 errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing shortName')
      return
   endif

   if (present(longName)) then
      call POP_FieldAttributeSet(field, 'longName', trim(longName), &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'longName', 'unknown', &
                                 errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing longName')
      return
   endif

   if (present(units)) then
      call POP_FieldAttributeSet(field, 'units', trim(units), &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'units', 'unitless', &
                                 errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing units')
      return
   endif

   if (present(numDims)) then
      call POP_FieldAttributeSet(field, 'numDims', numDims, &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'numDims', 0, &
                                 errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing number of dimensions')
      return
   endif

   if (present(horzLoc)) then
      call POP_FieldAttributeSet(field, 'horzLoc', horzLoc, &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'horzLoc', &
                                 POP_gridHorzLocUnknown, errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing horizontal field loc')
      return
   endif

   if (present(vertLoc)) then
      call POP_FieldAttributeSet(field, 'vertLoc', vertLoc, &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'vertLoc', &
                                 POP_gridVertLocUnknown, errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing vertical field loc')
      return
   endif

   if (present(dataType)) then
      call POP_FieldAttributeSet(field, 'dataType', dataType, &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'dataType', &
                                 POP_fieldDataTypeUnknown, errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing field data type')
      return
   endif

   if (present(fieldKind)) then
      call POP_FieldAttributeSet(field, 'fieldKind', fieldKind, &
                                 errorCode)
   else
      call POP_FieldAttributeSet(field, 'fieldKind', &
                                 POP_fieldKindUnknown, errorCode)
   endif
  
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldCreate: error initializing field kind')
      return
   endif

   if (present(fieldDims)) then
      if (present(numDims))  then
         nDims = numDims
      else
         nDims = size(fieldDims)
      endif
     
      allocate(field%fieldDims(nDims), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                        'POP_FieldCreate: error allocating fieldDims')
         return
      else
         field%fieldDims(:) = fieldDims(:)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_FieldCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldDestroy
! !INTERFACE:

   subroutine POP_FieldDestroy(field, errorCode)

! !DESCRIPTION:
!  Destroys a field type and deallocates all memory associated with
!  the field.  Field data is not deallocated or destroyed, but the
!  pointer to that data is nullified.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to be destroyed

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
      istat             ! status flag for allocates

!-----------------------------------------------------------------------
!
!  Clear all components of data type.
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (associated(field%fieldDims)) &
      deallocate(field%fieldDims, stat=istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
                      'POP_FieldDestroy: error deallocating fieldDims')
      return
   else
      nullify(field%fieldDims)
   endif

   if (associated(field%attribNameChar)) then
      deallocate(field%attribNameChar, &
                 field%attribValChar, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                 'POP_FieldDestroy: error deallocating attribChar')
         return
      else
         nullify(field%attribNameChar)
         nullify(field%attribValChar)
      endif
   endif

   if (associated(field%attribNameLog)) then
      deallocate(field%attribNameLog, &
                 field%attribValLog, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                 'POP_FieldDestroy: error deallocating attribLog')
         return
      else
         nullify(field%attribNameLog)
         nullify(field%attribValLog)
      endif
   endif

   if (associated(field%attribNameI4)) then
      deallocate(field%attribNameI4, &
                 field%attribValI4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                 'POP_FieldDestroy: error deallocating attribI4')
         return
      else
         nullify(field%attribNameI4)
         nullify(field%attribValI4)
      endif
   endif


   if (associated(field%attribNameR4)) then
      deallocate(field%attribNameR4, &
                 field%attribValR4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                 'POP_FieldDestroy: error deallocating attribR4')
         return
      else
         nullify(field%attribNameR4)
         nullify(field%attribValR4)
      endif
   endif

   if (associated(field%attribNameR8)) then
      deallocate(field%attribNameR8, &
                 field%attribValR8, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
                 'POP_FieldDestroy: error deallocating attribR8')
         return
      else
         nullify(field%attribNameR8)
         nullify(field%attribValR8)
      endif
   endif

   nullify(field%data2DLog)
   nullify(field%data3DLog)
   nullify(field%data4DLog)
   nullify(field%data2DI4)
   nullify(field%data3DI4)
   nullify(field%data4DI4)
   nullify(field%data2DR4)
   nullify(field%data3DR4)
   nullify(field%data4DR4)
   nullify(field%data2DR8)
   nullify(field%data3DR8)
   nullify(field%data4DR8)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldDestroy

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetNumAttributes
! !INTERFACE:

   subroutine POP_FieldGetNumAttributes(field,    errorCode,           &
                                        nAtts,    nAttsChar, nAttsLog, &
                                        nAttsI4,  nAttsR4  , nAttsR8)

! !DESCRIPTION:
!  Retrieves the number of attributes in a field.  Can also retrieve
!  the number of attributes of any given type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field  ! field type to be queried

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   integer (POP_i4), intent(out), optional :: &
      nAtts,       &! total number of attributes
      nAttsChar,   &! number of character attributes
      nAttsLog,    &! number of logical   attributes
      nAttsI4,     &! number of integer   attributes
      nAttsR4,     &! number of real      attributes
      nAttsR8       ! number of double    attributes
   
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  get relevant attribute counters
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(nAtts)) then
      nAtts = field%nAtts
   endif
   
   if (present(nAttsChar)) then
      nAttsChar = field%nAttsChar
   endif
   
   if (present(nAttsLog)) then
      nAttsLog  = field%nAttsLog
   endif
   
   if (present(nAttsI4)) then
      nAttsI4   = field%nAttsI4
   endif
   
   if (present(nAttsR4  )) then
      nAttsR4   = field%nAttsR4  
   endif
   
   if (present(nAttsR8)) then
      nAttsR8   = field%nAttsR8
   endif
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetNumAttributes

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetDimSize
! !INTERFACE:

   subroutine POP_FieldGetDimSize(field, dimIndex, dimSize, errorCode)

! !DESCRIPTION:
!  Retrieves the size of each dimension in a field. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field  ! field type to be queried

   integer (POP_i4), intent(in) :: &
      dimIndex       ! index of dimension for which size requested

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode,   &! returned error code
      dimSize       ! size of the requested dimension
   
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_Logical) :: &
      foundData              ! flag for finding attached data

!-----------------------------------------------------------------------
!
!  initialize
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   foundData = .false.

   if (dimIndex < 1 .or. dimIndex > 7) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldGetDimSize: bad dimension index')
      return
   endif

!-----------------------------------------------------------------------
!
!  if grid dimensions have been defined, check them directly
!
!-----------------------------------------------------------------------

   if (associated(field%fieldDims)) then

      if (dimIndex > size(field%fieldDims)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldGetDimSize: dimension index > defined dims')
         return
      endif
      
      call POP_GridDimGet(field%fieldDims(dimIndex), errorCode, &
                          length = dimSize)
                          
      if (errorCode /= POP_Success) then
      	 call POP_ErrorSet(errorCode, &
      	    'POP_FieldGetDimSize: error retrieving size from grid dim')
      	 return
      endif

      foundData = .true.
      
!-----------------------------------------------------------------------
!
!  otherwise try to determine dimension size from attached data
!
!-----------------------------------------------------------------------

   else

      if (associated(field%data2DLog)) then

         if (dimIndex > 2) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data2DLog, dim=dimIndex)

      else if (associated(field%data3DLog)) then

         if (dimIndex > 3) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data3DLog, dim=dimIndex)

      else if (associated(field%data4DLog)) then

         if (dimIndex > 4) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data4DLog, dim=dimIndex)

      else if (associated(field%data2DI4)) then

         if (dimIndex > 2) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data2DI4, dim=dimIndex)

      else if (associated(field%data3DI4)) then

         if (dimIndex > 3) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data3DI4, dim=dimIndex)

      else if (associated(field%data4DI4)) then

         if (dimIndex > 4) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data4DI4, dim=dimIndex)

      else if (associated(field%data2DR4)) then

         if (dimIndex > 2) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data2DR4, dim=dimIndex)

      else if (associated(field%data3DR4)) then

         if (dimIndex > 3) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data3DR4, dim=dimIndex)

      else if (associated(field%data4DR4)) then

         if (dimIndex > 4) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data4DR4, dim=dimIndex)

      else if (associated(field%data2DR8)) then

         if (dimIndex > 2) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data2DR8, dim=dimIndex)

      else if (associated(field%data3DR8)) then

         if (dimIndex > 3) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data3DR8, dim=dimIndex)

      else if (associated(field%data4DR8)) then

         if (dimIndex > 4) then
            call POP_ErrorSet(errorCode, &
               'POP_FieldGetDimSize: dimension index > defined dims')
            return
         endif

         dimSize = size(field%data4DR8, dim=dimIndex)

      endif
   endif

!-----------------------------------------------------------------------
!
!  check for success
!
!-----------------------------------------------------------------------

   if (.not. foundData) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldGetDimSize: could not determine dimension size')
      return
   endif
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetDimSize

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldDetachData
! !INTERFACE:

   subroutine POP_FieldDetachData(field, errorCode)

! !DESCRIPTION:
!  Detaches a data array from the field type and associated metadata.
!  The field data itself is not deallocated.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field from which to detach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  disassociate all data pointers
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   nullify(field%data2DLog)
   nullify(field%data3DLog)
   nullify(field%data4DLog)
   nullify(field%data2DI4)
   nullify(field%data3DI4)
   nullify(field%data4DI4)
   nullify(field%data2DR4)
   nullify(field%data3DR4)
   nullify(field%data4DR4)
   nullify(field%data2DR8)
   nullify(field%data3DR8)
   nullify(field%data4DR8)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldDetachData

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData2DLog
! !INTERFACE:

   subroutine POP_FieldAttachData2DLog(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (POP_Logical), dimension(:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data2DLog => data

   call POP_FieldAttributeSet(field, 'numDims', 2, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeLogical, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData2DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData3DLog
! !INTERFACE:

   subroutine POP_FieldAttachData3DLog(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (POP_Logical), dimension(:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data3DLog => data

   call POP_FieldAttributeSet(field, 'numDims', 3, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeLogical, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData3DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData4DLog
! !INTERFACE:

   subroutine POP_FieldAttachData4DLog(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (POP_Logical), dimension(:,:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data4DLog => data

   call POP_FieldAttributeSet(field, 'numDims', 4, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeLogical, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData4DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData2DI4
! !INTERFACE:

   subroutine POP_FieldAttachData2DI4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data2DI4 => data

   call POP_FieldAttributeSet(field, 'numDims', 2, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeI4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData3DI4
! !INTERFACE:

   subroutine POP_FieldAttachData3DI4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data3DI4 => data

   call POP_FieldAttributeSet(field, 'numDims', 3, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeI4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData3DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData4DI4
! !INTERFACE:

   subroutine POP_FieldAttachData4DI4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data4DI4 => data

   call POP_FieldAttributeSet(field, 'numDims', 4, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeI4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData4DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData2DR4
! !INTERFACE:

   subroutine POP_FieldAttachData2DR4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data2DR4 => data

   call POP_FieldAttributeSet(field, 'numDims', 2, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData3DR4
! !INTERFACE:

   subroutine POP_FieldAttachData3DR4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data3DR4 => data

   call POP_FieldAttributeSet(field, 'numDims', 3, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData3DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData4DR4
! !INTERFACE:

   subroutine POP_FieldAttachData4DR4(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data4DR4 => data

   call POP_FieldAttributeSet(field, 'numDims', 4, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR4, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData4DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData2DR8
! !INTERFACE:

   subroutine POP_FieldAttachData2DR8(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data2DR8 => data

   call POP_FieldAttributeSet(field, 'numDims', 2, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR8, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData3DR8
! !INTERFACE:

   subroutine POP_FieldAttachData3DR8(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data3DR8 => data

   call POP_FieldAttributeSet(field, 'numDims', 3, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR8, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData3DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttachData4DR8
! !INTERFACE:

   subroutine POP_FieldAttachData4DR8(field, data, errorCode)

! !DESCRIPTION:
!  Attaches a data array to a field data type for use in routines
!  requiring field data.  This routine is a specific interface for
!  the generic POP\_FieldAttachData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:), target, intent(in) :: &
      data   ! data array to attach to field

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type to attach data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   field%data4DR8 => data

   call POP_FieldAttributeSet(field, 'numDims', 4, errorCode)
   call POP_FieldAttributeSet(field, 'dataType', &
                              POP_fieldDataTypeR8, errorCode)
   
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttachData4DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeSetChar
! !INTERFACE:

 subroutine POP_FieldAttributeSetChar(field, attName, attValue, &
                                      errorCode)

! !DESCRIPTION:
!  This routine sets an attribute in an existing field.  If the
!  attribute already exists, the value is reset.  If the attribute
!  does not exist, it is added to the field structure and the value
!  is set to the input value.  This is a specific interface for the
!  generic POP\_FieldAttributeSet interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName,       &! name  of attribute to be added
      attValue        ! value of attribute to be added

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field to which attribute is added

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      istat,           &! status flag for internal function calls
      numAttribs        ! current number of attributes defined

   character (POP_CharLength), dimension(:), allocatable :: &
      nameTmp,         &! temp space for resizing attrib name  array
      valTmp            ! temp space for resizing attrib value array

   logical (POP_Logical) :: &
      attExists         ! attribute already defined

!-----------------------------------------------------------------------
!
!  if this is the first defined attribute, allocate space and 
!  set the attribute name, value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (field%nAttsChar == 0) then

      !--- update attribute counters
      
      field%nAtts = field%nAtts + 1
      field%nAttsChar = 1

      allocate(field%attribValChar(1), &
               field%attribNameChar(1), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetChar: error allocating attr array')
         return
      endif

      field%attribNameChar (1) = attName
      field%attribValChar(1) = attValue

!-----------------------------------------------------------------------
!
!  if not the first user attribute, see if an attribute by this name
!   exists and overwrite value
!  if does not exist, resize the attribute array and store the 
!  attributes
!
!-----------------------------------------------------------------------

   else

      !--- search for attribute

      attExists = .false.
      numAttribs = field%nAttsChar

      attSearch: do n=1,numAttribs
         if (trim(field%attribNameChar(n)) == attName) then
            !--- reset value if attribute exists
            field%attribValChar(n) = trim(attValue)
            attExists = .true.
            exit attSearch
         endif
      end do attSearch

      if (.not. attExists) then

         !--- does not exist - resize attribute array to make room

         allocate(nameTmp(numAttribs), valTmp(numAttribs), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetChar: error allocating temp array')
            return
         endif

         nameTmp(:) = field%attribNameChar(:)
         valTmp (:) = field%attribValChar (:)

         deallocate(field%attribNameChar )
         deallocate(field%attribValChar)

         numAttribs = numAttribs + 1
         field%nAtts     = field%nAtts + 1
         field%nAttsChar = numAttribs

         allocate(field%attribNameChar (numAttribs), &
                  field%attribValChar(numAttribs), &
                  stat = istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetChar: error allocating new att array')
            return
         endif

         field%attribNameChar (1:numAttribs-1) = nameTmp
         field%attribValChar(1:numAttribs-1) =  valTmp
         field%attribNameChar (numAttribs) = attName
         field%attribValChar(numAttribs) = attValue

         deallocate(nameTmp, valTmp, stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetChar: error deallocating tmp array')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeSetChar

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeSetLog
! !INTERFACE:

 subroutine POP_FieldAttributeSetLog(field, attName, attValue, &
                                     errorCode)

! !DESCRIPTION:
!  This routine sets an attribute in an existing field.  If the
!  attribute already exists, the value is reset.  If the attribute
!  does not exist, it is added to the field structure and the value
!  is set to the input value.  This is a specific interface for the
!  generic POP\_FieldAttributeSet interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName         ! name  of attribute to be added

   logical (POP_Logical), intent(in) :: &
      attValue        ! value of attribute to be added

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field to which attribute is added

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      istat,           &! status flag for internal function calls
      numAttribs        ! current number of attributes defined

   character (POP_CharLength), dimension(:), allocatable :: &
      nameTmp           ! temp space for resizing attrib name  array

   logical (POP_Logical), dimension(:), allocatable :: &
      valTmp            ! temp space for resizing attrib value array

   logical (POP_Logical) :: &
      attExists         ! attribute already defined

!-----------------------------------------------------------------------
!
!  if this is the first user-defined attribute, allocate space and 
!  set the attribute name, value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (field%nAttsLog == 0) then

      !--- update attribute counters
      
      field%nAtts = field%nAtts + 1
      field%nAttsLog = 1

      !--- allocate and fill first attribute

      allocate(field%attribValLog(1), &
               field%attribNameLog(1), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetLog: error allocating attr array')
         return
      endif

      field%attribNameLog (1) = attName
      field%attribValLog(1) = attValue

!-----------------------------------------------------------------------
!
!  if not the first user attribute, see if an attribute by this name
!   exists and overwrite value
!  if does not exist, resize the attribute array and store the 
!  attributes
!
!-----------------------------------------------------------------------

   else

      !--- search for attribute

      attExists = .false.
      numAttribs = size(field%attribValLog(:))

      attSearch: do n=1,numAttribs
         if (trim(field%attribNameLog(n)) == attName) then
            !--- reset value if attribute exists
            field%attribValLog(n) = attValue
            attExists = .true.
            exit attSearch
         endif
      end do attSearch

      if (.not. attExists) then

         !--- does not exist - resize attribute array to make room

         allocate(nameTmp(numAttribs), valTmp(numAttribs), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetLog: error allocating temp array')
            return
         endif

         nameTmp(:) = field%attribNameLog (:)
         valTmp (:) = field%attribValLog(:)

         deallocate(field%attribNameLog )
         deallocate(field%attribValLog)

         numAttribs         = numAttribs + 1
         field%nAtts    = field%nAtts + 1
         field%nAttsLog = numAttribs

         allocate(field%attribNameLog (numAttribs), &
                  field%attribValLog(numAttribs), &
                  stat = istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetLog: error allocating new att array')
            return
         endif

         field%attribNameLog (1:numAttribs-1) = nameTmp
         field%attribValLog(1:numAttribs-1) =  valTmp
         field%attribNameLog (numAttribs) = attName
         field%attribValLog(numAttribs) = attValue

         deallocate(nameTmp,valTmp, stat=istat)
         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetLog: error deallocating tmp array')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeSetLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeSetI4
! !INTERFACE:

 subroutine POP_FieldAttributeSetI4(field, attName, attValue, &
                                      errorCode)

! !DESCRIPTION:
!  This routine sets an attribute in an existing field.  If the
!  attribute already exists, the value is reset.  If the attribute
!  does not exist, it is added to the field structure and the value
!  is set to the input value.  This is a specific interface for the
!  generic POP\_FieldAttributeSet interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName         ! name  of attribute to be added

   integer (POP_i4), intent(in) :: &
      attValue        ! value of attribute to be added

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field to which attribute is added

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      istat,           &! status flag for internal function calls
      numAttribs        ! current number of attributes defined

   character (POP_CharLength), dimension(:), allocatable :: &
      nameTmp           ! temp space for resizing attrib name  array

   integer (POP_i4), dimension(:), allocatable :: &
      valTmp            ! temp space for resizing attrib value array

   logical (POP_Logical) :: &
      attExists         ! attribute already defined

!-----------------------------------------------------------------------
!
!  if this is the first user-defined attribute, allocate space and 
!  set the attribute name, value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (field%nAttsI4 == 0) then

      !--- update attribute counters
      
      field%nAtts = field%nAtts + 1
      field%nAttsI4  = 1

      !--- allocate and fill first attribute

      allocate(field%attribValI4(1), &
               field%attribNameI4(1), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetI4: error allocating attr array')
         return
      endif

      field%attribNameI4 (1) = attName
      field%attribValI4(1) = attValue

!-----------------------------------------------------------------------
!
!  if not the first user attribute, see if an attribute by this name
!   exists and overwrite value
!  if does not exist, resize the attribute array and store the 
!  attributes
!
!-----------------------------------------------------------------------

   else

      !--- search for attribute

      attExists = .false.
      numAttribs = size(field%attribValI4(:))

      attSearch: do n=1,numAttribs
         if (trim(field%attribNameI4(n)) == attName) then
            !--- reset value if attribute exists
            field%attribValI4(n) = attValue
            attExists = .true.
            exit attSearch
         endif
      end do attSearch

      if (.not. attExists) then

         !--- does not exist - resize attribute array to make room

         allocate(nameTmp(numAttribs), valTmp(numAttribs), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetI4: error allocating temp array')
            return
         endif

         nameTmp(:) = field%attribNameI4 (:)
         valTmp (:) = field%attribValI4(:)

         deallocate(field%attribNameI4 )
         deallocate(field%attribValI4)

         numAttribs        = numAttribs + 1
         field%nAtts   = field%nAtts + 1
         field%nAttsI4 = numAttribs

         allocate(field%attribNameI4 (numAttribs), &
                  field%attribValI4(numAttribs), &
                  stat = istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetI4: error allocating new att array')
            return
         endif

         field%attribNameI4 (1:numAttribs-1) = nameTmp
         field%attribValI4(1:numAttribs-1) =  valTmp
         field%attribNameI4 (numAttribs) = attName
         field%attribValI4(numAttribs) = attValue

         deallocate(nameTmp,valTmp, stat = istat)
         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetI4: error deallocating tmp array')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeSetI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeSetR4
! !INTERFACE:

 subroutine POP_FieldAttributeSetR4(field, attName, attValue, &
                                    errorCode)

! !DESCRIPTION:
!  This routine sets an attribute in an existing field.  If the
!  attribute already exists, the value is reset.  If the attribute
!  does not exist, it is added to the field structure and the value
!  is set to the input value.  This is a specific interface for the
!  generic POP\_FieldAttributeSet interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName         ! name  of attribute to be added

   real (POP_r4), intent(in) :: &
      attValue        ! value of attribute to be added

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field to which attribute is added

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      istat,           &! status flag for internal function calls
      numAttribs        ! current number of attributes defined

   character (POP_CharLength), dimension(:), allocatable :: &
      nameTmp           ! temp space for resizing attrib name  array

   real (POP_r4), dimension(:), allocatable :: &
      valTmp            ! temp space for resizing attrib value array

   logical (POP_Logical) :: &
      attExists         ! attribute already defined

!-----------------------------------------------------------------------
!
!  if this is the first user-defined attribute, allocate space and 
!  set the attribute name, value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (field%nAttsR4 == 0) then

      !--- update attribute counters
      
      field%nAtts   = field%nAtts + 1
      field%nAttsR4 = 1

      !--- allocate and fill first attribute

      allocate(field%attribValR4(1), &
               field%attribNameR4(1), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetR4: error allocating attr array')
         return
      endif

      field%attribNameR4 (1) = attName
      field%attribValR4(1) = attValue

!-----------------------------------------------------------------------
!
!  if not the first user attribute, see if an attribute by this name
!   exists and overwrite value
!  if does not exist, resize the attribute array and store the 
!  attributes
!
!-----------------------------------------------------------------------

   else

      !--- search for attribute

      attExists = .false.
      numAttribs = size(field%attribValR4(:))

      attSearch: do n=1,numAttribs
         if (trim(field%attribNameR4(n)) == attName) then
            !--- reset value if attribute exists
            field%attribValR4(n) = attValue
            attExists = .true.
            exit attSearch
         endif
      end do attSearch

      if (.not. attExists) then

         !--- does not exist - resize attribute array to make room

         allocate(nameTmp(numAttribs), valTmp(numAttribs), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetR4: error allocating temp array')
            return
         endif

         nameTmp(:) = field%attribNameR4(:)
         valTmp (:) = field%attribValR4 (:)

         deallocate(field%attribNameR4 )
         deallocate(field%attribValR4)

         numAttribs        = numAttribs + 1
         field%nAtts   = field%nAtts + 1
         field%nAttsR4 = numAttribs

         allocate(field%attribNameR4 (numAttribs), &
                  field%attribValR4(numAttribs), &
                  stat = istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetR4: error allocating new att array')
            return
         endif

         field%attribNameR4 (1:numAttribs-1) = nameTmp
         field%attribValR4(1:numAttribs-1) =  valTmp
         field%attribNameR4 (numAttribs) = attName
         field%attribValR4(numAttribs) = attValue

         deallocate(nameTmp,valTmp, stat = istat)
         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetR4: error deallocating tmp array')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeSetR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeSetR8
! !INTERFACE:

 subroutine POP_FieldAttributeSetR8(field, attName, attValue, &
                                    errorCode)

! !DESCRIPTION:
!  This routine sets an attribute in an existing field.  If the
!  attribute already exists, the value is reset.  If the attribute
!  does not exist, it is added to the field structure and the value
!  is set to the input value.  This is a specific interface for the
!  generic POP\_FieldAttributeSet interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName         ! name  of attribute to be added

   real (POP_r8), intent(in) :: &
      attValue        ! value of attribute to be added

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field to which attribute is added

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      istat,           &! status flag for internal function calls
      numAttribs        ! current number of attributes defined

   character (POP_CharLength), dimension(:), allocatable :: &
      nameTmp           ! temp space for resizing attrib name  array

   real (POP_r8), dimension(:), allocatable :: &
      valTmp            ! temp space for resizing attrib value array

   logical (POP_Logical) :: &
      attExists         ! attribute already defined

!-----------------------------------------------------------------------
!
!  if this is the first user-defined attribute, allocate space and 
!  set the attribute name, value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (field%nAttsR8 == 0) then

      !--- update attribute counters
      
      field%nAtts   = field%nAtts + 1
      field%nAttsR8 = 1

      !--- allocate and fill first attribute

      allocate(field%attribValR8(1), &
               field%attribNameR8(1), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetR8: error allocating attr array')
         return
      endif

      field%attribNameR8 (1) = attName
      field%attribValR8(1) = attValue

!-----------------------------------------------------------------------
!
!  if not the first user attribute, see if an attribute by this name
!   exists and overwrite value
!  if does not exist, resize the attribute array and store the 
!  attributes
!
!-----------------------------------------------------------------------

   else

      !--- search for attribute

      attExists = .false.
      numAttribs = size(field%attribValR8(:))

      attSearch: do n=1,numAttribs
         if (trim(field%attribNameR8(n)) == attName) then
            !--- reset value if attribute exists
            field%attribValR8(n) = attValue
            attExists = .true.
            exit attSearch
         endif
      end do attSearch

      if (.not. attExists) then

         !--- does not exist - resize attribute array to make room

         allocate(nameTmp(numAttribs), valTmp(numAttribs), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
              'POP_FieldAttributeSetR8: error allocating temp array')
            return
         endif

         nameTmp(:) = field%attribNameR8(:)
         valTmp (:) = field%attribValR8 (:)

         deallocate(field%attribNameR8 )
         deallocate(field%attribValR8)

         numAttribs        = numAttribs + 1
         field%nAtts   = field%nAtts + 1
         field%nAttsR8 = numAttribs

         allocate(field%attribNameR8 (numAttribs), &
                  field%attribValR8(numAttribs), &
                  stat = istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetR8: error allocating new att array')
            return
         endif

         field%attribNameR8 (1:numAttribs-1) = nameTmp
         field%attribValR8(1:numAttribs-1) =  valTmp
         field%attribNameR8 (numAttribs) = attName
         field%attribValR8(numAttribs) = attValue

         deallocate(nameTmp,valTmp, stat=istat)
         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeSetR8: error deallocating tmp array')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeSetR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetChar
! !INTERFACE:

 subroutine POP_FieldAttributeGetChar(field, attValue, errorCode, &
                                      attName, attIndex)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  It can 
!  retrieve the attribute by name or by index number.  The latter
!  case is useful for querying a field for all available attributes
!  and both the name (if requested) and value are returned. 
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field           ! field from which attribute is retrieved

   integer (POP_i4), intent(in), optional :: &
      attIndex        ! location of attribute in att array

! !INPUT/OUTPUT PARAMETERS:

   character (POP_CharLength), intent(inout), optional :: &
      attName         ! on input: name of attribute to be retrieved
                      ! if attIndex is supplied the attName will be
                      !    returned as an output

! !OUTPUT PARAMETERS:

   character (POP_CharLength), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      numAttribs        ! current number of attributes defined

   logical (POP_Logical) :: &
      attFound         ! attribute found

!-----------------------------------------------------------------------
!
!  set defaults
!
!-----------------------------------------------------------------------

   errorCode  = POP_Success
   attFound   = .false.
   numAttribs = field%nAttsChar

!-----------------------------------------------------------------------
!
!  if attribute is requested by name (index not supplied),
!  search for attribute name.  first check standard attributes.
!
!-----------------------------------------------------------------------

   if (.not. present(attIndex)) then

      if (.not. present(attName)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetChar: name or index must be supplied')
         return
      endif
      
      if (.not. attFound) then

         !--- search for attribute

         attSearchName: do n=1,numAttribs
            if (trim(field%attribNameChar(n)) == trim(attName)) then
               attFound = .true.
               attValue = trim(field%attribValChar(n))
               exit attSearchName
            endif
         end do attSearchName

      endif

!-----------------------------------------------------------------------
!
!  if attIndex supplied, return both name (if requested) and value
!
!-----------------------------------------------------------------------

   else  ! attIndex present

      !--- check for bad attIndex
      
      if (attIndex < 1 .or. attIndex > numAttribs) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetChar: attIndex out of range')
         return
      endif

      !--- grab the attribute name and value
      
      attFound = .true.
      attValue = trim(field%attribValChar(attIndex))
      if (present(attName)) &
         attName = trim(field%attribNameChar(attIndex))

   endif
  
!-----------------------------------------------------------------------
!
!  return error if attribute not found
!
!-----------------------------------------------------------------------

   if (.not. attFound) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldAttributeGetChar: attribute not found')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetChar

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetLog
! !INTERFACE:

 subroutine POP_FieldAttributeGetLog(field, attValue, errorCode, &
                                            attName, attIndex)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  It can 
!  retrieve the attribute by name or by index number.  The latter
!  case is useful for querying a field for all available attributes
!  and both the name (if requested) and value are returned. 
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field           ! field from which attribute is retrieved

   integer (POP_i4), intent(in), optional :: &
      attIndex        ! location of attribute in att array

! !INPUT/OUTPUT PARAMETERS:

   character (POP_CharLength), intent(inout), optional :: &
      attName         ! on input: name of attribute to be retrieved
                      ! if attIndex is supplied the attName will be
                      !    returned as an output

! !OUTPUT PARAMETERS:

   logical (POP_Logical), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      numAttribs        ! current number of attributes defined

   logical (POP_Logical) :: &
      attFound         ! attribute found

!-----------------------------------------------------------------------
!
!  set defaults
!
!-----------------------------------------------------------------------

   errorCode  = POP_Success
   attFound   = .false.
   numAttribs = field%nAttsLog

!-----------------------------------------------------------------------
!
!  if attribute is requested by name (index not supplied),
!  search for attribute name.  first check standard attributes.
!
!-----------------------------------------------------------------------

   if (.not. present(attIndex)) then

      if (.not. present(attName)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetLog: name or index must be supplied')
         return
      endif
      
      if (.not. attFound) then

         !--- search for attribute

         attSearchName: do n=1,numAttribs
            if (trim(field%attribNameLog(n)) == trim(attName)) then
               attFound = .true.
               attValue = field%attribValLog(n)
               exit attSearchName
            endif
         end do attSearchName

      endif

!-----------------------------------------------------------------------
!
!  if attIndex supplied, return both name (if requested) and value
!
!-----------------------------------------------------------------------

   else  ! attIndex present

      !--- check for bad attIndex
      
     if (attIndex < 1 .or. attIndex > numAttribs) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetLog: attIndex out of range')
         return
      endif

      !--- grab the attribute name and value
      
      attFound = .true.
      attValue = field%attribValLog(attIndex)
      if (present(attName)) &
         attName = trim(field%attribNameLog(attIndex))

   endif
  
!-----------------------------------------------------------------------
!
!  return error if attribute not found
!
!-----------------------------------------------------------------------

   if (.not. attFound) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldAttributeGetLog: attribute not found')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetI4
! !INTERFACE:

 subroutine POP_FieldAttributeGetI4(field, attValue, errorCode, &
                                           attName, attIndex)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  It can 
!  retrieve the attribute by name or by index number.  The latter
!  case is useful for querying a field for all available attributes
!  and both the name (if requested) and value are returned. 
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field           ! field from which attribute is retrieved

   integer (POP_i4), intent(in), optional :: &
      attIndex        ! location of attribute in att array

! !INPUT/OUTPUT PARAMETERS:

   character (POP_CharLength), intent(inout), optional :: &
      attName         ! on input: name of attribute to be retrieved
                      ! if attIndex is supplied the attName will be
                      !    returned as an output

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      numAttribs        ! current number of attributes defined

   logical (POP_Logical) :: &
      attFound         ! attribute found

!-----------------------------------------------------------------------
!
!  set defaults
!
!-----------------------------------------------------------------------

   errorCode  = POP_Success
   attFound   = .false.
   numAttribs = field%nAttsI4

!-----------------------------------------------------------------------
!
!  if attribute is requested by name (index not supplied),
!  search for attribute name.  first check standard attributes.
!
!-----------------------------------------------------------------------

   if (.not. present(attIndex)) then

      if (.not. present(attName)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetI4: name or index must be supplied')
         return
      endif
      
      if (.not. attFound) then

         !--- search for attribute

         attSearchName: do n=1,numAttribs
            if (trim(field%attribNameI4(n)) == trim(attName)) then
               attFound = .true.
               attValue = field%attribValI4(n)
               exit attSearchName
            endif
         end do attSearchName

      endif

!-----------------------------------------------------------------------
!
!  if attIndex supplied, return both name (if requested) and value
!
!-----------------------------------------------------------------------

   else  ! attIndex present

      !--- check for bad attIndex
      
      if (attIndex < 1 .or. attIndex > numAttribs) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetI4: attIndex out of range')
         return
      endif

      !--- grab the attribute name and value
      
      attFound = .true.
      attValue = field%attribValI4(attIndex)
      if (present(attName)) &
         attName = trim(field%attribNameI4(attIndex))

   endif

!-----------------------------------------------------------------------
!
!  return error if attribute not found
!
!-----------------------------------------------------------------------

   if (.not. attFound) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldAttributeGetI4: attribute not found')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetR4
! !INTERFACE:

 subroutine POP_FieldAttributeGetR4(field, attValue, errorCode, &
                                           attName, attIndex)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  It can 
!  retrieve the attribute by name or by index number.  The latter
!  case is useful for querying a field for all available attributes
!  and both the name (if requested) and value are returned. 
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field           ! field from which attribute is retrieved

   integer (POP_i4), intent(in), optional :: &
      attIndex        ! location of attribute in att array

! !INPUT/OUTPUT PARAMETERS:

   character (POP_CharLength), intent(inout), optional :: &
      attName         ! on input: name of attribute to be retrieved
                      ! if attIndex is supplied the attName will be
                      !    returned as an output

! !OUTPUT PARAMETERS:

   real (POP_r4), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      numAttribs        ! current number of attributes defined

   logical (POP_Logical) :: &
      attFound         ! attribute found

!-----------------------------------------------------------------------
!
!  set defaults
!
!-----------------------------------------------------------------------

   errorCode  = POP_Success
   attFound   = .false.
   numAttribs = field%nAttsR4

!-----------------------------------------------------------------------
!
!  if attribute is requested by name (index not supplied),
!  search for attribute name.  first check standard attributes.
!
!-----------------------------------------------------------------------

   if (.not. present(attIndex)) then

      if (.not. present(attName)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetR4: name or index must be supplied')
         return
      endif
      
      if (.not. attFound) then

         !--- search for attribute

         attSearchName: do n=1,numAttribs
            if (trim(field%attribNameR4(n)) == trim(attName)) then
               attFound = .true.
               attValue = field%attribValR4(n)
               exit attSearchName
            endif
         end do attSearchName

      endif

!-----------------------------------------------------------------------
!
!  if attIndex supplied, return both name (if requested) and value
!
!-----------------------------------------------------------------------

   else  ! attIndex present

      !--- check for bad attIndex
      
      if (attIndex < 1 .or. attIndex > numAttribs) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetR4: attIndex out of range')
         return
      endif

      !--- grab the attribute name and value
      
      attFound = .true.
      attValue = field%attribValR4(attIndex)
      if (present(attName)) &
         attName = trim(field%attribNameR4(attIndex))

   endif
  
!-----------------------------------------------------------------------
!
!  return error if attribute not found
!
!-----------------------------------------------------------------------

   if (.not. attFound) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldAttributeGetR4: attribute not found')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetR8
! !INTERFACE:

 subroutine POP_FieldAttributeGetR8(field, attValue, errorCode, &
                                           attName, attIndex)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  It can 
!  retrieve the attribute by name or by index number.  The latter
!  case is useful for querying a field for all available attributes
!  and both the name (if requested) and value are returned. 
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(in) :: &
      field           ! field from which attribute is retrieved

   integer (POP_i4), intent(in), optional :: &
      attIndex        ! location of attribute in att array

! !INPUT/OUTPUT PARAMETERS:

   character (POP_CharLength), intent(inout), optional :: &
      attName         ! on input: name of attribute to be retrieved
                      ! if attIndex is supplied the attName will be
                      !    returned as an output

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! loop index
      numAttribs        ! current number of attributes defined

   logical (POP_Logical) :: &
      attFound         ! attribute found

!-----------------------------------------------------------------------
!
!  set defaults
!
!-----------------------------------------------------------------------

   errorCode  = POP_Success
   attFound   = .false.
   numAttribs = field%nAttsR8

!-----------------------------------------------------------------------
!
!  if attribute is requested by name (index not supplied),
!  search for attribute name.  first check standard attributes.
!
!-----------------------------------------------------------------------

   if (.not. present(attIndex)) then

      if (.not. present(attName)) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetR8: name or index must be supplied')
         return
      endif
      
      if (.not. attFound) then

         !--- search for attribute

         attSearchName: do n=1,numAttribs
            if (trim(field%attribNameR8(n)) == trim(attName)) then
               attFound = .true.
               attValue = field%attribValR8(n)
               exit attSearchName
            endif
         end do attSearchName

      endif

!-----------------------------------------------------------------------
!
!  if attIndex supplied, return both name (if requested) and value
!
!-----------------------------------------------------------------------

   else  ! attIndex present

      !--- check for bad attIndex
      
      if (attIndex < 1 .or. attIndex > numAttribs) then
         call POP_ErrorSet(errorCode, &
            'POP_FieldAttributeGetR8: attIndex out of range')
         return
      endif

      !--- grab the attribute name and value
      
      attFound = .true.
      attValue = field%attribValR8(attIndex)
      if (present(attName)) &
         attName = trim(field%attribNameR8(attIndex))

   endif
  
!-----------------------------------------------------------------------
!
!  return error if attribute not found
!
!-----------------------------------------------------------------------

   if (.not. attFound) then
      call POP_ErrorSet(errorCode, &
         'POP_FieldAttributeGetR8: attribute not found')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldAttributeGetDims
! !INTERFACE:

 subroutine POP_FieldAttributeGetDims(field, attName, attValue, &
                                      errorCode)

! !DESCRIPTION:
!  This routine gets an attribute from an existing field.  If the 
!  attribute name does not exist, it returns an error.
!  This is a specific interface for the generic POP\_FieldAttributeGet 
!  interface corresponding to a grid dimension array.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      attName         ! name  of attribute to be retrieved

! !INPUT/OUTPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field           ! field from which attribute is retrieved

! !OUTPUT PARAMETERS:

   type (POP_GridDim), dimension(:), intent(out) :: &
      attValue        ! value of attribute to be retrieved

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if this is one of the required attributes, reset the value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success


   if (associated(field%fieldDims)) then
      attValue = field%fieldDims

   else

      call POP_ErrorSet(errorCode, &
                        'POP_FieldAttributeGetDims: no dims defined')

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldAttributeGetDims

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData2DLog
! !INTERFACE:

   subroutine POP_FieldGetData2DLog(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   logical (POP_Logical), dimension(:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data2DLog)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData2DLog: data pointer not associated')
      return
   endif

   data => field%data2DLog

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData2DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData3DLog
! !INTERFACE:

   subroutine POP_FieldGetData3DLog(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   logical (POP_Logical), dimension(:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data3DLog)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData3DLog: data pointer not associated')
      return
   endif

   data => field%data3DLog

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData3DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData4DLog
! !INTERFACE:

   subroutine POP_FieldGetData4DLog(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   logical (POP_Logical), dimension(:,:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data4DLog)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData4DLog: data pointer not associated')
      return
   endif

   data => field%data4DLog

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData4DLog

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData2DI4
! !INTERFACE:

   subroutine POP_FieldGetData2DI4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   integer (POP_i4), dimension(:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data2DI4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData2DI4: data pointer not associated')
      return
   endif

   data => field%data2DI4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData3DI4
! !INTERFACE:

   subroutine POP_FieldGetData3DI4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   integer (POP_i4), dimension(:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data3DI4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData3DI4: data pointer not associated')
      return
   endif

   data => field%data3DI4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData3DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData4DI4
! !INTERFACE:

   subroutine POP_FieldGetData4DI4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   integer (POP_i4), dimension(:,:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data4DI4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData4DI4: data pointer not associated')
      return
   endif

   data => field%data4DI4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData4DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData2DR4
! !INTERFACE:

   subroutine POP_FieldGetData2DR4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r4), dimension(:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data2DR4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData2DR4: data pointer not associated')
      return
   endif

   data => field%data2DR4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData3DR4
! !INTERFACE:

   subroutine POP_FieldGetData3DR4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r4), dimension(:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data3DR4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData3DR4: data pointer not associated')
      return
   endif

   data => field%data3DR4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData3DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData4DR4
! !INTERFACE:

   subroutine POP_FieldGetData4DR4(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r4), dimension(:,:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data4DR4)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData4DR4: data pointer not associated')
      return
   endif

   data => field%data4DR4

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData4DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData2DR8
! !INTERFACE:

   subroutine POP_FieldGetData2DR8(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r8), dimension(:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data2DR8)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData2DR8: data pointer not associated')
      return
   endif

   data => field%data2DR8

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData3DR8
! !INTERFACE:

   subroutine POP_FieldGetData3DR8(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r8), dimension(:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data3DR8)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData3DR8: data pointer not associated')
      return
   endif

   data => field%data3DR8

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData3DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_FieldGetData4DR8
! !INTERFACE:

   subroutine POP_FieldGetData4DR8(field, data, errorCode)

! !DESCRIPTION:
!  Returns a pointer to the data attached to a field.  Only the
!  pointer is returned and not a duplicate copy. The user is
!  responsible for making any local copies of the data. 
!  This routine is a specific interface for
!  the generic POP\_FieldGetData interface.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_Field), intent(inout) :: &
      field  ! field type from which to retrieve data

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode ! returned error code

   real (POP_r8), dimension(:,:,:,:,:), pointer :: &
      data   ! pointer to data array attached to this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  set pointer to data array
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. associated(field%data4DR8)) then
      call POP_ErrorSet(errorCode, &
          'POP_FieldGetData4DR8: data pointer not associated')
      return
   endif

   data => field%data4DR8

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_FieldGetData4DR8

!***********************************************************************

 end module POP_FieldMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

