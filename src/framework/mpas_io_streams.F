! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io_streams

#define COMMA ,
#define STREAM_DEBUG_WRITE(M) write(stderrUnit,*) M

#ifdef SINGLE_PRECISION
#define REAL_IO_TYPE MPAS_IO_REAL
#else
#define REAL_IO_TYPE MPAS_IO_DOUBLE
#endif

   use mpas_attlist
   use mpas_grid_types
   use mpas_timekeeping
   use mpas_io
   use mpas_io_units

   type field_list_type
      integer :: field_type = -999
      logical :: isDecomposed
      integer :: totalDimSize     ! Total size of outer dimension across all blocks for decomposed fields
      logical, dimension(:), pointer :: isAvailable => null() ! Used for reading var-arrays where one or more 
                                                              !   constitutent arrays may not be present in the input file
      type (field0dInteger), pointer :: int0dField => null()
      type (field1dInteger), pointer :: int1dField => null()
      type (field2dInteger), pointer :: int2dField => null()
      type (field3dInteger), pointer :: int3dField => null()
      type (field0dReal), pointer :: real0dField => null()
      type (field1dReal), pointer :: real1dField => null()
      type (field2dReal), pointer :: real2dField => null()
      type (field3dReal), pointer :: real3dField => null()
      type (field4dReal), pointer :: real4dField => null()
      type (field5dReal), pointer :: real5dField => null()
      type (field0dChar), pointer :: char0dField => null()
      type (field1dChar), pointer :: char1dField => null()
      type (field_list_type), pointer :: next => null()
   end type field_list_type

   type MPAS_Stream_type
      logical :: isInitialized = .false.
      integer :: ioFormat
      integer :: ioDirection
      integer :: framesPerFile
      integer :: defaultPrecision = MPAS_IO_NATIVE_PRECISION
      logical :: clobberRecords = .false.
      type (MPAS_IO_Handle_type) :: fileHandle
      type (att_list_type), pointer :: attList => null()
      type (field_list_type), pointer :: fieldList => null()
   end type MPAS_Stream_type


   interface MPAS_streamAddField
      module procedure MPAS_streamAddField_0dInteger
      module procedure MPAS_streamAddField_1dInteger
      module procedure MPAS_streamAddField_2dInteger
      module procedure MPAS_streamAddField_3dInteger
      module procedure MPAS_streamAddField_0dReal
      module procedure MPAS_streamAddField_1dReal
      module procedure MPAS_streamAddField_2dReal
      module procedure MPAS_streamAddField_3dReal
      module procedure MPAS_streamAddField_4dReal
      module procedure MPAS_streamAddField_5dReal
      module procedure MPAS_streamAddField_0dChar
   end interface MPAS_streamAddField

   interface MPAS_streamUpdateField
      module procedure MPAS_streamUpdateField_0dInteger
      module procedure MPAS_streamUpdateField_1dInteger
      module procedure MPAS_streamUpdateField_2dInteger
      module procedure MPAS_streamUpdateField_3dInteger
      module procedure MPAS_streamUpdateField_0dReal
      module procedure MPAS_streamUpdateField_1dReal
      module procedure MPAS_streamUpdateField_2dReal
      module procedure MPAS_streamUpdateField_3dReal
      module procedure MPAS_streamUpdateField_4dReal
      module procedure MPAS_streamUpdateField_5dReal
      module procedure MPAS_streamUpdateField_0dChar
   end interface MPAS_streamUpdateField

   interface MPAS_readStreamAtt
      module procedure MPAS_readStreamAtt_0dInteger
      module procedure MPAS_readStreamAtt_1dInteger
      module procedure MPAS_readStreamAtt_0dReal
      module procedure MPAS_readStreamAtt_1dReal
      module procedure MPAS_readStreamAtt_text
   end interface MPAS_readStreamAtt

   interface MPAS_writeStreamAtt
      module procedure MPAS_writeStreamAtt_0dInteger
      module procedure MPAS_writeStreamAtt_1dInteger
      module procedure MPAS_writeStreamAtt_0dReal
      module procedure MPAS_writeStreamAtt_1dReal
      module procedure MPAS_writeStreamAtt_text
   end interface MPAS_writeStreamAtt

   integer, parameter :: MPAS_STREAM_EXACT_TIME              = 100, &
                         MPAS_STREAM_NEAREST                 = 101, &
                         MPAS_STREAM_LATEST_BEFORE           = 102, &
                         MPAS_STREAM_EARLIEST_AFTER          = 103, &
                         MPAS_STREAM_LATEST_STRICTLY_BEFORE  = 104, &
                         MPAS_STREAM_EARLIEST_STRICTLY_AFTER = 105


   ! Error codes
   integer, parameter :: MPAS_STREAM_NOERR              =  0, &
                         MPAS_STREAM_NOT_INITIALIZED    = -1, &
                         MPAS_STREAM_FIELD_NOT_FOUND    = -2, &
                         MPAS_STREAM_CLOBBER_FILE       = -3, &
                         MPAS_STREAM_CLOBBER_RECORD     = -4, &
                         MPAS_IO_ERR                    = -5

   integer, parameter :: FIELD_0D_INT   =  1, &
                         FIELD_1D_INT   =  2, &
                         FIELD_2D_INT   =  3, &
                         FIELD_3D_INT   =  4, &
                         FIELD_0D_REAL  =  5, &
                         FIELD_1D_REAL  =  6, &
                         FIELD_2D_REAL  =  7, &
                         FIELD_3D_REAL  =  8, &
                         FIELD_4D_REAL  =  9, &
                         FIELD_5D_REAL  =  10, &
                         FIELD_0D_CHAR  =  11, &
                         FIELD_1D_CHAR  =  12

   private mergeArrays


   contains


   subroutine MPAS_createStream(stream, fileName, ioFormat, ioDirection, framesPerFile, precision, &
                                clobberRecords, clobberFiles, truncateFiles, ierr)

      implicit none

      type (MPAS_Stream_type), intent(out) :: stream
      character (len=*), intent(in) :: fileName
      integer, intent(in) :: ioFormat
      integer, intent(in) :: ioDirection
      integer, intent(in) :: framesPerFile
      integer, intent(in), optional :: precision
      logical, intent(in), optional :: clobberRecords
      logical, intent(in), optional :: clobberFiles
      logical, intent(in), optional :: truncateFiles
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR


      stream % fileHandle = MPAS_io_open(fileName, ioDirection, ioFormat, clobber_file=clobberFiles, truncate_file=truncateFiles, &
                                         ierr=io_err) 
      !
      ! Catch a few special errors
      !
      if (io_err == MPAS_IO_ERR_WOULD_CLOBBER) then
         if (present(ierr)) ierr = MPAS_STREAM_CLOBBER_FILE
         return
      end if

      ! General error
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      stream % ioDirection = ioDirection
      stream % ioFormat = ioFormat
      stream % framesPerFile = framesPerFile
      if (present(clobberRecords)) then
         stream % clobberRecords = clobberRecords
      else
         stream % clobberRecords = .false.
      end if
      if (present(precision)) stream % defaultPrecision = precision

      stream % isInitialized = .true.

   end subroutine MPAS_createStream


   integer function MPAS_seekStream(stream, seekTime, seekPosition, actualTime, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: seekTime
      integer, intent(in) :: seekPosition
      character (len=*), intent(out) :: actualTime
      integer, intent(out), optional :: ierr
      
      integer :: io_err
      integer :: i
      integer :: timeDim
      character (len=StrKIND), dimension(:), pointer :: xtimes
      character (len=StrKIND) :: strTemp
      type (MPAS_Time_type) :: sliceTime, startTime
      type (MPAS_TimeInterval_type) :: timeDiff, minTimeDiff

!      write(stderrUnit,*) 'Called MPAS_seekStream'

      !
      ! Initialize output arguments
      !
      if (present(ierr)) ierr = MPAS_STREAM_NOERR
      MPAS_seekStream = 0


      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      call MPAS_io_inq_dim(stream % fileHandle, 'Time', timeDim, io_err)
      if (timeDim <= 0 .or. io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

!write(stderrUnit,*) 'Found ', timeDim, ' times in file' 

      call MPAS_io_inq_var(stream % fileHandle, 'xtime', ierr=io_err)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

!write(stderrUnit,*) 'Found xtime variable'

      allocate(xtimes(timeDim))

      do i=1,timeDim
         call MPAS_io_set_frame(stream % fileHandle, i, io_err)
         call MPAS_io_get_var(stream % fileHandle, 'xtime', xtimes(i), io_err)
!write(stderrUnit,*) '... just read in xtime='//xtimes(i)
      end do

      if (len(seekTime) > 32) then
         write(strTemp, '(a)') seekTime(1:32)
      else
         write(strTemp, '(a)') trim(seekTime)
      end if
      call mpas_set_timeInterval(interval=minTimeDiff, DD=10000)
      call mpas_set_time(curr_time=startTime, dateTimeString=strTemp)

      do i=1,timeDim
         write(strTemp, '(a)') trim(xtimes(i)(1:32))
!write(stderrUnit,*) '... converted strTemp='//strTemp
         call mpas_set_time(curr_time=sliceTime, dateTimeString=strTemp)
         if (seekPosition == MPAS_STREAM_EXACT_TIME) then
            if (sliceTime == startTime) then
               minTimeDiff = timeDiff
               MPAS_seekStream = i
            end if
         else if (seekPosition == MPAS_STREAM_NEAREST) then
            timeDiff = abs(sliceTime - startTime)
            if (timeDiff < minTimeDiff) then
               minTimeDiff = timeDiff
               MPAS_seekStream = i
            end if
         else if (seekPosition == MPAS_STREAM_LATEST_BEFORE) then
            if (sliceTime <= startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  MPAS_seekStream = i
               end if
            end if
         else if (seekPosition == MPAS_STREAM_EARLIEST_AFTER) then
            if (sliceTime >= startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  MPAS_seekStream = i
               end if
            end if
         else
            write(stderrUnit,*) 'Error in MPAS_seekStream: unrecognized seekPosition'
            deallocate(xtimes)
            if (present(ierr)) ierr = MPAS_IO_ERR
            return
         end if
      end do

      if (MPAS_seekStream == 0) then
         deallocate(xtimes)
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      write(actualTime, '(a)') trim(xtimes(MPAS_seekStream)(1:32))

      deallocate(xtimes)

   end function MPAS_seekStream


   subroutine MPAS_streamTime(stream, frame, frameValidTime, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      character (len=*), intent(out) :: frameValidTime
      integer, intent(out), optional :: ierr
      
      integer :: io_err
      integer :: timeDim


!      write(stderrUnit,*) 'Called MPAS_streamTime'

      !
      ! Initialize output arguments
      !
      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      call MPAS_io_inq_var(stream % fileHandle, 'xtime', ierr=io_err)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

!write(stderrUnit,*) 'Found xtime variable'

      call MPAS_io_set_frame(stream % fileHandle, frame, io_err)
      call MPAS_io_get_var(stream % fileHandle, 'xtime', frameValidTime, io_err)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

   end subroutine MPAS_streamTime


   subroutine MPAS_streamAddField_0dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field0dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = 0

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      isDecomposed = .false.
      globalDimSize = 0
      totalDimSize = 0

      call MPAS_streamAddField_generic(stream, trim(field % fieldName), MPAS_IO_INT, dimNames, dimSizes, &
                                       field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= MPAS_STREAM_NOERR) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_INT
      new_field_list_node % int0dField => field

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_0dInteger


   subroutine MPAS_streamAddField_1dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field1DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field1dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), MPAS_IO_INT, dimNames0, &
                                             dimSizes0, field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), MPAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_INT
      new_field_list_node % int1dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_1dInteger


   subroutine MPAS_streamAddField_2dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field2DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field2dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), MPAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), MPAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_INT
      new_field_list_node % int2dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_2dInteger


   subroutine MPAS_streamAddField_3dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field3DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field3dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), MPAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), MPAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_INT
      new_field_list_node % int3dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_3dInteger


   subroutine MPAS_streamAddField_0dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field0dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      integer :: local_precision

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = 0

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      isDecomposed = .false.
      globalDimSize = 0
      totalDimSize = 0

      call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , dimNames, dimSizes, &
                                       field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= MPAS_STREAM_NOERR) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_REAL
      new_field_list_node % real0dField => field

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_0dReal


   subroutine MPAS_streamAddField_1dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field1DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field1dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , dimNames0, &
                                             dimSizes0, field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, precision=local_precision, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_REAL
      new_field_list_node % real1dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_1dReal


   subroutine MPAS_streamAddField_2dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field2DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field2dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, precision=local_precision, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_REAL
      new_field_list_node % real2dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_2dReal


   subroutine MPAS_streamAddField_3dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field3DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field3dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
!write(stderrUnit,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^ we are adding a var-array'
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, precision=local_precision, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_REAL
      new_field_list_node % real3dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'
!write(stderrUnit,*) 'DEBUGGING : Finished adding 3d real field '//trim(field % fieldName)

   end subroutine MPAS_streamAddField_3dReal


   subroutine MPAS_streamAddField_4dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field4DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field4dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
!write(stderrUnit,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^ we are adding a var-array'
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, precision=local_precision, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_4D_REAL
      new_field_list_node % real4dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'
!write(stderrUnit,*) 'DEBUGGING : Finished adding 4d real field '//trim(field % fieldName)

   end subroutine MPAS_streamAddField_4dReal


   subroutine MPAS_streamAddField_5dReal(stream, field, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field5DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field5dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (mpas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = size(field % dimSizes)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
!write(stderrUnit,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^ we are adding a var-array'
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, precision=local_precision, ierr=io_err)
            if (io_err == MPAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == MPAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_5D_REAL
      new_field_list_node % real5dField => field
      new_field_list_node % isAvailable => isAvailable

!write(stderrUnit,*) '... done adding field'
!write(stderrUnit,*) 'DEBUGGING : Finished adding 5d real field '//trim(field % fieldName)

   end subroutine MPAS_streamAddField_5dReal


   subroutine MPAS_streamAddField_0dChar(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      logical :: isDecomposed
      integer :: ndims
      type (field0dChar), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizes
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

!write(stderrUnit,*) '... Adding field '//trim(field % fieldName)//' to stream'

      ndims = 1

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(1))
      dimSizes(1) = 64
      dimNames(1) = 'StrLen'
      isDecomposed = .false.
      globalDimSize = 64
      totalDimSize = 64

      
      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call MPAS_streamAddField_generic(stream, trim(field % constituentNames(i)), MPAS_IO_CHAR, dimNames(1:1), &
                                             dimSizes, field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
         end do
      else
         call MPAS_streamAddField_generic(stream, trim(field % fieldName), MPAS_IO_CHAR, dimNames(1:1), dimSizes, &
                                          field % hasTimeDimension, isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
      end if

      deallocate(indices)
      deallocate(dimSizes)
      if (io_err /= MPAS_STREAM_NOERR) then
          if (present(ierr)) ierr = MPAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_CHAR
      new_field_list_node % char0dField => field

!write(stderrUnit,*) '... done adding field'

   end subroutine MPAS_streamAddField_0dChar


   subroutine MPAS_streamAddField_generic(stream, fieldName, fieldType, dimNames, dimSizes, hasTimeDimension, isDecomposed, &
                                          totalDimSize, globalDimSize, indices, precision, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: fieldName
      integer, intent(in) :: fieldType
      character (len=StrKIND), dimension(:), intent(in) :: dimNames
      integer, dimension(:), intent(in) :: dimSizes
      logical, intent(in) :: hasTimeDimension
      logical, intent(in) :: isDecomposed
      integer, intent(in) :: totalDimSize
      integer, intent(in) :: globalDimSize
      integer, dimension(:), intent(in) :: indices
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: ndims
      integer :: dimTemp
      character (len=StrKIND), dimension(5) :: dimNamesLocal
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizesInq
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      ndims = size(dimNames)

!write(stderrUnit,*) '... field has ', ndims, ' dimensions'

      allocate(new_field_list_node)
      nullify(new_field_list_node % next)

      if (stream % ioDirection == MPAS_IO_WRITE) then
!write(stderrUnit,*) '... defining field'

         !
         ! Define inner dimensions
         !
         do idim = 1, ndims-1
!write(stderrUnit,*) '... defining dimension ', trim(dimNames(idim)), dimSizes(idim)
            write(dimNamesLocal(idim),'(a)') dimNames(idim)
            call MPAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), dimSizes(idim), io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
               if (present(ierr)) ierr = MPAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end do

         !
         ! Define outer-most dimension
         !
         idim = ndims
         if (idim > 0) write(dimNamesLocal(idim),'(a)') dimNames(idim)

         if (isDecomposed) then
            new_field_list_node % totalDimSize = totalDimSize
         else
            new_field_list_node % totalDimSize = globalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed
  
         if (ndims > 0) then
!write(stderrUnit,*) '... defining dimension ', trim(dimNames(idim)), globalDimSize
            call MPAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), globalDimSize, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
               if (present(ierr)) ierr = MPAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end if

         !
         ! Define time dimension if necessary
         !
         if (hasTimeDimension) then
!write(stderrUnit,*) '... defining Time dimension '
            call MPAS_io_def_dim(stream % fileHandle, 'Time', MPAS_IO_UNLIMITED_DIM, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
               if (present(ierr)) ierr = MPAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
            ndims = ndims + 1
            write(dimNamesLocal(ndims),'(a)') 'Time'
         end if

         !
         ! Define variable to low-level interface
         !
!write(stderrUnit,*) '... defining var to low-level interface with ndims ', ndims

         call MPAS_io_def_var(stream % fileHandle, trim(fieldName), fieldType, dimNamesLocal(1:ndims), precision=precision, ierr=io_err)
         call MPAS_io_err_mesg(io_err, .false.)
         if (io_err /= MPAS_IO_NOERR) then
            if (present(ierr)) ierr = MPAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
        
      else if (stream % ioDirection == MPAS_IO_READ) then
!write(stderrUnit,*) '... inquiring about'

         call MPAS_io_inq_var(stream % fileHandle, trim(fieldName), dimnames=dimNamesInq, dimsizes=dimSizesInq, ierr=io_err)
         ! If the field does not exist in the input file, we should handle this situation gracefully at higher levels
         !   without printing disconcerting error messages
         !call MPAS_io_err_mesg(io_err, .false.)
         if (io_err /= MPAS_IO_NOERR) then
            if (present(ierr)) ierr = MPAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if

! Here, we should probably do a check to make sure the file agrees with what MPAS expects for the field
         do i=1,ndims
!write(stderrUnit,*) 'Comparing '//trim(dimNames(i))//' '//trim(dimNamesInq(i))
            if (trim(dimNames(i)) /= trim(dimNamesInq(i))) then
!write(stderrUnit,*) 'Mismatched dimension name in field'
               if (present(ierr)) ierr = MPAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
            if (i < ndims) then
               dimTemp = dimSizes(i)
            else
               if (trim(dimNamesInq(i)) == 'nCells' .or. &
                   trim(dimNamesInq(i)) == 'nEdges' .or. &
                   trim(dimNamesInq(i)) == 'nVertices'   &
                  ) then
                  dimTemp = globalDimSize
               else
                  dimTemp = dimSizes(i)
               end if
            end if
!write(stderrUnit,*) 'Comparing ', dimTemp, ' ', dimSizesInq(i)
            if (dimTemp /= dimSizesInq(i)) then
!write(stderrUnit,*) 'Mismatched dimension size in field'
               if (present(ierr)) ierr = MPAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
         end do 

         ! Set outer dimension sizes depending on whether the field is decomposed
         if (isDecomposed) then
            new_field_list_node % totalDimSize = totalDimSize
         else
            new_field_list_node % totalDimSize = globalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed

         deallocate(dimNamesInq)
         deallocate(dimSizesInq)

      end if

      
      !
      ! Set variable indices
      !
      if (ndims > 0 .and. isDecomposed) then
         call MPAS_io_set_var_indices(stream % fileHandle, trim(fieldName), indices, io_err)
         call MPAS_io_err_mesg(io_err, .false.)
         if (io_err /= MPAS_IO_NOERR) then
            if (present(ierr)) ierr = MPAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
      end if


      !
      ! Add field pointer to the list
      !
      if (.not. associated(stream % fieldList)) then
!write(stderrUnit,*) 'Adding field to the head of the list'
         stream % fieldList => new_field_list_node
      else
!write(stderrUnit,*) 'Adding field to the tail of the list'
         field_list_cursor => stream % fieldList
         do while (associated(field_list_cursor % next))
            field_list_cursor => field_list_cursor % next
         end do
         field_list_cursor % next => new_field_list_node
      end if

   end subroutine MPAS_streamAddField_generic


   subroutine MPAS_streamUpdateField_0dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_INT) then
            if (field_cursor % int0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d integer named '//trim(field_cursor % int0dField % fieldname) )
               field_cursor % int0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_0dInteger


   subroutine MPAS_streamUpdateField_1dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field1DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 1d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_1D_INT) then
            if (field_cursor % int1dField % fieldname == field % fieldname .and. &
                field_cursor % int1dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int1dField % dimNames(1) == field % dimNames(1)) then
               STREAM_DEBUG_WRITE( '... found 1d integer named '//trim(field_cursor % int1dField % fieldname) )
               field_cursor % int1dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 1d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_1dInteger


   subroutine MPAS_streamUpdateField_2dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field2DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 2d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_2D_INT) then
            if (field_cursor % int2dField % fieldname == field % fieldname .and. &
                field_cursor % int2dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int2dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % int2dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % int2dField % dimNames(2) == field % dimNames(2)) then
               STREAM_DEBUG_WRITE( '... found 2d integer named '//trim(field_cursor % int2dField % fieldname) )
               field_cursor % int2dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 2d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_2dInteger


   subroutine MPAS_streamUpdateField_3dInteger(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field3DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 3d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_3D_INT) then
            if (field_cursor % int3dField % fieldname == field % fieldname .and. &
                field_cursor % int3dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int3dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % int3dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % int3dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % int3dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % int3dField % dimNames(3) == field % dimNames(3)) then
               STREAM_DEBUG_WRITE( '... found 3d integer named '//trim(field_cursor % int3dField % fieldname) )
               field_cursor % int3dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 3d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_3dInteger


   subroutine MPAS_streamUpdateField_0dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_REAL) then
            if (field_cursor % real0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d real named '//trim(field_cursor % real0dField % fieldname) )
               field_cursor % real0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_0dReal


   subroutine MPAS_streamUpdateField_1dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field1DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 1d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_1D_REAL) then
            if (field_cursor % real1dField % fieldname == field % fieldname .and. &
                field_cursor % real1dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real1dField % dimNames(1) == field % dimNames(1)) then
               STREAM_DEBUG_WRITE( '... found 1d real named '//trim(field_cursor % real1dField % fieldname) )
               field_cursor % real1dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 1d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_1dReal


   subroutine MPAS_streamUpdateField_2dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field2DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 2d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_2D_REAL) then
            if (field_cursor % real2dField % fieldname == field % fieldname .and. &
                field_cursor % real2dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real2dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real2dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real2dField % dimNames(2) == field % dimNames(2)) then
               STREAM_DEBUG_WRITE( '... found 2d real named '//trim(field_cursor % real2dField % fieldname) )
               field_cursor % real2dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 2d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_2dReal


   subroutine MPAS_streamUpdateField_3dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field3DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 3d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_3D_REAL) then
            if (field_cursor % real3dField % fieldname == field % fieldname .and. &
                field_cursor % real3dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real3dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real3dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real3dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real3dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real3dField % dimNames(3) == field % dimNames(3)) then
               STREAM_DEBUG_WRITE( '... found 3d real named '//trim(field_cursor % real3dField % fieldname) )
               field_cursor % real3dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 3d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_3dReal


   subroutine MPAS_streamUpdateField_4dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field4DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 4d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_4D_REAL) then
            if (field_cursor % real4dField % fieldname == field % fieldname .and. &
                field_cursor % real4dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real4dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real4dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real4dField % dimSizes(4) == field % dimSizes(4) .and. &
                field_cursor % real4dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real4dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real4dField % dimNames(3) == field % dimNames(3) .and. &
                field_cursor % real4dField % dimNames(4) == field % dimNames(4)) then
               STREAM_DEBUG_WRITE( '... found 4d real named '//trim(field_cursor % real4dField % fieldname) )
               field_cursor % real4dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 4d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_4dReal


   subroutine MPAS_streamUpdateField_5dReal(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field5DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 5d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_5D_REAL) then
            if (field_cursor % real5dField % fieldname == field % fieldname .and. &
                field_cursor % real5dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real5dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real5dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real5dField % dimSizes(4) == field % dimSizes(4) .and. &
                field_cursor % real5dField % dimSizes(5) == field % dimSizes(5) .and. &
                field_cursor % real5dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real5dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real5dField % dimNames(3) == field % dimNames(3) .and. &
                field_cursor % real5dField % dimNames(4) == field % dimNames(4) .and. &
                field_cursor % real5dField % dimNames(5) == field % dimNames(5)) then
               STREAM_DEBUG_WRITE( '... found 5d real named '//trim(field_cursor % real5dField % fieldname) )
               field_cursor % real5dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 5d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_5dReal


   subroutine MPAS_streamUpdateField_0dChar(stream, field, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      type (field0DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d char field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_CHAR) then
            if (field_cursor % char0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d char named '//trim(field_cursor % char0dField % fieldname) )
               field_cursor % char0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d char field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = MPAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine MPAS_streamUpdateField_0dChar


   subroutine MPAS_readStream(stream, frame, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i, j
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                            :: real0d_temp
      real (kind=RKIND), dimension(:),     pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),   pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:), pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:), pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if


      !
      ! Set time frame to real
      !
      call MPAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if


      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_INT) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % int0dField % fieldName)
!write(stderrUnit,*) 'Reading in field '//trim(field_cursor % int0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'MGD calling MPAS_io_get_var now...'
            call MPAS_io_get_var(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
                if (present(ierr)) ierr = MPAS_IO_ERR
                return
            end if

!write(stderrUnit,*) 'Distributing and Copying field to other blocks'

            call mpas_dmpar_bcast_int(field_cursor % int0dField % block % domain % dminfo, int0d_temp)
            field_0dint_ptr => field_cursor % int0dField
            do while (associated(field_0dint_ptr))
               field_0dint_ptr % scalar = int0d_temp
               field_0dint_ptr => field_0dint_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_INT) then
!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % int1dField % fieldName)

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (.not. field_cursor % int1dField % isVarArray) then
                      deallocate(int1d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_1dint_ptr % dimSizes(1)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
                        field_1dint_ptr % array(j) = int0d_temp
                     else
                        field_1dint_ptr % array(1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dint_ptr => field_1dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int1dField % isVarArray) then
                     call mpas_dmpar_bcast_int(field_cursor % int1dField % block % domain % dminfo, int0d_temp)
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(j) = int0d_temp
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_ints(field_cursor % int1dField % block % domain % dminfo, size(int1d_temp), int1d_temp(:))
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(:) = int1d_temp(:)
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % int2dField % fieldName)
            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % int2dField % isVarArray) then
                      deallocate(int1d_temp)
                   else
                      deallocate(int2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_2dint_ptr % dimSizes(2)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        field_2dint_ptr % array(j,1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     else
                        field_2dint_ptr % array(:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dint_ptr => field_2dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int2dField % isVarArray) then
                     call mpas_dmpar_bcast_ints(field_cursor % int2dField % block % domain % dminfo, size(int1d_temp), int1d_temp(:))
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(j,:) = int1d_temp(:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_ints(field_cursor % int2dField % block % domain % dminfo, size(int2d_temp), int2d_temp(:,1))
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(:,:) = int2d_temp(:,:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % int3dField % fieldName)
            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % int3dField % isVarArray) then
                      deallocate(int2d_temp)
                   else
                      deallocate(int3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_3dint_ptr % dimSizes(3)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        field_3dint_ptr % array(j,:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dint_ptr % array(:,:,1:ownedSize) = int3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dint_ptr => field_3dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int3dField % isVarArray) then
                     call mpas_dmpar_bcast_ints(field_cursor % int3dField % block % domain % dminfo, size(int2d_temp), int2d_temp(:,1))
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(j,:,:) = int2d_temp(:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_ints(field_cursor % int3dField % block % domain % dminfo, size(int3d_temp), int3d_temp(:,1,1))
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(:,:,:) = int3d_temp(:,:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real0dField % fieldName)
!write(stderrUnit,*) 'Reading in field '//trim(field_cursor % real0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'MGD calling MPAS_io_get_var now...'
            call MPAS_io_get_var(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
                if (present(ierr)) ierr = MPAS_IO_ERR
                return
            end if

!write(stderrUnit,*) 'Distributing and Copying field to other blocks'

            call mpas_dmpar_bcast_real(field_cursor % real0dField % block % domain % dminfo, real0d_temp)
            field_0dreal_ptr => field_cursor % real0dField
            do while (associated(field_0dreal_ptr))
               field_0dreal_ptr % scalar = real0d_temp
               field_0dreal_ptr => field_0dreal_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real1dField % fieldName)
            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (.not. field_cursor % real1dField % isVarArray) then
                      deallocate(real1d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1

                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_1dreal_ptr % dimSizes(1)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
                        field_1dreal_ptr % array(j) = real0d_temp
                     else
                        field_1dreal_ptr % array(1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dreal_ptr => field_1dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real1dField % isVarArray) then
                     call mpas_dmpar_bcast_real(field_cursor % real1dField % block % domain % dminfo, real0d_temp)
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(j) = real0d_temp
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_reals(field_cursor % real1dField % block % domain % dminfo, size(real1d_temp), real1d_temp(:))
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(:) = real1d_temp(:)
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real2dField % fieldName)
            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % real2dField % isVarArray) then
                      deallocate(real1d_temp)
                   else
                      deallocate(real2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_2dreal_ptr % dimSizes(2)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        field_2dreal_ptr % array(j,1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     else
                        field_2dreal_ptr % array(:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dreal_ptr => field_2dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real2dField % isVarArray) then
                     call mpas_dmpar_bcast_reals(field_cursor % real2dField % block % domain % dminfo, size(real1d_temp), real1d_temp(:))
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(j,:) = real1d_temp(:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_reals(field_cursor % real2dField % block % domain % dminfo, size(real2d_temp), real2d_temp(:,1))
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(:,:) = real2d_temp(:,:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName)
!write(stderrUnit,*) 'DEBUGGING : reading a 3d real array'
            if (field_cursor % real3dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : reading a 3d real var-array'
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!write(stderrUnit,*) 'DEBUGGING : calling get_var for a constitutent'
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % real3dField % isVarArray) then
                      deallocate(real2d_temp)
                   else
                      deallocate(real3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_3dreal_ptr % dimSizes(3)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : copying the temporary array'
                        field_3dreal_ptr % array(j,:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dreal_ptr % array(:,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dreal_ptr => field_3dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real3dField % isVarArray) then
                     call mpas_dmpar_bcast_reals(field_cursor % real3dField % block % domain % dminfo, size(real2d_temp), real2d_temp(:,1))
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(j,:,:) = real2d_temp(:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_reals(field_cursor % real3dField % block % domain % dminfo, size(real3d_temp), real3d_temp(:,1,1))
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(:,:,:) = real3d_temp(:,:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if
         else if (field_cursor % field_type == FIELD_4D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName)
!write(stderrUnit,*) 'DEBUGGING : reading a 4d real array'
            if (field_cursor % real4dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : reading a 4d real var-array'
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real4dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!write(stderrUnit,*) 'DEBUGGING : calling get_var for a constitutent'
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % real4dField % isVarArray) then
                      deallocate(real3d_temp)
                   else
                      deallocate(real4d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_4dreal_ptr % dimSizes(4)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : copying the temporary array'
                        field_4dreal_ptr % array(j, :,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     else
                        field_4dreal_ptr % array(:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_4dreal_ptr => field_4dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real3dField % isVarArray) then
                     call mpas_dmpar_bcast_reals(field_cursor % real4dField % block % domain % dminfo, size(real3d_temp), real3d_temp(:,1,1))
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(j,:,:,:) = real3d_temp(:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_reals(field_cursor % real4dField % block % domain % dminfo, size(real4d_temp), real4d_temp(:,1,1,1))
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(:,:,:,:) = real4d_temp(:,:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

         else if (field_cursor % field_type == FIELD_5D_REAL) then

!write(stderrUnit,*) 'DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName)
!write(stderrUnit,*) 'DEBUGGING : reading a 4d real array'
            if (field_cursor % real5dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : reading a 4d real var-array'
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real5dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!write(stderrUnit,*) 'DEBUGGING : calling get_var for a constitutent'
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call MPAS_io_get_var(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR) then
                   if (present(ierr)) ierr = MPAS_IO_ERR
                   if (field_cursor % real5dField % isVarArray) then
                      deallocate(real4d_temp)
                   else
                      deallocate(real5d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_5dreal_ptr % dimSizes(5)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
!write(stderrUnit,*) 'DEBUGGING : copying the temporary array'
                        field_5dreal_ptr % array(j,:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     else
                        field_5dreal_ptr % array(:,:,:,:,1:ownedSize) = real5d_temp(:,:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_5dreal_ptr => field_5dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real5dField % isVarArray) then
                     call mpas_dmpar_bcast_reals(field_cursor % real5dField % block % domain % dminfo, size(real4d_temp), real4d_temp(:,1,1,1))
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(j,:,:,:,:) = real4d_temp(:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  else
                     call mpas_dmpar_bcast_reals(field_cursor % real5dField % block % domain % dminfo, size(real5d_temp), real5d_temp(:,1,1,1,1))
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(:,:,:,:,:) = real5d_temp(:,:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if


         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!write(stderrUnit,*) 'Reading in field '//trim(field_cursor % char0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'MGD calling MPAS_io_get_var now...'
            call MPAS_io_get_var(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR) then
                if (present(ierr)) ierr = MPAS_IO_ERR
                return
            end if

!write(stderrUnit,*) 'Distributing and Copying field to other blocks'

            call mpas_dmpar_bcast_char(field_cursor % char0dField % block % domain % dminfo, field_cursor % char0dField % scalar)
            field_0dchar_ptr => field_cursor % char0dField
            do while (associated(field_0dchar_ptr))
               field_0dchar_ptr % scalar = field_cursor % char0dField % scalar
               field_0dchar_ptr => field_0dchar_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_CHAR) then
         end if
         field_cursor => field_cursor % next
      end do
!write(stderrUnit,*) 'Finished fieldlist loop...'

   end subroutine MPAS_readStream


   subroutine MPAS_writeStream(stream, frame, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i, j
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                                :: real0d_temp
      real (kind=RKIND), dimension(:),         pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),       pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:),     pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:),   pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      !
      ! Set time frame to write
      !
      call MPAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      !
      ! Check whether we will clobber any records
      !
      if (.not. stream % clobberRecords) then
         if (MPAS_io_would_clobber_records(stream % fileHandle, io_err)) then
            if (present(ierr)) ierr = MPAS_STREAM_CLOBBER_RECORD
            return
         end if
      end if

      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))

         if (field_cursor % field_type == FIELD_0D_INT) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % int0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
            int0d_temp = field_cursor % int0dField % scalar

!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_INT) then

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_1dint_ptr % dimSizes(1)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        int0d_temp = field_1dint_ptr % array(j)
                     else
                        int1d_temp(i:i+ownedSize-1) = field_1dint_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_1dint_ptr => field_1dint_ptr % next
                  end do
               else
                  if (field_cursor % int1dField % isVarArray) then
                     int0d_temp = field_cursor % int1dField % array(j)
                  else
                     int1d_temp(:) = field_cursor % int1dField % array(:)
                  end if
               end if

               if (field_cursor % int1dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_2dint_ptr % dimSizes(2)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        int1d_temp(i:i+ownedSize-1) = field_2dint_ptr % array(j,1:ownedSize)
                     else
                        int2d_temp(:,i:i+ownedSize-1) = field_2dint_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_2dint_ptr => field_2dint_ptr % next
                  end do
               else
                  if (field_cursor % int2dField % isVarArray) then
                     int1d_temp(:) = field_cursor % int2dField % array(j,:)
                  else
                     int2d_temp(:,:) = field_cursor % int2dField % array(:,:)
                  end if
               end if

               if (field_cursor % int2dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_3dint_ptr % dimSizes(3)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        int2d_temp(:,i:i+ownedSize-1) = field_3dint_ptr % array(j,:,1:ownedSize)
                     else
                        int3d_temp(:,:,i:i+ownedSize-1) = field_3dint_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_3dint_ptr => field_3dint_ptr % next
                  end do
               else
                  if (field_cursor % int3dField % isVarArray) then
                     int2d_temp(:,:) = field_cursor % int3dField % array(j,:,:)
                  else
                     int3d_temp(:,:,:) = field_cursor % int3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % int3dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % real0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
            real0d_temp = field_cursor % real0dField % scalar

!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_REAL) then

            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1
                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_1dreal_ptr % dimSizes(1)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        real0d_temp = field_1dreal_ptr % array(j)
                     else
                        real1d_temp(i:i+ownedSize-1) = field_1dreal_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_1dreal_ptr => field_1dreal_ptr % next
                  end do
               else
                  if (field_cursor % real1dField % isVarArray) then
                     real0d_temp = field_cursor % real1dField % array(j)
                  else
                     real1d_temp(:) = field_cursor % real1dField % array(:)
                  end if
               end if

               if (field_cursor % real1dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_2dreal_ptr % dimSizes(2)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        real1d_temp(i:i+ownedSize-1) = field_2dreal_ptr % array(j,1:ownedSize)
                     else
                        real2d_temp(:,i:i+ownedSize-1) = field_2dreal_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_2dreal_ptr => field_2dreal_ptr % next
                  end do
               else
                  if (field_cursor % real2dField % isVarArray) then
                     real1d_temp(:) = field_cursor % real2dField % array(j,:)
                  else
                     real2d_temp(:,:) = field_cursor % real2dField % array(:,:)
                  end if
               end if

               if (field_cursor % real2dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_REAL) then

            if (field_cursor % real3dField % isVarArray) then
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_3dreal_ptr % dimSizes(3)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
                        real2d_temp(:,i:i+ownedSize-1) = field_3dreal_ptr % array(j,:,1:ownedSize)
                     else
                        real3d_temp(:,:,i:i+ownedSize-1) = field_3dreal_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_3dreal_ptr => field_3dreal_ptr % next
                  end do
               else
                  if (field_cursor % real3dField % isVarArray) then
                     real2d_temp(:,:) = field_cursor % real3dField % array(j,:,:)
                  else
                     real3d_temp(:,:,:) = field_cursor % real3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % real3dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_4D_REAL) then

            if (field_cursor % real4dField % isVarArray) then
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_4dreal_ptr % dimSizes(4)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
                        real3d_temp(:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(j,:,:,1:ownedSize)
                     else
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_4dreal_ptr => field_4dreal_ptr % next
                  end do
               else
                  if (field_cursor % real4dField % isVarArray) then
                     real3d_temp(:,:,:) = field_cursor % real4dField % array(j,:,:,:)
                  else
                     real4d_temp(:,:,:,:) = field_cursor % real4dField % array(:,:,:,:)
                  end if
               end if

               if (field_cursor % real4dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

         else if (field_cursor % field_type == FIELD_5D_REAL) then

            if (field_cursor % real5dField % isVarArray) then
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        ownedSize = field_5dreal_ptr % dimSizes(5)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(j,:,:,:,1:ownedSize)
                     else
                        real5d_temp(:,:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(:,:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_5dreal_ptr => field_5dreal_ptr % next
                  end do
               else
                  if (field_cursor % real5dField % isVarArray) then
                     real4d_temp(:,:,:,:) = field_cursor % real5dField % array(j,:,:,:,:)
                  else
                     real5d_temp(:,:,:,:,:) = field_cursor % real5dField % array(:,:,:,:,:)
                  end if
               end if

               if (field_cursor % real5dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % char0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_CHAR) then
         end if
         field_cursor => field_cursor % next
      end do


      !
      ! Sync all fields with disk
      !
      call MPAS_io_sync(stream % fileHandle)

   end subroutine MPAS_writeStream


   subroutine MPAS_readStreamAtt_0dInteger(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, intent(out) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_readStreamAtt_0dInteger


   subroutine MPAS_readStreamAtt_1dInteger(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, dimension(:), pointer :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_readStreamAtt_1dInteger


   subroutine MPAS_readStreamAtt_0dReal(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(out) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_readStreamAtt_0dReal


   subroutine MPAS_readStreamAtt_1dReal(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), pointer :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_readStreamAtt_1dReal


   subroutine MPAS_readStreamAtt_text(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      character (len=*), intent(out) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_readStreamAtt_text


   subroutine MPAS_writeStreamAtt_0dInteger(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, intent(in) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_put_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_writeStreamAtt_0dInteger


   subroutine MPAS_writeStreamAtt_1dInteger(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, dimension(:), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_put_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_writeStreamAtt_1dInteger


   subroutine MPAS_writeStreamAtt_0dReal(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_put_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_writeStreamAtt_0dReal


   subroutine MPAS_writeStreamAtt_1dReal(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_put_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_writeStreamAtt_1dReal


   subroutine MPAS_writeStreamAtt_text(stream, attName, attValue, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      character (len=*), intent(in) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_put_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

   end subroutine MPAS_writeStreamAtt_text


   subroutine MPAS_closeStream(stream, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      integer, intent(out), optional :: ierr

      integer :: io_err
      type (field_list_type), pointer :: field_cursor

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      call MPAS_io_close(stream % fileHandle, io_err)
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

!write(stderrUnit,*) 'Deallocating global attribute list'
      call mpas_deallocate_attlist(stream % attList)

!write(stderrUnit,*) 'Deallocating field list'
      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (associated(field_cursor % isAvailable)) then
            deallocate(field_cursor % isAvailable)
!write(stderrUnit,*) 'Deallocating isAvailable array'
         end if
         stream % fieldList => stream % fieldList % next
         deallocate(field_cursor)
         field_cursor => stream % fieldList
      end do

      stream % isInitialized = .false.

   end subroutine MPAS_closeStream


   subroutine mergeArrays(array1, array2)

      implicit none

      integer, dimension(:), pointer :: array1
      integer, dimension(:), intent(in) :: array2

      integer :: n1, n2
      integer, dimension(:), pointer :: newArray

      n1 = size(array1)
      n2 = size(array2)

      allocate(newArray(n1+n2))

      newArray(1:n1) = array1(:)
      newArray(n1+1:n1+n2) = array2(:)

      deallocate(array1)
      array1 => newArray

   end subroutine mergeArrays


   subroutine put_get_field_atts(handle, ioDirection, fieldname, attList)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(in) :: ioDirection
      character (len=*), intent(in) :: fieldname
      type (att_list_type), pointer :: attList

      type (att_list_type), pointer :: att_cursor

      if (.not. associated(attList)) return

      att_cursor => attList
      if (ioDirection == MPAS_IO_WRITE) then
         do while (associated(att_cursor))
            select case (att_cursor % attType)
               case (ATT_INT)
                  call MPAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueInt, fieldname)
               case (ATT_INTA)
                  call MPAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueIntA, fieldname)
               case (ATT_REAL)
                  call MPAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueReal, fieldname)
               case (ATT_REALA)
                  call MPAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueRealA, fieldname)
               case (ATT_TEXT)
                  call MPAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueText, fieldname)
            end select 
            att_cursor => att_cursor % next
         end do
      else
         do while (associated(att_cursor))
            select case (att_cursor % attType)
               case (ATT_INT)
                  call MPAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueInt, fieldname)
               case (ATT_INTA)
                  call MPAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueIntA, fieldname)
               case (ATT_REAL)
                  call MPAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueReal, fieldname)
               case (ATT_REALA)
                  call MPAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueRealA, fieldname)
               case (ATT_TEXT)
                  call MPAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueText, fieldname)
            end select 
            att_cursor => att_cursor % next
         end do
      end if

   end subroutine put_get_field_atts

end module mpas_io_streams
