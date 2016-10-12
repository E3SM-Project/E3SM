! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_pointwise_stats
!
!> \brief MPAS-seaice analysis mode member: pointwise_stats
!> \author Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date   Jan 2016
!> \details
!>  MPAS-seaice analysis mode member: pointwise_stats
!>
!-----------------------------------------------------------------------

#define COMMA ,
#define PWAM_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define PWAM_WARNING_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define PWAM_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

module cice_pointwise_stats

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   !use cice_constants
   !use cice_diagnostics_routines

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: cice_bootstrap_pointwise_stats, &
             cice_init_pointwise_stats, &
             cice_precompute_pointwise_stats, &
             cice_compute_pointwise_stats, &
             cice_restart_pointwise_stats, &
             cice_finalize_pointwise_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=*), parameter :: AMPoolSuffix = 'FieldMapping'
   character (len=*), parameter :: AMFieldNameSuffix = 'PointStats'

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_pointwise_stats
!
!> \brief   Bootstrap pointwise stats AM
!> \author  Doug Jacobsen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    02/03/2016
!> \details
!>  This routine builds out the fields needed for the pointwise stats AM.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_pointwise_stats(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character (len=*), parameter :: AMName = 'pointwiseStats'
      character (len=StrKIND), pointer :: config_AM_output_stream

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: AMInputFields

      type (field0DInteger), pointer :: ptr0DInt, ptr0DIntNew
      type (field1DInteger), pointer :: ptr1DInt, ptr1DIntNew
      type (field2DInteger), pointer :: ptr2DInt, ptr2DIntNew
      type (field3DInteger), pointer :: ptr3DInt, ptr3DIntNew
      type (field0DReal), pointer :: ptr0DReal, ptr0DRealNew
      type (field1DReal), pointer :: ptr1DReal, ptr1DRealNew
      type (field2DReal), pointer :: ptr2DReal, ptr2DRealNew
      type (field3DReal), pointer :: ptr3DReal, ptr3DRealNew
      type (field4DReal), pointer :: ptr4DReal, ptr4DRealNew
      type (field5DReal), pointer :: ptr5DReal, ptr5DRealNew

      type (mpas_pool_iterator_type) :: poolItr
      type (mpas_pool_field_info_type) :: fieldInfo
      character (len=StrKIND), pointer :: charPtr
      character (len=StrKIND) :: fieldName
      logical :: fieldActive

      logical :: keepField
      integer :: iConst, iDim, iDim2

      integer, pointer :: nPoints

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_AM_' // trim(AMName) // '_output_stream', config_AM_output_stream)

      call mpas_pool_create_pool(AMInputFields)

      if ( config_AM_output_stream /= 'none' ) then
         call mpas_stream_mgr_begin_iteration(domain % streamManager, config_AM_output_stream)
         do while ( mpas_stream_mgr_get_next_field(domain % streamManager, config_AM_output_stream, fieldNAme, fieldActive) )
            keepField = .false.
            if ( fieldActive ) then
               call mpas_pool_get_field_info(domain % blocklist % allFields, fieldName, fieldInfo)

               if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                  if ( fieldInfo % nDims == 0 ) then
                     ! Can't compute stats for 0D real...
                  else if ( fieldInfo % nDims == 1 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr1DReal)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr1DReal % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr1DReal % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr1DRealNew)
                           call mpas_duplicate_field(ptr1DReal, ptr1DRealNew)

                           ptr1DRealNew % fieldName = trim(ptr1DReal % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr1DRealNew % fieldName))

                           ! Change constituent names, if it's a var array
                           if ( ptr1DRealNew % isVarArray ) then
                              do iConst = 1, size(ptr1DRealNew % constituentNames)
                                 ptr1DRealNew % constituentNames(iConst) = trim(ptr1DReal % constituentNames(iConst)) &
                                                                         // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr1DReal % fieldName), trim(ptr1DRealNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr1DRealNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr1DRealNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr1DRealNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr1DRealNew % dimNames(iDim) = 'nPoints'
                                 ptr1DRealNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr1DRealNew % array(ptr1DRealNew % dimSizes(1)))
                           ptr1DRealNew % array(:) = 0.0_RKIND

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr1DRealNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr1DRealNew % fieldName, ptr1DReal)

                              ptr1DReal % next => ptr1DRealNew
                              ptr1DRealNew % prev => ptr1DReal
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr1DRealNew % fieldName, ptr1DRealNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 2 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr2DReal)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr2DReal % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr2DReal % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr2DRealNew)
                           call mpas_duplicate_field(ptr2DReal, ptr2DRealNew)

                           ptr2DRealNew % fieldName = trim(ptr2DReal % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr2DRealNew % fieldName))

                           if ( ptr2DRealNew % isVarArray ) then
                              do iConst = 1, size(ptr2DRealNew % constituentNames)
                                 ptr2DRealNew % constituentNames(iConst) = trim(ptr2DReal % constituentNames(iConst)) &
                                                                         // trim(AMFieldNameSuffix)
                              end do
                           end if


                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr2DReal % fieldName), trim(ptr2DRealNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr2DRealNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr2DRealNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr2DRealNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr2DRealNew % dimNames(iDim) = 'nPoints'
                                 ptr2DRealNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr2DRealNew % array(ptr2DRealNew % dimSizes(1), ptr2DRealNew % dimSizes(2)))
                           ptr2DRealNew % array(:, :) = 0.0_RKIND

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr2DRealNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr2DRealNew % fieldName, ptr2DReal)

                              ptr2DReal % next => ptr2DRealNew
                              ptr2DRealNew % prev => ptr2DReal
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr2DRealNew % fieldName, ptr2DRealNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 3 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr3DReal)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr3DReal % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr3DReal % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr3DRealNew)
                           call mpas_duplicate_field(ptr3DReal, ptr3DRealNew)

                           ptr3DRealNew % fieldName = trim(ptr3DReal % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr3DRealNew % fieldName))

                           if ( ptr3DRealNew % isVarArray ) then
                              do iConst = 1, size(ptr3DRealNew % constituentNames)
                                 ptr3DRealNew % constituentNames(iConst) = trim(ptr3DReal % constituentNames(iConst)) &
                                                                         // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr3DReal % fieldName), trim(ptr3DRealNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr3DRealNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr3DRealNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr3DRealNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr3DRealNew % dimNames(iDim) = 'nPoints'
                                 ptr3DRealNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr3DRealNew % array(ptr3DRealNew % dimSizes(1), ptr3DRealNew % dimSizes(2), &
                                                         ptr3DRealNew % dimSizes(3)))
                           ptr3DRealNew % array(:, :, :) = 0.0_RKIND

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr3DRealNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr3DRealNew % fieldName, ptr3DReal)

                              ptr3DReal % next => ptr3DRealNew
                              ptr3DRealNew % prev => ptr3DReal
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr3DRealNew % fieldName, ptr3DRealNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 4 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr4DReal)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr4DReal % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr4DReal % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr4DRealNew)
                           call mpas_duplicate_field(ptr4DReal, ptr4DRealNew)

                           ptr4DRealNew % fieldName = trim(ptr4DReal % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr4DRealNew % fieldName))

                           if ( ptr4DRealNew % isVarArray ) then
                              do iConst = 1, size(ptr4DRealNew % constituentNames)
                                 ptr4DRealNew % constituentNames(iConst) = trim(ptr4DReal % constituentNames(iConst)) &
                                                                         // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr4DReal % fieldName), trim(ptr4DRealNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr4DRealNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr4DRealNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr4DRealNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr4DRealNew % dimNames(iDim) = 'nPoints'
                                 ptr4DRealNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr4DRealNew % array(ptr4DRealNew % dimSizes(1), ptr4DRealNew % dimSizes(2), &
                                                         ptr4DRealNew % dimSizes(3), ptr4DRealNew % dimSizes(4)))
                           ptr4DRealNew % array(:, :, :, :) = 0.0_RKIND

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr4DRealNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr4DRealNew % fieldName, ptr4DReal)

                              ptr4DReal % next => ptr4DRealNew
                              ptr4DRealNew % prev => ptr4DReal
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr4DRealNew % fieldName, ptr4DRealNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 5 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr5DReal)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr5DReal % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr5DReal % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr5DRealNew)
                           call mpas_duplicate_field(ptr5DReal, ptr5DRealNew)

                           ptr5DRealNew % fieldName = trim(ptr5DReal % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr5DRealNew % fieldName))

                           if ( ptr5DRealNew % isVarArray ) then
                              do iConst = 1, size(ptr5DRealNew % constituentNames)
                                 ptr5DRealNew % constituentNames(iConst) = trim(ptr5DReal % constituentNames(iConst)) &
                                                                         // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr5DReal % fieldName), trim(ptr5DRealNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr5DRealNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr5DRealNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr5DRealNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr5DRealNew % dimNames(iDim) = 'nPoints'
                                 ptr5DRealNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr5DRealNew % array( ptr5DRealNew % dimSizes(1), ptr5DRealNew % dimSizes(2), &
                                                          ptr5DRealNew % dimSizes(3), ptr5DRealNew % dimSizes(4), &
                                                          ptr5DRealNew % dimSizes(5)))
                           ptr5DRealNew % array(:, :, :, :, :) = 0.0_RKIND

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr5DRealNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr5DRealNew % fieldName, ptr5DReal)

                              ptr5DReal % next => ptr5DRealNew
                              ptr5DRealNew % prev => ptr5DReal
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr5DRealNew % fieldName, ptr5DRealNew)
                        end if
                        block => block % next
                     end do
                  end if
               else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                  if ( fieldInfo % nDims == 0 ) then
                     ! Can't compute stats for 0D integer...
                  else if ( fieldInfo % nDims == 1 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr1DInt)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr1DInt % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr1DInt % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr1DIntNew)
                           call mpas_duplicate_field(ptr1DInt, ptr1DIntNew)

                           ptr1DIntNew % fieldName = trim(ptr1DInt % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr1DIntNew % fieldName))

                           if ( ptr1DIntNew % isVarArray ) then
                              do iConst = 1, size(ptr1DIntNew % constituentNames)
                                 ptr1DIntNew % constituentNames(iConst) = trim(ptr1DInt % constituentNames(iConst)) &
                                                                        // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr1DInt % fieldName), trim(ptr1DIntNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr1DIntNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr1DIntNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr1DIntNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr1DIntNew % dimNames(iDim) = 'nPoints'
                                 ptr1DIntNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr1DIntNew % array(ptr1DIntNew % dimSizes(1)))
                           ptr1DIntNew % array(:) = 0

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr1DIntNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr1DIntNew % fieldName, ptr1DInt)

                              ptr1DInt % next => ptr1DIntNew
                              ptr1DIntNew % prev => ptr1DInt
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr1DIntNew % fieldName, ptr1DIntNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 2 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr2DInt)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr2DInt % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr2DInt % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr2DIntNew)
                           call mpas_duplicate_field(ptr2DInt, ptr2DIntNew)

                           ptr2DIntNew % fieldName = trim(ptr2DInt % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr2DIntNew % fieldName))

                           if ( ptr2DIntNew % isVarArray ) then
                              do iConst = 1, size(ptr2DIntNew % constituentNames)
                                 ptr2DIntNew % constituentNames(iConst) = trim(ptr2DInt % constituentNames(iConst)) &
                                                                        // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr2DInt % fieldName), trim(ptr2DIntNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr2DIntNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr2DIntNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr2DIntNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr2DIntNew % dimNames(iDim) = 'nPoints'
                                 ptr2DIntNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr2DIntNew % array(ptr2DIntNew % dimSizes(1), ptr2DIntNew % dimSizes(2)))
                           ptr2DIntNew % array(:, :) = 0

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr2DIntNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr2DIntNew % fieldName, ptr2DInt)

                              ptr2DInt % next => ptr2DIntNew
                              ptr2DIntNew % prev => ptr2DInt
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr2DIntNew % fieldName, ptr2DIntNew)
                        end if
                        block => block % next
                     end do
                  else if ( fieldInfo % nDims == 3 ) then
                     block => domain % blocklist
                     do while ( associated(block) )
                        call mpas_pool_get_field(block % allFields, fieldName, ptr3DInt)
                        call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)

                        ! Check for nCells as a dimension, since that's all we can use for now...
                        do iDim = 1, fieldInfo % nDims
                           if ( trim(ptr3DInt % dimNames(iDim)) == 'nCells' .or. &
                                trim(ptr3DInt % dimNames(iDim)) == 'nVertices') then
                              keepField = .true.
                           end if
                        end do

                        if ( keepField ) then
                           ! Need to create a copy of the field, And change nCells with nPoints
                           allocate(ptr3DIntNew)
                           call mpas_duplicate_field(ptr3DInt, ptr3DIntNew)

                           ptr3DIntNew % fieldName = trim(ptr3DInt % fieldName) // trim(AMFieldNameSuffix)
                           PWAM_DEBUG_WRITE(' -- New Name: ' // trim(ptr3DIntNew % fieldName))

                           if ( ptr3DIntNew % isVarArray ) then
                              do iConst = 1, size(ptr3DIntNew % constituentNames)
                                 ptr3DIntNew % constituentNames(iConst) = trim(ptr3DInt % constituentNames(iConst)) &
                                                                        // trim(AMFieldNameSuffix)
                              end do
                           end if

                           ! Add input name to input name pool, with a value of the new name
                           call mpas_pool_add_config(AMInputFields, trim(ptr3DInt % fieldName), trim(ptr3DIntNew % fieldName))

                           ! Deallocate array, so we can resize it.
                           deallocate(ptr3DIntNew % array)

                           ! Swap out nCells with nPoints
                           do iDim = 1, fieldInfo % nDims
                              if ( trim(ptr3DIntNew % dimNames(iDim)) == 'nCells' .or. &
                                   trim(ptr3DIntNew % dimNames(iDim)) == 'nVertices') then
                                 PWAM_DEBUG_WRITE(' -- Changing nCells/nVertices to nPoints')
                                 ptr3DIntNew % dimNames(iDim) = 'nPoints'
                                 ptr3DIntNew % dimSizes(iDim) = nPoints
                              end if
                           end do

                           ! Allocate new array size
                           allocate(ptr3DIntNew % array(ptr3DIntNew % dimSizes(1), ptr3DIntNew % dimSizes(2), &
                                    ptr3DIntNew % dimSizes(3)))
                           ptr3DIntNew % array(:, :, :) = 0

                           ! Mark the field as non-decomposed, since nPoints is not a decomposed dimension
                           ptr3DIntNew % isDecomposed = .false.

                           ! Link with previous block
                           if ( associated(block % prev) ) then
                              call mpas_pool_get_field(block % prev % allFields, ptr3DIntNew % fieldName, ptr3DInt)

                              ptr3DInt % next => ptr3DIntNew
                              ptr3DIntNew % prev => ptr3DInt
                           end if

                           ! Add field to allFields pool
                           call mpas_pool_add_field(block % allFields, ptr3DIntNew % fieldName, ptr3DIntNew)
                        end if
                        block => block % next
                     end do
                  end if
               end if
            end if
         end do
      end if

      ! Swap fields in the stream
      call mpas_pool_begin_iteration(AMInputFields)
      do while ( mpas_pool_get_next_member(AMInputFields, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            if ( poolItr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(AMInputFields, poolItr % memberName, charPtr)

               call mpas_stream_mgr_remove_field(domain % streamManager, config_AM_output_stream, poolItr % memberName)
               call mpas_stream_mgr_add_field(domain % streamManager, config_AM_output_stream, charPtr)
            end if
         end if
      end do

      call mpas_pool_add_subpool(domain % blocklist % structs, trim(AMName) // trim(AMPoolSuffix), AMInputFields)
      nullify(AMInputFields)

   end subroutine cice_bootstrap_pointwise_stats!}}}

!***********************************************************************
!
!  routine cice_init_pointwise_stats
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    Jan 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_pointwise_stats(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: pointPool

      integer, pointer :: nCells, nCellsSolve, nVertices, nVerticesSolve, nPoints
      integer :: iCell, iVertex, iPoint, i
      integer, dimension(:), pointer :: indexToCellID,pointCellGlobalID, pointCellLocalID, indexToPointCellLocalID
      integer, dimension(:), pointer :: indexToVertexID,pointVertexGlobalID, pointVertexLocalID, indexToPointVertexLocalID

      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nPoints', nPoints)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_subpool(block % structs, 'pointLocations', pointPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         call mpas_pool_get_array(pointPool, 'pointCellGlobalID', pointCellGlobalID)
         call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointCellLocalID)
         call mpas_pool_get_array(pointPool, 'indexToPointCellLocalID', indexToPointCellLocalID)

         call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)
         call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
         call mpas_pool_get_array(pointPool, 'pointVertexGlobalID', pointVertexGlobalID)
         call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointVertexLocalID)
         call mpas_pool_get_array(pointPool, 'indexToPointVertexLocalID', indexToPointVertexLocalID)

         pointCellLocalID = nCells + 1
         indexToPointCellLocalID = 0

         ! Initialize nCells index arrays to record pointwise data.
         i = 0
         do iCell = 1,nCellsSolve
            do iPoint = 1,nPoints
               if (indexToCellID(iCell).eq.pointCellGlobalID(iPoint)) then
                  i = i + 1
                  indexToPointCellLocalID(i) = iPoint
                  pointCellLocalID(iPoint) = iCell
               endif
            end do
         end do

         pointVertexLocalID = nVertices + 1
         indexToPointVertexLocalID = 0

         ! Initialize nVertices index arrays to record pointwise data.
         i = 0
         do iVertex = 1,nVerticesSolve
            do iPoint = 1,nPoints
               if (indexToVertexID(iVertex).eq.pointVertexGlobalID(iPoint)) then
                  i = i + 1
                  indexToPointVertexLocalID(i) = iPoint
                  pointVertexLocalID(iPoint) = iVertex
               endif
            end do
         end do

         block => block % next
      end do

   end subroutine cice_init_pointwise_stats!}}}

!***********************************************************************
!
!  routine cice_precompute_pointwise_stats
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    Jan 2016
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_pointwise_stats(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

    end subroutine cice_precompute_pointwise_stats

!***********************************************************************
!
!  routine cice_compute_pointwise_stats
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    Jan 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_pointwise_stats(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      integer, pointer :: nPoints, nCells, maxDim
      integer :: iPoint, i
      integer, dimension(:), pointer :: pointCellGlobalID, pointCellLocalID, pointLocalElementID

      character (len=*), parameter :: AMName = 'pointwiseStats'

      type(field0DReal), pointer :: ptr0DReal1
      type(field1DReal), pointer :: ptr1DReal1
      type(field2DReal), pointer :: ptr2DReal1
      type(field3DReal), pointer :: ptr3DReal1
      type(field4DReal), pointer :: ptr4DReal1
      type(field5DReal), pointer :: ptr5DReal1

      real (kind=RKIND), pointer :: ptr0DReal2
      real (kind=RKIND), dimension(:), pointer :: ptr1DReal2
      real (kind=RKIND), dimension(:, :), pointer :: ptr2DReal2
      real (kind=RKIND), dimension(:, :, :), pointer :: ptr3DReal2
      real (kind=RKIND), dimension(:, :, :, :), pointer :: ptr4DReal2
      real (kind=RKIND), dimension(:, :, :, :, :), pointer :: ptr5DReal2

      type(field0DInteger), pointer :: ptr0DInt1
      type(field1DInteger), pointer :: ptr1DInt1
      type(field2DInteger), pointer :: ptr2DInt1
      type(field3DInteger), pointer :: ptr3DInt1

      integer, pointer :: ptr0DInt2
      integer, dimension(:), pointer :: ptr1DInt2
      integer, dimension(:, :), pointer :: ptr2DInt2
      integer, dimension(:, :, :), pointer :: ptr3DInt2

      character (len=StrKIND), pointer :: mappedName
      type (mpas_pool_type), pointer :: AMFieldMapping, pointPool
      type (mpas_pool_iterator_type) :: poolItr
      type (mpas_pool_field_info_type) :: fieldInfo
      integer :: nElements
      integer, dimension(:), pointer :: arrShape

      real (kind=RKIND), dimension(:), pointer :: tempRealArrLocal, tempRealArrGlobal
      integer, dimension(:), pointer :: tempIntArrLocal, tempIntArrGlobal
      integer :: iDim, iElement, iDim1, iDim2, iDim3, iDim4, iDim5
      integer :: fieldTimeLevel

      err = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, trim(AMName) // trim(AMPoolSuffix), AMFieldMapping)

      call mpas_pool_begin_iteration(AMFieldMapping)

      do while ( mpas_pool_get_next_member(AMFieldMapping, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_CONFIG ) then
            call mpas_pool_get_config(AMFieldMapping, poolItr % memberName, mappedName)
            call mpas_pool_get_field_info(domain % blocklist % allFields, poolItr % memberName, fieldInfo)
            call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nPoints', nPoints)
            nElements = 0

            ! Set field time level, for retrieving the correct time level of the field later
            if ( fieldInfo % nTimeLevels < timeLevel ) then
               fieldTimeLevel = 1
            else
               fieldTimeLevel = timeLevel
            end if

            PWAM_DEBUG_WRITE(' -- Building pointer')
            ! Get pointer to field in the first block, to store processor sum in
            if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
               if ( fieldInfo % nDims == 0 ) then
               else if ( fieldInfo % nDims == 1 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr1DReal2)
                  ptr1DReal2(:) = 0.0_RKIND
                  nElements = size(ptr1DReal2)
               else if ( fieldInfo % nDims == 2 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr2DReal2)
                  ptr2DReal2(:, :) = 0.0_RKIND
                  nElements = size(ptr2DReal2)
               else if ( fieldInfo % nDims == 3 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr3DReal2)
                  ptr3DReal2(:, :, :) = 0.0_RKIND
                  nElements = size(ptr3DReal2)
               else if ( fieldInfo % nDims == 4 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr4DReal2)
                  ptr4DReal2(:, :, :, :) = 0.0_RKIND
                  nElements = size(ptr4DReal2)
               else if ( fieldInfo % nDims == 5 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr5DReal2)
                  ptr5DReal2(:, :, :, :, :) = 0.0_RKIND
                  nElements = size(ptr5DReal2)
               end if
            else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
               if ( fieldInfo % nDims == 0 ) then
               else if ( fieldInfo % nDims == 1 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr1DInt2)
                  ptr1DInt2(:) = 0
                  nElements = size(ptr1DInt2)
               else if ( fieldInfo % nDims == 2 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr2DInt2)
                  ptr2DInt2(:, :) = 0
                  nElements = size(ptr2DInt2)
               else if ( fieldInfo % nDims == 3 ) then
                  call mpas_pool_get_array(domain % blocklist % allFields, mappedName, ptr3DInt2)
                  ptr3DInt2(:, :, :) = 0
                  nElements = size(ptr3DInt2)
               end if
            end if

            ! Accumulate point data into the first block's field
            PWAM_DEBUG_WRITE(' -- Accumulating field pointer ' // trim(poolItr % memberName))
            block => domain % blocklist
            do while ( associated(block) )
               call mpas_pool_get_subpool(block % structs, 'pointLocations', pointPool)

               if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                  if ( fieldInfo % nDims == 0 ) then
                     ! Can't do 0D reals currently...
                  else if ( fieldInfo % nDims == 1 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr1DReal1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr1DReal1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr1DReal1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr1DReal2(iPoint) = ptr1DReal1 % array( pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 2 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr2DReal1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr2DReal1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr2DReal1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr2DReal2(:, iPoint) = ptr2DReal1 % array(:, pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 3 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr3DReal1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr3DReal1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr3DReal1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr3DReal2(:, :, iPoint) = ptr3DReal1 % array(:, :, pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 4 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr4DReal1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr4DReal1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr4DReal1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr4DReal2(:, :, :, iPoint) = ptr4DReal1 % array(:, :, :, pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 5 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr5DReal1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr5DReal1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr5DReal1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr5DReal2(:, :, :, :, iPoint) = ptr5DReal1 % array(:, :, :, :, pointLocalElementID(iPoint) )
                        end if
                     end do

                  end if
               else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                  if ( fieldInfo % nDims == 0 ) then
                     ! Can't do 0D ints currently...
                  else if ( fieldInfo % nDims == 1 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr1DInt1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr1DInt1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr1DInt1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr1DInt2(iPoint) = ptr1DInt1 % array( pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 2 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr2DInt1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr2DInt1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr2DInt1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr2DInt2(:, iPoint) = ptr2DInt1 % array(:, pointLocalElementID(iPoint) )
                        end if
                     end do

                  else if ( fieldInfo % nDims == 3 ) then

                     call mpas_pool_get_field(block % allFields, poolItr % memberName, ptr3DInt1, fieldTimeLevel)

                     do iDim = 1, fieldInfo % nDims
                        if (trim(ptr3DInt1 % dimNames(iDim)) == "nCells") then
                           call mpas_pool_get_dimension(block % dimensions, 'nCells', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointCellLocalID', pointLocalElementID)
                        else if (trim(ptr3DInt1 % dimNames(iDim)) == "nVertices") then
                           call mpas_pool_get_dimension(block % dimensions, 'nVertices', maxDim)
                           call mpas_pool_get_array(pointPool, 'pointVertexLocalID', pointLocalElementID)
                        end if
                     enddo ! iDim

                     do iPoint = 1, nPoints
                        if ( pointLocalElementID(iPoint) < maxDim + 1 ) then
                           ptr3DInt2(:, :, iPoint) = ptr3DInt1 % array(:, :, pointLocalElementID(iPoint) )
                        end if
                     end do

                  end if
               end if
               block => block % next
            end do

            PWAM_DEBUG_WRITE(' -- Reducing field ' // trim(poolItr % memberName) // ' with nElements = ' COMMA nElements)
            ! Need to sum field across processors
            if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
               allocate( tempRealArrLocal(nElements) )
               allocate( tempRealArrGlobal(nElements) )
               if ( fieldInfo % nDims == 0 ) then
               else if ( fieldInfo % nDims == 1 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr1DReal2, dim=1)
                     tempRealArrLocal(iElement) = ptr1DReal2(iDim1)

                     iElement = iElement + 1
                  end do

                  call mpas_dmpar_sum_real_array(domain % dminfo, nElements, tempRealArrLocal, tempRealArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr1DReal2, dim=1)
                     ptr1DReal2(iDim1) = tempRealArrGlobal(iElement)
                     iElement = iElement + 1
                  end do
               else if ( fieldInfo % nDims == 2 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr2DReal2, dim=2)
                     do iDim2 = 1, size(ptr2DReal2, dim=1)
                        tempRealArrLocal(iElement) = ptr2DReal2(iDim2, iDim1)

                        iElement = iElement + 1
                     end do
                  end do

                  call mpas_dmpar_sum_real_array(domain % dminfo, nElements, tempRealArrLocal, tempRealArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr2DReal2, dim=2)
                     do iDim2 = 1, size(ptr2DReal2, dim=1)
                        ptr2DReal2(iDim2, iDim1) = tempRealArrGlobal(iElement)
                        iElement = iElement + 1
                     end do
                  end do
               else if ( fieldInfo % nDims == 3 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr3DReal2, dim=3)
                     do iDim2 = 1, size(ptr3DReal2, dim=2)
                        do iDim3 = 1, size(ptr3DReal2, dim=1)
                           tempRealArrLocal(iElement) = ptr3DReal2(iDim3, iDim2, iDim1)

                           iElement = iElement + 1
                        end do
                     end do
                  end do

                  call mpas_dmpar_sum_real_array(domain % dminfo, nElements, tempRealArrLocal, tempRealArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr3DReal2, dim=3)
                     do iDim2 = 1, size(ptr3DReal2, dim=2)
                        do iDim3 = 1, size(ptr3DReal2, dim=1)
                           ptr3DReal2(iDim3, iDim2, iDim1) = tempRealArrGlobal(iElement)
                           iElement = iElement + 1
                        end do
                     end do
                  end do
               else if ( fieldInfo % nDims == 4 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr4DReal2, dim=4)
                     do iDim2 = 1, size(ptr4DReal2, dim=3)
                        do iDim3 = 1, size(ptr4DReal2, dim=2)
                           do iDim4 = 1, size(ptr4DReal2, dim=1)
                              tempRealArrLocal(iElement) = ptr4DReal2(iDim4, iDim3, iDim2, iDim1)

                              iElement = iElement + 1
                           end do
                        end do
                     end do
                  end do

                  call mpas_dmpar_sum_real_array(domain % dminfo, nElements, tempRealArrLocal, tempRealArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr4DReal2, dim=4)
                     do iDim2 = 1, size(ptr4DReal2, dim=3)
                        do iDim3 = 1, size(ptr4DReal2, dim=2)
                           do iDim4 = 1, size(ptr4DReal2, dim=1)
                              ptr4DReal2(iDim4, iDim3, iDim2, iDim1) = tempRealArrGlobal(iElement)
                              iElement = iElement + 1
                           end do
                        end do
                     end do
                  end do
               else if ( fieldInfo % nDims == 5 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr5DReal2, dim=5)
                     do iDim2 = 1, size(ptr5DReal2, dim=4)
                        do iDim3 = 1, size(ptr5DReal2, dim=3)
                           do iDim4 = 1, size(ptr5DReal2, dim=2)
                              do iDim5 = 1, size(ptr5DReal2, dim=1)
                                 tempRealArrLocal(iElement) = ptr5DReal2(iDim5, iDim4, iDim3, iDim2, iDim1)

                                 iElement = iElement + 1
                              end do
                           end do
                        end do
                     end do
                  end do

                  call mpas_dmpar_sum_real_array(domain % dminfo, nElements, tempRealArrLocal, tempRealArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr5DReal2, dim=5)
                     do iDim2 = 1, size(ptr5DReal2, dim=4)
                        do iDim3 = 1, size(ptr5DReal2, dim=3)
                           do iDim4 = 1, size(ptr5DReal2, dim=2)
                              do iDim5 = 1, size(ptr5DReal2, dim=1)
                                 ptr5DReal2(iDim5, iDim4, iDim3, iDim2, iDim1) = tempRealArrGlobal(iElement)
                                 iElement = iElement + 1
                              end do
                           end do
                        end do
                     end do
                  end do
               end if
               deallocate( tempRealArrLocal )
               deallocate( tempRealArrGlobal )
            else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
               allocate( tempIntArrLocal(nElements) )
               allocate( tempIntArrGlobal(nElements) )
               if ( fieldInfo % nDims == 0 ) then
               else if ( fieldInfo % nDims == 1 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr1DInt2, dim=1)
                     tempIntArrLocal(iElement) = ptr1DInt2(iDim1)

                     iElement = iElement + 1
                  end do

                  call mpas_dmpar_sum_int_array(domain % dminfo, nElements, tempIntArrLocal, tempIntArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr1DInt2, dim=1)
                     ptr1DInt2(iDim1) = tempIntArrGlobal(iElement)
                     iElement = iElement + 1
                  end do
               else if ( fieldInfo % nDims == 2 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr2DInt2, dim=2)
                     do iDim2 = 1, size(ptr2DInt2, dim=1)
                        tempIntArrLocal(iElement) = ptr2DInt2(iDim2, iDim1)

                        iElement = iElement + 1
                     end do
                  end do

                  call mpas_dmpar_sum_int_array(domain % dminfo, nElements, tempIntArrLocal, tempIntArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr2DInt2, dim=2)
                     do iDim2 = 1, size(ptr2DInt2, dim=1)
                        ptr2DInt2(iDim2, iDim1) = tempIntArrGlobal(iElement)
                        iElement = iElement + 1
                     end do
                  end do
               else if ( fieldInfo % nDims == 3 ) then
                  ! Pack local array
                  iElement = 1
                  do iDim1 = 1, size(ptr3DInt2, dim=3)
                     do iDim2 = 1, size(ptr3DInt2, dim=2)
                        do iDim3 = 1, size(ptr3DInt2, dim=1)
                           tempIntArrLocal(iElement) = ptr3DInt2(iDim3, iDim2, iDim1)

                           iElement = iElement + 1
                        end do
                     end do
                  end do

                  call mpas_dmpar_sum_int_array(domain % dminfo, nElements, tempIntArrLocal, tempIntArrGlobal)

                  ! Unpack global array
                  iElement = 1
                  do iDim1 = 1, size(ptr3DInt2, dim=3)
                     do iDim2 = 1, size(ptr3DInt2, dim=2)
                        do iDim3 = 1, size(ptr3DInt2, dim=1)
                           ptr3DInt2(iDim3, iDim2, iDim1) = tempIntArrGlobal(iElement)
                           iElement = iElement + 1
                        end do
                     end do
                  end do
               end if
               deallocate( tempIntArrLocal )
               deallocate( tempIntArrGlobal )
            end if
            PWAM_DEBUG_WRITE(' -- Completed field ' // trim(poolItr % memberName))
         end if
      end do

    end subroutine cice_compute_pointwise_stats!}}}

!***********************************************************************
!
!  routine cice_restart_pointwise_stats
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    Jan 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_pointwise_stats(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine cice_restart_pointwise_stats!}}}

!***********************************************************************
!
!  routine cice_finalize_pointwise_stats
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Mark Petersen (adapted for MPAS-seaice by Adrian K. Turner)
!> \date    Jan 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_pointwise_stats(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine cice_finalize_pointwise_stats!}}}

end module cice_pointwise_stats

! vim: foldmethod=marker
