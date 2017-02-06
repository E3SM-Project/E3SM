! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_deactivate_unneeded_fields
!
!> \brief
!> \author Adrian K. Turner
!> \date   6th December 2016
!> \details
!>
!-----------------------------------------------------------------------

module cice_deactivate_unneeded_fields

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   implicit none
   private
   save

   public :: &
        cice_deactivate_unneeded_stream_fields

contains

!***********************************************************************
!
!  routine cice_deactivate_unneeded_stream_fields
!
!> \brief
!> \author  Adrian K. Turner
!> \date    6th December 2016
!> \details
!
!-----------------------------------------------------------------------

  subroutine cice_deactivate_unneeded_stream_fields(domain, poolName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         poolName

    type(MPAS_pool_type), pointer :: &
         pool

    type(MPAS_pool_iterator_type) :: &
         poolItr

    type(MPAS_pool_field_info_type) :: &
           fieldInfo

    logical :: &
         fieldInStream, &
         activeProperty, &
         isActive

    integer :: &
         ierrStream, &
         ierrField, &
         directionProperty

    character(len=strKIND) :: &
         streamID, &
         fieldName

    type(field1DReal), pointer :: &
         variable1DReal

    type(field2DReal), pointer :: &
         variable2DReal

    ! loop over area variables
    call MPAS_pool_get_subpool(domain % blocklist % structs, trim(poolName), pool)

    call mpas_pool_begin_iteration(pool)
    do while(mpas_pool_get_next_member(pool, poolItr))

       ! check for the actual fields
       if (poolItr % memberType == MPAS_POOL_FIELD) then

          fieldInStream = .false.

          ! start the stream iteration
          call MPAS_stream_mgr_begin_iteration(domain % streamManager, ierr=ierrStream)
          do while (MPAS_stream_mgr_get_next_stream(domain % streamManager, streamID, directionProperty, activeProperty))

             if (activeProperty .and. directionProperty == MPAS_STREAM_OUTPUT) then

                ! iterate over fields in stream
                call MPAS_stream_mgr_begin_iteration(domain % streamManager, streamID, ierr=ierrField)
                do while (MPAS_stream_mgr_get_next_field(domain % streamManager, streamID, fieldName, isActive))

                   ! see if area field is in stream
                   if (trim(fieldName) == trim(poolItr % memberName) .or. &
                      index(trim(fieldName), "_"//trim(poolItr % memberName)//"_") /= 0) then

                      fieldInStream = .true.

                   endif

                enddo ! stream field iteration

             endif ! active stream

          enddo ! stream iteration

          ! make field inactive if not used
          if (.not. fieldInStream) then

             call MPAS_pool_get_field_info(pool, trim(poolItr % memberName), fieldInfo)

             if (fieldInfo % nDims == 1) then

                call MPAS_pool_get_field(pool, trim(poolItr % memberName), variable1DReal)
                variable1DReal % isActive = .false.
                deallocate(variable1DReal % array)

             else if (fieldInfo % nDims == 2) then

                call MPAS_pool_get_field(pool, trim(poolItr % memberName), variable2DReal)
                variable2DReal % isActive = .false.
                deallocate(variable2DReal % array)

             endif

          endif

       endif ! fields only

    enddo ! area field iteration

  end subroutine cice_deactivate_unneeded_stream_fields

  !-----------------------------------------------------------------------

end module cice_deactivate_unneeded_fields
