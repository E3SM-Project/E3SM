! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

#define COMMA ,
#define ACCESSOR_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

!***********************************************************************
!
!  mpas_field_accessor
!
!> \brief Module providing quick access to members of fields by name
!> \author  Michael Duda, Doug Jacobsen
!> \date    28 March 2016
!> \details 
!> This module provides routines for accessing members of field types
!> (e.g., missingValue) given only the name of the field and a pool 
!> in which the field may be found.
!
!-----------------------------------------------------------------------
module mpas_field_accessor

    use mpas_derived_types, only : mpas_pool_type, mpas_pool_field_info_type, &
                                   MPAS_POOL_REAL, MPAS_POOL_INTEGER, MPAS_POOL_CHARACTER, MPAS_POOL_LOGICAL, &
                                   field0DReal, field1DReal, field2DReal, field3DReal, field4DReal, field5DReal, &
                                   field0DInteger, field1DInteger, field2DInteger, field3DInteger, &
                                   field0DChar, field1DChar, &
                                   field0DLogical
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_pool_routines, only : mpas_pool_get_field_info, mpas_pool_get_field
    use mpas_io_units, only : stderrUnit

    interface mpas_field_access_missing_value
        module procedure mpas_field_access_msgval_real
        module procedure mpas_field_access_msgval_int
        module procedure mpas_field_access_msgval_char
        module procedure mpas_field_access_msgval_logical
    end interface mpas_field_access_missing_value


    contains


    !-----------------------------------------------------------------------
    !  subroutine mpas_field_access_missing_value
    !
    !> \brief Accesses the 'missingValue' member for a field given the field name
    !> \author Doug Jacobsen, Michael Duda
    !> \date   28 March 2016
    !> \details
    !> This routine returns the value of the 'missingValue' member from the field type
    !> for the specified field. The named field must exist in the specified pool; 
    !> if it does not, an error message will be printed.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_field_access_msgval_real(fieldPool, fieldName, missingValue)

        implicit none

        type (mpas_pool_type), intent(in) :: fieldPool
        character(len=*), intent(in) :: fieldName
        real(kind=RKIND), intent(out) :: missingValue

        type (mpas_pool_field_info_type) :: fieldInfo
        type (field0DReal), pointer :: r0 => null()
        type (field1DReal), pointer :: r1 => null()
        type (field2DReal), pointer :: r2 => null()
        type (field3DReal), pointer :: r3 => null()
        type (field4DReal), pointer :: r4 => null()
        type (field5DReal), pointer :: r5 => null()


        ! Initialize fieldType so we can detect whether returned info is valid
        fieldInfo % fieldType = MPAS_POOL_REAL - 1
        call mpas_pool_get_field_info(fieldPool, trim(fieldName), fieldInfo)

        if (fieldInfo % fieldType /= MPAS_POOL_REAL) then
            ACCESSOR_ERROR_WRITE('Cannot return missingValue for field '//trim(fieldName))
            ACCESSOR_ERROR_WRITE('Either '//trim(fieldName)//' was not found in the specified pool or is not a real-type field')
            return
        end if

        ! At this point, we know that the field exists in the pool and is a real-valued field, 
        ! so we should not need extensive error checking below...

        select case(fieldInfo % nDims)
            case (0)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r0)
                missingValue = r0 % missingValue
            case (1)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r1)
                missingValue = r1 % missingValue
            case (2)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r2)
                missingValue = r2 % missingValue
            case (3)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r3)
                missingValue = r3 % missingValue
            case (4)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r4)
                missingValue = r4 % missingValue
            case (5)
                call mpas_pool_get_field(fieldPool, trim(fieldName), r5)
                missingValue = r5 % missingValue
            case default
                ACCESSOR_ERROR_WRITE('Unhandled dimensionality (6-d or more) in mpas_field_access_msgval_real')
        end select

    end subroutine mpas_field_access_msgval_real


    !-----------------------------------------------------------------------
    !  subroutine mpas_field_access_missing_value
    !
    !> \brief Accesses the 'missingValue' member for a field given the field name
    !> \author Doug Jacobsen, Michael Duda
    !> \date   28 March 2016
    !> \details
    !> This routine returns the value of the 'missingValue' member from the field type
    !> for the specified field. The named field must exist in the specified pool; 
    !> if it does not, an error message will be printed.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_field_access_msgval_int(fieldPool, fieldName, missingValue)

        implicit none

        type (mpas_pool_type), intent(in) :: fieldPool
        character(len=*), intent(in) :: fieldName
        integer, intent(out) :: missingValue

        type (mpas_pool_field_info_type) :: fieldInfo
        type (field0DInteger), pointer :: i0 => null()
        type (field1DInteger), pointer :: i1 => null()
        type (field2DInteger), pointer :: i2 => null()
        type (field3DInteger), pointer :: i3 => null()


        ! Initialize fieldType so we can detect whether returned info is valid
        fieldInfo % fieldType = MPAS_POOL_INTEGER - 1
        call mpas_pool_get_field_info(fieldPool, trim(fieldName), fieldInfo)

        if (fieldInfo % fieldType /= MPAS_POOL_INTEGER) then
            ACCESSOR_ERROR_WRITE('Cannot return missingValue for field '//trim(fieldName))
            ACCESSOR_ERROR_WRITE('Either '//trim(fieldName)//' was not found in the specified pool or is not an integer-type field')
            return
        end if

        ! At this point, we know that the field exists in the pool and is an integer-valued field, 
        ! so we should not need extensive error checking below...

        select case(fieldInfo % nDims)
            case (0)
                call mpas_pool_get_field(fieldPool, trim(fieldName), i0)
                missingValue = i0 % missingValue
            case (1)
                call mpas_pool_get_field(fieldPool, trim(fieldName), i1)
                missingValue = i1 % missingValue
            case (2)
                call mpas_pool_get_field(fieldPool, trim(fieldName), i2)
                missingValue = i2 % missingValue
            case (3)
                call mpas_pool_get_field(fieldPool, trim(fieldName), i3)
                missingValue = i3 % missingValue
            case default
                ACCESSOR_ERROR_WRITE('Unhandled dimensionality (4-d or more) in mpas_field_access_msgval_int')
        end select

    end subroutine mpas_field_access_msgval_int


    !-----------------------------------------------------------------------
    !  subroutine mpas_field_access_missing_value
    !
    !> \brief Accesses the 'missingValue' member for a field given the field name
    !> \author Doug Jacobsen, Michael Duda
    !> \date   28 March 2016
    !> \details
    !> This routine returns the value of the 'missingValue' member from the field type
    !> for the specified field. The named field must exist in the specified pool; 
    !> if it does not, an error message will be printed.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_field_access_msgval_char(fieldPool, fieldName, missingValue)

        implicit none

        type (mpas_pool_type), intent(in) :: fieldPool
        character(len=*), intent(in) :: fieldName
        character(len=*), intent(out) :: missingValue

        type (mpas_pool_field_info_type) :: fieldInfo
        type (field0DChar), pointer :: c0 => null()
        type (field1DChar), pointer :: c1 => null()


        ! Initialize fieldType so we can detect whether returned info is valid
        fieldInfo % fieldType = MPAS_POOL_CHARACTER - 1
        call mpas_pool_get_field_info(fieldPool, trim(fieldName), fieldInfo)

        if (fieldInfo % fieldType /= MPAS_POOL_CHARACTER) then
            ACCESSOR_ERROR_WRITE('Cannot return missingValue for field '//trim(fieldName))
            ACCESSOR_ERROR_WRITE('Either '//trim(fieldName)//' was not found in the specified pool or is not a char-type field')
            return
        end if

        ! At this point, we know that the field exists in the pool and is a character-valued field, 
        ! so we should not need extensive error checking below...

        select case(fieldInfo % nDims)
            case (0)
                call mpas_pool_get_field(fieldPool, trim(fieldName), c0)
                if (len(missingValue) < len_trim(c0 % missingValue)) then
                    ACCESSOR_ERROR_WRITE('Truncating missingValue for field '//trim(fieldName))
                    ACCESSOR_ERROR_WRITE('Actual argument for missingValue is too short')
                    missingValue = c0 % missingValue(1:len(missingValue))
                else
                    missingValue = trim(c0 % missingValue)
                end if
            case (1)
                call mpas_pool_get_field(fieldPool, trim(fieldName), c1)
                if (len(missingValue) < len_trim(c1 % missingValue)) then
                    ACCESSOR_ERROR_WRITE('Truncating missingValue for field '//trim(fieldName))
                    ACCESSOR_ERROR_WRITE('Actual argument for missingValue is too short')
                    missingValue = c1 % missingValue(1:len(missingValue))
                else
                    missingValue = trim(c1 % missingValue)
                end if
            case default
                ACCESSOR_ERROR_WRITE('Unhandled dimensionality (2-d or more) in mpas_field_access_msgval_char')
        end select

    end subroutine mpas_field_access_msgval_char


    !-----------------------------------------------------------------------
    !  subroutine mpas_field_access_missing_value
    !
    !> \brief Accesses the 'missingValue' member for a field given the field name
    !> \author Doug Jacobsen, Michael Duda
    !> \date   28 March 2016
    !> \details
    !> This routine returns the value of the 'missingValue' member from the field type
    !> for the specified field. The named field must exist in the specified pool; 
    !> if it does not, an error message will be printed.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_field_access_msgval_logical(fieldPool, fieldName, missingValue)

        implicit none

        type (mpas_pool_type), intent(in) :: fieldPool
        character(len=*), intent(in) :: fieldName
        logical, intent(out) :: missingValue

        type (mpas_pool_field_info_type) :: fieldInfo
        type (field0DLogical), pointer :: l0 => null()


#ifdef POOL_LOGICAL_FIELD_SUPPORT
        ! Initialize fieldType so we can detect whether returned info is valid
        fieldInfo % fieldType = MPAS_POOL_LOGICAL - 1
        call mpas_pool_get_field_info(fieldPool, trim(fieldName), fieldInfo)

        if (fieldInfo % fieldType /= MPAS_POOL_LOGICAL) then
            ACCESSOR_ERROR_WRITE('Cannot return missingValue for field '//trim(fieldName))
            ACCESSOR_ERROR_WRITE('Either '//trim(fieldName)//' was not found in the specified pool or is not a logical-type field')
            return
        end if

        ! At this point, we know that the field exists in the pool and is a logical-valued field, 
        ! so we should not need extensive error checking below...

        select case(fieldInfo % nDims)
            case (0)
                call mpas_pool_get_field(fieldPool, trim(fieldName), l0)
                missingValue = l0 % missingValue
            case default
                ACCESSOR_ERROR_WRITE('Unhandled dimensionality (1-d or more) in mpas_field_access_msgval_logical')
        end select
#else
         ACCESSOR_ERROR_WRITE('Support for accessing missingValue for logical fields is not implemented')
#endif

    end subroutine mpas_field_access_msgval_logical


end module mpas_field_accessor
