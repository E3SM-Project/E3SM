module avect_wrapper_mod
  ! This module supports building attribute vectors for use in unit tests, as well as
  ! performing other operations on attribute vectors.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mct_mod

  implicit none
  private
  save


  ! The following two routines are the same, except for the meaning of the two dimensions
  ! of the 'data' array
  public :: create_aVect_with_data_rows_are_points    ! creates an attribute vector with a given set of real-valued fields, and fills it with the given data
  public :: create_aVect_with_data_rows_are_fields    ! creates an attribute vector with a given set of real-valued fields, and fills it with the given data

  public :: create_aVect_without_data ! creates an attribute vector with a given set of real-valued fields
  public :: aVect_importRattr         ! wrapper to mct_aVect_importRattr which doesn't require a pointer input
  public :: aVect_exportRattr         ! wrapper to mct_aVect_exportRattr which doesn't require pointer management for the output

contains

  !-----------------------------------------------------------------------
  subroutine create_aVect_with_data_rows_are_points(av, attr_tags, data)
    !
    ! !DESCRIPTION:
    ! Creates an attribute vector with a given set of fields, which are all assumed to be
    ! real-valued. Then fills it with the given data.
    !
    ! The data should be given as a 2-d array, [point, field]. So the second dimension
    ! should be the same size as the attr_tags array, with data(:,i) being used to fill
    ! the attr_tags(i) variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: av
    character(len=*), intent(in) :: attr_tags(:)
    real(r8), intent(in) :: data(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: nfields
    integer :: npoints
    integer :: field_index

    character(len=*), parameter :: subname = 'create_aVect_with_data_rows_are_points'
    !-----------------------------------------------------------------------

    npoints = size(data, 1)
    nfields = size(data, 2)

    if (size(attr_tags) /= nfields) then
       print *, subname, ' ERROR: dimensionality mismatch between attr_tags and data'
       stop
    end if

    call create_aVect_without_data(av, attr_tags, npoints)

    do field_index = 1, nfields
       call aVect_importRattr(av, trim(attr_tags(field_index)), data(:,field_index))
    end do

  end subroutine create_aVect_with_data_rows_are_points

  !-----------------------------------------------------------------------
  subroutine create_aVect_with_data_rows_are_fields(av, attr_tags, data)
    !
    ! !DESCRIPTION:
    ! Creates an attribute vector with a given set of fields, which are all assumed to be
    ! real-valued. Then fills it with the given data.
    !
    ! The data should be given as a 2-d array, [field, point]. So the first dimension
    ! should be the same size as the attr_tags array, with data(i,:) being used to fill
    ! the attr_tags(i) variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: av
    character(len=*), intent(in) :: attr_tags(:)
    real(r8), intent(in) :: data(:,:)
    !-----------------------------------------------------------------------

    call create_aVect_with_data_rows_are_points(av, attr_tags, transpose(data))

  end subroutine create_aVect_with_data_rows_are_fields

  !-----------------------------------------------------------------------
  subroutine create_aVect_without_data(av, attr_tags, lsize)
    !
    ! !DESCRIPTION:
    ! Creates an attribute vector with a given set of fields, with space for the given
    ! number of points in each field.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: av
    character(len=*), intent(in) :: attr_tags(:)
    integer, intent(in) :: lsize
    !
    ! !LOCAL VARIABLES:
    integer :: nfields
    integer :: field_index
    integer :: list_length
    character(len=:), allocatable :: attr_list

    character(len=*), parameter :: subname = 'create_aVect_without_data'
    !-----------------------------------------------------------------------

    nfields = size(attr_tags)
    list_length = nfields * (len(attr_tags) + 1)
    allocate(character(len=list_length) :: attr_list)

    attr_list = trim(attr_tags(1))
    do field_index = 2, nfields
       attr_list = trim(attr_list) // ":" // trim(attr_tags(field_index))
    end do

    call mct_aVect_init(av, rList = attr_list, lsize = lsize)

  end subroutine create_aVect_without_data

  !-----------------------------------------------------------------------
  subroutine aVect_importRattr(av, attr_tag, data)
    !
    ! !DESCRIPTION:
    ! This routine is similar to mct_aVect_importRattr, but it doesn't require a pointer
    ! input - so it is often more convenient than calling mct_aVect_importRattr directly.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: av
    character(len=*), intent(in)   :: attr_tag
    real(r8), intent(in)           :: data(:)
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: data_ptr(:)

    character(len=*), parameter :: subname = 'aVect_importRattr'
    !-----------------------------------------------------------------------

    allocate(data_ptr(size(data)))
    data_ptr(:) = data(:)
    call mct_aVect_importRattr(av, trim(attr_tag), data_ptr)
    deallocate(data_ptr)

  end subroutine aVect_importRattr

  !-----------------------------------------------------------------------
  function aVect_exportRattr(av, attr_tag) result(data)
    !
    ! !DESCRIPTION:
    ! This function is similar to mct_aVect_exportRattr, but (1) it is a function rather
    ! than a subroutine (so that it can be included inline in other statements), and (2)
    ! it handles the pointer management for you, so that the caller doesn't have to.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: data(:)  ! function result
    type(mct_aVect), intent(in) :: av
    character(len=*), intent(in) :: attr_tag
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: data_ptr(:)

    character(len=*), parameter :: subname = 'aVect_exportRattr'
    !-----------------------------------------------------------------------

    nullify(data_ptr)
    call mct_aVect_exportRattr(av, trim(attr_tag), data_ptr)
    data = data_ptr
    deallocate(data_ptr)
  end function aVect_exportRattr


end module avect_wrapper_mod
