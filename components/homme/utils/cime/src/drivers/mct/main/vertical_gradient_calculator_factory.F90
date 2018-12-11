module vertical_gradient_calculator_factory
  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module creates vertical gradient objects

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use vertical_gradient_calculator_2nd_order, only : vertical_gradient_calculator_2nd_order_type
  use mct_mod

  implicit none
  private

  public :: create_vertical_gradient_calculator_2nd_order

  ! The following routines are public just to support unit testing, and shouldn't be
  ! called from production code
  public :: extract_data_from_attr_vect

contains

  !-----------------------------------------------------------------------
  function create_vertical_gradient_calculator_2nd_order( &
       attr_vect, fieldname, toponame, elevclass_names, elevclass_bounds) &
       result(calculator)
    !
    ! !DESCRIPTION:
    ! Creates and returns a vertical_gradient_calculator_2nd_order_type object.
    !
    ! The attribute vector is assumed to have fields named fieldname //
    ! elevclass_names(1), toponame // elevclass_names(1), etc.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_2nd_order_type) :: calculator  ! function result
    type(mct_aVect)  , intent(in) :: attr_vect ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname ! base name of the field of interest
    character(len=*) , intent(in) :: toponame  ! base name of the topographic field
    character(len=*) , intent(in) :: elevclass_names(:) ! strings corresponding to each elevation class
    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8)         , intent(in) :: elevclass_bounds(0:)
    !
    ! !LOCAL VARIABLES:
    integer :: nelev
    real(r8), allocatable :: field(:,:)
    real(r8), allocatable :: topo(:,:)

    character(len=*), parameter :: subname = 'create_vertical_gradient_calculator_2nd_order'
    !-----------------------------------------------------------------------

    nelev = size(elevclass_names)
    SHR_ASSERT_ALL_FL((ubound(elevclass_bounds) == (/nelev/)), __FILE__, __LINE__)

    call extract_data_from_attr_vect(attr_vect, fieldname, toponame, elevclass_names, &
         field, topo)

    calculator = vertical_gradient_calculator_2nd_order_type( &
         field = field, topo = topo, elevclass_bounds = elevclass_bounds)

  end function create_vertical_gradient_calculator_2nd_order

  !-----------------------------------------------------------------------
  subroutine extract_data_from_attr_vect(attr_vect, fieldname, toponame, elevclass_names, &
       field_extracted, topo_extracted)
    !
    ! !DESCRIPTION:
    ! Extract topo and data from attribute vector.
    !
    ! Allocates and sets topo_extracted and data_extracted
    !
    ! The attribute vector is assumed to have fields named fieldname //
    ! elevclass_names(1), toponame // elevclass_names(1), etc.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_aVect)  , intent(in) :: attr_vect ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname ! base name of the field of interest
    character(len=*) , intent(in) :: toponame  ! base name of the topographic field
    character(len=*) , intent(in) :: elevclass_names(:) ! strings corresponding to each elevation class

    ! field_extracted(i,j) is point i, elevation class j; same for topo_extracted
    ! these are both allocated here
    real(r8), intent(out), allocatable :: field_extracted(:,:)
    real(r8), intent(out), allocatable :: topo_extracted(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: npts
    integer :: nelev
    integer :: ec
    character(len=:), allocatable :: fieldname_ec
    character(len=:), allocatable :: toponame_ec

    ! The following temporary array is needed because mct wants pointers
    real(r8), pointer :: temp(:)
    
    character(len=*), parameter :: subname = 'extract_data_from_attr_vect'
    !-----------------------------------------------------------------------

    nelev = size(elevclass_names)
    npts = mct_aVect_lsize(attr_vect)

    allocate(field_extracted(npts, nelev))
    allocate(topo_extracted(npts, nelev))
    allocate(temp(npts))

    do ec = 1, nelev
       fieldname_ec = trim(fieldname) // trim(elevclass_names(ec))
       call mct_aVect_exportRattr(attr_vect, fieldname_ec, temp)
       field_extracted(:,ec) = temp(:)

       toponame_ec = trim(toponame) // trim(elevclass_names(ec))
       call mct_aVect_exportRattr(attr_vect, toponame_ec, temp)
       topo_extracted(:,ec) = temp(:)
    end do

    deallocate(temp)

  end subroutine extract_data_from_attr_vect


end module vertical_gradient_calculator_factory
