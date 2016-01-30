module vertical_gradient_calculator_2nd_order

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type for
  ! computing vertical gradients using a second-order centered difference.
  !
  ! Currently, this code assumes that it is working with GLC elevation classes.

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use glc_elevclass_mod, only : glc_elevclass_as_string
  use mct_mod
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use shr_sys_mod, only : shr_sys_abort
  
  implicit none
  private

  public :: vertical_gradient_calculator_2nd_order_type
  
  type, extends(vertical_gradient_calculator_base_type) :: &
       vertical_gradient_calculator_2nd_order_type
     private

     integer :: min_elevation_class
     integer :: max_elevation_class
     integer :: num_points
     real(r8), allocatable :: field(:,:)   ! field(i,j) is point i, elevation class j
     real(r8), allocatable :: topo(:,:)    ! topo(i,j) is point i, elevation class j

   contains
     procedure :: calc_vertical_gradient
     
     procedure, private :: set_data_from_attr_vect ! extract data from an attribute vector

  end type vertical_gradient_calculator_2nd_order_type

  interface vertical_gradient_calculator_2nd_order_type
     module procedure constructor
  end interface vertical_gradient_calculator_2nd_order_type
  
contains

  !-----------------------------------------------------------------------
  function constructor(attr_vect, fieldname, toponame, &
       min_elevation_class, max_elevation_class) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates a vertical_gradient_calculator_2nd_order_type object by reading the
    ! necessary data from the provided attribute vector
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vertical_gradient_calculator_2nd_order_type) :: this  ! function result
    type(mct_aVect)  , intent(in) :: attr_vect           ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname           ! base name of the field of interest
    character(len=*) , intent(in) :: toponame            ! base name of the topographic field
    integer          , intent(in) :: min_elevation_class ! first elevation class index
    integer          , intent(in) :: max_elevation_class ! last elevation class index
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%min_elevation_class = min_elevation_class
    this%max_elevation_class = max_elevation_class
    call this%set_data_from_attr_vect(attr_vect, fieldname, toponame)
    
  end function constructor


  !-----------------------------------------------------------------------
  subroutine calc_vertical_gradient(this, elevation_class, vertical_gradient)
    !
    ! !DESCRIPTION:
    ! Calculates the vertical gradient for all points, at a given elevation class.
    !
    ! If the topo values are nearly equal across the gradient (i.e., denominator is near
    ! 0), returns a gradient of 0.
    !
    ! If there is only one elevation class, returns a gradient of 0.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    integer, intent(in) :: elevation_class

    ! vertical_gradient should already be allocated to the appropriate size
    real(r8), intent(out) :: vertical_gradient(:)
    !
    ! !LOCAL VARIABLES:
    ! Tolerance for considering two topo values to be nearly equal
    real(r8), parameter :: topo_equality_tolerance = 1.e-13_r8

    integer :: ec_low   ! elevation class index to use as the lower bound of the gradient
    integer :: ec_high  ! elevation class index to use as the upper bound of the gradient
    
    character(len=*), parameter :: subname = 'calc_vertical_gradient'
    !-----------------------------------------------------------------------

    ! Assert pre-conditions
    
    SHR_ASSERT((size(vertical_gradient) == this%num_points), errMsg(__FILE__, __LINE__))

    if (elevation_class < this%min_elevation_class .or. &
         elevation_class > this%max_elevation_class) then
       write(logunit,*) subname, ': ERROR: elevation class out of bounds: ', &
            elevation_class, this%min_elevation_class, this%max_elevation_class
       call shr_sys_abort(subname//': ERROR: elevation class out of bounds')
    end if

    ! Do the calculations

    if (this%min_elevation_class == this%max_elevation_class) then
       vertical_gradient(:) = 0._r8

    else
       
       if (elevation_class == this%min_elevation_class) then
          ec_low = elevation_class
          ec_high = elevation_class + 1
       else if (elevation_class == this%max_elevation_class) then
          ec_low = elevation_class - 1
          ec_high = elevation_class
       else
          ec_low = elevation_class - 1
          ec_high = elevation_class + 1
       end if

       where(abs(this%topo(:, ec_high) - this%topo(:, ec_low)) < topo_equality_tolerance)
          vertical_gradient = 0._r8
       elsewhere
          vertical_gradient = &
               (this%field(:, ec_high) - this%field(:, ec_low)) / &
               (this%topo (:, ec_high) - this%topo (:, ec_low))
       end where

    end if
    
  end subroutine calc_vertical_gradient


  !-----------------------------------------------------------------------
  subroutine set_data_from_attr_vect(this, attr_vect, fieldname, toponame)
    !
    ! !DESCRIPTION:
    ! Extract data from an attribute vector.
    !
    ! Sets this%num_points, and allocates and sets this%field and this%topo.
    !
    ! !USES:
    use mct_mod
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(inout) :: this
    type(mct_aVect)  , intent(in) :: attr_vect ! attribute vector in which we can find the data
    character(len=*) , intent(in) :: fieldname ! base name of the field of interest
    character(len=*) , intent(in) :: toponame  ! base name of the topographic field
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
    character(len=:), allocatable :: elevclass_as_string
    character(len=:), allocatable :: fieldname_ec
    character(len=:), allocatable :: toponame_ec

    ! The following temporary array is needed because mct wants pointers
    real(r8), pointer :: temp(:)
    
    character(len=*), parameter :: subname = 'set_data_from_attr_vect'
    !-----------------------------------------------------------------------

    this%num_points = mct_aVect_lsize(attr_vect)

    allocate(this%field(this%num_points, this%min_elevation_class:this%max_elevation_class))
    allocate(this%topo(this%num_points, this%min_elevation_class:this%max_elevation_class))
    allocate(temp(this%num_points))
    
    do elevclass = this%min_elevation_class, this%max_elevation_class
       elevclass_as_string = glc_elevclass_as_string(elevclass)
       fieldname_ec = trim(fieldname) // elevclass_as_string
       call mct_aVect_exportRattr(attr_vect, fieldname_ec, temp)
       this%field(:,elevclass) = temp(:)

       elevclass_as_string = glc_elevclass_as_string(elevclass)
       toponame_ec = trim(toponame) // elevclass_as_string
       call mct_aVect_exportRattr(attr_vect, toponame_ec, temp)
       this%topo(:,elevclass) = temp(:)
    end do

    deallocate(temp)
    
  end subroutine set_data_from_attr_vect

end module vertical_gradient_calculator_2nd_order

