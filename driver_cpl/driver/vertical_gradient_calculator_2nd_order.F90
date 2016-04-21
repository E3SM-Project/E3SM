module vertical_gradient_calculator_2nd_order

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module defines a subclass of vertical_gradient_calculator_base_type for
  ! computing vertical gradients using a second-order centered difference.

#include "shr_assert.h"
  use seq_comm_mct, only : logunit
  use vertical_gradient_calculator_base, only : vertical_gradient_calculator_base_type
  use shr_kind_mod, only : r8 => shr_kind_r8
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

     ! Bounds of each elevation class. This array has one more element than the number of
     ! elevation classes, since it contains lower and upper bounds for each elevation
     ! class. The indices go (min_elevation_class-1):max_elevation_class
     real(r8), allocatable :: elevclass_bounds(:)

   contains
     procedure :: calc_vertical_gradient
     
     procedure, private :: set_data_from_attr_vect ! extract data from an attribute vector
     procedure, private :: check_topo ! check topographic heights
     procedure, private :: limit_gradient

  end type vertical_gradient_calculator_2nd_order_type

  interface vertical_gradient_calculator_2nd_order_type
     module procedure constructor
  end interface vertical_gradient_calculator_2nd_order_type
  
contains

  !-----------------------------------------------------------------------
  function constructor(attr_vect, fieldname, toponame, &
       min_elevation_class, max_elevation_class, elevclass_names, &
       elevclass_bounds) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates a vertical_gradient_calculator_2nd_order_type object by reading the
    ! necessary data from the provided attribute vector.
    !
    ! Pre-condition: Topographic heights in the attribute vector must all lie inside the
    ! bounds of their respective elevation class (given by elevclass_bounds), with the
    ! possible exception of the lowest elevation class (topographic heights can lie below
    ! the arbitrary lower bound of the elevation class) and the highest elevation class
    ! (topographic heights can lie above the arbitrary upper bound of the elevation
    ! class). (This pre-condition is mainly important for the sake of calculating the
    ! limiter.)
    !
    ! The attribute vector is assumed to have fields named fieldname //
    ! elevclass_names(1), toponame // elevclass_names(1), etc.
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

    ! strings corresponding to each elevation class
    character(len=*) , intent(in) :: elevclass_names( min_elevation_class: )

    ! bounds of each elevation class; this array should have one more element than the
    ! number of elevation classes, since it contains lower and upper bounds for each
    ! elevation class
    real(r8)         , intent(in) :: elevclass_bounds( min_elevation_class-1 : )
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(elevclass_names) == (/max_elevation_class/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(elevclass_bounds) == (/max_elevation_class/)), errMsg(__FILE__, __LINE__))

    this%min_elevation_class = min_elevation_class
    this%max_elevation_class = max_elevation_class
    allocate(this%elevclass_bounds((min_elevation_class-1):max_elevation_class))
    this%elevclass_bounds(:) = elevclass_bounds(:)
    call this%set_data_from_attr_vect(attr_vect, fieldname, toponame, elevclass_names)

    call this%check_topo()

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

    integer :: i
    integer :: ec_low   ! elevation class index to use as the lower bound of the gradient
    integer :: ec_high  ! elevation class index to use as the upper bound of the gradient
    logical :: two_sided ! true if we're estimating the gradient with a two-sided difference

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

    ! Start by assuming we're doing a two-sided difference; we'll set this to false if we aren't
    two_sided = .true.

    if (this%min_elevation_class == this%max_elevation_class) then
       vertical_gradient(:) = 0._r8
       two_sided = .false.

    else
       
       if (elevation_class == this%min_elevation_class) then
          ec_low = elevation_class
          ec_high = elevation_class + 1
          two_sided = .false.
       else if (elevation_class == this%max_elevation_class) then
          ec_low = elevation_class - 1
          ec_high = elevation_class
          two_sided = .false.
       else
          ec_low = elevation_class - 1
          ec_high = elevation_class + 1
       end if

       do i = 1, this%num_points
          if (abs(this%topo(i, ec_high) - this%topo(i, ec_low)) < topo_equality_tolerance) then
             vertical_gradient(i) = 0._r8
          else
             vertical_gradient(i) = &
                  (this%field(i, ec_high) - this%field(i, ec_low)) / &
                  (this%topo (i, ec_high) - this%topo (i, ec_low))
          end if
       end do

       if (two_sided) then
          call this%limit_gradient(elevation_class, ec_low, ec_high, vertical_gradient)
       end if

    end if
    
  end subroutine calc_vertical_gradient


  !-----------------------------------------------------------------------
  subroutine set_data_from_attr_vect(this, attr_vect, fieldname, toponame, elevclass_names)
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
    character(len=*) , intent(in) :: elevclass_names( this%min_elevation_class: ) ! strings corresponding to each elevation class
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
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
       fieldname_ec = trim(fieldname) // trim(elevclass_names(elevclass))
       call mct_aVect_exportRattr(attr_vect, fieldname_ec, temp)
       this%field(:,elevclass) = temp(:)

       toponame_ec = trim(toponame) // trim(elevclass_names(elevclass))
       call mct_aVect_exportRattr(attr_vect, toponame_ec, temp)
       this%topo(:,elevclass) = temp(:)
    end do

    deallocate(temp)
    
  end subroutine set_data_from_attr_vect

  !-----------------------------------------------------------------------
  subroutine check_topo(this)
    !
    ! !DESCRIPTION:
    ! Check topographic heights; abort if there is a problem
    !
    ! Topographic heights in the attribute vector must all lie inside the bounds of their
    ! respective elevation class (given by elevclass_bounds), with the possible exception
    ! of the lowest elevation class (topographic heights can lie below the arbitrary lower
    ! bound of the elevation class) and the highest elevation class (topographic heights
    ! can lie above the arbitrary upper bound of the elevation class)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: elevclass
    integer :: i

    ! Absolute tolerance for error checks. This is chosen so that it allows for
    ! double-precision roundoff-level errors on values of order 10,000.
    real(r8), parameter :: tol = 1.e-10_r8

    character(len=*), parameter :: subname = 'check_topo'
    !-----------------------------------------------------------------------

    do elevclass = this%min_elevation_class, this%max_elevation_class
       if (elevclass > this%min_elevation_class) then
          do i = 1, this%num_points
             if (this%topo(i,elevclass) - this%elevclass_bounds(elevclass-1) < -tol) then
                write(logunit,*) subname, ': ERROR: topo lower than lower bound of elevation class:'
                write(logunit,*) 'i, elevclass, topo, lower_bound = ', &
                     i, elevclass, this%topo(i,elevclass), this%elevclass_bounds(elevclass-1)
                call shr_sys_abort(subname//': ERROR: topo lower than lower bound of elevation class')
             end if
          end do
       end if

       if (elevclass < this%max_elevation_class) then
          do i = 1, this%num_points
             if (this%topo(i,elevclass) - this%elevclass_bounds(elevclass) > tol) then
                write(logunit,*) subname, ': ERROR: topo higher than upper bound of elevation class:'
                write(logunit,*) 'i, elevclass, topo, upper_bound = ', &
                     i, elevclass, this%topo(i,elevclass), this%elevclass_bounds(elevclass)
                call shr_sys_abort(subname//': ERROR: topo higher than upper bound of elevation class')
             end if
          end do
       end if
    end do

  end subroutine check_topo

  !-----------------------------------------------------------------------
  subroutine limit_gradient(this, k, ec_low, ec_high, vertical_gradient)
    !
    ! !DESCRIPTION:
    ! Limit the gradient: Ensure that the interface values lie inside the range defined
    ! by the max and min of the mean values in this class and its 2 adjacent neighbors.
    !
    ! Should only be called for two-sided differences (ec_low < k < ec_high)
    !
    ! !ARGUMENTS:
    class(vertical_gradient_calculator_2nd_order_type), intent(in) :: this
    integer , intent(in) :: k       ! elevation class index
    integer , intent(in) :: ec_low  ! elevation class index used as the lower bound of the gradient
    integer , intent(in) :: ec_high ! elevation class index used as the upper bound of the gradient
    real(r8), intent(inout) :: vertical_gradient(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    real(r8) :: deviation_high
    real(r8) :: deviation_low
    real(r8) :: deviation_max
    real(r8) :: deviation_min
    real(r8) :: diff_max
    real(r8) :: diff_min
    real(r8) :: factor1
    real(r8) :: factor2
    real(r8) :: limiting_factor

    character(len=*), parameter :: subname = 'limit_gradient'
    !-----------------------------------------------------------------------

    ! Basic idea: In 1D with a linear reconstruction, the extreme values of the data will
    ! lie at the interfaces between adjacent elevation classes. The interface values
    ! should not lie outside the range defined by the max and min of the mean values in
    ! this class and its 2 adjacent neighbors.

    ! This code only works correctly if we're doing a two-sided difference (otherwise,
    ! one of diff_min or diff_max will be 0, leading to 0 gradient - when in fact we
    ! don't want to do any limiting for a one-sided difference).
    SHR_ASSERT(ec_low < k, subname//': Only works for two-sided difference: must have ec_low < k')
    SHR_ASSERT(ec_high > k, subname//': Only works for two-sided difference: must have ec_high > k')

    do i = 1, this%num_points
       ! First compute the max and min values of the deviation of the data from its mean
       ! value. With a linear gradient, the max differences must lie at the adjacent
       ! interfaces.
       deviation_high = vertical_gradient(i) * (this%elevclass_bounds(k) - this%topo(i,k))
       deviation_low = vertical_gradient(i) *  (this%elevclass_bounds(k-1) - this%topo(i,k))
       deviation_max = max(deviation_high, deviation_low)
       deviation_min = min(deviation_high, deviation_low)

       ! Now compute the max and min of the data in the cell and its nearest neighbors.
       ! (Actually, the difference between this max/min value and the mean value in the
       ! current class.)
       diff_max = max(this%field(i,ec_high), this%field(i,k), this%field(i,ec_low)) - this%field(i,k)
       diff_min = min(this%field(i,ec_high), this%field(i,k), this%field(i,ec_low)) - this%field(i,k)

       ! Now limit the gradient using the information computed above.

       if (abs(deviation_min) > 0._r8) then
          factor1 = max(0._r8, diff_min/deviation_min)
       else
          factor1 = 1._r8
       endif
       
       if (abs(deviation_max) > 0._r8) then
          factor2 = max(0._r8, diff_max/deviation_max)
       else
          factor2 = 1._r8
       endif

       ! limiting factor will lie between 0 and 1
       limiting_factor = min(1._r8, factor1, factor2)
       vertical_gradient(i) = vertical_gradient(i) * limiting_factor
    end do

  end subroutine limit_gradient


end module vertical_gradient_calculator_2nd_order

