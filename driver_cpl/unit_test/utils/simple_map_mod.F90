module simple_map_mod

  ! This module defines a class for holding data describing a mapping between two grids

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod, only : errMsg => shr_log_errMsg
  
  implicit none
  private

  
  type, public :: simple_map_type
     private
     integer :: n_overlaps ! number of overlaps between the source & destination grid (size of sparse matrix)
     integer, allocatable :: source_indices(:)
     integer, allocatable :: dest_indices(:)
     real(r8), allocatable :: overlap_weights(:)
   contains
     procedure, public :: get_n_overlaps      ! get number of overlaps (size of sparse matrix)
     procedure, public :: get_n_source_points ! get number of source points
     procedure, public :: get_n_dest_points   ! get number of destination points
     procedure, public :: get_source_indices  ! get source indices in the sparse matrix
     procedure, public :: get_dest_indices    ! get dest indices in the sparse matrix
     procedure, public :: get_overlap_weights ! get overlap weights (the values in the sparse matrix)

     procedure, private :: check_okay  ! check if the data in this object are valid
     procedure, private :: check_for_duplicate_overlaps
     procedure, private :: check_for_nonpositive_weights
  end type simple_map_type

  interface simple_map_type
     module procedure constructor
  end interface simple_map_type

  ! Note: This could be written as a constructor, but instead is made a module-level
  ! routine so that it can be called with a more meaningful name.
  public :: create_simple_map_with_one_source  ! create a simple_map_type instance with a single source cell
  
contains

  ! ========================================================================
  ! Constructors and creation methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(source_indices, dest_indices, overlap_weights) result(this)
    !
    ! !DESCRIPTION:
    ! Create a simple_map_type instance.
    !
    ! The sizes of source_indices, dest_indices and overlap_weights must all be the same.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(simple_map_type) :: this  ! function result
    integer, intent(in) :: source_indices(:)
    integer, intent(in) :: dest_indices(:)
    real(r8), intent(in) :: overlap_weights(:)
    !
    ! !LOCAL VARIABLES:
    integer :: n_overlaps
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    n_overlaps = size(overlap_weights)
    call shr_assert(size(source_indices) == n_overlaps, errMsg(__FILE__, __LINE__))
    call shr_assert(size(dest_indices) == n_overlaps, errMsg(__FILE__, __LINE__))

    this%n_overlaps = n_overlaps
    this%source_indices = source_indices
    this%dest_indices = dest_indices
    this%overlap_weights = overlap_weights

    ! Perform some error-checking
    call this%check_okay()
  end function constructor

  !-----------------------------------------------------------------------
  function create_simple_map_with_one_source(ndest) result(simple_map)
    !
    ! !DESCRIPTION:
    ! Create a simple_map_type instance with a single source cell.
    !
    ! Assumes that all destination cells are fully contained within this single source cell.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(simple_map_type) :: simple_map ! function return value
    integer, intent(in) :: ndest        ! number of destination cells
    !
    ! !LOCAL VARIABLES:
    integer :: dest_index
    integer :: source_indices(ndest)
    integer :: dest_indices(ndest)
    real(r8) :: overlap_weights(ndest)
    
    character(len=*), parameter :: subname = 'create_simple_map_with_one_source'
    !-----------------------------------------------------------------------

    source_indices(:) = 1
    dest_indices = [(dest_index, dest_index = 1, ndest)]
    overlap_weights(:) = 1._r8
    simple_map = simple_map_type(source_indices=source_indices, dest_indices=dest_indices,&
         overlap_weights=overlap_weights)
    
  end function create_simple_map_with_one_source


  ! ========================================================================
  ! Class methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  function get_n_overlaps(this) result(n_overlaps)
    !
    ! !DESCRIPTION:
    ! Get number of overlaps (size of sparse matrix)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: n_overlaps  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_n_overlaps'
    !-----------------------------------------------------------------------

    n_overlaps = this%n_overlaps
  end function get_n_overlaps

  !-----------------------------------------------------------------------
  function get_n_source_points(this) result(n_source_points)
    !
    ! !DESCRIPTION:
    ! Get number of source points
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: n_source_points  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_n_source_points'
    !-----------------------------------------------------------------------

    n_source_points = maxval(this%source_indices)
  end function get_n_source_points

  !-----------------------------------------------------------------------
  function get_n_dest_points(this) result(n_dest_points)
    !
    ! !DESCRIPTION:
    ! Get number of destination points
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: n_dest_points  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_n_dest_points'
    !-----------------------------------------------------------------------

    n_dest_points = maxval(this%dest_indices)
  end function get_n_dest_points

  !-----------------------------------------------------------------------
  function get_source_indices(this) result(source_indices)
    !
    ! !DESCRIPTION:
    ! Get source indices in the sparse matrix
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, allocatable, dimension(:) :: source_indices  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_source_indices'
    !-----------------------------------------------------------------------

    source_indices = this%source_indices
  end function get_source_indices

    !-----------------------------------------------------------------------
  function get_dest_indices(this) result(dest_indices)
    !
    ! !DESCRIPTION:
    ! Get dest indices in the sparse matrix
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, allocatable, dimension(:) :: dest_indices  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_dest_indices'
    !-----------------------------------------------------------------------

    dest_indices = this%dest_indices
  end function get_dest_indices

  !-----------------------------------------------------------------------
  function get_overlap_weights(this) result(overlap_weights)
    !
    ! !DESCRIPTION:
    ! Get overlap weights (the values in the sparse matrix)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), allocatable, dimension(:) :: overlap_weights  ! function result
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_overlap_weights'
    !-----------------------------------------------------------------------

    overlap_weights = this%overlap_weights
  end function get_overlap_weights

  !-----------------------------------------------------------------------
  subroutine check_okay(this)
    !
    ! !DESCRIPTION:
    ! Makes sure that the data in this object are valid
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'check_okay'
    !-----------------------------------------------------------------------

    call this%check_for_duplicate_overlaps()
    call this%check_for_nonpositive_weights()
  end subroutine check_okay


  
  !-----------------------------------------------------------------------
  subroutine check_for_duplicate_overlaps(this)
    !
    ! !DESCRIPTION:
    ! Confirms that there are not multiple overlaps with the same source and destination
    ! indices.
    !
    ! Aborts if any duplicates are found.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    logical, allocatable :: overlap_found(:,:)
    integer :: overlap_index
    integer :: source_index
    integer :: dest_index
    
    character(len=*), parameter :: subname = 'check_for_duplicate_overlaps'
    !-----------------------------------------------------------------------

    allocate(overlap_found(this%get_n_source_points(), this%get_n_dest_points()), source=.false.)

    do overlap_index = 1, this%get_n_overlaps()
       source_index = this%source_indices(overlap_index)
       dest_index = this%dest_indices(overlap_index)
       if (overlap_found(source_index, dest_index)) then
          print *, subname, ' ERROR: duplicate found at: ', overlap_index, source_index, &
               dest_index
          stop
       end if
       overlap_found(source_index, dest_index) = .true.
    end do
       
  end subroutine check_for_duplicate_overlaps

  !-----------------------------------------------------------------------
  subroutine check_for_nonpositive_weights(this)
    !
    ! !DESCRIPTION:
    ! Confirms that all weights are positive.
    !
    ! Aborts if any zero or negative weights are found.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(simple_map_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: overlap_index
    
    character(len=*), parameter :: subname = 'check_for_nonpositive_weights'
    !-----------------------------------------------------------------------

    do overlap_index = 1, this%get_n_overlaps()
       if (this%overlap_weights(overlap_index) <= 0) then
          print *, subname, ' ERROR: non-positive weight found at: ', overlap_index, &
               this%overlap_weights(overlap_index)
          stop
       end if
    end do
  end subroutine check_for_nonpositive_weights

  
end module simple_map_mod
