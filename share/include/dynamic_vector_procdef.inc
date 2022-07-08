! Type-bound procedures for a dynamic vector.

#ifdef USE_PURE
#define PURE pure
#else
#define PURE
#endif

! Construct an empty vector.
PURE function new_vector_default() result(new_vec)
  ! Create an empty vector
  type( VECTOR_NAME ) :: new_vec

  ! Currently, this does nothing. But some compilers may do weird things if
  ! you don't "define" new_vec somehow, and clearing the vector is safe.
  call new_vec%clear()

end function new_vector_default

! Construct a vector from another vector.
PURE function new_vector_copy(vec) result(new_vec)
  ! Create a vector from a pre-existing array.
  type( VECTOR_NAME ), intent(in) :: vec
  type( VECTOR_NAME ) :: new_vec

  new_vec = vec

end function new_vector_copy

! Construct a vector from an array.
PURE function new_vector_array(array) result(new_vec)
  ! Create a vector from a pre-existing array.
  TYPE_NAME, intent(in) :: array(:)
  type( VECTOR_NAME ) :: new_vec

  new_vec = array

end function new_vector_array

! Query if the vector is empty.
PURE function empty_vec(self) result(is_empty)
  class( VECTOR_NAME ), intent(in) :: self
  logical :: is_empty

  is_empty = (self%vec_size == 0)

end function empty_vec

! Get size of the vector.
PURE function size_vec(self) result(vec_size)
  class( VECTOR_NAME ), intent(in) :: self
  integer :: vec_size

  vec_size = self%vec_size

end function size_vec

! Get maximum size the vector can have.
PURE function max_size_vec(self) result(max_size)
  class( VECTOR_NAME ), intent(in) :: self
  integer :: max_size

  ! The only theoretical limitation that can be determined without a system
  ! call is the maximum size of an integer.
  max_size = huge(self%vec_size)

end function max_size_vec

! Query current memory capacity of vector.
PURE function capacity_vec(self) result(capacity)
  class( VECTOR_NAME ), intent(in) :: self
  integer :: capacity

  if (allocated(self%data)) then
     capacity = size(self%data)
  else
     capacity = 0
  end if

end function capacity_vec

! Get one item based on an index.
PURE function get_single_vec(self, index) result(item)
  class( VECTOR_NAME ), intent(in) :: self
  integer, intent(in) :: index
  TYPE_NAME, allocatable :: item

  if (index > self%vec_size .or. index < 1) then
     THROW(OOBMsg("get", [1, self%vec_size], index))
     ! Purely to satisfy uninitialized data checks.
     allocate(item)
     return
  end if

  allocate(item, source=self%data(index))

end function get_single_vec

! Get items within a certain range.
PURE function get_range_vec(self, begin, end, stride) result(items)
  class( VECTOR_NAME ), intent(in) :: self
  integer, intent(in) :: begin
  integer, intent(in) :: end
  integer, intent(in), optional :: stride

  ! Have to use an allocatable, because we have to check if stride is
  ! present before we know what the size should be.
  TYPE_NAME, allocatable :: items(:)

  ! An oddity: since in Fortran function results must be "defined", we have
  ! to allocate "items" to portably avoid a segfault and allow the user to
  ! recover from an error. This is true regardless of what the function
  ! result is assigned to.
  if (end > self%vec_size) then
     allocate(items(0))
     THROW(OOBMsg("get", [1, self%vec_size], end))
     return
  end if
  if (begin < 1) then
     allocate(items(0))
     THROW(OOBMsg("get", [1, self%vec_size], begin))
     return
  end if

  if (present(stride)) then
     ! For strided access, the number of elements is the size over the
     ! stride, but rounded up, rather than down as in typical integer
     ! division. The shortcut for this is that (x-1)/y + 1 is the same as
     ! x/y rounded up.
     allocate(items((end-begin)/stride + 1))
     items = self%data(begin:end:stride)
  else
     allocate(items(end+1-begin))
     items = self%data(begin:end)
  end if

end function get_range_vec

! Get an array containing a copy of the vector's elements.
! If array is not allocated, returns a size zero array.
PURE function get_array_vec(self) result(array)
  class( VECTOR_NAME ), intent(in) :: self
  TYPE_NAME :: array(self%vec_size)

  if (allocated(self%data)) then
     array = self%data(:self%vec_size)
  end if

end function get_array_vec

! Get first item in the array
PURE function front_vec(self) result(item)
  class( VECTOR_NAME ), intent(in) :: self
  TYPE_NAME :: item

  item = self%get(1)

end function front_vec

! Get last item in the array
PURE function back_vec(self) result(item)
  class( VECTOR_NAME ), intent(in) :: self
  TYPE_NAME :: item

  item = self%get(self%vec_size)

end function back_vec

! Declare the vector to have zero size.
! Does not change vector capacity.
PURE subroutine clear_vec(self)
  class( VECTOR_NAME ), intent(inout) :: self

  call self%resize(0)

end subroutine clear_vec

! Declare the vector to have different size.
! Does not reduce vector capacity, but will enforce size <= capacity by
! growing array if necessary.
! Resizing to negative value is equivalent to resizing to 0.
PURE subroutine resize_vec(self, new_size, fill_value)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: new_size
  TYPE_NAME, intent(in), optional :: fill_value

  integer :: request_capacity

  ! If not big enough, request capacity twice as big
  ! as we have now (or 4 or 8 or... times, if necessary).
  if (new_size > self%capacity()) then
     request_capacity = max(self%capacity(),1)

     do while (request_capacity < new_size)
        request_capacity = request_capacity * 2
     end do

     call self%reserve(request_capacity)
  end if

  if (present(fill_value)) then
     self%data((self%vec_size+1):new_size) = fill_value
  end if

  self%vec_size = max(new_size,0)

end subroutine resize_vec

! Set one item based on an index.
PURE subroutine set_single_vec(self, item, index)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: item
  integer, intent(in) :: index

  if (index > self%vec_size .or. index < 1) then
     THROW(OOBMsg("set", [1, self%vec_size], index))
     return
  end if

  self%data(index) = item

end subroutine set_single_vec

! Set range in array.
PURE subroutine set_range_vec(self, array, begin, end, stride)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: array(:)
  integer, intent(in) :: begin
  integer, intent(in) :: end
  integer, intent(in), optional :: stride

  if (end > self%vec_size) then
     THROW(OOBMsg("set", [1, self%vec_size], end))
     return
  end if
  if (begin < 1) then
     THROW(OOBMsg("set", [1, self%vec_size], begin))
     return
  end if

  if (present(stride)) then
     self%data(begin:end:stride) = array
  else
     self%data(begin:end) = array
  end if

end subroutine set_range_vec

! Set range in array with a fill value.
PURE subroutine set_range_fill_vec(self, fill_value, begin, end, stride)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: fill_value
  integer, intent(in) :: begin
  integer, intent(in) :: end
  integer, intent(in), optional :: stride

  if (end > self%vec_size) then
     THROW(OOBMsg("set", [1, self%vec_size], end))
     return
  end if
  if (begin < 1) then
     THROW(OOBMsg("set", [1, self%vec_size], begin))
     return
  end if

  if (present(stride)) then
     self%data(begin:end:stride) = fill_value
  else
     self%data(begin:end) = fill_value
  end if

end subroutine set_range_fill_vec

! Set array from an array.
PURE subroutine set_array_vec(self, array)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: array(:)

  if (size(array) /= self%vec_size) then
     THROW("Input array is not the same size as the vector it sets.")
  end if

  if (self%vec_size > 0) then
     self%data(:self%vec_size) = array(:self%vec_size)
  end if

end subroutine set_array_vec

! Set array from a fill value.
! Bounds-checking unnecessary; empty arrays are left empty.
PURE subroutine set_fill_vec(self, fill_value)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: fill_value

  if (allocated(self%data)) then
     self%data(:self%vec_size) = fill_value
  end if

end subroutine set_fill_vec

! Add new object as last element.
PURE subroutine push_back_vec(self, item)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: item

  call self%resize(self%vec_size+1)

  call self%set(item, self%vec_size)

end subroutine push_back_vec

! Remove last element.
PURE subroutine pop_back_vec(self)
  class( VECTOR_NAME ), intent(inout) :: self

  if (self%empty()) then
     THROW("Attempted to pop an element from an empty vector.")
  end if

  call self%resize(self%vec_size-1)

end subroutine pop_back_vec

! Insert element
! Valid values are 1 to self%vec_size+1.
! Inserting at self%vec_size+1 is equivalent to push_back.
PURE subroutine insert_single_vec(self, index, item)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: index
  TYPE_NAME, intent(in) :: item

  if (index > self%vec_size+1 .or. index < 1) then
     THROW(OOBMsg("insert", [1, self%vec_size], index))
     return
  end if

  call self%resize(self%vec_size+1)

  ! Move everything forward
  self%data(index+1:self%vec_size) = &
       self%data(index:self%vec_size-1)

  call self%set(item, index)

end subroutine insert_single_vec

! Insert array
PURE subroutine insert_array_vec(self, index, items)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: index
  TYPE_NAME, intent(in) :: items(:)

  if (index > self%vec_size+1 .or. index < 1) then
     THROW(OOBMsg("insert", [1, self%vec_size], index))
     return
  end if

  call self%resize(self%vec_size+size(items))

  ! Move everything forward
  self%data(index+size(items):self%vec_size) = &
       self%data(index:self%vec_size-size(items))

  call self%set(items, index, index+size(items)-1)

end subroutine insert_array_vec

! Insert repeated value
PURE subroutine insert_repeat_vec(self, index, item, repeats)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: index
  TYPE_NAME, intent(in) :: item
  integer, intent(in) :: repeats

  if (index > self%vec_size+1 .or. index < 1) then
     THROW(OOBMsg("insert", [1, self%vec_size], index))
     return
  end if

  call self%resize(self%vec_size+repeats)

  ! Move everything forward
  self%data(index+repeats:self%vec_size) = &
       self%data(index:self%vec_size-repeats)

  call self%set(item, index, index+repeats-1)

end subroutine insert_repeat_vec

! Erase element
PURE subroutine erase_single_vec(self, index)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: index

  if (index > self%vec_size .or. index < 1) then
     THROW(OOBMsg("erase", [1, self%vec_size], index))
     return
  end if

  ! Move everything back
  self%data(index:(self%vec_size-1)) = self%data((index+1):self%vec_size)

  call self%pop_back()

end subroutine erase_single_vec

! Erase "repeats" elements at index.
PURE subroutine erase_range_vec(self, begin, end)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: begin
  integer, intent(in) :: end

  if (end > self%vec_size) then
     THROW(OOBMsg("erase", [1, self%vec_size], end))
     return
  end if
  if (begin < 1) then
     THROW(OOBMsg("erase", [1, self%vec_size], begin))
     return
  end if

  ! Move everything back
  self%data(begin:self%vec_size-end+begin-1) = &
       self%data(end+1:self%vec_size)

  call self%resize(self%vec_size - end + begin-1)

end subroutine erase_range_vec

! Shrink vector to minimum size necessary to hold all elements.
PURE subroutine shrink_to_fit_vec(self)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, allocatable :: tmp_array(:)

  ! Don't do anything unless we have to.
  if (self%vec_size < self%capacity()) then
     ! If size is zero, just deallocate array.
     if (self%vec_size == 0) then
        if (allocated(self%data)) deallocate(self%data)
     else
        ! Allocate temporary at minimum size
        allocate(tmp_array(self%vec_size))
        tmp_array = self%data(:self%vec_size)

        deallocate(self%data)
        call move_alloc(tmp_array, self%data)
     end if
  end if

end subroutine shrink_to_fit_vec

! Reserve a certain size, if vector is not already that big.
PURE subroutine reserve_vec(self, capacity)
  class( VECTOR_NAME ), intent(inout) :: self
  integer, intent(in) :: capacity

  TYPE_NAME, allocatable :: tmp_array(:)

  ! Only do anything if we need to get bigger.
  if (capacity > self%capacity()) then

     if (self%empty()) then
        ! No data to copy
        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(capacity))
     else
        ! Allocate new size
        allocate(tmp_array(capacity))
        ! Copy data
        tmp_array(:self%vec_size) = self%data(:self%vec_size)

        ! Replace array with new copy.
        deallocate(self%data)
        call move_alloc(tmp_array, self%data)
     end if
  end if

end subroutine reserve_vec

! Move allocatable array into self
! Note: Declaring self as intent(out) automatically empties the vector the
!       moment we enter this procedure!
PURE subroutine move_in_vec(self, array)
  class( VECTOR_NAME ), intent(out) :: self
  TYPE_NAME, allocatable, intent(inout) :: array(:)

  if (allocated(array)) then
     call move_alloc(array, self%data)
     self%vec_size = size(self%data)
  end if

end subroutine move_in_vec

! Move self into output allocatable array.
! For empty vector, do not allocate output.
PURE subroutine move_out_vec(self, array)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, allocatable, intent(out) :: array(:)

  call self%shrink_to_fit()

  if (.not. self%empty()) then
     call move_alloc(self%data, array)
  end if

  call self%clear()

end subroutine move_out_vec

! Efficient swapping (no de/reallocation)
PURE subroutine swap_vec(self, other)
  class( VECTOR_NAME ), intent(inout) :: self
  class( VECTOR_NAME ), intent(inout) :: other

  integer :: tmp_size
  TYPE_NAME, allocatable :: tmp_array(:)

  ! The following order is designed to work even if self and other are the
  ! same vector.
  if (allocated(other%data)) then
     call move_alloc(other%data, tmp_array)
  end if

  if (allocated(self%data)) then
#ifndef CPRPGI
     call move_alloc(self%data, other%data)
#else
     ! The above should work, but a PGI bug forces us to copy and
     ! deallocate.
     allocate(other%data, source=self%data)
     deallocate(self%data)
#endif
  end if

  if (allocated(tmp_array)) then
#ifndef CPRPGI
     call move_alloc(tmp_array, self%data)
#else
     ! The above should work, but a PGI bug forces us to copy and
     ! deallocate.
     allocate(self%data, source=tmp_array)
     deallocate(tmp_array)
#endif
  end if

  tmp_size = other%vec_size
  other%vec_size = self%vec_size
  self%vec_size = tmp_size

end subroutine swap_vec

! Assign self from an array
PURE subroutine array_assign_vec(self, array)
  class( VECTOR_NAME ), intent(inout) :: self
  TYPE_NAME, intent(in) :: array(:)

  call self%resize(size(array))

  call self%set(array)

end subroutine array_assign_vec

! Assign self from another vector.
! Copy-and-swap is used to ensure that at most one copy of the array is
! performed.
! This would allow assignment to self in other languages, but Fortran 2003
! is vague about whether this should work, since "other" must be
! "intent(in)" for an assignment, and this routine would modify it if it is
! the same as "self".
! Use of the "target" attribute is intended to mitigate the risk of a
! problem, warning the compiler that the two objects may overlap with other
! variables.
PURE subroutine vector_assign_vec(self, other)
  class( VECTOR_NAME ), intent(inout), target :: self
  class( VECTOR_NAME ), intent(in), target :: other

  class( VECTOR_NAME ), allocatable :: temp

  allocate(temp, source=other)

  call self%swap(temp)

  deallocate(temp)

end subroutine vector_assign_vec

#undef PURE
