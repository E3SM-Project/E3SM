! Test template for dynamic vector template.
!
! Include in module after "contains" with these prerequisites:
!  - Use the vector from the module where it is defined.
!  - Use shr_log_OOBMsg from shr_log_mod.
!  - Define TYPE_NAME as the name of the type held by the vector.
!  - Define TYPE_DECL as a declaration of the type held by the vector (e.g.
!    "real(r8)" or "type(foo)".
!  - Define VECTOR_NAME as the name of the vector type.
!  - Define test arrays global to the module:
!     - One named "test_array" with 3 elements
!     - One named "test_array_2" with 5 elements
!  - Define a scalar value of the type held by the vector,
!    named "new_val".
!

! We want to prefix with the type name so that it's possible to tell which test
! failed when pFUnit reports the test routine's name.
! This concatenation likely requires a real cpp for concatenation, rather than
! an fpp that may not implement "##".
#define CONCAT_UNDERSCORE(x, y) x ## _ ## y
#define CONCAT_EXPAND_UNDERSCORE(x, y) CONCAT_UNDERSCORE(x, y)
#define ADD_PREFIX(name) CONCAT_EXPAND_UNDERSCORE(TYPE_NAME, name)

@Test
subroutine ADD_PREFIX(empty)()

  type( VECTOR_NAME ) :: vec

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(empty)

@Test
subroutine ADD_PREFIX(alloc_empty)()

  type( VECTOR_NAME ) :: vec

  call vec%reserve(5)

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(alloc_empty)

@Test
subroutine ADD_PREFIX(non_empty)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  @assertFalse(vec%empty())

end subroutine ADD_PREFIX(non_empty)

@Test
subroutine ADD_PREFIX(default_constructor)()

  type( VECTOR_NAME ) :: vec

  vec = VECTOR_NAME ()

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(default_constructor)

@Test
subroutine ADD_PREFIX(assign_empty_array)()

  type( VECTOR_NAME ) :: vec
  TYPE_DECL :: new_array(0)

  vec = new_array

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(assign_empty_array)

@Test
subroutine ADD_PREFIX(size)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  @assertEqual(size(test_array), vec%vsize())

end subroutine ADD_PREFIX(size)

@Test
subroutine ADD_PREFIX(size_empty)()

  type( VECTOR_NAME ) :: vec

  @assertEqual(0, vec%vsize())

end subroutine ADD_PREFIX(size_empty)

@Test
subroutine ADD_PREFIX(max_size)()

  type( VECTOR_NAME ) :: vec

  ! No idea how you could verify this. Just
  ! make sure it's more than 0.
  @assertLessThan(0, vec%max_size())

end subroutine ADD_PREFIX(max_size)

@Test
subroutine ADD_PREFIX(capacity)()

  type( VECTOR_NAME ) :: vec

  call vec%reserve(5)

  @assertEqual(5, vec%capacity())

end subroutine ADD_PREFIX(capacity)

@Test
subroutine ADD_PREFIX(empty_capacity)()

  type( VECTOR_NAME ) :: vec

  @assertEqual(0, vec%capacity())

end subroutine ADD_PREFIX(empty_capacity)

@Test
subroutine ADD_PREFIX(get)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL :: compare_array(size(test_array))

  integer :: i

  vec = test_array

  do i = 1, size(test_array)
     compare_array(i) = vec%get(i)
  end do

  @assertEqual(test_array, compare_array)

end subroutine ADD_PREFIX(get)

@Test
subroutine ADD_PREFIX(get_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL :: throw_away

  integer :: i

  vec = test_array

  throw_away = vec%get(0)

  call assertExceptionRaised(OOBMsg("get", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(get_bnd_chk_low)

@Test
subroutine ADD_PREFIX(get_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL :: throw_away

  integer :: i

  vec = test_array

  throw_away = vec%get(size(test_array)+1)

  call assertExceptionRaised(OOBMsg("get", [1, vec%vsize()], size(test_array)+1))

end subroutine ADD_PREFIX(get_bnd_chk_high)

@Test
subroutine ADD_PREFIX(get_range)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 2
  integer, parameter :: end = 4

  vec = test_array_2

  @assertEqual(test_array_2(begin:end), vec%get(begin, end))

end subroutine ADD_PREFIX(get_range)

@Test
subroutine ADD_PREFIX(get_range_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: throw_away(:)

  vec = test_array

  throw_away = vec%get(0, size(test_array))

  call assertExceptionRaised(OOBMsg("get", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(get_range_bnd_chk_low)

@Test
subroutine ADD_PREFIX(get_range_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: throw_away(:)

  vec = test_array

  throw_away = vec%get(1, size(test_array)+1)

  call assertExceptionRaised(OOBMsg("get", [1, vec%vsize()], size(test_array)+1))

end subroutine ADD_PREFIX(get_range_bnd_chk_high)

@Test
subroutine ADD_PREFIX(get_range_stride)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 1
  integer, parameter :: end = 5
  integer, parameter :: stride = 2

  vec = test_array_2

  @assertEqual(test_array_2(begin:end:stride), vec%get(begin, end, stride))

end subroutine ADD_PREFIX(get_range_stride)

@Test
subroutine ADD_PREFIX(get_array)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  @assertEqual(test_array, vec%get())

end subroutine ADD_PREFIX(get_array)

@Test
subroutine ADD_PREFIX(array_constructor)()

  type( VECTOR_NAME ) :: vec

  vec = VECTOR_NAME (test_array)

  @assertEqual(test_array, vec%get())

end subroutine ADD_PREFIX(array_constructor)

@Test
subroutine ADD_PREFIX(copy_constructor)()

  type( VECTOR_NAME ) :: vec1
  type( VECTOR_NAME ) :: vec2

  vec1 = test_array
  vec2 = VECTOR_NAME (vec1)

  @assertEqual(test_array, vec2%get())

end subroutine ADD_PREFIX(copy_constructor)

@Test
subroutine ADD_PREFIX(front)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  @assertEqual(test_array(1), vec%front())

end subroutine ADD_PREFIX(front)

@Test
subroutine ADD_PREFIX(front_bnd_chk)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL :: throw_away

  throw_away = vec%front()

  call assertExceptionRaised(OOBMsg("get", [1, 0], 1))

end subroutine ADD_PREFIX(front_bnd_chk)

@Test
subroutine ADD_PREFIX(back)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  @assertEqual(test_array(size(test_array)), vec%back())

end subroutine ADD_PREFIX(back)

@Test
subroutine ADD_PREFIX(back_bnd_chk)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL :: throw_away

  throw_away = vec%back()

  call assertExceptionRaised(OOBMsg("get", [1, 0], 0))

end subroutine ADD_PREFIX(back_bnd_chk)

@Test
subroutine ADD_PREFIX(clear)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%clear()

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(clear)

@Test
subroutine ADD_PREFIX(resize)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: new_size = 5

  vec = test_array

  call vec%resize(new_size)

  @assertEqual(new_size, vec%vsize())
  @assertLessThanOrEqual(new_size, vec%capacity())

end subroutine ADD_PREFIX(resize)

@Test
subroutine ADD_PREFIX(resize_fill)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: new_size = 5

  vec = test_array

  call vec%resize(new_size, new_val)

  @assertEqual(new_size, vec%vsize())
  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(new_val, vec%get(size(test_array)+1, vec%vsize()))

end subroutine ADD_PREFIX(resize_fill)

@Test
subroutine ADD_PREFIX(resize_negative)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: new_size = -1

  vec = test_array

  call vec%resize(new_size)

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(resize_negative)

@Test
subroutine ADD_PREFIX(set)()

  type( VECTOR_NAME ) :: vec

  integer :: i

  vec = test_array_2

  do i = 1, size(test_array)
     call vec%set(test_array(i), i)
  end do

  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(test_array_2(size(test_array)+1:), vec%get(size(test_array)+1, vec%vsize()))

end subroutine ADD_PREFIX(set)

@Test
subroutine ADD_PREFIX(set_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(new_val, 0)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(set_bnd_chk_low)

@Test
subroutine ADD_PREFIX(set_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(new_val, size(test_array)+1)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], size(test_array)+1))

end subroutine ADD_PREFIX(set_bnd_chk_high)

@Test
subroutine ADD_PREFIX(set_range)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 2
  integer, parameter :: end = 4

  TYPE_DECL :: compare_array(size(test_array_2))

  vec = test_array_2

  call vec%set(test_array, begin, end)

  compare_array = test_array_2
  compare_array(begin:end) = test_array

  @assertEqual(compare_array, vec%get())

end subroutine ADD_PREFIX(set_range)

@Test
subroutine ADD_PREFIX(set_range_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array_2

  call vec%set(test_array, 0, size(test_array)-1)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(set_range_bnd_chk_low)

@Test
subroutine ADD_PREFIX(set_range_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array_2

  call vec%set(test_array, &
       size(test_array_2)-size(test_array)+1, size(test_array_2)+1)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], size(test_array_2)+1))

end subroutine ADD_PREFIX(set_range_bnd_chk_high)

@Test
subroutine ADD_PREFIX(set_range_stride)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 1
  integer, parameter :: end = 5
  integer, parameter :: stride = 2

  TYPE_DECL :: compare_array(size(test_array_2))

  vec = test_array_2

  call vec%set(test_array, begin, end, stride)

  compare_array = test_array_2
  compare_array(begin:end:stride) = test_array

  @assertEqual(compare_array, vec%get())

end subroutine ADD_PREFIX(set_range_stride)

@Test
subroutine ADD_PREFIX(set_range_fill)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 2
  integer, parameter :: end = 4

  TYPE_DECL :: compare_array(size(test_array_2))

  vec = test_array_2

  call vec%set(new_val, begin, end)

  compare_array = test_array_2
  compare_array(begin:end) = new_val

  @assertEqual(compare_array, vec%get())

end subroutine ADD_PREFIX(set_range_fill)

@Test
subroutine ADD_PREFIX(set_range_fill_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(new_val, 0, 1)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(set_range_fill_bnd_chk_low)

@Test
subroutine ADD_PREFIX(set_range_fill_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(new_val, 1, size(test_array)+1)

  call assertExceptionRaised(OOBMsg("set", [1, vec%vsize()], size(test_array)+1))

end subroutine ADD_PREFIX(set_range_fill_bnd_chk_high)

@Test
subroutine ADD_PREFIX(set_range_stride_fill)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 1
  integer, parameter :: end = 5
  integer, parameter :: stride = 2

  TYPE_DECL :: compare_array(size(test_array_2))

  vec = test_array_2

  call vec%set(new_val, begin, end, stride)

  compare_array = test_array_2
  compare_array(begin:end:stride) = new_val

  @assertEqual(compare_array, vec%get())

end subroutine ADD_PREFIX(set_range_stride_fill)

@Test
subroutine ADD_PREFIX(set_array)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  ! Note: This does NOT resize the vector.
  call vec%set(test_array_2(:size(test_array)))

  @assertEqual(test_array_2(:size(test_array)), vec%get())

end subroutine ADD_PREFIX(set_array)

@Test
subroutine ADD_PREFIX(set_array_chk)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(test_array_2)

  call assertExceptionRaised("Input array is not the same size as the vector it sets.")

end subroutine ADD_PREFIX(set_array_chk)

@Test
subroutine ADD_PREFIX(set_fill)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%set(new_val)

  @assertEqual(new_val, vec%get())

end subroutine ADD_PREFIX(set_fill)

@Test
subroutine ADD_PREFIX(push_back)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%push_back(new_val)

  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(new_val, vec%get(size(test_array)+1))

end subroutine ADD_PREFIX(push_back)

@Test
subroutine ADD_PREFIX(pop_back)()

  type( VECTOR_NAME ) :: vec

  vec = test_array_2

  call vec%pop_back()

  @assertEqual(test_array_2(:(size(test_array_2)-1)), vec%get())

end subroutine ADD_PREFIX(pop_back)

@Test
subroutine ADD_PREFIX(pop_back_throw)()

  type( VECTOR_NAME ) :: vec

  call vec%pop_back()

  call assertExceptionRaised("Attempted to pop an element from an empty vector.")

end subroutine ADD_PREFIX(pop_back_throw)

@Test
subroutine ADD_PREFIX(insert)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: ins = 2

  vec = test_array

  call vec%insert(ins, new_val)

  @assertEqual(test_array(:ins-1), vec%get(1, ins-1))
  @assertEqual(new_val, vec%get(ins))
  @assertEqual(test_array(ins:), vec%get(ins+1, size(test_array)+1))

end subroutine ADD_PREFIX(insert)

@Test
subroutine ADD_PREFIX(insert_at_end)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(size(test_array)+1, new_val)

  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(new_val, vec%get(size(test_array)+1))

end subroutine ADD_PREFIX(insert_at_end)

@Test
subroutine ADD_PREFIX(insert_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(0, new_val)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(insert_bnd_chk_low)

@Test
subroutine ADD_PREFIX(insert_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(size(test_array)+2, new_val)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], size(test_array)+2))

end subroutine ADD_PREFIX(insert_bnd_chk_high)

@Test
subroutine ADD_PREFIX(insert_array)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: ins = 2

  vec = test_array

  call vec%insert(ins, test_array_2)

  @assertEqual(test_array(:ins-1), vec%get(1, ins-1))
  @assertEqual(test_array_2, vec%get(ins, ins-1+size(test_array_2)))
  @assertEqual(test_array(ins:), vec%get(ins+size(test_array_2), size(test_array)+size(test_array_2)))

end subroutine ADD_PREFIX(insert_array)

@Test
subroutine ADD_PREFIX(insert_array_at_end)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(size(test_array)+1, test_array_2)

  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(test_array_2, vec%get(size(test_array)+1, size(test_array)+size(test_array_2)))

end subroutine ADD_PREFIX(insert_array_at_end)

@Test
subroutine ADD_PREFIX(insert_array_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(0, test_array_2)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(insert_array_bnd_chk_low)

@Test
subroutine ADD_PREFIX(insert_array_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(size(test_array)+2, test_array_2)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], size(test_array)+2))

end subroutine ADD_PREFIX(insert_array_bnd_chk_high)

@Test
subroutine ADD_PREFIX(insert_repeat)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: ins = 2

  integer, parameter :: repeat_num = 3

  vec = test_array

  call vec%insert(ins, new_val, repeat_num)

  @assertEqual(test_array(:ins-1), vec%get(1, ins-1))
  @assertEqual(new_val, vec%get(ins, ins-1+repeat_num))
  @assertEqual(test_array(ins:), vec%get(ins+repeat_num, size(test_array)+repeat_num))

end subroutine ADD_PREFIX(insert_repeat)

@Test
subroutine ADD_PREFIX(insert_repeat_at_end)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: repeat_num = 3

  vec = test_array

  call vec%insert(size(test_array)+1, new_val, repeat_num)

  @assertEqual(test_array, vec%get(1, size(test_array)))
  @assertEqual(new_val, vec%get(size(test_array)+1, size(test_array)+repeat_num))

end subroutine ADD_PREFIX(insert_repeat_at_end)

@Test
subroutine ADD_PREFIX(insert_repeat_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(0, new_val, 5)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(insert_repeat_bnd_chk_low)

@Test
subroutine ADD_PREFIX(insert_repeat_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%insert(size(test_array)+2, new_val, 5)

  call assertExceptionRaised(OOBMsg("insert", [1, vec%vsize()], size(test_array)+2))

end subroutine ADD_PREFIX(insert_repeat_bnd_chk_high)

@Test
subroutine ADD_PREFIX(erase)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: pos = 2

  vec = test_array

  call vec%erase(pos)

  @assertEqual(test_array(:pos-1), vec%get(1, pos-1))
  @assertEqual(test_array(pos+1:), vec%get(pos, size(test_array)-1))

end subroutine ADD_PREFIX(erase)

@Test
subroutine ADD_PREFIX(erase_at_end)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%erase(size(test_array))

  @assertEqual(test_array(:size(test_array)-1), vec%get())

end subroutine ADD_PREFIX(erase_at_end)

@Test
subroutine ADD_PREFIX(erase_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%erase(0)

  call assertExceptionRaised(OOBMsg("erase", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(erase_bnd_chk_low)

@Test
subroutine ADD_PREFIX(erase_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%erase(size(test_array)+1)

  call assertExceptionRaised(OOBMsg("erase", [1, vec%vsize()], size(test_array)+1))

end subroutine ADD_PREFIX(erase_bnd_chk_high)

@Test
subroutine ADD_PREFIX(erase_range)()

  type( VECTOR_NAME ) :: vec

  integer, parameter :: begin = 2
  integer, parameter :: end = 3

  vec = test_array_2

  call vec%erase(begin, end)

  @assertEqual(test_array_2(:begin-1), vec%get(1, begin-1))
  @assertEqual(test_array_2(end+1:), vec%get(begin, size(test_array_2)-end+begin-1))

end subroutine ADD_PREFIX(erase_range)


@Test
subroutine ADD_PREFIX(erase_range_at_end)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%erase(size(test_array)-1, size(test_array))

  @assertEqual(test_array(:size(test_array)-2), vec%get())

end subroutine ADD_PREFIX(erase_range_at_end)
@Test
subroutine ADD_PREFIX(erase_range_bnd_chk_low)()

  type( VECTOR_NAME ) :: vec

  vec = test_array_2

  call vec%erase(0, 1)

  call assertExceptionRaised(OOBMsg("erase", [1, vec%vsize()], 0))

end subroutine ADD_PREFIX(erase_range_bnd_chk_low)

@Test
subroutine ADD_PREFIX(erase_range_bnd_chk_high)()

  type( VECTOR_NAME ) :: vec

  vec = test_array_2

  call vec%erase(1, size(test_array_2)+1)

  call assertExceptionRaised(OOBMsg("erase", [1, vec%vsize()], size(test_array_2)+1))

end subroutine ADD_PREFIX(erase_range_bnd_chk_high)

@Test
subroutine ADD_PREFIX(shrink_capacity_empty)()

  type( VECTOR_NAME ) :: vec

  call vec%shrink_to_fit()

  @assertEqual(0, vec%capacity())

end subroutine ADD_PREFIX(shrink_capacity_empty)

@Test
subroutine ADD_PREFIX(shrink_capacity_alloc)()

  type( VECTOR_NAME ) :: vec

  call vec%reserve(5)

  call vec%shrink_to_fit()

  @assertEqual(0, vec%capacity())

end subroutine ADD_PREFIX(shrink_capacity_alloc)

@Test
subroutine ADD_PREFIX(shrink_capacity_full)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%shrink_to_fit()

  @assertEqual(size(test_array), vec%capacity())

end subroutine ADD_PREFIX(shrink_capacity_full)

@Test
subroutine ADD_PREFIX(reserve)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%reserve(size(test_array) + 1)

  @assertLessThanOrEqual(size(test_array) + 1, vec%capacity())

end subroutine ADD_PREFIX(reserve)

@Test
subroutine ADD_PREFIX(reserve_noop)()

  type( VECTOR_NAME ) :: vec

  call vec%reserve(5)

  call vec%reserve(3)

  @assertLessThanOrEqual(5, vec%capacity())

end subroutine ADD_PREFIX(reserve_noop)

@Test
subroutine ADD_PREFIX(shrink_capacity_fit)()

  type( VECTOR_NAME ) :: vec

  vec = test_array

  call vec%reserve(5)

  call vec%shrink_to_fit()

  @assertEqual(size(test_array), vec%capacity())

end subroutine ADD_PREFIX(shrink_capacity_fit)

@Test
subroutine ADD_PREFIX(move_in)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: new_array(:)

  new_array = test_array

  call vec%move_in(new_array)

  @assertEqual(test_array, vec%get())

end subroutine ADD_PREFIX(move_in)

@Test
subroutine ADD_PREFIX(move_in_empty)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: new_array(:)

  call vec%move_in(new_array)

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(move_in_empty)

@Test
subroutine ADD_PREFIX(move_in_overwrite_with_empty)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: new_array(:)

  vec = test_array

  call vec%move_in(new_array)

  @assertTrue(vec%empty())

end subroutine ADD_PREFIX(move_in_overwrite_with_empty)

@Test
subroutine ADD_PREFIX(move_out)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: new_array(:)

  vec = test_array

  call vec%move_out(new_array)

  @assertEqual(test_array, new_array)

end subroutine ADD_PREFIX(move_out)

@Test
subroutine ADD_PREFIX(move_out_empty)()

  type( VECTOR_NAME ) :: vec

  TYPE_DECL, allocatable :: new_array(:)

  call vec%move_out(new_array)

  @assertFalse(allocated(new_array))

end subroutine ADD_PREFIX(move_out_empty)

@Test
subroutine ADD_PREFIX(swap)()

  type( VECTOR_NAME ) :: vec1
  type( VECTOR_NAME ) :: vec2

  vec1 = test_array
  vec2 = test_array_2

  call vec1%swap(vec2)

  @assertEqual(test_array_2, vec1%get())
  @assertEqual(test_array, vec2%get())

end subroutine ADD_PREFIX(swap)

@Test
subroutine ADD_PREFIX(swap_1st_empty)()

  type( VECTOR_NAME ) :: vec1
  type( VECTOR_NAME ) :: vec2

  vec2 = test_array_2

  call vec1%swap(vec2)

  @assertEqual(test_array_2, vec1%get())
  @assertTrue(vec2%empty())

end subroutine ADD_PREFIX(swap_1st_empty)

@Test
subroutine ADD_PREFIX(swap_2nd_empty)()

  type( VECTOR_NAME ) :: vec1
  type( VECTOR_NAME ) :: vec2

  vec1 = test_array

  call vec1%swap(vec2)

  @assertTrue(vec1%empty())
  @assertEqual(test_array, vec2%get())

end subroutine ADD_PREFIX(swap_2nd_empty)

@Test
subroutine ADD_PREFIX(swap_both_empty)()

  type( VECTOR_NAME ) :: vec1
  type( VECTOR_NAME ) :: vec2

  call vec1%swap(vec2)

  @assertTrue(vec1%empty())
  @assertTrue(vec2%empty())

end subroutine ADD_PREFIX(swap_both_empty)

@Test
subroutine ADD_PREFIX(swap_with_self)()

  type( VECTOR_NAME ) :: vec1

  vec1 = test_array

  call vec1%swap(vec1)

  @assertEqual(test_array, vec1%get())

end subroutine ADD_PREFIX(swap_with_self)
