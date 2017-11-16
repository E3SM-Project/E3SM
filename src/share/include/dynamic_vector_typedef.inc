!
! Clone of C++ standard library vectors
!
! This type is a wrapper for an allocatable array, which provides
! efficient utilities for dynamic array operations, such as appending new
! elements, truncation, and reserving/retaining memory independently from
! changes to the array's apparent size.
!
! Dynamic arrays allocate a somewhat larger buffer of contiguous memory
! (the "capacity") than is actually being used at any given time (the
! "size"). This allow elements to be efficiently added to one end, with the
! object automatically reallocating a new buffer as necessary whenever the
! current capacity is exhausted. The capacity increases geometrically,
! wasting O(N) space, but requiring only O(1) time (amortized) to add each
! element.
!
! One downside is that this wrapper class does not support many of
! Fortran's intrinsic array operations. For instance, if you have a
! vector of reals, and you want to take the sine, you have to either
! iterate in a loop (slow), or set the upper bound yourself (without the
! safety of bounds checking). The latter looks like this:
!
!   x = sin(vec%data(:vec%size()))
!
! Because of this, it's probably preferable to use a standard array instead
! of a vector of reals for numerical work.
!
! Because this type uses an allocatable instead of a pointer, it should not
! cause a memory leak. However, deallocation can be forced by using "clear"
! followed by "shrink_to_fit", or by explicit deallocation of the data
! component.
!
! How to create a vector type:
! ----------------------------
!
! Define VECTOR_NAME and TYPE_NAME in a module, then include this file
! to create the type. Include this file before "contains" in the module,
! and the "procdef" file afterward.
!
! There must be a function in scope called OOBMsg (or a function macro of
! this name). This must accept a string representing the operation, a size
! 2 integer array representing the bounds of the array, and an integer
! representing an index into the array. It should return a string
! representing an error message for out-of-bounds access.
!
! Finally, define the function macro THROW to an error handling mechanism.
! THROW accepts one argument, a string representing an error message.
!
! Some tips:
! ----------
!
!  - Do not directly use the "data" component, unless it's unavoidable to
!    get decent efficiency.
!
!  - The data is assumed to always have lower bound 1.
!
!  - If you are finished with adding/removing elements, you can convert
!    this type into a standard allocatable array with the "move_out"
!    method. (You can do the reverse conversion cheaply with "move_in".)
!
!  - Don't include these files twice in the same module, as this will cause
!    name clashes.
!
!  - Don't use this type if you need pointers into the array to remain
!    valid as you add and remove elements. As with the C++ type, the array
!    is often reallocated if you are adding elements, and this invalidates
!    pointers into it.
!
! Advanced features:
! ------------------
!
!  - Define the macro "USE_PURE" if you need to mark all methods as pure.
!    This effectively requires errors to be silent (because THROW cannot do
!    anything useful if it has no side effects).
!
! Developer's notes:
! ------------------
!
! 1) The main difference from the C++ types is that we use Fortran array
!    indexing conventions:
!     - Indexing starts at 1, not 0.
!     - The last element of a range is included in the range. E.g. using
!       "vec%erase(2,3)" erases two elements, not just one.
!     - When an array size would be negative, it is treated as size 0.
!     - When an operation's ending index is smaller than the beginning
!       index, it is a no-op (unless a negative stride is provided).
!
! 2) We could have iterator types like in C++, but they don't really give
!    you anything more than integer indices. Other types, like linked
!    lists, will likely require companion iterator types.
!
! 3) For access with bounds-checking, C++ uses the "at" method to provide
!    references to individual elements. To avoid working with pointers, and
!    to provide an interface somewhat closer to Fortran array conventions,
!    this type uses set/get methods instead.
!
!    These methods are likely to produce extra copies, which may negatively
!    impact performance compared to direct access of the underlying data.
!    This is one reason why the data component is public, not private.
!
! 4) All vectors with vec_size = 0 are valid empty vectors, regardless of
!    whether or not "data" is allocated, and regardless of its size. This
!    slightly complicates some of the methods. However, it means that the
!    user does not have to initialize vectors, or treat empty vectors
!    differently depending on how they became empty.
!
! 5) Dynamic arrays have a time/space tradeoff parameter, which is the
!    factor by which the array's capacity grows whenever it is
!    automatically reallocated to hold more elements. In this code, the
!    factor is 2, which is a common, simple, and reasonably fast choice.
!
!    If there is too much wasted memory over a wide range of use cases,
!    however, it may be reasonable to consider using 1.5 or even lower
!    (with appropriate attention given to rounding issues).

type VECTOR_NAME

   TYPE_NAME, allocatable :: data(:)

   integer, private :: vec_size = 0
 contains

   !------------------------
   ! Query functions
   !------------------------

   ! Test whether there are any elements present
   procedure, pass(self) :: empty => empty_vec
   ! Return current size
   procedure, pass(self) :: vsize => size_vec
   ! Estimate maximum possible size
   procedure, pass(self) :: max_size => max_size_vec
   ! Return maximum number of elements that can be held before the data
   ! array will be reallocated to a larger size.
   procedure, pass(self) :: capacity => capacity_vec

   !------------------------
   ! Retrieving data
   !------------------------

   ! Get the value of the element at a particular index
   procedure, pass(self), private :: get_single_vec
   ! Get an array of values of all the elements within a range
   procedure, pass(self), private :: get_range_vec
   ! Get a copy of all the data
   procedure, pass(self), private :: get_array_vec
   ! Generic for all of the above.
   generic :: get => get_single_vec, get_range_vec, get_array_vec

   ! Get the value of the first element
   procedure, pass(self) :: front => front_vec
   ! Get the value of the last element
   procedure, pass(self) :: back => back_vec

   !------------------------
   ! Modifying data
   !------------------------

   ! Reset the vector to size 0 (without changing capacity)
   procedure, pass(self) :: clear => clear_vec

   ! Resize the vector (will not reduce capacity)
   ! Resizing to a larger size than the capacity causes reallocation of the
   ! data array.
   procedure, pass(self) :: resize => resize_vec

   ! None of the "set" routines below will grow the array. Setting elements
   ! past the end of the vector will result in an out-of-bounds error; use
   ! "insert", "push_back", or explicit resizing to add elements.

   ! Set the element at a particular index
   procedure, pass(self), private :: set_single_vec
   ! Set the elements in a range from an array
   procedure, pass(self), private :: set_range_vec
   ! Fill all the elements in a range with a scalar value
   procedure, pass(self), private :: set_range_fill_vec
   ! Set the data to a copy of some array
   procedure, pass(self), private :: set_array_vec
   ! Fill the data will a scalar value
   procedure, pass(self), private :: set_fill_vec
   ! Generic for all of the above.
   generic :: set => set_single_vec, set_range_vec, set_range_fill_vec, &
        set_array_vec, set_fill_vec

   ! Add an element to the back of the vector
   procedure, pass(self) :: push_back => push_back_vec
   ! Remove the element at the back of the vector
   procedure, pass(self) :: pop_back => pop_back_vec

   ! All of the insert routines add elements; the vector will be expanded
   ! and data shuffled to ensure that this is non-destructive. For a vector
   ! of size n, new elements can be inserted anywhere from 1 to n+1.
   ! Inserting at point n+1 is equivalent to adding the new elements one-
   ! by-one with push_back.

   ! Insert one element at a particular point
   procedure, pass(self), private :: insert_single_vec
   ! Insert all elements from an array at a particular point
   procedure, pass(self), private :: insert_array_vec
   ! Insert multiple copies of the same value at a particular point.
   procedure, pass(self), private :: insert_repeat_vec
   ! Generic for all of the above.
   generic :: insert => insert_single_vec, insert_array_vec, insert_repeat_vec

   ! Erase the element at a particular point
   procedure, pass(self), private :: erase_single_vec
   ! Erase all the elements between two points (inclusive)
   procedure, pass(self), private :: erase_range_vec
   ! Generic for all of the above.
   generic :: erase => erase_single_vec, erase_range_vec

   !------------------------
   ! Adjusting capacity
   !------------------------

   ! Shrink the vector's capacity to fit its size, releasing unneeded
   ! memory
   procedure, pass(self) :: shrink_to_fit => shrink_to_fit_vec

   ! Expand the vector to have at least as much capacity as requested
   ! Mostly useful to avoid unnecessary reallocation when you know that the
   ! data is unlikely to exceed some upper bound on its size.
   procedure, pass(self) :: reserve => reserve_vec

   !------------------------
   ! Move operations
   !------------------------

   ! Convert an allocatable array into a dynamic vector
   ! No copies or reallocations are performed, but afterward the array is
   ! no longer allocated.
   procedure, pass(self) :: move_in => move_in_vec

   ! Convert a dynamic vector to an allocatable array
   ! An empty vector is converted to an unallocated array. A reallocation
   ! and copy is often performed otherwise. Afterward the vector is empty.
   procedure, pass(self) :: move_out => move_out_vec

   ! Swap the contents of this vector with another one
   ! No copies or reallocations are performed.
   procedure, pass(self) :: swap => swap_vec

   !------------------------
   ! Copy/assignment
   !------------------------

   ! Overwrite contents of this vector with those of an array
   procedure, pass(self), private :: array_assign_vec
   ! Overwrite contents of this vector with those of another vector
   procedure, pass(self), private :: vector_assign_vec
   generic :: assignment(=) => array_assign_vec, vector_assign_vec

end type VECTOR_NAME

!------------------------
! Constructors
!------------------------

interface VECTOR_NAME
   ! Construct empty vector
   module procedure new_vector_default
   ! Construct vector as a copy of another vector
   module procedure new_vector_copy
   ! Construct vector with contents from an array
   module procedure new_vector_array
end interface
