module harness_handle_mod
   use iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_char, &
                            c_null_char, c_associated
   use api_functions_mod
   implicit none
   private

   public :: EMULATOR_REAL64, EMULATOR_INT32

   enum, bind(C)
      enumerator :: EMULATOR_REAL64 = 0
      enumerator :: EMULATOR_INT32 = 1
   end enum

   type, public :: harness_handle
      private
      type(c_ptr) :: handle = c_null_ptr
   contains
      procedure :: init
      procedure :: register_field
      final :: finalize
   end type

contains

   subroutine init(self)
      class(harness_handle), intent(inout) :: self

      if (c_associated(self%handle)) error stop "Harness already created"

      self%handle = c_harness_create()
      if (.not. c_associated(self%handle)) error stop "Harness creation failed"
   end subroutine

   subroutine register_field(self, name, extent, dtype, ptr)
      class(harness_handle), intent(in) :: self
      character(len=*), intent(in) :: name
      integer(c_int), intent(in) :: extent, dtype
      type(c_ptr), intent(in) :: ptr

      if (.not. c_associated(self%handle)) error stop "Harness not created"

      call c_harness_register_field(self%handle, trim(name)//c_null_char, &
                                    extent, dtype, ptr)
   end subroutine

   subroutine finalize(self)
      type(harness_handle), intent(inout) :: self

      if (c_associated(self%handle)) then
         call c_harness_destroy(self%handle)
         self%handle = c_null_ptr
      end if
   end subroutine

end module harness_handle_mod

