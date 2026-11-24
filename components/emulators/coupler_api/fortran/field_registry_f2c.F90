module field_registry_f2c
   use coupler_types, only: registered_field_desc
   use iso_c_binding, only: c_ptr
   implicit none
   private

   interface
      subroutine register_field(handle, desc) bind(c)
         import :: registered_field_desc, c_ptr
         type(c_ptr), intent(in) :: handle
         type(registered_field_desc), intent(in) :: desc
      end subroutine register_field
   end interface

end module field_registry_f2c
