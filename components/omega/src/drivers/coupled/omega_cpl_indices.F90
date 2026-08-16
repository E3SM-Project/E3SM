module omega_cpl_indices

   use, intrinsic :: iso_c_binding, only: c_int, c_char

   implicit none
   private

   integer, parameter, public :: num_omega_imports = 2
   integer, parameter, public :: num_omega_exports = 5
   integer, public :: num_coupler_imports, num_coupler_exports

   ! Names of import/export fields as defined by seq_flds_mod
   character(len=32, kind=c_char), public, target :: &
      import_field_names(num_omega_imports), &
      export_field_names(num_omega_exports)

   ! Coupler indices for import/export fields
   integer(kind=c_int), public, target :: &
      import_field_indices(num_omega_imports), &
      export_field_indices(num_omega_exports)

   ! Full CIME x2o/o2x field-name lists (colon-separated), needed by the
   ! MOAB bridge to define coupler tags covering every field CIME expects
   ! on this app, not just the fields Omega currently imports/exports.
   character(len=:, kind=c_char), public, allocatable, target :: &
      cpl_x2o_field_names, &
      cpl_o2x_field_names

   public :: omega_set_cpl_indices

contains

   subroutine omega_set_cpl_indices()

      use, intrinsic :: iso_c_binding, only: c_null_char
      use mct_mod, only: mct_aVect, mct_aVect_init, mct_aVect_clean
      use seq_flds_mod, only: seq_flds_o2x_fields, seq_flds_x2o_fields

      type(mct_aVect) :: x2o, o2x

      ! Create dummy mct_aVect objects, so that we can inquire about field
      ! indices prior to the omega instance, and therefore the HorzMesh
      ! instance which has the real lsize value, is created.
      call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
      call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

      ! full CIME field lists, for the MOAB bridge's tag definitions
      cpl_x2o_field_names = trim(seq_flds_x2o_fields)//c_null_char
      cpl_o2x_field_names = trim(seq_flds_o2x_fields)//c_null_char

      ! total number of import/export fields in coupler data arrays
      num_coupler_imports = size(x2o%rAttr, 1)
      num_coupler_exports = size(o2x%rAttr, 1)

      ! Import (x2o) Coupler field names
      import_field_names(1) = "Foxx_taux"
      import_field_names(2) = "Foxx_tauy"

      ! get mct_avect_index value for each import field name
      call get_indices_from_names( &
         import_field_indices, import_field_names, x2o &
         )

      ! append null_char to each import field name for C interoperability
      call append_null_char_to_names(import_field_names)

      ! Export (o2x) Coupler field names
      export_field_names(1) = "So_t"
      export_field_names(2) = "So_s"
      export_field_names(3) = "So_u"
      export_field_names(4) = "So_v"
      export_field_names(5) = "So_ssh"

      ! get mct_avect_index value for each export field name
      call get_indices_from_names( &
         export_field_indices, export_field_names, o2x &
         )

      ! append null_char to each import field name for C interoperability
      call append_null_char_to_names(export_field_names)

      ! clean up the dummy mct_aVect objects
      call mct_aVect_clean(x2o)
      call mct_aVect_clean(o2x)

   end subroutine omega_set_cpl_indices

   subroutine get_indices_from_names(field_indices, field_names, data)

      use, intrinsic :: iso_c_binding, only: c_char, c_int
      use mct_mod, only: mct_aVect, mct_avect_indexra

      integer(kind=c_int), intent(inout) :: field_indices(:)
      character(len=32, kind=c_char), intent(in) :: field_names(:)
      type(mct_aVect), intent(in) :: data

      integer :: i

      do i = 1, size(field_names)
         field_indices(i) = mct_avect_indexra(data, trim(field_names(i)))
      end do

   end subroutine get_indices_from_names

   subroutine append_null_char_to_names(field_names)

      use, intrinsic :: iso_c_binding, only: c_char, c_null_char

      character(len=32, kind=c_char), intent(inout) :: field_names(:)
      integer :: i

      do i = 1, size(field_names)
         field_names(i) = trim(field_names(i))//c_null_char
      end do

   end subroutine append_null_char_to_names

end module omega_cpl_indices
