!> Fortran side of cpl7_harness_c_api.cpp: lets a CPL7 cap drive one emulator
!> through the Harness. Errors come back as text; the wrappers return it so the
!> cap can abort with shr_sys_abort and the C++ message.
module emulator_cpl7_mod

   use iso_c_binding
   implicit none
   private

   public :: emulator_cpl7_type
   public :: emulator_cpl7_init, emulator_cpl7_step, emulator_cpl7_final
   public :: emulator_cpl7_write_restart, emulator_cpl7_read_restart

   integer, parameter :: errlen = 4096

   type :: emulator_cpl7_type
      type(c_ptr) :: handle = c_null_ptr
   end type emulator_cpl7_type

   interface
      integer(c_int) function c_create(f_comm, config_path, plan_path, lsize, gids, gsize, &
            coupling_dt, x2c_fields, x2c, x2c_nflds, c2x_fields, c2x, c2x_nflds, handle,   &
            errmsg, elen) bind(C, name="emulator_cpl7_create")
         import :: c_int, c_char, c_double, c_ptr
         integer(c_int), value :: f_comm, lsize, gsize, coupling_dt, x2c_nflds, c2x_nflds, elen
         character(kind=c_char), intent(in) :: config_path(*), plan_path(*)
         character(kind=c_char), intent(in) :: x2c_fields(*), c2x_fields(*)
         integer(c_int), intent(in) :: gids(*)
         type(c_ptr), value :: x2c, c2x
         type(c_ptr), intent(out) :: handle
         character(kind=c_char), intent(out) :: errmsg(*)
      end function c_create

      integer(c_int) function c_run(handle, ymd, tod, errmsg, elen) bind(C, name="emulator_cpl7_run")
         import :: c_int, c_char, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: ymd, tod, elen
         character(kind=c_char), intent(out) :: errmsg(*)
      end function c_run

      integer(c_int) function c_restart(handle, path, errmsg, elen) bind(C, name="emulator_cpl7_write_restart")
         import :: c_int, c_char, c_ptr
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
         integer(c_int), value :: elen
         character(kind=c_char), intent(out) :: errmsg(*)
      end function c_restart

      integer(c_int) function c_read_restart(handle, path, errmsg, elen) bind(C, name="emulator_cpl7_read_restart")
         import :: c_int, c_char, c_ptr
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
         integer(c_int), value :: elen
         character(kind=c_char), intent(out) :: errmsg(*)
      end function c_read_restart

      subroutine c_summary(handle, out, len) bind(C, name="emulator_cpl7_summary")
         import :: c_int, c_char, c_ptr
         type(c_ptr), value :: handle
         character(kind=c_char), intent(out) :: out(*)
         integer(c_int), value :: len
      end subroutine c_summary

      subroutine c_finalize(handle) bind(C, name="emulator_cpl7_finalize")
         import :: c_ptr
         type(c_ptr), value :: handle
      end subroutine c_finalize
   end interface

contains

   !> Creates the harness. x2c and c2x are the first element of the attribute
   !> vectors' rAttr (nflds x lsize); they must stay allocated for the run.
   subroutine emulator_cpl7_init(emu, comm, config_path, plan_path, gids, gsize, coupling_dt, &
         x2c_fields, x2c, c2x_fields, c2x, summary, errmsg, rc)
      type(emulator_cpl7_type), intent(inout) :: emu
      integer, intent(in) :: comm, gsize, coupling_dt
      character(len=*), intent(in) :: config_path, plan_path, x2c_fields, c2x_fields
      integer, intent(in) :: gids(:)
      real(c_double), target, intent(inout) :: x2c(:,:), c2x(:,:)
      character(len=*), intent(out) :: summary, errmsg
      integer, intent(out) :: rc

      character(kind=c_char) :: cerr(errlen), csum(errlen)
      integer(c_int), allocatable :: cgids(:)
      type(c_ptr) :: px2c, pc2x

      allocate(cgids(size(gids)))
      cgids = int(gids, c_int)
      px2c = c_null_ptr
      pc2x = c_null_ptr
      if (size(x2c) > 0) px2c = c_loc(x2c(1,1))
      if (size(c2x) > 0) pc2x = c_loc(c2x(1,1))
      rc = c_create(int(comm, c_int), cstr(config_path), cstr(plan_path), int(size(gids), c_int), &
         cgids, int(gsize, c_int), int(coupling_dt, c_int), cstr(x2c_fields), px2c,               &
         int(size(x2c, 1), c_int), cstr(c2x_fields), pc2x, int(size(c2x, 1), c_int), emu%handle,  &
         cerr, int(errlen, c_int))
      errmsg = ''
      summary = ''
      if (rc /= 0) then
         errmsg = fstr(cerr)
         return
      end if
      call c_summary(emu%handle, csum, int(errlen, c_int))
      summary = fstr(csum)
   end subroutine emulator_cpl7_init

   !> One coupling interval that starts at (ymd, tod).
   subroutine emulator_cpl7_step(emu, ymd, tod, errmsg, rc)
      type(emulator_cpl7_type), intent(inout) :: emu
      integer, intent(in) :: ymd, tod
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: rc
      character(kind=c_char) :: cerr(errlen)
      rc = c_run(emu%handle, int(ymd, c_int), int(tod, c_int), cerr, int(errlen, c_int))
      errmsg = ''
      if (rc /= 0) errmsg = fstr(cerr)
   end subroutine emulator_cpl7_step

   !> Saves the emulator's state, for an exact restart.
   subroutine emulator_cpl7_write_restart(emu, path, errmsg, rc)
      type(emulator_cpl7_type), intent(inout) :: emu
      character(len=*), intent(in) :: path
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: rc
      character(kind=c_char) :: cerr(errlen)
      rc = c_restart(emu%handle, cstr(path), cerr, int(errlen, c_int))
      errmsg = ''
      if (rc /= 0) errmsg = fstr(cerr)
   end subroutine emulator_cpl7_write_restart

   !> Restores what emulator_cpl7_write_restart saved.
   subroutine emulator_cpl7_read_restart(emu, path, errmsg, rc)
      type(emulator_cpl7_type), intent(inout) :: emu
      character(len=*), intent(in) :: path
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: rc
      character(kind=c_char) :: cerr(errlen)
      rc = c_read_restart(emu%handle, cstr(path), cerr, int(errlen, c_int))
      errmsg = ''
      if (rc /= 0) errmsg = fstr(cerr)
   end subroutine emulator_cpl7_read_restart

   subroutine emulator_cpl7_final(emu)
      type(emulator_cpl7_type), intent(inout) :: emu
      if (c_associated(emu%handle)) call c_finalize(emu%handle)
      emu%handle = c_null_ptr
   end subroutine emulator_cpl7_final

   function cstr(s) result(c)
      character(len=*), intent(in) :: s
      character(kind=c_char) :: c(len_trim(s) + 1)
      integer :: i
      do i = 1, len_trim(s)
         c(i) = s(i:i)
      end do
      c(len_trim(s) + 1) = c_null_char
   end function cstr

   function fstr(c) result(s)
      character(kind=c_char), intent(in) :: c(:)
      character(len=size(c)) :: s
      integer :: i
      s = ''
      do i = 1, size(c)
         if (c(i) == c_null_char) exit
         s(i:i) = c(i)
      end do
   end function fstr

end module emulator_cpl7_mod
