module shoc_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to shoc fortran.
!

contains
  subroutine append_precision(string, prefix)

    character(kind=c_char, len=128), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

  subroutine shoc_init_c(lookup_file_dir_c, info) bind(c)
    use array_io_mod
    use micro_p3, only: p3_init_a, p3_init_b, p3_set_tables, p3_get_tables

    type(c_ptr), intent(in) :: lookup_file_dir_c
    integer(kind=c_int), intent(out) :: info

    real(kind=c_real), dimension(150), target :: mu_r_table
    real(kind=c_real), dimension(300,10), target :: vn_table, vm_table, revap_table

    character(len=256), pointer :: lookup_file_dir
    character(kind=c_char, len=128) :: mu_r_filename, revap_filename, vn_filename, vm_filename
    integer :: len
    logical :: ok
    character(len=16) :: p3_version="4"  ! TODO: Change to be dependent on table version and path specified in p3_functions.hpp

    call c_f_pointer(lookup_file_dir_c, lookup_file_dir)
    len = index(lookup_file_dir, C_NULL_CHAR) - 1
    call p3_init_a(lookup_file_dir(1:len),p3_version)

    info = 0
    ok = .false.

    call append_precision(mu_r_filename, c_char_"mu_r_table.dat")
    call append_precision(revap_filename, c_char_"revap_table.dat")
    call append_precision(vn_filename, c_char_"vn_table.dat")
    call append_precision(vm_filename, c_char_"vm_table.dat")
    ok = array_io_file_exists(mu_r_filename) .and. &
         array_io_file_exists(revap_filename) .and. &
         array_io_file_exists(vn_filename) .and. &
         array_io_file_exists(vm_filename)
    if (ok) then
       ok = array_io_read(mu_r_filename, c_loc(mu_r_table), size(mu_r_table)) .and. &
            array_io_read(revap_filename, c_loc(revap_table), size(revap_table)) .and. &
            array_io_read(vn_filename, c_loc(vn_table), size(vn_table)) .and. &
            array_io_read(vm_filename, c_loc(vm_table), size(vm_table))
       if (.not. ok) then
          print *, 'shoc_iso_c::shoc_init: One or more table files exists but gave a read error.'
          info = -1
       end if
    end if

    if (ok) then
       call p3_set_tables(mu_r_table, revap_table, vn_table, vm_table)
    else
       call p3_init_b()
       call p3_get_tables(mu_r_table, revap_table, vn_table, vm_table)
       ok = array_io_write(mu_r_filename, c_loc(mu_r_table), size(mu_r_table)) .and. &
            array_io_write(revap_filename, c_loc(revap_table), size(revap_table)) .and. &
            array_io_write(vn_filename, c_loc(vn_table), size(vn_table)) .and. &
            array_io_write(vm_filename, c_loc(vm_table), size(vm_table))
       if (.not. ok) then
          print *, 'shoc_iso_c::shoc_init: Error when writing table files.'
          info = -1
       end if
    end if

  end subroutine shoc_init_c

  subroutine shoc_main_c(shcol,nlev,nlevi,dtime,nadv,host_dx, host_dy,thv,   &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner,phis, host_dse, tke, thetal,  &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt) bind(C)
    use shoc, only : shoc_main

    integer(kind=c_int), value, intent(in) :: shcol, nlev, nlevi, num_qtracers, nadv
    real(kind=c_real), intent(in), dimension(shcol) :: host_dx, host_dy
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: zt_grid
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: zi_grid
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pres
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: presi
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pdel, thv, w_field
    real(kind=c_real), intent(in), dimension(shcol) :: wthl_sfc, wqw_sfc, uw_sfc, vw_sfc
    real(kind=c_real), intent(in), dimension(shcol, num_qtracers) :: wtracer_sfc
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: exner
    real(kind=c_real), intent(in), dimension(shcol) :: phis

    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: host_dse, tke, &
       thetal, qw, u_wind, v_wind, wthv_sec
    real(kind=c_real), intent(inout), dimension(shcol, nlev, num_qtracers) :: qtracers
    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: tk, tkh

    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_cldfrac, shoc_ql
    real(kind=c_real), intent(out), dimension(shcol) :: pblh
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_mix, w_sec
    real(kind=c_real), intent(out), dimension(shcol, nlevi) :: thl_sec, qw_sec, &
       qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: wqls_sec, isotropy, brunt

    call shoc_main(shcol,nlev,nlevi,dtime,nadv,host_dx, host_dy,thv,   &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner,phis, host_dse, tke, thetal,  &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt)
  end subroutine shoc_main_c

  subroutine shoc_use_cxx_c(arg_use_cxx) bind(C)
    use shoc, only: use_cxx

    logical(kind=c_bool), value, intent(in) :: arg_use_cxx

    use_cxx = arg_use_cxx
  end subroutine shoc_use_cxx_c

end module shoc_iso_c
