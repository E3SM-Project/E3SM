module p3_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to micro_p3 fortran.
!

contains
  subroutine init_tables_from_f90_c(vn_table_vals_c, vm_table_vals_c, revap_table_vals_c, mu_table_c) bind(C)
    use micro_p3, only: p3_get_tables

    real(kind=c_real), intent(inout), dimension(300,10) :: vn_table_vals_c, vm_table_vals_c, revap_table_vals_c
    real(kind=c_real), intent(inout), dimension(150)    :: mu_table_c

    real(kind=c_real), dimension(150), target :: mu_table_f
    real(kind=c_real), dimension(300,10), target :: vn_table_vals_f, vm_table_vals_f, revap_table_vals_f

    call p3_get_tables(mu_table_f, revap_table_vals_f, vn_table_vals_f, vm_table_vals_f)
    vn_table_vals_c(:,:)    = vn_table_vals_f(:,:)
    vm_table_vals_c(:,:)    = vm_table_vals_f(:,:)
    revap_table_vals_c(:,:) = revap_table_vals_f(:,:)
    mu_table_c(:)      = mu_table_f(:)

  end subroutine init_tables_from_f90_c

  subroutine p3_init_c(lookup_file_dir_c, info, write_tables) bind(c)
    use ekat_array_io_mod, only: array_io_file_exists
#ifdef SCREAM_DOUBLE_PRECISION
    use ekat_array_io_mod, only: array_io_read=>array_io_read_double, array_io_write=>array_io_write_double
#else
    use ekat_array_io_mod, only: array_io_read=>array_io_read_float, array_io_write=>array_io_write_float
#endif
    use micro_p3, only: p3_init_a, p3_init_b, p3_set_tables, p3_get_tables

    type(c_ptr), intent(in) :: lookup_file_dir_c
    integer(kind=c_int), intent(out) :: info
    logical(kind=c_bool), intent(in) :: write_tables

    real(kind=c_real), dimension(150), target :: mu_r_table_vals
    real(kind=c_real), dimension(300,10), target :: vn_table_vals, vm_table_vals, revap_table_vals

    character(len=256), pointer :: lookup_file_dir
    character(kind=c_char, len=512) :: mu_r_filename, revap_filename, vn_filename, vm_filename
    integer :: len
    logical :: ok
    character(len=16) :: p3_version="4.1.1"  ! TODO: Change to be dependent on table version and path specified in p3_functions.hpp

    call c_f_pointer(lookup_file_dir_c, lookup_file_dir)
    len = index(lookup_file_dir, C_NULL_CHAR) - 1
    call p3_init_a(lookup_file_dir(1:len),p3_version)

    info = 0
    ok = .false.

#ifdef SCREAM_DOUBLE_PRECISION
    mu_r_filename  = lookup_file_dir(1:len)//'/mu_r_table_vals.dat8'//C_NULL_CHAR
    revap_filename = lookup_file_dir(1:len)//'/revap_table_vals.dat8'//C_NULL_CHAR
    vn_filename    = lookup_file_dir(1:len)//'/vn_table_vals.dat8'//C_NULL_CHAR
    vm_filename    = lookup_file_dir(1:len)//'/vm_table_vals.dat8'//C_NULL_CHAR
#else
    mu_r_filename  = lookup_file_dir(1:len)//'/mu_r_table_vals.dat4'//C_NULL_CHAR
    revap_filename = lookup_file_dir(1:len)//'/revap_table_vals.dat4'//C_NULL_CHAR
    vn_filename    = lookup_file_dir(1:len)//'/vn_table_vals.dat4'//C_NULL_CHAR
    vm_filename    = lookup_file_dir(1:len)//'/vm_table_vals.dat4'//C_NULL_CHAR
#endif

    if (write_tables) then
       call p3_init_b()
       call p3_get_tables(mu_r_table_vals, revap_table_vals, vn_table_vals, vm_table_vals)
       ok = array_io_write(mu_r_filename, c_loc(mu_r_table_vals), size(mu_r_table_vals)) .and. &
            array_io_write(revap_filename, c_loc(revap_table_vals), size(revap_table_vals)) .and. &
            array_io_write(vn_filename, c_loc(vn_table_vals), size(vn_table_vals)) .and. &
            array_io_write(vm_filename, c_loc(vm_table_vals), size(vm_table_vals))
       if (.not. ok) then
          print *, 'p3_iso_c::p3_init: Error when writing table files.'
          info = -1
       end if
    else
      ! Check table files exist
      ok = array_io_file_exists(mu_r_filename) .and. &
           array_io_file_exists(revap_filename) .and. &
           array_io_file_exists(vn_filename) .and. &
           array_io_file_exists(vm_filename)
      if (.not. ok) then
        print *, 'p3_iso_c::p3_init: One or more table files does not exist'
        info = -2
        return
      endif

      ! Read files
      if (.not. array_io_read(mu_r_filename, c_loc(mu_r_table_vals), size(mu_r_table_vals))) then
         print *, "p3_iso_c::p3_init: error reading mu_r table from file "//mu_r_filename
         info = -3
         return
      elseif (.not. array_io_read(revap_filename, c_loc(revap_table_vals), size(revap_table_vals))) then
         print *, "p3_iso_c::p3_init: error reading revap table from file "//revap_filename
         info = -4
         return

      elseif (.not. array_io_read(vn_filename, c_loc(vn_table_vals), size(vn_table_vals))) then
         print *, "p3_iso_c::p3_init: error reading vn table from file "//vn_filename
         info = -5
         return
      elseif (.not. array_io_read(vm_filename, c_loc(vm_table_vals), size(vm_table_vals))) then
         print *, "p3_iso_c::p3_init: error reading vm table from file "//vm_filename
         info = -6
         return
      endif

      call p3_set_tables(mu_r_table_vals, revap_table_vals, vn_table_vals, vm_table_vals)
    end if

  end subroutine p3_init_c

  subroutine micro_p3_utils_init_c(Cpair, Rair, RH2O, RHO_H2O, &
                 MWH2O, MWdry, gravit, LatVap, LatIce,        &
                 CpLiq, Tmelt, Pi, masterproc) bind(C)

    use micro_p3_utils, only: micro_p3_utils_init
    use iso_fortran_env, only: OUTPUT_UNIT
    real(kind=c_real), value, intent(in) :: Cpair
    real(kind=c_real), value, intent(in) :: Rair
    real(kind=c_real), value, intent(in) :: RH2O
    real(kind=c_real), value, intent(in) :: RHO_H2O
    real(kind=c_real), value, intent(in) :: MWH2O
    real(kind=c_real), value, intent(in) :: MWdry
    real(kind=c_real), value, intent(in) :: gravit
    real(kind=c_real), value, intent(in) :: LatVap
    real(kind=c_real), value, intent(in) :: LatIce
    real(kind=c_real), value, intent(in) :: CpLiq
    real(kind=c_real), value, intent(in) :: Tmelt
    real(kind=c_real), value, intent(in) :: Pi
    logical(kind=c_bool), value, intent(in)  :: masterproc

    call micro_p3_utils_init(Cpair,Rair,RH2O,RHO_H2O,MWH2O,MWdry,gravit,LatVap,LatIce, &
                   CpLiq,Tmelt,Pi,OUTPUT_UNIT,masterproc)
  end subroutine micro_p3_utils_init_c

  subroutine p3_init_a_c(ice_table_vals_c, collect_table_vals_c) bind(C)
    use micro_p3, only: ice_table_vals, collect_table_vals
    use micro_p3_utils, only: densize,rimsize,isize,ice_table_size,rcollsize,collect_table_size

    real(kind=c_real), intent(out), dimension(densize,rimsize,isize,ice_table_size) :: ice_table_vals_c
    real(kind=c_real), intent(out), dimension(densize,rimsize,isize,rcollsize,collect_table_size) :: collect_table_vals_c

    ice_table_vals_c(:,:,:,:)       = ice_table_vals(:,:,:,:)
    collect_table_vals_c(:,:,:,:,:) = collect_table_vals(:,:,:,:,:)
  end subroutine p3_init_a_c

end module p3_iso_c
