module micro_p3_iso_c
  use iso_c_binding
  implicit none

#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

contains
  subroutine append_precision(string, prefix)
    use iso_c_binding

    character(kind=c_char, len=128), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

  subroutine p3_init_c(lookup_file_dir_c, ncat, info) bind(c)
    use array_io_mod
    use micro_p3, only: p3_init_a, p3_init_b, p3_set_tables, p3_get_tables

    type(c_ptr), intent(in) :: lookup_file_dir_c
    integer(kind=c_int), value, intent(in) :: ncat
    integer(kind=c_int), intent(out) :: info

    real(kind=c_real), dimension(150), target :: mu_r_table
    real(kind=c_real), dimension(300,10), target :: vn_table, vm_table, revap_table

    character(len=256), pointer :: lookup_file_dir
    character(kind=c_char, len=128) :: mu_r_filename, revap_filename, vn_filename, vm_filename
    integer :: len
    logical :: ok

    call c_f_pointer(lookup_file_dir_c, lookup_file_dir)
    len = index(lookup_file_dir, C_NULL_CHAR) - 1
    call p3_init_a(lookup_file_dir(1:len), ncat)

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
          print *, 'micro_p3_iso_c::p3_init: One or more table files exists but gave a read error.'
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
          print *, 'micro_p3_iso_c::p3_init: Error when writing table files.'
          info = -1
       end if
    end if
  end subroutine p3_init_c

  subroutine p3_main_c(qc,nc,qr,nr,th_old,th,qv_old,qv,dt,qitot,qirim,nitot,birim,ssat,uzpl,   &
       pres,dzq,it,prt_liq,prt_sol,its,ite,kts,kte,nCat,diag_ze,diag_effc,     &
       diag_effi,diag_vmi,diag_di,diag_rhoi,n_diag_2d,diag_2d,n_diag_3d,       &
       diag_3d,log_predictNc_in,typeDiags_ON_in,model_in,prt_drzl,prt_rain,prt_crys,    &
       prt_snow,prt_grpl,prt_pell,prt_hail,prt_sndp) bind(C)
    use micro_p3, only : p3_main

    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qc, nc, qr, nr, ssat, qv, th, th_old, qv_old
    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte,nCat) :: qitot, qirim, nitot, birim
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: uzpl, pres, dzq
    real(kind=c_real), value, intent(in) :: dt
    real(kind=c_real), intent(out), dimension(its:ite) :: prt_liq, prt_sol
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_ze, diag_effc
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte,nCat) :: diag_effi, diag_vmi, diag_di, diag_rhoi
    real(kind=c_real), intent(out), dimension(its:ite,n_diag_2d) :: diag_2d
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte,n_diag_3d) :: diag_3d
    integer(kind=c_int), value, intent(in) :: its,ite, kts,kte, it, nCat, n_diag_2d, n_diag_3d
    logical(kind=c_bool), value, intent(in) :: log_predictNc_in, typeDiags_ON_in
    type(c_ptr), intent(in) :: model_in
    real(kind=c_real), intent(out), dimension(its:ite), optional :: &
         prt_drzl, prt_rain, prt_crys, prt_snow, prt_grpl, prt_pell, prt_hail, prt_sndp

    character(len=64), pointer :: model
    logical :: log_predictNc, typeDiags_ON
    integer :: len

    log_predictNc = log_predictNc_in
    typeDiags_ON = typeDiags_ON_in

    call c_f_pointer(model_in, model)
    len = index(model, C_NULL_CHAR) - 1

    call p3_main(qc,nc,qr,nr,th_old,th,qv_old,qv,dt,qitot,qirim,nitot,birim,ssat,uzpl,   &
         pres,dzq,it,prt_liq,prt_sol,its,ite,kts,kte,nCat,diag_ze,diag_effc,     &
         diag_effi,diag_vmi,diag_di,diag_rhoi,n_diag_2d,diag_2d,n_diag_3d,       &
         diag_3d,log_predictNc,typeDiags_ON,model(1:len),prt_drzl,prt_rain,prt_crys,    &
         prt_snow,prt_grpl,prt_pell,prt_hail,prt_sndp)
  end subroutine p3_main_c
end module micro_p3_iso_c
