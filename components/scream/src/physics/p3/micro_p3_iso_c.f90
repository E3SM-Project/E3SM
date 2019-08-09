module micro_p3_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
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

  subroutine p3_init_c(lookup_file_dir_c, info) bind(c)
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
    character(len=16) :: p3_version="2.8.2"

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

  subroutine p3_main_c(qc,nc,qr,nr,th,qv,dt,qitot,qirim,nitot,birim,   &
       pres,dzq,npccn,naai,it,prt_liq,prt_sol,its,ite,kts,kte,diag_ze,diag_effc,     &
       diag_effi,diag_vmi,diag_di,diag_rhoi,log_predictNc_in, &
       pdel,exner,cmeiout,prain,nevapr,prer_evap,rflx,sflx,rcldm,lcldm,icldm, &
       pratot,prctot,p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange, &
       vap_ice_exchange, vap_cld_exchange) bind(C)
    use micro_p3, only : p3_main

    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qc, nc, qr, nr, qv, th
    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qitot, qirim, nitot, birim
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: pres, dzq
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: npccn,naai
    real(kind=c_real), value, intent(in) :: dt
    real(kind=c_real), intent(out), dimension(its:ite) :: prt_liq, prt_sol
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_ze, diag_effc
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_effi, diag_vmi, diag_di, diag_rhoi
    integer(kind=c_int), value, intent(in) :: its,ite, kts,kte, it
    logical(kind=c_bool), value, intent(in) :: log_predictNc_in

    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: pdel
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: exner
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: cmeiout
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: prain
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: nevapr
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: prer_evap
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte+1)    :: rflx
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte+1)    :: sflx
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: icldm, lcldm, rcldm
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: pratot,prctot
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte,49)   :: p3_tend_out
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: mu_c,lamc
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: liq_ice_exchange
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: vap_liq_exchange
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: vap_ice_exchange
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: vap_cld_exchange

    logical :: log_predictNc

    log_predictNc = log_predictNc_in

    call p3_main(qc,nc,qr,nr,th,qv,dt,qitot,qirim,nitot,birim,   &
         pres,dzq,npccn,naai,it,prt_liq,prt_sol,its,ite,kts,kte,diag_ze,diag_effc,     &
         diag_effi,diag_vmi,diag_di,diag_rhoi,log_predictNc, &
         pdel,exner,cmeiout,prain,nevapr,prer_evap,rflx,sflx,rcldm,lcldm,icldm, &
         pratot,prctot,p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange, &
         vap_ice_exchange, vap_cld_exchange)
  end subroutine p3_main_c

   
  subroutine micro_p3_utils_init_c(Cpair, Rair, RH2O, RhoH2O, &
                 MWH2O, MWdry, gravit, LatVap, LatIce,        &
                 CpLiq, Tmelt, Pi, iulog_in, masterproc_in) bind(C)

    use micro_p3_utils, only: micro_p3_utils_init
    real(kind=c_real), value, intent(in) :: Cpair
    real(kind=c_real), value, intent(in) :: Rair
    real(kind=c_real), value, intent(in) :: RH2O
    real(kind=c_real), value, intent(in) :: RhoH2O
    real(kind=c_real), value, intent(in) :: MWH2O
    real(kind=c_real), value, intent(in) :: MWdry
    real(kind=c_real), value, intent(in) :: gravit
    real(kind=c_real), value, intent(in) :: LatVap
    real(kind=c_real), value, intent(in) :: LatIce
    real(kind=c_real), value, intent(in) :: CpLiq
    real(kind=c_real), value, intent(in) :: Tmelt
    real(kind=c_real), value, intent(in) :: Pi
    integer(kind=c_int), value, intent(in)   :: iulog_in
    logical(kind=c_bool), value, intent(in)  :: masterproc_in

    logical :: masterproc
    integer :: iulog

    masterproc = masterproc_in
    iulog = iulog_in

    call micro_p3_utils_init(Cpair,Rair,RH2O,RhoH2O,MWH2O,MWdry,gravit,LatVap,LatIce, &
                   CpLiq,Tmelt,Pi,iulog,masterproc)
  end subroutine micro_p3_utils_init_c

end module micro_p3_iso_c
