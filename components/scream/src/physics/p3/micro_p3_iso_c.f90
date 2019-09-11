module micro_p3_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to  micro_p3 fortran.
!

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

  subroutine p3_use_cxx_c(arg_use_cxx) bind(C)
    use micro_p3, only: use_cxx

    logical(kind=c_bool), value, intent(in) :: arg_use_cxx

    use_cxx = arg_use_cxx
  end subroutine p3_use_cxx_c

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

  subroutine p3_init_a_c(itab_c, itabcoll_c) bind(C)
    use micro_p3, only: itab, itabcoll
    use micro_p3_utils, only: densize,rimsize,isize,tabsize,rcollsize,colltabsize

    real(kind=c_real), intent(out), dimension(densize,rimsize,isize,tabsize) :: itab_c
    real(kind=c_real), intent(out), dimension(densize,rimsize,isize,rcollsize,colltabsize) :: itabcoll_c

    itab_c(:,:,:,:) = itab(:,:,:,:)
    itabcoll_c(:,:,:,:,:) = itabcoll(:,:,:,:,:)
  end subroutine p3_init_a_c

  subroutine find_lookuptable_indices_1a_c(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qitot,nitot,qirim,rhop) bind(C)
    use micro_p3, only: find_lookupTable_indices_1a
    use micro_p3_utils, only: densize,rimsize,isize

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real), value, intent(in)  :: qitot,nitot,qirim,rhop

    call find_lookupTable_indices_1a(dumi, dumjj, dumii, dumzz, dum1, dum4, dum5, dum6,      &
         isize, rimsize, densize, qitot, nitot, qirim, rhop)
  end subroutine find_lookuptable_indices_1a_c

  subroutine find_lookuptable_indices_1b_c(dumj,dum3,qr,nr) bind(C)
    use micro_p3, only: find_lookupTable_indices_1b
    use micro_p3_utils, only: rcollsize

    integer(kind=c_int), intent(out) :: dumj
    real(kind=c_real), intent(out) :: dum3
    real(kind=c_real), value, intent(in) :: qr, nr

    call find_lookupTable_indices_1b(dumj, dum3, rcollsize, qr, nr)
  end subroutine find_lookupTable_indices_1b_c

  subroutine access_lookup_table_c(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc) bind(C)
    use micro_p3, only: access_lookup_table

    integer(kind=c_int), value, intent(in) :: dumjj, dumii, dumi, index
    real(kind=c_real), value, intent(in) :: dum1, dum4, dum5
    real(kind=c_real), intent(out) :: proc

    call access_lookup_table(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc)
  end subroutine access_lookup_table_c

  subroutine access_lookup_table_coll_c(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc) bind(C)
    use micro_p3, only: access_lookup_table_coll

    integer(kind=c_int), value, intent(in) :: dumjj,dumii,dumj,dumi,index
    real(kind=c_real), value, intent(in) :: dum1,dum3,dum4,dum5
    real(kind=c_real), intent(out) :: proc

    call access_lookup_table_coll(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc)
  end subroutine access_lookup_table_coll_c

  subroutine get_cloud_dsd2_c(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1,lcldm) bind(C)
    use micro_p3, only: get_cloud_dsd2
    use micro_p3_utils, only: dnu

    !arguments:
    real(kind=c_real), value, intent(in)        :: qc,rho,lcldm
    real(kind=c_real), intent(inout)            :: nc
    real(kind=c_real), intent(out)              :: mu_c,nu,lamc,cdist,cdist1

    call get_cloud_dsd2(qc,nc,mu_c,rho,nu,dnu,lamc,cdist,cdist1,lcldm)
  end subroutine get_cloud_dsd2_c

  subroutine get_rain_dsd2_c(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm) bind(C)
    use micro_p3, only: get_rain_dsd2

    !arguments:
    real(kind=c_real), value, intent(in) :: qr,rcldm
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r

    call get_rain_dsd2(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm)
  end subroutine get_rain_dsd2_c

end module micro_p3_iso_c
