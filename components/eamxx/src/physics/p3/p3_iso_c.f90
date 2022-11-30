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
  subroutine append_precision(string, prefix)

    character(kind=c_char, len=256), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

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
    character(kind=c_char, len=256) :: mu_r_filename, revap_filename, vn_filename, vm_filename
    integer :: len
    logical :: ok
    character(len=16) :: p3_version="4.1.1"  ! TODO: Change to be dependent on table version and path specified in p3_functions.hpp

    call c_f_pointer(lookup_file_dir_c, lookup_file_dir)
    len = index(lookup_file_dir, C_NULL_CHAR) - 1
    call p3_init_a(lookup_file_dir(1:len),p3_version)

    info = 0
    ok = .false.

    call append_precision(mu_r_filename, SCREAM_DATA_DIR//"/tables/mu_r_table_vals.dat")
    call append_precision(revap_filename, SCREAM_DATA_DIR//"/tables/revap_table_vals.dat")
    call append_precision(vn_filename, SCREAM_DATA_DIR//"/tables/vn_table_vals.dat")
    call append_precision(vm_filename, SCREAM_DATA_DIR//"/tables/vm_table_vals.dat")

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

  subroutine p3_main_c(qc,nc,qr,nr,th_atm,qv,dt,qi,qm,ni,bm,   &
       pres,dz,nc_nuceat_tend,nccn_prescribed,ni_activated,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_eff_radius_qc,     &
       diag_eff_radius_qi,rho_qi,do_predict_nc,do_prescribed_CCN,dpres,inv_exner,qv2qi_depos_tend, &
       precip_liq_flux,precip_ice_flux,cld_frac_r,cld_frac_l,cld_frac_i,liq_ice_exchange, &
       vap_liq_exchange, vap_ice_exchange, qv_prev, t_prev, elapsed_s) bind(C)
    use micro_p3, only : p3_main

    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qc, nc, qr, nr, qv, th_atm
    real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qi, qm, ni, bm
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: pres, dz
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: nc_nuceat_tend,nccn_prescribed,ni_activated
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: inv_qc_relvar
    real(kind=c_real), value, intent(in) :: dt
    real(kind=c_real), intent(out), dimension(its:ite) :: precip_liq_surf, precip_ice_surf
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_eff_radius_qc
    real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_eff_radius_qi, rho_qi
    integer(kind=c_int), value, intent(in) :: its,ite, kts,kte, it
    logical(kind=c_bool), value, intent(in) :: do_predict_nc,do_prescribed_CCN

    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: dpres
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: inv_exner
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: qv2qi_depos_tend
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_liq_flux
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte+1)    :: precip_ice_flux
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: cld_frac_i, cld_frac_l, cld_frac_r
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: liq_ice_exchange
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: vap_liq_exchange
    real(kind=c_real), intent(out),   dimension(its:ite,kts:kte)      :: vap_ice_exchange
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: qv_prev
    real(kind=c_real), intent(in),    dimension(its:ite,kts:kte)      :: t_prev

    real(kind=c_real), intent(out) :: elapsed_s

    real(kind=c_real), dimension(its:ite,kts:kte,49)   :: p3_tend_out
    real(kind=c_real), dimension(its:ite,3) :: col_location
    real(kind=c_real), dimension(its:ite,kts:kte)      :: mu_c, lamc
    real(kind=c_real), dimension(its:ite,kts:kte)      :: precip_total_tend
    real(kind=c_real), dimension(its:ite,kts:kte)      :: nevapr
    real(kind=c_real), dimension(its:ite,kts:kte)      :: qr_evap_tend
    integer :: i
    do i = its,ite
      col_location(i,:) = real(i)
    end do

    call p3_main(qc,nc,qr,nr,th_atm,qv,dt,qi,qm,ni,bm,   &
         pres,dz,nc_nuceat_tend,nccn_prescribed,ni_activated,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_eff_radius_qc, &
         diag_eff_radius_qi,rho_qi,do_predict_nc,do_prescribed_CCN,dpres,inv_exner,qv2qi_depos_tend,precip_total_tend,nevapr, &
         qr_evap_tend,precip_liq_flux,precip_ice_flux,cld_frac_r,cld_frac_l,cld_frac_i,p3_tend_out,mu_c,lamc,liq_ice_exchange,&
         vap_liq_exchange,vap_ice_exchange,qv_prev,t_prev,col_location,elapsed_s)
  end subroutine p3_main_c

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

  subroutine find_lookuptable_indices_1a_c(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qi,ni,qm,rhop) bind(C)
    use micro_p3, only: find_lookupTable_indices_1a
    use micro_p3_utils, only: densize,rimsize,isize

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real), value, intent(in)  :: qi,ni,qm,rhop

    call find_lookupTable_indices_1a(dumi, dumjj, dumii, dumzz, dum1, dum4, dum5, dum6,      &
         isize, rimsize, densize, qi, ni, qm, rhop)
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

  subroutine back_to_cell_average_c(cld_frac_l,cld_frac_r,cld_frac_i, qc2qr_accret_tend,qr2qv_evap_tend,qc2qr_autoconv_tend,&
    nc_accret_tend,nc_selfcollect_tend,nc2nr_autoconv_tend,nr_selfcollect_tend,nr_evap_tend,ncautr,qi2qv_sublim_tend,nr_ice_shed_tend,qc2qi_hetero_freeze_tend,&
    qr2qi_collect_tend,qc2qr_ice_shed_tend,qi2qr_melt_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,ni2nr_melt_tend,nc_collect_tend,ncshdc,nc2ni_immers_freeze_tend,nr_collect_tend,ni_selfcollect_tend,&
    qv2qi_vapdep_tend,nr2ni_immers_freeze_tend,ni_sublim_tend,qv2qi_nucleat_tend,ni_nucleat_tend,qc2qi_berg_tend) bind(C)

    use micro_p3, only: back_to_cell_average
    real(kind=c_real), value, intent(in) :: cld_frac_l, cld_frac_r, cld_frac_i

    real(kind=c_real), intent(inout) :: qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend,  &
                                        nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend,  &
                                        nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, &
                                        qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend,&
                                        ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend,  &
                                        qc2qi_berg_tend

    call back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,&
      nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend,&
      qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend,&
      qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend, qc2qi_berg_tend)
  end subroutine back_to_cell_average_c

  subroutine cloud_water_conservation_c(qc,dt,qc2qr_autoconv_tend,qc2qr_accret_tend,qc2qi_collect_tend,qc2qi_hetero_freeze_tend,qc2qr_ice_shed_tend,     &
    qc2qi_berg_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend) bind(C)
    use micro_p3, only: cloud_water_conservation

    real(kind=c_real), value, intent(in) :: qc, dt
    real(kind=c_real), intent(inout) :: qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend

    call cloud_water_conservation(qc,dt,qc2qr_autoconv_tend,qc2qr_accret_tend,qc2qi_collect_tend,qc2qi_hetero_freeze_tend,qc2qr_ice_shed_tend,qc2qi_berg_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend)
  end subroutine cloud_water_conservation_c

  subroutine rain_water_conservation_c(qr,qc2qr_autoconv_tend,qc2qr_accret_tend,qi2qr_melt_tend,qc2qr_ice_shed_tend,dt,    &
    qr2qv_evap_tend,qr2qi_collect_tend,qr2qi_immers_freeze_tend) bind(C)
    use micro_p3, only: rain_water_conservation

    real(kind=c_real), value, intent(in) :: qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt
    real(kind=c_real), intent(inout) :: qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend

    call rain_water_conservation(qr,qc2qr_autoconv_tend,qc2qr_accret_tend,qi2qr_melt_tend,qc2qr_ice_shed_tend,dt,qr2qv_evap_tend,qr2qi_collect_tend,qr2qi_immers_freeze_tend)
  end subroutine rain_water_conservation_c

  subroutine rain_self_collection_c(rho, qr_incld, nr_incld, nr_selfcollect_tend) bind(C)
    use micro_p3, only: rain_self_collection

    real(kind=c_real), value, intent(in) :: rho, qr_incld, nr_incld
    real(kind=c_real), intent(out) :: nr_selfcollect_tend

    call rain_self_collection(rho, qr_incld, nr_incld, nr_selfcollect_tend)
  end subroutine rain_self_collection_c

  subroutine ice_water_conservation_c(qi,qv2qi_vapdep_tend,qv2qi_nucleat_tend,qc2qi_berg_tend,qr2qi_collect_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,qc2qi_hetero_freeze_tend,dt,    &
    qi2qv_sublim_tend,qi2qr_melt_tend) bind(C)
    use micro_p3, only: ice_water_conservation

    real(kind=c_real), value, intent(in) :: qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, qc2qi_berg_tend, dt
    real(kind=c_real), intent(inout) :: qi2qv_sublim_tend, qi2qr_melt_tend

    call ice_water_conservation(qi,qv2qi_vapdep_tend,qv2qi_nucleat_tend,qr2qi_collect_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,qc2qi_hetero_freeze_tend,qc2qi_berg_tend,dt,qi2qv_sublim_tend,qi2qr_melt_tend)
  end subroutine ice_water_conservation_c

  subroutine get_cloud_dsd2_c(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1) bind(C)
    use micro_p3, only: get_cloud_dsd2
    use micro_p3_utils, only: dnu

    !arguments:
    real(kind=c_real), value, intent(in)        :: qc,rho
    real(kind=c_real), intent(inout)            :: nc
    real(kind=c_real), intent(out)              :: mu_c,nu,lamc,cdist,cdist1

    call get_cloud_dsd2(qc,nc,mu_c,rho,nu,dnu,lamc,cdist,cdist1)
  end subroutine get_cloud_dsd2_c

  subroutine get_rain_dsd2_c(qr,nr,mu_r,lamr,cdistr,logn0r) bind(C)
    use micro_p3, only: get_rain_dsd2

    !arguments:
    real(kind=c_real), value, intent(in) :: qr
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r

    call get_rain_dsd2(qr,nr,mu_r,lamr,cdistr,logn0r)
  end subroutine get_rain_dsd2_c

  subroutine calc_rime_density_c(T_atm,rhofaci,table_val_qi_fallspd,acn,lamc,mu_c,qc_incld,qc2qi_collect_tend, &
                                 vtrmi1,rho_qm_cloud) bind(C)

      use micro_p3, only: calc_rime_density
      real(kind=c_real), value, intent(in) :: T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend
      real(kind=c_real), intent(out) :: vtrmi1, rho_qm_cloud

      call calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend, vtrmi1, rho_qm_cloud)
  end subroutine calc_rime_density_c

  subroutine cldliq_immersion_freezing_c(T_atm,lamc,mu_c,cdist1,qc_incld,inv_qc_relvar,qc2qi_hetero_freeze_tend,nc2ni_immers_freeze_tend) bind(C)

      use micro_p3, only: cldliq_immersion_freezing
      real(kind=c_real), value, intent(in) :: T_atm, lamc, mu_c, cdist1, qc_incld,inv_qc_relvar
      real(kind=c_real), intent(out) :: qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend

      call cldliq_immersion_freezing(T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar, qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend)
  end subroutine cldliq_immersion_freezing_c

  subroutine rain_immersion_freezing_c(T_atm,lamr,mu_r,cdistr,qr_incld,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend) bind(C)

      use micro_p3, only: rain_immersion_freezing
      real(kind=c_real), value, intent(in) :: T_atm, lamr, mu_r, cdistr, qr_incld
      real(kind=c_real), intent(out) :: qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend

      call rain_immersion_freezing(T_atm, lamr, mu_r, cdistr, qr_incld, qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend)
  end subroutine rain_immersion_freezing_c

  subroutine droplet_self_collection_c(rho,inv_rho,qc_incld,mu_c,nu,nc2nr_autoconv_tend,nc_selfcollect_tend) bind(C)

      use micro_p3, only: droplet_self_collection
      real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend
      real(kind=c_real), intent(out) :: nc_selfcollect_tend

      call droplet_self_collection(rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend, nc_selfcollect_tend)
  end subroutine droplet_self_collection_c

  subroutine cloud_rain_accretion_c(rho,inv_rho,qc_incld,nc_incld,qr_incld,inv_qc_relvar,qc2qr_accret_tend,nc_accret_tend) bind(C)

      use micro_p3, only: cloud_rain_accretion
      real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, nc_incld, qr_incld,inv_qc_relvar
      real(kind=c_real), intent(out) :: qc2qr_accret_tend, nc_accret_tend

      call cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar, qc2qr_accret_tend, nc_accret_tend)
  end subroutine cloud_rain_accretion_c

  subroutine cloud_water_autoconversion_c(rho,qc_incld,nc_incld,inv_qc_relvar,qc2qr_autoconv_tend,nc2nr_autoconv_tend,ncautr) bind(C)

      use micro_p3, only: cloud_water_autoconversion
      real(kind=c_real), value, intent(in) :: rho, qc_incld, nc_incld,inv_qc_relvar
      real(kind=c_real), intent(inout) :: qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr

      call cloud_water_autoconversion(rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr)
  end subroutine cloud_water_autoconversion_c

  subroutine impose_max_total_ni_c(ni_local, max_total_ni, inv_rho_local) bind(C)
    use micro_p3, only: impose_max_total_ni

    real(kind=c_real), intent(inout) :: ni_local
    real(kind=c_real), value, intent(in) :: max_total_ni, inv_rho_local

    call impose_max_total_ni(ni_local, max_total_ni, inv_rho_local)
  end subroutine impose_max_total_ni_c

  subroutine calc_first_order_upwind_step_c(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, num_arrays, fluxes, vs, qnx) bind(C)
    use micro_p3, only: calc_first_order_upwind_step, realptr

    !arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
    real(kind=c_real), value, intent(in) :: dt_sub
    real(kind=c_real), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dz
    type(c_ptr), intent(in), dimension(num_arrays) :: fluxes, vs, qnx

    type(realptr), dimension(num_arrays) :: fluxes_f, vs_f, qnx_f
    integer :: i

    do i = 1, num_arrays
       call c_f_pointer(fluxes(i), fluxes_f(i)%p, [(kte-kts)+1])
       call c_f_pointer(vs(i),     vs_f(i)%p,     [(kte-kts)+1])
       call c_f_pointer(qnx(i),    qnx_f(i)%p ,   [(kte-kts)+1])
    end do

    call calc_first_order_upwind_step(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, num_arrays, fluxes_f, vs_f, qnx_f)

  end subroutine calc_first_order_upwind_step_c

  subroutine generalized_sedimentation_c(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnx) bind(C)
    use micro_p3, only: generalized_sedimentation, realptr

    ! arguments
    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
    integer(kind=c_int), intent(inout) :: k_qxbot
    real(kind=c_real), value, intent(in) :: Co_max
    real(kind=c_real), intent(inout) :: dt_left, prt_accum
    real(kind=c_real), dimension(kts:kte), intent(in) :: inv_dz, inv_rho, rho
    type(c_ptr), intent(in), dimension(num_arrays) :: vs, fluxes, qnx

    type(realptr), dimension(num_arrays) :: fluxes_f, vs_f, qnx_f
    integer :: i

    do i = 1, num_arrays
       call c_f_pointer(fluxes(i), fluxes_f(i)%p, [(kte-kts)+1])
       call c_f_pointer(vs(i),     vs_f(i)%p,     [(kte-kts)+1])
       call c_f_pointer(qnx(i),    qnx_f(i)%p ,   [(kte-kts)+1])
    end do

    call generalized_sedimentation(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dz, inv_rho, rho, num_arrays, vs_f, fluxes_f, qnx_f)

  end subroutine generalized_sedimentation_c

  subroutine cloud_sedimentation_c(kts,kte,ktop,kbot,kdir,   &
       qc_incld,rho,inv_rho,cld_frac_l,acn,inv_dz,&
       dt,inv_dt,do_predict_nc, &
       qc, nc, nc_incld,mu_c,lamc,precip_liq_surf,qc_tend,nc_tend) bind(C)
    use micro_p3, only: cloud_sedimentation, dnu

    ! arguments
    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: cld_frac_l
    real(kind=c_real), intent(in), dimension(kts:kte) :: acn
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dz

    real(kind=c_real),    value, intent(in) :: dt
    real(kind=c_real),    value, intent(in) :: inv_dt
    logical(kind=c_bool), value, intent(in) :: do_predict_nc

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_c
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamc
    real(kind=c_real), intent(inout) :: precip_liq_surf
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_tend

    call cloud_sedimentation(kts,kte,ktop,kbot,kdir,   &
         qc_incld,rho,inv_rho,cld_frac_l,acn,inv_dz,&
         dt,inv_dt,dnu,do_predict_nc, &
         qc, nc, nc_incld,mu_c,lamc,precip_liq_surf,qc_tend,nc_tend)

  end subroutine cloud_sedimentation_c

  subroutine ice_sedimentation_c(kts,kte,ktop,kbot,kdir,    &
       rho,inv_rho,rhofaci,cld_frac_i,inv_dz,dt,inv_dt,  &
       qi,qi_incld,ni,qm,qm_incld,bm,bm_incld,ni_incld,precip_ice_surf,qi_tend,ni_tend) bind(C)
    use micro_p3, only: ice_sedimentation

    ! arguments
    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofaci
    real(kind=c_real), intent(in), dimension(kts:kte) :: cld_frac_i
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dz
    real(kind=c_real), value, intent(in) :: dt, inv_dt

    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qi
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: ni
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qm_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: bm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: bm_incld

    real(kind=c_real), intent(inout) :: precip_ice_surf
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni_tend

    call ice_sedimentation(kts,kte,ktop,kbot,kdir,    &
         rho,inv_rho,rhofaci,cld_frac_i,inv_dz,dt,inv_dt,  &
         qi,qi_incld,ni,qm,qm_incld,bm,bm_incld,ni_incld,precip_ice_surf,qi_tend,ni_tend)

  end subroutine ice_sedimentation_c

  subroutine rain_sedimentation_c(kts,kte,ktop,kbot,kdir,   &
       qr_incld,rho,inv_rho,rhofacr,cld_frac_r,inv_dz,dt,inv_dt,  &
       qr,nr,nr_incld,mu_r,lamr,precip_liq_surf,precip_liq_flux,qr_tend,nr_tend) bind(C)
    use micro_p3, only: rain_sedimentation

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofacr
    real(kind=c_real), intent(in), dimension(kts:kte) :: cld_frac_r
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dz
    real(kind=c_real), value, intent(in) :: dt, inv_dt

    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: nr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_r
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamr
    real(kind=c_real), intent(inout) :: precip_liq_surf
    real(kind=c_real), intent(inout), dimension(kts:kte+1) :: precip_liq_flux
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_tend

    call rain_sedimentation(kts,kte,ktop,kbot,kdir,   &
         qr_incld,rho,inv_rho,rhofacr,cld_frac_r,inv_dz,dt,inv_dt,  &
         qr,nr,nr_incld,mu_r,lamr,precip_liq_surf,precip_liq_flux,qr_tend,nr_tend)

  end subroutine rain_sedimentation_c

  subroutine calc_bulk_rho_rime_c(qi_tot, qi_rim, bi_rim, rho_rime) bind(C)
    use micro_p3, only: calc_bulkRhoRime

    ! arguments:
    real(kind=c_real),   value, intent(in)  :: qi_tot
    real(kind=c_real),   intent(inout) :: qi_rim, bi_rim
    real(kind=c_real),   intent(out) :: rho_rime

    call calc_bulkRhoRime(qi_tot, qi_rim, bi_rim, rho_rime)
  end subroutine calc_bulk_rho_rime_c

  subroutine homogeneous_freezing_c(kts,kte,ktop,kbot,kdir,T_atm,inv_exner,latent_heat_fusion,    &
   qc,nc,qr,nr,qi,ni,qm,bm,th_atm) bind(C)
    use micro_p3, only: homogeneous_freezing

    ! arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir
    real(kind=c_real), intent(in), dimension(kts:kte) :: T_atm
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_exner
    real(kind=c_real), intent(in), dimension(kts:kte) :: latent_heat_fusion

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: bm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: th_atm

    call homogeneous_freezing(kts,kte,ktop,kbot,kdir,T_atm,inv_exner,latent_heat_fusion,    &
         qc,nc,qr,nr,qi,ni,qm,bm,th_atm)
  end subroutine homogeneous_freezing_c

  subroutine compute_rain_fall_velocity_c(qr_incld, rhofacr, nr_incld, mu_r, lamr, V_qr, V_nr) bind(C)
    use micro_p3, only: compute_rain_fall_velocity

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, rhofacr
    real(kind=c_real), intent(inout) :: nr_incld
    real(kind=c_real), intent(out) :: mu_r, lamr, V_qr, V_nr

    call compute_rain_fall_velocity(qr_incld, rhofacr, nr_incld, mu_r, lamr, V_qr, V_nr)
  end subroutine compute_rain_fall_velocity_c

subroutine  update_prognostic_ice_c(qc2qi_hetero_freeze_tend,qc2qi_collect_tend,qc2qr_ice_shed_tend,nc_collect_tend,nc2ni_immers_freeze_tend,ncshdc,qr2qi_collect_tend,nr_collect_tend,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend,nr_ice_shed_tend, &
       qi2qr_melt_tend,ni2nr_melt_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend,qv2qi_nucleat_tend,ni_nucleat_tend,ni_selfcollect_tend,ni_sublim_tend,qc2qi_berg_tend,inv_exner,latent_heat_sublim,latent_heat_fusion,do_predict_nc,log_wetgrowth, &
       dt,nmltratio,rho_qm_cloud,th_atm,qv,qi,ni,qm,bm,qc,nc,qr,nr) bind(C)
    use micro_p3, only: update_prognostic_ice

    ! arguments
    real(kind=c_real), value, intent(in) :: qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend, &
         qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend, ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, inv_exner, &
         latent_heat_fusion, latent_heat_sublim, dt, nmltratio, rho_qm_cloud

    logical(kind=c_bool), value, intent(in) :: do_predict_nc, log_wetgrowth

    real(kind=c_real), intent(inout) :: th_atm, qv, qc, nc, qr, nr, qi, ni, qm, bm

    call update_prognostic_ice(qc2qi_hetero_freeze_tend,qc2qi_collect_tend,qc2qr_ice_shed_tend,nc_collect_tend,nc2ni_immers_freeze_tend,ncshdc,qr2qi_collect_tend,nr_collect_tend,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend,nr_ice_shed_tend, &
         qi2qr_melt_tend,ni2nr_melt_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend,qv2qi_nucleat_tend,ni_nucleat_tend,ni_selfcollect_tend,ni_sublim_tend,qc2qi_berg_tend,inv_exner,latent_heat_sublim,latent_heat_fusion,do_predict_nc,log_wetgrowth, &
         dt,nmltratio,rho_qm_cloud,th_atm,qv,qi,ni,qm,bm,qc,nc,qr,nr)

  end subroutine update_prognostic_ice_c

  subroutine get_time_space_phys_variables_c(T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt, &
                                             ab, abi, kap, eii) bind(C)
    use micro_p3, only: get_time_space_phys_variables

    !arguments
    real(kind=c_real), value, intent(in) :: T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i
    real(kind=c_real), intent(out) :: mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii

    call get_time_space_phys_variables(T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt, &
                                       ab, abi, kap, eii)
  end subroutine get_time_space_phys_variables_c

  subroutine ice_cldliq_collection_c(rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld, ni_incld, &
                                     nc_incld, qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc) bind(C)
    use micro_p3, only: ice_cldliq_collection

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, temp, rhofaci, table_val_qc2qi_collect
    real(kind=c_real), value, intent(in) :: qi_incld, qc_incld, ni_incld, nc_incld
    real(kind=c_real), intent(out) :: qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc

    call ice_cldliq_collection(rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld, ni_incld, &
                               nc_incld, qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc)
  end subroutine ice_cldliq_collection_c

  subroutine ice_rain_collection_c(rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect, &
                                   qi_incld, ni_incld, qr_incld, qr2qi_collect_tend, nr_collect_tend) bind(C)
    use micro_p3, only: ice_rain_collection

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect
    real(kind=c_real), value, intent(in) :: qi_incld, ni_incld, qr_incld
    real(kind=c_real), intent(out) :: qr2qi_collect_tend, nr_collect_tend

    call ice_rain_collection(rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect,  &
                             qi_incld, ni_incld, qr_incld, qr2qi_collect_tend, nr_collect_tend)
  end subroutine ice_rain_collection_c

  subroutine ice_self_collection_c(rho, rhofaci, table_val_ni_self_collect, eii, qm_incld, &
                                   qi_incld, ni_incld, ni_selfcollect_tend) bind(C)
    use micro_p3, only: ice_self_collection

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, rhofaci, table_val_ni_self_collect, eii, qm_incld
    real(kind=c_real), value, intent(in) :: qi_incld, ni_incld
    real(kind=c_real), intent(out) :: ni_selfcollect_tend

    call ice_self_collection(rho, rhofaci, table_val_ni_self_collect, eii, qm_incld, &
                             qi_incld, ni_incld, ni_selfcollect_tend)
  end subroutine ice_self_collection_c

  subroutine evaporate_rain_c(qr_incld,qc_incld,nr_incld,qi_incld, &
       cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i, &
       ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt,&
       qr2qv_evap_tend,nr_evap_tend) bind(C)
    use micro_p3, only: evaporate_rain

    ! arguments
    real(kind=c_real), value, intent(in) :: qr_incld,qc_incld,nr_incld,qi_incld, &
         cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i, &
         ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt

    real(kind=c_real), intent(out) :: qr2qv_evap_tend, nr_evap_tend

    call evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld, &
       cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i, &
       ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt,&
       qr2qv_evap_tend,nr_evap_tend)
  end subroutine evaporate_rain_c

  subroutine  update_prognostic_liquid_c(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend,nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, &
       qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend, do_predict_nc, do_prescribed_CCN, inv_rho, inv_exner, latent_heat_vapor, dt, th_atm, qv, qc, nc, qr, nr) bind(C)
    use micro_p3, only: update_prognostic_liquid

    ! arguments
    real(kind=c_real), value, intent(in) :: qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, &
         qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend

    logical(kind=c_bool), value, intent(in) :: do_predict_nc
    logical(kind=c_bool), value, intent(in) :: do_prescribed_CCN

    real(kind=c_real), value, intent(in) :: inv_rho, inv_exner, latent_heat_vapor, dt

    real(kind=c_real), intent(inout) :: th_atm, qv, qc, nc, qr, nr

    call update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend,nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, &
       qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend, do_predict_nc, do_prescribed_CCN, inv_rho, inv_exner, latent_heat_vapor, dt, th_atm, qv, qc, nc, qr, nr)

  end subroutine update_prognostic_liquid_c

  subroutine ice_deposition_sublimation_c(qi_incld, ni_incld, t_atm, qv_sat_l, qv_sat_i, epsi, abi, qv, inv_dt, qidep, qi2qv_sublim_tend, ni_sublim_tend, qiberg) bind(C)
    use micro_p3, only : ice_deposition_sublimation

    real(kind=c_real) , value, intent(in) :: qi_incld, ni_incld, t_atm, qv_sat_l, qv_sat_i, epsi, abi, qv, inv_dt
    real(kind=c_real) , intent(out) :: qidep, qi2qv_sublim_tend, ni_sublim_tend, qiberg

    call ice_deposition_sublimation(qi_incld, ni_incld, t_atm, qv_sat_l, qv_sat_i, epsi, abi, qv, inv_dt, qidep, qi2qv_sublim_tend, ni_sublim_tend, qiberg)
  end subroutine ice_deposition_sublimation_c

  subroutine ice_relaxation_timescale_c(rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt,   &
                                        dv, mu, sc, qi_incld, ni_incld, &
                                        epsi, epsi_tot) bind(C)
    use micro_p3, only: calc_ice_relaxation_timescale

    ! arguments
    real(kind=c_real), value, intent(in) :: rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, &
                                            dv, mu, sc, qi_incld, ni_incld
    real(kind=c_real), intent(out)   :: epsi
    real(kind=c_real), intent(inout) :: epsi_tot

    call calc_ice_relaxation_timescale(rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt,   &
                                       dv, mu, sc, qi_incld, ni_incld, &
                                       epsi, epsi_tot)
  end subroutine ice_relaxation_timescale_c

  subroutine calc_liq_relaxation_timescale_c(rho, f1r, f2r, dv, mu, sc, mu_r, &
                                             lamr, cdistr, cdist, qr_incld,   &
                                             qc_incld, epsr, epsc) bind(C)
    use micro_p3, only: calc_liq_relaxation_timescale

    ! arguments
    real(kind=c_real), value, intent(in) :: rho,f1r,f2r,dv,mu,sc,mu_r,lamr, &
                                            cdistr,cdist,qr_incld,qc_incld
    real(kind=c_real), intent(out) :: epsr
    real(kind=c_real), intent(out) :: epsc

    call calc_liq_relaxation_timescale(rho,f1r,f2r,dv,mu,sc,mu_r,lamr,      &
                                       cdistr,cdist,qr_incld,qc_incld,epsr, &
                                       epsc)
  end subroutine calc_liq_relaxation_timescale_c

  subroutine ice_nucleation_c(temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt, &
                              do_predict_nc, do_prescribed_CCN, qv2qi_nucleat_tend, ni_nucleat_tend) bind(C)
    use micro_p3, only: ice_nucleation

    ! arguments
    real(kind=c_real), value, intent(in) :: temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt
    logical(c_bool), value, intent(in) :: do_predict_nc
    logical(c_bool), value, intent(in) :: do_prescribed_CCN

    real(kind=c_real), intent(inout) :: qv2qi_nucleat_tend, ni_nucleat_tend

    call ice_nucleation(temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt, &
                        do_predict_nc, do_prescribed_CCN, qv2qi_nucleat_tend, ni_nucleat_tend)
 end subroutine ice_nucleation_c

 subroutine ice_melting_c(rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion, &
    dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend) bind(C)
    use micro_p3, only: ice_melting

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld
    real(kind=c_real), intent(out) :: qi2qr_melt_tend,ni2nr_melt_tend

    call ice_melting(rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend)

  end subroutine ice_melting_c

 subroutine ice_cldliq_wet_growth_c(rho, temp, pres, rhofaci, table_val_qi2qr_melting, &
                                    table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld,  &
                                    qi_incld, ni_incld, qr_incld, &
                                    log_wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend) bind(C)
   use micro_p3, only: ice_cldliq_wet_growth

   ! argmens
   real(kind=c_real), value, intent(in) :: rho, temp ,pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, &
                                           kap, mu, sc, qv, qc_incld, qi_incld, ni_incld,qr_incld
   logical(kind=c_bool), intent(inout) :: log_wetgrowth
   real(kind=c_real), intent(inout) :: qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend

   call ice_cldliq_wet_growth(rho, temp, pres, rhofaci, table_val_qi2qr_melting, &
                              table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld, &
                              qi_incld, ni_incld, qr_incld, &
                              log_wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend)
 end subroutine ice_cldliq_wet_growth_c

 subroutine get_latent_heat_c(its,ite,kts,kte,v,s,f) bind(C)
   use micro_p3, only: get_latent_heat

   ! arguments
   integer(kind=c_int), intent(in), value :: its, ite, kts, kte
   real(kind=c_real), dimension(its:ite, kts:kte), intent(out) :: v, s, f

   call get_latent_heat(its,ite,kts,kte,v,s,f)
 end subroutine get_latent_heat_c

 function subgrid_variance_scaling_c(relvar,expon) result(res) bind(C)
   use micro_p3, only: subgrid_variance_scaling

   ! arguments
   real(kind=c_real), value, intent(in) :: relvar,expon
   real(kind=c_real) :: res

   res = subgrid_variance_scaling(relvar,expon)
   return
 end function subgrid_variance_scaling_c

 subroutine check_values_c(qv, temp, kts, kte, timestepcount, &
                           force_abort, source_ind, col_loc) bind(C)
   use micro_p3, only: check_values

   ! argmens
   real(kind=c_real), intent(in) :: qv(kts:kte), temp(kts:kte), col_loc(3)
   integer(kind=c_int), value, intent(in) :: kts, kte, timestepcount, source_ind
   logical(kind=c_bool), value, intent(in) :: force_abort

   call check_values(qv,Temp,kts,kte,timestepcount,force_abort,source_ind,col_loc)
 end subroutine check_values_c

 subroutine calculate_incloud_mixingratios_c(qc, qr, qi, qm, nc, nr, ni, bm,   &
                                             inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,              &
                                             qc_incld, qr_incld, qi_incld, qm_incld, &
                                             nc_incld, nr_incld, ni_incld, bm_incld) bind(C)
   use micro_p3, only: calculate_incloud_mixingratios

   ! argumens
   real(kind=c_real), value, intent(in) :: qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
   real(kind=c_real), intent(inout) :: qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld

   call calculate_incloud_mixingratios(qc, qr, qi, qm, nc, nr, ni, bm,   &
                                       inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,              &
                                       qc_incld, qr_incld, qi_incld, qm_incld, &
                                       nc_incld, nr_incld, ni_incld, bm_incld)
 end subroutine calculate_incloud_mixingratios_c

 subroutine p3_main_part1_c(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, &
       pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, nccn_prescribed, &
       T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm, qc, nc, qr, nr, &
       qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, qm_incld, &
       nc_incld, nr_incld, ni_incld, bm_incld, is_nucleat_possible, is_hydromet_present) bind(C)

   use micro_p3, only: p3_main_part1

   ! arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: do_predict_nc, do_prescribed_CCN
   real(kind=c_real), value, intent(in) :: dt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, &
        inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, nccn_prescribed

   real(kind=c_real), intent(inout), dimension(kts:kte) :: T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, &
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, &
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld

   logical(kind=c_bool), intent(out) :: is_nucleat_possible, is_hydromet_present

   call p3_main_part1(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, &
        pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, nccn_prescribed, &
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm, qc, nc, qr, nr, &
        qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, qm_incld, &
        nc_incld, nr_incld, ni_incld, bm_incld, is_nucleat_possible, is_hydromet_present)

 end subroutine p3_main_part1_c

 subroutine p3_main_part2_c(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, &
       pres, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, &
       ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev, &
       T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofaci, acn, qv, th_atm, qc, nc, qr, nr, qi, ni, &
       qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, &
       ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, &
       nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
       prctot, is_hydromet_present) bind(C)

   use micro_p3, only: p3_main_part2

   !arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: do_predict_nc, do_prescribed_CCN
   real(kind=c_real), value, intent(in) :: dt, inv_dt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, inv_exner, inv_cld_frac_l, inv_cld_frac_i, &
        inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev

   real(kind=c_real), intent(inout), dimension(kts:kte) :: T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofaci, acn, &
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, &
        qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, &
        cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, vap_liq_exchange, &
        vap_ice_exchange, liq_ice_exchange, pratot, prctot

   logical(kind=c_bool), intent(out) :: is_hydromet_present

   ! throwaway
   real(kind=c_real), dimension(kts:kte,49) :: p3_tend_out

   call p3_main_part2(kts, kte, kbot, ktop, kdir, do_predict_nc, do_prescribed_CCN, dt, inv_dt, &
        pres, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev, t_prev, &
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofaci, acn, qv, th_atm, qc, nc, qr, nr, qi, ni, &
        qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, &
        ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, &
        nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
        prctot, p3_tend_out, is_hydromet_present)

 end subroutine p3_main_part2_c

 subroutine p3_main_part3_c(kts, kte, kbot, ktop, kdir, &
      inv_exner, cld_frac_l, cld_frac_r, cld_frac_i, &
      rho, inv_rho, rhofaci, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, &
      ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc) bind(C)

   use micro_p3, only: p3_main_part3

   ! args

   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   real(kind=c_real), intent(in), dimension(kts:kte) :: inv_exner, cld_frac_l, cld_frac_r, cld_frac_i
   real(kind=c_real), intent(inout), dimension(kts:kte) :: rho, inv_rho, rhofaci, &
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
        mu_c, nu, lamc, mu_r, &
        lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc

   call p3_main_part3(kts, kte, kbot, ktop, kdir, &
        inv_exner, cld_frac_l, cld_frac_r, cld_frac_i, &
        rho, inv_rho, rhofaci, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
        mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc)

 end subroutine p3_main_part3_c

  subroutine ice_supersat_conservation_c(qidep, qinuc, cld_frac_i, qv, qv_sat_i, latent_heat_sublim, t_atm, dt, qi2qv_sublim_tend, qr2qv_evap_tend) bind(C)
    use micro_p3, only : ice_supersat_conservation

    real(kind=c_real) , intent(inout) :: qidep, qinuc
    real(kind=c_real) , value, intent(in) :: cld_frac_i, qv, qv_sat_i, latent_heat_sublim, t_atm, dt, qi2qv_sublim_tend, qr2qv_evap_tend

    call ice_supersat_conservation(qidep, qinuc, cld_frac_i, qv, qv_sat_i, latent_heat_sublim, t_atm, dt, qi2qv_sublim_tend, qr2qv_evap_tend)
  end subroutine ice_supersat_conservation_c
  subroutine nc_conservation_c(nc, nc_selfcollect_tend, dt, nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend) bind(C)
    use micro_p3, only : nc_conservation

    real(kind=c_real) , value, intent(in) :: nc, nc_selfcollect_tend, dt
    real(kind=c_real) , intent(inout) :: nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend

    call nc_conservation(nc, nc_selfcollect_tend, dt, nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend)
  end subroutine nc_conservation_c
  subroutine nr_conservation_c(nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, dt, nmltratio, nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend) bind(C)
    use micro_p3, only : nr_conservation

    real(kind=c_real) , value, intent(in) :: nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, dt, nmltratio
    real(kind=c_real) , intent(inout) :: nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend

    call nr_conservation(nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, dt, nmltratio, nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend)
  end subroutine nr_conservation_c
  subroutine ni_conservation_c(ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, dt, ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend) bind(C)
    use micro_p3, only : ni_conservation

    real(kind=c_real) , value, intent(in) :: ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, dt
    real(kind=c_real) , intent(inout) :: ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend

    call ni_conservation(ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, dt, ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend)
  end subroutine ni_conservation_c
  subroutine prevent_liq_supersaturation_c(pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend) bind(C)
    use micro_p3, only : prevent_liq_supersaturation

    real(kind=c_real) , value, intent(in) :: pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc
    real(kind=c_real) , intent(inout) :: qi2qv_sublim_tend, qr2qv_evap_tend

    call prevent_liq_supersaturation(pres, t_atm, qv, latent_heat_vapor, latent_heat_sublim, dt, qidep, qinuc, qi2qv_sublim_tend, qr2qv_evap_tend)
  end subroutine prevent_liq_supersaturation_c
end module p3_iso_c
