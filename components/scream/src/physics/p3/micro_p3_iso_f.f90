module micro_p3_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from micro_p3 fortran to scream c++.
!

interface

  subroutine find_lookuptable_indices_1a_f(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qi,ni,qm,rhop) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real),   value, intent(in)  :: qi,ni,qm,rhop
  end subroutine find_lookuptable_indices_1a_f

  subroutine find_lookuptable_indices_1b_f(dumj,dum3,qr,nr) bind(C)
    use iso_c_binding

    integer(kind=c_int), intent(out) :: dumj
    real(kind=c_real),   intent(out) :: dum3
    real(kind=c_real),   value, intent(in) :: qr, nr
  end subroutine find_lookuptable_indices_1b_f

  subroutine access_lookup_table_f(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj, dumii, dumi, index
    real(kind=c_real),   value, intent(in) :: dum1, dum4, dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_f

  subroutine access_lookup_table_coll_f(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj,dumii,dumj,dumi,index
    real(kind=c_real),   value, intent(in) :: dum1,dum3,dum4,dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_coll_f

  subroutine back_to_cell_average_f(cld_frac_l,cld_frac_r,cld_frac_i, qc2qr_accret_tend,qr2qv_evap_tend,qc2qr_autoconv_tend,&
    nc_accret_tend,nc_selfcollect_tend,nc2nr_autoconv_tend,nr_selfcollect_tend,nr_evap_tend,ncautr,qi2qv_sublim_tend,nr_ice_shed_tend,qc2qi_hetero_freeze_tend,&
    qr2qi_collect_tend,qc2qr_ice_shed_tend,qi2qr_melt_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,ni2nr_melt_tend, &
    nc_collect_tend,ncshdc,nc2ni_immers_freeze_tend,nr_collect_tend,ni_selfcollect_tend,&
    qv2qi_vapdep_tend,nr2ni_immers_freeze_tend,ni_sublim_tend,qv2qi_nucleat_tend,ni_nucleat_tend,qc2qi_berg_tend) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: cld_frac_l, cld_frac_r, cld_frac_i
    real(kind=c_real), intent(inout) :: qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend, nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend,  &
                                        nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend,  &
                                        nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, &
                                        qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend, ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend,&
                                        ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend,  &
                                        qc2qi_berg_tend
  end subroutine back_to_cell_average_f

  subroutine prevent_ice_overdepletion_f(pres,t,qv,latent_heat_sublim,inv_dt,    &
     qv2qi_vapdep_tend,qi2qv_sublim_tend) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: pres, t, qv, latent_heat_sublim, inv_dt
    real(kind=c_real), intent(inout) :: qv2qi_vapdep_tend, qi2qv_sublim_tend

  end subroutine prevent_ice_overdepletion_f

  subroutine cloud_water_conservation_f(qc,dt,qc2qr_autoconv_tend,qc2qr_accret_tend,qc2qi_collect_tend,qc2qi_hetero_freeze_tend,qc2qr_ice_shed_tend,     &
    qc2qi_berg_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qc, dt
    real(kind=c_real), intent(inout) :: qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, &
        qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend
  end subroutine cloud_water_conservation_f

  subroutine rain_water_conservation_f(qr,qc2qr_autoconv_tend,qc2qr_accret_tend,qi2qr_melt_tend,qc2qr_ice_shed_tend,dt,    &
    qr2qv_evap_tend,qr2qi_collect_tend,qr2qi_immers_freeze_tend) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt
    real(kind=c_real), intent(inout) :: qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend
  end subroutine rain_water_conservation_f

  subroutine ice_water_conservation_f(qi,qv2qi_vapdep_tend,qv2qi_nucleat_tend,qc2qi_berg_tend,qr2qi_collect_tend,qc2qi_collect_tend,qr2qi_immers_freeze_tend,qc2qi_hetero_freeze_tend,dt,    &
    qi2qv_sublim_tend,qi2qr_melt_tend) bind(C)
    use iso_c_binding

    real(kind=c_real), value, intent(in) :: qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qr2qi_collect_tend, qc2qi_collect_tend, &
        qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, qc2qi_berg_tend, dt
    real(kind=c_real), intent(inout) :: qi2qv_sublim_tend, qi2qr_melt_tend

  end subroutine ice_water_conservation_f

  subroutine get_cloud_dsd2_f(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1,cld_frac_l) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)         :: qc,rho,cld_frac_l
    real(kind=c_real), intent(inout)             :: nc
    real(kind=c_real), intent(out)               :: mu_c,nu,lamc,cdist,cdist1
  end subroutine get_cloud_dsd2_f

  subroutine get_rain_dsd2_f(qr,nr,mu_r,lamr,cdistr,logn0r,cld_frac_r) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: qr,cld_frac_r
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r
  end subroutine get_rain_dsd2_f

  subroutine calc_rime_density_f(t,rhofaci,table_val_qi_fallspd,acn,lamc,mu_c,qc_incld,qc2qi_collect_tend,vtrmi1,rho_qm_cloud) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend
    real(kind=c_real), intent(out) :: vtrmi1, rho_qm_cloud
  end subroutine calc_rime_density_f

  subroutine cldliq_immersion_freezing_f(t,lamc,mu_c,cdist1,qc_incld,inv_qc_relvar,qc2qi_hetero_freeze_tend,nc2ni_immers_freeze_tend) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar
    real(kind=c_real), intent(out) :: qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend
  end subroutine cldliq_immersion_freezing_f

  subroutine rain_immersion_freezing_f(t,lamr,mu_r,cdistr,qr_incld,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: t, lamr, mu_r, cdistr, qr_incld
    real(kind=c_real), intent(out) :: qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend
  end subroutine rain_immersion_freezing_f

  subroutine droplet_self_collection_f(rho,inv_rho,qc_incld,mu_c,nu,nc2nr_autoconv_tend,nc_selfcollect_tend) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend
    real(kind=c_real), intent(out) :: nc_selfcollect_tend
  end subroutine droplet_self_collection_f

  subroutine cloud_rain_accretion_f(rho,inv_rho,qc_incld,nc_incld,qr_incld,inv_qc_relvar,qc2qr_accret_tend,nc_accret_tend) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar
    real(kind=c_real), intent(out) :: qc2qr_accret_tend, nc_accret_tend
  end subroutine cloud_rain_accretion_f

  subroutine cloud_water_autoconversion_f(rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, qc_incld, nc_incld, inv_qc_relvar
    real(kind=c_real), intent(inout) :: qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr
  end subroutine cloud_water_autoconversion_f

  subroutine rain_self_collection_f(rho, qr_incld, nr_incld, nr_selfcollect_tend) bind(C)
    use iso_c_binding

    !arguments;
    real(kind=c_real), value, intent(in) :: rho, qr_incld, nr_incld
    real(kind=c_real), intent(out) :: nr_selfcollect_tend

  end subroutine rain_self_collection_f

  subroutine ice_melting_f(rho,t,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt, &
             latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho,t,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt, &
        latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld
    real(kind=c_real), intent(out) :: qi2qr_melt_tend,ni2nr_melt_tend
  end subroutine ice_melting_f

  subroutine impose_max_total_ni_f(ni_local, max_total_Ni, inv_rho_local) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), intent(inout) :: ni_local
    real(kind=c_real), value, intent(in) :: max_total_Ni, inv_rho_local
  end subroutine impose_max_total_ni_f

  subroutine calc_first_order_upwind_step_f(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, num_arrays, fluxes, vs, qnx) bind(C)
    use iso_c_binding

    !arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
    real(kind=c_real), value, intent(in) :: dt_sub
    real(kind=c_real), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dz
    type(c_ptr), intent(in), dimension(num_arrays) :: fluxes, vs, qnx
  end subroutine calc_first_order_upwind_step_f

  subroutine generalized_sedimentation_f(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dz, inv_rho, rho, num_arrays, vs, fluxes, qnx) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
    integer(kind=c_int), intent(inout) :: k_qxbot
    real(kind=c_real), value, intent(in) :: Co_max
    real(kind=c_real), intent(inout) :: dt_left, prt_accum
    real(kind=c_real), dimension(kts:kte), intent(in) :: inv_dz, inv_rho, rho

    type(c_ptr), intent(in), dimension(num_arrays) :: vs, fluxes, qnx
  end subroutine generalized_sedimentation_f

  subroutine cloud_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qc_incld,rho,inv_rho,cld_frac_l,acn,inv_dz,&
       dt,inv_dt,do_predict_nc, &
       qc, nc, nc_incld,mu_c,lamc,precip_liq_surf,qc_tend,nc_tend) bind(C)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qc_incld
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
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_c
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamc
    real(kind=c_real), intent(inout) :: precip_liq_surf
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_tend
  end subroutine cloud_sedimentation_f

  subroutine ice_sedimentation_f(kts,kte,ktop,kbot,kdir,    &
       rho,inv_rho,rhofaci,cld_frac_i,inv_dz,dt,inv_dt,  &
       qi,qi_incld,ni,qm,qm_incld,bm,bm_incld,ni_incld,precip_ice_surf,qi_tend,ni_tend) bind(C)

    use iso_c_binding

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
  end subroutine ice_sedimentation_f

  subroutine rain_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qr_incld,rho,inv_rho,rhofacr,cld_frac_r,inv_dz,dt,inv_dt,  &
       qr,nr,nr_incld,mu_r,lamr,precip_liq_surf,precip_liq_flux,qr_tend,nr_tend) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qr_incld

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofacr
    real(kind=c_real), intent(in), dimension(kts:kte) :: cld_frac_r
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dz
    real(kind=c_real), value, intent(in) :: dt, inv_dt

    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: nr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_r
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamr
    real(kind=c_real), intent(inout) :: precip_liq_surf
    real(kind=c_real), intent(inout), dimension(kts:kte+1) :: precip_liq_flux
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_tend

  end subroutine rain_sedimentation_f

  subroutine calc_bulk_rho_rime_f(qi_tot, qi_rim, bi_rim, rho_rime) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real),   value, intent(in)  :: qi_tot
    real(kind=c_real),   intent(inout) :: qi_rim, bi_rim
    real(kind=c_real),   intent(out) :: rho_rime
  end subroutine calc_bulk_rho_rime_f

  subroutine homogeneous_freezing_f(kts,kte,ktop,kbot,kdir,t,exner,latent_heat_fusion,    &
   qc,nc,qr,nr,qi,ni,qm,bm,th) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir
    real(kind=c_real), intent(in), dimension(kts:kte) :: t
    real(kind=c_real), intent(in), dimension(kts:kte) :: exner
    real(kind=c_real), intent(in), dimension(kts:kte) :: latent_heat_fusion

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: bm
    real(kind=c_real), intent(inout), dimension(kts:kte) :: th
  end subroutine homogeneous_freezing_f

  subroutine compute_rain_fall_velocity_f(qr_incld, cld_frac_r, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, cld_frac_r, rhofacr
    real(kind=c_real), intent(inout) :: nr, nr_incld
    real(kind=c_real), intent(out) :: mu_r, lamr, V_qr, V_nr
  end subroutine compute_rain_fall_velocity_f

  subroutine get_time_space_phys_variables_f(t, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt, &
       ab, abi, kap, eii) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: t, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i
    real(kind=c_real), intent(out) :: mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii
  end subroutine get_time_space_phys_variables_f

subroutine  update_prognostic_ice_f(qc2qi_hetero_freeze_tend,qc2qi_collect_tend,qc2qr_ice_shed_tend,nc_collect_tend,nc2ni_immers_freeze_tend,ncshdc, &
       qr2qi_collect_tend,nr_collect_tend,qr2qi_immers_freeze_tend,nr2ni_immers_freeze_tend,nr_ice_shed_tend, &
       qi2qr_melt_tend,ni2nr_melt_tend,qi2qv_sublim_tend,qv2qi_vapdep_tend,qv2qi_nucleat_tend,ni_nucleat_tend,ni_selfcollect_tend,ni_sublim_tend, &
       qc2qi_berg_tend,exner,latent_heat_sublim,latent_heat_fusion,do_predict_nc,log_wetgrowth, &
       dt,nmltratio,rho_qm_cloud,th,qv,qi,ni,qm,bm,qc,nc,qr,nr) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, &
         nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend, &
         qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend, ni2nr_melt_tend, &
         qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, qc2qi_berg_tend, exner, &
         latent_heat_fusion, latent_heat_sublim, dt, nmltratio, rho_qm_cloud

    logical(kind=c_bool), value, intent(in) :: do_predict_nc
    logical(kind=c_bool), value, intent(in) :: log_wetgrowth

    real(kind=c_real), intent(inout) :: th, qv, qc, nc, qr, nr, qi, ni, qm, bm

  end subroutine update_prognostic_ice_f

  subroutine ice_cldliq_collection_f(rho, t, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld, ni_incld, &
                                     nc_incld, qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, t, rhofaci, table_val_qc2qi_collect
    real(kind=c_real), value, intent(in) :: qi_incld, qc_incld, ni_incld, nc_incld
    real(kind=c_real), intent(out) :: qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc
  end subroutine ice_cldliq_collection_f

  subroutine ice_rain_collection_f(rho, t, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect, &
                                   qi_incld, ni_incld, qr_incld, qr2qi_collect_tend, nr_collect_tend) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, t, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect
    real(kind=c_real), value, intent(in) :: qi_incld, ni_incld, qr_incld
    real(kind=c_real), intent(out) :: qr2qi_collect_tend, nr_collect_tend
  end subroutine ice_rain_collection_f

  subroutine ice_self_collection_f(rho, rhofaci, table_val_ni_self_collect, eii, qm_incld, qi_incld, ni_incld, ni_selfcollect_tend) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, rhofaci, table_val_ni_self_collect, eii, qm_incld, qi_incld, ni_incld
    real(kind=c_real), intent(out) :: ni_selfcollect_tend
  end subroutine ice_self_collection_f

  subroutine evaporate_sublimate_precip_f(qr_incld, qc_incld, nr_incld, qi_incld,  cld_frac_l, cld_frac_r, qv_sat_l, ab, &
       epsr, qv, qr2qv_evap_tend, nr_evap_tend) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, qc_incld, nr_incld, qi_incld,  cld_frac_l, cld_frac_r, qv_sat_l, ab, &
         epsr, qv
    real(kind=c_real), intent(out) :: qr2qv_evap_tend, nr_evap_tend

  end subroutine evaporate_sublimate_precip_f

  subroutine update_prognostic_liquid_f(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend,nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, &
       qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend, do_predict_nc, inv_rho, exner, latent_heat_vapor, dt, th, qv, qc, nc, qr, nr) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, &
         qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend

    logical(kind=c_bool), value, intent(in) :: do_predict_nc

    real(kind=c_real), value, intent(in) :: inv_rho, exner, latent_heat_vapor, dt

    real(kind=c_real), intent(inout) :: th, qv, qc, nc, qr, nr

  end subroutine update_prognostic_liquid_f

  subroutine ice_deposition_sublimation_f(qi_incld, ni_incld, t, qv_sat_l, qv_sat_i, epsi, abi, qv, &
       qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend) bind(C)

    use iso_c_binding
    !arguments

    real(kind=c_real), value, intent(in) :: qi_incld, ni_incld, t, qv_sat_l, qv_sat_i, epsi, abi, qv
    real(kind=c_real), intent(out) :: qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend

  end subroutine ice_deposition_sublimation_f

  subroutine ice_relaxation_timescale_f(rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt,   &
                                        dv, mu, sc, qi_incld, ni_incld, &
                                        epsi, epsi_tot) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, &
                                            dv, mu, sc, qi_incld, ni_incld
    real(kind=c_real), intent(out) :: epsi
    real(kind=c_real), intent(inout) :: epsi_tot
  end subroutine ice_relaxation_timescale_f

  subroutine calc_liq_relaxation_timescale_f(rho, f1r, f2r, dv, mu, sc, mu_r, &
                                             lamr, cdistr, cdist, qr_incld,   &
                                             qc_incld, epsr, epsc) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: rho,f1r,f2r,dv,mu,sc,mu_r,lamr, &
                                            cdistr,cdist,qr_incld,qc_incld
    real(kind=c_real), intent(out) :: epsr
    real(kind=c_real), intent(out) :: epsc

  end subroutine calc_liq_relaxation_timescale_f

  subroutine ice_nucleation_f(temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt, &
                              do_predict_nc, qv2qi_nucleat_tend, ni_nucleat_tend) bind(C)
    use iso_c_binding

    ! arguments
    real(kind=c_real), value, intent(in) :: temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt
    logical(kind=c_bool), value, intent(in) :: do_predict_nc
    real(kind=c_real), intent(inout) :: qv2qi_nucleat_tend, ni_nucleat_tend
 end subroutine ice_nucleation_f

 subroutine ice_cldliq_wet_growth_f(rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, &
                                    mu, sc, qv, qc_incld,qi_incld, ni_incld, qr_incld,     &
                                    log_wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend) bind(C)
   use iso_c_binding

   ! argmens
   real(kind=c_real), value, intent(in) :: rho, temp ,pres, rhofaci, table_val_qi2qr_melting,table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, &
                                           kap, mu, sc, qv, qc_incld, qi_incld, ni_incld,qr_incld
   logical(kind=c_bool), intent(inout) :: log_wetgrowth
   real(kind=c_real), intent(inout) :: qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend
 end subroutine ice_cldliq_wet_growth_f

 subroutine get_latent_heat_f(its, ite, kts, kte, v, s, f) bind(C)
   use iso_c_binding

   ! arguments
   integer(kind=c_int), value, intent(in) :: its, ite, kts, kte
   real(kind=c_real), dimension(its:ite, kts:kte), intent(out) :: v, s, f
 end subroutine get_latent_heat_f

 real(kind=c_real) function subgrid_variance_scaling_f(relvar,expon) bind(C)
   use iso_c_binding

   ! arguments
   real(kind=c_real), value, intent(in) :: relvar,expon
   ! return
   !real(kind=c_real) :: res

 end function subgrid_variance_scaling_f

 subroutine check_values_f(qv, temp, kts, kte, timestepcount, force_abort, source_ind, col_loc) bind(C)
   use iso_c_binding

   integer(kind=c_int), value, intent(in) :: kts, kte, timestepcount, source_ind
   logical(kind=c_bool), value, intent(in) :: force_abort

   ! arguments
   real(kind=c_real), intent(in) :: qv(kts:kte)
   real(kind=c_real), intent(in) :: temp(kts:kte)
   real(kind=c_real), intent(in) :: col_loc(3)

 end subroutine check_values_f

 subroutine calculate_incloud_mixingratios_f(qc, qr, qi, qm, nc, nr, ni, bm,  &
                                             inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, &
                                             qc_incld, qr_incld, qi_incld, qm_incld, &
                                             nc_incld, nr_incld, ni_incld, bm_incld) bind(C)
   use iso_c_binding

   ! argumens
   real(kind=c_real), value, intent(in) :: qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
   real(kind=c_real), intent(out) :: qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld
 end subroutine calculate_incloud_mixingratios_f

 subroutine p3_main_part1_f(kts, kte, kbot, ktop, kdir, do_predict_nc, dt, &
       pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, &
       t, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, &
       qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, qm_incld, &
       nc_incld, nr_incld, ni_incld, bm_incld, is_nucleat_possible, is_hydromet_present) bind(C)

   use iso_c_binding

   ! arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: do_predict_nc
   real(kind=c_real), value, intent(in) :: dt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, &
        inv_cld_frac_r, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion

   real(kind=c_real), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, &
        acn, qv, th, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, &
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld

   logical(kind=c_bool), intent(out) :: is_nucleat_possible, is_hydromet_present
 end subroutine p3_main_part1_f

 subroutine p3_main_part2_f(kts, kte, kbot, ktop, kdir, do_predict_nc, dt, inv_dt, &
       pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r,&
       t, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th, qc, nc, qr, nr, qi, ni, &
       qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, &
       ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, cmeiout, precip_total_tend, &
       nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange, pratot, &
       prctot, is_hydromet_present) bind(C)

   use iso_c_binding

   !arguments
   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   logical(kind=c_bool), value, intent(in) :: do_predict_nc
   real(kind=c_real), value, intent(in) :: dt, inv_dt

   real(kind=c_real), intent(in), dimension(kts:kte) :: pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, &
        inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r

   real(kind=c_real), intent(inout), dimension(kts:kte) :: t, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, &
        qv, th, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, &
        qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, &
        cdistr, mu_r, lamr, logn0r, cmeiout, precip_total_tend, nevapr, qr_evap_tend, vap_liq_exchange, &
        vap_ice_exchange, liq_ice_exchange, pratot, prctot

   logical(kind=c_bool), intent(out) :: is_hydromet_present

 end subroutine p3_main_part2_f

 subroutine p3_main_part3_f(kts, kte, kbot, ktop, kdir, &
      exner, cld_frac_l, cld_frac_r, &
      rho, inv_rho, rhofaci, qv, th, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, &
      ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, rho_qi, diag_ze, diag_effc) bind(C)

   use iso_c_binding

   ! args

   integer(kind=c_int), value, intent(in) :: kts, kte, kbot, ktop, kdir
   real(kind=c_real), intent(in), dimension(kts:kte) :: exner, cld_frac_l, cld_frac_r
   real(kind=c_real), intent(inout), dimension(kts:kte) :: rho, inv_rho, rhofaci, &
        qv, th, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, &
        mu_c, nu, lamc, mu_r, &
        lamr, vap_liq_exchange, &
        ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, rho_qi, diag_ze, diag_effc

 end subroutine p3_main_part3_f

 subroutine p3_main_f(qc,nc,qr,nr,th,qv,dt,qi,qm,ni,bm,   &
      pres,dz,nc_nuceat_tend,ni_activated,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_ze,diag_effc,     &
      diag_effi,diag_vmi,diag_di,rho_qi,do_predict_nc, &
      dpres,exner,cmeiout,precip_total_tend,nevapr,qr_evap_tend,precip_liq_flux,precip_ice_flux,cld_frac_r,cld_frac_l,cld_frac_i,  &
      pratot,prctot,mu_c,lamc,liq_ice_exchange,vap_liq_exchange, vap_ice_exchange) bind(C)

   use iso_c_binding

   ! args

   real(kind=c_real), intent(inout), dimension(its:ite,kts:kte) :: qc, nc, qr, nr, qi, qm, ni, bm, qv, th
   real(kind=c_real), intent(in),  dimension(its:ite,kts:kte) :: pres, dz, nc_nuceat_tend, ni_activated, dpres, exner, cld_frac_i, cld_frac_l, cld_frac_r, inv_qc_relvar
   real(kind=c_real), intent(out), dimension(its:ite,kts:kte) :: diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, rho_qi, mu_c, &
        lamc, cmeiout, precip_total_tend, nevapr, qr_evap_tend, pratot, prctot, liq_ice_exchange, vap_liq_exchange, vap_ice_exchange
   real(kind=c_real), intent(out), dimension(its:ite,kts:kte+1) :: precip_liq_flux, precip_ice_flux
   real(kind=c_real), intent(out), dimension(its:ite) :: precip_liq_surf, precip_ice_surf

   integer(kind=c_int), value, intent(in)  :: its, ite, kts, kte, it
   logical(kind=c_bool), value, intent(in) :: do_predict_nc
   real(kind=c_real), value, intent(in)    :: dt

 end subroutine p3_main_f

end interface

end module micro_p3_iso_f
