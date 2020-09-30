#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to micro_p3 fortran calls and vice versa
//

extern "C" {

void p3_init_a_c(Real* ice_table_vals, Real* collect_table_vals);

void find_lookuptable_indices_1a_c(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qi, Real ni, Real qm, Real rhop);

void find_lookuptable_indices_1b_c(Int* dumj, Real* dum3, Real qr, Real nr);

void access_lookup_table_c(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc);

void access_lookup_table_coll_c(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc);

void back_to_cell_average_c(Real cld_frac_l_, Real cld_frac_r_, Real cld_frac_i_,
                            Real* qc2qr_accret_tend_, Real* qr2qv_evap_tend_, Real* qc2qr_autoconv_tend_,
                            Real* nc_accret_tend_, Real* nc_selfcollect_tend_, Real* nc2nr_autoconv_tend_,
                            Real* nr_selfcollect_tend_, Real* nr_evap_tend_, Real* ncautr_,
                            Real* qi2qv_sublim_tend_,
                            Real* nr_ice_shed_tend_, Real* qc2qi_hetero_freeze_tend_, Real* qr2qi_collect_tend_,
                            Real* qc2qr_ice_shed_tend_, Real* qi2qr_melt_tend_, Real* qc2qi_collect_tend_,
                            Real* qr2qi_immers_freeze_tend_, Real* ni2nr_melt_tend_, Real* nc_collect_tend_,
                            Real* ncshdc_, Real* nc2ni_immers_freeze_tend_, Real* nr_collect_tend_,
                            Real* ni_selfcollect_tend_, Real* qv2qi_vapdep_tend_, Real* nr2ni_immers_freeze_tend_,
                            Real* ni_sublim_tend_, Real* qv2qi_nucleat_tend_, Real* ni_nucleat_tend_,
                            Real* qc2qi_berg_tend_);

void prevent_ice_overdepletion_c(Real pres, Real T_atm, Real qv, Real latent_heat_sublim,
                                 Real inv_dt, Real* qv2qi_vapdep_tend, Real* qi2qv_sublim_tend);

void cloud_water_conservation_c(Real qc, Real dt, Real* qc2qr_autoconv_tend, Real* qc2qr_accret_tend, Real* qc2qi_collect_tend,
  Real* qc2qi_hetero_freeze_tend, Real* qc2qr_ice_shed_tend, Real* qc2qi_berg_tend, Real* qi2qv_sublim_tend, Real* qv2qi_vapdep_tend);

void rain_water_conservation_c(Real qr, Real qc2qr_autoconv_tend, Real qc2qr_accret_tend, Real qi2qr_melt_tend, Real qc2qr_ice_shed_tend,
  Real dt, Real* qr2qv_evap_tend, Real* qr2qi_collect_tend, Real* qr2qi_immers_freeze_tend);

void ice_water_conservation_c(Real qi, Real qv2qi_vapdep_tend, Real qv2qi_nucleat_tend, Real qc2qi_berg_tend, Real qr2qi_collect_tend, Real qc2qi_collect_tend,
  Real qr2qi_immers_freeze_tend, Real qc2qi_hetero_freeze_tend, Real dt, Real* qi2qv_sublim_tend, Real* qi2qr_melt_tend);

void get_cloud_dsd2_c(Real qc, Real* nc, Real* mu_c, Real rho, Real* nu, Real* lamc,
                      Real* cdist, Real* cdist1);

void get_rain_dsd2_c(Real qr, Real* nr, Real* mu_r, Real* lamr, Real* cdistr, Real* logn0r);

void calc_rime_density_c(Real T_atm, Real rhofaci, Real table_val_qi_fallspd, Real acn,
                         Real lamc, Real mu_c, Real qc_incld, Real qc2qi_collect_tend,
                         Real* vtrmi1, Real* rho_qm_cloud);

void cldliq_immersion_freezing_c(Real T_atm, Real lamc, Real mu_c, Real cdist1,
                                 Real qc_incld, Real inv_qc_relvar, Real* qc2qi_hetero_freeze_tend, Real* nc2ni_immers_freeze_tend);

void rain_immersion_freezing_c(Real T_atm, Real lamr, Real mu_r, Real cdistr,
                               Real qr_incld, Real* qr2qi_immers_freeze_tend, Real* nr2ni_immers_freeze_tend);

void droplet_self_collection_c(Real rho, Real inv_rho, Real qc_incld, Real mu_c,
                               Real nu, Real nc2nr_autoconv_tend, Real* nc_accret_tend);

void cloud_rain_accretion_c(Real rho, Real inv_rho, Real qc_incld, Real nc_incld,
                            Real qr_incld, Real inv_qc_relvar, Real* qc2qr_accret_tend, Real* nc_accret_tend);

void cloud_water_autoconversion_c(Real rho, Real qc_incld, Real nc_incld, Real inv_qc_relvar, Real* qc2qr_autoconv_tend, Real* nc2nr_autoconv_tend, Real* ncautr);

void rain_self_collection_c(Real rho, Real qr_incld, Real nr_incld, Real* nr_selfcollect_tend);

void impose_max_total_ni_c(Real* ni_local, Real max_total_ni, Real inv_rho_local);

void ice_melting_c(Real rho,Real T_atm,Real pres,Real rhofaci,Real table_val_qi2qr_melting,Real table_val_qi2qr_vent_melt,
                   Real latent_heat_vapor,Real latent_heat_fusion,Real dv,Real sc,Real mu,Real kap,Real qv,Real qi_incld,
                   Real ni_incld,Real* qi2qr_melt_tend,Real* ni2nr_melt_tend);

void calc_first_order_upwind_step_c(Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho,
                                    Real* inv_rho, Real* inv_dz, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

void generalized_sedimentation_c(Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);
void cloud_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend);

void ice_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend);

void rain_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend);

void calc_bulk_rho_rime_c(Real qi_tot, Real* qi_rim, Real* bi_rim, Real* rho_rime);

void homogeneous_freezing_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* T_atm, Real* exner, Real* latent_heat_fusion,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm);

void get_time_space_phys_variables_c(Real T_atm, Real pres, Real rho, Real latent_heat_vapor, Real latent_heat_sublim, Real qv_sat_l, Real qv_sat_i,
  Real* mu, Real* dv, Real* sc, Real* dqsdt, Real* dqsidt, Real* ab, Real* abi, Real* kap, Real* eii);

void  update_prognostic_ice_c(
  Real qc2qi_hetero_freeze_tend, Real qc2qi_collect_tend, Real qc2qr_ice_shed_tend,  Real nc_collect_tend,  Real nc2ni_immers_freeze_tend, Real ncshdc,
  Real qr2qi_collect_tend,  Real nr_collect_tend, Real qr2qi_immers_freeze_tend, Real nr2ni_immers_freeze_tend, Real nr_ice_shed_tend,
  Real qi2qr_melt_tend, Real ni2nr_melt_tend, Real qi2qv_sublim_tend, Real qv2qi_vapdep_tend, Real qv2qi_nucleat_tend, Real ni_nucleat_tend,
  Real ni_selfcollect_tend, Real ni_sublim_tend, Real qc2qi_berg_tend, Real exner, Real latent_heat_sublim, Real latent_heat_fusion,
  bool do_predict_nc, bool log_wetgrowth, Real dt, Real nmltratio,
  Real rho_qm_cloud, Real* th_atm, Real* qv, Real* qi, Real* ni, Real* qm,
  Real* bm, Real* qc, Real* nc, Real* qr, Real* nr);

void evaporate_rain_c( Real qr_incld, Real qc_incld, Real nr_incld, Real qi_incld,
  Real cld_frac_l, Real cld_frac_r, Real qv, Real qv_prev,
  Real qv_sat_l, Real qv_sat_i, Real ab, Real abi,
  Real epsr, Real epsi_tot, Real t, Real t_prev,
  Real latent_heat_sublim, Real dqsdt, Real dt,
  Real* qr2qv_evap_tend, Real* nr_evap_tend);

void update_prognostic_liquid_c(
  Real qc2qr_accret_tend, Real nc_accret_tend, Real qc2qr_autoconv_tend, Real nc2nr_autoconv_tend, Real ncautr,
  Real nc_selfcollect_tend, Real  qr2qv_evap_tend, Real nr_evap_tend, Real nr_selfcollect_tend , bool do_predict_nc,
  Real inv_rho, Real exner, Real latent_heat_vapor, Real dt, Real* th_atm, Real* qv,
  Real* qc, Real* nc, Real* qr, Real* nr);

void ice_deposition_sublimation_c(
  Real qi_incld, Real ni_incld, Real T_atm, Real qv_sat_l, Real qv_sat_i, Real epsi,
  Real abi, Real qv, Real* qv2qi_vapdep_tend, Real* qi2qv_sublim_tend, Real* ni_sublim_tend, Real* qc2qi_berg_tend);

void compute_rain_fall_velocity_c(Real qr_incld, Real rhofacr,
                                  Real* nr_incld, Real* mu_r, Real* lamr, Real* V_qr, Real* V_nr);

void ice_cldliq_collection_c(Real rho, Real temp, Real rhofaci, Real table_val_qc2qi_collect,
                             Real qi_incld,Real qc_incld, Real ni_incld, Real nc_incld,
                             Real* qc2qi_collect_tend, Real* nc_collect_tend, Real* qc2qr_ice_shed_tend, Real* ncshdc);

void ice_rain_collection_c(Real rho, Real temp, Real rhofaci, Real logn0r, Real table_val_nr_collect, Real table_val_qr2qi_collect,
                           Real qi_incld, Real ni_incld, Real qr_incld, Real* qr2qi_collect_tend, Real* nr_collect_tend);


void ice_self_collection_c(Real rho, Real rhofaci, Real table_val_ni_self_collect, Real eii,
                           Real qm_incld, Real qi_incld, Real ni_incld, Real* ni_selfcollect_tend);

void ice_relaxation_timescale_c(Real rho, Real temp, Real rhofaci, Real table_val_qi2qr_melting, Real table_val_qi2qr_vent_melt,
                                Real dv, Real mu, Real sc, Real qi_incld, Real ni_incld,
                                Real* epsi, Real* epsi_tot);

void calc_liq_relaxation_timescale_c(Real rho, Real f1r, Real f2r, Real dv,
                                     Real mu, Real sc, Real mu_r, Real lamr,
                                     Real cdistr, Real cdist, Real qr_incld,
                                     Real qc_incld, Real* epsr, Real* epsc);

void ice_nucleation_c(Real temp, Real inv_rho, Real ni, Real ni_activated,
                      Real qv_supersat_i, Real inv_dt, bool do_predict_nc,
                      Real* qv2qi_nucleat_tend, Real* ni_nucleat_tend);

void ice_cldliq_wet_growth_c(Real rho, Real temp, Real pres, Real rhofaci, Real table_val_qi2qr_melting,
                             Real table_val_qi2qr_vent_melt, Real latent_heat_vapor, Real latent_heat_fusion, Real dv,
                             Real kap, Real mu, Real sc, Real qv, Real qc_incld,
                             Real qi_incld, Real ni_incld, Real qr_incld, bool* log_wetgrowth,
                             Real* qr2qi_collect_tend, Real* qc2qi_collect_tend, Real* qc_growth_rate, Real* nr_ice_shed_tend, Real* qc2qr_ice_shed_tend);

void get_latent_heat_c(Int its, Int ite, Int kts, Int kte, Real* s, Real* v, Real* f);

Real subgrid_variance_scaling_c(Real relvar, Real expon);

void check_values_c(Real* qv, Real* temp, Int kts, Int kte, Int timestepcount,
                    Int force_abort, Int source_ind, Real* col_loc);

void calculate_incloud_mixingratios_c(Real qc, Real qr, Real qi, Real qm, Real nc, Real nr, Real ni, Real bm,
                                      Real inv_cld_frac_l, Real inv_cld_frac_i, Real inv_cld_frac_r,
                                      Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld,
                                      Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld);

void p3_main_part1_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present);

void p3_main_part2_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, Real dt, Real inv_dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci, Real* acn,
  Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion, Real* qc_incld,
  Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1,
  Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* is_hydromet_present);

void p3_main_part3_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr,
  Real* qi, Real* ni, Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi, Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc);

void p3_main_c(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* rho_qi, bool do_predict_nc, Real* dpres, Real* exner,
  Real* qv2qi_depos_tend, Real* precip_total_tend, Real* nevapr, Real* qr_evap_tend, Real* precip_liq_flux,
  Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i, Real* mu_c, Real* lamc,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* qv_prev, Real* t_prev);

} // extern "C" : end _c decls

namespace scream {
namespace p3 {

//
// In all C++ -> Fortran bridge functions you should see p3_init(). P3 needs
// to be initialized since most of its function depend on global tables to be
// populated. The 'true' argument is to set p3 to use its fortran implementations
// instead of calling back to C++. We want this behavior since it doesn't make much
// sense for C++ to bridge over to fortran only to have fortran bridge back to C++.
// If the client wanted the C++ implementation, they should just call it directly.
//

void p3_init_a(P3InitAFortranData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  p3_init_a_c(d.ice_table_vals.data(), d.collect_table_vals.data());
}

void find_lookuptable_indices_1a(LookupIceData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  find_lookuptable_indices_1a_c(&d.dumi, &d.dumjj, &d.dumii, &d.dumzz,
                                &d.dum1, &d.dum4, &d.dum5, &d.dum6,
                                d.qi, d.ni, d.qm, d.rhop);
}

void find_lookuptable_indices_1b(LookupIceDataB& d)
{
  p3_init();
  find_lookuptable_indices_1b_c(&d.dumj, &d.dum3, d.qr, d.nr);
}

void access_lookup_table(AccessLookupTableData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  access_lookup_table_c(d.lid.dumjj, d.lid.dumii, d.lid.dumi, d.index,
                        d.lid.dum1, d.lid.dum4, d.lid.dum5, &d.proc);
}

void access_lookup_table_coll(AccessLookupTableCollData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  access_lookup_table_coll_c(d.lid.dumjj, d.lid.dumii, d.lidb.dumj, d.lid.dumi, d.index,
                             d.lid.dum1, d.lidb.dum3, d.lid.dum4, d.lid.dum5, &d.proc);
}

void BackToCellAverageData::randomize()
{
  // Populate the struct with numbers between 0 and 1.
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  cld_frac_l = data_dist(generator);
  cld_frac_r = data_dist(generator);
  cld_frac_i = data_dist(generator);
  qc2qr_accret_tend = data_dist(generator);
  qr2qv_evap_tend = data_dist(generator);
  qc2qr_autoconv_tend = data_dist(generator);
  nc_accret_tend = data_dist(generator);
  nc_selfcollect_tend = data_dist(generator);
  nc2nr_autoconv_tend = data_dist(generator);
  nr_selfcollect_tend = data_dist(generator);
  nr_evap_tend = data_dist(generator);
  ncautr = data_dist(generator);
  qcnuc = data_dist(generator);
  nc_nuceat_tend = data_dist(generator);
  qi2qv_sublim_tend = data_dist(generator);
  nr_ice_shed_tend = data_dist(generator);
  qc2qi_hetero_freeze_tend = data_dist(generator);
  qr2qi_collect_tend = data_dist(generator);
  qc2qr_ice_shed_tend = data_dist(generator);
  qi2qr_melt_tend = data_dist(generator);
  qc2qi_collect_tend = data_dist(generator);
  qr2qi_immers_freeze_tend = data_dist(generator);
  ni2nr_melt_tend = data_dist(generator);
  nc_collect_tend = data_dist(generator);
  ncshdc = data_dist(generator);
  nc2ni_immers_freeze_tend = data_dist(generator);
  nr_collect_tend = data_dist(generator);
  ni_selfcollect_tend = data_dist(generator);
  qv2qi_vapdep_tend = data_dist(generator);
  nr2ni_immers_freeze_tend = data_dist(generator);
  ni_sublim_tend = data_dist(generator);
  qv2qi_nucleat_tend = data_dist(generator);
  ni_nucleat_tend = data_dist(generator);
  qc2qi_berg_tend = data_dist(generator);
}

void back_to_cell_average(BackToCellAverageData& d)
{
  p3_init();
  back_to_cell_average_c(d.cld_frac_l, d.cld_frac_r, d.cld_frac_i, &d.qc2qr_accret_tend, &d.qr2qv_evap_tend,
    &d.qc2qr_autoconv_tend, &d.nc_accret_tend, &d.nc_selfcollect_tend, &d.nc2nr_autoconv_tend, &d.nr_selfcollect_tend, &d.nr_evap_tend, &d.ncautr,
    &d.qi2qv_sublim_tend, &d.nr_ice_shed_tend, &d.qc2qi_hetero_freeze_tend, &d.qr2qi_collect_tend, &d.qc2qr_ice_shed_tend,
    &d.qi2qr_melt_tend, &d.qc2qi_collect_tend, &d.qr2qi_immers_freeze_tend, &d.ni2nr_melt_tend, &d.nc_collect_tend, &d.ncshdc, &d.nc2ni_immers_freeze_tend,
    &d.nr_collect_tend, &d.ni_selfcollect_tend, &d.qv2qi_vapdep_tend, &d.nr2ni_immers_freeze_tend, &d.ni_sublim_tend, &d.qv2qi_nucleat_tend, &d.ni_nucleat_tend,
    &d.qc2qi_berg_tend);
}

void prevent_ice_overdepletion(PreventIceOverdepletionData& d)
{
  p3_init();
  prevent_ice_overdepletion_c(d.pres, d.T_atm, d.qv, d.latent_heat_sublim, d.inv_dt, &d.qv2qi_vapdep_tend,
                              &d.qi2qv_sublim_tend);
}

void calc_rime_density(CalcRimeDensityData& d)
{
  p3_init();
  calc_rime_density_c(d.T_atm, d.rhofaci, d.table_val_qi_fallspd, d.acn, d.lamc, d.mu_c,
                      d.qc_incld, d.qc2qi_collect_tend, &d.vtrmi1, &d.rho_qm_cloud);
}

void cldliq_immersion_freezing(CldliqImmersionFreezingData& d)
{
  p3_init();
  cldliq_immersion_freezing_c(d.T_atm, d.lamc, d.mu_c, d.cdist1, d.qc_incld, d.inv_qc_relvar,
                              &d.qc2qi_hetero_freeze_tend, &d.nc2ni_immers_freeze_tend);
}

LatentHeatData::LatentHeatData(Int kts_, Int kte_, Int its_, Int ite_) :
  PhysicsTestData((ite_ - its_) + 1, (kte_ - kts_) + 1, {&v, &s, &f}),
  its(its_), ite(ite_), kts(kts_), kte(kte_)
{}

void get_latent_heat(LatentHeatData& d)
{
  p3_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  get_latent_heat_c(d.its, d.ite, d.kts, d.kte, d.v, d.s, d.f);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void droplet_self_collection(DropletSelfCollectionData& d)
{
  p3_init();
  droplet_self_collection_c(d.rho, d.inv_rho, d.qc_incld, d.mu_c, d.nu, d.nc2nr_autoconv_tend,
                            &d.nc_selfcollect_tend);
}

void rain_immersion_freezing(RainImmersionFreezingData& d)
{
  p3_init();
  rain_immersion_freezing_c(d.T_atm, d.lamr, d.mu_r, d.cdistr, d.qr_incld,
                            &d.qr2qi_immers_freeze_tend, &d.nr2ni_immers_freeze_tend);
}

void cloud_rain_accretion(CloudRainAccretionData& d)
{
  p3_init();
  cloud_rain_accretion_c(d.rho, d.inv_rho, d.qc_incld, d.nc_incld, d.qr_incld, d.inv_qc_relvar,
                         &d.qc2qr_accret_tend, &d.nc_accret_tend);
}

void cloud_water_conservation(CloudWaterConservationData& d){
  p3_init();
  cloud_water_conservation_c(d.qc, d.dt, &d.qc2qr_autoconv_tend, &d.qc2qr_accret_tend, &d.qc2qi_collect_tend, &d.qc2qi_hetero_freeze_tend,
  &d.qc2qr_ice_shed_tend, &d.qc2qi_berg_tend, &d.qi2qv_sublim_tend, &d.qv2qi_vapdep_tend);
}

void rain_water_conservation(RainWaterConservationData& d){
  p3_init();
  rain_water_conservation_c(d.qr, d.qc2qr_autoconv_tend, d.qc2qr_accret_tend, d.qi2qr_melt_tend, d.qc2qr_ice_shed_tend,
                            d.dt, &d.qr2qv_evap_tend, &d.qr2qi_collect_tend, &d.qr2qi_immers_freeze_tend);
}

void ice_water_conservation(IceWaterConservationData& d){
  p3_init();
  ice_water_conservation_c(d.qi, d.qv2qi_vapdep_tend, d.qv2qi_nucleat_tend, d.qc2qi_berg_tend, d.qr2qi_collect_tend, d.qc2qi_collect_tend, d.qr2qi_immers_freeze_tend,
    d.qc2qi_hetero_freeze_tend, d.dt, &d.qi2qv_sublim_tend, &d.qi2qr_melt_tend);
}

void cloud_water_autoconversion(CloudWaterAutoconversionData& d){
  p3_init();
  cloud_water_autoconversion_c(d.rho, d.qc_incld, d.nc_incld, d.inv_qc_relvar,
    &d.qc2qr_autoconv_tend, &d.nc2nr_autoconv_tend, &d.ncautr);
}

void rain_self_collection(RainSelfCollectionData& d){
  p3_init();
  rain_self_collection_c(d.rho, d.qr_incld, d.nr_incld, &d.nr_selfcollect_tend);
}

void impose_max_total_ni(ImposeMaxTotalNiData& d){
  p3_init();
  impose_max_total_ni_c(&d.ni_local, d.max_total_ni, d.inv_rho_local);
}

void get_cloud_dsd2(GetCloudDsd2Data& d)
{
  p3_init();
  Real nc_in = d.nc_in;
  get_cloud_dsd2_c(d.qc, &nc_in, &d.mu_c, d.rho, &d.nu, &d.lamc, &d.cdist, &d.cdist1);
  d.nc_out = nc_in;
}

void get_rain_dsd2(GetRainDsd2Data& d)
{
  p3_init();
  Real nr_in = d.nr_in;
  get_rain_dsd2_c(d.qr, &nr_in, &d.mu_r, &d.lamr, &d.cdistr, &d.logn0r);
  d.nr_out = nr_in;
}

void ice_cldliq_collection(IceCldliqCollectionData& d)
{
  p3_init();
  ice_cldliq_collection_c(d.rho, d.temp, d.rhofaci, d.table_val_qc2qi_collect,
                          d.qi_incld, d.qc_incld, d.ni_incld, d.nc_incld,
                          &d.qc2qi_collect_tend, &d.nc_collect_tend, &d.qc2qr_ice_shed_tend, &d.ncshdc);
}

void ice_rain_collection(IceRainCollectionData& d)
{
  p3_init();
  ice_rain_collection_c(d.rho, d.temp, d.rhofaci, d.logn0r, d.table_val_nr_collect, d.table_val_qr2qi_collect,
                        d.qi_incld, d.ni_incld, d.qr_incld,
                        &d.qr2qi_collect_tend, &d.nr_collect_tend);
}

void ice_self_collection(IceSelfCollectionData& d)
{
  p3_init();
  ice_self_collection_c(d.rho, d.rhofaci, d.table_val_ni_self_collect, d.eii, d.qm_incld,
                        d.qi_incld, d.ni_incld,
                        &d.ni_selfcollect_tend);
}

void get_time_space_phys_variables(GetTimeSpacePhysVarsData& d)
{
  p3_init();
  get_time_space_phys_variables_c(d.T_atm, d.pres, d.rho, d.latent_heat_vapor, d.latent_heat_sublim, d.qv_sat_l, d.qv_sat_i, &d.mu, &d.dv,
				  &d.sc, &d.dqsdt, &d.dqsidt, &d.ab, &d.abi, &d.kap, &d.eii);
}

void ice_relaxation_timescale(IceRelaxationData& d)
{
  p3_init();
  ice_relaxation_timescale_c(d.rho, d.temp, d.rhofaci, d.table_val_qi2qr_melting, d.table_val_qi2qr_vent_melt,
                             d.dv, d.mu, d.sc, d.qi_incld, d.ni_incld,
                             &d.epsi, &d.epsi_tot);
}

void CalcLiqRelaxationData::randomize()
{
  // Populate the struct's input fields with numbers between 0 and 1.
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  rho = data_dist(generator);
  f1r = data_dist(generator);
  f2r = data_dist(generator);
  dv = data_dist(generator);
  mu = data_dist(generator);
  sc = data_dist(generator);
  mu_r = data_dist(generator);
  lamr = data_dist(generator);
  cdistr = data_dist(generator);
  cdist = data_dist(generator);
  qr_incld = data_dist(generator);
  qc_incld = data_dist(generator);
}

void calc_liq_relaxation_timescale(CalcLiqRelaxationData& d)
{
  p3_init();
  calc_liq_relaxation_timescale_c(d.rho, d.f1r, d.f2r, d.dv, d.mu, d.sc, d.mu_r,
    d.lamr, d.cdistr, d.cdist, d.qr_incld, d.qc_incld, &d.epsr, &d.epsc);
}

void ice_nucleation(IceNucleationData& d)
{
  p3_init();
  ice_nucleation_c(d.temp, d.inv_rho, d.ni, d.ni_activated,
                   d.qv_supersat_i, d.inv_dt, d.do_predict_nc,&d.qv2qi_nucleat_tend, &d.ni_nucleat_tend);
}

void ice_cldliq_wet_growth(IceWetGrowthData& d)
{
  p3_init();

  ice_cldliq_wet_growth_c(d.rho, d.temp, d.pres, d.rhofaci, d.table_val_qi2qr_melting,
                          d.table_val_qi2qr_vent_melt, d.latent_heat_vapor, d.latent_heat_fusion, d.dv,
                          d.kap, d.mu, d.sc, d.qv, d.qc_incld,
                          d.qi_incld, d.ni_incld, d.qr_incld, &d.log_wetgrowth,
                          &d.qr2qi_collect_tend, &d.qc2qi_collect_tend, &d.qc_growth_rate, &d.nr_ice_shed_tend, &d.qc2qr_ice_shed_tend);
}

CheckValuesData::CheckValuesData(
  Int kts_, Int kte_, Int timestepcount_, Int source_ind_, bool force_abort_) :
  PhysicsTestData((kte_-kts_)+1, {&qv, &temp, &col_loc}),
  kts(kts_), kte(kte_), timestepcount(timestepcount_), source_ind(source_ind_), force_abort(force_abort_)
{
  EKAT_REQUIRE_MSG(nk() >= 3 || (kte == 0 && kts == 0), "nk too small to use for col_loc");
}

void check_values(CheckValuesData& d)
{
  p3_init();
  check_values_c(d.qv, d.temp, d.kts, d.kte, d.timestepcount,
                 d.force_abort, d.source_ind, d.col_loc);
}

void calculate_incloud_mixingratios(IncloudMixingData& d)
{
  p3_init();

  calculate_incloud_mixingratios_c(d.qc, d.qr, d.qi, d.qm, d.nc, d.nr, d.ni, d.bm, d.inv_cld_frac_l, d.inv_cld_frac_i, d.inv_cld_frac_r,
                                   &d.qc_incld, &d.qr_incld, &d.qi_incld, &d.qm_incld,
                                   &d.nc_incld, &d.nr_incld, &d.ni_incld, &d.bm_incld);

}

void update_prognostic_ice(P3UpdatePrognosticIceData& d){
  p3_init();
  update_prognostic_ice_c(d.qc2qi_hetero_freeze_tend, d.qc2qi_collect_tend, d.qc2qr_ice_shed_tend,  d.nc_collect_tend,  d.nc2ni_immers_freeze_tend, d.ncshdc,
                          d.qr2qi_collect_tend,  d.nr_collect_tend, d.qr2qi_immers_freeze_tend, d.nr2ni_immers_freeze_tend, d.nr_ice_shed_tend,
                          d.qi2qr_melt_tend,  d.ni2nr_melt_tend, d.qi2qv_sublim_tend,  d.qv2qi_vapdep_tend,  d.qv2qi_nucleat_tend,  d.ni_nucleat_tend,
                          d.ni_selfcollect_tend,  d.ni_sublim_tend, d.qc2qi_berg_tend, d.exner,  d.latent_heat_sublim,   d.latent_heat_fusion,
                          d.do_predict_nc,  d.log_wetgrowth,    d.dt,     d.nmltratio,
                          d.rho_qm_cloud,      &d.th_atm,    &d.qv,    &d.qi, &d.ni, &d.qm,
                          &d.bm,         &d.qc,    &d.nc,    &d.qr, &d.nr);
}

void evaporate_rain(EvapRainData& d)
{
  p3_init();
  evaporate_rain_c(d.qr_incld,d.qc_incld,d.nr_incld,d.qi_incld,
		   d.cld_frac_l,d.cld_frac_r,d.qv,d.qv_prev,d.qv_sat_l,d.qv_sat_i,
		   d.ab,d.abi,d.epsr,d.epsi_tot,d.t,d.t_prev,d.latent_heat_sublim,d.dqsdt,d.dt,
		   &d.qr2qv_evap_tend,&d.nr_evap_tend);
}

void  update_prognostic_liquid(P3UpdatePrognosticLiqData& d){
  p3_init();
  update_prognostic_liquid_c(d.qc2qr_accret_tend, d.nc_accret_tend, d.qc2qr_autoconv_tend, d.nc2nr_autoconv_tend, d.ncautr,
			      d.nc_selfcollect_tend, d. qr2qv_evap_tend, d.nr_evap_tend, d.nr_selfcollect_tend , d.do_predict_nc,
			      d.inv_rho, d.exner, d.latent_heat_vapor, d.dt, &d.th_atm, &d.qv,
			      &d.qc, &d.nc, &d.qr, &d.nr);
}

void ice_deposition_sublimation(IceDepSublimationData& d){
  p3_init();
  ice_deposition_sublimation_c(d.qi_incld, d.ni_incld, d.T_atm, d.qv_sat_l, d.qv_sat_i, d.epsi, d.abi,
			       d.qv, &d.qv2qi_vapdep_tend, &d.qi2qv_sublim_tend, &d.ni_sublim_tend, &d.qc2qi_berg_tend);
}

CalcUpwindData::CalcUpwindData(
  Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_) :
  PhysicsTestData((kte_ - kts_) + 1, num_arrays_, {&vs, &qnx, &fluxes}, {&rho, &inv_rho, &inv_dz}),
  kts(kts_), kte(kte_), kdir(kdir_), kbot(kbot_), k_qxtop(k_qxtop_), num_arrays(num_arrays_), dt_sub(dt_sub_)
{}

void CalcUpwindData::convert_to_ptr_arr(std::vector<Real*>& mem_space, Real**& fluxes_, Real**& vs_, Real**& qnx_)
{
  mem_space.resize(num_arrays*3);
  for (Int i = 0; i < num_arrays; ++i) {
    mem_space[i]              = fluxes + (i*nk());
    mem_space[i+num_arrays]   = vs     + (i*nk());
    mem_space[i+num_arrays*2] = qnx    + (i*nk());
  }
  fluxes_ = mem_space.data();
  vs_     = mem_space.data() + num_arrays;
  qnx_    = mem_space.data() + num_arrays*2;
}

void calc_first_order_upwind_step(CalcUpwindData& d)
{
  p3_init();
  std::vector<Real*> tmp;
  Real** fluxes, **vs, **qnx;
  d.convert_to_ptr_arr(tmp, fluxes, vs, qnx);
  calc_first_order_upwind_step_c(d.kts, d.kte, d.kdir, d.kbot, d.k_qxtop, d.dt_sub, d.rho, d.inv_rho, d.inv_dz, d.num_arrays, fluxes, vs, qnx);
}

GenSedData::GenSedData(
  Int kts_, Int kte_, Int kdir_, Int k_qxtop_, Int k_qxbot_, Int kbot_, Real Co_max_, Real dt_left_,
  Real prt_accum_, Int num_arrays_) :
  CalcUpwindData(kts_, kte_, kdir_, kbot_, k_qxtop_, num_arrays_, 0.0),
  Co_max(Co_max_), k_qxbot(k_qxbot_), dt_left(dt_left_), prt_accum(prt_accum_)
{ }

void generalized_sedimentation(GenSedData& d)
{
  p3_init();
  std::vector<Real*> tmp;
  Real** fluxes, **vs, **qnx;
  d.convert_to_ptr_arr(tmp, fluxes, vs, qnx);
  generalized_sedimentation_c(d.kts, d.kte, d.kdir, d.k_qxtop, &d.k_qxbot, d.kbot, d.Co_max,
                              &d.dt_left, &d.prt_accum, d.inv_dz, d.inv_rho, d.rho,
                              d.num_arrays, fluxes, vs, qnx);
}

CloudSedData::CloudSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, bool do_predict_nc_, Real precip_liq_surf_) :
  PhysicsTestData((kte_ - kts_) + 1, {&qc_incld, &rho, &inv_rho, &cld_frac_l, &acn, &inv_dz, &qc, &nc, &nc_incld, &mu_c, &lamc, &qc_tend, &nc_tend}),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), do_predict_nc(do_predict_nc_), precip_liq_surf(precip_liq_surf_)
{}

void cloud_sedimentation(CloudSedData& d)
{
  p3_init();
  cloud_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                        d.qc_incld, d.rho, d.inv_rho, d.cld_frac_l, d.acn, d.inv_dz,
                        d.dt, d.inv_dt, d.do_predict_nc,
                        d.qc, d.nc, d.nc_incld, d.mu_c, d.lamc, &d.precip_liq_surf, d.qc_tend, d.nc_tend);
}

IceSedData::IceSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, Real precip_ice_surf_) :
  PhysicsTestData((kte_ - kts_) + 1, {&rho, &inv_rho, &rhofaci, &cld_frac_i, &inv_dz, &qi, &qi_incld, &ni, &ni_incld, &qm, &qm_incld, &bm, &bm_incld, &qi_tend, &ni_tend}),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), precip_ice_surf(precip_ice_surf_)
{}

void ice_sedimentation(IceSedData& d)
{
  p3_init();
  ice_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                      d.rho, d.inv_rho, d.rhofaci, d.cld_frac_i, d.inv_dz, d.dt, d.inv_dt,
                      d.qi, d.qi_incld, d.ni, d.qm, d.qm_incld, d.bm, d.bm_incld, d.ni_incld,
                      &d.precip_ice_surf, d.qi_tend, d.ni_tend);
}

RainSedData::RainSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, Real precip_liq_surf_) :
  PhysicsTestData((kte_ - kts_) + 2, // extra real at end for precip_liq_flux, so just add 1 to all
                  {&rho, &inv_rho, &rhofacr, &cld_frac_r, &inv_dz, &qr_incld, &qr, &nr, &nr_incld, &mu_r, &lamr, &qr_tend, &nr_tend, &precip_liq_flux}),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), precip_liq_surf(precip_liq_surf_)
{}

void rain_sedimentation(RainSedData& d)
{
  p3_init();
  rain_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                       d.qr_incld, d.rho, d.inv_rho, d.rhofacr, d.cld_frac_r, d.inv_dz,
                       d.dt, d.inv_dt,
                       d.qr, d.nr, d.nr_incld, d.mu_r, d.lamr, &d.precip_liq_surf, d.precip_liq_flux, d.qr_tend, d.nr_tend);
}

void calc_bulk_rho_rime(CalcBulkRhoRimeData& d)
{
  p3_init();
  calc_bulk_rho_rime_c(d.qi_tot, &d.qi_rim, &d.bi_rim, &d.rho_rime);
}

HomogeneousFreezingData::HomogeneousFreezingData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_) :
  PhysicsTestData((kte_ - kts_) + 1, {&T_atm, &exner, &latent_heat_fusion, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &th_atm}),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_)
{}

void homogeneous_freezing(HomogeneousFreezingData& d)
{
  p3_init();
  homogeneous_freezing_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.T_atm, d.exner, d.latent_heat_fusion,
                         d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.th_atm);
}

void ice_melting(IceMeltingData& d){
  p3_init();
  ice_melting_c(d.rho,d.T_atm,d.pres,d.rhofaci,d.table_val_qi2qr_melting,d.table_val_qi2qr_vent_melt,
		d.latent_heat_vapor,d.latent_heat_fusion,d.dv,d.sc,d.mu,d.kap,
		d.qv,d.qi_incld,d.ni_incld,&d.qi2qr_melt_tend,&d.ni2nr_melt_tend);
}

Real subgrid_variance_scaling(SubgridVarianceScalingData& d){
  p3_init();
  return subgrid_variance_scaling_c(d.relvar,d.expon);
}

void compute_rain_fall_velocity(ComputeRainFallVelocityData& d)
{
  p3_init();
  compute_rain_fall_velocity_c(d.qr_incld, d.rhofacr,
                               &d.nr_incld, &d.mu_r, &d.lamr, &d.V_qr, &d.V_nr);
}

P3MainPart1Data::P3MainPart1Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool do_predict_nc_, Real dt_) :
  PhysicsTestData((kte_ - kts_) + 1, {
    &pres, &dpres, &dz, &nc_nuceat_tend, &exner, &inv_exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci,
    &acn, &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &qc_incld, &qr_incld, &qi_incld,
    &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld}),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), dt(dt_)
{}

void p3_main_part1(P3MainPart1Data& d)
{
  p3_init();
  p3_main_part1_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.do_predict_nc,
    d.dt,
    d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.exner, d.inv_exner, d.inv_cld_frac_l, d.inv_cld_frac_i, d.inv_cld_frac_r, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion,
    d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci,
    d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.qc_incld, d.qr_incld, d.qi_incld,
    d.qm_incld, d.nc_incld, d.nr_incld, d.ni_incld, d.bm_incld,
    &d.is_nucleat_possible, &d.is_hydromet_present);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart2Data::P3MainPart2Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool do_predict_nc_, Real dt_) :
  PhysicsTestData((kte_ - kts_) + 1, {
    &pres, &dpres, &dz, &nc_nuceat_tend, &exner, &inv_exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &ni_activated, &inv_qc_relvar, &cld_frac_i, &cld_frac_l, &cld_frac_r, &qv_prev, &t_prev,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci, &acn,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion, &qc_incld, &qr_incld,
    &qi_incld, &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld, &mu_c, &nu, &lamc, &cdist, &cdist1,
    &cdistr, &mu_r, &lamr, &logn0r, &qv2qi_depos_tend, &precip_total_tend, &nevapr, &qr_evap_tend, &vap_liq_exchange,
    &vap_ice_exchange, &liq_ice_exchange, &pratot, &prctot}),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), dt(dt_), inv_dt(1 / dt)
{}

void p3_main_part2(P3MainPart2Data& d)
{
  p3_init();
  p3_main_part2_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir, d.do_predict_nc, d.dt, d.inv_dt,
    d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.exner, d.inv_exner, d.inv_cld_frac_l, d.inv_cld_frac_i, d.inv_cld_frac_r, d.ni_activated, d.inv_qc_relvar, d.cld_frac_i, d.cld_frac_l, d.cld_frac_r, d.qv_prev, d.t_prev,
    d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci, d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni,
    d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion, d.qc_incld, d.qr_incld, d.qi_incld, d.qm_incld, d.nc_incld, d.nr_incld,
    d.ni_incld, d.bm_incld, d.mu_c, d.nu, d.lamc, d.cdist, d.cdist1, d.cdistr, d.mu_r, d.lamr, d.logn0r, d.qv2qi_depos_tend, d.precip_total_tend,
    d.nevapr, d.qr_evap_tend, d.vap_liq_exchange, d.vap_ice_exchange, d.liq_ice_exchange, d.pratot,
    d.prctot, &d.is_hydromet_present);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart3Data::P3MainPart3Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_) :
  PhysicsTestData((kte_ - kts_) + 1, {
    &exner, &cld_frac_l, &cld_frac_r, &cld_frac_i,
    &rho, &inv_rho, &rhofaci,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim,
    &mu_c, &nu, &lamc, &mu_r,
    &lamr, &vap_liq_exchange,
    &ze_rain, &ze_ice, &diag_vm_qi, &diag_eff_radius_qi, &diag_diam_qi, &rho_qi, &diag_equiv_reflectivity, &diag_eff_radius_qc}),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_)
{}

void p3_main_part3(P3MainPart3Data& d)
{
  p3_init();
  p3_main_part3_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.exner, d.cld_frac_l, d.cld_frac_r, d.cld_frac_i, 
    d.rho, d.inv_rho, d.rhofaci, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim,
    d.mu_c, d.nu, d.lamc, d.mu_r, d.lamr, d.vap_liq_exchange,
    d. ze_rain, d.ze_ice, d.diag_vm_qi, d.diag_eff_radius_qi, d.diag_diam_qi, d.rho_qi, d.diag_equiv_reflectivity, d.diag_eff_radius_qc);
}

///////////////////////////////////////////////////////////////////////////////

P3MainData::P3MainData(
  Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_) :
  PhysicsTestData( (ite_ - its_) + 1, (kte_ - kts_) + 1, (kte_ - kts_) + 2, {
    &pres, &dz, &nc_nuceat_tend, &ni_activated, &dpres, &exner, &cld_frac_i, &cld_frac_l, &cld_frac_r,
    &inv_qc_relvar, &qc, &nc, &qr, &nr, &qi, &qm, &ni, &bm, &qv, &th_atm, &qv_prev, &t_prev,
    &diag_eff_radius_qc, &diag_eff_radius_qi, &rho_qi, &mu_c, &lamc, &qv2qi_depos_tend, &precip_total_tend, &nevapr,
    &qr_evap_tend, &liq_ice_exchange, &vap_liq_exchange, &vap_ice_exchange, &precip_liq_flux,
    &precip_ice_flux},
    {&precip_liq_surf, &precip_ice_surf}), // these two are (ni, nk+1)
  its(its_), ite(ite_), kts(kts_), kte(kte_), it(it_), dt(dt_), do_predict_nc(do_predict_nc_)
{}

//This is the variable ordering from micro_p3.F90
void p3_main(P3MainData& d)
{
  p3_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  p3_main_c(
    d.qc, d.nc, d.qr, d.nr, d.th_atm, d.qv, d.dt, d.qi, d.qm, d.ni,
    d.bm, d.pres, d.dz, d.nc_nuceat_tend, d.ni_activated, d.inv_qc_relvar, d.it, d.precip_liq_surf,
    d.precip_ice_surf, d.its, d.ite, d.kts, d.kte, d.diag_eff_radius_qc, d.diag_eff_radius_qi,
    d.rho_qi, d.do_predict_nc, d.dpres, d.exner, d.qv2qi_depos_tend, d.precip_total_tend, d.nevapr,
    d.qr_evap_tend, d.precip_liq_flux, d.precip_ice_flux, d.cld_frac_r, d.cld_frac_l, d.cld_frac_i, d.mu_c, d.lamc,
    d.liq_ice_exchange, d.vap_liq_exchange, d.vap_ice_exchange, d.qv_prev, d.t_prev);
  d.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

///////////////////////////////////////////////////////////////////////////////

std::shared_ptr<P3GlobalForFortran::Views> P3GlobalForFortran::s_views;

const P3GlobalForFortran::Views& P3GlobalForFortran::get()
{
  if (!P3GlobalForFortran::s_views) {
    P3GlobalForFortran::s_views = std::make_shared<Views>();
    P3F::init_kokkos_ice_lookup_tables(s_views->m_ice_table_vals, s_views->m_collect_table_vals);
    P3F::init_kokkos_tables(s_views->m_vn_table_vals, s_views->m_vm_table_vals,
      s_views->m_revap_table_vals, s_views->m_mu_r_table_vals, s_views->m_dnu);
  }
  return *P3GlobalForFortran::s_views;
}

void P3GlobalForFortran::deinit()
{
  P3GlobalForFortran::s_views = nullptr;
}

//
// _f function definitions
//

void find_lookuptable_indices_1a_f(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qi_, Real ni_, Real qm_, Real rhop_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using TableIce = typename P3F::TableIce;

  typename P3F::Spack qi(qi_), ni(ni_), qm(qm_), rhop(rhop_);
  typename P3F::view_1d<TableIce> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    P3F::lookup_ice(qi, ni, qm, rhop, t_d(0));
  });
  Kokkos::deep_copy(t_h, t_d);
  auto& tab = t_h(0);

  // adjust for 1-based indexing
  *dumi  = tab.dumi[0]  + 1;
  *dumjj = tab.dumjj[0] + 1;
  *dumii = tab.dumii[0] + 1;
  *dumzz = tab.dumzz[0] + 1;

  *dum1 = tab.dum1[0];
  *dum4 = tab.dum4[0];
  *dum5 = tab.dum5[0];
  *dum6 = tab.dum6[0];
}

void find_lookuptable_indices_1b_f(Int* dumj, Real* dum3, Real qr_, Real nr_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using TableRain = typename P3F::TableRain;

  typename P3F::Spack qr(qr_), nr(nr_);
  typename P3F::view_1d<TableRain> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    P3F::lookup_rain(qr, nr, t_d(0));
  });
  Kokkos::deep_copy(t_h, t_d);
  auto& tab = t_h(0);

  // adjust for 1-based indexing
  *dumj = tab.dumj[0] + 1;

  *dum3 = tab.dum3[0];
}

void access_lookup_table_f(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::TableIce tab;

  // Adjust for 0-based indexing
  tab.dumi  = dumi  - 1;
  tab.dumjj = dumjj - 1;
  tab.dumii = dumii - 1;

  int adjusted_index = index - 1;

  tab.dum1 = dum1;
  tab.dum4 = dum4;
  tab.dum5 = dum5;

  auto ice_table_vals = P3GlobalForFortran::ice_table_vals();
  Real result;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Real& value) {
    value = P3F::apply_table_ice(adjusted_index, ice_table_vals, tab)[0];
  }, result);
  *proc = result;
}

void access_lookup_table_coll_f(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::TableIce ti;
  typename P3F::TableRain tr;

  // Adjust for 0-based indexing
  ti.dumi  = dumi  - 1;
  ti.dumjj = dumjj - 1;
  ti.dumii = dumii - 1;
  tr.dumj  = dumj  - 1;

  int adjusted_index = index - 1;

  ti.dum1 = dum1;
  ti.dum4 = dum4;
  ti.dum5 = dum5;
  tr.dum3 = dum3;

  auto collect_table_vals = P3GlobalForFortran::collect_table_vals();
  Real result;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Real& value) {
    value = P3F::apply_table_coll(adjusted_index, collect_table_vals, ti, tr)[0];
  }, result);
  *proc = result;
}

void get_cloud_dsd2_f(Real qc_, Real* nc_, Real* mu_c_, Real rho_, Real* nu_, Real* lamc_,
                      Real* cdist_, Real* cdist1_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_d", 6);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nc = *nc_;
  const auto dnu = P3GlobalForFortran::dnu();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qc(qc_), nc(local_nc), rho(rho_);
    typename P3F::Spack mu_c, nu, lamc, cdist, cdist1;

    P3F::get_cloud_dsd2(qc, nc, mu_c, rho, nu, dnu, lamc, cdist, cdist1);

    t_d(0) = nc[0];
    t_d(1) = mu_c[0];
    t_d(2) = nu[0];
    t_d(3) = lamc[0];
    t_d(4) = cdist[0];
    t_d(5) = cdist1[0];
  });
  Kokkos::deep_copy(t_h, t_d);

  *nc_     = t_h(0);
  *mu_c_   = t_h(1);
  *nu_     = t_h(2);
  *lamc_   = t_h(3);
  *cdist_  = t_h(4);
  *cdist1_ = t_h(5);
}

void get_rain_dsd2_f(Real qr_, Real* nr_, Real* mu_r_, Real* lamr_, Real* cdistr_, Real* logn0r_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_d", 5);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_nr = *nr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qr(qr_), nr(local_nr);
    typename P3F::Spack lamr, mu_r, cdistr, logn0r;

    P3F::get_rain_dsd2(qr, nr, mu_r, lamr, cdistr, logn0r);

    t_d(0) = nr[0];
    t_d(1) = mu_r[0];
    t_d(2) = lamr[0];
    t_d(3) = cdistr[0];
    t_d(4) = logn0r[0];
  });
  Kokkos::deep_copy(t_h, t_d);

  *nr_     = t_h(0);
  *mu_r_   = t_h(1);
  *lamr_   = t_h(2);
  *cdistr_ = t_h(3);
  *logn0r_ = t_h(4);
}

void get_time_space_phys_variables_f(Real T_atm_, Real pres_, Real rho_, Real latent_heat_vapor_, Real latent_heat_sublim_, Real qv_sat_l_, Real qv_sat_i_,
				     Real* mu_, Real* dv_, Real* sc_, Real* dqsdt_, Real* dqsidt_, Real* ab_,
				     Real* abi_, Real* kap_, Real* eii_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 9);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack T_atm(T_atm_), pres(pres_), rho(rho_), latent_heat_vapor(latent_heat_vapor_), latent_heat_sublim(latent_heat_sublim_), qv_sat_l(qv_sat_l_), qv_sat_i(qv_sat_i_);
      typename P3F::Spack mu, dv, sc, dqsdt,dqsidt, ab, abi, kap, eii;

      P3F::get_time_space_phys_variables(T_atm, pres, rho, latent_heat_vapor, latent_heat_sublim, qv_sat_l, qv_sat_i, mu, dv, sc, dqsdt, dqsidt,
					 ab, abi, kap, eii);

      t_d(0) = mu[0];
      t_d(1) = dv[0];
      t_d(2) = sc[0];
      t_d(3) = dqsdt[0];
      t_d(4) = dqsidt[0];
      t_d(5) = ab[0];
      t_d(6) = abi[0];
      t_d(7) = kap[0];
      t_d(8) = eii[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *mu_     = t_h(0);
  *dv_     = t_h(1);
  *sc_     = t_h(2);
  *dqsdt_  = t_h(3);
  *dqsidt_ = t_h(4);
  *ab_     = t_h(5);
  *abi_    = t_h(6);
  *kap_    = t_h(7);
  *eii_    = t_h(8);
}

void update_prognostic_ice_f( Real qc2qi_hetero_freeze_tend_, Real qc2qi_collect_tend_, Real qc2qr_ice_shed_tend_,  Real nc_collect_tend_,  Real nc2ni_immers_freeze_tend_, Real ncshdc_,
                              Real qr2qi_collect_tend_,  Real nr_collect_tend_, Real qr2qi_immers_freeze_tend_, Real nr2ni_immers_freeze_tend_, Real nr_ice_shed_tend_,
                              Real qi2qr_melt_tend_, Real ni2nr_melt_tend_, Real qi2qv_sublim_tend_, Real qv2qi_vapdep_tend_, Real qv2qi_nucleat_tend_, Real ni_nucleat_tend_,
                              Real ni_selfcollect_tend_, Real ni_sublim_tend_, Real qc2qi_berg_tend_, Real exner_, Real latent_heat_sublim_, Real latent_heat_fusion_,
                              bool do_predict_nc_, bool log_wetgrowth_, Real dt_, Real nmltratio_,
                              Real rho_qm_cloud_, Real* th_atm_, Real* qv_, Real* qi_, Real* ni_, Real* qm_,
                              Real* bm_, Real* qc_, Real* nc_, Real* qr_, Real* nr_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 10);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_th    = *th_atm_;
  Real local_qv	   = *qv_;
  Real local_qc	   = *qc_;
  Real local_nc	   = *nc_;
  Real local_qr	   = *qr_;
  Real local_nr	   = *nr_;
  Real local_qi = *qi_;
  Real local_ni = *ni_;
  Real local_qm = *qm_;
  Real local_bm = *bm_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qc2qi_hetero_freeze_tend(qc2qi_hetero_freeze_tend_), qc2qi_collect_tend(qc2qi_collect_tend_),qc2qr_ice_shed_tend(qc2qr_ice_shed_tend_),  nc_collect_tend(nc_collect_tend_),
	nc2ni_immers_freeze_tend(nc2ni_immers_freeze_tend_),  ncshdc(ncshdc_),  qr2qi_collect_tend(qr2qi_collect_tend_),  nr_collect_tend(nr_collect_tend_),  qr2qi_immers_freeze_tend(qr2qi_immers_freeze_tend_),
	nr2ni_immers_freeze_tend(nr2ni_immers_freeze_tend_),  nr_ice_shed_tend(nr_ice_shed_tend_),  qi2qr_melt_tend(qi2qr_melt_tend_),  ni2nr_melt_tend(ni2nr_melt_tend_),  qi2qv_sublim_tend(qi2qv_sublim_tend_),
	qv2qi_vapdep_tend(qv2qi_vapdep_tend_),  qv2qi_nucleat_tend(qv2qi_nucleat_tend_),  ni_nucleat_tend(ni_nucleat_tend_),  ni_selfcollect_tend(ni_selfcollect_tend_),  ni_sublim_tend(ni_sublim_tend_),
	qc2qi_berg_tend(qc2qi_berg_tend_),  exner(exner_),  latent_heat_fusion(latent_heat_fusion_),  latent_heat_sublim(latent_heat_sublim_),
	rho_qm_cloud(rho_qm_cloud_);
      bool do_predict_nc(do_predict_nc_);
      typename P3F::Smask log_wetgrowth(log_wetgrowth_);
      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th_atm(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr),
	nr(local_nr), qi(local_qi), ni(local_ni), qm(local_qm), bm(local_bm);

      P3F::update_prognostic_ice(qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend,ncshdc,
				 qr2qi_collect_tend,   nr_collect_tend,  qr2qi_immers_freeze_tend,  nr2ni_immers_freeze_tend,  nr_ice_shed_tend,
				 qi2qr_melt_tend,  ni2nr_melt_tend,  qi2qv_sublim_tend,  qv2qi_vapdep_tend,  qv2qi_nucleat_tend,  ni_nucleat_tend,
				 ni_selfcollect_tend,  ni_sublim_tend,  qc2qi_berg_tend,  exner,  latent_heat_sublim,  latent_heat_fusion,
				 do_predict_nc, log_wetgrowth,  dt,  nmltratio_,
				 rho_qm_cloud, th_atm, qv, qi, ni, qm,
				 bm, qc, nc, qr, nr);


      t_d(0) = th_atm[0];
      t_d(1) = qv[0];
      t_d(2) = qi[0];
      t_d(3) = ni[0];
      t_d(4) = qm[0];
      t_d(5) = bm[0];
      t_d(6) = qc[0];
      t_d(7) = nc[0];
      t_d(8) = qr[0];
      t_d(9) = nr[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *th_atm_    = t_h(0);
  *qv_    = t_h(1);
  *qi_ = t_h(2);
  *ni_ = t_h(3);
  *qm_ = t_h(4);
  *bm_ = t_h(5);
  *qc_    = t_h(6);
  *nc_    = t_h(7);
  *qr_    = t_h(8);
  *nr_    = t_h(9);
}

void evaporate_rain_f(Real qr_incld_, Real qc_incld_, Real nr_incld_, Real qi_incld_,
		      Real cld_frac_l_, Real cld_frac_r_, Real qv_, Real qv_prev_,
		      Real qv_sat_l_, Real qv_sat_i_, Real ab_, Real abi_,
		      Real epsr_, Real epsi_tot_, Real t_, Real t_prev_,
		      Real latent_heat_sublim_, Real dqsdt_, Real dt_,
		      Real* qr2qv_evap_tend_, Real* nr_evap_tend_)

{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qr2qv_evap_tend = *qr2qv_evap_tend_;
  Real local_nr_evap_tend = *nr_evap_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qr_incld(qr_incld_), qc_incld(qc_incld_), nr_incld(nr_incld_), qi_incld(qi_incld_),
	cld_frac_l(cld_frac_l_), cld_frac_r(cld_frac_r_), qv(qv_), qv_prev(qv_prev_), qv_sat_l(qv_sat_l_), qv_sat_i(qv_sat_i_),
        ab(ab_), abi(abi_), epsr(epsr_), epsi_tot(epsi_tot_), t(t_), t_prev(t_prev_), latent_heat_sublim(latent_heat_sublim_),
        dqsdt(dqsdt_);
      
      typename P3F::Scalar dt(dt_);
      
      typename P3F::Spack qr2qv_evap_tend(local_qr2qv_evap_tend), nr_evap_tend(local_nr_evap_tend);

      P3F::evaporate_rain(qr_incld,qc_incld,nr_incld,qi_incld,
			  cld_frac_l,cld_frac_r,qv,qv_prev,qv_sat_l,qv_sat_i,
			  ab,abi,epsr,epsi_tot,t,t_prev,latent_heat_sublim,dqsdt,dt,
			  qr2qv_evap_tend,nr_evap_tend);

      t_d(0) = qr2qv_evap_tend[0];
      t_d(1) = nr_evap_tend[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qr2qv_evap_tend_ = t_h(0);
  *nr_evap_tend_ = t_h(1);
}

void update_prognostic_liquid_f(Real qc2qr_accret_tend_, Real nc_accret_tend_, Real qc2qr_autoconv_tend_, Real nc2nr_autoconv_tend_, Real ncautr_,
				Real nc_selfcollect_tend_, Real  qr2qv_evap_tend_, Real nr_evap_tend_, Real nr_selfcollect_tend_, bool do_predict_nc_,
				Real inv_rho_, Real exner_, Real latent_heat_vapor_, Real dt_, Real* th_atm_, Real* qv_,
				Real* qc_, Real* nc_, Real* qr_, Real* nr_)

{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 6);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_th = *th_atm_;
  Real local_qv = *qv_;
  Real local_qc = *qc_;
  Real local_nc = *nc_;
  Real local_qr = *qr_;
  Real local_nr = *nr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qc2qr_accret_tend(qc2qr_accret_tend_), nc_accret_tend(nc_accret_tend_), qc2qr_autoconv_tend(qc2qr_autoconv_tend_), nc2nr_autoconv_tend(nc2nr_autoconv_tend_),
	ncautr(ncautr_), nc_selfcollect_tend(nc_selfcollect_tend_),  qr2qv_evap_tend( qr2qv_evap_tend_), nr_evap_tend(nr_evap_tend_), nr_selfcollect_tend(nr_selfcollect_tend_), inv_rho(inv_rho_),
	exner(exner_), latent_heat_vapor(latent_heat_vapor_);

      bool do_predict_nc(do_predict_nc_);

      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th_atm(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr), nr(local_nr);

      P3F::update_prognostic_liquid(qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr,
				    nc_selfcollect_tend,  qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend , do_predict_nc,
				    inv_rho, exner, latent_heat_vapor, dt, th_atm, qv,
				    qc, nc, qr, nr);

      t_d(0) = th_atm[0];
      t_d(1) = qv[0];
      t_d(2) = qc[0];
      t_d(3) = nc[0];
      t_d(4) = qr[0];
      t_d(5) = nr[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *th_atm_    = t_h(0);
  *qv_    = t_h(1);
  *qc_    = t_h(2);
  *nc_    = t_h(3);
  *qr_    = t_h(4);
  *nr_    = t_h(5);
}

void ice_deposition_sublimation_f(Real qi_incld_, Real ni_incld_, Real T_atm_, Real qv_sat_l_,
				  Real qv_sat_i_, Real epsi_, Real abi_, Real qv_,
				  Real* qv2qi_vapdep_tend_, Real* qi2qv_sublim_tend_, Real* ni_sublim_tend_, Real* qc2qi_berg_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 4);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qv2qi_vapdep_tend  = *qv2qi_vapdep_tend_;
  Real local_qi2qv_sublim_tend  = *qi2qv_sublim_tend_;
  Real local_ni_sublim_tend  = *ni_sublim_tend_;
  Real local_qc2qi_berg_tend = *qc2qi_berg_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qi_incld(qi_incld_), ni_incld(ni_incld_), T_atm(T_atm_), qv_sat_l(qv_sat_l_), qv_sat_i(qv_sat_i_),
	epsi(epsi_), abi(abi_), qv(qv_);

      typename P3F::Spack qv2qi_vapdep_tend(local_qv2qi_vapdep_tend), qi2qv_sublim_tend(local_qi2qv_sublim_tend), ni_sublim_tend(local_ni_sublim_tend), qc2qi_berg_tend(local_qc2qi_berg_tend);

      P3F::ice_deposition_sublimation(qi_incld, ni_incld, T_atm, qv_sat_l, qv_sat_i, epsi, abi, qv,
				      qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend);

      t_d(0) = qv2qi_vapdep_tend[0];
      t_d(1) = qi2qv_sublim_tend[0];
      t_d(2) = ni_sublim_tend[0];
      t_d(3) = qc2qi_berg_tend[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qv2qi_vapdep_tend_  = t_h(0);
  *qi2qv_sublim_tend_  = t_h(1);
  *ni_sublim_tend_  = t_h(2);
  *qc2qi_berg_tend_ = t_h(3);
}

template <int N, typename T>
Kokkos::Array<T*, N> ptr_to_arr(T** data)
{
  Kokkos::Array<T*, N> result;
  for (int i = 0; i < N; ++i) result[i] = data[i];

  return result;
}

template <int N>
void calc_first_order_upwind_step_f_impl(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub,
  Real* rho, Real* inv_rho, Real* inv_dz,
  Real** fluxes, Real** vs, Real** qnx)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Spack>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using view_1d_ptr_array = typename P3F::view_1d_ptr_array<Spack, N>;
  using uview_1d = typename P3F::uview_1d<Spack>;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Setup views
  Kokkos::Array<view_1d, 3> temp_d;
  Kokkos::Array<view_1d, N> fluxes_d, vs_d, qnx_d;

  ekat::host_to_device({rho, inv_rho, inv_dz}, nk, temp_d);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dz_d(temp_d[2]);

  ekat::host_to_device(ptr_to_arr<N>((const Real**)fluxes), nk, fluxes_d);
  ekat::host_to_device(ptr_to_arr<N>((const Real**)vs)    , nk, vs_d);
  ekat::host_to_device(ptr_to_arr<N>((const Real**)qnx)   , nk, qnx_d);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_d[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_d[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_d[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dz_d(inv_dz_d);
    P3F::calc_first_order_upwind_step<N>(urho_d, uinv_rho_d, uinv_dz_d, team, nk, kbot, k_qxtop, kdir, dt_sub, fluxes_ptr, vs_ptr, qnx_ptr);
  });

  // Sync back to host
  ekat::device_to_host(ptr_to_arr<N>(fluxes), nk, fluxes_d);
  ekat::device_to_host(ptr_to_arr<N>(qnx), nk, qnx_d);
}

template <int N>
void generalized_sedimentation_f_impl(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
  Real** vs, Real** fluxes, Real** qnx)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using Singlep = typename ekat::Pack<Real, 1>;
  using view_1d = typename P3F::view_1d<Spack>;
  using view_1ds = typename P3F::view_1d<Singlep>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using view_1d_ptr_array = typename P3F::view_1d_ptr_array<Spack, N>;
  using uview_1d = typename P3F::uview_1d<Spack>;
  using ekat::host_to_device;
  using ekat::device_to_host;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;
  *k_qxbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, 3> temp_d;
  Kokkos::Array<view_1d, N> fluxes_d, vs_d, qnx_d;
  Kokkos::Array<view_1ds, 1> scalar_temp;
  std::vector<Real> scalars = {*prt_accum, *dt_left, static_cast<Real>(*k_qxbot)};

  host_to_device({rho, inv_rho, inv_dz}, nk, temp_d);
  host_to_device({scalars.data()}, scalars.size(), scalar_temp);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dz_d(temp_d[2]);
  view_1ds scalars_d(scalar_temp[0]);

  host_to_device(ptr_to_arr<N>((const Real**)fluxes), nk, fluxes_d);
  host_to_device(ptr_to_arr<N>((const Real**)vs)    , nk, vs_d);
  host_to_device(ptr_to_arr<N>((const Real**)qnx)   , nk, qnx_d);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_d[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_d[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_d[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dz_d(inv_dz_d);

    // Each thread needs their own copy, like we expect in the main program, or else we will hit
    // data race issues
    Real prt_accum_k = scalars_d(0)[0];
    Real dt_left_k   = scalars_d(1)[0];
    Int k_qxbot_k    = static_cast<int>(scalars_d(2)[0]);

    P3F::generalized_sedimentation<N>(urho_d, uinv_rho_d, uinv_dz_d, team, nk, k_qxtop, k_qxbot_k, kbot, kdir, Co_max, dt_left_k, prt_accum_k, fluxes_ptr, vs_ptr, qnx_ptr);

    scalars_d(0)[0] = prt_accum_k;
    scalars_d(1)[0] = dt_left_k;
    scalars_d(2)[0] = k_qxbot_k;
  });

  // Sync back to host
  device_to_host(ptr_to_arr<N>(fluxes), nk, fluxes_d);
  device_to_host(ptr_to_arr<N>(qnx), nk, qnx_d);
  device_to_host({scalars.data()}, scalars.size(), scalar_temp);

  // Set scalars
  *prt_accum = scalars[0];
  *dt_left   = scalars[1];
  *k_qxbot   = scalars[2] + 1;
}

void calc_first_order_upwind_step_f(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub,
  Real* rho, Real* inv_rho, Real* inv_dz,
  Int num_arrays, Real** fluxes, Real** vs, Real** qnx)
{
  if (num_arrays == 1) {
    calc_first_order_upwind_step_f_impl<1>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else if (num_arrays == 2) {
    calc_first_order_upwind_step_f_impl<2>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else if (num_arrays == 4) {
    calc_first_order_upwind_step_f_impl<4>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else {
    EKAT_REQUIRE_MSG(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void generalized_sedimentation_f(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
  Int num_arrays, Real** vs, Real** fluxes, Real** qnx)
{
  if (num_arrays == 1) {
    generalized_sedimentation_f_impl<1>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 2) {
    generalized_sedimentation_f_impl<2>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 4) {
    generalized_sedimentation_f_impl<4>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else {
    EKAT_REQUIRE_MSG(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void cloud_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Spack>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  const auto dnu = P3GlobalForFortran::dnu();

  Kokkos::Array<view_1d, CloudSedData::NUM_ARRAYS> temp_d;

  ekat::host_to_device({qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz, qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend},
                       nk, temp_d);

  view_1d
    qc_incld_d(temp_d[0]),
    rho_d     (temp_d[1]),
    inv_rho_d (temp_d[2]),
    cld_frac_l_d   (temp_d[3]),
    acn_d     (temp_d[4]),
    inv_dz_d (temp_d[5]),
    qc_d      (temp_d[6]),
    nc_d      (temp_d[7]),
    nc_incld_d(temp_d[8]),
    mu_c_d    (temp_d[9]),
    lamc_d    (temp_d[10]),
    qc_tend_d (temp_d[11]),
    nc_tend_d (temp_d[12]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_liq_surf_k) {

    P3F::cloud_sedimentation(
      qc_incld_d, rho_d, inv_rho_d, cld_frac_l_d, acn_d, inv_dz_d, dnu,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, inv_dt, do_predict_nc,
      qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d,
      precip_liq_surf_k);

  }, *precip_liq_surf);

  // Sync back to host
  Kokkos::Array<view_1d, 7> inout_views = {qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d};
  ekat::device_to_host({qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend}, nk, inout_views);
}

void ice_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, IceSedData::NUM_ARRAYS> temp_d;

  ekat::host_to_device({rho, inv_rho, rhofaci, cld_frac_i, inv_dz, qi, qi_incld, ni, qm, qm_incld, bm, bm_incld, ni_incld, qi_tend, ni_tend},
                       nk, temp_d);

  view_1d
    rho_d        (temp_d[0]),
    inv_rho_d    (temp_d[1]),
    rhofaci_d    (temp_d[2]),
    cld_frac_i_d (temp_d[3]),
    inv_dz_d     (temp_d[4]),
    qi_d         (temp_d[5]),
    qi_incld_d   (temp_d[6]),
    ni_d         (temp_d[7]),
    qm_d         (temp_d[8]),
    qm_incld_d   (temp_d[9]),
    bm_d         (temp_d[10]),
    bm_incld_d   (temp_d[11]),
    ni_incld_d   (temp_d[12]),
    qi_tend_d    (temp_d[13]),
    ni_tend_d    (temp_d[14]);

  // Call core function from kernel
  auto ice_table_vals = P3GlobalForFortran::ice_table_vals();
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 6, policy);
  Real my_precip_ice_surf = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_ice_surf_k) {

    P3F::ice_sedimentation(
      rho_d, inv_rho_d, rhofaci_d, cld_frac_i_d, inv_dz_d,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, inv_dt,
      qi_d, qi_incld_d, ni_d, ni_incld_d, qm_d, qm_incld_d, bm_d, bm_incld_d,
      qi_tend_d, ni_tend_d, ice_table_vals,
      precip_ice_surf_k);

  }, my_precip_ice_surf);
  *precip_ice_surf += my_precip_ice_surf;

  // Sync back to host
  Kokkos::Array<view_1d, 10> inout_views = {qi_d, qi_incld_d, ni_d, ni_incld_d, qm_d, qm_incld_d,
                                            bm_d, bm_incld_d, qi_tend_d, ni_tend_d};
  ekat::device_to_host({qi, qi_incld, ni, ni_incld, qm, qm_incld, bm, bm_incld, qi_tend, ni_tend}, nk, inout_views);
}

void rain_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, RainSedData::NUM_ARRAYS> temp_d;
  Kokkos::Array<size_t, RainSedData::NUM_ARRAYS> sizes;
  for (size_t i = 0; i < RainSedData::NUM_ARRAYS; ++i) sizes[i] = nk;
  sizes[RainSedData::NUM_ARRAYS - 1] = nk+1;

  ekat::host_to_device({qr_incld, rho, inv_rho, rhofacr, cld_frac_r, inv_dz, qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, precip_liq_flux},
                       sizes, temp_d);

  view_1d
    qr_incld_d       (temp_d[0]),
    rho_d            (temp_d[1]),
    inv_rho_d        (temp_d[2]),
    rhofacr_d        (temp_d[3]),
    cld_frac_r_d     (temp_d[4]),
    inv_dz_d         (temp_d[5]),
    qr_d             (temp_d[6]),
    nr_d             (temp_d[7]),
    nr_incld_d       (temp_d[8]),
    mu_r_d           (temp_d[9]),
    lamr_d           (temp_d[10]),
    qr_tend_d        (temp_d[11]),
    nr_tend_d        (temp_d[12]),
    precip_liq_flux_d(temp_d[13]);

  // Call core function from kernel
  auto vn_table_vals = P3GlobalForFortran::vn_table_vals();
  auto vm_table_vals = P3GlobalForFortran::vm_table_vals();
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Real my_precip_liq_surf = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_liq_surf_k) {

    P3F::rain_sedimentation(
      rho_d, inv_rho_d, rhofacr_d, cld_frac_r_d, inv_dz_d, qr_incld_d,
      team, wsm.get_workspace(team), vn_table_vals, vm_table_vals,
      nk, ktop, kbot, kdir, dt, inv_dt,
      qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, precip_liq_flux_d, qr_tend_d, nr_tend_d,
      precip_liq_surf_k);

  }, my_precip_liq_surf);
  *precip_liq_surf += my_precip_liq_surf;

  // Sync back to host
  Kokkos::Array<size_t, 8> sizes_out;
  for (int i = 0; i < 8; ++i) sizes_out[i] = nk;
  sizes_out[7] = nk+1;

  Kokkos::Array<view_1d, 8> inout_views = {qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, qr_tend_d, nr_tend_d, precip_liq_flux_d};
  ekat::device_to_host({qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, precip_liq_flux}, sizes_out, inout_views);
}

void back_to_cell_average_f(Real cld_frac_l_, Real cld_frac_r_, Real cld_frac_i_,
                            Real* qc2qr_accret_tend_, Real* qr2qv_evap_tend_, Real* qc2qr_autoconv_tend_,
                            Real* nc_accret_tend_, Real* nc_selfcollect_tend_, Real* nc2nr_autoconv_tend_,
                            Real* nr_selfcollect_tend_, Real* nr_evap_tend_, Real* ncautr_,
                            Real* qi2qv_sublim_tend_,
                            Real* nr_ice_shed_tend_, Real* qc2qi_hetero_freeze_tend_, Real* qr2qi_collect_tend_,
                            Real* qc2qr_ice_shed_tend_, Real* qi2qr_melt_tend_, Real* qc2qi_collect_tend_,
                            Real* qr2qi_immers_freeze_tend_, Real* ni2nr_melt_tend_, Real* nc_collect_tend_,
                            Real* ncshdc_, Real* nc2ni_immers_freeze_tend_, Real* nr_collect_tend_,
                            Real* ni_selfcollect_tend_, Real* qv2qi_vapdep_tend_, Real* nr2ni_immers_freeze_tend_,
                            Real* ni_sublim_tend_, Real* qv2qi_nucleat_tend_, Real* ni_nucleat_tend_,
                            Real* qc2qi_berg_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 29);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qc2qr_accret_tend = *qc2qr_accret_tend_;
  Real local_qr2qv_evap_tend = *qr2qv_evap_tend_;
  Real local_qc2qr_autoconv_tend = *qc2qr_autoconv_tend_;
  Real local_nc_accret_tend = *nc_accret_tend_;
  Real local_nc_selfcollect_tend = *nc_selfcollect_tend_;
  Real local_nc2nr_autoconv_tend = *nc2nr_autoconv_tend_;
  Real local_nr_selfcollect_tend = *nr_selfcollect_tend_;
  Real local_nr_evap_tend = *nr_evap_tend_;
  Real local_ncautr = *ncautr_;
  Real local_qi2qv_sublim_tend = *qi2qv_sublim_tend_;
  Real local_nr_ice_shed_tend = *nr_ice_shed_tend_;
  Real local_qc2qi_hetero_freeze_tend = *qc2qi_hetero_freeze_tend_;
  Real local_qr2qi_collect_tend = *qr2qi_collect_tend_;
  Real local_qc2qr_ice_shed_tend = *qc2qr_ice_shed_tend_;
  Real local_qi2qr_melt_tend = *qi2qr_melt_tend_;
  Real local_qc2qi_collect_tend = *qc2qi_collect_tend_;
  Real local_qr2qi_immers_freeze_tend = *qr2qi_immers_freeze_tend_;
  Real local_ni2nr_melt_tend = *ni2nr_melt_tend_;
  Real local_nc_collect_tend = *nc_collect_tend_;
  Real local_ncshdc = *ncshdc_;
  Real local_nc2ni_immers_freeze_tend = *nc2ni_immers_freeze_tend_;
  Real local_nr_collect_tend = *nr_collect_tend_;
  Real local_ni_selfcollect_tend = *ni_selfcollect_tend_;
  Real local_qv2qi_vapdep_tend = *qv2qi_vapdep_tend_;
  Real local_nr2ni_immers_freeze_tend = *nr2ni_immers_freeze_tend_;
  Real local_ni_sublim_tend = *ni_sublim_tend_;
  Real local_qv2qi_nucleat_tend = *qv2qi_nucleat_tend_;
  Real local_ni_nucleat_tend = *ni_nucleat_tend_;
  Real local_qc2qi_berg_tend = *qc2qi_berg_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack cld_frac_l(cld_frac_l_), cld_frac_r(cld_frac_r_), cld_frac_i(cld_frac_i_),
      qc2qr_accret_tend(local_qc2qr_accret_tend), qr2qv_evap_tend(local_qr2qv_evap_tend), qc2qr_autoconv_tend(local_qc2qr_autoconv_tend), nc_accret_tend(local_nc_accret_tend),
      nc_selfcollect_tend(local_nc_selfcollect_tend), nc2nr_autoconv_tend(local_nc2nr_autoconv_tend), nr_selfcollect_tend(local_nr_selfcollect_tend), nr_evap_tend(local_nr_evap_tend),
      ncautr(local_ncautr), qi2qv_sublim_tend(local_qi2qv_sublim_tend),
      nr_ice_shed_tend(local_nr_ice_shed_tend), qc2qi_hetero_freeze_tend(local_qc2qi_hetero_freeze_tend), qr2qi_collect_tend(local_qr2qi_collect_tend), qc2qr_ice_shed_tend(local_qc2qr_ice_shed_tend),
      qi2qr_melt_tend(local_qi2qr_melt_tend), qc2qi_collect_tend(local_qc2qi_collect_tend), qr2qi_immers_freeze_tend(local_qr2qi_immers_freeze_tend), ni2nr_melt_tend(local_ni2nr_melt_tend),
      nc_collect_tend(local_nc_collect_tend), ncshdc(local_ncshdc), nc2ni_immers_freeze_tend(local_nc2ni_immers_freeze_tend), nr_collect_tend(local_nr_collect_tend),
      ni_selfcollect_tend(local_ni_selfcollect_tend), qv2qi_vapdep_tend(local_qv2qi_vapdep_tend), nr2ni_immers_freeze_tend(local_nr2ni_immers_freeze_tend), ni_sublim_tend(local_ni_sublim_tend),
      qv2qi_nucleat_tend(local_qv2qi_nucleat_tend), ni_nucleat_tend(local_ni_nucleat_tend), qc2qi_berg_tend(local_qc2qi_berg_tend);

    P3F::back_to_cell_average(cld_frac_l, cld_frac_r, cld_frac_i, qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,
      nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend,
      nr_ice_shed_tend, qc2qi_hetero_freeze_tend, qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend,
      ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend, qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend,
      qc2qi_berg_tend);

    t_d(0) = qc2qr_accret_tend[0];
    t_d(1) = qr2qv_evap_tend[0];
    t_d(2) = qc2qr_autoconv_tend[0];
    t_d(3) = nc_accret_tend[0];
    t_d(4) = nc_selfcollect_tend[0];
    t_d(5) = nc2nr_autoconv_tend[0];
    t_d(6) = nr_selfcollect_tend[0];
    t_d(7) = nr_evap_tend[0];
    t_d(8) = ncautr[0];
    t_d(9) = qi2qv_sublim_tend[0];
    t_d(10) = nr_ice_shed_tend[0];
    t_d(11) = qc2qi_hetero_freeze_tend[0];
    t_d(12) = qr2qi_collect_tend[0];
    t_d(13) = qc2qr_ice_shed_tend[0];
    t_d(14) = qi2qr_melt_tend[0];
    t_d(15) = qc2qi_collect_tend[0];
    t_d(16) = qr2qi_immers_freeze_tend[0];
    t_d(17) = ni2nr_melt_tend[0];
    t_d(18) = nc_collect_tend[0];
    t_d(19) = ncshdc[0];
    t_d(20) = nc2ni_immers_freeze_tend[0];
    t_d(21) = nr_collect_tend[0];
    t_d(22) = ni_selfcollect_tend[0];
    t_d(23) = qv2qi_vapdep_tend[0];
    t_d(24) = nr2ni_immers_freeze_tend[0];
    t_d(25) = ni_sublim_tend[0];
    t_d(26) = qv2qi_nucleat_tend[0];
    t_d(27) = ni_nucleat_tend[0];
    t_d(28) = qc2qi_berg_tend[0];

  });
  Kokkos::deep_copy(t_h, t_d);

  *qc2qr_accret_tend_        = t_h(0);
  *qr2qv_evap_tend_          = t_h(1);
  *qc2qr_autoconv_tend_      = t_h(2);
  *nc_accret_tend_           = t_h(3);
  *nc_selfcollect_tend_      = t_h(4);
  *nc2nr_autoconv_tend_      = t_h(5);
  *nr_selfcollect_tend_      = t_h(6);
  *nr_evap_tend_             = t_h(7);
  *ncautr_                   = t_h(8);
  *qi2qv_sublim_tend_        = t_h(9);
  *nr_ice_shed_tend_         = t_h(10);
  *qc2qi_hetero_freeze_tend_ = t_h(11);
  *qr2qi_collect_tend_       = t_h(12);
  *qc2qr_ice_shed_tend_      = t_h(13);
  *qi2qr_melt_tend_          = t_h(14);
  *qc2qi_collect_tend_       = t_h(15);
  *qr2qi_immers_freeze_tend_ = t_h(16);
  *ni2nr_melt_tend_          = t_h(17);
  *nc_collect_tend_          = t_h(18);
  *ncshdc_                   = t_h(19);
  *nc2ni_immers_freeze_tend_ = t_h(20);
  *nr_collect_tend_          = t_h(21);
  *ni_selfcollect_tend_      = t_h(22);
  *qv2qi_vapdep_tend_        = t_h(23);
  *nr2ni_immers_freeze_tend_ = t_h(24);
  *ni_sublim_tend_           = t_h(25);
  *qv2qi_nucleat_tend_       = t_h(26);
  *ni_nucleat_tend_          = t_h(27);
  *qc2qi_berg_tend_          = t_h(28);
}

void prevent_ice_overdepletion_f(
  Real pres_, Real T_atm_, Real qv_, Real latent_heat_sublim_, Real inv_dt_, Real* qv2qi_vapdep_tend_,
  Real* qi2qv_sublim_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qv2qi_vapdep_tend = *qv2qi_vapdep_tend_;
  Real local_qi2qv_sublim_tend = *qi2qv_sublim_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack pres(pres_), T_atm(T_atm_), qv(qv_), latent_heat_sublim(latent_heat_sublim_),
      qv2qi_vapdep_tend(local_qv2qi_vapdep_tend), qi2qv_sublim_tend(local_qi2qv_sublim_tend);
    P3F::prevent_ice_overdepletion(pres, T_atm, qv, latent_heat_sublim, inv_dt_, qv2qi_vapdep_tend, qi2qv_sublim_tend);

    t_d(0) = qv2qi_vapdep_tend[0];
    t_d(1) = qi2qv_sublim_tend[0];

  });
  Kokkos::deep_copy(t_h, t_d);

  *qv2qi_vapdep_tend_ = t_h(0);
  *qi2qv_sublim_tend_ = t_h(1);
}

void calc_rime_density_f(
  Real T_atm_, Real rhofaci_, Real table_val_qi_fallspd_, Real acn_, Real lamc_, Real mu_c_,
  Real qc_incld_, Real qc2qi_collect_tend_, Real* vtrmi1_, Real* rho_qm_cloud_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_vtrmi1 = *vtrmi1_;
  Real local_rho_qm_cloud = *rho_qm_cloud_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack T_atm(T_atm_), rhofaci(rhofaci_), table_val_qi_fallspd(table_val_qi_fallspd_), acn(acn_),
                          lamc(lamc_), mu_c(mu_c_), qc_incld(qc_incld_),
                          qc2qi_collect_tend(qc2qi_collect_tend_), vtrmi1(local_vtrmi1),
                          rho_qm_cloud(local_rho_qm_cloud);
      P3F::calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld,
                             qc2qi_collect_tend, vtrmi1, rho_qm_cloud);

      t_d(0) = vtrmi1[0];
      t_d(1) = rho_qm_cloud[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *vtrmi1_       = t_h(0);
  *rho_qm_cloud_ = t_h(1);
}

void cldliq_immersion_freezing_f(
  Real T_atm_, Real lamc_, Real mu_c_, Real cdist1_, Real qc_incld_, Real inv_qc_relvar_,
  Real* qc2qi_hetero_freeze_tend_, Real* nc2ni_immers_freeze_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qc2qi_hetero_freeze_tend = *qc2qi_hetero_freeze_tend_;
  Real local_nc2ni_immers_freeze_tend = *nc2ni_immers_freeze_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack T_atm(T_atm_), lamc(lamc_), mu_c(mu_c_), cdist1(cdist1_),qc_incld(qc_incld_),
	                  inv_qc_relvar(inv_qc_relvar_),qc2qi_hetero_freeze_tend(local_qc2qi_hetero_freeze_tend), nc2ni_immers_freeze_tend(local_nc2ni_immers_freeze_tend);
      P3F::cldliq_immersion_freezing(T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar,
                                     qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend);

      t_d(0) = qc2qi_hetero_freeze_tend[0];
      t_d(1) = nc2ni_immers_freeze_tend[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qc2qi_hetero_freeze_tend_ = t_h(0);
  *nc2ni_immers_freeze_tend_ = t_h(1);
}

void rain_immersion_freezing_f(
  Real T_atm_, Real lamr_, Real mu_r_, Real cdistr_, Real qr_incld_,
  Real* qr2qi_immers_freeze_tend_, Real* nr2ni_immers_freeze_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qr2qi_immers_freeze_tend = *qr2qi_immers_freeze_tend_;
  Real local_nr2ni_immers_freeze_tend = *nr2ni_immers_freeze_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack T_atm(T_atm_), lamr(lamr_), mu_r(mu_r_),
                          cdistr(cdistr_), qr_incld(qr_incld_),
                          qr2qi_immers_freeze_tend(local_qr2qi_immers_freeze_tend), nr2ni_immers_freeze_tend(local_nr2ni_immers_freeze_tend);
      P3F::rain_immersion_freezing(T_atm, lamr, mu_r, cdistr, qr_incld,
                                   qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend);

      t_d(0) = qr2qi_immers_freeze_tend[0];
      t_d(1) = nr2ni_immers_freeze_tend[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qr2qi_immers_freeze_tend_ = t_h(0);
  *nr2ni_immers_freeze_tend_ = t_h(1);
}

void droplet_self_collection_f(
  Real rho_, Real inv_rho_, Real qc_incld_, Real mu_c_, Real nu_,
  Real nc2nr_autoconv_tend_, Real* nc_selfcollect_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_nc_selfcollect_tend = *nc_selfcollect_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), inv_rho(inv_rho_), qc_incld(qc_incld_),
                          mu_c(mu_c_), nu(nu_), nc2nr_autoconv_tend(nc2nr_autoconv_tend_),
                          nc_selfcollect_tend(local_nc_selfcollect_tend);
      P3F::droplet_self_collection(rho, inv_rho, qc_incld, mu_c, nu, nc2nr_autoconv_tend,
                                   nc_selfcollect_tend);

      t_d(0) = nc_selfcollect_tend[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *nc_selfcollect_tend_ = t_h(0);
}

void cloud_rain_accretion_f(
  Real rho_, Real inv_rho_, Real qc_incld_, Real nc_incld_, Real qr_incld_, Real inv_qc_relvar_,
  Real* qc2qr_accret_tend_, Real* nc_accret_tend_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qc2qr_accret_tend = *qc2qr_accret_tend_;
  Real local_nc_accret_tend = *nc_accret_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), inv_rho(inv_rho_), qc_incld(qc_incld_),
                          nc_incld(nc_incld_), qr_incld(qr_incld_),
	                  qc2qr_accret_tend(local_qc2qr_accret_tend), nc_accret_tend(local_nc_accret_tend), inv_qc_relvar(inv_qc_relvar_);
      P3F::cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld, inv_qc_relvar,
                                qc2qr_accret_tend, nc_accret_tend);

      t_d(0) = qc2qr_accret_tend[0];
      t_d(1) = nc_accret_tend[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qc2qr_accret_tend_ = t_h(0);
  *nc_accret_tend_    = t_h(1);
}

void cloud_water_autoconversion_f(
     Real rho_, Real qc_incld_, Real nc_incld_, Real inv_qc_relvar_,
     Real* qc2qr_autoconv_tend_, Real* nc2nr_autoconv_tend_, Real* ncautr_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 3);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qc2qr_autoconv_tend = *qc2qr_autoconv_tend_;
  Real local_nc2nr_autoconv_tend = *nc2nr_autoconv_tend_;
  Real local_ncautr              = *ncautr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), qc_incld(qc_incld_), nc_incld(nc_incld_), qc2qr_autoconv_tend(local_qc2qr_autoconv_tend),
	nc2nr_autoconv_tend(local_nc2nr_autoconv_tend), ncautr(local_ncautr), inv_qc_relvar(inv_qc_relvar_);
      P3F::cloud_water_autoconversion(rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr);

      t_d(0) = qc2qr_autoconv_tend[0];
      t_d(1) = nc2nr_autoconv_tend[0];
      t_d(2) = ncautr[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qc2qr_autoconv_tend_ = t_h(0);
  *nc2nr_autoconv_tend_ = t_h(1);
  *ncautr_              = t_h(2);
}

void rain_self_collection_f(Real rho_, Real qr_incld_, Real nr_incld_, Real* nr_selfcollect_tend_){
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nr_selfcollect_tend = *nr_selfcollect_tend_;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), qr_incld(qr_incld_), nr_incld(nr_incld_),  nr_selfcollect_tend(local_nr_selfcollect_tend);
      P3F::rain_self_collection(rho, qr_incld, nr_incld, nr_selfcollect_tend);

      t_d(0) = nr_selfcollect_tend[0];

    });

  Kokkos::deep_copy(t_h, t_d);
  *nr_selfcollect_tend_ = t_h(0);
}

  void ice_melting_f(Real rho_,Real T_atm_,Real pres_,Real rhofaci_,Real table_val_qi2qr_melting_,Real table_val_qi2qr_vent_melt_,Real latent_heat_vapor_,Real latent_heat_fusion_,Real dv_,Real sc_,Real mu_,Real kap_,Real qv_,Real qi_incld_,Real ni_incld_,Real* qi2qr_melt_tend_,Real* ni2nr_melt_tend_){
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qi2qr_melt_tend = *qi2qr_melt_tend_;
  Real local_ni2nr_melt_tend = *ni2nr_melt_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), T_atm(T_atm_), pres(pres_), rhofaci(rhofaci_),table_val_qi2qr_melting(table_val_qi2qr_melting_), table_val_qi2qr_vent_melt(table_val_qi2qr_vent_melt_), latent_heat_vapor(latent_heat_vapor_), latent_heat_fusion(latent_heat_fusion_),dv(dv_), sc(sc_), mu(mu_), kap(kap_),qv(qv_), qi_incld(qi_incld_), ni_incld(ni_incld_), qi2qr_melt_tend(local_qi2qr_melt_tend), ni2nr_melt_tend(local_ni2nr_melt_tend);
      P3F::ice_melting(rho,T_atm,pres,rhofaci,table_val_qi2qr_melting,table_val_qi2qr_vent_melt,latent_heat_vapor,latent_heat_fusion,dv,sc,mu,kap,qv,qi_incld,ni_incld,qi2qr_melt_tend,ni2nr_melt_tend);

      t_d(0) = qi2qr_melt_tend[0];
      t_d(1) = ni2nr_melt_tend[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qi2qr_melt_tend_ = t_h(0);
  *ni2nr_melt_tend_ = t_h(1);
}

void impose_max_total_ni_f(Real* ni_local_, Real max_total_ni_, Real inv_rho_local_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;


  view_1d t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_ni_local = *ni_local_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack ni_local(local_ni_local);
    Spack inv_rho_local(inv_rho_local_);

    P3F::impose_max_total_ni(ni_local, max_total_ni_, inv_rho_local);
    t_d(0) = ni_local[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *ni_local_ = t_h(0);
}

void calc_bulk_rho_rime_f(Real qi_tot_, Real* qi_rim_, Real* bi_rim_, Real* rho_rime_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  Real local_qi_rim = *qi_rim_, local_bi_rim = *bi_rim_;
  view_1d t_d("t_d", 3);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack qi_tot(qi_tot_), qi_rim(local_qi_rim), bi_rim(local_bi_rim);

    const auto result = P3F::calc_bulk_rho_rime(qi_tot, qi_rim, bi_rim);
    t_d(0) = qi_rim[0];
    t_d(1) = bi_rim[0];
    t_d(2) = result[0];
  });
  Kokkos::deep_copy(t_h, t_d);

  *qi_rim_   = t_h(0);
  *bi_rim_   = t_h(1);
  *rho_rime_ = t_h(2);
}

void homogeneous_freezing_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* T_atm, Real* exner, Real* latent_heat_fusion,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, HomogeneousFreezingData::NUM_ARRAYS> temp_d;

  ekat::host_to_device({T_atm, exner, latent_heat_fusion, qc, nc, qr, nr, qi, ni, qm, bm, th_atm},
                       nk, temp_d);

  view_1d
    t_d                   (temp_d[0]),
    exner_d               (temp_d[1]),
    latent_heat_fusion_d  (temp_d[2]),
    qc_d                  (temp_d[3]),
    nc_d                  (temp_d[4]),
    qr_d                  (temp_d[5]),
    nr_d                  (temp_d[6]),
    qi_d                  (temp_d[7]),
    ni_d                  (temp_d[8]),
    qm_d                  (temp_d[9]),
    bm_d                  (temp_d[10]),
    th_atm_d              (temp_d[11]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::homogeneous_freezing(
      t_d, exner_d, latent_heat_fusion_d,
      team,
      nk, ktop, kbot, kdir,
      qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, th_atm_d);
  });

  // Sync back to host
  Kokkos::Array<view_1d, 9> inout_views = {qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, th_atm_d};

  ekat::device_to_host({qc, nc, qr, nr, qi, ni, qm, bm, th_atm}, nk, inout_views);
}

void compute_rain_fall_velocity_f(Real qr_incld_, Real rhofacr_,
                                  Real* nr_incld_, Real* mu_r_, Real* lamr_, Real* V_qr_, Real* V_nr_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  Real local_nr_incld = *nr_incld_;
  view_1d t_d("t_d", 5);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  const auto vn_table_vals = P3GlobalForFortran::vn_table_vals();
  const auto vm_table_vals = P3GlobalForFortran::vm_table_vals();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack qr_incld(qr_incld_), rhofacr(rhofacr_), nr_incld(local_nr_incld),
      mu_r, lamr, V_qr, V_nr;

    P3F::compute_rain_fall_velocity(vn_table_vals, vm_table_vals,
                                    qr_incld, rhofacr, nr_incld, mu_r, lamr, V_qr, V_nr);
    t_d(0) = nr_incld[0];
    t_d(1) = mu_r[0];
    t_d(2) = lamr[0];
    t_d(3) = V_qr[0];
    t_d(4) = V_nr[0];
  });
  Kokkos::deep_copy(t_h, t_d);

  *nr_incld_ = t_h(0);
  *mu_r_     = t_h(1);
  *lamr_     = t_h(2);
  *V_qr_     = t_h(3);
  *V_nr_     = t_h(4);
}

void ice_cldliq_collection_f(Real rho_, Real temp_, Real rhofaci_, Real table_val_qc2qi_collect_,
                             Real qi_incld_,Real qc_incld_, Real ni_incld_, Real nc_incld_,
                             Real* qc2qi_collect_tend_, Real* nc_collect_tend_, Real* qc2qr_ice_shed_tend_, Real* ncshdc_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 4);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, table_val_qc2qi_collect{table_val_qc2qi_collect_}, qi_incld{qi_incld_},
          qc_incld{qc_incld_}, ni_incld{ni_incld_}, nc_incld{nc_incld_};
    Spack qc2qi_collect_tend{0.}, nc_collect_tend{0.}, qc2qr_ice_shed_tend{0.}, ncshdc{0.};

    P3F::ice_cldliq_collection(rho, temp, rhofaci, table_val_qc2qi_collect, qi_incld, qc_incld, ni_incld, nc_incld,
                               qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc);

    t_d(0) = qc2qi_collect_tend[0];
    t_d(1) = nc_collect_tend[0];
    t_d(2) = qc2qr_ice_shed_tend[0];
    t_d(3) = ncshdc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qc2qi_collect_tend_     = t_h(0);
  *nc_collect_tend_        = t_h(1);
  *qc2qr_ice_shed_tend_    = t_h(2);
  *ncshdc_                 = t_h(3);
}

void ice_rain_collection_f(Real rho_, Real temp_, Real rhofaci_, Real logn0r_, Real table_val_nr_collect_, Real table_val_qr2qi_collect_,
                           Real qi_incld_, Real ni_incld_, Real qr_incld_, Real* qr2qi_collect_tend_, Real* nr_collect_tend_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, logn0r{logn0r_}, table_val_nr_collect{table_val_nr_collect_}, table_val_qr2qi_collect{table_val_qr2qi_collect_},
          qi_incld{qi_incld_}, qr_incld{qr_incld_}, ni_incld{ni_incld_};
    Spack qr2qi_collect_tend{0.}, nr_collect_tend{0.};

    P3F::ice_rain_collection(rho, temp, rhofaci, logn0r, table_val_nr_collect, table_val_qr2qi_collect,
                             qi_incld, ni_incld, qr_incld,
                             qr2qi_collect_tend, nr_collect_tend);

    t_d(0) = qr2qi_collect_tend[0];
    t_d(1) = nr_collect_tend[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qr2qi_collect_tend_  = t_h(0);
  *nr_collect_tend_     = t_h(1);
}


void ice_self_collection_f(Real rho_, Real rhofaci_, Real table_val_ni_self_collect_, Real eii_,
                           Real qm_incld_, Real qi_incld_, Real ni_incld_, Real* ni_selfcollect_tend_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, rhofaci{rhofaci_}, table_val_ni_self_collect{table_val_ni_self_collect_}, eii{eii_}, qm_incld{qm_incld_},
          qi_incld{qi_incld_}, ni_incld{ni_incld_};
    Spack ni_selfcollect_tend{0.};

    P3F::ice_self_collection(rho, rhofaci, table_val_ni_self_collect, eii, qm_incld, qi_incld, ni_incld,
                             ni_selfcollect_tend);

    t_d(0) = ni_selfcollect_tend[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *ni_selfcollect_tend_     = t_h(0);

}


void ice_relaxation_timescale_f(Real rho_, Real temp_, Real rhofaci_, Real table_val_qi2qr_melting_, Real table_val_qi2qr_vent_melt_,
                                Real dv_, Real mu_, Real sc_, Real qi_incld_, Real ni_incld_,
                                Real* epsi_, Real* epsi_tot_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, table_val_qi2qr_melting{table_val_qi2qr_melting_}, table_val_qi2qr_vent_melt{table_val_qi2qr_vent_melt_}, dv{dv_},
          mu{mu_}, sc{sc_}, qi_incld{qi_incld_}, ni_incld{ni_incld_};

    Spack epsi{0.0}, epsi_tot{0.0};

    P3F::ice_relaxation_timescale(rho, temp, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld, ni_incld,
                                  epsi, epsi_tot);

    t_d(0) = epsi[0];
    t_d(1) = epsi_tot[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *epsi_      = t_h(0);
  *epsi_tot_  = t_h(1);
}

void calc_liq_relaxation_timescale_f(Real rho_, Real f1r_, Real f2r_, Real dv_,
                                     Real mu_, Real sc_, Real mu_r_, Real lamr_,
                                     Real cdistr_, Real cdist_, Real qr_incld_,
                                     Real qc_incld_, Real* epsr_, Real* epsc_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);
  auto revap_table_vals = P3GlobalForFortran::revap_table_vals();

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, dv{dv_},
          mu{mu_}, sc{sc_}, mu_r{mu_r_}, lamr{lamr_}, cdistr{cdistr_},
          cdist{cdist_}, qr_incld{qr_incld_}, qc_incld{qc_incld_};

    Spack epsr{0.0}, epsc{0.0};

    P3F::calc_liq_relaxation_timescale(revap_table_vals, rho, f1r_, f2r_, dv, mu, sc,
      mu_r, lamr, cdistr, cdist, qr_incld, qc_incld, epsr, epsc);

    t_d(0) = epsr[0];
    t_d(1) = epsc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *epsr_ = t_h(0);
  *epsc_ = t_h(1);
}

void ice_nucleation_f(Real temp_, Real inv_rho_, Real ni_, Real ni_activated_,
                      Real qv_supersat_i_, Real inv_dt_, bool do_predict_nc_,
                      Real* qv2qi_nucleat_tend_, Real* ni_nucleat_tend_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using view_1d      = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack temp{temp_}, inv_rho{inv_rho_}, ni{ni_}, ni_activated{ni_activated_}, qv_supersat_i{qv_supersat_i_};
    Spack qv2qi_nucleat_tend{0.0}, ni_nucleat_tend{0.0};

    P3F::ice_nucleation(temp, inv_rho, ni, ni_activated, qv_supersat_i, inv_dt_, do_predict_nc_,
                        qv2qi_nucleat_tend, ni_nucleat_tend);

    t_d(0) = qv2qi_nucleat_tend[0];
    t_d(1) = ni_nucleat_tend[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qv2qi_nucleat_tend_         = t_h(0);
  *ni_nucleat_tend_         = t_h(1);
}

void ice_cldliq_wet_growth_f(Real rho_, Real temp_, Real pres_, Real rhofaci_, Real table_val_qi2qr_melting_,
                             Real table_val_qi2qr_vent_melt_, Real latent_heat_vapor_, Real latent_heat_fusion_, Real dv_,
                             Real kap_, Real mu_, Real sc_, Real qv_, Real qc_incld_,
                             Real qi_incld_, Real ni_incld_, Real qr_incld_, bool* log_wetgrowth_,
                             Real* qr2qi_collect_tend_, Real* qc2qi_collect_tend_, Real* qc_growth_rate_, Real* nr_ice_shed_tend_, Real* qc2qr_ice_shed_tend_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using Smask        = typename P3F::Smask;
  using view_1d      = typename P3F::view_1d<Real>;
  using bool_view_1d = typename P3F::view_1d<bool>;

  bool_view_1d b_d("b_d", 1);
  view_1d t_d("t_d", 5);
  const auto b_h = Kokkos::create_mirror_view(b_d);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  const bool log_wetgrowth_local = *log_wetgrowth_;
  Real local_qr2qi_collect_tend = *qr2qi_collect_tend_, local_qc2qi_collect_tend = *qc2qi_collect_tend_, local_qc_growth_rate = *qc_growth_rate_, 
       local_nr_ice_shed_tend = *nr_ice_shed_tend_, local_qc2qr_ice_shed_tend = *qc2qr_ice_shed_tend_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, pres{pres_}, rhofaci{rhofaci_}, table_val_qi2qr_melting{table_val_qi2qr_melting_}, 
          table_val_qi2qr_vent_melt{table_val_qi2qr_vent_melt_}, latent_heat_vapor{latent_heat_vapor_},
          latent_heat_fusion{latent_heat_fusion_}, dv{dv_}, kap{kap_}, mu{mu_}, sc{sc_}, qv{qv_}, qc_incld{qc_incld_}, qi_incld{qi_incld_},
          ni_incld{ni_incld_}, qr_incld{qr_incld_};

    Smask log_wetgrowth{log_wetgrowth_local};

    Spack qr2qi_collect_tend{local_qr2qi_collect_tend}, qc2qi_collect_tend{local_qc2qi_collect_tend}, qc_growth_rate{local_qc_growth_rate}, 
          nr_ice_shed_tend{local_nr_ice_shed_tend}, qc2qr_ice_shed_tend{local_qc2qr_ice_shed_tend};

    P3F::ice_cldliq_wet_growth(rho, temp, pres, rhofaci, table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor, latent_heat_fusion, dv, kap, mu, sc, qv, qc_incld,
                              qi_incld, ni_incld, qr_incld, log_wetgrowth,
                              qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend);

    b_d(0) = log_wetgrowth[0];
    t_d(0) = qr2qi_collect_tend[0];
    t_d(1) = qc2qi_collect_tend[0];
    t_d(2) = qc_growth_rate[0];
    t_d(3) = nr_ice_shed_tend[0];
    t_d(4) = qc2qr_ice_shed_tend[0];
  });

  Kokkos::deep_copy(t_h, t_d);
  Kokkos::deep_copy(b_h, b_d);

  *log_wetgrowth_         = b_h(0);
  *qr2qi_collect_tend_    = t_h(0);
  *qc2qi_collect_tend_    = t_h(1);
  *qc_growth_rate_        = t_h(2);
  *nr_ice_shed_tend_      = t_h(3);
  *qc2qr_ice_shed_tend_   = t_h(4);
}

void get_latent_heat_f(Int its, Int ite, Int kts, Int kte, Real* v, Real* s, Real* f)
{
  using P3F        = Functions<Real, DefaultDevice>;
  using Spack      = typename P3F::Spack;
  using view_2d    = typename P3F::view_2d<Spack>;

  EKAT_REQUIRE_MSG(kte >= kts,
                     "kte must be >= kts, kts=" << kts << " kte=" << kte);

  EKAT_REQUIRE_MSG(ite >= its,
                     "ite must be >= its, its=" << its << " ite=" << ite);

  kts -= 1;
  kte -= 1;
  its -= 1;
  ite -= 1;

  Int nk = (kte - kts) + 1;
  Int nj = (ite - its) + 1;

  // Set up views
  view_2d v_d("v_d", nj, nk),
    s_d("s_d", nj, nk),
    f_d("f_d", nj, nk);

  P3F::get_latent_heat(nj, nk, v_d, s_d, f_d);

  Kokkos::Array<view_2d, 3> out_views = {v_d, s_d, f_d};
  ekat::device_to_host({v, s, f}, nj, nk, out_views, true);
}

Real subgrid_variance_scaling_f(Real relvar_, Real expon_)
{
  //The fortran version calling this function operates on scalar inputs
  //and expects scalar output. The C++ version expects relvar to be a Spack
  //and expon to be a scalar and returns a Spack.

  using P3F = Functions<Real, DefaultDevice>;
  using Spack = typename P3F::Spack;
  using Scalar = typename P3F::Scalar;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      Spack relvar(relvar_);
      Scalar expon(expon_);
      Spack out;

      out=P3F::subgrid_variance_scaling(relvar,expon);
      t_d(0) = out[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  return t_h[0];
}

void check_values_f(Real* qv, Real* temp, Int kstart, Int kend,
                    Int timestepcount, bool force_abort, Int source_ind, Real* col_loc)
{
  using P3F        = Functions<Real, DefaultDevice>;
  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using suview_1d  = typename P3F::uview_1d<Real>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kend > kstart,
                    "ktop must be larger than kstart, kstart, kend " << kend << kstart);

  kstart -= 1;
  kend -= 1;
  const auto nk = (kend - kstart) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::Array<view_1d, CheckValuesData::NUM_ARRAYS+1> cvd_d;

  ekat::host_to_device<Int,3>({qv, temp, col_loc}, {nk, nk, 3}, cvd_d);

  view_1d qv_d(cvd_d[0]), temp_d(cvd_d[1]), col_loc_d(cvd_d[2]);
  suview_1d ucol_loc_d(reinterpret_cast<Real*>(col_loc_d.data()), 3);

  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::check_values(qv_d, temp_d, kstart, kend, timestepcount, force_abort, source_ind, team,
                      ucol_loc_d);
  });
}

void calculate_incloud_mixingratios_f(Real qc_, Real qr_, Real qi_, Real qm_, Real nc_, Real nr_, Real ni_, Real bm_,
                                      Real inv_cld_frac_l_, Real inv_cld_frac_i_, Real inv_cld_frac_r_,
                                      Real* qc_incld_, Real* qr_incld_, Real* qi_incld_, Real* qm_incld_,
                                      Real* nc_incld_, Real* nr_incld_, Real* ni_incld_, Real* bm_incld_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using view_1d      = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 8);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack qc{qc_}, qr{qr_}, qi{qi_}, qm{qm_}, nc{nc_}, nr{nr_}, ni{ni_},
          bm{bm_}, inv_cld_frac_l{inv_cld_frac_l_}, inv_cld_frac_i{inv_cld_frac_i_}, inv_cld_frac_r{inv_cld_frac_r_};

    Spack qc_incld{0.}, qr_incld{0.}, qi_incld{0.}, qm_incld{0.},
          nc_incld{0.}, nr_incld{0.}, ni_incld{0.}, bm_incld{0.};

    P3F::calculate_incloud_mixingratios(qc, qr, qi, qm, nc, nr, ni, bm, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
                           qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld);

    t_d(0) = qc_incld[0];
    t_d(1) = qr_incld[0];
    t_d(2) = qi_incld[0];
    t_d(3) = qm_incld[0];
    t_d(4) = nc_incld[0];
    t_d(5) = nr_incld[0];
    t_d(6) = ni_incld[0];
    t_d(7) = bm_incld[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qc_incld_  = t_h(0);
  *qr_incld_  = t_h(1);
  *qi_incld_  = t_h(2);
  *qm_incld_  = t_h(3);
  *nc_incld_  = t_h(4);
  *nr_incld_  = t_h(5);
  *ni_incld_  = t_h(6);
  *bm_incld_  = t_h(7);
}

// Cuda implementations of std math routines are not necessarily BFB
// with the host.
template <typename ScalarT, typename DeviceT>
struct CudaWrap
{
  using Scalar = ScalarT;

  static Scalar cxx_pow(Scalar base, Scalar exp)
  {
    Scalar result;
    Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Scalar& value) {
        value = std::pow(base, exp);
    }, result);

    return result;
  }

#define cuda_wrap_single_arg(wrap_name, func_call)      \
static Scalar wrap_name(Scalar input) {                 \
  Scalar result;                                        \
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Scalar& value) { \
    value = func_call(input);                                         \
  }, result);                                                         \
  return result;                                                      \
}

  cuda_wrap_single_arg(cxx_gamma, std::tgamma)
  cuda_wrap_single_arg(cxx_sqrt, std::sqrt)
  cuda_wrap_single_arg(cxx_cbrt, std::cbrt)
  cuda_wrap_single_arg(cxx_log, std::log)
  cuda_wrap_single_arg(cxx_log10, std::log10)
  cuda_wrap_single_arg(cxx_exp, std::exp)
  cuda_wrap_single_arg(cxx_expm1, std::expm1)

#undef cuda_wrap_single_arg
};

Real cxx_pow(Real base, Real exp)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_pow(base, exp);
#else
  return std::pow(base, exp);
#endif
}

Real cxx_gamma(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_gamma(input);
#else
  return std::tgamma(input);
#endif
}

Real cxx_cbrt(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_cbrt(input);
#else
  return std::cbrt(input);
#endif
}

Real cxx_sqrt(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_sqrt(input);
#else
  return std::sqrt(input);
#endif
}

Real cxx_log(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_log(input);
#else
  return std::log(input);
#endif
}

Real cxx_log10(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_log10(input);
#else
  return std::log10(input);
#endif
}

Real cxx_exp(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_exp(input);
#else
  return std::exp(input);
#endif
}

void cloud_water_conservation_f(Real qc_, Real dt, Real* qc2qr_autoconv_tend_, Real* qc2qr_accret_tend_, Real* qc2qi_collect_tend_,
  Real* qc2qi_hetero_freeze_tend_, Real* qc2qr_ice_shed_tend_, Real* qc2qi_berg_tend_, Real* qi2qv_sublim_tend_, Real* qv2qi_vapdep_tend_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qc(qc_), qc2qr_autoconv_tend(*qc2qr_autoconv_tend_), qc2qr_accret_tend(*qc2qr_accret_tend_), qc2qi_collect_tend(*qc2qi_collect_tend_), qc2qi_hetero_freeze_tend(*qc2qi_hetero_freeze_tend_);
  Spack qc2qr_ice_shed_tend(*qc2qr_ice_shed_tend_), qc2qi_berg_tend(*qc2qi_berg_tend_), qi2qv_sublim_tend(*qi2qv_sublim_tend_), qv2qi_vapdep_tend(*qv2qi_vapdep_tend_);

  P3F::cloud_water_conservation(qc, dt, qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend);
  *qc2qr_autoconv_tend_ = qc2qr_autoconv_tend[0];
  *qc2qr_accret_tend_ = qc2qr_accret_tend[0];
  *qc2qi_collect_tend_ = qc2qi_collect_tend[0];
  *qc2qi_hetero_freeze_tend_ = qc2qi_hetero_freeze_tend[0];
  *qc2qr_ice_shed_tend_ = qc2qr_ice_shed_tend[0];
  *qc2qi_berg_tend_ = qc2qi_berg_tend[0];
  *qi2qv_sublim_tend_ = qi2qv_sublim_tend[0];
  *qv2qi_vapdep_tend_ = qv2qi_vapdep_tend[0];
}

void rain_water_conservation_f(Real qr_, Real qc2qr_autoconv_tend_, Real qc2qr_accret_tend_, Real qi2qr_melt_tend_, Real qc2qr_ice_shed_tend_,
  Real dt, Real* qr2qv_evap_tend_, Real* qr2qi_collect_tend_, Real* qr2qi_immers_freeze_tend_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qr(qr_), qc2qr_autoconv_tend(qc2qr_autoconv_tend_), qc2qr_accret_tend(qc2qr_accret_tend_), qi2qr_melt_tend(qi2qr_melt_tend_),
        qc2qr_ice_shed_tend(qc2qr_ice_shed_tend_), qr2qv_evap_tend(*qr2qv_evap_tend_);
  Spack qr2qi_collect_tend(*qr2qi_collect_tend_), qr2qi_immers_freeze_tend(*qr2qi_immers_freeze_tend_);

  P3F::rain_water_conservation(qr, qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt, qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend);
  *qr2qv_evap_tend_ = qr2qv_evap_tend[0];
  *qr2qi_collect_tend_ = qr2qi_collect_tend[0];
  *qr2qi_immers_freeze_tend_ = qr2qi_immers_freeze_tend[0];
}

void ice_water_conservation_f(Real qi_, Real qv2qi_vapdep_tend_, Real qv2qi_nucleat_tend_, Real qc2qi_berg_tend_, Real qr2qi_collect_tend_, Real qc2qi_collect_tend_,
  Real qr2qi_immers_freeze_tend_, Real qc2qi_hetero_freeze_tend_, Real dt, Real* qi2qv_sublim_tend_, Real* qi2qr_melt_tend_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qi(qi_), qv2qi_vapdep_tend(qv2qi_vapdep_tend_), qv2qi_nucleat_tend(qv2qi_nucleat_tend_), qc2qi_berg_tend(qc2qi_berg_tend_),
        qr2qi_collect_tend(qr2qi_collect_tend_), qc2qi_collect_tend(qc2qi_collect_tend_);
  Spack qr2qi_immers_freeze_tend(qr2qi_immers_freeze_tend_), qc2qi_hetero_freeze_tend(qc2qi_hetero_freeze_tend_),
        qi2qv_sublim_tend(*qi2qv_sublim_tend_), qi2qr_melt_tend(*qi2qr_melt_tend_);

  P3F::ice_water_conservation(qi, qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend,
       qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt, qi2qv_sublim_tend, qi2qr_melt_tend);
  *qi2qv_sublim_tend_ = qi2qv_sublim_tend[0];
  *qi2qr_melt_tend_ = qi2qr_melt_tend[0];
}

void p3_main_part1_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart1Data::NUM_ARRAYS> temp_d;

  ekat::host_to_device({pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci,
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld, qi_incld,
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld},
    nk, temp_d);

  view_1d
    pres_d               (temp_d[0]),
    dpres_d              (temp_d[1]),
    dz_d                 (temp_d[2]),
    nc_nuceat_tend_d     (temp_d[3]),
    exner_d              (temp_d[4]),
    inv_exner_d          (temp_d[5]),
    inv_cld_frac_l_d     (temp_d[6]),
    inv_cld_frac_i_d     (temp_d[7]),
    inv_cld_frac_r_d     (temp_d[8]),
    t_d                  (temp_d[9]),
    rho_d                (temp_d[10]),
    inv_rho_d            (temp_d[11]),
    qv_sat_l_d           (temp_d[12]),
    qv_sat_i_d           (temp_d[13]),
    qv_supersat_i_d      (temp_d[14]),
    rhofacr_d            (temp_d[15]),
    rhofaci_d            (temp_d[16]),
    acn_d                (temp_d[17]),
    qv_d                 (temp_d[18]),
    th_atm_d             (temp_d[19]),
    qc_d                 (temp_d[20]),
    nc_d                 (temp_d[21]),
    qr_d                 (temp_d[22]),
    nr_d                 (temp_d[23]),
    qi_d                 (temp_d[24]),
    ni_d                 (temp_d[25]),
    qm_d                 (temp_d[26]),
    bm_d                 (temp_d[27]),
    latent_heat_vapor_d  (temp_d[28]),
    latent_heat_sublim_d (temp_d[29]),
    latent_heat_fusion_d (temp_d[30]),
    qc_incld_d           (temp_d[31]),
    qr_incld_d           (temp_d[32]),
    qi_incld_d           (temp_d[33]),
    qm_incld_d           (temp_d[34]),
    nc_incld_d           (temp_d[35]),
    nr_incld_d           (temp_d[36]),
    ni_incld_d           (temp_d[37]),
    bm_incld_d           (temp_d[38]);

  // Call core function from kernel
  bview_1d bools_d("bools", 2);
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part1(
      team, nk, do_predict_nc, dt,
      pres_d, dpres_d, dz_d, nc_nuceat_tend_d, exner_d, inv_exner_d, inv_cld_frac_l_d, inv_cld_frac_i_d,
      inv_cld_frac_r_d, latent_heat_vapor_d, latent_heat_sublim_d, latent_heat_fusion_d,
      t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d,
      acn_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, qc_incld_d, qr_incld_d, qi_incld_d,
      qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d,
      bools_d(0), bools_d(1));
  });

  // Sync back to host
  Kokkos::Array<view_1d, 28> inout_views = {
    t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d,
    acn_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, qc_incld_d, qr_incld_d, qi_incld_d,
    qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d};

  ekat::device_to_host({T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci,
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld,
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *is_nucleat_possible  = bools_h(0);
  *is_hydromet_present = bools_h(1);
}

void p3_main_part2_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, Real dt, Real inv_dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* exner, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion, Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* is_hydromet_present)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart2Data::NUM_ARRAYS> temp_d;

  ekat::host_to_device({pres, dpres, dz, nc_nuceat_tend, exner, inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r,
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld,
        qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1,
        cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, vap_liq_exchange,
        vap_ice_exchange, liq_ice_exchange, pratot, prctot, qv_prev, t_prev
        },
    nk, temp_d);

  view_1d
    pres_d              (temp_d[0]),
    dpres_d             (temp_d[1]),
    dz_d                (temp_d[2]),
    nc_nuceat_tend_d    (temp_d[3]),
    exner_d             (temp_d[4]),
    inv_exner_d         (temp_d[5]),
    inv_cld_frac_l_d    (temp_d[6]),
    inv_cld_frac_i_d    (temp_d[7]),
    inv_cld_frac_r_d    (temp_d[8]),
    ni_activated_d      (temp_d[9]),
    inv_qc_relvar_d     (temp_d[10]),
    cld_frac_i_d        (temp_d[11]),
    cld_frac_l_d        (temp_d[12]),
    cld_frac_r_d        (temp_d[13]),
    t_d                 (temp_d[14]),
    rho_d               (temp_d[15]),
    inv_rho_d           (temp_d[16]),
    qv_sat_l_d          (temp_d[17]),
    qv_sat_i_d          (temp_d[18]),
    qv_supersat_i_d     (temp_d[19]),
    rhofacr_d           (temp_d[20]),
    rhofaci_d           (temp_d[21]),
    acn_d               (temp_d[22]),
    qv_d                (temp_d[23]),
    th_atm_d            (temp_d[24]),
    qc_d                (temp_d[25]),
    nc_d                (temp_d[26]),
    qr_d                (temp_d[27]),
    nr_d                (temp_d[28]),
    qi_d                (temp_d[29]),
    ni_d                (temp_d[30]),
    qm_d                (temp_d[31]),
    bm_d                (temp_d[32]),
    latent_heat_vapor_d (temp_d[33]),
    latent_heat_sublim_d(temp_d[34]),
    latent_heat_fusion_d(temp_d[35]),
    qc_incld_d          (temp_d[36]),
    qr_incld_d          (temp_d[37]),
    qi_incld_d          (temp_d[38]),
    qm_incld_d          (temp_d[39]),
    nc_incld_d          (temp_d[40]),
    nr_incld_d          (temp_d[41]),
    ni_incld_d          (temp_d[42]),
    bm_incld_d          (temp_d[43]),
    mu_c_d              (temp_d[44]),
    nu_d                (temp_d[45]),
    lamc_d              (temp_d[46]),
    cdist_d             (temp_d[47]),
    cdist1_d            (temp_d[48]),
    cdistr_d            (temp_d[49]),
    mu_r_d              (temp_d[50]),
    lamr_d              (temp_d[51]),
    logn0r_d            (temp_d[52]),
    qv2qi_depos_tend_d           (temp_d[53]),
    precip_total_tend_d (temp_d[54]),
    nevapr_d            (temp_d[55]),
    qr_evap_tend_d      (temp_d[56]),
    vap_liq_exchange_d  (temp_d[57]),
    vap_ice_exchange_d  (temp_d[58]),
    liq_ice_exchange_d  (temp_d[59]),
    pratot_d            (temp_d[60]),
    prctot_d            (temp_d[61]),
    qv_prev_d           (temp_d[62]),
    t_prev_d            (temp_d[63]);

  // Call core function from kernel
  const auto dnu         = P3GlobalForFortran::dnu();
  const auto ice_table_vals        = P3GlobalForFortran::ice_table_vals();
  const auto collect_table_vals     = P3GlobalForFortran::collect_table_vals();
  const auto revap_table_vals = P3GlobalForFortran::revap_table_vals();
  bview_1d bools_d("bools", 1);
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part2(
      team, nk_pack, do_predict_nc, dt, inv_dt, dnu, ice_table_vals, collect_table_vals, revap_table_vals,
      pres_d, dpres_d, dz_d, nc_nuceat_tend_d, exner_d, inv_exner_d, inv_cld_frac_l_d,
      inv_cld_frac_i_d, inv_cld_frac_r_d, ni_activated_d, inv_qc_relvar_d, cld_frac_i_d, cld_frac_l_d, cld_frac_r_d,
      qv_prev_d, t_prev_d, t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d, acn_d,
      qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d,
      latent_heat_vapor_d, latent_heat_sublim_d, latent_heat_fusion_d, qc_incld_d, qr_incld_d, qi_incld_d,
      qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d,
      mu_c_d, nu_d, lamc_d, cdist_d, cdist1_d, cdistr_d, mu_r_d, lamr_d,
      logn0r_d, qv2qi_depos_tend_d, precip_total_tend_d, nevapr_d, qr_evap_tend_d, vap_liq_exchange_d,
      vap_ice_exchange_d, liq_ice_exchange_d, pratot_d, prctot_d, bools_d(0));
  });

  // Sync back to host. Skip intent in variables.
  Kokkos::Array<view_1d, 48> inout_views = {
    t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d, acn_d,
    qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d,
    latent_heat_vapor_d, latent_heat_sublim_d, latent_heat_fusion_d, qc_incld_d, qr_incld_d, qi_incld_d, qm_incld_d,
    nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d, mu_c_d, nu_d, lamc_d,
    cdist_d, cdist1_d, cdistr_d, mu_r_d, lamr_d, logn0r_d, qv2qi_depos_tend_d, precip_total_tend_d,
    nevapr_d, qr_evap_tend_d, vap_liq_exchange_d, vap_ice_exchange_d,
    liq_ice_exchange_d, pratot_d, prctot_d
  };

  ekat::device_to_host({
      T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm, qc, nc,
      qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion, qc_incld, qr_incld,
      qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld,
      mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend,
      nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange,
      pratot, prctot},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *is_hydromet_present = bools_h(0);
}

void p3_main_part3_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc,
  Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm,
  Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* mu_c, Real* nu, Real* lamc,
  Real* mu_r, Real* lamr, Real* vap_liq_exchange, Real* ze_rain, Real* ze_ice,
  Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi,
  Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart3Data::NUM_ARRAYS> temp_d;

  ekat::host_to_device({
      exner, cld_frac_l, cld_frac_r, cld_frac_i, rho, inv_rho, rhofaci, qv, th_atm, qc,
      nc, qr, nr, qi, ni, qm, bm, latent_heat_vapor, latent_heat_sublim, mu_c, nu, lamc, mu_r,
      lamr, vap_liq_exchange, ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi,
      rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc},
    nk, temp_d);

  view_1d
    exner_d                    (temp_d[0]),
    cld_frac_l_d               (temp_d[1]),
    cld_frac_r_d               (temp_d[2]),
    cld_frac_i_d               (temp_d[3]),
    rho_d                      (temp_d[4]),
    inv_rho_d                  (temp_d[5]),
    rhofaci_d                  (temp_d[6]),
    qv_d                       (temp_d[7]),
    th_atm_d                   (temp_d[8]),
    qc_d                       (temp_d[9]),
    nc_d                       (temp_d[10]),
    qr_d                       (temp_d[11]),
    nr_d                       (temp_d[12]),
    qi_d                       (temp_d[13]),
    ni_d                       (temp_d[14]),
    qm_d                       (temp_d[15]),
    bm_d                       (temp_d[16]),
    latent_heat_vapor_d        (temp_d[17]),
    latent_heat_sublim_d       (temp_d[18]),
    mu_c_d                     (temp_d[19]),
    nu_d                       (temp_d[20]),
    lamc_d                     (temp_d[21]),
    mu_r_d                     (temp_d[22]),
    lamr_d                     (temp_d[23]),
    vap_liq_exchange_d         (temp_d[24]),
    ze_rain_d                  (temp_d[25]),
    ze_ice_d                   (temp_d[26]),
    diag_vm_qi_d               (temp_d[27]),
    diag_eff_radius_qi_d          (temp_d[28]),
    diag_diam_qi_d             (temp_d[29]),
    rho_qi_d                   (temp_d[30]),
    diag_equiv_reflectivity_d  (temp_d[31]),
    diag_eff_radius_qc_d          (temp_d[32]);

  // Call core function from kernel
  const auto dnu            = P3GlobalForFortran::dnu();
  const auto ice_table_vals = P3GlobalForFortran::ice_table_vals();
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part3(team, nk_pack, dnu, ice_table_vals,
                       exner_d, cld_frac_l_d, cld_frac_r_d, cld_frac_i_d, rho_d, inv_rho_d,
                       rhofaci_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d,
                       qi_d, ni_d, qm_d, bm_d, latent_heat_vapor_d,
                       latent_heat_sublim_d, mu_c_d, nu_d, lamc_d, mu_r_d, lamr_d,
                       vap_liq_exchange_d, ze_rain_d, ze_ice_d,
                       diag_vm_qi_d, diag_eff_radius_qi_d, diag_diam_qi_d, rho_qi_d,
                       diag_equiv_reflectivity_d, diag_eff_radius_qc_d);
  });

  // Sync back to host
  Kokkos::Array<view_1d, 29> inout_views = {
    rho_d, inv_rho_d, rhofaci_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d,
    ni_d, qm_d, bm_d, latent_heat_vapor_d, latent_heat_sublim_d, mu_c_d, nu_d, lamc_d, mu_r_d,
    lamr_d, vap_liq_exchange_d, ze_rain_d, ze_ice_d, diag_vm_qi_d, diag_eff_radius_qi_d,
    diag_diam_qi_d, rho_qi_d, diag_equiv_reflectivity_d, diag_eff_radius_qc_d
  };

  ekat::device_to_host({
      rho, inv_rho, rhofaci, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm,
      latent_heat_vapor, latent_heat_sublim, mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, ze_rain, ze_ice,
      diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc
    },
    nk, inout_views);
}

void p3_main_f(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* rho_qi, bool do_predict_nc, Real* dpres, Real* exner,
  Real* qv2qi_depos_tend, Real* precip_total_tend, Real* nevapr, Real* qr_evap_tend, Real* precip_liq_flux,
  Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i, Real* mu_c, Real* lamc,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* qv_prev, Real* t_prev)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_2d    = typename P3F::view_2d<Spack>;
  using sview_1d   = typename P3F::view_1d<Real>;
  using sview_2d   = typename P3F::view_2d<Real>;

  EKAT_REQUIRE_MSG(its == 1, "its must be 1, got " << its);
  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  its  -= 1;
  ite  -= 1;
  kts  -= 1;
  kte  -= 1;

  const Int nj    = (ite - its) + 1;
  const Int nk    = (kte - kts) + 1;

  // Set up views, pretend all views are input views for the sake of initializing kokkos views
  Kokkos::Array<view_2d, P3MainData::NUM_ARRAYS> temp_d;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS> dim1_sizes;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS> dim2_sizes;
  Kokkos::Array<const Real*, P3MainData::NUM_ARRAYS> ptr_array = {
    pres, dz, nc_nuceat_tend, ni_activated, dpres, exner, cld_frac_i, cld_frac_l, cld_frac_r, inv_qc_relvar,
    qc, nc, qr, nr, qi, qm, ni, bm, qv, th_atm, qv_prev, t_prev, diag_eff_radius_qc, diag_eff_radius_qi,
    rho_qi, mu_c, lamc, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, liq_ice_exchange,
    vap_liq_exchange, vap_ice_exchange, precip_liq_flux, precip_ice_flux, precip_liq_surf, precip_ice_surf
  };

  for (size_t i = 0; i < P3MainData::NUM_ARRAYS; ++i) dim1_sizes[i] = nj;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS; ++i) dim2_sizes[i] = nk;

  //PMC - hardcoding the index for each variable is very brittle :-(.
  dim2_sizes[34] = nk+1; // precip_liq_flux
  dim2_sizes[35] = nk+1; // precip_ice_flux
  dim1_sizes[36] = 1; dim2_sizes[36] = nj; // precip_liq_surf
  dim1_sizes[37] = 1; dim2_sizes[37] = nj; // precip_ice_surf

  // Initialize outputs to avoid uninitialized read warnings in memory checkers
  for (size_t i = P3MainData::NUM_INPUT_ARRAYS; i < P3MainData::NUM_ARRAYS; ++i) {
    for (size_t j = 0; j < dim1_sizes[i]*dim2_sizes[i]; ++j) {
      const_cast<Real*>(ptr_array[i])[j] = 0;
    }
  }

  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

  int counter = 0;
  view_2d
    pres_d                 (temp_d[counter++]), //0
    dz_d                   (temp_d[counter++]),
    nc_nuceat_tend_d       (temp_d[counter++]),
    ni_activated_d         (temp_d[counter++]),
    dpres_d                (temp_d[counter++]), 
    exner_d                (temp_d[counter++]), //5
    cld_frac_i_d           (temp_d[counter++]),
    cld_frac_l_d           (temp_d[counter++]),
    cld_frac_r_d           (temp_d[counter++]),
    inv_qc_relvar_d        (temp_d[counter++]), 
    qc_d                   (temp_d[counter++]), //10
    nc_d                   (temp_d[counter++]),
    qr_d                   (temp_d[counter++]),
    nr_d                   (temp_d[counter++]),
    qi_d                   (temp_d[counter++]), 
    qm_d                   (temp_d[counter++]), //15
    ni_d                   (temp_d[counter++]),
    bm_d                   (temp_d[counter++]),
    qv_d                   (temp_d[counter++]),
    th_atm_d               (temp_d[counter++]),
    qv_prev_d              (temp_d[counter++]), //20
    t_prev_d               (temp_d[counter++]),
    diag_eff_radius_qc_d      (temp_d[counter++]),
    diag_eff_radius_qi_d      (temp_d[counter++]),
    rho_qi_d               (temp_d[counter++]),
    mu_c_d                 (temp_d[counter++]),   // 25
    lamc_d                 (temp_d[counter++]),
    qv2qi_depos_tend_d     (temp_d[counter++]),
    precip_total_tend_d    (temp_d[counter++]),
    nevapr_d               (temp_d[counter++]), 
    qr_evap_tend_d         (temp_d[counter++]), //30
    liq_ice_exchange_d     (temp_d[counter++]),
    vap_liq_exchange_d     (temp_d[counter++]),
    vap_ice_exchange_d     (temp_d[counter++]),
    precip_liq_flux_d      (temp_d[counter++]), 
    precip_ice_flux_d      (temp_d[counter++]), //35
    precip_liq_surf_temp_d (temp_d[counter++]),
    precip_ice_surf_temp_d (temp_d[counter++]); //37

  // Special cases: precip_liq_surf=1d<scalar>(ni), precip_ice_surf=1d<scalar>(ni), col_location=2d<scalar>(nj, 3)
  sview_1d precip_liq_surf_d("precip_liq_surf_d", nj), precip_ice_surf_d("precip_ice_surf_d", nj);
  sview_2d col_location_d("col_location_d", nj, 3);

  Kokkos::parallel_for(nj, KOKKOS_LAMBDA(const Int& i) {
    precip_liq_surf_d(i) = precip_liq_surf_temp_d(0, i / Spack::n)[i % Spack::n];
    precip_ice_surf_d(i) = precip_ice_surf_temp_d(0, i / Spack::n)[i % Spack::n];

    for (int j = 0; j < 3; ++j) {
      col_location_d(i, j) = i+1;
    }
  });
  
  // Pack our data into structs and ship it off to p3_main.
  P3F::P3PrognosticState prog_state{qc_d, nc_d, qr_d, nr_d, qi_d, qm_d,
                                    ni_d, bm_d, qv_d, th_atm_d};
  P3F::P3DiagnosticInputs diag_inputs{nc_nuceat_tend_d, ni_activated_d, inv_qc_relvar_d, cld_frac_i_d,
                                      cld_frac_l_d, cld_frac_r_d, pres_d, dz_d, dpres_d,
                                      exner_d, qv_prev_d, t_prev_d};
  P3F::P3DiagnosticOutputs diag_outputs{mu_c_d, lamc_d, qv2qi_depos_tend_d, precip_liq_surf_d,
                                        precip_ice_surf_d, diag_eff_radius_qc_d, diag_eff_radius_qi_d,
                                        rho_qi_d, precip_total_tend_d, nevapr_d,
                                        qr_evap_tend_d, precip_liq_flux_d, precip_ice_flux_d};
  P3F::P3Infrastructure infrastructure{dt, it, its, ite, kts, kte,
                                       do_predict_nc, col_location_d};
  P3F::P3HistoryOnly history_only{liq_ice_exchange_d, vap_liq_exchange_d,
                                  vap_ice_exchange_d};
  
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, nj, nk);
  
  Kokkos::parallel_for(nj, KOKKOS_LAMBDA(const Int& i) {
    precip_liq_surf_temp_d(0, i / Spack::n)[i % Spack::n] = precip_liq_surf_d(i);
    precip_ice_surf_temp_d(0, i / Spack::n)[i % Spack::n] = precip_ice_surf_d(i);
  });

  // Sync back to host
  Kokkos::Array<view_2d, P3MainData::NUM_ARRAYS - 12> inout_views = {
    qc_d, nc_d, qr_d, nr_d, qi_d, qm_d, ni_d, bm_d, qv_d, th_atm_d,
    diag_eff_radius_qc_d, diag_eff_radius_qi_d, rho_qi_d, mu_c_d, lamc_d, qv2qi_depos_tend_d, precip_total_tend_d,
    nevapr_d, qr_evap_tend_d, liq_ice_exchange_d, vap_liq_exchange_d,
    vap_ice_exchange_d, precip_liq_flux_d, precip_ice_flux_d, precip_liq_surf_temp_d, precip_ice_surf_temp_d
  };
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS - 12> dim1_sizes_out;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS - 12> dim2_sizes_out;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS - 12; ++i) dim1_sizes_out[i] = nj;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS - 12; ++i) dim2_sizes_out[i] = nk;

  dim2_sizes_out[22] = nk+1; // precip_liq_flux
  dim2_sizes_out[23] = nk+1; // precip_ice_flux
  dim1_sizes_out[24] = 1; dim2_sizes_out[24] = nj; // precip_liq_surf
  dim1_sizes_out[25] = 1; dim2_sizes_out[25] = nj; // precip_ice_surf
  
  ekat::device_to_host({
      qc, nc, qr, nr, qi, qm, ni, bm, qv, th_atm, diag_eff_radius_qc, diag_eff_radius_qi,
      rho_qi, mu_c, lamc, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, liq_ice_exchange,
      vap_liq_exchange, vap_ice_exchange, precip_liq_flux, precip_ice_flux, precip_liq_surf, precip_ice_surf
    },
    dim1_sizes_out, dim2_sizes_out, inout_views, true);
}

} // namespace p3
} // namespace scream
