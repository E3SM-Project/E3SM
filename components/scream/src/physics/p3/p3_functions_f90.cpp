#include "p3_functions_f90.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
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
  Real* T_atm, Real* inv_exner, Real* latent_heat_fusion,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm);

void get_time_space_phys_variables_c(Real T_atm, Real pres, Real rho, Real latent_heat_vapor, Real latent_heat_sublim, Real qv_sat_l, Real qv_sat_i,
  Real* mu, Real* dv, Real* sc, Real* dqsdt, Real* dqsidt, Real* ab, Real* abi, Real* kap, Real* eii);

void  update_prognostic_ice_c(
  Real qc2qi_hetero_freeze_tend, Real qc2qi_collect_tend, Real qc2qr_ice_shed_tend,  Real nc_collect_tend,  Real nc2ni_immers_freeze_tend, Real ncshdc,
  Real qr2qi_collect_tend,  Real nr_collect_tend, Real qr2qi_immers_freeze_tend, Real nr2ni_immers_freeze_tend, Real nr_ice_shed_tend,
  Real qi2qr_melt_tend, Real ni2nr_melt_tend, Real qi2qv_sublim_tend, Real qv2qi_vapdep_tend, Real qv2qi_nucleat_tend, Real ni_nucleat_tend,
  Real ni_selfcollect_tend, Real ni_sublim_tend, Real qc2qi_berg_tend, Real inv_exner, Real latent_heat_sublim, Real latent_heat_fusion,
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
  Real nc_selfcollect_tend, Real  qr2qv_evap_tend, Real nr_evap_tend, Real nr_selfcollect_tend , bool do_predict_nc, bool do_prescribed_CCN,
  Real inv_rho, Real inv_exner, Real latent_heat_vapor, Real dt, Real* th_atm, Real* qv,
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
                      Real qv_supersat_i, Real inv_dt, bool do_predict_nc, bool do_prescribed_CCN,
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
  bool do_predict_nc, bool do_prescribed_CCN,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* nccn_prescribed, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present);

void p3_main_part2_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, bool do_prescribed_CCN, Real dt, Real inv_dt,
  Real* pres, Real* inv_exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofaci, Real* acn,
  Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim, Real* latent_heat_fusion, Real* qc_incld,
  Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1,
  Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* is_hydromet_present);

void p3_main_part3_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* inv_exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr,
  Real* qi, Real* ni, Real* qm, Real* bm, Real* latent_heat_vapor, Real* latent_heat_sublim,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi, Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc);

void p3_main_c(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* nccn_prescribed, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* rho_qi, bool do_predict_nc, bool do_prescribed, Real* dpres, Real* inv_exner,
  Real* qv2qi_depos_tend, Real* precip_liq_flux, Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* qv_prev, Real* t_prev, Real* elapsed_s);

void ice_supersat_conservation_c(Real* qidep, Real* qinuc, Real cld_frac_i, Real qv, Real qv_sat_i, Real latent_heat_sublim, Real t_atm, Real dt, Real qi2qv_sublim_tend, Real qr2qv_evap_tend);
void nc_conservation_c(Real nc, Real nc_selfcollect_tend, Real dt, Real* nc_collect_tend, Real* nc2ni_immers_freeze_tend, Real* nc_accret_tend, Real* nc2nr_autoconv_tend);
void nr_conservation_c(Real nr, Real ni2nr_melt_tend, Real nr_ice_shed_tend, Real ncshdc, Real nc2nr_autoconv_tend, Real dt, Real nmltratio, Real* nr_collect_tend, Real* nr2ni_immers_freeze_tend, Real* nr_selfcollect_tend, Real* nr_evap_tend);
void ni_conservation_c(Real ni, Real ni_nucleat_tend, Real nr2ni_immers_freeze_tend, Real nc2ni_immers_freeze_tend, Real dt, Real* ni2nr_melt_tend, Real* ni_sublim_tend, Real* ni_selfcollect_tend);
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

void BackToCellAverageData::randomize(std::mt19937_64& engine)
{
  // Populate the struct with numbers between 0 and 1.
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  cld_frac_l = data_dist(engine);
  cld_frac_r = data_dist(engine);
  cld_frac_i = data_dist(engine);
  qc2qr_accret_tend = data_dist(engine);
  qr2qv_evap_tend = data_dist(engine);
  qc2qr_autoconv_tend = data_dist(engine);
  nc_accret_tend = data_dist(engine);
  nc_selfcollect_tend = data_dist(engine);
  nc2nr_autoconv_tend = data_dist(engine);
  nr_selfcollect_tend = data_dist(engine);
  nr_evap_tend = data_dist(engine);
  ncautr = data_dist(engine);
  qcnuc = data_dist(engine);
  nc_nuceat_tend = data_dist(engine);
  qi2qv_sublim_tend = data_dist(engine);
  nr_ice_shed_tend = data_dist(engine);
  qc2qi_hetero_freeze_tend = data_dist(engine);
  qr2qi_collect_tend = data_dist(engine);
  qc2qr_ice_shed_tend = data_dist(engine);
  qi2qr_melt_tend = data_dist(engine);
  qc2qi_collect_tend = data_dist(engine);
  qr2qi_immers_freeze_tend = data_dist(engine);
  ni2nr_melt_tend = data_dist(engine);
  nc_collect_tend = data_dist(engine);
  ncshdc = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  nr_collect_tend = data_dist(engine);
  ni_selfcollect_tend = data_dist(engine);
  qv2qi_vapdep_tend = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  ni_sublim_tend = data_dist(engine);
  qv2qi_nucleat_tend = data_dist(engine);
  ni_nucleat_tend = data_dist(engine);
  qc2qi_berg_tend = data_dist(engine);
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
  PhysicsTestData( { {(ite_ - its_) + 1, (kte_ - kts_) + 1} },
                   { {&v, &s, &f} }),
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

void CalcLiqRelaxationData::randomize(std::mt19937_64& engine)
{
  // Populate the struct's input fields with numbers between 0 and 1.
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  rho = data_dist(engine);
  f1r = data_dist(engine);
  f2r = data_dist(engine);
  dv = data_dist(engine);
  mu = data_dist(engine);
  sc = data_dist(engine);
  mu_r = data_dist(engine);
  lamr = data_dist(engine);
  cdistr = data_dist(engine);
  cdist = data_dist(engine);
  qr_incld = data_dist(engine);
  qc_incld = data_dist(engine);
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
                   d.qv_supersat_i, d.inv_dt, d.do_predict_nc, d.do_prescribed_CCN, &d.qv2qi_nucleat_tend, &d.ni_nucleat_tend);
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
  PhysicsTestData( { {(kte_-kts_)+1} },
                   { {&qv, &temp, &col_loc} }),
  kts(kts_), kte(kte_), timestepcount(timestepcount_), source_ind(source_ind_), force_abort(force_abort_)
{
  EKAT_REQUIRE_MSG(nk() >= 3 || (kte == 1 && kts == 1), "nk too small to use for col_loc");
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
                          d.ni_selfcollect_tend,  d.ni_sublim_tend, d.qc2qi_berg_tend, d.inv_exner,  d.latent_heat_sublim,   d.latent_heat_fusion,
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
			      d.nc_selfcollect_tend, d. qr2qv_evap_tend, d.nr_evap_tend, d.nr_selfcollect_tend , d.do_predict_nc, d.do_prescribed_CCN,
			      d.inv_rho, d.inv_exner, d.latent_heat_vapor, d.dt, &d.th_atm, &d.qv,
			      &d.qc, &d.nc, &d.qr, &d.nr);
}

void ice_deposition_sublimation(IceDepSublimationData& d){
  p3_init();
  ice_deposition_sublimation_c(d.qi_incld, d.ni_incld, d.T_atm, d.qv_sat_l, d.qv_sat_i, d.epsi, d.abi,
			       d.qv, &d.qv2qi_vapdep_tend, &d.qi2qv_sublim_tend, &d.ni_sublim_tend, &d.qc2qi_berg_tend);
}

CalcUpwindData::CalcUpwindData(
  Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_) :
  PhysicsTestData({ {(kte_ - kts_)+1, num_arrays_}, {(kte_ - kts_)+1} },
                  { {&vs, &qnx, &fluxes},           {&rho, &inv_rho, &inv_dz} }),
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
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&qc_incld, &rho, &inv_rho, &cld_frac_l, &acn, &inv_dz, &qc, &nc, &nc_incld, &mu_c, &lamc, &qc_tend, &nc_tend} }),
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
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&rho, &inv_rho, &rhofaci, &cld_frac_i, &inv_dz, &qi, &qi_incld, &ni, &ni_incld, &qm, &qm_incld, &bm, &bm_incld, &qi_tend, &ni_tend} }),
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
  PhysicsTestData({ {(kte_ - kts_) + 2} }, // extra real at end for precip_liq_flux, so just add 1 to all
                  { {&rho, &inv_rho, &rhofacr, &cld_frac_r, &inv_dz, &qr_incld, &qr, &nr, &nr_incld, &mu_r, &lamr, &qr_tend, &nr_tend, &precip_liq_flux} }),
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
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&T_atm, &inv_exner, &latent_heat_fusion, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &th_atm} }),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_)
{}

void homogeneous_freezing(HomogeneousFreezingData& d)
{
  p3_init();
  homogeneous_freezing_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.T_atm, d.inv_exner, d.latent_heat_fusion,
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
  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &pres, &dpres, &dz, &nc_nuceat_tend, &inv_exner, &exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion, &nccn_prescribed,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci,
    &acn, &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &qc_incld, &qr_incld, &qi_incld,
    &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_), dt(dt_)
{}

void p3_main_part1(P3MainPart1Data& d)
{
  p3_init();
  p3_main_part1_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.do_predict_nc, d.do_prescribed_CCN,
    d.dt,
    d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.nccn_prescribed, d.inv_exner, d.exner, d.inv_cld_frac_l, d.inv_cld_frac_i, d.inv_cld_frac_r, d.latent_heat_vapor,
    d.latent_heat_sublim, d.latent_heat_fusion,
    d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci,
    d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.qc_incld, d.qr_incld, d.qi_incld,
    d.qm_incld, d.nc_incld, d.nr_incld, d.ni_incld, d.bm_incld,
    &d.is_nucleat_possible, &d.is_hydromet_present);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart2Data::P3MainPart2Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &pres, &dpres, &dz, &nc_nuceat_tend, &inv_exner, &exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &ni_activated, &inv_qc_relvar, &cld_frac_i, &cld_frac_l, &cld_frac_r, &qv_prev, &t_prev,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci, &acn,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion, &qc_incld, &qr_incld,
    &qi_incld, &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld, &mu_c, &nu, &lamc, &cdist, &cdist1,
    &cdistr, &mu_r, &lamr, &logn0r, &qv2qi_depos_tend, &precip_total_tend, &nevapr, &qr_evap_tend, &vap_liq_exchange,
    &vap_ice_exchange, &liq_ice_exchange, &pratot, &prctot} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_), dt(dt_), inv_dt(1 / dt)
{}

void p3_main_part2(P3MainPart2Data& d)
{
  p3_init();
  p3_main_part2_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir, d.do_predict_nc, d.do_prescribed_CCN, d.dt, d.inv_dt,
    d.pres, d.inv_exner, d.inv_cld_frac_l, d.inv_cld_frac_i, d.inv_cld_frac_r, d.ni_activated, d.inv_qc_relvar,
    d.cld_frac_i, d.cld_frac_l, d.cld_frac_r, d.qv_prev, d.t_prev,
    d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofaci, d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni,
    d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion, d.qc_incld, d.qr_incld, d.qi_incld, d.qm_incld, d.nc_incld, d.nr_incld,
    d.ni_incld, d.bm_incld, d.mu_c, d.nu, d.lamc, d.cdist, d.cdist1, d.cdistr, d.mu_r, d.lamr, d.logn0r, d.qv2qi_depos_tend, d.precip_total_tend,
    d.nevapr, d.qr_evap_tend, d.vap_liq_exchange, d.vap_ice_exchange, d.liq_ice_exchange, d.pratot,
    d.prctot, &d.is_hydromet_present);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart3Data::P3MainPart3Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &inv_exner, &cld_frac_l, &cld_frac_r, &cld_frac_i,
    &rho, &inv_rho, &rhofaci,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim,
    &mu_c, &nu, &lamc, &mu_r,
    &lamr, &vap_liq_exchange,
    &ze_rain, &ze_ice, &diag_vm_qi, &diag_eff_radius_qi, &diag_diam_qi, &rho_qi, &diag_equiv_reflectivity, &diag_eff_radius_qc} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_)
{}

void p3_main_part3(P3MainPart3Data& d)
{
  p3_init();
  p3_main_part3_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.inv_exner, d.cld_frac_l, d.cld_frac_r, d.cld_frac_i,
    d.rho, d.inv_rho, d.rhofaci, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim,
    d.mu_c, d.nu, d.lamc, d.mu_r, d.lamr, d.vap_liq_exchange,
    d. ze_rain, d.ze_ice, d.diag_vm_qi, d.diag_eff_radius_qi, d.diag_diam_qi, d.rho_qi, d.diag_equiv_reflectivity, d.diag_eff_radius_qc);
}

///////////////////////////////////////////////////////////////////////////////

P3MainData::P3MainData(
  Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_, bool do_prescribed_CCN_) :
  PhysicsTestData( { {(ite_ - its_) + 1, (kte_ - kts_) + 1}, {(ite_ - its_) + 1, (kte_ - kts_) + 2} }, { {
    &pres, &dz, &nc_nuceat_tend, &nccn_prescribed, &ni_activated, &dpres, &inv_exner, &cld_frac_i, &cld_frac_l, &cld_frac_r,
    &inv_qc_relvar, &qc, &nc, &qr, &nr, &qi, &qm, &ni, &bm, &qv, &th_atm, &qv_prev, &t_prev,
    &diag_eff_radius_qc, &diag_eff_radius_qi, &rho_qi, &mu_c, &lamc, &qv2qi_depos_tend, &precip_total_tend, &nevapr,
    &qr_evap_tend, &liq_ice_exchange, &vap_liq_exchange, &vap_ice_exchange, &precip_liq_flux,
    &precip_ice_flux},
    {&precip_liq_surf, &precip_ice_surf} }), // these two are (ni, nk+1)
  its(its_), ite(ite_), kts(kts_), kte(kte_), it(it_), dt(dt_), do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_)
{}

//This is the variable ordering from micro_p3.F90
void p3_main(P3MainData& d)
{
  p3_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  p3_main_c(
    d.qc, d.nc, d.qr, d.nr, d.th_atm, d.qv, d.dt, d.qi, d.qm, d.ni,
    d.bm, d.pres, d.dz, d.nc_nuceat_tend, d.nccn_prescribed, d.ni_activated, d.inv_qc_relvar, d.it, d.precip_liq_surf,
    d.precip_ice_surf, d.its, d.ite, d.kts, d.kte, d.diag_eff_radius_qc, d.diag_eff_radius_qi,
    d.rho_qi, d.do_predict_nc, d.do_prescribed_CCN, d.dpres, d.inv_exner, d.qv2qi_depos_tend,
    d.precip_liq_flux, d.precip_ice_flux, d.cld_frac_r, d.cld_frac_l, d.cld_frac_i,
    d.liq_ice_exchange, d.vap_liq_exchange, d.vap_ice_exchange, d.qv_prev, d.t_prev, &d.elapsed_s);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void ice_supersat_conservation(IceSupersatConservationData& d)
{
  p3_init();
  ice_supersat_conservation_c(&d.qidep, &d.qinuc, d.cld_frac_i, d.qv, d.qv_sat_i, d.latent_heat_sublim, d.t_atm, d.dt, d.qi2qv_sublim_tend, d.qr2qv_evap_tend);
}

void nc_conservation(NcConservationData& d)
{
  p3_init();
  nc_conservation_c(d.nc, d.nc_selfcollect_tend, d.dt, &d.nc_collect_tend, &d.nc2ni_immers_freeze_tend, &d.nc_accret_tend, &d.nc2nr_autoconv_tend);
}

void nr_conservation(NrConservationData& d)
{
  p3_init();
  nr_conservation_c(d.nr, d.ni2nr_melt_tend, d.nr_ice_shed_tend, d.ncshdc, d.nc2nr_autoconv_tend, d.dt, d.nmltratio, &d.nr_collect_tend, &d.nr2ni_immers_freeze_tend, &d.nr_selfcollect_tend, &d.nr_evap_tend);
}

void ni_conservation(NiConservationData& d)
{
  p3_init();
  ni_conservation_c(d.ni, d.ni_nucleat_tend, d.nr2ni_immers_freeze_tend, d.nc2ni_immers_freeze_tend, d.dt, &d.ni2nr_melt_tend, &d.ni_sublim_tend, &d.ni_selfcollect_tend);
}

void IceSupersatConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  cld_frac_i         = data_dist(engine);
  qv                 = data_dist(engine);
  qv_sat_i           = data_dist(engine);
  latent_heat_sublim = data_dist(engine);
  t_atm              = data_dist(engine);
  dt                 = data_dist(engine);
  qi2qv_sublim_tend  = data_dist(engine);
  qr2qv_evap_tend    = data_dist(engine);
  qidep              = data_dist(engine);
  qinuc              = data_dist(engine);
}

void NcConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  nc                       = data_dist(engine);
  nc_selfcollect_tend      = data_dist(engine);
  dt                       = data_dist(engine);
  nc_collect_tend          = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  nc_accret_tend           = data_dist(engine);
  nc2nr_autoconv_tend      = data_dist(engine);
}

void NrConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  nr                       = data_dist(engine);
  ni2nr_melt_tend          = data_dist(engine);
  nr_ice_shed_tend         = data_dist(engine);
  ncshdc                   = data_dist(engine);
  nc2nr_autoconv_tend      = data_dist(engine);
  dt                       = data_dist(engine);
  nmltratio                = data_dist(engine);
  nr_collect_tend          = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  nr_selfcollect_tend      = data_dist(engine);
  nr_evap_tend             = data_dist(engine);
}

void NiConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  ni                       = data_dist(engine);
  ni_nucleat_tend          = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  dt                       = data_dist(engine);
  ni2nr_melt_tend          = data_dist(engine);
  ni_sublim_tend           = data_dist(engine);
  ni_selfcollect_tend      = data_dist(engine);
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

} // namespace p3
} // namespace scream
