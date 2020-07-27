#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/scream_pack_kokkos.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to micro_p3 fortran calls and vice versa
//

extern "C" {

void p3_init_a_c(Real* itab, Real* itabcol);

void find_lookuptable_indices_1a_c(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qitot, Real nitot, Real qirim, Real rhop);

void find_lookuptable_indices_1b_c(Int* dumj, Real* dum3, Real qr, Real nr);

void access_lookup_table_c(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc);

void access_lookup_table_coll_c(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc);

void back_to_cell_average_c(Real lcldm_, Real rcldm_, Real icldm_,
                            Real* qcacc_, Real* qrevp_, Real* qcaut_,
                            Real* ncacc_, Real* ncslf_, Real* ncautc_,
                            Real* nrslf_, Real* nrevp_, Real* ncautr_,
                            Real* qisub_,
                            Real* nrshdr_, Real* qcheti_, Real* qrcol_,
                            Real* qcshd_, Real* qimlt_, Real* qccol_,
                            Real* qrheti_, Real* nimlt_, Real* nccol_,
                            Real* ncshdc_, Real* ncheti_, Real* nrcol_,
                            Real* nislf_, Real* qidep_, Real* nrheti_,
                            Real* nisub_, Real* qinuc_, Real* ninuc_,
                            Real* qiberg_);

void prevent_ice_overdepletion_c(Real pres, Real t, Real qv, Real xxls,
                                 Real odt, Real* qidep, Real* qisub);

void cloud_water_conservation_c(Real qc, Real dt, Real* qcaut, Real* qcacc, Real* qccol,
  Real* qcheti, Real* qcshd, Real* qiberg, Real* qisub, Real* qidep);

void rain_water_conservation_c(Real qr, Real qcaut, Real qcacc, Real qimlt, Real qcshd,
  Real dt, Real* qrevp, Real* qrcol, Real* qrheti);

void ice_water_conservation_c(Real qitot, Real qidep, Real qinuc, Real qiberg, Real qrcol, Real qccol,
  Real qrheti, Real qcheti, Real dt, Real* qisub, Real* qimlt);

void get_cloud_dsd2_c(Real qc, Real* nc, Real* mu_c, Real rho, Real* nu, Real* lamc,
                      Real* cdist, Real* cdist1, Real lcldm);

void get_rain_dsd2_c(Real qr, Real* nr, Real* mu_r, Real* lamr, Real* cdistr, Real* logn0r, Real rcldm);

void calc_rime_density_c(Real t, Real rhofaci, Real f1pr02, Real acn,
                         Real lamc, Real mu_c, Real qc_incld, Real qccol,
                         Real* vtrmi1, Real* rhorime_c);

void cldliq_immersion_freezing_c(Real t, Real lamc, Real mu_c, Real cdist1,
                                 Real qc_incld, Real qc_relvar, Real* qcheti, Real* ncheti);

void rain_immersion_freezing_c(Real t, Real lamr, Real mu_r, Real cdistr,
                               Real qr_incld, Real* qrheti, Real* nrheti);

void droplet_self_collection_c(Real rho, Real inv_rho, Real qc_incld, Real mu_c,
                               Real nu, Real ncautc, Real* ncacc);

void cloud_rain_accretion_c(Real rho, Real inv_rho, Real qc_incld, Real nc_incld,
                            Real qr_incld, Real qc_relvar, Real* qcacc, Real* ncacc);

void cloud_water_autoconversion_c(Real rho, Real qc_incld, Real nc_incld, Real qc_relvar, Real* qcaut, Real* ncautc, Real* ncautr);

void rain_self_collection_c(Real rho, Real qr_incld, Real nr_incld, Real* nrslf);

void impose_max_total_ni_c(Real* nitot_local, Real max_total_Ni, Real inv_rho_local);

void ice_melting_c(Real rho,Real t,Real pres,Real rhofaci,Real f1pr05,Real f1pr14,Real xxlv,Real xlf,Real dv,Real sc,Real mu,Real kap,Real qv,Real qitot_incld, Real nitot_incld,Real* qimlt,Real* nimlt);

void calc_first_order_upwind_step_c(Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub, Real* rho, Real* inv_rho, Real* inv_dzq, Int num_arrays, Real** fluxes, Real** vs, Real** qnx);

void generalized_sedimentation_c(Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
                                 Real* dt_left, Real* prt_accum, Real* inv_dzq, Real* inv_rho, Real* rho,
                                 Int num_arrays, Real** vs, Real** fluxes, Real** qnx);
void cloud_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* lcldm, Real* acn, Real* inv_dzq,
  Real dt, Real odt, bool log_predictNc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* prt_liq, Real* qc_tend, Real* nc_tend);

void ice_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* icldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qitot, Real* qitot_incld, Real* nitot, Real* qirim, Real* qirim_incld, Real* birim, Real* birim_incld,
  Real* nitot_incld, Real* prt_sol, Real* qi_tend, Real* ni_tend);

void rain_sedimentation_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* rcldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* prt_liq, Real* rflx, Real* qr_tend, Real* nr_tend);

void calc_bulk_rho_rime_c(Real qi_tot, Real* qi_rim, Real* bi_rim, Real* rho_rime);

void homogeneous_freezing_c(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* t, Real* exner, Real* xlf,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* th);

void get_time_space_phys_variables_c(Real t, Real pres, Real rho, Real xxlv, Real xxls, Real qvs, Real qvi,
  Real* mu, Real* dv, Real* sc, Real* dqsdt, Real* dqsidt, Real* ab, Real* abi, Real* kap, Real* eii);

void  update_prognostic_ice_c(
  Real qcheti, Real qccol, Real qcshd,  Real nccol,  Real ncheti, Real ncshdc,
  Real qrcol,  Real nrcol, Real qrheti, Real nrheti, Real nrshdr,
  Real qimlt, Real nimlt, Real qisub, Real qidep, Real qinuc, Real ninuc,
  Real nislf, Real nisub, Real qiberg, Real exner, Real xxls, Real xlf,
  bool log_predictNc, bool log_wetgrowth, Real dt, Real nmltratio,
  Real rhorime_c, Real* th, Real* qv, Real* qitot, Real* nitot, Real* qirim,
  Real* birim, Real* qc, Real* nc, Real* qr, Real* nr);

void evaporate_sublimate_precip_c(Real qr_incld, Real qc_incld, Real nr_incld, Real qitot_incld, Real lcldm,
  Real rcldm, Real qvs, Real ab, Real epsr, Real qv, Real* qrevp, Real* nrevp);

void update_prognostic_liquid_c(
  Real qcacc, Real ncacc, Real qcaut, Real ncautc, Real ncautr,
  Real ncslf, Real  qrevp, Real nrevp, Real nrslf , bool log_predictNc,
  Real inv_rho, Real exner, Real xxlv, Real dt, Real* th, Real* qv,
  Real* qc, Real* nc, Real* qr, Real* nr);

void ice_deposition_sublimation_c(
  Real qitot_incld, Real nitot_incld, Real t, Real qvs, Real qvi, Real epsi,
  Real abi, Real qv, Real* qidep, Real* qisub, Real* nisub, Real* qiberg);

void compute_rain_fall_velocity_c(Real qr_incld, Real rcldm, Real rhofacr,
                                  Real* nr_incld, Real* mu_r, Real* lamr, Real* V_qr, Real* V_nr);

void ice_cldliq_collection_c(Real rho, Real temp, Real rhofaci, Real f1pr04,
                             Real qitot_incld,Real qc_incld, Real nitot_incld, Real nc_incld,
                             Real* qccol, Real* nccol, Real* qcshd, Real* ncshdc);

void ice_rain_collection_c(Real rho, Real temp, Real rhofaci, Real logn0r, Real f1pr07, Real f1pr08,
                           Real qitot_incld, Real nitot_incld, Real qr_incld, Real* qrcol, Real* nrcol);


void ice_self_collection_c(Real rho, Real rhofaci, Real f1pr03, Real eii,
                           Real qirim_incld, Real qitot_incld, Real nitot_incld, Real* nislf);

void ice_relaxation_timescale_c(Real rho, Real temp, Real rhofaci, Real f1pr05, Real f1pr14,
                                Real dv, Real mu, Real sc, Real qitot_incld, Real nitot_incld,
                                Real* epsi, Real* epsi_tot);

void calc_liq_relaxation_timescale_c(Real rho, Real f1r, Real f2r, Real dv,
                                     Real mu, Real sc, Real mu_r, Real lamr,
                                     Real cdistr, Real cdist, Real qr_incld,
                                     Real qc_incld, Real* epsr, Real* epsc);

void ice_nucleation_c(Real temp, Real inv_rho, Real nitot, Real naai,
                      Real supi, Real odt, bool log_predictNc,
                      Real* qinuc, Real* ninuc);

void ice_cldliq_wet_growth_c(Real rho, Real temp, Real pres, Real rhofaci, Real f1pr05,
                             Real f1pr14, Real xxlv, Real xlf, Real dv,
                             Real kap, Real mu, Real sc, Real qv, Real qc_incld,
                             Real qitot_incld, Real nitot_incld, Real qr_incld, bool* log_wetgrowth,
                             Real* qrcol, Real* qccol, Real* qwgrth, Real* nrshdr, Real* qcshd);

void get_latent_heat_c(Int its, Int ite, Int kts, Int kte, Real* s, Real* v, Real* f);

Real subgrid_variance_scaling_c(Real relvar, Real expon);

void check_values_c(Real* qv, Real* temp, Int kts, Int kte, Int timestepcount,
                    Int force_abort, Int source_ind, Real* col_loc);

void calculate_incloud_mixingratios_c(Real qc, Real qr, Real qitot, Real qirim, Real nc, Real nr, Real nitot, Real birim,
                                      Real inv_lcldm, Real inv_icldm, Real inv_rcldm,
                                      Real* qc_incld, Real* qr_incld, Real* qitot_incld, Real* qirim_incld,
                                      Real* nc_incld, Real* nr_incld, Real* nitot_incld, Real* birim_incld);

void p3_main_part1_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool log_predictNc,
  Real dt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* xxlv, Real* xxls, Real* xlf,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* qc_incld, Real* qr_incld, Real* qitot_incld,
  Real* qirim_incld, Real* nc_incld, Real* nr_incld, Real* nitot_incld, Real* birim_incld,
  bool* log_nucleationPossible, bool* log_hydrometeorsPresent);

void p3_main_part2_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool log_predictNc, Real dt, Real odt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* naai, Real* qc_relvar, Real* icldm, Real* lcldm, Real* rcldm,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot,
  Real* qirim, Real* birim, Real* xxlv, Real* xxls, Real* xlf, Real* qc_incld, Real* qr_incld, Real* qitot_incld, Real* qirim_incld, Real* nc_incld, Real* nr_incld,
  Real* nitot_incld, Real* birim_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* cmeiout, Real* prain,
  Real* nevapr, Real* prer_evap, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* log_hydrometeorsPresent);

void p3_main_part3_c(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* lcldm, Real* rcldm,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* xxlv, Real* xxls,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vmi, Real* diag_effi, Real* diag_di, Real* diag_rhoi, Real* diag_ze, Real* diag_effc);

void p3_main_c(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th, Real* qv, Real dt, Real* qitot, Real* qirim, Real* nitot, Real* birim,
  Real* pres, Real* dzq, Real* ncnuc, Real* naai, Real* qc_relvar, Int it, Real* prt_liq, Real* prt_sol, Int its, Int ite, Int kts, Int kte, Real* diag_ze, Real* diag_effc,
  Real* diag_effi, Real* diag_vmi, Real* diag_di, Real* diag_rhoi, bool log_predictNc,
  Real* pdel, Real* exner, Real* cmeiout, Real* prain, Real* nevapr, Real* prer_evap, Real* rflx, Real* sflx, Real* rcldm, Real* lcldm, Real* icldm,
  Real* pratot, Real* prctot, Real* mu_c, Real* lamc, Real* liq_ice_exchange, Real* vap_liq_exchange,
  Real* vap_ice_exchange);

}

namespace scream {
namespace p3 {

// helper functions
namespace {

template <size_t N, size_t M>
void gen_random_data(const std::array<std::pair<Real, Real>, N>& ranges,
                     const std::array<Real**, M>& ptrs,
                     Real* data, Int nk)
{
  // You can provide more ptrs than ranges to initialize non-input data
  static_assert(N <= M, "Require at least as many ptrs as ranges");

  Int offset = 0;
  std::default_random_engine generator;

  for (size_t i = 0; i < N; ++i) {
    std::uniform_real_distribution<Real> data_dist(ranges[i].first, ranges[i].second);
    *ptrs[i] = data + offset;
    offset += nk;
    for(Int k = 0; k < nk; ++k) {
      (*ptrs[i])[k] = data_dist(generator);
    }
  }

  for (size_t i = N; i < M; ++i) {
    *ptrs[i] = data + offset;
    offset += nk;
  }
}

}

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
  p3_init_a_c(d.itab.data(), d.itabcol.data());
}

void find_lookuptable_indices_1a(LookupIceData& d)
{
  p3_init(); // need to initialize p3 first so that tables are loaded
  find_lookuptable_indices_1a_c(&d.dumi, &d.dumjj, &d.dumii, &d.dumzz,
                                &d.dum1, &d.dum4, &d.dum5, &d.dum6,
                                d.qitot, d.nitot, d.qirim, d.rhop);
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
  lcldm = data_dist(generator);
  rcldm = data_dist(generator);
  icldm = data_dist(generator);
  qcacc = data_dist(generator);
  qrevp = data_dist(generator);
  qcaut = data_dist(generator);
  ncacc = data_dist(generator);
  ncslf = data_dist(generator);
  ncautc = data_dist(generator);
  nrslf = data_dist(generator);
  nrevp = data_dist(generator);
  ncautr = data_dist(generator);
  qcnuc = data_dist(generator);
  ncnuc = data_dist(generator);
  qisub = data_dist(generator);
  nrshdr = data_dist(generator);
  qcheti = data_dist(generator);
  qrcol = data_dist(generator);
  qcshd = data_dist(generator);
  qimlt = data_dist(generator);
  qccol = data_dist(generator);
  qrheti = data_dist(generator);
  nimlt = data_dist(generator);
  nccol = data_dist(generator);
  ncshdc = data_dist(generator);
  ncheti = data_dist(generator);
  nrcol = data_dist(generator);
  nislf = data_dist(generator);
  qidep = data_dist(generator);
  nrheti = data_dist(generator);
  nisub = data_dist(generator);
  qinuc = data_dist(generator);
  ninuc = data_dist(generator);
  qiberg = data_dist(generator);
}

void back_to_cell_average(BackToCellAverageData& d)
{
  p3_init();
  back_to_cell_average_c(d.lcldm, d.rcldm, d.icldm, &d.qcacc, &d.qrevp,
    &d.qcaut, &d.ncacc, &d.ncslf, &d.ncautc, &d.nrslf, &d.nrevp, &d.ncautr,
    &d.qisub, &d.nrshdr, &d.qcheti, &d.qrcol, &d.qcshd,
    &d.qimlt, &d.qccol, &d.qrheti, &d.nimlt, &d.nccol, &d.ncshdc, &d.ncheti,
    &d.nrcol, &d.nislf, &d.qidep, &d.nrheti, &d.nisub, &d.qinuc, &d.ninuc,
    &d.qiberg);
}

void prevent_ice_overdepletion(PreventIceOverdepletionData& d)
{
  p3_init();
  prevent_ice_overdepletion_c(d.pres, d.t, d.qv, d.xxls, d.odt, &d.qidep,
                              &d.qisub);
}

void calc_rime_density(CalcRimeDensityData& d)
{
  p3_init();
  calc_rime_density_c(d.t, d.rhofaci, d.f1pr02, d.acn, d.lamc, d.mu_c,
                      d.qc_incld, d.qccol, &d.vtrmi1, &d.rhorime_c);
}

void cldliq_immersion_freezing(CldliqImmersionFreezingData& d)
{
  p3_init();
  cldliq_immersion_freezing_c(d.t, d.lamc, d.mu_c, d.cdist1, d.qc_incld, d.qc_relvar,
                              &d.qcheti, &d.ncheti);
}

LatentHeatData::LatentHeatData(Int kts_, Int kte_, Int its_, Int ite_) :
  its(its_), ite(ite_), kts(kts_), kte(kte_),
  m_ni((ite_ - its_) + 1), m_nk((kte_ - kts_) + 1), m_total(m_ni*m_nk),
  m_data( NUM_ARRAYS * m_total, 0.0)
{
  init_ptrs();
}

LatentHeatData::LatentHeatData(const LatentHeatData& rhs) :
  its(rhs.its), ite(rhs.ite), kts(rhs.kts), kte(rhs.kte),
  m_ni(rhs.m_ni), m_nk(rhs.m_nk), m_total(rhs.m_total),
  m_data(rhs.m_data)
{
  init_ptrs();
}

LatentHeatData& LatentHeatData::operator=(const LatentHeatData& rhs)
{
  its     = rhs.its;
  ite     = rhs.ite;
  kts     = rhs.kts;
  kte     = rhs.kte;
  m_ni    = rhs.m_ni;
  m_nk    = rhs.m_nk;
  m_total = rhs.m_total;
  m_data  = rhs.m_data; // copy

  init_ptrs();

  return *this;
}

void LatentHeatData::transpose()
{
  LatentHeatData d_trans(*this);
  util::transpose<util::TransposeDirection::f2c>(v, d_trans.v, m_ni, m_nk);
  util::transpose<util::TransposeDirection::f2c>(s, d_trans.s, m_ni, m_nk);
  util::transpose<util::TransposeDirection::f2c>(f, d_trans.f, m_ni, m_nk);

  *this = d_trans;
}

void LatentHeatData::init_ptrs()
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  std::array<Real**, NUM_ARRAYS> ptrs = {&v, &s, &f};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_total;
  }
}

void get_latent_heat(LatentHeatData& d)
{
  p3_init();
  get_latent_heat_c(d.its, d.ite, d.kts, d.kte, d.v, d.s, d.f);
  d.transpose();
}

void droplet_self_collection(DropletSelfCollectionData& d)
{
  p3_init();
  droplet_self_collection_c(d.rho, d.inv_rho, d.qc_incld, d.mu_c, d.nu, d.ncautc,
                            &d.ncslf);
}

void rain_immersion_freezing(RainImmersionFreezingData& d)
{
  p3_init();
  rain_immersion_freezing_c(d.t, d.lamr, d.mu_r, d.cdistr, d.qr_incld,
                            &d.qrheti, &d.nrheti);
}

void cloud_rain_accretion(CloudRainAccretionData& d)
{
  p3_init();
  cloud_rain_accretion_c(d.rho, d.inv_rho, d.qc_incld, d.nc_incld, d.qr_incld, d.qc_relvar,
                         &d.qcacc, &d.ncacc);
}

void cloud_water_conservation(CloudWaterConservationData& d){
  p3_init();
  cloud_water_conservation_c(d.qc, d.dt, &d.qcaut, &d.qcacc, &d.qccol, &d.qcheti,
  &d.qcshd, &d.qiberg, &d.qisub, &d.qidep);
}

void rain_water_conservation(RainWaterConservationData& d){
  p3_init();
  rain_water_conservation_c(d.qr, d.qcaut, d.qcacc, d.qimlt, d.qcshd, d.dt, &d.qrevp, &d.qrcol, &d.qrheti);
}

void ice_water_conservation(IceWaterConservationData& d){
  p3_init();
  ice_water_conservation_c(d.qitot, d.qidep, d.qinuc, d.qiberg, d.qrcol, d.qccol, d.qrheti,
    d.qcheti, d.dt, &d.qisub, &d.qimlt);
}

void cloud_water_autoconversion(CloudWaterAutoconversionData& d){
  p3_init();
  cloud_water_autoconversion_c(d.rho, d.qc_incld, d.nc_incld, d.qc_relvar,
    &d.qcaut, &d.ncautc, &d.ncautr);
}

void rain_self_collection(RainSelfCollectionData& d){
  p3_init();
  rain_self_collection_c(d.rho, d.qr_incld, d.nr_incld, &d.nrslf);
}

void impose_max_total_Ni(ImposeMaxTotalNiData& d){
  p3_init();
  impose_max_total_ni_c(&d.nitot_local, d.max_total_Ni, d.inv_rho_local);
}

void get_cloud_dsd2(GetCloudDsd2Data& d)
{
  p3_init();
  Real nc_in = d.nc_in;
  get_cloud_dsd2_c(d.qc, &nc_in, &d.mu_c, d.rho, &d.nu, &d.lamc, &d.cdist, &d.cdist1, d.lcldm);
  d.nc_out = nc_in;
}

void get_rain_dsd2(GetRainDsd2Data& d)
{
  p3_init();
  Real nr_in = d.nr_in;
  get_rain_dsd2_c(d.qr, &nr_in, &d.mu_r, &d.lamr, &d.cdistr, &d.logn0r, d.rcldm);
  d.nr_out = nr_in;
}

void ice_cldliq_collection(IceCldliqCollectionData& d)
{
  p3_init();
  ice_cldliq_collection_c(d.rho, d.temp, d.rhofaci, d.f1pr04,
                          d.qitot_incld, d.qc_incld, d.nitot_incld, d.nc_incld,
                          &d.qccol, &d.nccol, &d.qcshd, &d.ncshdc);
}

void ice_rain_collection(IceRainCollectionData& d)
{
  p3_init();
  ice_rain_collection_c(d.rho, d.temp, d.rhofaci, d.logn0r, d.f1pr07, d.f1pr08,
                        d.qitot_incld, d.nitot_incld, d.qr_incld,
                        &d.qrcol, &d.nrcol);
}

void ice_self_collection(IceSelfCollectionData& d)
{
  p3_init();
  ice_self_collection_c(d.rho, d.rhofaci, d.f1pr03, d.eii, d.qirim_incld,
                        d.qitot_incld, d.nitot_incld,
                        &d.nislf);
}

void get_time_space_phys_variables(GetTimeSpacePhysVarsData& d)
{
  p3_init();
  get_time_space_phys_variables_c(d.t, d.pres, d.rho, d.xxlv, d.xxls, d.qvs, d.qvi, &d.mu, &d.dv,
				  &d.sc, &d.dqsdt, &d.dqsidt, &d.ab, &d.abi, &d.kap, &d.eii);
}

void ice_relaxation_timescale(IceRelaxationData& d)
{
  p3_init();
  ice_relaxation_timescale_c(d.rho, d.temp, d.rhofaci, d.f1pr05, d.f1pr14,
                             d.dv, d.mu, d.sc, d.qitot_incld, d.nitot_incld,
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
  ice_nucleation_c(d.temp, d.inv_rho, d.nitot, d.naai,
                   d.supi, d.odt, d.log_predictNc,&d.qinuc, &d.ninuc);
}

void ice_cldliq_wet_growth(IceWetGrowthData& d)
{
  p3_init();

  ice_cldliq_wet_growth_c(d.rho, d.temp, d.pres, d.rhofaci, d.f1pr05,
                          d.f1pr14, d.xxlv, d.xlf, d.dv,
                          d.kap, d.mu, d.sc, d.qv, d.qc_incld,
                          d.qitot_incld, d.nitot_incld, d.qr_incld, &d.log_wetgrowth,
                          &d.qrcol, &d.qccol, &d.qwgrth, &d.nrshdr, &d.qcshd);
}

CheckValuesData::CheckValuesData(
  Int kts_, Int kte_, Int timestepcount_, Int source_ind_, bool force_abort_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), timestepcount(timestepcount_), source_ind(source_ind_), force_abort(force_abort_),
  m_nk((kte_-kts_)+1),
  m_data(NUM_ARRAYS*m_nk+3, 1.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {&qv, &temp};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  Real* data_begin = m_data.data();
  Int offset = m_nk*NUM_ARRAYS;
  col_loc = data_begin + offset;
  col_loc[1] = 2.0;
  col_loc[2] = 3.0;
}

CheckValuesData::CheckValuesData(const CheckValuesData& rhs) :
  kts(rhs.kts), kte(rhs.kte), timestepcount(rhs.timestepcount), source_ind(rhs.source_ind),
  force_abort(rhs.force_abort),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS+1] = {&qv, &temp, &col_loc};

  for (size_t i = 0; i < NUM_ARRAYS+1; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
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

  calculate_incloud_mixingratios_c(d.qc, d.qr, d.qitot, d.qirim, d.nc, d.nr, d.nitot, d.birim, d.inv_lcldm, d.inv_icldm, d.inv_rcldm,
                                   &d.qc_incld, &d.qr_incld, &d.qitot_incld, &d.qirim_incld,
                                   &d.nc_incld, &d.nr_incld, &d.nitot_incld, &d.birim_incld);

}

void update_prognostic_ice(P3UpdatePrognosticIceData& d){
  p3_init();
  update_prognostic_ice_c(d.qcheti, d.qccol, d.qcshd,  d.nccol,  d.ncheti, d.ncshdc,
                          d.qrcol,  d.nrcol, d.qrheti, d.nrheti, d.nrshdr,
                          d.qimlt,  d.nimlt, d.qisub,  d.qidep,  d.qinuc,  d.ninuc,
                          d.nislf,  d.nisub, d.qiberg, d.exner,  d.xxls,   d.xlf,
                          d.log_predictNc,  d.log_wetgrowth,    d.dt,     d.nmltratio,
                          d.rhorime_c,      &d.th,    &d.qv,    &d.qitot, &d.nitot, &d.qirim,
                          &d.birim,         &d.qc,    &d.nc,    &d.qr, &d.nr);
}

void evaporate_sublimate_precip(EvapSublimatePrecipData& d)
{
  p3_init();
  evaporate_sublimate_precip_c(d.qr_incld, d.qc_incld, d.nr_incld, d.qitot_incld,
			       d.lcldm, d.rcldm, d.qvs, d.ab, d.epsr, d.qv,
			       &d.qrevp, &d.nrevp);
}

void  update_prognostic_liquid(P3UpdatePrognosticLiqData& d){
  p3_init();
  update_prognostic_liquid_c(d.qcacc, d.ncacc, d.qcaut, d.ncautc, d.ncautr,
			      d.ncslf, d. qrevp, d.nrevp, d.nrslf , d.log_predictNc,
			      d.inv_rho, d.exner, d.xxlv, d.dt, &d.th, &d.qv,
			      &d.qc, &d.nc, &d.qr, &d.nr);
  }

void ice_deposition_sublimation(IceDepSublimationData& d){
  p3_init();
  ice_deposition_sublimation_c(d.qitot_incld, d.nitot_incld, d.t, d.qvs, d.qvi, d.epsi, d.abi,
			       d.qv, &d.qidep, &d.qisub, &d.nisub, &d.qiberg);
  }

CalcUpwindData::CalcUpwindData(
  Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_,
  std::pair<Real, Real> rho_range, std::pair<Real, Real> inv_dzq_range,
  std::pair<Real, Real> vs_range, std::pair<Real, Real> qnx_range) :
  kts(kts_), kte(kte_), kdir(kdir_), kbot(kbot_), k_qxtop(k_qxtop_), num_arrays(num_arrays_), dt_sub(dt_sub_),
  m_nk((kte_ - kts_) + 1),
  m_data( (3 + num_arrays_*3) * m_nk, 0.0),
  m_ptr_data(num_arrays_*3)
{
  Int offset = 0;

  rho     = m_data.data();
  inv_rho = rho + (offset+=m_nk);
  inv_dzq = rho + (offset+=m_nk);

  fluxes = m_ptr_data.data();
  vs     = fluxes + num_arrays;
  qnx    = vs + num_arrays;

  for (Int i = 0; i < num_arrays; ++i) {
    fluxes[i]  = rho + (offset+=m_nk);
    vs[i]      = rho + (offset+=m_nk);
    qnx[i]     = rho + (offset+=m_nk);
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real>
    rho_dist(rho_range.first, rho_range.second),
    inv_dzq_dist(inv_dzq_range.first, inv_dzq_range.second),
    vs_dist(vs_range.first, vs_range.second),
    qnx_dist(qnx_range.first, qnx_range.second);

  for (Int k = 0; k < m_nk; ++k) {
    rho[k]     = rho_dist(generator);
    inv_rho[k] = 1 / rho[k];
    inv_dzq[k] = inv_dzq_dist(generator);

    for (Int i = 0; i < num_arrays; ++i) {
      vs    [i][k] = vs_dist(generator);
      qnx   [i][k] = qnx_dist(generator);
    }
  }
}

CalcUpwindData::CalcUpwindData(const CalcUpwindData& rhs) :
  kts(rhs.kts), kte(rhs.kte), kdir(rhs.kdir), kbot(rhs.kbot), k_qxtop(rhs.k_qxtop), num_arrays(rhs.num_arrays), dt_sub(rhs.dt_sub),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data),
  m_ptr_data(rhs.m_ptr_data.size())
{
  Int offset = 0;

  rho     = m_data.data();
  inv_rho = rho + (offset+=m_nk);
  inv_dzq = rho + (offset+=m_nk);

  fluxes = m_ptr_data.data();
  vs     = fluxes + num_arrays;
  qnx    = vs + num_arrays;

  for (Int i = 0; i < num_arrays; ++i) {
    fluxes[i] = rho + (offset+=m_nk);
    vs[i]     = rho + (offset+=m_nk);
    qnx[i]    = rho + (offset+=m_nk);
  }
}

void calc_first_order_upwind_step(CalcUpwindData& d)
{
  p3_init();
  calc_first_order_upwind_step_c(d.kts, d.kte, d.kdir, d.kbot, d.k_qxtop, d.dt_sub, d.rho, d.inv_rho, d.inv_dzq, d.num_arrays, d.fluxes, d.vs, d.qnx);
}

GenSedData::GenSedData(
  Int kts_, Int kte_, Int kdir_, Int k_qxtop_, Int k_qxbot_, Int kbot_, Real Co_max_, Real dt_left_,
  Real prt_accum_, Int num_arrays_,
  std::pair<Real, Real> rho_range, std::pair<Real, Real> inv_dzq_range,
  std::pair<Real, Real> vs_range, std::pair<Real, Real> qnx_range) :
  CalcUpwindData(kts_, kte_, kdir_, kbot_, k_qxtop_, num_arrays_, 0.0, rho_range, inv_dzq_range, vs_range, qnx_range),
  Co_max(Co_max_), k_qxbot(k_qxbot_), dt_left(dt_left_), prt_accum(prt_accum_)
{ }

void generalized_sedimentation(GenSedData& d)
{
  p3_init();
  generalized_sedimentation_c(d.kts, d.kte, d.kdir, d.k_qxtop, &d.k_qxbot, d.kbot, d.Co_max,
                              &d.dt_left, &d.prt_accum, d.inv_dzq, d.inv_rho, d.rho,
                              d.num_arrays, d.vs, d.fluxes, d.qnx);
}

CloudSedData::CloudSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real odt_, bool log_predictNc_, Real prt_liq_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), odt(odt_), log_predictNc(log_predictNc_), prt_liq(prt_liq_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs =
    {&qc_incld, &rho, &inv_rho, &lcldm, &acn, &inv_dzq, &qc, &nc, &nc_incld, &mu_c, &lamc, &qc_tend, &nc_tend};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite inv_rho
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k] = 1 / rho[k];
  }
}

CloudSedData::CloudSedData(const CloudSedData& rhs) :
  kts(rhs.kts), kte(rhs.kte), ktop(rhs.ktop), kbot(rhs.kbot), kdir(rhs.kdir),
  dt(rhs.dt), odt(rhs.odt), log_predictNc(rhs.log_predictNc), prt_liq(rhs.prt_liq),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS] =
    {&qc_incld, &rho, &inv_rho, &lcldm, &acn, &inv_dzq, &qc, &nc, &nc_incld, &mu_c, &lamc, &qc_tend, &nc_tend};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void cloud_sedimentation(CloudSedData& d)
{
  p3_init();
  cloud_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                        d.qc_incld, d.rho, d.inv_rho, d.lcldm, d.acn, d.inv_dzq,
                        d.dt, d.odt, d.log_predictNc,
                        d.qc, d.nc, d.nc_incld, d.mu_c, d.lamc, &d.prt_liq, d.qc_tend, d.nc_tend);
}

IceSedData::IceSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real odt_, Real prt_sol_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), odt(odt_), prt_sol(prt_sol_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs =
    {&rho, &inv_rho, &rhofaci, &icldm, &inv_dzq, &qitot, &qitot_incld, &nitot, &nitot_incld, &qirim, &qirim_incld,
     &birim, &birim_incld, &qi_tend, &ni_tend};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite inv_rho
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k] = 1 / rho[k];
  }
}

IceSedData::IceSedData(const IceSedData& rhs) :
  kts(rhs.kts), kte(rhs.kte), ktop(rhs.ktop), kbot(rhs.kbot), kdir(rhs.kdir),
  dt(rhs.dt), odt(rhs.odt), prt_sol(rhs.prt_sol),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS] =
    {&rho, &inv_rho, &rhofaci, &icldm, &inv_dzq, &qitot, &qitot_incld, &nitot, &nitot_incld, &qirim, &qirim_incld,
     &birim, &birim_incld, &qi_tend, &ni_tend};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void ice_sedimentation(IceSedData& d)
{
  p3_init();
  ice_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                      d.rho, d.inv_rho, d.rhofaci, d.icldm, d.inv_dzq, d.dt, d.odt,
                      d.qitot, d.qitot_incld, d.nitot, d.qirim, d.qirim_incld, d.birim, d.birim_incld, d.nitot_incld,
                      &d.prt_sol, d.qi_tend, d.ni_tend);
}

RainSedData::RainSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real odt_, Real prt_liq_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), odt(odt_), prt_liq(prt_liq_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk + 1 /*extra real at end for rflx*/, 0.0)
{
  // harmless to leave last rflx value set to 0.0
  std::array<Real**, NUM_ARRAYS> ptrs =
    {&rho, &inv_rho, &rhofacr, &rcldm, &inv_dzq, &qr_incld,
     &qr, &nr, &nr_incld, &mu_r, &lamr, &qr_tend, &nr_tend, &rflx};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite inv_rho
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k] = 1 / rho[k];
  }
}

RainSedData::RainSedData(const RainSedData& rhs) :
  kts(rhs.kts), kte(rhs.kte), ktop(rhs.ktop), kbot(rhs.kbot), kdir(rhs.kdir),
  dt(rhs.dt), odt(rhs.odt), prt_liq(rhs.prt_liq),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS] =
    {&rho, &inv_rho, &rhofacr, &rcldm, &inv_dzq, &qr_incld,
     &qr, &nr, &nr_incld, &mu_r, &lamr, &qr_tend, &nr_tend, &rflx};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void rain_sedimentation(RainSedData& d)
{
  p3_init();
  rain_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                       d.qr_incld, d.rho, d.inv_rho, d.rhofacr, d.rcldm, d.inv_dzq,
                       d.dt, d.odt,
                       d.qr, d.nr, d.nr_incld, d.mu_r, d.lamr, &d.prt_liq, d.rflx, d.qr_tend, d.nr_tend);
}

void calc_bulk_rho_rime(CalcBulkRhoRimeData& d)
{
  p3_init();
  calc_bulk_rho_rime_c(d.qi_tot, &d.qi_rim, &d.bi_rim, &d.rho_rime);
}

HomogeneousFreezingData::HomogeneousFreezingData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs =
    {&t, &exner, &xlf, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &th};
  gen_random_data(ranges, ptrs, m_data.data(), m_nk);
}

HomogeneousFreezingData::HomogeneousFreezingData(const HomogeneousFreezingData& rhs) :
  kts(rhs.kts), kte(rhs.kte), ktop(rhs.ktop), kbot(rhs.kbot), kdir(rhs.kdir),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  Real** ptrs[NUM_ARRAYS] =
    {&t, &exner, &xlf, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &th};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void homogeneous_freezing(HomogeneousFreezingData& d)
{
  p3_init();
  homogeneous_freezing_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.t, d.exner, d.xlf,
                         d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.th);
}

void ice_melting(IceMeltingData& d){
  p3_init();
  ice_melting_c(d.rho,d.t,d.pres,d.rhofaci,d.f1pr05,d.f1pr14,
		d.xxlv,d.xlf,d.dv,d.sc,d.mu,d.kap,
		d.qv,d.qitot_incld,d.nitot_incld,&d.qimlt,&d.nimlt);
}

Real subgrid_variance_scaling(SubgridVarianceScalingData& d){
  p3_init();
  return subgrid_variance_scaling_c(d.relvar,d.expon);
}

void compute_rain_fall_velocity(ComputeRainFallVelocityData& d)
{
  p3_init();
  compute_rain_fall_velocity_c(d.qr_incld, d.rcldm, d.rhofacr,
                               &d.nr_incld, &d.mu_r, &d.lamr, &d.V_qr, &d.V_nr);
}

P3MainPart1Data::P3MainPart1Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool log_predictNc_, Real dt_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  log_predictNc(log_predictNc_), dt(dt_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &pdel, &dzq, &ncnuc, &exner, &inv_exner, &inv_lcldm, &inv_icldm, &inv_rcldm, &xxlv, &xxls, &xlf,
    &t, &rho, &inv_rho, &qvs, &qvi, &supi, &rhofacr, &rhofaci,
    &acn, &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &qc_incld, &qr_incld, &qitot_incld,
    &qirim_incld, &nc_incld, &nr_incld, &nitot_incld, &birim_incld};

  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite invs
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k] = 1 / rho[k];
    inv_exner[k] = 1 / exner[k];
  }
}

P3MainPart1Data::P3MainPart1Data(const P3MainPart1Data& rhs) :
  kts(rhs.kts), kte(rhs.kte), kbot(rhs.kbot), ktop(rhs.ktop), kdir(rhs.kdir),
  log_predictNc(rhs.log_predictNc), dt(rhs.dt),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &pdel, &dzq, &ncnuc, &exner, &inv_exner, &inv_lcldm, &inv_icldm, &inv_rcldm, &xxlv, &xxls, &xlf,
    &t, &rho, &inv_rho, &qvs, &qvi, &supi, &rhofacr, &rhofaci,
    &acn, &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &qc_incld, &qr_incld, &qitot_incld,
    &qirim_incld, &nc_incld, &nr_incld, &nitot_incld, &birim_incld};

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void p3_main_part1(P3MainPart1Data& d)
{
  p3_init();
  p3_main_part1_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.log_predictNc,
    d.dt,
    d.pres, d.pdel, d.dzq, d.ncnuc, d.exner, d.inv_exner, d.inv_lcldm, d.inv_icldm, d.inv_rcldm, d.xxlv, d.xxls, d.xlf,
    d.t, d.rho, d.inv_rho, d.qvs, d.qvi, d.supi, d.rhofacr, d.rhofaci,
    d.acn, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.qc_incld, d.qr_incld, d.qitot_incld,
    d.qirim_incld, d.nc_incld, d.nr_incld, d.nitot_incld, d.birim_incld,
    &d.log_nucleationPossible, &d.log_hydrometeorsPresent);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart2Data::P3MainPart2Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool log_predictNc_, Real dt_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  log_predictNc(log_predictNc_), dt(dt_), odt(1 / dt),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &pdel, &dzq, &ncnuc, &exner, &inv_exner, &inv_lcldm, &inv_icldm, &inv_rcldm, &naai, &qc_relvar, &icldm, &lcldm, &rcldm,
    &t, &rho, &inv_rho, &qvs, &qvi, &supi, &rhofacr, &rhofaci, &acn,
    &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &xxlv, &xxls, &xlf, &qc_incld, &qr_incld,
    &qitot_incld, &qirim_incld, &nc_incld, &nr_incld, &nitot_incld, &birim_incld, &mu_c, &nu, &lamc, &cdist, &cdist1,
    &cdistr, &mu_r, &lamr, &logn0r, &cmeiout, &prain, &nevapr, &prer_evap, &vap_liq_exchange,
    &vap_ice_exchange, &liq_ice_exchange, &pratot, &prctot
  };

  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite invs
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k]   = 1 / rho[k];
    inv_exner[k] = 1 / exner[k];
    inv_lcldm[k] = 1 / lcldm[k];
    inv_icldm[k] = 1 / icldm[k];
    inv_rcldm[k] = 1 / rcldm[k];
  }
}

P3MainPart2Data::P3MainPart2Data(const P3MainPart2Data& rhs) :
  kts(rhs.kts), kte(rhs.kte), kbot(rhs.kbot), ktop(rhs.ktop), kdir(rhs.kdir),
  log_predictNc(rhs.log_predictNc), dt(rhs.dt), odt(rhs.odt),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &pdel, &dzq, &ncnuc, &exner, &inv_exner, &inv_lcldm, &inv_icldm, &inv_rcldm, &naai, &qc_relvar, &icldm, &lcldm, &rcldm,
    &t, &rho, &inv_rho, &qvs, &qvi, &supi, &rhofacr, &rhofaci, &acn,
    &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &xxlv, &xxls, &xlf, &qc_incld, &qr_incld,
    &qitot_incld, &qirim_incld, &nc_incld, &nr_incld, &nitot_incld, &birim_incld, &mu_c, &nu, &lamc, &cdist, &cdist1,
    &cdistr, &mu_r, &lamr, &logn0r, &cmeiout, &prain, &nevapr, &prer_evap, &vap_liq_exchange,
    &vap_ice_exchange, &liq_ice_exchange, &pratot, &prctot
  };

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void p3_main_part2(P3MainPart2Data& d)
{
  p3_init();
  p3_main_part2_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir, d.log_predictNc, d.dt, d.odt,
    d.pres, d.pdel, d.dzq, d.ncnuc, d.exner, d.inv_exner, d.inv_lcldm, d.inv_icldm, d.inv_rcldm, d.naai, d.qc_relvar, d.icldm, d.lcldm, d.rcldm,
    d.t, d.rho, d.inv_rho, d.qvs, d.qvi, d.supi, d.rhofacr, d.rhofaci, d.acn, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot,
    d.qirim, d.birim, d.xxlv, d.xxls, d.xlf, d.qc_incld, d.qr_incld, d.qitot_incld, d.qirim_incld, d.nc_incld, d.nr_incld,
    d.nitot_incld, d.birim_incld, d.mu_c, d.nu, d.lamc, d.cdist, d.cdist1, d.cdistr, d.mu_r, d.lamr, d.logn0r, d.cmeiout, d.prain,
    d.nevapr, d.prer_evap, d.vap_liq_exchange, d.vap_ice_exchange, d.liq_ice_exchange, d.pratot,
    d.prctot, &d.log_hydrometeorsPresent);
}

///////////////////////////////////////////////////////////////////////////////

P3MainPart3Data::P3MainPart3Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  const std::array< std::pair<Real, Real>, NUM_ARRAYS >& ranges) :
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  m_nk((kte_ - kts_) + 1),
  m_data( NUM_ARRAYS * m_nk, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {
    &exner, &lcldm, &rcldm,
    &rho, &inv_rho, &rhofaci,
    &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &xxlv, &xxls,
    &mu_c, &nu, &lamc, &mu_r,
    &lamr, &vap_liq_exchange,
    &ze_rain, &ze_ice, &diag_vmi, &diag_effi, &diag_di, &diag_rhoi, &diag_ze, &diag_effc
  };

  gen_random_data(ranges, ptrs, m_data.data(), m_nk);

  // overwrite invs
  for (Int k = 0; k < m_nk; ++k) {
    inv_rho[k]   = 1 / rho[k];
  }
}

P3MainPart3Data::P3MainPart3Data(const P3MainPart3Data& rhs) :
  kts(rhs.kts), kte(rhs.kte), kbot(rhs.kbot), ktop(rhs.ktop), kdir(rhs.kdir),
  m_nk(rhs.m_nk),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  std::array<Real**, NUM_ARRAYS> ptrs = {
    &exner, &lcldm, &rcldm,
    &rho, &inv_rho, &rhofaci,
    &qv, &th, &qc, &nc, &qr, &nr, &qitot, &nitot, &qirim, &birim, &xxlv, &xxls,
    &mu_c, &nu, &lamc, &mu_r,
    &lamr, &vap_liq_exchange,
    &ze_rain, &ze_ice, &diag_vmi, &diag_effi, &diag_di, &diag_rhoi, &diag_ze, &diag_effc
  };

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nk;
  }
}

void p3_main_part3(P3MainPart3Data& d)
{
  p3_init();
  p3_main_part3_c(
    d.kts, d.kte, d.kbot, d.ktop, d.kdir,
    d.exner, d.lcldm, d.rcldm,
    d.rho, d.inv_rho, d.rhofaci, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.xxlv, d.xxls,
    d.mu_c, d.nu, d.lamc, d.mu_r, d.lamr, d.vap_liq_exchange,
    d. ze_rain, d.ze_ice, d.diag_vmi, d.diag_effi, d.diag_di, d.diag_rhoi, d.diag_ze, d.diag_effc);
}

///////////////////////////////////////////////////////////////////////////////

P3MainData::P3MainData(
  Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool log_predictNc_,
  const std::array< std::pair<Real, Real>, NUM_INPUT_ARRAYS >& ranges) :
  its(its_), ite(ite_), kts(kts_), kte(kte_), it(it_), dt(dt_), log_predictNc(log_predictNc_),
  m_ni((ite_ - its_) + 1), m_nk((kte_ - kts_) + 1),
  m_nt(m_ni * (m_nk + 1)), // overprovision since a couple data blocks are bigger than (ni, nk)
  m_data( NUM_ARRAYS * m_nt, 0.0)
{
  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &dzq, &ncnuc, &naai, &pdel, &exner, &icldm, &lcldm, &rcldm, &qc_relvar,
    &qc, &nc, &qr, &nr, &qitot, &qirim, &nitot, &birim, &qv, &th,
    &diag_ze, &diag_effc, &diag_effi, &diag_vmi, &diag_di, &diag_rhoi, &mu_c, &lamc, &cmeiout, &prain, &nevapr, &prer_evap, &pratot, &prctot, &liq_ice_exchange, &vap_liq_exchange, &vap_ice_exchange, &rflx, &sflx, &prt_liq, &prt_sol
  };

  gen_random_data(ranges, ptrs, m_data.data(), m_nt);
}

P3MainData::P3MainData(const P3MainData& rhs) :
  its(rhs.its), ite(rhs.ite), kts(rhs.kts), kte(rhs.kte), it(rhs.it), dt(rhs.dt), log_predictNc(rhs.log_predictNc),
  m_ni(rhs.m_ni), m_nk(rhs.m_nk), m_nt(rhs.m_nt),
  m_data(rhs.m_data)
{
  Int offset = 0;
  Real* data_begin = m_data.data();

  std::array<Real**, NUM_ARRAYS> ptrs = {
    &pres, &dzq, &ncnuc, &naai, &pdel, &exner, &icldm, &lcldm, &rcldm, &qc_relvar,
    &qc, &nc, &qr, &nr, &qitot, &qirim, &nitot, &birim, &qv, &th,
    &diag_ze, &diag_effc, &diag_effi, &diag_vmi, &diag_di, &diag_rhoi, &mu_c, &lamc, &cmeiout, &prain, &nevapr, &prer_evap, &pratot, &prctot, &liq_ice_exchange, &vap_liq_exchange, &vap_ice_exchange, &rflx, &sflx, &prt_liq, &prt_sol
  };

  for (size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_nt;
  }
}

void p3_main(P3MainData& d)
{
  p3_init();
  d.transpose<util::TransposeDirection::c2f>();
  p3_main_c(
    d.qc, d.nc, d.qr, d.nr, d.th, d.qv, d.dt, d.qitot, d.qirim, d.nitot, d.birim,
    d.pres, d.dzq, d.ncnuc, d.naai, d.qc_relvar, d.it, d.prt_liq, d.prt_sol, d.its, d.ite, d.kts, d.kte, d.diag_ze, d.diag_effc,
    d.diag_effi, d.diag_vmi, d.diag_di, d.diag_rhoi, d.log_predictNc,
    d.pdel, d.exner, d.cmeiout, d.prain, d.nevapr, d.prer_evap, d.rflx, d.sflx, d.rcldm, d.lcldm, d.icldm,
    d.pratot, d.prctot, d.mu_c, d.lamc, d.liq_ice_exchange, d.vap_liq_exchange,
    d.vap_ice_exchange);
  d.transpose<util::TransposeDirection::f2c>();
}

///////////////////////////////////////////////////////////////////////////////

std::shared_ptr<P3GlobalForFortran::Views> P3GlobalForFortran::s_views;

const P3GlobalForFortran::Views& P3GlobalForFortran::get()
{
  if (!P3GlobalForFortran::s_views) {
    P3GlobalForFortran::s_views = std::make_shared<Views>();
    P3F::init_kokkos_ice_lookup_tables(s_views->m_itab, s_views->m_itabcol);
    P3F::init_kokkos_tables(s_views->m_vn_table, s_views->m_vm_table,
      s_views->m_revap_table, s_views->m_mu_r_table, s_views->m_dnu);
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
                                   Real qitot_, Real nitot_, Real qirim_, Real rhop_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using TableIce = typename P3F::TableIce;

  typename P3F::Spack qitot(qitot_), nitot(nitot_), qirim(qirim_), rhop(rhop_);
  typename P3F::view_1d<TableIce> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    P3F::lookup_ice(qitot, nitot, qirim, rhop, t_d(0));
  });
  Kokkos::deep_copy(t_h, t_d);
  auto& t = t_h(0);

  // adjust for 1-based indexing
  *dumi  = t.dumi[0]  + 1;
  *dumjj = t.dumjj[0] + 1;
  *dumii = t.dumii[0] + 1;
  *dumzz = t.dumzz[0] + 1;

  *dum1 = t.dum1[0];
  *dum4 = t.dum4[0];
  *dum5 = t.dum5[0];
  *dum6 = t.dum6[0];
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
  auto& t = t_h(0);

  // adjust for 1-based indexing
  *dumj = t.dumj[0] + 1;

  *dum3 = t.dum3[0];
}

void access_lookup_table_f(Int dumjj, Int dumii, Int dumi, Int index,
                           Real dum1, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::TableIce t;

  // Adjust for 0-based indexing
  t.dumi  = dumi  - 1;
  t.dumjj = dumjj - 1;
  t.dumii = dumii - 1;

  int adjusted_index = index - 1;

  t.dum1 = dum1;
  t.dum4 = dum4;
  t.dum5 = dum5;

  auto itab = P3GlobalForFortran::itab();
  Real result;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Real& value) {
    value = P3F::apply_table_ice(adjusted_index, itab, t)[0];
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

  auto itabcol = P3GlobalForFortran::itabcol();
  Real result;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Real& value) {
    value = P3F::apply_table_coll(adjusted_index, itabcol, ti, tr)[0];
  }, result);
  *proc = result;
}

void get_cloud_dsd2_f(Real qc_, Real* nc_, Real* mu_c_, Real rho_, Real* nu_, Real* lamc_,
                      Real* cdist_, Real* cdist1_, Real lcldm_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_d", 6);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nc = *nc_;
  const auto dnu = P3GlobalForFortran::dnu();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qc(qc_), nc(local_nc), rho(rho_), lcldm(lcldm_);
    typename P3F::Spack mu_c, nu, lamc, cdist, cdist1;

    P3F::get_cloud_dsd2(qc, nc, mu_c, rho, nu, dnu, lamc, cdist, cdist1, lcldm);

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

void get_rain_dsd2_f(Real qr_, Real* nr_, Real* mu_r_, Real* lamr_, Real* cdistr_, Real* logn0r_, Real rcldm_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_d", 5);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_nr = *nr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qr(qr_), rcldm(rcldm_), nr(local_nr);
    typename P3F::Spack lamr, mu_r, cdistr, logn0r;

    P3F::get_rain_dsd2(qr, nr, mu_r, lamr, cdistr, logn0r, rcldm);

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

void get_time_space_phys_variables_f(Real t_, Real pres_, Real rho_, Real xxlv_, Real xxls_, Real qvs_, Real qvi_,
				     Real* mu_, Real* dv_, Real* sc_, Real* dqsdt_, Real* dqsidt_, Real* ab_,
				     Real* abi_, Real* kap_, Real* eii_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 9);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack t(t_), pres(pres_), rho(rho_), xxlv(xxlv_), xxls(xxls_), qvs(qvs_), qvi(qvi_);
      typename P3F::Spack mu, dv, sc, dqsdt,dqsidt, ab, abi, kap, eii;

      P3F::get_time_space_phys_variables(t, pres, rho, xxlv, xxls, qvs, qvi, mu, dv, sc, dqsdt, dqsidt,
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

void update_prognostic_ice_f( Real qcheti_, Real qccol_, Real qcshd_,  Real nccol_,  Real ncheti_, Real ncshdc_,
                              Real qrcol_,  Real nrcol_, Real qrheti_, Real nrheti_, Real nrshdr_,
                              Real qimlt_, Real nimlt_, Real qisub_, Real qidep_, Real qinuc_, Real ninuc_,
                              Real nislf_, Real nisub_, Real qiberg_, Real exner_, Real xxls_, Real xlf_,
                              bool log_predictNc_, bool log_wetgrowth_, Real dt_, Real nmltratio_,
                              Real rhorime_c_, Real* th_, Real* qv_, Real* qitot_, Real* nitot_, Real* qirim_,
                              Real* birim_, Real* qc_, Real* nc_, Real* qr_, Real* nr_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 10);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_th    = *th_;
  Real local_qv	   = *qv_;
  Real local_qc	   = *qc_;
  Real local_nc	   = *nc_;
  Real local_qr	   = *qr_;
  Real local_nr	   = *nr_;
  Real local_qitot = *qitot_;
  Real local_nitot = *nitot_;
  Real local_qirim = *qirim_;
  Real local_birim = *birim_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qcheti(qcheti_), qccol(qccol_),qcshd(qcshd_),  nccol(nccol_),
	ncheti(ncheti_),  ncshdc(ncshdc_),  qrcol(qrcol_),  nrcol(nrcol_),  qrheti(qrheti_),
	nrheti(nrheti_),  nrshdr(nrshdr_),  qimlt(qimlt_),  nimlt(nimlt_),  qisub(qisub_),
	qidep(qidep_),  qinuc(qinuc_),  ninuc(ninuc_),  nislf(nislf_),  nisub(nisub_),
	qiberg(qiberg_),  exner(exner_),  xlf(xlf_),  xxls(xxls_),
	rhorime_c(rhorime_c_);
      bool log_predictNc(log_predictNc_);
      typename P3F::Smask log_wetgrowth(log_wetgrowth_);
      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr),
	nr(local_nr), qitot(local_qitot), nitot(local_nitot), qirim(local_qirim), birim(local_birim);

      P3F::update_prognostic_ice(qcheti, qccol, qcshd, nccol, ncheti,ncshdc,
				 qrcol,   nrcol,  qrheti,  nrheti,  nrshdr,
				 qimlt,  nimlt,  qisub,  qidep,  qinuc,  ninuc,
				 nislf,  nisub,  qiberg,  exner,  xxls,  xlf,
				 log_predictNc, log_wetgrowth,  dt,  nmltratio_,
				 rhorime_c, th, qv, qitot, nitot, qirim,
				 birim, qc, nc, qr, nr);


      t_d(0) = th[0];
      t_d(1) = qv[0];
      t_d(2) = qitot[0];
      t_d(3) = nitot[0];
      t_d(4) = qirim[0];
      t_d(5) = birim[0];
      t_d(6) = qc[0];
      t_d(7) = nc[0];
      t_d(8) = qr[0];
      t_d(9) = nr[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *th_    = t_h(0);
  *qv_    = t_h(1);
  *qitot_ = t_h(2);
  *nitot_ = t_h(3);
  *qirim_ = t_h(4);
  *birim_ = t_h(5);
  *qc_    = t_h(6);
  *nc_    = t_h(7);
  *qr_    = t_h(8);
  *nr_    = t_h(9);
}

void evaporate_sublimate_precip_f(Real qr_incld_, Real qc_incld_, Real nr_incld_, Real qitot_incld_, Real lcldm_,
				  Real rcldm_, Real qvs_, Real ab_, Real epsr_, Real qv_,
				  Real* qrevp_, Real* nrevp_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qrevp = *qrevp_;
  Real local_nrevp = *nrevp_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qr_incld(qr_incld_), qc_incld(qc_incld_), nr_incld(nr_incld_), qitot_incld(qitot_incld_),
	lcldm(lcldm_), rcldm(rcldm_), qvs(qvs_), ab(ab_), epsr(epsr_), qv(qv_);

      typename P3F::Spack qrevp(local_qrevp), nrevp(local_nrevp);

      P3F::evaporate_sublimate_precip(qr_incld, qc_incld, nr_incld, qitot_incld,  lcldm, rcldm, qvs, ab,
      epsr, qv, qrevp, nrevp);

      t_d(0) = qrevp[0];
      t_d(1) = nrevp[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qrevp_ = t_h(0);
  *nrevp_ = t_h(1);
}

void update_prognostic_liquid_f(Real qcacc_, Real ncacc_, Real qcaut_, Real ncautc_, Real ncautr_,
				Real ncslf_, Real  qrevp_, Real nrevp_, Real nrslf_, bool log_predictNc_,
				Real inv_rho_, Real exner_, Real xxlv_, Real dt_, Real* th_, Real* qv_,
				Real* qc_, Real* nc_, Real* qr_, Real* nr_)

{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 6);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_th = *th_;
  Real local_qv = *qv_;
  Real local_qc = *qc_;
  Real local_nc = *nc_;
  Real local_qr = *qr_;
  Real local_nr = *nr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qcacc(qcacc_), ncacc(ncacc_), qcaut(qcaut_), ncautc(ncautc_),
	ncautr(ncautr_), ncslf(ncslf_),  qrevp( qrevp_), nrevp(nrevp_), nrslf(nrslf_), inv_rho(inv_rho_),
	exner(exner_), xxlv(xxlv_);

      bool log_predictNc(log_predictNc_);

      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr), nr(local_nr);

      P3F::update_prognostic_liquid(qcacc, ncacc, qcaut, ncautc, ncautr,
				    ncslf,  qrevp, nrevp, nrslf , log_predictNc,
				    inv_rho, exner, xxlv, dt, th, qv,
				    qc, nc, qr, nr);

      t_d(0) = th[0];
      t_d(1) = qv[0];
      t_d(2) = qc[0];
      t_d(3) = nc[0];
      t_d(4) = qr[0];
      t_d(5) = nr[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *th_    = t_h(0);
  *qv_    = t_h(1);
  *qc_    = t_h(2);
  *nc_    = t_h(3);
  *qr_    = t_h(4);
  *nr_    = t_h(5);
}

void ice_deposition_sublimation_f(Real qitot_incld_, Real nitot_incld_, Real t_, Real qvs_,
				  Real qvi_, Real epsi_, Real abi_, Real qv_,
				  Real* qidep_, Real* qisub_, Real* nisub_, Real* qiberg_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 4);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qidep  = *qidep_;
  Real local_qisub  = *qisub_;
  Real local_nisub  = *nisub_;
  Real local_qiberg = *qiberg_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack qitot_incld(qitot_incld_), nitot_incld(nitot_incld_), t(t_), qvs(qvs_), qvi(qvi_),
	epsi(epsi_), abi(abi_), qv(qv_);

      typename P3F::Spack qidep(local_qidep), qisub(local_qisub), nisub(local_nisub), qiberg(local_qiberg);

      P3F::ice_deposition_sublimation(qitot_incld, nitot_incld, t, qvs, qvi, epsi, abi, qv,
				      qidep, qisub, nisub, qiberg);

      t_d(0) = qidep[0];
      t_d(1) = qisub[0];
      t_d(2) = nisub[0];
      t_d(3) = qiberg[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qidep_  = t_h(0);
  *qisub_  = t_h(1);
  *nisub_  = t_h(2);
  *qiberg_ = t_h(3);
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
  Real* rho, Real* inv_rho, Real* inv_dzq,
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

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Setup views
  Kokkos::Array<view_1d, 3> temp_d;
  Kokkos::Array<view_1d, N> fluxes_d, vs_d, qnx_d;

  pack::host_to_device({rho, inv_rho, inv_dzq}, nk, temp_d);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dzq_d(temp_d[2]);

  pack::host_to_device(ptr_to_arr<N>((const Real**)fluxes), nk, fluxes_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)vs)    , nk, vs_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)qnx)   , nk, qnx_d);

  // Call core function from kernel
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_d[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_d[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_d[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dzq_d(inv_dzq_d);
    P3F::calc_first_order_upwind_step<N>(urho_d, uinv_rho_d, uinv_dzq_d, team, nk, kbot, k_qxtop, kdir, dt_sub, fluxes_ptr, vs_ptr, qnx_ptr);
  });

  // Sync back to host
  pack::device_to_host(ptr_to_arr<N>(fluxes), nk, fluxes_d);
  pack::device_to_host(ptr_to_arr<N>(qnx), nk, qnx_d);
}

template <int N>
void generalized_sedimentation_f_impl(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dzq, Real* inv_rho, Real* rho,
  Real** vs, Real** fluxes, Real** qnx)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using Singlep = typename pack::Pack<Real, 1>;
  using view_1d = typename P3F::view_1d<Spack>;
  using view_1ds = typename P3F::view_1d<Singlep>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using view_1d_ptr_array = typename P3F::view_1d_ptr_array<Spack, N>;
  using uview_1d = typename P3F::uview_1d<Spack>;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;
  *k_qxbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, 3> temp_d;
  Kokkos::Array<view_1d, N> fluxes_d, vs_d, qnx_d;
  Kokkos::Array<view_1ds, 1> scalar_temp;
  std::vector<Real> scalars = {*prt_accum, *dt_left, static_cast<Real>(*k_qxbot)};

  pack::host_to_device({rho, inv_rho, inv_dzq}, nk, temp_d);
  pack::host_to_device({scalars.data()}, scalars.size(), scalar_temp);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dzq_d(temp_d[2]);
  view_1ds scalars_d(scalar_temp[0]);

  pack::host_to_device(ptr_to_arr<N>((const Real**)fluxes), nk, fluxes_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)vs)    , nk, vs_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)qnx)   , nk, qnx_d);

  // Call core function from kernel
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_d[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_d[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_d[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dzq_d(inv_dzq_d);

    // Each thread needs their own copy, like we expect in the main program, or else we will hit
    // data race issues
    Real prt_accum_k = scalars_d(0)[0];
    Real dt_left_k   = scalars_d(1)[0];
    Int k_qxbot_k    = static_cast<int>(scalars_d(2)[0]);

    P3F::generalized_sedimentation<N>(urho_d, uinv_rho_d, uinv_dzq_d, team, nk, k_qxtop, k_qxbot_k, kbot, kdir, Co_max, dt_left_k, prt_accum_k, fluxes_ptr, vs_ptr, qnx_ptr);

    scalars_d(0)[0] = prt_accum_k;
    scalars_d(1)[0] = dt_left_k;
    scalars_d(2)[0] = k_qxbot_k;
  });

  // Sync back to host
  pack::device_to_host(ptr_to_arr<N>(fluxes), nk, fluxes_d);
  pack::device_to_host(ptr_to_arr<N>(qnx), nk, qnx_d);
  pack::device_to_host({scalars.data()}, scalars.size(), scalar_temp);

  // Set scalars
  *prt_accum = scalars[0];
  *dt_left   = scalars[1];
  *k_qxbot   = scalars[2] + 1;
}

void calc_first_order_upwind_step_f(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub,
  Real* rho, Real* inv_rho, Real* inv_dzq,
  Int num_arrays, Real** fluxes, Real** vs, Real** qnx)
{
  if (num_arrays == 1) {
    calc_first_order_upwind_step_f_impl<1>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, fluxes, vs, qnx);
  }
  else if (num_arrays == 2) {
    calc_first_order_upwind_step_f_impl<2>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, fluxes, vs, qnx);
  }
  else if (num_arrays == 4) {
    calc_first_order_upwind_step_f_impl<4>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, fluxes, vs, qnx);
  }
  else {
    scream_require_msg(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void generalized_sedimentation_f(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dzq, Real* inv_rho, Real* rho,
  Int num_arrays, Real** vs, Real** fluxes, Real** qnx)
{
  if (num_arrays == 1) {
    generalized_sedimentation_f_impl<1>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dzq, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 2) {
    generalized_sedimentation_f_impl<2>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dzq, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 4) {
    generalized_sedimentation_f_impl<4>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dzq, inv_rho, rho, vs, fluxes, qnx);
  }
  else {
    scream_require_msg(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void cloud_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* lcldm, Real* acn, Real* inv_dzq,
  Real dt, Real odt, bool log_predictNc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* prt_liq, Real* qc_tend, Real* nc_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Spack>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  const auto dnu = P3GlobalForFortran::dnu();

  Kokkos::Array<view_1d, CloudSedData::NUM_ARRAYS> temp_d;

  pack::host_to_device({qc_incld, rho, inv_rho, lcldm, acn, inv_dzq, qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend},
                       nk, temp_d);

  view_1d
    qc_incld_d(temp_d[0]),
    rho_d     (temp_d[1]),
    inv_rho_d (temp_d[2]),
    lcldm_d   (temp_d[3]),
    acn_d     (temp_d[4]),
    inv_dzq_d (temp_d[5]),
    qc_d      (temp_d[6]),
    nc_d      (temp_d[7]),
    nc_incld_d(temp_d[8]),
    mu_c_d    (temp_d[9]),
    lamc_d    (temp_d[10]),
    qc_tend_d (temp_d[11]),
    nc_tend_d (temp_d[12]);

  // Call core function from kernel
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_liq_k) {

    P3F::cloud_sedimentation(
      qc_incld_d, rho_d, inv_rho_d, lcldm_d, acn_d, inv_dzq_d, dnu,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, odt, log_predictNc,
      qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d,
      prt_liq_k);

  }, *prt_liq);

  // Sync back to host
  Kokkos::Array<view_1d, 7> inout_views = {qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d};
  pack::device_to_host({qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend}, nk, inout_views);
}

void ice_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* icldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qitot, Real* qitot_incld, Real* nitot, Real* qirim, Real* qirim_incld, Real* birim, Real* birim_incld,
  Real* nitot_incld, Real* prt_sol, Real* qi_tend, Real* ni_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, IceSedData::NUM_ARRAYS> temp_d;

  pack::host_to_device({rho, inv_rho, rhofaci, icldm, inv_dzq, qitot, qitot_incld, nitot, qirim, qirim_incld, birim, birim_incld, nitot_incld, qi_tend, ni_tend},
                       nk, temp_d);

  view_1d
    rho_d        (temp_d[0]),
    inv_rho_d    (temp_d[1]),
    rhofaci_d    (temp_d[2]),
    icldm_d      (temp_d[3]),
    inv_dzq_d    (temp_d[4]),
    qitot_d      (temp_d[5]),
    qitot_incld_d(temp_d[6]),
    nitot_d      (temp_d[7]),
    qirim_d      (temp_d[8]),
    qirim_incld_d(temp_d[9]),
    birim_d      (temp_d[10]),
    birim_incld_d(temp_d[11]),
    nitot_incld_d(temp_d[12]),
    qi_tend_d    (temp_d[13]),
    ni_tend_d    (temp_d[14]);

  // Call core function from kernel
  auto itab = P3GlobalForFortran::itab();
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 6, policy);
  Real my_prt_sol = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_sol_k) {

    P3F::ice_sedimentation(
      rho_d, inv_rho_d, rhofaci_d, icldm_d, inv_dzq_d,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, odt,
      qitot_d, qitot_incld_d, nitot_d, nitot_incld_d, qirim_d, qirim_incld_d, birim_d, birim_incld_d,
      qi_tend_d, ni_tend_d, itab,
      prt_sol_k);

  }, my_prt_sol);
  *prt_sol += my_prt_sol;

  // Sync back to host
  Kokkos::Array<view_1d, 10> inout_views = {qitot_d, qitot_incld_d, nitot_d, nitot_incld_d, qirim_d, qirim_incld_d,
                                            birim_d, birim_incld_d, qi_tend_d, ni_tend_d};
  pack::device_to_host({qitot, qitot_incld, nitot, nitot_incld, qirim, qirim_incld, birim, birim_incld, qi_tend, ni_tend}, nk, inout_views);
}

void rain_sedimentation_f(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* rcldm, Real* inv_dzq,
  Real dt, Real odt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* prt_liq, Real* rflx, Real* qr_tend, Real* nr_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, RainSedData::NUM_ARRAYS> temp_d;
  Kokkos::Array<size_t, RainSedData::NUM_ARRAYS> sizes;
  for (size_t i = 0; i < RainSedData::NUM_ARRAYS; ++i) sizes[i] = nk;
  sizes[RainSedData::NUM_ARRAYS - 1] = nk+1;

  pack::host_to_device({qr_incld, rho, inv_rho, rhofacr, rcldm, inv_dzq, qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, rflx},
                       sizes, temp_d);

  view_1d
    qr_incld_d   (temp_d[0]),
    rho_d        (temp_d[1]),
    inv_rho_d    (temp_d[2]),
    rhofacr_d    (temp_d[3]),
    rcldm_d      (temp_d[4]),
    inv_dzq_d    (temp_d[5]),
    qr_d         (temp_d[6]),
    nr_d         (temp_d[7]),
    nr_incld_d   (temp_d[8]),
    mu_r_d       (temp_d[9]),
    lamr_d       (temp_d[10]),
    qr_tend_d    (temp_d[11]),
    nr_tend_d    (temp_d[12]),
    rflx_d       (temp_d[13]);

  // Call core function from kernel
  auto vn_table = P3GlobalForFortran::vn_table();
  auto vm_table = P3GlobalForFortran::vm_table();
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Real my_prt_liq = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_liq_k) {

    P3F::rain_sedimentation(
      rho_d, inv_rho_d, rhofacr_d, rcldm_d, inv_dzq_d, qr_incld_d,
      team, wsm.get_workspace(team), vn_table, vm_table,
      nk, ktop, kbot, kdir, dt, odt,
      qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, rflx_d, qr_tend_d, nr_tend_d,
      prt_liq_k);

  }, my_prt_liq);
  *prt_liq += my_prt_liq;

  // Sync back to host
  Kokkos::Array<size_t, 8> sizes_out;
  for (int i = 0; i < 8; ++i) sizes_out[i] = nk;
  sizes_out[7] = nk+1;

  Kokkos::Array<view_1d, 8> inout_views = {qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, qr_tend_d, nr_tend_d, rflx_d};
  pack::device_to_host({qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, rflx}, sizes_out, inout_views);
}

void back_to_cell_average_f(Real lcldm_, Real rcldm_, Real icldm_,
                            Real* qcacc_, Real* qrevp_, Real* qcaut_,
                            Real* ncacc_, Real* ncslf_, Real* ncautc_,
                            Real* nrslf_, Real* nrevp_, Real* ncautr_,
                            Real* qisub_,
                            Real* nrshdr_, Real* qcheti_, Real* qrcol_,
                            Real* qcshd_, Real* qimlt_, Real* qccol_,
                            Real* qrheti_, Real* nimlt_, Real* nccol_,
                            Real* ncshdc_, Real* ncheti_, Real* nrcol_,
                            Real* nislf_, Real* qidep_, Real* nrheti_,
                            Real* nisub_, Real* qinuc_, Real* ninuc_,
                            Real* qiberg_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 29);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_qcacc = *qcacc_;
  Real local_qrevp = *qrevp_;
  Real local_qcaut = *qcaut_;
  Real local_ncacc = *ncacc_;
  Real local_ncslf = *ncslf_;
  Real local_ncautc = *ncautc_;
  Real local_nrslf = *nrslf_;
  Real local_nrevp = *nrevp_;
  Real local_ncautr = *ncautr_;
  Real local_qisub = *qisub_;
  Real local_nrshdr = *nrshdr_;
  Real local_qcheti = *qcheti_;
  Real local_qrcol = *qrcol_;
  Real local_qcshd = *qcshd_;
  Real local_qimlt = *qimlt_;
  Real local_qccol = *qccol_;
  Real local_qrheti = *qrheti_;
  Real local_nimlt = *nimlt_;
  Real local_nccol = *nccol_;
  Real local_ncshdc = *ncshdc_;
  Real local_ncheti = *ncheti_;
  Real local_nrcol = *nrcol_;
  Real local_nislf = *nislf_;
  Real local_qidep = *qidep_;
  Real local_nrheti = *nrheti_;
  Real local_nisub = *nisub_;
  Real local_qinuc = *qinuc_;
  Real local_ninuc = *ninuc_;
  Real local_qiberg = *qiberg_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack lcldm(lcldm_), rcldm(rcldm_), icldm(icldm_),
      qcacc(local_qcacc), qrevp(local_qrevp), qcaut(local_qcaut), ncacc(local_ncacc),
      ncslf(local_ncslf), ncautc(local_ncautc), nrslf(local_nrslf), nrevp(local_nrevp),
      ncautr(local_ncautr), qisub(local_qisub),
      nrshdr(local_nrshdr), qcheti(local_qcheti), qrcol(local_qrcol), qcshd(local_qcshd),
      qimlt(local_qimlt), qccol(local_qccol), qrheti(local_qrheti), nimlt(local_nimlt),
      nccol(local_nccol), ncshdc(local_ncshdc), ncheti(local_ncheti), nrcol(local_nrcol),
      nislf(local_nislf), qidep(local_qidep), nrheti(local_nrheti), nisub(local_nisub),
      qinuc(local_qinuc), ninuc(local_ninuc), qiberg(local_qiberg);

    P3F::back_to_cell_average(lcldm, rcldm, icldm, qcacc, qrevp, qcaut,
      ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qisub,
      nrshdr, qcheti, qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol,
      ncshdc, ncheti, nrcol, nislf, qidep, nrheti, nisub, qinuc, ninuc,
      qiberg);

    t_d(0) = qcacc[0];
    t_d(1) = qrevp[0];
    t_d(2) = qcaut[0];
    t_d(3) = ncacc[0];
    t_d(4) = ncslf[0];
    t_d(5) = ncautc[0];
    t_d(6) = nrslf[0];
    t_d(7) = nrevp[0];
    t_d(8) = ncautr[0];
    t_d(9) = qisub[0];
    t_d(10) = nrshdr[0];
    t_d(11) = qcheti[0];
    t_d(12) = qrcol[0];
    t_d(13) = qcshd[0];
    t_d(14) = qimlt[0];
    t_d(15) = qccol[0];
    t_d(16) = qrheti[0];
    t_d(17) = nimlt[0];
    t_d(18) = nccol[0];
    t_d(19) = ncshdc[0];
    t_d(20) = ncheti[0];
    t_d(21) = nrcol[0];
    t_d(22) = nislf[0];
    t_d(23) = qidep[0];
    t_d(24) = nrheti[0];
    t_d(25) = nisub[0];
    t_d(26) = qinuc[0];
    t_d(27) = ninuc[0];
    t_d(28) = qiberg[0];

  });
  Kokkos::deep_copy(t_h, t_d);

  *qcacc_ = t_h(0);
  *qrevp_ = t_h(1);
  *qcaut_ = t_h(2);
  *ncacc_ = t_h(3);
  *ncslf_ = t_h(4);
  *ncautc_ = t_h(5);
  *nrslf_ = t_h(6);
  *nrevp_ = t_h(7);
  *ncautr_ = t_h(8);
  *qisub_ = t_h(9);
  *nrshdr_ = t_h(10);
  *qcheti_ = t_h(11);
  *qrcol_ = t_h(12);
  *qcshd_ = t_h(13);
  *qimlt_ = t_h(14);
  *qccol_ = t_h(15);
  *qrheti_ = t_h(16);
  *nimlt_ = t_h(17);
  *nccol_ = t_h(18);
  *ncshdc_ = t_h(19);
  *ncheti_ = t_h(20);
  *nrcol_ = t_h(21);
  *nislf_ = t_h(22);
  *qidep_ = t_h(23);
  *nrheti_ = t_h(24);
  *nisub_ = t_h(25);
  *qinuc_ = t_h(26);
  *ninuc_ = t_h(27);
  *qiberg_ = t_h(28);
}

void prevent_ice_overdepletion_f(
  Real pres_, Real t_, Real qv_, Real xxls_, Real odt_, Real* qidep_,
  Real* qisub_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qidep = *qidep_;
  Real local_qisub = *qisub_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack pres(pres_), t(t_), qv(qv_), xxls(xxls_),
      qidep(local_qidep), qisub(local_qisub);
    P3F::prevent_ice_overdepletion(pres, t, qv, xxls, odt_, qidep, qisub);

    t_d(0) = qidep[0];
    t_d(1) = qisub[0];

  });
  Kokkos::deep_copy(t_h, t_d);

  *qidep_ = t_h(0);
  *qisub_ = t_h(1);
}

void calc_rime_density_f(
  Real t_, Real rhofaci_, Real f1pr02_, Real acn_, Real lamc_, Real mu_c_,
  Real qc_incld_, Real qccol_, Real* vtrmi1_, Real* rhorime_c_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_vtrmi1 = *vtrmi1_;
  Real local_rhorime_c = *rhorime_c_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack t(t_), rhofaci(rhofaci_), f1pr02(f1pr02_), acn(acn_),
                          lamc(lamc_), mu_c(mu_c_), qc_incld(qc_incld_),
                          qccol(qccol_), vtrmi1(local_vtrmi1),
                          rhorime_c(local_rhorime_c);
      P3F::calc_rime_density(t, rhofaci, f1pr02, acn, lamc, mu_c, qc_incld,
                             qccol, vtrmi1, rhorime_c);

      t_d(0) = vtrmi1[0];
      t_d(1) = rhorime_c[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *vtrmi1_ = t_h(0);
  *rhorime_c_ = t_h(1);
}

void cldliq_immersion_freezing_f(
  Real t_, Real lamc_, Real mu_c_, Real cdist1_, Real qc_incld_, Real qc_relvar_,
  Real* qcheti_, Real* ncheti_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qcheti = *qcheti_;
  Real local_ncheti = *ncheti_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack t(t_), lamc(lamc_), mu_c(mu_c_), cdist1(cdist1_),qc_incld(qc_incld_),
	                  qc_relvar(qc_relvar_),qcheti(local_qcheti), ncheti(local_ncheti);
      P3F::cldliq_immersion_freezing(t, lamc, mu_c, cdist1, qc_incld, qc_relvar,
                                     qcheti, ncheti);

      t_d(0) = qcheti[0];
      t_d(1) = ncheti[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qcheti_ = t_h(0);
  *ncheti_ = t_h(1);
}

void rain_immersion_freezing_f(
  Real t_, Real lamr_, Real mu_r_, Real cdistr_, Real qr_incld_,
  Real* qrheti_, Real* nrheti_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qrheti = *qrheti_;
  Real local_nrheti = *nrheti_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack t(t_), lamr(lamr_), mu_r(mu_r_),
                          cdistr(cdistr_), qr_incld(qr_incld_),
                          qrheti(local_qrheti), nrheti(local_nrheti);
      P3F::rain_immersion_freezing(t, lamr, mu_r, cdistr, qr_incld,
                                   qrheti, nrheti);

      t_d(0) = qrheti[0];
      t_d(1) = nrheti[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *qrheti_ = t_h(0);
  *nrheti_ = t_h(1);
}

void droplet_self_collection_f(
  Real rho_, Real inv_rho_, Real qc_incld_, Real mu_c_, Real nu_,
  Real ncautc_, Real* ncslf_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_ncslf = *ncslf_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), inv_rho(inv_rho_), qc_incld(qc_incld_),
                          mu_c(mu_c_), nu(nu_), ncautc(ncautc_),
                          ncslf(local_ncslf);
      P3F::droplet_self_collection(rho, inv_rho, qc_incld, mu_c, nu, ncautc,
                                   ncslf);

      t_d(0) = ncslf[0];
    });
  Kokkos::deep_copy(t_h, t_d);

  *ncslf_ = t_h(0);
}

void cloud_rain_accretion_f(
  Real rho_, Real inv_rho_, Real qc_incld_, Real nc_incld_, Real qr_incld_, Real qc_relvar_,
  Real* qcacc_, Real* ncacc_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qcacc = *qcacc_;
  Real local_ncacc = *ncacc_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), inv_rho(inv_rho_), qc_incld(qc_incld_),
                          nc_incld(nc_incld_), qr_incld(qr_incld_),
	                  qcacc(local_qcacc), ncacc(local_ncacc), qc_relvar(qc_relvar_);
      P3F::cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld, qc_relvar,
                                qcacc, ncacc);

      t_d(0) = qcacc[0];
      t_d(1) = ncacc[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qcacc_ = t_h(0);
  *ncacc_ = t_h(1);
}

void cloud_water_autoconversion_f(
     Real rho_, Real qc_incld_, Real nc_incld_, Real qc_relvar_,
     Real* qcaut_, Real* ncautc_, Real* ncautr_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 3);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qcaut = *qcaut_;
  Real local_ncautc = *ncautc_;
  Real local_ncautr = *ncautr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), qc_incld(qc_incld_), nc_incld(nc_incld_), qcaut(local_qcaut),
	ncautc(local_ncautc), ncautr(local_ncautr), qc_relvar(qc_relvar_);
      P3F::cloud_water_autoconversion(rho, qc_incld, nc_incld, qc_relvar, qcaut, ncautc, ncautr);

      t_d(0) = qcaut[0];
      t_d(1) = ncautc[0];
      t_d(2) = ncautr[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qcaut_ = t_h(0);
  *ncautc_ = t_h(1);
  *ncautr_ = t_h(2);
}

void rain_self_collection_f(Real rho_, Real qr_incld_, Real nr_incld_, Real* nrslf_){
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nrslf = *nrslf_;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), qr_incld(qr_incld_), nr_incld(nr_incld_),  nrslf(local_nrslf);
      P3F::rain_self_collection(rho, qr_incld, nr_incld, nrslf);

      t_d(0) = nrslf[0];

    });

  Kokkos::deep_copy(t_h, t_d);
  *nrslf_ = t_h(0);
}

  void ice_melting_f(Real rho_,Real t_,Real pres_,Real rhofaci_,Real f1pr05_,Real f1pr14_,Real xxlv_,Real xlf_,Real dv_,Real sc_,Real mu_,Real kap_,Real qv_,Real qitot_incld_,Real nitot_incld_,Real* qimlt_,Real* nimlt_){
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qimlt = *qimlt_;
  Real local_nimlt = *nimlt_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), t(t_), pres(pres_), rhofaci(rhofaci_),f1pr05(f1pr05_), f1pr14(f1pr14_), xxlv(xxlv_), xlf(xlf_),dv(dv_), sc(sc_), mu(mu_), kap(kap_),qv(qv_), qitot_incld(qitot_incld_), nitot_incld(nitot_incld_), qimlt(local_qimlt), nimlt(local_nimlt);
      P3F::ice_melting(rho,t,pres,rhofaci,f1pr05,f1pr14,xxlv,xlf,dv,sc,mu,kap,qv,qitot_incld,nitot_incld,qimlt,nimlt);

      t_d(0) = qimlt[0];
      t_d(1) = nimlt[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qimlt_ = t_h(0);
  *nimlt_ = t_h(1);
}

void impose_max_total_ni_f(Real* nitot_local_, Real max_total_Ni_, Real inv_rho_local_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;


  view_1d t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nitot_local = *nitot_local_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack nitot_local(local_nitot_local);
    Spack inv_rho_local(inv_rho_local_);

    P3F::impose_max_total_Ni(nitot_local, max_total_Ni_, inv_rho_local);
    t_d(0) = nitot_local[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *nitot_local_ = t_h(0);
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
  Real* t, Real* exner, Real* xlf,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* th)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, HomogeneousFreezingData::NUM_ARRAYS> temp_d;

  pack::host_to_device({t, exner, xlf, qc, nc, qr, nr, qitot, nitot, qirim, birim, th},
                       nk, temp_d);

  view_1d
    t_d    (temp_d[0]),
    exner_d(temp_d[1]),
    xlf_d  (temp_d[2]),
    qc_d   (temp_d[3]),
    nc_d   (temp_d[4]),
    qr_d   (temp_d[5]),
    nr_d   (temp_d[6]),
    qitot_d(temp_d[7]),
    nitot_d(temp_d[8]),
    qirim_d(temp_d[9]),
    birim_d(temp_d[10]),
    th_d   (temp_d[11]);

  // Call core function from kernel
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::homogeneous_freezing(
      t_d, exner_d, xlf_d,
      team,
      nk, ktop, kbot, kdir,
      qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, th_d);
  });

  // Sync back to host
  Kokkos::Array<view_1d, 9> inout_views = {qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, th_d};

  pack::device_to_host({qc, nc, qr, nr, qitot, nitot, qirim, birim, th}, nk, inout_views);
}

void compute_rain_fall_velocity_f(Real qr_incld_, Real rcldm_, Real rhofacr_,
                                  Real* nr_incld_, Real* mu_r_, Real* lamr_, Real* V_qr_, Real* V_nr_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  Real local_nr_incld = *nr_incld_;
  view_1d t_d("t_d", 5);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  const auto vn_table = P3GlobalForFortran::vn_table();
  const auto vm_table = P3GlobalForFortran::vm_table();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack qr_incld(qr_incld_), rcldm(rcldm_), rhofacr(rhofacr_), nr_incld(local_nr_incld),
      mu_r, lamr, V_qr, V_nr;

    P3F::compute_rain_fall_velocity(vn_table, vm_table,
                                    qr_incld, rcldm, rhofacr, nr_incld, mu_r, lamr, V_qr, V_nr);
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

void ice_cldliq_collection_f(Real rho_, Real temp_, Real rhofaci_, Real f1pr04_,
                             Real qitot_incld_,Real qc_incld_, Real nitot_incld_, Real nc_incld_,
                             Real* qccol_, Real* nccol_, Real* qcshd_, Real* ncshdc_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 4);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, f1pr04{f1pr04_}, qitot_incld{qitot_incld_},
          qc_incld{qc_incld_}, nitot_incld{nitot_incld_}, nc_incld{nc_incld_};
    Spack qccol{0.}, nccol{0.}, qcshd{0.}, ncshdc{0.};

    P3F::ice_cldliq_collection(rho, temp, rhofaci, f1pr04, qitot_incld, qc_incld, nitot_incld, nc_incld,
                               qccol, nccol, qcshd, ncshdc);

    t_d(0) = qccol[0];
    t_d(1) = nccol[0];
    t_d(2) = qcshd[0];
    t_d(3) = ncshdc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qccol_     = t_h(0);
  *nccol_     = t_h(1);
  *qcshd_     = t_h(2);
  *ncshdc_    = t_h(3);
}

void ice_rain_collection_f(Real rho_, Real temp_, Real rhofaci_, Real logn0r_, Real f1pr07_, Real f1pr08_,
                           Real qitot_incld_, Real nitot_incld_, Real qr_incld_, Real* qrcol_, Real* nrcol_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, logn0r{logn0r_}, f1pr07{f1pr07_}, f1pr08{f1pr08_},
          qitot_incld{qitot_incld_}, qr_incld{qr_incld_}, nitot_incld{nitot_incld_};
    Spack qrcol{0.}, nrcol{0.};

    P3F::ice_rain_collection(rho, temp, rhofaci, logn0r, f1pr07, f1pr08,
                             qitot_incld, nitot_incld, qr_incld,
                             qrcol, nrcol);

    t_d(0) = qrcol[0];
    t_d(1) = nrcol[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qrcol_     = t_h(0);
  *nrcol_     = t_h(1);
}


void ice_self_collection_f(Real rho_, Real rhofaci_, Real f1pr03_, Real eii_,
                           Real qirim_incld_, Real qitot_incld_, Real nitot_incld_, Real* nislf_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 1);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, rhofaci{rhofaci_}, f1pr03{f1pr03_}, eii{eii_}, qirim_incld{qirim_incld_},
          qitot_incld{qitot_incld_}, nitot_incld{nitot_incld_};
    Spack nislf{0.};

    P3F::ice_self_collection(rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld,
                             nislf);

    t_d(0) = nislf[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *nislf_     = t_h(0);

}


void ice_relaxation_timescale_f(Real rho_, Real temp_, Real rhofaci_, Real f1pr05_, Real f1pr14_,
                                Real dv_, Real mu_, Real sc_, Real qitot_incld_, Real nitot_incld_,
                                Real* epsi_, Real* epsi_tot_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, rhofaci{rhofaci_}, f1pr05{f1pr05_}, f1pr14{f1pr14_}, dv{dv_},
          mu{mu_}, sc{sc_}, qitot_incld{qitot_incld_}, nitot_incld{nitot_incld_};

    Spack epsi{0.0}, epsi_tot{0.0};

    P3F::ice_relaxation_timescale(rho, temp, rhofaci, f1pr05, f1pr14, dv, mu, sc, qitot_incld, nitot_incld,
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
  auto revap_table = P3GlobalForFortran::revap_table();

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, dv{dv_},
          mu{mu_}, sc{sc_}, mu_r{mu_r_}, lamr{lamr_}, cdistr{cdistr_},
          cdist{cdist_}, qr_incld{qr_incld_}, qc_incld{qc_incld_};

    Spack epsr{0.0}, epsc{0.0};

    P3F::calc_liq_relaxation_timescale(revap_table, rho, f1r_, f2r_, dv, mu, sc,
      mu_r, lamr, cdistr, cdist, qr_incld, qc_incld, epsr, epsc);

    t_d(0) = epsr[0];
    t_d(1) = epsc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *epsr_ = t_h(0);
  *epsc_ = t_h(1);
}

void ice_nucleation_f(Real temp_, Real inv_rho_, Real nitot_, Real naai_,
                      Real supi_, Real odt_, bool log_predictNc_,
                      Real* qinuc_, Real* ninuc_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using view_1d      = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack temp{temp_}, inv_rho{inv_rho_}, nitot{nitot_}, naai{naai_}, supi{supi_};
    Spack qinuc{0.0}, ninuc{0.0};

    P3F::ice_nucleation(temp, inv_rho, nitot, naai, supi, odt_, log_predictNc_,
                        qinuc, ninuc);

    t_d(0) = qinuc[0];
    t_d(1) = ninuc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qinuc_         = t_h(0);
  *ninuc_         = t_h(1);
}

void ice_cldliq_wet_growth_f(Real rho_, Real temp_, Real pres_, Real rhofaci_, Real f1pr05_,
                             Real f1pr14_, Real xxlv_, Real xlf_, Real dv_,
                             Real kap_, Real mu_, Real sc_, Real qv_, Real qc_incld_,
                             Real qitot_incld_, Real nitot_incld_, Real qr_incld_, bool* log_wetgrowth_,
                             Real* qrcol_, Real* qccol_, Real* qwgrth_, Real* nrshdr_, Real* qcshd_)
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
  Real local_qrcol = *qrcol_, local_qccol = *qccol_, local_qwgrth = *qwgrth_, local_nrshdr = *nrshdr_, local_qcshd = *qcshd_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack rho{rho_}, temp{temp_}, pres{pres_}, rhofaci{rhofaci_}, f1pr05{f1pr05_}, f1pr14{f1pr14_}, xxlv{xxlv_},
          xlf{xlf_}, dv{dv_}, kap{kap_}, mu{mu_}, sc{sc_}, qv{qv_}, qc_incld{qc_incld_}, qitot_incld{qitot_incld_},
          nitot_incld{nitot_incld_}, qr_incld{qr_incld_};

    Smask log_wetgrowth{log_wetgrowth_local};

    Spack qrcol{local_qrcol}, qccol{local_qccol}, qwgrth{local_qwgrth}, nrshdr{local_nrshdr}, qcshd{local_qcshd};

    P3F::ice_cldliq_wet_growth(rho, temp, pres, rhofaci, f1pr05, f1pr14, xxlv, xlf, dv, kap, mu, sc, qv, qc_incld,
                              qitot_incld, nitot_incld, qr_incld, log_wetgrowth,
                              qrcol, qccol, qwgrth, nrshdr, qcshd);

    b_d(0) = log_wetgrowth[0];
    t_d(0) = qrcol[0];
    t_d(1) = qccol[0];
    t_d(2) = qwgrth[0];
    t_d(3) = nrshdr[0];
    t_d(4) = qcshd[0];
  });

  Kokkos::deep_copy(t_h, t_d);
  Kokkos::deep_copy(b_h, b_d);

  *log_wetgrowth_ = b_h(0);
  *qrcol_         = t_h(0);
  *qccol_         = t_h(1);
  *qwgrth_        = t_h(2);
  *nrshdr_        = t_h(3);
  *qcshd_         = t_h(4);
}

void get_latent_heat_f(Int its, Int ite, Int kts, Int kte, Real* v, Real* s, Real* f)
{
  using P3F        = Functions<Real, DefaultDevice>;
  using Spack      = typename P3F::Spack;
  using uview_1d   = typename P3F::uview_1d<Spack>;
  using view_2d    = typename P3F::view_2d<Spack>;

  scream_require_msg(kte >= kts,
                     "kte must be >= kts, kts=" << kts << " kte=" << kte);

  scream_require_msg(ite >= its,
                     "ite must be >= its, its=" << its << " ite=" << ite);

  kts -= 1;
  kte -= 1;
  its -= 1;
  ite -= 1;

  Int nk = (kte - kts) + 1;
  Int ni = (ite - its) + 1;
  Int total = ni*nk;

  // Set up views
  view_2d v_d("v_d", ni, nk),
    s_d("s_d", ni, nk),
    f_d("f_d", ni, nk);

  P3F::get_latent_heat(ni, nk, v_d, s_d, f_d);

  Kokkos::Array<view_2d, 3> out_views = {v_d, s_d, f_d};
  pack::device_to_host({v, s, f}, ni, nk, out_views, true);
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

  scream_require_msg(kend > kstart,
                    "ktop must be larger than kstart, kstart, kend " << kend << kstart);

  kstart -= 1;
  kend -= 1;
  const unsigned long nk = (unsigned long)((kend - kstart) + 1);
  const Int nk_pack = scream::pack::npack<Spack>(nk);
  Kokkos::Array<view_1d, CheckValuesData::NUM_ARRAYS+1> cvd_d;

  pack::host_to_device({qv, temp, col_loc}, {nk, nk, 3}, cvd_d);

  view_1d qv_d(cvd_d[0]), temp_d(cvd_d[1]), col_loc_d(cvd_d[2]);
  suview_1d ucol_loc_d(reinterpret_cast<Real*>(col_loc_d.data()), 3);

  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::check_values(qv_d, temp_d, kstart, kend, timestepcount, force_abort, source_ind, team,
                      ucol_loc_d);
  });
}

void calculate_incloud_mixingratios_f(Real qc_, Real qr_, Real qitot_, Real qirim_, Real nc_, Real nr_, Real nitot_, Real birim_,
                                      Real inv_lcldm_, Real inv_icldm_, Real inv_rcldm_,
                                      Real* qc_incld_, Real* qr_incld_, Real* qitot_incld_, Real* qirim_incld_,
                                      Real* nc_incld_, Real* nr_incld_, Real* nitot_incld_, Real* birim_incld_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using view_1d      = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 8);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack qc{qc_}, qr{qr_}, qitot{qitot_}, qirim{qirim_}, nc{nc_}, nr{nr_}, nitot{nitot_},
          birim{birim_}, inv_lcldm{inv_lcldm_}, inv_icldm{inv_icldm_}, inv_rcldm{inv_rcldm_};

    Spack qc_incld{0.}, qr_incld{0.}, qitot_incld{0.}, qirim_incld{0.},
          nc_incld{0.}, nr_incld{0.}, nitot_incld{0.}, birim_incld{0.};

    P3F::calculate_incloud_mixingratios(qc, qr, qitot, qirim, nc, nr, nitot, birim, inv_lcldm, inv_icldm, inv_rcldm,
                           qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld);

    t_d(0) = qc_incld[0];
    t_d(1) = qr_incld[0];
    t_d(2) = qitot_incld[0];
    t_d(3) = qirim_incld[0];
    t_d(4) = nc_incld[0];
    t_d(5) = nr_incld[0];
    t_d(6) = nitot_incld[0];
    t_d(7) = birim_incld[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qc_incld_         = t_h(0);
  *qr_incld_         = t_h(1);
  *qitot_incld_      = t_h(2);
  *qirim_incld_      = t_h(3);
  *nc_incld_         = t_h(4);
  *nr_incld_         = t_h(5);
  *nitot_incld_      = t_h(6);
  *birim_incld_      = t_h(7);
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

void cloud_water_conservation_f(Real qc_, Real dt, Real* qcaut_, Real* qcacc_, Real* qccol_,
  Real* qcheti_, Real* qcshd_, Real* qiberg_, Real* qisub_, Real* qidep_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qc(qc_), qcaut(*qcaut_), qcacc(*qcacc_), qccol(*qccol_), qcheti(*qcheti_);
  Spack qcshd(*qcshd_), qiberg(*qiberg_), qisub(*qisub_), qidep(*qidep_);

  P3F::cloud_water_conservation(qc, dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);
  *qcaut_ = qcaut[0];
  *qcacc_ = qcacc[0];
  *qccol_ = qccol[0];
  *qcheti_ = qcheti[0];
  *qcshd_ = qcshd[0];
  *qiberg_ = qiberg[0];
  *qisub_ = qisub[0];
  *qidep_ = qidep[0];
}

void rain_water_conservation_f(Real qr_, Real qcaut_, Real qcacc_, Real qimlt_, Real qcshd_,
  Real dt, Real* qrevp_, Real* qrcol_, Real* qrheti_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qr(qr_), qcaut(qcaut_), qcacc(qcacc_), qimlt(qimlt_), qcshd(qcshd_), qrevp(*qrevp_);
  Spack qrcol(*qrcol_), qrheti(*qrheti_);

  P3F::rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, dt, qrevp, qrcol, qrheti);
  *qrevp_ = qrevp[0];
  *qrcol_ = qrcol[0];
  *qrheti_ = qrheti[0];
}

void ice_water_conservation_f(Real qitot_, Real qidep_, Real qinuc_, Real qiberg_, Real qrcol_, Real qccol_,
  Real qrheti_, Real qcheti_, Real dt, Real* qisub_, Real* qimlt_)
{
  using P3F = Functions<Real, HostDevice>;
  using Spack   = typename P3F::Spack;

  Spack qitot(qitot_), qidep(qidep_), qinuc(qinuc_), qiberg(qiberg_), qrcol(qrcol_), qccol(qccol_);
  Spack qrheti(qrheti_), qcheti(qcheti_), qisub(*qisub_), qimlt(*qimlt_);

  P3F::ice_water_conservation(qitot, qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt, qisub, qimlt);
  *qisub_ = qisub[0];
  *qimlt_ = qimlt[0];
}

void p3_main_part1_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool log_predictNc,
  Real dt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* xxlv, Real* xxls, Real* xlf,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* qc_incld, Real* qr_incld, Real* qitot_incld,
  Real* qirim_incld, Real* nc_incld, Real* nr_incld, Real* nitot_incld, Real* birim_incld,
  bool* log_nucleationPossible, bool* log_hydrometeorsPresent)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart1Data::NUM_ARRAYS> temp_d;

  pack::host_to_device({pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm,
        t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci,
        acn, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld, qitot_incld,
        qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld},
    nk, temp_d);

  view_1d
    pres_d        (temp_d[0]),
    pdel_d        (temp_d[1]),
    dzq_d         (temp_d[2]),
    ncnuc_d       (temp_d[3]),
    exner_d       (temp_d[4]),
    inv_exner_d   (temp_d[5]),
    inv_lcldm_d   (temp_d[6]),
    inv_icldm_d   (temp_d[7]),
    inv_rcldm_d   (temp_d[8]),
    t_d           (temp_d[9]),
    rho_d         (temp_d[10]),
    inv_rho_d     (temp_d[11]),
    qvs_d         (temp_d[12]),
    qvi_d         (temp_d[13]),
    supi_d        (temp_d[14]),
    rhofacr_d     (temp_d[15]),
    rhofaci_d     (temp_d[16]),
    acn_d         (temp_d[17]),
    qv_d          (temp_d[18]),
    th_d          (temp_d[19]),
    qc_d          (temp_d[20]),
    nc_d          (temp_d[21]),
    qr_d          (temp_d[22]),
    nr_d          (temp_d[23]),
    qitot_d       (temp_d[24]),
    nitot_d       (temp_d[25]),
    qirim_d       (temp_d[26]),
    birim_d       (temp_d[27]),
    xxlv_d        (temp_d[28]),
    xxls_d        (temp_d[29]),
    xlf_d         (temp_d[30]),
    qc_incld_d    (temp_d[31]),
    qr_incld_d    (temp_d[32]),
    qitot_incld_d (temp_d[33]),
    qirim_incld_d (temp_d[34]),
    nc_incld_d    (temp_d[35]),
    nr_incld_d    (temp_d[36]),
    nitot_incld_d (temp_d[37]),
    birim_incld_d (temp_d[38]);

  // Call core function from kernel
  bview_1d bools_d("bools", 2);
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part1(
      team, nk, log_predictNc, dt,
      pres_d, pdel_d, dzq_d, ncnuc_d, exner_d, inv_exner_d, inv_lcldm_d, inv_icldm_d, inv_rcldm_d, xxlv_d, xxls_d, xlf_d,
      t_d, rho_d, inv_rho_d, qvs_d, qvi_d, supi_d, rhofacr_d, rhofaci_d,
      acn_d, qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, qc_incld_d, qr_incld_d, qitot_incld_d,
      qirim_incld_d, nc_incld_d, nr_incld_d, nitot_incld_d, birim_incld_d,
      bools_d(0), bools_d(1));
  });

  // Sync back to host
  Kokkos::Array<view_1d, 28> inout_views = {
    t_d, rho_d, inv_rho_d, qvs_d, qvi_d, supi_d, rhofacr_d, rhofaci_d,
    acn_d, qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, qc_incld_d, qr_incld_d, qitot_incld_d,
    qirim_incld_d, nc_incld_d, nr_incld_d, nitot_incld_d, birim_incld_d};

  pack::device_to_host({t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci,
        acn, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, qc_incld, qr_incld, qitot_incld,
        qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *log_nucleationPossible  = bools_h(0);
  *log_hydrometeorsPresent = bools_h(1);
}

void p3_main_part2_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool log_predictNc, Real dt, Real odt,
  Real* pres, Real* pdel, Real* dzq, Real* ncnuc, Real* exner, Real* inv_exner, Real* inv_lcldm, Real* inv_icldm, Real* inv_rcldm, Real* naai, Real* qc_relvar, Real* icldm, Real* lcldm, Real* rcldm,
  Real* t, Real* rho, Real* inv_rho, Real* qvs, Real* qvi, Real* supi, Real* rhofacr, Real* rhofaci, Real* acn, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot,
  Real* qirim, Real* birim, Real* xxlv, Real* xxls, Real* xlf, Real* qc_incld, Real* qr_incld, Real* qitot_incld, Real* qirim_incld, Real* nc_incld, Real* nr_incld,
  Real* nitot_incld, Real* birim_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* cmeiout, Real* prain,
  Real* nevapr, Real* prer_evap, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange, Real* pratot,
  Real* prctot, bool* log_hydrometeorsPresent)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart2Data::NUM_ARRAYS> temp_d;

  pack::host_to_device({pres, pdel, dzq, ncnuc, exner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, naai, qc_relvar, icldm, lcldm, rcldm,
        t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn,
        qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld,
        qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1,
        cdistr, mu_r, lamr, logn0r, cmeiout, prain, nevapr, prer_evap, vap_liq_exchange,
        vap_ice_exchange, liq_ice_exchange, pratot, prctot
        },
    nk, temp_d);

  view_1d
    pres_d             (temp_d[0]),
    pdel_d             (temp_d[1]),
    dzq_d              (temp_d[2]),
    ncnuc_d            (temp_d[3]),
    exner_d            (temp_d[4]),
    inv_exner_d        (temp_d[5]),
    inv_lcldm_d        (temp_d[6]),
    inv_icldm_d        (temp_d[7]),
    inv_rcldm_d        (temp_d[8]),
    naai_d             (temp_d[9]),
    qc_relvar_d        (temp_d[10]),
    icldm_d            (temp_d[11]),
    lcldm_d            (temp_d[12]),
    rcldm_d            (temp_d[13]),
    t_d                (temp_d[14]),
    rho_d              (temp_d[15]),
    inv_rho_d          (temp_d[16]),
    qvs_d              (temp_d[17]),
    qvi_d              (temp_d[18]),
    supi_d             (temp_d[19]),
    rhofacr_d          (temp_d[20]),
    rhofaci_d          (temp_d[21]),
    acn_d              (temp_d[22]),
    qv_d               (temp_d[23]),
    th_d               (temp_d[24]),
    qc_d               (temp_d[25]),
    nc_d               (temp_d[26]),
    qr_d               (temp_d[27]),
    nr_d               (temp_d[28]),
    qitot_d            (temp_d[29]),
    nitot_d            (temp_d[30]),
    qirim_d            (temp_d[31]),
    birim_d            (temp_d[32]),
    xxlv_d             (temp_d[33]),
    xxls_d             (temp_d[34]),
    xlf_d              (temp_d[35]),
    qc_incld_d         (temp_d[36]),
    qr_incld_d         (temp_d[37]),
    qitot_incld_d      (temp_d[38]),
    qirim_incld_d      (temp_d[39]),
    nc_incld_d         (temp_d[40]),
    nr_incld_d         (temp_d[41]),
    nitot_incld_d      (temp_d[42]),
    birim_incld_d      (temp_d[43]),
    mu_c_d             (temp_d[44]),
    nu_d               (temp_d[45]),
    lamc_d             (temp_d[46]),
    cdist_d            (temp_d[47]),
    cdist1_d           (temp_d[48]),
    cdistr_d           (temp_d[49]),
    mu_r_d             (temp_d[50]),
    lamr_d             (temp_d[51]),
    logn0r_d           (temp_d[52]),
    cmeiout_d          (temp_d[53]),
    prain_d            (temp_d[54]),
    nevapr_d           (temp_d[55]),
    prer_evap_d        (temp_d[56]),
    vap_liq_exchange_d (temp_d[57]),
    vap_ice_exchange_d (temp_d[58]),
    liq_ice_exchange_d (temp_d[59]),
    pratot_d           (temp_d[60]),
    prctot_d           (temp_d[61]);

  // Call core function from kernel
  const auto dnu         = P3GlobalForFortran::dnu();
  const auto itab        = P3GlobalForFortran::itab();
  const auto itabcol     = P3GlobalForFortran::itabcol();
  const auto revap_table = P3GlobalForFortran::revap_table();
  bview_1d bools_d("bools", 1);
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part2(
      team, nk_pack, log_predictNc, dt, odt, dnu, itab, itabcol, revap_table,
      pres_d, pdel_d, dzq_d, ncnuc_d, exner_d, inv_exner_d, inv_lcldm_d, inv_icldm_d, inv_rcldm_d, naai_d, qc_relvar_d, icldm_d, lcldm_d, rcldm_d,
      t_d, rho_d, inv_rho_d, qvs_d, qvi_d, supi_d, rhofacr_d, rhofaci_d, acn_d,
      qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, xxlv_d, xxls_d, xlf_d, qc_incld_d, qr_incld_d,
      qitot_incld_d, qirim_incld_d, nc_incld_d, nr_incld_d, nitot_incld_d, birim_incld_d, mu_c_d, nu_d, lamc_d, cdist_d, cdist1_d,
      cdistr_d, mu_r_d, lamr_d, logn0r_d, cmeiout_d, prain_d, nevapr_d, prer_evap_d, vap_liq_exchange_d,
      vap_ice_exchange_d, liq_ice_exchange_d, pratot_d, prctot_d,
      bools_d(0));
  });

  // Sync back to host
  Kokkos::Array<view_1d, 48> inout_views = {
    t_d, rho_d, inv_rho_d, qvs_d, qvi_d, supi_d, rhofacr_d, rhofaci_d, acn_d,
    qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, xxlv_d, xxls_d, xlf_d, qc_incld_d, qr_incld_d,
    qitot_incld_d, qirim_incld_d, nc_incld_d, nr_incld_d, nitot_incld_d, birim_incld_d, mu_c_d, nu_d, lamc_d, cdist_d, cdist1_d,
    cdistr_d, mu_r_d, lamr_d, logn0r_d, cmeiout_d, prain_d, nevapr_d, prer_evap_d, vap_liq_exchange_d,
    vap_ice_exchange_d, liq_ice_exchange_d, pratot_d, prctot_d
  };

  pack::device_to_host({t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn,
        qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls, xlf, qc_incld, qr_incld,
        qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, mu_c, nu, lamc, cdist, cdist1,
        cdistr, mu_r, lamr, logn0r, cmeiout, prain, nevapr, prer_evap, vap_liq_exchange,
        vap_ice_exchange, liq_ice_exchange, pratot, prctot},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *log_hydrometeorsPresent = bools_h(0);
}

void p3_main_part3_f(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* exner, Real* lcldm, Real* rcldm,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th, Real* qc, Real* nc, Real* qr, Real* nr, Real* qitot, Real* nitot, Real* qirim, Real* birim, Real* xxlv, Real* xxls,
  Real* mu_c, Real* nu, Real* lamc, Real* mu_r, Real* lamr, Real* vap_liq_exchange,
  Real*  ze_rain, Real* ze_ice, Real* diag_vmi, Real* diag_effi, Real* diag_di, Real* diag_rhoi, Real* diag_ze, Real* diag_effc)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = scream::pack::npack<Spack>(nk);

  // Set up views
  Kokkos::Array<view_1d, P3MainPart3Data::NUM_ARRAYS> temp_d;

  pack::host_to_device({
      exner, lcldm, rcldm,
      rho, inv_rho, rhofaci,
      qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls,
      mu_c, nu, lamc, mu_r,
      lamr, vap_liq_exchange,
      ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc
    },
    nk, temp_d);

  view_1d
    exner_d            (temp_d[0]),
    lcldm_d            (temp_d[1]),
    rcldm_d            (temp_d[2]),
    rho_d              (temp_d[3]),
    inv_rho_d          (temp_d[4]),
    rhofaci_d          (temp_d[5]),
    qv_d               (temp_d[6]),
    th_d               (temp_d[7]),
    qc_d               (temp_d[8]),
    nc_d               (temp_d[9]),
    qr_d               (temp_d[10]),
    nr_d               (temp_d[11]),
    qitot_d            (temp_d[12]),
    nitot_d            (temp_d[13]),
    qirim_d            (temp_d[14]),
    birim_d            (temp_d[15]),
    xxlv_d             (temp_d[16]),
    xxls_d             (temp_d[17]),
    mu_c_d             (temp_d[18]),
    nu_d               (temp_d[19]),
    lamc_d             (temp_d[20]),
    mu_r_d             (temp_d[21]),
    lamr_d             (temp_d[22]),
    vap_liq_exchange_d (temp_d[23]),
    ze_rain_d          (temp_d[24]),
    ze_ice_d           (temp_d[25]),
    diag_vmi_d         (temp_d[26]),
    diag_effi_d        (temp_d[27]),
    diag_di_d          (temp_d[28]),
    diag_rhoi_d        (temp_d[29]),
    diag_ze_d          (temp_d[30]),
    diag_effc_d        (temp_d[31]);

  // Call core function from kernel
  const auto dnu         = P3GlobalForFortran::dnu();
  const auto itab        = P3GlobalForFortran::itab();
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part3(team, nk_pack, dnu, itab,
                                exner_d, lcldm_d, rcldm_d,
                                rho_d, inv_rho_d, rhofaci_d, qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, xxlv_d, xxls_d,
                                mu_c_d, nu_d, lamc_d, mu_r_d, lamr_d, vap_liq_exchange_d,
                                ze_rain_d, ze_ice_d, diag_vmi_d, diag_effi_d, diag_di_d, diag_rhoi_d, diag_ze_d, diag_effc_d);
  });

  // Sync back to host
  Kokkos::Array<view_1d, 29> inout_views = {
    rho_d, inv_rho_d, rhofaci_d, qv_d, th_d, qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, xxlv_d, xxls_d,
    mu_c_d, nu_d, lamc_d, mu_r_d, lamr_d, vap_liq_exchange_d,
    ze_rain_d, ze_ice_d, diag_vmi_d, diag_effi_d, diag_di_d, diag_rhoi_d, diag_ze_d, diag_effc_d
  };

  pack::device_to_host({
      rho, inv_rho, rhofaci, qv, th, qc, nc, qr, nr, qitot, nitot, qirim, birim, xxlv, xxls,
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange,
      ze_rain, ze_ice, diag_vmi, diag_effi, diag_di, diag_rhoi, diag_ze, diag_effc
    },
    nk, inout_views);
}

void p3_main_f(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th, Real* qv, Real dt, Real* qitot, Real* qirim, Real* nitot, Real* birim,
  Real* pres, Real* dzq, Real* ncnuc, Real* naai, Real* qc_relvar, Int it, Real* prt_liq, Real* prt_sol, Int its, Int ite, Int kts, Int kte, Real* diag_ze, Real* diag_effc,
  Real* diag_effi, Real* diag_vmi, Real* diag_di, Real* diag_rhoi, bool log_predictNc,
  Real* pdel, Real* exner, Real* cmeiout, Real* prain, Real* nevapr, Real* prer_evap, Real* rflx, Real* sflx, Real* rcldm, Real* lcldm, Real* icldm,
  Real* pratot, Real* prctot, Real* mu_c, Real* lamc, Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_2d    = typename P3F::view_2d<Spack>;
  using sview_1d   = typename P3F::view_1d<Real>;
  using sview_2d   = typename P3F::view_2d<Real>;
  using KT         = typename P3F::KT;

  scream_require_msg(its == 1, "its must be 1, got " << its);
  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  its  -= 1;
  ite  -= 1;
  kts  -= 1;
  kte  -= 1;

  const Int ni    = (ite - its) + 1;
  const Int nk    = (kte - kts) + 1;

  // Set up views, pretend all views are input views for the sake of initializing kokkos views
  Kokkos::Array<view_2d, P3MainData::NUM_ARRAYS> temp_d;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS> dim1_sizes;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS> dim2_sizes;
  Kokkos::Array<const Real*, P3MainData::NUM_ARRAYS> ptr_array = {
    pres, dzq, ncnuc, naai, pdel, exner, icldm, lcldm, rcldm, qc_relvar,
    qc, nc, qr, nr, qitot, qirim, nitot, birim, qv, th,
    diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, diag_rhoi, mu_c, lamc, cmeiout, prain, nevapr, prer_evap, pratot, prctot, liq_ice_exchange, vap_liq_exchange, vap_ice_exchange, rflx, sflx, prt_liq, prt_sol
  };

  for (size_t i = 0; i < P3MainData::NUM_ARRAYS; ++i) dim1_sizes[i] = ni;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS; ++i) dim2_sizes[i] = nk;

  dim2_sizes[37] = nk+1; // rflx
  dim2_sizes[38] = nk+1; // sflx
  dim1_sizes[39] = 1; dim2_sizes[39] = ni; // prt_liq
  dim1_sizes[40] = 1; dim2_sizes[40] = ni; // prt_sol

  // Initialize outputs to avoid uninitialized read warnings in memory checkers
  for (size_t i = P3MainData::NUM_INPUT_ARRAYS; i < P3MainData::NUM_ARRAYS; ++i) {
    for (size_t j = 0; j < dim1_sizes[i]*dim2_sizes[i]; ++j) {
      const_cast<Real*>(ptr_array[i])[j] = 0;
    }
  }

  pack::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

  int counter = 0;
  view_2d
    pres_d             (temp_d[counter++]),
    dzq_d              (temp_d[counter++]),
    ncnuc_d            (temp_d[counter++]),
    naai_d             (temp_d[counter++]),
    pdel_d             (temp_d[counter++]),
    exner_d            (temp_d[counter++]),
    icldm_d            (temp_d[counter++]),
    lcldm_d            (temp_d[counter++]),
    rcldm_d            (temp_d[counter++]),
    qc_relvar_d        (temp_d[counter++]),
    qc_d               (temp_d[counter++]),
    nc_d               (temp_d[counter++]),
    qr_d               (temp_d[counter++]),
    nr_d               (temp_d[counter++]),
    qitot_d            (temp_d[counter++]),
    qirim_d            (temp_d[counter++]),
    nitot_d            (temp_d[counter++]),
    birim_d            (temp_d[counter++]),
    qv_d               (temp_d[counter++]),
    th_d               (temp_d[counter++]),
    diag_ze_d          (temp_d[counter++]),
    diag_effc_d        (temp_d[counter++]),
    diag_effi_d        (temp_d[counter++]),
    diag_vmi_d         (temp_d[counter++]),
    diag_di_d          (temp_d[counter++]),
    diag_rhoi_d        (temp_d[counter++]),
    mu_c_d             (temp_d[counter++]),
    lamc_d             (temp_d[counter++]),
    cmeiout_d          (temp_d[counter++]),
    prain_d            (temp_d[counter++]),
    nevapr_d           (temp_d[counter++]),
    prer_evap_d        (temp_d[counter++]),
    pratot_d           (temp_d[counter++]),
    prctot_d           (temp_d[counter++]),
    liq_ice_exchange_d (temp_d[counter++]),
    vap_liq_exchange_d (temp_d[counter++]),
    vap_ice_exchange_d (temp_d[counter++]),
    rflx_d             (temp_d[counter++]),
    sflx_d             (temp_d[counter++]),
    prt_liq_temp_d     (temp_d[counter++]),
    prt_sol_temp_d     (temp_d[counter++]);

  // Special cases: prt_liq=1d<scalar>(ni), prt_sol=1d<scalar>(ni), col_location=2d<scalar>(ni, 3)
  sview_1d prt_liq_d("prt_liq_d", ni), prt_sol_d("prt_sol_d", ni);
  sview_2d col_location_d("col_location_d", ni, 3);

  Kokkos::parallel_for(ni, KOKKOS_LAMBDA(const Int& i) {
    prt_liq_d(i) = prt_liq_temp_d(0, i / Spack::n)[i % Spack::n];
    prt_sol_d(i) = prt_sol_temp_d(0, i / Spack::n)[i % Spack::n];

    for (int j = 0; j < 3; ++j) {
      col_location_d(i, j) = i+1;
    }
  });

  P3F::p3_main(pres_d, dzq_d, ncnuc_d, naai_d, qc_relvar_d, dt, ni, nk, it, log_predictNc, pdel_d, exner_d,
               icldm_d, lcldm_d, rcldm_d, col_location_d, qc_d, nc_d, qr_d, nr_d, qitot_d, qirim_d, nitot_d,
               birim_d, qv_d, th_d, prt_liq_d, prt_sol_d, diag_ze_d, diag_effc_d, diag_effi_d, diag_vmi_d, diag_di_d,
               diag_rhoi_d, mu_c_d, lamc_d, cmeiout_d, prain_d, nevapr_d, prer_evap_d, rflx_d, sflx_d, pratot_d,
               prctot_d, liq_ice_exchange_d, vap_liq_exchange_d, vap_ice_exchange_d);

  Kokkos::parallel_for(ni, KOKKOS_LAMBDA(const Int& i) {
    prt_liq_temp_d(0, i / Spack::n)[i % Spack::n] = prt_liq_d(i);
    prt_sol_temp_d(0, i / Spack::n)[i % Spack::n] = prt_sol_d(i);
  });

  // Sync back to host
  Kokkos::Array<view_2d, P3MainData::NUM_ARRAYS - 10> inout_views = {
    qc_d, nc_d, qr_d, nr_d, qitot_d, qirim_d, nitot_d, birim_d, qv_d, th_d,
    diag_ze_d, diag_effc_d, diag_effi_d, diag_vmi_d, diag_di_d, diag_rhoi_d, mu_c_d, lamc_d, cmeiout_d, prain_d, nevapr_d, prer_evap_d, pratot_d, prctot_d, liq_ice_exchange_d, vap_liq_exchange_d, vap_ice_exchange_d, rflx_d, sflx_d, prt_liq_temp_d, prt_sol_temp_d
  };
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS - 10> dim1_sizes_out;
  Kokkos::Array<size_t,  P3MainData::NUM_ARRAYS - 10> dim2_sizes_out;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS - 10; ++i) dim1_sizes_out[i] = ni;
  for (size_t i = 0; i < P3MainData::NUM_ARRAYS - 10; ++i) dim2_sizes_out[i] = nk;


  dim2_sizes_out[27] = nk+1; // rflx
  dim2_sizes_out[28] = nk+1; // sflx
  dim1_sizes_out[29] = 1; dim2_sizes_out[29] = ni; // prt_liq
  dim1_sizes_out[30] = 1; dim2_sizes_out[30] = ni; // prt_sol

  pack::device_to_host({
      qc, nc, qr, nr, qitot, qirim, nitot, birim, qv, th,
      diag_ze, diag_effc, diag_effi, diag_vmi, diag_di, diag_rhoi, mu_c, lamc, cmeiout, prain, nevapr, prer_evap, pratot, prctot, liq_ice_exchange, vap_liq_exchange, vap_ice_exchange, rflx, sflx, prt_liq, prt_sol
    },
    dim1_sizes_out, dim2_sizes_out, inout_views, true);
}

} // namespace p3
} // namespace scream
