#include "p3_functions_f90.hpp"

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/scream_kokkos_utils.hpp"
#include "share/scream_pack_kokkos.hpp"
#include "p3_f90.hpp"

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
                            Real* qcnuc_, Real* ncnuc_, Real* qisub_,
                            Real* nrshdr_, Real* qcheti_, Real* qrcol_,
                            Real* qcshd_, Real* qimlt_, Real* qccol_,
                            Real* qrheti_, Real* nimlt_, Real* nccol_,
                            Real* ncshdc_, Real* ncheti_, Real* nrcol_,
                            Real* nislf_, Real* qidep_, Real* nrheti_,
                            Real* nisub_, Real* qinuc_, Real* ninuc_,
                            Real* qiberg_);

void prevent_ice_overdepletion_c(Real pres, Real t, Real qv, Real xxls,
                                 Real odt, Real* qidep, Real* qisub);

void cloud_water_conservation_c(Real qc, Real qcnuc, Real dt, Real* qcaut, Real* qcacc, Real* qccol,
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
                                 Real qc_incld, Real* qcheti, Real* ncheti);

void rain_immersion_freezing_c(Real t, Real lamr, Real mu_r, Real cdistr,
                               Real qr_incld, Real* qrheti, Real* nrheti);

void droplet_self_collection_c(Real rho, Real inv_rho, Real qc_incld, Real mu_c,
                               Real nu, Real ncautc, Real* ncacc);

void cloud_rain_accretion_c(Real rho, Real inv_rho, Real qc_incld, Real nc_incld,
                            Real qr_incld, Real* qcacc, Real* ncacc);

void cloud_water_autoconversion_c(Real rho, Real qc_incld, Real nc_incld, Real* qcaut, Real* ncautc, Real* ncautr);

void rain_self_collection_c(Real rho, Real qr_incld, Real nr_incld, Real* nrslf);

void impose_max_total_ni_c(Real* nitot_local, Real max_total_Ni, Real inv_rho_local);

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
  Real qcacc, Real ncacc, Real qcaut, Real ncautc, Real qcnuc, Real ncautr,
  Real ncslf, Real  qrevp, Real nrevp, Real nrslf , bool log_predictNc,
  Real inv_rho, Real exner, Real xxlv, Real dt, Real* th, Real* qv,
  Real* qc, Real* nc, Real* qr, Real* nr);

void ice_deposition_sublimation_c(
  Real qitot_incld, Real nitot_incld, Real t, Real qvs, Real qvi, Real epsi,
  Real abi, Real qv, Real* qidep, Real* qisub, Real* nisub, Real* qiberg);

void compute_rain_fall_velocity_c(Real qr_incld, Real rcldm, Real rhofacr,
                                  Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* V_qr, Real* V_nr);

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

void droplet_activation_c(Real temp, Real pres, Real qv, Real qc,
                          Real inv_rho, Real sup, Real xxlv, Real npccn,
                          bool log_predictNc, Real odt,
                          Real* qcnuc, Real* ncnuc);

void ice_cldliq_wet_growth_c(Real rho, Real temp, Real pres, Real rhofaci, Real f1pr05,
                             Real f1pr14, Real xxlv, Real xlf, Real dv,
                             Real kap, Real mu, Real sc, Real qv, Real qc_incld,
                             Real qitot_incld, Real nitot_incld, Real qr_incld, bool* log_wetgrowth,
                             Real* qrcol, Real* qccol, Real* qwgrth, Real* nrshdr, Real* qcshd);

}

namespace scream {
namespace p3 {

// helper functions
namespace {

template <size_t N>
void gen_random_data(const std::array<std::pair<Real, Real>, N>& ranges,
                     const std::array<Real**, N>& ptrs,
                     Real* data, Int nk)
{
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
}

}

//
// In all C++ -> Fortran bridge functions you should see p3_init(true). P3 needs
// to be initialized since most of its function depend on global tables to be
// populated. The 'true' argument is to set p3 to use its fortran implementations
// instead of calling back to C++. We want this behavior since it doesn't make much
// sense for C++ to bridge over to fortran only to have fortran bridge back to C++.
// If the client wanted the C++ implementation, they should just call it directly.
//

void p3_init_a(P3InitAFortranData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  p3_init_a_c(d.itab.data(), d.itabcol.data());
}

void find_lookuptable_indices_1a(LookupIceData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  find_lookuptable_indices_1a_c(&d.dumi, &d.dumjj, &d.dumii, &d.dumzz,
                                &d.dum1, &d.dum4, &d.dum5, &d.dum6,
                                d.qitot, d.nitot, d.qirim, d.rhop);
}

void find_lookuptable_indices_1b(LookupIceDataB& d)
{
  p3_init(true);
  find_lookuptable_indices_1b_c(&d.dumj, &d.dum3, d.qr, d.nr);
}

void access_lookup_table(AccessLookupTableData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
  access_lookup_table_c(d.lid.dumjj, d.lid.dumii, d.lid.dumi, d.index,
                        d.lid.dum1, d.lid.dum4, d.lid.dum5, &d.proc);
}

void access_lookup_table_coll(AccessLookupTableCollData& d)
{
  p3_init(true); // need to initialize p3 first so that tables are loaded
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
  p3_init(true);
  back_to_cell_average_c(d.lcldm, d.rcldm, d.icldm, &d.qcacc, &d.qrevp,
    &d.qcaut, &d.ncacc, &d.ncslf, &d.ncautc, &d.nrslf, &d.nrevp, &d.ncautr,
    &d.qcnuc, &d.ncnuc, &d.qisub, &d.nrshdr, &d.qcheti, &d.qrcol, &d.qcshd,
    &d.qimlt, &d.qccol, &d.qrheti, &d.nimlt, &d.nccol, &d.ncshdc, &d.ncheti,
    &d.nrcol, &d.nislf, &d.qidep, &d.nrheti, &d.nisub, &d.qinuc, &d.ninuc,
    &d.qiberg);
}

void prevent_ice_overdepletion(PreventIceOverdepletionData& d)
{
  p3_init(true);
  prevent_ice_overdepletion_c(d.pres, d.t, d.qv, d.xxls, d.odt, &d.qidep,
                              &d.qisub);
}

void calc_rime_density(CalcRimeDensityData& d)
{
  p3_init(true);
  calc_rime_density_c(d.t, d.rhofaci, d.f1pr02, d.acn, d.lamc, d.mu_c,
                      d.qc_incld, d.qccol, &d.vtrmi1, &d.rhorime_c);
}

void cldliq_immersion_freezing(CldliqImmersionFreezingData& d)
{
  p3_init(true);
  cldliq_immersion_freezing_c(d.t, d.lamc, d.mu_c, d.cdist1, d.qc_incld,
                              &d.qcheti, &d.ncheti);
}

void droplet_self_collection(DropletSelfCollectionData& d)
{
  p3_init(true);
  droplet_self_collection_c(d.rho, d.inv_rho, d.qc_incld, d.mu_c, d.nu, d.ncautc,
                            &d.ncslf);
}

void rain_immersion_freezing(RainImmersionFreezingData& d)
{
  p3_init(true);
  rain_immersion_freezing_c(d.t, d.lamr, d.mu_r, d.cdistr, d.qr_incld,
                            &d.qrheti, &d.nrheti);
}

void cloud_rain_accretion(CloudRainAccretionData& d)
{
  p3_init(true);
  cloud_rain_accretion_c(d.rho, d.inv_rho, d.qc_incld, d.nc_incld, d.qr_incld,
                         &d.qcacc, &d.ncacc);
}

void cloud_water_conservation(CloudWaterConservationData& d){
  p3_init(true);
  cloud_water_conservation_c(d.qc, d.qcnuc, d.dt, &d.qcaut, &d.qcacc, &d.qccol, &d.qcheti,
  &d.qcshd, &d.qiberg, &d.qisub, &d.qidep);
}

void rain_water_conservation(RainWaterConservationData& d){
  p3_init(true);
  rain_water_conservation_c(d.qr, d.qcaut, d.qcacc, d.qimlt, d.qcshd, d.dt, &d.qrevp, &d.qrcol, &d.qrheti);
}

void ice_water_conservation(IceWaterConservationData& d){
  p3_init(true);
  ice_water_conservation_c(d.qitot, d.qidep, d.qinuc, d.qiberg, d.qrcol, d.qccol, d.qrheti,
    d.qcheti, d.dt, &d.qisub, &d.qimlt);
}

void cloud_water_autoconversion(CloudWaterAutoconversionData& d){
  p3_init(true);
  cloud_water_autoconversion_c(d.rho, d.qc_incld, d.nc_incld, &d.qcaut, &d.ncautc, &d.ncautr);
}

void rain_self_collection(RainSelfCollectionData& d){
  p3_init(true);
  rain_self_collection_c(d.rho, d.qr_incld, d.nr_incld, &d.nrslf);
}

void impose_max_total_Ni(ImposeMaxTotalNiData& d){
  p3_init(true);
  impose_max_total_ni_c(&d.nitot_local, d.max_total_Ni, d.inv_rho_local);
}

void get_cloud_dsd2(GetCloudDsd2Data& d)
{
  p3_init(true);
  Real nc_in = d.nc_in;
  get_cloud_dsd2_c(d.qc, &nc_in, &d.mu_c, d.rho, &d.nu, &d.lamc, &d.cdist, &d.cdist1, d.lcldm);
  d.nc_out = nc_in;
}

void get_rain_dsd2(GetRainDsd2Data& d)
{
  p3_init(true);
  Real nr_in = d.nr_in;
  get_rain_dsd2_c(d.qr, &nr_in, &d.mu_r, &d.lamr, &d.cdistr, &d.logn0r, d.rcldm);
  d.nr_out = nr_in;
}

void ice_cldliq_collection(IceCldliqCollectionData& d)
{
  p3_init(true);
  ice_cldliq_collection_c(d.rho, d.temp, d.rhofaci, d.f1pr04,
                          d.qitot_incld, d.qc_incld, d.nitot_incld, d.nc_incld,
                          &d.qccol, &d.nccol, &d.qcshd, &d.ncshdc);
}

void ice_rain_collection(IceRainCollectionData& d)
{
  p3_init(true);
  ice_rain_collection_c(d.rho, d.temp, d.rhofaci, d.logn0r, d.f1pr07, d.f1pr08,
                        d.qitot_incld, d.nitot_incld, d.qr_incld,
                        &d.qrcol, &d.nrcol);
}

void ice_self_collection(IceSelfCollectionData& d)
{
  p3_init(true);
  ice_self_collection_c(d.rho, d.rhofaci, d.f1pr03, d.eii, d.qirim_incld,
                        d.qitot_incld, d.nitot_incld,
                        &d.nislf);
}


void ice_relaxation_timescale(IceRelaxationData& d)
{
  p3_init(true);
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
  p3_init(true);
  calc_liq_relaxation_timescale_c(d.rho, d.f1r, d.f2r, d.dv, d.mu, d.sc, d.mu_r,
    d.lamr, d.cdistr, d.cdist, d.qr_incld, d.qc_incld, &d.epsr, &d.epsc);
}

void ice_nucleation(IceNucleationData& d)
{
  p3_init(true);
  ice_nucleation_c(d.temp, d.inv_rho, d.nitot, d.naai,
                   d.supi, d.odt, d.log_predictNc,&d.qinuc, &d.ninuc);
}

void droplet_activation(DropletActivationData& d)
{
  p3_init(true);

  droplet_activation_c(d.temp, d.pres, d.qv, d.qc,
                       d.inv_rho, d.sup, d.xxlv, d.npccn,
                       d.log_predictNc, d.odt,
                       &d.qcnuc, &d.ncnuc);

}

void ice_cldliq_wet_growth(IceWetGrowthData& d)
{
  p3_init(true);

  ice_cldliq_wet_growth_c(d.rho, d.temp, d.pres, d.rhofaci, d.f1pr05,
                          d.f1pr14, d.xxlv, d.xlf, d.dv,
                          d.kap, d.mu, d.sc, d.qv, d.qc_incld,
                          d.qitot_incld, d.nitot_incld, d.qr_incld, &d.log_wetgrowth,
                          &d.qrcol, &d.qccol, &d.qwgrth, &d.nrshdr, &d.qcshd);
}

  void  update_prognostic_ice(P3UpdatePrognosticIceData& d){
    p3_init(true);
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
  p3_init(true);
  evaporate_sublimate_precip_c(d.qr_incld, d.qc_incld, d.nr_incld, d.qitot_incld,
			       d.lcldm, d.rcldm, d.qvs, d.ab, d.epsr, d.qv,
			       &d.qrevp, &d.nrevp);
}

void  update_prognostic_liquid(P3UpdatePrognosticLiqData& d){
  p3_init(true);
  update_prognostic_liquid_c(d.qcacc, d.ncacc, d.qcaut, d.ncautc, d.qcnuc, d.ncautr,
			      d.ncslf, d. qrevp, d.nrevp, d.nrslf , d.log_predictNc,
			      d.inv_rho, d.exner, d.xxlv, d.dt, &d.th, &d.qv,
			      &d.qc, &d.nc, &d.qr, &d.nr);
  }

void ice_deposition_sublimation(IceDepSublimationData& d){
  p3_init(true);
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
  p3_init(true);
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
  p3_init(true);
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
  p3_init(true);
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
  p3_init(true);
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
  p3_init(true);
  rain_sedimentation_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                       d.qr_incld, d.rho, d.inv_rho, d.rhofacr, d.rcldm, d.inv_dzq,
                       d.dt, d.odt,
                       d.qr, d.nr, d.nr_incld, d.mu_r, d.lamr, &d.prt_liq, d.rflx, d.qr_tend, d.nr_tend);
}

void calc_bulk_rho_rime(CalcBulkRhoRimeData& d)
{
  p3_init(true);
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
  p3_init(true);
  homogeneous_freezing_c(d.kts, d.kte, d.ktop, d.kbot, d.kdir,
                         d.t, d.exner, d.xlf,
                         d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.th);
}

void compute_rain_fall_velocity(ComputeRainFallVelocityData& d)
{
  p3_init(true);
  compute_rain_fall_velocity_c(d.qr_incld, d.rcldm, d.rhofacr,
                               &d.nr, &d.nr_incld, &d.mu_r, &d.lamr, &d.V_qr, &d.V_nr);
}

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

void find_lookuptable_indices_1a_f(Int* dumi, Int* dumjj, Int* dumii, Int* dumzz,
                                   Real* dum1, Real* dum4, Real* dum5, Real* dum6,
                                   Real qitot_, Real nitot_, Real qirim_, Real rhop_)
{
  using P3F = Functions<Real, DefaultDevice>;
  using TableIce = typename P3F::TableIce;

  typename P3F::Smask qiti_gt_small(qitot_ > P3F::C::QSMALL);
  typename P3F::Spack qitot(qitot_), nitot(nitot_), qirim(qirim_), rhop(rhop_);
  typename P3F::view_1d<TableIce> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    P3F::lookup_ice(qiti_gt_small, qitot, nitot, qirim, rhop, t_d(0));
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

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);

  typename P3F::Spack qr(qr_), nr(nr_);
  typename P3F::view_1d<TableRain> t_d("t_h", 1);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    P3F::lookup_rain(qiti_gt_small, qr, nr, t_d(0));
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

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);
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
    value = P3F::apply_table_ice(qiti_gt_small, adjusted_index, itab, t)[0];
  }, result);
  *proc = result;
}

void access_lookup_table_coll_f(Int dumjj, Int dumii, Int dumj, Int dumi, Int index,
                                Real dum1, Real dum3, Real dum4, Real dum5, Real* proc)
{
  using P3F = Functions<Real, DefaultDevice>;

  // we can assume fortran would not be calling this routine if qiti_gt_small was not true
  typename P3F::Smask qiti_gt_small(true);

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
    value = P3F::apply_table_coll(qiti_gt_small, adjusted_index, itabcol, ti, tr)[0];
  }, result);
  *proc = result;
}

void get_cloud_dsd2_f(Real qc_, Real* nc_, Real* mu_c_, Real rho_, Real* nu_, Real* lamc_,
                      Real* cdist_, Real* cdist1_, Real lcldm_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::Smask qc_gt_small(qc_ > P3F::C::QSMALL);
  typename P3F::view_1d<Real> t_d("t_d", 6);
  auto t_h = Kokkos::create_mirror_view(t_d);

  Real local_nc = *nc_;
  const auto dnu = P3GlobalForFortran::dnu();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qc(qc_), nc(local_nc), rho(rho_), lcldm(lcldm_);
    typename P3F::Spack mu_c, nu, lamc, cdist, cdist1;

    P3F::get_cloud_dsd2(qc_gt_small, qc, nc, mu_c, rho, nu, dnu, lamc, cdist, cdist1, lcldm);

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

  typename P3F::Smask qr_gt_small(qr_ > P3F::C::QSMALL);
  typename P3F::view_1d<Real> t_d("t_d", 5);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_nr = *nr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    typename P3F::Spack qr(qr_), rcldm(rcldm_), nr(local_nr);
    typename P3F::Spack lamr, mu_r, cdistr, logn0r;

    P3F::get_rain_dsd2(qr_gt_small, qr, nr, mu_r, lamr, cdistr, logn0r, rcldm);

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
	qiberg(qiberg_),  exner(exner_),  xlf(xlf_),  xxls(xxls_),  nmltratio(nmltratio_),
	rhorime_c(rhorime_c_);
      bool log_predictNc(log_predictNc_), log_wetgrowth(log_wetgrowth_);
      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr),
	nr(local_nr), qitot(local_qitot), nitot(local_nitot), qirim(local_qirim), birim(local_birim);

      P3F::update_prognostic_ice(qcheti, qccol, qcshd, nccol, ncheti,ncshdc,
				 qrcol,   nrcol,  qrheti,  nrheti,  nrshdr,
				 qimlt,  nimlt,  qisub,  qidep,  qinuc,  ninuc,
				 nislf,  nisub,  qiberg,  exner,  xxls,  xlf,
				 log_predictNc, log_wetgrowth,  dt,  nmltratio,
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

void update_prognostic_liquid_f(Real qcacc_, Real ncacc_, Real qcaut_, Real ncautc_, Real qcnuc_, Real ncautr_,
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
      typename P3F::Spack qcacc(qcacc_), ncacc(ncacc_), qcaut(qcaut_), ncautc(ncautc_), qcnuc(qcnuc_),
	ncautr(ncautr_), ncslf(ncslf_),  qrevp( qrevp_), nrevp(nrevp_), nrslf(nrslf_), inv_rho(inv_rho_),
	exner(exner_), xxlv(xxlv_);

      bool log_predictNc(log_predictNc_);

      typename P3F::Scalar dt(dt_);

      typename P3F::Spack th(local_th), qv(local_qv), qc(local_qc), nc(local_nc), qr(local_qr), nr(local_nr);

      P3F::update_prognostic_liquid(qcacc, ncacc, qcaut, ncautc, qcnuc, ncautr,
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

  // Setup views
  Kokkos::Array<view_1d, 3> temp_d;
  Kokkos::Array<view_1d, N> fluxes_d, vs_d, qnx_d;

  pack::host_to_device({rho, inv_rho, inv_dzq}, nk, temp_d);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dzq_d(temp_d[2]);

  pack::host_to_device(ptr_to_arr<N>((const Real**)fluxes), nk, fluxes_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)vs)    , nk, vs_d);
  pack::host_to_device(ptr_to_arr<N>((const Real**)qnx)   , nk, qnx_d);

  // Call core function from kernel
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
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
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
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
  using uview_1d = typename P3F::uview_1d<Spack>;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;

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
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_liq_k) {

    uview_1d
      uqc_incld_d(temp_d[0]),
      urho_d     (temp_d[1]),
      uinv_rho_d (temp_d[2]),
      ulcldm_d   (temp_d[3]),
      uacn_d     (temp_d[4]),
      uinv_dzq_d (temp_d[5]),
      uqc_d      (temp_d[6]),
      unc_d      (temp_d[7]),
      unc_incld_d(temp_d[8]),
      umu_c_d    (temp_d[9]),
      ulamc_d    (temp_d[10]),
      uqc_tend_d (temp_d[11]),
      unc_tend_d (temp_d[12]);

    P3F::cloud_sedimentation(
      uqc_incld_d, urho_d, uinv_rho_d, ulcldm_d, uacn_d, uinv_dzq_d, dnu,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, odt, log_predictNc,
      uqc_d, unc_d, unc_incld_d, umu_c_d, ulamc_d, uqc_tend_d, unc_tend_d,
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
  using uview_1d   = typename P3F::uview_1d<Spack>;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;

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
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 6, policy);
  Real my_prt_sol = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_sol_k) {

    uview_1d
      urho_d        (temp_d[0]),
      uinv_rho_d    (temp_d[1]),
      urhofaci_d    (temp_d[2]),
      uicldm_d      (temp_d[3]),
      uinv_dzq_d    (temp_d[4]),
      uqitot_d      (temp_d[5]),
      uqitot_incld_d(temp_d[6]),
      unitot_d      (temp_d[7]),
      uqirim_d      (temp_d[8]),
      uqirim_incld_d(temp_d[9]),
      ubirim_d      (temp_d[10]),
      ubirim_incld_d(temp_d[11]),
      unitot_incld_d(temp_d[12]),
      uqi_tend_d    (temp_d[13]),
      uni_tend_d    (temp_d[14]);

    P3F::ice_sedimentation(
      urho_d, uinv_rho_d, urhofaci_d, uicldm_d, uinv_dzq_d,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, odt,
      uqitot_d, uqitot_incld_d, unitot_d, unitot_incld_d, uqirim_d, uqirim_incld_d, ubirim_d, ubirim_incld_d,
      uqi_tend_d, uni_tend_d, itab,
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
  using uview_1d   = typename P3F::uview_1d<Spack>;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;

  // Set up views
  Kokkos::Array<view_1d, RainSedData::NUM_ARRAYS> temp_d;
  Kokkos::Array<size_t, RainSedData::NUM_ARRAYS> sizes;
  for (int i = 0; i < RainSedData::NUM_ARRAYS; ++i) sizes[i] = nk;
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
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
  WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Real my_prt_liq = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& prt_liq_k) {

    uview_1d
      uqr_incld_d   (temp_d[0]),
      urho_d        (temp_d[1]),
      uinv_rho_d    (temp_d[2]),
      urhofacr_d    (temp_d[3]),
      urcldm_d      (temp_d[4]),
      uinv_dzq_d    (temp_d[5]),
      uqr_d         (temp_d[6]),
      unr_d         (temp_d[7]),
      unr_incld_d   (temp_d[8]),
      umu_r_d       (temp_d[9]),
      ulamr_d       (temp_d[10]),
      uqr_tend_d    (temp_d[11]),
      unr_tend_d    (temp_d[12]),
      urflx_d       (temp_d[13]);

    P3F::rain_sedimentation(
      urho_d, uinv_rho_d, urhofacr_d, urcldm_d, uinv_dzq_d, uqr_incld_d,
      team, wsm.get_workspace(team), vn_table, vm_table,
      nk, ktop, kbot, kdir, dt, odt,
      uqr_d, unr_d, unr_incld_d, umu_r_d, ulamr_d, urflx_d, uqr_tend_d, unr_tend_d,
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
                            Real* qcnuc_, Real* ncnuc_, Real* qisub_,
                            Real* nrshdr_, Real* qcheti_, Real* qrcol_,
                            Real* qcshd_, Real* qimlt_, Real* qccol_,
                            Real* qrheti_, Real* nimlt_, Real* nccol_,
                            Real* ncshdc_, Real* ncheti_, Real* nrcol_,
                            Real* nislf_, Real* qidep_, Real* nrheti_,
                            Real* nisub_, Real* qinuc_, Real* ninuc_,
                            Real* qiberg_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 31);
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
  Real local_qcnuc = *qcnuc_;
  Real local_ncnuc = *ncnuc_;
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
      ncautr(local_ncautr), qcnuc(local_qcnuc), ncnuc(local_ncnuc), qisub(local_qisub),
      nrshdr(local_nrshdr), qcheti(local_qcheti), qrcol(local_qrcol), qcshd(local_qcshd),
      qimlt(local_qimlt), qccol(local_qccol), qrheti(local_qrheti), nimlt(local_nimlt),
      nccol(local_nccol), ncshdc(local_ncshdc), ncheti(local_ncheti), nrcol(local_nrcol),
      nislf(local_nislf), qidep(local_qidep), nrheti(local_nrheti), nisub(local_nisub),
      qinuc(local_qinuc), ninuc(local_ninuc), qiberg(local_qiberg);

    P3F::back_to_cell_average(lcldm, rcldm, icldm, qcacc, qrevp, qcaut,
      ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qcnuc, ncnuc, qisub,
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
    t_d(9) = qcnuc[0];
    t_d(10) = ncnuc[0];
    t_d(11) = qisub[0];
    t_d(12) = nrshdr[0];
    t_d(13) = qcheti[0];
    t_d(14) = qrcol[0];
    t_d(15) = qcshd[0];
    t_d(16) = qimlt[0];
    t_d(17) = qccol[0];
    t_d(18) = qrheti[0];
    t_d(19) = nimlt[0];
    t_d(20) = nccol[0];
    t_d(21) = ncshdc[0];
    t_d(22) = ncheti[0];
    t_d(23) = nrcol[0];
    t_d(24) = nislf[0];
    t_d(25) = qidep[0];
    t_d(26) = nrheti[0];
    t_d(27) = nisub[0];
    t_d(28) = qinuc[0];
    t_d(29) = ninuc[0];
    t_d(30) = qiberg[0];

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
  *qcnuc_ = t_h(9);
  *ncnuc_ = t_h(10);
  *qisub_ = t_h(11);
  *nrshdr_ = t_h(12);
  *qcheti_ = t_h(13);
  *qrcol_ = t_h(14);
  *qcshd_ = t_h(15);
  *qimlt_ = t_h(16);
  *qccol_ = t_h(17);
  *qrheti_ = t_h(18);
  *nimlt_ = t_h(19);
  *nccol_ = t_h(20);
  *ncshdc_ = t_h(21);
  *ncheti_ = t_h(22);
  *nrcol_ = t_h(23);
  *nislf_ = t_h(24);
  *qidep_ = t_h(25);
  *nrheti_ = t_h(26);
  *nisub_ = t_h(27);
  *qinuc_ = t_h(28);
  *ninuc_ = t_h(29);
  *qiberg_ = t_h(30);
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
      typename P3F::Spack pres(pres_), t(t_), qv(qv_), xxls(xxls_), odt(odt_),
                          qidep(local_qidep), qisub(local_qisub);
      P3F::prevent_ice_overdepletion(pres, t, qv, xxls, odt, qidep, qisub);

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
  Real t_, Real lamc_, Real mu_c_, Real cdist1_, Real qc_incld_,
  Real* qcheti_, Real* ncheti_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 2);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qcheti = *qcheti_;
  Real local_ncheti = *ncheti_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack t(t_), lamc(lamc_), mu_c(mu_c_), cdist1(cdist1_),
                          qc_incld(qc_incld_), qcheti(local_qcheti), ncheti(local_ncheti);
      P3F::cldliq_immersion_freezing(t, lamc, mu_c, cdist1, qc_incld,
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
  Real rho_, Real inv_rho_, Real qc_incld_, Real nc_incld_, Real qr_incld_,
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
                          qcacc(local_qcacc), ncacc(local_ncacc);
      P3F::cloud_rain_accretion(rho, inv_rho, qc_incld, nc_incld, qr_incld,
                                qcacc, ncacc);

      t_d(0) = qcacc[0];
      t_d(1) = ncacc[0];

    });
  Kokkos::deep_copy(t_h, t_d);

  *qcacc_ = t_h(0);
  *ncacc_ = t_h(1);
}

void cloud_water_autoconversion_f(
  Real rho_, Real qc_incld_, Real nc_incld_, Real* qcaut_, Real* ncautc_, Real* ncautr_)
{
  using P3F = Functions<Real, DefaultDevice>;

  typename P3F::view_1d<Real> t_d("t_h", 3);
  auto t_h = Kokkos::create_mirror_view(t_d);
  Real local_qcaut = *qcaut_;
  Real local_ncautc = *ncautc_;
  Real local_ncautr = *ncautr_;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      typename P3F::Spack rho(rho_), qc_incld(qc_incld_), nc_incld(nc_incld_), qcaut(local_qcaut), ncautc(local_ncautc), ncautr(local_ncautr);
      P3F::cloud_water_autoconversion(rho, qc_incld, nc_incld, qcaut, ncautc, ncautr);

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
    Spack max_total_Ni(max_total_Ni_); 
    Spack inv_rho_local(inv_rho_local_);

    P3F::impose_max_total_Ni(nitot_local, max_total_Ni, inv_rho_local);
    t_d(0) = nitot_local[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *nitot_local_ = t_h(0); 
}

void calc_bulk_rho_rime_f(Real qi_tot_, Real* qi_rim_, Real* bi_rim_, Real* rho_rime_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using Smask   = typename P3F::Smask;
  using view_1d = typename P3F::view_1d<Real>;

  Real local_qi_rim = *qi_rim_, local_bi_rim = *bi_rim_;
  view_1d t_d("t_d", 3);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Smask qi_gt_small(qi_tot_ > P3F::C::QSMALL);
    Spack qi_tot(qi_tot_), qi_rim(local_qi_rim), bi_rim(local_bi_rim);

    const auto result = P3F::calc_bulk_rho_rime(qi_gt_small, qi_tot, qi_rim, bi_rim);
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
  using uview_1d   = typename P3F::uview_1d<Spack>;

  scream_require_msg(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;

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
  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    uview_1d
      ut_d    (temp_d[0]),
      uexner_d(temp_d[1]),
      uxlf_d  (temp_d[2]),
      uqc_d   (temp_d[3]),
      unc_d   (temp_d[4]),
      uqr_d   (temp_d[5]),
      unr_d   (temp_d[6]),
      uqitot_d(temp_d[7]),
      unitot_d(temp_d[8]),
      uqirim_d(temp_d[9]),
      ubirim_d(temp_d[10]),
      uth_d   (temp_d[11]);

    P3F::homogeneous_freezing(
      ut_d, uexner_d, uxlf_d,
      team,
      nk, ktop, kbot, kdir,
      uqc_d, unc_d, uqr_d, unr_d, uqitot_d, unitot_d, uqirim_d, ubirim_d, uth_d);
  });

  // Sync back to host
  Kokkos::Array<view_1d, 9> inout_views = {qc_d, nc_d, qr_d, nr_d, qitot_d, nitot_d, qirim_d, birim_d, th_d};

  pack::device_to_host({qc, nc, qr, nr, qitot, nitot, qirim, birim, th}, nk, inout_views);
}

void compute_rain_fall_velocity_f(Real qr_incld_, Real rcldm_, Real rhofacr_,
                                  Real* nr_, Real* nr_incld_, Real* mu_r_, Real* lamr_, Real* V_qr_, Real* V_nr_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack   = typename P3F::Spack;
  using Smask   = typename P3F::Smask;
  using view_1d = typename P3F::view_1d<Real>;

  Real local_nr = *nr_, local_nr_incld = *nr_incld_;
  view_1d t_d("t_d", 6);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  const auto vn_table = P3GlobalForFortran::vn_table();
  const auto vm_table = P3GlobalForFortran::vm_table();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Smask qr_gt_small(qr_incld_ > P3F::C::QSMALL);
    Spack qr_incld(qr_incld_), rcldm(rcldm_), rhofacr(rhofacr_), nr(local_nr), nr_incld(local_nr_incld),
      mu_r, lamr, V_qr, V_nr;

    P3F::compute_rain_fall_velocity(qr_gt_small, vn_table, vm_table,
                                    qr_incld, rcldm, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr);
    t_d(0) = nr[0];
    t_d(1) = nr_incld[0];
    t_d(2) = mu_r[0];
    t_d(3) = lamr[0];
    t_d(4) = V_qr[0];
    t_d(5) = V_nr[0];
  });
  Kokkos::deep_copy(t_h, t_d);

  *nr_       = t_h(0);
  *nr_incld_ = t_h(1);
  *mu_r_     = t_h(2);
  *lamr_     = t_h(3);
  *V_qr_     = t_h(4);
  *V_nr_     = t_h(5);
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

    Spack rho{rho_}, f1r{f1r_}, f2r{f2r_}, dv{dv_},
          mu{mu_}, sc{sc_}, mu_r{mu_r_}, lamr{lamr_}, cdistr{cdistr_},
          cdist{cdist_}, qr_incld{qr_incld_}, qc_incld{qc_incld_};

    Spack epsr{0.0}, epsc{0.0};

    P3F::calc_liq_relaxation_timescale(revap_table, rho, f1r, f2r, dv, mu, sc,
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
  using Smask        = typename P3F::Smask;
  using view_1d      = typename P3F::view_1d<Real>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Spack temp{temp_}, inv_rho{inv_rho_}, nitot{nitot_}, naai{naai_}, supi{supi_}, odt{odt_};
    Smask log_predictNc{log_predictNc_};
    Spack qinuc{0.0}, ninuc{0.0};

    P3F::ice_nucleation(temp, inv_rho, nitot, naai, supi, odt, log_predictNc,
                        qinuc, ninuc);

    t_d(0) = qinuc[0];
    t_d(1) = ninuc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qinuc_         = t_h(0);
  *ninuc_         = t_h(1);
}

void droplet_activation_f(Real temp_, Real pres_, Real qv_, Real qc_,
                          Real inv_rho_, Real sup_, Real xxlv_, Real npccn_,
                          bool log_predictNc_, Real odt_,
                          Real* qcnuc_, Real* ncnuc_)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack        = typename P3F::Spack;
  using Smask        = typename P3F::Smask;
  using view_1d      = typename P3F::view_1d<Real>;
  using bool_view_1d = typename P3F::view_1d<bool>;

  view_1d t_d("t_d", 2);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {

    Real qcnuc_loc{*qcnuc_}, ncnuc_loc{*ncnuc_};

    Spack temp{temp_}, pres{pres_}, qv{qv_}, qc{qc_}, inv_rho{inv_rho_}, sup{sup_}, xxlv{xxlv_},
          npccn{npccn_}, odt{odt_};

    Smask log_predictNc{log_predictNc_};

    Spack qcnuc{qcnuc_loc}, ncnuc{ncnuc_loc};

    P3F::droplet_activation(temp, pres, qv, qc, inv_rho, sup, xxlv, npccn, log_predictNc, odt, qcnuc, ncnuc);

    t_d(0) = qcnuc[0];
    t_d(1) = ncnuc[0];
  });

  Kokkos::deep_copy(t_h, t_d);

  *qcnuc_  = t_h(0);
  *ncnuc_  = t_h(1);
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

void cloud_water_conservation_f(Real qc_, Real qcnuc_, Real dt, Real* qcaut_, Real* qcacc_, Real* qccol_,
  Real* qcheti_, Real* qcshd_, Real* qiberg_, Real* qisub_, Real* qidep_)
  {
    using P3F = Functions<Real, HostDevice>;
    using Spack   = typename P3F::Spack;

    Spack qc(qc_), qcnuc(qcnuc_), qcaut(*qcaut_), qcacc(*qcacc_), qccol(*qccol_), qcheti(*qcheti_);
    Spack qcshd(*qcshd_), qiberg(*qiberg_), qisub(*qisub_), qidep(*qidep_);

    P3F::cloud_water_conservation(qc, qcnuc, dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);
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

} // namespace p3
} // namespace scream
