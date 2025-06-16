#include "gw_test_data.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to gw fortran calls and vice versa
//

namespace scream {
namespace gw {

using GWF = Functions<Real, HostDevice>;
using GWC = typename GWF::C;

extern "C" {

void gwd_compute_tendencies_from_stress_divergence_c(Int ncol, Int ngwv, bool do_taper, Real dt, Real effgw, Int* tend_level, Real* lat, Real* dpm, Real* rdpm, Real* c, Real* ubm, Real* t, Real* nm, Real* xv, Real* yv, Real* tau, Real* gwut, Real* utgw, Real* vtgw);

void gw_init_c(Int pver_in, Int pgwv_in, Real dc_in, Real* cref_in, bool orographic_only, bool do_molec_diff_in, bool tau_0_ubc_in, Int nbot_molec_in, Int ktop_in, Int kbotbg_in, Real fcrit2_in, Real kwv_in, Real gravit_in, Real rair_in, Real* alpha_in);

void gw_prof_c(Int ncol, Real cpair, Real* t, Real* pmid, Real* pint, Real* rhoi, Real* ti, Real* nm, Real* ni);

void momentum_energy_conservation_c(Int ncol, Int* tend_level, Real dt, Real* taucd, Real* pint, Real* pdel, Real* u, Real* v, Real* dudt, Real* dvdt, Real* dsdt, Real* utgw, Real* vtgw, Real* ttgw);

void gwd_compute_stress_profiles_and_diffusivities_c(Int ncol, Int ngwv, Int* src_level, Real* ubi, Real* c, Real* rhoi, Real* ni, Real* kvtt, Real* t, Real* ti, Real* piln, Real* tau);

void gwd_project_tau_c(Int ncol, Int ngwv, Int* tend_level, Real* tau, Real* ubi, Real* c, Real* xv, Real* yv, Real* taucd);

void gwd_precalc_rhoi_c(Int pcnst, Int ncol, Int ngwv, Real dt, Int* tend_level, Real* pmid, Real* pint, Real* t, Real* gwut, Real* ubm, Real* nm, Real* rdpm, Real* c, Real* q, Real* dse, Real* egwdffi, Real* qtgw, Real* dttdf, Real* dttke, Real* ttgw);

void gw_drag_prof_c(Int pcnst, Int ncol, Int ngwv, Int* src_level, Int* tend_level, bool do_taper, Real dt, Real* lat, Real* t, Real* ti, Real* pmid, Real* pint, Real* dpm, Real* rdpm, Real* piln, Real* rhoi, Real* nm, Real* ni, Real* ubm, Real* ubi, Real* xv, Real* yv, Real effgw, Real* c, Real* kvtt, Real* q, Real* dse, Real* tau, Real* utgw, Real* vtgw, Real* ttgw, Real* qtgw, Real* taucd, Real* egwdffi, Real* gwut, Real* dttdf, Real* dttke);

void gw_front_init_c(Real taubgnd, Real frontgfc_in, Int kfront_in);

void gw_front_project_winds_c(Int ncol, Int kbot, Real* u, Real* v, Real* xv, Real* yv, Real* ubm, Real* ubi);

void gw_front_gw_sources_c(Int ncol, Int ngwv, Int kbot, Real* frontgf, Real* tau);

void gw_cm_src_c(Int ncol, Int ngwv, Int kbot, Real* u, Real* v, Real* frontgf, Int* src_level, Int* tend_level, Real* tau, Real* ubm, Real* ubi, Real* xv, Real* yv, Real* c);

} // extern "C" : end _c decls

// Wrapper around gw_init
void gw_init(GwInit& init)
{
  gw_init_c(init.pver, init.pgwv, init.dc, init.cref, init.orographic_only, init.do_molec_diff, init.tau_0_ubc, init.nbot_molec, init.ktop, init.kbotbg, init.fcrit2, init.kwv, GWC::gravit, GWC::Rair, init.alpha);
}

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_compute_tendencies_from_stress_divergence_c(d.ncol, d.ngwv, d.do_taper, d.dt, d.effgw, d.tend_level, d.lat, d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.xv, d.yv, d.tau, d.gwut, d.utgw, d.vtgw);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gw_prof(GwProfData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gw_prof_c(d.ncol, d.cpair, d.t, d.pmid, d.pint, d.rhoi, d.ti, d.nm, d.ni);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void momentum_energy_conservation(MomentumEnergyConservationData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  momentum_energy_conservation_c(d.ncol, d.tend_level, d.dt, d.taucd, d.pint, d.pdel, d.u, d.v, d.dudt, d.dvdt, d.dsdt, d.utgw, d.vtgw, d.ttgw);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gwd_compute_stress_profiles_and_diffusivities(GwdComputeStressProfilesAndDiffusivitiesData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_compute_stress_profiles_and_diffusivities_c(d.ncol, d.ngwv, d.src_level, d.ubi, d.c, d.rhoi, d.ni, d.kvtt, d.t, d.ti, d.piln, d.tau);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gwd_project_tau(GwdProjectTauData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_project_tau_c(d.ncol, d.ngwv, d.tend_level, d.tau, d.ubi, d.c, d.xv, d.yv, d.taucd);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gwd_precalc_rhoi(GwdPrecalcRhoiData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gwd_precalc_rhoi_c(d.pcnst, d.ncol, d.ngwv, d.dt, d.tend_level, d.pmid, d.pint, d.t, d.gwut, d.ubm, d.nm, d.rdpm, d.c, d.q, d.dse, d.egwdffi, d.qtgw, d.dttdf, d.dttke, d.ttgw);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gw_drag_prof(GwDragProfData& d)
{
  gw_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gw_drag_prof_c(d.pcnst, d.ncol, d.ngwv, d.src_level, d.tend_level, d.do_taper, d.dt, d.lat, d.t, d.ti, d.pmid, d.pint, d.dpm, d.rdpm, d.piln, d.rhoi, d.nm, d.ni, d.ubm, d.ubi, d.xv, d.yv, d.effgw, d.c, d.kvtt, d.q, d.dse, d.tau, d.utgw, d.vtgw, d.ttgw, d.qtgw, d.taucd, d.egwdffi, d.gwut, d.dttdf, d.dttke);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gw_front_init(GwFrontInitData& d)
{
  gw_init(d.init);
  gw_front_init_c(d.taubgnd, d.frontgfc_in, d.kfront_in);
}

void gw_front_project_winds(GwFrontProjectWindsData& d)
{
  gw_front_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gw_front_project_winds_c(d.ncol, d.kbot, d.u, d.v, d.xv, d.yv, d.ubm, d.ubi);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gw_front_gw_sources(GwFrontGwSourcesData& d)
{
  gw_front_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gw_front_gw_sources_c(d.ncol, d.ngwv, d.kbot, d.frontgf, d.tau);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void gw_cm_src(GwCmSrcData& d)
{
  gw_front_init(d.init);
  d.transpose<ekat::TransposeDirection::c2f>();
  gw_cm_src_c(d.ncol, d.ngwv, d.kbot, d.u, d.v, d.frontgf, d.src_level, d.tend_level, d.tau, d.ubm, d.ubi, d.xv, d.yv, d.c);
  d.transpose<ekat::TransposeDirection::f2c>();
}

// end _c impls

} // namespace gw
} // namespace scream
