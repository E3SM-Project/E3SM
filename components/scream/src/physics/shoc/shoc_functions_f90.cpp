#include "shoc_functions_f90.hpp"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/scream_pack_kokkos.hpp"
#include "shoc_f90.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A C interface to SHOC fortran calls. The stubs below will link to fortran definitions in shoc_iso_c.f90
//
extern "C" {

void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                 Real zvir, Real latvap, Real latice, Real karman,
                 Real* pref_mid, int nbot_shoc, int ntop_shoc);

void shoc_grid_c(int shcol, int nlev, int nlevi, Real *zt_grid, Real *zi_grid,
                 Real *pdel, Real *dz_zt, Real *dzi_zi, Real *rho_zt);

void update_host_dse_c(Int shcol, Int nlev, Real *thlm, Real *shoc_ql,
                       Real *exner, Real *zt_grid, Real *phis, Real *host_dse);

void shoc_energy_integrals_c(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int);
			     
void shoc_energy_total_fixer(Int shcol, Int nlev, Int nlevi, Real dtime, Real nadv,
                             Real *zt_grid, Real *zi_grid,
                             Real *se_b, Real *ke_b, Real *wv_b, Real *wl_b,
                             Real *se_a, Real *ke_a, Real *wv_a, Real *wl_a,
                             Real *wthl_sfc, Real *wqw_sfc, Real *rho_zt,
                             Real *te_a, Real *te_b)			     

void calc_shoc_varorcovar_c(Int shcol, Int nlev, Int nlevi,  Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
			    Real *invar1, Real *invar2, Real *varorcovar);

void integ_column_stability_c(Int nlev, Int shcol, Real *dz_zt, Real *pres,
			      Real *brunt, Real *brunt_int);

void compute_shr_prod_c(Int nlevi, Int nlev, Int shcol, Real *dz_zi,
                        Real *u_wind, Real *v_wind, Real *sterm);

void isotropic_ts_c(Int nlev, Int shcol, Real *brunt_int, Real *tke,
                    Real *a_diss, Real *brunt, Real *isotropy);

void adv_sgs_tke_c(Int nlev, Int shcol, Real dtime, Real *shoc_mix,
                   Real *wthv_sec, Real *sterm_zt, Real *tk,
                   Real *tke, Real *a_diss);

void eddy_diffusivities_c(Int nlev, Int shcol, Real *obklen, Real *pblh,
                          Real *zt_grid, Real *shoc_mix, Real *sterm_zt,
                          Real *isotropy, Real *tke, Real *tkh, Real *tk);

void calc_shoc_vertflux_c(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
}

namespace scream {
namespace shoc {

//
// Data struct
//

SHOCDataBase::SHOCDataBase(Int shcol_, Int nlev_, Int nlevi_,
                           const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i) :
  shcol(shcol_),
  nlev(nlev_),
  nlevi(nlevi_),
  m_total(shcol_ * nlev_),
  m_totali(shcol_ * nlevi_),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_data(m_ptrs.size() * m_total + m_ptrs_i.size() * m_totali, 0)
{
  init_ptrs();
}

SHOCDataBase::SHOCDataBase(const SHOCDataBase &rhs, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i) :
  shcol(rhs.shcol),
  nlev(rhs.nlev),
  nlevi(rhs.nlevi),
  m_total(rhs.m_total),
  m_totali(rhs.m_totali),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_data(rhs.m_data)
{
  init_ptrs();
}

SHOCDataBase& SHOCDataBase::operator=(const SHOCDataBase& rhs)
{
  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  nlevi    = rhs.nlevi;
  m_total  = rhs.m_total;
  m_totali = rhs.m_totali;
  m_data   = rhs.m_data; // Copy

  init_ptrs();

  return *this;
}

void SHOCDataBase::init_ptrs()
{
  Int offset         = 0;
  Real *data_begin   = m_data.data();

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    *(m_ptrs[i]) = data_begin + offset;
    offset += m_total;
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    *(m_ptrs_i[i]) = data_begin + offset;
    offset += m_totali;
  }
}

void SHOCDataBase::randomize()
{
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    for (size_t j = 0; j < m_total; ++j) {
      (*(m_ptrs[i]))[j] = data_dist(generator);
    }
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    for (size_t j = 0; j < m_totali; ++j) {
      (*(m_ptrs_i[i]))[j] = data_dist(generator);
    }
  }
}

//
// Glue functions to call fortran from from C++ with the Data struct
//
// In all C++ -> Fortran bridge functions you should see shoc_init(nlev, true).
// We are provisionally following P3 here in case SHOC uses global data. The
// 'true' argument is to set shoc to use its fortran implementations instead of
// calling back to C++. We want this behavior since it doesn't make much sense
// for C++ to bridge over to fortran only to have fortran bridge back to C++.
// Anyone who wants the C++ implementation should call it directly. We need
// need to be aware of data layout since f90 is different from cxx. All these
// functions will expect incoming data to be C layout. They will transpose to f90
// before calling fortran and then back to C before returning.
//

void calc_shoc_varorcovar(SHOCVarorcovarData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi,
                         d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<util::TransposeDirection::f2c>();
}

void shoc_grid(SHOCGridData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<util::TransposeDirection::f2c>();
}

void update_host_dse(SHOCEnergydseData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  update_host_dse_c(d.shcol, d.nlev, d.thlm, d.shoc_ql, d.exner, 
                    d.zt_grid, d.phis, d.host_dse);
  d.transpose<util::TransposeDirection::f2c>();
}

void shoc_energy_integrals(SHOCEnergyintData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  shoc_energy_integrals_c(d.shcol, d.nlev, d.host_dse, d.pdel,
                          d.rtm, d.rcm, d.u_wind, d.v_wind,
                          d.se_int, d.ke_int, d.wv_int, d.wl_int);
  d.transpose<util::TransposeDirection::f2c>();
}

void shoc_energy_total_fixer(SHOCEnergytotData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  shoc_energy_total_fixer_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv,
                            d.zt_grid, d.zi_grid,       
                            d.se_b, d.ke_b, d.wv_b, d.wl_b,
                            d.se_a, d.ke_a, d.wv_a, d.wl_a,
                            d.wthl_sfc, d.wqw_sfc, d.rho_zt,
                            d.te_a, d.te_b);
  d.transpose<util::TransposeDirection::f2c>();
}

void calc_shoc_vertflux(SHOCVertfluxData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<util::TransposeDirection::f2c>();
}

void integ_column_stability(SHOCColstabData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  integ_column_stability_c(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
  d.transpose<util::TransposeDirection::f2c>();
}

void compute_shr_prod(SHOCTkeshearData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  compute_shr_prod_c(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind,
                       d.v_wind, d.sterm);
  d.transpose<util::TransposeDirection::f2c>();
}

void isotropic_ts(SHOCIsotropicData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  isotropic_ts_c(d.nlev, d.shcol, d.brunt_int, d.tke, d.a_diss,
                 d.brunt, d.isotropy);
  d.transpose<util::TransposeDirection::f2c>();
}

void adv_sgs_tke(SHOCAdvsgstkeData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  adv_sgs_tke_c(d.nlev, d.shcol, d.dtime, d.shoc_mix, d.wthv_sec,
                d.sterm_zt, d.tk, d.tke, d.a_diss);
  d.transpose<util::TransposeDirection::f2c>();
}

void eddy_diffusivities(SHOCEddydiffData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  eddy_diffusivities_c(d.nlev, d.shcol, d.obklen, d.pblh, d.zt_grid,
     d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
  d.transpose<util::TransposeDirection::f2c>();
}

//
// _f function definitions. These expect data in C layout
//

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = 4;

  Kokkos::Array<view_2d, num_arrays> temp_d;
  Kokkos::Array<size_t, num_arrays> dim1_sizes     = {shcol,  shcol, shcol, shcol};
  Kokkos::Array<size_t, num_arrays> dim2_sizes     = {nlevi,  nlevi, nlev,  nlevi};
  Kokkos::Array<const Real*, num_arrays> ptr_array = {tkh_zi, dz_zi, invar, vertflux};

  // Sync to device
  pack::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

  view_2d
    tkh_zi_d  (temp_d[0]),
    dz_zi_d   (temp_d[1]),
    invar_d   (temp_d[2]),
    vertflux_d(temp_d[3]);

  const Int nk_pack = scream::pack::npack<Spack>(nlev);
  const auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto otkh_zi_d   = util::subview(tkh_zi_d, i);
    const auto odz_zi_d    = util::subview(dz_zi_d, i);
    const auto oinvar_d    = util::subview(invar_d, i);
    const auto overtflux_d = util::subview(vertflux_d, i);

    SHF::calc_shoc_vertflux(team, nlev, otkh_zi_d, odz_zi_d, oinvar_d, overtflux_d);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {vertflux_d};
  pack::device_to_host({vertflux}, {shcol}, {nlevi}, inout_views, true);
}

} // namespace shoc
} // namespace scream
