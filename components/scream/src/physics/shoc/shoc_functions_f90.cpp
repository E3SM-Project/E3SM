#include "shoc_functions_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
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

void shoc_energy_total_fixer_c(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv,
                               Real *zt_grid, Real *zi_grid,
                               Real *se_b, Real *ke_b, Real *wv_b, Real *wl_b,
                               Real *se_a, Real *ke_a, Real *wv_a, Real *wl_a,
                               Real *wthl_sfc, Real *wqw_sfc, Real *rho_zt,
                               Real *te_a, Real *te_b);

void shoc_energy_threshold_fixer_c(Int shcol, Int nlev, Int nlevi,
                             Real *pint, Real *tke, Real *te_a, Real *te_b,
                             Real *se_dis, Int *shoctop);

void shoc_energy_dse_fixer_c(Int shcol, Int nlev,
                             Real *se_dis, Int *shoctop,
                             Real *host_dse);

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

void compute_brunt_shoc_length_c(Int nlev, Int nlevi, Int shcol ,Real *dz_zt,
                                 Real *thv, Real *thv_zi, Real *brunt);

void compute_l_inf_shoc_length_c(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);

void compute_conv_vel_shoc_length_c(Int nlev, Int shcol, Real *pblh, Real *zt_grid,
                                    Real *dz_zt, Real *thv, Real *wthv_sec,
				    Real *conv_vel);

void compute_conv_time_shoc_length_c(Int shcol, Real *pblh, Real *conv_vel,
                                     Real *tscale);

void compute_shoc_mix_shoc_length_c(Int nlev, Int shcol, Real *tke, Real* brunt,
                                    Real *tscale, Real *zt_grid, Real *l_inf,
				    Real *shoc_mix);

void check_length_scale_shoc_length_c(Int nlev, Int shcol, Real *host_dx,
                                    Real *host_dy, Real *shoc_mix);

void shoc_diag_second_moments_srf_c(Int shcol, Real* wthl, Real* uw, Real* vw,
                                   Real* ustar2, Real* wstar);

void linear_interp_c(Real *x1, Real *x2, Real *y1, Real *y2, Int km1,
                     Int km2, Int ncol, Real minthresh);

}

namespace scream {
namespace shoc {

//
// Data struct
//

SHOCDataBase::SHOCDataBase(Int shcol_, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  SHOCDataBase(shcol_, 0, 0, {}, {}, ptrs_c, idx_c)
{}

SHOCDataBase::SHOCDataBase(Int shcol_, Int nlev_, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  SHOCDataBase(shcol_, nlev_, 0, ptrs, {}, ptrs_c, idx_c) {}

SHOCDataBase::SHOCDataBase(Int shcol_, Int nlev_, Int nlevi_,
                           const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i,
			   const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  shcol(shcol_),
  nlev(nlev_),
  nlevi(nlevi_),
  m_total(shcol_ * nlev_),
  m_totali(shcol_ * nlevi_),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_ptrs_c(ptrs_c),
  m_indices_c(idx_c),
  m_data(m_ptrs.size() * m_total + m_ptrs_i.size() * m_totali + m_ptrs_c.size() * shcol, 0),
  m_idx_data(idx_c.size() * shcol, 0)
{
  init_ptrs();
}

SHOCDataBase::SHOCDataBase(const SHOCDataBase &rhs,
                           const std::vector<Real**>& real_ptrs, const std::vector<Int**>& int_ptrs) :
  shcol(rhs.shcol),
  nlev(rhs.nlev),
  nlevi(rhs.nlevi),
  m_total(rhs.m_total),
  m_totali(rhs.m_totali),
  m_ptrs(rhs.m_ptrs.size()),
  m_ptrs_i(rhs.m_ptrs_i.size()),
  m_ptrs_c(rhs.m_ptrs_c.size()),
  m_indices_c(rhs.m_indices_c.size()),
  m_data(rhs.m_data),
  m_idx_data(rhs.m_idx_data)
{
  EKAT_ASSERT(real_ptrs.size() == (m_ptrs.size() + m_ptrs_i.size() + m_ptrs_c.size()));
  EKAT_ASSERT_MSG(int_ptrs.size() == m_indices_c.size(), int_ptrs.size() << " != " << m_indices_c.size());

  size_t real_offset = 0;
  for (auto& item : {&m_ptrs, &m_ptrs_i, &m_ptrs_c}) {
    std::copy(real_ptrs.begin() + real_offset, real_ptrs.begin() + real_offset + item->size(), item->begin());
    real_offset += item->size();
  }

  std::copy(int_ptrs.begin(), int_ptrs.end(), m_indices_c.begin());

  init_ptrs();
}

SHOCDataBase& SHOCDataBase::operator=(const SHOCDataBase& rhs)
{
  shcol      = rhs.shcol;
  nlev       = rhs.nlev;
  nlevi      = rhs.nlevi;
  m_total    = rhs.m_total;
  m_totali   = rhs.m_totali;
  m_data     = rhs.m_data;      // Copy
  m_idx_data = rhs.m_idx_data;  // Copy

  init_ptrs();

  return *this;
}

void SHOCDataBase::init_ptrs()
{
  Int offset       = 0;
  Real *data_begin = m_data.data();

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    *(m_ptrs[i]) = data_begin + offset;
    offset += m_total;
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    *(m_ptrs_i[i]) = data_begin + offset;
    offset += m_totali;
  }

  for (size_t i = 0; i < m_ptrs_c.size(); ++i) {
    *(m_ptrs_c[i]) = data_begin + offset;
    offset += shcol;
  }

  for (size_t i = 0; i < m_indices_c.size(); ++i) {
    *(m_indices_c[i]) = m_idx_data.data() + shcol*i;
  }
}

void SHOCDataBase::randomize(const std::vector<std::pair<Real, Real> >& ranges,
                             const std::vector<std::pair<Real, Real> >& ranges_i,
                             const std::vector<std::pair<Real, Real> >& ranges_c,
                             const std::vector<std::pair<Int, Int> >&   ranges_idx)
{
  std::default_random_engine generator;

  EKAT_ASSERT_MSG(ranges.size() <= m_ptrs.size(), "Provided more ranges than data items");
  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    std::uniform_real_distribution<Real> data_dist(i < ranges.size() ? ranges[i].first  : 0.0,
                                                   i < ranges.size() ? ranges[i].second : 1.0);
    for (int j = 0; j < m_total; ++j) {
      (*(m_ptrs[i]))[j] = data_dist(generator);
    }
  }

  EKAT_ASSERT_MSG(ranges_i.size() <= m_ptrs_i.size(), "Provided more ranges_i than data items");
  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    std::uniform_real_distribution<Real> data_dist(i < ranges_i.size() ? ranges_i[i].first  : 0.0,
                                                   i < ranges_i.size() ? ranges_i[i].second : 1.0);

    for (int j = 0; j < m_totali; ++j) {
      (*(m_ptrs_i[i]))[j] = data_dist(generator);
    }
  }

  EKAT_ASSERT_MSG(ranges_i.size() <= m_ptrs_i.size(), "Provided more ranges_c than data items");
  for (size_t i = 0; i < m_ptrs_c.size(); ++i) {
    std::uniform_real_distribution<Real> data_dist(i < ranges_c.size() ? ranges_c[i].first  : 0.0,
                                                   i < ranges_c.size() ? ranges_c[i].second : 1.0);
    for (int j = 0; j < shcol; ++j) {
      (*(m_ptrs_c[i]))[j] = data_dist(generator);
    }
  }

  EKAT_ASSERT_MSG(ranges_idx.size() == m_indices_c.size(), "Must provide ranges_idx for index data");
  for (size_t i = 0; i < m_indices_c.size(); ++i) {
    std::uniform_int_distribution<Int> data_dist(ranges_idx[i].first, ranges_idx[i].second);
    for (int j = 0; j < shcol; ++j) {
      (*(m_indices_c[i]))[j] = data_dist(generator);
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
  d.transpose<ekat::util::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi,
                         d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_grid(SHOCGridData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void update_host_dse(SHOCEnergydseData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  update_host_dse_c(d.shcol, d.nlev, d.thlm, d.shoc_ql, d.exner,
                    d.zt_grid, d.phis, d.host_dse);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_energy_integrals(SHOCEnergyintData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_energy_integrals_c(d.shcol, d.nlev, d.host_dse, d.pdel,
                          d.rtm, d.rcm, d.u_wind, d.v_wind,
                          d.se_int, d.ke_int, d.wv_int, d.wl_int);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_energy_total_fixer(SHOCEnergytotData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_energy_total_fixer_c(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv,
                            d.zt_grid, d.zi_grid,
                            d.se_b, d.ke_b, d.wv_b, d.wl_b,
                            d.se_a, d.ke_a, d.wv_a, d.wl_a,
                            d.wthl_sfc, d.wqw_sfc, d.rho_zt,
                            d.te_a, d.te_b);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_energy_threshold_fixer(SHOCEnergythreshfixerData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_energy_threshold_fixer_c(d.shcol, d.nlev, d.nlevi,
                          d.pint, d.tke, d.te_a, d.te_b,
			  d.se_dis, d.shoctop);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_energy_dse_fixer(SHOCEnergydsefixerData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_energy_dse_fixer_c(d.shcol, d.nlev,
                          d.se_dis, d.shoctop, d.host_dse);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void calc_shoc_vertflux(SHOCVertfluxData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void integ_column_stability(SHOCColstabData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  integ_column_stability_c(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_shr_prod(SHOCTkeshearData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_shr_prod_c(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind,
                       d.v_wind, d.sterm);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void isotropic_ts(SHOCIsotropicData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  isotropic_ts_c(d.nlev, d.shcol, d.brunt_int, d.tke, d.a_diss,
                 d.brunt, d.isotropy);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void adv_sgs_tke(SHOCAdvsgstkeData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  adv_sgs_tke_c(d.nlev, d.shcol, d.dtime, d.shoc_mix, d.wthv_sec,
                d.sterm_zt, d.tk, d.tke, d.a_diss);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void eddy_diffusivities(SHOCEddydiffData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  eddy_diffusivities_c(d.nlev, d.shcol, d.obklen, d.pblh, d.zt_grid,
     d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_brunt_shoc_length(SHOCBruntlengthData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_brunt_shoc_length_c(d.nlev,d.nlevi,d.shcol,d.dz_zt,d.thv,d.thv_zi,d.brunt);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_l_inf_shoc_length(SHOCInflengthData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_l_inf_shoc_length_c(d.nlev,d.shcol,d.zt_grid,d.dz_zt,d.tke,d.l_inf);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_conv_vel_shoc_length(SHOCConvvelData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_conv_vel_shoc_length_c(d.nlev,d.shcol,d.pblh,d.zt_grid,
                                 d.dz_zt,d.thv,d.wthv_sec,d.conv_vel);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_conv_time_shoc_length(SHOCConvtimeData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_conv_time_shoc_length_c(d.shcol,d.pblh,d.conv_vel,d.tscale);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void compute_shoc_mix_shoc_length(SHOCMixlengthData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  compute_shoc_mix_shoc_length_c(d.nlev,d.shcol,d.tke,d.brunt,d.tscale,
                                 d.zt_grid,d.l_inf,d.shoc_mix);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void check_length_scale_shoc_length(SHOCMixcheckData &d) {
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  check_length_scale_shoc_length_c(d.nlev,d.shcol,d.host_dx,d.host_dy,d.shoc_mix);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void shoc_diag_second_moments_srf(SHOCSecondMomentSrfData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  shoc_diag_second_moments_srf_c(d.shcol, d.wthl, d.uw, d.vw, d.ustar2, d.wstar);
  d.transpose<ekat::util::TransposeDirection::f2c>();
}

void linear_interp(SHOCLinearintData& d)
{
  shoc_init(d.nlev, true);
  d.transpose<ekat::util::TransposeDirection::c2f>();
  linear_interp_c(d.x1,d.x2,d.y1,d.y2,d.nlev,d.nlevi,d.shcol,d.minthresh);
  d.transpose<ekat::util::TransposeDirection::f2c>();
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
  ekat::pack::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d, true);

  view_2d
    tkh_zi_d  (temp_d[0]),
    dz_zi_d   (temp_d[1]),
    invar_d   (temp_d[2]),
    vertflux_d(temp_d[3]);

  const Int nk_pack = ekat::pack::npack<Spack>(nlev);
  const auto policy = ekat::util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const auto otkh_zi_d   = ekat::util::subview(tkh_zi_d, i);
    const auto odz_zi_d    = ekat::util::subview(dz_zi_d, i);
    const auto oinvar_d    = ekat::util::subview(invar_d, i);
    const auto overtflux_d = ekat::util::subview(vertflux_d, i);

    SHF::calc_shoc_vertflux(team, nlev, otkh_zi_d, odz_zi_d, oinvar_d, overtflux_d);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {vertflux_d};
  ekat::pack::device_to_host({vertflux}, {shcol}, {nlevi}, inout_views, true);
}

void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw, Real* ustar2, Real* wstar)
{
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Pack1      = typename ekat::pack::Pack<Real, 1>;
  using view_1d    = typename SHOC::view_1d<Pack1>;

  Kokkos::Array<view_1d, 3> temp_d;
  ekat::pack::host_to_device({wthl, uw, vw}, shcol, temp_d);

  // inputs
  view_1d
    wthl_d (temp_d[0]),
    uw_d   (temp_d[1]),
    vw_d   (temp_d[2]);

  // outputs
  view_1d ustar2_d("ustar2", shcol),
          wstar_d ("wstar", shcol);

  Kokkos::parallel_for("parallel_moments_srf", shcol, KOKKOS_LAMBDA (const int& i) {

     Scalar wthl_s{wthl_d(i)[0]};
     Scalar uw_s{uw_d(i)[0]};
     Scalar vw_s{vw_d(i)[0]};

     Scalar ustar2_s{0};
     Scalar wstar_s{0};

     SHOC::shoc_diag_second_moments_srf(wthl_s, uw_s, vw_s, ustar2_s, wstar_s);

     ustar2_d(i)[0] = ustar2_s;
     wstar_d(i)[0]  = wstar_s;
   });

  Kokkos::Array<view_1d, 2> out_views = {ustar2_d, wstar_d};
  ekat::pack::device_to_host({ustar2, wstar}, shcol, out_views);
}

} // namespace shoc
} // namespace scream
