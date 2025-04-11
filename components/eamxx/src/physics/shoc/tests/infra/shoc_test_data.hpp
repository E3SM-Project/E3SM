#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "share/eamxx_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "shoc_functions.hpp"
#include "physics_constants.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

struct ShocTestGridDataBase : public PhysicsTestData
{
  Real *zt_grid, *zi_grid;

  ShocTestGridDataBase(
    const std::vector<std::vector<Int> >& dims,
    const std::vector<std::vector<Real**> >& reals,
    const std::vector<std::vector<Int**> >& ints = {},
    const std::vector<std::vector<bool**> >& bools = {}) :
    PhysicsTestData(dims, reals, ints, bools)
  {}

  template <typename Engine>
  void randomize(Engine& engine, const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges = {})
  {
    PhysicsTestData::randomize(engine, ranges);

    const auto shcol = dim(zt_grid, 0);
    const auto nlev  = dim(zt_grid, 1);
    const auto nlevi = dim(zi_grid, 1);

    EKAT_REQUIRE_MSG(shcol == dim(zi_grid, 0), "Mismatched shcol dim for zt_grid and zi_grid");
    EKAT_REQUIRE_MSG(nlev == nlevi-1, "Mismatched lev dim for zt_grid and zi_grid");

    // Don't want true randomness in the grid data, need interleaved grid points with some minimum separation
    for (auto i = decltype(shcol){0}; i < shcol; ++i) {
      Real* this_col_zi = zi_grid + nlevi*i;
      std::sort(this_col_zi, this_col_zi + nlevi);
      const auto min = this_col_zi[0];
      const auto max = this_col_zi[nlevi-1];

      const auto avg_jump = (max-min) / (nlevi-1);
      const auto avg_jump_d = avg_jump / 2.5;

      for (auto k = decltype(nlevi){1}; k < (nlevi-1); ++k) {
        std::uniform_real_distribution<Real> x2_dist(min + k*avg_jump - avg_jump_d, min + k*avg_jump + avg_jump_d);
        this_col_zi[k] = x2_dist(engine);
      }

      std::sort(this_col_zi, this_col_zi + nlevi, std::greater<Real>());
      for (auto k = decltype(nlev){0}; k < nlev; ++k) {
        std::uniform_real_distribution<Real> x2_dist(zi_grid[nlevi*i + k], zi_grid[nlevi*i + k+1]);
        zt_grid[nlev*i + k] = x2_dist(engine);
      }
    }
  }
};

struct ShocGridData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *pdel;

  // Outputs
  Real *dz_zt, *dz_zi, *rho_zt;

  ShocGridData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &zt_grid, &pdel, &dz_zt, &rho_zt }, { &zi_grid, &dz_zi }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ShocGridData, 3, shcol, nlev, nlevi);
};

struct ShocDiagObklenData : public PhysicsTestData {
  // Inputs
  Int shcol;
  Real *uw_sfc, *vw_sfc, *wthl_sfc, *wqw_sfc, *thl_sfc, *cldliq_sfc, *qv_sfc;

  // Outputs
  Real *ustar, *kbfs, *obklen;

  ShocDiagObklenData(Int shcol_) :
    PhysicsTestData({{ shcol_ }}, {{ &uw_sfc, &vw_sfc, &wthl_sfc, &wqw_sfc, &thl_sfc, &cldliq_sfc, &qv_sfc, &ustar, &kbfs, &obklen }}), shcol(shcol_) {}

  PTD_STD_DEF(ShocDiagObklenData, 1, shcol);
};

struct UpdateHostDseData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *thlm, *shoc_ql, *inv_exner, *zt_grid, *phis;

  // Outputs
  Real *host_dse;

  UpdateHostDseData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }}, {{ &thlm, &shoc_ql, &inv_exner, &zt_grid, &host_dse }, { &phis }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(UpdateHostDseData, 2, shcol, nlev);
};

struct ShocEnergyFixerData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi, nadv;
  Real dtime;
  Real *se_b, *ke_b, *wv_b, *wl_b, *se_a, *ke_a, *wv_a, *wl_a, *wthl_sfc, *wqw_sfc, *rho_zt, *tke, *pint;

  // Inputs/Outputs
  Real *host_dse;

  ShocEnergyFixerData(Int shcol_, Int nlev_, Int nlevi_, Real dtime_, Int nadv_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_ }}, {{ &zt_grid, &rho_zt, &tke, &host_dse }, { &zi_grid, &pint }, { &se_b, &ke_b, &wv_b, &wl_b, &se_a, &ke_a, &wv_a, &wl_a, &wthl_sfc, &wqw_sfc }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), nadv(nadv_), dtime(dtime_) {}

  PTD_STD_DEF(ShocEnergyFixerData, 5, shcol, nlev, nlevi, dtime, nadv);
};

struct ShocEnergyIntegralsData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *host_dse, *pdel, *rtm, *rcm, *u_wind, *v_wind;

  // Outputs
  Real *se_int, *ke_int, *wv_int, *wl_int;

  ShocEnergyIntegralsData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }}, {{ &host_dse, &pdel, &rtm, &rcm, &u_wind, &v_wind }, { &se_int, &ke_int, &wv_int, &wl_int }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(ShocEnergyIntegralsData, 2, shcol, nlev);
};

struct CalcShocVertfluxData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *tkh_zi, *dz_zi, *invar;

  // Outputs
  Real *vertflux;

  CalcShocVertfluxData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlevi_ }, { shcol_, nlev_ }}, {{ &tkh_zi, &dz_zi, &vertflux }, { &invar }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(CalcShocVertfluxData, 3, shcol, nlev, nlevi);
};

struct CalcShocVarorcovarData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real tunefac;
  Real *isotropy_zi, *tkh_zi, *dz_zi, *invar1, *invar2;

  // Inputs/Outputs
  Real *varorcovar;

  CalcShocVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_) :
    PhysicsTestData({{ shcol_, nlevi_ }, { shcol_, nlev_ }}, {{ &isotropy_zi, &tkh_zi, &dz_zi, &varorcovar }, { &invar1, &invar2 }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), tunefac(tunefac_) {}

  PTD_STD_DEF(CalcShocVarorcovarData, 4, shcol, nlev, nlevi, tunefac);
};

struct ComputeTmpiData : public PhysicsTestData {
  // Inputs
  Int shcol, nlevi;
  Real dtime;
  Real *rho_zi, *dz_zi;

  // Outputs
  Real *tmpi;

  ComputeTmpiData(Int shcol_, Int nlevi_, Real dtime_) :
    PhysicsTestData({{ shcol_, nlevi_ }}, {{ &rho_zi, &dz_zi, &tmpi }}), shcol(shcol_), nlevi(nlevi_), dtime(dtime_) {}

  PTD_STD_DEF(ComputeTmpiData, 3, shcol, nlevi, dtime);
};

struct DpInverseData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *rho_zt, *dz_zt;

  // Outputs
  Real *rdp_zt;

  DpInverseData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &rho_zt, &dz_zt, &rdp_zt }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(DpInverseData, 2, shcol, nlev);
};

struct IntegColumnStabilityData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *dz_zt, *pres, *brunt;

  // Outputs
  Real *brunt_int;

  IntegColumnStabilityData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }}, {{ &dz_zt, &pres, &brunt }, { &brunt_int }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(IntegColumnStabilityData, 2, shcol, nlev);
};

struct CheckTkeData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;

  // Inputs/Outputs
  Real *tke;

  CheckTkeData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &tke }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(CheckTkeData, 2, shcol, nlev);
};

struct ShocTkeData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real dtime;
  Real *wthv_sec, *shoc_mix, *dz_zi, *dz_zt, *pres, *tabs, *u_wind, *v_wind, *brunt, *pblh;

  // Inputs/Outputs
  Real *tke, *tk, *tkh;

  // Outputs
  Real *isotropy;

  ShocTkeData(Int shcol_, Int nlev_, Int nlevi_, Real dtime_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_ }},
                         {{ &wthv_sec, &shoc_mix, &dz_zt, &pres, &tabs, &u_wind, &v_wind, &brunt, &zt_grid, &tke, &tk, &tkh, &isotropy },
                          { &dz_zi, &zi_grid },
                          { &pblh }}),
                        shcol(shcol_), nlev(nlev_), nlevi(nlevi_), dtime(dtime_) {}

  PTD_STD_DEF(ShocTkeData, 4, shcol, nlev, nlevi, dtime);
};

struct ComputeShrProdData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *dz_zi, *u_wind, *v_wind;

  // Outputs
  Real *sterm;

  ComputeShrProdData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlevi_ }, { shcol_, nlev_ }}, {{ &dz_zi, &sterm }, { &u_wind, &v_wind }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ComputeShrProdData, 3, shcol, nlev, nlevi);
};

struct IsotropicTsData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *brunt_int, *tke, *a_diss, *brunt;

  // Outputs
  Real *isotropy;

  IsotropicTsData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_ }, { shcol_, nlev_ }}, {{ &brunt_int }, { &tke, &a_diss, &brunt, &isotropy }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(IsotropicTsData, 2, shcol, nlev);
};

struct AdvSgsTkeData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real dtime;
  Real *shoc_mix, *wthv_sec, *sterm_zt, *tk;

  // Inputs/Outputs
  Real *tke;

  // Outputs
  Real *a_diss;

  AdvSgsTkeData(Int shcol_, Int nlev_, Real dtime_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss }}), shcol(shcol_), nlev(nlev_), dtime(dtime_) {}

  PTD_STD_DEF(AdvSgsTkeData, 3, shcol, nlev, dtime);
};

struct EddyDiffusivitiesData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *pblh, *zt_grid, *tabs, *shoc_mix, *sterm_zt, *isotropy, *tke;

  // Outputs
  Real *tkh, *tk;

  EddyDiffusivitiesData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_ }, { shcol_, nlev_ }}, {{ &pblh }, { &zt_grid, &tabs, &shoc_mix, &sterm_zt, &isotropy, &tke, &tkh, &tk }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(EddyDiffusivitiesData, 2, shcol, nlev);
};

struct ShocLengthData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *host_dx, *host_dy, *tke, *dz_zt, *thv;

  // Outputs
  Real *brunt, *shoc_mix;

  ShocLengthData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_ }, { shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &host_dx, &host_dy }, { &zt_grid, &dz_zt, &tke, &thv, &brunt, &shoc_mix }, { &zi_grid }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ShocLengthData, 3, shcol, nlev, nlevi);
};

struct ComputeBruntShocLengthData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *dz_zt, *thv, *thv_zi;

  // Outputs
  Real *brunt;

  ComputeBruntShocLengthData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &dz_zt, &thv, &brunt }, { &thv_zi }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ComputeBruntShocLengthData, 3, shcol, nlev, nlevi);
};

struct ComputeLInfShocLengthData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *zt_grid, *dz_zt, *tke;

  // Inputs/Outputs
  Real *l_inf;

  ComputeLInfShocLengthData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }}, {{ &zt_grid, &dz_zt, &tke }, { &l_inf }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(ComputeLInfShocLengthData, 2, shcol, nlev);
};

struct ComputeConvTimeShocLengthData : public PhysicsTestData {
  // Inputs
  Int shcol;
  Real *pblh;

  // Inputs/Outputs
  Real *conv_vel, *tscale;

  ComputeConvTimeShocLengthData(Int shcol_) :
    PhysicsTestData({{ shcol_ }}, {{ &pblh, &conv_vel, &tscale }}), shcol(shcol_) {}

  PTD_STD_DEF(ComputeConvTimeShocLengthData, 1, shcol);
};

struct ComputeShocMixShocLengthData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *tke, *brunt, *tscale, *zt_grid, *l_inf;

  // Outputs
  Real *shoc_mix;

  ComputeShocMixShocLengthData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }}, {{ &tke, &brunt, &zt_grid, &shoc_mix }, { &tscale, &l_inf }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(ComputeShocMixShocLengthData, 2, shcol, nlev);
};

struct CheckLengthScaleShocLengthData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *host_dx, *host_dy;

  // Inputs/Outputs
  Real *shoc_mix;

  CheckLengthScaleShocLengthData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_ }, { shcol_, nlev_ }}, {{ &host_dx, &host_dy }, { &shoc_mix }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(CheckLengthScaleShocLengthData, 2, shcol, nlev);
};

struct ClippingDiagThirdShocMomentsData : public PhysicsTestData {
  // Inputs
  Int shcol, nlevi;
  Real *w_sec_zi;

  // Inputs/Outputs
  Real *w3;

  ClippingDiagThirdShocMomentsData(Int shcol_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlevi_ }}, {{ &w_sec_zi, &w3 }}), shcol(shcol_), nlevi(nlevi_) {}

  PTD_STD_DEF(ClippingDiagThirdShocMomentsData, 2, shcol, nlevi);
};

struct DiagSecondMomentsSrfData : public PhysicsTestData {
  // Inputs
  Int shcol;
  Real *wthl_sfc, *uw_sfc, *vw_sfc;

  // Outputs
  Real *ustar2, *wstar;

  DiagSecondMomentsSrfData(Int shcol_) :
    PhysicsTestData({{ shcol_ }}, {{ &wthl_sfc, &uw_sfc, &vw_sfc, &ustar2, &wstar }}), shcol(shcol_) {}

  PTD_STD_DEF(DiagSecondMomentsSrfData, 1, shcol);
};

struct LinearInterpData : public PhysicsTestData {
  // Inputs
  Real *x1, *x2, *y1;
  Int ncol, km1, km2;
  Real minthresh;

  // Outputs
  Real *y2;

  LinearInterpData(Int ncol_, Int km1_, Int km2_, Real minthresh_) :
    PhysicsTestData({{ ncol_, km1_ }, { ncol_, km2_ }}, {{ &x1, &y1 }, { &x2, &y2 }}), ncol(ncol_), km1(km1_), km2(km2_), minthresh(minthresh_) {}

  PTD_STD_DEF(LinearInterpData, 4, ncol, km1, km2, minthresh);
};

struct DiagThirdShocMomentsData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *w_sec, *thl_sec, *wthl_sec, *isotropy, *brunt, *thetal, *tke, *dz_zt, *dz_zi;

  // Outputs
  Real *w3;

  DiagThirdShocMomentsData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &w_sec, &isotropy, &brunt, &thetal, &tke, &dz_zt, &zt_grid }, { &thl_sec, &wthl_sec, &dz_zi, &zi_grid, &w3 }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(DiagThirdShocMomentsData, 3, shcol, nlev, nlevi);
};

struct ComputeDiagThirdShocMomentData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *w_sec, *thl_sec, *wthl_sec, *tke, *dz_zt, *dz_zi, *isotropy_zi, *brunt_zi, *w_sec_zi, *thetal_zi;

  // Outputs
  Real *w3;

  ComputeDiagThirdShocMomentData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &w_sec, &tke, &dz_zt }, { &thl_sec, &wthl_sec, &dz_zi, &isotropy_zi, &brunt_zi, &w_sec_zi, &thetal_zi, &w3 }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ComputeDiagThirdShocMomentData, 3, shcol, nlev, nlevi);
};

struct ShocAssumedPdfData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *thetal, *qw, *w_field, *thl_sec, *qw_sec, *wthl_sec, *w_sec, *wqw_sec, *qwthl_sec, *w3, *pres;

  // Outputs
  Real *shoc_cldfrac, *shoc_ql, *wqls, *wthv_sec, *shoc_ql2;

  ShocAssumedPdfData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &thetal, &qw, &w_field, &w_sec, &pres, &zt_grid, &shoc_cldfrac, &shoc_ql, &wqls, &wthv_sec, &shoc_ql2 }, { &thl_sec, &qw_sec, &wthl_sec, &wqw_sec, &qwthl_sec, &w3, &zi_grid }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(ShocAssumedPdfData, 3, shcol, nlev, nlevi);
};

struct ShocAssumedPdfTildeToRealData {
  // Inputs
  Real w_first, sqrtw2;

  // Inputs/Outputs
  Real w1;
};

struct ShocAssumedPdfVvParametersData {
  // Inputs
  Real w_first, w_sec, w3var, w_tol_sqd;

  // Outputs
  Real skew_w, w1_1, w1_2, w2_1, w2_2, a;
};

struct ShocAssumedPdfThlParametersData {
  // Inputs
  Real wthlsec, sqrtw2, sqrtthl, thlsec, thl_first, w1_1, w1_2, skew_w, a, thl_tol, w_thresh;

  // Outputs
  Real thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2;
};

struct ShocAssumedPdfQwParametersData {
  // Inputs
  Real wqwsec, sqrtw2, skew_w, sqrtqt, qwsec, w1_2, w1_1, qw_first, a, rt_tol, w_thresh;

  // Outputs
  Real qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2;
};

struct ShocAssumedPdfInplumeCorrelationsData {
  // Inputs
  Real sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2, sqrtthl2_2, qwthlsec, qw1_1, qw_first, thl1_1, thl_first, qw1_2, thl1_2;

  // Outputs
  Real r_qwthl_1;
};

struct ShocAssumedPdfComputeTemperatureData {
  // Inputs
  Real thl1, pval;

  // Outputs
  Real tl1;
};

struct ShocAssumedPdfComputeQsData {
  // Inputs
  Real tl1_1, tl1_2, pval;

  // Outputs
  Real qs1, beta1, qs2, beta2;
};

struct ShocAssumedPdfComputeSData {
  // Inputs
  Real qw1, qs1, beta, pval, thl2, qw2, sqrtthl2, sqrtqw2, r_qwthl;

  // Outputs
  Real s, std_s, qn, c;
};

struct ShocAssumedPdfComputeSgsLiquidData {
  // Inputs
  Real a, ql1, ql2;

  // Outputs
  Real shoc_ql;
};

struct ShocAssumedPdfComputeCloudLiquidVarianceData {
  // Inputs
  Real a, s1, ql1, c1, std_s1, s2, ql2, c2, std_s2, shoc_ql;

  // Outputs
  Real shoc_ql2;
};

struct ShocAssumedPdfComputeLiquidWaterFluxData {
  // Inputs
  Real a, w1_1, w_first, ql1, w1_2, ql2;

  // Outputs
  Real wqls;
};

struct ShocAssumedPdfComputeBuoyancyFluxData {
  // Inputs
  Real wthlsec, wqwsec, pval, wqls;

  // Outputs
  Real wthv_sec;
};

struct DiagSecondMomentsUbycondData : public PhysicsTestData {
  // Inputs
  Int shcol;

  // Outputs
  Real *thl_sec, *qw_sec, *wthl_sec, *wqw_sec, *qwthl_sec, *uw_sec, *vw_sec, *wtke_sec;

  DiagSecondMomentsUbycondData(Int shcol_) :
    PhysicsTestData({{ shcol_ }}, {{ &thl_sec, &qw_sec, &wthl_sec, &wqw_sec, &qwthl_sec, &uw_sec, &vw_sec, &wtke_sec }}), shcol(shcol_) {}

  PTD_STD_DEF(DiagSecondMomentsUbycondData, 1, shcol);
};

struct PblintdInitPotData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *thl, *ql, *q;

  // Outputs
  Real *thv;

  PblintdInitPotData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &thl, &ql, &q, &thv }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(PblintdInitPotData, 2, shcol, nlev);
};

struct PblintdCldcheckData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *zi, *cldn;

  // Inputs/Outputs
  Real *pblh;

  PblintdCldcheckData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlevi_ }, { shcol_, nlev_ }, { shcol_ }}, {{ &zi }, { &cldn }, { &pblh }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(PblintdCldcheckData, 3, shcol, nlev, nlevi);
};

struct DiagSecondMomentsLbycondData : public PhysicsTestData {
  // Inputs
  Int shcol;
  Real *wthl_sfc, *wqw_sfc, *uw_sfc, *vw_sfc, *ustar2, *wstar;

  // Outputs
  Real *wthl_sec, *wqw_sec, *uw_sec, *vw_sec, *wtke_sec, *thl_sec, *qw_sec, *qwthl_sec;

  DiagSecondMomentsLbycondData(Int shcol_) :
    PhysicsTestData({{ shcol_ }}, {{ &wthl_sfc, &wqw_sfc, &uw_sfc, &vw_sfc, &ustar2, &wstar, &wthl_sec, &wqw_sec, &uw_sec, &vw_sec, &wtke_sec, &thl_sec, &qw_sec, &qwthl_sec }}), shcol(shcol_) {}

  PTD_STD_DEF(DiagSecondMomentsLbycondData, 1, shcol);
};

struct DiagSecondMomentsData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *thetal, *qw, *u_wind, *v_wind, *tke, *isotropy, *tkh, *tk, *dz_zi, *shoc_mix;

  // Inputs/Outputs
  Real *thl_sec, *qw_sec, *wthl_sec, *wqw_sec, *qwthl_sec, *uw_sec, *vw_sec, *wtke_sec;

  // Outputs
  Real *w_sec;

  DiagSecondMomentsData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }}, {{ &thetal, &qw, &u_wind, &v_wind, &tke, &isotropy, &tkh, &tk, &zt_grid, &shoc_mix, &w_sec }, { &dz_zi, &zi_grid, &thl_sec, &qw_sec, &wthl_sec, &wqw_sec, &qwthl_sec, &uw_sec, &vw_sec, &wtke_sec }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(DiagSecondMomentsData, 3, shcol, nlev, nlevi);
};

struct DiagSecondShocMomentsData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *thetal, *qw, *u_wind, *v_wind, *tke, *isotropy, *tkh, *tk, *dz_zi, *shoc_mix, *wthl_sfc, *wqw_sfc, *uw_sfc, *vw_sfc;

  // Outputs
  Real *thl_sec, *qw_sec, *wthl_sec, *wqw_sec, *qwthl_sec, *uw_sec, *vw_sec, *wtke_sec, *w_sec;

  DiagSecondShocMomentsData(Int shcol_, Int nlev_, Int nlevi_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_ }}, {{ &thetal, &qw, &u_wind, &v_wind, &tke, &isotropy, &tkh, &tk, &zt_grid, &shoc_mix, &w_sec }, { &dz_zi, &zi_grid, &thl_sec, &qw_sec, &wthl_sec, &wqw_sec, &qwthl_sec, &uw_sec, &vw_sec, &wtke_sec }, { &wthl_sfc, &wqw_sfc, &uw_sfc, &vw_sfc }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(DiagSecondShocMomentsData, 3, shcol, nlev, nlevi);
};

struct ComputeShocVaporData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *qw, *ql;

  // Outputs
  Real *qv;

  ComputeShocVaporData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &qw, &ql, &qv }}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(ComputeShocVaporData, 2, shcol, nlev);
};

struct UpdatePrognosticsImplicitData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi, num_tracer;
  Real dtime;
  Real *dz_zt, *dz_zi, *rho_zt, *tk, *tkh, *uw_sfc, *vw_sfc, *wthl_sfc, *wqw_sfc, *wtracer_sfc;

  // Inputs/Outputs
  Real *thetal, *qw, *tracer, *tke, *u_wind, *v_wind;

  UpdatePrognosticsImplicitData(Int shcol_, Int nlev_, Int nlevi_, Int num_tracer_, Real dtime_) :
    ShocTestGridDataBase({{ shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_ }, { shcol_, num_tracer_ }, { shcol_, nlev_, num_tracer_ }}, {{ &dz_zt, &rho_zt, &zt_grid, &tk, &tkh, &thetal, &qw, &tke, &u_wind, &v_wind }, { &dz_zi, &zi_grid }, { &uw_sfc, &vw_sfc, &wthl_sfc, &wqw_sfc }, { &wtracer_sfc }, { &tracer }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), num_tracer(num_tracer_), dtime(dtime_) {}

  PTD_STD_DEF(UpdatePrognosticsImplicitData, 5, shcol, nlev, nlevi, num_tracer, dtime);
};

struct ShocMainData : public ShocTestGridDataBase {
  // Inputs
  Int shcol, nlev, nlevi, nadv, num_qtracers;
  Real dtime;
  Real *host_dx, *host_dy, *thv, *pres, *presi, *pdel, *wthl_sfc, *wqw_sfc, *uw_sfc, *vw_sfc, *wtracer_sfc, *w_field, *inv_exner, *phis;

  // Inputs for shoc_init
  Int nbot_shoc, ntop_shoc;
  Real *pref_mid;

  // Inputs/Outputs
  Real *host_dse, *tke, *thetal, *qw, *u_wind, *v_wind, *qtracers, *wthv_sec, *tkh, *tk, *shoc_ql, *shoc_cldfrac;

  // Outputs
  Real *pblh, *shoc_mix, *isotropy, *w_sec, *thl_sec, *qw_sec, *qwthl_sec, *wthl_sec, *wqw_sec, *wtke_sec, *uw_sec, *vw_sec, *w3, *wqls_sec, *brunt, *shoc_ql2;

  Real elapsed_s;

  ShocMainData(Int shcol_, Int nlev_, Int nlevi_, Int num_qtracers_, Real dtime_, Int nadv_, Int nbot_shoc_, Int ntop_shoc_) :
    ShocTestGridDataBase({{ shcol_ }, { shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_, num_qtracers_ }, { shcol_, nlev_, num_qtracers_ }, { nlev_ }},
                    {{ &host_dx, &host_dy, &wthl_sfc, &wqw_sfc, &uw_sfc, &vw_sfc, &phis, &pblh },
                     { &thv, &zt_grid, &pres, &pdel, &w_field, &inv_exner, &host_dse, &tke, &thetal,
                       &qw, &u_wind, &v_wind, &wthv_sec, &tkh, &tk, &shoc_ql, &shoc_cldfrac,
                       &shoc_mix, &isotropy, &w_sec, &wqls_sec, &brunt, &shoc_ql2 },
                     { &zi_grid, &presi, &thl_sec, &qw_sec, &qwthl_sec, &wthl_sec,
                       &wqw_sec, &wtke_sec, &uw_sec, &vw_sec, &w3 },
                     { &wtracer_sfc },
                     { &qtracers },
                     { &pref_mid }}),
                    shcol(shcol_), nlev(nlev_), nlevi(nlevi_), nadv(nadv_),
                    num_qtracers(num_qtracers_), dtime(dtime_), nbot_shoc(nbot_shoc_), ntop_shoc(ntop_shoc_) {}

  PTD_STD_DEF(ShocMainData, 8, shcol, nlev, nlevi, num_qtracers, dtime, nadv, nbot_shoc, ntop_shoc);

  template <size_t N>
  Real interpolate_data(const std::array<Real, N>& ref_elevations,
                        const std::array<Real, N>& ref_data,
                        Real z)
  {
    auto pos = std::lower_bound(ref_elevations.begin(), ref_elevations.end(), z);
    Int index = pos - ref_elevations.begin();
    if (index == 0)
      return ref_data[0];
    else if (index < (Int)N) {
      const Real a = (z - ref_elevations[index-1]) /
        (ref_elevations[index] - ref_elevations[index-1]);
      return (1.0 - a) * ref_data[index-1] + a * ref_data[index];
    }
    else {
      // Don't extrapolate off the end of the table.
      return ref_data[N-1];
    }
  }

  void compute_column_pressure(Int shcol, Int nlev, const Real* z,
                               Real* pres) {
    using consts = scream::physics::Constants<Real>;
    const Real k = consts::Rair / consts::Cpair;
    const Real c = -consts::gravit * pow(consts::P0, k) / consts::Rair;
    const Real p_s = 1015e2;

    const std::array<Real, 5> z_ref = {0.0, 520.0, 1480.0, 2000.0, 3000.0};
    const std::array<Real, 5> theta_ref = {299.7, 298.7, 302.4, 308.2, 312.85};

    // Move up the column, computing the pressures at each elevation.
    for (Int i = 0; i < shcol; ++i) {
      const Int offset = nlev*i;
      for (Int j = 0; j < nlev; ++j) {
        Real z0 = (j == 0) ? 0.0 : z[offset + j-1];
        Real z1 = z[offset+ j];
        Real th0 = interpolate_data(z_ref, theta_ref, z0);
        Real th1 = interpolate_data(z_ref, theta_ref, z1);
        Real p0 = (j == 0) ? p_s : pres[offset + j-1];
        if (std::abs(th0 - th1) < 1e-14 * th0) {
          pres[offset + j] = pow(pow(p0, k) + k*c*(z1 - z0)/th0, 1.0/k);
        }
        else {
          Real ra = (z1 - z0)/(th1 - th0);
          pres[offset + j] = pow(pow(p0, k) + k*c*ra*log(th1/th0), 1.0/k);
        }
      }
    }
  }

  template <typename Engine>
  void randomize(Engine& engine, const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges = {})
  {
    using consts = scream::physics::Constants<Real>;

    ShocTestGridDataBase::randomize(engine, ranges);

    const auto shcol = dim(zt_grid, 0);
    const auto nlev  = dim(zt_grid, 1);
    const auto nlevi = dim(zi_grid, 1);

    EKAT_REQUIRE_MSG(shcol == dim(zi_grid, 0), "Mismatched shcol dim for zt_grid and zi_grid");
    EKAT_REQUIRE_MSG(nlev == nlevi-1, "Mismatched lev dim for zt_grid and zi_grid");

    compute_column_pressure(shcol, nlev, zt_grid, pres);
    compute_column_pressure(shcol, nlevi, zi_grid, presi);

    const Real pot_temp = 300; // keep pot_temp fixed for now

    for (auto i = decltype(shcol){0}; i < shcol; ++i) {
      const auto nlev_offset = i * nlev;
      const auto nlevi_offset = i * nlevi;
      for (auto k = decltype(nlev){0}; k < nlev; ++k) {
        pdel[nlev_offset + k] = std::abs(presi[nlevi_offset + k] - presi[nlevi_offset + k+1]);
        inv_exner[nlev_offset + k] = pow(pres[nlev_offset + k]/consts::P0, consts::Rair/consts::Cpair);
        host_dse[nlev_offset + k] = consts::Cpair * inv_exner[nlev_offset + k] * thv[nlev_offset + k] +
          consts::gravit * zt_grid[nlev_offset + k];

        const Real qv = qw[nlev_offset+k] - shoc_ql[nlev_offset+k];
        thetal[nlev_offset+k] = pot_temp - (consts::LatVap/consts::Cpair)*shoc_ql[nlev_offset+k];
        thv[nlev_offset+k] = pot_temp * (1 + 0.61*qv - shoc_ql[nlev_offset+k]);
        inv_exner[nlev_offset+k] = 1/std::pow(pres[nlev_offset+k]/consts::P0,consts::Rair/consts::Cpair);
      }
    }

    // 3 types of pref_mid ranges, pick one randomly
    std::uniform_int_distribution<Int> int_dist(0, 2);
    const auto pref_type = int_dist(engine);
    std::pair<Real,Real> pref_mid_range;
    if (pref_type == 0) {
      // all(pref_mid) >= pblmaxp
      pref_mid_range.first  = 1e5;
      pref_mid_range.second = 8e5;
    }
    else if (pref_type == 1) {
      // all(pref_mid) < pblmaxp
      pref_mid_range.first  = 1e3;
      pref_mid_range.second = 8e3;
    }
    else {
      // both pref_mid >= pblmaxp and pref_mid < pblmaxp values
      pref_mid_range.first  = 1e4;
      pref_mid_range.second = 8e4;
    }

    // pref_mid must be monotonically increasing
    Real upper = pref_mid_range.first;
    Real lower = pref_mid_range.second;
    for (Int k=0; k<nlev; ++k) {
      pref_mid[k] = upper - k*(upper-lower)/(nlev-1);
    }

    // TKH and TK get the same values on purpose
    std::copy(tkh, tkh + shcol*nlev, tk);
  }
};

struct PblintdHeightData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, npbl;
  Real *z, *u, *v, *ustar, *thv, *thv_ref;

  // Inputs/Outputs
  Real *rino;
  bool *check;

  // Outputs
  Real *pblh;

  PblintdHeightData(Int shcol_, Int nlev_, Int npbl_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }, { shcol_ }}, {{ &z, &u, &v, &thv, &rino }, { &ustar, &thv_ref, &pblh }}, {}, {{ &check }}), shcol(shcol_), nlev(nlev_), npbl(npbl_) {}

  PTD_STD_DEF(PblintdHeightData, 3, shcol, nlev, npbl);
};

struct VdShocDecompandSolveData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi, n_rhs;
  Real dtime;
  Real *kv_term, *tmpi, *rdp_zt, *flux;

  // Inputs/Outputs
  Real *var;

  VdShocDecompandSolveData(Int shcol_, Int nlev_, Int nlevi_, Real dtime_, Int n_rhs_) :
    PhysicsTestData({{shcol_}, {shcol_, nlev_},               {shcol_, nlevi_},  {shcol_, nlev_, n_rhs_}},
                    {{&flux}, {&rdp_zt}, {&kv_term, &tmpi}, {&var }                }, {}),
                    shcol(shcol_), nlev(nlev_), nlevi(nlevi_), n_rhs(n_rhs_), dtime(dtime_) {}

  PTD_STD_DEF(VdShocDecompandSolveData, 5, shcol, nlev, nlevi, dtime, n_rhs);
};

struct PblintdSurfTempData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *z, *ustar, *obklen, *kbfs, *thv;

  // Inputs/Outputs
  Real *pblh, *rino;
  bool *check;

  // Outputs
  Real *tlv;

  PblintdSurfTempData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }, { shcol_ }}, {{ &z, &thv, &rino }, { &ustar, &obklen, &kbfs, &tlv, &pblh }}, {}, {{ &check }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(PblintdSurfTempData, 3, shcol, nlev, nlevi);
};

struct PblintdCheckPblhData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi;
  Real *z, *ustar;
  bool *check;

  // Outputs
  Real *pblh;

  PblintdCheckPblhData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_ }, { shcol_ }}, {{ &z }, { &ustar, &pblh }}, {}, {{ &check }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_) {}

  PTD_STD_DEF(PblintdCheckPblhData, 3, shcol, nlev, nlevi);
};

struct PblintdData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev, nlevi, npbl;
  Real *z, *zi, *thl, *ql, *q, *u, *v, *ustar, *obklen, *kbfs, *cldn;

  // Outputs
  Real *pblh;

  PblintdData(Int shcol_, Int nlev_, Int nlevi_, Int npbl_) :
    PhysicsTestData({{ shcol_, nlev_ }, { shcol_, nlevi_ }, { shcol_ }}, {{ &z, &thl, &ql, &q, &u, &v, &cldn }, { &zi }, { &ustar, &obklen, &kbfs, &pblh }}), shcol(shcol_), nlev(nlev_), nlevi(nlevi_), npbl(npbl_) {}

  PTD_STD_DEF(PblintdData, 4, shcol, nlev, nlevi, npbl);
};

struct ComputeShocTempData : public PhysicsTestData {
  // Inputs
  Int shcol, nlev;
  Real *thetal, *ql, *inv_exner;

  // Outputs
  Real *tabs;

  ComputeShocTempData(Int shcol_, Int nlev_) :
    PhysicsTestData({{ shcol_, nlev_ }}, {{ &thetal, &ql, &inv_exner, &tabs}}), shcol(shcol_), nlev(nlev_) {}

  PTD_STD_DEF(ComputeShocTempData, 2, shcol, nlev);
};

// Glue functions to call from host with the Data struct

void shoc_grid                                      (ShocGridData& d);
void shoc_diag_obklen                               (ShocDiagObklenData& d);
void update_host_dse                                (UpdateHostDseData& d);
void shoc_energy_fixer                              (ShocEnergyFixerData& d);
void shoc_energy_integrals                          (ShocEnergyIntegralsData& d);
void calc_shoc_vertflux                             (CalcShocVertfluxData& d);
void calc_shoc_varorcovar                           (CalcShocVarorcovarData& d);
void compute_tmpi                                   (ComputeTmpiData& d);
void dp_inverse                                     (DpInverseData& d);
void integ_column_stability                         (IntegColumnStabilityData& d);
void check_tke                                      (CheckTkeData& d);
void shoc_tke                                       (ShocTkeData& d);
void compute_shr_prod                               (ComputeShrProdData& d);
void isotropic_ts                                   (IsotropicTsData& d);
void adv_sgs_tke                                    (AdvSgsTkeData& d);
void eddy_diffusivities                             (EddyDiffusivitiesData& d);
void shoc_length                                    (ShocLengthData& d);
void compute_brunt_shoc_length                      (ComputeBruntShocLengthData& d);
void compute_l_inf_shoc_length                      (ComputeLInfShocLengthData& d);
void compute_shoc_mix_shoc_length                   (ComputeShocMixShocLengthData& d);
void check_length_scale_shoc_length                 (CheckLengthScaleShocLengthData& d);
void clipping_diag_third_shoc_moments               (ClippingDiagThirdShocMomentsData& d);
void diag_second_moments_srf                        (DiagSecondMomentsSrfData& d);
void linear_interp                                  (LinearInterpData& d);
void diag_third_shoc_moments                        (DiagThirdShocMomentsData& d);
void compute_diag_third_shoc_moment                 (ComputeDiagThirdShocMomentData& d);
void shoc_assumed_pdf                               (ShocAssumedPdfData& d);
void shoc_assumed_pdf_tilde_to_real                 (ShocAssumedPdfTildeToRealData& d);
void shoc_assumed_pdf_vv_parameters                 (ShocAssumedPdfVvParametersData& d);
void shoc_assumed_pdf_thl_parameters                (ShocAssumedPdfThlParametersData& d);
void shoc_assumed_pdf_qw_parameters                 (ShocAssumedPdfQwParametersData& d);
void shoc_assumed_pdf_inplume_correlations          (ShocAssumedPdfInplumeCorrelationsData& d);
void shoc_assumed_pdf_compute_temperature           (ShocAssumedPdfComputeTemperatureData& d);
void shoc_assumed_pdf_compute_qs                    (ShocAssumedPdfComputeQsData& d);
void shoc_assumed_pdf_compute_s                     (ShocAssumedPdfComputeSData& d);
void shoc_assumed_pdf_compute_sgs_liquid            (ShocAssumedPdfComputeSgsLiquidData& d);
void shoc_assumed_pdf_compute_cloud_liquid_variance (ShocAssumedPdfComputeCloudLiquidVarianceData& d);
void shoc_assumed_pdf_compute_liquid_water_flux     (ShocAssumedPdfComputeLiquidWaterFluxData& d);
void shoc_assumed_pdf_compute_buoyancy_flux         (ShocAssumedPdfComputeBuoyancyFluxData& d);
void diag_second_moments_ubycond                    (DiagSecondMomentsUbycondData& d);
void pblintd_init_pot                               (PblintdInitPotData& d);
void pblintd_cldcheck                               (PblintdCldcheckData& d);
void diag_second_moments_lbycond                    (DiagSecondMomentsLbycondData& d);
void diag_second_moments                            (DiagSecondMomentsData& d);
void diag_second_shoc_moments                       (DiagSecondShocMomentsData& d);
void compute_shoc_vapor                             (ComputeShocVaporData& d);
void update_prognostics_implicit                    (UpdatePrognosticsImplicitData& d);
void shoc_main                                      (ShocMainData& d);
void pblintd_height                                 (PblintdHeightData& d);
void vd_shoc_decomp_and_solve                       (VdShocDecompandSolveData& d);
void pblintd_surf_temp(PblintdSurfTempData& d);
void pblintd_check_pblh(PblintdCheckPblhData& d);
void pblintd(PblintdData& d);
void compute_shoc_temperature(ComputeShocTempData& d);

// Call from host

void calc_shoc_varorcovar_host(Int shcol, Int nlev, Int nlevi, Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
                            Real *invar1, Real *invar2, Real *varorcovar);
void calc_shoc_vertflux_host(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
void shoc_diag_second_moments_srf_host(Int shcol, Real* wthl, Real* uw, Real* vw,
                          Real* ustar2, Real* wstar);
void shoc_diag_second_moments_ubycond_host(Int shcol, Real* thl, Real* qw, Real* wthl,
                          Real* wqw, Real* qwthl, Real* uw, Real* vw, Real* wtke);
void update_host_dse_host(Int shcol, Int nlev, Real* thlm, Real* shoc_ql, Real* inv_exner, Real* zt_grid,
                       Real* phis, Real* host_dse);
void compute_diag_third_shoc_moment_host(Int shcol, Int nlev, Int nlevi, Real* w_sec,
                                      Real* thl_sec, Real* wthl_sec, Real* tke,
                                      Real* dz_zt, Real* dz_zi, Real* isotropy_zi,
                                      Real* brunt_zi, Real* w_sec_zi, Real* thetal_zi,
                                      Real* w3);
void shoc_pblintd_init_pot_host(Int shcol, Int nlev, Real* thl, Real* ql, Real* q, Real* thv);
void compute_shoc_mix_shoc_length_host(Int nlev, Int shcol, Real* tke, Real* brunt,
                                    Real* zt_grid, Real* l_inf, Real* shoc_mix);
void check_tke_host(Int shcol, Int nlev, Real* tke);
void linear_interp_host(Real* x1, Real* x2, Real* y1, Real* y2, Int km1, Int km2, Int ncol, Real minthresh);
void clipping_diag_third_shoc_moments_host(Int nlevi, Int shcol, Real *w_sec_zi,
                                        Real *w3);
void shoc_energy_integrals_host(Int shcol, Int nlev, Real *host_dse, Real *pdel,
                             Real *rtm, Real *rcm, Real *u_wind, Real *v_wind,
                             Real *se_int, Real *ke_int, Real *wv_int, Real *wl_int);
void compute_brunt_shoc_length_host(Int nlev, Int nlevi, Int shcol, Real* dz_zt, Real* thv,
                                 Real* thv_zi, Real* brunt);
void compute_l_inf_shoc_length_host(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);
void check_length_scale_shoc_length_host(Int nlev, Int shcol, Real* host_dx, Real* host_dy,
                                      Real* shoc_mix);
void diag_second_moments_lbycond_host(Int shcol, Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* ustar2, Real* wstar,
                                  Real* wthl_sec, Real* wqw_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* thl_sec,
                                  Real* qw_sec, Real* qwthl_sec);
void diag_second_moments_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke, Real* isotropy,
                          Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix, Real* thl_sec, Real* qw_sec,
                          Real* wthl_sec, Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec);
void diag_second_shoc_moments_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* u_wind, Real* v_wind, Real* tke,
                                Real* isotropy, Real* tkh, Real* tk, Real* dz_zi, Real* zt_grid, Real* zi_grid, Real* shoc_mix,
                                Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc, Real* thl_sec, Real* qw_sec, Real* wthl_sec,
                                Real* wqw_sec, Real* qwthl_sec, Real* uw_sec, Real* vw_sec, Real* wtke_sec, Real* w_sec);
void shoc_diag_obklen_host(Int shcol, Real* uw_sfc, Real* vw_sfc, Real* wthl_sfc, Real* wqw_sfc,
                        Real* thl_sfc, Real* cldliq_sfc, Real* qv_sfc, Real* ustar, Real* kbfs, Real* obklen);
void shoc_pblintd_cldcheck_host(Int shcol, Int nlev, Int nlevi, Real* zi, Real* cldn, Real* pblh);
void compute_shr_prod_host(Int nlevi, Int nlev, Int shcol, Real* dz_zi, Real* u_wind, Real* v_wind, Real* sterm);
void shoc_length_host(Int shcol, Int nlev, Int nlevi, Real* host_dx, Real* host_dy,
                   Real* zt_grid, Real* zi_grid, Real*dz_zt, Real* tke,
                   Real* thv, Real*brunt, Real* shoc_mix);
void shoc_energy_fixer_host(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* zt_grid,
                         Real* zi_grid, Real* se_b, Real* ke_b, Real* wv_b, Real* wl_b,
                         Real* se_a, Real* ke_a, Real* wv_a, Real* wl_a, Real* wthl_sfc,
                         Real* wqw_sfc, Real* rho_zt, Real* tke, Real* pint,
                         Real* host_dse);
void compute_shoc_vapor_host(Int shcol, Int nlev, Real* qw, Real* ql, Real* qv);
void update_prognostics_implicit_host(Int shcol, Int nlev, Int nlevi, Int num_tracer, Real dtime,
                                   Real* dz_zt, Real* dz_zi, Real* rho_zt, Real* zt_grid,
                                   Real* zi_grid, Real* tk, Real* tkh, Real* uw_sfc, Real* vw_sfc,
                                   Real* wthl_sfc, Real* wqw_sfc, Real* wtracer_sfc, Real* thetal,
                                   Real* qw, Real* tracer, Real* tke, Real* u_wind, Real* v_wind);
void diag_third_shoc_moments_host(Int shcol, Int nlev, Int nlevi, Real* w_sec, Real* thl_sec,
                               Real* wthl_sec, Real* isotropy, Real* brunt, Real* thetal,
                               Real* tke, Real* dz_zt, Real* dz_zi, Real* zt_grid, Real* zi_grid,
                               Real* w3);
void adv_sgs_tke_host(Int nlev, Int shcol, Real dtime, Real* shoc_mix, Real* wthv_sec, Real* sterm_zt,
                   Real* tk, Real* tke, Real* a_diss);
void shoc_assumed_pdf_host(Int shcol, Int nlev, Int nlevi, Real* thetal, Real* qw, Real* w_field,
                        Real* thl_sec, Real* qw_sec, Real* wthl_sec, Real* w_sec, Real* wqw_sec,
                        Real* qwthl_sec, Real* w3, Real* pres, Real* zt_grid, Real* zi_grid,
                        Real* shoc_cldfrac, Real* shoc_ql, Real* wqls, Real* wthv_sec, Real* shoc_ql2);
void compute_tmpi_host(Int nlevi, Int shcol, Real dtime, Real *rho_zi, Real *dz_zi, Real *tmpi);
void integ_column_stability_host(Int nlev, Int shcol, Real *dz_zt,
                              Real *pres, Real *brunt, Real *brunt_int);
void isotropic_ts_host(Int nlev, Int shcol, Real* brunt_int, Real* tke,
                    Real* a_diss, Real* brunt, Real* isotropy);
void dp_inverse_host(Int nlev, Int shcol, Real *rho_zt, Real *dz_zt, Real *rdp_zt);

int shoc_init_host(Int nlev, Real* pref_mid, Int nbot_shoc, Int ntop_shoc);
Int shoc_main_host(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Int npbl, Real* host_dx, Real* host_dy, Real* thv,
                Real* zt_grid, Real* zi_grid, Real* pres, Real* presi, Real* pdel, Real* wthl_sfc, Real* wqw_sfc,
                Real* uw_sfc, Real* vw_sfc, Real* wtracer_sfc, Int num_qtracers, Real* w_field, Real* inv_exner,
                Real* phis, Real* host_dse, Real* tke, Real* thetal, Real* qw, Real* u_wind, Real* v_wind,
                Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk, Real* shoc_ql, Real* shoc_cldfrac, Real* pblh,
                Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec, Real* qw_sec, Real* qwthl_sec,
                Real* wthl_sec, Real* wqw_sec, Real* wtke_sec, Real* uw_sec, Real* vw_sec, Real* w3, Real* wqls_sec,
                Real* brunt, Real* shoc_ql2);

void pblintd_height_host(Int shcol, Int nlev, Int npbl, Real* z, Real* u, Real* v, Real* ustar, Real* thv, Real* thv_ref, Real* pblh, Real* rino, bool* check);

void vd_shoc_decomp_and_solve_host(Int shcol, Int nlev, Int nlevi, Int num_rhs, Real dtime, Real* kv_term, Real* tmpi, Real* rdp_zt, Real* flux, Real* var);


void pblintd_surf_temp_host(Int shcol, Int nlev, Int nlevi, Real* z, Real* ustar, Real* obklen, Real* kbfs, Real* thv, Real* tlv, Real* pblh, bool* check, Real* rino);

void pblintd_check_pblh_host(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* ustar, bool* check, Real* pblh);

void pblintd_host(Int shcol, Int nlev, Int nlevi, Int npbl, Real* z, Real* zi, Real* thl, Real* ql, Real* q, Real* u, Real* v, Real* ustar, Real* obklen, Real* kbfs, Real* cldn, Real* pblh);
void shoc_grid_host(Int shcol, Int nlev, Int nlevi, Real* zt_grid, Real* zi_grid, Real* pdel, Real* dz_zt, Real* dz_zi, Real* rho_zt);
void eddy_diffusivities_host(Int nlev, Int shcol, Real* pblh, Real* zt_grid, Real* tabs, Real* shoc_mix, Real* sterm_zt, Real* isotropy,
                          Real* tke, Real* tkh, Real* tk);
void shoc_tke_host(Int shcol, Int nlev, Int nlevi, Real dtime, Real* wthv_sec, Real* shoc_mix, Real* dz_zi, Real* dz_zt, Real* pres,
                Real* u_wind, Real* v_wind, Real* brunt, Real* obklen, Real* zt_grid, Real* zi_grid, Real* pblh, Real* tke,
                Real* tk, Real* tkh, Real* isotropy);
void compute_shoc_temperature_host(Int shcol, Int nlev, Real* thetal, Real* ql, Real* inv_exner, Real* tabs);

void shoc_energy_total_fixer_host(Int shcol, Int nlev, Int nlevi, Real dtime, Int nadv, Real* zt_grid, Real* zi_grid, Real* se_b, Real* ke_b, Real* wv_b, Real* wl_b, Real* se_a, Real* ke_a, Real* wv_a, Real* wl_a, Real* wthl_sfc, Real* wqw_sfc, Real* rho_zt, Real* pint, Real* te_a, Real* te_b);

void shoc_assumed_pdf_tilde_to_real_host(Real w_first, Real sqrtw2, Real* w1);

void shoc_assumed_pdf_vv_parameters_host(Real w_first, Real w_sec, Real w3var, Real w_tol_sqd, Real* skew_w, Real* w1_1, Real* w1_2, Real* w2_1, Real* w2_2, Real* a);

void shoc_assumed_pdf_thl_parameters_host(Real wthlsec, Real sqrtw2, Real sqrtthl, Real thlsec, Real thl_first, Real w1_1, Real w1_2, Real skew_w, Real a, Real thl_tol, Real w_thresh, Real* thl1_1, Real* thl1_2, Real* thl2_1, Real* thl2_2, Real* sqrtthl2_1, Real* sqrtthl2_2);

void shoc_assumed_pdf_qw_parameters_host(Real wqwsec, Real sqrtw2, Real skew_w, Real sqrtqt, Real qwsec, Real w1_2, Real w1_1, Real qw_first, Real a, Real rt_tol, Real w_thresh, Real* qw1_1, Real* qw1_2, Real* qw2_1, Real* qw2_2, Real* sqrtqw2_1, Real* sqrtqw2_2);

void shoc_assumed_pdf_inplume_correlations_host(Real sqrtqw2_1, Real sqrtthl2_1, Real a, Real sqrtqw2_2, Real sqrtthl2_2, Real qwthlsec, Real qw1_1, Real qw_first, Real thl1_1, Real thl_first, Real qw1_2, Real thl1_2, Real* r_qwthl_1);

void shoc_assumed_pdf_compute_temperature_host(Real thl1, Real pval, Real* tl1);

void shoc_assumed_pdf_compute_qs_host(Real tl1_1, Real tl1_2, Real pval, Real* qs1, Real* beta1, Real* qs2, Real* beta2);

void shoc_assumed_pdf_compute_s_host(Real qw1, Real qs1, Real beta, Real pval, Real thl2, Real qw2, Real sqrtthl2, Real sqrtqw2, Real r_qwthl, Real* s, Real* std_s, Real* qn, Real* c);

void shoc_assumed_pdf_compute_sgs_liquid_host(Real a, Real ql1, Real ql2, Real* shoc_ql);

void shoc_assumed_pdf_compute_cloud_liquid_variance_host(Real a, Real s1, Real ql1, Real c1, Real std_s1, Real s2, Real ql2, Real c2, Real std_s2, Real shoc_ql, Real* shoc_ql2);

void shoc_assumed_pdf_compute_liquid_water_flux_host(Real a, Real w1_1, Real w_first, Real ql1, Real w1_2, Real ql2, Real* wqls);

void shoc_assumed_pdf_compute_buoyancy_flux_host(Real wthlsec, Real wqwsec, Real pval, Real wqls, Real* wthv_sec);

// end _host function decls

}  // namespace shoc
}  // namespace scream

#endif // SCREAM_SHOC_FUNCTIONS_F90_HPP
