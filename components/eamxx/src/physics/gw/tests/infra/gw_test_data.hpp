#ifndef SCREAM_GW_FUNCTIONS_F90_HPP
#define SCREAM_GW_FUNCTIONS_F90_HPP

#include "physics/gw/gw_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

//
// Bridge functions to call fortran version of gw functions from C++
//

namespace scream {
namespace gw {

// The Data struct is special; it is used to do gw initialization, which
// must be called before any gw function.
struct GwInit : public PhysicsTestData {
  // Inputs
  Int pver, pgwv;
  Real dc;
  bool orographic_only, do_molec_diff, tau_0_ubc;
  Int nbot_molec, ktop, kbotbg;
  Real fcrit2, kwv;
  Real *cref, *alpha;

  GwInit(Int pver_, Int pgwv_, Real dc_, bool orographic_only_, bool do_molec_diff_, bool tau_0_ubc_, Int nbot_molec_, Int ktop_, Int kbotbg_, Real fcrit2_, Real kwv_) :
    PhysicsTestData({
      {pgwv_*2 + 1},
      {pver_ + 1}
    },
    {
      {&cref},
      {&alpha}
    }),
    pver(pver_), pgwv(pgwv_), dc(dc_), orographic_only(orographic_only_), do_molec_diff(do_molec_diff_), tau_0_ubc(tau_0_ubc_), nbot_molec(nbot_molec_), ktop(ktop_), kbotbg(kbotbg_), fcrit2(fcrit2_), kwv(kwv_)
  {
    // Assert valid init data?
    assert(kbotbg <= pver);
    assert(ktop >= 0);
    assert(kbotbg >= ktop);
    assert(pgwv > 0);
    assert(nbot_molec >= 0);
    assert(nbot_molec >= ktop);
  }

  PTD_STD_DEF(GwInit, 11, pver, pgwv, dc, orographic_only, do_molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    shift_int_scalar<D>(ktop);
    shift_int_scalar<D>(kbotbg);
  }
};

struct GwdComputeTendenciesFromStressDivergenceData : public PhysicsTestData {
  // Inputs
  Int ncol;
  bool do_taper;
  Real dt, effgw;
  Int *tend_level;
  Real *lat, *dpm, *rdpm, *c, *ubm, *t, *nm, *xv, *yv;
  GwInit init;

  // Inputs/Outputs
  Real *tau;

  // Outputs
  Real *gwut, *utgw, *vtgw;

  GwdComputeTendenciesFromStressDivergenceData(Int ncol_, bool do_taper_, Real dt_, Real effgw_, GwInit init_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, init_.pver},
      {ncol_, 2*init_.pgwv + 1},
      {ncol_, 2*init_.pgwv + 1, init_.pver + 1},
      {ncol_, init_.pver, 2*init_.pgwv + 1},
      {ncol_}
    },
    {
      {&lat, &xv, &yv},
      {&dpm, &rdpm, &ubm, &t, &nm, &utgw, &vtgw},
      {&c},
      {&tau},
      {&gwut}
    },
    {
      {&tend_level}
    }),
    ncol(ncol_), do_taper(do_taper_), dt(dt_), effgw(effgw_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwdComputeTendenciesFromStressDivergenceData, 4, ncol, do_taper, dt, effgw);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }

};

struct GwProfData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real cpair;
  Real *t, *pmid, *pint;
  GwInit init;

  // Outputs
  Real *rhoi, *ti, *nm, *ni;

  GwProfData(Int ncol_, Real cpair_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver},
      {ncol_, init_.pver + 1}
    },
    {
      {&t, &pmid, &nm},
      {&pint, &rhoi, &ti, &ni}
    }),
    ncol(ncol_), cpair(cpair_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwProfData, 2, ncol, cpair);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }

};

struct MomentumEnergyConservationData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Int *tend_level;
  Real dt;
  Real *taucd, *pint, *pdel, *u, *v;
  GwInit init;

  // Inputs/Outputs
  Real *dudt, *dvdt, *dsdt, *utgw, *vtgw, *ttgw;

  MomentumEnergyConservationData(Int ncol_, Real dt_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver + 1, 4},
      {ncol_, init_.pver + 1},
      {ncol_, init_.pver},
      {ncol_}
    },
    {
      {&taucd},
      {&pint},
      {&pdel, &u, &v, &utgw, &vtgw, &ttgw, &dudt, &dvdt, &dsdt}
    },
    {
      {&tend_level}
    }),
    ncol(ncol_), dt(dt_), init(init_)
  {}

  PTD_STD_DEF_INIT(MomentumEnergyConservationData, 2, ncol, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwdComputeStressProfilesAndDiffusivitiesData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Int *src_level;
  Real *ubi, *c, *rhoi, *ni, *kvtt, *t, *ti, *piln;
  GwInit init;

  // Inputs/Outputs
  Real *tau;

  GwdComputeStressProfilesAndDiffusivitiesData(Int ncol_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver + 1},
      {ncol_, init_.pgwv*2 + 1},
      {ncol_, init_.pver},
      {ncol_, init_.pgwv*2 + 1, init_.pver + 1},
      {ncol_}
    },
    {
      {&ubi, &rhoi, &ni, &kvtt, &ti, &piln},
      {&c},
      {&t},
      {&tau}
    },
    {
      {&src_level}
    }),
    ncol(ncol_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwdComputeStressProfilesAndDiffusivitiesData, 1, ncol);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwdProjectTauData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Int *tend_level;
  Real *tau, *ubi, *c, *xv, *yv;
  GwInit init;

  // Outputs
  Real *taucd;

  GwdProjectTauData(Int ncol_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pgwv*2 + 1, init_.pver + 1},
      {ncol_, init_.pver + 1},
      {ncol_, init_.pgwv*2 + 1},
      {ncol_},
      {ncol_, init_.pver + 1, 4},
      {ncol_}
    },
    {
      {&tau},
      {&ubi},
      {&c},
      {&xv, &yv},
      {&taucd}
    },
    {
      {&tend_level}
    }),
    ncol(ncol_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwdProjectTauData, 1, ncol);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwdPrecalcRhoiData : public PhysicsTestData {
  // Inputs
  Int pcnst, ncol;
  Real dt;
  Int *tend_level;
  Real *pmid, *pint, *t, *gwut, *ubm, *nm, *rdpm, *c, *q, *dse;
  GwInit init;

  // Outputs
  Real *egwdffi, *qtgw, *dttdf, *dttke, *ttgw;

  GwdPrecalcRhoiData(Int pcnst_, Int ncol_, Real dt_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver},
      {ncol_, init_.pver + 1},
      {ncol_, init_.pver, init_.pgwv*2 + 1},
      {ncol_, init_.pgwv*2 + 1},
      {ncol_, init_.pver, pcnst_},
      {ncol_}
    },
    {
      {&pmid, &t, &ubm, &nm, &rdpm, &dse, &dttdf, &dttke, &ttgw},
      {&pint, &egwdffi},
      {&gwut},
      {&c},
      {&q, &qtgw}
    },
    {
      {&tend_level}
    }),
    pcnst(pcnst_), ncol(ncol_), dt(dt_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwdPrecalcRhoiData, 3, pcnst, ncol, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwDragProfData : public PhysicsTestData {
  // Inputs
  Int pcnst, ncol;
  Int *src_level, *tend_level;
  bool do_taper;
  Real dt, effgw;
  Real *lat, *t, *ti, *pmid, *pint, *dpm, *rdpm, *piln, *rhoi, *nm, *ni, *ubm, *ubi, *xv, *yv, *c, *kvtt, *q, *dse;
  GwInit init;

  // Inputs/Outputs
  Real *tau;

  // Outputs
  Real *utgw, *vtgw, *ttgw, *qtgw, *taucd, *egwdffi, *gwut, *dttdf, *dttke;

  GwDragProfData(Int pcnst_, Int ncol_, bool do_taper_, Real dt_, Real effgw_, GwInit init_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, init_.pver},
      {ncol_, init_.pver + 1},
      {ncol_, init_.pgwv*2 + 1},
      {ncol_, init_.pver, pcnst_},
      {ncol_, init_.pgwv*2 + 1, init_.pver + 1},
      {ncol_, init_.pver + 1, 4},
      {ncol_, init_.pver, init_.pgwv*2 + 1},
      {ncol_}
    },
    {
      {&lat, &xv, &yv},
      {&t, &pmid, &dpm, &rdpm, &nm, &ubm, &dse, &utgw, &vtgw, &ttgw, &dttdf, &dttke},
      {&ti, &pint, &piln, &rhoi, &ni, &ubi, &kvtt, &egwdffi},
      {&c},
      {&q, &qtgw},
      {&tau},
      {&taucd},
      {&gwut}
    },
    {
      {&src_level, &tend_level}
    }),
    pcnst(pcnst_), ncol(ncol_), do_taper(do_taper_), dt(dt_), effgw(effgw_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwDragProfData, 5, pcnst, ncol, do_taper, dt, effgw);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwFrontInitData : public PhysicsTestData{
  // Inputs
  Real taubgnd, frontgfc_in;
  Int kfront_in;
  GwInit init;

  GwFrontInitData(Real taubgnd_, Real frontgfc_in_, Int kfront_in_, GwInit init_) :
    PhysicsTestData({}, {}, {}),
    taubgnd(taubgnd_),
    frontgfc_in(frontgfc_in_),
    kfront_in(kfront_in_),
    init(init_)
  {}

  PTD_STD_DEF_INIT(GwFrontInitData, 3, taubgnd, frontgfc_in, kfront_in);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();

    shift_int_scalar<D>(kfront_in);
  }
};

struct GwFrontProjectWindsData : public PhysicsTestData {
  // Inputs
  Int ncol, kbot;
  Real *u, *v;
  GwFrontInitData init;

  // Outputs
  Real *xv, *yv, *ubm, *ubi;

  GwFrontProjectWindsData(Int ncol_, Int kbot_, GwFrontInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_},
      {ncol_, init_.init.pver + 1}
    },
    {
      {&u, &v, &ubm},
      {&xv, &yv},
      {&ubi}
    }),
    ncol(ncol_), kbot(kbot_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwFrontProjectWindsData, 2, ncol, kbot);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwFrontGwSourcesData : public PhysicsTestData {
  // Inputs
  Int ncol, kbot;
  Real *frontgf;
  GwFrontInitData init;

  // Outputs
  Real *tau;

  GwFrontGwSourcesData(Int ncol_, Int kbot_, GwFrontInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_, init_.init.pgwv*2 + 1, init_.init.pver + 1}
    },
    {
      {&frontgf},
      {&tau}
    }),
    ncol(ncol_), kbot(kbot_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwFrontGwSourcesData, 2, ncol, kbot);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwCmSrcData : public PhysicsTestData {
  // Inputs
  Int ncol, kbot;
  Real *u, *v, *frontgf;
  GwFrontInitData init;

  // Outputs
  Int *src_level, *tend_level;
  Real *tau, *ubm, *ubi, *xv, *yv, *c;

  GwCmSrcData(Int ncol_, Int kbot_, GwFrontInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_, init_.init.pgwv*2 + 1, init_.init.pver + 1},
      {ncol_, init_.init.pver + 1},
      {ncol_},
      {ncol_, init_.init.pgwv*2 + 1},
      {ncol_}
    },
    {
      {&u, &v, &ubm, &frontgf},
      {&tau},
      {&ubi},
      {&xv, &yv},
      {&c}
    },
    {
      {&src_level, &tend_level}
    }),
    ncol(ncol_), kbot(kbot_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwCmSrcData, 2, ncol, kbot);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();

    shift_int_scalar<D>(kbot);
  }
};

struct GwConvectInitData : public PhysicsTestData{
  // Inputs
  Int maxh, maxuh;
  Real plev_src_wind;
  Real *mfcc_in;
  GwInit init;

  GwConvectInitData(Int maxh_, Int maxuh_, Real plev_src_wind_, GwInit init_) :
    PhysicsTestData({
      {maxh_, maxuh_*2 + 1, init_.pgwv*2 + 1}
    },
    {
      {&mfcc_in}
    }),
    maxh(maxh_), maxuh(maxuh_), plev_src_wind(plev_src_wind_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwConvectInitData, 3, maxh, maxuh, plev_src_wind);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwConvectProjectWindsData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real *u, *v;
  GwConvectInitData init;

  // Outputs
  Real *xv, *yv, *ubm, *ubi;

  GwConvectProjectWindsData(Int ncol_, GwConvectInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_},
      {ncol_, init_.init.pver + 1}
    },
    {
      {&u, &v, &ubm},
      {&xv, &yv},
      {&ubi}
    }),
    ncol(ncol_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwConvectProjectWindsData, 1, ncol);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwHeatingDepthData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real maxq0_conversion_factor, hdepth_scaling_factor;
  bool use_gw_convect_old;
  Real *zm, *netdt;
  GwConvectInitData init;

  // Outputs
  Int *mini, *maxi;
  Real *hdepth, *maxq0_out, *maxq0;

  GwHeatingDepthData(Int ncol_, Real maxq0_conversion_factor_, Real hdepth_scaling_factor_, bool use_gw_convect_old_, GwConvectInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_},
      {ncol_}
    },
    {
      {&zm, &netdt},
      {&hdepth, &maxq0_out, &maxq0}
    },
    {
      {&mini, &maxi}
    }),
    ncol(ncol_), maxq0_conversion_factor(maxq0_conversion_factor_), hdepth_scaling_factor(hdepth_scaling_factor_), use_gw_convect_old(use_gw_convect_old_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwHeatingDepthData, 4, ncol, maxq0_conversion_factor, hdepth_scaling_factor, use_gw_convect_old);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwStormSpeedData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real storm_speed_min;
  Real *ubm;
  Int *mini, *maxi;
  GwConvectInitData init;

  // Outputs
  Int *storm_speed;
  Real *uh, *umin, *umax;

  GwStormSpeedData(Int ncol_, Real storm_speed_min_, GwConvectInitData init_) :
    PhysicsTestData({
      {ncol_, init_.init.pver},
      {ncol_},
      {ncol_}
    },
    {
      {&ubm},
      {&uh, &umin, &umax}
    },
    {
      {&mini, &maxi, &storm_speed}
    }),
    ncol(ncol_), storm_speed_min(storm_speed_min_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwStormSpeedData, 2, ncol, storm_speed_min);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwConvectGwSourcesData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real *lat, *hdepth, *netdt, *uh, *maxq0, *umin, *umax;
  Real hdepth_min;
  Int *mini, *maxi, *storm_speed;
  GwConvectInitData init;

  // Outputs
  Real *tau;

  GwConvectGwSourcesData(Int ncol_, Real hdepth_min_, GwConvectInitData init_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, init_.init.pver},
      {ncol_, init_.init.pgwv*2 + 1, init_.init.pver + 1},
      {ncol_}
    },
    {
      {&lat, &hdepth, &uh, &maxq0, &umin, &umax},
      {&netdt},
      {&tau}
    },
    {
      {&mini, &maxi, &storm_speed}
    }),
    ncol(ncol_), hdepth_min(hdepth_min_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwConvectGwSourcesData, 2, ncol, hdepth_min);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwBeresSrcData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real *lat, *u, *v, *netdt, *zm;
  Real maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min;
  bool use_gw_convect_old;
  GwConvectInitData init;

  // Outputs
  Int *src_level, *tend_level;
  Real *tau, *ubm, *ubi, *xv, *yv, *c, *hdepth, *maxq0_out;

  GwBeresSrcData(Int ncol_, Real maxq0_conversion_factor_, Real hdepth_scaling_factor_, Real hdepth_min_, Real storm_speed_min_, bool use_gw_convect_old_, GwConvectInitData init_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, init_.init.pver},
      {ncol_, init_.init.pgwv*2 + 1, init_.init.pver + 1},
      {ncol_, init_.init.pver + 1},
      {ncol_, init_.init.pgwv*2 + 1},
      {ncol_}
    },
    {
      {&lat, &xv, &yv, &hdepth, &maxq0_out},
      {&u, &v, &zm, &ubm, &netdt},
      {&tau},
      {&ubi},
      {&c}
    },
    {
      {&src_level, &tend_level}
    }),
    ncol(ncol_), maxq0_conversion_factor(maxq0_conversion_factor_), hdepth_scaling_factor(hdepth_scaling_factor_), hdepth_min(hdepth_min_), storm_speed_min(storm_speed_min_), use_gw_convect_old(use_gw_convect_old_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwBeresSrcData, 6, ncol, maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min, use_gw_convect_old);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

struct GwEdiffData : public PhysicsTestData {
  // Inputs
  Int ncol, kbot, ktop;
  Int *tend_level;
  Real *gwut, *ubm, *nm, *rho, *pmid, *rdpm, *c;
  Real dt;
  GwInit init;

  // Outputs
  Real *egwdffi;
  Real *decomp_ca, *decomp_cc, *decomp_dnom, *decomp_ze;

  GwEdiffData(Int ncol_, Int kbot_, Int ktop_, Real dt_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver, 2*init_.pgwv + 1},
      {ncol_, init_.pver},
      {ncol_, init_.pver + 1},
      {ncol_, 2*init_.pgwv + 1},
      {ncol_}
    },
    {
      {&gwut},
      {&ubm, &nm, &pmid, &rdpm, &decomp_ca, &decomp_cc, &decomp_dnom, &decomp_ze},
      {&rho, &egwdffi},
      {&c}
    },
    {
      {&tend_level}
    }),
    ncol(ncol_), kbot(kbot_), ktop(ktop_), dt(dt_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwEdiffData, 4, ncol, kbot, ktop, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();

    shift_int_scalar<D>(kbot);
    shift_int_scalar<D>(ktop);
  }
};

struct GwDiffTendData : public PhysicsTestData {
  // Inputs
  Int ncol, kbot, ktop;
  Real *q;
  Real dt;
  Real *decomp_ca, *decomp_cc, *decomp_dnom, *decomp_ze;
  GwInit init;

  // Outputs
  Real *dq;

  GwDiffTendData(Int ncol_, Int kbot_, Int ktop_, Real dt_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver}
    },
    {
      {&q, &dq, &decomp_ca, &decomp_cc, &decomp_dnom, &decomp_ze}
    }),
    ncol(ncol_), kbot(kbot_), ktop(ktop_), dt(dt_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwDiffTendData, 4, ncol, kbot, ktop, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();

    shift_int_scalar<D>(kbot);
    shift_int_scalar<D>(ktop);
  }
};

struct GwOroSrcData : public PhysicsTestData {
  // Inputs
  Int ncol;
  Real *u, *v, *t, *sgh, *pmid, *pint, *dpm, *zm, *nm;
  GwInit init;

  // Outputs
  Int *src_level, *tend_level;
  Real *tau, *ubm, *ubi, *xv, *yv, *c;

  GwOroSrcData(Int ncol_, GwInit init_) :
    PhysicsTestData({
      {ncol_, init_.pver},
      {ncol_},
      {ncol_, init_.pver + 1},
      {ncol_, init_.pgwv*2 + 1, init_.pver + 1},
      {ncol_, init_.pgwv*2 + 1},
      {ncol_}
    },
    {
      {&u, &v, &t, &pmid, &dpm, &zm, &nm, &ubm},
      {&sgh, &xv, &yv},
      {&pint, &ubi},
      {&tau},
      {&c}
    },
    {
      {&src_level, &tend_level}
    }),
    ncol(ncol_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwOroSrcData, 1, ncol);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();

    init.transition<D>();
  }
};

// Glue functions to call fortran from from C++ with the Data struct
void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d);
void gwd_compute_tendencies_from_stress_divergence_f(GwdComputeTendenciesFromStressDivergenceData& d);
void gw_prof(GwProfData& d);
void momentum_energy_conservation_f(MomentumEnergyConservationData& d);
void momentum_energy_conservation(MomentumEnergyConservationData& d);
void gwd_compute_stress_profiles_and_diffusivities_f(GwdComputeStressProfilesAndDiffusivitiesData& d);
void gwd_compute_stress_profiles_and_diffusivities(GwdComputeStressProfilesAndDiffusivitiesData& d);
void gwd_project_tau(GwdProjectTauData& d);
void gwd_precalc_rhoi(GwdPrecalcRhoiData& d);
void gw_drag_prof(GwDragProfData& d);
void gw_front_project_winds(GwFrontProjectWindsData& d);
void gw_front_gw_sources(GwFrontGwSourcesData& d);
void gw_cm_src(GwCmSrcData& d);
void gw_convect_project_winds(GwConvectProjectWindsData& d);
void gw_heating_depth(GwHeatingDepthData& d);
void gw_storm_speed(GwStormSpeedData& d);
void gw_convect_gw_sources(GwConvectGwSourcesData& d);
void gw_beres_src(GwBeresSrcData& d);
void gw_ediff(GwEdiffData& d);
void gw_diff_tend(GwDiffTendData& d);
void gw_oro_src(GwOroSrcData& d);

extern "C" { // _f function decls
}

}  // namespace gw
}  // namespace scream

#endif
