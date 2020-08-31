#ifndef SCREAM_SHOC_FUNCTIONS_F90_HPP
#define SCREAM_SHOC_FUNCTIONS_F90_HPP

#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include "shoc_functions.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of shoc functions from C++
//

namespace scream {
namespace shoc {

// Base class for common SHOC data setup
struct SHOCDataBase
{
  Int shcol, nlev, nlevi;

  SHOCDataBase(Int shcol_, Int nlev_, Int nlevi_,
               const std::vector<Real**>& ptrs = {}, const std::vector<Real**>& ptrs_i = {},
	       const std::vector<Real**>& ptrs_c = {}, const std::vector<Int**>& idx_c = {});

  SHOCDataBase(const SHOCDataBase &rhs) = delete;

  SHOCDataBase(const SHOCDataBase &rhs,
               const std::vector<Real**>& ptrs = {}, const std::vector<Real**>& ptrs_i = {},
               const std::vector<Real**>& ptrs_c = {}, const std::vector<Int**>& idx_c = {});

  SHOCDataBase &operator=(const SHOCDataBase &rhs);

  // Since we are also preparing index data, this function is doing more than transposing. It's shifting the
  // format of all data from one language to another
  template <util::TransposeDirection::Enum D>
  void transpose() {
    std::vector<Real> data(m_data.size());

    // Transpose on the zt grid
    for (size_t i = 0; i < m_ptrs.size(); ++i) {
      util::transpose<D>(*(m_ptrs[i]), data.data() + (m_total*i) , shcol, nlev);
    }

    // Transpose on the zi grid
    for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
      util::transpose<D>(*(m_ptrs_i[i]), data.data() + (m_ptrs.size()*m_total) + (m_totali*i), shcol, nlevi);
    }

    // Copy the column only grid
    const Int c_start_offset = m_ptrs.size()*m_total + m_ptrs_i.size()*m_totali;
    std::copy(m_data.begin() + c_start_offset, m_data.end(), data.begin() + c_start_offset);

    m_data = data;

    // Shift the indices
    for (size_t i = 0; i < m_idx_data.size(); ++i) {
      m_idx_data[i] += (D == util::TransposeDirection::c2f ? 1 : -1);
      scream_assert_msg(m_idx_data[i] >= 0, "Bad index: " << m_idx_data[i]);
    }
  }

  void randomize(const std::vector<std::pair<Real, Real> >& ranges = {},
                 const std::vector<std::pair<Real, Real> >& ranges_i = {},
                 const std::vector<std::pair<Real, Real> >& ranges_c = {},
                 const std::vector<std::pair<Int, Int> >&   ranges_idx = {});

  Int total() const { return m_total; }
  Int totali() const { return m_totali; }

 private:
  void init_ptrs();

  // Internals
  Int m_total, m_totali;
  std::vector<Real**> m_ptrs, m_ptrs_i, m_ptrs_c;
  std::vector<Int**> m_indices_c;
  std::vector<Real> m_data;
  std::vector<Int> m_idx_data;
};

//Create data structure to hold data for shoc_grid
struct SHOCGridData : public SHOCDataBase {
  // Inputs
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&zt_grid, &dz_zt, &pdel, &rho_zt}, {&zi_grid, &dz_zi}) {}

  SHOCGridData(const SHOCGridData &rhs) : SHOCDataBase(rhs, {&zt_grid, &dz_zt, &pdel, &rho_zt}, {&zi_grid, &dz_zi}) {}

  SHOCGridData &operator=(const SHOCGridData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};

//Create data structure to hold data for integ_column_stability
struct SHOCColstabData : public SHOCDataBase {
  // Inputs
  Real *dz_zt, *pres, *brunt;

  // Output
  Real *brunt_int;

  SHOCColstabData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&dz_zt, &pres, &brunt}, {}, {&brunt_int}) {}
  SHOCColstabData(const SHOCColstabData &rhs) : SHOCDataBase(rhs, {&dz_zt, &pres, &brunt}, {}, {&brunt_int}) {}
  SHOCColstabData &operator=(const SHOCColstabData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCColstabData

//Create data structure to hold data for compute_shr_prod
struct SHOCTkeshearData : public SHOCDataBase {
  // Inputs
  Real *dz_zi, *u_wind, *v_wind;

  // In/out
  Real *sterm;

  //functions to initialize data
  SHOCTkeshearData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&u_wind, &v_wind}, {&dz_zi, &sterm}) {}
  SHOCTkeshearData(const SHOCTkeshearData &rhs) : SHOCDataBase(rhs, {&u_wind, &v_wind}, {&dz_zi, &sterm}) {}
  SHOCTkeshearData &operator=(const SHOCTkeshearData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCTkeshearData

//Create data structure to hold data for isotropic_ts
struct SHOCIsotropicData : public SHOCDataBase {
  // Inputs
  Real *tke, *a_diss, *brunt, *brunt_int;

  // Output
  Real *isotropy;

  //functions to initialize data
  SHOCIsotropicData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&tke, &a_diss, &brunt, &isotropy}, {}, {&brunt_int}) {}
  SHOCIsotropicData(const SHOCIsotropicData &rhs) : SHOCDataBase(rhs, {&tke, &a_diss, &brunt, &isotropy}, {}, {&brunt_int}) {}
  SHOCIsotropicData &operator=(const SHOCIsotropicData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCIsotropicData

//Create data structure to hold data for adv_sgs_tke
struct SHOCAdvsgstkeData : public SHOCDataBase {
  // Inputs
  Real dtime;
  Real *shoc_mix, *wthv_sec, *sterm_zt, *tk;

  // In/out
  Real *tke;

  // Outputs
  Real *a_diss;

  //functions to initialize data
  SHOCAdvsgstkeData(Int shcol_, Int nlev_, Real dtime_) :
    SHOCDataBase(shcol_, nlev_, 0, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}), dtime(dtime_) {}
  SHOCAdvsgstkeData(const SHOCAdvsgstkeData &rhs) : SHOCDataBase(rhs, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}), dtime(rhs.dtime) {}
  SHOCAdvsgstkeData &operator=(const SHOCAdvsgstkeData &rhs)
  { SHOCDataBase::operator=(rhs); dtime = rhs.dtime; return *this; }
};//SHOCAdvsgstkeData

//Create data structure to hold data for eddy_diffusivities
struct SHOCEddydiffData : public SHOCDataBase {
  // Inputs
  Real *pblh, *obklen, *zt_grid, *shoc_mix, *sterm_zt,
        *isotropy, *tke;

  // Output
  Real *tk, *tkh;

  //functions to initialize data
  SHOCEddydiffData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {}, {&obklen, &pblh}) {}
  SHOCEddydiffData(const SHOCEddydiffData &rhs) : SHOCDataBase(rhs, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {}, {&obklen, &pblh}) {}
  SHOCEddydiffData &operator=(const SHOCEddydiffData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEddydiffData


//create data structure for update_host_dse
struct SHOCEnergydseData : public SHOCDataBase {
  // Inputs
  Real *thlm, *shoc_ql, *exner, *zt_grid, *phis;

  // Output
  Real *host_dse;

  //functions to initialize data
  SHOCEnergydseData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&thlm, &shoc_ql, &exner, &zt_grid, &host_dse}, {}, {&phis}) {}
  SHOCEnergydseData(const SHOCEnergydseData &rhs) : SHOCDataBase(rhs, {&thlm, &shoc_ql, &exner, &zt_grid, &host_dse}, {}, {&phis}) {}
  SHOCEnergydseData &operator=(const SHOCEnergydseData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEnergydseData

//create data structure for shoc_energy_integrals
struct SHOCEnergyintData : public SHOCDataBase {
  // Inputs
  Real *host_dse, *pdel, *rtm, *rcm, *u_wind, *v_wind;

  // Output
  Real *se_int, *ke_int, *wv_int, *wl_int;

  //functions to initialize data
  SHOCEnergyintData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&host_dse, &pdel, &rtm, &rcm, &u_wind, &v_wind}, {}, {&se_int, &ke_int, &wv_int, &wl_int}) {}
  SHOCEnergyintData(const SHOCEnergyintData &rhs) : SHOCDataBase(rhs, {&host_dse, &pdel, &rtm, &rcm, &u_wind, &v_wind}, {}, {&se_int, &ke_int, &wv_int, &wl_int}) {}
  SHOCEnergyintData &operator=(const SHOCEnergyintData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEnergyintData

//Create data structure for shoc_energy_total_fixer
struct SHOCEnergytotData : public SHOCDataBase {
  // Inputs
  Int nadv;
  Real dtime;
  Real *zt_grid, *zi_grid, *se_b, *ke_b, *wv_b, *wl_b, *se_a;
  Real *ke_a, *wv_a, *wl_a, *wthl_sfc, *wqw_sfc, *rho_zt;

  // Output
  Real *te_a, *te_b;

  //functions to initialize data for shoc_energy_total_fixer
  SHOCEnergytotData(Int shcol_, Int nlev_, Int nlevi_, Real dtime_, Int nadv_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&zt_grid, &rho_zt}, {&zi_grid}, {&se_b, &ke_b, &wv_b, &wl_b, &se_a, &ke_a, &wv_a, &wl_a, &wthl_sfc, &wqw_sfc, &te_a, &te_b}), dtime(dtime_), nadv(nadv_) {}
  SHOCEnergytotData(const SHOCEnergytotData &rhs) : SHOCDataBase(rhs, {&zt_grid, &rho_zt}, {&zi_grid}, {&se_b, &ke_b, &wv_b, &wl_b, &se_a, &ke_a, &wv_a, &wl_a, &wthl_sfc, &wqw_sfc, &te_a, &te_b}) {}
  SHOCEnergytotData &operator=(const SHOCEnergytotData &rhs) { SHOCDataBase::operator=(rhs); dtime = rhs.dtime; nadv = rhs.nadv; return *this; }
};//SHOCEnergytotData

//create data structure for shoc_energy_threshold_fixer
struct SHOCEnergythreshfixerData : public SHOCDataBase {
  // Inputs
  Real *pint, *tke, *te_a, *te_b;

  // In/out
  Real *se_dis;
  Int *shoctop;

  //functions to initialize data
  SHOCEnergythreshfixerData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&tke}, {&pint}, {&se_dis, &te_a, &te_b}, {&shoctop}) {}
  SHOCEnergythreshfixerData(const SHOCEnergythreshfixerData &rhs) : SHOCDataBase(rhs, {&tke}, {&pint}, {&se_dis, &te_a, &te_b}, {&shoctop}) {}
  SHOCEnergythreshfixerData &operator=(const SHOCEnergythreshfixerData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEnergythreshfixerData

//create data structure for shoc_energy_dse_fixer
struct SHOCEnergydsefixerData : public SHOCDataBase {
  // Inputs
  Real *se_dis;
  Int *shoctop;

  // In/out
  Real *host_dse;

  //functions to initialize data
  SHOCEnergydsefixerData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&host_dse}, {}, {&se_dis}, {&shoctop}) {}
  SHOCEnergydsefixerData(const SHOCEnergydsefixerData &rhs) : SHOCDataBase(rhs, {&host_dse}, {}, {&se_dis}, {&shoctop}) {}
  SHOCEnergydsefixerData &operator=(const SHOCEnergydsefixerData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEnergydsefixerData

//Create data structure to hold data for calc_shoc_vertflux
struct SHOCVertfluxData : public SHOCDataBase {
  // Inputs
  Real *tkh_zi, *dz_zi, *invar;

  // In/out
  Real *vertflux;

  SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData(const SHOCVertfluxData &rhs) : SHOCDataBase(rhs, {&invar}, {&tkh_zi, &dz_zi, &vertflux}) {}
  SHOCVertfluxData &operator=(const SHOCVertfluxData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
}; //SHOCVertfluxData

//Create data structure to hold data for calc_shoc_varorcovar
struct SHOCVarorcovarData : public SHOCDataBase {
  // Inputs
  Real tunefac;
  Real *tkh_zi, *dz_zi, *isotropy_zi, *invar1, *invar2;

  // In/out
  Real *varorcovar;

  SHOCVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&invar1, &invar2}, {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(tunefac_) {}
  SHOCVarorcovarData(const SHOCVarorcovarData &rhs) :
    SHOCDataBase(rhs, {&invar1, &invar2}, {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar}), tunefac(rhs.tunefac) {}
  SHOCVarorcovarData &operator=(const SHOCVarorcovarData &rhs)
  { SHOCDataBase::operator=(rhs); tunefac = rhs.tunefac; return *this; }
};//SHOCVarorcovarData

//Create data structure to hold data for compute_brunt_shoc_length
struct SHOCBruntlengthData : public SHOCDataBase {
  // Inputs
  Real *dz_zt, *thv, *thv_zi;

  // In/out
  Real *brunt;

  SHOCBruntlengthData(Int shcol_, Int nlev_, Int nlevi_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&dz_zt, &thv, &brunt}, {&thv_zi}) {}
  SHOCBruntlengthData(const SHOCBruntlengthData &rhs) :
    SHOCDataBase(rhs, {&dz_zt, &thv, &brunt}, {&thv_zi}) {}
  SHOCBruntlengthData &operator=(const SHOCBruntlengthData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCBruntlengthData

//Create data structure to hold data for compute_l_inf_shoc_length
struct SHOCInflengthData : public SHOCDataBase {
  // Inputs
  Real *zt_grid, *dz_zt, *tke;

  // In/out
  Real *l_inf;

  SHOCInflengthData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&zt_grid, &dz_zt, &tke}, {}, {&l_inf}) {}
  SHOCInflengthData(const SHOCInflengthData &rhs) :
    SHOCDataBase(rhs, {&zt_grid, &dz_zt, &tke}, {}, {&l_inf}) {}
  SHOCInflengthData &operator=(const SHOCInflengthData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCInflengthData

//Create data structure to hold data for compute_vel_shoc_length
struct SHOCConvvelData : public SHOCDataBase {
  // Inputs
  Real *pblh, *zt_grid, *dz_zt, *thv, *wthv_sec;

  // In/out
  Real *conv_vel;

  SHOCConvvelData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {}, {&conv_vel, &pblh}) {}
  SHOCConvvelData(const SHOCConvvelData &rhs) :
    SHOCDataBase(rhs, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {}, {&conv_vel, &pblh}) {}
  SHOCConvvelData &operator=(const SHOCConvvelData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCConvvelData

//Create data structure to hold data for compute_conv_time_shoc_length
struct SHOCConvtimeData : public SHOCDataBase {
  // Inputs
  Real *pblh, *conv_vel;

  // In/out
  Real *tscale;

  SHOCConvtimeData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {}, {}, {&conv_vel, &pblh, &tscale}) {}
  SHOCConvtimeData(const SHOCConvtimeData &rhs) :
    SHOCDataBase(rhs, {}, {}, {&conv_vel, &pblh, &tscale}) {}
  SHOCConvtimeData &operator=(const SHOCConvtimeData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCConvtimeData

//Create data structure to hold data for compute_shoc_mix_shoc_length
struct SHOCMixlengthData : public SHOCDataBase {
  // Inputs
  Real *tke, *brunt, *tscale, *zt_grid, *l_inf;

  // In/out
  Real *shoc_mix;

  SHOCMixlengthData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&tke, &brunt, &zt_grid, &shoc_mix}, {}, {&l_inf, &tscale}) {}
  SHOCMixlengthData(const SHOCMixlengthData &rhs) :
    SHOCDataBase(rhs, {&tke, &brunt, &zt_grid, &shoc_mix}, {}, {&l_inf, &tscale}) {}
  SHOCMixlengthData &operator=(const SHOCMixlengthData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCMixlengthData

//Create data structure to hold data for check_length_scale_shoc_length
struct SHOCMixcheckData : public SHOCDataBase {
  // Inputs
  Real *host_dx, *host_dy;

  // In/out
  Real *shoc_mix;

  SHOCMixcheckData(Int shcol_, Int nlev_) :
    SHOCDataBase(shcol_, nlev_, 0, {&shoc_mix}, {}, {&host_dx, &host_dy}) {}
  SHOCMixcheckData(const SHOCMixcheckData &rhs) :
    SHOCDataBase(rhs, {&shoc_mix}, {}, {&host_dx, &host_dy}) {}
  SHOCMixcheckData &operator=(const SHOCMixcheckData &rhs)
  { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCMixcheckData

struct SHOCSecondMomentSrfData : public SHOCDataBase {
  // Inputs
  Real *wthl, *uw, *vw;

  // out
  Real *ustar2, *wstar;

  SHOCSecondMomentSrfData(Int shcol_) :
  SHOCDataBase(shcol_, 1, 1, {&wthl, &uw, &vw, &ustar2, &wstar}, {}) {}
  SHOCSecondMomentSrfData(const SHOCSecondMomentSrfData &rhs) : SHOCDataBase(rhs, {&wthl, &uw, &vw, &ustar2, &wstar}, {}) {}
  SHOCSecondMomentSrfData &operator=(const SHOCSecondMomentSrfData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};

//Create data structure to hold data for linear_interp
struct SHOCLinearintData : public SHOCDataBase {
  // Inputs
  Real minthresh;
  Real *x1, *x2, *y1;

  // In/out
  Real *y2;

  SHOCLinearintData(Int shcol_, Int nlev_, Int nlevi_, Real minthresh_) :
    SHOCDataBase(shcol_, nlev_, nlevi_, {&x1, &y1}, {&x2, &y2}, {}), minthresh(minthresh_) {}
  SHOCLinearintData(const SHOCLinearintData &rhs) :
    SHOCDataBase(rhs, {&x1, &y1}, {&x2, &y2}, {}), minthresh(rhs.minthresh) {}
  SHOCLinearintData &operator=(const SHOCLinearintData &rhs)
  { SHOCDataBase::operator=(rhs); minthresh = rhs.minthresh; return *this; }
};//SHOCLinearintData

//Create data structure to hold data for shoc_assumed_pdf_tilda_to_real
struct SHOCPDFtildaData
{
  // inputs
  Real w_first, sqrtw2;
  
  // outputs
  Real w1;

};

// Create data structure to hold data for shoc_assumed_pdf_vv_parameters
struct SHOCPDFvvparamData
{
  // inputs
  Real w_first, w_sec, w3var;
  
  // outputs
  Real Skew_w, w1_1, w1_2, w2_1, w2_2, a;
};

// Create data structure to hold data for shoc_assumed_pdf_thl_parameters
struct SHOCPDFthlparamData
{
  // inputs
  Real wthlsec, sqrtw2, sqrtthl, thlsec, thl_first, w1_1, w1_2, Skew_w, a;
  bool dothetal_skew;
  
  // outputs
  Real thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2, Skew_thl;
};

// Create data structure to hold data for shoc_assumed_pdf_qw_parameters
struct SHOCPDFqwparamData
{
  // inputs
  Real wqwsec, qwsec, sqrtw2, sqrtqt, qw_first, w1_1, w1_2, Skew_w, a;
  
  // outputs
  Real qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2, Skew_qw;
};

// Create data structure to hold data for shoc_assumed_pdf_inplume_correlations
struct SHOCPDFinplumeData
{
  // inputs
  Real sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2;
  Real qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2;
  
  // outputs
  Real r_qwthl_1;
};

//Create data structure to hold data for shoc_assumed_pdf_compute_temperature
struct SHOCPDFcomptempData
{
  // inputs
  Real thl1, basepres, pval;
  
  // outputs
  Real Tl1;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_qs
struct SHOCPDFcompqsData
{
  // inputs
  Real Tl1_1, Tl1_2, pval;
  
  // outputs
  Real qs1, beta1, qs2, beta2;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_s
struct SHOCPDFcompsData
{
  // inputs
  Real qw1, qs1, beta, pval, thl2, qw2, sqrtthl2, sqrtqw2, r_qwthl;
  
  // outputs
  Real s, cthl, cqt, std_s, qn, C;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_sgs_liquid
struct SHOCPDFcompsgsliqData
{
  // inputs
  Real a, ql1, ql2;
  
  // outputs
  Real shoc_ql;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_cloud_liquid_variance
struct SHOCPDFcompcloudvarData
{
  // inputs
  Real a, s1, ql1, C1, std_s1, s2, ql2, C2, std_s2, shoc_ql;
  
  // outputs
  Real shoc_ql2;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_liquid_water_flux
struct SHOCPDFcompliqfluxData
{
  // inputs
  Real a, w1_1, w_first, ql1, w1_2, ql2;
  
  // outputs
  Real wqls;

};

//Create data structure to hold data for shoc_assumed_pdf_compute_buoyancy_flux
struct SHOCPDFcompbuoyfluxData
{
  // inputs
  Real wthlsec, epsterm, wqwsec, pval, wqls;
  
  // outputs
  Real wthv_sec;

};

//
// Glue functions to call fortran from from C++ with the Data struct
//

// This function initialzes the grid used by shoc. Given the
// locations of the cell center (location of thermodynaics quantities), cell
// interfaces, and pressure gradient the functon returns dz_zi, dz_zt,
// and density.
void shoc_grid(SHOCGridData &d);
void update_host_dse(SHOCEnergydseData &d);
void shoc_energy_integrals(SHOCEnergyintData &d);
void shoc_energy_total_fixer(SHOCEnergytotData &d);
void shoc_energy_threshold_fixer(SHOCEnergythreshfixerData &d);
void shoc_energy_dse_fixer(SHOCEnergydsefixerData &d);
void calc_shoc_vertflux(SHOCVertfluxData &d);
void calc_shoc_varorcovar(SHOCVarorcovarData &d);
void integ_column_stability(SHOCColstabData &d);
void compute_shr_prod(SHOCTkeshearData &d);
void isotropic_ts(SHOCIsotropicData &d);
void adv_sgs_tke(SHOCAdvsgstkeData &d);
void eddy_diffusivities(SHOCEddydiffData &d);
void compute_brunt_shoc_length(SHOCBruntlengthData &d);
void compute_l_inf_shoc_length(SHOCInflengthData &d);
void compute_conv_vel_shoc_length(SHOCConvvelData &d);
void compute_conv_time_shoc_length(SHOCConvtimeData &d);
void compute_shoc_mix_shoc_length(SHOCMixlengthData &d);
void check_length_scale_shoc_length(SHOCMixcheckData &d);
void shoc_diag_second_moments_srf(SHOCSecondMomentSrfData& d);
void linear_interp(SHOCLinearintData &d);
void shoc_assumed_pdf_tilda_to_real(SHOCPDFtildaData &d);
void shoc_assumed_pdf_vv_parameters(SHOCPDFvvparamData &d);
void shoc_assumed_pdf_thl_parameters(SHOCPDFthlparamData &d);
void shoc_assumed_pdf_qw_parameters(SHOCPDFqwparamData &d);
void shoc_assumed_pdf_inplume_correlations(SHOCPDFinplumeData &d);
void shoc_assumed_pdf_compute_temperature(SHOCPDFcomptempData &d);
void shoc_assumed_pdf_compute_qs(SHOCPDFcompqsData &d);
void shoc_assumed_pdf_compute_sgs_liquid(SHOCPDFcompsgsliqData &d);
void shoc_assumed_pdf_compute_cloud_liquid_variance(SHOCPDFcompcloudvarData &d);
void shoc_assumed_pdf_compute_liquid_water_flux(SHOCPDFcompliqfluxData &d);
void shoc_assumed_pdf_compute_buoyancy_flux(SHOCPDFcompbuoyfluxData &d);

//
// _f functions decls
//
extern "C" {

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
void shoc_diag_second_moments_srf_f(Int shcol, Real* wthl, Real* uw, Real* vw,
                         Real* ustar2, Real* wstar);

}

}  // namespace shoc
}  // namespace scream

#endif
