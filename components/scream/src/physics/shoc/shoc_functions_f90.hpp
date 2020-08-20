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
               const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i);

  SHOCDataBase(const SHOCDataBase &rhs) = delete;
  SHOCDataBase(const SHOCDataBase &rhs, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i);
  SHOCDataBase &operator=(const SHOCDataBase &rhs);

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

    m_data = data;
  }

  void randomize();

  Int total() const { return m_total; }
  Int totali() const { return m_totali; }

 private:
  void init_ptrs();

  // Internals
  Int m_total, m_totali;
  std::vector<Real**> m_ptrs, m_ptrs_i;
  std::vector<Real> m_data;
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
    SHOCDataBase(shcol_, nlev_, 1, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}
  SHOCColstabData(const SHOCColstabData &rhs) : SHOCDataBase(rhs, {&dz_zt, &pres, &brunt}, {&brunt_int}) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {&tke, &a_diss, &brunt, &isotropy}, {&brunt_int}) {}
  SHOCIsotropicData(const SHOCIsotropicData &rhs) : SHOCDataBase(rhs, {&tke, &a_diss, &brunt, &isotropy}, {&brunt_int}) {}
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
    SHOCDataBase(shcol_, nlev_, 0, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}, {}), dtime(dtime_) {}
  SHOCAdvsgstkeData(const SHOCAdvsgstkeData &rhs) : SHOCDataBase(rhs, {&shoc_mix, &wthv_sec, &sterm_zt, &tk, &tke, &a_diss}, {}), dtime(rhs.dtime) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {&obklen, &pblh}) {}
  SHOCEddydiffData(const SHOCEddydiffData &rhs) : SHOCDataBase(rhs, {&zt_grid, &shoc_mix, &isotropy, &tke, &tk, &tkh, &sterm_zt}, {&obklen, &pblh}) {}
  SHOCEddydiffData &operator=(const SHOCEddydiffData &rhs) { SHOCDataBase::operator=(rhs); return *this; }
};//SHOCEddydiffData

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
    SHOCDataBase(shcol_, nlev_, 1, {&zt_grid, &dz_zt, &tke}, {&l_inf}) {}
  SHOCInflengthData(const SHOCInflengthData &rhs) :
    SHOCDataBase(rhs, {&zt_grid, &dz_zt, &tke}, {&l_inf}) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {&conv_vel, &pblh}) {}
  SHOCConvvelData(const SHOCConvvelData &rhs) :
    SHOCDataBase(rhs, {&zt_grid, &dz_zt, &thv, &wthv_sec}, {&conv_vel, &pblh}) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {}, {&conv_vel, &pblh, &tscale}) {}
  SHOCConvtimeData(const SHOCConvtimeData &rhs) :
    SHOCDataBase(rhs, {}, {&conv_vel, &pblh, &tscale}) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {&tke, &brunt, &zt_grid, &shoc_mix}, {&l_inf, &tscale}) {}
  SHOCMixlengthData(const SHOCMixlengthData &rhs) :
    SHOCDataBase(rhs, {&tke, &brunt, &zt_grid, &shoc_mix}, {&l_inf, &tscale}) {}
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
    SHOCDataBase(shcol_, nlev_, 1, {&shoc_mix}, {&host_dx, &host_dy}) {}
  SHOCMixcheckData(const SHOCMixcheckData &rhs) :
    SHOCDataBase(rhs, {&shoc_mix}, {&host_dx, &host_dy}) {}
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

//
// Glue functions to call fortran from from C++ with the Data struct
//

// This function initialzes the grid used by shoc. Given the
// locations of the cell center (location of thermodynaics quantities), cell
// interfaces, and pressure gradient the functon returns dz_zi, dz_zt,
// and density.
void shoc_grid(SHOCGridData &d);
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
