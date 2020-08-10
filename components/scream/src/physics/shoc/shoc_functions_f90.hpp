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

///////////////////////////////////////////////////////////////////////////////
// Converted subroutine helpers go here.
struct SHOCGridData {
  static constexpr size_t NUM_ARRAYS   = 4;
  static constexpr size_t NUM_ARRAYS_i = 2;

  // Inputs
  Int shcol, nlev, nlevi;
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int shcol_, Int nevl_, Int nlevi_);
  SHOCGridData(const SHOCGridData &rhs);
  SHOCGridData &operator=(const SHOCGridData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevi, m_total, m_totali;
  std::vector<Real> m_data;
  std::vector<Real> m_datai;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCGridData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(zt_grid, d_trans.zt_grid, shcol, nlev);
    util::transpose<D>(dz_zt, d_trans.dz_zt, shcol, nlev);
    util::transpose<D>(pdel, d_trans.pdel, shcol, nlev);
    util::transpose<D>(rho_zt, d_trans.rho_zt, shcol, nlev);

    // Transpose on the zi grid
    util::transpose<D>(zi_grid, d_trans.zi_grid, shcol, nlevi);
    util::transpose<D>(dz_zi, d_trans.dz_zi, shcol, nlevi);

    *this = std::move(d_trans);
  }
};

// This function initialzes the grid used by shoc. Given the 
// locations of the cell center (location of thermodynaics quantities), cell 
// interfaces, and pressure gradient the functon returns dz_zi, dz_zt,
// and density. 
void shoc_grid(Int nlev, SHOCGridData &d);


//Create data structure to hold data for integ_column_stability
struct SHOCColstabData {
  static constexpr size_t NUM_ARRAYS   = 3; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 1; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *dz_zt, *pres, *brunt;

  // Output
  Real *brunt_int;

  //functions to initialize data
  SHOCColstabData(Int shcol_, Int nlev_);
  SHOCColstabData(const SHOCColstabData &rhs);
  SHOCColstabData &operator=(const SHOCColstabData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCColstabData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(dz_zt, d_trans.dz_zt, shcol, nlev);
    util::transpose<D>(pres, d_trans.pres, shcol, nlev);
    util::transpose<D>(brunt, d_trans.brunt, shcol, nlev);
    
    // Transpose on the column only grid
    util::transpose<D>(brunt_int, d_trans.brunt_int, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCColstabData

void integ_column_stability(Int nlev, SHOCColstabData &d);

//Create data structure to hold data for calc_shoc_vertflux
struct SHOCTkeshearData {
  static constexpr size_t NUM_ARRAYS   = 2; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_i = 2; //# of arrays with values at interface centers (zi grid)

  // Inputs
  Int   shcol, nlev, nlevi;
  Real *dz_zi, *u_wind, *v_wind;

  // In/out
  Real *sterm;

  //functions to initialize data
  SHOCTkeshearData(Int shcol_, Int nlev_, Int nlevi_);
  SHOCTkeshearData(const SHOCTkeshearData &rhs);
  SHOCTkeshearData &operator=(const SHOCTkeshearData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevi, m_total, m_totali;
  std::vector<Real> m_data;
  std::vector<Real> m_datai;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCTkeshearData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(u_wind, d_trans.u_wind, shcol, nlev);
    util::transpose<D>(v_wind, d_trans.v_wind, shcol, nlev);

    // Transpose on the zi grid
    util::transpose<D>(dz_zi, d_trans.dz_zi, shcol, nlevi);
    util::transpose<D>(sterm, d_trans.sterm, shcol, nlevi);

    *this = std::move(d_trans);
  }
};//SHOCTkeshearData

void compute_shr_prod(Int nlev, SHOCTkeshearData &d);

//Create data structure to hold data for integ_column_stability
struct SHOCIsotropicData {
  static constexpr size_t NUM_ARRAYS   = 4; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 1; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *tke, *a_diss, *brunt, *brunt_int;

  // Output
  Real *isotropy;

  //functions to initialize data
  SHOCIsotropicData(Int shcol_, Int nlev_);
  SHOCIsotropicData(const SHOCIsotropicData &rhs);
  SHOCIsotropicData &operator=(const SHOCIsotropicData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCIsotropicData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(tke, d_trans.tke, shcol, nlev);
    util::transpose<D>(a_diss, d_trans.a_diss, shcol, nlev);
    util::transpose<D>(brunt, d_trans.brunt, shcol, nlev);
    util::transpose<D>(isotropy, d_trans.isotropy, shcol, nlev);
    
    // Transpose on the column only grid
    util::transpose<D>(brunt_int, d_trans.brunt_int, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCIsotropicData

void isotropic_ts(Int nlev, SHOCIsotropicData &d);

//Create data structure to hold data for adv_sgs_tke
struct SHOCAdvsgstkeData {
  static constexpr size_t NUM_ARRAYS   = 7; //# of arrays with values at cell centers (zt grid)

  // Inputs
  Int   shcol, nlev;
  Real dtime;
  Real *shoc_mix, *wthv_sec, *sterm_zt, *tk;

  // In/out
  Real *tke;
  
  // Outputs
  Real *a_diss;

  //functions to initialize data
  SHOCAdvsgstkeData(Int shcol_, Int nlev_, Real dtime_);
  SHOCAdvsgstkeData(const SHOCAdvsgstkeData &rhs);
  SHOCAdvsgstkeData &operator=(const SHOCAdvsgstkeData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_total;
  std::vector<Real> m_data;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCAdvsgstkeData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(shoc_mix, d_trans.shoc_mix, shcol, nlev);
    util::transpose<D>(wthv_sec, d_trans.wthv_sec, shcol, nlev);
    util::transpose<D>(sterm_zt, d_trans.sterm_zt, shcol, nlev);
    util::transpose<D>(tk, d_trans.tk, shcol, nlev);
    util::transpose<D>(tke, d_trans.tke, shcol, nlev);
    util::transpose<D>(a_diss, d_trans.a_diss, shcol, nlev);

    *this = std::move(d_trans);
  }
};//SHOCAdvsgstkeData

void adv_sgs_tke(Int nlev, SHOCAdvsgstkeData &d);

//Create data structure to hold data for eddy_diffusivities
struct SHOCEddydiffData {
  static constexpr size_t NUM_ARRAYS   = 7; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 2; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *pblh, *obklen, *zt_grid, *shoc_mix, *sterm_zt, 
        *isotropy, *tke;

  // Output
  Real *tk, *tkh;

  //functions to initialize data
  SHOCEddydiffData(Int shcol_, Int nlev_);
  SHOCEddydiffData(const SHOCEddydiffData &rhs);
  SHOCEddydiffData &operator=(const SHOCEddydiffData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCEddydiffData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(zt_grid, d_trans.zt_grid, shcol, nlev);
    util::transpose<D>(shoc_mix, d_trans.shoc_mix, shcol, nlev);
    util::transpose<D>(isotropy, d_trans.isotropy, shcol, nlev);
    util::transpose<D>(tke, d_trans.tke, shcol, nlev);
    util::transpose<D>(tk, d_trans.tk, shcol, nlev);
    util::transpose<D>(tkh, d_trans.tkh, shcol, nlev);
    
    // Transpose on the column only grid
    util::transpose<D>(obklen, d_trans.obklen, shcol, 1);
    util::transpose<D>(pblh, d_trans.pblh, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCEddydiffData

void eddy_diffusivities(Int nlev, SHOCEddydiffData &d);


//Create data structure to hold data for calc_shoc_vertflux
struct SHOCVertfluxData {
  static constexpr size_t NUM_ARRAYS   = 1; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_i = 3; //# of arrays with values at interface centers (zi grid)

  // Inputs
  Int   shcol, nlev, nlevi;
  Real *tkh_zi, *dz_zi, *invar;

  // In/out
  Real *vertflux;

  //functions to initialize data
  SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_);
  SHOCVertfluxData(const SHOCVertfluxData &rhs);
  SHOCVertfluxData &operator=(const SHOCVertfluxData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevi, m_total, m_totali;
  std::vector<Real> m_data;
  std::vector<Real> m_datai;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCVertfluxData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(invar, d_trans.invar, shcol, nlev);

    // Transpose on the zi grid
    util::transpose<D>(tkh_zi, d_trans.tkh_zi, shcol, nlevi);
    util::transpose<D>(dz_zi, d_trans.dz_zi, shcol, nlevi);
    util::transpose<D>(vertflux, d_trans.vertflux, shcol, nlevi);

    *this = std::move(d_trans);
  }
};//SHOCVertfluxData

void calc_shoc_vertflux(Int nlev, SHOCVertfluxData &d);

struct SHOCVarorcovarData {
  static constexpr size_t NUM_ARRAYS   = 2;
  static constexpr size_t NUM_ARRAYS_i = 4;

  // Inputs
  Int   shcol, nlev, nlevi;
  Real tunefac;
  Real *tkh_zi, *dz_zi, *isotropy_zi, *invar1, *invar2;

  // In/out
  Real *varorcovar;

  SHOCVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_);
  SHOCVarorcovarData(const SHOCVarorcovarData &rhs);
  SHOCVarorcovarData &operator=(const SHOCVarorcovarData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevi, m_total, m_totali;
  std::vector<Real> m_data;
  std::vector<Real> m_datai;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCVarorcovarData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(invar1, d_trans.invar1, shcol, nlev);
    util::transpose<D>(invar2, d_trans.invar2, shcol, nlev);

    // Transpose on the zi grid
    util::transpose<D>(tkh_zi, d_trans.tkh_zi, shcol, nlevi);
    util::transpose<D>(dz_zi, d_trans.dz_zi, shcol, nlevi);
    util::transpose<D>(isotropy_zi, d_trans.isotropy_zi, shcol, nlevi);
    util::transpose<D>(varorcovar, d_trans.varorcovar, shcol, nlevi);

    *this = std::move(d_trans);
  }
};//SHOCVarorcovarData

void calc_shoc_varorcovar(Int nlev, SHOCVarorcovarData &d);

//Create data structure to hold data for compute_brunt_shoc_length
struct SHOCBruntlengthData {
  static constexpr size_t NUM_ARRAYS   = 3; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_i = 1; //# of arrays with values at interface centers (zi grid)

  // Inputs
  Int   shcol, nlev, nlevi;
  Real *dz_zt, *thv, *thv_zi;

  // Outputs
  Real *brunt;

  //functions to initialize data
  SHOCBruntlengthData(Int shcol_, Int nlev_, Int nlevi_);
  SHOCBruntlengthData(const SHOCBruntlengthData &rhs);
  SHOCBruntlengthData &operator=(const SHOCBruntlengthData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevi, m_total, m_totali;
  std::vector<Real> m_data;
  std::vector<Real> m_datai;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCBruntlengthData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(dz_zt, d_trans.dz_zt, shcol, nlev);
    util::transpose<D>(thv, d_trans.thv, shcol, nlev);
    util::transpose<D>(brunt, d_trans.brunt, shcol, nlev);

    // Transpose on the zi grid
    util::transpose<D>(thv_zi, d_trans.thv_zi, shcol, nlevi);

    *this = std::move(d_trans);
  }
};//SHOCBruntlengthData

void compute_brunt_shoc_length(Int nlev, SHOCBruntlengthData &d);

//Create data structure to hold data for compute_l_inf_shoc_length
struct SHOCInflengthData {
  static constexpr size_t NUM_ARRAYS   = 3; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 1; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *zt_grid, *dz_zt, *tke;

  // Outputs
  Real *l_inf;

  //functions to initialize data
  SHOCInflengthData(Int shcol_, Int nlev_);
  SHOCInflengthData(const SHOCInflengthData &rhs);
  SHOCInflengthData &operator=(const SHOCInflengthData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevc, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCInflengthData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(zt_grid, d_trans.zt_grid, shcol, nlev);
    util::transpose<D>(dz_zt, d_trans.dz_zt, shcol, nlev);
    util::transpose<D>(tke, d_trans.tke, shcol, nlev);

    // Transpose on the column only grid
    util::transpose<D>(l_inf, d_trans.l_inf, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCInflengthData

void compute_l_inf_shoc_length(Int nlev, SHOCInflengthData &d);

//Create data structure to hold data for compute_conv_vel_shoc_length
struct SHOCConvvelData {
  static constexpr size_t NUM_ARRAYS   = 4; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 2; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *pblh, *zt_grid, *dz_zt, *thv, *wthv_sec;

  // Outputs
  Real *conv_vel;

  //functions to initialize data
  SHOCConvvelData(Int shcol_, Int nlev_);
  SHOCConvvelData(const SHOCConvvelData &rhs);
  SHOCConvvelData &operator=(const SHOCConvvelData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevc, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCConvvelData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(zt_grid, d_trans.zt_grid, shcol, nlev);
    util::transpose<D>(dz_zt, d_trans.dz_zt, shcol, nlev);
    util::transpose<D>(thv, d_trans.thv, shcol, nlev);
    util::transpose<D>(wthv_sec, d_trans.wthv_sec, shcol, nlev);

    // Transpose on the column only grid
    util::transpose<D>(pblh, d_trans.pblh, shcol, 1);
    util::transpose<D>(conv_vel, d_trans.conv_vel, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCConvvelData

void compute_conv_vel_shoc_length(Int nlev, SHOCConvvelData &d);

//Create data structure to hold data for compute_conv_time_shoc_length
struct SHOCConvtimeData {
  static constexpr size_t NUM_ARRAYS_c = 3; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *pblh, *conv_vel;

  // Outputs
  Real *tscale;

  //functions to initialize data
  SHOCConvtimeData(Int shcol_);
  SHOCConvtimeData(const SHOCConvtimeData &rhs);
  SHOCConvtimeData &operator=(const SHOCConvtimeData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_totalc;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCConvtimeData d_trans(*this);

    // Transpose on the column only grid
    util::transpose<D>(pblh, d_trans.pblh, shcol, 1);
    util::transpose<D>(conv_vel, d_trans.conv_vel, shcol, 1);
    util::transpose<D>(tscale, d_trans.tscale, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCConvtimeData

void compute_conv_time_shoc_length(Int nlev, SHOCConvtimeData &d);

//Create data structure to hold data for compute_shoc_mix_shoc_length
struct SHOCMixlengthData {
  static constexpr size_t NUM_ARRAYS   = 4; //# of arrays with values at cell centers (zt grid)
  static constexpr size_t NUM_ARRAYS_c = 2; //# of arrays with column only dimensions

  // Inputs
  Int   shcol, nlev;
  Real *tke, *brunt, *tscale, *zt_grid, *l_inf;

  // Outputs
  Real *shoc_mix;

  //functions to initialize data
  SHOCMixlengthData(Int shcol_, Int nlev_);
  SHOCMixlengthData(const SHOCMixlengthData &rhs);
  SHOCMixlengthData &operator=(const SHOCMixlengthData &rhs);

  void init_ptrs();

  // Internals
  Int m_shcol, m_nlev, m_nlevc, m_total, m_totalc;
  std::vector<Real> m_data;
  std::vector<Real> m_datac;

  template <util::TransposeDirection::Enum D>
  void transpose() {
    SHOCMixlengthData d_trans(*this);

    // Transpose on the zt grid
    util::transpose<D>(tke, d_trans.tke, shcol, nlev);
    util::transpose<D>(brunt, d_trans.brunt, shcol, nlev);
    util::transpose<D>(zt_grid, d_trans.zt_grid, shcol, nlev);
    util::transpose<D>(shoc_mix, d_trans.shoc_mix, shcol, nlev);

    // Transpose on the column only grid
    util::transpose<D>(l_inf, d_trans.l_inf, shcol, 1);
    util::transpose<D>(tscale, d_trans.tscale, shcol, 1);

    *this = std::move(d_trans);
  }
};//SHOCMixlengthData

void compute_shoc_mix_shoc_length(Int nlev, SHOCMixlengthData &d);

}  // namespace shoc
}  // namespace scream

#endif
