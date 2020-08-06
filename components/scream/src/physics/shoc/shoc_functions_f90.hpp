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

}  // namespace shoc
}  // namespace scream

#endif
