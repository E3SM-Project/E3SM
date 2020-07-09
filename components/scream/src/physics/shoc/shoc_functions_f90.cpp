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
// A C++ interface to SHOC fortran calls and vice versa
//

extern "C" {

void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                 Real zvir, Real latvap, Real latice, Real karman,
                 Real* pref_mid, int nbot_shoc, int ntop_shoc);

void shoc_grid_c(int shcol, int nlev, int nlevi, Real *zt_grid, Real *zi_grid,
                 Real *pdel, Real *dz_zt, Real *dzi_zi, Real *rho_zt);
}

namespace scream {
namespace shoc {

// helper functions

// In all C++ -> Fortran bridge functions you should see shoc_init(nlev, true).
// We are provisionally following P3 here in case SHOC uses global data. The
// 'true' argument is to set shoc to use its fortran implementations instead of
// calling back to C++. We want this behavior since it doesn't make much sense
// for C++ to bridge over to fortran only to have fortran bridge back to C++.
// Anyone who wants the C++ implementation should call it directly.

struct SHOCSubroutineData  // example data struct
{
  // In
  Real in1, in2, in3;

  // Out
  Real out1, out2, out3;
};

void shoc_subroutine(SHOCSubroutineData &d)  // example wrapper function
{
  Int nlev = 128;
  shoc_init(nlev, true);
  // shoc_subroutine_c(d.in1, d.in2, d.in3, &d.out1, &d.out2, &d.out3);
}

SHOCGridData::SHOCGridData(Int shcol_, Int nlev_, Int nlevi_)
    : shcol(shcol_),
      nlev(nlev_),
      nlevi(nlevi_),
      m_total(shcol_ * nlev_),
      m_totali(shcol_ * nlevi_),
      m_data(NUM_ARRAYS * m_total, 0),
      m_datai(NUM_ARRAYS_i * m_totali, 0) {
  init_ptrs();
}

SHOCGridData::SHOCGridData(const SHOCGridData &rhs)
    : shcol(rhs.shcol),
      nlev(rhs.nlev),
      nlevi(rhs.nlevi),
      m_total(rhs.m_total),
      m_totali(rhs.m_totali),
      m_data(rhs.m_data),
      m_datai(rhs.m_datai) {
  init_ptrs();
}

SHOCGridData &SHOCGridData::operator=(const SHOCGridData &rhs) {
  init_ptrs();

  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  nlevi    = rhs.nlevi;
  m_total  = rhs.m_total;
  m_totali = rhs.m_totali;
  m_data   = rhs.m_data;
  m_datai  = rhs.m_datai;  // Copy

  return *this;
}

void SHOCGridData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_i = m_datai.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&zt_grid, &dz_zt, &pdel, &rho_zt};
  std::array<Real **, NUM_ARRAYS_i> ptrs_i = {&zi_grid, &dz_zi};

  for(size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_total;
  }

  offset = 0;
  for(size_t i = 0; i < NUM_ARRAYS_i; ++i) {
    *ptrs_i[i] = data_begin_i + offset;
    offset += m_totali;
  }
}

void shoc_grid(Int nlev, SHOCGridData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  shoc_grid_c(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt,
              d.dz_zi, d.rho_zt);
  d.transpose<util::TransposeDirection::f2c>();
}

} // namespace shoc
} // namespace scream
