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

void calc_shoc_varorcovar_c(Int shcol, Int nlev, Int nlevi,  Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
			    Real *invar1, Real *invar2, Real *varorcovar);

void calc_shoc_vertflux_c(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
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


//Initialize data for calc_shoc_vertflux function
SHOCVertfluxData::SHOCVertfluxData(Int shcol_, Int nlev_, Int nlevi_)
  : shcol(shcol_),
    nlev(nlev_),
    nlevi(nlevi_),
    m_total(shcol_ * nlev_),
    m_totali(shcol_ * nlevi_),
    m_data(NUM_ARRAYS * m_total, 0),
    m_datai(NUM_ARRAYS_i * m_totali, 0) {
  init_ptrs();
}

SHOCVertfluxData::SHOCVertfluxData(const SHOCVertfluxData &rhs)
  : shcol(rhs.shcol),
    nlev(rhs.nlev),
    nlevi(rhs.nlevi),
    m_total(rhs.m_total),
    m_totali(rhs.m_totali),
    m_data(rhs.m_data),
    m_datai(rhs.m_datai) {
  init_ptrs();
}


SHOCVertfluxData  &SHOCVertfluxData::operator=(const SHOCVertfluxData &rhs) {
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


void SHOCVertfluxData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_i = m_datai.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&invar};
  std::array<Real **, NUM_ARRAYS_i> ptrs_i = {&tkh_zi, &dz_zi, &vertflux};

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

SHOCVarorcovarData::SHOCVarorcovarData(Int shcol_, Int nlev_, Int nlevi_, Real tunefac_)
  : shcol(shcol_),
    nlev(nlev_),
    nlevi(nlevi_),
    tunefac(tunefac_),
    m_total(shcol_ * nlev_),
    m_totali(shcol_ * nlevi_),
    m_data(NUM_ARRAYS * m_total, 0),
    m_datai(NUM_ARRAYS_i * m_totali, 0) {
  init_ptrs();
}

SHOCVarorcovarData::SHOCVarorcovarData(const SHOCVarorcovarData &rhs)
  : shcol(rhs.shcol),
    nlev(rhs.nlev),
    nlevi(rhs.nlevi),
    tunefac(rhs.tunefac),
    m_total(rhs.m_total),
    m_totali(rhs.m_totali),
    m_data(rhs.m_data),
    m_datai(rhs.m_datai) {
  init_ptrs();
}


SHOCVarorcovarData  &SHOCVarorcovarData::operator=(const SHOCVarorcovarData &rhs) {
  init_ptrs();

  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  nlevi    = rhs.nlevi;
  tunefac  = rhs.tunefac;
  m_total  = rhs.m_total;
  m_totali = rhs.m_totali;
  m_data   = rhs.m_data;
  m_datai  = rhs.m_datai;  // Copy

  return *this;
}


void SHOCVarorcovarData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_i = m_datai.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&invar1, &invar2};
  std::array<Real **, NUM_ARRAYS_i> ptrs_i = {&tkh_zi, &dz_zi, &isotropy_zi, &varorcovar};

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

void calc_shoc_varorcovar(Int nlev, SHOCVarorcovarData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_varorcovar_c(d.shcol, d.nlev, d.nlevi, d.tunefac, d.isotropy_zi, d.tkh_zi,
                         d.dz_zi, d.invar1, d.invar2, d.varorcovar);
  d.transpose<util::TransposeDirection::f2c>();
}

//
// Glue functions to call fortran from from C++ with the Data struct
//

//Initialize shoc parameterization, trnaspose data from c to fortran,
//call calc_shoc_vertflux fortran subroutine and transpose data back to c
void calc_shoc_vertflux(SHOCVertfluxData &d) {
  shoc_init(d.nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<util::TransposeDirection::f2c>();
}

//
// _f function definitions
//
extern "C" {

void calc_shoc_vertflux_f(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux)
{
  using SHF = Functions<Real, DefaultDevice>;

  using Spack      = typename SHF::Spack;
  using view_2d    = typename SHF::view_2d<Spack>;
  using KT         = typename SHF::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename SHF::MemberType;

  static constexpr Int num_arrays = SHOCVertfluxData::NUM_ARRAYS + SHOCVertfluxData::NUM_ARRAYS_i;

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

  auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    SHF::calc_shoc_vertflux(team, shcol, nlev, nlevi, tkh_zi_d, dz_zi_d, invar_d, vertflux_d);
  });

  // Sync back to host
  Kokkos::Array<view_2d, 1> inout_views = {vertflux_d};
  pack::device_to_host({vertflux}, {shcol}, {nlevi}, inout_views, true);
}

}

} // namespace shoc
} // namespace scream
