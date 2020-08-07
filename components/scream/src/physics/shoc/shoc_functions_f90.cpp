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
		 
void calc_shoc_varorcovar_c(Int shcol, Int nlev, Int nlevi,  Real tunefac,
                            Real *isotropy_zi, Real *tkh_zi, Real *dz_zi,
			    Real *invar1, Real *invar2, Real *varorcovar);		 

void calc_shoc_vertflux_c(Int shcol, Int nlev, Int nlevi, Real *tkh_zi,
			  Real *dz_zi, Real *invar, Real *vertflux);
			  
void compute_brunt_shoc_length_c(Int nlev, Int nlevi, Int shcol ,Real *dz_zt,
                                 Real *thv, Real *thv_zi, Real *brunt);
				 
void compute_l_inf_shoc_length_c(Int nlev, Int shcol, Real *zt_grid, Real *dz_zt,
                                 Real *tke, Real *l_inf);	
				 
void compute_conv_vel_shoc_length_c(Int nlev, Int shcol, Real *pblh, Real *zt_grid,
                                    Real *dz_zt, Real *thv, Real *wthv_sec,
				    Real *conv_vel);
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

//Initialize shoc parameterization, trnaspose data from c to fortran,
//call calc_shoc_vertflux fortran subroutine and transpose data back to c
void calc_shoc_vertflux(Int nlev, SHOCVertfluxData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  calc_shoc_vertflux_c(d.shcol, d.nlev, d.nlevi, d.tkh_zi, d.dz_zi, d.invar,
		       d.vertflux);
  d.transpose<util::TransposeDirection::f2c>();
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

//Initialize data for compute_brunt_shoc_length function
SHOCBruntlengthData::SHOCBruntlengthData(Int shcol_, Int nlev_, Int nlevi_)
  : shcol(shcol_),
    nlev(nlev_),
    nlevi(nlevi_),
    m_total(shcol_ * nlev_),
    m_totali(shcol_ * nlevi_),
    m_data(NUM_ARRAYS * m_total, 0),
    m_datai(NUM_ARRAYS_i * m_totali, 0) {
  init_ptrs();
}

SHOCBruntlengthData::SHOCBruntlengthData(const SHOCBruntlengthData &rhs)
  : shcol(rhs.shcol),
    nlev(rhs.nlev),
    nlevi(rhs.nlevi),
    m_total(rhs.m_total),
    m_totali(rhs.m_totali),
    m_data(rhs.m_data),
    m_datai(rhs.m_datai) {
  init_ptrs();
}


SHOCBruntlengthData  &SHOCBruntlengthData::operator=(const SHOCBruntlengthData &rhs) {
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


void SHOCBruntlengthData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_i = m_datai.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&dz_zt, &thv, &brunt};
  std::array<Real **, NUM_ARRAYS_i> ptrs_i = {&thv_zi};

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

//Initialize shoc parameterization, trnaspose data from c to fortran,
//call compute_brunt_shoc_length fortran subroutine and transpose data back to c
void compute_brunt_shoc_length(Int nlev, SHOCBruntlengthData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  compute_brunt_shoc_length_c(d.nlev,d.nlevi,d.shcol,d.dz_zt,d.thv,d.thv_zi,d.brunt);
  d.transpose<util::TransposeDirection::f2c>();
}


//Initialize data for compute_l_inf_shoc_length function
SHOCInflengthData::SHOCInflengthData(Int shcol_, Int nlev_)
  : shcol(shcol_),
    nlev(nlev_),
    m_total(shcol_ * nlev_),
    m_totalc(shcol_),
    m_data(NUM_ARRAYS * m_total, 0),
    m_datac(NUM_ARRAYS_c * m_totalc,0) {
  init_ptrs();
}

SHOCInflengthData::SHOCInflengthData(const SHOCInflengthData &rhs)
  : shcol(rhs.shcol),
    nlev(rhs.nlev),
    m_total(rhs.m_total),
    m_totalc(rhs.m_totalc),
    m_data(rhs.m_data),
    m_datac(rhs.m_datac) {
  init_ptrs();
}


SHOCInflengthData  &SHOCInflengthData::operator=(const SHOCInflengthData &rhs) {
  init_ptrs();

  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  m_total  = rhs.m_total;
  m_totalc = rhs.m_totalc;
  m_data   = rhs.m_data;
  m_datac  = rhs.m_datac;  // Copy

  return *this;
}


void SHOCInflengthData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_c = m_datac.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&zt_grid, &dz_zt, &tke};
  std::array<Real **, NUM_ARRAYS_c> ptrs_c = {&l_inf};

  for(size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_total;
  }

  offset = 0;
  for(size_t i = 0; i < NUM_ARRAYS_c; ++i) {
    *ptrs_c[i] = data_begin_c + offset;
    offset += m_totalc;
  }
}

//Initialize shoc parameterization, trnaspose data from c to fortran,
//call compute_l_inf_shoc_length fortran subroutine and transpose data back to c
void compute_l_inf_shoc_length(Int nlev, SHOCInflengthData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  compute_l_inf_shoc_length_c(d.nlev,d.shcol,d.zt_grid,d.dz_zt,d.tke,d.l_inf);
  d.transpose<util::TransposeDirection::f2c>();
}

//Initialize data for compute_l_inf_shoc_length function
SHOCConvvelData::SHOCConvvelData(Int shcol_, Int nlev_)
  : shcol(shcol_),
    nlev(nlev_),
    m_total(shcol_ * nlev_),
    m_totalc(shcol_),
    m_data(NUM_ARRAYS * m_total, 0),
    m_datac(NUM_ARRAYS_c * m_totalc,0) {
  init_ptrs();
}

SHOCConvvelData::SHOCConvvelData(const SHOCConvvelData &rhs)
  : shcol(rhs.shcol),
    nlev(rhs.nlev),
    m_total(rhs.m_total),
    m_totalc(rhs.m_totalc),
    m_data(rhs.m_data),
    m_datac(rhs.m_datac) {
  init_ptrs();
}


SHOCConvvelData  &SHOCConvvelData::operator=(const SHOCConvvelData &rhs) {
  init_ptrs();

  shcol    = rhs.shcol;
  nlev     = rhs.nlev;
  m_total  = rhs.m_total;
  m_totalc = rhs.m_totalc;
  m_data   = rhs.m_data;
  m_datac  = rhs.m_datac;  // Copy

  return *this;
}


void SHOCConvvelData::init_ptrs() {
  Int offset         = 0;
  Real *data_begin   = m_data.data();
  Real *data_begin_c = m_datac.data();

  std::array<Real **, NUM_ARRAYS> ptrs     = {&zt_grid, &dz_zt, &thv, &wthv_sec};
  std::array<Real **, NUM_ARRAYS_c> ptrs_c = {&pblh, &conv_vel};

  for(size_t i = 0; i < NUM_ARRAYS; ++i) {
    *ptrs[i] = data_begin + offset;
    offset += m_total;
  }

  offset = 0;
  for(size_t i = 0; i < NUM_ARRAYS_c; ++i) {
    *ptrs_c[i] = data_begin_c + offset;
    offset += m_totalc;
  }
}

//Initialize shoc parameterization, trnaspose data from c to fortran,
//call compute_conv_vel_shoc_length fortran subroutine and transpose data back to c
void compute_conv_vel_shoc_length(Int nlev, SHOCConvvelData &d) {
  shoc_init(nlev, true);
  d.transpose<util::TransposeDirection::c2f>();
  compute_conv_vel_shoc_length_c(d.nlev,d.shcol,d.pblh,d.zt_grid,
                                 d.dz_zt,d.thv,d.wthv_sec,d.conv_vel);
  d.transpose<util::TransposeDirection::f2c>();
}

} // namespace shoc
} // namespace scream
