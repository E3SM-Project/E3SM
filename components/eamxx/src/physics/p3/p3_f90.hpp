#ifndef SCREAM_P3_F90_HPP
#define SCREAM_P3_F90_HPP

#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace p3 {

// Data format we can use to communicate with Fortran version.
struct FortranData {
  typedef std::shared_ptr<FortranData> Ptr;

  using KT     = KokkosTypes<HostDevice>;
  using Scalar = Real;

  using Array1 = typename KT::template lview<Scalar*>;
  using Array2 = typename KT::template lview<Scalar**>;
  using Array3 = typename KT::template lview<Scalar***>;

  bool do_predict_nc;
  bool do_prescribed_CCN;
  const Int ncol, nlev;

  // In
  Real dt;
  Int it;
  Array2 qv, th_atm, pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated, inv_qc_relvar, qc, nc, qr, nr,  qi,
    ni, qm, bm, dpres, inv_exner, qv_prev, t_prev;
  // Out
  Array1 precip_liq_surf, precip_ice_surf;
  Array2 diag_eff_radius_qc, diag_eff_radius_qi, rho_qi, qv2qi_depos_tend,
         precip_liq_flux, precip_ice_flux, cld_frac_r, cld_frac_l, cld_frac_i;
  Array3 p3_tend_out;
  Array2 liq_ice_exchange,vap_liq_exchange,vap_ice_exchange;

  FortranData(Int ncol, Int nlev);
};

// Iterate over a FortranData's arrays. For examples, see Baseline::write, read.
struct FortranDataIterator {
  struct RawArray {
    std::string name;
    Int dim;
    Int extent[3];
    FortranData::Scalar* data;
    FortranData::Array1::size_type size;
  };

  explicit FortranDataIterator(const FortranData::Ptr& d);

  Int nfield () const { return fields_.size(); }
  const RawArray& getfield(Int i) const;

private:
  FortranData::Ptr d_;
  std::vector<RawArray> fields_;

  void init(const FortranData::Ptr& d);
};

void p3_init(const bool write_tables = false,
             const bool masterproc = false);

// We will likely want to remove these checks in the future, as we're not tied
// to the exact implementation or arithmetic in P3. For now, these checks are
// here to establish that the initial regression-testing code gives results that
// match the python f2py tester, without needing a data file.
Int check_against_python(const FortranData& d);

int test_FortranData();

}  // namespace p3
}  // namespace scream

#endif
