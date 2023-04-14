#ifndef SCREAM_SHOC_F90_HPP
#define SCREAM_SHOC_F90_HPP

#include "share/scream_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace shoc {

// Data format we can use to communicate with Fortran version.
struct FortranData {
  typedef std::shared_ptr<FortranData> Ptr;

  using KT     = KokkosTypes<HostDevice>;
  using Scalar = Real;

  using Array1 = typename KT::template lview<Scalar*>;
  using Array2 = typename KT::template lview<Scalar**>;
  using Array3 = typename KT::template lview<Scalar***>;

  Int shcol, nlev, nlevi, num_qtracers, nadv;

  // In
  Real dtime;
  Array1 host_dx, host_dy;
  Array2 zt_grid, zi_grid, pres, presi, pdel, thv, w_field;
  Array1 wthl_sfc, wqw_sfc, uw_sfc, vw_sfc;
  Array2 wtracer_sfc, inv_exner;
  Array1 phis;

  // In-out
  Array2 host_dse, tke, thetal, qw, u_wind, v_wind, wthv_sec;
  Array3 qtracers;
  Array2 tk, tkh;

  // Out
  Array2 shoc_cldfrac, shoc_ql;
  Array1 pblh;
  Array2 shoc_mix, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec,
    wtke_sec, uw_sec, vw_sec, w3, wqls_sec, isotropy, brunt, shoc_ql2;

  FortranData() = delete;
  FortranData(Int shcol, Int nlev, Int nlevi, Int num_qtracers);
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

// Initialize SHOC with the given number of levels.
void shoc_init(Int nlev, bool use_fortran=false, bool force_reinit=false);

// We will likely want to remove these checks in the future, as we're not tied
// to the exact implementation or arithmetic in SHOC. For now, these checks are
// here to establish that the initial regression-testing code gives results that
// match the python f2py tester, without needing a data file.
Int check_against_python(const FortranData& d);

int test_FortranData();
int test_shoc_init(bool use_fortran);

}  // namespace shoc
}  // namespace scream

#endif
