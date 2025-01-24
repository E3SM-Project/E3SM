#ifndef SCREAM_P3_DATA_HPP
#define SCREAM_P3_DATA_HPP

#include "share/eamxx_types.hpp"

#include <memory>
#include <vector>

namespace scream {
namespace p3 {

// Data format we can use to store (and read/write) data for a full P3 run.
struct P3Data {
  typedef std::shared_ptr<P3Data> Ptr;

  using KT     = KokkosTypes<HostDevice>;
  using Scalar = Real;

  using Array1 = typename KT::template view_1d<Scalar>;
  using Array2 = typename KT::template view_2d<Scalar>;
  using Array3 = typename KT::template view_3d<Scalar>;

  bool do_predict_nc;
  bool do_prescribed_CCN;
  bool use_hetfrz_classnuc;
  const Int ncol, nlev;

  // In
  Real dt;
  Int it;
  Array2 qv, th_atm, pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated, inv_qc_relvar, qc, nc, qr, nr,  qi,
    ni, qm, bm, dpres, inv_exner, qv_prev, t_prev;
  Array2 hetfrz_immersion_nucleation_tend, hetfrz_contact_nucleation_tend, hetfrz_deposition_nucleation_tend;
  // Out
  Array1 precip_liq_surf, precip_ice_surf;
  Array2 diag_eff_radius_qc, diag_eff_radius_qi, diag_eff_radius_qr, rho_qi, qv2qi_depos_tend,
         precip_liq_flux, precip_ice_flux, cld_frac_r, cld_frac_l, cld_frac_i;
  Array3 p3_tend_out;
  Array2 liq_ice_exchange,vap_liq_exchange,vap_ice_exchange;
  Array2 P3_qr2qv_evap, P3_qi2qv_sublim, P3_qc2qr_accret, P3_qc2qr_autoconv, P3_qv2qi_vapdep, P3_qc2qi_berg, P3_qc2qr_ice_shed, P3_qc2qi_collect, P3_qr2qi_collect, P3_qc2qi_hetero_freeze, P3_qr2qi_immers_freeze, P3_qi2qr_melt, P3_qr_sed, P3_qc_sed, P3_qi_sed;

  P3Data(Int ncol, Int nlev);
};

// Iterate over a P3Data's arrays. For examples, see Baseline::write, read.
struct P3DataIterator {
  struct RawArray {
    std::string name;
    Int dim;
    Int extent[3];
    P3Data::Scalar* data;
    P3Data::Array1::size_type size;
  };

  explicit P3DataIterator(const P3Data::Ptr& d);

  Int nfield () const { return fields_.size(); }
  const RawArray& getfield(Int i) const;

private:
  P3Data::Ptr d_;
  std::vector<RawArray> fields_;

  void init(const P3Data::Ptr& d);
};

// We will likely want to remove these checks in the future, as we're not tied
// to the exact implementation or arithmetic in P3. For now, these checks are
// here to establish that the initial regression-testing code gives results that
// match the python f2py tester, without needing a data file.
Int check_against_python(const P3Data& d);

int test_P3Data();

}  // namespace p3
}  // namespace scream

#endif
