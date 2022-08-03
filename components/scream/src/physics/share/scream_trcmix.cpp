#include "scream_trcmix.hpp"
#include "physics_constants.hpp"

#include "ekat/ekat_pack_kokkos.hpp"

#include <vector>

using scream::Real;
using scream::Int;

extern "C" {

void trcmix_c2f(const char** name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q,
                const Real mwdry, const Real mwco2, const Real mwn2o, const Real mwch4, const Real mwf11, const Real mwf12,
                const Real o2mmr, const Real co2vmr_rad, const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr);

}

namespace scream {
namespace physics {

void trcmix(
  const char* name, Int ncol, Int pcols, Int pver,
  trcmix_view1d<const Real> const& clat,  // latitude for columns in degrees
  trcmix_view2d<const Real> const& pmid,  // model pressures
  trcmix_view2d<Real>            & q,     // constituent mass mixing ratio (output)
  const Real co2vmr_rad, const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr)
{
  using C = Constants<Real>;
  using P = ekat::Pack<Real, 1>;

  //
  // View bounds checking
  //
  EKAT_REQUIRE(clat.extent(0) == static_cast<size_t>(pcols));
  EKAT_REQUIRE(pmid.extent(0) == static_cast<size_t>(pcols) && pmid.extent(1) >= static_cast<size_t>(pver));
  EKAT_REQUIRE(q.extent(0) == static_cast<size_t>(pcols) && q.extent(1) >= static_cast<size_t>(pver));

  // Transpose and copy to host. We reinterpret views as views of Pack<Real, 1> so
  // that we can use the ekat::device_to_host, host_to_device API
  std::vector<Real> hclat(pcols), hpmid(pcols*pver), hq(pcols*pver);

  {
    std::vector<trcmix_view1d<const Real>> stuff = {clat};
    ekat::device_to_host(
      {hclat.data()},
      pcols,
      reinterpret_cast<std::vector<trcmix_view1d<const P>>&>(stuff));
  }
  {
    std::vector<trcmix_view2d<const Real>> stuff = {pmid};
    ekat::device_to_host(
      {hpmid.data()},
      pcols, pver,
      reinterpret_cast<std::vector<trcmix_view2d<const P>>&>(stuff), true);
  }


  // Convert to hclat to radians
  for (auto& item : hclat) {
    item *= C::Pi/180.0;
  }

  trcmix_c2f(
    &name, ncol, pcols, pver, hclat.data(), hpmid.data(), hq.data(),
    C::MWdry,
    C::get_gas_mol_weight("co2"),
    C::get_gas_mol_weight("n2o"),
    C::get_gas_mol_weight("ch4"),
    C::get_gas_mol_weight("cfc11"),
    C::get_gas_mol_weight("cfc12"),
    C::o2mmr,
    co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr);

  {
    std::vector<trcmix_view2d<Real>> stuff = {q};
    ekat::host_to_device(
      {hq.data()},
      pcols, pver,
      reinterpret_cast<std::vector<trcmix_view2d<P>>&>(stuff), true);
  }
}

#undef CHECK_TRCMIX_VIEW

} // namespace physics
} // namespace scream
