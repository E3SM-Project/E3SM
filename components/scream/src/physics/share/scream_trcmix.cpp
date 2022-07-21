#include "scream_trcmix.hpp"
#include "physics_constants.hpp"

#include "ekat/ekat_pack_kokkos.hpp"

#include <vector>

using scream::Real;
using scream::Int;

extern "C" {

void trcmix_cf2(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q,
                const Real mwdry, const Real mwco2, const Real mwn2o, const Real mwch4, const Real mwf11, const Real mwf12,
                const Real o2mmr, const Real co2vmr_rad, const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr);

}

namespace scream {
namespace physics {

#define CHECK_TRCMIX_VIEW(vtype, ranks)           \
static_assert(Kokkos::is_view<vtype>::value, "trcmix: " #vtype " is not a KokkosView"); \
static_assert(std::is_same<typename vtype::non_const_value_type, Real>::value, "trcmix: " #vtype " is not View of Reals"); \
static_assert(vtype::rank == ranks, "trcmix: " #vtype " does not have expected rank: " #ranks);

template <typename ClatView, typename PmidView, typename QView>
void trcmix(
  const char* name, Int ncol, Int pcols, Int pver,
  const ClatView& clat, const PmidView& pmid, QView& q,
  const Real co2vmr_rad, const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr)
{
  using C = Constants<Real>;

  //
  // View type checking
  //
  CHECK_TRCMIX_VIEW(ClatView, 1);
  CHECK_TRCMIX_VIEW(PmidView, 2);
  CHECK_TRCMIX_VIEW(QView,    2);

  EKAT_REQUIRE(clat.extent(0) == pcols);
  EKAT_REQUIRE(pmid.extent(0) == pcols && pmid.extent(1) == pver);
  EKAT_REQUIRE(q.extent(0) == pcols && q.extent(1) == pver);

  // Transpose and copy to host
  std::vector<Real> hclat(pcols), hpmid(pcols*pver), hq(pcols*pver);
  ekat::device_to_host(
    {hclat.data(), hpmid.data(), hq.data()},
    {pcols, pcols, pcols},
    {1, pver, pver},
    {clat, pmid, q}, true);

  // Convert to hclat to radians
  for (auto& item : hclat) {
    item *= C::Pi/180.0;
  }

  trcmix_cf2(
    name, ncol, pcols, pver, hclat.data(), hpmid.data(), hq.data(),
    C::MWdry,
    C::get_gas_mol_weight("co2"),
    C::get_gas_mol_weight("n2o"),
    C::get_gas_mol_weight("ch4"),
    C::get_gas_mol_weight("cfc11"),
    C::get_gas_mol_weight("cfc12"),
    C::o2mmr,
    co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr);

  ekat::host_to_device(
    {hq.data()},
    {pcols},
    {pver},
    {q}, true);
}

#undef CHECK_TRCMIX_VIEW

} // namespace physics
} // namespace scream
