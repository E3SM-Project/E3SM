#include "scream_trcmix.hpp"
#include "physics_constants.hpp"

using scream::Real;
using scream::Int;

extern "C" {

void trcmix_cf2(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q,
                const Real mwdry, const Real mwco2, const Real mwn2o, const Real mwch4, const Real mwf11, const Real mwf12,
                const Real o2mmr, const Real co2vmr_rad, const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr);

}

namespace scream {
namespace physics {

void trcmix(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q)
{
  using C = Constants<Real>;

  Real o2mmr, co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr;

  //d.transpose<ekat::TransposeDirection::c2f>();
  trcmix_cf2(
    name, ncol, pcols, pver, clat, pmid, q,
    C::MWdry,
    C::get_gas_mol_weight("co2"),
    C::get_gas_mol_weight("n2o"),
    C::get_gas_mol_weight("ch4"),
    C::get_gas_mol_weight("cfc11"),
    C::get_gas_mol_weight("cfc12"),
    o2mmr, co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr);
  //d.transpose<ekat::TransposeDirection::f2c>();
}

} // namespace physics
} // namespace scream
