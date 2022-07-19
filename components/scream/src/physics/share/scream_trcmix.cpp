#include "scream_trcmix.hpp"
#include "physics_constants.hpp"

using scream::Real;
using scream::Int;

extern "C" {

void trcmix_cf2(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q,
                const Real mwdry, const Real mwco2, const Real mwn2o, const Real mwch4, const Real mwf11, const Real mwf12);

}

namespace scream {
namespace physics {

void trcmix(const char* name, Int ncol, Int pcols, Int pver, Real* clat, Real* pmid, Real* q)
{
  using C = Constants<Real>;

  //d.transpose<ekat::TransposeDirection::c2f>();
  trcmix_cf2(
    name, ncol, pcols, pver, clat, pmid, q,
    C::MWdry,
    C::get_gas_mol_weight("co2"),
    C::get_gas_mol_weight("n2o"),
    C::get_gas_mol_weight("ch4"),
    C::get_gas_mol_weight("cfc11"),
    C::get_gas_mol_weight("cfc12"));
  //d.transpose<ekat::TransposeDirection::f2c>();
}

} // namespace physics
} // namespace scream
