#include "scream_trcmix.hpp"
#include "physics_constants.hpp"

#include "ekat/util/ekat_math_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <cmath>

using scream::Real;
using scream::Int;

namespace scream {
namespace physics {

void trcmix(
  const std::string& name,
  trcmix_view1d<const Real> const& clat,  // latitude for columns in degrees
  trcmix_view2d<const Real> const& pmid,  // model pressures
  trcmix_view2d<Real>            & q,     // constituent mass mixing ratio (output)
  const Real co2vmr, const Real n2ovmr, const Real ch4vmr, const Real f11vmr, const Real f12vmr)
{
  using C = Constants<Real>;
  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = KT::ExeSpace;
  using MemberType = KT::MemberType;

  const auto ncols = clat.extent(0);
  const auto nlevs = pmid.extent(1);

  //
  // View bounds checking. pmid and q might be padded if simd vector size > 1
  //
  EKAT_REQUIRE(pmid.extent(0) == ncols);
  EKAT_REQUIRE(q.extent(0) == ncols && q.extent(1) == nlevs);

  const auto mwn2o = C::get_gas_mol_weight("n2o");
  const auto mwch4 = C::get_gas_mol_weight("ch4");
  const auto mwf11 = C::get_gas_mol_weight("cfc11");
  const auto mwf12 = C::get_gas_mol_weight("cfc12");
  const auto mwco2 = C::get_gas_mol_weight("co2");

  const auto rmwn2o = mwn2o/C::MWdry;
  const auto rmwch4 = mwch4/C::MWdry;
  const auto rmwf11 = mwf11/C::MWdry;
  const auto rmwf12 = mwf12/C::MWdry;
  const auto rmwco2 = mwco2/C::MWdry;

  // Constants map: gas_name -> [trop_mmr, scale1_base, scale1_fact, scale2_base, scale2_fact]
  std::map<std::string, std::vector<Real> > const_map = {
    {"ch4",  {rmwch4*ch4vmr, 0.2353, 0.     , 0.2353, 0.0225489}},
    {"n2o",  {rmwn2o*n2ovmr, 0.3478, 0.00116, 0.4000, 0.013333 }},
    {"cfc11",{rmwf11*f11vmr, 0.7273, 0.00606, 1.    , 0.013333 }},
    {"cfc12",{rmwf12*f12vmr, 0.4000, 0.00222, 0.5   , 0.024444 }}
  };

  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols, nlevs);

  if (name == "o2" || name == "co2") {
    const auto val = name == "o2" ? C::o2mmr : rmwco2 * co2vmr;
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const Int& k) {
        q(i, k) = val;
      });
    });
  }
  else {
    const auto it = const_map.find(name);
    EKAT_REQUIRE_MSG(it != const_map.end(), "trcmix: Unhandled gas_name: " + name);

    const auto trop_mmr    = it->second[0];
    const auto scale1_base = it->second[1];
    const auto scale1_fact = it->second[2];
    const auto scale2_base = it->second[3];
    const auto scale2_fact = it->second[4];

    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const Int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const Int& k) {
        // set stratospheric scale height factor for gases. Need to convert
        // clat to radians
        const auto clat_r = clat(i) * C::Pi/180.0;
        const auto dlat = std::abs(57.2958 * clat_r);
        const auto scale = dlat <= 45.0
          ? scale1_base + scale1_fact * dlat
          : scale2_base + scale2_fact * (dlat-45);

        // pressure of tropopause
        const auto ptrop = 250.0e2 - 150.0e2*std::pow(std::cos(clat_r), 2);

        // determine output mass mixing ratios
        if (pmid(i,k) >= ptrop) {
          q(i,k) = trop_mmr;
        }
        else {
          const auto pratio = pmid(i,k)/ptrop;
          q(i,k) = trop_mmr * std::pow(pratio, scale);
        }
      });
    });
  }
}

} // namespace physics
} // namespace scream
