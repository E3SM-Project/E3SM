#ifndef EAMXX_WATER_ISOTOPES_HPP
#define EAMXX_WATER_ISOTOPES_HPP

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_process_interface.hpp"
#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"
#include "eamxx_water_isotopes_fractionation.hpp"
#include "eamxx_water_isotopes_parameters.hpp"
#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{
namespace wiso
{
/*
 * Upward (positive) ocean evaporation flux of one water isotopologue at the
 * bottom of a column [kg/m2/s], with equilibrium fractionation only:
 *
 *   f <= 0 (not a valid ocean point): F_i = 0
 *   E_ocn > 0 (evaporation):          F_i = bulk_ocean_flux(f, rho, u, r, s*R_ocn/alpha, q_i)
 *   E_ocn <= 0 (dew):                 F_i = wtrc_ratio(q_i, q) * E_ocn
 *
 * alpha = R_liquid/R_vapor at the SST, so s*R_ocn/alpha is the isotopologue's
 * saturation specific humidity over the ocean (the delq term of iCAM's
 * wiso_flxoce). For H216O alpha == 1 exactly, and the default normalized
 * ratios give R_ocn == 1, so q_i == q yields F_i == E_ocn bit-for-bit.
 *
 * Inputs, all at the bottom level: f = ocnfrac (the caller zeroes it where the
 * ocean state is invalid), rho = air density, u = ocean friction velocity,
 * r = moisture transfer coefficient (re), s = ocean saturation specific
 * humidity (ssq), q = qv, q_i = this tracer, e_ocn = bulk_ocean_flux(f, rho,
 * u, r, s, q), t_s = SST [K], r_ocn = the species' ocean ratio.
 *
 * alpha is evaluated only on evaporating lanes. Elsewhere the SST may be
 * missing or garbage, and passing it to alpha's temperature guard would abort
 * the run. ScalarT is ekat::Pack<Real,N>; WaterIsotopes calls it with N=1.
 */
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT isotope_ocean_flux (
    const ScalarT& f, const ScalarT& rho, const ScalarT& u, const ScalarT& r,
    const ScalarT& s, const ScalarT& q, const ScalarT& q_i, const ScalarT& e_ocn,
    const ScalarT& t_s,
    const typename ekat::ScalarTraits<ScalarT>::scalar_type r_ocn,
    const WaterIsotopologues species,
    const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const auto ocean = f > RealT(0);
  const auto evap  = ocean && (e_ocn > RealT(0));
  const auto dew   = ocean && !(e_ocn > RealT(0));

  ScalarT flux(0);
  if (evap.any()) {
    // Lanes outside evap get an in-range placeholder temperature inside
    // alpha_liquid_vapor; their result is discarded by the masked set.
    const ScalarT alpha = WaterIsotopeFractionation::alpha_liquid_vapor(
        t_s, species, CondensedOverVapor, iso_params, evap);
    flux.set(evap, water_tracers::bulk_ocean_flux(f, rho, u, r, s*r_ocn/alpha, q_i));
  }
  if (dew.any()) {
    flux.set(dew, water_tracers::wtrc_ratio(q_i, q)*e_ocn);
  }
  return flux;
}

} // namespace wiso

/*
 * The water_isotopes process: WaterTracers plus water isotopologue tracers.
 *
 * In addition to everything WaterTracers handles, it accepts `type: isotope`
 * entries in the shared tracers_file, for species H216O, HDO, H218O, H217O or
 * HTO. Each isotope's ocean ratio, R_ocn = std_iso_ratio * mean_ocean_enrichment
 * for the selected runtime formulations, is resolved on the host at
 * construction. It is also the isotope's default `init: ratio` value.
 *
 * The only isotope physics so far is the diagnostic ocean evaporation flux
 * <name>_surf_flux, with equilibrium liquid/vapor fractionation at the SST
 * (wiso::isotope_ocean_flux). Nothing is applied to qv or to any tracer.
 */

class WaterIsotopes : public WaterTracers
{
public:
  // Constructors
  WaterIsotopes (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Override the name to distinguish from base WaterTracers
  std::string name () const override { return "water_isotopes"; }

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  // Writes <name>_surf_flux for every isotope tracer. WaterTracers::run_impl
  // calls this after filling the bottom-level bulk state (m_sfc_*) and before
  // the property-tracer fluxes, so an isotope can be a property carrier.
  void compute_isotope_fluxes () override;

private:

  struct IsotopeTracer {
    std::string              name;
    wiso::WaterIsotopologues species;
    Real                     r_ocn;
  };

  // Runtime configuration options for isotope formulations
  wiso::WaterIsotopeRuntimeOptions m_runtime_options;

  // Built once on the host from m_runtime_options and captured by value in
  // kernels. Device code then reads only the resolved coefficient rows, never
  // the constexpr formulation tables (PR #8697 B4).
  wiso::WaterIsotopeParameters<Real> m_iso_params;

  std::vector<IsotopeTracer> m_isotopes;

}; // class WaterIsotopes

} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_HPP
