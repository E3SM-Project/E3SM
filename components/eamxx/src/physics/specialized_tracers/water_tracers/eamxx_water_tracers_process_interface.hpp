#ifndef EAMXX_WATER_TRACERS_HPP
#define EAMXX_WATER_TRACERS_HPP

#include "eamxx_water_tracers_registry.hpp"
#include "../water_isotopes/eamxx_water_isotopes.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{
/*
 * Water tracers: extra water species that follow qv through the model.
 *
 * The tracers are defined in a YAML file (the tracers_file parameter; see
 * eamxx_water_tracers_registry.hpp). Each one is registered as an advected
 * tracer <name> [kg/kg], plus a diagnostic <name>_surf_flux [kg/m2/s], the
 * upward ocean evaporation flux of that tracer. The flux is computed every
 * step from the ocean bulk state imported from the coupler (ocn_ssq, ocn_re,
 * ocn_ustar) and the bottom-level atmosphere, and is NOT applied to anything:
 * this process never changes qv or a tracer after initialization.
 *
 * Tracer kinds: generic (bulk_proportional or bulk_formula), property
 * tracers carried by a mass tracer, and isotopes. Isotope entries are only
 * accepted when EAMXX_HAS_WATER_ISOTOPES is defined; isotope-specific state
 * and math live in the composed wiso::WaterIsotopes member, m_isotopes.
 */

class WaterTracers : public AtmosphereProcess
{
public:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  std::string name () const override { return "water_tracers"; }

  void create_requests () override;

protected:

#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif

  void initialize_impl (const RunType run_type) override;
  void run_impl        (const double dt) override;
  void finalize_impl   () override {}

  // Steps of run_impl
  void compute_bulk_state ();
  void compute_generic_fluxes ();
  void compute_isotope_fluxes ();
  void compute_property_fluxes ();

  // Sets q = r_init*q_ref for every init: ratio tracer on the given grid
  void apply_ratio_initial_conditions (const std::string& grid_name);

protected:

  int m_num_cols;
  int m_num_levs;

  std::shared_ptr<const AbstractGrid> m_grid;

  // One entry per tracer in tracers_file, sorted by name. r_ocn/r_init for
  // isotope specs are resolved by m_isotopes.resolve() in the constructor,
  // before initialize/run time reads them.
  std::vector<water_tracers::TracerSpec> m_specs;

  // Isotope-specific state (species table, runtime formulation options,
  // resolved isotope tracers). Always constructed; empty/inert unless
  // EAMXX_HAS_WATER_ISOTOPES is defined and the tracers file defines
  // isotope entries. See compute_isotope_fluxes().
  wiso::WaterIsotopes m_isotopes;

  // Per-column bottom-level bulk state, filled at the start of run_impl:
  //  - m_sfc_f: ocnfrac, set to 0 where the ocean state is not valid
  //  - m_sfc_ustar, m_sfc_re, m_sfc_ssq: ocn_ustar, ocn_re, ocn_ssq, set to
  //    0 where the ocean state is not valid (so they never hold fill values)
  //  - m_sfc_rho: bottom-level air density [kg/m3]
  //  - m_sfc_qv: bottom-level qv [kg/kg]
  //  - m_sfc_evap: E_ocn = bulk_ocean_flux(f,rho,u,r,s,q) [kg/m2/s], upward
  //    positive, exactly 0 where the ocean state is not valid
  //  - m_sfc_sst: sst [K], as imported
  KT::view_1d<Real> m_sfc_f, m_sfc_rho, m_sfc_ustar, m_sfc_re, m_sfc_ssq,
                    m_sfc_qv, m_sfc_evap, m_sfc_sst;

private:

  // The grid on which the driver sets initial conditions, if it differs from
  // m_grid (physics_gll in a pgN run), else null. See create_requests.
  std::shared_ptr<const AbstractGrid> m_ic_grid;

}; // class WaterTracers

} // namespace scream

#endif // EAMXX_WATER_TRACERS_HPP
