#ifndef SCREAM_MASS_ENERGY_COLUMN_CONSERVATION_CHECK_HPP
#define SCREAM_MASS_ENERGY_COLUMN_CONSERVATION_CHECK_HPP

#include "share/property_checks/property_check.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"

namespace scream {

// This property check ensures that energy has been conserved.
// It is a column-local check meant only for column independant processes.
class MassAndEnergyColumnConservationCheck: public PropertyCheck {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpaceUtils = ekat::ExeSpaceUtils<KT::ExeSpace>;

  template<typename ScalarT>
  using view_1d = typename KT::template view_1d<ScalarT>;
  template<typename ScalarT>
  using view_2d = typename KT::template view_2d<ScalarT>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

public:

  // Constructor
  MassAndEnergyColumnConservationCheck (const std::shared_ptr<const AbstractGrid>& grid,
                                        const Real                                 mass_error_tolerance,
                                        const Real                                 energy_error_tolerance,
                                        const std::shared_ptr<const Field>&        pseudo_density_ptr,
                                        const std::shared_ptr<const Field>&        ps_ptr,
                                        const std::shared_ptr<const Field>&        phis_ptr,
                                        const std::shared_ptr<const Field>&        horiz_winds_ptr,
                                        const std::shared_ptr<const Field>&        T_mid_ptr,
                                        const std::shared_ptr<const Field>&        qv_ptr,
                                        const std::shared_ptr<const Field>&        qc_ptr,
                                        const std::shared_ptr<const Field>&        qr_ptr,
                                        const std::shared_ptr<const Field>&        qi_ptr,
                                        const std::shared_ptr<const Field>&        vapor_flux_ptr,
                                        const std::shared_ptr<const Field>&        water_flux_ptr,
                                        const std::shared_ptr<const Field>&        ice_flux_ptr,
                                        const std::shared_ptr<const Field>&        heat_flux_ptr);

  // The name of the property check
  std::string name () const override { return "Mass and energy column conservation check"; }

  PropertyType type () const override { return PropertyType::ColumnWise; }

  // Computes mass and energy and tests against a tolerance.
  ResultAndMsg check () const override;

  std::shared_ptr<const AbstractGrid> get_grid () const { return m_grid; }

  // Set the timestep for the process running the check. This
  // should be the timestep of the current subcycle, equiv to
  // dt = model_dt/num_subcycles.
  void set_dt (const int dt) { m_dt = dt; }

  // Compute total mass and store into m_current_mass.
  // Each process that calls this checker needs to
  // call this function before updating any fields
  // in m_fields.
  void compute_current_mass ();

  // Compute total energy and store into m_current_energy.
  // Each process that calls this checker needs to
  // call this function before updating any fields
  // in m_fields.
  void compute_current_energy ();

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifndef KOKKOS_ENABLE_CUDA
  protected:
#endif

  KOKKOS_INLINE_FUNCTION
  static Real compute_total_mass_on_column (const KT::MemberType&       team,
                                            const int                   nlevs,
                                            const uview_1d<const Real>& pseudo_density,
                                            const uview_1d<const Real>& qv,
                                            const uview_1d<const Real>& qc,
                                            const uview_1d<const Real>& qi,
                                            const uview_1d<const Real>& qr);

  KOKKOS_INLINE_FUNCTION
  static Real compute_mass_boundary_flux_on_column (const Real vapor_flux,
                                                    const Real water_flux);

  KOKKOS_INLINE_FUNCTION
  static Real compute_total_energy_on_column (const KT::MemberType&       team,
                                              const int                   nlevs,
                                              const uview_1d<const Real>& pseudo_density,
                                              const uview_1d<const Real>& T_mid,
                                              const uview_2d<const Real>& horiz_winds,
                                              const uview_1d<const Real>& qv,
                                              const uview_1d<const Real>& qc,
                                              const uview_1d<const Real>& qr,
                                              const Real                  ps,
                                              const Real                  phis);

  KOKKOS_INLINE_FUNCTION
  static Real compute_energy_boundary_flux_on_column (const Real vapor_flux,
                                                      const Real water_flux,
                                                      const Real ice_flux,
                                                      const Real heat_flux);

protected:

  std::shared_ptr<const AbstractGrid>                 m_grid;
  std::map<std::string, std::shared_ptr<const Field>> m_fields;

  int m_num_cols;
  int m_num_levs;
  Real m_dt;
  Real m_mass_tol;
  Real m_energy_tol;

  // Current value for total energy. These values
  // should be updated before a process is run.
  view_1d<Real> m_current_energy;
  view_1d<Real> m_current_mass;
}; // class EnergyConservationCheck

} // namespace scream

#endif //SCREAM_MASS_ENERGY_COLUMN_CONSERVATION_CHECK_HPP
