#ifndef EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
#define EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_vertical_interpolation.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class FieldAtPressureLevel : public AtmosphereDiagnostic
{
public:
  using mPack = ekat::Pack<Real,1>;

  using KT = KokkosTypes<DefaultDevice>;
  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  // Constructors
  FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return m_field_name + " @ pressure " + std::to_string(m_pressure_level) + " hPa"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:



  // Keep track of field dimensions
  std::string         m_field_name;
  FieldLayout         m_field_layout;
  ekat::units::Units  m_field_units;
  std::string         m_pres_name;

  view_1d<mPack>      m_p_tgt;
  Field               m_mask_field;
  Real                m_pressure_level;
  int                 m_num_levs;
  int                 m_num_cols;
  Real                m_mask_val;

}; // class FieldAtPressureLevel

} //namespace scream

#endif // EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
