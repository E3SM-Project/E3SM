#ifndef EAMXX_VERTICALLY_REMAPPED_FIELD_HPP
#define EAMXX_VERTICALLY_REMAPPED_FIELD_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_vertical_interpolation.hpp"

namespace scream
{

class VerticallyRemappedField : public AtmosphereDiagnostic
{
public:
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT = KokkosTypes<DefaultDevice>;
  using view_1d_const = typename KT::template view_1d<const Pack>;
  using view_1d = typename KT::template view_1d<Pack>;

  // Constructors
  VerticallyRemappedField (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return m_field_name; }

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

  int                 m_num_levs;
  int                 m_num_cols;
  int                 m_tgt_num_levs;
  const view_1d_const m_tgt_pres_levs;

}; // class VerticallyRemappedField

} //namespace scream

#endif // EAMXX_VERTICALLY_REMAPPED_FIELD_HPP
