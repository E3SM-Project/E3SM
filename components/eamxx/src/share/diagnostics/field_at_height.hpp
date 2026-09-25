#ifndef EAMXX_FIELD_AT_HEIGHT_HPP
#define EAMXX_FIELD_AT_HEIGHT_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given height (above surface)
 */

class FieldAtHeight : public AbstractDiagnostic
{
public:

  // Constructors
  FieldAtHeight (const ekat::Comm& comm, const ekat::ParameterList& params,
                 const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtHeight"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
protected:
  void initialize_impl ();

  std::string    m_diag_name;
  std::string    m_z_name;
  std::string    m_z_suffix;
  std::string    m_field_name;

  // Names of the (mid/int) level-index diagnostics we could depend on;
  // only one of the two is actually used, decided in initialize_impl,
  // once we know whether the input field lives on LEV or ILEV.
  std::string    m_index_name_mid;
  std::string    m_index_name_int;
  std::string    m_index_name;

  Real           m_z;
};

} //namespace scream

#endif // EAMXX_FIELD_AT_HEIGHT_HPP
