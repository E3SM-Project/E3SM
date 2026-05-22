#ifndef EAMXX_WATER_PATH_HPP
#define EAMXX_WATER_PATH_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class WaterPath : public AbstractDiagnostic
{
public:
  // Constructors
  WaterPath (const ekat::Comm& comm, const ekat::ParameterList& params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "WaterPath"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl () override;
protected:

  std::string m_qname;
  std::string m_kind;
};

} //namespace scream

#endif // EAMXX_WATER_PATH_HPP
