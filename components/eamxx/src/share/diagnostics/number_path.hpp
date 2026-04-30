#ifndef EAMXX_NUMBER_PATH_HPP
#define EAMXX_NUMBER_PATH_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will produce the "number" path.
 */

class NumberPath : public AbstractDiagnostic {
 public:
  // Constructors
  NumberPath(const ekat::Comm &comm,
             const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "NumberPath"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl() override;

 protected:
  std::string m_qname;
  std::string m_nname;
  std::string m_kind;
};

}  // namespace scream

#endif  // EAMXX_NUMBER_PATH_HPP
