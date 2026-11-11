#ifndef EAMXX_NUMBER_PATH_HPP
#define EAMXX_NUMBER_PATH_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will produce the "number" paths of liquid, ice and rain,
 * in the fields LiqNumberPath, IceNumberPath and RainNumberPath.
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
  Field m_liq;
  Field m_ice;
  Field m_rain;
};

}  // namespace scream

#endif  // EAMXX_NUMBER_PATH_HPP
