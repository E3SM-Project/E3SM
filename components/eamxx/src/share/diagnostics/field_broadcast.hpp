#ifndef EAMXX_FIELD_BROADCAST_DIAG_HPP
#define EAMXX_FIELD_BROADCAST_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic broadcasts a field to the layout of another field.
 *
 * The layout of the field to broadcast must be an ordered subset of the target layout
 * (see Field::broadcast). The target field is only used to get the layout:
 * its data is never read.
 *
 * Required parameters:
 *  - field_name:  the name of the field to broadcast
 *  - target_name: the name of the field whose layout is the broadcast target
 *
 * In an expression, this is X.broadcast_like(Y).
 */

class FieldBroadcast : public AbstractDiagnostic {
 public:
  FieldBroadcast (const ekat::Comm &comm, const ekat::ParameterList &params,
                  const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "FieldBroadcast"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_name;
  std::string m_tgt_name;

  // A (read-only) view of the input, broadcasted to the target layout
  Field m_broadcasted;
};

}  // namespace scream

#endif  // EAMXX_FIELD_BROADCAST_DIAG_HPP
