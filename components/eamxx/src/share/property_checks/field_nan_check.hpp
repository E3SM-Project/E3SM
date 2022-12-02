#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/property_checks/property_check.hpp"
#include "share/grid/abstract_grid.hpp"

namespace scream
{

// Inspect whether a field contains any NaN values.
// No repair allowed. If we find NaN's, we should crash.
class FieldNaNCheck: public PropertyCheck {
public:
  FieldNaNCheck (const Field& f,
                 const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the field check
  std::string name () const override {
    return "NaN check for field " + fields().front().name();
  }

  PropertyType type () const override { return PropertyType::PointWise; }

  ResultAndMsg check() const override;

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifndef EAMXX_ENABLE_GPU
protected:
#endif
  template<typename ST>
  ResultAndMsg check_impl() const;

private:

  std::shared_ptr<const AbstractGrid> m_grid;
};

} // namespace scream

#endif // SCREAM_FIELD_NAN_CHECK_HPP
