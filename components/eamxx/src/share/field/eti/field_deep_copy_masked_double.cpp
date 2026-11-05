#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

// NOTE: deep_copy_masked is actually just deep_copy with a mask parameter
// It uses the same deep_copy_impl template, so we instantiate it here for completeness
template void Field::deep_copy_impl<double> (const double) const;

} // namespace scream
