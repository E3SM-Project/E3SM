#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

// Only CombineMode::Replace is needed: deep_copy with allow_narrowing=true is the
// only intentional path for float←double (write-time precision conversion in IO).
template void Field::update_impl<CombineMode::Replace, false, float, double>(const Field&, const float, const float);

} // namespace scream

