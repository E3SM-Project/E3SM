#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

// No fill-aware float←double instantiations needed:
// float←double conversion is only used via deep_copy(src, allow_narrowing=true),
// which uses CombineMode::Replace (non-fill-aware).

namespace scream {
} // namespace scream

