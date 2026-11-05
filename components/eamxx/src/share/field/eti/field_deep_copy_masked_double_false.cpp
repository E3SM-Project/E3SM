#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_masked<false,double> (const double, const Field&) const;

} // namespace scream
