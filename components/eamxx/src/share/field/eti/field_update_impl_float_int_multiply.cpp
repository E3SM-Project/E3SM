#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Multiply, float, int>(const Field&, const float, const float, const float) const;

} // namespace scream
