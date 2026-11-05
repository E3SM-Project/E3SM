#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update, int, float>(const Field&, const int, const int, const int) const;

} // namespace scream
