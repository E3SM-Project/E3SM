#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, int, int>(const Field&, const int, const int) const;
template void Field::update_impl<CombineMode::Multiply, true, int, int>(const Field&, const int, const int) const;
template void Field::update_impl<CombineMode::Divide,   true, int, int>(const Field&, const int, const int) const;
template void Field::update_impl<CombineMode::Max,      true, int, int>(const Field&, const int, const int) const;
template void Field::update_impl<CombineMode::Min,      true, int, int>(const Field&, const int, const int) const;

} // namespace scream
