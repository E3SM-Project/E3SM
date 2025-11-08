#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   false, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply, false, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,   false, int, int>(const Field&, const int, const int);

} // namespace scream
