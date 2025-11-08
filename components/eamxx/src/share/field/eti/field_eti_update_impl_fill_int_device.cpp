#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,   true, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply, true, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,   true, int, int>(const Field&, const int, const int);

} // namespace scream
