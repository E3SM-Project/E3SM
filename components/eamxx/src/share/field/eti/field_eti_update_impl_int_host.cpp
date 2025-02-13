#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Host, false, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply,Host, false, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,  Host, false, int>(const Field&, const int, const int);

} // namespace scream
