#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Host, true, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply,Host, true, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,  Host, true, int>(const Field&, const int, const int);

} // namespace scream
