#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, false, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply,Device, false, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,  Device, false, int>(const Field&, const int, const int);

} // namespace scream
