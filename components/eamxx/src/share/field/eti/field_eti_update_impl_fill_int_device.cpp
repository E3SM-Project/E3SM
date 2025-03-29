#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Device, true, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply,Device, true, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,  Device, true, int, int>(const Field&, const int, const int);

} // namespace scream
