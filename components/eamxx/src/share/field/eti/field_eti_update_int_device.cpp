#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Device, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Multiply,Device, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Divide,  Device, int>(const Field&, const int, const int);

} // namespace scream
