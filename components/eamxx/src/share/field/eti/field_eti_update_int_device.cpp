#include "share/field/field.hpp"

namespace scream {

template void Field::update<Device, CombineMode::Update,  int>(const Field&, const int, const int);
template void Field::update<Device, CombineMode::Multiply,int>(const Field&, const int, const int);
template void Field::update<Device, CombineMode::Divide,  int>(const Field&, const int, const int);

} // namespace scream
