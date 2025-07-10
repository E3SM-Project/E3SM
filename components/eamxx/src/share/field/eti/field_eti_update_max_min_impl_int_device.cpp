#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, Device, false, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Min, Device, false, int, int>(const Field&, const int, const int);

} // namespace scream
