#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Max, false, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Min, false, int, int>(const Field&, const int, const int);

} // namespace scream
