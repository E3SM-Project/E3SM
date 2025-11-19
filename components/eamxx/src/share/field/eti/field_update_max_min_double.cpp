#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Min, double>(const Field&, const double, const double);

} // namespace scream
