#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, Host, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Min, Host, double>(const Field&, const double, const double);

} // namespace scream
