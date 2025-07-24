#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Max, Host, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Min, Host, int>(const Field&, const int, const int);

} // namespace scream
