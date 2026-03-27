#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  int, int>(const Field&, const int, const int);
template void Field::update_fill_aware<CombineMode::Update,   int, int>(const Field&, const int, const int);
template void Field::update_fill_aware<CombineMode::Multiply, int, int>(const Field&, const int, const int);
template void Field::update_fill_aware<CombineMode::Divide,   int, int>(const Field&, const int, const int);
template void Field::update_fill_aware<CombineMode::Max,      int, int>(const Field&, const int, const int);
template void Field::update_fill_aware<CombineMode::Min,      int, int>(const Field&, const int, const int);

} // namespace scream
