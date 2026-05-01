#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  int, int>(const Field&, const int, const int, const int, const Field&);
template void Field::update_masked<CombineMode::Update,   int, int>(const Field&, const int, const int, const int, const Field&);
template void Field::update_masked<CombineMode::Multiply, int, int>(const Field&, const int, const int, const int, const Field&);
template void Field::update_masked<CombineMode::Divide,   int, int>(const Field&, const int, const int, const int, const Field&);
template void Field::update_masked<CombineMode::Max,      int, int>(const Field&, const int, const int, const int, const Field&);
template void Field::update_masked<CombineMode::Min,      int, int>(const Field&, const int, const int, const int, const Field&);

} // namespace scream
