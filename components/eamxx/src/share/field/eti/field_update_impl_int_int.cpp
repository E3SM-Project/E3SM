#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace,  int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Update,   int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Multiply, int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Divide,   int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Max,      int, int>(const Field&, const int, const int);
template void Field::update_impl<CombineMode::Min,      int, int>(const Field&, const int, const int);

} // namespace scream
