#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<true, CombineMode::Replace,  int, int>(const Field&, const int, const int, const Field*);
template void Field::update_impl<true, CombineMode::Update,   int, int>(const Field&, const int, const int, const Field*);
template void Field::update_impl<true, CombineMode::Multiply, int, int>(const Field&, const int, const int, const Field*);
template void Field::update_impl<true, CombineMode::Divide,   int, int>(const Field&, const int, const int, const Field*);
template void Field::update_impl<true, CombineMode::Max,      int, int>(const Field&, const int, const int, const Field*);
template void Field::update_impl<true, CombineMode::Min,      int, int>(const Field&, const int, const int, const Field*);

} // namespace scream
