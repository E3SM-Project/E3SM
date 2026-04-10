#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  double, float>(const Field&, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Update,   double, float>(const Field&, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Multiply, double, float>(const Field&, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Divide,   double, float>(const Field&, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Max,      double, float>(const Field&, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Min,      double, float>(const Field&, const double, const double, const Field&);

} // namespace scream
