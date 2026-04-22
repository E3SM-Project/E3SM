#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  double, int>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Update,   double, int>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Multiply, double, int>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Divide,   double, int>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Max,      double, int>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Min,      double, int>(const Field&, const double, const double, const double, const Field&);

} // namespace scream
