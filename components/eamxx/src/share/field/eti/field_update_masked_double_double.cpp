#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace,  double, double>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Update,   double, double>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Multiply, double, double>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Divide,   double, double>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Max,      double, double>(const Field&, const double, const double, const double, const Field&);
template void Field::update_masked<CombineMode::Min,      double, double>(const Field&, const double, const double, const double, const Field&);

} // namespace scream
