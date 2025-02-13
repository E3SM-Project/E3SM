#include "share/field/field.hpp"

namespace scream {

template void Field::update<Host, CombineMode::Update,  double>(const Field&, const double, const double);
template void Field::update<Host, CombineMode::Multiply,double>(const Field&, const double, const double);
template void Field::update<Host, CombineMode::Divide,  double>(const Field&, const double, const double);

} // namespace scream
