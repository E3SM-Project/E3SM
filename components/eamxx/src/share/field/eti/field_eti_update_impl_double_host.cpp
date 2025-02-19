#include "share/field/field.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Update,  Host, false, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Multiply,Host, false, double>(const Field&, const double, const double);
template void Field::update_impl<CombineMode::Divide,  Host, false, double>(const Field&, const double, const double);

} // namespace scream
