#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Host, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Multiply,Host, double>(const Field&, const double, const double);
template void Field::update<CombineMode::Divide,  Host, double>(const Field&, const double, const double);

} // namespace scream
