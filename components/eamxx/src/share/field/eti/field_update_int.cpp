#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,   int>(const Field&, const int, const int);
template void Field::update<CombineMode::Multiply, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Divide,   int>(const Field&, const int, const int);

} // namespace scream
