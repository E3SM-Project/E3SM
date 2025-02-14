#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update,  Host, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Multiply,Host, int>(const Field&, const int, const int);
template void Field::update<CombineMode::Divide,  Host, int>(const Field&, const int, const int);

} // namespace scream
