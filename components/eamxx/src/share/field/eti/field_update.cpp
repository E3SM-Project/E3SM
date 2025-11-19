#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update>(const Field&, const ScalarWrapper, const ScalarWrapper);

} // namespace scream
