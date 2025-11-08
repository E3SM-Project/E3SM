#include "share/field/field.hpp"

namespace scream {

template void Field::update<CombineMode::Update>(const Field&, const ScalarWrapper, const ScalarWrapper);
template void Field::update<CombineMode::Multiply>(const Field&, const ScalarWrapper, const ScalarWrapper);
template void Field::update<CombineMode::Divide>(const Field&, const ScalarWrapper, const ScalarWrapper);
template void Field::update<CombineMode::Max>(const Field&, const ScalarWrapper, const ScalarWrapper);
template void Field::update<CombineMode::Min>(const Field&, const ScalarWrapper, const ScalarWrapper);

} // namespace scream
