#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<true,float>(const float, const Field&);
template void Field::deep_copy_impl<false,float>(const float, const Field&);

} // namespace scream
