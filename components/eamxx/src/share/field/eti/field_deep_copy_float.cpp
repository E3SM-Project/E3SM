#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,float>(const float, const Field&);
template void Field::deep_copy_impl<false,float>(const float, const Field&);

} // namespace scream
