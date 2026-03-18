#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,true,float>(const float, const Field*);
template void Field::deep_copy_impl<true,false,float>(const float, const Field*);
template void Field::deep_copy_impl<false,false,float>(const float, const Field*);

} // namespace scream
