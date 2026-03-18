#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,true,int>(const int, const Field*);
template void Field::deep_copy_impl<true,false,int>(const int, const Field*);
template void Field::deep_copy_impl<false,false,int>(const int, const Field*);

} // namespace scream
