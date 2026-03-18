#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,true,double>(const double, const Field*);
template void Field::deep_copy_impl<true,false,double>(const double, const Field*);
template void Field::deep_copy_impl<false,false,double>(const double, const Field*);

} // namespace scream
