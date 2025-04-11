#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Device,true,float>(const float, const Field&);
template void Field::deep_copy_impl<Device,false,float>(const float, const Field&);

} // namespace scream
