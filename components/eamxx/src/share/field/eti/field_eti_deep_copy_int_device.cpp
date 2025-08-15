#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Device,true,int>(const int, const Field&);
template void Field::deep_copy_impl<Device,false,int>(const int, const Field&);

} // namespace scream
