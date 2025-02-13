#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Device,float>(const float);

} // namespace scream
