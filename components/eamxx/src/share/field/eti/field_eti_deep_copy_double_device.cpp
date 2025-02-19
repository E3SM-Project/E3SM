#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Device,double>(const double);

} // namespace scream
