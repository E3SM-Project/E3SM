#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Host,true,float>(const float, const Field&);
template void Field::deep_copy_impl<Host,false,float>(const float, const Field&);

} // namespace scream
