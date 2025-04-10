#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Host,true,int>(const int, const Field&);
template void Field::deep_copy_impl<Host,false,int>(const int, const Field&);

} // namespace scream
