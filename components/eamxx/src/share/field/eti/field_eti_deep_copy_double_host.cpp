#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<Host,true,double>(const double, const Field&);
template void Field::deep_copy_impl<Host,false,double>(const double, const Field&);

} // namespace scream
