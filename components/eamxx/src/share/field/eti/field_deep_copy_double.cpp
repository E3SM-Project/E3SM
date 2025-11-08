#include "share/field/field.hpp"

namespace scream {

template void Field::deep_copy_impl<true,double>(const double, const Field&);
template void Field::deep_copy_impl<false,double>(const double, const Field&);

} // namespace scream
