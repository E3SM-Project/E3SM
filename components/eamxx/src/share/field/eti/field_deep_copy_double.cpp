#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,double>(const double, const Field&) const;
template void Field::deep_copy_impl<false,double>(const double, const Field&) const;

} // namespace scream
