#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::deep_copy_impl<true,int>(const int, const Field&) const;
template void Field::deep_copy_impl<false,int>(const int, const Field&) const;

} // namespace scream
