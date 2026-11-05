// Unity build file for field_update_impl (device)
// Combines multiple ETI instantiations into a single translation unit
// to reduce header parsing overhead.

#include "eti/field_update_impl_double_double.cpp"
#include "eti/field_update_impl_double_float.cpp"
#include "eti/field_update_impl_double_int.cpp"
#include "eti/field_update_impl_float_float.cpp"
#include "eti/field_update_impl_float_int.cpp"
#include "eti/field_update_impl_int_int.cpp"
