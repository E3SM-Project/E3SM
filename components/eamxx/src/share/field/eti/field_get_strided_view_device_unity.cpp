// Unity build file for field_get_strided_view (device)
// Combines multiple ETI instantiations into a single translation unit
// to reduce header parsing overhead.

#include "eti/field_get_strided_view_double_device.cpp"
#include "eti/field_get_strided_view_float_device.cpp"
#include "eti/field_get_strided_view_int_device.cpp"
