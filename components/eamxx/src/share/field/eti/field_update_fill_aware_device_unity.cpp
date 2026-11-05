// Unity build file for field_update_fill_aware (device)
// Combines multiple ETI instantiations into a single translation unit
// to reduce header parsing overhead.

#include "eti/field_update_fill_aware_double_double.cpp"
#include "eti/field_update_fill_aware_double_float.cpp"
#include "eti/field_update_fill_aware_double_int.cpp"
#include "eti/field_update_fill_aware_float_float.cpp"
#include "eti/field_update_fill_aware_float_int.cpp"
#include "eti/field_update_fill_aware_int_int.cpp"
