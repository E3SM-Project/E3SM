// Unity build file for field_update_masked (device)
// Combines multiple ETI instantiations into a single translation unit
// to reduce header parsing overhead.

#include "eti/field_update_masked_double_double.cpp"
#include "eti/field_update_masked_double_float.cpp"
#include "eti/field_update_masked_double_int.cpp"
#include "eti/field_update_masked_float_float.cpp"
#include "eti/field_update_masked_float_int.cpp"
#include "eti/field_update_masked_int_int.cpp"
