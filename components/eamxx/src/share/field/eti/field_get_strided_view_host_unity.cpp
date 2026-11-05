// Unity build file for field_get_strided_view (host, GPU-only)
// Combines multiple ETI instantiations into a single translation unit
// to reduce header parsing overhead. Only compiled when EAMXX_ENABLE_GPU=ON.

#include "eti/field_get_strided_view_double_host.cpp"
#include "eti/field_get_strided_view_float_host.cpp"
#include "eti/field_get_strided_view_int_host.cpp"
