#include "share/scream_types.hpp"

namespace scream
{

extern "C"
{

// Set ptrs from F90 into C views

// I believe there's no standard portable way of passing a string from F90 to C. So we cannot do
// void init_field_manager(const char*& name1, Real*& field1, int& rank1, int*& dims1, [other_options_for_field_1, ]
//                         const char*& name2, Real*& field2, int& rank2, int*& dims2, [other_options_for_field_2, ]
//                         ...);
// I _think_ our only option is multiple setup functions (eeeewwww). E.g., :
// void init_co2_field(Real*& co2_field, int& rank, int*& dims [, other_options_for_field_1 ]); // not sure if we need rank/dims or if we can assume to know how each field looks
// void init_qv_field(Real*& qv_field, int& rank, int*& dims [, other_options_for_field_1 ]);   // not sure if we need rank/dims or if we can assume to know how each field looks
//...

// Sync fields that are shared with other components
// void sync_f90_to_cxx ();
// void sync_cxx_to_f90 ();
} // extern "C"

} // namespace scream
