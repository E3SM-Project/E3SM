/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_INTERNAL_DIAGNOSTICS_HPP
#define HOMMEXX_INTERNAL_DIAGNOSTICS_HPP

#include <string>

namespace Homme {

// Print a hash of ElementsState and Tracer values on rank 0.
void print_global_state_hash(const std::string& label);

} // Homme

#endif // INTERNAL_DIAGNOSTICS_HPP
