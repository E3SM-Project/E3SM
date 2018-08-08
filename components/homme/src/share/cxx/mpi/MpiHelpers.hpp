/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_MPI_HELPERS_HPP
#define HOMMEXX_MPI_HELPERS_HPP

namespace Homme
{

// The tags to assign to messages. We use tags close to MPI_TAG_UB to avoid clashes
// with tags from Fortran (this is relevant only while refactoring is ongoing).
// NOTE: you are NOT allowed to overlap two exchanges with the same tag

enum ExchangeType : short int {
  MPI_EXCHANGE         = 1000,
  MPI_EXCHANGE_MIN_MAX = 2000
};

// For min/max exchange, we store the two values in a single array, and often need to access it
// to retrieve the max or the min. Instead of hard-coding 0/1, we use a wordy name
enum MinMaxId : int {
  MIN_ID = 0,
  MAX_ID = 1
};

} // namespace Homme

#endif // HOMMEXX_MPI_HELPERS_HPP
