/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include <iostream>

#include <mpi.h>

#include "ErrorDefs.hpp"
#include "Hommexx_Session.hpp"

namespace Homme {
namespace Errors {

class HommeException : public std::exception
{
public:
  HommeException (const std::string& msg)
   : m_msg (msg)
  {
    // Nothing to do here
  }
  ~HommeException () = default;

  const char* what () const noexcept override {
    return m_msg.c_str();
  }
private:

  const std::string m_msg;

};

void runtime_check(bool cond, const std::string& message, int code) {
  if (!cond) {
    runtime_abort(message,code);
  }
}

void runtime_abort(const std::string& message, int code) {
  if (Session::m_throw_instead_of_abort) {
    throw HommeException(message);
  } else {
    std::cerr << message << std::endl << "Exiting..." << std::endl;
    finalize_hommexx_session();
    MPI_Abort(MPI_COMM_WORLD, code);
  }
}

} // namespace Homme
} // namespace Errors
