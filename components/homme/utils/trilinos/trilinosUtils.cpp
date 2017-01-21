#include "trilinosUtils.hpp"

// ==============================================================================
// check function return value
//
// opt == 0 means check if returned NULL pointer
// opt == 1 means check if flag /= 0
// ==============================================================================

void check_flag(void *flagvalue, const std::string message, int opt, int myid)
{

  if (opt == 0 && flagvalue == NULL) {

    // Error, function returned NULL pointer, halt run
    std::cout << "\n" << myid << " ERROR: " 
              << message << ", returned NULL pointer\n\n";
    std::cerr << "\n" << myid << " ERROR: " 
              << message << ", returned NULL pointer\n\n";
    abortrun();
 
  } else if (opt == 1) {

    int *errflag;
    errflag = (int *) flagvalue;

    // Warning, positive return value, continue run    
    if (*errflag > 0) {
      std::cerr << "\n" << myid << " WARNING: " 
                << message << ", returned with flag = " << *errflag << "\n\n";
    }

    // Error, negative return value, halt run  
    if (*errflag < 0) {
      std::cout << "\n" << myid << " ERROR: " 
                << message << ", returned with flag = " << *errflag << "\n\n";
      std::cerr << "\n" << myid << " ERROR: " 
                << message << ", returned with flag = " << *errflag << "\n\n";
      abortrun();
    }
  }
}
