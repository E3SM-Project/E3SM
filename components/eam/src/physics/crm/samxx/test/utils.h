
#include <string>
#include <iostream>
#include <cmath>
#include "abcoefs.h"



// Computes a relative difference between v1 and v2, where v1 is the normalizer.
// Compares the relative difference against the tolerance
//    If rel diff <= tol, returns 0
//    Otherwise returns -1
// If v1 is zero, then the absolute difference is compared against tolerance to avoid
// division by zero
inline template <class T> int compareScalar( std::string const &label , T v1 , T v2 , T tol ) {
  double diff = std::abs( (double) v1 - (double) v2);
  if (v1 != 0) { diff /= (double) v1; }
  bool pass = diff <= tol;
  std::cout << std::scientific;
  std::cout << label << ": " << (pass ? "PASS" : "FAIL") << std::endl;
  if (!pass) {
    std::cout << "Rel Diff:  " << diff << std::endl;
    std::cout << "v1:        " << v1   << std::endl;
    std::cout << "v2:        " << v2   << std::endl;
    std::cout << "Tolerance: " << tol  << std::endl;
  }
  std::cout << std::endl;
  return pass ? 0 : -1;
}



inline void compareArray( std::string const &label , real1d &v1 , real1d &v2 ) {}

inline void compareArray( std::string const &label , real2d &v1 , real2d &v2 ) {}

inline void compareArray( std::string const &label , real3d &v1 , real3d &v2 ) {}

inline void compareArray( std::string const &label , real4d &v1 , real4d &v2 ) {}

inline void compareArray( std::string const &label , real5d &v1 , real5d &v2 ) {}

inline void compareArray( std::string const &label , real6d &v1 , real6d &v2 ) {}



inline void initDomain( Domain &dom ) {}



template <class T> inline void loadScalar( std::string const &label , T &data ) {}



inline void loadArray( std::string const &label , real1d &data ) {}

inline void loadArray( std::string const &label , real2d &data ) {}

inline void loadArray( std::string const &label , real3d &data ) {}

inline void loadArray( std::string const &label , real4d &data ) {}

inline void loadArray( std::string const &label , real5d &data ) {}

inline void loadArray( std::string const &label , real6d &data ) {}





