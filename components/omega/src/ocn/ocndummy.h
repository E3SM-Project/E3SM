// OCN dummy header file

#include "YAKL.h"
#include <iostream>

using yakl::Array;
using yakl::memDevice;
using yakl::memHost;
using yakl::styleC;
using yakl::c::Bounds;
using yakl::c::parallel_for;

typedef float real;
typedef Array<real, 1, memDevice, styleC> real1d;
typedef Array<real, 2, memDevice, styleC> real2d;

void die(std::string msg) { yakl::yakl_throw(msg.c_str()); }
