#ifndef INCLUDE_SCREAM_TYPES
#define INCLUDE_SCREAM_TYPES

namespace scream {

#ifdef SCREAM_DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

typedef int Int;

} // namespace scream

#endif
