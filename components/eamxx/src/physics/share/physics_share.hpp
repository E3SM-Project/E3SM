#ifndef SCREAM_PHYSICS_SHARE_HPP
#define SCREAM_PHYSICS_SHARE_HPP

#include "share/scream_types.hpp"

namespace scream {

extern "C" {

Real cxx_pow(Real base, Real exp);
Real cxx_sqrt(Real base);
Real cxx_cbrt(Real base);
Real cxx_gamma(Real input);
Real cxx_log(Real input);
Real cxx_log10(Real input);
Real cxx_exp(Real input);
Real cxx_expm1(Real input);
Real cxx_tanh(Real input);
Real cxx_erf(Real input);

}

} // namespace scream

#endif // SCREAM_PHYSICS_SHARE_HPP
