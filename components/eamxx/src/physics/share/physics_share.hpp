#ifndef SCREAM_PHYSICS_SHARE_HPP
#define SCREAM_PHYSICS_SHARE_HPP

#include "share/eamxx_types.hpp"

namespace scream {

extern "C" {

Real scream_pow(Real base, Real exp);
Real scream_sqrt(Real base);
Real scream_cbrt(Real base);
Real scream_gamma(Real input);
Real scream_log(Real input);
Real scream_log10(Real input);
Real scream_exp(Real input);
Real scream_expm1(Real input);
Real scream_tanh(Real input);
Real scream_erf(Real input);

}

} // namespace scream

#endif // SCREAM_PHYSICS_SHARE_HPP
