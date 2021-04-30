#ifndef RRTMGP_HEATING_RATE_HPP
#define RRTMGP_HEATING_RATE_HPP
#include "physics/share/physics_constants.hpp"
#include "YAKL/YAKL.h"
#include "cpp/const.h"
namespace scream {
    namespace rrtmgp {
        // Provide a routine to compute heating due to radiative fluxes. This is
        // computed as net flux into a layer, converted to a heating rate. It is
        // the responsibility of the user to ensure fields are passed with the
        // proper units. I.e., pressure at level interfaces should be in Pa,
        // fluxes in W m-2, Cpair in J kg-1 K-1, gravit in m s-2. This will give
        // heating in units of K s-1.
        // TODO: we should probably update this to use the pseudo-density dp instead
        // of approximating dp by differencing the level interface pressures.
        // We are leaving this for the time being for consistency with SCREAMv0,
        // from which this code was directly ported.
        template <class T, int myMem, int myStyle> void compute_heating_rate (
                yakl::Array<T,2,myMem,myStyle> const &flux_up, yakl::Array<T,2,myMem,myStyle> const &flux_dn, 
                yakl::Array<T,2,myMem,myStyle> const &dp     , yakl::Array<T,2,myMem,myStyle> &heating_rate
            ) {
            using physconst = scream::physics::Constants<Real>;
            auto ncol = flux_up.dimension[0];
            auto nlay = flux_up.dimension[1]-1;
            parallel_for(Bounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
                heating_rate(icol,ilay) = (
                    flux_up(icol,ilay+1) - flux_up(icol,ilay) - 
                    flux_dn(icol,ilay+1) + flux_dn(icol,ilay)
                ) * physconst::gravit / (physconst::Cpair * dp(icol,ilay));
            });
        }
    }
}
#endif
