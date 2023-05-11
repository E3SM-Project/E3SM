#ifndef RRTMGP_UTILS_HPP
#define RRTMGP_UTILS_HPP

#include "physics/share/physics_constants.hpp"
#include "cpp/rrtmgp_const.h"
#include "YAKL.h"
#include "YAKL_Bounds_fortran.h"

namespace scream {
    namespace rrtmgp {

        // Things we need from YAKL
        using yakl::intrinsics::maxval;
        using yakl::intrinsics::minval;
        using yakl::intrinsics::count;
        using yakl::intrinsics::sum;

        // Provide a routine to compute heating due to radiative fluxes. This is
        // computed as net flux into a layer, converted to a heating rate. It is
        // the responsibility of the user to ensure fields are passed with the
        // proper units. I.e., pressure at level interfaces should be in Pa,
        // fluxes in W m-2, Cpair in J kg-1 K-1, gravit in m s-2. This will give
        // heating in units of K s-1.
        // TODO: we should probably update this to use the pseudo-density pdel instead
        // of approximating pdel by differencing the level interface pressures.
        // We are leaving this for the time being for consistency with SCREAMv0,
        // from which this code was directly ported.
        template <class T, int myMem, int myStyle> void compute_heating_rate (
                yakl::Array<T,2,myMem,myStyle> const &flux_up, yakl::Array<T,2,myMem,myStyle> const &flux_dn, 
                yakl::Array<T,2,myMem,myStyle> const &pdel   , yakl::Array<T,2,myMem,myStyle> &heating_rate
            ) {
            using physconst = scream::physics::Constants<Real>;
            auto ncol = flux_up.dimension[0];
            auto nlay = flux_up.dimension[1]-1;
            yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
                heating_rate(icol,ilay) = (
                    flux_up(icol,ilay+1) - flux_up(icol,ilay) - 
                    flux_dn(icol,ilay+1) + flux_dn(icol,ilay)
                ) * physconst::gravit / (physconst::Cpair * pdel(icol,ilay));
            });
        }

        inline bool radiation_do(const int irad, const int nstep) {
            // If irad == 0, then never do radiation;
            // Otherwise, we always call radiation at the first step,
            // and afterwards we do radiation if the timestep is divisible
            // by irad
            if (irad == 0) {
                return false;
            } else {
                return ( (nstep == 0) || (nstep % irad == 0) );
            }
        }


        // Verify that array only contains values within valid range, and if not
        // report min and max of array
        template <class T> bool check_range(T x, Real xmin, Real xmax, std::string msg, std::ostream& out=std::cout) {
            bool pass = true;
            auto _xmin = minval(x);
            auto _xmax = maxval(x);
            if (_xmin < xmin or _xmax > xmax) {
                // How many outside range?
                auto bad_mask = x.createDeviceCopy();
                memset(bad_mask, 0);
                yakl::c::parallel_for(yakl::c::SimpleBounds<1>(x.totElems()), YAKL_LAMBDA (int i) {
                  if (x.data()[i] < xmin or x.data()[i] > xmax) {
                      bad_mask.data()[i] = 1;
                  }
                });
                auto num_bad = sum(bad_mask);
                pass = false;
                out << msg << ": "
                    << num_bad << " values outside range "
                    << "[" << xmin << "," << xmax << "]"
                    << "; minval = " << _xmin 
                    << "; maxval = " << _xmax << "\n";
            }
            return pass;
        }

    }

}
#endif
