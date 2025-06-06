#ifndef RRTMGP_UTILS_HPP
#define RRTMGP_UTILS_HPP

#include "physics/share/physics_constants.hpp"
#include "cpp/rrtmgp_const.h"
#include "cpp/rrtmgp_conversion.h"

namespace scream {
namespace rrtmgp {

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
template<class View1, class View2, class View3, class View4>
void compute_heating_rate (
  View1 const &flux_up,
  View2 const &flux_dn,
  View3 const &pdel   ,
  View4 const &heating_rate)
{
  using physconst = scream::physics::Constants<Real>;
  using LayoutT = typename View1::array_layout;
  const int ncol = (int)flux_up.extent(0);
  const int nlay = (int)flux_up.extent(1)-1;
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
    heating_rate(icol,ilay) = (
      flux_up(icol,ilay+1) - flux_up(icol,ilay) -
      flux_dn(icol,ilay+1) + flux_dn(icol,ilay)
                               ) * physconst::gravit / (physconst::Cpair * pdel(icol,ilay));
                                  ));
}

inline bool radiation_do(const int rad_freq, const int nstep) {
  // If rad_freq == 0, then never do radiation;
  // Otherwise, we always call radiation at the first step,
  // and afterwards we do radiation if the timestep is divisible
  // by rad_freq
  if (rad_freq == 0) {
    return false;
  } else {
    return nstep % rad_freq == 0;
  }
}

// Verify that array only contains values within valid range, and if not
// report min and max of array
template <class T, typename std::enable_if<T::rank == 1>::type* dummy = nullptr>
bool check_range_k(T x, typename T::const_value_type xmin, typename T::const_value_type xmax,
                   std::string msg, std::ostream& out=std::cout) {
  bool pass = true;
  auto _xmin = conv::minval(x);
  auto _xmax = conv::maxval(x);
  if (_xmin < xmin or _xmax > xmax) {
    // How many outside range?
    Kokkos::View<bool*> bad_mask("bad_mask", x.extent(0));
    Kokkos::parallel_for(x.extent(0), KOKKOS_LAMBDA (int i) {
      if (x(i) < xmin or x(i) > xmax) {
        bad_mask(i) = true;
      }
    });
    auto num_bad = conv::sum(bad_mask);
    pass = false;
    out << msg << ": "
        << num_bad << " values outside range "
        << "[" << xmin << "," << xmax << "]"
        << "; minval = " << _xmin
        << "; maxval = " << _xmax << "\n";
  }
  return pass;
}

template <class T, typename std::enable_if<T::rank == 2>::type* dummy = nullptr>
bool check_range_k(T x, typename T::const_value_type xmin, typename T::const_value_type xmax,
                   std::string msg, std::ostream& out=std::cout) {
  bool pass = true;
  auto _xmin = conv::minval(x);
  auto _xmax = conv::maxval(x);
  if (_xmin < xmin or _xmax > xmax) {
    // How many outside range?
    Kokkos::View<bool**> bad_mask("bad_mask", x.extent(0), x.extent(1));
    Kokkos::parallel_for(x.extent(0), KOKKOS_LAMBDA (int i) {
      for (size_t j = 0; j < x.extent(1); ++j) {
        if (x(i, j) < xmin or x(i, j) > xmax) {
          bad_mask(i, j) = true;
        }
      }
    });
    auto num_bad = conv::sum(bad_mask);
    if (num_bad > 0) {
      pass = false;
      out << msg << ": "
          << num_bad << " values outside range "
          << "[" << xmin << "," << xmax << "]"
          << "; minval = " << _xmin
          << "; maxval = " << _xmax << "\n";
    }
  }
  return pass;
}

template <class T, typename std::enable_if<T::rank == 3>::type* dummy = nullptr>
bool check_range_k(T x, typename T::const_value_type xmin, typename T::const_value_type xmax,
                   std::string msg, std::ostream& out=std::cout) {
  bool pass = true;
  auto _xmin = conv::minval(x);
  auto _xmax = conv::maxval(x);
  if (_xmin < xmin or _xmax > xmax) {
    // How many outside range?
    Kokkos::View<bool***> bad_mask("bad_mask", x.extent(0), x.extent(1), x.extent(2));
    Kokkos::parallel_for(x.extent(0), KOKKOS_LAMBDA (int i) {
      for (size_t j = 0; j < x.extent(1); ++j) {
        for (size_t k = 0; k < x.extent(2); ++k) {
          if (x(i, j, k) < xmin or x(i, j, k) > xmax) {
            bad_mask(i, j, k) = true;
          }
        }
      }
    });
    auto num_bad = conv::sum(bad_mask);
    if (num_bad > 0) {
      pass = false;
      out << msg << ": "
          << num_bad << " values outside range "
          << "[" << xmin << "," << xmax << "]"
          << "; minval = " << _xmin
          << "; maxval = " << _xmax << "\n";
    }
  }
  return pass;
}

} // namespace rrtmgp
} // namespace scream

#endif
