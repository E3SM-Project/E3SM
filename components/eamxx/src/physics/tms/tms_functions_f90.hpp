#ifndef SCREAM_TMS_FUNCTIONS_F90_HPP
#define SCREAM_TMS_FUNCTIONS_F90_HPP

#include "share/eamxx_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "tms_functions.hpp"
#include "physics_constants.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of tms functions from C++
//

namespace scream {
namespace tms {

struct ComputeTMSData : public PhysicsTestData
{
  // Input
  int ncols, nlevs;
  Real *u_wind, *v_wind, *t_mid, *p_mid, *exner, *z_mid, *sgh, *landfrac;

  // Output
  Real *ksrf, *taux, *tauy;

  ComputeTMSData(int ncols_, int nlevs_)
   : PhysicsTestData({ {ncols_}, {ncols_, nlevs_} },
                     { {&sgh, &landfrac, &ksrf, &taux, &tauy},
                       {&u_wind, &v_wind, &t_mid, &p_mid, &exner, &z_mid} }),
    ncols(ncols_), nlevs(nlevs_)
  {}

  PTD_STD_DEF(ComputeTMSData, 2, ncols, nlevs);
};

// Glue functions to call fortran from from C++ with the Data struct
void compute_tms(ComputeTMSData& d);

// _f function decls
extern "C" {
void compute_tms_f(int ncols, int nlevs,
                   Real *u_wind, Real *v_wind, Real *t_mid, Real *p_mid, Real *exner, Real *z_mid,
                   Real *sgh, Real *landfrac, Real *ksrf, Real *taux, Real *tauy);
} // end _f function decls

}  // namespace tms
}  // namespace scream

#endif // SCREAM_TMS_FUNCTIONS_F90_HPP
