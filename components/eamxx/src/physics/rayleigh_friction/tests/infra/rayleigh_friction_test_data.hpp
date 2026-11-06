#ifndef SCREAM_RAYLEIGH_FRICTION_TEST_DATA_HPP
#define SCREAM_RAYLEIGH_FRICTION_TEST_DATA_HPP

#include "share/core/eamxx_types.hpp"
#include "share/physics/physics_test_data.hpp"

#include "rayleigh_friction_functions.hpp"

//
// Data structs and glue functions to call the C++ Rayleigh friction
// functions with host data in C layout
//

namespace scream {
namespace rayleigh_friction {

struct RayleighFrictionInitData : public PhysicsTestData
{
  // Input
  int nlevs, rayk0;
  Real raykrange, raytau0;

  // Output
  Real *otau;

  RayleighFrictionInitData(int nlevs_, int rayk0_, Real raykrange_, Real raytau0_)
   : PhysicsTestData({ {nlevs_} }, { {&otau} }),
    nlevs(nlevs_), rayk0(rayk0_), raykrange(raykrange_), raytau0(raytau0_)
  {}

  PTD_STD_DEF(RayleighFrictionInitData, 4, nlevs, rayk0, raykrange, raytau0);
};

struct RayleighFrictionTendData : public PhysicsTestData
{
  // Input
  int ncols, nlevs;
  Real dt;
  Real *otau;

  // Input/Output
  Real *u_wind, *v_wind, *t_mid;

  RayleighFrictionTendData(int ncols_, int nlevs_, Real dt_)
   : PhysicsTestData({ {nlevs_}, {ncols_, nlevs_} },
                     { {&otau}, {&u_wind, &v_wind, &t_mid} }),
    ncols(ncols_), nlevs(nlevs_), dt(dt_)
  {}

  PTD_STD_DEF(RayleighFrictionTendData, 3, ncols, nlevs, dt);
};

// Glue functions to call the C++ implementation with the Data structs
void rayleigh_friction_init(RayleighFrictionInitData& d);
void rayleigh_friction_tend(RayleighFrictionTendData& d);

}  // namespace rayleigh_friction
}  // namespace scream

#endif // SCREAM_RAYLEIGH_FRICTION_TEST_DATA_HPP
