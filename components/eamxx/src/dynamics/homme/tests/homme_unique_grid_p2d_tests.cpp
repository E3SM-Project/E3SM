#include <catch2/catch.hpp>

#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/interface/eamxx_homme_interface.hpp"
#include "share/remap/unique_grid_remapper.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"

#include "SimulationParams.hpp"

#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_grids_manager.hpp"

#include <ekat_comm.hpp>
#include <ekat_test_utils.hpp>

#include <chrono>
#include <cstdlib>
#include <memory>
#include <string>

extern "C" {
// These are specific C/F calls for these tests (i.e., not part of eamxx_homme_interface.hpp)
void init_test_params_f90 ();
void cleanup_test_f90 ();
}

namespace {

/*
 * Compare the two ways of doing a phys->dyn remap on the same (real) Homme grids:
 *  - PhysicsDynamicsRemapper (fwd): uses Homme's BoundaryExchange to sync the shared dofs
 *  - UniqueGridRemapper (bwd): uses a p2p comm schedule built from the gids
 * The two must give the exact same result. We also print the time spent in setup and in the
 * remap, as a very rough comparison (set the env var HOMME_P2D_NE to run on a larger grid,
 * and HOMME_P2D_ITERS to change the number of remap calls that are timed).
 */

int get_env_int (const char* name, const int def) {
  const char* v = std::getenv(name);
  return v ? std::atoi(v) : def;
}

using clock_type = std::chrono::steady_clock;
double seconds_since (const clock_type::time_point& t0) {
  return std::chrono::duration<double>(clock_type::now()-t0).count();
}

TEST_CASE("homme_unique_grid_p2d", "") {
  using namespace scream;
  using namespace ShortFieldTagsNames;

  constexpr int pg_gll = 0;
  constexpr int NQ  = HOMMEXX_QSIZE_D;
  constexpr int PackSize = HOMMEXX_VECTOR_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  auto& c = Homme::Context::singleton();
  auto& sp = c.create<Homme::SimulationParams>();
  sp.qsize = NQ;

  const int ne    = get_env_int("HOMME_P2D_NE",2);
  const int iters = get_env_int("HOMME_P2D_ITERS",3);
  set_homme_param("ne",ne);

  ekat::ParameterList params;
  params.set<std::string>("physics_grid_type","gll");
  params.set<std::string>("vertical_coordinate_filename","NONE");
  HommeGridsManager gm(comm,params);
  gm.build_grids();

  auto phys_grid = gm.get_grid("physics_gll");
  auto dyn_grid  = gm.get_grid("dynamics");
  EKAT_REQUIRE_MSG(phys_grid->get_num_local_dofs()>0, "Internal test error! Fix homme_unique_grid_p2d, please.\n");

  const auto units = ekat::units::m;
  auto make_field = [&](const std::string& name, const FieldLayout& l, const std::string& gn) {
    Field f(FieldIdentifier(name,l,units,gn));
    // Only 3d fields are packed (padding a 2d field would break the static extents that Homme expects)
    if (l.has_tag(LEV)) {
      f.get_header().get_alloc_properties().request_allocation(PackSize);
    }
    f.allocate_view();
    return f;
  };

  // Fields like the ones remapped in HommeDynamics (2d/3d, scalar/vector, tracers)
  const std::vector<std::string> names = {"s2d","v2d","s3d","v3d","tr3d"};
  auto get_layout = [&](const AbstractGrid& g, const std::string& n) {
    if (n=="s2d")  return g.get_2d_scalar_layout();
    if (n=="v2d")  return g.get_2d_vector_layout(2,"vec");
    if (n=="s3d")  return g.get_3d_scalar_layout(LEV);
    if (n=="v3d")  return g.get_3d_vector_layout(LEV,2,"vec");
    return g.get_3d_vector_layout(LEV,NQ,"tracers");
  };

  // Phys fields (source), and two sets of dyn fields (one per remapper)
  std::vector<Field> phys, dyn_pd, dyn_ur;
  for (const auto& n : names) {
    phys.push_back  (make_field(n,get_layout(*phys_grid,n),phys_grid->name()));
    dyn_pd.push_back(make_field(n,get_layout(*dyn_grid,n),dyn_grid->name()));
    dyn_ur.push_back(make_field(n,get_layout(*dyn_grid,n),dyn_grid->name()));
  }

  // Fill phys fields with random numbers (the same values are given to both remappers)
  for (size_t i=0; i<phys.size(); ++i) {
    randomize_uniform(phys[i],comm.rank()*100+i+1);
  }

  // --- PhysicsDynamicsRemapper (BoundaryExchange) --- //
  auto t0 = clock_type::now();
  auto pd = std::make_shared<PhysicsDynamicsRemapper>(phys_grid,dyn_grid);
  for (size_t i=0; i<phys.size(); ++i) {
    pd->register_field(phys[i],dyn_pd[i]);
  }
  pd->registration_ends();
  double t_pd_setup = seconds_since(t0);

  // --- UniqueGridRemapper (p2p) --- //
  t0 = clock_type::now();
  auto ur = std::make_shared<UniqueGridRemapper>(dyn_grid);
  for (size_t i=0; i<phys.size(); ++i) {
    ur->register_field(dyn_ur[i],phys[i]);
  }
  ur->registration_ends();
  double t_up_setup = seconds_since(t0);

  // The first call of remap_bwd builds the comm schedule
  comm.barrier();
  t0 = clock_type::now();
  ur->remap_bwd();
  Kokkos::fence();
  double t_up_first = seconds_since(t0);

  // The first call of the PD remapper (for a fair comparison)
  comm.barrier();
  t0 = clock_type::now();
  pd->remap_fwd();
  Kokkos::fence();
  double t_pd_first = seconds_since(t0);

  // Timed calls
  double t_pd = 0, t_up = 0;
  for (int it=0; it<iters; ++it) {
    comm.barrier();
    t0 = clock_type::now();
    pd->remap_fwd();
    Kokkos::fence();
    t_pd += seconds_since(t0);

    comm.barrier();
    t0 = clock_type::now();
    ur->remap_bwd();
    Kokkos::fence();
    t_up += seconds_since(t0);
  }
  t_pd /= iters;
  t_up /= iters;

  // Same answer (exactly)
  for (size_t i=0; i<names.size(); ++i) {
    CHECK (views_are_equal(dyn_pd[i],dyn_ur[i]));
  }

  // Report the max time over ranks
  double times[6] = {t_pd_setup, t_up_setup, t_pd_first, t_up_first, t_pd, t_up};
  double max_times[6];
  comm.all_reduce(times,max_times,6,MPI_MAX);
  if (comm.am_i_root()) {
    printf("\n  p2d remap, ne=%d, ranks=%d, global cols=%d, iters=%d (max over ranks, in sec)\n",
           ne,comm.size(),(int)phys_grid->get_num_global_dofs(),iters);
    printf("               setup      1st call   avg call\n");
    printf("   PD (BEX)   %9.3e  %9.3e  %9.3e\n",max_times[0],max_times[2],max_times[4]);
    printf("   unique->grid p2p %9.3e  %9.3e  %9.3e\n\n",max_times[1],max_times[3],max_times[5]);
  }

  // Delete remappers before finalizing the mpi context, since PD has some MPI stuff in it
  pd = nullptr;
  ur = nullptr;

  Homme::Context::finalize_singleton();
  cleanup_test_f90();
}

} // anonymous namespace
