#include <catch2/catch.hpp>

#include "share/ekat_pack.hpp"
#include "share/field/field.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "share/util/scream_test_utils.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "mpi/Comm.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/MpiContext.hpp"

#include <random>
#include <numeric>

namespace {

TEST_CASE("remap", "") {

  // +---------------------------------------------+
  // |  We assume four elements, and the global id |
  // |  of columns is as follows                   |
  // +---------------------------------------------+
  //
  //        Elem 0               Elem 1
  //
  //         (E)                  (S)
  //     15-16-17-12          12-13-14-15
  //      |  |  |  |           |  |  |  |
  //     28-29-30-31          31-34-35-28
  //  (N) |  |  |  | (S)   (E) |  |  |  | (W)
  //     24-25-26-27          27-32-33-24
  //      |  |  |  |           |  |  |  |
  //      3-22-23--0           0--1--2--3
  //         (W)                  (N)
  //  
  //         (E)                  (N)
  //      3 22-23--0           0--1--2--3
  //      |  |  |  |           |  |  |  |
  //      7-20-21--4           4--5--6--7
  //  (N) |  |  |  | (S)   (W) |  |  |  | (E)
  //     11-18-19--8           8--9-10-11
  //      |  |  |  |           |  |  |  |
  //     15-16-17-12          12-13-14-15
  //         (W)                  (S)
  //
  //        Elem 3               Elem 2
  //
  // so that, e.g., the W edge of elem 2 is shared with the S edge of elem 3
  // Note: This is a 'torus' configuration, where moving along the 'long' circle
  //       you see the element pairs 0-1 and 2-3, while moving along the 'short'
  //       circle you see the element pairs 0-3 and 1-2

  using namespace scream;
  using namespace ShortFieldTagsNames;

  // Some type defs
  using Device = DefaultDevice;
  using PackType = ekat::pack::Pack<Homme::Real,HOMMEXX_VECTOR_SIZE>;
  using Remapper = PhysicsDynamicsRemapper<Homme::Real,Device>;
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  std::random_device rd;
  const int seed = rd();
  std::mt19937_64 engine(seed);
  RPDF pdf(-1.0,1.0);

  // Get or create a comm
  Homme::Comm& comm = Homme::MpiContext::singleton().get_comm();
  comm.reset_mpi_comm(MPI_COMM_WORLD);

  const int world_size = comm.size();
  const int world_rank = comm.rank();
  if (world_size!=1 && world_size!=2 && world_size!=4) {
    // The test runs only with 1 or 2 ranks
    REQUIRE (false);
  }

  // Global counters
  constexpr int num_elems = 4;

  // Local counters
  const int num_local_elems = num_elems / world_size;
  int num_local_columns = -1;
  switch(world_size) {
    case 1:
      num_local_columns = 36;
      break;
    case 2:
      num_local_columns = (world_rank==0 ? 24 : 12);
      break;
    case 4:
      switch (world_rank) {
        case 0:
          num_local_columns = 16;
          break;
        case 1: // Fallthrough
        case 2:
          num_local_columns = 8;
          break;
        case 3:
          num_local_columns = 4;
          break;
      }
      break;
  }
  EKAT_REQUIRE_MSG(num_local_columns>0, "Internal test error! Fix homme_pd_remap_tests, please.\n");

  auto get_elem_pid = [&](const int elem_gid)->int {
    switch (world_size) {
      case 1: return 0;
      case 2: return elem_gid / 2;
      case 4: return elem_gid;
    }
    return -1;
  };
  auto get_elem_lid = [&](const int elem_gid)->int {
    return elem_gid - get_elem_pid(elem_gid)*num_local_elems;
  };

  // Create the columns global id mappings
  typename SEGrid::dofs_list_type p_dofs("p dofs",num_local_columns);
  typename SEGrid::dofs_list_type d_dofs("d dofs",num_local_elems*NP*NP);
  typename SEGrid::dofs_map_type p2d("p2d",num_local_columns);
  typename SEGrid::dofs_map_type d2p("d2p",num_local_elems*NP*NP);

  auto h_p2d = Kokkos::create_mirror_view(p2d);
  auto h_d2p = Kokkos::create_mirror_view(d2p);
  auto h_p_dofs = Kokkos::create_mirror_view(p_dofs);
  auto h_d_dofs = Kokkos::create_mirror_view(d_dofs);

  int col_ids[4][16] = { 
                         { 0, 27, 31, 12, 23, 26, 30, 17, 22, 25, 29, 16,  3, 24, 28, 15},
                         {15, 14, 13, 12, 28, 35, 34, 31, 24, 33, 32, 27,  3,  2,  1,  0},
                         {12, 13, 14, 15,  8,  9, 10, 11,  4,  5,  6,  7,  0,  1,  2,  3},
                         {12,  8,  4,  0, 17, 19, 21, 23, 16, 18, 20, 22, 15, 11,  7,  3}
                       };

  std::vector<int> my_elems_gids (num_local_elems);
  for (int i=0; i<num_local_elems; ++i) {
    my_elems_gids[i] = world_rank*num_local_elems + i;
  }

  for (int ie=0,col_lid=0; ie<num_local_elems; ++ie) {
    int elem_gid = ie + num_local_elems*world_rank;
    for (int i=0; i<NP; ++i) {
      for (int j=0; j<NP; ++j) {
        const int idof = ie*NP*NP+i*NP+j;
        const long col_gid = col_ids[elem_gid][i*NP+j];

        h_d2p(idof,0) = ie;
        h_d2p(idof,1) = i;
        h_d2p(idof,2) = j;
        h_d_dofs(idof) = col_gid;
        
        // Check if another element already owns the column
        bool found = false;
        for (int prev_elgid=0; prev_elgid<elem_gid; ++prev_elgid) {
          auto start = &col_ids[prev_elgid][0];
          auto end   = start+16;
          auto it = std::find(start,end,col_gid);
          if (it!=end) {
            found = true;
            break;
          }
        }
        if (!found) {
          // I own this column
          h_p2d(col_lid,0) = ie;
          h_p2d(col_lid,1) = i;
          h_p2d(col_lid,2) = j;
          h_p_dofs(col_lid) = col_gid;
          ++col_lid;
        }
      }
    }
  }

  Kokkos::deep_copy(p_dofs,h_p_dofs);
  Kokkos::deep_copy(d_dofs,h_d_dofs);
  Kokkos::deep_copy(p2d,h_p2d);
  Kokkos::deep_copy(d2p,h_d2p);

  // Create the physics and dynamics grids
  auto phys_grid = std::make_shared<SEGrid>(p2d,p_dofs,"Physics",GridType::SE_NodeBased);
  auto dyn_grid  = std::make_shared<SEGrid>(d2p,d_dofs,"Dynamics",GridType::SE_CellBased);

  // Create connectivity, and add the only connection
  // Note: catch2 runs this routine several times, but the connectivity can only be init-ed once, so check first.
  std::shared_ptr<Homme::Connectivity> connectivity = Homme::MpiContext::singleton().get_connectivity();

  constexpr auto NORTH = Homme::etoi(Homme::ConnectionName::NORTH);
  constexpr auto SOUTH = Homme::etoi(Homme::ConnectionName::SOUTH);
  constexpr auto EAST  = Homme::etoi(Homme::ConnectionName::EAST);
  constexpr auto WEST  = Homme::etoi(Homme::ConnectionName::WEST);
  constexpr auto NEAST = Homme::etoi(Homme::ConnectionName::NEAST);
  constexpr auto NWEST = Homme::etoi(Homme::ConnectionName::NWEST);
  constexpr auto SEAST = Homme::etoi(Homme::ConnectionName::SEAST);
  constexpr auto SWEST = Homme::etoi(Homme::ConnectionName::SWEST);

  if (!connectivity->is_finalized()) {
    connectivity->set_comm(comm);
    connectivity->set_num_elements(num_local_elems);
    if (get_elem_pid(0)==world_rank) {
      connectivity->add_connection(get_elem_lid(0),0,SOUTH,get_elem_pid(0),
                                   get_elem_lid(1),1,EAST ,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(0),0,NORTH,get_elem_pid(0),
                                   get_elem_lid(1),1,WEST ,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(0),0,WEST ,get_elem_pid(0),
                                   get_elem_lid(3),3,EAST ,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(0),0,EAST ,get_elem_pid(0),
                                   get_elem_lid(3),3,WEST ,get_elem_pid(3));

      connectivity->add_connection(get_elem_lid(0),0,SEAST,get_elem_pid(0),
                                   get_elem_lid(2),2,SWEST,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,NWEST,get_elem_pid(0),
                                   get_elem_lid(2),2,NEAST,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,SWEST,get_elem_pid(0),
                                   get_elem_lid(2),2,NWEST,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,NEAST,get_elem_pid(0),
                                   get_elem_lid(2),2,SEAST,get_elem_pid(2));
    }
    if (get_elem_pid(1)==world_rank) {
      connectivity->add_connection(get_elem_lid(1),1,EAST ,get_elem_pid(1),
                                   get_elem_lid(0),0,SOUTH,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(1),1,WEST ,get_elem_pid(1),
                                   get_elem_lid(0),0,NORTH,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(1),1,NORTH,get_elem_pid(1),
                                   get_elem_lid(2),2,NORTH,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(1),1,SOUTH,get_elem_pid(1),
                                   get_elem_lid(2),2,SOUTH,get_elem_pid(2));

      connectivity->add_connection(get_elem_lid(1),1,SEAST,get_elem_pid(1),
                                   get_elem_lid(3),3,SWEST,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,NWEST,get_elem_pid(1),
                                   get_elem_lid(3),3,NEAST,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,NEAST,get_elem_pid(1),
                                   get_elem_lid(3),3,SEAST,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,SWEST,get_elem_pid(1),
                                   get_elem_lid(3),3,NWEST,get_elem_pid(3));
    }
    if (get_elem_pid(2)==world_rank) {
      connectivity->add_connection(get_elem_lid(2),2,NORTH,get_elem_pid(2),
                                   get_elem_lid(1),1,NORTH,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(2),2,SOUTH,get_elem_pid(2),
                                   get_elem_lid(1),1,SOUTH,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(2),2,WEST ,get_elem_pid(2),
                                   get_elem_lid(3),3,SOUTH,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(2),2,EAST ,get_elem_pid(2),
                                   get_elem_lid(3),3,NORTH,get_elem_pid(3));

      connectivity->add_connection(get_elem_lid(2),2,NWEST,get_elem_pid(2),
                                   get_elem_lid(0),0,SWEST,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,SEAST,get_elem_pid(2),
                                   get_elem_lid(0),0,NEAST,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,NEAST,get_elem_pid(2),
                                   get_elem_lid(0),0,NWEST,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,SWEST,get_elem_pid(2),
                                   get_elem_lid(0),0,SEAST,get_elem_pid(0));
    }
    if (get_elem_pid(3)==world_rank) {
      connectivity->add_connection(get_elem_lid(3),3,EAST ,get_elem_pid(3),
                                   get_elem_lid(0),0,WEST ,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(3),3,WEST ,get_elem_pid(3),
                                   get_elem_lid(0),0,EAST ,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(3),3,SOUTH,get_elem_pid(3),
                                   get_elem_lid(2),2,WEST ,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(3),3,NORTH,get_elem_pid(3),
                                   get_elem_lid(2),2,EAST ,get_elem_pid(2));

      connectivity->add_connection(get_elem_lid(3),3,SEAST,get_elem_pid(3),
                                   get_elem_lid(1),1,NEAST,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,NWEST,get_elem_pid(3),
                                   get_elem_lid(1),1,SWEST,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,NEAST,get_elem_pid(3),
                                   get_elem_lid(1),1,NWEST,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,SWEST,get_elem_pid(3),
                                   get_elem_lid(1),1,SEAST,get_elem_pid(1));
    }
    connectivity->finalize(/* sanity check = */ false);
  }

  constexpr int NVL = HOMMEXX_NUM_PHYSICAL_LEV;
  constexpr int NTL = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int NTQ = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int NQ  = HOMMEXX_QSIZE_D;
  const int ne = num_local_elems;
  const int nc = num_local_columns;
  const auto units = units::m;  // Placeholder units (we don't care about units here)

  // Create tags and dimensions
  std::vector<FieldTag> scalar_2d_dyn_tags        = {EL,          GP, GP    };
  std::vector<FieldTag> vector_2d_dyn_tags        = {EL,     CMP, GP, GP    };
  std::vector<FieldTag> scalar_3d_dyn_tags        = {EL,          GP, GP, VL};
  std::vector<FieldTag> vector_3d_dyn_tags        = {EL,     CMP, GP, GP, VL};
  std::vector<FieldTag> scalar_state_3d_dyn_tags  = {EL, TL,      GP, GP, VL};
  std::vector<FieldTag> vector_state_3d_dyn_tags  = {EL, TL, CMP, GP, GP, VL};

  std::vector<FieldTag> scalar_2d_phys_tags       = {COL         };
  std::vector<FieldTag> vector_2d_phys_tags       = {COL, CMP    };
  std::vector<FieldTag> scalar_3d_phys_tags       = {COL,      VL};
  std::vector<FieldTag> vector_3d_phys_tags       = {COL, CMP, VL};
  std::vector<FieldTag> vector_state_3d_phys_tags = {COL, CMP, VL};
  std::vector<FieldTag> scalar_state_3d_phys_tags = {COL,      VL};

  std::vector<int> scalar_2d_dyn_dims        = {ne,          NP, NP     };
  std::vector<int> vector_2d_dyn_dims        = {ne,       2, NP, NP     };
  std::vector<int> scalar_3d_dyn_dims        = {ne,          NP, NP, NVL};
  std::vector<int> vector_3d_dyn_dims        = {ne,       2, NP, NP, NVL};
  std::vector<int> scalar_state_3d_dyn_dims  = {ne, NTL,     NP, NP, NVL};
  std::vector<int> vector_state_3d_dyn_dims  = {ne, NTL,  2, NP, NP, NVL};
  std::vector<int> tracer_state_3d_dyn_dims  = {ne, NTQ, NQ, NP, NP, NVL};

  std::vector<int> scalar_2d_phys_dims       = {nc              };
  std::vector<int> vector_2d_phys_dims       = {nc,       2     };
  std::vector<int> scalar_3d_phys_dims       = {nc,          NVL};
  std::vector<int> vector_3d_phys_dims       = {nc,       2, NVL};
  std::vector<int> scalar_state_3d_phys_dims = {nc,     NVL};
  std::vector<int> vector_state_3d_phys_dims = {nc,  2, NVL};
  std::vector<int> tracer_state_3d_phys_dims = {nc, NQ, NVL};

  // Create identifiers
  FieldIdentifier scalar_2d_dyn_fid  ("scalar_2d_dynamics", FieldLayout(scalar_2d_dyn_tags,  scalar_2d_dyn_dims), units, dyn_grid->name());
  FieldIdentifier vector_2d_dyn_fid  ("vector_2d_dynamics", FieldLayout(vector_2d_dyn_tags,  vector_2d_dyn_dims), units, dyn_grid->name());
  FieldIdentifier scalar_3d_dyn_fid  ("scalar_3d_dynamics", FieldLayout(scalar_3d_dyn_tags,  scalar_3d_dyn_dims), units, dyn_grid->name());
  FieldIdentifier vector_3d_dyn_fid  ("vector_3d_dynamics", FieldLayout(vector_3d_dyn_tags,  vector_3d_dyn_dims), units, dyn_grid->name());
  FieldIdentifier scalar_state_3d_dyn_fid ("scalar_state_3d_dynamics", FieldLayout(scalar_state_3d_dyn_tags, scalar_state_3d_dyn_dims),units, dyn_grid->name());
  FieldIdentifier vector_state_3d_dyn_fid ("vector_state_3d_dynamics", FieldLayout(vector_state_3d_dyn_tags, vector_state_3d_dyn_dims),units, dyn_grid->name());
  FieldIdentifier tracer_state_3d_dyn_fid ("tracer_state_3d_dynamics", FieldLayout(vector_state_3d_dyn_tags, tracer_state_3d_dyn_dims),units, dyn_grid->name());

  FieldIdentifier scalar_2d_phys_fid ("scalar_2d_physics",  FieldLayout(scalar_2d_phys_tags, scalar_2d_phys_dims),units, phys_grid->name());
  FieldIdentifier vector_2d_phys_fid ("vector_2d_physics",  FieldLayout(vector_2d_phys_tags, vector_2d_phys_dims),units, phys_grid->name());
  FieldIdentifier scalar_3d_phys_fid ("scalar_3d_physics",  FieldLayout(scalar_3d_phys_tags, scalar_3d_phys_dims),units, phys_grid->name());
  FieldIdentifier vector_3d_phys_fid ("vector_3d_physics",  FieldLayout(vector_3d_phys_tags, vector_3d_phys_dims),units, phys_grid->name());
  FieldIdentifier scalar_state_3d_phys_fid ("scalar_state_3d_physics", FieldLayout(scalar_state_3d_phys_tags, scalar_state_3d_phys_dims),units, phys_grid->name());
  FieldIdentifier vector_state_3d_phys_fid ("vector_state_3d_physics", FieldLayout(vector_state_3d_phys_tags, vector_state_3d_phys_dims),units, phys_grid->name());
  FieldIdentifier tracer_state_3d_phys_fid ("tracer_state_3d_physics", FieldLayout(vector_state_3d_phys_tags, tracer_state_3d_phys_dims),units, phys_grid->name());

  // Create fields
  Field<Real,Device> scalar_2d_field_phys (scalar_2d_phys_fid);
  Field<Real,Device> vector_2d_field_phys (vector_2d_phys_fid);
  Field<Real,Device> scalar_3d_field_phys (scalar_3d_phys_fid);
  Field<Real,Device> vector_3d_field_phys (vector_3d_phys_fid);
  Field<Real,Device> scalar_state_3d_field_phys (scalar_state_3d_phys_fid);
  Field<Real,Device> vector_state_3d_field_phys (vector_state_3d_phys_fid);
  Field<Real,Device> tracer_state_3d_field_phys (tracer_state_3d_phys_fid);

  Field<Real,Device> scalar_2d_field_dyn(scalar_2d_dyn_fid);
  Field<Real,Device> vector_2d_field_dyn(vector_2d_dyn_fid);
  Field<Real,Device> scalar_3d_field_dyn(scalar_3d_dyn_fid);
  Field<Real,Device> vector_3d_field_dyn(vector_3d_dyn_fid);
  Field<Real,Device> scalar_state_3d_field_dyn (scalar_state_3d_dyn_fid);
  Field<Real,Device> vector_state_3d_field_dyn (vector_state_3d_dyn_fid);
  Field<Real,Device> tracer_state_3d_field_dyn (tracer_state_3d_dyn_fid);

  // Request allocation to fit packs of reals
  scalar_2d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_2d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  scalar_3d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_3d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  scalar_state_3d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_state_3d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  tracer_state_3d_field_phys.get_header().get_alloc_properties().request_value_type_allocation<PackType>();

  scalar_2d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_2d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  scalar_3d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_3d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  scalar_state_3d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  vector_state_3d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  tracer_state_3d_field_dyn.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
  tracer_state_3d_field_dyn.get_header().set_extra_data("Is Tracer State", true);

  // Allocate view
  scalar_2d_field_phys.allocate_view();
  vector_2d_field_phys.allocate_view();
  scalar_3d_field_phys.allocate_view();
  vector_3d_field_phys.allocate_view();
  scalar_state_3d_field_phys.allocate_view();
  vector_state_3d_field_phys.allocate_view();
  tracer_state_3d_field_phys.allocate_view();

  scalar_2d_field_dyn.allocate_view();
  vector_2d_field_dyn.allocate_view();
  scalar_3d_field_dyn.allocate_view();
  vector_3d_field_dyn.allocate_view();
  scalar_state_3d_field_dyn.allocate_view();
  vector_state_3d_field_dyn.allocate_view();
  tracer_state_3d_field_dyn.allocate_view();

  // Build the remapper, and register the fields
  std::shared_ptr<Remapper> remapper(new Remapper(phys_grid,dyn_grid));
  remapper->registration_begins();
  remapper->register_field(scalar_2d_field_phys, scalar_2d_field_dyn);
  remapper->register_field(vector_2d_field_phys, vector_2d_field_dyn);
  remapper->register_field(scalar_3d_field_phys, scalar_3d_field_dyn);
  remapper->register_field(vector_3d_field_phys, vector_3d_field_dyn);
  remapper->register_field(scalar_state_3d_field_phys, scalar_state_3d_field_dyn);
  remapper->register_field(vector_state_3d_field_phys, vector_state_3d_field_dyn);
  remapper->register_field(tracer_state_3d_field_phys, tracer_state_3d_field_dyn);
  remapper->registration_ends();

  const int np1 = IPDF(0,NTL-1)(engine);
  const int np1_qdp = IPDF(0,NTQ-1)(engine);
  {
    auto& tl = Homme::Context::singleton().get_time_level();
    tl.np1 = np1;
    tl.nm1 = (np1+1) % NTL;
    tl.n0  = (np1+2) % NTL;
    tl.np1_qdp = np1_qdp;
    tl.n0_qdp  = (np1_qdp+1) % NTQ;
  }

  SECTION ("remap") {

    for (bool fwd : {true, false}) {

      // Generate random numbers
      if (fwd) {
        ekat::util::genRandArray(scalar_2d_field_phys,  engine, pdf);
        ekat::util::genRandArray(vector_2d_field_phys,  engine, pdf);
        ekat::util::genRandArray(scalar_3d_field_phys,  engine, pdf);
        ekat::util::genRandArray(vector_3d_field_phys,  engine, pdf);

        ekat::util::genRandArray(scalar_state_3d_field_phys,  engine, pdf);
        ekat::util::genRandArray(vector_state_3d_field_phys,  engine, pdf);
        ekat::util::genRandArray(tracer_state_3d_field_phys,  engine, pdf);
      } else {
        // Note: for the dyn->phys test to run correctly, the dynamics input vector must be synced,
        //       meaning that the values at the interface between two elements must match.
        //       To do this, we initialize each entry in the dynamic vector with the id
        //       of the corresponding column.
        auto scalar_2d_view = scalar_2d_field_dyn.get_reshaped_view<Homme::Real*   [NP][NP]     >();
        auto vector_2d_view = vector_2d_field_dyn.get_reshaped_view<Homme::Real*[2][NP][NP]     >();
        auto scalar_3d_view = scalar_3d_field_dyn.get_reshaped_view<Homme::Real*   [NP][NP][NVL]>();
        auto vector_3d_view = vector_3d_field_dyn.get_reshaped_view<Homme::Real*[2][NP][NP][NVL]>();
        auto scalar_state_3d_view = scalar_state_3d_field_dyn.get_reshaped_view<Homme::Real*[NTL]    [NP][NP][NVL]>();
        auto vector_state_3d_view = vector_state_3d_field_dyn.get_reshaped_view<Homme::Real*[NTL] [2][NP][NP][NVL]>();
        auto tracer_state_3d_view = tracer_state_3d_field_dyn.get_reshaped_view<Homme::Real*[NTQ][NQ][NP][NP][NVL]>();

        auto h_scalar_2d_view = Kokkos::create_mirror_view(scalar_2d_view);
        auto h_vector_2d_view = Kokkos::create_mirror_view(vector_2d_view);
        auto h_scalar_3d_view = Kokkos::create_mirror_view(scalar_3d_view);
        auto h_vector_3d_view = Kokkos::create_mirror_view(vector_3d_view);
        auto h_scalar_state_3d_view = Kokkos::create_mirror_view(scalar_state_3d_view);
        auto h_vector_state_3d_view = Kokkos::create_mirror_view(vector_state_3d_view);
        auto h_tracer_state_3d_view = Kokkos::create_mirror_view(tracer_state_3d_view);
        for (int ie=0; ie<num_local_elems; ++ie) {
          for (int ip=0; ip<NP; ++ip) {
            for (int jp=0; jp<NP; ++jp) {
              const int idof = ie*NP*NP + ip*NP + jp;
              h_scalar_2d_view(ie,ip,jp) = h_d_dofs(idof);
              h_vector_2d_view(ie,0,ip,jp) = h_d_dofs(idof);
              h_vector_2d_view(ie,1,ip,jp) = h_d_dofs(idof);
              for (int il=0; il<NVL; ++ il) {
                h_scalar_3d_view(ie,ip,jp,il) = h_d_dofs(idof);
                h_vector_3d_view(ie,0,ip,jp,il) = h_d_dofs(idof);
                h_vector_3d_view(ie,1,ip,jp,il) = h_d_dofs(idof);

                for (int itl=0; itl<NTL; ++itl) {
                  h_scalar_state_3d_view(ie,itl,ip,jp,il) = h_d_dofs(idof);
                  h_vector_state_3d_view(ie,itl,0,ip,jp,il) = h_d_dofs(idof);
                  h_vector_state_3d_view(ie,itl,1,ip,jp,il) = h_d_dofs(idof);
                }
                for (int itl=0; itl<NTQ; ++itl) {
                  for (int iq=0; iq<NTQ; ++iq) {
                    h_tracer_state_3d_view(ie,itl,iq,ip,jp,il) = h_d_dofs(idof);
                  }
                }
              }
            }
          }
        }
        Kokkos::deep_copy(scalar_2d_view,h_scalar_2d_view);
        Kokkos::deep_copy(vector_2d_view,h_vector_2d_view);
        Kokkos::deep_copy(scalar_3d_view,h_scalar_3d_view);
        Kokkos::deep_copy(vector_3d_view,h_vector_3d_view);

        Kokkos::deep_copy(scalar_state_3d_view,h_scalar_state_3d_view);
        Kokkos::deep_copy(vector_state_3d_view,h_vector_state_3d_view);
        Kokkos::deep_copy(tracer_state_3d_view,h_tracer_state_3d_view);
      }

      // Remap
      remapper->remap(fwd);

      // Check
      {
        // 2d scalar
        auto phys = Kokkos::create_mirror_view(scalar_2d_field_phys.template get_reshaped_view<Homme::Real*>());
        auto dyn = Kokkos::create_mirror_view(scalar_2d_field_dyn.template get_reshaped_view<Homme::Real***>());
        Kokkos::deep_copy(phys,scalar_2d_field_phys.template get_reshaped_view<Homme::Real*>());
        Kokkos::deep_copy(dyn,scalar_2d_field_dyn.template get_reshaped_view<Homme::Real***>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          if (phys(icol)!=dyn(ie,ip,jp)) {
              printf("p_in(%d) = %2.16f\n",icol,phys(icol));
              printf("d_out(%d,%d,%d) = %2.16f\n",ie,ip,jp,dyn(ie,ip,jp));
          }
          REQUIRE (phys(icol)==dyn(ie,ip,jp));
        }
      }

      {
        // 2d vector
        auto phys = Kokkos::create_mirror_view(vector_2d_field_phys.template get_reshaped_view<Homme::Real**>());
        auto dyn = Kokkos::create_mirror_view(vector_2d_field_dyn.template get_reshaped_view<Homme::Real****>());
        Kokkos::deep_copy(phys,vector_2d_field_phys.template get_reshaped_view<Homme::Real**>());
        Kokkos::deep_copy(dyn,vector_2d_field_dyn.template get_reshaped_view<Homme::Real****>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int icomp=0; icomp<2; ++icomp) {
            if (phys(icol,icomp)!=dyn(ie,icomp,ip,jp)) {
                printf("p_in(%d,%d) = %2.16f\n",icol,icomp,phys(icol,icomp));
                printf("d_out(%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,dyn(ie,icomp,ip,jp));
            }
            REQUIRE (phys(icol,icomp)==dyn(ie,icomp,ip,jp));
          }
        }
      }

      {
        // 3d scalar
        auto phys = Kokkos::create_mirror_view(scalar_3d_field_phys.template get_reshaped_view<Homme::Real**>());
        auto dyn = Kokkos::create_mirror_view(scalar_3d_field_dyn.template get_reshaped_view<Homme::Real****>());
        Kokkos::deep_copy(phys,scalar_3d_field_phys.template get_reshaped_view<Homme::Real**>());
        Kokkos::deep_copy(dyn,scalar_3d_field_dyn.template get_reshaped_view<Homme::Real****>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int ilev=0; ilev<NVL; ++ilev) {
            if (phys(icol,ilev)!=dyn(ie,ip,jp,ilev)) {
                printf("p_in(%d,%d) = %2.16f\n",icol,ilev,phys(icol,ilev));
                printf("d_out(%d,%d,%d,%d) = %2.16f\n",ie,ip,jp,ilev,dyn(ie,ip,jp,ilev));
            }
            REQUIRE (phys(icol,ilev)==dyn(ie,ip,jp,ilev));
          }
        }
      }

      {
        // 3d vector
        auto phys = Kokkos::create_mirror_view(vector_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        auto dyn = Kokkos::create_mirror_view(vector_3d_field_dyn.template get_reshaped_view<Homme::Real*****>());
        Kokkos::deep_copy(phys,vector_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        Kokkos::deep_copy(dyn,vector_3d_field_dyn.template get_reshaped_view<Homme::Real*****>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int icomp=0; icomp<2; ++icomp) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(icol,icomp,ilev)!=dyn(ie,icomp,ip,jp,ilev)) {
                  printf("p_in(%d,%d,%d) = %2.16f\n",icol,icomp,ilev,phys(icol,icomp,ilev));
                  printf("d_out(%d,%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,icomp,ip,jp,ilev));
              }
              REQUIRE (phys(icol,icomp,ilev)==dyn(ie,icomp,ip,jp,ilev));
            }
          }
        }
      }

      {
        // 3d scalar state
        const int itl = Homme::Context::singleton().get_time_level().np1;
        auto phys = Kokkos::create_mirror_view(scalar_state_3d_field_phys.template get_reshaped_view<Homme::Real**>());
        auto dyn = Kokkos::create_mirror_view(scalar_state_3d_field_dyn.template get_reshaped_view<Homme::Real*****>());
        Kokkos::deep_copy(phys,scalar_state_3d_field_phys.template get_reshaped_view<Homme::Real**>());
        Kokkos::deep_copy(dyn,scalar_state_3d_field_dyn.template get_reshaped_view<Homme::Real*****>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int ilev=0; ilev<NVL; ++ilev) {
            if (phys(icol,ilev)!=dyn(ie,itl,ip,jp,ilev)) {
                printf("p_in(%d,%d) = %2.16f\n",icol,ilev,phys(icol,ilev));
                printf("d_out(%d,%d,%d,%d) = %2.16f\n",ie,ip,jp,ilev,dyn(ie,itl,ip,jp,ilev));
            }
            REQUIRE (phys(icol,ilev)==dyn(ie,itl,ip,jp,ilev));
          }
        }
      }

      {
        // 3d vector state
        const int itl = Homme::Context::singleton().get_time_level().np1;
        auto phys = Kokkos::create_mirror_view(vector_state_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        auto dyn = Kokkos::create_mirror_view(vector_state_3d_field_dyn.template get_reshaped_view<Homme::Real******>());
        Kokkos::deep_copy(phys,vector_state_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        Kokkos::deep_copy(dyn,vector_state_3d_field_dyn.template get_reshaped_view<Homme::Real******>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int icomp=0; icomp<2; ++icomp) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(icol,icomp,ilev)!=dyn(ie,itl,icomp,ip,jp,ilev)) {
                  printf("p_in(%d,%d,%d) = %2.16f\n",icol,icomp,ilev,phys(icol,icomp,ilev));
                  printf("d_out(%d,%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,ilev,dyn(ie,itl,icomp,ip,jp,ilev));
              }
              REQUIRE (phys(icol,icomp,ilev)==dyn(ie,itl,icomp,ip,jp,ilev));
            }
          }
        }
      }

      {
        // 3d tracer state
        const int itl = Homme::Context::singleton().get_time_level().np1_qdp;
        auto phys = Kokkos::create_mirror_view(tracer_state_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        auto dyn = Kokkos::create_mirror_view(tracer_state_3d_field_dyn.template get_reshaped_view<Homme::Real******>());
        Kokkos::deep_copy(phys,tracer_state_3d_field_phys.template get_reshaped_view<Homme::Real***>());
        Kokkos::deep_copy(dyn,tracer_state_3d_field_dyn.template get_reshaped_view<Homme::Real******>());
        for (int icol=0; icol<num_local_columns; ++icol) {
          auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
          const int ie = elgp[0];
          const int ip = elgp[1];
          const int jp = elgp[2];
          for (int iq=0; iq<2; ++iq) {
            for (int ilev=0; ilev<NVL; ++ilev) {
              if (phys(icol,iq,ilev)!=dyn(ie,itl,iq,ip,jp,ilev)) {
                  printf("p_in(%d,%d,%d) = %2.16f\n",icol,iq,ilev,phys(icol,iq,ilev));
                  printf("d_out(%d,%d,%d,%d,%d) = %2.16f\n",ie,iq,ip,jp,ilev,dyn(ie,itl,iq,ip,jp,ilev));
              }
              REQUIRE (phys(icol,iq,ilev)==dyn(ie,itl,iq,ip,jp,ilev));
            }
          }
        }
      }
    }
  }

  // Delete remapper before finalizing the mpi context, since the remapper has some MPI stuff in it
  remapper = nullptr;

  // Finalize Homme::MpiContext (deletes buffers manager)
  Homme::MpiContext::finalize_singleton();
}

} // anonymous namespace
