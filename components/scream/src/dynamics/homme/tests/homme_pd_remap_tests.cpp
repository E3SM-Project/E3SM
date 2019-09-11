#include <catch2/catch.hpp>

#include "share/scream_pack.hpp"
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

  // Some type defs
  using Device = DefaultDevice;
  using PackType = pack::Pack<Homme::Real,HOMMEXX_VECTOR_SIZE>;
  using Remapper = PhysicsDynamicsRemapper<Homme::Real,Device>;
  using RPDF = std::uniform_real_distribution<Real>;

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
  scream_require_msg(num_local_columns>0, "Internal test error! Fix homme_pd_remap_tests, please.\n");

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
  typename SEGrid::dofs_map_type p2d("p2d",num_local_columns);
  typename SEGrid::dofs_map_type d2p("d2p",num_local_elems*NP*NP);

  auto h_p2d = Kokkos::create_mirror_view(p2d);
  auto h_d2p = Kokkos::create_mirror_view(d2p);

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
        auto col_gid = col_ids[elem_gid][i*NP+j];
        auto elgp = Kokkos::subview(h_d2p,ie*NP*NP+i*NP+j,Kokkos::ALL());
        elgp(0) = ie;
        elgp(1) = i;
        elgp(2) = j;
        elgp(3) = col_gid;
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
          h_p2d(col_lid,3) = col_gid;
          ++col_lid;
        }
      }
    }
  }

  Kokkos::deep_copy(p2d,h_p2d);
  Kokkos::deep_copy(d2p,h_d2p);

  // Create the physics and dynamics grids
  auto phys_grid = std::make_shared<SEGrid>(p2d,"Physics",GridType::SE_NodeBased);
  auto dyn_grid  = std::make_shared<SEGrid>(d2p,"Dynamics",GridType::SE_CellBased);

  // Create connectivity, and add the only connection
  // Note: catch2 runs this routine several times, but the connectivity can only be init-ed once, so check first.
  std::shared_ptr<Homme::Connectivity> connectivity = Homme::MpiContext::singleton().get_connectivity();

  constexpr auto NORTH = Homme::ConnectionName::NORTH;
  constexpr auto SOUTH = Homme::ConnectionName::SOUTH;
  constexpr auto EAST  = Homme::ConnectionName::EAST;
  constexpr auto WEST  = Homme::ConnectionName::WEST;
  constexpr auto NEAST = Homme::ConnectionName::NEAST;
  constexpr auto NWEST = Homme::ConnectionName::NWEST;
  constexpr auto SEAST = Homme::ConnectionName::SEAST;
  constexpr auto SWEST = Homme::ConnectionName::SWEST;

  if (!connectivity->is_finalized()) {
    connectivity->set_comm(comm);
    connectivity->set_num_elements(num_local_elems);
    if (get_elem_pid(0)==world_rank) {
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(SOUTH),get_elem_pid(0),
                                   get_elem_lid(1),1,Homme::etoi(EAST) ,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(NORTH),get_elem_pid(0),
                                   get_elem_lid(1),1,Homme::etoi(WEST) ,get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(WEST) ,get_elem_pid(0),
                                   get_elem_lid(3),3,Homme::etoi(EAST) ,get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(EAST) ,get_elem_pid(0),
                                   get_elem_lid(3),3,Homme::etoi(WEST) ,get_elem_pid(3));

      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(SEAST),get_elem_pid(0),
                                   get_elem_lid(2),2,Homme::etoi(SWEST),get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(NWEST),get_elem_pid(0),
                                   get_elem_lid(2),2,Homme::etoi(NEAST),get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(SWEST),get_elem_pid(0),
                                   get_elem_lid(2),2,Homme::etoi(NWEST),get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(0),0,Homme::etoi(NEAST),get_elem_pid(0),
                                   get_elem_lid(2),2,Homme::etoi(SEAST),get_elem_pid(2));
    }
    if (get_elem_pid(1)==world_rank) {
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(EAST) ,get_elem_pid(1),
                                   get_elem_lid(0),0,Homme::etoi(SOUTH),get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(WEST) ,get_elem_pid(1),
                                   get_elem_lid(0),0,Homme::etoi(NORTH),get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(NORTH),get_elem_pid(1),
                                   get_elem_lid(2),2,Homme::etoi(NORTH),get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(SOUTH),get_elem_pid(1),
                                   get_elem_lid(2),2,Homme::etoi(SOUTH),get_elem_pid(2));

      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(SEAST),get_elem_pid(1),
                                   get_elem_lid(3),3,Homme::etoi(SWEST),get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(NWEST),get_elem_pid(1),
                                   get_elem_lid(3),3,Homme::etoi(NEAST),get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(NEAST),get_elem_pid(1),
                                   get_elem_lid(3),3,Homme::etoi(SEAST),get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(1),1,Homme::etoi(SWEST),get_elem_pid(1),
                                   get_elem_lid(3),3,Homme::etoi(NWEST),get_elem_pid(3));
    }
    if (get_elem_pid(2)==world_rank) {
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(NORTH),get_elem_pid(2),
                                   get_elem_lid(1),1,Homme::etoi(NORTH),get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(SOUTH),get_elem_pid(2),
                                   get_elem_lid(1),1,Homme::etoi(SOUTH),get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(WEST) ,get_elem_pid(2),
                                   get_elem_lid(3),3,Homme::etoi(SOUTH),get_elem_pid(3));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(EAST) ,get_elem_pid(2),
                                   get_elem_lid(3),3,Homme::etoi(NORTH),get_elem_pid(3));

      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(NWEST),get_elem_pid(2),
                                   get_elem_lid(0),0,Homme::etoi(SWEST),get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(SEAST),get_elem_pid(2),
                                   get_elem_lid(0),0,Homme::etoi(NEAST),get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(NEAST),get_elem_pid(2),
                                   get_elem_lid(0),0,Homme::etoi(NWEST),get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(2),2,Homme::etoi(SWEST),get_elem_pid(2),
                                   get_elem_lid(0),0,Homme::etoi(SEAST),get_elem_pid(0));
    }
    if (get_elem_pid(3)==world_rank) {
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(EAST) ,get_elem_pid(3),
                                   get_elem_lid(0),0,Homme::etoi(WEST) ,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(WEST) ,get_elem_pid(3),
                                   get_elem_lid(0),0,Homme::etoi(EAST) ,get_elem_pid(0));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(SOUTH),get_elem_pid(3),
                                   get_elem_lid(2),2,Homme::etoi(WEST) ,get_elem_pid(2));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(NORTH),get_elem_pid(3),
                                   get_elem_lid(2),2,Homme::etoi(EAST) ,get_elem_pid(2));

      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(SEAST),get_elem_pid(3),
                                   get_elem_lid(1),1,Homme::etoi(NEAST),get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(NWEST),get_elem_pid(3),
                                   get_elem_lid(1),1,Homme::etoi(SWEST),get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(NEAST),get_elem_pid(3),
                                   get_elem_lid(1),1,Homme::etoi(NWEST),get_elem_pid(1));
      connectivity->add_connection(get_elem_lid(3),3,Homme::etoi(SWEST),get_elem_pid(3),
                                   get_elem_lid(1),1,Homme::etoi(SEAST),get_elem_pid(1));
    }
    connectivity->finalize(/* sanity check = */ false);
  }

  SECTION ("remap fwd") {

    // Create tags and dimensions
    std::vector<FieldTag> scalar_2d_dyn_tags  = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
    std::vector<FieldTag> scalar_2d_phys_tags = {FieldTag::Column};
    std::vector<FieldTag> vector_2d_dyn_tags  = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint};
    std::vector<FieldTag> vector_2d_phys_tags = {FieldTag::Column, FieldTag::Component};
    std::vector<FieldTag> scalar_3d_dyn_tags  = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
    std::vector<FieldTag> scalar_3d_phys_tags = {FieldTag::Column, FieldTag::VerticalLevel};
    std::vector<FieldTag> vector_3d_dyn_tags  = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
    std::vector<FieldTag> vector_3d_phys_tags = {FieldTag::Column, FieldTag::Component, FieldTag::VerticalLevel};

    std::vector<int> scalar_2d_dyn_dims  = {num_local_elems, NP, NP};
    std::vector<int> vector_2d_dyn_dims  = {num_local_elems, 2, NP, NP};
    std::vector<int> scalar_3d_dyn_dims  = {num_local_elems,    NP, NP, HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> vector_3d_dyn_dims  = {num_local_elems, 2, NP, NP, HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> scalar_2d_phys_dims = {num_local_columns};
    std::vector<int> vector_2d_phys_dims = {num_local_columns, 2};
    std::vector<int> scalar_3d_phys_dims = {num_local_columns,    HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> vector_3d_phys_dims = {num_local_columns, 2, HOMMEXX_NUM_PHYSICAL_LEV};

    // Create identifiers
    const auto units = units::m;  // Placeholder units
    FieldIdentifier scalar_2d_dyn_fid  ("scalar_2d_dynamics", FieldLayout(scalar_2d_dyn_tags,  scalar_2d_dyn_dims), units, dyn_grid->name());
    FieldIdentifier scalar_2d_phys_fid ("scalar_2d_physics",  FieldLayout(scalar_2d_phys_tags, scalar_2d_phys_dims),units,  phys_grid->name());
    FieldIdentifier vector_2d_dyn_fid  ("vector_2d_dynamics", FieldLayout(vector_2d_dyn_tags,  vector_2d_dyn_dims), units, dyn_grid->name());
    FieldIdentifier vector_2d_phys_fid ("vector_2d_physics",  FieldLayout(vector_2d_phys_tags, vector_2d_phys_dims),units,  phys_grid->name());
    FieldIdentifier scalar_3d_dyn_fid  ("scalar_3d_dynamics", FieldLayout(scalar_3d_dyn_tags,  scalar_3d_dyn_dims), units, dyn_grid->name());
    FieldIdentifier scalar_3d_phys_fid ("scalar_3d_physics",  FieldLayout(scalar_3d_phys_tags, scalar_3d_phys_dims),units,  phys_grid->name());
    FieldIdentifier vector_3d_dyn_fid  ("vector_3d_dynamics", FieldLayout(vector_3d_dyn_tags,  vector_3d_dyn_dims), units, dyn_grid->name());
    FieldIdentifier vector_3d_phys_fid ("vector_3d_physics",  FieldLayout(vector_3d_phys_tags, vector_3d_phys_dims),units,  phys_grid->name());

    // Create fields
    Field<Real,Device> scalar_2d_field_in (scalar_2d_phys_fid);
    Field<Real,Device> scalar_2d_field_out(scalar_2d_dyn_fid);
    Field<Real,Device> vector_2d_field_in (vector_2d_phys_fid);
    Field<Real,Device> vector_2d_field_out(vector_2d_dyn_fid);
    Field<Real,Device> scalar_3d_field_in (scalar_3d_phys_fid);
    Field<Real,Device> scalar_3d_field_out(scalar_3d_dyn_fid);
    Field<Real,Device> vector_3d_field_in (vector_3d_phys_fid);
    Field<Real,Device> vector_3d_field_out(vector_3d_dyn_fid);

    // Request allocation to fit packs of reals
    scalar_2d_field_in.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    scalar_2d_field_out.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    vector_2d_field_in.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    vector_2d_field_out.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    scalar_3d_field_in.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    scalar_3d_field_out.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    vector_3d_field_in.get_header().get_alloc_properties().request_value_type_allocation<PackType>();
    vector_3d_field_out.get_header().get_alloc_properties().request_value_type_allocation<PackType>();

    // Allocate view
    scalar_2d_field_in.allocate_view();
    scalar_2d_field_out.allocate_view();
    vector_2d_field_in.allocate_view();
    vector_2d_field_out.allocate_view();
    scalar_3d_field_in.allocate_view();
    scalar_3d_field_out.allocate_view();
    vector_3d_field_in.allocate_view();
    vector_3d_field_out.allocate_view();

    // Build the remapper, and register the fields
    std::unique_ptr<Remapper> remapper(new Remapper(phys_grid,dyn_grid));
    remapper->registration_begins();
    remapper->register_field(scalar_2d_field_in, scalar_2d_field_out);
    remapper->register_field(vector_2d_field_in, vector_2d_field_out);
    remapper->register_field(scalar_3d_field_in, scalar_3d_field_out);
    remapper->register_field(vector_3d_field_in, vector_3d_field_out);
    remapper->registration_complete();

    // Generate random numbers
    util::genRandArray(scalar_2d_field_in,  engine, pdf);
    util::genRandArray(vector_2d_field_in,  engine, pdf);
    util::genRandArray(scalar_3d_field_in,  engine, pdf);
    util::genRandArray(vector_3d_field_in,  engine, pdf);

    // Remap
    remapper->remap(/* forward = */ true);

    // Check
    {
      auto phys_in = Kokkos::create_mirror_view(scalar_2d_field_in.template get_reshaped_view<Homme::Real*>());
      auto dyn_out = Kokkos::create_mirror_view(scalar_2d_field_out.template get_reshaped_view<Homme::Real*[NP][NP]>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        if (phys_in(icol)!=dyn_out(ie,ip,jp)) {
            printf("p_in(%d) = %2.16f\n",icol,phys_in(icol));
            printf("d_out(%d,%d,%d) = %2.16f\n",ie,ip,jp,dyn_out(ie,ip,jp));
        }
        REQUIRE (phys_in(icol)==dyn_out(ie,ip,jp));
      }
    }

    {
      auto phys_in = Kokkos::create_mirror_view(vector_2d_field_in.template get_reshaped_view<Homme::Real**>());
      auto dyn_out = Kokkos::create_mirror_view(vector_2d_field_out.template get_reshaped_view<Homme::Real**[NP][NP]>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int icomp=0; icomp<2; ++icomp) {
          if (phys_in(icol,icomp)!=dyn_out(ie,icomp,ip,jp)) {
              printf("p_in(%d,%d) = %2.16f\n",icol,icomp,phys_in(icol,icomp));
              printf("d_out(%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,dyn_out(ie,icomp,ip,jp));
          }
          REQUIRE (phys_in(icol,icomp)==dyn_out(ie,icomp,ip,jp));
        }
      }
    }

    {
      auto phys_in = Kokkos::create_mirror_view(scalar_3d_field_in.template get_reshaped_view<Homme::Real**>());
      auto dyn_out = Kokkos::create_mirror_view(scalar_3d_field_out.template get_reshaped_view<Homme::Real*[NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int ilev=0; ilev<HOMMEXX_NUM_PHYSICAL_LEV; ++ilev) {
          if (phys_in(icol,ilev)!=dyn_out(ie,ip,jp,ilev)) {
              printf("p_in(%d,%d) = %2.16f\n",icol,ilev,phys_in(icol,ilev));
              printf("d_out(%d,%d,%d,%d) = %2.16f\n",ie,ip,jp,ilev,dyn_out(ie,ip,jp,ilev));
          }
          REQUIRE (phys_in(icol,ilev)==dyn_out(ie,ip,jp,ilev));
        }
      }
    }

    {
      auto phys_in = Kokkos::create_mirror_view(vector_3d_field_in.template get_reshaped_view<Homme::Real***>());
      auto dyn_out = Kokkos::create_mirror_view(vector_3d_field_out.template get_reshaped_view<Homme::Real**[NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int icomp=0; icomp<2; ++icomp) {
          for (int ilev=0; ilev<HOMMEXX_NUM_PHYSICAL_LEV; ++ilev) {
            if (phys_in(icol,icomp,ilev)!=dyn_out(ie,icomp,ip,jp,ilev)) {
                printf("p_in(%d,%d,%d) = %2.16f\n",icol,icomp,ilev,phys_in(icol,icomp,ilev));
                printf("d_out(%d,%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,ilev,dyn_out(ie,icomp,ip,jp,ilev));
            }
            REQUIRE (phys_in(icol,icomp,ilev)==dyn_out(ie,icomp,ip,jp,ilev));
          }
        }
      }
    }
  }

  SECTION ("remap bwd") {

    // Create tags and dimensions
    std::vector<FieldTag> scalar_2d_dyn_tags  = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint};
    std::vector<FieldTag> scalar_2d_phys_tags = {FieldTag::Column};
    std::vector<FieldTag> vector_2d_dyn_tags  = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint};
    std::vector<FieldTag> vector_2d_phys_tags = {FieldTag::Column, FieldTag::Component};
    std::vector<FieldTag> scalar_3d_dyn_tags  = {FieldTag::Element, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
    std::vector<FieldTag> scalar_3d_phys_tags = {FieldTag::Column, FieldTag::VerticalLevel};
    std::vector<FieldTag> vector_3d_dyn_tags  = {FieldTag::Element, FieldTag::Component, FieldTag::GaussPoint, FieldTag::GaussPoint, FieldTag::VerticalLevel};
    std::vector<FieldTag> vector_3d_phys_tags = {FieldTag::Column, FieldTag::Component, FieldTag::VerticalLevel};

    std::vector<int> scalar_2d_dyn_dims  = {num_local_elems,    NP, NP};
    std::vector<int> vector_2d_dyn_dims  = {num_local_elems, 2, NP, NP};
    std::vector<int> scalar_3d_dyn_dims  = {num_local_elems,    NP, NP, HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> vector_3d_dyn_dims  = {num_local_elems, 2, NP, NP, HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> scalar_2d_phys_dims = {num_local_columns};
    std::vector<int> vector_2d_phys_dims = {num_local_columns, 2};
    std::vector<int> scalar_3d_phys_dims = {num_local_columns,    HOMMEXX_NUM_PHYSICAL_LEV};
    std::vector<int> vector_3d_phys_dims = {num_local_columns, 2, HOMMEXX_NUM_PHYSICAL_LEV};

    // Create identifiers
    const auto units = units::m;  // Placeholder units
    FieldIdentifier scalar_2d_dyn_fid  ("scalar_2d_dynamics", FieldLayout(scalar_2d_dyn_tags,  scalar_2d_dyn_dims),  units, dyn_grid->name());
    FieldIdentifier scalar_2d_phys_fid ("scalar_2d_physics",  FieldLayout(scalar_2d_phys_tags, scalar_2d_phys_dims), units, phys_grid->name());
    FieldIdentifier vector_2d_dyn_fid  ("vector_2d_dynamics", FieldLayout(vector_2d_dyn_tags,  vector_2d_dyn_dims),  units, dyn_grid->name());
    FieldIdentifier vector_2d_phys_fid ("vector_2d_physics",  FieldLayout(vector_2d_phys_tags, vector_2d_phys_dims), units, phys_grid->name());
    FieldIdentifier scalar_3d_dyn_fid  ("scalar_3d_dynamics", FieldLayout(scalar_3d_dyn_tags,  scalar_3d_dyn_dims),  units, dyn_grid->name());
    FieldIdentifier scalar_3d_phys_fid ("scalar_3d_physics",  FieldLayout(scalar_3d_phys_tags, scalar_3d_phys_dims), units, phys_grid->name());
    FieldIdentifier vector_3d_dyn_fid  ("vector_3d_dynamics", FieldLayout(vector_3d_dyn_tags,  vector_3d_dyn_dims),  units, dyn_grid->name());
    FieldIdentifier vector_3d_phys_fid ("vector_3d_physics",  FieldLayout(vector_3d_phys_tags, vector_3d_phys_dims), units, phys_grid->name());

    // Create fields
    Field<Real,Device> scalar_2d_field_in (scalar_2d_dyn_fid);
    Field<Real,Device> scalar_2d_field_out(scalar_2d_phys_fid);
    Field<Real,Device> vector_2d_field_in (vector_2d_dyn_fid);
    Field<Real,Device> vector_2d_field_out(vector_2d_phys_fid);
    Field<Real,Device> scalar_3d_field_in (scalar_3d_dyn_fid);
    Field<Real,Device> scalar_3d_field_out(scalar_3d_phys_fid);
    Field<Real,Device> vector_3d_field_in (vector_3d_dyn_fid);
    Field<Real,Device> vector_3d_field_out(vector_3d_phys_fid);

    // Request allocation to fit packs of reals
    scalar_2d_field_in.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    scalar_2d_field_out.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    vector_2d_field_in.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    vector_2d_field_out.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    scalar_3d_field_in.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    scalar_3d_field_out.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    vector_3d_field_in.get_header().get_alloc_properties().request_value_type_allocation<Real>();
    vector_3d_field_out.get_header().get_alloc_properties().request_value_type_allocation<Real>();

    // Allocate view
    scalar_2d_field_in.allocate_view();
    scalar_2d_field_out.allocate_view();
    vector_2d_field_in.allocate_view();
    vector_2d_field_out.allocate_view();
    scalar_3d_field_in.allocate_view();
    scalar_3d_field_out.allocate_view();
    vector_3d_field_in.allocate_view();
    vector_3d_field_out.allocate_view();

    // Build the remapper, and register the fields
    std::unique_ptr<Remapper> remapper(new Remapper(phys_grid,dyn_grid));
    remapper->registration_begins();
    remapper->register_field(scalar_2d_field_out, scalar_2d_field_in);
    remapper->register_field(vector_2d_field_out, vector_2d_field_in);
    remapper->register_field(scalar_3d_field_out, scalar_3d_field_in);
    remapper->register_field(vector_3d_field_out, vector_3d_field_in);
    remapper->registration_complete();

    // Generate random numbers
    // Note: for the test to run correctly, the dynamics input vector must be synced,
    //       meaning that the values at the interface between two elements must match.
    //       To do this, we initialize each entry in the dynamic vector with the id
    //       of the corresponding column.
    auto scalar_2d_view = scalar_2d_field_in.get_reshaped_view<Homme::Real*[NP][NP]>();
    auto vector_2d_view = vector_2d_field_in.get_reshaped_view<Homme::Real*[2][NP][NP]>();
    auto scalar_3d_view = scalar_3d_field_in.get_reshaped_view<Homme::Real*[NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
    auto vector_3d_view = vector_3d_field_in.get_reshaped_view<Homme::Real*[2][NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>();
    auto h_scalar_2d_view = Kokkos::create_mirror_view(scalar_2d_view);
    auto h_vector_2d_view = Kokkos::create_mirror_view(vector_2d_view);
    auto h_scalar_3d_view = Kokkos::create_mirror_view(scalar_3d_view);
    auto h_vector_3d_view = Kokkos::create_mirror_view(vector_3d_view);
    for (int ie=0; ie<num_local_elems; ++ie) {
      for (int ip=0; ip<NP; ++ip) {
        for (int jp=0; jp<NP; ++jp) {
          h_scalar_2d_view(ie,ip,jp) = h_d2p(ie*NP*NP + ip*NP + jp,3);
        }
      }
    }
    Kokkos::deep_copy(scalar_2d_view,h_scalar_2d_view);
    Kokkos::deep_copy(vector_2d_view,h_vector_2d_view);
    Kokkos::deep_copy(scalar_3d_view,h_scalar_3d_view);
    Kokkos::deep_copy(vector_3d_view,h_vector_3d_view);

    // Remap
    remapper->remap(/* forward = */ false);

    // Check
    {
      auto dyn_in   = Kokkos::create_mirror_view(scalar_2d_field_in.template get_reshaped_view<Homme::Real*[NP][NP]>());
      auto phys_out = Kokkos::create_mirror_view(scalar_2d_field_out.template get_reshaped_view<Homme::Real*>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        if (phys_out(icol)!=dyn_in(ie,ip,jp)) {
            printf("d_in(%d,%d,%d) = %2.16f\n",ie,ip,jp,dyn_in(ie,ip,jp));
            printf("p_out(%d) = %2.16f\n",icol,phys_out(icol));
        }
        REQUIRE (phys_out(icol)==dyn_in(ie,ip,jp));
      }
    }

    {
      auto dyn_in   = Kokkos::create_mirror_view(vector_2d_field_in.template get_reshaped_view<Homme::Real**[NP][NP]>());
      auto phys_out = Kokkos::create_mirror_view(vector_2d_field_out.template get_reshaped_view<Homme::Real**>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int icomp=0; icomp<2; ++icomp) {
          if (phys_out(icol,icomp)!=dyn_in(ie,icomp,ip,jp)) {
              printf("d_in(%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,dyn_in(ie,icomp,ip,jp));
              printf("p_out(%d,%d) = %2.16f\n",icol,icomp,phys_out(icol,icomp));
          }
          REQUIRE (phys_out(icol,icomp)==dyn_in(ie,icomp,ip,jp));
        }
      }
    }

    {
      auto dyn_in   = Kokkos::create_mirror_view(scalar_3d_field_in.template get_reshaped_view<Homme::Real*[NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>());
      auto phys_out = Kokkos::create_mirror_view(scalar_3d_field_out.template get_reshaped_view<Homme::Real**>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int ilev=0; ilev<2; ++ilev) {
          if (phys_out(icol,ilev)!=dyn_in(ie,ip,jp,ilev)) {
              printf("d_in(%d,%d,%d,%d) = %2.16f\n",ie,ip,jp,ilev,dyn_in(ie,ip,jp,ilev));
              printf("p_out(%d,%d) = %2.16f\n",icol,ilev,phys_out(icol,ilev));
          }
          REQUIRE (phys_out(icol,ilev)==dyn_in(ie,ip,jp,ilev));
        }
      }
    }

    {
      auto dyn_in   = Kokkos::create_mirror_view(vector_3d_field_in.template get_reshaped_view<Homme::Real**[NP][NP][HOMMEXX_NUM_PHYSICAL_LEV]>());
      auto phys_out = Kokkos::create_mirror_view(vector_3d_field_out.template get_reshaped_view<Homme::Real***>());
      for (int icol=0; icol<num_local_columns; ++icol) {
        auto elgp = Kokkos::subview(h_p2d,icol,Kokkos::ALL());
        const int ie = elgp[0];
        const int ip = elgp[1];
        const int jp = elgp[2];
        for (int icomp=0; icomp<2; ++icomp) {
          for (int ilev=0; ilev<2; ++ilev) {
            if (phys_out(icol,icomp,ilev)!=dyn_in(ie,icomp,ip,jp,ilev)) {
                printf("d_in(%d,%d,%d,%d,%d) = %2.16f\n",ie,icomp,ip,jp,ilev,dyn_in(ie,icomp,ip,jp,ilev));
                printf("p_out(%d,%d,%d) = %2.16f\n",icol,icomp,ilev,phys_out(icol,icomp,ilev));
            }
            REQUIRE (phys_out(icol,icomp,ilev)==dyn_in(ie,icomp,ip,jp,ilev));
          }
        }
      }
    }
  }

  // Finalize Homme::MpiContext (deletes buffers manager)
  Homme::MpiContext::finalize_singleton();
}

} // anonymous namespace
