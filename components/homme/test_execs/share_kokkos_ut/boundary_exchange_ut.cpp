#include <catch2/catch.hpp>

#include "mpi/MpiContext.hpp"
#include "mpi/BuffersManager.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/Connectivity.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "Types.hpp"

#include <random>
#include <iomanip>

using namespace Homme;

extern "C" {

void initmp_f90 ();
void init_cube_geometry_f90 (const int& ne);
void init_connectivity_f90 (const int& num_min_max_fields_1d, const int& num_scalar_fields_2d,
                            const int& num_scalar_fields_3d,  const int& num_vector_fields_3d,
                            const int& vector_dim);
void cleanup_f90 ();
void boundary_exchange_test_f90 (F90Ptr& field_min_1d_ptr, F90Ptr& field_max_1d_ptr,
                                 F90Ptr& field_2d_ptr, F90Ptr& field_3d_ptr, F90Ptr& field_4d_ptr,
                                 const int& inner_dim_4d, const int& num_time_levels,
                                 const int& idim_2d, const int& idim_3d, const int& idim_4d,
                                 const int& minmax_split);

} // extern "C"

// =========================== TESTS ============================ //

TEST_CASE ("Boundary Exchange", "Testing the boundary exchange framework")
{
  //std::random_device rd;
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  rngAlg engine(rd());
  std::uniform_real_distribution<Real> dreal(-1.0, 1.0);
  // neighbor_minmax imposes a hard positivity cutoff, so we can't test the bdy
  // exchange alone (for min) with numbers < 0.
  std::uniform_real_distribution<Real> dreal_minmax(0.0, 1.0);
  std::uniform_int_distribution<int>   dint(0,1);

  constexpr int ne        = 2;
  constexpr int num_tests = 1;
  constexpr int DIM       = 2;
  constexpr double test_tolerance = 1e-13;
  constexpr int num_min_max_fields_1d = 1; // Count min and max of a field as 1, does not count the x2 due to min and max
  constexpr int num_scalar_fields_2d  = 1;
  constexpr int num_scalar_fields_3d  = 1;
  constexpr int num_vector_fields_3d  = 1;
  constexpr int field_2d_idim = 0;
  constexpr int field_3d_idim = 1;
  constexpr int field_4d_outer_idim = 2;

  // Initialize f90 mpi stuff
  initmp_f90();

  // Create cube geometry
  init_cube_geometry_f90(ne);

  // Create connectivity
  init_connectivity_f90(num_min_max_fields_1d,num_scalar_fields_2d, num_scalar_fields_3d, num_vector_fields_3d, DIM);
  std::shared_ptr<Connectivity> connectivity = MpiContext::singleton().get_connectivity();

  // Retrieve local number of elements
  int num_elements = connectivity->get_num_local_elements();
  int rank = connectivity->get_comm().rank();

  // Create input data arrays
  HostViewManaged<Real*[num_min_max_fields_1d][NUM_PHYSICAL_LEV]> field_min_1d_f90("", num_elements);
  HostViewManaged<Real*[num_min_max_fields_1d][NUM_PHYSICAL_LEV]> field_max_1d_f90("", num_elements);
  ExecViewManaged<Scalar*[num_min_max_fields_1d][2][NUM_LEV]>     field_1d_cxx("", num_elements);
  auto field_1d_cxx_host = Kokkos::create_mirror_view(field_1d_cxx);

  HostViewManaged<Real*[NUM_TIME_LEVELS][NP][NP]> field_2d_f90("", num_elements);
  ExecViewManaged<Real*[NUM_TIME_LEVELS][NP][NP]> field_2d_cxx("", num_elements);
  ExecViewManaged<Real*[NUM_TIME_LEVELS][NP][NP]>::HostMirror field_2d_cxx_host;
  field_2d_cxx_host = Kokkos::create_mirror_view(field_2d_cxx);

  HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]> field_3d_f90("", num_elements);
  ExecViewManaged<Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]> field_3d_cxx ("", num_elements);
  ExecViewManaged<Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]>::HostMirror field_3d_cxx_host;
  field_3d_cxx_host = Kokkos::create_mirror_view(field_3d_cxx);

  HostViewManaged<Real*[NUM_TIME_LEVELS][DIM][NUM_PHYSICAL_LEV][NP][NP]> field_4d_f90 ("", num_elements);
  ExecViewManaged<Scalar*[NUM_TIME_LEVELS][DIM][NP][NP][NUM_LEV]> field_4d_cxx ("", num_elements);
  ExecViewManaged<Scalar*[NUM_TIME_LEVELS][DIM][NP][NP][NUM_LEV]>::HostMirror field_4d_cxx_host;
  field_4d_cxx_host = Kokkos::create_mirror_view(field_4d_cxx);

  // Get the buffers manager
  std::shared_ptr<BuffersManager> buffers_manager = MpiContext::singleton().get_buffers_manager(MPI_EXCHANGE);
  std::shared_ptr<BuffersManager> buffers_manager_min_max = MpiContext::singleton().get_buffers_manager(MPI_EXCHANGE_MIN_MAX);

  // Create boundary exchanges
  std::shared_ptr<BoundaryExchange> be1 = std::make_shared<BoundaryExchange>(connectivity,buffers_manager);
  std::shared_ptr<BoundaryExchange> be2 = std::make_shared<BoundaryExchange>(connectivity,buffers_manager);
  std::shared_ptr<BoundaryExchange> be3 = std::make_shared<BoundaryExchange>(connectivity,buffers_manager_min_max);

  // Setup the be objects
  be1->set_num_fields(0,num_scalar_fields_2d,DIM*num_vector_fields_3d);
  be1->register_field(field_2d_cxx,1,field_2d_idim);
  be1->register_field(field_4d_cxx,  field_4d_outer_idim,DIM,0);
  be1->registration_completed();

  be2->set_num_fields(0,0,num_scalar_fields_3d);
  be2->register_field(field_3d_cxx,1,field_3d_idim);
  be2->registration_completed();

  be3->set_num_fields(num_min_max_fields_1d,0,0);
  be3->register_min_max_fields(field_1d_cxx,num_min_max_fields_1d,0);
  be3->registration_completed();

  for (int itest=0; itest<num_tests; ++itest)
  {
    // Whether the neighbor min/max should be done as a whole or with two separate calls (start/pack_and_send and finish/recv_and_unpack)
    int minmax_split = dint(engine);

    // Initialize input data to random values
    genRandArray(field_min_1d_f90,engine,dreal_minmax);
    genRandArray(field_max_1d_f90,engine,dreal_minmax);
    for (int ie=0; ie<num_elements; ++ie) {
      for (int ifield=0; ifield<num_min_max_fields_1d; ++ifield) {
        for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
          const int ilev = level / VECTOR_SIZE;
          const int ivec = level % VECTOR_SIZE;
          if (field_min_1d_f90(ie,ifield,level) > field_max_1d_f90(ie,ifield,level)) {
            std::swap(field_min_1d_f90(ie,ifield,level), field_max_1d_f90(ie,ifield,level));
          }
          field_1d_cxx_host(ie,ifield,MIN_ID,ilev)[ivec] = field_min_1d_f90(ie,ifield,level);
          field_1d_cxx_host(ie,ifield,MAX_ID,ilev)[ivec] = field_max_1d_f90(ie,ifield,level);
    }}}
    Kokkos::deep_copy(field_1d_cxx, field_1d_cxx_host);

    genRandArray(field_2d_f90,engine,dreal);
    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++jgp) {
            field_2d_cxx_host(ie,itl,igp,jgp) = field_2d_f90(ie,itl,igp,jgp);
    }}}}
    Kokkos::deep_copy(field_2d_cxx, field_2d_cxx_host);

    genRandArray(field_3d_f90,engine,dreal);
    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
          const int ilev = level / VECTOR_SIZE;
          const int ivec = level % VECTOR_SIZE;
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              field_3d_cxx_host(ie,itl,igp,jgp,ilev)[ivec] = field_3d_f90(ie,itl,level,igp,jgp);
    }}}}}
    Kokkos::deep_copy(field_3d_cxx, field_3d_cxx_host);

    genRandArray(field_4d_f90,engine,dreal);
    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int idim=0; idim<DIM; ++idim) {
          for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
            const int ilev = level / VECTOR_SIZE;
            const int ivec = level % VECTOR_SIZE;
            for (int igp=0; igp<NP; ++igp) {
              for (int jgp=0; jgp<NP; ++jgp) {
                field_4d_cxx_host(ie,itl,idim,igp,jgp,ilev)[ivec] = field_4d_f90(ie,itl,idim,level,igp,jgp);
    }}}}}}
    Kokkos::deep_copy(field_4d_cxx, field_4d_cxx_host);

    // Perform boundary exchange
    boundary_exchange_test_f90(field_min_1d_f90.data(), field_max_1d_f90.data(),
                               field_2d_f90.data(), field_3d_f90.data(), field_4d_f90.data(),
                               DIM, NUM_TIME_LEVELS, field_2d_idim+1, field_3d_idim+1, field_4d_outer_idim+1, minmax_split);
    minmax_split = 1;
    if (minmax_split==0) {
      be1->exchange();
      be2->exchange();
      be3->exchange_min_max();
    } else {
      be3->pack_and_send_min_max();
      be1->pack_and_send();
      be1->recv_and_unpack();
      be2->pack_and_send();
      be2->recv_and_unpack();
      be3->recv_and_unpack_min_max();
    }
    Kokkos::deep_copy(field_1d_cxx_host, field_1d_cxx);
    Kokkos::deep_copy(field_2d_cxx_host,     field_2d_cxx);
    Kokkos::deep_copy(field_3d_cxx_host,     field_3d_cxx);
    Kokkos::deep_copy(field_4d_cxx_host,     field_4d_cxx);

    // Compare answers
    for (int ie=0; ie<num_elements; ++ie) {
      for (int ifield=0; ifield<num_min_max_fields_1d; ++ifield) {
        for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
          const int ilev = level / VECTOR_SIZE;
          const int ivec = level % VECTOR_SIZE;
          REQUIRE(compare_answers(field_min_1d_f90(ie,ifield,level),field_1d_cxx_host(ie,ifield,MIN_ID,ilev)[ivec]) < test_tolerance);
          if(compare_answers(field_min_1d_f90(ie,ifield,level),field_1d_cxx_host(ie,ifield,MIN_ID,ilev)[ivec]) >= test_tolerance) {
            std::cout << std::setprecision(17) << "rank,ie,ifield,ilev,iv: " << rank << ", " << ie << ", " << ifield << ", " << ilev << ", " << ivec << "\n";
            std::cout << std::setprecision(17) << "f90: " << field_min_1d_f90(ie,ifield,level) << "\n";
            std::cout << std::setprecision(17) << "cxx: " << field_1d_cxx_host(ie,ifield,MIN_ID,ilev)[ivec] << "\n";
          }
          REQUIRE(compare_answers(field_max_1d_f90(ie,ifield,level),field_1d_cxx_host(ie,ifield,MAX_ID,ilev)[ivec]) < test_tolerance);
          if(compare_answers(field_max_1d_f90(ie,ifield,level),field_1d_cxx_host(ie,ifield,MAX_ID,ilev)[ivec]) >= test_tolerance) {
            std::cout << std::setprecision(17) << "rank,ie,ifield,ilev,iv: " << rank << ", " << ie << ", " << ifield << ", " << ilev << ", " << ivec << "\n";
            std::cout << std::setprecision(17) << "f90: " << field_max_1d_f90(ie,ifield,level) << "\n";
            std::cout << std::setprecision(17) << "cxx: " << field_1d_cxx_host(ie,ifield,MAX_ID,ilev)[ivec] << "\n";
          }
    }}}

    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++jgp) {
            if(compare_answers(field_2d_f90(ie,itl,igp,jgp),field_2d_cxx_host(ie,itl,igp,jgp)) >= test_tolerance) {
              std::cout << "rank,ie,itl,igp,jgp: " << rank << ", " << ie << ", " << itl << ", " << igp << ", " << jgp << "\n";
              std::cout << "f90: " << field_2d_f90(ie,itl,igp,jgp) << "\n";
              std::cout << "cxx: " << field_2d_cxx_host(ie,itl,igp,jgp) << "\n";
            }
            REQUIRE(compare_answers(field_2d_f90(ie,itl,igp,jgp),field_2d_cxx_host(ie,itl,igp,jgp)) < test_tolerance);
    }}}}

    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
          const int ilev = level / VECTOR_SIZE;
          const int ivec = level % VECTOR_SIZE;
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              if(compare_answers(field_3d_f90(ie,itl,level,igp,jgp),field_3d_cxx_host(ie,itl,igp,jgp,ilev)[ivec]) >= test_tolerance) {
                std::cout << std::setprecision(17) << "rank,ie,itl,igp,jgp,ilev,iv: " << rank << ", " << ie << ", " << itl << ", " << igp << ", " << jgp << ", " << ilev << ", " << ivec << "\n";
                std::cout << std::setprecision(17) << "f90: " << field_3d_f90(ie,itl,level,igp,jgp) << "\n";
                std::cout << std::setprecision(17) << "cxx: " << field_3d_cxx_host(ie,itl,igp,jgp,ilev)[ivec] << "\n";
              }
              REQUIRE(compare_answers(field_3d_f90(ie,itl,level,igp,jgp),field_3d_cxx_host(ie,itl,igp,jgp,ilev)[ivec]) < test_tolerance);
    }}}}}

    for (int ie=0; ie<num_elements; ++ie) {
      for (int itl=0; itl<NUM_TIME_LEVELS; ++itl) {
        for (int idim=0; idim<DIM; ++idim) {
          for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
            const int ilev = level / VECTOR_SIZE;
            const int ivec = level % VECTOR_SIZE;
            for (int igp=0; igp<NP; ++igp) {
              for (int jgp=0; jgp<NP; ++jgp) {
                if(compare_answers(field_4d_f90(ie,itl,idim,level,igp,jgp),field_4d_cxx_host(ie,itl,idim,igp,jgp,ilev)[ivec]) >= test_tolerance) {
                  std::cout << std::setprecision(17) << "rank,ie,itl,idim,igp,jgp,ilev,iv: " << rank << ", " << ie << ", " << itl << ", " << idim << ", " << igp << ", " << jgp << ", " << ilev << ", " << ivec << "\n";
                  std::cout << std::setprecision(17) << "f90: " << field_4d_f90(ie,itl,idim,level,igp,jgp) << "\n";
                  std::cout << std::setprecision(17) << "cxx: " << field_4d_cxx_host(ie,itl,idim,igp,jgp,ilev)[ivec] << "\n";
                }
                REQUIRE(compare_answers(field_4d_f90(ie,itl,idim,level,igp,jgp),field_4d_cxx_host(ie,itl,idim,igp,jgp,ilev)[ivec]) < test_tolerance);
    }}}}}}
  }

  // Cleanup
  cleanup_f90();  // Deallocate stuff in the F90 module
  be1->clean_up();
  be2->clean_up();
  be3->clean_up();
}
