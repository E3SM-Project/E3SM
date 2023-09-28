#include <catch2/catch.hpp>

#include <limits>

#include "Dimensions.hpp"
#include "KernelVariables.hpp"
#include "SphereOperators.hpp"
#include "Types.hpp"
#include "PhysicalConstants.hpp"
#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"

#include <assert.h>
#include <stdio.h>
#include <random>

using namespace Homme;

using rngAlg = std::mt19937_64;

class compute_sphere_operator_test_ml {
 public:

  compute_sphere_operator_test_ml(int num_elems)
      : _num_elems(num_elems),
        scalar_input_d("scalar input", num_elems),
        scalar_input_COPY2_d("scalar input 2", num_elems),
        vector_input_d("vector input", num_elems),
        scalar_output_d("scalar output", num_elems),
        vector_output_d("vector output", num_elems),

        // Make certain Kokkos doesn't use the same arrays
        // as used on the device when they are mutable
        scalar_input_host("scalar input host", num_elems),
        scalar_input_COPY2_host("scalar input host copy", num_elems),
        vector_input_host("vector input host", num_elems),
        scalar_output_host(
            Kokkos::create_mirror_view(scalar_output_d)),
        vector_output_host(
            Kokkos::create_mirror_view(vector_output_d)),
        sphere_ops(PhysicalConstants::rearth0, 1/PhysicalConstants::rearth0)
  {
    std::random_device rd;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);

    // Although scalar_output_d is an output, it can be used as input by divergence_sphere_update
    genRandArray(
        scalar_output_host, engine,
        std::uniform_real_distribution<Real>(-1000.0,
                                             1000.0));
    Kokkos::deep_copy(scalar_output_host, scalar_output_d);

    genRandArray(
        scalar_input_host, engine,
        std::uniform_real_distribution<Real>(-1000.0,
                                             1000.0));
    Kokkos::deep_copy(scalar_input_d, scalar_input_host);
    Kokkos::deep_copy(scalar_input_COPY2_d, scalar_input_host);
    Kokkos::deep_copy(scalar_input_COPY2_host, scalar_input_host);

    genRandArray(
        vector_input_host, engine,
        std::uniform_real_distribution<Real>(-1000.0,
                                             1000.0));
    Kokkos::deep_copy(vector_input_d, vector_input_host);

    // D
    ExecViewManaged<Real * [2][2][NP][NP]> d_d("",num_elems);
    d_host = Kokkos::create_mirror_view(d_d);
    genRandArray(d_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(d_d, d_host);

    // Dinv
    ExecViewManaged<Real * [2][2][NP][NP]> dinv_d("",num_elems);
    dinv_host = Kokkos::create_mirror_view(dinv_d);
    genRandArray(dinv_host,
                 engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(dinv_d, dinv_host);

    // metinv
    ExecViewManaged<Real * [2][2][NP][NP]> metinv_d("",num_elems);
    metinv_host = Kokkos::create_mirror_view(metinv_d);
    genRandArray(metinv_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(metinv_d, metinv_host);

    // metdet
    ExecViewManaged<Real * [NP][NP]> metdet_d("",num_elems);
    metdet_host = Kokkos::create_mirror_view(metdet_d);
    genRandArray(metdet_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(metdet_d, metdet_host);

    // spheremp
    ExecViewManaged<Real * [NP][NP]> spheremp_d("",num_elems);
    spheremp_host = Kokkos::create_mirror_view(spheremp_d);
    genRandArray(spheremp_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(spheremp_d, spheremp_host);

    // mp
    ExecViewManaged<Real [NP][NP]> mp_d("");
    mp_host = Kokkos::create_mirror_view(mp_d);
    genRandArray(mp_host,
                 engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(mp_d, mp_host);

    // dvv
    ExecViewManaged<Real[NP][NP]> dvv_d("");
    dvv_host = Kokkos::create_mirror_view(dvv_d);
    genRandArray(dvv_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(dvv_d, dvv_host);

    genRandArray(&alpha, 1, engine,
                 std::uniform_real_distribution<Real>(
                     -30.0, 30.0));

    sphere_ops.set_views(dvv_d,d_d,dinv_d,metinv_d,metdet_d,spheremp_d,mp_d);
  }  // end of constructor

  int _num_elems;  // league size, serves as ie index

  // device
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>
      scalar_input_d;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>
      scalar_input_COPY2_d;
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>
      vector_input_d;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>
      scalar_output_d;
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>
      vector_output_d;
  // host
  // rely on fact NUM_PHYSICAL_LEV=NUM_LEV*VECTOR_SIZE
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>::HostMirror
      scalar_input_host,
      scalar_input_COPY2_host;
  const int scalar_input_len =
      NUM_PHYSICAL_LEV * NP * NP;  // temp code

  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>::HostMirror
      vector_input_host;
  const int vector_input_len =
      NUM_PHYSICAL_LEV * 2 * NP * NP;

  ExecViewManaged<Real * [2][2][NP][NP]>::HostMirror d_host;
  const int d_len = 2 * 2 * NP * NP;  // temp code

  ExecViewManaged<Real * [2][2][NP][NP]>::HostMirror
      dinv_host;
  const int dinv_len = 2 * 2 * NP * NP;  // temp code

  ExecViewManaged<Real * [2][2][NP][NP]>::HostMirror
      metinv_host;
  const int metinv_len = 2 * 2 * NP * NP;  // temp code

  ExecViewManaged<Real * [NP][NP]>::HostMirror metdet_host;
  const int metdet_len = NP * NP;

  ExecViewManaged<Real * [NP][NP]>::HostMirror
      spheremp_host;
  const int spheremp_len = NP * NP;

  ExecViewManaged<Real [NP][NP]>::HostMirror mp_host;
  const int mp_len = NP * NP;

  ExecViewManaged<Real[NP][NP]>::HostMirror dvv_host;
  const int dvv_len = NP * NP;

  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>::HostMirror
      scalar_output_host;
  ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]>::HostMirror
      vector_output_host;

  SphereOperators     sphere_ops;

  Real alpha;

  // tag for divergence_sphere_wk
  struct TagDivergenceSphereWkML {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagDivergenceSphereWkML &,
                  const TeamMember& team) const {
    KernelVariables kv(team);

    sphere_ops.divergence_sphere_wk(kv,
                         Homme::subview(vector_input_d, kv.ie),
                         Homme::subview(scalar_output_d,kv.ie));
  }  // end of op() for divergence_sphere_wk_ml

  void run_functor_divergence_sphere_wk() {
    auto policy = Homme::get_default_team_policy<ExecSpace, TagDivergenceSphereWkML>(_num_elems);
    sphere_ops.allocate_buffers(policy);

    //for (int i = 0; i < 1000000; i++)   
    Kokkos::parallel_for(policy, *this);

    Kokkos::fence();
    Kokkos::deep_copy(scalar_output_host, scalar_output_d);
  };

};  // end of class def compute_sphere_op_test_ml


TEST_CASE("divergence_sphere_wk",
          "divergence_sphere_wk") {
  constexpr const int elements = 1000000;

  compute_sphere_operator_test_ml testing_div_ml(elements);
  testing_div_ml.run_functor_divergence_sphere_wk();

  std::cout << "test div_wk multilevel finished. \n";

}  // end of test div_sphere_wk_ml

