#include <catch2/catch.hpp>

#include <limits>

#include "Elements.hpp"
#include "Dimensions.hpp"
#include "KernelVariables.hpp"
#include "SphereOperators.hpp"
#include "Types.hpp"
#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"

#include <assert.h>
#include <stdio.h>
#include <random>

using namespace Homme;

using rngAlg = std::mt19937_64;

extern "C" {

void laplace_simple_c_callable(const Real *input,
                               const Real *dvv,
                               const Real *dinv,
                               const Real *metdet,
                               Real *output);

void gradient_sphere_c_callable(const Real *input,
                                const Real *dvv,
                                const Real *dinv,
                                Real *output);

void divergence_sphere_wk_c_callable(const Real *input,
                                     const Real *dvv,
                                     const Real *spheremp,
                                     const Real *dinv,
                                     Real *output);

}  // extern C

class compute_sphere_operator_test {
 public:
  compute_sphere_operator_test(int num_elems)
      : _num_elems(num_elems),
        scalar_input_d("scalar input", num_elems),
        vector_input_d("vector input", num_elems),
        scalar_output_d("scalar output", num_elems),
        vector_output_d("vector output", num_elems),
        scalar_input_host(
            Kokkos::create_mirror_view(scalar_input_d)),
        vector_input_host(
            Kokkos::create_mirror_view(vector_input_d)),
        scalar_output_host(
            Kokkos::create_mirror_view(scalar_output_d)),
        vector_output_host(
            Kokkos::create_mirror_view(vector_output_d))
  {
    // constructor's body
    // init randonly

    std::random_device rd;
    rngAlg engine(rd());

    // check singularities? divergence_wk uses both D and
    // Dinv, does it matter if they are not inverses of
    // each other?
    genRandArray(scalar_input_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 100.0));
    Kokkos::deep_copy(scalar_input_d, scalar_input_host);

    genRandArray(vector_input_host, engine,
                 std::uniform_real_distribution<Real>(
                     -100.0, 100.0));
    Kokkos::deep_copy(vector_input_d, vector_input_host);

    // D
    ExecViewManaged<Real * [2][2][NP][NP]> d_d("",num_elems);
    d_host = Kokkos::create_mirror_view(d_d);
    genRandArray(d_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 1.0));
    Kokkos::deep_copy(d_d, d_host);

    // Dinv
    ExecViewManaged<Real * [2][2][NP][NP]> dinv_d("",num_elems);
    dinv_host = Kokkos::create_mirror_view(dinv_d);
    genRandArray(dinv_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 1.0));
    Kokkos::deep_copy(dinv_d, dinv_host);

    // metdet
    ExecViewManaged<Real * [NP][NP]> metdet_d("",num_elems);
    metdet_host = Kokkos::create_mirror_view(metdet_d);
    genRandArray(metdet_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 1.0));
    Kokkos::deep_copy(metdet_d, metdet_host);

    // spheremp
    ExecViewManaged<Real * [NP][NP]> spheremp_d ("",num_elems);
    spheremp_host = Kokkos::create_mirror_view(spheremp_d);
    genRandArray(spheremp_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 1.0));
    Kokkos::deep_copy(spheremp_d, spheremp_host);

    // dvv
    ExecViewManaged<Real[NP][NP]> dvv_d("");
    dvv_host = Kokkos::create_mirror_view(dvv_d);
    genRandArray(dvv_host, engine,
                 std::uniform_real_distribution<Real>(
                     0, 1.0));
    Kokkos::deep_copy(dvv_d, dvv_host);

    // Set device views in SphereOperators
    ExecViewManaged<Real*      [NP][NP]> mp_d("",num_elems);  // Unused by this test, but needed by sphere_ops
    ExecViewManaged<Real*[2][2][NP][NP]> metinv_d("",num_elems);  // Unused by this test, but needed by sphere_ops
    sphere_ops.set_views(dvv_d,d_d,dinv_d,metinv_d,metdet_d,spheremp_d,mp_d);
  }

  const int _num_elems;  // league size, serves as ie index

  // device views
  ExecViewManaged<Real * [NP][NP]> scalar_input_d;
  ExecViewManaged<Real * [2][NP][NP]> vector_input_d;
  ExecViewManaged<Real * [NP][NP]> scalar_output_d;
  ExecViewManaged<Real * [2][NP][NP]> vector_output_d;

  // host views, one dim is num_elems. Spherical operators
  // do not take ie or nlev fields,  but to make it a more
  // reasonable test and to have parallel_for we ise another
  // dimension.
  ExecViewManaged<Real * [NP][NP]>::HostMirror
      scalar_input_host;
  // how to get total length of view? use dim0*dim1*...till
  // dim7
  const int scalar_input_len = NP * NP;  // temp code
  ExecViewManaged<Real * [2][NP][NP]>::HostMirror
      vector_input_host;
  // how to get total length of view? use dim0*dim1*...till
  // dim7
  const int vector_input_len = 2 * NP * NP;  // temp code
  ExecViewManaged<Real * [2][2][NP][NP]>::HostMirror d_host;
  const int d_len = 2 * 2 * NP * NP;  // temp code
  ExecViewManaged<Real * [2][2][NP][NP]>::HostMirror
      dinv_host;
  const int dinv_len = 2 * 2 * NP * NP;  // temp code
  ExecViewManaged<Real * [NP][NP]>::HostMirror metdet_host;
  const int metdet_len = NP * NP;
  ExecViewManaged<Real * [NP][NP]>::HostMirror
      spheremp_host;
  const int spheremp_len = NP * NP;
  ExecViewManaged<Real[NP][NP]>::HostMirror dvv_host;
  const int dvv_len = NP * NP;
  ExecViewManaged<Real * [NP][NP]>::HostMirror
      scalar_output_host;
  ExecViewManaged<Real * [2][NP][NP]>::HostMirror
      vector_output_host;

  SphereOperators     sphere_ops;

  // tag for laplace_simple()
  struct TagSimpleLaplace {};
  // tag for gradient_sphere()
  struct TagGradientSphere {};
  // tag for divergence_sphere_wk
  struct TagDivergenceSphereWk {};
  // tag for default, a dummy
  struct TagDefault {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagDefault &,
                  TeamMember team) const {
      // do nothing or print a message
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagSimpleLaplace &,
                  TeamMember team) const {
    KernelVariables kv(team);

    sphere_ops.laplace_wk_sl(kv,
                  Homme::subview(scalar_input_d,kv.ie),
                  Homme::subview(scalar_output_d,kv.ie));
  };  // end of op() for laplace_simple

  /*
   * A comment on how these tests work:
   * Consider 160 threads available, with _num_elems=10;
   *
   * Below are lines
   *     Kokkos::TeamPolicy<ExecSpace, TagSimpleLaplace>
   * policy(_num_elems, 16); Kokkos::parallel_for(policy,
   * *this); this one will call operator() with, say,
   * weak divergence tag. Then first 160 threads will
   * be clustered into 10 leagues (as many as _num_elems).
   * Each league will contain 16 threads (16 as the second
   * argument in policy().) Each league will have its
   * separate subview input and subview output (subview of
   * global arrays based on team.league_rank), so, league 2
   * will have input from local_vector_input_d(2,:,:), etc.
   * When divergence_sphere_wk is called, it will be
   * executed by 16 threads in league, each sharing input
   * and output. So, it is not a perfect situation and not a
   * perfect test, because 16 team threads run the same code
   * and OVERWRITE the same output. A better test should
   * have another level of parallelism, a loop with
   * TeamThreadRange. Also, one should note that
   * divergence_sphere_wk as well as other SphereOperators
   * should be called from loop with aeamThreadRange.
   */

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagDivergenceSphereWk &,
                  TeamMember team) const {
    KernelVariables kv(team);

    sphere_ops.divergence_sphere_wk_sl(team,
                            Homme::subview(vector_input_d, kv.ie),
                            Homme::subview(scalar_output_d,kv.ie));
  };  // end of op() for divergence_sphere_wk

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagGradientSphere &,
                  TeamMember team) const {
    KernelVariables kv(team);

    sphere_ops.gradient_sphere_sl(team,
                       Homme::subview(scalar_input_d, kv.ie),
                       Homme::subview(vector_output_d,kv.ie));
  };

  // this could be even nicer,
  // put in a param in run_functor(param) to only branch
  // policy type
  void run_functor_simple_laplace() {
    auto policy = Homme::get_default_team_policy<ExecSpace,TagSimpleLaplace>(_num_elems);
    sphere_ops.allocate_buffers(policy);
    Kokkos::parallel_for(policy, *this);
    ExecSpace::fence();
    // TO FROM
    Kokkos::deep_copy(scalar_output_host, scalar_output_d);
  };

  void run_functor_gradient_sphere() {
    auto policy = Homme::get_default_team_policy<ExecSpace,TagGradientSphere>(_num_elems);
    sphere_ops.allocate_buffers(policy);
    Kokkos::parallel_for(policy, *this);
    ExecSpace::fence();
    // TO FROM
    Kokkos::deep_copy(vector_output_host, vector_output_d);
  };

  void run_functor_div_wk() {
    auto policy = Homme::get_default_team_policy<ExecSpace,TagDivergenceSphereWk>(_num_elems);
    sphere_ops.allocate_buffers(policy);
    Kokkos::parallel_for(policy, *this);
    ExecSpace::fence();
    // TO FROM
    // remember to copy correct output
    Kokkos::deep_copy(scalar_output_host, scalar_output_d);
  };

};  // end of definition of compute_sphere_operator_test()

TEST_CASE("testing_laplace_simple_sl",
          "laplace_simple_sl") {
  constexpr const int elements = 10;

  compute_sphere_operator_test testing_laplace(elements);

  testing_laplace.run_functor_simple_laplace();

  HostViewManaged<Real[NP][NP]> local_fortran_output("fortran results");

  for(int ie = 0; ie < elements; ie++) {
    HostViewUnmanaged<Real[NP][NP]> local_scalar_input =
        Homme::subview(testing_laplace.scalar_input_host,
                        ie);

    HostViewUnmanaged<Real[2][2][NP][NP]> local_dinv =
        Homme::subview(testing_laplace.dinv_host, ie);

    HostViewUnmanaged<Real[NP][NP]> local_spheremp =
        Homme::subview(testing_laplace.spheremp_host,
                        ie);

    // run F code
    laplace_simple_c_callable(
        local_scalar_input.data(),
        testing_laplace.dvv_host.data(),
        local_dinv.data(), local_spheremp.data(),
        local_fortran_output.data());

    // compare answers
    for(int igp = 0; igp < NP; ++igp) {
      for(int jgp = 0; jgp < NP; ++jgp) {
        REQUIRE(
            !std::isnan(local_fortran_output(igp, jgp)));
        REQUIRE(
            !std::isnan(testing_laplace.scalar_output_host(
                ie, igp, jgp)));
        Real rel_error = compare_answers(
                    local_fortran_output(igp, jgp),
                    testing_laplace.scalar_output_host(
                        ie, igp, jgp));
        REQUIRE(local_fortran_output(igp, jgp) ==
                    testing_laplace.scalar_output_host(
                        ie, igp, jgp));
      }  // jgp
    }    // igp
  }      // end of for loop for elements

  std::cout
      << "simple_laplace_sl single level test finished.\n";

};  // end of TEST_CASE(..., "simple laplace")

TEST_CASE("Testing div_wk_sl()", "div_wk_sl") {
  constexpr const int elements = 1;

  compute_sphere_operator_test testing_divwk(elements);

  testing_divwk.run_functor_div_wk();

  for(int ie = 0; ie < elements; ie++) {
    Real local_fortran_output[NP][NP];

    HostViewUnmanaged<Real[2][NP][NP]> local_vector_input =
        Homme::subview(testing_divwk.vector_input_host,
                        ie);

    HostViewUnmanaged<Real[2][2][NP][NP]> local_dinv =
        Homme::subview(testing_divwk.dinv_host, ie);

    HostViewUnmanaged<Real[NP][NP]> local_spheremp =
        Homme::subview(testing_divwk.spheremp_host, ie);

    Real vf[2][NP][NP];
    Real dvvf[NP][NP];
    Real dinvf[2][2][NP][NP];
    Real sphf[NP][NP];

    for(int _i = 0; _i < NP; _i++)
      for(int _j = 0; _j < NP; _j++) {
        dvvf[_i][_j] = testing_divwk.dvv_host(_i, _j);
        sphf[_i][_j] = local_spheremp(_i, _j);
        for(int _d1 = 0; _d1 < 2; _d1++) {
          vf[_d1][_i][_j] = local_vector_input(_d1, _i, _j);
          for(int _d2 = 0; _d2 < 2; _d2++)
            dinvf[_d1][_d2][_i][_j] =
                local_dinv(_d1, _d2, _i, _j);
        }
      }

    divergence_sphere_wk_c_callable(
        &(vf[0][0][0]), &(dvvf[0][0]), &(sphf[0][0]),
        &(dinvf[0][0][0][0]),
        &(local_fortran_output[0][0]));

    for(int igp = 0; igp < NP; ++igp) {
      for(int jgp = 0; jgp < NP; ++jgp) {
        REQUIRE(
            !std::isnan(local_fortran_output[igp][jgp]));
        REQUIRE(
            !std::isnan(testing_divwk.scalar_output_host(
                ie, igp, jgp)));
        REQUIRE(local_fortran_output[igp][jgp] ==
                    testing_divwk.scalar_output_host(
                        ie, igp, jgp));
      }  // jgp
    }    // igp
  };     // end of elements loop

  std::cout << "div_wk single level test finished.\n";

}  // end of TEST_CASE(...,"divergence_sphere_wk")

TEST_CASE("Testing gradient_sphere_sl()",
          "gradient_sphere") {
  constexpr const int elements = 10;

  compute_sphere_operator_test testing_grad(elements);

  // running kokkos version of operator
  testing_grad.run_functor_gradient_sphere();

  for(int ie = 0; ie < elements; ie++) {
    Real local_fortran_output[2][NP][NP];

    HostViewUnmanaged<Real[NP][NP]> local_scalar_input =
        Homme::subview(testing_grad.scalar_input_host,
                        ie);

    HostViewUnmanaged<Real[2][2][NP][NP]> local_dinv =
        Homme::subview(testing_grad.dinv_host, ie);

    Real sf[NP][NP];
    Real dvvf[NP][NP];
    Real dinvf[2][2][NP][NP];

    // flipping arrays -- WRITE DOWN HOW THIS SHOULD BE DONE
    for(int _i = 0; _i < NP; _i++)
      for(int _j = 0; _j < NP; _j++) {
        sf[_i][_j] = local_scalar_input(_i, _j);
        dvvf[_i][_j] = testing_grad.dvv_host(_i, _j);
        for(int _d1 = 0; _d1 < 2; _d1++)
          for(int _d2 = 0; _d2 < 2; _d2++)
            dinvf[_d1][_d2][_i][_j] =
                local_dinv(_d1, _d2, _i, _j);
      }

    // running F version of operator
    gradient_sphere_c_callable(
        &(sf[0][0]), &(dvvf[0][0]), &(dinvf[0][0][0][0]),
        &(local_fortran_output[0][0][0]));

    // comparing answers from kokkos and F
    for(int igp = 0; igp < NP; ++igp)
      for(int jgp = 0; jgp < NP; ++jgp)
        for(int _d = 0; _d < 2; ++_d) {
          REQUIRE(!std::isnan(
              local_fortran_output[_d][igp][jgp]));
          REQUIRE(
              !std::isnan(testing_grad.vector_output_host(
                  ie, _d, igp, jgp)));
          REQUIRE(local_fortran_output[_d][igp][jgp] ==
                      testing_grad.vector_output_host(
                          ie, _d, igp, jgp));
        }  // end of comparing answers

  }  // end of loop for elements

  std::cout << "grad single level test finished.\n";

};  // end of TEST_CASE(..., "gradient_sphere")

// SHMEM ????
