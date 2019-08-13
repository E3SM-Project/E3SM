/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "HybridVCoord.hpp"
#include "ErrorDefs.hpp"
#include "utilities/TestUtils.hpp"

#include <random>

namespace Homme
{

void HybridVCoord::init(const Real ps0_in,
                   CRCPtr hybrid_am_ptr,
                   CRCPtr hybrid_ai_ptr,
                   CRCPtr hybrid_bm_ptr,
                   CRCPtr hybrid_bi_ptr)
{
  ps0 = ps0_in;

  // hybrid_am = ExecViewManaged<Real[NUM_PHYSICAL_LEV]>(
  //    "Hybrid coordinates; coefficient A_midpoints");
  hybrid_ai = ExecViewManaged<Real[NUM_INTERFACE_LEV]>(
      "Hybrid coordinates; coefficient A_interfaces");
  hybrid_bi = ExecViewManaged<Real[NUM_INTERFACE_LEV]>(
      "Hybrid coordinates; coefficient B_interfaces");
  // hybrid_bm = ExecViewManaged<Real[NUM_PHYSICAL_LEV]>(
  //    "Hybrid coordinates; coefficient B_midpoints");

  // HostViewUnmanaged<const Real[NUM_PHYSICAL_LEV]>
  // host_hybrid_am(hybrid_am_ptr);
  // HostViewUnmanaged<const Real[NUM_PHYSICAL_LEV]>
  // host_hybrid_bm(hybrid_bm_ptr);
  HostViewUnmanaged<const Real[NUM_INTERFACE_LEV]> host_hybrid_ai(
      hybrid_ai_ptr);
  Kokkos::deep_copy(hybrid_ai, host_hybrid_ai);
  HostViewUnmanaged<const Real[NUM_INTERFACE_LEV]> host_hybrid_bi(
      hybrid_bi_ptr);
  Kokkos::deep_copy(hybrid_bi, host_hybrid_bi);

  // i don't think this saves us much now
  {
    // Only hybrid_ai(0) is needed.
    hybrid_ai0 = hybrid_ai_ptr[0];
  }
  // this is not in master anymore?
  //  assert(hybrid_ai_ptr != nullptr);
  //  assert(hybrid_bi_ptr != nullptr);

  compute_deltas();
}

void HybridVCoord::random_init(int seed) {
  hybrid_ai = ExecViewManaged<Real[NUM_INTERFACE_LEV]>(
      "Hybrid a_interface coefs");
  hybrid_bi = ExecViewManaged<Real[NUM_INTERFACE_LEV]>(
      "Hybrid b_interface coefs");

  std::mt19937_64 engine(seed);
  ps0 = 1.0;

  using urd = std::uniform_real_distribution<Real>;

  // eta=hybrid_a+hybrid_b increasing ensures p is monotonic (sufficient, not necessary cond).
  // Also, hybrid_a and hybrid_b must be between 0 and 1.
  // So easy way: generate uniform eta=a+b in [eta_min,1], with eta_min>0.
  // That means eta(i) = i*delta_eta, with delta_eta = (1-eta_min)/NUM_PHYSICAL_LEV.
  // Then pick a = c*b, with c<1 (say c=0.1). Then eta = b+c*b=(c+1)b,
  // so b(i) = i*delta_eta/(c+1) and a(i) = c*i*delta_eta/(c+1).
  // To add randomness, perturb sligthly: a(i)=a(i)+da and b(i)=b(i)+db,
  // without violating eta(i+1)>eta(i). For instance, let da+db<delta_eta/2,
  // which is true if da<delta_eta/4 and db<delta_eta/4.

  decltype(hybrid_ai)::HostMirror host_hybrid_ai("Host hybrid ai coefs");
  decltype(hybrid_bi)::HostMirror host_hybrid_bi("Host hybrid bi coefs");
  const Real eta_min = 0.001;

  // The proportionality constant between a and b.
  Real c;
  genRandArray(&c,1,engine,urd(0.0001,0.25));

  // The eta uniform grid spacing
  Real deta = (1-eta_min)/NUM_PHYSICAL_LEV;
  Real perturbation = deta/5;

  for (int i=0; i<NUM_INTERFACE_LEV; ++i) {
    // Compute eta
    Real eta_i = eta_min+deta*i;

    // Create b and a
    host_hybrid_bi(i) = eta_i/(c+1);
    host_hybrid_ai(i) = c*host_hybrid_bi(i);

    // Add some randomness to a and b
    Real dab[2];
    Real perturbation_adjusted = std::min(std::min(perturbation,eta_i),std::min(perturbation,1-eta_i));;

    auto perturbation_dist = urd(-perturbation_adjusted,perturbation_adjusted);
    genRandArray(dab,2,engine,perturbation_dist);
    host_hybrid_ai(i) += dab[0];
    host_hybrid_bi(i) += dab[1];

    // Enforce 0<a,b<1
    if (host_hybrid_ai(i)<0) {
      host_hybrid_ai(i) = 0;
    }
    if (host_hybrid_bi(i)<0) {
      host_hybrid_bi(i) = 0;
    }
    if (host_hybrid_ai(i)>1) {
      host_hybrid_ai(i) = 1;
    }
    if (host_hybrid_bi(i)>1) {
      host_hybrid_bi(i) = 1;
    }
  }

  // Enforce b(0)=0, a(0)=eta_min, and b(end)=1, a(end)=0;
  host_hybrid_ai(0) = eta_min;
  host_hybrid_bi(0) = 0;
  host_hybrid_ai(NUM_PHYSICAL_LEV) = 0;
  host_hybrid_bi(NUM_PHYSICAL_LEV) = 1;
  
  // Safety check: a+b should be monotone here
  for (int i=1; i<NUM_INTERFACE_LEV; ++i) {
    Real curr = host_hybrid_ai(i) + host_hybrid_bi(i);
    Real prev = host_hybrid_ai(i-1) + host_hybrid_bi(i-1);   

    Errors::runtime_check(curr>prev,"Error! hybrid_a+hybrid_b is not increasing.\n", -1);
  }

  Kokkos::deep_copy(hybrid_ai, host_hybrid_ai);
  Kokkos::deep_copy(hybrid_bi, host_hybrid_bi);

  hybrid_ai0 = host_hybrid_ai(0);

  compute_deltas();
}

void HybridVCoord::compute_deltas ()
{
  const auto host_hybrid_ai = Kokkos::create_mirror_view(hybrid_ai);
  const auto host_hybrid_bi = Kokkos::create_mirror_view(hybrid_bi);
  Kokkos::deep_copy(host_hybrid_ai, hybrid_ai);
  Kokkos::deep_copy(host_hybrid_bi, hybrid_bi);

  hybrid_ai_delta = ExecViewManaged<Scalar[NUM_LEV]>(
      "Difference in Hybrid a coordinates between consecutive interfaces");
  hybrid_bi_delta = ExecViewManaged<Scalar[NUM_LEV]>(
      "Difference in Hybrid b coordinates between consecutive interfaces");

  decltype(hybrid_ai_delta)::HostMirror host_hybrid_ai_delta =
      Kokkos::create_mirror_view(hybrid_ai_delta);
  decltype(hybrid_bi_delta)::HostMirror host_hybrid_bi_delta =
      Kokkos::create_mirror_view(hybrid_bi_delta);
  for (int level = 0; level < NUM_PHYSICAL_LEV; ++level) {
    const int ilev = level / VECTOR_SIZE;
    const int ivec = level % VECTOR_SIZE;

    host_hybrid_ai_delta(ilev)[ivec] =
        host_hybrid_ai(level + 1) - host_hybrid_ai(level);
    host_hybrid_bi_delta(ilev)[ivec] =
        host_hybrid_bi(level + 1) - host_hybrid_bi(level);
  }
  for(int level = NUM_PHYSICAL_LEV; level < NUM_LEV * VECTOR_SIZE; ++level) {
    const int ilev = level / VECTOR_SIZE;
    const int ivec = level % VECTOR_SIZE;

    host_hybrid_ai_delta(ilev)[ivec] = std::numeric_limits<Real>::quiet_NaN();
    host_hybrid_bi_delta(ilev)[ivec] = std::numeric_limits<Real>::quiet_NaN();
  }
  Kokkos::deep_copy(hybrid_ai_delta, host_hybrid_ai_delta);
  Kokkos::deep_copy(hybrid_bi_delta, host_hybrid_bi_delta);
  {
    dp0 = ExecViewManaged<Scalar[NUM_LEV]>("dp0");
    const auto hdp0 = Kokkos::create_mirror_view(dp0);
    for (int ilev = 0; ilev < NUM_LEV; ++ilev) {
      // BFB way of writing it.
      hdp0(ilev) =
          host_hybrid_ai_delta(ilev) * ps0 + host_hybrid_bi_delta(ilev) * ps0;
    }
    Kokkos::deep_copy(dp0, hdp0);
  }
}

} // namespace Homme
