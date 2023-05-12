#ifndef HOMMEXX_DIAGNOSTICS_THETA_HPP
#define HOMMEXX_DIAGNOSTICS_THETA_HPP

#include "Types.hpp"

#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ColumnOps.hpp"
#include "EquationOfState.hpp"
#include "ElementOps.hpp"
#include "HybridVCoord.hpp"
#include "Tracers.hpp"
#include "KernelVariables.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/ViewUtils.hpp"

#if ! defined(NDEBUG)
#define RESOLVE_ISSUE_WITH_ASSERTS
#endif

namespace Homme
{

struct FunctorsBuffersManager;

class Diagnostics
{
private:
  struct Buffers {
    static constexpr int num_3d_scalar_mid_buf    = 4;
    static constexpr int num_3d_scalar_int_buf_hy = 2;
    static constexpr int num_3d_scalar_int_buf_nh = 1;

    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    pnh;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    exner;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    phi;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    dp_ref;

    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV_P]>  phi_i;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV_P]>  dpnh_dp_i;
  };

//In debug regime with asserts default team sizes are too big.
#if defined RESOLVE_ISSUE_WITH_ASSERTS
  template <typename FunctorTag>
  typename std::enable_if<OnGpu<ExecSpace>::value == false,
                          Kokkos::TeamPolicy<ExecSpace, FunctorTag> >::type
  d_team_policy(const int num_exec) {
    return Homme::get_default_team_policy<ExecSpace, FunctorTag>(num_exec);
  }

  template <typename FunctorTag>
  typename std::enable_if<OnGpu<ExecSpace>::value == true,
                          Kokkos::TeamPolicy<ExecSpace, FunctorTag> >::type
  d_team_policy(const int num_exec) {
    ThreadPreferences tp;
    tp.max_threads_usable = 8;  //16
    tp.max_vectors_usable = 32; //32
    tp.prefer_larger_team = true;
    return Homme::get_default_team_policy<ExecSpace, FunctorTag>(num_exec, tp);
  }
#endif



public:


  Diagnostics (const int num_elems, const int num_tracers, const bool theta_hydrostatic_mode) :
#if ! defined(RESOLVE_ISSUE_WITH_ASSERTS)
    m_policy_diag_scalars(Homme::get_default_team_policy<ExecSpace,DiagScalarsTag>(num_elems*num_tracers)),
    m_policy_energy_halftimes(Homme::get_default_team_policy<ExecSpace,EnergyHalfTimesTag>(num_elems)),
#else
    m_policy_diag_scalars(d_team_policy<DiagScalarsTag>(num_elems*num_tracers)),
    m_policy_energy_halftimes(d_team_policy<EnergyHalfTimesTag>(num_elems)),
#endif
    m_tu(m_policy_energy_halftimes),
    m_num_elems(num_elems),
    m_num_tracers(num_tracers),
    m_theta_hydrostatic_mode(theta_hydrostatic_mode)
  {}

  void init (const ElementsState& state, const ElementsGeometry& geometry,
             const HybridVCoord& hvcoord,  const Tracers& tracers,
             F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,
             F90Ptr& elem_accum_qmass_ptr, F90Ptr& elem_accum_q1mass_ptr,
             F90Ptr& elem_accum_iener_ptr, F90Ptr& elem_accum_kener_ptr,
             F90Ptr& elem_accum_pener_ptr);

  int requested_buffer_size () const;
  void init_buffers (const FunctorsBuffersManager& fbm);

  void sync_diagnostics_to_host ();

  void run_diagnostics (const bool before_advance, const int ivar);

  void prim_diag_scalars (const bool before_advance, const int ivar);
  void prim_energy_halftimes (const bool before_advance, const int ivar);

  struct DiagScalarsTag {};
  struct EnergyHalfTimesTag {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const DiagScalarsTag&, const TeamMember& team) const {
    KernelVariables kv(team, m_num_tracers);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;

      const auto Q   = viewAsReal(Homme::subview(m_tracers.Q,   kv.ie,         kv.iq, igp, jgp));
      const auto qdp = viewAsReal(Homme::subview(m_tracers.qdp, kv.ie, t2_qdp, kv.iq, igp, jgp));

      auto& Qvar =   d_Qvar  (kv.ie,m_ivar,kv.iq,igp,jgp);
      auto& Qmass =  d_Qmass (kv.ie,m_ivar,kv.iq,igp,jgp);
      auto& Q1mass = d_Q1mass(kv.ie,       kv.iq,igp,jgp);

      // Compute Qvar
      Qvar = 0.0;
      Dispatch<>::parallel_reduce(kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [&](const int ilev, Real& accumulator){
        accumulator += qdp(ilev)*Q(ilev);
      }, Qvar);

      // Compute Qmass
      Qmass = 0.0;
      Dispatch<>::parallel_reduce(kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [&](const int ilev, Real& accumulator){
        accumulator += qdp(ilev);
      }, Qmass);

      // Compute Q1mass
      Q1mass = 0.0;
      Dispatch<>::parallel_reduce(kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [&](const int ilev, Real& accumulator){
        accumulator += qdp(ilev);
      }, Q1mass);
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const EnergyHalfTimesTag&, const TeamMember& team) const {
    // This checks that the buffers were init-ed (debug only)
    assert (m_buffers.phi.size()>0);
    assert (m_buffers.exner.size()>0);
    assert (m_buffers.pnh.size()>0);
    assert (m_buffers.dpnh_dp_i.size()>0);

    KernelVariables kv(team, m_tu);

    // Subview inputs/outputs
    auto vtheta_dp = Homme::subview(m_state.m_vtheta_dp, kv.ie,t1);
    auto dpt1      = Homme::subview(m_state.m_dp3d,      kv.ie,t1);
    auto phi_i     = m_theta_hydrostatic_mode
                   ? Homme::subview(m_buffers.phi_i,   kv.team_idx)
                   : Homme::subview(m_state.m_phinh_i, kv.ie,t1);

    auto phi       = Homme::subview(m_buffers.phi,       kv.team_idx);
    auto exner     = Homme::subview(m_buffers.exner,     kv.team_idx);
    auto pnh       = Homme::subview(m_buffers.pnh,       kv.team_idx);
    auto dpnh_dp_i = Homme::subview(m_buffers.dpnh_dp_i, kv.team_idx);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;

      // Reinterpret everything as Real* instead of Scalar*.
      // Give up some vectorization on CPU, but diagnostics are not a performance sensitive part
      auto u              = viewAsReal(Homme::subview(m_state.m_v,kv.ie,t1,0,igp,jgp));
      auto v              = viewAsReal(Homme::subview(m_state.m_v,kv.ie,t1,1,igp,jgp));
      auto pnh_real       = viewAsReal(Homme::subview(pnh,igp,jgp));
      auto phi_real       = viewAsReal(Homme::subview(phi,igp,jgp));
      auto dpt1_real      = viewAsReal(Homme::subview(dpt1,igp,jgp));
      auto phi_i_real     = viewAsReal(Homme::subview(phi_i,igp,jgp));
      auto vtheta_dp_real = viewAsReal(Homme::subview(vtheta_dp,igp,jgp));
      auto exner_real     = viewAsReal(Homme::subview(exner,igp,jgp));

      // Compute exner and pnh
      if (m_theta_hydrostatic_mode) {
        // Use phi_i to store the temporary p_i
        m_elem_ops.compute_hydrostatic_p(kv,Homme::subview(dpt1,igp,jgp),
                                            Homme::subview(phi_i,igp,jgp),
                                            Homme::subview(pnh,igp,jgp));
        m_eos.compute_exner(kv,Homme::subview(pnh,igp,jgp),
                               Homme::subview(exner,igp,jgp));

        m_eos.compute_phi_i(kv,m_geometry.m_phis(kv.ie,igp,jgp),
                               Homme::subview(vtheta_dp,igp,jgp),
                               Homme::subview(exner,igp,jgp),
                               Homme::subview(pnh,igp,jgp),
                               Homme::subview(phi_i,igp,jgp));
      } else {
        m_eos.compute_pnh_and_exner(kv,Homme::subview(vtheta_dp,igp,jgp),
                                       Homme::subview(phi_i,igp,jgp),
                                       Homme::subview(pnh,igp,jgp),
                                       Homme::subview(exner,igp,jgp));
      }

      // Compute phi at midpoints
      ColumnOps::compute_midpoint_values(kv,Homme::subview(phi_i,igp,jgp),
                                            Homme::subview(phi,  igp,jgp));

      auto& KEner = d_KEner(kv.ie,m_ivar,igp,jgp);
      auto& PEner = d_PEner(kv.ie,m_ivar,igp,jgp);
      auto& IEner = d_IEner(kv.ie,m_ivar,igp,jgp);

      // Compute KEner
      KEner = 0.0;
      Dispatch<>::parallel_reduce(kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [&](const int ilev, Real& accumulator){
        accumulator += ((u(ilev)*u(ilev) + v(ilev)*v(ilev))/2.0) * dpt1_real(ilev);
      }, KEner);

      if (!m_theta_hydrostatic_mode) {
        Real sum = 0.0;
        auto w_i = viewAsReal(Homme::subview(m_state.m_w_i,kv.ie,t1,igp,jgp));
        Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                                    [&](const int ilev, Real& accumulator){
          accumulator += (w_i(ilev)*w_i(ilev) + w_i(ilev+1)*w_i(ilev+1))/4.0 *dpt1_real(ilev);
        },sum);

        // Only one thread can update KEner
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          KEner += sum;
        });
      }

      // Compute PEner
      PEner = 0.0;
      Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                [&](const int ilev, Real& accumulator){
        accumulator += phi_real(ilev)*dpt1_real(ilev);
      },PEner);

      // Compute IEner
      IEner = 0.0;
      Kokkos::Real2 sum;
      Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [&](const int ilev, Kokkos::Real2& accumulator){
        accumulator.v[0] += PhysicalConstants::cp*vtheta_dp_real(ilev)*exner_real(ilev);
        accumulator.v[1] += (phi_i_real(ilev+1)-phi_i_real(ilev))*pnh_real(ilev);
      },sum);

      IEner = sum.v[0] + sum.v[1] + pnh_real(0)*pnh_real(0);
    });
  }

private:

  static constexpr int NUM_DIAG_TIMES = 6;

  HostViewUnmanaged<Real*[NUM_DIAG_TIMES][NP][NP]> h_IEner;
  HostViewUnmanaged<Real*[NUM_DIAG_TIMES][NP][NP]> h_KEner;
  HostViewUnmanaged<Real*[NUM_DIAG_TIMES][NP][NP]> h_PEner;

  ExecViewManaged<Real*[NUM_DIAG_TIMES][NP][NP]> d_IEner;
  ExecViewManaged<Real*[NUM_DIAG_TIMES][NP][NP]> d_KEner;
  ExecViewManaged<Real*[NUM_DIAG_TIMES][NP][NP]> d_PEner;

  HostViewUnmanaged<Real*[NUM_DIAG_TIMES][QSIZE_D][NP][NP]> h_Qvar;
  HostViewUnmanaged<Real*[NUM_DIAG_TIMES][QSIZE_D][NP][NP]> h_Qmass;
  HostViewUnmanaged<Real*                [QSIZE_D][NP][NP]> h_Q1mass;

  ExecViewManaged<Real*[NUM_DIAG_TIMES][QSIZE_D][NP][NP]> d_Qvar;
  ExecViewManaged<Real*[NUM_DIAG_TIMES][QSIZE_D][NP][NP]> d_Qmass;
  ExecViewManaged<Real*                [QSIZE_D][NP][NP]> d_Q1mass;

  HybridVCoord      m_hvcoord;
  EquationOfState   m_eos;
  ElementOps        m_elem_ops;
  ElementsState     m_state;
  ElementsGeometry  m_geometry;
  Tracers           m_tracers;
  Buffers           m_buffers;

  Kokkos::TeamPolicy<ExecSpace, DiagScalarsTag>     m_policy_diag_scalars;
  Kokkos::TeamPolicy<ExecSpace, EnergyHalfTimesTag> m_policy_energy_halftimes;
  TeamUtils<ExecSpace> m_tu;

  int t1,t1_qdp,t2_qdp;

  int m_ivar;
  int m_num_elems;
  int m_num_tracers;

  bool m_theta_hydrostatic_mode;
};

} // namespace Homme

#endif // HOMMEXX_DIAGNOSTICS_THETA_HPP

