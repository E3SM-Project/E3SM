#ifndef HOMMEXX_DIAGNOSTICS_THETA_HPP
#define HOMMEXX_DIAGNOSTICS_THETA_HPP

#include "Types.hpp"

#include "DiagnosticsBase.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ColumnOps.hpp"
#include "EquationOfState.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/ViewUtils.hpp"

namespace Homme
{

struct FunctorsBuffersManager;

class DiagnosticsTheta : public DiagnosticsBase
{
private:
  struct Buffers {
    static constexpr int num_3d_scalar_mid_buf = 4;
    static constexpr int num_3d_scalar_int_buf = 1;

    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    pnh;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    exner;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    phi;
    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV]>    dp_ref;

    ExecViewUnmanaged<Scalar *[NP][NP][NUM_LEV_P]>  dpnh_dp_i;
  };

public:
  void init (const ElementsState& state, const ElementsGeometry& geometry,
             const HybridVCoord& hvcoord, const bool theta_hydrostatic_mode,
             F90Ptr& elem_state_q_ptr,
             F90Ptr& elem_accum_qvar_ptr,  F90Ptr& elem_accum_qmass_ptr,
             F90Ptr& elem_accum_q1mass_ptr,F90Ptr& elem_accum_iener_ptr,
             F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr);

  int requested_buffer_size () const;
  void init_buffers (const FunctorsBuffersManager& fbm);

  void prim_energy_halftimes (const bool before_advance, const int ivar);

  HostViewUnmanaged<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> h_Q;

  struct EnergyHalfTimesTag {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const EnergyHalfTimesTag&, const TeamMember& team) const {
    KernelVariables kv(team);

    // Subview inputs/outputs
    auto vtheta_dp = Homme::subview(m_state.m_vtheta_dp, kv.ie,t1);
    auto ps_v      = Homme::subview(m_state.m_ps_v,      kv.ie,t1);

    auto phi       = Homme::subview(m_buffers.phi,       kv.team_idx);
    auto exner     = Homme::subview(m_buffers.exner,     kv.team_idx);
    auto pnh       = Homme::subview(m_buffers.pnh,       kv.team_idx);
    auto dpt1      = Homme::subview(m_buffers.dp_ref,    kv.team_idx);
    auto dpnh_dp_i = Homme::subview(m_buffers.dpnh_dp_i, kv.team_idx);

    ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV_P]> phi_i;
    if (m_theta_hydrostatic_mode) {
      phi_i = Homme::subview(m_buffers.dpnh_dp_i, kv.team_idx);
    } else {
      phi_i = Homme::subview(m_state.m_phinh_i,kv.ie,t1);
    }

    // Compute delta_hyai*ps0+delta_hybi*ps_v
    m_hvcoord.compute_dp_ref(kv,ps_v,dpt1);

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
        m_eos.compute_hydrostatic_p(kv,Homme::subview(dpt1,igp,jgp),
                                       Homme::subview(phi_i,igp,jgp),
                                       Homme::subview(pnh,igp,jgp));
        m_eos.compute_exner(kv,Homme::subview(pnh,igp,jgp),
                               Homme::subview(exner,igp,jgp));

        m_eos.compute_phi_i(kv.team,
                            Homme::subview(m_geometry.m_phis,kv.ie),
                            vtheta_dp, exner, pnh, phi_i);
      } else {
        m_eos.compute_pnh_and_exner(kv,Homme::subview(vtheta_dp,igp,jgp),
                                       Homme::subview(phi_i,igp,jgp),
                                       Homme::subview(pnh,igp,jgp),
                                       Homme::subview(exner,igp,jgp));
      }

      // Compute phi at midpoints
      m_col_ops.compute_midpoint_values(kv,Homme::subview(phi_i,igp,jgp),
                                            Homme::subview(phi,  igp,jgp));

      auto& KEner = m_KEner(kv.ie,m_ivar,igp,jgp);
      auto& PEner = m_PEner(kv.ie,m_ivar,igp,jgp);
      auto& IEner = m_IEner(kv.ie,m_ivar,igp,jgp);

      // Compute KEner
      KEner = 0.0;
      Dispatch<>::parallel_reduce(kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [=](const int ilev, Real& accumulator){
        accumulator += (u(ilev)*u(ilev) + v(ilev)*v(ilev))/2.0 * dpt1_real(ilev);
      }, KEner);

      if (!m_theta_hydrostatic_mode) {
        Real sum = 0.0;
        auto w_i = viewAsReal(Homme::subview(m_state.m_w_i,kv.ie,t1,igp,jgp));
        Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                                    [=](const int ilev, Real& accumulator){
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
                                  [=](const int ilev, Real& accumulator){
        accumulator += phi_real(ilev)*dpt1_real(ilev);
      },PEner);

      // Compute IEner
      IEner = 0.0;
      Kokkos::Real2 sum;
      Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                                  [=](const int ilev, Kokkos::Real2& accumulator){
        accumulator.v[0] += PhysicalConstants::cp*vtheta_dp_real(ilev)*exner_real(ilev);
        accumulator.v[1] += (phi_i_real(ilev+1)-phi_i_real(ilev))*pnh_real(ilev);
      },sum);

      IEner = sum.v[0] + sum.v[1] + pnh_real(0)*pnh_real(0);
    });
  }

private:

  HostViewUnmanaged<Real*[4][NP][NP]> h_IEner;
  HostViewUnmanaged<Real*[4][NP][NP]> h_KEner;
  HostViewUnmanaged<Real*[4][NP][NP]> h_PEner;

  ExecViewManaged<Real*[4][NP][NP]>   m_IEner;
  ExecViewManaged<Real*[4][NP][NP]>   m_KEner;
  ExecViewManaged<Real*[4][NP][NP]>   m_PEner;

  HybridVCoord      m_hvcoord;
  EquationOfState   m_eos;
  ElementsState     m_state;
  ElementsGeometry  m_geometry;
  ColumnOps         m_col_ops;
  Buffers           m_buffers;

  int t1,t1_qdp;

  int m_ivar;

  bool m_theta_hydrostatic_mode;
};

} // namespace Homme

#endif // HOMMEXX_DIAGNOSTICS_THETA_HPP

