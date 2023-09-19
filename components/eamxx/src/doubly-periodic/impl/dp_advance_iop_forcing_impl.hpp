#ifndef DP_ADVANCE_IOP_FORCING_IMPL_HPP
#define DP_ADVANCE_IOP_FORCING_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp advance_iop_forcing. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::plevs0(
  // Input arguments
  const Int& nver,
  const Scalar& ps,
  const uview_1d<const Spack>& hyai,
  const uview_1d<const Spack>& hyam,
  const uview_1d<const Spack>& hybi,
  const uview_1d<const Spack>& hybm,
  // Kokkos stuff
  const MemberType& team,
  // Output arguments
  const uview_1d<Spack>& pint,
  const uview_1d<Spack>& pmid,
  const uview_1d<Spack>& pdel)
{
  const auto ps0 = C::P0;
  const Int nver_pack = ekat::npack<Spack>(nver);

  // Set interface pressures
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nver_pack), [&] (Int k) {
      pint(k) = hyai(k)*ps0 + hybi(k)*ps;
      pmid(k) = hyam(k)*ps0 + hybm(k)*ps;
  });

  // Set midpoint pressures and layer thicknesses
  const auto pint_s = scalarize(pint);
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nver_pack), [&] (Int k) {
      Spack spint, spint_1;
      IntSmallPack range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
      auto range_pack2_p1_safe = range_pack1;
      range_pack2_p1_safe.set(range_pack1 > nver-1, nver-1);
      ekat::index_and_shift<1>(pint_s, range_pack2_p1_safe, spint, spint_1);
      pdel(k) = spint_1 - spint;
  });
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::advance_iop_forcing(
  // Input arguments
  const Int& plev,
  const Int& pcnst,
  const bool& have_u,
  const bool& have_v,
  const bool& dp_crm,
  const bool& use_3dfrc,
  const Scalar& scm_dt,
  const Scalar& ps_in,
  const uview_1d<const Spack>& u_in,
  const uview_1d<const Spack>& v_in,
  const uview_1d<const Spack>& t_in,
  const uview_2d<const Spack>& q_in,
  const uview_1d<const Spack>& t_phys_frc,
  const uview_1d<const Spack>& divt3d,
  const uview_2d<const Spack>& divq3d,
  const uview_1d<const Spack>& divt,
  const uview_2d<const Spack>& divq,
  const uview_1d<const Spack>& wfld,
  const uview_1d<const Spack>& uobs,
  const uview_1d<const Spack>& vobs,
  const uview_1d<const Spack>& hyai,
  const uview_1d<const Spack>& hyam,
  const uview_1d<const Spack>& hybi,
  const uview_1d<const Spack>& hybm,
  // Kokkos stuff
  const MemberType& team,
  const Workspace& workspace,
  // Output arguments
  const uview_1d<Spack>& u_update,
  const uview_1d<Spack>& v_update,
  const uview_1d<Spack>& t_update,
  const uview_2d<Spack>& q_update)
{
  // Local variables
  uview_1d<Spack>
    pmidm1, // pressure at model levels
    pintm1, // pressure at model interfaces (dim=plev+1)
    pdelm1; // pdel(k)   = pint  (k+1)-pint  (k)
  workspace.template take_many_contiguous_unsafe<3>(
    {"pmidm1", "pintm1", "pdelm1"},
    {&pmidm1, &pintm1, &pdelm1});

  // Get vertical level profiles
  plevs0(plev, ps_in, hyai, hyam, hybi, hybm, team, pintm1, pmidm1, pdelm1);

  ////////////////////////////////////////////////////////////
  //  Advance T and Q due to large scale forcing

  uview_1d<const Spack> t_lsf; // storage for temperature large scale forcing
  uview_2d<const Spack> q_lsf; // storage for moisture large scale forcing

  if (use_3dfrc) {
    t_lsf = divt3d;
    q_lsf = divq3d;
  }
  else {
    t_lsf = divt;
    q_lsf = divq;
  }

  const Int plev_pack = ekat::npack<Spack>(plev);
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
      // Initialize thermal expansion term to zero.  This term is only
      //  considered if using the preq-x dycore and if three dimensional
      //  forcing is not provided by IOP forcing file.
      Spack t_expan = 0;
      t_update(k) = t_in(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k));
      for (Int m = 0; m < pcnst; ++m) {
        q_update(m, k) = q_in(m, k) + scm_dt*q_lsf(m, k);
      }
  });

  ////////////////////////////////////////////////////////////
  //  Set U and V fields

  uview_1d<const Spack> u_src, v_src;

  if (have_v && have_u && !dp_crm) {
    u_src = uobs;
    v_src = vobs;
  }
  else {
    u_src = u_in;
    v_src = v_in;
  }

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev_pack), [&] (Int k) {
      u_update(k) = u_src(k);
      v_update(k) = v_src(k);
  });

  workspace.template release_many_contiguous<3>(
    {&pmidm1, &pintm1, &pdelm1});
}

} // namespace dp
} // namespace scream

#endif
