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
  const Int& ncol,
  const Int& ncold,
  const Int& nver,
  const Spack& ps,
  // Kokkos stuff
  const MemberType& team,
  // Output arguments
  const uview_1d<Spack>& pint,
  const uview_1d<Spack>& pmid,
  const uview_1d<Spack>& pdel)
{
// !-----------------------------------------------------------------------
//   integer , intent(in)  :: ncol               ! Longitude dimension
//   integer , intent(in)  :: ncold              ! Declared longitude dimension
//   integer , intent(in)  :: nver               ! vertical dimension
//   real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
//   real(r8), intent(out) :: pint(ncold,nver+1) ! Pressure at model interfaces
//   real(r8), intent(out) :: pmid(ncold,nver)   ! Pressure at model levels
//   real(r8), intent(out) :: pdel(ncold,nver)   ! Layer thickness (pint(k+1) - pint(k))
// !-----------------------------------------------------------------------

// !---------------------------Local workspace-----------------------------
//   integer i,k             ! Longitude, level indices
// !-----------------------------------------------------------------------
// !
// ! Set interface pressures
// !
//   !$OMP PARALLEL DO PRIVATE (K, I)
//   do k=1,nver+1
//        do i=1,ncol
//             pint(i,k) = hyai(k)*ps0 + hybi(k)*ps(i)
//      end do
//   end do
// !
// ! Set midpoint pressures and layer thicknesses
// !
//   !$OMP PARALLEL DO PRIVATE (K, I)
//   do k=1,nver
//        do i=1,ncol
//             pmid(i,k) = hyam(k)*ps0 + hybm(k)*ps(i)
//             pdel(i,k) = pint(i,k+1) - pint(i,k)
//      end do
//   end do
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::advance_iop_forcing(
  // Input arguments
  const Int& plev,
  const Int& plon,
  const Int& pcnst,
  const bool& use_3dfrc,
  const Spack& scm_dt,
  const Spack& ps_in,
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
  uview_1d<Spack> pmidm1, pintm1, pdelm1, t_lsf;
  uview_2d<Spack> q_lsf;
  workspace.template take_many_contiguous_unsafe<3>(
    {"pmidm1", "pintm1", "pdelm1"},
    {&pmidm1, &pintm1, &pdelm1});

  // real(r8) pmidm1(plev)  ! pressure at model levels
  // real(r8) pintm1(plevp) ! pressure at model interfaces
  // real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
  // real(r8) t_lsf(plev)       ! storage for temperature large scale forcing
  // real(r8) q_lsf(plev,pcnst) ! storage for moisture large scale forcing
  // real(r8) fac, t_expan

  // integer k,m           ! longitude, level, constituent indices

  // Get vertical level profiles
  const Int nlon = 1; // number of columns for plevs0 routine
  plevs0(nlon, plon, plev, ps_in, team, pintm1, pmidm1, pdelm1);

  constexpr Scalar rair = C::Rair;
  constexpr Scalar cpair = C::Cpair;

  ////////////////////////////////////////////////////////////
  //  Advance T and Q due to large scale forcing

  // if (use_3dfrc) {
  //   t_lsf = divt3d;
  //   q_lsf = divq3d;
  // }
  // else {
  //   t_lsf = divt;
  //   q_lsf = divq;
  // }

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, plev), [&] (Int k) {
      // Initialize thermal expansion term to zero.  This term is only
      //  considered if using the preq-x dycore and if three dimensional
      //  forcing is not provided by IOP forcing file.
      Spack t_expan = 0;
#ifndef MODEL_THETA_L
      // this term is already taken into account through
      //  LS vertical advection in theta-l dycore
      if (!use_3dfrc) {
        t_expan = scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k));
      }
#endif
      t_update(k) = t_in(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k));
      for (Int m = 0; m < pcnst; ++m) {
        q_update(k,m) = q_in(k,m) + scm_dt*q_lsf(k,m);
      }
  });

  ////////////////////////////////////////////////////////////
  //  Set U and V fields

  Kokkos::deep_copy(u_update, u_in);
  Kokkos::deep_copy(v_update, v_in);

  workspace.template release_many_contiguous<3>(
    {&pmidm1, &pintm1, &pdelm1});
}

} // namespace dp
} // namespace scream

#endif
