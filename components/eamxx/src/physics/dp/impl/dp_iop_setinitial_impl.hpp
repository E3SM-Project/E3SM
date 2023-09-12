#ifndef DP_IOP_SETINITIAL_IMPL_HPP
#define DP_IOP_SETINITIAL_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_setinitial. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::iop_setinitial(
    const Int& plev,
    const Int& pcnst,
    const Int& nelemd,
    const Int& np,
    const Int& nstep,
    const bool& use_replay,
    const bool& dynproc,
    const bool& have_t,
    const bool& have_q,
    const bool& have_ps,
    const bool& have_u,
    const bool& have_v,
    const bool& have_numliq,
    const bool& have_cldliq,
    const bool& have_numice,
    const bool& have_cldice,
    const bool& scm_zero_non_iop_tracers,
    const bool& is_first_restart_step,
    const uview_1d<const Scalar>& qmin,
    const uview_1d<const Spack>& uobs,
    const uview_1d<const Spack>& vobs,
    const uview_1d<const Spack>& numliqobs,
    const uview_1d<const Spack>& numiceobs,
    const uview_1d<const Spack>& cldliqobs,
    const uview_1d<const Spack>& cldiceobs,
    const Scalar& psobs,
    const uview_1d<const Scalar>& dx_short,
    Scalar& dyn_dx_size,
    tracer_t& tracers,
    element_t& elem,
    const uview_1d<Spack>& tobs,
    const uview_1d<Spack>& qobs)
{
  // Made these up

  //FieldGroup g = ...;
  //const auto& names = g.m_info->m_field_names;
  constexpr Int inumliq = 1;
  constexpr Int inumice = 2;
  constexpr Int icldliq = 3;
  constexpr Int icldice = 4;

  const Int plev_packs = ekat::npack<Spack>(plev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nelemd, plev_packs);

  if (!use_replay && nstep == 0 && dynproc) {

    Kokkos::parallel_for(
      "iop_setinitial loop",
      policy,
      KOKKOS_LAMBDA(const MemberType& team) {

      const Int ie = team.league_rank();
        for (Int j = 0; j < np; ++j) {
          for (Int i = 0; i < np; ++i) {

          // Find level where tobs is no longer zero
          Int thelev=0;
          Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(team, plev_packs), [&] (int plev_pack, Int& pmin) {
            auto zmask = tobs(plev_pack) != 0;
            if (zmask.any()) {
              pmin = plev_pack;
            }
          }, Kokkos::Min<Int>(thelev));

          if (nstep <= 1) {
            Kokkos::parallel_for(
              Kokkos::TeamVectorRange(team, thelev+1), [&] (int k) {
              auto zmask = k < thelev ? Smask(true) : tobs(k) == 0;
              vector_simd for (Int p = 0; p < Spack::n; ++p) {
                if (zmask[p]) {
                  tobs(k)[p] = elem.m_forcing.m_ft(ie,i,j,k)[p];
                  qobs(k)[p] = tracers.Q(ie,0,i,j,k)[p]; // Tracer index 0 is qobs
                }
              }
            });
          }
          else {
            Kokkos::parallel_for(
              Kokkos::TeamVectorRange(team, plev_packs), [&] (int k) {
              // Ekat packs don't know how to assign themselves to KokkosKernels::Batched::Experimental::SIMD
              vector_simd for (Int p = 0; p < Spack::n; ++p) {
                tobs(k)[p] = elem.m_forcing.m_ft(ie,i,j,k)[p];
                qobs(k)[p] = tracers.Q(ie,0,i,j,k)[p];
              }
            });
          }

          if (nstep == 0) {
            if (scm_zero_non_iop_tracers) {
              for (Int cix = 0; cix < pcnst; ++cix) {
                Kokkos::parallel_for(
                  Kokkos::TeamVectorRange(team, plev_packs), [&] (int k) {
                  tracers.Q(ie, cix, i, j, k) = qmin(cix);
                });
              }
            }

            Kokkos::parallel_for(
              Kokkos::TeamVectorRange(team, thelev, plev_packs), [&] (int k) {
              auto zmask = k > thelev ? Smask(true) : tobs(k) != 0;
              if (have_t) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  if (zmask[p]) {
                    elem.m_forcing.m_ft(ie,i,j,k)[p] = tobs(k)[p];
                  }
                }
              }
              if (have_q) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  if (zmask[p]) {
                    tracers.Q(ie,0,i,j,k)[p] = qobs(k)[p];
                  }
                }
              }
            });

            if (have_ps) {
              elem.m_state.m_ps_v(ie, 0, i, j) = psobs; // what is [NUM_TIME_LEVELS]?
            }

            Kokkos::parallel_for(
              Kokkos::TeamVectorRange(team, plev_packs), [&] (int k) {
              if (have_u) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  elem.m_state.m_v(ie, k, 0, i, j, k)[p] = uobs(k)[p]; // [2] is 0?
                }
              }
              if (have_v) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  elem.m_state.m_v(ie, k, 1, i, j, k)[p] = vobs(k)[p]; // [2] is 1?
                }
              }
              if (have_numliq) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  tracers.Q(ie, inumliq, i, j, k)[p] = numliqobs(k)[p];
                }
              }
              if (have_cldliq) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  tracers.Q(ie, icldliq, i, j, k)[p] = cldliqobs(k)[p];
                }
              }
              if (have_numice) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  tracers.Q(ie, inumice, i, j, k)[p] = numiceobs(k)[p];
                }
              }
              if (have_cldice) {
                vector_simd for (Int p = 0; p < Spack::n; ++p) {
                  tracers.Q(ie, icldice, i, j, k)[p] = cldiceobs(k)[p];
                }
              }

              // If DP-CRM mode we do NOT want to write over the dy-core vertical
              // velocity with the large-scale one.  wfld is used in forecast.F90
              // for the compuation of the large-scale subsidence.
              elem.m_derived.m_omega_p(ie, i, j, k) = 0;
            });
          }
        }
      }
    });
  }

  // If DP-CRM mode then SHOC/CLUBB needs to know about grid
  // length size.  The calculations of this based on a sphere in the
  // SHOC and CLUBB interefaces are not valid for a planar grid, thus
  // save the grid length from the dycore. Note that planar dycore
  // only supports uniform grids, thus we only save one value.
  // Set this if it is the first time step or the first restart step
  if ( (nstep == 0 || is_first_restart_step) && dynproc) {
    // for (Int ie = 0; ie < nelemd; ++ie) {
    //   dyn_dx_size = dx_short(ie) * 1000; // why not just grab the last ie?
    // }
    Kokkos::parallel_reduce(
      "iop_setinitial loop",
      policy,
      KOKKOS_LAMBDA(const MemberType& team, Scalar& dyn) {
        const Int ie = team.league_rank();
        if (ie == nelemd-1) {
          dyn = dx_short(nelemd-1) * 1000;
        }
      }, dyn_dx_size);
  }
}

} // namespace dp
} // namespace scream

#endif
