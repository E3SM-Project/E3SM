#ifndef EAMXX_MAM_DRY_DEPOSITION_FUNCTIONS_HPP
#define EAMXX_MAM_DRY_DEPOSITION_FUNCTIONS_HPP

#include <ekat/kokkos/ekat_subview_utils.hpp>
#include <mam4xx/aero_config.hpp>
#include <mam4xx/convproc.hpp>
#include <mam4xx/mam4.hpp>

namespace scream {

namespace {
void compute_tendencies(
    // inputs
    const int ncol, const int nlev, const double dt,
    const MAMDryDep::const_view_1d obklen,
    const MAMDryDep::const_view_1d surfric,
    const MAMDryDep::const_view_1d landfrac,
    const MAMDryDep::const_view_1d icefrac,
    const MAMDryDep::const_view_1d ocnfrac,
    const MAMDryDep::const_view_1d friction_velocity,
    const MAMDryDep::const_view_1d aerodynamical_resistance,
    const MAMDryDep::const_view_2d fraction_landuse_,
    const MAMDryDep::const_view_3d dgncur_awet_,
    const MAMDryDep::const_view_3d wet_dens_,
    const mam_coupling::DryAtmosphere dry_atm,
    const mam_coupling::AerosolState dry_aero,

    // input-outputs
    MAMDryDep::view_3d qqcw_,

    // outputs
    MAMDryDep::view_3d ptend_q, MAMDryDep::view_2d aerdepdrycw,
    MAMDryDep::view_2d aerdepdryis,

    // work arrays
    MAMDryDep::view_2d rho_, MAMDryDep::view_4d vlc_dry_,
    MAMDryDep::view_3d vlc_trb_, MAMDryDep::view_4d vlc_grv_,
    MAMDryDep::view_3d dqdt_tmp_, MAMDryDep::view_3d qtracers) {
  static constexpr int num_aero_modes = mam_coupling::num_aero_modes();
  const auto policy =
      ekat::ExeSpaceUtils<MAMDryDep::KT::ExeSpace>::get_default_team_policy(
          ncol, nlev);

  // Parallel loop over all the columns
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MAMDryDep::KT::MemberType &team) {
        static constexpr int num_aero_species =
            mam_coupling::num_aero_species();

        const int icol = team.league_rank();
        // Parallel loop over all the levels to populate qtracers array using
        // dry_aero
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev), [&](const int lev) {
              for(int mode = 0; mode < num_aero_modes; ++mode) {
                int icnst = mam4::ConvProc::numptrcw_amode(mode);
                qtracers(icol, lev, icnst) =
                    dry_aero.int_aero_nmr[mode](icol, lev);
                for(int species = 0; species < num_aero_species; ++species) {
                  icnst = mam4::ConvProc::lmassptrcw_amode(species, mode);
                  if(-1 < icnst) {
                    qtracers(icol, lev, icnst) =
                        dry_aero.int_aero_mmr[mode][species](icol, lev);
                  }
                }
              }
            });  // parallel_for for nlevs
        team.team_barrier();

        // Create atm and progs objects
        mam4::Atmosphere atm    = atmosphere_for_column(dry_atm, icol);
        mam4::Prognostics progs = aerosols_for_column(dry_aero, icol);

        // Extract column data (or 1d view) from 2d views of data
        mam4::ConstColumnView dgncur_awet[num_aero_modes];
        mam4::ConstColumnView wet_dens[num_aero_modes];

        for(int i = 0; i < num_aero_modes; ++i) {
          dgncur_awet[i] = ekat::subview(dgncur_awet_, icol, i);
          wet_dens[i]    = ekat::subview(wet_dens_, icol, i);
        }

        mam4::ColumnView rho;
        rho = ekat::subview(rho_, icol);

        static constexpr int n_land_type = MAMDryDep::n_land_type;
        Real fraction_landuse[n_land_type];
        for(int i = 0; i < n_land_type; ++i) {
          fraction_landuse[i] = fraction_landuse_(icol, i);
        }

        static constexpr int nmodes = mam4::AeroConfig::num_modes();
        mam4::ColumnView vlc_dry[nmodes][MAMDryDep::aerosol_categories_];
        mam4::ColumnView vlc_grv[nmodes][MAMDryDep::aerosol_categories_];
        Real vlc_trb[nmodes][MAMDryDep::aerosol_categories_];

        for(int i = 0; i < nmodes; ++i) {
          for(int j = 0; j < MAMDryDep::aerosol_categories_; ++j) {
            vlc_dry[i][j] = ekat::subview(vlc_dry_, i, j, icol);
            vlc_trb[i][j] = vlc_trb_(i, j, icol);
            vlc_grv[i][j] = ekat::subview(vlc_grv_, i, j, icol);
          }
        }
        static constexpr int pcnst = mam4::aero_model::pcnst;
        mam4::ColumnView qqcw[pcnst];
        mam4::ColumnView dqdt_tmp[pcnst];
        for(int i = 0; i < pcnst; ++i) {
          qqcw[i]     = ekat::subview(qqcw_, i, icol);
          dqdt_tmp[i] = ekat::subview(dqdt_tmp_, i, icol);
        }
        // Extract qqcw from Prognostics
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
          for(int m = 0; m < nmodes; ++m) {
            qqcw[mam4::ConvProc::numptrcw_amode(m)][kk] = progs.n_mode_c[m][kk];
            for(int a = 0; a < mam4::AeroConfig::num_aerosol_ids(); ++a)
              if(-1 < mam4::ConvProc::lmassptrcw_amode(a, m))
                qqcw[mam4::ConvProc::lmassptrcw_amode(a, m)][kk] =
                    progs.q_aero_c[m][a][kk];
          }
        });  // parallel_for nlevs
        team.team_barrier();
        bool ptend_lq[pcnst];  // currently unused
        mam4::aero_model_drydep(
            // inputs
            team, fraction_landuse, atm.temperature, atm.pressure,
            atm.interface_pressure, atm.hydrostatic_dp,
            ekat::subview(qtracers, icol), dgncur_awet, wet_dens, obklen[icol],
            surfric[icol], landfrac[icol], icefrac[icol], ocnfrac[icol],
            friction_velocity[icol], aerodynamical_resistance[icol], dt,
            // input-outputs
            qqcw,
            // outputs
            ekat::subview(ptend_q, icol), ptend_lq,
            ekat::subview(aerdepdrycw, icol), ekat::subview(aerdepdryis, icol),
            // work arrays
            rho, vlc_dry, vlc_trb, vlc_grv, dqdt_tmp);
      });  // parallel_for for ncols
}  // Compute_tendencies ends

// Update interstitial aerosols using ptend_q tendencies
void update_interstitial_mmrs(const MAMDryDep::view_3d ptend_q, const double dt,
                              const int ncol, const int nlev,
                              // output
                              const mam_coupling::AerosolState dry_aero) {
  const auto policy =
      ekat::ExeSpaceUtils<MAMDryDep::KT::ExeSpace>::get_default_team_policy(
          ncol, nlev);
  static constexpr int nmodes = mam4::AeroConfig::num_modes();
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MAMDryDep::KT::MemberType &team) {
        const int icol = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
          for(int m = 0; m < nmodes; ++m) {
            dry_aero.int_aero_nmr[m](icol, kk) +=
                ptend_q(icol, kk, mam4::ConvProc::numptrcw_amode(m)) * dt;
            for(int a = 0; a < mam4::AeroConfig::num_aerosol_ids(); ++a)
              if(-1 < mam4::ConvProc::lmassptrcw_amode(a, m))
                dry_aero.int_aero_mmr[m][a](icol, kk) +=
                    ptend_q(icol, kk, mam4::ConvProc::lmassptrcw_amode(a, m)) *
                    dt;
          }
        });  // parallel_for nlevs
      });    // parallel_for icol
}  // Update interstitial aerosols ends

// Update cloud borne aerosols using qqcw
void update_cloudborne_mmrs(const MAMDryDep::view_3d qqcw, const double dt,
                            const int nlev_,
                            // output
                            const mam_coupling::AerosolState dry_aero) {
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    Kokkos::deep_copy(dry_aero.cld_aero_nmr[m],
                      ekat::subview(qqcw, mam4::ConvProc::numptrcw_amode(m)));
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      if(dry_aero.cld_aero_mmr[m][a].data()) {
        Kokkos::deep_copy(
            dry_aero.cld_aero_mmr[m][a],
            ekat::subview(qqcw, mam4::ConvProc::lmassptrcw_amode(a, m)));
      }
    }
  }
}  // Update cloud borne aerosols ends

}  // namespace
}  // namespace scream

#endif
