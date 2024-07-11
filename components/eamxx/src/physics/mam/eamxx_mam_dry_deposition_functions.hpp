#ifndef EAMXX_MAM_DRY_DEPOSITION_FUNCTIONS_HPP
#define EAMXX_MAM_DRY_DEPOSITION_FUNCTIONS_HPP

#include <ekat/kokkos/ekat_subview_utils.hpp>
#include <mam4xx/mam4.hpp>

namespace scream {

namespace {
void compute_tendencies(
    const int ncol, const int nlev, const mam4::DryDeposition dry_deposition,
    const double dt, const MAMDryDep::const_view_1d obklen,
    const MAMDryDep::const_view_1d surfric,
    const MAMDryDep::const_view_1d landfrac,
    const MAMDryDep::const_view_1d icefrac,
    const MAMDryDep::const_view_1d ocnfrac,
    const MAMDryDep::const_view_1d friction_velocity,
    const MAMDryDep::const_view_1d aerodynamical_resistance,
    MAMDryDep::view_3d qtracers, MAMDryDep::view_3d d_qtracers_dt,
    const MAMDryDep::view_3d dgncur_awet_, const MAMDryDep::view_3d wet_dens_,
    const mam_coupling::DryAtmosphere dry_atm,
    const mam_coupling::AerosolState dry_aero,
    const mam_coupling::AerosolState wet_aero, MAMDryDep::view_2d aerdepdrycw,
    MAMDryDep::view_2d aerdepdryis, MAMDryDep::view_3d tendencies) {
  const auto policy =
      ekat::ExeSpaceUtils<MAMDryDep::KT::ExeSpace>::get_default_team_policy(
          ncol, nlev);
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MAMDryDep::KT::MemberType &team) {
        const int num_aero_modes   = mam_coupling::num_aero_modes();
        const int num_aero_species = mam_coupling::num_aero_species();
        const int icol             = team.league_rank();
        const Real t               = 0;

        compute_wet_mixing_ratios(team, dry_atm, dry_aero, wet_aero, icol);
        team.team_barrier();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev), [&](const int lev) {
              for(int mode = 0; mode < num_aero_modes; ++mode) {
                int icnst = mam4::ConvProc::numptrcw_amode(mode);
                qtracers(icol, lev, icnst) =
                    wet_aero.int_aero_nmr[mode](icol, lev);
                for(int species = 0; species < num_aero_species; ++species) {
                  icnst = mam4::ConvProc::lmassptrcw_amode(species, mode);
                  if(-1 < icnst) {
                    qtracers(icol, lev, icnst) =
                        wet_aero.int_aero_mmr[mode][species](icol, lev);
                  }
                }
              }
            });
        team.team_barrier();

        mam4::Atmosphere atm    = atmosphere_for_column(dry_atm, icol);
        mam4::Prognostics progs = aerosols_for_column(dry_aero, icol);
        mam4::Surface surf;
        mam4::Diagnostics diags;
        mam4::Tendencies tends;

        for(int i = 0; i < num_aero_modes; ++i) {
          diags.wet_geometric_mean_diameter_i[i] =
              ekat::subview(dgncur_awet_, i, icol);
          diags.wet_density[i] = ekat::subview(wet_dens_, i, icol);
        }
        diags.tracer_mixing_ratio      = ekat::subview(qtracers, icol);
        diags.d_tracer_mixing_ratio_dt = ekat::subview(d_qtracers_dt, icol);
        diags.deposition_flux_of_cloud_borne_aerosols =
            ekat::subview(aerdepdrycw, icol);
        diags.deposition_flux_of_interstitial_aerosols =
            ekat::subview(aerdepdryis, icol);

        diags.Obukhov_length           = obklen[icol];
        diags.surface_friction_velocty = surfric[icol];
        diags.land_fraction            = landfrac[icol];
        diags.ice_fraction             = icefrac[icol];
        diags.ocean_fraction           = ocnfrac[icol];
        diags.friction_velocity        = friction_velocity[icol];
        diags.aerodynamical_resistance = aerodynamical_resistance[icol];

        // Fill Tendency views
        for(int m = 0; m < num_aero_modes; ++m) {
          int iconv         = mam4::ConvProc::numptrcw_amode(m);
          tends.n_mode_c[m] = ekat::subview(tendencies, icol, iconv);
          for(int a = 0; a < num_aero_species; ++a) {
            iconv = mam4::ConvProc::lmassptrcw_amode(a, m);
            if(-1 < iconv)
              tends.q_aero_c[m][a] = ekat::subview(tendencies, icol, iconv);
          }
        }

        const mam4::AeroConfig aero_config;
        dry_deposition.compute_tendencies(aero_config, team, t, dt, atm, surf,
                                          progs, diags, tends);
      });
}
}  // namespace
}  // namespace scream

#endif