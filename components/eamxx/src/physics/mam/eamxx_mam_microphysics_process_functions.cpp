#include <mam4xx/mam4.hpp>
#include <physics/mam/eamxx_mam_microphysics_process_interface.hpp>

namespace scream {

void MAMMicrophysics::run_small_kernels_microphysics(const double dt, const double eccf)
{
  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  const int nlev = nlev_;
  //
  // set external forcing
  constexpr int extcnt = mam4::gas_chemistry::extcnt;
  const auto& extfrc = extfrc_;
  Kokkos::deep_copy(extfrc,0.0);
  for (int mm = 0; mm < extcnt; ++mm) {
    // Fortran to C++ indexing
    const int nn = forcings_[mm].frc_ndx - 1;
    for (int isec = 0; isec < forcings_[mm].nsectors; ++isec) {
    const auto& field = forcings_[mm].fields[isec];
    if (forcings_[mm].file_alt_data) {
      Kokkos::parallel_for(
        "MAMMicrophysics::run_impl::forcing", policy,
        KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol     = team.league_rank();   // column index
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev),
         [&](int kk) {
          extfrc(icol,kk,nn) += field(icol,nlev - 1 - kk);
        });
      });

   } else {
     Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::forcing", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
        extfrc(icol,kk,nn) += field(icol,kk);
      });
    });
   }
  } // isec
  }   // end mm
  // set external forcing ends

  // set invariants
  const auto& invariants = invariants_;
  const auto& cnst_offline=cnst_offline_;
  const auto& dry_atm = dry_atm_;
  Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::setinv", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      // fetch column-specific atmosphere state data
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
      const auto invariants_icol = ekat::subview(invariants, icol);
      view_1d cnst_offline_icol[mam4::mo_setinv::num_tracer_cnst];
        for(int i = 0; i < mam4::mo_setinv::num_tracer_cnst; ++i) {
          cnst_offline_icol[i] = ekat::subview(cnst_offline[i], icol);
        }
      mam4::mo_setinv::setinv(team,                                    // in
                          invariants_icol,                         // out
                          atm.temperature, atm.vapor_mixing_ratio, // in
                          cnst_offline_icol, atm.pressure);        // in

  });

  // set invariants ends

  // set o3_col_dens
  // set extract_stateq
  const int offset_aerosol = mam4::utils::gasses_start_ind();
  const auto& dry_aero = dry_aero_;
  const int pcnst              = mam4::pcnst;
  const auto& state_q = state_q_;
  const auto& qqcw_pcnst = qqcw_pcnst_;
  const auto& qq = qq_;
  const auto& qqcw = qqcw_;
  const auto& vmr = vmr_;
  const auto& vmrcw = vmrcw_;
  constexpr int num_gas_aerosol_constituents = mam_coupling::gas_pcnst();
  Real adv_mass_kg_per_moles[num_gas_aerosol_constituents];
  for(int i = 0; i < num_gas_aerosol_constituents; ++i) {
    // NOTE: state_q is kg/kg-dry-air; adv_mass is in g/mole.
    // Convert adv_mass to kg/mole as vmr_from_mmr function uses
    // molec_weight_dry_air with kg/mole units
    adv_mass_kg_per_moles[i] = mam4::gas_chemistry::adv_mass[i] / 1e3;
  }

  Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::extract_stateq", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
  // extract aerosol state variables into "working arrays" (mass
    // mixing ratios) (in EAM, this is done in the gas_phase_chemdr
    // subroutine defined within
    //  mozart/mo_gas_phase_chemdr.F90)
    const int icol     = team.league_rank();   // column index
    const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
    mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

    const auto state_q_icol = ekat::subview(state_q,icol);
    const auto qqcw_pcnst_icol = ekat::subview(qqcw_pcnst,icol);
    const auto qq_icol = ekat::subview(qq,icol);
    const auto qqcw_icol = ekat::subview(qqcw,icol);
    const auto vmr_icol = ekat::subview(vmr,icol);
    const auto vmrcw_icol = ekat::subview(vmrcw,icol);
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, nlev),
      [&](const int kk) {
        const auto state_q_kk = ekat::subview(state_q_icol,kk);
        const auto qqcw_pcnst_kk = ekat::subview(qqcw_pcnst_icol,kk);
        const auto qq_kk = ekat::subview(qq_icol,kk);
        const auto qqcw_kk = ekat::subview(qqcw_icol,kk);
        const auto vmr_kk = ekat::subview(vmr_icol,kk);
        const auto vmrcw_kk = ekat::subview(vmrcw_icol,kk);
        // output (state_q)
        mam4::utils::extract_stateq_from_prognostics(progs, atm, state_q_kk, kk);
        // output (qqcw_pcnst)
        mam4::utils::extract_qqcw_from_prognostics(progs, qqcw_pcnst_kk, kk);
        for (int i = offset_aerosol; i < pcnst; ++i) {
          qq_kk[i - offset_aerosol] = state_q_kk[i];
          qqcw_kk[i - offset_aerosol] = qqcw_pcnst_kk[i];
        }
        // convert mass mixing ratios to volume mixing ratios (VMR),
        // equivalent to tracer mixing ratios (TMR)
        // output (vmr)
        mam4::microphysics::mmr2vmr(qq_kk.data(), adv_mass_kg_per_moles, vmr_kk.data());
        // output (vmrcw)
        mam4::microphysics::mmr2vmr(qqcw_kk.data(), adv_mass_kg_per_moles, vmrcw_kk.data());
    });
  });


  const auto& o3_col_dens = buffer_.scratch[8];
  const auto& work_set_het  = work_set_het_;
  Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::compute_o3_column_density", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    // calculate o3 column densities (first component of col_dens in Fortran
    // code)
    auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
    const auto& invariants_icol = ekat::subview(invariants, icol);
    const auto& atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
    const auto& vmr_icol = ekat::subview(vmr,icol);
    const auto work_set_het_icol = ekat::subview(work_set_het, icol);
    auto work_set_het_ptr = (Real *)work_set_het_icol.data();
    const auto& o3_col_deltas  = view_1d(work_set_het_ptr, mam4::nlev + 1);
    // NOTE: if we need o2 column densities, set_ub_col and setcol must be changed
    const int nlev = mam4::nlev;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](const int kk) {
      const Real pdel = atm.hydrostatic_dp(kk);
      const auto vmr_kk = ekat::subview(vmr_icol,kk);
      const auto invariants_k = ekat::subview(invariants_icol,kk);
      // compute the change in o3 density for this column above its neighbor
      mam4::mo_photo::set_ub_col(o3_col_deltas(kk + 1),     // out
                               vmr_kk.data(), invariants_k.data(), pdel); // out
    });
    team.team_barrier();
    // sum the o3 column deltas to densities
    mam4::mo_photo::setcol(team, o3_col_deltas.data(), // in
                         o3_col_dens_i);        // out
  });
  // set o3_col_dens ends

  // set photo table
  const auto &photo_rates                     = photo_rates_;
  const auto& d_sfc_alb_dir_vis = d_sfc_alb_dir_vis_;
  Kokkos::deep_copy(photo_rates_,0.0);
  const auto& work_photo_table =work_photo_table_;
  const auto& zenith_angle = acos_cosine_zenith_;
  const auto& photo_table=photo_table_;

  Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::photo_table", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    const auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
    const auto &work_photo_table_icol =
            ekat::subview(work_photo_table, icol);
    const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
    const auto &photo_rates_icol = ekat::subview(photo_rates, icol);
        // set up photolysis work arrays for this column.
    mam4::mo_photo::PhotoTableWorkArrays photo_work_arrays_icol;

    //  set work view using 1D photo_work_arrays_icol
    // Note: We are not allocating views here.
    mam4::mo_photo::set_photo_table_work_arrays(photo_table,
                                              work_photo_table_icol,   // in
                                              photo_work_arrays_icol); // out

    team.team_barrier();
    mam4::mo_photo::table_photo(team, photo_rates_icol,                    // out
                              atm.pressure, atm.hydrostatic_dp,          // in
                              atm.temperature, o3_col_dens_i,            // in
                              zenith_angle(icol), d_sfc_alb_dir_vis(icol), // in
                              atm.liquid_mixing_ratio, atm.cloud_fraction, // in
                              eccf, photo_table,                           // in
                              photo_work_arrays_icol); // out
    });

    // set photo table ends
    const auto& col_latitudes = col_latitudes_;
    const auto& het_rates =het_rates_;
    const auto &cmfdqr       = cmfdqr_;

    // Stratiform rain production rate [kg/kg/s]
    const auto &prain =
      get_field_in("precip_total_tend").get_view<const Real **>();
    // Evaporation from stratiform rain [kg/kg/s]
    const auto &nevapr = get_field_in("nevapr").get_view<const Real **>();

    // set sethet
    Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::sethet", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
    const auto phis= dry_atm.phis(icol);
    const Real col_lat = col_latitudes(icol);  // column latitude (degrees?)
    // convert column latitude to radians
    const Real rlats = col_lat * M_PI / 180.0;
    const auto prain_icol        = ekat::subview(prain, icol);
    const auto nevapr_icol       = ekat::subview(nevapr, icol);
    const auto invariants_icol = ekat::subview(invariants, icol);
    const auto het_rates_icol  = ekat::subview(het_rates, icol);
    const auto vmr_icol = ekat::subview(vmr, icol);

    const auto work_set_het_icol = ekat::subview(work_set_het, icol);
    auto work_set_het_ptr = (Real *)work_set_het_icol.data();
    const int sethet_work_len = mam4::mo_sethet::get_work_len_sethet();
    const auto work_sethet_call = view_1d(work_set_het_ptr, sethet_work_len);
    work_set_het_ptr += sethet_work_len;
    mam4::mo_sethet::sethet(team, atm, het_rates_icol, rlats, phis, cmfdqr, prain_icol,
                          nevapr_icol, dt, invariants_icol, vmr_icol,
                          work_sethet_call);
  });

  // set_het end

  // set drydep_xactive

  const auto& index_season_lai=index_season_lai_;

  // U wind component [m/s]
  const const_view_2d u_wind =
      get_field_in("horiz_winds").get_component(0).get_view<const Real **>();

  // V wind component [m/s]
  const const_view_2d v_wind =
      get_field_in("horiz_winds").get_component(1).get_view<const Real **>();

  // Liquid precip [kg/m2]
  const const_view_1d precip_liq_surf_mass =
      get_field_in("precip_liq_surf_mass").get_view<const Real *>();

  // drydep_xactive
  // Snow depth on land [m]
  const const_view_1d snow_depth_land =
      get_field_in("snow_depth_land").get_view<const Real *>();

  // Fractional land use [fraction]
  const const_view_2d fraction_landuse =
      get_field_in("fraction_landuse").get_view<const Real **>();

    // Downwelling solar flux at the surface [w/m2]
  const const_view_2d sw_flux_dn =
      get_field_in("SW_flux_dn").get_view<const Real **>();

  const mam4::seq_drydep::Data drydep_data =
      mam4::seq_drydep::set_gas_drydep_data();

  const int month              = start_of_step_ts().get_month();  // 1-based

  // Surface temperature [K]
  const const_view_1d sfc_temperature =
      get_field_in("surf_radiative_T").get_view<const Real *>();

  // Surface pressure [Pa]
  const const_view_1d sfc_pressure =
      get_field_in("ps").get_view<const Real *>();

  // Constituent fluxes of gas and aerosol species
  view_2d constituent_fluxes =
      get_field_out("constituent_fluxes").get_view<Real **>();

  // Ice precip [kg/m2]
  const const_view_1d precip_ice_surf_mass =
      get_field_in("precip_ice_surf_mass").get_view<const Real *>();

  Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::drydep_xactive", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);

    const auto work_set_het_icol = ekat::subview(work_set_het, icol);
    auto work_set_het_ptr = (Real *)work_set_het_icol.data();
    // deposition velocity [1/cm/s]
    const auto& dflx_col = view_1d(work_set_het_ptr, num_gas_aerosol_constituents);
    work_set_het_ptr += num_gas_aerosol_constituents;
    // deposition flux [1/cm^2/s]
    const auto& dvel_col = view_1d(work_set_het_ptr, num_gas_aerosol_constituents);
    work_set_het_ptr += num_gas_aerosol_constituents;

    // Snow depth on land [m]
    const Real snow_height = snow_depth_land(icol);

    Real fraction_landuse_icol[mam4::mo_drydep::n_land_type];
    for(int i = 0; i < mam4::mo_drydep::n_land_type; ++i) {
      fraction_landuse_icol[i] = fraction_landuse(icol, i);
    }

    int index_season[mam4::mo_drydep::n_land_type];
    {
      //-------------------------------------------------------------------------------------
      // define which season (relative to Northern hemisphere climate)
      //-------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------
      // define season index based on fixed LAI
      //-------------------------------------------------------------------------------------
      for(int lt = 0; lt < mam4::mo_drydep::n_land_type; ++lt) {
            index_season[lt] = index_season_lai(icol, month - 1);
      }

      //-------------------------------------------------------------------------------------
      // special case for snow covered terrain
      //-------------------------------------------------------------------------------------
      if(snow_height > 0.01) {  // BAD_CONSTANT
        for(int lt = 0; lt < mam4::mo_drydep::n_land_type; ++lt) {
          index_season[lt] = 3;
        }
        }
      }
      const int surface_lev = nlev - 1; // Surface level
      // specific humidity [kg/kg]
      const Real spec_hum = atm.vapor_mixing_ratio(surface_lev);
      // surface air temperature [K]
      const Real air_temp = atm.temperature(surface_lev);
      // potential temperature [K] *(temp*(1+vapor_mixing_ratio))
      //(FIXME: We followed Fortran, compare it with MAM4xx's potential temp
      // func)
      const Real tv = air_temp * (1.0 + spec_hum);
      // 10-meter pressure [Pa]
      // Surface pressure at 10m (Followed the fortran code)
      const Real pressure_10m = atm.pressure(surface_lev);

      // Wind speed at the surface
      const Real wind_speed =
            haero::sqrt(u_wind(icol, surface_lev) * u_wind(icol, surface_lev) +
                        v_wind(icol, surface_lev) * v_wind(icol, surface_lev));
      // Total rain at the surface
      const Real rain =
            precip_liq_surf_mass(icol) + precip_ice_surf_mass(icol);
      // Downwelling solar flux at the surface (value at interface) [w/m2]
      const Real solar_flux = sw_flux_dn(icol, surface_lev + 1);
      const auto qq_icol = ekat::subview(qq,icol);
      const auto qq_sfc = ekat::subview(qq_icol,surface_lev);



      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
       mam4::mo_drydep::drydep_xactive(
        drydep_data,
        fraction_landuse_icol, // fraction of land use for column by land type
        index_season,     // column-specific mapping of month indices to
                          // seasonal land-type indices [-]
        sfc_temperature(icol),         // surface temperature [K]
        air_temp,         // surface air temperature [K]
        tv,               // potential temperature [K]
        sfc_pressure(icol),     // surface pressure [Pa]
        pressure_10m,     // 10-meter pressure [Pa]
        spec_hum,         // specific humidity [kg/kg]
        wind_speed,       // 10-meter wind spped [m/s]
        rain,             // rain content [??]
        solar_flux,       // direct shortwave surface radiation [W/m^2]
        qq_sfc.data(),               // constituent MMRs [kg/kg]
        dvel_col.data(),             // deposition velocity [1/cm/s]
        dflx_col.data()              // deposition flux [1/cm^2/s]
      );
      });
      // Update constituent fluxes with gas drydep fluxes (dflx)
        // FIXME: Possible units mismatch (dflx is in kg/cm2/s but
        // constituent_fluxes is kg/m2/s) (Following mimics Fortran code
        // behavior but we should look into it)
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, offset_aerosol, pcnst), [&](int ispc) {
          constituent_fluxes(icol, ispc) -= dflx_col(ispc - offset_aerosol);
        });
    });


    // set drydep_xactive ends

    // gas_phase_chemistry

    // Store mixing ratios before gas chemistry changes the mixing ratios
    const auto& vmr0 = vmr0_;
    Kokkos::deep_copy(vmr0,vmr);
    // NOTE: Making copies of clsmap_4 and permute_4 to fix undefined arrays on
    // the device.
    int clsmap_4[num_gas_aerosol_constituents], permute_4[num_gas_aerosol_constituents];
    for(int i = 0; i < num_gas_aerosol_constituents; ++i) {
      clsmap_4[i]              = mam4::gas_chemistry::clsmap_4[i];
      permute_4[i]             = mam4::gas_chemistry::permute_4[i];
    }

    view_3d gas_phase_chemistry_dvmrdt;
    if (extra_mam4_aero_microphys_diags_) {
      gas_phase_chemistry_dvmrdt = get_field_out("mam4_microphysics_tendency_gas_phase_chemistry").get_view<Real ***>();
    }

    Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::gas_phase_chemistry", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
      const auto &photo_rates_icol = ekat::subview(photo_rates, icol);
      const auto invariants_icol = ekat::subview(invariants, icol);
      const auto extfrc_icol = ekat::subview(extfrc, icol);
      const auto het_rates_icol = ekat::subview(het_rates, icol);
      const auto& vmr_icol = ekat::subview(vmr, icol);

      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        const auto &extfrc_k = ekat::subview(extfrc_icol, kk);
        const auto &invariants_k = ekat::subview(invariants_icol, kk);
        const auto &photo_rates_k = ekat::subview(photo_rates_icol, kk);
        const auto &het_rates_k = ekat::subview(het_rates_icol, kk);
        // extract atm state variables (input)
        const Real temperature = atm.temperature(kk);
        const auto &vmr_kk = ekat::subview(vmr_icol, kk);
        mam4::microphysics::gas_phase_chemistry(
        // in
        temperature, dt, photo_rates_k.data(), extfrc_k.data(), invariants_k.data(),
        clsmap_4, permute_4, het_rates_k.data(),
        // out
        vmr_kk);
      });

    });

    if (gas_phase_chemistry_dvmrdt.size()) {

    Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::gas_phase_chemistry_dvmrdt", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        Real pdel = atm.hydrostatic_dp(kk);
        const Real mbar = haero::Constants::molec_weight_dry_air;
      const Real gravit = haero::Constants::gravity;
      const Real x = 1.0 / mbar * pdel / gravit;
      for (int m = 0; m < num_gas_aerosol_constituents; ++m)
        gas_phase_chemistry_dvmrdt(icol, m, kk) =
            x * adv_mass_kg_per_moles[m] * (vmr(icol,kk,m) - vmr0(icol,kk,m)) / dt;
       });

    });
  }
  // gas_phase_chemistry ends

  // setsox_single_level
  const Config &config                        = config_;
  const auto& vmr_pregas =vmr_pregas_;
  const auto& vmr_precld=vmr_precld_;
  Kokkos::deep_copy(vmr_pregas, vmr);
  Kokkos::deep_copy(vmr_precld, vmrcw );

  const auto& vmr_bef_aq_chem = vmr_pregas;
  const auto& config_setsox = config.setsox;
  const auto& dqdt_aqso4 = dqdt_aqso4_;
  const auto& dqdt_aqh2so4 = dqdt_aqh2so4_;

  Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::setsox_single_level", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {

    const int icol     = team.league_rank();   // column index
    const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
    //----------------------
    // Aerosol microphysics
    //----------------------
    // the logic below is taken from the aero_model_gasaerexch
    // subroutine in eam/src/chemistry/modal_aero/aero_model.F90
    const auto dqdt_aqso4_icol = ekat::subview(dqdt_aqso4,icol);
    const auto dqdt_aqh2so4_icol =ekat::subview(dqdt_aqh2so4,icol);

    const auto invariants_icol = ekat::subview(invariants, icol);
    const auto & vmrcw_icol = ekat::subview(vmrcw,icol);
    const auto & vmr_icol = ekat::subview(vmr,icol);

    Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        // extract atm state variables (input)
        Real temp = atm.temperature(kk);
        Real pmid = atm.pressure(kk);
        Real pdel = atm.hydrostatic_dp(kk);
        Real lwc = atm.liquid_mixing_ratio(kk);
        Real cldfrac = atm.cloud_fraction(kk);
        Real cldnum = atm.cloud_liquid_number_mixing_ratio(kk);
        const auto &invariants_k = ekat::subview(invariants_icol, kk);
        // aqueous chemistry ...
       constexpr Real mbar = haero::Constants::molec_weight_dry_air;
       constexpr int indexm = mam4::gas_chemistry::indexm;
       const auto &dqdt_aqso4_k = ekat::subview(dqdt_aqso4_icol, kk);
       const auto &dqdt_aqh2so4_k = ekat::subview(dqdt_aqh2so4_icol, kk);
       const auto & vmrcw_k = ekat::subview(vmrcw_icol,kk);
       const auto & vmr_k = ekat::subview(vmr_icol,kk);

    mam4::mo_setsox::setsox_single_level(
        // in
        offset_aerosol, dt, pmid, pdel, temp, mbar, lwc, cldfrac, cldnum,
        invariants_k(indexm), config_setsox,
        // out
        dqdt_aqso4_k.data(), dqdt_aqh2so4_k.data(), vmrcw_k.data(), vmr_k.data());
    });
    });

    view_3d aqueous_chemistry_dvmrdt;
    if (extra_mam4_aero_microphys_diags_) {
      aqueous_chemistry_dvmrdt = get_field_out("mam4_microphysics_tendency_aqueous_chemistry").get_view<Real ***>();
    }

    if (aqueous_chemistry_dvmrdt.size()) {
      Kokkos::parallel_for(
      "MAMMicrophysics::run_impl::aqueous_chemistry_dvmrdt", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol     = team.league_rank();   // column index
        const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
        Kokkos::parallel_for(
         Kokkos::TeamVectorRange(team, nlev),
         [&](const int kk) {
         Real pdel = atm.hydrostatic_dp(kk);
         const Real mbar = haero::Constants::molec_weight_dry_air;
         const Real gravit = haero::Constants::gravity;
         const Real x = 1.0 / mbar * pdel / gravit;
         for (int m = 0; m < num_gas_aerosol_constituents; ++m)
          aqueous_chemistry_dvmrdt(icol, m, kk) =
            x * adv_mass_kg_per_moles[m] * (vmr(icol,kk,m) - vmr_bef_aq_chem(icol,kk,m)) / dt;
         });

      });
    }
    //setsox_single_level ends

    // modal_aero_amicphys_intr
    const auto wet_diameter =
      get_field_in("dgnumwet").get_view<const Real ***>();
    const auto dry_diameter =
      get_field_in("dgnum").get_view<const Real ***>();
    const auto wetdens = get_field_in("wetdens").get_view<const Real ***>();

    auto& dgncur_awet_loc = dgncur_awet_;
    auto& dgncur_a_loc = dgncur_a_;
    auto& wetdens_loc = wetdens_;
    constexpr int nmodes = mam4::AeroConfig::num_modes();
     Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::modal_aero_amicphys_intr_precompute", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
    Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        for (int imode = 0; imode < nmodes; imode++) {
         dgncur_awet_loc(icol, kk, imode) = wet_diameter(icol, imode, kk);
         dgncur_a_loc(icol, kk, imode) = dry_diameter(icol, imode, kk);
         wetdens_loc(icol, kk, imode) = wetdens(icol, imode, kk);
        }
      });
    });

    view_3d gas_aero_exchange_condensation, gas_aero_exchange_renaming,
          gas_aero_exchange_nucleation, gas_aero_exchange_coagulation,
          gas_aero_exchange_renaming_cloud_borne;

    if (extra_mam4_aero_microphys_diags_) {
      gas_aero_exchange_condensation = get_field_out("mam4_microphysics_tendency_condensation").get_view<Real***>();
      gas_aero_exchange_renaming = get_field_out("mam4_microphysics_tendency_renaming").get_view<Real***>();
      gas_aero_exchange_nucleation = get_field_out("mam4_microphysics_tendency_nucleation").get_view<Real***>();
      gas_aero_exchange_coagulation = get_field_out("mam4_microphysics_tendency_coagulation").get_view<Real***>();
      gas_aero_exchange_renaming_cloud_borne = get_field_out("mam4_microphysics_tendency_renaming_cloud_borne").get_view<Real***>();
    }

    const bool extra_mam4_aero_microphys_diags  = extra_mam4_aero_microphys_diags_;

    const auto& config_amicphys = config.amicphys;
     Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::modal_aero_amicphys_intr", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {

      const int icol     = team.league_rank();   // column index
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);

      mam4::MicrophysDiagnosticArrays diag_arrays;
      if (extra_mam4_aero_microphys_diags) {
          diag_arrays.gas_aero_exchange_condensation = ekat::subview(gas_aero_exchange_condensation, icol);
          diag_arrays.gas_aero_exchange_renaming = ekat::subview(gas_aero_exchange_renaming, icol);
          diag_arrays.gas_aero_exchange_nucleation = ekat::subview(gas_aero_exchange_nucleation, icol);
          diag_arrays.gas_aero_exchange_coagulation = ekat::subview(gas_aero_exchange_coagulation, icol);
          diag_arrays.gas_aero_exchange_renaming_cloud_borne = ekat::subview(gas_aero_exchange_renaming_cloud_borne, icol);
      }
      auto vmrcw_icol = ekat::subview(vmrcw,icol);
      auto vmr_icol = ekat::subview(vmr,icol);
      const auto & vmr0_icol = ekat::subview(vmr0,icol);
      const auto & vmr_pregas_icol = ekat::subview(vmr_pregas,icol);
      const auto & vmr_precld_icol = ekat::subview(vmr_precld,icol);

      const auto& dgncur_awet_icol = ekat::subview(dgncur_awet_loc,icol);
      const auto& dgncur_a_icol = ekat::subview(dgncur_a_loc,icol);
      const auto& wetdens_icol = ekat::subview(wetdens_loc,icol);

      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        // calculate aerosol water content using water uptake treatment
        // * dry and wet diameters [m]
        // * wet densities [kg/m3]
        // * aerosol water mass mixing ratio [kg/kg]
        const Real temp = atm.temperature(kk);
        const Real pmid = atm.pressure(kk);
        const Real pdel = atm.hydrostatic_dp(kk);
        const Real zm = atm.height(kk);
        const Real pblh = atm.planetary_boundary_layer_height;
        const Real qv = atm.vapor_mixing_ratio(kk);
        const Real cldfrac = atm.cloud_fraction(kk);
        const auto& dgncur_awet_kk = ekat::subview(dgncur_awet_icol,kk);
        const auto& dgncur_a_kk = ekat::subview(dgncur_a_icol,kk);
        const auto& wetdens_kk = ekat::subview(wetdens_icol,kk);

        auto vmr_kk = ekat::subview(vmr_icol,kk);
        auto vmrcw_kk = ekat::subview(vmrcw_icol,kk);
        const auto & vmr0_kk = ekat::subview(vmr0_icol,kk);
        const auto & vmr_pregas_kk = ekat::subview(vmr_pregas_icol,kk);
        const auto & vmr_precld_kk = ekat::subview(vmr_precld_icol,kk);
    // Perform aerosol microphysics (gas-aerosol exchange, nucleation,
    // coagulation)
    mam4::microphysics::modal_aero_amicphys_intr(
        // in
        config_amicphys, dt, temp, pmid, pdel, zm, pblh, qv, cldfrac,
        // out
        vmr_kk, vmrcw_kk,
        // diagnostics (out)
        kk, diag_arrays.gas_aero_exchange_condensation,
        diag_arrays.gas_aero_exchange_renaming,
        diag_arrays.gas_aero_exchange_nucleation,
        diag_arrays.gas_aero_exchange_coagulation,
        diag_arrays.gas_aero_exchange_renaming_cloud_borne,
        // in
        vmr0_kk, vmr_pregas_kk, vmr_precld_kk, dgncur_a_kk, dgncur_awet_kk, wetdens_kk);
      });
    });

    // modal_aero_amicphys_intr ends

    // vmr2mmr_cw
    Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::vmr2mmr_cw", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      auto vmrcw_icol = ekat::subview(vmrcw,icol);
      auto qqcw_icol = ekat::subview(qqcw,icol);
      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        auto vmrcw_kk = ekat::subview(vmrcw_icol,kk);
        auto qqcw_kk = ekat::subview(qqcw_icol,kk);
        mam4::microphysics::vmr2mmr(vmrcw_kk.data(),
           adv_mass_kg_per_moles, qqcw_kk.data());
       });
    });
    // vmr2mmr_cw ends

    // linoz
    if (config.linoz.compute) {

      // climatology data for linear stratospheric chemistry
      // ozone (climatology) [vmr]
      view_2d linoz_o3_clim =  buffer_.scratch[0];
      // column o3 above box (climatology) [Dobson Units (DU)]
      view_2d linoz_o3col_clim = buffer_.scratch[1];
      // temperature (climatology) [K]
      view_2d linoz_t_clim = buffer_.scratch[2];
      // P minus L (climatology) [vmr/s]
      view_2d linoz_PmL_clim = buffer_.scratch[3];
      // sensitivity of P minus L to O3 [1/s]
      view_2d linoz_dPmL_dO3 = buffer_.scratch[4];
      // sensitivity of P minus L to T3 [K]
      view_2d linoz_dPmL_dT = buffer_.scratch[5];
      // sensitivity of P minus L to overhead O3 column [vmr/DU]
      view_2d linoz_dPmL_dO3col = buffer_.scratch[6];
      // Cariolle parameter for PSC loss of ozone [1/s]
      view_2d linoz_cariolle_pscs = buffer_.scratch[7];
      const auto& linoz_conf=config.linoz;
      const int o3_ndx = static_cast<int>(mam4::GasId::O3);
      Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::linoz", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const Real col_lat = col_latitudes(icol);  // column latitude (degrees?)
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
      const auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
      // convert column latitude to radians
      const Real rlats = col_lat * M_PI / 180.0;
      mam4::microphysics::LinozData linoz_data;
      if (config.linoz.compute) {
          linoz_data.linoz_o3_clim_icol = ekat::subview(linoz_o3_clim, icol);
          linoz_data.linoz_t_clim_icol  = ekat::subview(linoz_t_clim, icol);
          linoz_data.linoz_o3col_clim_icol =
            ekat::subview(linoz_o3col_clim, icol);
          linoz_data.linoz_PmL_clim_icol = ekat::subview(linoz_PmL_clim, icol);
          linoz_data.linoz_dPmL_dO3_icol = ekat::subview(linoz_dPmL_dO3, icol);
          linoz_data.linoz_dPmL_dT_icol  = ekat::subview(linoz_dPmL_dT, icol);
          linoz_data.linoz_dPmL_dO3col_icol =
            ekat::subview(linoz_dPmL_dO3col, icol);
          linoz_data.linoz_cariolle_pscs_icol =
            ekat::subview(linoz_cariolle_pscs, icol);
      }
      const auto& vmr_icol = ekat::subview(vmr,icol);

      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
      //-----------------
      // LINOZ chemistry
      //-----------------
      const Real temp = atm.temperature(kk);
      const Real pmid = atm.pressure(kk);
      const Real pdel = atm.hydrostatic_dp(kk);

      // the following things are diagnostics, which we're not
      // including in the first rev
      Real do3_linoz = 0, do3_linoz_psc = 0, ss_o3 = 0, o3col_du_diag = 0,
           o3clim_linoz_diag = 0, zenith_angle_degrees = 0;


      const auto& vmr_kk = ekat::subview(vmr_icol,kk);

      // index of "O3" in solsym array (in EAM)
      mam4::lin_strat_chem::lin_strat_chem_solve_kk(
          // in
          o3_col_dens_i(kk), temp, zenith_angle(icol), pmid, dt, rlats,
          linoz_data.linoz_o3_clim_icol(kk), linoz_data.linoz_t_clim_icol(kk),
          linoz_data.linoz_o3col_clim_icol(kk),
          linoz_data.linoz_PmL_clim_icol(kk),
          linoz_data.linoz_dPmL_dO3_icol(kk), linoz_data.linoz_dPmL_dT_icol(kk),
          linoz_data.linoz_dPmL_dO3col_icol(kk),
          linoz_data.linoz_cariolle_pscs_icol(kk), linoz_conf.chlorine_loading,
          linoz_conf.psc_T,
          // out
          vmr_kk[o3_ndx],
          // outputs that are not used
          do3_linoz, do3_linoz_psc, ss_o3, o3col_du_diag, o3clim_linoz_diag,
          zenith_angle_degrees);

      // Update source terms above the ozone decay threshold
      if (kk >= nlev - linoz_conf.o3_lbl) {
        const Real o3l_vmr_old = vmr_kk(o3_ndx);
        Real do3mass = 0;
        const Real o3l_vmr_new =
            mam4::lin_strat_chem::lin_strat_sfcsink_kk(dt, pdel,          // in
                                                       o3l_vmr_old,       // in
                                                       linoz_conf.o3_sfc, // in
                                                       linoz_conf.o3_tau, // in
                                                       do3mass);          // out
        // Update the mixing ratio (vmr) for O3
        vmr_kk(o3_ndx) = o3l_vmr_new;
      }
        });
        });

    }
    // linoz ends
    Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::inject_to_progs", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const auto& vmr_icol = ekat::subview(vmr,icol);
      const auto& qq_icol = ekat::subview(qq,icol);
      const auto& state_q_icol = ekat::subview(state_q,icol);
      const auto& qqcw_pcnst_icol = ekat::subview(qqcw_pcnst,icol);
      // fetch column-specific subviews into aerosol prognostics
      mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);
     Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
      // Check for negative values and reset to zero
    for (int i = 0; i < num_gas_aerosol_constituents; ++i) {
      if (vmr(icol,kk,i) < 0.0)
        vmr(icol,kk,i) = 0.0;
    }
    const auto& vmr_kk = ekat::subview(vmr_icol,kk);
    const auto& qq_kk = ekat::subview(qq_icol,kk);

    mam4::microphysics::vmr2mmr(vmr_kk.data(), adv_mass_kg_per_moles, qq_kk.data());
    for (int i = offset_aerosol; i < pcnst; ++i) {
      state_q(icol, kk,i) = qq(icol, kk,i - offset_aerosol);
      qqcw_pcnst(icol, kk,i) = qqcw(icol, kk,i - offset_aerosol);
    }
    const auto& state_q_kk = ekat::subview(state_q_icol,kk);
    const auto& qqcw_pcnst_kk = ekat::subview(qqcw_pcnst_icol,kk);
    mam4::utils::inject_stateq_to_prognostics(state_q_kk, progs, kk);
    mam4::utils::inject_qqcw_to_prognostics(qqcw_pcnst_kk, progs, kk);
    });
    });

    // diagnostics
    // - dvmr/dt: Tendencies for mixing ratios  [kg/kg/s]
    view_2d dqdt_so4_aqueous_chemistry, dqdt_h2so4_uptake;
    view_3d aqso4_incloud_mmr_tendency, aqh2so4_incloud_mmr_tendency;
    if (extra_mam4_aero_microphys_diags_) {
      dqdt_so4_aqueous_chemistry = get_field_out("dqdt_so4_aqueous_chemistry").get_view<Real **>();
      dqdt_h2so4_uptake = get_field_out("dqdt_h2so4_uptake").get_view<Real **>();
      aqso4_incloud_mmr_tendency   = get_field_out("mam4_microphysics_tendency_aqso4").get_view<Real ***>();
      aqh2so4_incloud_mmr_tendency = get_field_out("mam4_microphysics_tendency_aqh2so4").get_view<Real ***>();
    }

    Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::diagnostics", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
      const auto dqdt_aqso4_icol = ekat::subview(dqdt_aqso4,icol);
      const auto dqdt_aqh2so4_icol =ekat::subview(dqdt_aqh2so4,icol);

      mam4::MicrophysDiagnosticArrays diag_arrays;
      if (extra_mam4_aero_microphys_diags) {
          diag_arrays.dqdt_so4_aqueous_chemistry = ekat::subview(dqdt_so4_aqueous_chemistry, icol);
          diag_arrays.dqdt_h2so4_uptake          = ekat::subview(dqdt_h2so4_uptake, icol);
          diag_arrays.aqh2so4_incloud_mmr_tendency = ekat::subview(aqh2so4_incloud_mmr_tendency, icol);
          diag_arrays.aqso4_incloud_mmr_tendency = ekat::subview(aqso4_incloud_mmr_tendency, icol);
      }

      // Diagnose the column-integrated flux (kg/m2/s) using
      // volume mixing ratios ( // kmol/kmol(air) )
      const auto &pdel = atm.hydrostatic_dp; // layer thickness (Pa)
      for (int m = 0; m < nmodes; ++m) {
      const int ll = config_setsox.lptr_so4_cw_amode[m] - offset_aerosol;
      if (0 <= ll) {
        const auto adv_mass = adv_mass_kg_per_moles[ll];
        Real vmr_so4 = 0.0;
        const Real gravit = mam4::Constants::gravity;
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team, nlev),
          [&](int kk, Real &lsum) { lsum += dqdt_aqso4_icol(kk,ll) * pdel(kk) / gravit; },
          vmr_so4);
        Real vmr_h2s = 0.0;
        Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(team, nlev),
          [&](int kk, Real &lsum) { lsum += dqdt_aqh2so4_icol(kk,ll) * pdel(kk) / gravit; },
          vmr_h2s);


      if (dqdt_so4_aqueous_chemistry.size())
        diag_arrays.dqdt_so4_aqueous_chemistry(m) =
            mam4::conversions::mmr_from_vmr(vmr_so4, adv_mass);
      if (dqdt_h2so4_uptake.size())
        diag_arrays.dqdt_h2so4_uptake(m) = mam4::conversions::mmr_from_vmr(vmr_h2s, adv_mass);
      if (aqso4_incloud_mmr_tendency.size()) {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev), [&](const int kk) {
              diag_arrays.aqso4_incloud_mmr_tendency(m, kk) =
                  mam4::conversions::mmr_from_vmr(dqdt_aqso4_icol(kk,ll), adv_mass);
            });
      }
      if (aqh2so4_incloud_mmr_tendency.size()) {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev), [&](const int kk) {
              diag_arrays.aqh2so4_incloud_mmr_tendency(m, kk) =
                  mam4::conversions::mmr_from_vmr(dqdt_aqh2so4_icol(kk,ll), adv_mass);
            });
      }
    } // (if 0 <= ll)
  } // for loop over num_modes


    });


}//run_small_kernels_microphysics

}  // namespace scream
