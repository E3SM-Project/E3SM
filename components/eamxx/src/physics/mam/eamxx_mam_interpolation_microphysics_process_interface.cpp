#include <mam4xx/mam4.hpp>
#include <physics/mam/eamxx_mam_interpolation_microphysics_process_interface.hpp>
// impl namespace for some driver level functions for microphysics

#include "share/util/eamxx_data_interpolation.hpp"
#include "readfiles/tracer_reader_utils.hpp"

namespace scream {

MAMInterpolationMicrophysics::MAMInterpolationMicrophysics(const ekat::Comm &comm,
                                 const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}
// ================================================================
//  SET_GRIDS
// ================================================================
void MAMInterpolationMicrophysics::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  // set grid for all the inputs and outputs
  // use physics grid
  grid_                 = grids_manager->get_grid("physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // number of levels per column
  const FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  using namespace ekat::units;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  add_tracers_wet_atm();
  add_fields_dry_atm();
  // cloud liquid number mixing ratio [1/kg]
  auto n_unit           = 1 / kg;   // units of number mixing ratios of tracers
  add_tracer<Required>("nc", grid_, n_unit);


  m_var_names_oxi = {"O3", "OH", "NO3", "HO2"};
  for(const auto &field_name : m_var_names_oxi) {
      add_field<Computed>(field_name, scalar3d_mid, nondim, grid_name);
  }

  // linoz
  m_var_names_linoz = {
        "o3_clim",  "o3col_clim", "t_clim",      "PmL_clim",
        "dPmL_dO3", "dPmL_dT",    "dPmL_dO3col", "cariolle_pscs"};

  for(const auto &field_name : m_var_names_linoz) {
      add_field<Computed>(field_name, scalar3d_mid, nondim, grid_name);
  }

  // This order corresponds to files in namelist e3smv2
  m_elevated_emis_var_names["so2"]    = {"BB", "ENE_ELEV", "IND_ELEV",
                                          "contvolc"};
  m_elevated_emis_var_names["so4_a1"] = {"BB", "ENE_ELEV", "IND_ELEV",
                                          "contvolc"};
  m_elevated_emis_var_names["so4_a2"] = {"contvolc"};
  m_elevated_emis_var_names["pom_a4"] = {"BB"};
  m_elevated_emis_var_names["bc_a4"]  = {"BB"};
  m_elevated_emis_var_names["num_a1"] = {
        "num_a1_SO4_ELEV_BB", "num_a1_SO4_ELEV_ENE", "num_a1_SO4_ELEV_IND",
        "num_a1_SO4_ELEV_contvolc"};
  m_elevated_emis_var_names["num_a2"] = {"num_a2_SO4_ELEV_contvolc"};
  // num_a4
  // FIXME: why the sectors in this files are num_a1;
  //  I guess this should be num_a4? Is this a bug in the orginal nc files?
  m_elevated_emis_var_names["num_a4"] = {"num_a1_BC_ELEV_BB",
                                          "num_a1_POM_ELEV_BB"};
  m_elevated_emis_var_names["soag"] = {"SOAbb_src", "SOAbg_src", "SOAff_src"};

  for(const auto &pair : m_elevated_emis_var_names) {
    const auto &var_name =pair.first;
    for(const auto &field_name : pair.second) {
      add_field<Computed>(field_name+"_"+var_name, scalar3d_mid, nondim, grid_name);
    }
  }

}  // set_grids

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by
// the above. Buffer type given the number of columns and vertical
// levels
size_t MAMInterpolationMicrophysics::requested_buffer_size_in_bytes() const {
return mam_coupling::buffer_size(ncol_, nlev_, 0, 0);;
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to
// store intermediate (dry) quantities on the given number of
// columns with the given number of vertical levels. Returns the
// number of bytes allocated.

void MAMInterpolationMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMMicrophysics."
                   " Used memory: "
                       << std::to_string(used_mem)
                       << "."
                          " Requested memory: "
                       << std::to_string(requested_buffer_size_in_bytes())
                       << ". \n");
}
// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMInterpolationMicrophysics::initialize_impl(const RunType run_type) {

  std::vector<Field> elevated_fields;

  const auto oxid_file_name = m_params.get<std::string>("mam4_oxid_file_name");
  const std::string oxid_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");

  for(const auto &field_name : m_var_names_oxi) {
      elevated_fields.push_back(get_field_out(field_name));
  }

  auto pmid = get_field_in("p_mid");

  // in format YYYYMMDD
  const int oxid_ymd = m_params.get<int>("mam4_oxid_ymd");
  util::TimeStamp ref_ts_oxid= mam_coupling::convert_date(oxid_ymd);
  m_data_interpolation = std::make_shared<DataInterpolation>(grid_,elevated_fields);
  m_data_interpolation->setup_time_database ({oxid_file_name},util::TimeLine::YearlyPeriodic, ref_ts_oxid);

  const std::string extfrc_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
  DataInterpolation::RemapData remap_data;
  remap_data.hremap_file = extfrc_map_file=="none" ? "" : extfrc_map_file;
  remap_data.vr_type = DataInterpolation::MAM4_PSRef;
  remap_data.pname = "PS";
  remap_data.pmid = pmid;
  m_data_interpolation->setup_remappers (remap_data);
  m_data_interpolation->init_data_interval (start_of_step_ts());
  // linoz

  // in format YYYYMMDD
  const int linoz_cyclical_ymd = m_params.get<int>("mam4_linoz_ymd");
  util::TimeStamp ref_ts_linoz = mam_coupling::convert_date(linoz_cyclical_ymd);
  // util::TimeStamp ref_ts_linoz (1,1,1,0,0,0); // Beg of any year, since we use yearly periodic timeline
  const auto m_linoz_file_name = m_params.get<std::string>("mam4_linoz_file_name");
  const std::string linoz_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");

  std::vector<Field> linoz_fields;
  for(const auto &field_name : m_var_names_linoz) {
      linoz_fields.push_back(get_field_out(field_name));
  }

  m_data_interpolation_linoz = std::make_shared<DataInterpolation>(grid_,linoz_fields);
  m_data_interpolation_linoz->setup_time_database ({m_linoz_file_name},util::TimeLine::YearlyPeriodic, ref_ts_linoz);
  DataInterpolation::RemapData remap_data_linoz;
  remap_data_linoz.hremap_file = linoz_map_file=="none" ? "" : linoz_map_file;
  remap_data_linoz.vr_type = DataInterpolation::MAM4_ZONAL;
  remap_data_linoz.pname = "lev";
  remap_data_linoz.pmid = pmid;
  m_data_interpolation_linoz->setup_remappers (remap_data_linoz);
  m_data_interpolation_linoz->init_data_interval (start_of_step_ts());

  // Populate the wet atmosphere state with views from fields
  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  //FIXME: units
  Field z_iface(FieldIdentifier("z_iface",grid_->get_3d_scalar_layout(false),ekat::units::Pa,grid_->name()));
  // z_iface.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  z_iface.allocate_view();
  const auto z_iface_view = z_iface.get_view<  Real **>();
  dry_atm_.z_iface=z_iface_view;

  // std::string var_name="so2";
  int elevated_emiss_cyclical_ymd = m_params.get<int>("elevated_emiss_ymd");

  for(const auto &pair : m_elevated_emis_var_names) {
    const auto& var_name=pair.first;
    std::string item_name = "mam4_" + var_name + "_elevated_emiss_file_name";
    const auto file_name  = m_params.get<std::string>(item_name);
    util::TimeStamp ref_ts_vertical = mam_coupling::convert_date(elevated_emiss_cyclical_ymd);
    std::vector<Field> vertical_fields;
    for(const auto &field_name :pair.second) {
      vertical_fields.push_back(get_field_out(field_name+"_"+var_name).alias(field_name));
    }
    std::shared_ptr<DataInterpolation> di_vertical = std::make_shared<DataInterpolation>(grid_,vertical_fields);
    di_vertical->setup_time_database ({file_name},util::TimeLine::YearlyPeriodic, ref_ts_vertical);
    DataInterpolation::RemapData remap_data_vertical;
    // FIXME linoz_map_file
    remap_data_vertical.hremap_file = linoz_map_file=="none" ? "" : linoz_map_file;
    remap_data_vertical.vr_type = DataInterpolation::MAM4_ELEVATED_EMISSIONS;
    remap_data_vertical.pname = "altitude_int";
    remap_data_vertical.pmid = z_iface;
    di_vertical->setup_remappers (remap_data_vertical);
    di_vertical->init_data_interval (start_of_step_ts());
    m_data_interpolation_vertical.push_back(di_vertical);
  }//end var_name


}  // initialize_impl

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMInterpolationMicrophysics::run_impl(const double dt) {
  m_data_interpolation->run(end_of_step_ts());
  m_data_interpolation_linoz->run(end_of_step_ts());

  const auto &wet_atm = wet_atm_;
  const auto &dry_atm = dry_atm_;

  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  Kokkos::parallel_for(
      scan_policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int i = team.league_rank();  // column index
        mam_coupling::compute_dry_mixing_ratios(team, wet_atm, dry_atm, i);
        team.team_barrier();
        // vertical heights has to be computed after computing dry mixing ratios
        // for atmosphere
        mam_coupling::compute_vertical_layer_heights(team, dry_atm, i);
      });

  for (size_t i = 0; i < m_elevated_emis_var_names.size(); ++i) {
    m_data_interpolation_vertical[i]->run(end_of_step_ts());
  }

}  // MAMInterpolationMicrophysics::run_impl

}  // namespace scream
