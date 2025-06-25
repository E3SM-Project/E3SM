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
  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);

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


}  // set_grids

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by
// the above. Buffer type given the number of columns and vertical
// levels
size_t MAMInterpolationMicrophysics::requested_buffer_size_in_bytes() const {
return 0;
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to
// store intermediate (dry) quantities on the given number of
// columns with the given number of vertical levels. Returns the
// number of bytes allocated.

void MAMInterpolationMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
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
}  // initialize_impl

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMInterpolationMicrophysics::run_impl(const double dt) {
  m_data_interpolation->run(end_of_step_ts());
  m_data_interpolation_linoz->run(end_of_step_ts());
}  // MAMInterpolationMicrophysics::run_impl

}  // namespace scream
