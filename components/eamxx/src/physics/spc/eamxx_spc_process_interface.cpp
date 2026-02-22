#include "eamxx_spc_process_interface.hpp"

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

SPC::SPC (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  EKAT_REQUIRE_MSG(m_params.isParameter("spc_data_file"),
      "ERROR: spc_data_file is missing from SPC parameter list.");
}

// =========================================================================================
void SPC::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  constexpr auto nondim = Units::nondimensional();
  constexpr int ps = SCREAM_PACK_SIZE;

  m_model_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_model_grid->name();
  
  // Define the different field layouts that will be used for this process
  auto scalar3d_mid    = m_model_grid->get_3d_scalar_layout(true);

  // Set of fields used strictly as input
  add_field<Required>("p_mid"      , scalar3d_mid, Pa,     grid_name, ps);

  // Set of fields used strictly as output
  add_field<Computed>("o3_volume_mix_ratio", scalar3d_mid,    mol/mol,   grid_name, ps);
}

// =========================================================================================
void SPC::initialize_impl (const RunType /* run_type */)
{
  using namespace ekat::units;

  // NOTE: SPC does not have an internal persistent state, so run_type is irrelevant

  std::vector<Field> spc_fields = {
      get_field_out("o3_volume_mix_ratio").alias("O3")
  };
  auto spc_data_file = m_params.get<std::string>("spc_data_file");
  auto spc_map_file  = m_params.get<std::string>("spc_remap_file","");

  auto pmid = get_field_in("p_mid");

  util::TimeStamp ref_ts (1,1,1,0,0,0); // Beg of any year, since we use yearly periodic timeline
  m_data_interpolation = std::make_shared<DataInterpolation>(m_model_grid,spc_fields);
  m_data_interpolation->setup_time_database ({spc_data_file},util::TimeLine::YearlyPeriodic, DataInterpolation::Linear, ref_ts);

  if (m_iop_data_manager!=nullptr) {
    // IOP cases cannot have a remap file. We will create a IOPRemapper as the horiz remapper
    EKAT_REQUIRE_MSG(spc_map_file == "" or spc_map_file=="none",
      "Error! Cannot define spc_remap_file for cases with an Intensive Observation Period defined. "
      "The IOP class defines it's own remap from file data -> model data.\n");

    // TODO: expose tgt lat/lon in IOPDataManager, to avoid injecting knowledge
    // of its param list structure in other places
    Real iop_lat = m_iop_data_manager->get_params().get<Real>("target_latitude");
    Real iop_lon = m_iop_data_manager->get_params().get<Real>("target_longitude");
    m_data_interpolation->create_horiz_remappers (iop_lat,iop_lon);
  } else {
    m_data_interpolation->create_horiz_remappers (spc_map_file=="none" ? "" : spc_map_file);
  }
  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = DataInterpolation::Dynamic3DRef;
  vremap_data.pname = "PS";
  vremap_data.pmid = pmid;
  m_data_interpolation->create_vert_remapper (vremap_data);
  m_data_interpolation->init_data_interval (start_of_step_ts());

  // Set property checks for fields in this process
  using FWI = FieldWithinIntervalCheck;
  const auto eps = std::numeric_limits<double>::epsilon();

  add_postcondition_check<FWI>(get_field_out("o3_volume_mix_ratio"),m_model_grid,1e-36,1e-2,true);
}

// =========================================================================================
void SPC::run_impl (const double /* dt */)
{
  m_data_interpolation->run(end_of_step_ts());
}

} // namespace scream
