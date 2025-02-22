#include "eamxx_spa_process_interface.hpp"

#include "share/util/eamxx_data_interpolation.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat/ekat_assert.hpp>
#include <ekat/util/ekat_units.hpp>

namespace scream
{

SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  EKAT_REQUIRE_MSG(m_params.isParameter("spa_data_file"),
      "ERROR: spa_data_file is missing from SPA parameter list.");
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  constexpr auto nondim = Units::nondimensional();
  constexpr int ps = SCREAM_PACK_SIZE;

  m_model_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_model_grid->name();

  // Get bands info from file, and log it
  const auto spa_data_file = m_params.get<std::string>("spa_data_file");
  const int nswbands = scorpio::get_dimlen(spa_data_file,"swband");
  const int nlwbands = scorpio::get_dimlen(spa_data_file,"lwband");
  this->log(LogLevel::info,
      "SPA data file bands dimensions:\n"
      "  - num sw bands: " + std::to_string(nswbands) + "\n"
      "  - num lw bands: " + std::to_string(nlwbands) + "\n");

  // Define the different field layouts that will be used for this process
  auto scalar3d_mid    = m_model_grid->get_3d_scalar_layout(true);
  auto scalar2d        = m_model_grid->get_2d_scalar_layout();
  auto scalar1d_mid    = m_model_grid->get_vertical_layout(true);
  auto scalar3d_swband = m_model_grid->get_3d_vector_layout(true,nswbands,"swband");
  auto scalar3d_lwband = m_model_grid->get_3d_vector_layout(true,nlwbands,"lwband");

  // Set of fields used strictly as input
  add_field<Required>("p_mid"      , scalar3d_mid, Pa,     grid_name, ps);

  // Set of fields used strictly as output
  add_field<Computed>("nccn",        scalar3d_mid,    1/kg,   grid_name, ps);
  add_field<Computed>("aero_g_sw",   scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_ssa_sw", scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_tau_sw", scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_tau_lw", scalar3d_lwband, nondim, grid_name, ps);
}

// =========================================================================================
void SPA::initialize_impl (const RunType /* run_type */)
{
  using namespace ekat::units;

  // NOTE: SPA does not have an internal persistent state, so run_type is irrelevant

  std::vector<Field> spa_fields = {
      get_field_out("nccn").alias("CCN3"),
      get_field_out("aero_g_sw").alias("AER_G_SW"),
      get_field_out("aero_ssa_sw").alias("AER_SSA_SW"),
      get_field_out("aero_tau_sw").alias("AER_TAU_SW"),
      get_field_out("aero_tau_lw").alias("AER_TAU_LW")
  };
  auto spa_data_file = m_params.get<std::string>("spa_data_file");
  auto spa_map_file  = m_params.get<std::string>("spa_remap_file","");

  // SPA doesn't really *need* pint, but DataInterpolation does. It's important to stress that
  // NO FIELD VALUES from p_int are accessed in the DataInterpolation we build, since we
  // don't remap any field on interfaces. But when the VerticalRemapper in the DataInterpolation
  // is setup, pint must store the correct layout of a field defined at interfaces; it does not
  // have to have COL dimension. It just need to have the ILEV tag in the layout AND have alloc
  // properties compatible with SCREAM_PACK_SIZE.
  // NOTE: we could just add p_int as a required field, but that would be misleading in the DAG
  auto pmid = get_field_in("p_mid");
  Field pint(FieldIdentifier("p_int",m_model_grid->get_vertical_layout(false),Pa,m_model_grid->name()));
  pint.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  pint.allocate_view();

  util::TimeStamp ref_ts (1,1,1,0,0,0); // Beg of any year, since we use yearly periodic timeline
  m_data_interpolation = std::make_shared<DataInterpolation>(m_model_grid,spa_fields);
  m_data_interpolation->setup_time_database ({spa_data_file},util::TimeLine::YearlyPeriodic, ref_ts);

  DataInterpolation::RemapData remap_data;
  remap_data.hremap_file = spa_map_file=="None" ? "" : spa_map_file;
  if (m_iop_data_manager!=nullptr) {
    // IOP cases cannot have a remap file. We will create a IOPRemapper as the horiz remapper
    EKAT_REQUIRE_MSG(spa_map_file == "" or spa_map_file=="None",
      "Error! Cannot define spa_remap_file for cases with an Intensive Observation Period defined. "
      "The IOP class defines it's own remap from file data -> model data.\n");

    // TODO: expose tgt lat/lon in IOPDataManager, to avoid injecting knowledge
    // of its param list structure in other places
    remap_data.iop_lat = m_iop_data_manager->get_params().get<Real>("target_latitude");
    remap_data.iop_lon = m_iop_data_manager->get_params().get<Real>("target_longitude");
    remap_data.has_iop = true;
  }
  remap_data.vr_type = DataInterpolation::Dynamic3DRef;
  remap_data.pname = "PS";
  remap_data.pmid = pmid;
  remap_data.pint = pint;
  m_data_interpolation->setup_remappers (remap_data);
  m_data_interpolation->init_data_interval (timestamp());

  // Set property checks for fields in this process
  using FWI = FieldWithinIntervalCheck;
  const auto eps = std::numeric_limits<double>::epsilon();

  add_postcondition_check<FWI>(get_field_out("nccn"),m_model_grid,0,1e11,true,-1e11*eps);

  // TODO: add an epslon to max possible upper bound of aero_ssa_sw?
  add_postcondition_check<FWI>(get_field_out("aero_g_sw"),m_model_grid,0.0,1.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_ssa_sw"),m_model_grid,0.0,1.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_tau_sw"),m_model_grid,0.0,1.0,true);
  add_postcondition_check<FWI>(get_field_out("aero_tau_lw"),m_model_grid,0.0,1.0,true);
}

// =========================================================================================
void SPA::run_impl (const double dt)
{
  m_data_interpolation->run(timestamp()+dt);
}

} // namespace scream
