#include "eamxx_spc_process_interface.hpp"

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>
#include <ekat_string_utils.hpp>

namespace scream
{

SPC::SPC (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  EKAT_REQUIRE_MSG(m_params.isParameter("spc_data_file"),
      "ERROR: spc_data_file is missing from SPC parameter list.");
  EKAT_REQUIRE_MSG(m_params.isParameter("gas_species"),
      "ERROR: gas_species is missing from SPC parameter list.");

  m_gas_species = m_params.get<std::vector<std::string>>("gas_species");
  EKAT_REQUIRE_MSG(m_gas_species.size()>0,
      "ERROR: gas_species list in SPC parameter list is empty.");
}

// =========================================================================================
void SPC::create_requests()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  constexpr int ps = SCREAM_PACK_SIZE;

  m_model_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_model_grid->name();
  
  // Define the different field layouts that will be used for this process
  auto scalar3d_mid = m_model_grid->get_3d_scalar_layout(LEV);

  // Set of fields used strictly as input
  add_field<Required>("p_mid"      , scalar3d_mid, Pa,     grid_name, ps);

  // Set of fields used strictly as output
  for (const auto& species : m_gas_species) {
    add_field<Computed>(species + "_volume_mix_ratio", scalar3d_mid, mol/mol, grid_name, ps);
  }
}

// =========================================================================================
void SPC::initialize_impl (const RunType /* run_type */)
{
  using namespace ekat::units;

  // NOTE: SPC does not have an internal persistent state, so run_type is irrelevant

  std::vector<Field> spc_fields;
  spc_fields.reserve(m_gas_species.size());
  for (const auto& species : m_gas_species) {
    spc_fields.push_back(get_field_out(species + "_volume_mix_ratio").alias(ekat::upper_case(species)));
  }
  auto spc_data_file = m_params.get<std::string>("spc_data_file");
  auto spc_map_file  = m_params.get<std::string>("spc_remap_file","");
  auto time_interpolation_method = m_params.get<std::string>("time_interpolation_method","yearly_periodic");

  auto pmid = get_field_in("p_mid");
  
  m_data_interpolation = std::make_shared<DataInterpolation>(m_model_grid,spc_fields);
  m_data_interpolation->setup_time_database({spc_data_file}, time_interpolation_method);
  m_data_interpolation->create_horiz_remappers(spc_map_file, m_iop_data_manager);
  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = DataInterpolation::Dynamic3DRef;
  vremap_data.pname = "PS";
  vremap_data.pmid = pmid;
  m_data_interpolation->create_vert_remapper (vremap_data);
  m_data_interpolation->init_time_interpolation (start_of_step_ts(), DataInterpolation::Linear);

  // Set property checks for fields in this process
  using FWI = FieldWithinIntervalCheck;

  for (const auto& species : m_gas_species) {
    add_postcondition_check<FWI>(get_field_out(species + "_volume_mix_ratio"),m_model_grid,1e-36,1e-2,true);
  }
}

// =========================================================================================
void SPC::run_impl (const double /* dt */)
{
  m_data_interpolation->run(end_of_step_ts());
}

} // namespace scream
