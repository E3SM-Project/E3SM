#include "diagnostics/field_at_pressure_level.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

namespace scream
{

// =========================================================================================
FieldAtPressureLevel::FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
  ,m_field_layout(m_params.get<FieldLayout>("Field Layout"))
  ,m_field_units(m_params.get<ekat::units::Units>("Field Units"))
  ,m_field_name(m_params.get<std::string>("Field Name"))
  ,m_pressure_level(m_params.get<Real>("Field Target Pressure"))
{
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (ekat::contains(std::vector<FieldTag>{LEV,ILEV},m_field_layout.tags().back()),
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");

  m_p_tgt = view_1d<mPack>("",1);
  Kokkos::deep_copy(m_p_tgt, m_pressure_level);

  m_mask_val = m_params.get<Real>("mask_value",Real(-99999));
}

// =========================================================================================
void FieldAtPressureLevel::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  const auto& gname  = m_params.get<std::string>("Grid Name");
  auto m_grid = grids_manager->get_grid(gname);
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();
  auto num_levs = m_grid->get_num_vertical_levels();

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);

  m_pres_name = m_field_layout.tags().back()==LEV ? "p_mid" : "p_int";
  add_field<Required>(m_pres_name, m_field_layout, Pa, gname);
  m_num_levs = m_field_layout.dims().back();

  FieldLayout diag_layout { {COL}, {m_num_cols} };
  FieldIdentifier fid (name(),diag_layout, m_field_units, gname);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void FieldAtPressureLevel::compute_diagnostic_impl()
{
  using namespace scream::vinterp;

  //This is 2D source pressure
  const Field& pressure_f = get_field_in(m_pres_name);
  const auto pressure = pressure_f.get_view<const mPack**>();

  //input field
  const Field& f = get_field_in(m_field_name);
  const auto f_data_src = f.get_view<const mPack**>();

  //output field on new grid
  auto d_data_tgt = m_diagnostic_output.get_view<Real*>();
  view_2d<mPack> data_tgt_tmp(reinterpret_cast<mPack*>(d_data_tgt.data()),m_num_cols,1);  // Note, vertical interp wants a 2D view, so we create a temporary one

  perform_vertical_interpolation(pressure,m_p_tgt,f_data_src,data_tgt_tmp,m_num_levs,1,m_mask_val);

}

} //namespace scream
