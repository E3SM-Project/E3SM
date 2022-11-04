#include "diagnostics/vertically_remapped_field.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

#include <iostream>
#include <fstream>

namespace scream
{

// =========================================================================================
VerticallyRemappedField::VerticallyRemappedField (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
  , m_field_layout(m_params.get<FieldLayout>("Field Layout"))
  , m_field_units(m_params.get<ekat::units::Units>("Field Units")) 
  , m_field_name(m_params.get<std::string>("Field Name"))
  , m_tgt_pres_levs(m_params.get<view_1d_const>("press_levels"))
  , m_tgt_num_levs(m_params.get<int>("tgt_num_levs"))
{
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (ekat::contains(std::vector<FieldTag>{LEV,ILEV},m_field_layout.tags().back()),
      "Error! VerticallyRemappedField diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");
}

// =========================================================================================
void VerticallyRemappedField::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  const auto& gname  = m_params.get<std::string>("Grid Name");
  auto m_grid = grids_manager->get_grid(gname);
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);

  m_pres_name = m_field_layout.tags().back()==LEV ? "p_mid" : "p_int";
  add_field<Required>(m_pres_name, m_field_layout, Pa, gname);
  m_num_levs = m_field_layout.dims().back();

  constexpr int ps = Pack::n;
  FieldLayout diag_layout { {COL,LEV}, {m_num_cols,m_tgt_num_levs} };
  FieldIdentifier fid (name(),diag_layout, m, gname);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void VerticallyRemappedField::compute_diagnostic_impl()
{
  using namespace scream::vinterp;

  //This is 2D source pressure
  const Field& pressure_f = get_field_in(m_pres_name);
  const auto pressure = pressure_f.get_view<const Pack**>();
 
  //input field
  const Field& f = get_field_in(m_field_name);
  const auto f_data_src = f.get_view<const Pack**>();
 
  //output field on new grid
  auto d_data_tgt = m_diagnostic_output.get_view<Pack**>();
  perform_vertical_interpolation(pressure,
                                 m_tgt_pres_levs,
                                 f_data_src, 
                                 d_data_tgt,
                                 m_num_levs,
                                 m_tgt_num_levs);

}

} //namespace scream
