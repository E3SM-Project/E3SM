#include "diagnostics/vertically_remapped_field.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

#include <iostream>
#include <fstream>

namespace scream
{

// =========================================================================================
VerticallyRemappedField::VerticallyRemappedField (const ekat::Comm& comm, 
                                            const ekat::ParameterList& params,
                                            const view_1d_const& m_tgt_pres_levs_, 
                                            int m_tgt_num_levs_)
  : AtmosphereDiagnostic(comm,params)
  , m_field_layout(m_params.get<FieldLayout>("Field Layout"))
  , m_field_units(m_params.get<ekat::units::Units>("Field Units")) 
  , m_field_name(m_params.get<std::string>("Field Name"))
  , m_tgt_pres_levs(m_tgt_pres_levs_)
  , m_tgt_num_levs(m_tgt_num_levs_)
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
  m_num_levs = m_grid->get_num_vertical_levels();

  constexpr int ps = Pack::n;
  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);
  //  if (ekat::contains(std::vector<FieldTag>{LEV},m_field_layout.tags().back())) {
  if (m_field_layout.tags().back()==LEV) {
    FieldLayout pres_layout { {COL,LEV}, {m_num_cols,m_num_levs} };
    m_pres_name = "p_mid";
    add_field<Required>(m_pres_name, pres_layout, Pa, gname);
  } else {
    FieldLayout pres_layout { {COL,ILEV}, {m_num_cols,m_num_levs+1} };
    m_pres_name = "p_int";
    add_field<Required>(m_pres_name, pres_layout, Pa, gname);
  }

  FieldLayout diag_layout { {COL,LEV}, {m_num_cols,m_tgt_num_levs} };
  FieldIdentifier fid (name(),diag_layout, m, gname);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  //C_ap.request_allocation(ps);
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
  perform_vertical_interpolation(pressure,m_tgt_pres_levs,f_data_src,d_data_tgt,m_num_levs,m_tgt_num_levs);

}

} //namespace scream
