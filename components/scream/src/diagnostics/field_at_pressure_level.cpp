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

  constexpr int ps = Pack::n;

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);
  if (m_field_layout.tags().back()==LEV) {
    m_pres_name = "p_mid";
    m_num_levs = num_levs;
    FieldLayout pres_layout { {COL,LEV}, {m_num_cols,m_num_levs} };
    add_field<Required>(m_pres_name, pres_layout, Pa, gname, ps);
  } else {
    m_pres_name = "p_int";
    m_num_levs = num_levs+1;
    FieldLayout pres_layout { {COL,ILEV}, {m_num_cols,m_num_levs} };
    add_field<Required>(m_pres_name, pres_layout, Pa, gname, ps);
  }

  FieldLayout diag_layout { {COL}, {m_num_cols} };
  FieldIdentifier fid (name(),diag_layout, m, gname);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void FieldAtPressureLevel::compute_diagnostic_impl()
{
  using namespace scream::vinterp;

  //This is 2D source pressure
  const Field& pressure_f = get_field_in(m_pres_name);
  const auto pressure = pressure_f.get_view<const Pack**>();

  //This is the 1D target pressure
  view_1d<Pack> p_tgt = view_1d<Pack>("",1);  // We only plan to map onto a pressure level
  Kokkos::deep_copy(p_tgt, m_pressure_level);

  //input field
  const Field& f = get_field_in(m_field_name);
  const auto f_data_src = f.get_view<const Pack**>();

  //output field on new grid
  auto d_data_tgt = m_diagnostic_output.get_view<Real*>();
  view_2d<Pack> data_tgt_tmp("",m_num_cols,1);  // Note, vertical interp wants a 2D view, so we create a temporary one

  perform_vertical_interpolation(pressure,p_tgt,f_data_src,data_tgt_tmp,m_num_levs,1);
  Kokkos::parallel_for("", m_num_cols, KOKKOS_LAMBDA (const int& icol) {
    d_data_tgt(icol) = data_tgt_tmp(icol,0)[0];
  });

}

} //namespace scream
