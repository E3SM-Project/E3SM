#include "diagnostics/field_at_single_pressure.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

namespace scream
{

// =========================================================================================
FieldAtSinglePressure::FieldAtSinglePressure (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
  , m_field_layout(m_params.get<FieldLayout>("Field Layout"))
  , m_field_units(m_params.get<ekat::units::Units>("Field Units"))
{
  m_field_name  = m_params.get<std::string>("Field Name");

  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (ekat::contains(std::vector<FieldTag>{LEV,ILEV},m_field_layout.tags().back()),
      "Error! FieldAtSinglePressure diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");
}

// =========================================================================================
void FieldAtSinglePressure::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  const auto& gname  = m_params.get<std::string>("Grid Name");
  auto m_grid = grids_manager->get_grid(gname);
  const auto& grid_name = m_grid->name();
  int ncol = m_grid->get_num_local_dofs();
  int nlev = m_grid->get_num_vertical_levels();

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);
  if (ekat::contains(std::vector<FieldTag>{LEV},m_field_layout.tags().back())) {
    FieldLayout pres_layout { {COL,LEV}, {ncol,nlev} };
    m_pres_name = "p_mid";
    add_field<Required>(m_pres_name, pres_layout, Pa, gname);
  } else {
    FieldLayout pres_layout { {COL,ILEV}, {ncol,nlev+1} };
    m_pres_name = "p_int";
    add_field<Required>(m_pres_name, pres_layout, Pa, gname);
  }
}
// =========================================================================================
void FieldAtSinglePressure::compute_diagnostic_impl()
{
  const auto& pressure = get_field_in(m_pres_name);
  auto p_data = pressure.get_view<const Real**>();

  const auto& f = get_field_in(m_field_name);
  auto f_data = f.get_internal_view_data<const Real>();

  auto d_data = m_diagnostic_output.get_internal_view_data<Real>();
  auto f_size = f.get_header().get_alloc_properties().get_num_scalars();
  auto stride = f.get_header().get_alloc_properties().get_last_extent();
  auto d_size = f_size/stride;
  auto pres_level  = m_pressure_level;

  // Here we use vertical interpolation to interpolate the field for each column onto a 
  // particular pressure coordinate.
  Kokkos::parallel_for(m_diagnostic_output.name(),
                       Kokkos::RangePolicy<>(0,d_size),
                       KOKKOS_LAMBDA(const int& idx) {
//ASD <this is the old code from FieldAtLevel>    d_data[idx] = f_data[idx*stride + pressure];  // Need to insert an interpolation here.
  });
  Kokkos::fence();
}

void FieldAtSinglePressure::set_required_field_impl (const Field& f)
{
  const auto& f_fid = f.get_header().get_identifier();

  // Now that we have the exact units of f, we can build the diagnostic field
  const auto& pres_str = m_params.get<std::string>("Field Pressure");
  auto is_int = [](const std::string& s) -> bool {
    return s.find_first_not_of("0123456789")==std::string::npos;
  };
  EKAT_REQUIRE_MSG (is_int(pres_str),
      "Error! Entry 'Field Pressure' must be a string representation of an integer.\n");

  m_pressure_level  = std::stoi(pres_str);
  EKAT_REQUIRE_MSG (m_pressure_level>=0,
      "Error! Invalid value for 'Field Pressure' in FieldAtSinglePressure diagnostic.\n"
      "  - field name  : " + f_fid.name() + "\n"
      "  - field layout: " + to_string(f_fid.get_layout()) + "\n"
      "  - field pressure : " + std::to_string(m_pressure_level) + "\n");

  FieldIdentifier d_fid (f_fid.name() + "_" + std::to_string(m_pressure_level) + "hPa",
                         m_field_layout.strip_dim(m_field_layout.rank()-1),
                         f_fid.get_units(),
                         f_fid.get_grid_name());

  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

} //namespace scream
