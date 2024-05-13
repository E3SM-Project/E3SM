#include "diagnostics/field_at_pressure_level.hpp"
#include "share/util/scream_vertical_interpolation.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

namespace scream
{

// =========================================================================================
FieldAtPressureLevel::
FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  m_field_name = m_params.get<std::string>("field_name");

  // Figure out the pressure value
  const auto& location = m_params.get<std::string>("vertical_location");
  auto chars_start = location.find_first_not_of("0123456789.");
  EKAT_REQUIRE_MSG (chars_start!=0 && chars_start!=std::string::npos,
      "Error! Invalid string for pressure value for FieldAtPressureLevel.\n"
      " - input string   : " + location + "\n"
      " - expected format: Nxyz, with N integer, and xyz='mb', 'hPa', or 'Pa'\n");
  const auto press_str = location.substr(0,chars_start);
  m_pressure_level = std::stod(press_str);

  const auto units = location.substr(chars_start);
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid string for pressure value for FieldAtPressureLevel.\n"
      " - input string   : " + location + "\n"
      " - expected format: Nxyz, with N integer, and xyz='mb', 'hPa', or 'Pa'\n");

  // Convert pressure level to Pa, the units of pressure in the simulation
  if (units=="mb" || units=="hPa") {
    m_pressure_level *= 100;
  }

  m_p_tgt = view_1d<Pack1>("",1);
  Kokkos::deep_copy(m_p_tgt, m_pressure_level);

  m_mask_val = m_params.get<double>("mask_value",Real(std::numeric_limits<float>::max()/10.0));

  m_diag_name = m_field_name + "_at_" + location;
}

void FieldAtPressureLevel::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  const auto& gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_field_name,gname);

  // We don't know yet which one we need
  add_field<Required>("p_mid",gname);
  add_field<Required>("p_int",gname);
}

void FieldAtPressureLevel::
initialize_impl (const RunType /*run_type*/)
{
  const auto& f = get_field_in(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=3,
      "Error! Field rank not supported by FieldAtPressureLevel.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  FieldIdentifier d_fid (m_diag_name,layout.clone().strip_dim(tag),fid.get_units(),fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  m_pressure_name = tag==LEV ? "p_mid" : "p_int";
  m_num_levs = layout.dims().back();
  auto num_cols = layout.dims().front();

  // Take care of mask tracking for this field, in case it is needed.  This has two steps:
  //   1.  We need to actually track the masked columns, so we create a 2d (COL only) field.
  //       NOTE: Here we assume that even a source field of rank 3+ will be masked the same
  //       across all components so the mask is represented by a column-wise slice.
  //   2.  We also need to create a helper field that may be used w/ the vertical remapper to set
  //       the mask.  This field is 3d (COLxLEV) to mimic the type of interpolation that is
  //       being conducted on the source field.

  // Add a field representing the mask as extra data to the diagnostic field.
  auto nondim = ekat::units::Units::nondimensional();
  const auto& gname = fid.get_grid_name();

  std::string mask_name = name() + " mask";
  FieldLayout mask_layout( {COL}, {num_cols});
  FieldIdentifier mask_fid (mask_name,mask_layout, nondim, gname);
  Field diag_mask(mask_fid);
  diag_mask.allocate_view();
  m_diagnostic_output.get_header().set_extra_data("mask_data",diag_mask);
  m_diagnostic_output.get_header().set_extra_data("mask_value",m_mask_val);

  // Allocate helper views
  FieldLayout mask_src_layout( {COL, LEV}, {num_cols, m_num_levs});
  FieldIdentifier mask_src_fid ("mask_tmp",mask_src_layout, nondim, gname);
  m_mask_field = Field(mask_src_fid);
  m_mask_field.allocate_view();

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src = get_fields_in().front();
  const auto& src_atts = src.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

// =========================================================================================
void FieldAtPressureLevel::compute_diagnostic_impl()
{
  using namespace scream::vinterp;

  //This is 2D source pressure
  const Field& pressure_f = get_field_in(m_pressure_name);
  const auto pressure = pressure_f.get_view<const Pack1**>();
  view_Nd<const Pack1,2> pres(pressure.data(),pressure.extent_int(0),pressure.extent_int(1));
  // Kokkos::deep_copy(pres,pressure);

  const Field& f = get_field_in(m_field_name);

  // The setup for interpolation varies depending on the rank of the input field:
  const int rank = f.rank();

  m_mask_field.deep_copy(1.0);
  auto mask_v_tmp = m_mask_field.get_view<Pack1**>();
  if (rank==2) {
    const auto f_data_src = f.get_view<const Pack1**>();
    //output field on new grid
    auto d_data_tgt = m_diagnostic_output.get_view<Pack1*>();
    view_Nd<Pack1,2> data_tgt_tmp(d_data_tgt.data(),d_data_tgt.extent_int(0),1);  // Note, vertical interp wants a 2D view, so we create a temporary one
    perform_vertical_interpolation<Real,1,2>(pres,m_p_tgt,f_data_src,data_tgt_tmp,m_num_levs,1,m_mask_val);

    // Track mask
    auto mask = m_diagnostic_output.get_header().get_extra_data<Field>("mask_data");
    auto d_mask_tgt = mask.get_view<Pack1*>();
    view_Nd<Pack1,2> mask_tgt_tmp(d_mask_tgt.data(),d_mask_tgt.extent_int(0),1);  
    perform_vertical_interpolation<Real,1,2>(pres,m_p_tgt,mask_v_tmp,mask_tgt_tmp,m_num_levs,1,0);
  } else if (rank==3) {
    const auto f_data_src = f.get_view<const Pack1***>();
    //output field on new grid
    auto d_data_tgt = m_diagnostic_output.get_view<Pack1**>();
    view_Nd<Pack1,3> data_tgt_tmp(d_data_tgt.data(),d_data_tgt.extent_int(0),d_data_tgt.extent_int(1),1);  

    perform_vertical_interpolation<Real,1,3>(pres,m_p_tgt,f_data_src,data_tgt_tmp,m_num_levs,1,m_mask_val);

    // Track mask
    auto mask = m_diagnostic_output.get_header().get_extra_data<Field>("mask_data");
    auto d_mask_tgt = mask.get_view<Pack1*>();
    view_Nd<Pack1,2> mask_tgt_tmp(d_mask_tgt.data(),d_mask_tgt.extent_int(0),1);  
    perform_vertical_interpolation<Real,1,2>(pres,m_p_tgt,mask_v_tmp,mask_tgt_tmp,m_num_levs,1,0);
  } else {
    EKAT_ERROR_MSG("Error! field at pressure level only supports fields ranks 2 and 3 \n");
  }

}

} //namespace scream
