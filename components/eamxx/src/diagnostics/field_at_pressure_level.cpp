#include "diagnostics/field_at_pressure_level.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

namespace scream
{

// =========================================================================================
FieldAtPressureLevel::
FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
 , m_field_name(m_params.get<std::string>("Field Name"))
 , m_field_layout(m_params.get<FieldLayout>("Field Layout"))
 , m_field_units(m_params.get<ekat::units::Units>("Field Units"))
 , m_pressure_level(m_params.get<Real>("Field Target Pressure"))
{
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (ekat::contains(std::vector<FieldTag>{LEV,ILEV},m_field_layout.tags().back()),
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");

  m_p_tgt = view_1d<mPack>("",1);
  Kokkos::deep_copy(m_p_tgt, m_pressure_level);

  m_mask_val = m_params.get<Real>("mask_value",Real(-999999.0));
}

// =========================================================================================
void FieldAtPressureLevel::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  const auto& gname  = m_params.get<std::string>("Grid Name");
  auto m_grid = grids_manager->get_grid(gname);
  m_num_cols = m_grid->get_num_local_dofs();

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);

  m_pres_name = m_field_layout.tags().back()==LEV ? "p_mid" : "p_int";
  add_field<Required>(m_pres_name, m_field_layout, Pa, gname);
  m_num_levs = m_field_layout.dims().back();

  // A field at a specific pressure level is just the same field w/ the LEV/ILEV dimension stripped from it.
  FieldLayout diag_layout = m_field_layout.strip_dim(m_field_layout.tags().back());
  FieldIdentifier fid (name(),diag_layout, m_field_units, gname);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();

  // Add a field representing the mask as extra data to the diagnostic field.
  Units nondim(0,0,0,0,0,0,0);
  std::string mask_name = name() + " mask";
  FieldIdentifier mask_fid (mask_name,diag_layout, nondim, gname);
  Field diag_mask(mask_fid);
  diag_mask.allocate_view();
  m_diagnostic_output.get_header().set_extra_data("mask_data",diag_mask);
  
}
// =========================================================================================
void FieldAtPressureLevel::compute_diagnostic_impl()
{
  using namespace scream::vinterp;

  //This is 2D source pressure
  const Field& pressure_f = get_field_in(m_pres_name);
  const auto pressure = pressure_f.get_view<const mPack**>();
  auto pres = view_Nd<mPack,2>("",pressure.extent_int(0),pressure.extent_int(1));
  Kokkos::deep_copy(pres,pressure);

  //input field
  const Field& f = get_field_in(m_field_name);

  // The setup for interpolation varies depending on the rank of the input field:
  const int rank = f.rank();

  if (rank==2) {
    const auto f_data_src = f.get_view<const mPack**>();
    auto fdat = view_Nd<mPack,2>("",f_data_src.extent_int(0),f_data_src.extent_int(1));
    Kokkos::deep_copy(fdat,f_data_src);
    //output field on new grid
    auto d_data_tgt = m_diagnostic_output.get_view<Real*>();
    view_Nd<mPack,2> data_tgt_tmp(reinterpret_cast<mPack*>(d_data_tgt.data()),d_data_tgt.extent_int(0),1);  // Note, vertical interp wants a 2D view, so we create a temporary one
    perform_vertical_interpolation<Real,1,2>(pres,m_p_tgt,fdat,data_tgt_tmp,m_num_levs,1,m_mask_val);

    // Track mask
    Kokkos::deep_copy(fdat,1.0); // reuse fdat
    auto extra_data = m_diagnostic_output.get_header().get_extra_data().at("mask_data");
    auto d_mask     = ekat::any_cast<Field>(extra_data);
    auto d_mask_tgt = d_mask.get_view<Real*>();
    view_Nd<mPack,2> mask_tgt_tmp(reinterpret_cast<mPack*>(d_mask_tgt.data()),d_mask_tgt.extent_int(0),1);  // Note, vertical interp wants a 2D view, so we create a temporary one
    perform_vertical_interpolation<Real,1,2>(pres,m_p_tgt,fdat,mask_tgt_tmp,m_num_levs,1,0);
  } else if (rank==3) {
    const auto f_data_src = f.get_view<const mPack***>();
    auto fdat = view_Nd<mPack,3>("",f_data_src.extent_int(0),f_data_src.extent_int(1),f_data_src.extent_int(2));
    Kokkos::deep_copy(fdat,f_data_src);
    //output field on new grid
    auto d_data_tgt = m_diagnostic_output.get_view<Real**>();
    view_Nd<mPack,3> data_tgt_tmp(reinterpret_cast<mPack*>(d_data_tgt.data()),d_data_tgt.extent_int(0),d_data_tgt.extent_int(1),1);  // Note, vertical interp wants a 2D view, so we create a temporary one

    perform_vertical_interpolation<Real,1,3>(pres,m_p_tgt,fdat,data_tgt_tmp,m_num_levs,1,m_mask_val);

    // Track mask
    Kokkos::deep_copy(fdat,1.0); // reuse fdat
    auto extra_data = m_diagnostic_output.get_header().get_extra_data().at("mask_data");
    auto d_mask     = ekat::any_cast<Field>(extra_data);
    auto d_mask_tgt = d_mask.get_view<Real*>();
    view_Nd<mPack,3> mask_tgt_tmp(reinterpret_cast<mPack*>(d_mask_tgt.data()),d_mask_tgt.extent_int(0),1);  // Note, vertical interp wants a 2D view, so we create a temporary one
    perform_vertical_interpolation<Real,1,3>(pres,m_p_tgt,fdat,mask_tgt_tmp,m_num_levs,1,0);
  } else {
    EKAT_ERROR_MSG("Error! field at pressure level only supports fields ranks 2 and 3 \n");
  }

}

} //namespace scream
