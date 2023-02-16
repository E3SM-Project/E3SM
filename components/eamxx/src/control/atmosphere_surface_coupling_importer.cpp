#include "atmosphere_surface_coupling_importer.hpp"

#include "share/property_checks/field_within_interval_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
// =========================================================================================
SurfaceCouplingImporter::SurfaceCouplingImporter (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{

}
// =========================================================================================
void SurfaceCouplingImporter::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();

  m_num_cols = m_grid->get_num_local_dofs();      // Number of columns on this rank
 
  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Qunit = kg/kg;
  Qunit.set_string("kg/kg");
  auto nondim = Units::nondimensional();
  auto Wm2 = W / m / m;
  Wm2.set_string("W/m2)");
  const auto m2 = m*m;

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  FieldLayout scalar2d_layout { {COL     }, {m_num_cols   } };
  FieldLayout vector2d_layout { {COL, CMP}, {m_num_cols, 2} };

  add_field<Computed>("sfc_alb_dir_vis",  scalar2d_layout, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dir_nir",  scalar2d_layout, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dif_vis",  scalar2d_layout, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dif_nir",  scalar2d_layout, nondim,  grid_name);
  add_field<Computed>("surf_lw_flux_up",  scalar2d_layout, W/m2,    grid_name);
  add_field<Computed>("surf_sens_flux",   scalar2d_layout, W/m2,    grid_name);
  add_field<Computed>("surf_evap",        scalar2d_layout, kg/m2/s, grid_name);
  add_field<Computed>("surf_mom_flux",    vector2d_layout, N/m2,    grid_name);
  add_field<Computed>("surf_radiative_T", scalar2d_layout, K,       grid_name);
  add_field<Computed>("T_2m",             scalar2d_layout, K,       grid_name);
  add_field<Computed>("qv_2m",            scalar2d_layout, Qunit,   grid_name);
  add_field<Computed>("wind_speed_10m",   scalar2d_layout, m/s,     grid_name);
  add_field<Computed>("snow_depth_land",  scalar2d_layout, m,       grid_name);
}
// =========================================================================================
  void SurfaceCouplingImporter::setup_surface_coupling_data(const SCDataManager &sc_data_manager)
{
  m_num_cpl_imports    = sc_data_manager.get_num_cpl_fields();
  m_num_scream_imports = sc_data_manager.get_num_scream_fields();

  EKAT_ASSERT_MSG(m_num_scream_imports <= m_num_cpl_imports,
                  "Error! More SCREAM imports than actual cpl imports.\n");
  EKAT_ASSERT_MSG(m_num_cols == sc_data_manager.get_field_size(),
                  "Error! Surface Coupling imports need to have size ncols.\n");

  // The import data is of size ncols,num_cpl_imports. All other data is of size num_scream_imports
  m_cpl_imports_view_h = decltype(m_cpl_imports_view_h) (sc_data_manager.get_field_data_ptr(),
                                                         m_num_cols, m_num_cpl_imports);
  m_cpl_imports_view_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(),
                                                             m_cpl_imports_view_h);
  m_import_field_names = new name_t[m_num_scream_imports];
  std::memcpy(m_import_field_names, sc_data_manager.get_field_name_ptr(), m_num_scream_imports*32*sizeof(char));

  m_cpl_indices_view =
    decltype(m_cpl_indices_view)          (sc_data_manager.get_field_cpl_indices_ptr(),
					 m_num_scream_imports);
  m_vector_components_view =
      decltype(m_vector_components_view)  (sc_data_manager.get_field_vector_components_ptr(),
                                           m_num_scream_imports);
  m_constant_multiple_view =
      decltype(m_constant_multiple_view)  (sc_data_manager.get_field_constant_multiple_ptr(),
                                           m_num_scream_imports);
  m_do_import_during_init_view =
    decltype(m_do_import_during_init_view)(sc_data_manager.get_field_transfer_during_init_ptr(),
					   m_num_scream_imports);

  m_column_info_d = decltype(m_column_info_d) ("m_info", m_num_scream_imports);
  m_column_info_h = Kokkos::create_mirror_view(m_column_info_d);
}
// =========================================================================================
void SurfaceCouplingImporter::initialize_impl (const RunType /* run_type */)
{
  bool any_initial_imports = false;

  for (int i=0; i<m_num_scream_imports; ++i) {

    std::string fname = m_import_field_names[i];
    Field field = get_field_out(fname);
    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");

    // Set view data ptr
    m_column_info_h(i).data = field.get_internal_view_data<Real>();

    // Get column info from field utility function
    get_col_info_for_surface_values(field.get_header_ptr(),
                                    m_vector_components_view(i),
                                    m_column_info_h(i).col_offset, m_column_info_h(i).col_stride);

    // Set constant multiple
    m_column_info_h(i).constant_multiple = m_constant_multiple_view(i);

    // Set whether or not this field should be imported during initialization
    m_column_info_h(i).transfer_during_initialization = m_do_import_during_init_view(i);
    if (m_do_import_during_init_view(i)) any_initial_imports = true;

    // Set index for referencing cpl data.
    m_column_info_h(i).cpl_indx = m_cpl_indices_view(i);
  }

  // Copy data to device for use in do_import()
  Kokkos::deep_copy(m_column_info_d, m_column_info_h);

  // Set property checks for fields in this proces
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dir_vis"),m_grid,0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dir_nir"),m_grid,0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dif_vis"),m_grid,0.0,1.0,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dif_nir"),m_grid,0.0,1.0,true);

  // Perform initial import (if any are marked for import during initialization) 
  if (any_initial_imports) do_import(true);
}
// =========================================================================================
void SurfaceCouplingImporter::run_impl (const double /* dt */)
{
  do_import();
}
// =========================================================================================
void SurfaceCouplingImporter::do_import(const bool called_during_initialization)
{
  using policy_type = KokkosTypes<DefaultDevice>::RangePolicy;

  // Local copies, to deal with CUDA's handling of *this
  const auto col_info           = m_column_info_d;
  const auto cpl_imports_view_d = m_cpl_imports_view_d;
  const int  num_cols           = m_num_cols;
  const int  num_imports        = m_num_scream_imports;

  // Deep copy cpl host array to devic
  Kokkos::deep_copy(m_cpl_imports_view_d,m_cpl_imports_view_h);

  // Unpack the fields
  auto unpack_policy = policy_type(0,num_imports*num_cols);
  Kokkos::parallel_for(unpack_policy, KOKKOS_LAMBDA(const int& i) {
    const int ifield = i / num_cols;
    const int icol   = i % num_cols;

    const auto& info = col_info(ifield);

    auto offset = icol*info.col_stride + info.col_offset;

    // if this is during initialization, check whether or not the field should be imported
    bool do_import = (not called_during_initialization || info.transfer_during_initialization);
    if (do_import) {
      info.data[offset] = cpl_imports_view_d(icol,info.cpl_indx)*info.constant_multiple;
    }
  });
}
// =========================================================================================
void SurfaceCouplingImporter::finalize_impl()
{

}
// =========================================================================================
} // namespace scream
