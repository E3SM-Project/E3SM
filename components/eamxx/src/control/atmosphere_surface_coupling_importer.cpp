#include "atmosphere_surface_coupling_importer.hpp"

#include "share/property_checks/field_within_interval_check.hpp"
#include "physics/share/physics_constants.hpp"

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
  constexpr auto nondim = Units::nondimensional();
  constexpr auto m2 = pow(m, 2);

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  const FieldLayout scalar2d = m_grid->get_2d_scalar_layout();
  const FieldLayout vector2d = m_grid->get_2d_vector_layout(2);
  const FieldLayout vector4d = m_grid->get_2d_vector_layout(4);

  add_field<Computed>("sfc_alb_dir_vis",  scalar2d, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dir_nir",  scalar2d, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dif_vis",  scalar2d, nondim,  grid_name);
  add_field<Computed>("sfc_alb_dif_nir",  scalar2d, nondim,  grid_name);
  add_field<Computed>("surf_lw_flux_up",  scalar2d, W/m2,    grid_name);
  add_field<Computed>("surf_sens_flux",   scalar2d, W/m2,    grid_name);
  add_field<Computed>("surf_evap",        scalar2d, kg/m2/s, grid_name);
  add_field<Computed>("surf_mom_flux",    vector2d, N/m2,    grid_name);
  add_field<Computed>("surf_radiative_T", scalar2d, K,       grid_name);
  add_field<Computed>("T_2m",             scalar2d, K,       grid_name);
  add_field<Computed>("qv_2m",            scalar2d, kg/kg,   grid_name);
  add_field<Computed>("wind_speed_10m",   scalar2d, m/s,     grid_name);
  add_field<Computed>("snow_depth_land",  scalar2d, m,       grid_name);
  add_field<Computed>("ocnfrac",          scalar2d, nondim,  grid_name);
  add_field<Computed>("landfrac",         scalar2d, nondim,  grid_name);
  add_field<Computed>("icefrac",          scalar2d, nondim,  grid_name);
  // Friction velocity [m/s]
  add_field<Computed>("fv",               scalar2d, m/s,     grid_name);
  // Aerodynamical resistance
  add_field<Computed>("ram1",             scalar2d, s/m,     grid_name);
  // Sea surface temperature [K]
  add_field<Computed>("sst",              scalar2d, K,       grid_name);
  //dust fluxes [kg/m^2/s]: Four flux values for eacch column
  add_field<Computed>("dstflx",           vector4d, kg/m2/s, grid_name);

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
#ifdef HAVE_MOAB
  // The import data is of size num_cpl_imports, ncol. All other data is of size num_scream_imports
  m_moab_cpl_imports_view_h = decltype(m_moab_cpl_imports_view_h) (sc_data_manager.get_field_data_moab_ptr(),
                                                         m_num_cpl_imports, m_num_cols);
  m_moab_cpl_imports_view_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(),
                                                             m_moab_cpl_imports_view_h);
#endif
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
#ifdef HAVE_MOAB
  // Deep copy cpl host array to device
  const auto moab_cpl_imports_view_d = m_moab_cpl_imports_view_d;
  Kokkos::deep_copy(m_moab_cpl_imports_view_d,m_moab_cpl_imports_view_h);
#endif


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

#ifdef HAVE_MOAB
  Kokkos::parallel_for(unpack_policy, KOKKOS_LAMBDA(const int& i) {

    const int icol   = i / num_imports;
    const int ifield = i % num_imports;

    const auto& info = col_info(ifield);

    auto offset = icol*info.col_stride + info.col_offset;

    // if this is during initialization, check whether or not the field should be imported
    bool do_import = (not called_during_initialization || info.transfer_during_initialization);
    if (do_import) {
      info.data[offset] = moab_cpl_imports_view_d(info.cpl_indx, icol)*info.constant_multiple;
    }
  });
#endif

  if (m_iop_data_manager) {
    if (m_iop_data_manager->get_params().get<bool>("iop_srf_prop")) {
      // Overwrite imports with data from IOP file
      overwrite_iop_imports(called_during_initialization);
    }
  }
}
// =========================================================================================
void SurfaceCouplingImporter::overwrite_iop_imports (const bool called_during_initialization)
{
  using policy_type = KokkosTypes<DefaultDevice>::RangePolicy;
  using C = physics::Constants<Real>;

  const auto has_lhflx = m_iop_data_manager->has_iop_field("lhflx");
  const auto has_shflx = m_iop_data_manager->has_iop_field("shflx");
  const auto has_Tg    = m_iop_data_manager->has_iop_field("Tg");

  // Read IOP file for current time step, if necessary
  // TODO: this is using the TS from the beg of the step. Should it use end_of_step_ts() instead?
  m_iop_data_manager->read_iop_file_data(start_of_step_ts());

  static constexpr Real latvap = C::LatVap;
  static constexpr Real stebol = C::stebol;

  const auto& col_info_h = m_column_info_h;
  const auto& col_info_d = m_column_info_d;

  for (int ifield=0; ifield<m_num_scream_imports; ++ifield) {
    const std::string fname = m_import_field_names[ifield];
    const auto& info_h = col_info_h(ifield);

    // If we are in initialization and field should not be imported, skip
    if (called_during_initialization && not info_h.transfer_during_initialization) {
      continue;
    }

    // Store IOP surf data into col_val
    Real col_val(std::nan(""));
    if (fname == "surf_evap" && has_lhflx) {
      const auto f = m_iop_data_manager->get_iop_field("lhflx");
      f.sync_to_host();
      col_val = f.get_view<Real, Host>()()/latvap;
    } else if (fname == "surf_sens_flux" && has_shflx) {
      const auto f = m_iop_data_manager->get_iop_field("shflx");
      f.sync_to_host();
      col_val = f.get_view<Real, Host>()();
    } else if (fname == "surf_radiative_T" && has_Tg) {
      const auto f = m_iop_data_manager->get_iop_field("Tg");
      f.sync_to_host();
      col_val = f.get_view<Real, Host>()();
    } else if (fname == "surf_lw_flux_up" && has_Tg) {
      const auto f = m_iop_data_manager->get_iop_field("Tg");
      f.sync_to_host();
      col_val = stebol*std::pow(f.get_view<Real, Host>()(), 4);
    } else {
      // If import field doesn't satisify above, skip
      continue;
    }

    // Overwrite iop imports with col_val for each column
    auto policy = policy_type(0, m_num_cols);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int& icol) {
      const auto& info_d = col_info_d(ifield);
      const auto offset = icol*info_d.col_stride + info_d.col_offset;
      info_d.data[offset] = col_val;
#ifdef HAVE_MOAB
   //  TODO
#endif
    });
  }
}
// =========================================================================================
void SurfaceCouplingImporter::finalize_impl()
{

}
// =========================================================================================
} // namespace scream
