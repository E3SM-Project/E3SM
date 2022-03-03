#include "atmosphere_prescribed_aerosol.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace spa;
// =========================================================================================
SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto nondim = Units::nondimensional();

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column
  m_dofs_gids = grid->get_dofs_gids();
  m_total_global_dofs = grid->get_num_global_dofs();
  m_min_global_dof    = grid->get_global_min_dof_gid();

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  FieldLayout scalar2d_layout     { {COL},     {m_num_cols} };
  FieldLayout scalar1d_layout_mid { {LEV},     {m_num_levs} };
  // Use VAR field tag for gases for now; consider adding a tag?
  FieldLayout scalar3d_swband_layout { {COL,SWBND, LEV}, {m_num_cols, m_nswbands, m_num_levs} }; 
  FieldLayout scalar3d_lwband_layout { {COL,LWBND, LEV}, {m_num_cols, m_nlwbands, m_num_levs} }; 

  // Set of fields used strictly as input
  constexpr int ps = Pack::n;
  add_field<Required>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Required>("hyam"       , scalar1d_layout_mid, nondim, grid_name, ps); // TODO: These fields should  be loaded from file and not registered with the field manager.
  add_field<Required>("hybm"       , scalar1d_layout_mid, nondim, grid_name, ps); // TODO: These fields should  be loaded from file and not registered with the field manager.

  // Set of fields used strictly as output
  add_field<Computed>("nc_activated",   scalar3d_layout_mid,    1/kg,   grid_name,ps);
  add_field<Computed>("aero_g_sw",      scalar3d_swband_layout, 1/kg,   grid_name,ps);
  add_field<Computed>("aero_ssa_sw",    scalar3d_swband_layout, 1/kg,   grid_name,ps);
  add_field<Computed>("aero_tau_sw",    scalar3d_swband_layout, 1/kg,   grid_name,ps);
  add_field<Computed>("aero_tau_lw",    scalar3d_lwband_layout, 1/kg,   grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}
// =========================================================================================
int SPA::requested_buffer_size_in_bytes() const
{
  const Int num_mid_packs    = ekat::npack<Spack>(m_num_levs);
  const Int num_int_packs = ekat::npack<Spack>(m_num_levs+1);

  // Number of Reals needed by local views in the interface
  const int interface_request =
      // 1d view scalar, size (ncol)
      Buffer::num_1d_scalar*m_num_cols*sizeof(Real) +
      // 2d view packed, size (ncol, nlev_packs)
      Buffer::num_2d_vector*m_num_cols*num_mid_packs*sizeof(Spack) +
      Buffer::num_2dp1_vector*m_num_cols*num_int_packs*sizeof(Spack);

  // Number of Reals needed by the WorkspaceManager
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, num_mid_packs);
  const int wsm_request   = WSM::get_total_bytes_needed(num_mid_packs, 3, policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void SPA::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  m_buffer.ps_src = decltype(m_buffer.ps_src)(mem, m_num_cols);
  mem += m_buffer.ps_src.size();

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2d packed views
  const Int num_mid_packs    = ekat::npack<Spack>(m_num_levs);

  m_buffer.p_mid_src = decltype(m_buffer.p_mid_src)(s_mem, m_num_cols, num_mid_packs);
  s_mem += m_buffer.p_mid_src.size();
  m_buffer.ccn3_src = decltype(m_buffer.ccn3_src)(s_mem, m_num_cols, num_mid_packs);
  s_mem += m_buffer.ccn3_src.size();

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy  = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, num_mid_packs);
  const int wsm_size = WSM::get_total_bytes_needed(num_mid_packs, 3, policy)/sizeof(Spack);
  s_mem += wsm_size;

  int used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SPA.");
}

// =========================================================================================
void SPA::initialize_impl (const RunType /* run_type */)
{
  // Initialize SPA pressure state stucture and set pointers for the SPA output data to
  // field managed variables.
  SPAPressureState.ncols         = m_num_cols;
  SPAPressureState.nlevs         = m_num_levs;
  SPAPressureState.hyam          = get_field_in("hyam").get_view<const Pack*>();
  SPAPressureState.hybm          = get_field_in("hybm").get_view<const Pack*>();
  SPAPressureState.pmid          = get_field_in("p_mid").get_view<const Pack**>();
  SPAData_out.CCN3               = get_field_out("nc_activated").get_view<Pack**>();
  SPAData_out.AER_G_SW           = get_field_out("aero_g_sw").get_view<Pack***>();
  SPAData_out.AER_SSA_SW         = get_field_out("aero_ssa_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_SW         = get_field_out("aero_tau_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_LW         = get_field_out("aero_tau_lw").get_view<Pack***>();

  // Retrieve the remap and data file locations from the parameter list:
  EKAT_REQUIRE_MSG(m_params.isParameter("SPA Remap File"),"ERROR: SPA Remap File is missing from SPA parameter list.");
  EKAT_REQUIRE_MSG(m_params.isParameter("SPA Data File"),"ERROR: SPA Data File is missing from SPA parameter list.");
  m_spa_remap_file = m_params.get<std::string>("SPA Remap File");
  m_spa_data_file = m_params.get<std::string>("SPA Data File");

  // Set the SPA remap weights.  
  // TODO: We may want to provide an option to calculate weights on-the-fly. 
  //       If so, then the EKAT_REQUIRE_MSG above will need to be removed and 
  //       we can have a default m_spa_data_file option that is online calculation.
  using ci_string = ekat::CaseInsensitiveString;
  ci_string no_filename = "none";
  if (m_spa_remap_file == no_filename) {
    printf("WARNING: SPA Remap File has been set to 'NONE', assuming that SPA data and simulation are on the same grid - skipping horizontal interpolation");
    SPAFunc::set_remap_weights_one_to_one(m_total_global_dofs,m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  } else {
    SPAFunc::get_remap_weights_from_file(m_spa_remap_file,m_total_global_dofs,m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  }
  // Note: only the number of levels associated with this data haven't been set.  We can
  //       take this information directly from the spa data file.
  scorpio::register_file(m_spa_data_file,scorpio::Read);
  const int source_data_nlevs = scorpio::get_dimlen_c2f(m_spa_data_file.c_str(),"lev")+2; // Add 2 for padding
  SPAHorizInterp.m_comm = m_comm;

  // Initialize the size of the SPAData structures:
  SPAData_start = SPAFunc::SPAData(m_dofs_gids.size(), source_data_nlevs, m_nswbands, m_nlwbands);
  SPAData_end   = SPAFunc::SPAData(m_dofs_gids.size(), source_data_nlevs, m_nswbands, m_nlwbands);

  // Update the local time state information and load the first set of SPA data for interpolation:
  auto ts = timestamp();
  SPATimeState.inited = false;
  SPATimeState.current_month = ts.get_month();
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);
}

// =========================================================================================
void SPA::run_impl (const int /* dt */)
{
  /* Gather time and state information for interpolation */
  auto ts = timestamp();
  /* Update time state and if the month has changed, update the data.*/
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // Call the main SPA routine to get interpolated aerosol forcings.
  SPAFunc::spa_main(SPATimeState, SPAPressureState,SPAData_start,SPAData_end,SPAData_out,m_num_cols,m_num_levs,m_nswbands,m_nlwbands);
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

} // namespace scream
