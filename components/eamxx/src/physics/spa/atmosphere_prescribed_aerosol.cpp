#include "atmosphere_prescribed_aerosol.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"

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

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column
  m_dofs_gids = m_grid->get_dofs_gids().get_view<const gid_type*>();
  m_min_global_dof    = m_grid->get_global_min_dof_gid();

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

  // Set of fields used strictly as output
  add_field<Computed>("nccn",   scalar3d_layout_mid,    1/kg,   grid_name,ps);
  add_field<Computed>("aero_g_sw",      scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_ssa_sw",    scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_tau_sw",    scalar3d_swband_layout, nondim, grid_name,ps);
  add_field<Computed>("aero_tau_lw",    scalar3d_lwband_layout, nondim, grid_name,ps);

  // Init output data structure
  SPAData_out.init(m_num_cols,m_num_levs,m_nswbands,m_nlwbands,false);

  // Note: only the number of levels associated with this data haven't been set.  We can
  //       take this information directly from the spa data file.
  m_spa_data_file = m_params.get<std::string>("spa_data_file");
  scorpio::register_file(m_spa_data_file,scorpio::Read);
  m_num_src_levs = scorpio::get_dimlen_c2f(m_spa_data_file.c_str(),"lev");
  scorpio::eam_pio_closefile(m_spa_data_file);
  SPAHorizInterp.m_comm = m_comm;

}
// =========================================================================================
size_t SPA::requested_buffer_size_in_bytes() const
{
  using PackInfo = ekat::PackInfo<Spack::n>;

  // Recall: the quantities in spa_temp defined over vlevs have 1 Real of
  //         padding in each column (at beginning and end).
  //         That's why we have m_num_levs+2
  const int nlevs = m_num_src_levs+2;
  const int num_mid_packs = PackInfo::num_packs(nlevs);
  const int nlevs_alloc = num_mid_packs*Spack::n;

  // We have
  //  - one (ncols) view (spa_temp's ps)
  //  - two (ncols,nlevs) mid view (p_mid_src, spa_temp's ccn)
  //  - three (ncols,nswbands,nlevs) views (spa_temp's aer_g_sw, aer_ssa_sw, aer_tau_sw)
  //  - one (ncols,nlwbands,nlevs) view (aer_tau_lw)
  const int num_reals = m_num_cols*(1+nlevs_alloc*(2 + 3*m_nswbands + m_nlwbands));

  return num_reals*sizeof(Real);
}

// =========================================================================================
void SPA::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Buffers size not sufficient.\n");

  using PackInfo = ekat::PackInfo<Spack::n>;

  // Short names make following rows fit on text editor screen
  // Recall: the quantities in spa_temp defined over vlevs have 1 Real of
  //         padding in each column (at beginning and end).
  //         That's why we have m_num_levs+2
  const int nlevs  = m_num_src_levs+2;
  const int npacks = PackInfo::num_packs(nlevs);
  const int ncols  = m_num_cols;
  const int nswb   = m_nswbands;
  const int nlwb   = m_nlwbands;

  Spack* mem = reinterpret_cast<Spack*>(buffer_manager.get_memory());

  // Source pressure levels
  m_buffer.p_mid_src = decltype(m_buffer.p_mid_src)(mem, ncols, npacks);
  mem += m_buffer.p_mid_src.size();

  // SPA temporaries
  auto& spa_data = m_buffer.spa_temp.data;
  spa_data.init(ncols,nlevs,nswb,nlwb,false);

  spa_data.CCN3 = decltype(spa_data.CCN3)(mem, ncols, npacks);
  mem += spa_data.CCN3.size();

  spa_data.AER_G_SW = decltype(spa_data.AER_G_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_G_SW.size();
  spa_data.AER_SSA_SW = decltype(spa_data.AER_SSA_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_SSA_SW.size();
  spa_data.AER_TAU_SW = decltype(spa_data.AER_TAU_SW)(mem, ncols, nswb, npacks);
  mem += spa_data.AER_TAU_SW.size();

  spa_data.AER_TAU_LW = decltype(spa_data.AER_TAU_LW)(mem, ncols, nlwb, npacks);
  mem += spa_data.AER_TAU_LW.size();

  Real* r_mem = reinterpret_cast<Real*>(mem);
  m_buffer.spa_temp.PS = decltype(m_buffer.spa_temp.PS)(r_mem,ncols);
  r_mem += m_buffer.spa_temp.PS.size();

  size_t used_mem = (r_mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for SPA.\n"
      "   - used mem     : " + std::to_string(used_mem) + "\n"
      "   - requested mem: " + std::to_string(requested_buffer_size_in_bytes()) + "\n");
}

// =========================================================================================
void SPA::initialize_impl (const RunType /* run_type */)
{
  // Initialize SPA pressure state stucture and set pointers for the SPA output data to
  // field managed variables.
  SPAData_out.CCN3               = get_field_out("nccn").get_view<Pack**>();
  SPAData_out.AER_G_SW           = get_field_out("aero_g_sw").get_view<Pack***>();
  SPAData_out.AER_SSA_SW         = get_field_out("aero_ssa_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_SW         = get_field_out("aero_tau_sw").get_view<Pack***>();
  SPAData_out.AER_TAU_LW         = get_field_out("aero_tau_lw").get_view<Pack***>();

  // Retrieve the remap and data file locations from the parameter list:
  EKAT_REQUIRE_MSG(m_params.isParameter("spa_remap_file"),"ERROR: spa_remap_file is missing from SPA parameter list.");
  EKAT_REQUIRE_MSG(m_params.isParameter("spa_data_file"),"ERROR: spa_data_file is missing from SPA parameter list.");
  m_spa_remap_file = m_params.get<std::string>("spa_remap_file");

  // Set the SPA remap weights.  
  // TODO: We may want to provide an option to calculate weights on-the-fly. 
  //       If so, then the EKAT_REQUIRE_MSG above will need to be removed and 
  //       we can have a default m_spa_data_file option that is online calculation.
  using ci_string = ekat::CaseInsensitiveString;
  ci_string no_filename = "none";
  if (m_spa_remap_file == no_filename) {
    if (m_comm.am_i_root()) {
      printf("WARNING: spa_remap_file has been set to 'NONE', assuming that SPA data and simulation are on the same grid - skipping horizontal interpolation\n");
    }
    SPAFunc::set_remap_weights_one_to_one(m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  } else {
    SPAFunc::get_remap_weights_from_file(m_spa_remap_file,m_min_global_dof,m_dofs_gids,SPAHorizInterp);
  }

  // Initialize the size of the SPAData structures:  add 2 to number of levels for padding
  SPAData_start = SPAFunc::SPAInput(m_dofs_gids.size(), m_num_src_levs+2, m_nswbands, m_nlwbands);
  SPAData_end   = SPAFunc::SPAInput(m_dofs_gids.size(), m_num_src_levs+2, m_nswbands, m_nlwbands);

  // Update the local time state information and load the first set of SPA data for interpolation:
  auto ts = timestamp();
  SPATimeState.inited = false;
  SPATimeState.current_month = ts.get_month();
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  const auto eps = std::numeric_limits<double>::epsilon();

  add_postcondition_check<Interval>(get_field_out("nccn"),m_grid,0,1e11,true,-1e11*eps);
  // TODO: add an epslon to max possible upper bound of aero_ssa_sw?
  add_postcondition_check<Interval>(get_field_out("aero_g_sw"),m_grid,0.0,1.0,true);
  add_postcondition_check<Interval>(get_field_out("aero_ssa_sw"),m_grid,0.0,1.0,true);
  add_postcondition_check<Interval>(get_field_out("aero_tau_sw"),m_grid,0.0,1.0,true);
  add_postcondition_check<Interval>(get_field_out("aero_tau_lw"),m_grid,0.0,1.0,true);
}

// =========================================================================================
void SPA::run_impl (const double dt)
{
  /* Gather time and state information for interpolation */
  auto ts = timestamp()+dt;
  /* Update the SPATimeState to reflect the current time, note the addition of dt */
  SPATimeState.t_now = ts.frac_of_year_in_days();
  /* Update time state and if the month has changed, update the data.*/
  SPAFunc::update_spa_timestate(m_spa_data_file,m_nswbands,m_nlwbands,ts,SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // Call the main SPA routine to get interpolated aerosol forcings.
  const auto& pmid_tgt = get_field_in("p_mid").get_view<const Pack**>();
  SPAFunc::spa_main(SPATimeState, pmid_tgt, m_buffer.p_mid_src,
                    SPAData_start,SPAData_end,m_buffer.spa_temp,SPAData_out);
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

} // namespace scream
