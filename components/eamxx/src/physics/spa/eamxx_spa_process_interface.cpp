#include "eamxx_spa_process_interface.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

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
  EKAT_REQUIRE_MSG(m_params.isParameter("spa_data_file"),
      "ERROR: spa_data_file is missing from SPA parameter list.");
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto nondim = Units::nondimensional();

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar1d_mid = m_grid->get_vertical_layout(true);
  // Use VAR field tag for gases for now; consider adding a tag?
  FieldLayout scalar3d_swband = m_grid->get_3d_vector_layout(true,m_nswbands,"swband");
  FieldLayout scalar3d_lwband = m_grid->get_3d_vector_layout(true,m_nlwbands,"lwband");

  // Set of fields used strictly as input
  constexpr int ps = Spack::n;
  add_field<Required>("p_mid"      , scalar3d_mid, Pa,     grid_name, ps);

  // Set of fields used strictly as output
  add_field<Computed>("nccn",        scalar3d_mid,    1/kg,   grid_name, ps);
  add_field<Computed>("aero_g_sw",   scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_ssa_sw", scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_tau_sw", scalar3d_swband, nondim, grid_name, ps);
  add_field<Computed>("aero_tau_lw", scalar3d_lwband, nondim, grid_name, ps);

  // We can already create some of the spa structures

  // 1. Create SPAHorizInterp remapper
  auto spa_data_file = m_params.get<std::string>("spa_data_file");
  auto spa_map_file  = m_params.get<std::string>("spa_remap_file","");

  // IOP cases cannot have a remap file. IOP file reader is itself a remaper,
  // where a single column of data corresponding to the closest lat/lon pair to
  // the IOP lat/lon parameters is read from file, and that column data is mapped
  // to all columns of the IdentityRemapper source fields.
  EKAT_REQUIRE_MSG(spa_map_file == "" or spa_map_file == "None" or not m_iop,
    "Error! Cannot define spa_remap_file for cases with an Intensive Observation Period defined. "
    "The IOP class defines it's own remap from file data -> model data.\n");

  SPAHorizInterp = SPAFunc::create_horiz_remapper (m_grid,spa_data_file,spa_map_file, m_iop!=nullptr);

  // Grab a sw and lw field from the horiz interp, and check sw/lw dim against what we hardcoded in this class
  auto nswbands_data = SPAHorizInterp->get_src_field(4).get_header().get_identifier().get_layout().dim("swband");
  auto nlwbands_data = SPAHorizInterp->get_src_field(5).get_header().get_identifier().get_layout().dim("lwband");
  EKAT_REQUIRE_MSG (nswbands_data==m_nswbands,
      "Error! Spa data file has a different number of sw bands than the model.\n"
      " - spa data swbands: " + std::to_string(nswbands_data) + "\n"
      " - model swbands   : " + std::to_string(m_nswbands) + "\n");
  EKAT_REQUIRE_MSG (nlwbands_data==m_nlwbands,
      "Error! Spa data file has a different number of lw bands than the model.\n"
      " - spa data lwbands: " + std::to_string(nlwbands_data) + "\n"
      " - model lwbands   : " + std::to_string(m_nlwbands) + "\n");

  const auto io_grid = SPAHorizInterp->get_src_grid();

  // 2. Initialize the size of the SPAData structures.
  // Note: add 2 to number of levels to allow extrapolation if model pressure is outside the data range
  m_num_src_levs = io_grid->get_num_vertical_levels();
  SPAData_start = SPAFunc::SPAInput(m_num_cols, m_num_src_levs+2, m_nswbands, m_nlwbands);
  SPAData_end   = SPAFunc::SPAInput(m_num_cols, m_num_src_levs+2, m_nswbands, m_nlwbands);
  SPAData_out.init(m_num_cols,m_num_levs,m_nswbands,m_nlwbands,false);

  // 3 Read in hyam/hybm in start/end data, and pad them
  Field hyam(FieldIdentifier("hyam",io_grid->get_vertical_layout(true),nondim,io_grid->name()));
  Field hybm(FieldIdentifier("hybm",io_grid->get_vertical_layout(true),nondim,io_grid->name()));
  hyam.allocate_view();
  hybm.allocate_view();

  AtmosphereInput hvcoord_reader(spa_data_file,io_grid,{hyam,hybm},true);
  hvcoord_reader.read_variables();
  hvcoord_reader.finalize();

  // Do the copy on host, cause it's just easier,
  // and this is just a cheap loop during init
  auto hyam_h = hyam.get_view<const Real*,Host>();
  auto hybm_h = hybm.get_view<const Real*,Host>();
  auto nlevs = io_grid->get_num_vertical_levels();
  for (auto data : {SPAData_start, SPAData_end} ) {
    auto spa_hyam = ekat::scalarize(data.hyam);
    auto spa_hybm = ekat::scalarize(data.hybm);
    auto spa_hyam_h = Kokkos::create_mirror_view(spa_hyam);
    auto spa_hybm_h = Kokkos::create_mirror_view(spa_hybm);
    for (int i=0; i<nlevs; ++i) {
      spa_hyam_h(i+1) = hyam_h(i);
      spa_hybm_h(i+1) = hybm_h(i);
    }
    spa_hyam_h(0) = 0;
    spa_hyam_h(nlevs+1) = 1e5;
    spa_hybm_h(0) = 0;
    spa_hybm_h(nlevs+1) = 0;

    Kokkos::deep_copy(spa_hyam,spa_hyam_h);
    Kokkos::deep_copy(spa_hybm,spa_hybm_h);
  }

  // 4. Create reader for spa data. The reader is either an
  //    AtmosphereInput object (for reading into standard
  //    grids) or a SpaFunctions::IOPReader (for reading into
  //    an IOP grid).
  if (m_iop) {
    SPAIOPDataReader = SPAFunc::create_spa_data_reader(m_iop,SPAHorizInterp,spa_data_file);
  } else {
    SPADataReader = SPAFunc::create_spa_data_reader(SPAHorizInterp,spa_data_file);
  }
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
  // Initialize SPAData_out with the views from the out fields
  SPAData_out.CCN3       = get_field_out("nccn").get_view<Spack**>();
  SPAData_out.AER_G_SW   = get_field_out("aero_g_sw").get_view<Spack***>();
  SPAData_out.AER_SSA_SW = get_field_out("aero_ssa_sw").get_view<Spack***>();
  SPAData_out.AER_TAU_SW = get_field_out("aero_tau_sw").get_view<Spack***>();
  SPAData_out.AER_TAU_LW = get_field_out("aero_tau_lw").get_view<Spack***>();

  // Load the first month into spa_end.
  // Note: At the first time step, the data will be moved into spa_beg,
  //       and spa_end will be reloaded from file with the new month.
  const int curr_month = timestamp().get_month()-1; // 0-based
  SPAFunc::update_spa_data_from_file(SPADataReader,SPAIOPDataReader,timestamp(),curr_month,*SPAHorizInterp,SPAData_end);

  // 6. Set property checks for fields in this process
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
    SPAFunc::update_spa_timestate(SPADataReader,SPAIOPDataReader,ts,*SPAHorizInterp,SPATimeState,SPAData_start,SPAData_end);

  // Call the main SPA routine to get interpolated aerosol forcings.
  const auto& pmid_tgt = get_field_in("p_mid").get_view<const Spack**>();
  SPAFunc::spa_main(SPATimeState, pmid_tgt, m_buffer.p_mid_src,
                    SPAData_start,SPAData_end,m_buffer.spa_temp,SPAData_out);
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

} // namespace scream
