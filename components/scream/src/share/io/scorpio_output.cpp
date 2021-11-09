#include "share/io/scorpio_output.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "share/io/scorpio_input.hpp"

#include "ekat/util/ekat_string_utils.hpp"

#include <numeric>
#include <fstream>

namespace scream
{

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const fm_type>& field_mgr,
                  const std::shared_ptr<const gm_type>& grids_mgr)
 : m_comm      (comm)
{
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  // Set grid, and possibly create a remapper
  set_grid(field_mgr->get_grid());

  // Parse the input parameter list
  set_params(params,field_mgr->get_grid()->name());

  // Sets the intermal field mgr, and possibly sets up the remapper
  set_field_manager(field_mgr,grids_mgr);

  // Setup I/O structures
  init ();
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::restart (const std::string& filename)
{

  // Create an input stream on the fly, and init averaging data
  ekat::ParameterList res_params("Input Parameters");
  res_params.set<std::string>("Filename",filename);
  res_params.set("Fields",m_fields_names);

  AtmosphereInput hist_restart (m_comm,res_params,m_grid,m_host_views_1d,m_layouts);
  hist_restart.read_variables();
  hist_restart.finalize();
}

void AtmosphereOutput::init()
{
  for (const auto& var_name : m_fields_names) {
    register_dimensions(var_name);
  }

  // Now that the fields have been gathered register the local views which will be used to determine output data to be written.
  register_views();

} // init
/*-----*/
void AtmosphereOutput::run (const std::string& filename, const bool is_write_step, const int nsteps_since_last_output)
{
  using namespace scream::scorpio;

  // If needed, remap fields from their grid to the unique grid, for I/O
  if (m_remapper) {
    m_remapper->remap(true);
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields_names) {
    // Get all the info for this field.
    const auto  field = m_field_mgr->get_field(name);
    const auto& layout = m_layouts.at(name);
    const auto& dims = layout.dims();
    const auto  rank = layout.rank();

    // Safety check: make sure that the field was written at least once before using it.
    EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Output field '" + name + "' has not been initialized yet\n.");

    // Make sure host data is up to date
    field.sync_to_host();

    // Manually update the 'running-tally' views with data from the field,
    // by combining new data with current avg values.
    // NOTE: the running-tally is not a tally for Instant avg_type.
    auto data = m_host_views_1d.at(name).data();
    switch (rank) {
      case 1:
      {
        auto new_view_1d = field.get_view<const Real*,Host>();
        auto avg_view_1d = view_Nd_host<1>(data,dims[0]);

        for (int i=0; i<dims[0]; ++i) {
          combine(new_view_1d(i), avg_view_1d(i),nsteps_since_last_output);
        }
        break;
      }
      case 2:
      {
        auto new_view_2d = field.get_view<const Real**,Host>();
        auto avg_view_2d = view_Nd_host<2>(data,dims[0],dims[1]);
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            combine(new_view_2d(i,j), avg_view_2d(i,j),nsteps_since_last_output);
        }}
        break;
      }
      case 3:
      {
        auto new_view_3d = field.get_view<const Real***,Host>();
        auto avg_view_3d = view_Nd_host<3>(data,dims[0],dims[1],dims[2]);
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              combine(new_view_3d(i,j,k), avg_view_3d(i,j,k),nsteps_since_last_output);
        }}}
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
    }

    if (is_write_step) {
      grid_write_data_array(filename,name,data);
    }
  }
} // run

/* ---------------------------------------------------------- */

void AtmosphereOutput::
set_params (const ekat::ParameterList& params, const std::string& grid_name)
{
  // Parse the parameters that controls this output instance.
  auto avg_type = params.get<std::string>("Averaging Type");
  m_avg_type = ekat::upper_case(avg_type);
  auto valid_avg_types = {"INSTANT", "MAX", "MIN", "AVERAGE"};
  EKAT_REQUIRE_MSG (ekat::contains(valid_avg_types,m_avg_type),
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

  // For each output field, ensure its dimensions are registered in the pio file
  using vos_t = std::vector<std::string>;
  if (params.isParameter("Fields")) {
    m_fields_names = params.get<vos_t>("Fields");
  } else {
    if (params.isSublist("Fields") && params.sublist("Fields").isParameter(grid_name)) {
      m_fields_names = params.sublist("Fields").get<vos_t>(grid_name);
    }
  }
}

void AtmosphereOutput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_field_mgr, "Error! Field manager was already set.\n");
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");

  m_field_mgr = field_mgr;

  if (not m_field_mgr->get_grid()->is_unique()) {
    // The grid is not unique. Build a unique field manager
    // and build the remapper, and we call this method again.
    build_remapper(grids_mgr);
  }
}

void AtmosphereOutput::
build_remapper (const std::shared_ptr<const gm_type>& grids_mgr)
{
  EKAT_REQUIRE_MSG (grids_mgr,
      "Error! Cannot build a remapper without a valid grids manager.\n");

  auto fm_grid = m_field_mgr->get_grid();
  auto unique_grid = fm_grid->get_unique_grid();

  // The fm is not on a unique grid. We need to create a helper fm,
  // with all the import fields defined on the unique grid,
  // set up scorpio with the helper fm, then create a remapper
  // from the input fm grid to the unique grid.
  EKAT_REQUIRE_MSG (unique_grid,
      "Error! The grid stored in the input field manager is not unique,\n"
      "       and is not storing a unique grid.\n"
      "    grid name: " + fm_grid->name());

  // Create the remapper and register fields from the fid's taken from the input fm
  m_remapper = grids_mgr->create_remapper(fm_grid,unique_grid);
  m_remapper->registration_begins();
  for (const auto& fname : m_fields_names) {
    auto f = m_field_mgr->get_field(fname);
    const auto& src_fid = f.get_header().get_identifier();

    m_remapper->register_field_from_src(src_fid);
  }
  m_remapper->registration_ends();

  // Now register the fields in the unique_fm. To generate their
  // field identifiers, use the tgt fids from the remapper
  auto unique_fm = std::make_shared<fm_type>(unique_grid);
  unique_fm->registration_begins();
  for (int i=0; i<m_remapper->get_num_fields(); ++i) {
    const auto& tgt_fid = m_remapper->get_tgt_field_id(i);
    auto f = m_field_mgr->get_field(tgt_fid.name());

    const int ps = f.get_header().get_alloc_properties().get_largest_pack_size();
    FieldRequest tgt_freq(tgt_fid,ps);
    unique_fm->register_field(tgt_freq);
  }
  unique_fm->registration_ends();

  // Finally, bind the src/tgt fields in the remapper
  // While at it, set the same time stamp of the new field
  // (otherwise an error is thrown during the run() call)
  for (const auto& fname : m_fields_names) {
    auto src = m_field_mgr->get_field(fname);
    auto tgt = unique_fm->get_field(fname);
    tgt.get_header().get_tracking().update_time_stamp(src.get_header().get_tracking().get_time_stamp());
    m_remapper->bind_field(src,tgt);
  }

  // This should never fail, but just in case
  EKAT_REQUIRE_MSG (m_remapper->get_num_fields()==m_remapper->get_num_bound_fields(),
      "Error! Something went wrong while building the scorpio input remapper.\n");

  // Invalidate field mgr ptr, so the checks in the following call don't' throw
  m_field_mgr = nullptr;

  // Replace the field mgr with the unique one
  // NOTE: this call to set_field_manager should not have to call build_remapper anymore.
  //       Therefore, pass nullptr as grids manager, in case there's some bug and this function
  //       gets called again (recursively). A nullptr will cause this function to crap out.
  set_field_manager(unique_fm,nullptr);
}

void AtmosphereOutput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (not m_grid, "Error! Grid pointer was already set.\n");
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");

  auto io_grid = grid;
  if (not grid->is_unique()) {
    io_grid = grid->get_unique_grid();
    EKAT_REQUIRE_MSG (io_grid,
        "Error! Input grid wasn't unique, and did not store a unique grid.\n");
  }
  EKAT_REQUIRE_MSG (io_grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");
  EKAT_REQUIRE_MSG (
      (io_grid->get_global_max_dof_gid()-io_grid->get_global_min_dof_gid()+1)==io_grid->get_num_global_dofs(),
      "Error! In order for IO to work, the grid must (globally) have dof gids in interval [gid_0,gid_0+num_global_dofs).\n");

  EKAT_REQUIRE_MSG(m_comm.size()<=io_grid->get_num_global_dofs(),
      "Error! PIO interface requires the size of the IO MPI group to be\n"
      "       no greater than the global number of columns.\n"
      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_grid = io_grid;
}

void AtmosphereOutput::register_dimensions(const std::string& name)
{
/* 
 * Checks that the dimensions associated with a specific variable will be registered with IO file.
 * INPUT:
 *   field_manager: is a pointer to the field_manager for this simulation.
 *   name: is a string name of the variable who is to be added to the list of variables in this IO stream.
 */
  using namespace scorpio;

  // Store the field layout
  const auto& fid = m_field_mgr->get_field(name).get_header().get_identifier();
  const auto& layout = fid.get_layout();
  m_layouts.emplace(name,layout);

  // Now check taht all the dims of this field are already set to be registered.
  for (int i=0; i<layout.rank(); ++i) {
    // check tag against m_dims map.  If not in there, then add it.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    const auto tag_name = get_nc_tag_name(tags[i],dims[i]);
    auto tag_loc = m_dims.find(tag_name);
    if (tag_loc == m_dims.end()) {
      int tag_len = 0;
      if(e2str(tags[i]) == "COL") {
        // Note: This is because only cols are decomposed over mpi ranks. 
        //       In this case, we need the GLOBAL number of cols.
        tag_len = m_grid->get_num_global_dofs();
      } else {
        tag_len = layout.dim(i);
      }
      m_dims.emplace(std::make_pair(get_nc_tag_name(tags[i],dims[i]),tag_len));
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag_name)==dims[i] or e2str(tags[i])=="COL",
        "Error! Dimension " + tag_name + " on field " + name + " has conflicting lengths");
    }
  }
} // register_dimensions
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_views()
{
  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);

    // These local views are really only needed if the averaging time is not 'Instant',
    // to store running tallies for the average operation. However, we create them
    // also for Instant avg_type, for simplicity later on.

    // If we have an 'Instant' avg type, we can alias the tmp views with the
    // host view of the field, provided that the field does not have padding,
    // and that it is not a subfield of another field (or else the view
    // would be strided).
    bool can_alias_field_view =
        m_avg_type=="Instant" &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent().expired();

    const auto size = m_layouts.at(name).size();
    if (can_alias_field_view) {
      // Alias field's data, to save storage.
      m_host_views_1d.emplace(name,view_1d_host(field.get_internal_view_data<Host>(),size));
    } else {
      // Create a local host view.
      m_host_views_1d.emplace(name,view_1d_host("",size));
    }
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_variables(const std::string& filename)
{
  using namespace scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    auto& fid  = field.get_header().get_identifier();
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    const auto& layout = fid.get_layout();
    for (int i=0; i<fid.get_layout().rank(); ++i) {
      const auto tag_name = get_nc_tag_name(layout.tag(i), layout.dim(i));
      io_decomp_tag += "-" + tag_name; // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }
    // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger? 
    // Should we register dimension variables (such as ncol and lat/lon) elsewhere
    // in the dimension registration?  These won't have time. 
    io_decomp_tag += "-time";
    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to fully C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.

     // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.
    register_variable(filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);
  }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<int> AtmosphereOutput::get_var_dof_offsets(const FieldLayout& layout)
{
  std::vector<int> var_dof(layout.size());

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // Since we order the global array based on dof gid, and we *assume* (we actually
  // check this during set_grid) that the grid global gids are in the interval
  // [gid_0, gid_0+num_global_dofs), the offset is simply given by
  // (dof_gid-gid_0)*column_size (for partitioned arrays).
  // NOTE: we allow gid_0!=0, so that we don't have to worry about 1-based numbering
  //       vs 0-based numbering. The key feature is that the global gids are a
  //       contiguous array. The starting point doesn't matter.
  if (layout.has_tag(ShortFieldTagsNames::COL)) {
    const int num_cols = m_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smalle.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    int col_size = layout.size() / num_cols;

    auto dofs = m_grid->get_dofs_gids();
    auto dofs_h = Kokkos::create_mirror_view(dofs);
    Kokkos::deep_copy(dofs_h,dofs);

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_grid->get_global_min_dof_gid();
    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      auto offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else {
    // This field is *not* defined over columns, so it is not partitioned.
    std::iota(var_dof.begin(),var_dof.end(),0);
  } 

  return var_dof; 
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_degrees_of_freedom(const std::string& filename)
{
  using namespace scorpio;
  using namespace ShortFieldTagsNames;

  // Cycle through all fields and set dof.
  for (auto const& name : m_fields_names) {
    auto field = m_field_mgr->get_field(name);
    const auto& fid  = field.get_header().get_identifier();
    auto var_dof = get_var_dof_offsets(fid.get_layout());
    set_dof(filename,name,var_dof.size(),var_dof.data());
    m_dofs.emplace(std::make_pair(name,var_dof.size()));
  }

  /* TODO: 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
void AtmosphereOutput::setup_output_file(const std::string& filename)
{
  using namespace scream::scorpio;

  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    register_dimension(filename,it.first,it.first,it.second);
  }

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename);

  // Set the offsets of the local dofs in the global vector.
  set_degrees_of_freedom(filename);
}

} // namespace scream
