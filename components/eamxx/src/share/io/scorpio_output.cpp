#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/scream_array_utils.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include "share/util/scream_timing.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

#include <numeric>
#include <fstream>

namespace scream
{

// This helper function updates the current output val with a new one,
// according to the "averaging" type, and according to the number of
// model time steps since the last output step.
KOKKOS_INLINE_FUNCTION
void combine (const Real& new_val, Real& curr_val, const OutputAvgType avg_type)
{
  switch (avg_type) {
    case OutputAvgType::Instant:
      curr_val = new_val;
      break;
    case OutputAvgType::Max:
      curr_val = ekat::impl::max(curr_val,new_val);
      break;
    case OutputAvgType::Min:
      curr_val = ekat::impl::min(curr_val,new_val);
      break;
    case OutputAvgType::Average:
      curr_val += new_val;
      break;
    default:
      EKAT_KERNEL_ERROR_MSG ("Unexpected value for m_avg_type. Please, contact developers.\n");
  }
}

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm,
                  const std::vector<Field>& fields,
                  const std::shared_ptr<const grid_type>& grid)
 : m_comm (comm)
{
  // This version of AtmosphereOutput is for quick output of fields
  m_avg_type = OutputAvgType::Instant;
  m_add_time_dim = false;

  // Create a FieldManager with the input fields
  auto fm = std::make_shared<FieldManager> (grid);
  fm->registration_begins();
  fm->registration_ends();
  for (auto f : fields) {
    fm->add_field(f);
  }

  set_field_manager (fm,"sim");

  for (auto f : fields) {
    m_fields_names.push_back(f.name());
  }

  set_grid (grid);
  set_field_manager (fm,"io");

  // Setup I/O structures
  init ();
}


AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const fm_type>& field_mgr,
                  const std::shared_ptr<const gm_type>& grids_mgr)
 : m_comm         (comm)
 , m_add_time_dim (true)
{
  using vos_t = std::vector<std::string>;

  // Figure out what kind of averaging is requested
  auto avg_type = params.get<std::string>("Averaging Type");
  m_avg_type = str2avg(avg_type);
  EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: Instant, Max, Min, Average. Case insensitive.\n");

  // Set all internal field managers to the simulation field manager to start with.  If
  // vertical remapping, horizontal remapping or both are used then those remapper will
  // set things accordingly.
  set_field_manager (field_mgr,{"sim","io","int"});

  // By default, IO is done directly on the field mgr grid
  m_grids_manager = grids_mgr;
  std::shared_ptr<const grid_type> fm_grid, io_grid;
  io_grid = fm_grid = field_mgr->get_grid();
  if (params.isParameter("Field Names")) {
    // This simple parameter list option does *not* allow to remap fields
    // to an io grid different from that of the field manager. In order to
    // use that functionality, you need the full syntax
    m_fields_names = params.get<vos_t>("Field Names");
  } else if (params.isSublist("Fields")){
    const auto& f_pl = params.sublist("Fields");
    const auto& io_grid_aliases = io_grid->aliases();
    bool grid_found = false;
    for (const auto& grid_name : io_grid_aliases) {
      if (f_pl.isSublist(grid_name)) {
        grid_found = true;
        const auto& pl = f_pl.sublist(grid_name);
        if (pl.isType<vos_t>("Field Names")) {
          m_fields_names = pl.get<vos_t>("Field Names");
        } else if (pl.isType<std::string>("Field Names")) {
          m_fields_names.resize(1, pl.get<std::string>("Field Names"));
          if (m_fields_names[0]=="NONE") {
            m_fields_names.clear();
          }
        }

        // Check if the user wants to remap fields on a different grid first
        if (pl.isParameter("IO Grid Name")) {
          io_grid = grids_mgr->get_grid(pl.get<std::string>("IO Grid Name"));
        }
        break;
      }
    }
    EKAT_REQUIRE_MSG (grid_found,
        "Error! Bad formatting of output yaml file. Missing 'Fields->$grid_name` sublist.\n");
  }

  // Check if remapping and if so create the appropriate remapper 
  // Note: We currently support three remappers
  //   - vertical remapping from file
  //   - horizontal remapping from file
  //   - online remapping which is setup using the create_remapper function
  const bool use_vertical_remap_from_file = params.isParameter("vertical_remap_file");
  const bool use_horiz_remap_from_file = params.isParameter("horiz_remap_file");
  const bool use_online_remapper = io_grid->name()!=fm_grid->name();  // TODO: QUESTION, Do we anticipate online remapping w/ horiz_remap_from file?
  // Check that we are not requesting online remapping w/ horiz and/or vertical remapping.  Which is not currently supported.
  if (use_online_remapper) {
    EKAT_REQUIRE_MSG(!use_vertical_remap_from_file and !use_horiz_remap_from_file,"ERROR: scorpio_output - online remapping not supported with vertical and/or horizontal remapping from file");
  }

  // Try to set the IO grid (checks will be performed)
  set_grid (io_grid);

  // Register any diagnostics needed by this output stream
  set_diagnostics();

  // Setup remappers - if needed
  if (use_vertical_remap_from_file) {  
    using namespace ShortFieldTagsNames;
    // We build a remapper, to remap fields from the fm grid to the io grid
    auto vert_remap_file   = params.get<std::string>("vertical_remap_file");
    auto f_lev = get_field("p_mid","sim");
    auto f_ilev = get_field("p_int","sim");
    m_vert_remapper = std::make_shared<VerticalRemapper>(io_grid,vert_remap_file,f_lev,f_ilev); // We use the default mask value 
    io_grid = m_vert_remapper->get_tgt_grid();
    set_grid(io_grid);

    // Register all output fields in the remapper.
    m_vert_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      auto f = get_field(fname,"sim");
      const auto& src_fid = f.get_header().get_identifier();
      EKAT_REQUIRE_MSG(src_fid.data_type()==DataType::RealType,
          "Error! I/O supports only Real data, for now.\n");
      m_vert_remapper->register_field_from_src(src_fid);
    }
    m_vert_remapper->registration_ends();

    // Now create a new FM on io grid, and create copies of output fields from FM.
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (int i=0; i<m_vert_remapper->get_num_fields(); ++i) {
      const auto& tgt_fid = m_vert_remapper->get_tgt_field_id(i);
      const auto  name    = tgt_fid.name();
      const auto  src     = get_field(name,"sim");
      const auto  packsize = src.get_header().get_alloc_properties().get_largest_pack_size();
      io_fm->register_field(FieldRequest(tgt_fid,packsize)); 
    }
    // Now end registration of variables in io_fm
    io_fm->registration_ends();

    // Now that fields have been allocated on the io grid, we can bind them in the remapper
    for (const auto& fname : m_fields_names) {
      auto src = get_field(fname,"sim");
      auto tgt = io_fm->get_field(src.name());
      m_vert_remapper->bind_field(src,tgt);
    }
    
    // Set the field manager for IO and add an entry "after vertical remapping"
    set_field_manager(io_fm,"io");
    set_field_manager(io_fm,"after_vertical_remap");

    // This should never fail, but just in case
    EKAT_REQUIRE_MSG (m_vert_remapper->get_num_fields()==m_vert_remapper->get_num_bound_fields(),
        "Error! Something went wrong while building the scorpio input remapper.\n");

  }

  // Online remapper and horizontal remapper follow a similar pattern so we check in the same conditional.
  if (use_online_remapper || use_horiz_remap_from_file) {

    if (use_vertical_remap_from_file) {
      // Then we set the horizontal remapper fm input to the one after vertical remapping
      const auto fm_tmp = get_field_manager("after_vertical_remap");
      set_field_manager(fm_tmp,"before_horizontal_remap");
    } else {
      // Then the input FM for horizontal remapping is the simulation remapper.
      const auto fm_tmp = get_field_manager("sim");
      set_field_manager(fm_tmp,"before_horizontal_remap");
    }
    // We build a remapper, to remap fields from the fm grid to the io grid
    // We may need to track masked fields in the horizontal remapper, so we
    // declare the need variables here.
    std::vector<Field>        mask_fields;
    if (use_horiz_remap_from_file) {
      std::map<std::string,int> mask_map;
      // First we check if we need to support masking
      for (const auto& fname : m_fields_names) {
        auto f = get_field(fname,"before_horizontal_remap");
        auto f_extra = f.get_header().get_extra_data();
        if (f_extra.count("mask_data")) {
          // Then this field has a mask attached to it.
          // We only want to tally up a unique set of masks, some fields
          // may share a mask.  We check if a specific mask has already
          // been added.
          auto f_mask = ekat::any_cast<Field>(f_extra.at("mask_data"));
          bool found  = false;
          for (int ii=0; ii<mask_fields.size(); ii++) {
            auto fi = mask_fields[ii];
            if (f_mask.equivalent(fi)) {
              // Then we have already registered an equivalent mask so don't need to do it again.
              mask_map.emplace(fname,ii);
              found = true;
              break;
            }
          }
          // If we haven't already add this mask, we do that now. 
          if (!found) {
            mask_fields.push_back(f_mask);
            mask_map.emplace(f.name(),mask_fields.size()-1);
          }
        } else {
          // No mask attached, give this field an index of -1 to flag it as unmasked
          mask_map.emplace(f.name(),-1);
        }
      }
      // Construct the coarsening remapper
      auto horiz_remap_file   = params.get<std::string>("horiz_remap_file");
      m_horiz_remapper = std::make_shared<CoarseningRemapper>(io_grid,horiz_remap_file,mask_fields,mask_map);
      io_grid = m_horiz_remapper->get_tgt_grid();
      set_grid(io_grid);
    } else {
      m_horiz_remapper = grids_mgr->create_remapper(fm_grid,io_grid);
    }

    // Register all output fields in the remapper.
    m_horiz_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      auto f = get_field(fname,"before_horizontal_remap");
      const auto& src_fid = f.get_header().get_identifier();
      EKAT_REQUIRE_MSG(src_fid.data_type()==DataType::RealType,
          "Error! I/O supports only Real data, for now.\n");
      m_horiz_remapper->register_field_from_src(src_fid);
    }
    m_horiz_remapper->registration_ends();

    // Now create a new FM on io grid, and create copies of output fields from FM.
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (int i=0; i<m_horiz_remapper->get_num_fields(); ++i) {
      const auto& tgt_fid = m_horiz_remapper->get_tgt_field_id(i);
      io_fm->register_field(FieldRequest(tgt_fid));
    }
    io_fm->registration_ends();

    // Now that fields have been allocated on the io grid, we can bind them in the remapper
    for (const auto& fname : m_fields_names) {
      auto src = get_field(fname,"before_horizontal_remap");
      auto tgt = io_fm->get_field(src.name());
      m_horiz_remapper->bind_field(src,tgt);
    }
    // Bind any masked fields, since they won't be listed in the `m_fields_names`
    for (const auto& src : mask_fields) {
      auto tgt  = io_fm->get_field(src.name());
      m_horiz_remapper->bind_field(src,tgt);
    }

    // This should never fail, but just in case
    EKAT_REQUIRE_MSG (m_horiz_remapper->get_num_fields()==m_horiz_remapper->get_num_bound_fields(),
        "Error! Something went wrong while building the scorpio input remapper.\n");

    // Reset the IO field manager
    set_field_manager(io_fm,"io");
  } 

  // Setup I/O structures
  init ();
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::restart (const std::string& filename)
{
  // Create an input stream on the fly, and init averaging data
  ekat::ParameterList res_params("Input Parameters");
  res_params.set<std::string>("Filename",filename);
  res_params.set("Field Names",m_fields_names);

  AtmosphereInput hist_restart (res_params,m_io_grid,m_host_views_1d,m_layouts);
  hist_restart.read_variables();
  hist_restart.finalize();
  for (auto& it : m_host_views_1d) {
    const auto& name = it.first;
    const auto& host = it.second;
    const auto& dev  = m_dev_views_1d.at(name);
    Kokkos::deep_copy(dev,host);
  }
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
  // If we do INSTANT output, but this is not an write step,
  // we can immediately return
  if (not is_write_step and m_avg_type==OutputAvgType::Instant) {
    return;
  }

  using namespace scream::scorpio;

  // Update all diagnostics, we need to do this before applying the remapper
  // to make sure that the remapped fields are the most up to date.
  // First we reset the diag computed map so that all diags are recomputed.
  m_diag_computed.clear();
  for (auto& it : m_diagnostics) {
    compute_diagnostic(it.first);
  }

  auto apply_remap = [&](const std::shared_ptr<AbstractRemapper> remapper)
  {
    remapper->remap(true);

    for (int i=0; i<remapper->get_num_fields(); ++i) {
      // Need to update the time stamp of the fields on the IO grid,
      // to avoid throwing an exception later
      auto src = remapper->get_src_field(i);
      auto tgt = remapper->get_tgt_field(i);

      auto src_t = src.get_header().get_tracking().get_time_stamp();
      tgt.get_header().get_tracking().update_time_stamp(src_t);
    }
  }; // end apply_remap

  // If needed, remap fields from their grid to the unique grid, for I/O
  if (m_vert_remapper) {
    start_timer("EAMxx::IO::vert_remap");
    apply_remap(m_vert_remapper);
    stop_timer("EAMxx::IO::vert_remap");
  }

  if (m_horiz_remapper) {
    start_timer("EAMxx::IO::horiz_remap");
    apply_remap(m_horiz_remapper);
    stop_timer("EAMxx::IO::horiz_remap");
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields_names) {
    // Get all the info for this field.
    const auto  field = get_field(name,"io");
    const auto& layout = m_layouts.at(name);
    const auto& dims = layout.dims();
    const auto  rank = layout.rank();

    // Safety check: make sure that the field was written at least once before using it.
    EKAT_REQUIRE_MSG (!m_add_time_dim || field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Time-dependent output field '" + name + "' has not been initialized yet\n.");

    const bool is_diagnostic = (m_diagnostics.find(name) != m_diagnostics.end());
    const bool is_aliasing_field_view =
        m_avg_type==OutputAvgType::Instant &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent().expired() &&
        not is_diagnostic;

    // Manually update the 'running-tally' views with data from the field,
    // by combining new data with current avg values.
    // NOTE: this is skipped for instant output, if IO view is aliasing Field view.
    auto view_dev = m_dev_views_1d.at(name);
    auto data = view_dev.data();
    KT::RangePolicy policy(0,layout.size());
    const auto extents = layout.extents();

    auto avg_type = m_avg_type;
    // If the dev_view_1d is aliasing the field device view (must be Instant output),
    // then there's no point in copying from the field's view to dev_view
    if (not is_aliasing_field_view) {
      switch (rank) {
        case 1:
        {
          auto new_view_1d = field.get_view<const Real*,Device>();
          auto avg_view_1d = view_Nd_dev<1>(data,dims[0]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
            combine(new_view_1d(i), avg_view_1d(i),avg_type);
          });
          break;
        }
        case 2:
        {
          auto new_view_2d = field.get_view<const Real**,Device>();
          auto avg_view_2d = view_Nd_dev<2>(data,dims[0],dims[1]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            combine(new_view_2d(i,j), avg_view_2d(i,j),avg_type);
          });
          break;
        }
        case 3:
        {
          auto new_view_3d = field.get_view<const Real***,Device>();
          auto avg_view_3d = view_Nd_dev<3>(data,dims[0],dims[1],dims[2]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            combine(new_view_3d(i,j,k), avg_view_3d(i,j,k),avg_type);
          });
          break;
        }
        case 4:
        {
          auto new_view_4d = field.get_view<const Real****,Device>();
          auto avg_view_4d = view_Nd_dev<4>(data,dims[0],dims[1],dims[2],dims[3]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            combine(new_view_4d(i,j,k,l), avg_view_4d(i,j,k,l),avg_type);
          });
          break;
        }
        case 5:
        {
          auto new_view_5d = field.get_view<const Real*****,Device>();
          auto avg_view_5d = view_Nd_dev<5>(data,dims[0],dims[1],dims[2],dims[3],dims[4]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            combine(new_view_5d(i,j,k,l,m), avg_view_5d(i,j,k,l,m),avg_type);
          });
          break;
        }
        case 6:
        {
          auto new_view_6d = field.get_view<const Real******,Device>();
          auto avg_view_6d = view_Nd_dev<6>(data,dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            combine(new_view_6d(i,j,k,l,m,n), avg_view_6d(i,j,k,l,m,n),avg_type);
          });
          break;
        }
        default:
          EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
      }
    }

    if (is_write_step) {
      if (avg_type==OutputAvgType::Average) {
        // Divide by steps count only when the summation is complete
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
          data[i] /= nsteps_since_last_output;
        });
      }
      // Bring data to host
      auto view_host = m_host_views_1d.at(name);
      Kokkos::deep_copy (view_host,view_dev);
      grid_write_data_array(filename,name,view_host.data(),view_host.size());
    }
  }
} // run

long long AtmosphereOutput::
res_dep_memory_footprint () const {
  long long rdmf = 0;
  const auto sim_field_mgr = get_field_manager("sim");

  // Cycle through all unique field managers in this output,
  // first make a copy so we can grab just unique pointers.
  std::set<std::string> grids;
  for (auto fm : m_field_mgrs) {
    auto field_mgr = fm.second;
    if (field_mgr != sim_field_mgr) {
      const auto& gn = field_mgr->get_grid()->name();
      if (grids.count(gn)>0) {
        continue; // Grid has already been parsed
      }
      grids.insert(gn);
      // This FM is done on a different grid than SIM and hasn't been included in
      // the memory calculation yet.  So we can safely add its footprint
      for (const auto& it : *field_mgr) {
        const auto& fap = it.second->get_header().get_alloc_properties();
        if (fap.is_subfield()) {
          continue;
        }
        rdmf += fap.get_alloc_size();
      }
    }
  }

  const auto io_field_mgr = get_field_manager("io");
  for (const auto& fn : m_fields_names) {
    bool is_diagnostic = (m_diagnostics.find(fn) != m_diagnostics.end());
    bool can_alias_field_view =
        m_avg_type==OutputAvgType::Instant && not is_diagnostic &&
        io_field_mgr->get_field(fn).get_header().get_alloc_properties().get_padding()==0 &&
        io_field_mgr->get_field(fn).get_header().get_parent().expired();

    if (not can_alias_field_view) {
      rdmf += m_dev_views_1d.size()*sizeof(Real);
    }
  }

  return rdmf;
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr, const std::vector<std::string>& modes)
{

  // Sanity checks
  EKAT_REQUIRE_MSG (field_mgr, "Error! Invalid field manager pointer.\n");
  EKAT_REQUIRE_MSG (field_mgr->get_grid(), "Error! Field manager stores an invalid grid pointer.\n");

  for (int ii=0; ii<modes.size(); ii++) {
    const auto mode = modes[ii];
    if (m_field_mgrs.count(mode)) {
      // We must redefine the field manager for this location in the map.
      m_field_mgrs.at(mode) = field_mgr;
    } else {
      m_field_mgrs.emplace(mode, field_mgr);
    }
  }

}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
set_field_manager (const std::shared_ptr<const fm_type>& field_mgr, const std::string& mode)
{
  const std::vector<std::string> modes = {mode};
  set_field_manager(field_mgr,modes);
}
/* ---------------------------------------------------------- */

void AtmosphereOutput::
set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (grid, "Error! Input grid pointer is invalid.\n");
  EKAT_REQUIRE_MSG (grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");
  EKAT_REQUIRE_MSG (
      (grid->get_global_max_dof_gid()-grid->get_global_min_dof_gid()+1)==grid->get_num_global_dofs(),
      "Error! In order for IO to work, the grid must (globally) have dof gids in interval [gid_0,gid_0+num_global_dofs).\n");

//ASD  IS THIS STILL TRUE
//ASD  EKAT_REQUIRE_MSG(m_comm.size()<=grid->get_num_global_dofs(),
//ASD      "Error! PIO interface requires the size of the IO MPI group to be\n"
//ASD      "       no greater than the global number of columns.\n"
//ASD      "       Consider decreasing the size of IO MPI group.\n");

  // The grid is good. Store it.
  m_io_grid = grid;
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
  const auto& fid = get_field(name,"io").get_header().get_identifier();
  const auto& layout = fid.get_layout();
  m_layouts.emplace(name,layout);

  // Now check taht all the dims of this field are already set to be registered.
  for (int i=0; i<layout.rank(); ++i) {
    // check tag against m_dims map.  If not in there, then add it.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    const auto tag_name = get_nc_tag_name(tags[i],dims[i]);
    auto tag_loc = m_dims.find(tag_name);
    auto is_partitioned = m_io_grid->get_partitioned_dim_tag()==tags[i];
    if (tag_loc == m_dims.end()) {
      int tag_len = 0;
      if(tags[i] == m_io_grid->get_partitioned_dim_tag()) {
        // This is the dimension that is partitioned across ranks.
        tag_len = m_io_grid->get_partitioned_dim_global_size();
      } else {
        tag_len = layout.dim(i);
      }
      m_dims[get_nc_tag_name(tags[i],dims[i])] = std::make_pair(tag_len,is_partitioned);
    } else {  
      EKAT_REQUIRE_MSG(m_dims.at(tag_name).first==dims[i] or is_partitioned,
        "Error! Dimension " + tag_name + " on field " + name + " has conflicting lengths");
    }
  }
} // register_dimensions
/* ---------------------------------------------------------- */
void AtmosphereOutput::register_views()
{
  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = get_field(name,"io");
    bool is_diagnostic = (m_diagnostics.find(name) != m_diagnostics.end());

    // These local views are really only needed if the averaging time is not 'Instant',
    // to store running tallies for the average operation. However, we create them
    // also for Instant avg_type, for simplicity later on.

    // If we have an 'Instant' avg type, we can alias the 1d views with the
    // views of the field, provided that the field does not have padding,
    // and that it is not a subfield of another field (or else the view
    // would be strided).
    //
    // We also don't want to alias to a diagnostic output since it could share memory
    // with another diagnostic.
    bool can_alias_field_view =
        m_avg_type==OutputAvgType::Instant &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent().expired() &&
        not is_diagnostic;

    const auto size = m_layouts.at(name).size();
    if (can_alias_field_view) {
      // Alias field's data, to save storage.
      m_dev_views_1d.emplace(name,view_1d_dev(field.get_internal_view_data<Real,Device>(),size));
      m_host_views_1d.emplace(name,view_1d_host(field.get_internal_view_data<Real,Host>(),size));
    } else {
      // Create a local view.
      m_dev_views_1d.emplace(name,view_1d_dev("",size));
      m_host_views_1d.emplace(name,Kokkos::create_mirror(m_dev_views_1d[name]));

    }
  }
  // Initialize the local views
  reset_dev_views();
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
reset_dev_views()
{
  // Reset the local device views depending on the averaging type
  // Init dev view with an "identity" for avg_type
  for (auto const& name : m_fields_names) {
    switch (m_avg_type) {
      case OutputAvgType::Instant:
        // No averaging
        break;
      case OutputAvgType::Max:
        Kokkos::deep_copy(m_dev_views_1d[name],-std::numeric_limits<Real>::infinity());
        break;
      case OutputAvgType::Min:
        Kokkos::deep_copy(m_dev_views_1d[name],std::numeric_limits<Real>::infinity());
        break;
      case OutputAvgType::Average:
        Kokkos::deep_copy(m_dev_views_1d[name],0);
        break;
      default:
        EKAT_ERROR_MSG ("Unrecognized averaging type.\n");
    }
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
register_variables(const std::string& filename,
                   const std::string& fp_precision)
{
  using namespace scorpio;

  // Cycle through all fields and register.
  for (auto const& name : m_fields_names) {
    auto field = get_field(name,"io");
    auto& fid  = field.get_header().get_identifier();
    // Make a unique tag for each decomposition. To reuse decomps successfully,
    // we must be careful to make the tags 1-1 with the intended decomp. Here we
    // use the I/O grid name and its global #DOFs, then append the local
    // dimension data.
    //   We use real here because the data type for the decomp is the one used
    // in the simulation and not the one used in the output file.
    std::string io_decomp_tag = (std::string("Real-") + m_io_grid->name() + "-" +
                                 std::to_string(m_io_grid->get_num_global_dofs()));
    std::vector<std::string> vec_of_dims;
    const auto& layout = fid.get_layout();
    std::string units = to_string(fid.get_units());
    for (int i=0; i<fid.get_layout().rank(); ++i) {
      const auto tag_name = get_nc_tag_name(layout.tag(i), layout.dim(i));
      // Concatenate the dimension string to the io-decomp string
      io_decomp_tag += "-" + tag_name;
      // If tag==CMP, we already attached the length to the tag name
      if (layout.tag(i)!=ShortFieldTagsNames::CMP) {
        io_decomp_tag += "_" + std::to_string(layout.dim(i));
      }
      vec_of_dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }

    // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO,
    // may need to delete this line when switching to fully C++/C implementation.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end());
    if (m_add_time_dim) {
      io_decomp_tag += "-time";
      vec_of_dims.push_back("time");  //TODO: See the above comment on time.
    } else {
      io_decomp_tag += "-notime";
    }

    // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.

    register_variable(filename, name, name, units, vec_of_dims,
                      "real",fp_precision, io_decomp_tag);
  }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<scorpio::offset_t>
AtmosphereOutput::get_var_dof_offsets(const FieldLayout& layout)
{
  std::vector<scorpio::offset_t> var_dof(layout.size());

  // Gather the offsets of the dofs of this variable w.r.t. the *global* array.
  // Since we order the global array based on dof gid, and we *assume* (we actually
  // check this during set_grid) that the grid global gids are in the interval
  // [gid_0, gid_0+num_global_dofs), the offset is simply given by
  // (dof_gid-gid_0)*column_size (for partitioned arrays).
  // NOTE: we allow gid_0!=0, so that we don't have to worry about 1-based numbering
  //       vs 0-based numbering. The key feature is that the global gids are a
  //       contiguous array. The starting point doesn't matter.
  // NOTE: a "dof" in the grid object is not the same as a "dof" in scorpio.
  //       For a SEGrid 3d vector field with (MPI local) layout (nelem,2,np,np,nlev),
  //       scorpio sees nelem*2*np*np*nlev dofs, while the SE grid sees nelem*np*np dofs.
  //       All we need to do in this routine is to compute the offset of all the entries
  //       of the MPI-local array w.r.t. the global array. So long as the offsets are in
  //       the same order as the corresponding entry in the data to be read/written, we're good.
  // NOTE: In the case of regional output this rank may have 0 columns to write, thus, var_dof
  //       should be empty, we check for this special case and return an empty var_dof.
  auto dofs_h = m_io_grid->get_dofs_gids().get_view<const AbstractGrid::gid_type*,Host>();
  if (layout.has_tag(ShortFieldTagsNames::COL)) {
    const int num_cols = m_io_grid->get_num_local_dofs();
    if (num_cols==0) {
      return var_dof;
    }

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      auto offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else if (layout.has_tag(ShortFieldTagsNames::EL)) {
    auto layout2d = m_io_grid->get_2d_scalar_layout();
    const int num_my_elems = layout2d.dim(0);
    const int ngp = layout2d.dim(1);
    const int num_cols = num_my_elems*ngp*ngp;
    if (num_cols==0) {
      return var_dof;
    }

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    // Precompute this *before* the loop, since it involves expensive collectives.
    // Besides, the loop might have different length on different ranks, so
    // computing it inside might cause deadlocks.
    auto min_gid = m_io_grid->get_global_min_dof_gid();
    for (int ie=0,icol=0; ie<num_my_elems; ++ie) {
      for (int igp=0; igp<ngp; ++igp) {
        for (int jgp=0; jgp<ngp; ++jgp,++icol) {
          // Get chunk of var_dof to fill
          auto start = var_dof.begin()+icol*col_size;
          auto end   = start+col_size;

          // Compute start of the column offset, then fill column adding 1 to each entry
          auto gid = dofs_h(icol);
          auto offset = (gid-min_gid)*col_size;
          std::iota(start,end,offset);
    }}}
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
    auto field = get_field(name,"io");
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
void AtmosphereOutput::
setup_output_file(const std::string& filename,
                  const std::string& fp_precision)
{
  using namespace scream::scorpio;

  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    register_dimension(filename,it.first,it.first,it.second.first,it.second.second);
  }

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename,fp_precision);

  // Set the offsets of the local dofs in the global vector.
  set_degrees_of_freedom(filename);
}
/* ---------------------------------------------------------- */
// This routine will evaluate the diagnostics stored in this
// output instance.
void AtmosphereOutput::compute_diagnostic(const std::string& name)
{
  auto skip_diag = m_diag_computed[name];
  if (skip_diag) {
    // Diagnostic already computed, just return
    return;
  }
  const auto& diag = m_diagnostics.at(name);
  // Check if the diagnostics has any dependencies, if so, evaluate
  // them as well.  Needed if a diagnostic relies on another
  // diagnostic.
  for (const auto& dep : m_diag_depends_on_diags.at(name)) {
    compute_diagnostic(dep);
  }
  diag->compute_diagnostic();
  m_diag_computed[name] = true;
}
/* ---------------------------------------------------------- */
// General get_field routine for output.
// This routine will first check if a field is in the local field
// manager.  If not it will next check to see if it is in the list
// of available diagnostics.  If neither of these two options it
// will throw an error.
Field AtmosphereOutput::get_field(const std::string& name, const std::string mode) const
{
  const auto field_mgr = get_field_manager(mode);
  const auto sim_field_mgr = get_field_manager("sim");
  if (field_mgr->has_field(name)) {
    return field_mgr->get_field(name);
  } else if (m_diagnostics.find(name) != m_diagnostics.end() && field_mgr==sim_field_mgr) {
    const auto& diag = m_diagnostics.at(name);
    return diag->get_diagnostic();
  } else if (m_fields_alt_name.find(name) != m_fields_alt_name.end()) {
    return get_field(m_fields_alt_name.at(name),mode);
  } else {
    EKAT_ERROR_MSG ("ERROR::AtmosphereOutput::get_field Field " + name + " not found in " + mode + " field manager or diagnostics list.");
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_diagnostics()
{
  const auto sim_field_mgr = get_field_manager("sim");
  // Create all diagnostics
  for (const auto& fname : m_fields_names) {
    if (!sim_field_mgr->has_field(fname)) {
      create_diagnostic(fname);
    }
  }

  // Set required fields for all diagnostics
  // NOTE: do this *after* creating all diags: in case the required
  //       field of certain diagnostics is itself a diagnostic,
  //       we want to make sure the required ones are all built.
  for (const auto& dd : m_diagnostics) {
    const auto& diag = dd.second;
    for (const auto& req : diag->get_required_field_requests()) {
      const auto& req_field = get_field(req.fid.name(),"sim");
      diag->set_required_field(req_field.get_const());
    }

    // Note: this inits with an invalid timestamp. If by any chance we try to
    //       output the diagnostic without computing it, we'll get an error.
    diag->initialize(util::TimeStamp(),RunType::Initial);
  }
}

void AtmosphereOutput::
create_diagnostic (const std::string& diag_field_name) {
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();

  // Construct a diagnostic by this name
  ekat::ParameterList params;
  std::string diag_name;

  // If the diagnostic is $field@lev$N/$field_bot/$field_top,
  // then we need to set some params
  auto tokens = ekat::split(diag_field_name,'@');
  auto last = tokens.back();

  // FieldAtLevel          follows convention variable@lev_N (where N is some integer)
  // FieldAtPressureLevel follows convention variable@999mb (where 999 is some integer)
  auto lev_and_idx = ekat::split(last,'_');
  auto pos = lev_and_idx[0].find_first_not_of("0123456789");
  auto lev_str = lev_and_idx[0].substr(pos);
  
  if (last=="tom" || last=="bot" || lev_str=="lev") {
    // Diagnostic is a horizontal slice at a specific level
    diag_name = "FieldAtLevel";
    tokens.pop_back();
    auto fname = ekat::join(tokens,"_");
    // If the field is itself a diagnostic, make sure it's built
    if (diag_factory.has_product(fname) and
        m_diagnostics.count(fname)==0) {
      create_diagnostic(fname);
      m_diag_depends_on_diags[diag_field_name].push_back(fname);
    } else {
      m_diag_depends_on_diags[diag_field_name].resize(0);
    }
    auto fid = get_field(fname,"sim").get_header().get_identifier();
    params.set("Field Name", fname);
    params.set("Grid Name",fid.get_grid_name());
    params.set("Field Layout",fid.get_layout());
    params.set("Field Units",fid.get_units());

    // If last is bot or top, will simply use that
    params.set("Field Level", lev_and_idx.back());
  } else if (lev_str=="mb" || lev_str=="hPa" || lev_str=="Pa") {
    // Diagnostic is a horizontal slice at a specific pressure level
    diag_name = "FieldAtPressureLevel";
    auto pres_str = lev_and_idx[0].substr(0,pos);
    auto pres_units = lev_and_idx[0].substr(pos);
    auto pres_level = std::stoi(pres_str);
    // Convert pressure level to Pa, the units of pressure in the simulation
    if (pres_units=="mb" || pres_units=="hPa") {
      pres_level *= 100;
    }
    tokens.pop_back();
    auto fname = ekat::join(tokens,"_");
    // If the field is itself a diagnostic, make sure it's built
    if (diag_factory.has_product(fname) and
        m_diagnostics.count(fname)==0) {
      create_diagnostic(fname);
      m_diag_depends_on_diags[diag_field_name].push_back(fname);
    } else {
      m_diag_depends_on_diags[diag_field_name].resize(0);
    }
    auto fid = get_field(fname,"sim").get_header().get_identifier();
    params.set("Field Name", fname);
    params.set("Grid Name",fid.get_grid_name());
    params.set("Field Layout",fid.get_layout());
    params.set("Field Units",fid.get_units());
    params.set<Real>("Field Target Pressure", pres_level);
  } else {
    diag_name = diag_field_name;
    m_diag_depends_on_diags[diag_field_name].resize(0);
  }

  // Create the diagnostic
  auto diag = diag_factory.create(diag_name,m_comm,params);
  diag->set_grids(m_grids_manager);
  m_diagnostics.emplace(diag_field_name,diag);
  // When using remappers with certain diagnostics the get_field command can be called with both the diagnostic
  // name as saved inside the diagnostic and with the name as it is given in the output control file.  If it is
  // the case that these names don't match we add their pairings to the alternate name map.
  if (diag->name() != diag_field_name) {
    m_fields_alt_name.emplace(diag->name(),diag_field_name);
    m_fields_alt_name.emplace(diag_field_name,diag->name());
  }
}

} // namespace scream
