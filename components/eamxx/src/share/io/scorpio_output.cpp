#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/eamxx_array_utils.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/field/field_utils.hpp"

#include "diagnostics/register_diagnostics.hpp"

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
// This one covers cases where a variable might be masked.
KOKKOS_INLINE_FUNCTION
void combine_and_fill (const Real& new_val, Real& curr_val, const OutputAvgType avg_type, const Real fill_value)
{
  const bool new_fill  = new_val  == fill_value;
  const bool curr_fill = curr_val == fill_value;
  if (curr_fill && new_fill) {
    // Then the value is already set to be filled and the new value doesn't change things.
    return;
  } else if (curr_fill) {
    // Then the current value is filled but the new value will replace that for all cases.
    curr_val = new_val;
  } else {
    switch (avg_type) {
      case OutputAvgType::Instant:
        curr_val = new_val;
        break;
      case OutputAvgType::Max:
        curr_val = new_fill ? curr_val : ekat::impl::max(curr_val,new_val);
        break;
      case OutputAvgType::Min:
        curr_val = new_fill ? curr_val : ekat::impl::min(curr_val,new_val);
        break;
      case OutputAvgType::Average:
        curr_val += (new_fill ? 0.0 : new_val);
        break;
      default:
        EKAT_KERNEL_ERROR_MSG ("Unexpected value for m_avg_type. Please, contact developers.\n");
    }
  }
}

// This helper function is used to make sure that the list of fields in
// m_fields_names is a list of unique strings, otherwise throw an error.
void sort_and_check(std::vector<std::string>& fields)
{
  std::sort(fields.begin(),fields.end());
  const bool hasDuplicates = std::adjacent_find(fields.begin(),fields.end()) != fields.end();
  EKAT_REQUIRE_MSG(!hasDuplicates,"ERROR!!! scorpio_output::check_for_duplicates - One of the output yaml files has duplicate field entries.  Please check");
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
  for (auto f : fields) {
    fm->add_field(f);
  }

  set_field_manager (fm,"sim");

  for (auto f : fields) {
    m_fields_names.push_back(f.name());
  }
  sort_and_check(m_fields_names);

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
  sort_and_check(m_fields_names);

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

  // Avg count only makes sense if we have
  //  - non-instant output
  //  - we have one between:
  //    - vertically remapped output
  //    - field_at_XhPa diagnostic
  // We already set m_track_avg_cnt to true if field_at_XhPa is found in set_diagnostics.
  // Hence, here we only check if vert remap is active

  if (params.isParameter("track_avg_cnt")) {
    // This is to be used for unit testing only, so that we can test avg cnt even
    // if there is no vert remap and no field_at_XhPa diagnostic in the stream
    m_track_avg_cnt = params.get<bool>("track_avg_cnt");
  }
  if (use_vertical_remap_from_file) {
    m_track_avg_cnt = true;
  }
  if (params.isParameter("fill_value")) {
    m_fill_value = static_cast<float>(params.get<double>("fill_value"));
  }
  if (params.isParameter("fill_threshold")) {
    m_avg_coeff_threshold = params.get<Real>("fill_threshold");
  }

  // Helper lambda, to copy io string attributes. This will be used if any
  // remapper is created, to ensure atts set by atm_procs are not lost
  auto transfer_io_str_atts = [&] (const Field& src, Field& tgt) {
    const std::string io_string_atts_key ="io: string attributes";
    using stratts_t = std::map<std::string,std::string>;
    const auto& src_atts = src.get_header().get_extra_data<stratts_t>(io_string_atts_key);
          auto& dst_atts = tgt.get_header().get_extra_data<stratts_t>(io_string_atts_key);
    for (const auto& [name,val] : src_atts) {
      dst_atts[name] = val;
    }
  };

  // Setup remappers - if needed
  if (use_vertical_remap_from_file) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    auto vert_remap_file   = params.get<std::string>("vertical_remap_file");
    auto p_mid = get_field("p_mid","sim");
    auto p_int = get_field("p_int","sim");
    auto vert_remapper = std::make_shared<VerticalRemapper>(io_grid,vert_remap_file);
    vert_remapper->set_source_pressure (p_mid,p_int);
    vert_remapper->set_mask_value(m_fill_value);
    vert_remapper->set_extrapolation_type(VerticalRemapper::Mask); // both Top AND Bot
    m_vert_remapper = vert_remapper;
    io_grid = m_vert_remapper->get_tgt_grid();
    set_grid(io_grid);

    // Now create a new FM on io grid, and create copies of output fields on that grid,
    // using the remapper to get the correct identifier on the tgt grid
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (const auto& fname : m_fields_names) {
      const auto src = get_field(fname,"sim");
      const auto tgt_fid = m_vert_remapper->create_tgt_fid(src.get_header().get_identifier());
      const auto packsize = src.get_header().get_alloc_properties().get_largest_pack_size();
      io_fm->register_field(FieldRequest(tgt_fid,packsize));
    }
    io_fm->registration_ends();
    for (const auto& fname : m_fields_names) {
      const auto& src = get_field(fname,"sim");
            auto& tgt = io_fm->get_field(fname);
      transfer_io_str_atts (src,tgt);
    }

    // Register all output fields in the remapper.
    m_vert_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      const auto src = get_field(fname,"sim");
      const auto tgt = io_fm->get_field(src.name());
      m_vert_remapper->register_field(src,tgt);
    }
    m_vert_remapper->registration_ends();

    // Reet the field manager for IO
    set_field_manager(io_fm,"io");

    // Store a handle to 'after-vremap' FM
    set_field_manager(io_fm,"after_vertical_remap");
  }

  // Online remapper and horizontal remapper follow a similar pattern so we check in the same conditional.
  if (use_online_remapper || use_horiz_remap_from_file) {

    // Whic FM is the one pre-horiz-remap depends on whether we did vert remap or not
    const auto fm_pre_hremap = use_vertical_remap_from_file
                             ? get_field_manager("after_vertical_remap")
                             : get_field_manager("sim");
    set_field_manager(fm_pre_hremap,"before_horizontal_remap");

    // We build a remapper, to remap fields from the fm grid to the io grid
    if (use_horiz_remap_from_file) {
      // Construct the coarsening remapper
      auto horiz_remap_file   = params.get<std::string>("horiz_remap_file");
      m_horiz_remapper = std::make_shared<CoarseningRemapper>(io_grid,horiz_remap_file,true);
      io_grid = m_horiz_remapper->get_tgt_grid();
      set_grid(io_grid);
    } else {
      // Construct a generic remapper (likely, SE->Point)
      m_horiz_remapper = grids_mgr->create_remapper(fm_grid,io_grid);
    }

    // Create a FM on the horiz remapper tgt grid, and register fields on it
    auto io_fm = std::make_shared<fm_type>(io_grid);
    io_fm->registration_begins();
    for (const auto& fname : m_fields_names) {
      const auto src = get_field(fname,"before_horizontal_remap");
      const auto tgt_fid = m_horiz_remapper->create_tgt_fid(src.get_header().get_identifier());
      const auto packsize = src.get_header().get_alloc_properties().get_largest_pack_size();
      io_fm->register_field(FieldRequest(tgt_fid,packsize));
    }
    io_fm->registration_ends();
    for (const auto& fname : m_fields_names) {
      const auto& src = get_field(fname,"before_horizontal_remap");
            auto& tgt = io_fm->get_field(fname);
      transfer_io_str_atts (src,tgt);
    }

    // Register all output fields in the remapper.
    m_horiz_remapper->registration_begins();
    for (const auto& fname : m_fields_names) {
      const auto src = get_field(fname,"before_horizontal_remap");
      const auto tgt = io_fm->get_field(src.name());
      EKAT_REQUIRE_MSG(src.data_type()==DataType::RealType,
          "Error! I/O supports only Real data, for now.\n");
      m_horiz_remapper->register_field(src,tgt);
    }
    m_horiz_remapper->registration_ends();

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
  std::vector<std::string> input_field_names = m_fields_names;
  input_field_names.insert(input_field_names.end(),m_avg_cnt_names.begin(),m_avg_cnt_names.end());
  res_params.set("Field Names",input_field_names);

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
}

void AtmosphereOutput::
init_timestep (const util::TimeStamp& start_of_step)
{
  for (auto& it : m_diagnostics) {
    it.second->init_timestep(start_of_step);
  }
}

void AtmosphereOutput::
run (const std::string& filename,
     const bool output_step, const bool checkpoint_step,
     const int nsteps_since_last_output,
     const bool allow_invalid_fields)
{
  // If we do INSTANT output, but this is not an write step,
  // we can immediately return
  const bool is_write_step = output_step or checkpoint_step;
  if (not is_write_step and m_avg_type==OutputAvgType::Instant) {
    return;
  }
  Real duration_write = 0.0;  // Record of time spent writing output
  if (is_write_step) {
    if (m_atm_logger) {
      m_atm_logger->info("[EAMxx::scorpio_output] Writing variables to file");
      m_atm_logger->info("  file name: " + filename);
    }
  }

  using namespace scream::scorpio;

  // Update all diagnostics, we need to do this before applying the remapper
  // to make sure that the remapped fields are the most up to date.
  // First we reset the diag computed map so that all diags are recomputed.
  m_diag_computed.clear();
  for (auto& it : m_diagnostics) {
    compute_diagnostic(it.first,allow_invalid_fields);
  }

  auto apply_remap = [&](const std::shared_ptr<AbstractRemapper> remapper)
  {
    remapper->remap_fwd();

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

  // Update all of the averaging count views (if needed)
  // The strategy is as follows:
  // For the update to the averaged value for this timestep we need to track if
  // a point in a specific layout is "filled" or not.  So we create a set of local
  // temporary views for each layout that are either 0 or 1 depending on if the
  // value is filled or unfilled.
  // We then use these values to update the overall average count views for that layout.
  if (m_track_avg_cnt) {
    // Note, we assume that all fields that share a layout are also masked/filled in the same
    // way. If we need to handle a case where only a subset of output variables are expected to
    // be masked/filled then the recommendation is to request those variables in a separate output
    // stream.
    // We cycle through all fields and we
    //  1. Find the avg_cnt view for this field.
    //  2. If we already processed the avg_cnt view, go to next field, and start from 1 again.
    //  3. Add 1 to all entries of avg_cnt where field!=fill_value
    std::set<std::string> avg_updated;
    for (const auto& name : m_fields_names) {
      auto avg_cnt_name = m_field_to_avg_cnt_map.at(name);
      if (avg_updated.count(avg_cnt_name)==1) {
        // We updated this avg_cnt by checking another field
        continue;
      }
      auto field = get_field(name,"io");
      update_avg_cnt_view(field,m_dev_views_1d.at(avg_cnt_name));

      // Make sure we don't double update this avg cnt
      avg_updated.insert(avg_cnt_name);
    }
  }

  // Take care of updating and possibly writing fields.
  // These are needed inside kernels, so crate local copies
  auto do_avg_cnt = m_track_avg_cnt;
  auto avg_type = m_avg_type;
  auto fill_value = m_fill_value;
  auto avg_coeff_threshold = m_avg_coeff_threshold;
  for (auto const& name : m_fields_names) {
    // Get all the info for this field.
          auto  field = get_field(name,"io");
    const auto& layout = m_layouts.at(field.name());
    const auto& dims = layout.dims();
    const auto  rank = layout.rank();

    if (not field.get_header().get_tracking().get_time_stamp().is_valid()) {
      // Safety check: make sure that the user is ok with this
      if (allow_invalid_fields) {
        field.deep_copy(m_fill_value);
      } else {
        EKAT_REQUIRE_MSG (!m_add_time_dim,
            "Error! Time-dependent output field '" + name + "' has not been initialized yet\n.");
      }
    }

    const bool is_diagnostic = (m_diagnostics.find(name) != m_diagnostics.end());
    const bool is_aliasing_field_view =
        m_avg_type==OutputAvgType::Instant &&
        field.get_header().get_alloc_properties().get_padding()==0 &&
        field.get_header().get_parent()==nullptr &&
        not is_diagnostic;

    // Manually update the 'running-tally' views with data from the field,
    // by combining new data with current avg values.
    // NOTE: this is skipped for instant output, if IO view is aliasing Field view.
    auto view_dev = m_dev_views_1d.at(name);
    auto data = view_dev.data();
    KT::RangePolicy policy(0,layout.size());
    const auto extents = layout.extents();

    // If the dev_view_1d is aliasing the field device view (must be Instant output),
    // then there's no point in copying from the field's view to dev_view
    if (not is_aliasing_field_view) {
      switch (rank) {
        case 0:
        {
          auto new_view_0d = field.get_view<const Real,Device>();
          auto avg_view_0d = view_Nd_dev<0>(data);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int) {
            if (do_avg_cnt) {
              combine_and_fill(new_view_0d(),avg_view_0d(),avg_type,fill_value);
            } else {
              combine(new_view_0d(),avg_view_0d(),avg_type);
            }
          });
          break;
        }
        case 1:
        {
          // For rank-1 views, we use strided layout, since it helps us
          // handling a few more scenarios
          auto new_view_1d = field.get_strided_view<const Real*,Device>();
          auto avg_view_1d = view_Nd_dev<1>(data,dims[0]);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
            if (do_avg_cnt) {
              combine_and_fill(new_view_1d(i),avg_view_1d(i),avg_type,fill_value);
            } else {
              combine(new_view_1d(i),avg_view_1d(i),avg_type);
            }
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
            if (do_avg_cnt) {
              combine_and_fill(new_view_2d(i,j),avg_view_2d(i,j),avg_type,fill_value);
            } else {
              combine(new_view_2d(i,j), avg_view_2d(i,j),avg_type);
            }
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
            if (do_avg_cnt) {
              combine_and_fill(new_view_3d(i,j,k),avg_view_3d(i,j,k),avg_type,fill_value);
            } else {
              combine(new_view_3d(i,j,k), avg_view_3d(i,j,k),avg_type);
            }
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
            if (do_avg_cnt) {
              combine_and_fill(new_view_4d(i,j,k,l), avg_view_4d(i,j,k,l),avg_type,fill_value);
            } else {
              combine(new_view_4d(i,j,k,l), avg_view_4d(i,j,k,l),avg_type);
            }
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
            if (do_avg_cnt) {
              combine_and_fill(new_view_5d(i,j,k,l,m), avg_view_5d(i,j,k,l,m),avg_type,fill_value);
            } else {
              combine(new_view_5d(i,j,k,l,m), avg_view_5d(i,j,k,l,m),avg_type);
            }
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
            if (do_avg_cnt) {
              combine_and_fill(new_view_6d(i,j,k,l,m,n), avg_view_6d(i,j,k,l,m,n), avg_type,fill_value);
            } else {
              combine(new_view_6d(i,j,k,l,m,n), avg_view_6d(i,j,k,l,m,n),avg_type);
            }
          });
          break;
        }
        default:
          EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by AtmosphereOutput.\n");
      }
    }

    if (is_write_step) {
      if (output_step and avg_type==OutputAvgType::Average) {
        if (do_avg_cnt) {
          const auto avg_cnt_lookup = m_field_to_avg_cnt_map.at(name);
          const auto avg_cnt_view = m_dev_views_1d.at(avg_cnt_lookup);
          const auto avg_nsteps = avg_cnt_view.data();
          // Divide by steps count only when the summation is complete
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
            Real coeff_percentage = Real(avg_nsteps[i])/nsteps_since_last_output;
            if (data[i] != fill_value && coeff_percentage > avg_coeff_threshold) {
              data[i] /= avg_nsteps[i];
            } else {
              data[i] = fill_value;
            }
          });
        } else {
          // Divide by steps count only when the summation is complete
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
            data[i] /= nsteps_since_last_output;
          });
        }
      }
      // Bring data to host
      auto view_host = m_host_views_1d.at(name);
      Kokkos::deep_copy (view_host,view_dev);
      auto func_start = std::chrono::steady_clock::now();
      scorpio::write_var(filename,name,view_host.data());
      auto func_finish = std::chrono::steady_clock::now();
      auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
      duration_write += duration_loc.count();
    }
  }
  // Handle writing the average count variables to file
  if (is_write_step) {
    for (const auto& name : m_avg_cnt_names) {
      auto& view_dev = m_dev_views_1d.at(name);
      // Bring data to host
      auto view_host = m_host_views_1d.at(name);
      Kokkos::deep_copy (view_host,view_dev);
      auto func_start = std::chrono::steady_clock::now();
      scorpio::write_var(filename,name,view_host.data());
      auto func_finish = std::chrono::steady_clock::now();
      auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
      duration_write += duration_loc.count();
    }
  }
  if (is_write_step) {
    if (m_atm_logger) {
      m_atm_logger->info("  Done! Elapsed time: " + std::to_string(duration_write/1000.0) +" seconds");
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
        io_field_mgr->get_field(fn).get_header().get_parent()==nullptr;

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

  for (unsigned ii=0; ii<modes.size(); ii++) {
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

std::shared_ptr<const FieldManager>
AtmosphereOutput::get_field_manager (const std::string& mode) const
{
  auto it = m_field_mgrs.find(mode);
  EKAT_REQUIRE_MSG (it!=m_field_mgrs.end(),
    "ERROR! AtmosphereOutput::get_field_manager FM for mode = " + mode + " not found in list of available field managers!.");
  return it->second;
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
  using namespace ShortFieldTagsNames;

  // Store the field layout
  const auto& fid = get_field(name,"io").get_header().get_identifier();
  const auto& layout = fid.get_layout();
  m_layouts.emplace(fid.name(),layout);

  // Now check taht all the dims of this field are already set to be registered.
  const auto& tags = layout.tags();
  const auto& dims = layout.dims();
  for (int i=0; i<layout.rank(); ++i) {
    // check tag against m_dims map.  If not in there, then add it.
    std::string tag_name = m_io_grid->has_special_tag_name(tags[i])
                         ? m_io_grid->get_special_tag_name(tags[i])
                         : layout.names()[i];

    // If t==CMP, and the name stored in the layout is the default ("dim"),
    // we append also the extent, to allow different vector dims in the file
    tag_name += tag_name=="dim" ? std::to_string(dims[i]) : "";

    auto is_partitioned = m_io_grid->get_partitioned_dim_tag()==tags[i];
    int dim_len = is_partitioned
                ? m_io_grid->get_partitioned_dim_global_size()
                : layout.dim(i);
    auto it_bool = m_dims.emplace(tag_name,dim_len);
    EKAT_REQUIRE_MSG(it_bool.second or it_bool.first->second==dim_len,
      "Error! Dimension " + tag_name + " on field " + name + " has conflicting lengths.\n"
      "  - old length: " + std::to_string(m_dims[tag_name]) + "\n"
      "  - new length: " + std::to_string(dim_len) + "\n"
      "If same name applies to different dims (e.g. PhysicsGLL and PhysicsPG2 define "
      "\"ncol\" at different lengths), reset tag name for one of the grids.\n");
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
        field.get_header().get_parent()==nullptr &&
        not is_diagnostic;

    const auto layout = m_layouts.at(field.name());
    const auto size = layout.size();
    if (can_alias_field_view) {
      // Alias field's data, to save storage.
      m_dev_views_1d.emplace(name,view_1d_dev(field.get_internal_view_data<Real,Device>(),size));
      m_host_views_1d.emplace(name,view_1d_host(field.get_internal_view_data<Real,Host>(),size));
    } else {
      // Create a local view.
      m_dev_views_1d.emplace(name,view_1d_dev("",size));
      m_host_views_1d.emplace(name,Kokkos::create_mirror(m_dev_views_1d[name]));
    }

    if (m_track_avg_cnt) {
      // Now create and store a dev view to track the averaging count for this layout (if we are tracking)
      // We don't need to track average counts for files that are not tracking the time dim
      set_avg_cnt_tracking(name,layout);
    }
  }

  // Initialize the local views
  reset_dev_views();
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_avg_cnt_tracking(const std::string& name, const FieldLayout& layout)
{
  // Make sure this field "name" hasn't already been registered with avg_cnt tracking.
  // Note, we check this because some diagnostics need to have their own tracking which
  // is created at the 'create_diagnostics' function.
  if (m_field_to_avg_cnt_map.count(name)>0) {
    return;
  }

  // If the field is not an output field, do not register the avg count. This can happen
  // if a diag depends on another diag. In this case, the inner diag is never outputed,
  // so we don't want to create an avg count for its layout, since it may contain dims
  // that are not in the list of registered dims (the dims of the output vars).
  // See issue https://github.com/E3SM-Project/scream/issues/2663
  if (not ekat::contains(m_fields_names,name)) {
    return;
  }

  // Now create and store a dev view to track the averaging count for this layout (if we are tracking)
  // We don't need to track average counts for files that are not tracking the time dim
  const auto& avg_cnt_suffix = m_field_to_avg_cnt_suffix[name];
  const auto size = layout.size();
  const auto tags = layout.tags();
  if (m_track_avg_cnt) {
    std::string avg_cnt_name = "avg_count" + avg_cnt_suffix;
    for (int i=0; i<layout.rank(); ++i) {
      const auto t = layout.tag(i);
      std::string tag_name = m_io_grid->has_special_tag_name(t)
                           ? m_io_grid->get_special_tag_name(t)
                           : layout.names()[i];

      // If t==CMP, and the name stored in the layout is the default ("dim"),
      // we append also the extent, to allow different vector dims in the file
      tag_name += tag_name=="dim" ? std::to_string(layout.dim(i)) : "";

      avg_cnt_name += "_" + tag_name;
    }
    if (std::find(m_avg_cnt_names.begin(),m_avg_cnt_names.end(),avg_cnt_name)==m_avg_cnt_names.end()) {
      m_avg_cnt_names.push_back(avg_cnt_name);
    }
    m_field_to_avg_cnt_map.emplace(name,avg_cnt_name);
    m_dev_views_1d.emplace(avg_cnt_name,view_1d_dev("",size));  // Note, emplace will only add a new key if one isn't already there
    m_host_views_1d.emplace(avg_cnt_name,Kokkos::create_mirror(m_dev_views_1d[avg_cnt_name]));
    m_layouts.emplace(avg_cnt_name,layout);
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
reset_dev_views()
{
  // Reset the local device views depending on the averaging type
  // Init dev view with an "identity" for avg_type
  const Real fill_for_average = m_track_avg_cnt ? m_fill_value : 0.0;
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
        Kokkos::deep_copy(m_dev_views_1d[name],fill_for_average);
        break;
      default:
        EKAT_ERROR_MSG ("Unrecognized averaging type.\n");
    }
  }
  // Reset all views for averaging count to 0
  for (auto const& name : m_avg_cnt_names) {
    Kokkos::deep_copy(m_dev_views_1d[name],0);
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
register_variables(const std::string& filename,
                   const std::string& fp_precision,
                   const scorpio::FileMode mode)
{
  using namespace ShortFieldTagsNames;
  using strvec_t = std::vector<std::string>;

  EKAT_REQUIRE_MSG (ekat::contains(strvec_t{"float","single","double","real"},fp_precision),
      "Error! Invalid/unsupported value for fp_precision.\n"
      "  - input value: " + fp_precision + "\n"
      "  - supported values: float, single, double, real\n");

  // Helper lambdas
  auto set_vec_of_dims = [&](const FieldLayout& layout) {
    std::vector<std::string> vec_of_dims;
    for (int i=0; i<layout.rank(); ++i) {
      const auto t = layout.tag(i);
      auto tag_name = m_io_grid->has_special_tag_name(t)
                    ? m_io_grid->get_special_tag_name(t)
                    : layout.names()[i];
      if (tag_name=="dim") {
        tag_name += std::to_string(layout.dim(i));
      }
      vec_of_dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }
    return vec_of_dims;
  };

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
    const auto& layout = fid.get_layout();
    auto vec_of_dims   = set_vec_of_dims(layout);
    std::string units = fid.get_units().to_string();

    // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.

    if (mode==scorpio::FileMode::Append) {
      // Simply check that the var is in the file, and has the right properties
      EKAT_REQUIRE_MSG (scorpio::has_var(filename,name),
          "Error! Cannot append, due to variable missing from the file.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n");
      const auto& var = scorpio::get_var(filename,name);
      EKAT_REQUIRE_MSG (var.dim_names()==vec_of_dims,
          "Error! Cannot append, due to variable dimensions mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n"
          "  - var dims : " + ekat::join(vec_of_dims,",") + "\n"
          "  - var dims from file: " + ekat::join(var.dim_names(),",") + "\n");
      EKAT_REQUIRE_MSG (var.units==units,
          "Error! Cannot append, due to variable units mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n"
          "  - var units: " + units + "\n"
          "  - var units from file: " + var.units + "\n");
      EKAT_REQUIRE_MSG (var.time_dep==m_add_time_dim,
          "Error! Cannot append, due to time dependency mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n"
          "  - var time dep: " + (m_add_time_dim ? "yes" : "no") + "\n"
          "  - var time dep from file: " + (var.time_dep ? "yes" : "no") + "\n");
    } else {
      scorpio::define_var (filename, name, units, vec_of_dims,
                            "real",fp_precision, m_add_time_dim);

      // Add FillValue as an attribute of each variable
      // FillValue is a protected metadata, do not add it if it already existed
      if (fp_precision=="double" or
          (fp_precision=="real" and std::is_same<Real,double>::value)) {
        double fill_value = m_fill_value;
        scorpio::set_attribute(filename, name, "_FillValue",fill_value);
      } else {
        float fill_value = m_fill_value;
        scorpio::set_attribute(filename, name, "_FillValue",fill_value);
      }

      // If this is has subfields, add list of its children
      const auto& children = field.get_header().get_children();
      if (children.size()>0) {
        // This field is a parent to a set of subfields
        std::string children_list;
        children_list += "[ ";
        for (const auto& ch_w : children) {
          auto child = ch_w.lock();
          children_list += child->get_identifier().name() + ", ";
        }
        // Replace last "," with "]"
        children_list.pop_back();
        children_list.pop_back();
        children_list += " ]";
        scorpio::set_attribute(filename,name,"sub_fields",children_list);
      }

      // If tracking average count variables then add the name of the tracking variable for this variable
      if (m_track_avg_cnt) {
        const auto lookup = m_field_to_avg_cnt_map.at(name);
        scorpio::set_attribute(filename,name,"averaging_count_tracker",lookup);
      }

      // Atm procs may have set some request for metadata.
      using stratts_t = std::map<std::string,std::string>;
      const auto& str_atts = field.get_header().get_extra_data<stratts_t>("io: string attributes");
      for (const auto& [att_name,att_val] : str_atts) {
        scorpio::set_attribute(filename,name,att_name,att_val);
      }

      // Gather longname (if not already in the io: string attributes)
      if (str_atts.count("long_name")==0) {
        auto longname = m_default_metadata.get_longname(name);
        scorpio::set_attribute(filename, name, "long_name", longname);
      }

      // Gather standard name, CF-Compliant (if not already in the io: string attributes)
      if (str_atts.count("standard_name")==0) {
        auto standardname = m_default_metadata.get_standardname(name);
        scorpio::set_attribute(filename, name, "standard_name", standardname);
      }
    }
  }
  // Now register the average count variables
  if (m_track_avg_cnt) {
    std::string unitless = "unitless";
    for (const auto& name : m_avg_cnt_names) {
      const auto layout = m_layouts.at(name);
      auto vec_of_dims   = set_vec_of_dims(layout);
      if (mode==scorpio::FileMode::Append) {
        // Similar to the regular fields above, check that the var is in the file, and has the right properties
        EKAT_REQUIRE_MSG (scorpio::has_var(filename,name),
            "Error! Cannot append, due to variable missing from the file.\n"
            "  - filename : " + filename + "\n"
            "  - varname  : " + name + "\n");
        const auto& var = scorpio::get_var(filename,name);
        EKAT_REQUIRE_MSG (var.dim_names()==vec_of_dims,
            "Error! Cannot append, due to variable dimensions mismatch.\n"
            "  - filename : " + filename + "\n"
            "  - varname  : " + name + "\n"
            "  - var dims : " + ekat::join(vec_of_dims,",") + "\n"
            "  - var dims from file: " + ekat::join(var.dim_names(),",") + "\n");
        EKAT_REQUIRE_MSG (var.units==unitless,
            "Error! Cannot append, due to variable units mismatch.\n"
            "  - filename : " + filename + "\n"
            "  - varname  : " + name + "\n"
            "  - var units: " + unitless + "\n"
            "  - var units from file: " + var.units + "\n");
        EKAT_REQUIRE_MSG (var.time_dep==m_add_time_dim,
            "Error! Cannot append, due to time dependency mismatch.\n"
            "  - filename : " + filename + "\n"
            "  - varname  : " + name + "\n"
            "  - var time dep: " + (m_add_time_dim ? "yes" : "no") + "\n"
            "  - var time dep from file: " + (var.time_dep ? "yes" : "no") + "\n");
      } else {
	// Note, unlike with regular output variables, for the average counting
	// variables we don't need to add all of the extra metadata.  So we simply
	// define the variable.
        scorpio::define_var(filename, name, unitless, vec_of_dims,
                            "real",fp_precision, m_add_time_dim);
      }
    }
  }
} // register_variables
/* ---------------------------------------------------------- */
std::vector<scorpio::offset_t>
AtmosphereOutput::get_var_dof_offsets(const FieldLayout& layout)
{
  using namespace ShortFieldTagsNames;

  // Precompute this *before* the early return, since it involves collectives.
  // If one rank owns zero cols, and returns prematurely, the others will be left waiting.
  AbstractGrid::gid_type min_gid = -1;
  if (layout.has_tag(COL) or layout.has_tag(EL)) {
    min_gid = m_io_grid->get_global_min_dof_gid();
  }

  // It may be that this MPI rank owns no chunk of the field
  if (layout.size()==0) {
    return {};
  }

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
  if (layout.has_tag(COL)) {
    const int num_cols = m_io_grid->get_num_local_dofs();

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

    for (int icol=0; icol<num_cols; ++icol) {
      // Get chunk of var_dof to fill
      auto start = var_dof.begin()+icol*col_size;
      auto end   = start+col_size;

      // Compute start of the column offset, then fill column adding 1 to each entry
      auto gid = dofs_h(icol);
      auto offset = (gid-min_gid)*col_size;
      std::iota(start,end,offset);
    }
  } else if (layout.has_tag(EL)) {
    auto layout2d = m_io_grid->get_2d_scalar_layout();
    const int num_my_elems = layout2d.dim(0);
    const int ngp = layout2d.dim(1);
    const int num_cols = num_my_elems*ngp*ngp;

    // Note: col_size might be *larger* than the number of vertical levels, or even smaller.
    //       E.g., (ncols,2,nlevs), or (ncols,2) respectively.
    scorpio::offset_t col_size = layout.size() / num_cols;

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

void AtmosphereOutput::set_decompositions(const std::string& filename)
{
  using namespace ShortFieldTagsNames;

  // First, check if any of the vars is indeed partitioned
  const auto decomp_tag  = m_io_grid->get_partitioned_dim_tag();

  bool has_decomposed_layouts = false;
  for (const auto& it : m_layouts) {
    if (it.second.has_tag(decomp_tag)) {
      has_decomposed_layouts = true;
      break;
    }
  }
  if (not has_decomposed_layouts) {
    // If none of the vars are decomposed on this grid,
    // then there's nothing to do here
    return;
  } 

  // Set the decomposition for the partitioned dimension
  const int local_dim = m_io_grid->get_partitioned_dim_local_size();
  std::string decomp_dim = m_io_grid->has_special_tag_name(decomp_tag)
                         ? m_io_grid->get_special_tag_name(decomp_tag)
                         : e2str(decomp_tag);
  auto gids_f = m_io_grid->get_partitioned_dim_gids();
  auto gids_h = gids_f.get_view<const AbstractGrid::gid_type*,Host>();
  auto min_gid = m_io_grid->get_global_min_partitioned_dim_gid();
  std::vector<scorpio::offset_t> offsets(local_dim);
  for (int idof=0; idof<local_dim; ++idof) {
    offsets[idof] = gids_h[idof] - min_gid;
  }
  scorpio::set_dim_decomp(filename,decomp_dim,offsets);
}

void AtmosphereOutput::
setup_output_file(const std::string& filename,
                  const std::string& fp_precision,
                  const scorpio::FileMode mode)
{
  // Register dimensions with netCDF file.
  for (auto it : m_dims) {
    if (mode==scorpio::FileMode::Append) {
      // Simply check that the dim is in the file, and has the right extent
      EKAT_REQUIRE_MSG (scorpio::has_dim(filename,it.first),
          "Error! Cannot append, due to missing dim in the file.\n"
          "  - filename: " + filename + "\n"
          "  - dimname : " + it.first + "\n");
      EKAT_REQUIRE_MSG (scorpio::get_dimlen(filename,it.first)==it.second,
          "Error! Cannot append, due to mismatch dim length.\n"
          "  - filename: " + filename + "\n"
          "  - dimname : " + it.first + "\n"
          "  - old len : " + std::to_string(scorpio::get_dimlen(filename,it.first)) + "\n"
          "  - new len : " + std::to_string(it.second) + "\n");
    } else {
      scorpio::define_dim(filename,it.first,it.second);
    }
  }

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename,fp_precision,mode);

  // Set the offsets of the local dofs in the global vector.
  set_decompositions(filename);
}
/* ---------------------------------------------------------- */
// This routine will evaluate the diagnostics stored in this
// output instance.
void AtmosphereOutput::
compute_diagnostic(const std::string& name, const bool allow_invalid_fields)
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
    compute_diagnostic(dep,allow_invalid_fields);
  }

  m_diag_computed[name] = true;
  if (allow_invalid_fields) {
    // If any input is invalid, fill the diagnostic with invalid data
    for (auto f : diag->get_fields_in()) {
      if (not f.get_header().get_tracking().get_time_stamp().is_valid()) {
        // Fill diag with invalid data and return
        diag->get_diagnostic().deep_copy(m_fill_value);
        return;
      }
    }
  }

  // Either allow_invalid_fields=false, or all inputs are valid. Proceed.
  diag->compute_diagnostic();

  // The diag may have failed to compute (e.g., t=0 output with a flux-like diag).
  // If we're allowing invalid fields, then we should simply set diag=m_fill_value
  if (allow_invalid_fields) {
    auto d = diag->get_diagnostic();
    if (not d.get_header().get_tracking().get_time_stamp().is_valid()) {
      d.deep_copy(m_fill_value);
    }
  }
}
/* ---------------------------------------------------------- */
// General get_field routine for output.
// This routine will first check if a field is in the local field
// manager.  If not it will next check to see if it is in the list
// of available diagnostics.  If neither of these two options it
// will throw an error.
Field AtmosphereOutput::
get_field(const std::string& name, const std::string& mode) const
{
  const auto field_mgr = get_field_manager(mode);
  const auto sim_field_mgr = get_field_manager("sim");
  const bool can_be_diag = field_mgr == sim_field_mgr;
  if (field_mgr->has_field(name)) {
    return field_mgr->get_field(name);
  } else if (m_diagnostics.find(name) != m_diagnostics.end() && can_be_diag) {
    const auto& diag = m_diagnostics.at(name);
    return diag->get_diagnostic();
  } else {
    EKAT_ERROR_MSG ("ERROR::AtmosphereOutput::get_field Field " + name + " not found in " + mode + " field manager or diagnostics list.");
  }
  static Field f;
  return f;
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::set_diagnostics()
{
  const auto sim_field_mgr = get_field_manager("sim");
  // Create all diagnostics
  for (const auto& fname : m_fields_names) {
    if (!sim_field_mgr->has_field(fname)) {
      m_diagnostics[fname] = create_diagnostic(fname);
    }
  }
}

/* ---------------------------------------------------------- */
std::shared_ptr<AtmosphereDiagnostic>
AtmosphereOutput::create_diagnostic (const std::string& diag_field_name)
{
  // We need scream scope resolution, since this->create_diagnostic is hiding it
  auto diag = scream::create_diagnostic(diag_field_name,get_field_manager("sim")->get_grid());

  // Some diags need some extra setup or trigger extra behaviors
  std::string diag_avg_cnt_name = "";
  auto& params = diag->get_params();
  if (diag->name()=="FieldAtPressureLevel") {
    params.set<double>("mask_value",m_fill_value);
    diag_avg_cnt_name = "_"
                      + params.get<std::string>("pressure_value")
                      + params.get<std::string>("pressure_units");
    m_track_avg_cnt = m_track_avg_cnt || m_avg_type!=OutputAvgType::Instant;
  } else if (diag->name()=="FieldAtHeight") {
    if (params.get<std::string>("surface_reference")=="sealevel") {
      diag_avg_cnt_name = "_"
                        + params.get<std::string>("height_value")
                        + params.get<std::string>("height_units") + "_above_sealevel";
      m_track_avg_cnt = m_track_avg_cnt || m_avg_type!=OutputAvgType::Instant;
    }
  }

  // Ensure there's an entry in the map for this diag, so .at(diag_name) always works
  auto& deps = m_diag_depends_on_diags[diag_field_name];

  // Initialize the diagnostic
  const auto sim_field_mgr = get_field_manager("sim");
  for (const auto& freq : diag->get_required_field_requests()) {
    const auto& fname = freq.fid.name();
    if (!sim_field_mgr->has_field(fname)) {
      // This diag depends on another diag. Create and init the dependency
      if (m_diagnostics.count(fname)==0) {
        m_diagnostics[fname] = create_diagnostic(fname);
      }
      deps.push_back(fname);
    }
    diag->set_required_field (get_field(fname,"sim"));
  }

  diag->initialize(util::TimeStamp(),RunType::Initial);

  // If specified, set avg_cnt tracking for this diagnostic.
  if (m_track_avg_cnt) {
    const auto diag_field = diag->get_diagnostic();
    const auto name       = diag_field.name();
    m_field_to_avg_cnt_suffix.emplace(name,diag_avg_cnt_name);
  }

  return diag;
}

// Helper function to mark filled points in a specific layout
void AtmosphereOutput::
update_avg_cnt_view(const Field& field, view_1d_dev& dev_view) {
  const auto& name   = field.name();
  const auto& layout = m_layouts.at(name);
  const auto& dims   = layout.dims();
        auto  data   = dev_view.data();
  const auto fill_value = m_fill_value;

  KT::RangePolicy policy(0,layout.size());
  const auto extents = layout.extents();
  switch (layout.rank()) {
    case 0:
    {
      auto src_view_0d = field.get_view<const Real,Device>();
      auto tgt_view_0d = view_Nd_dev<0>(data);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int) {
        if (src_view_0d()!=fill_value) {
          tgt_view_0d() += 1;
        }
      });
      break;
    }
    case 1:
    {
      // For rank-1 views, we use strided layout, since it helps us
      // handling a few more scenarios
      auto src_view_1d = field.get_strided_view<const Real*,Device>();
      auto tgt_view_1d = view_Nd_dev<1>(data,dims[0]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int i) {
        if (src_view_1d(i)!=fill_value) {
          tgt_view_1d(i) += 1;
        }
      });
      break;
    }
    case 2:
    {
      auto src_view_2d = field.get_view<const Real**,Device>();
      auto tgt_view_2d = view_Nd_dev<2>(data,dims[0],dims[1]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
        int i,j;
        unflatten_idx(idx,extents,i,j);
        if (src_view_2d(i,j)!=fill_value) {
          tgt_view_2d(i,j) += 1;
        }
      });
      break;
    }
    case 3:
    {
      auto src_view_3d = field.get_view<const Real***,Device>();
      auto tgt_view_3d = view_Nd_dev<3>(data,dims[0],dims[1],dims[2]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
        int i,j,k;
        unflatten_idx(idx,extents,i,j,k);
        if (src_view_3d(i,j,k)!=fill_value) {
          tgt_view_3d(i,j,k) += 1;
        }
      });
      break;
    }
    case 4:
    {
      auto src_view_4d = field.get_view<const Real****,Device>();
      auto tgt_view_4d = view_Nd_dev<4>(data,dims[0],dims[1],dims[2],dims[3]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
        int i,j,k,l;
        unflatten_idx(idx,extents,i,j,k,l);
        if (src_view_4d(i,j,k,l)!=fill_value) {
          tgt_view_4d(i,j,k,l) += 1;
        }
      });
      break;
    }
    case 5:
    {
      auto src_view_5d = field.get_view<const Real*****,Device>();
      auto tgt_view_5d = view_Nd_dev<5>(data,dims[0],dims[1],dims[2],dims[3],dims[4]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
        int i,j,k,l,m;
        unflatten_idx(idx,extents,i,j,k,l,m);
        if (src_view_5d(i,j,k,l,m)!=fill_value) {
          tgt_view_5d(i,j,k,l,m) += 1;
        }
      });
      break;
    }
    case 6:
    {
      auto src_view_6d = field.get_view<const Real******,Device>();
      auto tgt_view_6d = view_Nd_dev<6>(data,dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int idx) {
        int i,j,k,l,m,n;
        unflatten_idx(idx,extents,i,j,k,l,m,n);
        if (src_view_6d(i,j,k,l,m,n)!=fill_value) {
          tgt_view_6d(i,j,k,l,m,n) += 1;
        }
      });
      break;
    }
    default:
      EKAT_ERROR_MSG (
            "Error! Field rank not not supported by AtmosphereOutput.\n"
          "  - field name:   " + field.name() + "\n"
          "  - field layout: " + layout.to_string() + "\n");
  }
}

} // namespace scream
