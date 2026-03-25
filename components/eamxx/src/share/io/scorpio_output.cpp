#include "share/io/scorpio_output.hpp"

#include "share/field/field_utils.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/remap/horizontal_remapper.hpp"
#include "share/remap/vertical_remapper.hpp"
#include "share/util/eamxx_timing.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_string_utils.hpp>
#include <ekat_units.hpp>

#include <numeric>

namespace
{
// Helper lambda, to copy extra data (io string attributes plus filled settings).
// This will be used if any remapper is created, to ensure atts set by atm_procs are not lost
void
transfer_extra_data(const scream::Field &src, scream::Field &tgt)
{
  if (src.is_aliasing(tgt))
    return;

  // Transfer io string attributes
  const std::string io_string_atts_key = "io: string attributes";
  using stratts_t                      = std::map<std::string, std::string>;
  const auto &src_atts = src.get_header().get_extra_data<stratts_t>(io_string_atts_key);
  auto &dst_atts       = tgt.get_header().get_extra_data<stratts_t>(io_string_atts_key);
  for (const auto &[name, val] : src_atts) {
    dst_atts[name] = val;
  }
};

// Helper function to get the name of a transposed helper field from a layout and data type
std::string
get_transposed_helper_name(const scream::FieldLayout& layout, const scream::DataType data_type)
{
  return "transposed_" + layout.transpose().to_string() + "_" + e2str(data_type);
}

// Note: this is also declared in eamxx_scorpio_interface.cpp. Move it somewhere else?
template <typename T>
std::string
print_map_keys(const std::map<std::string, T> &map)
{
  std::string s;
  for (const auto &it : map) {
    s += it.first + ",";
  }
  s.pop_back();
  return s;
}
} // anonymous namespace

namespace scream
{

template <typename T>
bool
has_duplicates(const std::vector<T> &c)
{
  std::set<T> s(c.begin(), c.end());
  return c.size() > s.size();
}

AtmosphereOutput::AtmosphereOutput(const ekat::Comm &comm, const std::vector<Field> &fields,
                                   const std::shared_ptr<const grid_type> &grid)
 : m_comm(comm)
{
  // This version of AtmosphereOutput is for quick output of fields (no remaps, no time dim)
  m_avg_type     = OutputAvgType::Instant;
  m_add_time_dim = false;

  // Create a FieldManager with the input fields
  auto fm = std::make_shared<FieldManager>(grid);
  for (auto f : fields) {
    fm->add_field(f);
    m_fields_names.push_back(f.name());
  }

  // No remaps: set all FM except the one for scorpio (created in init())
  m_field_mgrs[FromModel] = m_field_mgrs[AfterVertRemap] = m_field_mgrs[AfterHorizRemap] = fm;

  // Setup I/O structures
  init();

  auto fname = [](const Field& f) { return f.name(); };
  m_stream_name = ekat::join(fields,fname,",");
}

AtmosphereOutput::AtmosphereOutput(const ekat::Comm &comm, const ekat::ParameterList &params,
                                   const std::shared_ptr<const fm_type> &field_mgr,
                                   const std::string &grid_name)
 : m_comm(comm),
   m_add_time_dim(true)
{
  using vos_t = std::vector<std::string>;

  // The param list name will be the name of this stream.
  // This works great for regular CIME cases, where the param list name is
  // the name of the yaml file where the options are read from.
  m_stream_name = params.name();

  // Is this output set to be transposed?
  if (params.isParameter("transpose")) {
    m_transpose = params.get<bool>("transpose");
  }

  auto gm = field_mgr->get_grids_manager();

  // Figure out what kind of averaging is requested
  auto avg_type = params.get<std::string>("averaging_type");
  m_avg_type    = str2avg(avg_type);
  EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: instant, Max, Min, Average. Case insensitive.\n");

  // By default, IO is done directly on the field mgr grid
  auto fm_grid = field_mgr->get_grids_manager()->get_grid(grid_name);

  std::string output_data_layout = "default";
  if (params.isParameter("field_names")) {
    // This simple parameter list option does *not* allow to remap fields
    // to an io grid different from that of the field manager. In order to
    // use that functionality, you need the full syntax
    m_fields_names = params.get<vos_t>("field_names");
    if (params.isParameter("output_data_layout"))
      output_data_layout = params.get<std::string>("output_data_layout");
  } else if (params.isSublist("fields")){
    const auto& f_pl = params.sublist("fields");
    bool grid_found = false;
    for (const auto& grid_name : fm_grid->aliases()) {
      if (f_pl.isSublist(grid_name)) {
        grid_found = true;
        const auto& pl = f_pl.sublist(grid_name);
        if (pl.isParameter("output_data_layout"))
          output_data_layout = pl.get<std::string>("output_data_layout");
        if (pl.isType<vos_t>("field_names")) {
          m_fields_names = pl.get<vos_t>("field_names");
        } else if (pl.isType<std::string>("field_names")) {
          m_fields_names.resize(1, pl.get<std::string>("field_names"));
          if (m_fields_names[0]=="NONE") {
            m_fields_names.clear();
          }
        }
      }
    }
    EKAT_REQUIRE_MSG (grid_found,
        "Error! Bad formatting of output yaml file. Missing 'fields->$grid_name` sublist.\n");
  }

  bool change_data_layout = false;
  if (not fm_grid->get_aux_grid(output_data_layout)) {
    EKAT_REQUIRE_MSG (output_data_layout=="native" or output_data_layout=="default",
        "Error! Grid for requested output_data_layout not found among this grid's aux grids.\n"
        " - grid name: " + fm_grid->name() + "\n"
        " - output_data_layout: " + output_data_layout + "\n");
  } else {
    change_data_layout = true;
  }

  // Check if remapping and if so create the appropriate remapper
  // Note: We currently support three remappers
  //   - vertical remapping from file
  //   - horizontal remapping from file
  //   - online remapping which is setup using the create_remapper function
  const bool use_vertical_remap_from_file = params.isParameter("vertical_remap_file");
  const bool use_horiz_remap_from_file = params.isParameter("horiz_remap_file");
  if (change_data_layout) {
    EKAT_REQUIRE_MSG(!use_vertical_remap_from_file and !use_horiz_remap_from_file,
        "[AtmosphereOutput] Error! Online Dyn->PhysGLL remapping not supported along with vertical and/or horizontal remapping from file");
  }

  auto& fm_model = m_field_mgrs[FromModel];
  auto& fm_after_vr = m_field_mgrs[AfterVertRemap];
  auto& fm_after_hr = m_field_mgrs[AfterHorizRemap];

  // For simplicity, we create a "copy" of the input fm, so we can stuff also diags in it
  fm_model = std::make_shared<FieldManager>(fm_grid,RepoState::Closed);

  // Add ALL field of the FM that are on the output grid
  for (const auto& [name,f_ptr] : field_mgr->get_repo(grid_name)) {
    fm_model->add_field(*f_ptr);
  }

  // Avg count only makes sense if we have
  //  - non-instant output
  //  - we have fields that have masks set. For instance:
  //    - vertically remapped output
  //    - field_at_XhPa diagnostic
  //    - cosp fields
  // NOTE: we are currently only doing this for AVERAGE output, but we should prob extend to MAX/MIN
  // NOTE 2: while we are in the process of transitioning to use Field's mask API, there are still
  //         some diags that use the "mask_data" extra data, so we must check for that too
  // NOTE 3: if in process_requested_fields we realize that NO field is masked, we'll reset it to false
  if (m_avg_type==OutputAvgType::Average) {
    m_track_avg_cnt = true;
    if (params.isParameter("fill_threshold")) {
      m_avg_coeff_threshold = params.get<Real>("fill_threshold");
    }
  }

  // Then we 1) create aliases, and b) create diagnostics, adding alias/diag fields to fm_model
  process_requested_fields ();

  // Setup remappers - if needed
  auto grid_after_vr = fm_grid;
  if (use_vertical_remap_from_file) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    auto vert_remap_file   = params.get<std::string>("vertical_remap_file");
    auto p_mid = fm_model->get_field("p_mid");
    auto p_int = fm_model->get_field("p_int");
    auto vert_remapper = std::make_shared<VerticalRemapper>(fm_model->get_grid(),vert_remap_file);
    vert_remapper->set_source_pressure (p_mid,p_int);
    vert_remapper->set_extrapolation_type(VerticalRemapper::Mask); // both Top AND Bot
    m_vert_remapper = vert_remapper;
    m_vert_remapper->set_name(m_stream_name + " VertRemap");
    if (params.isParameter("enable_fine_grain_timers")) {
      m_vert_remapper->toggle_timers(params.get<bool>("enable_fine_grain_timers"));
    }

    grid_after_vr = m_vert_remapper->get_tgt_grid();
    fm_after_vr = std::make_shared<FieldManager>(grid_after_vr,RepoState::Closed);

    for (const auto& fname : m_fields_names) {
      auto src = fm_model->get_field(fname,fm_grid->name());
      auto tgt = m_vert_remapper->register_field_from_src(src);
      transfer_extra_data (src,tgt);
      fm_after_vr->add_field(tgt);
    }
    m_vert_remapper->registration_ends();
  } else {
    // No vert remap. Simply alias the fm from the model
    fm_after_vr = fm_model;
  }

  // Online remapper and horizontal remapper follow a similar pattern so we check in the same conditional.
  auto grid_after_hr = grid_after_vr;
  if (change_data_layout || use_horiz_remap_from_file) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    if (use_horiz_remap_from_file) {
      // Construct the coarsening remapper
      auto horiz_remap_file   = params.get<std::string>("horiz_remap_file");
      m_horiz_remapper = std::make_shared<HorizontalRemapper>(grid_after_vr,horiz_remap_file,true);
    } else {
      // Construct a generic remapper (likely, Dyn->PhysicsGLL)
      grid_after_hr = fm_grid->get_aux_grid(output_data_layout);
      m_horiz_remapper = gm->create_remapper(grid_after_vr,grid_after_hr);
    }
    m_horiz_remapper->set_name(m_stream_name + " HorizRemap");
    if (params.isParameter("enable_fine_grain_timers")) {
      m_horiz_remapper->toggle_timers(params.get<bool>("enable_fine_grain_timers"));
    }

    grid_after_hr = m_horiz_remapper->get_tgt_grid();
    fm_after_hr = std::make_shared<FieldManager>(grid_after_hr,RepoState::Closed);

    for (const auto& fname : m_fields_names) {
      auto src = fm_after_vr->get_field(fname,grid_after_vr->name());
      auto tgt = m_horiz_remapper->register_field_from_src(src);
      transfer_extra_data (src,tgt);
      fm_after_hr->add_field(tgt);
    }
    m_horiz_remapper->registration_ends();
  } else {
    // No vert remap. Simply alias the fm after vr
    fm_after_hr = fm_after_vr;
  }

  // Setup I/O structures (including the scorpio FM)
  init ();
}

AtmosphereOutput::
~AtmosphereOutput()
{
  // NOTE: yes, we remove ALL the diags in the static var, even if some other
  //       output stream may still use it. But by the time this destructor runs,
  //       we are likely at finalization. This static var is only needed at init,
  //       to prevent two output streams creating the same diag. So when this
  //       destructor runs, it's fine to clean up this static var
  for (auto d : m_diagnostics) {
    const auto& name = d->get_diagnostic().name();
    m_diag_repo.erase(name);
  }
}

void AtmosphereOutput::
set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
  EKAT_REQUIRE_MSG (atm_logger, "Error! Invalid logger pointer.\n");
  m_atm_logger = atm_logger;
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::restart (const std::string& filename)
{
  // Create an input stream on the fly, and init averaging data
  const auto& fm = m_field_mgrs[Scorpio];
  std::vector<Field> fields;
  for (const auto& [name,f_ptr] : fm->get_repo()) {
    fields.push_back(*f_ptr);
  }
  for (const auto& f : m_avg_counts) {
    fields.push_back(f);
  }

  AtmosphereInput hist_restart (filename, fm->get_grid(), fields);
  hist_restart.read_variables();
}

void AtmosphereOutput::init()
{
  auto fm_after_hr = m_field_mgrs[AfterHorizRemap];
  m_io_grid  = fm_after_hr->get_grid();
  m_latlon_output = m_io_grid->has_geometry_data("lat_idx");

  EKAT_REQUIRE_MSG (m_io_grid->is_unique(),
      "Error! I/O only supports grids which are 'unique', meaning that the\n"
      "       map dof_gid->proc_id is well defined.\n");
  EKAT_REQUIRE_MSG (
      (m_io_grid->get_global_max_dof_gid()-m_io_grid->get_global_min_dof_gid()+1)==m_io_grid->get_num_global_dofs(),
      "Error! In order for IO to work, the grid must (globally) have dof gids in interval [gid_0,gid_0+num_global_dofs).\n");

  // Create FM for scorpio. The fields in this FM are guaranteed to NOT have parents/padding
  auto fm_scorpio = m_field_mgrs[Scorpio] = std::make_shared<FieldManager>(fm_after_hr->get_grid(),RepoState::Closed);
  for (size_t i = 0; i < m_fields_names.size(); ++i) {
    const auto& fname = m_fields_names[i];
    const auto& f = fm_after_hr->get_field(fname);
    const auto& fh = f.get_header();
    const auto& fid = fh.get_identifier();

    // Check if the field for scorpio can alias the field after hremap.
    // It can do so only for Instant output, and if the field is NOT a subfield ant NOT padded
    // Also, if we track avg cnt, we MUST add the fill_value extra data, to trigger fill-value logic
    // when calling Field's update methods
    auto can_alias = m_avg_type==OutputAvgType::Instant and
                     fh.get_alloc_properties().get_padding()==0 and
                     fh.get_parent()==nullptr;
    auto f_for_scorpio = can_alias ? f : Field(fid,true);
    transfer_extra_data (f,f_for_scorpio);
    fm_scorpio->add_field(f_for_scorpio);

    // Store the field layout, so that calls to setup_output_file are easier
    const auto& layout = fid.get_layout();
    m_vars_dims[fname] = get_var_dimnames(m_transpose ? layout.transpose() : layout);

    // Initialize a helper_field for each unique layout.  This can be used for operations
    // such as writing transposed output.
    if (m_transpose) {
      const auto helper_layout = layout.transpose();
      const auto data_type = fid.data_type();
      // Note: helper name is based on the ORIGINAL layout (not transposed) and data type, so that
      // when we look up the helper during write, we use the field's original layout and data type
      const std::string helper_name = get_transposed_helper_name(layout, data_type);
      if (m_helper_fields.find(helper_name) == m_helper_fields.end()) {
        // We can add a new helper field for this layout and data type
        // Use Units::invalid() since this helper is reused for fields with same layout but different units
        using namespace ekat::units;
        auto fid_helper = fid.clone(helper_name).reset_units(Units::invalid()).reset_layout(helper_layout);
        Field helper(fid_helper);
        helper.get_header().get_alloc_properties().request_allocation();
        helper.allocate_view();
        m_helper_fields[helper_name] = helper;
      }
    }

    // Now check that all the dims of this field are already set to be registered.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    for (int j=0; j<layout.rank(); ++j) {
      if (tags[j]==FieldTag::Column and m_latlon_output) {
        // We need to make sure we are registering lat and lon as dimensions
        auto lat = m_io_grid->get_geometry_data("lat");
        auto lon = m_io_grid->get_geometry_data("lon");

        m_dims_len.emplace("lat",lat.get_header().get_identifier().get_layout().size());
        m_dims_len.emplace("lon",lon.get_header().get_identifier().get_layout().size());

        continue;
      }
      // check tag against m_dims_len map.  If not in there, then add it.
      std::string dimname = m_io_grid->has_special_tag_name(tags[j])
                          ? m_io_grid->get_special_tag_name(tags[j])
                          : layout.names()[j];

      // If t==CMP, and the name stored in the layout is "dim" (the default) or "bin",
      // we append also the extent, to allow different vector dims in the file
      // TODO: generalize this to all tags, for now hardcoding to dim and bin only
      dimname += (dimname=="dim" or dimname=="bin") ? std::to_string(dims[j]) : "";

      auto is_partitioned = m_io_grid->get_partitioned_dim_tag()==tags[j];
      int dimlen = is_partitioned
                  ? m_io_grid->get_partitioned_dim_global_size()
                  : layout.dim(j);
      auto it_bool = m_dims_len.emplace(dimname,dimlen);
      EKAT_REQUIRE_MSG(it_bool.second or it_bool.first->second==dimlen,
        "Error! Dimension " + dimname + " on field " + fname + " has conflicting lengths.\n"
        "  - old length: " + std::to_string(it_bool.first->second) + "\n"
        "  - new length: " + std::to_string(dimlen) + "\n"
        "If same name applies to different dims (e.g. PhysicsGLL and PhysicsPG2 define "
        "\"ncol\" at different lengths), reset tag name for one of the grids.\n");

      if (is_partitioned) {
        EKAT_REQUIRE_MSG (m_decomp_dimname=="" or m_decomp_dimname==dimname,
            "Error! Decomposed dimension name was already set for a different dimension.\n"
            " - old name: " + m_decomp_dimname + "\n"
            " - new name: " + dimname + "\n");
        m_decomp_dimname = dimname;
      }
    }

    if (m_track_avg_cnt and f.has_mask()) {
      // Create and store a Field to track the averaging count for this layout
      set_avg_cnt_tracking(f);
    }
  }

  // For non-instantaneous output, ensure scorpio fields are
  // inited with correct value for accumulation
  if (m_avg_type!=OutputAvgType::Instant)
    reset_scorpio_fields();
}

void AtmosphereOutput::
init_timestep (const util::TimeStamp& start_of_step)
{
  for (auto diag : m_diagnostics) {
    diag->init_timestep(start_of_step);
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
    m_atm_logger->info("[EAMxx::scorpio_output] Writing variables to file");
    m_atm_logger->info("  file name: " + filename);
  }

  // Update all diagnostics, we need to do this before applying the remapper
  // to make sure that the remapped fields are the most up to date.
  compute_diagnostics(allow_invalid_fields);

  auto apply_remap = [&](AbstractRemapper& remapper)
  {
    remapper.remap_fwd();

    for (int i=0; i<remapper.get_num_fields(); ++i) {
      // Need to update the time stamp of the fields on the IO grid,
      // to avoid throwing an exception later
      auto src = remapper.get_src_field(i);
      auto tgt = remapper.get_tgt_field(i);

      auto src_t = src.get_header().get_tracking().get_time_stamp();
      tgt.get_header().get_tracking().update_time_stamp(src_t);
    }
  }; // end apply_remap

  // If needed, remap fields from their grid to the unique grid, for I/O
  if (m_vert_remapper) {
    start_timer("EAMxx::IO::vert_remap");
    apply_remap(*m_vert_remapper);
    stop_timer("EAMxx::IO::vert_remap");
  }

  if (m_horiz_remapper) {
    start_timer("EAMxx::IO::horiz_remap");
    apply_remap(*m_horiz_remapper);
    stop_timer("EAMxx::IO::horiz_remap");
  }

  auto fm_scorpio = m_field_mgrs[Scorpio];
  auto fm_after_hr = m_field_mgrs[AfterHorizRemap];

  // If tracking avg count, update the count at each field location separately.
  // We do count++ only where the fields are NOT equal to the fill value.
  // Note, we assume that all fields that share a layout are also masked/filled in the same way.
  if (m_track_avg_cnt) {
    // Since 2+ fields may have same avg count, make sure we update the counts only ONCE.
    std::set<std::string> updated;

    for (auto& [fname, count] : m_field_to_avg_count) {
      auto count_updated = updated.count(count.name())>0;
      if (count_updated) {
        // was already in the map (hence, already updated)
        continue;
      }

      auto field = fm_after_hr->get_field(fname);
      auto mask = field.get_mask();

      // mask=1 for "good" entries, and mask=0 otherwise.
      count.update(mask,1,1);

      // Handle writing the average count variables to file
      if (is_write_step) {
        // Bring data to host
        count.sync_to_host();

        auto func_start = std::chrono::steady_clock::now();
        if (m_transpose) {
          const auto& id = count.get_header().get_identifier();
          const auto& layout = id.get_layout();
          const auto data_type = id.data_type();
          const std::string helper_name = get_transposed_helper_name(layout, data_type);
          auto& temp = m_helper_fields.at(helper_name);
          transpose(count,temp);
          temp.sync_to_host();
          scorpio::write_var(filename,count.name(),temp.get_internal_view_data<int,Host>());
        } else {
          scorpio::write_var(filename,count.name(),count.get_internal_view_data<int,Host>());
        }
        auto func_finish = std::chrono::steady_clock::now();
        auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
        duration_write += duration_loc.count();

        // If it's an output step, for Avg we need to ensure count>threshold.
        // Compute the valid/invalid mask fields for this count field
        if (output_step and m_avg_type==OutputAvgType::Average) {
          int min_count = static_cast<int>(std::floor(m_avg_coeff_threshold*nsteps_since_last_output));

          // Recycle mask to find where count<thresh
          auto& valid_count = count.get_mask();
          compute_mask(count,min_count,Comparison::GT,valid_count);
        }
      }
      updated.insert(count.name());
    }
  }

  // Take care of updating and possibly writing fields.
  for (size_t i = 0; i < m_fields_names.size(); ++i) {
    const auto& field_name = m_fields_names[i];
    
    // Get all the info for this field.
    const auto& f_in  = fm_after_hr->get_field(field_name);
          auto& f_out = fm_scorpio->get_field(field_name);

    // Safety check: if we are computing an Average and we are tracking the count,
    // we must have created an avg-count tracking field; otherwise division by the raw
    // number of steps would bias the result wherever fill values occurred.
    if (m_avg_type==OutputAvgType::Average && m_track_avg_cnt && f_in.has_mask()) {
      EKAT_REQUIRE_MSG(m_field_to_avg_count.count(field_name),
        "[AtmosphereOutput::run] Error! No avg-count tracking field for this field.\n"
        " - field name : " + field_name + "\n"
        "Please, contact developers.\n");
    }

    const bool masked = f_in.has_mask();
    switch (m_avg_type) {
      case OutputAvgType::Instant:
        // Note: if f_in aliases f_out, this is a no-op
        f_out.deep_copy(f_in);
        if (masked)
          // We must f=FV wherever the mask is not valid
          f_out.deep_copy(constants::fill_value<Real>,f_in.get_mask(),true);
        break;
      case OutputAvgType::Max:
        if (masked)
          f_out.max(f_in,f_in.get_mask());
        else 
          f_out.max(f_in);
        break;
      case OutputAvgType::Min:
        if (masked)
          f_out.min(f_in,f_in.get_mask());
        else 
          f_out.min(f_in);
        break;
      case OutputAvgType::Average:
        if (masked)
          f_out.update(f_in,1,1,f_in.get_mask());
        else 
          f_out.update(f_in,1,1);
        break;
      default:
        EKAT_ERROR_MSG ("Unexpected/unsupported averaging type.\n");
    }

    // NOTE: we don't divide by the avg cnt when writing rhist files, so only do this for OUTPUT steps
    if (output_step and m_avg_type==OutputAvgType::Average) {
      if (masked) {
        const auto& avg_count = m_field_to_avg_count.at(field_name);
        const auto& valid_count = avg_count.get_mask();

        // Divide by avg_count where large enough, set to fill_value elsewhere
        f_out.scale_inv(avg_count,valid_count);
        f_out.deep_copy(constants::fill_value<Real>,valid_count,true);
      } else {
        // Simply divide by steps count
        f_out.scale(Real(1.0) / nsteps_since_last_output);
      }
    }

    if (is_write_step) {
      // Write to file
      auto func_start = std::chrono::steady_clock::now();
      if (m_transpose) {
        const auto& id = f_out.get_header().get_identifier();
        const auto& layout = id.get_layout();
        const auto data_type = id.data_type();
        const std::string helper_name = get_transposed_helper_name(layout, data_type);
        auto& temp = m_helper_fields.at(helper_name);
        transpose(f_out,temp);
        temp.sync_to_host();
        scorpio::write_var(filename,field_name,temp.get_internal_view_data<Real,Host>());
      } else {
        // Bring data to host (only needed for non-transposed output)
        f_out.sync_to_host();
        scorpio::write_var(filename,field_name,f_out.get_internal_view_data<Real,Host>());
      }
      auto func_finish = std::chrono::steady_clock::now();
      auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
      duration_write += duration_loc.count();
    }
  }

  if (is_write_step) {
    m_atm_logger->info("  Done! Elapsed time: " + std::to_string(duration_write/1000.0) +" seconds");
  }
} // run

long long AtmosphereOutput::
res_dep_memory_footprint () const
{
  long long rdmf = 0;

  // Loop over ALL field mgr, and ALL fields in each of them.
  // Keep track of Field obj we parse, so we don't count them twice
  std::set<Field*> fields;

  for (const auto& [phase,fm] : m_field_mgrs) {
    auto is_diag = [&](const Field& f) {
      const auto& groups = f.get_header().get_tracking().get_groups_names();
      return ekat::contains(groups,"diagnostic");
    };
    for (const auto& [fname,f_ptr] : fm->get_repo()) {
      auto [it, inserted] = fields.insert(f_ptr.get());
      if (not inserted) {
        // We already parsed this field
        continue;
      }
      const auto& fap = f_ptr->get_header().get_alloc_properties();
      if (fap.is_subfield()) {
        // We don't count subfields
        continue;
      }
      if (phase==FromModel and not is_diag(*f_ptr)) {
        // We don't count fields from the model, as we only shallow-copied them
        continue;
      }
      rdmf += fap.get_alloc_size();
    }
  }

  return rdmf;
}

void AtmosphereOutput::set_avg_cnt_tracking(const Field& f)
{
  const auto& mask = f.get_mask();
  const auto avg_cnt_name = "avg_cnt_" + mask.name();

  // Look for an avg count field with the right name
  for (const auto& count : m_avg_counts) {
    if (count.name()==avg_cnt_name) {
      // We already created this avg count field
      m_field_to_avg_count[f.name()] = count;
      return;
    }
  }

  // We have not created this avg count field yet.
  const auto& layout = mask.get_header().get_identifier().get_layout();
  m_vars_dims[avg_cnt_name] = get_var_dimnames(m_transpose ? layout.transpose() : layout);

  // The count is basically a mask with values in [0,max_int)
  auto count = mask.clone(avg_cnt_name);

  // For Average, we also check if count is larger than a threshold, so we need a mask
  // for count to establish if the count is "valid"
  count.create_mask();

  m_avg_counts.push_back(count);
  m_field_to_avg_count[f.name()] = count;
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
reset_scorpio_fields()
{
  // Reset the fields for scorpio to whatever is the proper accumulation value (if avg!=Instant)
  Real value;
  switch (m_avg_type) {
    case OutputAvgType::Max:
      value = -std::numeric_limits<Real>::infinity(); break;
    case OutputAvgType::Min:
      value =  std::numeric_limits<Real>::infinity(); break;
    case OutputAvgType::Average:
      value =  0.0;                                   break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unexpected averaging type.\n");
  }

  auto fm = m_field_mgrs[Scorpio];
  for (const auto& name : m_fields_names) {
    fm->get_field(name).deep_copy(value);
  }
  for (auto& count : m_avg_counts) {
    count.deep_copy(0);
  }
}
/* ---------------------------------------------------------- */
void AtmosphereOutput::
register_variables(const std::string& filename,
                   const std::string& fp_precision,
                   const scorpio::FileMode mode)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (ekat::contains(strvec_t{"float","single","double","real"},fp_precision),
      "Error! Invalid/unsupported value for fp_precision.\n"
      "  - input value: " + fp_precision + "\n"
      "  - supported values: float, single, double, real\n");

  // Cycle through all fields and register variables
  for (const auto& field_name : m_fields_names) {
    const auto& f = m_field_mgrs[Scorpio]->get_field(field_name);
    const auto& fid  = f.get_header().get_identifier();
    const auto& dimnames = m_vars_dims.at(field_name);
    std::string units = fid.get_units().to_string();

    // TODO  Need to change dtype to allow for other variables.
    // Currently the field_manager only stores Real variables so it is not an issue,
    // but in the future if non-Real variables are added we will want to accomodate that.

    if (mode==scorpio::FileMode::Append) {
      // Simply check that the var is in the file, and has the right properties
      EKAT_REQUIRE_MSG (scorpio::has_var(filename,field_name),
          "Error! Cannot append, due to variable missing from the file.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + field_name + "\n");
      const auto& var = scorpio::get_var(filename,field_name);
      EKAT_REQUIRE_MSG (var.dim_names()==dimnames,
          "Error! Cannot append, due to variable dimensions mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + field_name + "\n"
          "  - var dims : " + ekat::join(dimnames,",") + "\n"
          "  - var dims from file: " + ekat::join(var.dim_names(),",") + "\n");
      EKAT_REQUIRE_MSG (var.units==units,
          "Error! Cannot append, due to variable units mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + field_name + "\n"
          "  - var units: " + units + "\n"
          "  - var units from file: " + var.units + "\n");
      EKAT_REQUIRE_MSG (var.time_dep==m_add_time_dim,
          "Error! Cannot append, due to time dependency mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + field_name + "\n"
          "  - var time dep: " + (m_add_time_dim ? "yes" : "no") + "\n"
          "  - var time dep from file: " + (var.time_dep ? "yes" : "no") + "\n");
    } else {
      scorpio::define_var (filename, field_name, units, dimnames,
                            "real",fp_precision, m_add_time_dim);

      // Add FillValue as an attribute of each variable
      // FillValue is a protected metadata, do not add it if it already existed
      if (fp_precision=="double" or
          (fp_precision=="real" and std::is_same<Real,double>::value)) {
        scorpio::set_attribute(filename, field_name, "_FillValue",constants::fill_value<double>);
      } else {
        scorpio::set_attribute(filename, field_name, "_FillValue",constants::fill_value<float>);
      }
      if (m_alias_to_orig.count(field_name)==1) {
        // Store what this field is the alias of
        scorpio::set_attribute(filename, field_name, "alias_of",m_alias_to_orig[field_name]);
      }

      // If this is has subfields, add list of its children
      const auto& children = f.get_header().get_children();
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
        scorpio::set_attribute(filename,field_name,"sub_fields",children_list);
      }

      // If tracking average count variables then add the name of the tracking variable for this variable
      if (m_field_to_avg_count.count(field_name)) {
        const auto& count = m_field_to_avg_count.at(field_name);
        scorpio::set_attribute(filename,field_name,"averaging_count_tracker",count.name());
      }

      // Atm procs may have set some request for metadata.
      using stratts_t = std::map<std::string,std::string>;
      const auto& str_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
      for (const auto& [att_name,att_val] : str_atts) {
        scorpio::set_attribute(filename,field_name,att_name,att_val);
      }

      // Gather longname (if not already in the io: string attributes)
      if (str_atts.count("long_name")==0) {
        auto longname = m_default_metadata.get_longname(field_name);
        scorpio::set_attribute(filename, field_name, "long_name", longname);
      }

      // Gather standard name, CF-Compliant (if not already in the io: string attributes)
      if (str_atts.count("standard_name")==0) {
        auto standardname = m_default_metadata.get_standardname(field_name);
        scorpio::set_attribute(filename, field_name, "standard_name", standardname);
      }
      
      // If output represents an statistic over a time range add a "cell methods"
      // attribute.
      switch (m_avg_type) {
        case OutputAvgType::Instant:
          scorpio::set_attribute(filename, field_name, "cell_methods", "time: point");
          break;  // Don't add the attribute
        case OutputAvgType::Max:
          scorpio::set_attribute(filename, field_name, "cell_methods", "time: maximum");
          break;
        case OutputAvgType::Min:
          scorpio::set_attribute(filename, field_name, "cell_methods", "time: minimum");
          break;
        case OutputAvgType::Average:
          scorpio::set_attribute(filename, field_name, "cell_methods", "time: mean");
          break;
        default:
          EKAT_ERROR_MSG ("Unexpected/unsupported averaging type.\n");
      }

      // If output contains the column dimension add a "coordinates" attribute.
      if (fid.get_layout().has_tag(COL)) {
        scorpio::set_attribute(filename, field_name, "coordinates", "lat lon");
      }

      // If this is transposed output, mark the variable
      if (m_transpose) {
        scorpio::set_attribute(filename, field_name, "transposed_output", "true");
      }
    }
  }

  // Now register the average count variables (if any)
  for (const auto& f : m_avg_counts) {
    const auto& name = f.name();
    const auto& dimnames = m_vars_dims.at(name);
    if (mode==scorpio::FileMode::Append) {
      // Similar to the regular fields above, check that the var is in the file, and has the right properties
      EKAT_REQUIRE_MSG (scorpio::has_var(filename,name),
          "Error! Cannot append, due to variable missing from the file.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n");
      const auto& var = scorpio::get_var(filename,name);
      EKAT_REQUIRE_MSG (var.dim_names()==dimnames,
          "Error! Cannot append, due to variable dimensions mismatch.\n"
          "  - filename : " + filename + "\n"
          "  - varname  : " + name + "\n"
          "  - var dims : " + ekat::join(dimnames,",") + "\n"
          "  - var dims from file: " + ekat::join(var.dim_names(),",") + "\n");
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
      scorpio::define_var(filename, name, dimnames, "int", m_add_time_dim);
    }
  }
} // register_variables

void AtmosphereOutput::set_decompositions(const std::string& filename)
{
  if (m_latlon_output) {
    // We need to find out which (lat,lon) offsets we own
    auto lat_idx_h = m_io_grid->get_geometry_data("lat_idx").get_view<const int*,Host>();
    auto lon_idx_h = m_io_grid->get_geometry_data("lon_idx").get_view<const int*,Host>();
    int ncols = m_io_grid->get_num_local_dofs();
    int nlon = m_io_grid->get_geometry_data("lon").get_header().get_identifier().get_layout().size();
    std::vector<scorpio::offset_t> offsets(ncols);
    for (int i=0; i<ncols; ++i) {
      offsets[i] = lat_idx_h(i)*nlon + lon_idx_h(i);
    }
    scorpio::set_dims_decomp(filename,{"lat","lon"},offsets);
  } else {
    if (m_decomp_dimname=="")
      return;

    // Set the decomposition for the partitioned dimension
    const int local_dim = m_io_grid->get_partitioned_dim_local_size();
    auto gids_f = m_io_grid->get_partitioned_dim_gids();
    auto gids_h = gids_f.get_view<const AbstractGrid::gid_type*,Host>();
    auto min_gid = m_io_grid->get_global_min_partitioned_dim_gid();
    std::vector<scorpio::offset_t> offsets(local_dim);
    for (int idof=0; idof<local_dim; ++idof) {
      offsets[idof] = gids_h[idof] - min_gid;
    }
    scorpio::set_dim_decomp(filename,m_decomp_dimname,offsets);
  }
}

void AtmosphereOutput::
setup_output_file(const std::string& filename,
                  const std::string& fp_precision,
                  const scorpio::FileMode mode)
{
  // Register dimensions with netCDF file.
  for (const auto& [dimname,dimlen] : m_dims_len) {
    if (mode==scorpio::FileMode::Append) {
      // Simply check that the dim is in the file, and has the right extent
      EKAT_REQUIRE_MSG (scorpio::has_dim(filename,dimname),
          "Error! Cannot append, due to missing dim in the file.\n"
          "  - filename: " + filename + "\n"
          "  - dimname : " + dimname + "\n");
      EKAT_REQUIRE_MSG (scorpio::get_dimlen(filename,dimname)==dimlen,
          "Error! Cannot append, due to mismatch dim length.\n"
          "  - filename: " + filename + "\n"
          "  - dimname : " + dimname + "\n"
          "  - old len : " + std::to_string(scorpio::get_dimlen(filename,dimname)) + "\n"
          "  - new len : " + std::to_string(dimlen) + "\n");
    } else {
      scorpio::define_dim(filename,dimname,dimlen);
    }
  }

  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables(filename,fp_precision,mode);

  // Set the offsets of the local dofs in the global vector.
  set_decompositions(filename);
}

void AtmosphereOutput::
compute_diagnostics(const bool allow_invalid_fields)
{
  for (auto diag : m_diagnostics) {
    // Check if all inputs are valid
    bool computable = true;
    for (const auto& f : diag->get_fields_in()) {
      computable &= f.get_header().get_tracking().get_time_stamp().is_valid();

      EKAT_REQUIRE_MSG (computable or allow_invalid_fields,
        "Error! Cannot compute a diagnostic. One dependency has an invalid timestamp.\n"
        " - stream name: " + m_stream_name + "\n"
        " - diag name: " + diag->get_diagnostic().name() + "\n"
        " - dep  name: " + f.name() + "\n");
    }

    auto d = diag->get_diagnostic();
    if (computable) {
      diag->compute_diagnostic();
    }

    bool computed = d.get_header().get_tracking().get_time_stamp().is_valid();

    EKAT_REQUIRE_MSG (computed or allow_invalid_fields,
      "Error! Failed to compute diagnostic.\n"
      " - diag name: " + diag->get_diagnostic().name() + "\n");

    if (not computed) {
      // The diag was either not computable or it may have failed to compute
      // (e.g., t=0 output with a flux-like diag).
      // If we're allowing invalid fields, then we should simply set diag=fill_value
      d.deep_copy(constants::fill_value<float>);
    }
  }
}

void AtmosphereOutput::
process_requested_fields()
{
  // So far, all fields (on the output grid) that ARE in the model FM have been added
  // to the FM stored in this class for the FromModel phase. Anything missing
  // must be either a diagnostic or an alias.

  auto fm_model = m_field_mgrs[FromModel];
  auto fm_grid = m_field_mgrs[FromModel]->get_grid();

  // First, find out which field names are just aliases
  for (auto& name : m_fields_names) {
    auto tokens = ekat::split(name,":=");
    EKAT_REQUIRE_MSG(tokens.size()==2 or tokens.size()==1,
        "Error! Invalid alias request. Should be 'alias:=original'.\n"
        " - request: " + name + "\n");
    if (tokens.size()==2) {
      EKAT_REQUIRE_MSG (m_alias_to_orig.count(tokens[0])==0,
          "Error! The same alias has been used multiple times.\n"
          " - stream name: " + m_stream_name + "\n"
          " - first alias: " + tokens[0] + ":=" + m_alias_to_orig[tokens[0]] + "\n"
          " - second alias: " + tokens[0] + ":=" + tokens[1] + "\n");
      m_alias_to_orig[tokens[0]] = tokens[1];
      name = tokens[0];
    }
  }

  // In case someone has an alias of an alias, we need to resolve the TRUE orig names.
  bool has_multiple_aliasing_layers = false;
  do {
    for (auto it : m_alias_to_orig) {
      if (m_alias_to_orig.count(it.second)>0) {
        it.second = m_alias_to_orig[it.second];
        has_multiple_aliasing_layers = true;
      }
    }
  } while (has_multiple_aliasing_layers);

  EKAT_REQUIRE_MSG (not has_duplicates(m_fields_names),
      "Error! The list of requested output fields contains duplicates.\n"
      " - stream name:  " + m_stream_name + "\n"
      " - fields names: " + ekat::join(m_fields_names,",") + "\n");

  // Helper lambda that initializes a diagnostic
  auto init_diag = [&](const std::shared_ptr<AtmosphereDiagnostic>& diag) {
    // Set inputs in the diag
    for (const auto& freq : diag->get_field_requests()) {
      const auto& dep_name = freq.fid.name();

      auto dep = fm_model->get_field(dep_name);
      diag->set_required_field(dep);
    }

    // Initialize the diag
    diag->initialize(util::TimeStamp(),RunType::Initial);
  };

  // Now process each requested field, if possible. We can process a field if either:
  //  - it is already in the model FM
  //  - it is an alias of a field added to the FM
  //  - it is a diag that ONLY depends on fields already added to the FM
  // We keep scanning the fields names and process a field, removing it from the list.
  // If a field cannot be processed yet (e.g., it's a diag whose deps have not been
  // parsed yet, or an alias of a field not yet parsed), then we simply add it back
  // to the end of the list, which is effectively a queue.
  // To avoid infinite loops, we need to ensure that we either add or remove items
  // to the list. A check on the size is not ok (if we have an alias of an unprocessed
  // field, we are simply moving it back to the end of the line). Instead, we keep
  // track of how many iters it's been since the last fully processed field. This number
  // should NEVER exceed the legnth of the remaining fields
  // There MUST be at least ONE field we can process at every iteration, but
  // some field MAY have to wait until another is processed. E.g., if we have
  // 'foo_at_900hPa' and 'foo:=horiz_winds', the latter (an alias) must be
  // processed first, as it adds the 'foo' alias field to the fm, which the
  // diag 'foo_at_900hPa' needs to be correctly initialized.
  // Notice that, thanks to how we create and add diags to m_diagnostics,
  // the order in which diags appear in m_diagnostics follows the evaluation
  // order, meaning that if diag A depends on diag B, then B appears *before* A.
  // This ensures we can evaluate diags in order at runtime
  std::list<std::string> remaining(m_fields_names.begin(),m_fields_names.end());
  for (const auto& it : m_alias_to_orig) {
    remaining.push_back(it.second);
  }

  bool any_masked_field = false;
  size_t iters_since_last_done = 0;
  for (auto it=remaining.begin(); it!=remaining.end(); ) {
    const auto name = *it;

    EKAT_REQUIRE_MSG (iters_since_last_done<remaining.size(),
        "Error! It seems we're stuck in an infinite loop...\n"
        " Field '" + name + "' seem to cause circular deps.");
    ++iters_since_last_done;

    if (fm_model->has_field(name)) {
      // This is a regular field, not a diagnostic nor an alias.
      iters_since_last_done = 0;
      any_masked_field |= fm_model->get_field(name).has_mask();
    } else if (m_alias_to_orig.count(name)==1) {
      // An alias. If the aliased field was already processed, we can
      // process the alias as well
      if (fm_model->has_field(m_alias_to_orig[name])) {
        const auto& orig = fm_model->get_field(m_alias_to_orig[name]);
        auto alias = orig.alias(name);
        fm_model->add_field(alias);
        iters_since_last_done = 0;

        // TODO: move all remaining diags to use mask field API in Field
        any_masked_field |= alias.has_mask() ||
                            alias.get_header().has_extra_data("mask_data");
      } else {
        // Put this field back at the end of the list, while we wait to process the aliased field
        remaining.push_back(*it);
      }
    } else {
      auto& diag = m_diag_repo[name];
      if (not diag) {
        // First time we run into this diag. Create it
        diag = create_diagnostic(name,fm_model->get_grid());
      }
      // Add its deps to the list of fields to process (if not already in fm_model)
      bool deps_met = true;
      for (const auto& req : diag->get_field_requests()) {
        if (not fm_model->has_field(req.fid.name())) {
          deps_met = false;
          remaining.push_back(req.fid.name());
        }
      }

      // If we are missing any dep, we DELAY adding this diag to m_diagnostics, so that
      // the order in which diags appear is compatible with the evaluation order
      if (deps_met) {
        // Check if already inited (perhaps by another stream)
        if (not diag->is_initialized()) {
          init_diag(diag);
        }
        diag->get_diagnostic().get_header().get_tracking().add_group("diagnostic");
        fm_model->add_field(diag->get_diagnostic());
        m_diagnostics.push_back(diag);
        auto diag_field = diag->get_diagnostic();

        iters_since_last_done = 0;

        // TODO: move all remaining diags to use mask field API in Field
        any_masked_field |= diag_field.has_mask() ||
                            diag_field.get_header().has_extra_data("mask_data");
      } else {
        // we'll come back to this after its deps, so put it back at the end
        remaining.push_back(name);
      }
    }

    // Whether we fully processed this field, or we put it back at the end, we need
    // to remove the current entry, and update the iterator
    it = remaining.erase(it);
  }

  m_track_avg_cnt &= any_masked_field;
}

std::vector<std::string> AtmosphereOutput::
get_var_dimnames (const FieldLayout& layout) const
{
  using namespace ShortFieldTagsNames;
  strvec_t dims;
  for (int i=0; i<layout.rank(); ++i) {
    const auto t = layout.tag(i);
    if (t==COL and m_latlon_output) {
      // Lat-Lon remapping uses a PointGrid target grid, so we replace the
      // single column dimension with separate latitude and longitude dimensions.
      dims.push_back("lat");
      dims.push_back("lon");
    } else {
      auto tag_name = m_io_grid->has_special_tag_name(t)
                    ? m_io_grid->get_special_tag_name(t)
                    : layout.names()[i];
      if (tag_name=="dim" or tag_name=="bin") {
        tag_name += std::to_string(layout.dim(i));
      }
      dims.push_back(tag_name); // Add dimensions string to vector of dims.
    }
  }
  return dims;
}

// Instantiate the static member var
AtmosphereOutput::strmap_t<AtmosphereOutput::diag_ptr_type>
AtmosphereOutput::m_diag_repo;

} // namespace scream
