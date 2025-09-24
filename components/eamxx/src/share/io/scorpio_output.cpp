#include "share/io/scorpio_output.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/scorpio_input.hpp"
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

  // Transfer io string attributes
  const std::string io_string_atts_key = "io: string attributes";
  using stratts_t                      = std::map<std::string, std::string>;
  const auto &src_atts = src.get_header().get_extra_data<stratts_t>(io_string_atts_key);
  auto &dst_atts       = tgt.get_header().get_extra_data<stratts_t>(io_string_atts_key);
  for (const auto &[name, val] : src_atts) {
    dst_atts[name] = val;
  }

  // Transfer whether or not this field MAY contain fill_value (to trigger usage of the
  // proper implementation of Field update methods
  if (src.get_header().may_be_filled()) {
    tgt.get_header().set_may_be_filled(true);
  }
};

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

  // Then we 1) create aliases, and b) create diagnostics, adding alias/diag fields to fm_model
  process_requested_fields ();

  // Avg count only makes sense if we have
  //  - non-instant output
  //  - we have one between:
  //    - vertically remapped output
  //    - field_at_XhPa diagnostic
  //    - fields that can contain invalid values (not yet supported, but RAD may want this at some point)
  // We already set m_track_avg_cnt to true if field_at_XhPa is found in process_requested_fields.
  // Hence, here we only check if vert remap is active

  if (m_avg_type!=OutputAvgType::Instant) {
    if (params.isParameter("track_avg_cnt")) {
      // This is to be used for unit testing only, so that we can test avg cnt even
      // if there is no vert remap and no field_at_XhPa diagnostic in the stream
      m_track_avg_cnt = params.get<bool>("track_avg_cnt");
    }
    if (use_vertical_remap_from_file) {
      m_track_avg_cnt = true;
    }
    if (params.isParameter("fill_threshold")) {
      m_avg_coeff_threshold = params.get<Real>("fill_threshold");
    }
  }

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
      m_horiz_remapper = std::make_shared<CoarseningRemapper>(grid_after_vr,horiz_remap_file,true);
    } else {
      // Construct a generic remapper (likely, Dyn->PhysicsGLL)
      grid_after_hr = fm_grid->get_aux_grid(output_data_layout);
      m_horiz_remapper = gm->create_remapper(grid_after_vr,grid_after_hr);
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
    if (m_avg_type!=OutputAvgType::Instant or
        fh.get_alloc_properties().get_padding()>0 or
        fh.get_parent()!=nullptr) {
      Field copy(fid);
      copy.allocate_view();
      transfer_extra_data (f,copy);
      fm_scorpio->add_field(copy);
    } else {
      fm_scorpio->add_field(f);
    }

    // Store the field layout, so that calls to setup_output_file are easier
    const auto& layout = fid.get_layout();
    m_vars_dims[fname] = get_var_dimnames(layout);

    // Now check that all the dims of this field are already set to be registered.
    const auto& tags = layout.tags();
    const auto& dims = layout.dims();
    for (int j=0; j<layout.rank(); ++j) {
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

    if (m_track_avg_cnt) {
      // Create and store a Field to track the averaging count for this layout
      set_avg_cnt_tracking(fname,layout);
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
    for (auto& [fname, count] : m_field_to_avg_count) {
      count.get_header().set_extra_data("updated",false);
    }

    for (auto& [fname, count] : m_field_to_avg_count) {
      if (count.get_header().get_extra_data<bool>("updated")) {
        continue;
      }

      //////////////////////// TODO ////////////////////////
      // 1. Make the diags/procs COMPUTE the mask for the field, and
      //    set the mask as extra data. We DON'T want to compute it here
      // 2. Create count with same layout as the mask
      // 3. In avg=tally/count, count MAY have smaller layout, hence
      //    be incompatible. E.g., mask came from FieldAtPressureLevel,
      //    and was (ncol), but the field is (ncol,2). In this case, we
      //    ASSUME same mask for all slices, so do scaling on all slices:
      //      avg.component(i).scale_inv(count)
      // 4. In FieldAtPressureLevel, make ALL instances at same Plev share
      //    the same mask field. How? Two ideas:
      //      - store a static map string->Field in class, and use a string
      //        key that encodes grid and press level. Only compute mask if
      //        timestamp of mask field is "old"
      //      - make another diag that computes the mask, make F@plev depend
      //        on that diag.
      // 5. FieldAtPressureLevel (and vremap) should only fill the mask field
      //    if we later need it. E.g, if no AvgCount AND no hremap, we don't need it.
      //////////////////////////////////////////////////////
      auto field = fm_after_hr->get_field(fname);
      auto mask  = count.get_header().get_extra_data<Field>("mask");

      // Find where the field is NOT equal to fill_value
      compute_mask<Comparison::NE>(field,constants::fill_value<Real>,mask);

      // mask=1 for "good" entries, and mask=0 otherwise.
      count.update(mask,1,1);

      // Handle writing the average count variables to file
      if (is_write_step) {
        // Bring data to host
        count.sync_to_host();

        auto func_start = std::chrono::steady_clock::now();
        scorpio::write_var(filename,count.name(),count.get_internal_view_data<int,Host>());
        auto func_finish = std::chrono::steady_clock::now();
        auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
        duration_write += duration_loc.count();

        // If it's an output step, for Avg we need to ensure count>threshold.
        // If count<=threshold, we set count=fill_value, so that fill_val propagates
        // to the output fields when we divide by count later
        if (output_step and m_avg_type==OutputAvgType::Average) {
          int min_count = static_cast<int>(std::floor(m_avg_coeff_threshold*nsteps_since_last_output));

          // Recycle mask to find where count<thresh
          compute_mask<Comparison::LE>(count,min_count,mask);

          // Later, we divide fields by count. By setting count=1 where count<thresholt,
          // we can later do
          //   f.scale_inv(count); // Requires count!=0 anywhere
          //   f.deep_copy(fill_val,mask)
          count.deep_copy(1,mask);
        }
      }
      count.get_header().set_extra_data("updated",true);
    }
  }

  // Take care of updating and possibly writing fields.
  for (size_t i = 0; i < m_fields_names.size(); ++i) {
    const auto& field_name = m_fields_names[i];
    
    // Get all the info for this field.
    const auto& f_in  = fm_after_hr->get_field(field_name);
          auto& f_out = fm_scorpio->get_field(field_name);

    // Safety check: if a field may contain fill values and we are computing an Average,
    // we must have created an avg-count tracking field; otherwise division by the raw
    // number of steps would bias the result wherever fill values occurred.
    if (m_avg_type==OutputAvgType::Average && f_in.get_header().may_be_filled()) {
      EKAT_REQUIRE_MSG(m_field_to_avg_count.count(field_name),
        "[AtmosphereOutput::run] Error! Averaging a fill-aware field without avg-count tracking.\n"
        " - field name : " + field_name + "\n"
        "This indicates the field was marked may_be_filled after output initialization or tracking logic missed it." );
    }

    switch (m_avg_type) {
      case OutputAvgType::Instant:
        f_out.deep_copy(f_in);  break; // Note: if f_in aliases f_out, this is a no-op
      case OutputAvgType::Max:
        f_out.max(f_in);        break;
      case OutputAvgType::Min:
        f_out.min(f_in);        break;
      case OutputAvgType::Average:
        f_out.update(f_in,1,1); break;
      default:
        EKAT_ERROR_MSG ("Unexpected/unsupported averaging type.\n");
    }

    if (is_write_step) {
      // NOTE: we don't divide by the avg cnt for checkpoint output
      if (output_step and m_avg_type==OutputAvgType::Average) {
        // Even if m_track_avg_cnt=true, this field may not need it
        if (m_track_avg_cnt) {
          auto avg_count = m_field_to_avg_count.at(field_name);

          f_out.scale_inv(avg_count);

          const auto& mask = avg_count.get_header().get_extra_data<Field>("mask");
          f_out.deep_copy(constants::fill_value<Real>,mask);
        } else {
          // Divide by steps count only when the summation is complete
          f_out.scale(Real(1.0) / nsteps_since_last_output);
        }
      }

      // Bring data to host
      f_out.sync_to_host();

      // Write to file
      auto func_start = std::chrono::steady_clock::now();
      scorpio::write_var(filename,field_name,f_out.get_internal_view_data<Real,Host>());
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

void AtmosphereOutput::set_avg_cnt_tracking(const std::string& name, const FieldLayout& layout)
{
  // Now create a Field to track the averaging count for this layout
  const auto& avg_cnt_suffix = m_field_to_avg_cnt_suffix[name];
  const auto tags = layout.tags();
  std::string avg_cnt_name = "avg_count" + avg_cnt_suffix;
  for (int i=0; i<layout.rank(); ++i) {
    const auto t = layout.tag(i);
    std::string tag_name = m_io_grid->has_special_tag_name(t)
                         ? m_io_grid->get_special_tag_name(t)
                         : layout.names()[i];

    // If t==CMP, and the name stored in the layout is "dim" (the default) or "bin",
    // we append also the extent, to allow different vector dims in the file
    // TODO: generalize this to all tags, for now hardcoding to dim and bin only
    tag_name += (tag_name=="dim" or tag_name=="bin") ? std::to_string(layout.dim(i)) : "";

    avg_cnt_name += "_" + tag_name;
  }

  // Look for an avg count field with the right name
  for (const auto& f : m_avg_counts) {
    if (f.name()==avg_cnt_name) {
      // We already created this avg count field
      m_field_to_avg_count[name] = f;
      return;
    }
  }

  // We have not created this avg count field yet.
  m_vars_dims[avg_cnt_name] = get_var_dimnames(layout);

  auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier count_id (avg_cnt_name,layout,nondim,m_io_grid->name(),DataType::IntType);
  Field count(count_id);
  count.allocate_view();

  // We will use a helper field for updating cnt, so store it inside the field header
  auto mask = count.clone(count.name()+"_mask");
  count.get_header().set_extra_data("mask",mask);

  m_avg_counts.push_back(count);
  m_field_to_avg_count[name] = count;
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

  // Helper lambda to check if this fm_model field should trigger avg count
  auto check_for_avg_cnt = [&](const Field& f) {
    // We need avg-count tracking for any averaged (non-instant) field that:
    //  - supplies explicit mask info (mask_data or mask_field), OR
    //  - is marked as potentially containing fill values (may_be_filled()).
    // Without this, fill-aware updates skip fill_value during accumulation (good)
    // but we would still divide by the raw nsteps, biasing the result low.
    if (m_avg_type!=OutputAvgType::Instant) {
      const bool has_mask = f.get_header().has_extra_data("mask_data") || f.get_header().has_extra_data("mask_field");
      const bool may_be_filled = f.get_header().may_be_filled();
      if (has_mask || may_be_filled) {
        m_track_avg_cnt = true;
        // Avoid duplicate insertion if already present (e.g., mask + filled both true)
        if (m_field_to_avg_cnt_suffix.count(f.name())==0) {
          m_field_to_avg_cnt_suffix.emplace(f.name(), "_" + f.name());
        }
      }
    }
  };

  // Helper lambda that initializes a diagnostic
  auto init_diag = [&](const std::shared_ptr<AtmosphereDiagnostic>& diag) {
    // Set inputs in the diag
    for (const auto& freq : diag->get_required_field_requests()) {
      const auto& dep_name = freq.fid.name();

      auto dep = fm_model->get_field(dep_name);
      diag->set_required_field(dep);
    }

    // Initialize the diag
    diag->initialize(util::TimeStamp(),RunType::Initial);
  };

  auto check_diag_avg_cnt = [&](const std::shared_ptr<AtmosphereDiagnostic>& diag) {
    // Set the diag field in the FM
    auto diag_field = diag->get_diagnostic();

    // Add the field to the diag group
    diag_field.get_header().get_tracking().add_group("diagnostic");

    // Some diags need some extra setup or trigger extra behaviors
    std::string diag_avg_cnt_name = "";
    auto& params = diag->get_params();
    if (diag->name()=="FieldAtPressureLevel") {
      diag_avg_cnt_name = "_"
                        + params.get<std::string>("pressure_value")
                        + params.get<std::string>("pressure_units");
      m_track_avg_cnt |= m_avg_type!=OutputAvgType::Instant;
    } else if (diag->name()=="FieldAtHeight") {
      if (params.get<std::string>("surface_reference")=="sealevel") {
        diag_avg_cnt_name = "_"
                          + params.get<std::string>("height_value")
                          + params.get<std::string>("height_units") + "_above_sealevel";
        m_track_avg_cnt |= m_avg_type!=OutputAvgType::Instant;
      }
    } else if (diag->name()=="AerosolOpticalDepth550nm") {
      m_track_avg_cnt = m_track_avg_cnt || m_avg_type!=OutputAvgType::Instant;
      diag_avg_cnt_name = "_" + diag->name();
    }
    else if (diag_field.get_header().has_extra_data("mask_data")) {
      m_track_avg_cnt = m_track_avg_cnt || m_avg_type!=OutputAvgType::Instant;
      diag_avg_cnt_name = "_" + diag_field.name();
    }

    // If specified, set avg_cnt tracking for this diagnostic.
    if (m_track_avg_cnt) {
      m_field_to_avg_cnt_suffix.emplace(diag_field.name(),diag_avg_cnt_name);
    }
  };
  
  // Now process each requested field, if possible. We can process a field if either:
  //  - it is already in the model FM
  //  - it is an alias of a field added to the FM
  //  - it is a diag that ONLY depends on fields already added to the FM
  // We keep scanning the fields names and process a field, removing it from the list.
  // There MUST be at least ONE field we can process at every iteration, but
  // some field MAY have to wait until another is processed. E.g., if we have
  // 'foo_at_900hPa' and 'foo:=horiz_winds', the latter (an alias) must be
  // processed first, as it adds the 'foo' alias field to the fm, which the
  // diag 'foo_at_900hPa' needs to be correctly initialized.
  // Notice that, thanks to how we create and add diags to m_diagnostics,
  // the order in which diags appear in m_diagnostics follows the evaluation
  // order, meaning that if diag A depends on diag B, then B appears *before* A.
  // This ensures we can evaluate diags in order at runtime
  bool done = false;
  std::set<std::string> remaining(m_fields_names.begin(),m_fields_names.end());
  for (const auto& it : m_alias_to_orig) {
    remaining.insert(it.second);
  }
  while (not done) {
    // We can't add-to/rm-form a std:;set while iterating on it, as that could
    // change the end iterator. Hence, keep track of what we add or remove,
    // and add/remove after the for loop ends
    std::set<std::string> remove_these, add_these;
    for (const auto& name : remaining) {
      if (fm_model->has_field(name)) {
        // This is a regular field, not a diagnostic nor an alias.
        // Still, we might need to do some extra setup, like for avg_count.
        const auto& f = m_field_mgrs[FromModel]->get_field(name);
        check_for_avg_cnt(f);
        remove_these.insert(name);
      } else if (m_alias_to_orig.count(name)==1) {
        // An alias. If the aliased field was already processed, we can
        // process the alias as well
        if (fm_model->has_field(m_alias_to_orig[name])) {
          const auto& orig = fm_model->get_field(m_alias_to_orig[name]);
          auto alias = orig.alias(name);
          fm_model->add_field(alias);
          remove_these.insert(name);
        }
      } else {
        auto& diag = m_diag_repo[name];
        if (not diag) {
          // First time we run into this diag. Create it
          diag = create_diagnostic(name,fm_model->get_grid());
        }
        // Add its deps to the list of fields to process (if not already in fm_model)
        bool deps_met = true;
        for (const auto& req : diag->get_required_field_requests()) {
          if (not fm_model->has_field(req.fid.name())) {
            deps_met = false;
            add_these.insert(req.fid.name());
          }
        }

        // If we are missing any dep, we DELAY adding this diag to m_diagnostics, so that
        // the order in which diags appear is compatible with the evaluation order
        if (deps_met) {
          // Check if already inited (perhaps by another stream)
          if (not diag->is_initialized()) {
            init_diag(diag);
          }
          check_diag_avg_cnt (diag);
          remove_these.insert(name);
          fm_model->add_field(diag->get_diagnostic());
          m_diagnostics.push_back(diag);
        }
      }
    }

    EKAT_REQUIRE_MSG (add_these.size()>0 or remove_these.size()>0,
        "Error! We're stuck in an endless loop while processing output fields.\n"
        " - stream name: " + m_stream_name + "\n");

    // Remove fields we added to the FM, and add new diags we need
    for (const auto& n : remove_these) {
      remaining.erase(n);
    }
    for (const auto& n : add_these) {
      remaining.insert(n);
    }

    done = remaining.size()==0;
  }
}

std::vector<std::string> AtmosphereOutput::
get_var_dimnames (const FieldLayout& layout) const
{
  strvec_t dims;
  for (int i=0; i<layout.rank(); ++i) {
    const auto t = layout.tag(i);
    auto tag_name = m_io_grid->has_special_tag_name(t)
                  ? m_io_grid->get_special_tag_name(t)
                  : layout.names()[i];
    if (tag_name=="dim" or tag_name=="bin") {
      tag_name += std::to_string(layout.dim(i));
    }
    dims.push_back(tag_name); // Add dimensions string to vector of dims.
  }
  return dims;
}

// Instantiate the static member var
AtmosphereOutput::strmap_t<AtmosphereOutput::diag_ptr_type>
AtmosphereOutput::m_diag_repo;

} // namespace scream
