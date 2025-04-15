#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/vertical_remapper.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/field/field_utils.hpp"

#include <ekat/util/ekat_units.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/std_meta/ekat_std_utils.hpp>

#include <numeric>

namespace {
  // Helper lambda, to copy io string attributes. This will be used if any
  // remapper is created, to ensure atts set by atm_procs are not lost
  void transfer_io_str_atts  (const scream::Field& src, scream::Field& tgt) {
    const std::string io_string_atts_key ="io: string attributes";
    using stratts_t = std::map<std::string,std::string>;
    const auto& src_atts = src.get_header().get_extra_data<stratts_t>(io_string_atts_key);
          auto& dst_atts = tgt.get_header().get_extra_data<stratts_t>(io_string_atts_key);
    for (const auto& [name,val] : src_atts) {
      dst_atts[name] = val;
    }
  };
}

namespace scream
{

template<typename T>
bool has_duplicates (const std::vector<T>& c)
{
  std::set<T> s(c.begin(),c.end());
  return c.size()>s.size();
}

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm,
                  const std::vector<Field>& fields,
                  const std::shared_ptr<const grid_type>& grid)
 : m_comm (comm)
{
  // This version of AtmosphereOutput is for quick output of fields (no remaps, no time dim)
  m_avg_type = OutputAvgType::Instant;
  m_add_time_dim = false;

  // Create a FieldManager with the input fields
  auto fm = std::make_shared<FieldManager> (grid);
  for (auto f : fields) {
    fm->add_field(f);
    m_fields_names.push_back(f.name());
  }

  // No remaps: set all FM except the one for scorpio (created in init())
  m_field_mgrs[FromModel] = m_field_mgrs[AfterVertRemap] = m_field_mgrs[AfterHorizRemap] = fm;

  // Setup I/O structures
  init ();
}

AtmosphereOutput::
AtmosphereOutput (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const fm_type>& field_mgr,
                  const std::string& grid_name)
 : m_comm           (comm)
 , m_add_time_dim   (true)
{
  using vos_t = std::vector<std::string>;

  auto gm = field_mgr->get_grids_manager();

  // Figure out what kind of averaging is requested
  auto avg_type = params.get<std::string>("averaging_type");
  m_avg_type = str2avg(avg_type);
  EKAT_REQUIRE_MSG (m_avg_type!=OutputAvgType::Invalid,
      "Error! Unsupported averaging type '" + avg_type + "'.\n"
      "       Valid options: instant, Max, Min, Average. Case insensitive.\n");

  // By default, IO is done directly on the field mgr grid
  auto fm_grid = field_mgr->get_grids_manager()->get_grid(grid_name);
  std::string io_grid_name = fm_grid->name();
  if (params.isParameter("field_names")) {
    // This simple parameter list option does *not* allow to remap fields
    // to an io grid different from that of the field manager. In order to
    // use that functionality, you need the full syntax
    m_fields_names = params.get<vos_t>("field_names");
  } else if (params.isSublist("fields")){
    const auto& f_pl = params.sublist("fields");
    bool grid_found = false;
    for (const auto& grid_name : fm_grid->aliases()) {
      if (f_pl.isSublist(grid_name)) {
        grid_found = true;
        const auto& pl = f_pl.sublist(grid_name);
        if (pl.isType<vos_t>("field_names")) {
          m_fields_names = pl.get<vos_t>("field_names");
        } else if (pl.isType<std::string>("field_names")) {
          m_fields_names.resize(1, pl.get<std::string>("field_names"));
          if (m_fields_names[0]=="NONE") {
            m_fields_names.clear();
          }
        }

        // Check if the user wants to remap fields on a different grid first
        if (pl.isParameter("io_grid_name")) {
          io_grid_name = pl.get<std::string>("io_grid_name");
        }
        break;
      }
    }
    EKAT_REQUIRE_MSG (grid_found,
        "Error! Bad formatting of output yaml file. Missing 'fields->$grid_name` sublist.\n");
  }

  EKAT_REQUIRE_MSG (not has_duplicates(m_fields_names),
      "[AtmosphereOutput] Error! One of the output yaml files has duplicate field entries.\n"
      " - yaml file: " + params.name() + "\n"
      " - fields names; " + ekat::join(m_fields_names,",") + "\n");

  // Check if remapping and if so create the appropriate remapper
  // Note: We currently support three remappers
  //   - vertical remapping from file
  //   - horizontal remapping from file
  //   - online remapping which is setup using the create_remapper function
  const bool use_vertical_remap_from_file = params.isParameter("vertical_remap_file");
  const bool use_horiz_remap_from_file = params.isParameter("horiz_remap_file");
  const bool use_online_remapper = io_grid_name!=fm_grid->name();
  if (use_online_remapper) {
    EKAT_REQUIRE_MSG(!use_vertical_remap_from_file and !use_horiz_remap_from_file,
        "[AtmosphereOutput] Error! Online Dyn->PhysGLL remapping not supported along with vertical and/or horizontal remapping from file");
  }

  auto& fm_model = m_field_mgrs[FromModel];
  auto& fm_after_vr = m_field_mgrs[AfterVertRemap];
  auto& fm_after_hr = m_field_mgrs[AfterHorizRemap];

  // For simplicity, we create a "copy" of the input fm, so we can stuff also diags in it
  fm_model = std::make_shared<FieldManager>(field_mgr->get_grid(),RepoState::Closed);
  for (auto& fname : m_fields_names) {
    if (field_mgr->has_field(fname)) {
      fm_model->add_field(field_mgr->get_field(fname));
    } else {
      // This must be a diagnostic field. Crate the diag, and set the diag_field
      // in the "FromModel" field manager
      auto diag = create_diagnostic(fname);
      const auto& diag_fname = diag->get_diagnostic().name();
      m_diagnostics[diag_fname] = diag;

      // Note: some diag fields have a name different from what was used
      //       in the yaml file, so update the name with the actual
      //       diagnostic field name.
      // TODO: we should change this. The name used in the yaml file should
      //       MATCH the diag field name. This is confusing.
      fname = diag_fname;

      fm_model->add_field(m_diagnostics.at(fname)->get_diagnostic());
    }
  }

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

  // Setup remappers - if needed
  auto grid_after_vr = fm_grid;
  if (use_vertical_remap_from_file) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    auto vert_remap_file   = params.get<std::string>("vertical_remap_file");
    auto p_mid = field_mgr->get_field("p_mid");
    auto p_int = field_mgr->get_field("p_int");
    auto vert_remapper = std::make_shared<VerticalRemapper>(fm_model->get_grid(),vert_remap_file);
    vert_remapper->set_source_pressure (p_mid,p_int);
    vert_remapper->set_mask_value(m_fill_value);
    vert_remapper->set_extrapolation_type(VerticalRemapper::Mask); // both Top AND Bot
    m_vert_remapper = vert_remapper;

    grid_after_vr = m_vert_remapper->get_tgt_grid();
    fm_after_vr = std::make_shared<FieldManager>(grid_after_vr,RepoState::Closed);

    for (const auto& fname : m_fields_names) {
      auto src = fm_model->get_field(fname,fm_grid->name());
      auto tgt = m_vert_remapper->register_field_from_src(src);
      transfer_io_str_atts (src,tgt);
      fm_after_vr->add_field(tgt);
    }
  } else {
    // No vert remap. Simply alias the fm from the model
    fm_after_vr = fm_model;
  }

  // Online remapper and horizontal remapper follow a similar pattern so we check in the same conditional.
  auto grid_after_hr = grid_after_vr;
  if (use_online_remapper || use_horiz_remap_from_file) {
    // We build a remapper, to remap fields from the fm grid to the io grid
    if (use_horiz_remap_from_file) {
      // Construct the coarsening remapper
      auto horiz_remap_file   = params.get<std::string>("horiz_remap_file");
      m_horiz_remapper = std::make_shared<CoarseningRemapper>(grid_after_vr,horiz_remap_file,true);
    } else {
      // Construct a generic remapper (likely, Dyn->PhysicsGLL)
      grid_after_hr = gm->get_grid(io_grid_name);
      m_horiz_remapper = gm->create_remapper(grid_after_vr,grid_after_hr);
    }

    grid_after_hr = m_horiz_remapper->get_tgt_grid();
    fm_after_hr = std::make_shared<FieldManager>(grid_after_hr,RepoState::Closed);

    for (const auto& fname : m_fields_names) {
      auto src = fm_after_vr->get_field(fname,grid_after_vr->name());
      auto tgt = m_vert_remapper->register_field_from_src(src);
      transfer_io_str_atts (src,tgt);
      fm_after_hr->add_field(tgt);
    }
  } else {
    // No vert remap. Simply alias the fm after vr
    fm_after_hr = fm_after_vr;
  }

  // Setup I/O structures (including the scorpio FM)
  init ();
}

/* ---------------------------------------------------------- */
void AtmosphereOutput::restart (const std::string& filename)
{
  // Create an input stream on the fly, and init averaging data
  ekat::ParameterList res_params("Input Parameters");
  res_params.set<std::string>("filename",filename);

  AtmosphereInput hist_restart (res_params,m_field_mgrs[Scorpio]);
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
  for (const auto& fname : m_fields_names) {
    const auto& f = fm_after_hr->get_field(fname);
    const auto& fh = f.get_header();
    const auto& fid = fh.get_identifier();

    // Check if the field for scorpio can alias the field after hremap.
    // It can do so only for Instant output, and if the field is NOT a subfield ant NOT padded
    // Also, if we track avg cnt, we MUST add the mask_value extra data, to trigger fill-value logic
    // when calling Field's update methods
    if (m_avg_type!=OutputAvgType::Instant or
        fh.get_alloc_properties().get_padding()>0 or
        m_track_avg_cnt) {
      Field copy(fid);
      copy.allocate_view();
      transfer_io_str_atts (f,copy);
      if (m_track_avg_cnt) {
        copy.get_header().set_extra_data("mask_value",Real(m_fill_value));
      }
      fm_scorpio->add_field(copy);
    } else {
      fm_scorpio->add_field(f);
    }

    // Store the field layout, so that calls to setup_output_file are easier
    const auto& layout = fid.get_layout();
    m_layouts.emplace(fname,layout);

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
        "Error! Dimension " + tag_name + " on field " + fname + " has conflicting lengths.\n"
        "  - old length: " + std::to_string(m_dims[tag_name]) + "\n"
        "  - new length: " + std::to_string(dim_len) + "\n"
        "If same name applies to different dims (e.g. PhysicsGLL and PhysicsPG2 define "
        "\"ncol\" at different lengths), reset tag name for one of the grids.\n");
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
  for (auto& [name,diag] : m_diagnostics) {
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
  for (auto& [name,diag] : m_diagnostics) {
    compute_diagnostic(name,allow_invalid_fields);
  }

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
    std::set<std::string> avg_updated;
    for (const auto& [fname, avg_cnt_name] : m_field_to_avg_cnt_map) {
      if (avg_updated.count(avg_cnt_name)==1) {
        continue;
      }
      avg_updated.insert(avg_cnt_name);

      auto field = fm_after_hr->get_field(fname);
      auto count = fm_scorpio->get_field(avg_cnt_name);
      auto mask  = count.get_header().get_extra_data<Field>("mask");

      // Find where the field is NOT equal to m_fill_value
      compute_mask<Comparison::NE>(field,m_fill_value,mask);

      // mask=1 for "good" entries, and mask=0 otherwise.
      count.update(mask,1,1);

      // Handle writing the average count variables to file
      if (is_write_step) {
        // Bring data to host
        count.sync_to_host();

        auto func_start = std::chrono::steady_clock::now();
        scorpio::write_var(filename,avg_cnt_name,count.get_internal_view_data<int,Host>());
        auto func_finish = std::chrono::steady_clock::now();
        auto duration_loc = std::chrono::duration_cast<std::chrono::milliseconds>(func_finish - func_start);
        duration_write += duration_loc.count();

        // If it's an output step, for Avg we need to ensure count>threshold.
        // If count<=threshold, we set count=fill_value, so that fill_val propagates
        // to the output fields when we divide by count later
        if (output_step and m_avg_type==OutputAvgType::Average) {
          int min_count = m_avg_coeff_threshold*nsteps_since_last_output;

          // Recycle mask to find where count<thresh
          compute_mask<Comparison::LT>(count,min_count,mask);

          // Later, we divide fields by count. By setting count=1 where count<thresholt,
          // we can later do
          //   f.scale_inv(count); // Requires count!=0 anywhere
          //   f.deep_copy(fill_val,mask)
          count.deep_copy(1,mask);
        }
      }
    }
  }

  // Take care of updating and possibly writing fields.
  for (auto const& name : m_fields_names) {
    // Get all the info for this field.
    const auto& f_in  = fm_after_hr->get_field(name);
          auto& f_out = fm_scorpio->get_field(name);

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
      if (output_step and m_avg_type==OutputAvgType::Average) {
        // NOTE: we don't divide by the avg cnt for checkpoint output
        if (m_track_avg_cnt) {
          const auto& avg_cnt_name = m_field_to_avg_cnt_map.at(name);
          auto avg_cnt = fm_scorpio->get_field(avg_cnt_name);

          f_out.scale_inv(avg_cnt);

          const auto& mask = avg_cnt.get_header().get_extra_data<Field>("mask");
          f_out.deep_copy(m_fill_value,mask);
        } else {
          // Divide by steps count only when the summation is complete
          f_out.scale(1.0 / nsteps_since_last_output);
        }
      }

      // Bring data to host
      f_out.sync_to_host();

      // Write
      auto func_start = std::chrono::steady_clock::now();
      scorpio::write_var(filename,name,f_out.get_internal_view_data<Real,Host>());
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
res_dep_memory_footprint () const
{
  long long rdmf = 0;

  // Loop over ALL field mgr, and ALL fields in each of them.
  // Keep track of Field obj we parse, so we don't count them twice
  std::set<Field*> fields;
  for (const auto& [phase,fm] : m_field_mgrs) {
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
      if (phase==FromModel and m_diagnostics.count(fname)==0) {
        // We don't count fields from the model
        continue;
      }
      rdmf += fap.get_alloc_size();
    }
  }

  return rdmf;
}

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

  // Now create a Field to track the averaging count for this layout
  const auto& avg_cnt_suffix = m_field_to_avg_cnt_suffix[name];
  const auto tags = layout.tags();
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

  auto [it, inserted] = m_field_to_avg_cnt_map.emplace(name,avg_cnt_name);
  if (not inserted) {
    // We already created this avg cnt field
    return;
  }
  m_layouts.emplace(avg_cnt_name,layout);
  m_avg_cnt_names.push_back(avg_cnt_name);

  auto fm_scorpio = m_field_mgrs[Scorpio];

  // NOTE: while it seems right to use IntType as data type for cnt, we later use
  //       it in arithmetic updates of the avg field, where we need both fields
  //       to have the same data type
  auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier cnt_id (avg_cnt_name,layout,nondim,fm_scorpio->get_grid()->name(),DataType::IntType);
  Field cnt (cnt_id);
  cnt.allocate_view();

  // We will use a helper field for updating cnt, so store it inside the field header
  auto mask = cnt.clone(cnt.name()+"_mask");
  mask.get_header().set_extra_data("true_value",int(1));
  mask.get_header().set_extra_data("false_value",int(0));
  cnt.get_header().set_extra_data("mask",mask);
  
  fm_scorpio->add_field(cnt);
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
      value =  m_track_avg_cnt ? m_fill_value : 0.0;  break;
      break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unexpected averaging type.\n");
  }

  auto fm = m_field_mgrs[Scorpio];
  for (const auto& name : m_fields_names) {
    fm->get_field(name).deep_copy(value);
  }
  for (const auto& name : m_avg_cnt_names) {
    fm->get_field(name).deep_copy(0);
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
    const auto& f = m_field_mgrs[Scorpio]->get_field(name);
    const auto& fid  = f.get_header().get_identifier();
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
        scorpio::set_attribute(filename,name,"sub_fields",children_list);
      }

      // If tracking average count variables then add the name of the tracking variable for this variable
      if (m_track_avg_cnt) {
        const auto& lookup = m_field_to_avg_cnt_map.at(name);
        scorpio::set_attribute(filename,name,"averaging_count_tracker",lookup);
      }

      // Atm procs may have set some request for metadata.
      using stratts_t = std::map<std::string,std::string>;
      const auto& str_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
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
void AtmosphereOutput::set_diagnostics()
{
  auto fm_model = m_field_mgrs[FromModel];
  // Create all diagnostics
  for (const auto& fname : m_fields_names) {
    if (fm_model->has_field(fname)) {
      m_diagnostics[fname] = create_diagnostic(fname);
    }
  }
}

/* ---------------------------------------------------------- */
std::shared_ptr<AtmosphereDiagnostic>
AtmosphereOutput::create_diagnostic (const std::string& diag_field_name)
{
  // We need scream scope resolution, since this->create_diagnostic is hiding it
  auto fm_model = m_field_mgrs[FromModel];
  auto diag = scream::create_diagnostic(diag_field_name,fm_model->get_grid());

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
  for (const auto& freq : diag->get_required_field_requests()) {
    const auto& fname = freq.fid.name();
    if (not fm_model->has_field(fname)) {
      // This diag depends on another diag. Create and init the dependency
      if (m_diagnostics.count(fname)==0) {
        m_diagnostics[fname] = create_diagnostic(fname);
      }
      deps.push_back(fname);
    }
    diag->set_required_field (fm_model->get_field(fname));
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

} // namespace scream
