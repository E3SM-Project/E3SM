#include "share/io/eamxx_data_interpolation.hpp"

#include "share/remap/identity_remapper.hpp"
#include "share/remap/vertical_remapper.hpp"
#include "share/remap/refining_remapper_p2p.hpp"
#include "share/remap/iop_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/physics/physics_constants.hpp"

#include <ekat_team_policy_utils.hpp>

#include <filesystem>
#include <fstream>
#include <regex>
#include <unordered_set>

namespace scream{

DataInterpolation::
DataInterpolation (const std::shared_ptr<const AbstractGrid>& model_grid,
                   const std::vector<Field>& fields)
 : m_model_grid (model_grid)
 , m_fields     (fields)
{
  EKAT_REQUIRE_MSG (model_grid!=nullptr,
      "[DataInterpolation] Error! Invalid grid pointer.\n");

  m_nfields = m_fields.size();
  m_comm = model_grid->get_comm();

  using namespace ShortFieldTagsNames;
  m_input_files_dimnames[COL] = "ncol";
  m_input_files_dimnames[LEV] = "lev";
  e2str(LEV);

  // Initialize logger with a console logger for warnings (logs on all ranks)
  m_logger = console_logger(ekat::logger::LogLevel::warn);
}

void DataInterpolation::
set_logger (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  EKAT_REQUIRE_MSG (logger, "Error! Invalid logger pointer.\n");
  m_logger = logger;
}

void DataInterpolation::run (const util::TimeStamp& ts)
{
  EKAT_REQUIRE_MSG (m_data_initialized,
      "[DataInterpolation] Error! You must call 'init_data_interval' before calling 'run'.\n");

  // If we went past the current interval end, we need to update the end state
  if (not m_data_interval.contains(ts)) {
    shift_data_interval ();
  }

  // Perform the time interpolation: f_out = f_beg*alpha + f_end*(1-alpha),
  // where alpha = (ts-t_beg) / (t_end-t_beg).
  // NOTE: pay attention to time strategy, since for YearlyPeriodic you may
  //       have t_beg>t_end
  util::TimeInterval beg_to_ts (m_data_interval.beg,ts,m_data_interval.timeline);
  double alpha = beg_to_ts.length / m_data_interval.length;
  EKAT_REQUIRE_MSG (alpha>=0 and alpha<=1,
    "[DataInterpolation] Error! Input timestamp is outside the current data time interval.\n"
    "  data interval beg  ; " + m_data_interval.beg.to_string() + "\n"
    "  data interval end  ; " + m_data_interval.end.to_string() + "\n"
    "  input timestamp    ; " + ts.to_string() + "\n"
    "  interval length    : " + std::to_string(m_data_interval.length) + "\n"
    "  interpolation coeff: " + std::to_string(alpha) + "\n");

  for (int i=0; i<m_nfields; ++i) {
    const auto& beg = m_horiz_remapper_beg->get_tgt_field(i);
    const auto& end = m_horiz_remapper_end->get_tgt_field(i);
          auto  out = m_vert_remapper->get_src_field(i);

    out.deep_copy(beg);
    out.update(end,alpha,1-alpha);
  }

  // For Dynamic3D/Dynamic3D profile we also need to compute the source pressure profile
  // NOTE: this can't be done in the loop above, since p_data is not a "remapped"
  //       field in the vertical remapper (also, we need to use ad different ptr)
  if (m_vr_type==Dynamic3D) {
    // The pressure field is THE LAST registered in the horiz remappers
    const auto p_beg = m_horiz_remapper_beg->get_tgt_field(m_nfields);
    const auto p_end = m_horiz_remapper_end->get_tgt_field(m_nfields);

    auto p = m_helper_pressure_fields["p_data"];
    p.deep_copy(p_beg);
    p.update(p_end,alpha,1-alpha);
  } else if (m_vr_type==Dynamic3DRef) {
    // The surface pressure field is THE LAST registered in the horiz remappers
    const auto ps_beg = m_horiz_remapper_beg->get_tgt_field(m_nfields);
    const auto ps_end = m_horiz_remapper_end->get_tgt_field(m_nfields);

    auto p  = m_helper_pressure_fields["p_data"];
    auto ps = m_helper_pressure_fields["p_file"];
    ps.deep_copy(ps_beg);
    ps.update(ps_end,alpha,1-alpha);

    // Reconstruct reference p from ps, hyam, and hybm
    using KT = KokkosTypes<DefaultDevice>;
    using ExeSpace = typename KT::ExeSpace;
    using MemberType = typename KT::MemberType;
    using TPF = ekat::TeamPolicyFactory<ExeSpace>;
    using C = scream::physics::Constants<Real>;
    using PT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

    auto ps_v = ps.get_view<const Real*>();
    auto p_v  = p.get_view<PT**>();
    auto hyam = m_vert_remapper->get_src_grid()->get_geometry_data("hyam").get_view<const PT*>();
    auto hybm = m_vert_remapper->get_src_grid()->get_geometry_data("hybm").get_view<const PT*>();

    constexpr auto P0 = C::P0;

    const int ncols = ps_v.extent(0);
    const int num_vert_packs = p_v.extent(1);
    const auto policy = TPF::get_default_team_policy(ncols, num_vert_packs);

    Kokkos::parallel_for("spa_compute_p_src_loop", policy,
      KOKKOS_LAMBDA (const MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_vert_packs),
                           [&](const int k) {
        p_v(icol,k) = ps_v(icol) * hybm(k)  + P0 * hyam(k);
      });
    });
  }

  m_vert_remapper->remap_fwd();
}

void DataInterpolation::shift_data_interval ()
{
  m_curr_interval_idx.first = m_curr_interval_idx.second;
  m_curr_interval_idx.second = m_time_database.get_next_idx(m_curr_interval_idx.first);

  m_data_interval.advance(m_time_database.slices[m_curr_interval_idx.second].time);
  std::swap (m_horiz_remapper_beg,m_horiz_remapper_end);
  update_end_fields ();
}

void DataInterpolation::
update_end_fields ()
{
  // First, set the correct fields in the reader
  std::vector<Field> fields;
  for (int i=0; i<m_nfields; ++i) {
    fields.push_back(m_horiz_remapper_end->get_src_field(i));
  }

  if (m_vr_type==Dynamic3D or m_vr_type==Dynamic3DRef) {
    // We also need to read the src pressure profile
    fields.push_back(m_horiz_remapper_end->get_src_field(m_nfields));
  }
  m_reader->set_fields(fields);

  // If we're also changing the file, must (re)init the scorpio structures
  const auto& slice_beg = m_time_database.slices[m_curr_interval_idx.first];
  const auto& slice_end = m_time_database.slices[m_curr_interval_idx.second];
  if (m_reader->get_filename()!=slice_end.filename) {
    m_reader->reset_filename(slice_end.filename);
  }

  // Read and interpolate fields
  m_logger->info("[DataInterpolation] Reading end of interval fields.");
  m_logger->info(" - interval: [" + slice_beg.time.to_string() + ", " + slice_end.time.to_string() + "]");
  m_logger->info(" - filename: " + slice_end.filename);
  m_logger->info(" - file time idx: " + std::to_string(slice_end.time_idx));
  m_reader->read_variables(slice_end.time_idx);
  m_horiz_remapper_end->remap_fwd();
}

void DataInterpolation::
init_data_interval (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG (m_vert_remapper!=nullptr,
      "[DataInterpolation] Error! Cannot call 'init_data_interval' until after remappers creation.\n");

  register_fields_in_remappers ();

  // Create a bare reader. Fields and filename are set inside the update_end_fields call
  strvec_t fnames;
  for (auto f : m_fields) {
    fnames.push_back(f.name());
  }

  m_reader = std::make_shared<AtmosphereInput>(fnames,m_horiz_remapper_beg->get_src_grid());

  // Loop over all stored time slices to find an interval that contains t0
  auto t0_interval = m_time_database.find_interval(t0);
  const auto& t_beg = m_time_database.slices[t0_interval].time;

  // We need to read in the beg/end fields for the initial interval. However, our generic
  // framework can only load the end slice (since that's what we need at runtime).
  // So, load end state for t=t_beg, then call shift_data_interval
  // NOTE: don't compute length now, since beg time point is invalid (we don't need length yet).
  m_data_interval = util::TimeInterval (util::TimeStamp(),t_beg,m_time_database.timeline,false);
  m_curr_interval_idx.second = t0_interval;
  update_end_fields ();
  shift_data_interval ();

  m_data_initialized = true;
}

void DataInterpolation::
setup_time_database (const strvec_t& input_files,
                     const util::TimeLine timeline,
                     const util::TimeStamp& ref_ts)
{
  // Log the final list of files, so the user know if something went wrong (e.g. a bad regex)
  m_logger->debug("Setting up DataInerpolation object. List of input files:");
  for (const auto& fname : input_files) {
    m_logger->debug("  - " + fname);
  }

  // Make sure there are no repetitions
  auto num_unique_files = std::unordered_set(input_files.begin(),input_files.end()).size();
  EKAT_REQUIRE_MSG (num_unique_files==input_files.size(),
      "[DataInterpolation] Error! The input files list contains duplicates.\n"
      " - input_files:\n     " + ekat::join(input_files,"\n     ") + "\n");

  // We perform a bunch of checks on the input files
  namespace fs = std::filesystem;

  auto file_readable = [] (const std::string& fileName) {
    std::ifstream file(fileName);
    return file.good(); // Check if the file can be opened
  };

  // Read what time stamps we have in each file
  auto ts2str = [](const util::TimeStamp& t) { return t.to_string(); };
  std::vector<std::vector<util::TimeStamp>> times;
  for (const auto& fname : input_files) {
    EKAT_REQUIRE_MSG (file_readable(input_files.back()),
        "Error! One of the input files is not readable.\n"
        " - file   : " + input_files.back() + "\n");

    scorpio::register_file(fname,scorpio::Read);

    if (not scorpio::has_time_dim(fname)) {
      EKAT_REQUIRE_MSG (scorpio::has_dim(fname,"time"),
        "[DataInterpolation] Error! Input file does not contain a 'time' dimension.\n"
        " - file name: " + fname + "\n");
      scorpio::mark_dim_as_time(fname,"time");
    }
    auto file_times = scorpio::get_all_times(fname);
    EKAT_REQUIRE_MSG (file_times.size()>0,
        "[DataInterpolation] Error! Input file contains no time variable.\n"
        " - file name: " + fname + "\n");

    auto t_ref = ref_ts.is_valid() ? ref_ts : read_timestamp (fname,"reference_time_stamp");

    times.emplace_back();
    for (const auto& t : file_times) {
      times.back().push_back(t_ref + t*constants::seconds_per_day);
    }
    scorpio::release_file(fname);

    // Ensure time slices are sorted (it would make code messy otherwise)
    EKAT_REQUIRE_MSG (std::is_sorted(times.back().begin(),times.back().end()),
        "[DataInterpolation] Error! One of the input files has time slices not sorted.\n"
        " - file name  : " + fname + "\n"
        " - time stamps: " + ekat::join(times.back(),ts2str,", ") + "\n");
  }

  // Sort the files based on start date
  auto fileCmp = [](const std::vector<util::TimeStamp>& times1,
                    const std::vector<util::TimeStamp>& times2)
  {
    return times1.front() < times2.front();
  };
  std::sort(times.begin(),times.end(),fileCmp);

  // Setup the time database
  m_time_database.timeline = timeline;
  m_time_database.files = input_files;

  int nfiles = input_files.size();
  for (int i=0; i<nfiles; ++i) {
    int time_idx = 0;
    for (const auto& t : times[i]) {
      auto& slice = m_time_database.slices.emplace_back();
      slice.time = t;
      slice.filename = input_files[i];
      slice.time_idx = time_idx;
      ++time_idx;
    }

    if (i>0) {
      // Ensure files don't overlap (it would be a mess)
      const auto& prev = times[i-1];
      const auto& next = times[i];
      EKAT_REQUIRE_MSG (prev.back() < next.front(),
          "[DataInterpolation] Error! The input files contain overlapping time slices.\n"
          " - file1 name : " + input_files[i-1] + "\n"
          " - file2 name : " + input_files[i] + "\n"
          " - file1 times: " + ekat::join(prev,ts2str,", ") + "\n"
          " - file2 times: " + ekat::join(next,ts2str,", ") + "\n");
    }
  }

  // To avoid trouble in our logic of handling time stamps relationshipc,
  // we must ensure we have 2+ time slices overall
  EKAT_REQUIRE_MSG (m_time_database.size()>=2,
      "[DataInterpolation] Error! Input file(s) only contain 1 time slice overall.\n");

  m_time_db_created = true;
}

int DataInterpolation::TimeDatabase::
get_next_idx (int prev) const
{
  int next = prev+1;
  if (next >= size()) {
    EKAT_REQUIRE_MSG (timeline==util::TimeLine::YearlyPeriodic,
        "[TimeDatabase::get_next_idx] Error! Requesting slice that is past the database end.\n");
    next = next % size();
  }
  return next;
}

int DataInterpolation::TimeDatabase::
find_interval (const util::TimeStamp& t) const
{
  EKAT_REQUIRE_MSG (size()>1,
      "[TimeDatabase::find_interval] Error! The database has not been initialized yet.\n");

  auto contains = [&](int beg, int end, const util::TimeStamp& t) {
    const auto& t_beg = slices[beg].time;
    const auto& t_end = slices[end].time;
    util::TimeInterval t_int (t_beg,t_end,timeline);
    return t_int.contains(t);
  };
  int beg=0;
  int end=1;
  while (end<size()) {
    if (contains(beg,end,t)) {
      return beg;
    }
    beg = end;
    end = get_next_idx(beg);
  };

  // If we got here, no interval [i,i+1] contains t. If the timeline is Linear,
  // this is an error (in fact, there should have been an error before!).
  // But if the timeline is YearlyPeriodic, then it must be the case that
  // t is in the interval (slices[N],slices[0]).
  EKAT_REQUIRE_MSG (timeline==util::TimeLine::YearlyPeriodic and contains(size(),0,t),
      "[TimeDatabase::find_interval] Error! Could not locate interval containing input timestamp.\n"
      "  - input time : " + t.to_string() + "\n"
      "  - timeline   : " + (timeline==util::TimeLine::Linear ? "linear" : "yearly_periodic") + "\n"
      "  - first slice: " + slices.front().time.to_string() + "\n"
      "  - last slice : " + slices.back().time.to_string() + "\n"
      "Did you mean to use YearlyPeriodic timeline?\n");

  // Since we are in YearlyPeriodic timeline, the interval starts at the last slice
  return size()-1;
}

int DataInterpolation::
get_input_files_dimlen (const std::string& dimname) const
{
  // Retrieve a dim len from input file.
  // Also check that all files agree on that dim len
  int dimlen = -1;
  for (const auto& fname : m_time_database.files) {
    scorpio::register_file(fname,scorpio::Read);

    EKAT_REQUIRE_MSG (scorpio::has_dim(fname,dimname),
        "Error! Input file is missing '" + dimname + "' dimension.\n"
        "  - input file: " + fname + "\n");

    auto this_file_dimlen = scorpio::get_dimlen(fname,dimname);
    EKAT_REQUIRE_MSG (dimlen==-1 or dimlen==this_file_dimlen,
        "Error! Input files do not agree on '" + dimname + "' dimension length.\n"
        "  - file1: " + m_time_database.files.front() + "\n"
        "  - file2: " + fname + "\n"
        "  - file1 dim len: " + std::to_string(dimlen) + "\n"
        "  - file2 dim len: " + std::to_string(this_file_dimlen) + "\n");
    scorpio::release_file(fname);

    dimlen = this_file_dimlen;
  }
  return dimlen;
}

void DataInterpolation::
create_horiz_remappers (const std::string& map_file)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (m_horiz_remapper_beg==nullptr,
      "[DataInterpolation] Error! Horizontal remappers were already setup.\n");

  // Create hremap tgt grid
  int nlevs_data = get_input_files_dimlen (m_input_files_dimnames[LEV]);
  int ncols_data = get_input_files_dimlen (m_input_files_dimnames[COL]);
  m_grid_after_hremap = m_model_grid->clone("after_hremap",true);
  m_grid_after_hremap->reset_num_vertical_lev(nlevs_data);

  int ncols_model = m_model_grid->get_num_global_dofs();
  if (map_file!="") {
    m_horiz_remapper_beg = std::make_shared<RefiningRemapperP2P>(m_grid_after_hremap,map_file);
    m_horiz_remapper_end = std::make_shared<RefiningRemapperP2P>(m_grid_after_hremap,map_file);

    int map_ncols_src = m_horiz_remapper_beg->get_src_grid()->get_num_global_dofs();
    int map_ncols_tgt = m_horiz_remapper_beg->get_tgt_grid()->get_num_global_dofs();

    // Ensure that the map file was compatible with data/model grids
    EKAT_REQUIRE_MSG (map_ncols_src==ncols_data,
        "[DataInterpolation] Error! Map file src grid is incompatible with the input data.\n"
        " - map file ncols src: " + std::to_string(map_ncols_src) + "\n"
        " - data file ncols   : " + std::to_string(ncols_data) + "\n");
    EKAT_REQUIRE_MSG (map_ncols_tgt==ncols_model,
        "[DataInterpolation] Error! Map file tgt grid is incompatible with the model grid.\n"
        " - map file ncols tgt: " + std::to_string(map_ncols_tgt) + "\n"
        " - model file ncols  : " + std::to_string(ncols_model) + "\n");
  } else {
    // No hremap: 'ncols' from the data must match the model grid (nlev can differ; vremap is not set yet)
    EKAT_REQUIRE_MSG (ncols_data==ncols_model,
        "[DataInterpolation] Error! No horiz remap requested, but the 'ncol' dim from data does not match the model grid one.\n"
        " - model grid num global cols: " + std::to_string(ncols_model) + "\n"
        " - input data num global cols: " + std::to_string(ncols_data) + "\n");

    using IDR = IdentityRemapper;
    constexpr auto SAT = IDR::SrcAliasTgt;

    m_horiz_remapper_beg = std::make_shared<IDR>(m_grid_after_hremap,SAT);
    m_horiz_remapper_end = std::make_shared<IDR>(m_grid_after_hremap,SAT);
  }
}

void DataInterpolation::
create_horiz_remappers (const Real iop_lat, const Real iop_lon)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (m_horiz_remapper_beg==nullptr,
      "[DataInterpolation] Error! Horizontal remappers were already setup.\n");

  EKAT_REQUIRE_MSG (not std::isnan(iop_lat) and not std::isnan(iop_lon),
      "[DataInterpolation] Error! At least one between iop_lat and iop_lon appears to be invalid.\n"
      "  - iop_lat: " << iop_lat << "\n"
      "  - iop_lon: " << iop_lon << "\n");

  int nlevs_data = get_input_files_dimlen (m_input_files_dimnames[LEV]);
  int ncols_data = get_input_files_dimlen (m_input_files_dimnames[COL]);

  // Create grid for IO and load lat/lon field in IO grid from any data file
  auto data_grid = create_point_grid("data",ncols_data,nlevs_data,m_model_grid->get_comm());
  std::vector<Field> latlon = {
    data_grid->create_geometry_data("lat",data_grid->get_2d_scalar_layout()),
    data_grid->create_geometry_data("lon",data_grid->get_2d_scalar_layout())
  };
  AtmosphereInput latlon_reader (m_time_database.files.front(),data_grid,latlon);
  latlon_reader.read_variables();

  // Create iop remap tgt grid
  m_grid_after_hremap = m_model_grid->clone("after_hremap",true);
  m_grid_after_hremap->reset_num_vertical_lev(nlevs_data);

  // Create IOP remappers
  m_horiz_remapper_beg = std::make_shared<IOPRemapper>(data_grid,m_grid_after_hremap,iop_lat,iop_lon);
  m_horiz_remapper_end = std::make_shared<IOPRemapper>(data_grid,m_grid_after_hremap,iop_lat,iop_lon);
}

void DataInterpolation::
create_vert_remapper ()
{
  create_vert_remapper(VertRemapData());
}

void DataInterpolation::
create_vert_remapper (const VertRemapData& data)
{
  EKAT_REQUIRE_MSG (m_vert_remapper==nullptr,
      "[DataInterpolation] Error! Vertical remapper was already setup.\n");
  EKAT_REQUIRE_MSG (m_horiz_remapper_beg!=nullptr,
      "[DataInterpolation] Error! You must call `create_horiz_remappers` before `create_vert_remapper`.\n");

  m_vr_type = data.vr_type;

  if (m_vr_type==VRemapType::None) {
    // Not much to do. Set up a do-nothing remapper and return
    using IDR = IdentityRemapper;
    constexpr auto SAT = IDR::SrcAliasTgt;

    // If no vert remap is requested, model_grid and grid_after_hremap MUST have same nlevs
    int model_nlevs = m_model_grid->get_num_vertical_levels();
    int data_nlevs  = m_grid_after_hremap->get_num_vertical_levels();
    EKAT_REQUIRE_MSG (model_nlevs==data_nlevs,
        "[DataInterpolation] Error! No vertical remap was requested, but the 'lev' dim from file does not match the model grid one.\n"
        " - model grid num vert levels: " + std::to_string(model_nlevs) + "\n"
        " - input data num vert levels: " + std::to_string(data_nlevs) + "\n");
    m_vert_remapper = std::make_shared<IDR>(m_grid_after_hremap,SAT);
    return;
  }

  if (data.custom_remapper) {
    // The user provided a remapper. If vr_type is not Custom, this MUST be a VerticalRemapper.
    EKAT_REQUIRE_MSG (m_vr_type==Custom or std::dynamic_pointer_cast<VerticalRemapper>(data.custom_remapper),
        "[DataInterpolation] Error! Input vertical remapper MUST be of type VerticalRemapper (or derive from it) if the remap type is not 'Custom'.\n");

    // Sanity checks on src/tgt grids
    auto src_ncols = data.custom_remapper->get_src_grid()->get_num_global_dofs();
    auto tgt_ncols = data.custom_remapper->get_tgt_grid()->get_num_global_dofs();
    auto src_nlevs = data.custom_remapper->get_src_grid()->get_num_vertical_levels();
    auto tgt_nlevs = data.custom_remapper->get_tgt_grid()->get_num_vertical_levels();

    EKAT_REQUIRE_MSG (src_ncols==m_grid_after_hremap->get_num_global_dofs(),
        "[DataInterpolation] Error! Custom vert remapper src grid incompatible with the grid after horiz remap.\n"
        " - grid after horiz remap num global cols: " + std::to_string(m_grid_after_hremap->get_num_global_dofs()) + "\n"
        " - custom remapper num global cols       : " + std::to_string(src_ncols) + "\n");
    EKAT_REQUIRE_MSG (tgt_ncols==m_model_grid->get_num_global_dofs(),
        "[DataInterpolation] Error! Custom vert remapper tgt grid incompatible with the model grid.\n"
        " - model grid num global cols      : " + std::to_string(m_model_grid->get_num_global_dofs()) + "\n"
        " - custom remapper num global cols : " + std::to_string(tgt_ncols) + "\n");
    EKAT_REQUIRE_MSG (src_nlevs==m_grid_after_hremap->get_num_vertical_levels(),
        "[DataInterpolation] Error! Custom vert remapper src grid incompatible with the grid after horiz remap.\n"
        " - grid after horiz remap num levels: " + std::to_string(m_grid_after_hremap->get_num_vertical_levels()) + "\n"
        " - custom remapper num levels       : " + std::to_string(src_nlevs) + "\n");
    EKAT_REQUIRE_MSG (tgt_nlevs==m_model_grid->get_num_vertical_levels(),
        "[DataInterpolation] Error! Custom vert remapper tgt grid incompatible with the model grid.\n"
        " - model grid num levels      : " + std::to_string(m_model_grid->get_num_vertical_levels()) + "\n"
        " - custom remapper num levels : " + std::to_string(tgt_nlevs) + "\n");

    m_vert_remapper = data.custom_remapper;
  } else {
    auto s2et = [](const std::string& s) {
      if (s=="P0") {
        return VerticalRemapper::P0;
      } else if (s=="Mask") {
        return VerticalRemapper::Mask;
      } else {
        EKAT_ERROR_MSG (
            "Error! Invalid/unsupported extrapolation type.\n"
            " - input value : " + s + "\n"
            " - valid values: P0, Mask\n");
        return static_cast<VerticalRemapper::ExtrapType>(-1);
      }
    };

    // We need to build a vert remapper based on the input data.
    // Note: on src grid, we don't distinguish midpoints from interfaces, while on tgt we do.
    auto vremap = std::make_shared<VerticalRemapper>(m_grid_after_hremap,m_model_grid,true,false);

    vremap->set_extrapolation_type(s2et(data.extrap_top),VerticalRemapper::Top);
    vremap->set_extrapolation_type(s2et(data.extrap_bot),VerticalRemapper::Bot);

    m_vert_remapper = vremap;
  }

  // If the vremap type is not CUSTOM, we need to setup the source pressure profile
  if (m_vr_type!=Custom) {
    // If the remapper is not "custom", we MUST be able to cast down to VerticalRemapper
    auto vremap = std::dynamic_pointer_cast<VerticalRemapper>(m_vert_remapper);

    // Setup vertical pressure profiles (which can add 1 extra field to hremap)
    // NOTES:
    //  - both Dynamic3D and Dynamic3DRef use a 3d profile for the data
    //  - p_data is the full 3d pressure where data is defined, while p_file is the field
    //    we read from file. For Static1D and Dynamic3D they are the same, but for
    //    Dynamic3DRef, p_file is the surf pressure (2d), while p_data is the full 3d pmid
    auto p_layout = m_vr_type==Static1D ? m_grid_after_hremap->get_vertical_layout(true)
                                        : m_grid_after_hremap->get_3d_scalar_layout(true);
    auto& p_data = m_helper_pressure_fields ["p_data"];
    p_data = Field (FieldIdentifier("p_data",p_layout,ekat::units::Pa,m_grid_after_hremap->name()));
    p_data.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    p_data.allocate_view();
    if (m_vr_type==VRemapType::Dynamic3D) {
      // We load a full 3d profile, so p_file IS p_data
      m_helper_pressure_fields ["p_file"] = p_data.alias(data.pname);
    } else if (m_vr_type==VRemapType::Dynamic3DRef) {
      // We load the surface pressure, and reconstruct p_data via p=ps*hybm(k) + p0*hyam(k)
      auto& ps = m_helper_pressure_fields ["p_file"];
      ps = Field(FieldIdentifier(data.pname,m_grid_after_hremap->get_2d_scalar_layout(),ekat::units::Pa,m_grid_after_hremap->name()));
      ps.allocate_view();

      // We need to reconstruct the 3d pressure from ps, hybm, and hyam.
      // We read and store hyam/hybm in the vremap src grid
      auto layout = m_grid_after_hremap->get_vertical_layout(true);
      auto nondim = ekat::units::Units::nondimensional();
      DataType real_t = DataType::RealType;
      std::vector<Field> fields = {
        m_grid_after_hremap->create_geometry_data("hyam",layout,nondim,real_t,SCREAM_PACK_SIZE),
        m_grid_after_hremap->create_geometry_data("hybm",layout,nondim,real_t,SCREAM_PACK_SIZE)
      };
      AtmosphereInput hvcoord_reader (m_time_database.files.front(),m_grid_after_hremap,fields,true);
      hvcoord_reader.read_variables();
    } else if (m_vr_type==Static1D) {
      // Can load p now, since it's static
      std::vector<Field> fields = {p_data.alias(data.pname)};
      AtmosphereInput p_data_reader (m_time_database.files.front(),m_grid_after_hremap,fields,true);
      p_data_reader.read_variables();
    }
    vremap->set_source_pressure (m_helper_pressure_fields["p_data"],VerticalRemapper::Both);

    if (data.pint.is_allocated()) {
      vremap->set_target_pressure(data.pint,VerticalRemapper::Interfaces);
    }
    if (data.pmid.is_allocated()) {
      vremap->set_target_pressure(data.pmid,VerticalRemapper::Midpoints);
    }
  }
}

void DataInterpolation::register_fields_in_remappers ()
{
  // Register fields in the remappers. Vertical first, since we only have model-grid fields
  for (int i=0; i<m_nfields; ++i) {
    m_vert_remapper->register_field_from_tgt(m_fields[i]);
  }
  m_vert_remapper->registration_ends();

  for (int i=0; i<m_nfields; ++i) {
    const auto& f = m_vert_remapper->get_src_field(i);
    m_horiz_remapper_beg->register_field_from_tgt(f.clone(f.name(), m_horiz_remapper_beg->get_src_grid()->name()));
    m_horiz_remapper_end->register_field_from_tgt(f.clone(f.name(), m_horiz_remapper_end->get_src_grid()->name()));
  }
  if (m_vr_type==Dynamic3D or m_vr_type==Dynamic3DRef) {
    const auto& data_p = m_helper_pressure_fields["p_file"];
    m_horiz_remapper_beg->register_field_from_tgt(data_p.clone(data_p.name(), m_horiz_remapper_beg->get_src_grid()->name()));
    m_horiz_remapper_end->register_field_from_tgt(data_p.clone(data_p.name(), m_horiz_remapper_end->get_src_grid()->name()));
  }
  m_horiz_remapper_beg->registration_ends();
  m_horiz_remapper_end->registration_ends();
}

} // namespace scream
