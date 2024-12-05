#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/iop/intensive_observation_period.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

#include <numeric>

// Extend ekat mpi type for <Real,int> pairs,
// used for reduction of type MPI_MINLOC.
namespace ekat {
#ifdef SCREAM_DOUBLE_PRECISION
  template<>
  MPI_Datatype get_mpi_type<std::pair<double, int>> () {
    return MPI_DOUBLE_INT;
  }
#else
  template<>
  MPI_Datatype get_mpi_type<std::pair<float, int>> () {
    return MPI_FLOAT_INT;
  }
#endif
}

namespace scream {
namespace control {

IntensiveObservationPeriod::
IntensiveObservationPeriod(const ekat::Comm& comm,
                           const ekat::ParameterList& params,
                           const util::TimeStamp& run_t0,
                           const int model_nlevs,
                           const Field& hyam,
                           const Field& hybm)
{
  m_comm = comm;
  m_params = params;
  EKAT_REQUIRE_MSG(m_params.get<bool>("doubly_periodic_mode", false),
                   "Error! Currently doubly_periodic_mode is the only use case for "
	           "intensive observation period files.\n");

  EKAT_REQUIRE_MSG(m_params.isParameter("target_latitude") && m_params.isParameter("target_longitude"),
                   "Error! Using intensive observation period files requires "
                   "target_latitude and target_longitude be gives as parameters in "
                   "\"iop_options\" in the input yaml file.\n");
  const auto target_lat = m_params.get<Real>("target_latitude");
  const auto target_lon = m_params.get<Real>("target_longitude");
  EKAT_REQUIRE_MSG(-90 <= target_lat and target_lat <= 90,
                   "Error! IOP target_lat="+std::to_string(target_lat)+" outside of expected range [-90, 90].\n");
  EKAT_REQUIRE_MSG(0 <= target_lon and target_lon <= 360,
	           "Error! IOP target_lat="+std::to_string(target_lon)+" outside of expected range [0, 360].\n");

  // Set defaults for some parameters
  if (not m_params.isParameter("iop_srf_prop"))         m_params.set<bool>("iop_srf_prop",         false);
  if (not m_params.isParameter("iop_dosubsidence"))     m_params.set<bool>("iop_dosubsidence",     false);
  if (not m_params.isParameter("iop_coriolis"))         m_params.set<bool>("iop_coriolis",         false);
  if (not m_params.isParameter("iop_nudge_tq"))         m_params.set<bool>("iop_nudge_tq",         false);
  if (not m_params.isParameter("iop_nudge_uv"))         m_params.set<bool>("iop_nudge_uv",         false);
  if (not m_params.isParameter("iop_nudge_tq_low"))     m_params.set<Real>("iop_nudge_tq_low",     1050);
  if (not m_params.isParameter("iop_nudge_tq_high"))    m_params.set<Real>("iop_nudge_tq_high",    0);
  if (not m_params.isParameter("iop_nudge_tscale"))     m_params.set<Real>("iop_nudge_tscale",     10800);
  if (not m_params.isParameter("zero_non_iop_tracers")) m_params.set<bool>("zero_non_iop_tracers", false);

  // Store hybrid coords in helper fields
  m_helper_fields.insert({"hyam", hyam});
  m_helper_fields.insert({"hybm", hybm});

  // Use IOP file to initialize parameters
  // and timestepping information
  initialize_iop_file(run_t0, model_nlevs);
}

IntensiveObservationPeriod::
~IntensiveObservationPeriod ()
{
  const auto iop_file = m_params.get<std::string>("iop_file");
  scorpio::release_file(iop_file);
}

void IntensiveObservationPeriod::
initialize_iop_file(const util::TimeStamp& run_t0,
                    int model_nlevs)
{
  EKAT_REQUIRE_MSG(m_params.isParameter("iop_file"),
                   "Error! Using IOP requires defining an iop_file parameter.\n");

  const auto iop_file = m_params.get<std::string>("iop_file");

  // All the scorpio::has_var call can open the file on the fly, but since there
  // are a lot of those calls, for performance reasons we just open it now.
  // All the calls to register_file made on-the-fly inside has_var will be no-op.
  scorpio::register_file(iop_file,scorpio::FileMode::Read);

  // Lambda for allocating space and storing information for potential iop fields.
  // Inputs:
  //   - varnames:    Vector of possible variable names in the iop file.
  //                  First entry will be the variable name used when accessing in class
  //   - fl:          IOP field layout (acceptable ranks: 0, 1)
  //   - srf_varname: Name of surface variable potentially in iop file associated with iop variable.
  auto setup_iop_field = [&, this] (const vos&         varnames,
                                    const FieldLayout& fl,
                                    const std::string& srf_varname = "none") {
    EKAT_REQUIRE_MSG(fl.rank() == 0 || fl.rank() == 1,
                     "Error! IOP fields must have rank 0 or 1. "
                     "Attempting to setup "+varnames[0]+" with rank "
                     +std::to_string(fl.rank())+".\n");

    // Check if var exists in IOP file. Some variables will
    // need to check alternate names.
    const auto iop_varname = varnames[0];
    bool has_var = false;
    std::string file_varname = "";
    for (auto varname : varnames) {
      if (scorpio::has_var(iop_file, varname)) {
        has_var = true;
        file_varname = varname;
        break;
      };
    }
    if (has_var) {
      // Store if iop file has a different varname than the iop field
      if (iop_varname != file_varname) m_iop_file_varnames.insert({iop_varname, file_varname});
      // Store if variable contains a surface value in iop file
      if (scorpio::has_var(iop_file, srf_varname)) {
        m_iop_field_surface_varnames.insert({iop_varname, srf_varname});
      }
      // Store that the IOP variable is found in the IOP file
      m_iop_field_type.insert({iop_varname, IOPFieldType::FromFile});

      // Allocate field for variable
      FieldIdentifier fid(iop_varname, fl, ekat::units::Units::nondimensional(), "");
      const auto field_rank = fl.rank();
      EKAT_REQUIRE_MSG(field_rank <= 1,
                       "Error! Unexpected field rank "+std::to_string(field_rank)+" for iop file fields.\n");
      Field field(fid);
      if (fl.has_tag(FieldTag::LevelMidPoint) or fl.has_tag(FieldTag::LevelInterface)) {
        // Request packsize allocation for level layout
        field.get_header().get_alloc_properties().request_allocation(Pack::n);
      }
      field.allocate_view();
      m_iop_fields.insert({iop_varname, field});
    }
  };

  // Check if the following variables exist in the iop file

  // Scalar data
  FieldLayout fl_scalar({},{}); // Zero dim fields used for iop file scalars
  setup_iop_field({"Ps"},          fl_scalar);
  setup_iop_field({"Tg"},          fl_scalar);
  setup_iop_field({"lhflx", "lh"}, fl_scalar);
  setup_iop_field({"shflx", "sh"}, fl_scalar);

  // Level data
  FieldLayout fl_vector({FieldTag::LevelMidPoint}, {model_nlevs});
  setup_iop_field({"T"},        fl_vector, "Tsair");
  setup_iop_field({"q"},        fl_vector, "qsrf");
  setup_iop_field({"cld"},      fl_vector);
  setup_iop_field({"clwp"},     fl_vector);
  setup_iop_field({"divq"},     fl_vector, "divqsrf");
  setup_iop_field({"vertdivq"}, fl_vector, "vertdivqsrf");
  setup_iop_field({"NUMLIQ"},   fl_vector);
  setup_iop_field({"CLDLIQ"},   fl_vector);
  setup_iop_field({"CLDICE"},   fl_vector);
  setup_iop_field({"NUMICE"},   fl_vector);
  setup_iop_field({"divu"},     fl_vector, "divusrf");
  setup_iop_field({"divv"},     fl_vector, "divvsrf");
  setup_iop_field({"divT"},     fl_vector, "divtsrf");
  setup_iop_field({"vertdivT"}, fl_vector, "vertdivTsrf");
  setup_iop_field({"divT3d"},   fl_vector, "divT3dsrf");
  setup_iop_field({"u"},        fl_vector, "usrf");
  setup_iop_field({"u_ls"},     fl_vector, "usrf");
  setup_iop_field({"v"},        fl_vector, "vsrf");
  setup_iop_field({"v_ls"},     fl_vector, "vsrf");
  setup_iop_field({"Q1"},       fl_vector);
  setup_iop_field({"Q2"},       fl_vector);
  setup_iop_field({"omega"},    fl_vector, "Ptend");

  // Certain fields are required from the iop file
  EKAT_REQUIRE_MSG(has_iop_field("Ps"),
                   "Error! IOP file required to contain variable \"Ps\".\n");
  EKAT_REQUIRE_MSG(has_iop_field("T"),
                   "Error! IOP file required to contain variable \"T\".\n");
  EKAT_REQUIRE_MSG(has_iop_field("q"),
                   "Error! IOP file required to contain variable \"q\".\n");
  EKAT_REQUIRE_MSG(has_iop_field("divT"),
                   "Error! IOP file required to contain variable \"divT\".\n");
  EKAT_REQUIRE_MSG(has_iop_field("divq"),
                   "Error! IOP file required to contain variable \"divq\".\n");

  // Check for large scale winds and enfore "all-or-nothing" for u and v component
  const bool both_ls = (has_iop_field("u_ls") and has_iop_field("v_ls"));
  const bool neither_ls = (not (has_iop_field("u_ls") or has_iop_field("v_ls")));
  EKAT_REQUIRE_MSG(both_ls or neither_ls,
    "Error! Either u_ls and v_ls both defined in IOP file, or neither.\n");
  m_params.set<bool>("use_large_scale_wind", both_ls);

  // Require large scale winds if using Coriolis forcing
  if (m_params.get<bool>("iop_coriolis")) {
    EKAT_REQUIRE_MSG(m_params.get<bool>("use_large_scale_wind"),
                     "Error! Large scale winds required for coriolis forcing.\n");
  }

  // If we have the vertical component of T/Q forcing, define 3d forcing as a computed field.
  if (has_iop_field("vertdivT")) {
    FieldIdentifier fid("divT3d", fl_vector, ekat::units::Units::nondimensional(), "");
    Field field(fid);
    field.get_header().get_alloc_properties().request_allocation(Pack::n);
    field.allocate_view();
    m_iop_fields.insert({"divT3d", field});
    m_iop_field_type.insert({"divT3d", IOPFieldType::Computed});
  }
  if (has_iop_field("vertdivq")) {
    FieldIdentifier fid("divq3d", fl_vector, ekat::units::Units::nondimensional(), "");
    Field field(fid);
    field.get_header().get_alloc_properties().request_allocation(Pack::n);
    field.allocate_view();
    m_iop_fields.insert({"divq3d", field});
    m_iop_field_type.insert({"divq3d", IOPFieldType::Computed});
  }

  // Enforce that 3D forcing is all-or-nothing for T and q.
  const bool both_3d_forcing = (has_iop_field("divT3d") and has_iop_field("divq3d"));
  const bool neither_3d_forcing = (not (has_iop_field("divT3d") or has_iop_field("divq3d")));
  EKAT_REQUIRE_MSG(both_3d_forcing or neither_3d_forcing,
    "Error! Either T and q both have 3d forcing, or neither have 3d forcing.\n");
  m_params.set<bool>("use_3d_forcing", both_3d_forcing);

  // Initialize time information
  int bdate;
  std::string bdate_name;
  if      (scorpio::has_var(iop_file, "bdate"))    bdate_name = "bdate";
  else if (scorpio::has_var(iop_file, "basedate")) bdate_name = "basedate";
  else if (scorpio::has_var(iop_file, "nbdate"))   bdate_name = "nbdate";
  else EKAT_ERROR_MSG("Error! No valid name for bdate in "+iop_file+".\n");

  scorpio::read_var(iop_file, bdate_name, &bdate);

  int yr=bdate/10000;
  int mo=(bdate/100) - yr*100;
  int day=bdate - (yr*10000+mo*100);
  m_time_info.iop_file_begin_time = util::TimeStamp(yr,mo,day,0,0,0);

  std::string time_dimname;
  if      (scorpio::has_dim(iop_file, "time")) time_dimname = "time";
  else if (scorpio::has_dim(iop_file, "tsec")) time_dimname = "tsec";
  else EKAT_ERROR_MSG("Error! No valid dimension for tsec in "+iop_file+".\n");

  // When we read vars, "time" must be treated as unlimited, to avoid issues
  if (not scorpio::is_dim_unlimited(iop_file,time_dimname)) {
    scorpio::pretend_dim_is_unlimited(iop_file,time_dimname);
  }

  const auto ntimes = scorpio::get_dimlen(iop_file, time_dimname);
  m_time_info.iop_file_times_in_sec = view_1d_host<int>("iop_file_times", ntimes);
  for (int t=0; t<ntimes; ++t) {
    scorpio::read_var(iop_file,"tsec",&m_time_info.iop_file_times_in_sec(t),t);
  }

  // Check that lat/lon from iop file match the targets in parameters. Note that
  // longitude may be negtive in the iop file, we convert to positive before checking.
  const auto nlats = scorpio::get_dimlen(iop_file, "lat");
  const auto nlons = scorpio::get_dimlen(iop_file, "lon");
  EKAT_REQUIRE_MSG(nlats==1 and nlons==1, "Error! IOP data file requires a single lat/lon pair.\n");
  Real iop_file_lat, iop_file_lon;

  scorpio::read_var(iop_file,"lat",&iop_file_lat);
  scorpio::read_var(iop_file,"lon",&iop_file_lon);

  const Real rel_lat_err = std::fabs(iop_file_lat - m_params.get<Real>("target_latitude"))/
                             std::max(std::fabs(m_params.get<Real>("target_latitude")),(Real)0.1);
  const Real rel_lon_err = std::fabs(std::fmod(iop_file_lon + 360.0, 360.0)-m_params.get<Real>("target_longitude"))/
                             std::max(m_params.get<Real>("target_longitude"),(Real)0.1);
  EKAT_REQUIRE_MSG(rel_lat_err < std::numeric_limits<float>::epsilon(),
                   "Error! IOP file variable \"lat\" does not match target_latitude from IOP parameters.\n");
  EKAT_REQUIRE_MSG(rel_lon_err < std::numeric_limits<float>::epsilon(),
                   "Error! IOP file variable \"lon\" does not match target_longitude from IOP parameters.\n");

  // Store iop file pressure as helper field with dimension lev+1.
  // Load the first lev entries from iop file, the lev+1 entry will
  // be set when reading iop data.
  EKAT_REQUIRE_MSG(scorpio::has_var(iop_file, "lev"),
                    "Error! Using IOP file requires variable \"lev\".\n");
  const auto file_levs = scorpio::get_dimlen(iop_file, "lev");
  FieldIdentifier fid("iop_file_pressure",
                      FieldLayout({FieldTag::LevelMidPoint}, {file_levs+1}),
                      ekat::units::Units::nondimensional(),
                      "");
  Field iop_file_pressure(fid);
  iop_file_pressure.get_header().get_alloc_properties().request_allocation(Pack::n);
  iop_file_pressure.allocate_view();
  auto data = iop_file_pressure.get_view<Real*, Host>().data();
  scorpio::read_var(iop_file,"lev",data);

  // Convert to pressure to millibar (file gives pressure in Pa)
  for (int ilev=0; ilev<file_levs; ++ilev) data[ilev] /= 100;
  iop_file_pressure.sync_to_dev();
  m_helper_fields.insert({"iop_file_pressure", iop_file_pressure});

  // Create model pressure helper field (values will be computed
  // in read_iop_file_data())
  FieldIdentifier model_pres_fid("model_pressure",
                                  fl_vector,
                                  ekat::units::Units::nondimensional(), "");
  Field model_pressure(model_pres_fid);
  model_pressure.get_header().get_alloc_properties().request_allocation(Pack::n);
  model_pressure.allocate_view();
  m_helper_fields.insert({"model_pressure", model_pressure});
}

void IntensiveObservationPeriod::
setup_io_info(const std::string& file_name,
              const grid_ptr& grid)
{
  const auto grid_name = grid->name();

  // Create io grid if doesn't exist
  if (m_io_grids.count(grid_name) == 0) {
    // IO grid needs to have ncol dimension equal to the IC/topo file
    const auto nc_file_ncols = scorpio::get_dimlen(file_name, "ncol");
    const auto nlevs = grid->get_num_vertical_levels();
    m_io_grids[grid_name] = create_point_grid(grid_name,
                                              nc_file_ncols,
                                              nlevs,
                                              m_comm);
  }

  // Store closest lat/lon info for this grid if doesn't exist
  if (m_lat_lon_info.count(grid_name) == 0) {
    const auto& io_grid = m_io_grids[grid_name];

    // Create lat/lon fields
    const auto ncols = io_grid->get_num_local_dofs();
    std::vector<Field> fields;

    FieldIdentifier lat_fid("lat",
                            FieldLayout({FieldTag::Column},{ncols}),
                            ekat::units::Units::nondimensional(),
                            grid_name);
    Field lat_f(lat_fid);
    lat_f.allocate_view();
    fields.push_back(lat_f);

    FieldIdentifier lon_fid("lon",
                            FieldLayout({FieldTag::Column},{ncols}),
                            ekat::units::Units::nondimensional(),
                            grid_name);
    Field lon_f(lon_fid);
    lon_f.allocate_view();
    fields.push_back(lon_f);

    // Read from file
    AtmosphereInput file_reader(file_name, io_grid, fields);
    file_reader.read_variables();
    file_reader.finalize();

    // Find column index of closest lat/lon to target_lat/lon params
    auto lat_v = fields[0].get_view<Real*>();
    auto lon_v = fields[1].get_view<Real*>();
    const auto target_lat = m_params.get<Real>("target_latitude");
    const auto target_lon = m_params.get<Real>("target_longitude");
    using minloc_t = Kokkos::MinLoc<Real,int>;
    using minloc_value_t = typename minloc_t::value_type;
    minloc_value_t minloc;
    Kokkos::parallel_reduce(ncols, KOKKOS_LAMBDA (int icol, minloc_value_t& result) {
      auto dist = std::abs(lat_v(icol)-target_lat)+std::abs(lon_v(icol)-target_lon);
      if(dist<result.val) {
        result.val = dist;
        result.loc = icol;
      }
    }, minloc_t(minloc));

    // Find processor with closest lat/lon match
    const auto my_rank = m_comm.rank();
    std::pair<Real, int> min_dist_and_rank = {minloc.val, my_rank};
    m_comm.all_reduce<std::pair<Real, int>>(&min_dist_and_rank, 1, MPI_MINLOC);

    // Broadcast closest lat/lon values to all ranks
    const auto lat_v_h = lat_f.get_view<Real*,Host>();
    const auto lon_v_h = lon_f.get_view<Real*,Host>();
    auto local_column_idx = minloc.loc;
    auto min_dist_rank = min_dist_and_rank.second;
    Real lat_lon_vals[2];
    if (my_rank == min_dist_rank) {
      lat_lon_vals[0] = lat_v_h(local_column_idx);
      lat_lon_vals[1] = lon_v_h(local_column_idx);
    }
    m_comm.broadcast(lat_lon_vals, 2, min_dist_rank);

    // Set local_column_idx=-1 for mpi ranks not containing minimum lat/lon distance
    if (my_rank != min_dist_rank) local_column_idx = -1;

    // Store closest lat/lon info for this grid, used later when reading ICs
    m_lat_lon_info[grid_name] = ClosestLatLonInfo{lat_lon_vals[0], lat_lon_vals[1], min_dist_rank, local_column_idx};
  }
}

void IntensiveObservationPeriod::
read_fields_from_file_for_iop (const std::string& file_name,
                               const vos& field_names_nc,
                               const vos& field_names_eamxx,
                               const util::TimeStamp& initial_ts,
                               const field_mgr_ptr field_mgr,
                               const int time_index)
{
  const auto dummy_units = ekat::units::Units::nondimensional();

  EKAT_REQUIRE_MSG(field_names_nc.size()==field_names_eamxx.size(),
                  "Error! Field name arrays must have same size.\n");

  if (field_names_nc.size()==0) {
    return;
  }

  const auto& grid_name = field_mgr->get_grid()->name();
  EKAT_REQUIRE_MSG(m_io_grids.count(grid_name) > 0,
                   "Error! Attempting to read IOP initial conditions on "
                   +grid_name+" grid, but m_io_grid entry has not been created.\n");
  EKAT_REQUIRE_MSG(m_lat_lon_info.count(grid_name) > 0,
                   "Error! Attempting to read IOP initial conditions on "
                   +grid_name+" grid, but m_lat_lon_info entry has not been created.\n");

  auto io_grid = m_io_grids[grid_name];
  if (grid_name=="Physics GLL" && scorpio::has_dim(file_name,"ncol_d")) {
    // If we are on GLL grid, and nc file contains "ncol_d" dimension,
    // we need to reset COL dim tag
    using namespace ShortFieldTagsNames;
    auto grid = io_grid->clone(io_grid->name(),true);
    grid->reset_field_tag_name(COL,"ncol_d");
    io_grid = grid;
  }

  // Create vector of fields with correct dimensions to read from file
  std::vector<Field> io_fields;
  for (size_t i=0; i<field_names_nc.size(); ++i) {
    const auto& nc_name    = field_names_nc[i];
    const auto& eamxx_name = field_names_eamxx[i];
    const auto& fm_field = eamxx_name!=nc_name
                           ?
                           field_mgr->get_field(eamxx_name).alias(nc_name)
                           :
                           field_mgr->get_field(eamxx_name);
    auto fm_fid = fm_field.get_header().get_identifier();
    EKAT_REQUIRE_MSG(fm_fid.get_layout().tag(0)==FieldTag::Column,
                     "Error! IOP inputs read from IC/topo file must have Column "
                     "as first dim tag.\n");

    // Set first dimension to match input file
    FieldLayout io_fl = fm_fid.get_layout();
    io_fl.reset_dim(0,io_grid->get_num_local_dofs());
    FieldIdentifier io_fid(fm_fid.name(), io_fl, fm_fid.get_units(), io_grid->name());
    Field io_field(io_fid);
    io_field.allocate_view();
    io_fields.push_back(io_field);
  }

  // Read data from file
  AtmosphereInput file_reader(file_name,io_grid,io_fields);
  file_reader.read_variables(time_index);
  file_reader.finalize();

  // For each field, broadcast data from closest lat/lon column to all processors
  // and copy data into each field's column in the field manager.
  for (size_t i=0; i<io_fields.size(); ++i) {
    const auto& io_field = io_fields[i];
    const auto& fname = field_names_eamxx[i];
    auto& fm_field = field_mgr->get_field(fname);

    // Create a temporary field to store the data from the
    // single column of the closest lat/lon pair
    const auto io_fid = io_field.get_header().get_identifier();
    FieldLayout col_data_fl = io_fid.get_layout().clone().strip_dim(0);
    FieldIdentifier col_data_fid("col_data", col_data_fl, dummy_units, "");
    Field col_data(col_data_fid);
    col_data.allocate_view();

    // MPI rank with closest column index store column data
    const auto mpi_rank_with_col = m_lat_lon_info[grid_name].mpi_rank_of_closest_column;
    if (m_comm.rank() == mpi_rank_with_col) {
      const auto col_idx_with_data = m_lat_lon_info[grid_name].local_column_index_of_closest_column;
      col_data.deep_copy<Host>(io_field.subfield(0,col_idx_with_data));
    }

    // Broadcast column data to all other ranks
    const auto col_size = col_data.get_header().get_identifier().get_layout().size();
    m_comm.broadcast(col_data.get_internal_view_data<Real,Host>(), col_size, mpi_rank_with_col);

    // Copy column data to all columns in field manager field
    const auto ncols = fm_field.get_header().get_identifier().get_layout().dim(0);
    for (auto icol=0; icol<ncols; ++icol) {
      fm_field.subfield(0,icol).deep_copy<Host>(col_data);
    }

    // Sync fields to device
    fm_field.sync_to_dev();

    // Set the initial time stamp on FM fields
    fm_field.get_header().get_tracking().update_time_stamp(initial_ts);
  }
}

void IntensiveObservationPeriod::
read_iop_file_data (const util::TimeStamp& current_ts)
{
  // Query to see if we need to load data from IOP file.
  // If we are still in the time interval as the previous
  // read from iop file, there is no need to reload data.
  const auto iop_file_time_idx = m_time_info.get_iop_file_time_idx(current_ts);
  EKAT_REQUIRE_MSG(iop_file_time_idx >= m_time_info.time_idx_of_current_data,
                   "Error! Attempting to read previous iop file data time index.\n");
  if (iop_file_time_idx == m_time_info.time_idx_of_current_data) return;

  const auto iop_file = m_params.get<std::string>("iop_file");
  const auto file_levs = scorpio::get_dimlen(iop_file, "lev");
  const auto iop_file_pressure = m_helper_fields["iop_file_pressure"];
  const auto model_pressure = m_helper_fields["model_pressure"];
  const auto surface_pressure = m_iop_fields["Ps"];

  // Loop through iop fields, if any rank 1 fields are loaded from file,
  // we need to gather information for vertical interpolation
  bool has_level_data = false;
  for (auto& it : m_iop_fields) {
    if (it.second.rank() == 1
        and
        m_iop_field_type.at(it.first)==IOPFieldType::FromFile) {
      has_level_data = true;
      break;
    }
  }

  // Compute values and indices associate with pressure for interpolating data (if necessary).
  int adjusted_file_levs;
  int iop_file_start;
  int iop_file_end;
  int model_start;
  int model_end;
  if (has_level_data) {
    // Load surface pressure (Ps) from iop file
    auto ps_data = surface_pressure.get_view<Real, Host>().data();

    scorpio::read_var(iop_file,"Ps",ps_data,iop_file_time_idx);
    surface_pressure.sync_to_dev();

    // Read in IOP lev data
    auto data = iop_file_pressure.get_view<Real*, Host>().data();
    scorpio::read_var(iop_file,"lev",data);

    // Convert to pressure to millibar (file gives pressure in Pa)
    for (int ilev=0; ilev<file_levs; ++ilev) data[ilev] /= 100;
    iop_file_pressure.sync_to_dev();

    // Pre-process file pressures, store number of file levels
    // where the last level is the first level equal to surface pressure.
    const auto iop_file_pres_v = iop_file_pressure.get_view<Real*>();
    // Sanity check
    EKAT_REQUIRE_MSG(file_levs+1 == iop_file_pressure.get_header().get_identifier().get_layout().dim(0),
                    "Error! Unexpected size for helper field \"iop_file_pressure\"\n");
    const auto& Ps = surface_pressure.get_view<const Real>();
    Kokkos::parallel_reduce(file_levs+1, KOKKOS_LAMBDA (const int ilev, int& lmin) {
      if (ilev == file_levs) {
        // Add surface pressure to last iop file pressure entry
        iop_file_pres_v(ilev) = Ps()/100;
      }
      if (iop_file_pres_v(ilev) > Ps()/100) {
        // Set upper bound on pressure values
        iop_file_pres_v(ilev) = Ps()/100;
      }
      if (iop_file_pres_v(ilev) == Ps()/100) {
        // Find minimum number of levels where the final
        // level would contain the largest value.
        if (ilev < lmin) lmin = ilev+1;
      }
    }, Kokkos::Min<int>(adjusted_file_levs));

    EKAT_REQUIRE_MSG(adjusted_file_levs > 1,
                     "Error! Pressures in iop file "+iop_file+" is are inccorrectly set. "
                     "Surface pressure \"Ps\" (converted to millibar) should be greater "
                     "than at least the 1st entry in midpoint pressures \"lev\".\n");

    // Compute model pressure levels
    const auto model_pres_v = model_pressure.get_view<Real*>();
    const auto model_nlevs = model_pressure.get_header().get_identifier().get_layout().dim(0);
    const auto hyam_v = m_helper_fields["hyam"].get_view<const Real*>();
    const auto hybm_v = m_helper_fields["hybm"].get_view<const Real*>();
    Kokkos::parallel_for(model_nlevs, KOKKOS_LAMBDA (const int ilev) {
      model_pres_v(ilev) = 1000*hyam_v(ilev) + Ps()*hybm_v(ilev)/100;
    });

    // Find file pressure levels just outside the range of model pressure levels
    Kokkos::parallel_reduce(adjusted_file_levs, KOKKOS_LAMBDA (const int& ilev, int& lmax, int& lmin) {
      if (iop_file_pres_v(ilev) <= model_pres_v(0) && ilev > lmax) {
        lmax = ilev;
      }
      if (iop_file_pres_v(ilev) >= model_pres_v(model_nlevs-1) && ilev+1 < lmin) {
        lmin = ilev+1;
      }
    },
    Kokkos::Max<int>(iop_file_start),
    Kokkos::Min<int>(iop_file_end));

    // If no file pressures are found outide the reference pressure range, set to file level endpoints
    if (iop_file_start == Kokkos::reduction_identity<int>::max()) iop_file_start = 0;
    if (iop_file_end   == Kokkos::reduction_identity<int>::min()) iop_file_end = adjusted_file_levs;

    // Find model pressure levels just inside range of file pressure levels
    Kokkos::parallel_reduce(model_nlevs, KOKKOS_LAMBDA (const int& ilev, int& lmin, int& lmax) {
      if (model_pres_v(ilev) >= iop_file_pres_v(iop_file_start) && ilev < lmin) {
        lmin = ilev;
      }
      if (model_pres_v(ilev) <= iop_file_pres_v(iop_file_end-1) && ilev+1 > lmax) {
        lmax = ilev+1;
      }
    },
    Kokkos::Min<int>(model_start),
    Kokkos::Max<int>(model_end));

    // If not reference pressures are found inside file pressures, set to model level endpoints
    if (model_start == Kokkos::reduction_identity<int>::min()) model_start = model_nlevs-1;
    if (model_end   == Kokkos::reduction_identity<int>::max()) model_end = 1;
  }

  // Loop through fields and store data from file
  for (auto& it : m_iop_fields) {
    auto fname = it.first;
    auto field = it.second;

    // If this is a computed field, do not attempt to load from file
    if (m_iop_field_type.at(fname)==IOPFieldType::Computed) continue;

    // File may use different varname than IOP class
    auto file_varname = (m_iop_file_varnames.count(fname) > 0) ? m_iop_file_varnames[fname] : fname;

    if (field.rank()==0) {
      // For scalar data, read iop file variable directly into field data
      auto data = field.get_view<Real, Host>().data();
      scorpio::read_var(iop_file,file_varname,data,iop_file_time_idx);
      field.sync_to_dev();
    } else if (field.rank()==1) {
      // Create temporary fields for reading iop file variables. We use
      // adjusted_file_levels (computed above) which contains an unset
      // value for surface.
      FieldIdentifier fid(file_varname+"_iop_file",
                          FieldLayout({FieldTag::LevelMidPoint},
                          {adjusted_file_levs}),
                          ekat::units::Units::nondimensional(),
                          "");
      Field iop_file_field(fid);
      iop_file_field.get_header().get_alloc_properties().request_allocation(Pack::n);
      iop_file_field.allocate_view();

      // Read data from iop file.
      std::vector<Real> data(file_levs);
      scorpio::read_var(iop_file,file_varname,data.data(),iop_file_time_idx);

      // Copy first adjusted_file_levs-1 values to field
      auto iop_file_v_h = iop_file_field.get_view<Real*,Host>();
      for (int ilev=0; ilev<adjusted_file_levs-1; ++ilev) iop_file_v_h(ilev) = data[ilev];

      // Set or compute surface value
      const auto has_srf = m_iop_field_surface_varnames.count(fname)>0;
      if (has_srf) {
        const auto srf_varname = m_iop_field_surface_varnames[fname];
        scorpio::read_var(iop_file,srf_varname,&iop_file_v_h(adjusted_file_levs-1),iop_file_time_idx);
      } else {
        // No surface value exists, compute surface value
        const auto dx = iop_file_v_h(adjusted_file_levs-2) - iop_file_v_h(adjusted_file_levs-3);
        if (dx == 0) iop_file_v_h(adjusted_file_levs-1) = iop_file_v_h(adjusted_file_levs-2);
        else {
          iop_file_pressure.sync_to_host();
          const auto iop_file_pres_v_h = iop_file_pressure.get_view<const Real*, Host>();
          const auto dy = iop_file_pres_v_h(adjusted_file_levs-2) - iop_file_pres_v_h(adjusted_file_levs-3);
          const auto scale = dy/dx;

          iop_file_v_h(adjusted_file_levs-1) =
            (iop_file_pres_v_h(adjusted_file_levs-1)-iop_file_pres_v_h(adjusted_file_levs-2))/scale
            + iop_file_v_h(adjusted_file_levs-2);
        }
      }
      iop_file_field.sync_to_dev();

      // Vertically interpolate iop file data to iop fields.
      // Note: ekat lininterp requires packs. Use 1d packs here
      // to easily mask out levels which we do not want to interpolate.
      const auto iop_file_pres_v = iop_file_pressure.get_view<const Pack1d*>();
      const auto model_pres_v = model_pressure.get_view<const Pack1d*>();
      const auto iop_file_v = iop_file_field.get_view<const Pack1d*>();
      auto iop_field_v = field.get_view<Pack1d*>();

      const auto nlevs_input = iop_file_end - iop_file_start;
      const auto nlevs_output = model_end - model_start;
      const auto total_nlevs = field.get_header().get_identifier().get_layout().dim(0);

      ekat::LinInterp<Real,Pack1d::n> vert_interp(1, nlevs_input, nlevs_output);
      const auto policy = ESU::get_default_team_policy(1, total_nlevs);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
        const auto x_src  = Kokkos::subview(iop_file_pres_v, Kokkos::pair<int,int>(iop_file_start,iop_file_end));
        const auto x_tgt  = Kokkos::subview(model_pres_v, Kokkos::pair<int,int>(model_start,model_end));
        const auto input  = Kokkos::subview(iop_file_v, Kokkos::pair<int,int>(iop_file_start,iop_file_end));
        const auto output = Kokkos::subview(iop_field_v, Kokkos::pair<int,int>(model_start,model_end));

        vert_interp.setup(team, x_src, x_tgt);
        vert_interp.lin_interp(team, x_src, x_tgt, input, output);
      });
      Kokkos::fence();

      // For certain fields we need to make sure to fill in the ends of
      // the interpolated region with the value at model_start/model_end
      if (fname == "T"    || fname == "q" || fname == "u" ||
          fname == "u_ls" || fname == "v" || fname == "v_ls") {
        Kokkos::parallel_for(Kokkos::RangePolicy<>(0, model_start+1),
			     KOKKOS_LAMBDA (const int ilev) {
			       iop_field_v(ilev) = iop_file_v(0);
			     });
        Kokkos::parallel_for(Kokkos::RangePolicy<>(model_end-1, total_nlevs),
			     KOKKOS_LAMBDA (const int ilev) {
			       iop_field_v(ilev) = iop_file_v(adjusted_file_levs-1);
			     });
      }
    }
  }

  // Calculate 3d forcing (if applicable).
  if (has_iop_field("divT3d")) {
    if (m_iop_field_type.at("divT3d")==IOPFieldType::Computed) {
      const auto divT = get_iop_field("divT").get_view<const Real*>();
      const auto vertdivT = get_iop_field("vertdivT").get_view<const Real*>();
      const auto divT3d = get_iop_field("divT3d").get_view<Real*>();
      const auto nlevs = get_iop_field("divT3d").get_header().get_identifier().get_layout().dim(0);
      Kokkos::parallel_for(nlevs, KOKKOS_LAMBDA (const int ilev) {
        divT3d(ilev) = divT(ilev) + vertdivT(ilev);
      });
    }
  }
  if (has_iop_field("divq3d")) {
    if (m_iop_field_type.at("divq3d")==IOPFieldType::Computed) {
      const auto divq = get_iop_field("divq").get_view<const Real*>();
      const auto vertdivq = get_iop_field("vertdivq").get_view<const Real*>();
      const auto divq3d = get_iop_field("divq3d").get_view<Real*>();
      const auto nlevs = get_iop_field("divq3d").get_header().get_identifier().get_layout().dim(0);
      Kokkos::parallel_for(nlevs, KOKKOS_LAMBDA (const int ilev) {
        divq3d(ilev) = divq(ilev) + vertdivq(ilev);
      });
    }
  }

  // Now that data is loaded, reset the index of the currently loaded data.
  m_time_info.time_idx_of_current_data = iop_file_time_idx;
}

void IntensiveObservationPeriod::
set_fields_from_iop_data(const field_mgr_ptr field_mgr)
{
  if (m_params.get<bool>("zero_non_iop_tracers") && field_mgr->has_group("tracers")) {
    // Zero out all tracers before setting iop tracers (if requested)
    field_mgr->get_field_group("tracers").m_bundle->deep_copy(0);
  }

  EKAT_REQUIRE_MSG(field_mgr->get_grid()->name() == "Physics GLL",
                   "Error! Attempting to set non-GLL fields using "
                   "data from the IOP file.\n");

  // Find which fields need to be written
  const bool set_ps            = field_mgr->has_field("ps") && has_iop_field("Ps");
  const bool set_T_mid         = field_mgr->has_field("T_mid") && has_iop_field("T");
  const bool set_horiz_winds_u = field_mgr->has_field("horiz_winds") && has_iop_field("u");
  const bool set_horiz_winds_v = field_mgr->has_field("horiz_winds") && has_iop_field("v");
  const bool set_qv            = field_mgr->has_field("qv") && has_iop_field("q");
  const bool set_nc            = field_mgr->has_field("nc") && has_iop_field("NUMLIQ");
  const bool set_qc            = field_mgr->has_field("qc") && has_iop_field("CLDLIQ");
  const bool set_qi            = field_mgr->has_field("qi") && has_iop_field("CLDICE");
  const bool set_ni            = field_mgr->has_field("ni") && has_iop_field("NUMICE");

  // Create views/scalars for these field's data
  view_1d<Real> ps;
  view_2d<Real> T_mid, qv, nc, qc, qi, ni;
  view_3d<Real> horiz_winds;

  Real ps_iop;
  view_1d<Real> t_iop, u_iop, v_iop, qv_iop, nc_iop, qc_iop, qi_iop, ni_iop;

  if (set_ps) {
    ps = field_mgr->get_field("ps").get_view<Real*>();
    get_iop_field("Ps").sync_to_host();
    ps_iop = get_iop_field("Ps").get_view<Real, Host>()();
  }
  if (set_T_mid) {
    T_mid = field_mgr->get_field("T_mid").get_view<Real**>();
    t_iop = get_iop_field("T").get_view<Real*>();
  }
  if (set_horiz_winds_u || set_horiz_winds_v) {
    horiz_winds = field_mgr->get_field("horiz_winds").get_view<Real***>();
    if (set_horiz_winds_u) u_iop = get_iop_field("u").get_view<Real*>();
    if (set_horiz_winds_v) v_iop = get_iop_field("v").get_view<Real*>();
  }
  if (set_qv) {
    qv = field_mgr->get_field("qv").get_view<Real**>();
    qv_iop = get_iop_field("q").get_view<Real*>();
  }
  if (set_nc) {
    nc = field_mgr->get_field("nc").get_view<Real**>();
    nc_iop = get_iop_field("NUMLIQ").get_view<Real*>();
  }
  if (set_qc) {
    qc = field_mgr->get_field("qc").get_view<Real**>();
    qc_iop = get_iop_field("CLDLIQ").get_view<Real*>();
  }
  if (set_qi) {
    qi = field_mgr->get_field("qi").get_view<Real**>();
    qi_iop = get_iop_field("CLDICE").get_view<Real*>();
  }
  if (set_ni) {
    ni = field_mgr->get_field("ni").get_view<Real**>();
    ni_iop = get_iop_field("NUMICE").get_view<Real*>();
  }

  // Check if t_iop has any 0 entires near the top of the model
  // and correct t_iop and q_iop accordingly.
  correct_temperature_and_water_vapor(field_mgr);

  // Loop over all columns and copy IOP field values to FM views
  const auto ncols = field_mgr->get_grid()->get_num_local_dofs();
  const auto nlevs = field_mgr->get_grid()->get_num_vertical_levels();
  const auto policy = ESU::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const auto icol = team.league_rank();

    if (set_ps) {
      ps(icol) = ps_iop;
    }
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const int ilev) {
      if (set_T_mid) {
        T_mid(icol, ilev) = t_iop(ilev);
      }
      if (set_horiz_winds_u) {
        horiz_winds(icol, 0, ilev) = u_iop(ilev);
      }
      if (set_horiz_winds_v) {
        horiz_winds(icol, 1, ilev) = v_iop(ilev);
      }
      if (set_qv) {
        qv(icol, ilev) = qv_iop(ilev);
      }
      if (set_nc) {
        nc(icol, ilev) = nc_iop(ilev);
      }
      if (set_qc) {
        qc(icol, ilev) = qc_iop(ilev);
      }
      if (set_qi) {
        qi(icol, ilev) = qi_iop(ilev);
      }
      if (set_ni) {
        ni(icol, ilev) = ni_iop(ilev);
      }
    });
  });
}

void IntensiveObservationPeriod::
correct_temperature_and_water_vapor(const field_mgr_ptr field_mgr)
{
  // Find the first valid level index for t_iop, i.e., first non-zero entry
  int first_valid_idx;
  const auto nlevs = field_mgr->get_grid()->get_num_vertical_levels();
  auto t_iop = get_iop_field("T").get_view<Real*>();
  Kokkos::parallel_reduce(nlevs, KOKKOS_LAMBDA (const int ilev, int& lmin) {
    if (t_iop(ilev) > 0 && ilev < lmin) lmin = ilev;
  }, Kokkos::Min<int>(first_valid_idx));

  // If first_valid_idx>0, we must correct IOP fields T and q corresponding to
  // levels 0,...,first_valid_idx-1
  if (first_valid_idx > 0) {
    // If we have values of T and q to correct, we must have both T_mid and qv as FM fields
    EKAT_REQUIRE_MSG(field_mgr->has_field("T_mid"), "Error! IOP requires FM to define T_mid.\n");
    EKAT_REQUIRE_MSG(field_mgr->has_field("qv"),    "Error! IOP requires FM to define qv.\n");

    // Replace values of T and q where t_iop contains zeros
    auto T_mid = field_mgr->get_field("T_mid").get_view<const Real**>();
    auto qv   = field_mgr->get_field("qv").get_view<const Real**>();
    auto q_iop = get_iop_field("q").get_view<Real*>();
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0, first_valid_idx), KOKKOS_LAMBDA (const int ilev) {
      t_iop(ilev) = T_mid(0, ilev);
      q_iop(ilev) = qv(0, ilev);
    });
  }
}

} // namespace control
} // namespace scream
