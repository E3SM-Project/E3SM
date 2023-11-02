#include "control/intensive_observation_period.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include "ekat/ekat_assert.hpp"

// Extend ekat mpi type for pair of double int, used
// in find_closest_lat_lon_index_and_rank()
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
	    	           const ekat::ParameterList& params)
{
  m_comm = comm;
  m_params = params;
  EKAT_REQUIRE_MSG(m_params.get<bool>("doubly_periodic_mode", false),
                   "Error! Currently doubly_periodic_mode is the only use case for "
	           "intensive observation period files.\n");

  EKAT_REQUIRE_MSG(m_params.isParameter("target_latitude") && m_params.isParameter("target_longitude"),
                   "Error! Using intensive observation period files requires "
                   "target_latitude and target_longitude be gives as parameters in "
                   "\"intensive_observation_period_options\" in the input yaml file.\n");
  const auto target_lat = m_params.get<Real>("target_latitude");
  const auto target_lon = m_params.get<Real>("target_longitude");
  EKAT_REQUIRE_MSG(-90 <= target_lat and target_lat <= 90,
                   "Error! IOP target_lat="+std::to_string(target_lat)+" outside of expected range [-90, 90].\n");
  EKAT_REQUIRE_MSG(0 <= target_lon and target_lon <= 360,
	           "Error! IOP target_lat="+std::to_string(target_lon)+" outside of expected range [0, 360].\n");
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
                               const field_mgr_ptr field_mgr)
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
    auto dims = fm_fid.get_layout().dims();
    dims[0] = io_grid->get_num_local_dofs();
    FieldLayout io_fl(fm_fid.get_layout().tags(), dims);
    FieldIdentifier io_fid(fm_fid.name(), io_fl, fm_fid.get_units(), io_grid->name());
    Field io_field(io_fid);
    io_field.allocate_view();
    io_fields.push_back(io_field);
  }

  // Read data from file
  AtmosphereInput file_reader(file_name,io_grid,io_fields);
  file_reader.read_variables();
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
    FieldLayout col_data_fl = io_fid.get_layout().strip_dim(0);
    FieldIdentifier col_data_fid("col_data", col_data_fl, dummy_units, "");
    Field col_data(col_data_fid);
    col_data.allocate_view();

    // MPI rank with closest column index store column data
    const auto mpi_rank_with_col = m_lat_lon_info[grid_name].mpi_rank_of_closest_column;
    if (m_comm.rank() == mpi_rank_with_col) {
      const auto col_indx_with_data = m_lat_lon_info[grid_name].local_column_index_of_closest_column;
      col_data.deep_copy(io_field.subfield(0,col_indx_with_data));
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

} // namespace control
} // namespace scream


