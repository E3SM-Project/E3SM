#include "control/intensive_observation_period.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include "ekat/ekat_assert.hpp"

#include <pio.h>

// Extend ekat mpi type for pair of double int, used
// in find_closest_lat_lon_index_and_rank()
namespace ekat {
  template<>
  MPI_Datatype get_mpi_type<std::pair<double, int>> () {
    return MPI_DOUBLE_INT;
  }
}

namespace scream {
namespace control {

// Helper functions for reading files
namespace {
// Currently EAMxx scorpio interface doesn't have a way to inquire about dimension names
bool file_has_dim(const std::string& filename, const std::string& dimname)
{
  int ncid, dimid, err;

  bool was_open = scorpio::is_file_open_c2f(filename.c_str(),-1);
  if (not was_open) {
    scorpio::register_file(filename,scorpio::FileMode::Read);
  }

  ncid = scorpio::get_file_ncid_c2f (filename.c_str());
  err = PIOc_inq_dimid(ncid,dimname.c_str(),&dimid);
  if (err==PIO_EBADDIM) return false;

  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while retrieving dimension id.\n"
      " - filename : " + filename + "\n"
      " - dimname  : " + dimname + "\n"
      " - pio error: " + std::to_string(err) + "\n");
  if (not was_open) {
    scorpio::eam_pio_closefile(filename);
  }
  return true;
}
}

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
                   "target_lat and target_lon be gives as parameters in "
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
    AtmosphereInput dp_ic_reader(file_name, io_grid, fields);
    dp_ic_reader.read_variables();
    dp_ic_reader.finalize();

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
  if (grid_name=="Physics GLL" && file_has_dim(file_name,"ncol_d")) {
    // If we are on GLL grid, and nc file contains "ncol_d" dimension,
    // we need to reset COL dim tag
    using namespace ShortFieldTagsNames;
    auto grid = io_grid->clone(io_grid->name(),true);
    grid->reset_field_tag_name(COL,"ncol_d");
    io_grid = grid;
  }

  // Create vector of fields with correct dimensions to read from file
  std::vector<Field> fields;
  for (size_t i=0; i<field_names_nc.size(); ++i) {
    const auto& nc_name    = field_names_nc[i];
    const auto& eamxx_name = field_names_eamxx[i];

    auto fm_field = eamxx_name!=nc_name
                    ?
                    field_mgr->get_field(eamxx_name).alias(nc_name)
                    :
                    field_mgr->get_field(eamxx_name);
    auto fm_fid = fm_field.get_header().get_identifier();
    auto dims = fm_fid.get_layout().dims();
    EKAT_REQUIRE_MSG(fm_fid.get_layout().tag(0)==FieldTag::Column,
                     "Error! IOP inputs read from IC/topo file must have Column "
                     "as first dim tag.\n");

    // Set first dimension to match input file
    dims[0] = io_grid->get_num_local_dofs();
    FieldLayout dp_fl(fm_fid.get_layout().tags(), dims);
    FieldIdentifier dp_fid(fm_fid.name(), dp_fl, fm_fid.get_units(), io_grid->name());
    Field dp_field(dp_fid);
    dp_field.allocate_view();
    fields.push_back(dp_field);
  }

  // Read data from file
  AtmosphereInput dp_ic_reader(file_name,io_grid,fields);
  dp_ic_reader.read_variables();
  dp_ic_reader.finalize();

  // For each field, broadcast data from closest lat/lon column to all processors
  // and copy data into each field's column.

  // Find the maximum size of transfer data over fields.
  int max_data_size = 0;
  for (auto& f : fields) {
    const auto fl = f.get_header().get_identifier().get_layout();
    if (fl.size()/fl.dim(0) > max_data_size) max_data_size = fl.size()/fl.dim(0);
  }

  Real transfer_data[max_data_size];
  for (size_t i=0; i<fields.size(); ++i) {
    const auto f = fields[i];
    const auto fl = f.get_header().get_identifier().get_layout();
    const auto field_rank = fl.rank();
    const auto fname = field_names_eamxx[i];

    // If I am the rank with the minimum lat/lon distance,
    // store data of values to broadcast.
    const auto my_rank = m_comm.rank();
    const auto col_mpi_rank = m_lat_lon_info[grid_name].mpi_rank_of_closest_column;
    if (my_rank == col_mpi_rank) {
      const auto col_idx = m_lat_lon_info[grid_name].local_column_index_of_closest_column;
      switch (field_rank) {
        case 1: {
          auto src_view = f.get_view<Real*, Host>();
          transfer_data[0] = src_view(col_idx);
          break;
        }
        case 2: {
          auto src_view = f.get_view<Real**, Host>();
          for (auto i1=0; i1<fl.dim(1); ++i1) {
            transfer_data[i1] = src_view(col_idx, i1);
          }
          break;
        }
        case 3: {
          auto src_view = f.get_view<Real***, Host>();
          for (auto i1=0; i1<fl.dim(1); ++i1) {
            for (auto i2=0; i2<fl.dim(2); ++i2) {
              transfer_data[i2+i1*fl.dim(2)] = src_view(col_idx, i1, i2);
          }}
          break;
        }
        default:
          // Currently only up to rank 3 inputs needed from IC file.
          // Additional ranks can easily be added if needed.
          EKAT_ERROR_MSG ("Error! Unexpected field rank (" + std::to_string(field_rank) + ") when transferring "
                          "IOP ICs for field "+fname+". If needed, this can be implemented by adding another case.\n");
      }
    }

    // Broadcast data to other processors
    const auto data_size_needed = fl.size()/fl.dim(0);
    m_comm.broadcast(transfer_data, data_size_needed, col_mpi_rank);

    // For each processor, copy data  into each column
    auto fm_field = field_mgr->get_field(fname);
    auto fm_field_rank = fm_field.get_header().get_identifier().get_layout().rank();
    // Sanity check
    EKAT_REQUIRE_MSG(field_rank==fm_field_rank,
                     "Error! Field ranks mismatch for "+fname+" when transfering IOP ICs.\n");
    switch (field_rank) {
      case 1: {
        auto dst_view = fm_field.get_view<Real*, Host>();
        for (size_t i0=0; i0<dst_view.extent(0); ++i0) {
          dst_view(i0) = transfer_data[0];
        }
        break;
      }
      case 2: {
        auto dst_view = fm_field.get_view<Real**, Host>();
        for (size_t i0=0; i0<dst_view.extent(0); ++i0) {
          for (size_t i1=0; i1<dst_view.extent(1); ++i1) {
            dst_view(i0, i1) = transfer_data[i1];
        }}
        break;
      }
      case 3: {
        auto dst_view = fm_field.get_view<Real***, Host>();
        for (size_t i0=0; i0<dst_view.extent(0); ++i0) {
          for (size_t i1=0; i1<dst_view.extent(1); ++i1) {
            for (size_t i2=0; i2<dst_view.extent(2); ++i2) {
              dst_view(i0, i1, i2) = transfer_data[i2+i1*dst_view.extent(2)];
        }}}
        break;
      }
      default:
        // Currently only up to rank 3 inputs needed from IC file.
        // Additional ranks can easily be added if needed.
        EKAT_ERROR_MSG ("Error! Unexpected field rank (" + std::to_string(field_rank) + ") when transferring "
                        "IOP ICs for field "+fname+". If needed, this can be implemented by adding another case.\n");
    }

    // Sync fields to device
    fm_field.sync_to_dev();
  }

  // Set the initial time stamp on FM fields
  for (const auto& fn : field_names_eamxx) {
    field_mgr->get_field(fn).get_header().get_tracking().update_time_stamp(initial_ts);
  }
}

} // namespace control
} // namespace scream


