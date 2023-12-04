#ifndef SCREAM_IOP_HPP
#define SCREAM_IOP_HPP

#include "share/scream_types.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

namespace scream {
namespace control {
/*
 * Class which provides functionality for running EAMxx with an intensive
 * observation period (IOP). Currently the only use case is the doubly
 * periodic model (DP-SCREAM).
 */
class IntensiveObservationPeriod
{
  using vos = std::vector<std::string>;
  using field_mgr_ptr = std::shared_ptr<FieldManager>;
  using grid_ptr = std::shared_ptr<const AbstractGrid>;

public:

  // Constructor
  IntensiveObservationPeriod(const ekat::Comm& comm,
                             const ekat::ParameterList& params);

  // Default destructor
  ~IntensiveObservationPeriod() = default;

  // Setup io grids for reading data from file and determine the closest lat/lon
  // pair in a IC/topo file to the target lat/lon params for a specific grid. This
  // should be called on each grid that loads field data from file before reading
  // data since the data file is not guarenteed to contain lat/lon for the correct
  // grid (e.g., loading PHIS_d from topography file which only contains lat/lon on
  // PG2 grid). EAMxx expects the ic file to contain lat/lon on GLL grid, and
  // topography file to contain lat/lon on PG2 grid.
  void setup_io_info (const std::string& file_name,
	               	    const grid_ptr& grid);

  // Read ICs from file for IOP cases. We set all columns in the
  // given fields to the values of the column in the file with the
  // closest lat,lon pair to the target lat,lon in the parameters.
  // The setup_io_info must be called for the correct grids before
  // this function can be called.
  // Input:
  //  - file_name: Name of the file used to load field data (IC or topo file)
  //  - field_names_nc: Field names used by the input file
  //  - field_names_eamxx: Field names used by eamxx
  //  - initial_ts: Inital timestamp
  // Input/output
  //  - field_mgr: Field manager containing fields that need data read from files
  void read_fields_from_file_for_iop(const std::string& file_name,
                                     const vos& field_names_nc,
                                     const vos& field_names_eamxx,
                                     const util::TimeStamp& initial_ts,
                                     const field_mgr_ptr field_mgr);

  // Version of above, but where nc and eamxx field names are identical
  void read_fields_from_file_for_iop(const std::string& file_name,
				     const vos& field_names,
				     const util::TimeStamp& initial_ts,
				     const field_mgr_ptr field_mgr)
  {
    read_fields_from_file_for_iop(file_name, field_names, field_names, initial_ts, field_mgr);
  }

private:

  // Helper struct for storing info related
  // to the closest lat,lon pair
  struct ClosestLatLonInfo {
    // Value for the closest lat/lon in file.
    Real closest_lat;
    Real closest_lon;
    // MPI rank which owns the columns whose
    // lat,lon pair is closest to target lat,
    // lon parameters.
    int mpi_rank_of_closest_column;
    // Local column index of closest lat,lon pair.
    // Should be set -1 on ranks not equal to the
    // one above.
    int local_column_index_of_closest_column;
  };

  ekat::Comm m_comm;
  ekat::ParameterList m_params;
  std::map<std::string,ClosestLatLonInfo> m_lat_lon_info;
  std::map<std::string,grid_ptr> m_io_grids;
}; // class IntensiveObservationPeriod

} // namespace control
} // namespace scream

#endif // #ifndef SCREAM_IOP_HPP

