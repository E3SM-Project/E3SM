#ifndef SCREAM_IOP_HPP
#define SCREAM_IOP_HPP

#include "share/scream_types.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

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

  using KT = ekat::KokkosTypes<DefaultDevice>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using Pack = ekat::Pack<Real, SCREAM_PACK_SIZE>;
  using Pack1d = ekat::Pack<Real, 1>;

  template<typename ScalarT>
  using view_1d = KT::template view_1d<ScalarT>;
  template<typename ScalarT>
  using view_2d = KT::template view_2d<ScalarT>;
  template<typename ScalarT>
  using view_3d = KT::template view_3d<ScalarT>;
  template<typename ScalarT>
  using view_1d_host = typename view_1d<ScalarT>::HostMirror;

public:

  // Constructor
  // Input:
  //   - comm: MPI communicator
  //   - params: Input yaml file needs intensive_observation_period_options sublist
  //   - run_t0: Initial timestamp for the simulation
  //   - model_nlevs: Number of vertical levels in the simulation. Needed since
  //                  the iop file contains a (potentially) different number of levels
  IntensiveObservationPeriod(const ekat::Comm& comm,
                             const ekat::ParameterList& params,
                             const util::TimeStamp& run_t0,
                             const int model_nlevs,
                             const Field& hyam,
                             const Field& hybm);

  // Default destructor
  ~IntensiveObservationPeriod() = default;

  // Read data from IOP file and store internally.
  void read_iop_file_data(const util::TimeStamp& current_ts);

  // Setup io grids for reading data from file and determine the closest lat/lon
  // pair in a IC/topo file to the target lat/lon params for a specific grid. This
  // should be called on each grid that loads field data from file before reading
  // data since the data file is not guarenteed to contain lat/lon for the correct
  // grid (e.g., loading PHIS_d from topography file which only contains lat/lon on
  // PG2 grid). EAMxx expects the ic file to contain lat/lon on GLL grid, and
  // topography file to contain lat/lon on PG2 grid.
  void setup_io_info (const std::string& file_name,
	               	    const grid_ptr& grid);

  // Read ICs and SPA data from file and remap to fields in field_mgr.
  // The remap is defined by setting all columns in the given fields to the
  // values of the column in the file with the closest lat,lon pair to
  // the target lat,lon in the parameters.
  // The function setup_io_info() must be called for the grids corresponding
  // to the file data before this function can be called.
  // Fields in the field_mgr must have the same number of levels as the file.
  // Inputs and outputs:
  //  - file_name: Name of the file used to load field data (IC or topo file)
  //  - field_names_nc: Field names used by the input file
  //  - field_names_eamxx: Field names used by eamxx
  //  - initial_ts: Inital timestamp. If initial_ts.is_valid()==false, then no timestamp is
  //                set for FM fields after storing data. Ex: passing util::TimeStamp().
  //  - field_mgr: Field manager containing fields that need data read from files
  //  - time_index: Time index of read. time_index=-1 will read the latest time in file.
  void read_fields_from_file_for_iop(const std::string& file_name,
                                     const vos& field_names_nc,
                                     const vos& field_names_eamxx,
                                     const util::TimeStamp& initial_ts,
                                     const field_mgr_ptr field_mgr,
                                     const int time_index = -1);

  // Version of above, but where nc and eamxx field names are identical
  void read_fields_from_file_for_iop(const std::string& file_name,
                                     const vos& field_names,
                                     const util::TimeStamp& initial_ts,
                                     const field_mgr_ptr field_mgr,
                                     const int time_index = -1)
  {
    read_fields_from_file_for_iop(file_name, field_names, field_names, initial_ts, field_mgr, time_index);
  }

  // Set fields using data loaded from the iop file
  void set_fields_from_iop_data(const field_mgr_ptr field_mgr);

  // The IOP file may contain temperature values that are
  // 0 at or above the surface. Correct these values using
  // the temperature T_mid (from field_mgr) where we
  // replace all values T_iop(k) == 0 with T_mid(0, k).
  // Likewise, at these k indices, we will replace q_iop(k)
  // with qv(0, k).
  // Note: We only need to use the first column because during
  //       the loading of ICs, every columns will have the same
  //       data.
  void correct_temperature_and_water_vapor(const field_mgr_ptr field_mgr);

  // Store grid spacing for use in SHOC ad interface
  void set_grid_spacing (const Real dx_short) {
    m_dynamics_dx_size = dx_short*1000;
  }

  Real get_dynamics_dx_size () { return m_dynamics_dx_size; }

  ekat::ParameterList& get_params() { return m_params; }

  bool has_iop_field(const std::string& fname) {
    return m_iop_fields.count(fname) > 0;
  }

  Field get_iop_field(const std::string& fname) {
    EKAT_REQUIRE_MSG(has_iop_field(fname), "Error! Requesting IOP field \""+fname+"\", but field is not stored in object.\n");
    return m_iop_fields[fname];
  }

private:

  // Struct for storing info related
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

  // Struct for storing relevant time information
  struct TimeInfo {
    util::TimeStamp iop_file_begin_time;
    view_1d_host<int> iop_file_times_in_sec;

    int time_idx_of_current_data = -1;

    int get_iop_file_time_idx (const util::TimeStamp& current_ts)
    {
      // Get iop file time index that the given timestamp falls between.
      // Note: the last time in iop file represents the non-inclusive
      //       upper bound of acceptable model times.
      const auto n_iop_times = iop_file_times_in_sec.extent(0);
      int time_idx=-1;
      for (size_t t=0; t<n_iop_times-1; ++t) {
        if (iop_file_begin_time + iop_file_times_in_sec(t) <= current_ts &&
            current_ts < iop_file_begin_time + iop_file_times_in_sec(t+1)) {
          time_idx = t;
        }
      }

      EKAT_REQUIRE_MSG(time_idx>=0,
                       "Error! Current model time ("+current_ts.to_string()+") is not within "
                       "IOP time period: ["+iop_file_begin_time.to_string()+", "+
                       (iop_file_begin_time+iop_file_times_in_sec(n_iop_times-1)).to_string()+").\n");
      return time_idx;
    }
  };

  void initialize_iop_file(const util::TimeStamp& run_t0,
                           int model_nlevs);

  ekat::Comm m_comm;
  ekat::ParameterList m_params;

  std::map<std::string,ClosestLatLonInfo> m_lat_lon_info;
  TimeInfo m_time_info;

  Real m_dynamics_dx_size;

  std::map<std::string,grid_ptr> m_io_grids;

  std::map<std::string, Field> m_iop_fields;
  std::map<std::string, Field> m_helper_fields;

  std::map<std::string, std::string> m_iop_file_varnames;
  std::map<std::string, std::string> m_iop_field_surface_varnames;
}; // class IntensiveObservationPeriod

} // namespace control
} // namespace scream

#endif // #ifndef SCREAM_IOP_HPP

