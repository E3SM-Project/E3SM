#ifndef SCREAM_IOP_HPP
#define SCREAM_IOP_HPP

#include "share/eamxx_types.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {
namespace control {
/*
 * Class which data for an intensive observation period (IOP).
 */
class IOPDataManager
{
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
  //   - params: Input yaml file needs iop_options sublist
  //   - run_t0: Initial timestamp for the simulation
  //   - model_nlevs: Number of vertical levels in the simulation. Needed since
  //                  the iop file contains a (potentially) different number of levels
  IOPDataManager(const ekat::Comm& comm,
                 const ekat::ParameterList& params,
                 const util::TimeStamp& run_t0,
                 const int model_nlevs,
                 const Field& hyam,
                 const Field& hybm);

  // Destructor
  ~IOPDataManager();

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


  void read_fields_from_file_for_iop(const std::string& file_name,
                                     const std::vector<Field>& fields,
                                     const std::shared_ptr<const AbstractGrid>& tgt_grid);

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

  ekat::ParameterList& get_params() { return m_params; }

  bool has_iop_field(const std::string& fname) {
    return m_iop_fields.count(fname) > 0;
  }

  Field get_iop_field(const std::string& fname) {
    EKAT_REQUIRE_MSG(has_iop_field(fname), "Error! Requesting IOP field \""+fname+"\", but field is not stored in object.\n");
    return m_iop_fields[fname];
  }

private:

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

  enum IOPFieldType {
    FromFile,
    Computed
  };

  void initialize_iop_file(const util::TimeStamp& run_t0,
                           int model_nlevs);

  ekat::Comm m_comm;
  ekat::ParameterList m_params;

  TimeInfo m_time_info;

  Real m_dynamics_dx_size;

  std::map<std::string,grid_ptr> m_io_grids;

  std::map<std::string, Field> m_iop_fields;
  std::map<std::string, Field> m_helper_fields;

  std::map<std::string, std::string> m_iop_file_varnames;
  std::map<std::string, std::string> m_iop_field_surface_varnames;
  std::map<std::string, IOPFieldType> m_iop_field_type;
}; // class IOPDataManager

} // namespace control
} // namespace scream

#endif // #ifndef SCREAM_IOP_HPP

