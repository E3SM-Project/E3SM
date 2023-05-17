#ifndef EAMXX_TIME_INTERPOLATION_HPP
#define EAMXX_TIME_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"

#include "share/util/scream_time_stamp.hpp"

#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

namespace scream{
namespace util {

class TimeInterpolation {
public:
   using grid_ptr_type = std::shared_ptr<const AbstractGrid>;
   using vos_type = std::vector<std::string>;
   using fm_type = std::shared_ptr<FieldManager>;

  // Constructors
  TimeInterpolation() = default;
  TimeInterpolation(const grid_ptr_type& grid);
//  TimeInterpolation(const grid_ptr_type& grid, const TimeStamp& ts, const vos_type& list_of_files);

  // Running the interpolation
  void initialize_timestamps(const TimeStamp& ts_in);
  void initialize_data_from_field(const Field& field_in);
  void update_data_from_field(const Field& field_in);
  void update_timestamp(const TimeStamp& ts_in);
  std::map<std::string,Field> perform_time_interpolation(const TimeStamp& time_in);

  // Build interpolator
  void add_field(const Field& field_in);

  // Informational
  void print();

protected:

  // Helper functions to shift data
  void shift_data();
  void shift_data(const std::string& name);

  // Local field managers used to store two time snaps of data for interpolation
  fm_type m_fm_time0;
  fm_type m_fm_time1;

  // Store the timestamps associated with the two time snaps
  TimeStamp m_time0;
  TimeStamp m_time1;

}; // class TimeInterpolation

} // namespace util
} // namespace scream

#endif // EAMXX_TIME_INTERPOLATION_HPP
