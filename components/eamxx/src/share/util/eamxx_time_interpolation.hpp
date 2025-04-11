#ifndef EAMXX_TIME_INTERPOLATION_HPP
#define EAMXX_TIME_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"

#include "share/util/eamxx_time_stamp.hpp"

#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/io/scorpio_input.hpp"

namespace scream{
namespace util {

class TimeInterpolation {
public:
   using grid_ptr_type = std::shared_ptr<const AbstractGrid>;
   using vos_type = std::vector<std::string>;
   using fm_type = std::shared_ptr<FieldManager>;

  // Constructors & Destructor
  TimeInterpolation() = default;
  TimeInterpolation(const grid_ptr_type& grid);
  TimeInterpolation(const grid_ptr_type& grid, const vos_type& list_of_files);
  ~TimeInterpolation () = default;

  // Running the interpolation
  void initialize_timestamps(const TimeStamp& ts_in);
  void initialize_data_from_field(const Field& field_in);
  void initialize_data_from_files();
  void update_data_from_field(const Field& field_in);
  void update_timestamp(const TimeStamp& ts_in);
  void perform_time_interpolation(const TimeStamp& time_in);
  void finalize();

  // Build interpolator
  void add_field(const Field& field_in, const bool store_shallow_copy=false);

  // Getters
  Field get_field(const std::string& name) {
    return m_interp_fields.at(name);
  };

  // Informational
  void print();

  // Option to add a logger
  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& logger,
                  const std::string& header) {
      m_logger = logger;
      m_header = header;
  }

protected:

  // Internal structure to store data source triplets (when using data from file)
  // For each timesnap of data we have access to this triplet stores the
  //  - filename
  //  - timestamp
  //  - time index in the file.
  // Note, in many cases we will have files with multiple snaps of data and we
  // need a good way to organize this information.
  struct DataFromFileTriplet {
  public:
    std::string filename;
    TimeStamp   timestamp;
    int         time_idx;
  };

  // Helper functions to shift data
  void shift_data();

  // For the case where forcing data comes from files
  void set_file_data_triplets(const vos_type& list_of_files);
  void read_data();
  void check_and_update_data(const TimeStamp& ts_in);

  // Local field managers used to store two time snaps of data for interpolation
  fm_type  m_fm_time0;
  fm_type  m_fm_time1;
  vos_type m_field_names;
  std::map<std::string,Field> m_interp_fields;

  // Store the timestamps associated with the two time snaps
  TimeStamp m_time0;
  TimeStamp m_time1;

  // Variables related to the case where we use data from file
  std::vector<DataFromFileTriplet>           m_file_data_triplets;
  int                                        m_triplet_idx;
  std::shared_ptr<AtmosphereInput>           m_file_data_atm_input;
  bool                                       m_is_data_from_file=false;

  std::shared_ptr<ekat::logger::LoggerBase>  m_logger;
  std::string                                m_header;
}; // class TimeInterpolation

} // namespace util
} // namespace scream

#endif // EAMXX_TIME_INTERPOLATION_HPP
