#include "share/util/eamxx_time_interpolation.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/eamxx_io_utils.hpp"

namespace scream{
namespace util {

/*-----------------------------------------------------------------------------------------------*/
// Constructors
TimeInterpolation::TimeInterpolation(
  const grid_ptr_type& grid 
)
{
  // Given the grid initialize field managers to store interpolation data
  m_fm_time0 = std::make_shared<FieldManager>(grid);
  m_fm_time1 = std::make_shared<FieldManager>(grid);
  m_fm_time0->registration_begins();
  m_fm_time0->registration_ends();
  m_fm_time1->registration_begins();
  m_fm_time1->registration_ends();
}
/*-----------------------------------------------------------------------------------------------*/
TimeInterpolation::TimeInterpolation(
  const grid_ptr_type& grid, 
  const vos_type& list_of_files
) : TimeInterpolation(grid)
{
  set_file_data_triplets(list_of_files);
  m_is_data_from_file = true;
}
/*-----------------------------------------------------------------------------------------------*/
void TimeInterpolation::finalize()
{
  if (m_is_data_from_file) {
    m_file_data_atm_input = nullptr;
    m_is_data_from_file = false;
  }
}
/*-----------------------------------------------------------------------------------------------*/
/* A function to perform time interpolation using data from all the fields stored in the local
 * field managers.
 * Conducts a simple linear interpolation between two points using
 *        y* = w*y0 + (1-w)*y1
 * where w = (t1-t*)/(t1-t0)
 * Input:
 *   time_in - A timestamp to interpolate onto.
 * Output
 *   A map of (field name) and (Field) pairs such that each field is the interpolated values and
 *   the map makes it possible to quickly query individual outputs.
 */
void TimeInterpolation::perform_time_interpolation(const TimeStamp& time_in)
{
  // Declare the output map
  std::map<std::string,Field> interpolated_fields;

  // If data is handled by files we need to check that the timestamps are still relevant
  if (m_file_data_triplets.size()>0) {
    check_and_update_data(time_in);
  }

  // Gather weights for interpolation.  Note, timestamp differences are integers and we need a
  // real defined weight.
  const Real w_num = m_time1 - time_in;
  const Real w_den = m_time1 - m_time0;
  const Real weight0 = w_num/w_den;
  const Real weight1 = 1.0-weight0;

  // Cycle through all stored fields and conduct the time interpolation
  for (auto name : m_field_names)
  {
    const auto& field0   = m_fm_time0->get_field(name);
    const auto& field1   = m_fm_time1->get_field(name);
          auto field_out = m_interp_fields.at(name);
    field_out.deep_copy(field0);
    field_out.update(field1,weight1,weight0);
  }
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which registers a field in the local field managers.
 * Input:
 *   field_in - Is a field with the appropriate dimensions and metadata to match the interpolation
 * Output:
 *   None
 */
void TimeInterpolation::add_field(const Field& field_in, const bool store_shallow_copy)
{
  // First check that we haven't already added a field with the same name.
  const std::string name = field_in.name();
  EKAT_REQUIRE_MSG(!m_fm_time0->has_field(name) and !m_fm_time1->has_field(name),
		  "Error!! TimeInterpolation:add_field, field + " << name << " has already been added." << "\n");
  EKAT_REQUIRE_MSG (field_in.data_type()==DataType::FloatType or field_in.data_type()==DataType::DoubleType,
      "[TimeInterpolation] Error! Input field must have floating-point data type.\n"
      " - field name: " + field_in.name() + "\n"
      " - data type : " + e2str(field_in.data_type()) + "\n");

  // Clone the field for each field manager to get all the metadata correct.
  auto field0 = field_in.clone();
  auto field1 = field_in.clone();
  m_fm_time0->add_field(field0);
  m_fm_time1->add_field(field1);
  if (store_shallow_copy) {
    // Then we want to store the actual field_in and override it when interpolating
    m_interp_fields.emplace(name,field_in);
  } else {
    // We want to store a copy of the field but not ovveride
    auto field_out = field_in.clone();
    m_interp_fields.emplace(name,field_out);
  }
  m_field_names.push_back(name);
}
/*-----------------------------------------------------------------------------------------------*/
/* Function to shift all data from time1 to time0, update timestamp for time0
 */
void TimeInterpolation::shift_data()
{
  for (auto name : m_field_names)
  {
    auto& field0 = m_fm_time0->get_field(name);
    auto& field1 = m_fm_time1->get_field(name);
    std::swap(field0,field1);
  }
  m_file_data_atm_input->set_field_manager(m_fm_time1);
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which will initialize the TimeStamps.
 * Input:
 *   ts_in - A timestamp to set both time0 and time1 to. 
 *
 * At initialization we assume that only the first timestep of data has been set.  Subsequent
 * timesteps of data are added using update_data and then a call to update_timestamp will 
 * shift time0 to time1 and update time1.
 */
void TimeInterpolation::initialize_timestamps(const TimeStamp& ts_in)
{
  m_time0 = ts_in;
  m_time1 = ts_in;
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which will initialize field data given an input field.
 * Input:
 *   field_in - A field with a name matching one of the fields in the interpolator.  Data will be
 *              shifted and copied.
 */
void TimeInterpolation::initialize_data_from_field(const Field& field_in)
{
  const auto name = field_in.name();
  auto field0 = m_fm_time0->get_field(name);
  auto field1 = m_fm_time1->get_field(name);
  field0.deep_copy(field_in);
  field1.deep_copy(field_in);
  auto ts = field_in.get_header().get_tracking().get_time_stamp();
  m_time0 = ts;
  m_time1 = ts;
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which will initialize the full set of fields given a list of files pointing to where
 * the data can be found.
 * Input:
 *   list_of_files: A vector of strings representing all the files where interpolation data can be
 *                  gathered.
 */
void TimeInterpolation::initialize_data_from_files()
{
  auto triplet_curr = m_file_data_triplets[m_triplet_idx];
  // Initialize the AtmosphereInput object that will be used to gather data
  ekat::ParameterList input_params;
  input_params.set("Field Names",m_field_names);
  input_params.set("Filename",triplet_curr.filename);
  m_file_data_atm_input = std::make_shared<AtmosphereInput>(input_params,m_fm_time1);
  m_file_data_atm_input->set_logger(m_logger);
  // Assign the mask value gathered from the FillValue found in the source file.
  // TODO: Should we make it possible to check if FillValue is in the metadata and only assign mask_value if it is?
  for (auto& name : m_field_names) {
    auto& field0 = m_fm_time0->get_field(name);
    auto& field1 = m_fm_time1->get_field(name);
    auto& field_out = m_interp_fields.at(name);

    auto set_fill_value = [&](const auto var_fill_value) {
      const auto dt = field_out.data_type();
      if (dt==DataType::FloatType) {
        field0.get_header().set_extra_data("mask_value",static_cast<float>(var_fill_value));
        field1.get_header().set_extra_data("mask_value",static_cast<float>(var_fill_value));
        field_out.get_header().set_extra_data("mask_value",static_cast<float>(var_fill_value));
      } else if (dt==DataType::DoubleType) {
        field0.get_header().set_extra_data("mask_value",static_cast<double>(var_fill_value));
        field1.get_header().set_extra_data("mask_value",static_cast<double>(var_fill_value));
        field_out.get_header().set_extra_data("mask_value",static_cast<double>(var_fill_value));
      } else {
        EKAT_ERROR_MSG (
            "[TimeInterpolation] Unexpected/unsupported field data type.\n"
            " - field name: " + field_out.name() + "\n"
            " - data type : " + e2str(dt) + "\n");
      }
    };

    const auto& pio_var = scorpio::get_var(triplet_curr.filename,name);
    if (scorpio::refine_dtype(pio_var.nc_dtype)=="float") {
      auto var_fill_value = scorpio::get_attribute<float>(triplet_curr.filename,name,"_FillValue");
      set_fill_value(var_fill_value);
    } else if (scorpio::refine_dtype(pio_var.nc_dtype)=="double") {
      auto var_fill_value = scorpio::get_attribute<double>(triplet_curr.filename,name,"_FillValue");
      set_fill_value(var_fill_value);
    } else {
      EKAT_ERROR_MSG (
          "Unrecognized/unsupported data type\n"
          " - filename: " + triplet_curr.filename + "\n"
          " - varname : " + name + "\n"
          " - dtype   : " + pio_var.dtype + "\n");
    }

  }
  // Read first snap of data and shift to time0
  read_data();
  shift_data();
  update_timestamp(triplet_curr.timestamp);
  // Advance the iterator and read the next set of data for time1
  ++m_triplet_idx;
  read_data();
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which will update the timestamps by shifting time1 to time0 and setting time1.
 * Input:
 *   ts_in - A timestamp for the most recent timestamp of the interpolation data.
 *
 * It is assumed that any updated time is meant to replace time1 and that the time1
 * should be shifted to time0
 */
void TimeInterpolation::update_timestamp(const TimeStamp& ts_in)
{
  m_time0 = m_time1;
  m_time1 = ts_in;
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which will update field data given an input field.  Useful for time interpolation
 * related to the EAMxx state.
 * Input:
 *   field_in - A field with a name matching one of the fields in the interpolator.  Data will be
 *              shifted and copied.
 *
 * It is assumed that any updated data is meant to replace the data at time1 and that the time1
 * data should be shifted to time0
 */
void TimeInterpolation::update_data_from_field(const Field& field_in)
{
  const auto name = field_in.name();
  auto& field0 = m_fm_time0->get_field(name);
  auto& field1 = m_fm_time1->get_field(name);
  std::swap(field0,field1);
  // Now that we have swapped field0 and field 1 we need to grab field 1 from the field manager again.
  // Alternatively we could just update `field0` which is now inside m_fm_time1, but choosing this
  // approach for code readability.
  auto& field1_new = m_fm_time1->get_field(name);
  field1_new.deep_copy(field_in);
}
/*-----------------------------------------------------------------------------------------------*/
void TimeInterpolation::print()
{
  printf("Settings for time interpolator...\n");
  printf("Time 0 = %s\n",m_time0.to_string().c_str());
  printf("Time 1 = %s\n",m_time1.to_string().c_str());
  printf("List of Fields in interpolator:\n");
  for (auto name : m_field_names)
  {
    printf("     -   %16s\n",name.c_str());
  }

}
/*-----------------------------------------------------------------------------------------------*/
/* Function to organize the data available from data files by time.  This is necessary to quickly
 * grab the appropriate data each time interpolation is requested.
 * Input:
 *   list_of_files - Is a vector of strings representing all files that can be used for data.
 *
 * We create a reference timestamp to use when sorting the data snaps.
 */
void TimeInterpolation::set_file_data_triplets(const vos_type& list_of_files) {
  EKAT_REQUIRE_MSG(list_of_files.size()>0,"ERROR! TimeInterpolation::set_file_data_triplets - the list of files is empty. Please check.");
  TimeStamp ts_ref;
  // The first step is to grab all relevant metadata for the DataFromFileTriplet objects.
  // We will store the times in a map and take advantage of maps natural sorting to organize the triplets
  // in chronological order.  This ensures that if the list of files does not represent the chronological
  // order to the data we will still have sorted data.
  vos_type               filenames_tmp;
  std::vector<TimeStamp> timestamps_tmp;
  std::vector<int>       time_idx_tmp;
  std::map<int,int>      map_of_times_to_vector_idx;
  int running_idx = 0;
  for (size_t ii=0; ii<list_of_files.size(); ii++) {
    const auto filename = list_of_files[ii];
    // Reference TimeStamp
    scorpio::register_file(filename,scorpio::FileMode::Read);
    auto ts_file_start = read_timestamp(filename,"case_t0");
    // Gather the units of time
    auto time_units = scorpio::get_attribute<std::string>(filename,"time","units");
    int time_mult;
    if (time_units.find("seconds") != std::string::npos) {
      time_mult = 1;
    } else if (time_units.find("minutes") != std::string::npos) {
      time_mult = 60;
    } else if (time_units.find("hours") != std::string::npos) {
      time_mult = 3600;
    } else if (time_units.find("days") != std::string::npos) {
      time_mult = 86400;
    } else {
      EKAT_ERROR_MSG("Error!! TimeInterpolation::set_file_triplets - unsupported units of time = (" << time_units << ") in source data file " << filename << ", supported units are: seconds, minutes, hours and days");
    }
    // Gather information about time in this file
    if (ii==0) {
      ts_ref = ts_file_start;
    }
    const int ntime = scorpio::get_dimlen(filename,"time");
    for (int tt=0; tt<ntime; tt++) {
      auto time_snap = scorpio::get_time(filename,tt);
      TimeStamp ts_snap = ts_file_start;
      if (time_snap>0) {
        ts_snap += (time_snap*time_mult);
      }
      auto time = ts_snap.seconds_from(ts_ref);
      // Sanity check that we don't have multiples of the same timesnap
      EKAT_REQUIRE_MSG(map_of_times_to_vector_idx.count(time)==0,"Error! TimeInterpolation::set_file_data_triplets - The same time step has been encountered more than once in the data files, please check\n"  
		      << "    TimeStamp: " << ts_snap.to_string() << "\n"  
		      << "     Filename: " << filename << "\n");
      map_of_times_to_vector_idx.emplace(time,running_idx);
      filenames_tmp.push_back(filename);
      timestamps_tmp.push_back(ts_snap);
      time_idx_tmp.push_back(tt);
      ++running_idx;
    }
    scorpio::release_file(filename);
  }
  // Now that we have gathered all of the timesnaps we can arrange them in order as DataFromFileTriplet objects.
  // Taking advantage of maps automatically self-sorting by the first arg.
  for (auto a : map_of_times_to_vector_idx) {
    auto idx = a.second;
    DataFromFileTriplet my_trip;
    my_trip.filename  = filenames_tmp[idx];
    my_trip.timestamp = timestamps_tmp[idx];
    my_trip.time_idx  = time_idx_tmp[idx];
    m_file_data_triplets.push_back(my_trip);
  }
  // Finally set the iterator to point to the first triplet.
  m_triplet_idx = 0;
}	
/*-----------------------------------------------------------------------------------------------*/
/* Function to read a new set of data from file using the current iterator pointing to the current
 * DataFromFileTriplet.
 */
void TimeInterpolation::read_data()
{
  const auto triplet_curr = m_file_data_triplets[m_triplet_idx];
  if (not m_file_data_atm_input or triplet_curr.filename != m_file_data_atm_input->get_filename()) {
    // Then we need to close this input stream and open a new one
    ekat::ParameterList input_params;
    input_params.set("Field Names",m_field_names);
    input_params.set("Filename",triplet_curr.filename);
    m_file_data_atm_input = std::make_shared<AtmosphereInput>(input_params,m_fm_time1);
    m_file_data_atm_input->set_logger(m_logger);
    // Also determine the FillValue, if used
    // TODO: Should we make it possible to check if FillValue is in the metadata and only assign mask_value if it is?
    for (auto& name : m_field_names) {
      auto& field = m_fm_time1->get_field(name);
      const auto dt = field.data_type();
      if (dt==DataType::FloatType) {
        auto var_fill_value = scorpio::get_attribute<float>(triplet_curr.filename,name,"_FillValue");
        field.get_header().set_extra_data("mask_value",var_fill_value);
      } else if (dt==DataType::DoubleType) {
        auto var_fill_value = scorpio::get_attribute<double>(triplet_curr.filename,name,"_FillValue");
        field.get_header().set_extra_data("mask_value",var_fill_value);
      } else {
        EKAT_ERROR_MSG (
            "[TimeInterpolation] Unexpected/unsupported field data type.\n"
            " - field name: " + field.name() + "\n"
            " - data type : " + e2str(dt) + "\n");
      }
    }
  }

  if (m_logger) {
    m_logger->info(m_header);
    m_logger->info("[EAMxx:time_interpolation] Reading data at time " + triplet_curr.timestamp.to_string());
  }
  m_file_data_atm_input->read_variables(triplet_curr.time_idx);
  m_time1 = triplet_curr.timestamp;
}
/*-----------------------------------------------------------------------------------------------*/
/* Function to check the current set of interpolation data against a timestamp and, if needed,
 * update the set of interpolation data to ensure the passed timestamp is within the bounds of
 * the interpolation data.
 * Input:
 *   ts_in - A timestamp that we intend to interpolate onto.
 */
void TimeInterpolation::check_and_update_data(const TimeStamp& ts_in)
{
  // First check if the passed timestamp is within the bounds of time0 and time1.
  EKAT_REQUIRE_MSG(ts_in.seconds_from(m_time0) >= 0, "ERROR!!! TimeInterpolation::check_and_update_data - "
		  << "Current timestamp of " << ts_in.to_string() << " is lower than the TimeInterpolation bounds of " << m_time0.to_string());
  if (m_time1.seconds_from(ts_in) < 0) {
    // The timestamp is out of bounds, need to load new data.
    // First cycle through the DataFromFileTriplet's to find a timestamp that is greater than this one.
    bool found = false;
    int step_cnt = 0; // Track how many triplets we passed to find one that worked. 
    while (m_triplet_idx < static_cast<int>(m_file_data_triplets.size())) {
      ++m_triplet_idx;
      ++step_cnt;
      auto ts_tmp = m_file_data_triplets[m_triplet_idx].timestamp;
      if (ts_tmp.seconds_from(ts_in) >= 0) {
        // This timestamp is greater than the input timestamp, we can use it
	found = true;
	break;
      }
    }
    EKAT_REQUIRE_MSG(found,"ERROR!! TimeInterpolation::check_and_update_data - timestamp " << ts_in.to_string() << "is outside the bounds of the set of data files." << "\n"
		   <<  "     TimeStamp time0: " << m_time0.to_string() << "\n"
		   <<  "     TimeStamp time1: " << m_time1.to_string() << "\n");
    // Now we need to make sure we didn't jump more than one triplet, if we did then the data at time0 is
    // incorrect.
    if (step_cnt>1) {
      // Then we need to populate data for time1 as the previous triplet before shifting data to time0
      --m_triplet_idx;
      read_data();
      ++m_triplet_idx;
    }
    // We shift the time1 data to time0 and read the new data.
    shift_data();
    update_timestamp(m_file_data_triplets[m_triplet_idx].timestamp);
    read_data();
    // Sanity Check
    bool current_data_check = (ts_in.seconds_from(m_time0) >= 0) and (m_time1.seconds_from(ts_in) >= 0);
    EKAT_REQUIRE_MSG(current_data_check,"ERROR!! TimeInterpolation::check_and_update_data - Something went wrong in updating data:\n"
		   <<  "      TimeStamp    IN: " << ts_in.to_string() << "\n"
		   <<  "      TimeStamp time0: " << m_time0.to_string() << "\n"
		   <<  "      TimeStamp time1: " << m_time1.to_string() << "\n");

  }
}
/*-----------------------------------------------------------------------------------------------*/

} // namespace util
} // namespace scream
