#include "share/util/eamxx_time_interpolation.hpp"

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
/* A function to perform time interpolation using data from all the fields stored in the local
 * field managers.
 * Conducts a simple linear interpolation between two points using
 *        y* = w*y0 + (1-w)*y1
 * where w = (t1-dt)/(t1-t0)
 * Input:
 *   time_in - A timestamp to interpolate onto.
 * Output
 *   A map of (field name) and (Field) pairs such that each field is the interpolated values and
 *   the map makes it possible to quickly query individual outputs.
 */
std::map<std::string,Field> TimeInterpolation::perform_time_interpolation(const TimeStamp& time_in)
{
  // Declare the output map
  std::map<std::string,Field> interpolated_fields;

  // Gather weights for interpolation.  Note, timestamp differences are integers and we need a
  // real defined weight.
  const Real w_num = m_time1 - time_in;
  const Real w_den = m_time1 - m_time0;
  const Real weight = w_num/w_den;

  // Cycle through all stored fields and conduct the time interpolation
  for (auto ff = m_fm_time0->begin(); ff != m_fm_time0->end(); ff++)
  {
    const auto name    = ff->first;
    const auto& field0 = m_fm_time0->get_field(name);
    const auto& field1 = m_fm_time1->get_field(name);
    Field field_out    = field0.clone();
    field_out.update(field1,1.0-weight,weight);
    interpolated_fields.emplace(name,field_out);
  }

  // Return map of interpolated fields
  return interpolated_fields;
}
/*-----------------------------------------------------------------------------------------------*/
/* Function which registers a field in the local field managers.
 * Input:
 *   field_in - Is a field with the appropriate dimensions and metadata to match the interpolation
 * Output:
 *   None
 */
void TimeInterpolation::add_field(const Field& field_in)
{
  // First check that we haven't already added a field with the same name.
  const auto name = field_in.name();
  EKAT_REQUIRE_MSG(!m_fm_time0->has_field(name) and !m_fm_time1->has_field(name),
		  "Error!! TimeInterpolation:add_field, field + " << name << " has already been added." << "\n");

  // Clone the field for each field manager to get all the metadata correct.
  auto field0 = field_in.clone();
  auto field1 = field_in.clone();
  m_fm_time0->add_field(field0);
  m_fm_time1->add_field(field1);
}
/*-----------------------------------------------------------------------------------------------*/
/* Function to shift the data of a single field from time1 to time0
 * Input:
 *   name - The name of the field to be shifted.
 */
void TimeInterpolation::shift_data(const std::string& name)
{
        auto  field0 = m_fm_time0->get_field(name);
  const auto& field1 = m_fm_time1->get_field(name);
  field0.deep_copy(field1);
}
/*-----------------------------------------------------------------------------------------------*/
/* Function to shift all data from time1 to time0, update timestamp for time0
 */
void TimeInterpolation::shift_data()
{
  for (auto ff = m_fm_time0->begin(); ff != m_fm_time0->end(); ff++)
  {
    const auto name = ff->first;
    shift_data(name);
  }
  m_time0 = m_time1;
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
  shift_data(name);
  auto field1 = m_fm_time1->get_field(name);
  field1.deep_copy(field_in);
}
/*-----------------------------------------------------------------------------------------------*/

} // namespace util
} // namespace scream
