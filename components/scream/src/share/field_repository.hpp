#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "scream_types.hpp"
#include "field.hpp"

#include <map>

// I have to decide where to store a field's providers and customers, that is,
// a list of parametrizations that compute or need a given field. Should it
// be done by the FieldRepository? Or should a Field store its own? I lean toward
// the latter.
// NOTE: if two parametrizations both compute/update a field, we need to ensure 
//       that the two parametrizations are run sequentially. We need to devise
//       a mechanism to ensure that, or at the very least, a test to check it.
//       Notice that some fields are needed/computed/updated only within Scream, so
//       we can store a pointer to the Parametrization, but other are computed/needed
//       in other components, so for those we can at best store a string that
//       identifies the component (at least informally).

namespace scream
{

// We template a field (and a field manager) over the view type. The primary goal is
// to have two managers: one for the device views, which are be Managed, and one
// for the host views, which are Unmanaged, and simply store a pointer to the array
// that is managed by Fortran. The latter, is basically a manager for storing all
// the fields that are passed to the atmosphere by the coupler, or that the atmosphere
// needs to pass back to the coupler.
// NOTE: we will NOT pass fields back and forth at each iteration. We will instead set
//       all the pointers in the manager sometimes during the atm_init_mct call.
template<typename ViewType>
class FieldRepository {
public:

  // Public types
  using field_type  = Field<ViewType>;
  using header_type = typename field_type::header_type;
  using view_type   = typename field_type::view_type;

  // Do we need this? What structures need to be initialized?
  void initialize () { /* impl */ }

  // Closes the field registration phase
  // Up for debate: we should not allow calls to get_field*** before registration is completed.
  void registration_complete() { /* impl */ }

  // Cleans up the fm. This is needed since this class will most likely be contained inside
  // some singleton with static storage, which will be destroyed only after exit from main.
  // However, Kokkos prohibits to keep view objects alive after the call to Kokkos::finalize(),
  // which will be right before returning from main.
  void clean_up () { /* impl */ }

  // Methods to add new fields to the manager
  void register_field (const header_type& header) { create_field(header); }
  void register_field (const header_type& header, ViewType view) { create_field(header); *m_fields.at[header.m_name].m_field = view; }

  // Methods to query the database
  const field_type& get_field        (const std::string& name) const { return m_fields.at(name);          }
  const header_type& get_field_header (const std::string& name) const { return get_field(name).m_header;   }

  // Do we need/want this? It is only a shortcut to *get_field(name).m_field
  view_type get_field_view (const std::string& name) const { return *(get_field(name).m_field); }
  
protected:

  void create_field (const header_type& header) { /* check field does not exist first */ m_fields[header.m_name].m_field = std::make_shared<field_type>(); }

  // A map name->field
  std::map<std::string,field_type>   m_fields;
  // You could do
  // std::map<FieldHeader,Field>  m_fields;
  // which would allow to have more fields with the same name
  // (they would be distinguished by their rank/layout)
};

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
