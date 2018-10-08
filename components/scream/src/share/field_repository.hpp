#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "scream_types.hpp"
#include "error_defs.hpp"
#include "field.hpp"

#include <map>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the memory space. The data type of the underlying view
  *  will *always* be Real*, so that we can store all the fields in one place (regardless of
  *  const/nonconst, and regardless of layout). The Field class offers some capabilities for
  *  reshaping a field (see field.hpp for details).
  */

template<typename MemSpace>
class FieldRepository {
public:

  // Public types
  using field_type  = Field<Real*,MemSpace,true>;
  using header_type = typename field_type::header_type;
  using view_type   = typename field_type::view_type;

  FieldRepository () : m_registration_completed(false) {}

  // Methods to add fields to the database
  field_type register_field (std::shared_ptr<header_type> header, const bool abort_if_existing = false);

  // Methods to query the database
  field_type get_field (std::shared_ptr<header_type> header) const;

  // Closes the field registration phase
  // Using this checkpoint, allows to confine all the fields registration in one execution phase,
  // allowing for better debugging.
  // Up for debate: we should not allow calls to get_field*** before registration is completed.
  void registration_complete();
  
  // Cleans up the repo. This is needed since this class will most likely be contained inside
  // some singleton with static storage, which will be destroyed only after exit from main.
  // However, Kokkos prohibits to keep view objects alive after the call to Kokkos::finalize(),
  // which will be right before returning from main.
  void clean_up ();

protected:

  // Create a new field
  field_type create_field (std::shared_ptr<header_type> header);

  // When true, no more fields are allowed to be added to the repo
  bool m_registration_completed;

  // A map identifier->field
  std::map<std::string,field_type>   m_fields;
};

// Explicit instantiation
extern template class FieldRepository<ExecMemSpace>;
extern template class FieldRepository<HostMemSpace>;

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
