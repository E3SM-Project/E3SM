#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include <share/scream_types.hpp>
#include <share/error_defs.hpp>
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
  using field_type  = Field<Real*,MemSpace,MemoryManaged>;
  using header_type = typename field_type::header_type;
  using identifier_type = typename header_type::identifier_type;
  using view_type   = typename field_type::view_type;

  // Constructor(s)
  FieldRepository () : m_registration_completed(false) {}

  // No copies
  FieldRepository (const FieldRepository<MemSpace>&) = delete;
  FieldRepository& operator= (const FieldRepository<MemSpace>&) = delete;

  // Methods to add fields to the database
  field_type register_field (const identifier_type& identifier, const bool abort_if_existing = false);

  // Methods to query the database
  field_type get_field (const identifier_type& identifier) const;

  // Closes the field registration phase.
  // Using this checkpoint, allows to confine all the fields registration in one confined phase,
  // allowing for better debugging.
  void registration_complete ();

  // Queries whether the database is open for registration
  bool is_registration_open () const { return !m_registration_completed; }
  
  // Cleans up the repo. This is needed since this class will most likely be contained inside
  // some singleton with static storage, which will be destroyed only after exit from main.
  // However, Kokkos prohibits to keep view objects alive after the call to Kokkos::finalize(),
  // which will be right before returning from main.
  void clean_up ();

protected:

  // Create a new field
  field_type create_field (const identifier_type& identifier);

  // When true, no more fields are allowed to be added to the repo
  bool m_registration_completed;

  // A map identifier->field
  std::map<identifier_type,field_type>   m_fields;
};

// Explicit instantiation (declaration)
extern template class FieldRepository<ExecMemSpace>;
extern template class FieldRepository<HostMemSpace>;

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
