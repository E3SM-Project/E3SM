#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "share/scream_types.hpp"
#include "share/scream_assert.hpp"
#include "share/field/field.hpp"

#include <map>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the field value type and over
  *  the memory space of the views.
  */

enum class RepoState {
  CLEAN,
  OPEN,
  CLOSED
};

template<typename ScalarType, typename Device>
class FieldRepository {
public:

  // Public types
  using device_type     = Device;
  using scalar_type     = ScalarType;
  using field_type      = Field<scalar_type,device_type>;
  using header_type     = typename field_type::header_type;
  using identifier_type = typename header_type::identifier_type;
  using map_type        = std::map<identifier_type,field_type>;

  // Constructor(s)
  FieldRepository ();

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Change the state of the database
  void registration_begins ();
  void registration_ends ();
  void clean_up ();

  // Deduce the pack size from the scalar type (which must be of type Pack<ScalarType,N>, for some int N>0, or ScalarType)
  template<typename RequestedValueType = scalar_type>
  void register_field (const identifier_type& identifier);

  // Methods to query the database
  int size () const { return m_fields.size(); }
  bool has_field (const identifier_type& identifier) const { return m_fields.find(identifier)!=m_fields.end(); }
  field_type get_field (const identifier_type& identifier) const;
  RepoState repository_state () const { return m_state; }

protected:

  // The state of the repository
  RepoState     m_state;

  // The actual repo
  map_type  m_fields;
};

// ============================== IMPLEMENTATION ============================= //

template<typename ScalarType, typename Device>
FieldRepository<ScalarType,Device>::FieldRepository ()
 : m_state (RepoState::CLEAN)
{
  // Nothing to be done here
}

template<typename ScalarType, typename Device>
template<typename RequestedValueType>
void FieldRepository<ScalarType,Device>::register_field (const identifier_type& id) {
  // Check that ScalarOrPackType is indeed ScalarType or Pack<ScalarType,N>, for some N>0.
  static_assert(std::is_same<ScalarType,RequestedValueType>::value ||
                std::is_same<ScalarType,typename util::ScalarProperties<RequestedValueType>::scalar_type>::value,
                "Error! The template argument 'RequestedValueType' of this function must either match "
                "the template argument 'ScalarType' of this class or be a Pack type based on ScalarType.\n");

  // Sanity checks
  error::runtime_check(m_state==RepoState::OPEN,"Error! Registration of new fields no longer allowed.\n");

  // Try to create the field. Allow case where it is already existing.
  auto it_bool = m_fields.emplace(id,field_type(id));

  // Make sure the field can accommodate the requested value type
  it_bool.first->second.get_header().get_alloc_properties().template request_value_type_allocation<RequestedValueType>();
}

template<typename ScalarType, typename Device>
typename FieldRepository<ScalarType,Device>::field_type
FieldRepository<ScalarType,Device>::get_field (const identifier_type& id) const {
  error::runtime_check(m_state==RepoState::CLOSED,"Error! You are not allowed to grab fields from the repo until after the registration phase is completed.\n");

  auto it = m_fields.find(id);
  error::runtime_check(it!=m_fields.end(), "Error! Field not found.\n");
  return it->second;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::registration_begins () {
  // Update the state of the repo
  m_state = RepoState::OPEN;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::registration_ends () {
  // Proceed to allocate fields
  for (auto& it : m_fields) {
    it.second.allocate_view();
  }

  // Prohibit further registration of fields
  m_state = RepoState::CLOSED;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::clean_up() {
  m_fields.clear();
  m_state = RepoState::CLEAN;
}


} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
