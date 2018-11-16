#include "field_repository.hpp"

namespace scream {

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::register_field (const identifier_type& id,
                                           const bool abort_if_existing) {
  error::runtime_check(!m_registration_completed,"Error! Registration of new fields no longer allowed.\n");
  auto it = m_fields.find(id);
  if (it==m_fields.end()) {
    return create_field(id);
  } else if (abort_if_existing) {
    error::runtime_abort("Error! Field already registered in the database.\n");
  }
  return it->second;
}

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::get_field (const identifier_type& id) const {
  auto it = m_fields.find(id);
  error::runtime_check(it!=m_fields.end(), "Error! Field not found.\n");
  return it->second;
}

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::create_field (const identifier_type& id) {
  error::runtime_check(id.dimensions_set(), "Error! You need to set dimensions in the FieldIdentifier *before* creating the field.\n");
  auto it_bool = m_fields.emplace(id,id);
  error::runtime_check(it_bool.second, "Error! Someting went wrong while creating the field with identifier '" + id.get_identifier() + "'\n.");
  return it_bool.first->second;
}

template<typename MemSpace>
void FieldRepository<MemSpace>::registration_complete () {
  m_registration_completed = true;
}

template<typename MemSpace>
void FieldRepository<MemSpace>::clean_up() {
  m_fields.clear();
  m_registration_completed = false;
}

} // namespace scream
