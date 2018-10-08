#include "field_repository.hpp"

namespace scream {

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::register_field (std::shared_ptr<header_type> header,
                                           const bool abort_if_existing) {
  error::runtime_check(!m_registration_completed,"Error! Registration of new fields no longer allowed.\n",-1);
  error::runtime_check(static_cast<bool>(header),"Error! Invalid header pointer.\n",-1);
  auto it = m_fields.find(header->identifier());
  if (it==m_fields.end()) {
    return create_field(header);
  } else if (abort_if_existing) {
    error::runtime_abort("Error! Field already registered in the database.\n", -1);
  }
  return it->second;
}

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::get_field (std::shared_ptr<header_type> header) const {
  error::runtime_check(static_cast<bool>(header),"Error! Invalid header pointer.\n",-1);
  auto it = m_fields.find(header->identifier());
  error::runtime_check(it!=m_fields.end(), "Error! Field not found.\n", -1);
  return it->second;
}

template<typename MemSpace>
typename FieldRepository<MemSpace>::field_type
FieldRepository<MemSpace>::create_field (std::shared_ptr<header_type> header) {
  error::runtime_check(static_cast<bool>(header),"Error! Invalid header pointer.\n",-1);
  auto it_bool = m_fields.emplace(header->identifier(),header);
  error::runtime_check(it_bool.second, "Error! Someting went wrong while creating the field with identifier '" + header->identifier() + "'\n.", -1);
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
