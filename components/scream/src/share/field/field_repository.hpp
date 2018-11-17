#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include <share/scream_types.hpp>
#include <share/error_defs.hpp>
#include <share/field/field.hpp>

#include <map>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the field value type and over
  *  the memory space of the views. The data type of the underlying view
  *  will *always* be scalar_type*, so that we can store all the fields
  *  in one place (regardless of const/nonconst, and regardless of layout).
  *  You can find some utilities for reshaping a field in field_utils.hpp.
  */

template<typename Scalar, typename MemSpace>
class FieldRepository {
public:

  // Public types
  using scalar_type     = Scalar;
  using field_type      = Field<scalar_type*,MemSpace,MemoryManaged>;
  using header_type     = typename field_type::header_type;
  using identifier_type = typename header_type::identifier_type;
  using view_type       = typename field_type::view_type;

  // Constructor(s)
  FieldRepository () : m_registration_completed(false) {}

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Deduce the pack size from the scalar type (which must be of type Pack<ScalarType,N>, for some int N>0, or ScalarType)
  template<typename ScalarOrPackType = scalar_type>
  void register_field (const identifier_type& identifier, const bool abort_if_existing = false);

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

  // When true, no more fields are allowed to be added to the repo
  bool                                  m_registration_completed;

  // A map identifier->field
  std::map<identifier_type,field_type>  m_fields;
};

// ============================== IMPLEMENTATION ============================= //

template<typename Scalar, typename MemSpace>
template<typename ScalarOrPackType>
void FieldRepository<Scalar,MemSpace>::register_field (const identifier_type& id,
                                                       const bool abort_if_existing) {
  // Check that ScalarOrPackType is indeed Scalar or Pack<Scalar,N>, for some N>0.
  static_assert(std::is_same<Scalar,ScalarOrPackType>::value ||
                std::is_same<Scalar,typename util::is_pack<ScalarOrPackType>::scalar_type>::value,
                "Error! The template argument 'ScalarOrPackType' of this function must either match "
                "the template argument 'Scalar' of this class or be a Pack<Scalar,N>, for some N>0.\n");
  
  // Sanity checks
  error::runtime_check(!m_registration_completed,"Error! Registration of new fields no longer allowed.\n");

  // Try to create the field. Abort if we fail and abort_if_existing is true
  auto it_bool = m_fields.emplace(id,id);
  if (!it_bool.second && abort_if_existing) {
    error::runtime_abort("Error! Field with identifier '" + id.get_identifier() + "' already registered in the database.\n");
  }

  // Detect pack size
  int pack_size = 1;
  if (util::is_pack<ScalarOrPackType>::value) {
    pack_size = util::is_pack<ScalarOrPackType>::size;
  }

  // Set the requested pack size in the field
  it_bool.first->second.request_pack_size(pack_size);
}

template<typename Scalar, typename MemSpace>
typename FieldRepository<Scalar,MemSpace>::field_type
FieldRepository<Scalar,MemSpace>::get_field (const identifier_type& id) const {
  error::runtime_check(m_registration_completed,"Error! You are not allowed to grab fields from the repo until after the registration phase is completed.\n");

  auto it = m_fields.find(id);
  error::runtime_check(it!=m_fields.end(), "Error! Field not found.\n");
  return it->second;
}

template<typename Scalar, typename MemSpace>
void FieldRepository<Scalar,MemSpace>::registration_complete () {
  // Proceed to allocate fields
  for (auto& it : m_fields) {
    it.second.allocate_view();
  }

  // Prohibit further registration of fields
  m_registration_completed = true;
}

template<typename Scalar, typename MemSpace>
void FieldRepository<Scalar,MemSpace>::clean_up() {
  m_fields.clear();
  m_registration_completed = false;
}


} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
