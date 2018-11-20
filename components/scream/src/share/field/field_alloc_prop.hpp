#ifndef SCREAM_FIELD_ALLOC_PROP_HPP
#define SCREAM_FIELD_ALLOC_PROP_HPP

#include <share/error_defs.hpp>
#include <share/field/field_identifier.hpp>

#include <share/scream_types.hpp>

#include <vector>

namespace scream
{


/*
 *  Small structure holding a few properties of the field allocation
 *
 *  This class allows to keep track of a field allocation properties.
 *  The reason for this is that the allocation could be different from
 *  what one would infer by simply looking at the field identifier.
 *  For instance, say you need a 2d field with dimensions (3,10).
 *  In the field repo, this would be stored as a 1d field: View<Real*>.
 *  Let's say you have one class that uses the field as View<Real**>,
 *  with dimensions (3,10), and another class that used the field in
 *  a packed way, as View<Pack<Real,4>**>. The second view will have
 *  dimensions (3,4), in terms of its value type (i.e., Pack<Real,4>).
 *  This means the field needs an allocation bigger than 30 Real's, namely
 *  36 Real's. This class keeps track of all the requested sizes, so
 *  that 1) the field is allocated with enough memory to accommodate all
 *  uses, and 2) customers of the field can check what the allocation
 *  was, so that they know whether there is padding in the field.
 */
class FieldAllocProp {
public:

  FieldAllocProp (const FieldIdentifier& id);
  
  // Request allocation able to accommodate the given ValueType
  template<typename ValueType>
  void request_value_type_allocation ();

  // Locks the properties, preventing furter value types requests
  void commit ();

  // ---- Getters ---- //
  bool is_committed   () const { return m_committed; }
  int  get_alloc_size () const;

  template<typename ValueType>
  bool is_allocation_compatible_with_value_type () const;

protected:

  const FieldIdentifier& m_fid;

  int   m_alloc_value_type_size;
  int   m_scalar_type_size;
  int   m_fast_index_alloc_size;
  int   m_alloc_size;

  bool  m_committed;
};

// ================================= IMPLEMENTATION ================================== //

template<typename ValueType>
void FieldAllocProp::request_value_type_allocation () {

  error::runtime_check(!m_committed, "Error! Cannot change allocation properties after they have been commited.\n");

  if (m_scalar_type_size==0) {
    m_scalar_type_size = sizeof(typename util::ScalarProperties<ValueType>::scalar_type);
    m_alloc_value_type_size = sizeof(ValueType);
  } else {
    // Make sure the scalar_type of the new requested type is at least compatible with the one already stored
    error::runtime_check(sizeof(typename util::ScalarProperties<ValueType>::scalar_type)==m_scalar_type_size,
                         "Error! There was already a value type request for this allocation, and the stored scalar_type_size does not match the one from the new request.\n");
  }

  // If the current value type is large than the currently stored one, update the stored one.
  constexpr int vts = sizeof(ValueType);
  if (vts>m_alloc_value_type_size) {
    // Safety check: old largest value size should divide the new one
    // This is because we are thinking about Pack<T,N> and T's, and we only allow
    // Pack's with N=2^n, for some n>0.
    error::runtime_check(vts%m_alloc_value_type_size == 0, "Error! You are requesting a ValueType whose size is not a multiple of the currently stored value_type_size.\n");

    m_alloc_value_type_size = vts;
  }
}

inline int FieldAllocProp::get_alloc_size () const {
  error::runtime_check(m_committed,"Error! You cannot query the allocation properties until they have been committed.");
  return m_alloc_size;
}

template<typename ValueType>
bool FieldAllocProp::is_allocation_compatible_with_value_type () const {
  constexpr int sts = sizeof(typename util::ScalarProperties<ValueType>::scalar_type);
  constexpr int vts = sizeof(ValueType);

  return sts==m_scalar_type_size && (m_fast_index_alloc_size%vts==0);
}

} // namespace scream

#endif // SCREAM_FIELD_ALLOC_PROP_HPP
