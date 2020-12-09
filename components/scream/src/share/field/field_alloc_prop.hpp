#ifndef SCREAM_FIELD_ALLOC_PROP_HPP
#define SCREAM_FIELD_ALLOC_PROP_HPP

#include "share/field/field_identifier.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_assert.hpp"

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
 *  36 Real's. This class does the book-keeping for the allocation size,
 *  so that 1) the field can be allocated with enough memory to accommodate
 *  all requests, 2) customers of the field can check what the allocation
 *  is, so that they know whether there is padding in the field, and
 *  3) query whether the allocation is compatible with a given value type.
 *
 *  Note: at every request for a new value_type, this class checks the
 *        underlying scalar_type. We ASSUME that the dimensions in the
 *        field identifier refer to that scalar_type.
 */

class FieldAllocProp {
public:

  FieldAllocProp ();
  
  // Request allocation able to accommodate the given ValueType
  template<typename ValueType>
  void request_value_type_allocation ();

  // Locks the properties, preventing furter value types requests
  void commit (const FieldLayout& layout);

  // ---- Getters ---- //
  bool is_committed   () const { return m_committed; }
  int  get_alloc_size () const;
  int  get_last_dim_alloc_size () const;

  template<typename ValueType>
  bool is_allocation_compatible_with_value_type () const;

  // This is here just in case we need it for debugging.
  const std::vector<int>& get_requested_value_types_sizes () const { return m_value_type_sizes; }

protected:

  std::vector<int>    m_value_type_sizes;

  int         m_scalar_type_size;
  std::string m_scalar_type_name;

  int   m_alloc_size;
  int   m_last_dim_alloc_size;

  bool  m_committed;
};

// ================================= IMPLEMENTATION ================================== //

template<typename ValueType>
void FieldAllocProp::request_value_type_allocation () {

  using ekat::ScalarTraits;
  using namespace ekat::error;

  ekat::error::runtime_check(!m_committed, "Error! Cannot change allocation properties after they have been commited.\n");
  
  constexpr int vts = sizeof(ValueType);
  if (m_scalar_type_size==0) {
    // This is the first time we receive a request. Set the scalar type properties
    m_scalar_type_size = sizeof(typename ScalarTraits<ValueType>::scalar_type);
    m_scalar_type_name = ScalarTraits<typename ScalarTraits<ValueType>::scalar_type>::name();
  } else {

    // Make sure the scalar_type of the new requested type coincides with the one already stored
    runtime_check(ScalarTraits<typename ScalarTraits<ValueType>::scalar_type>::name()==m_scalar_type_name,
                  "Error! There was already a value type request for this allocation, and the stored scalar_type name (" +
                  m_scalar_type_name + ") does not match the one from the new request (" + ScalarTraits<ValueType>::name() + ").\n");
    runtime_check(sizeof(typename ScalarTraits<ValueType>::scalar_type)==m_scalar_type_size,
                  "Error! There was already a value type request for this allocation, and the size of the stored scalar_type (" +
                  std::to_string(m_scalar_type_size) + ") does not match the size of the scalar_type of the new request (" +
                  std::to_string(sizeof(typename ScalarTraits<ValueType>::scalar_type)) + ").\n");

    // Furthermore, if value type is not the same as scalar type (by comparing name and size), ValueType *must* be a pack
    const std::string vtn = ScalarTraits<ValueType>::name();
    runtime_check( (m_scalar_type_name==vtn && m_scalar_type_size==vts) || ScalarTraits<ValueType>::is_simd,
                    "Error! Template argument ValueType must be either the ScalarType of this allocation, or a pack type.\n");
  }

  // If ValueType only contains N scalar_type inside, this will pass. If not, then it's a bad ValueType,
  // or you have a wrong scalar_type in the specialization of ScalarTraits<ValueType>. Either way, error out.
  EKAT_REQUIRE_MSG(vts % m_scalar_type_size == 0,
                     "Error! The size of the scalar_type (" + std::to_string(m_scalar_type_size) +
                     ") does not divide the size of the given value_type (" + std::to_string(vts) + ").\n");

  // Store the size of the value type.
  m_value_type_sizes.push_back(vts);
}

inline int FieldAllocProp::get_alloc_size () const {
  ekat::error::runtime_check(m_committed,"Error! You cannot query the allocation properties until they have been committed.");
  return m_alloc_size;
}

inline int FieldAllocProp::get_last_dim_alloc_size () const {
  ekat::error::runtime_check(m_committed,"Error! You cannot query the allocation properties until they have been committed.");
  return m_last_dim_alloc_size;
}

template<typename ValueType>
bool FieldAllocProp::is_allocation_compatible_with_value_type () const {
  using NonConstValueType = typename std::remove_const<ValueType>::type;
  using ekat::ScalarTraits;

  constexpr int  sts = sizeof(typename ScalarTraits<ValueType>::scalar_type);
  constexpr int  vts = sizeof(ValueType);
  const auto stn = ScalarTraits<typename ScalarTraits<NonConstValueType>::scalar_type>::name();

  return stn==m_scalar_type_name
      && sts==m_scalar_type_size && (m_last_dim_alloc_size%vts==0);
}

} // namespace scream

#endif // SCREAM_FIELD_ALLOC_PROP_HPP
