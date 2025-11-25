#ifndef EAMXX_FIELD_BROADCAST_HPP
#define EAMXX_FIELD_BROADCAST_HPP

#include "field.hpp"

namespace scream {

/*
 * A specialization of a Field that represent a mask
 *
 * It is just like a field, but:
 * - data type is hard-coded to IntType
 * - units are hard-coded to nondimensional
 * - overload the update method to align with the idea of combining masks
 */

class FieldBroadcast : public Field
{
public:
  FieldBroadcast (const Field& f, const FieldLayout& tgt_layout);

  FieldBroadcast (const FieldBroadcast& src) = default;
  FieldBroadcast& operator= (const FieldBroadcast& src) = default;

  template<typename DT, HostOrDevice HD = Device>
  ekat::ViewBroadcast<get_view_type<DT,HD>>
  get_view () const {
    using ret_t = ekat::ViewBroadcast<get_view_type<DT,HD>>;
    using value_t = typename ret_t::view_type::traits::value_type;
    switch (m_src_field.rank()) {
      case 1:
      {
        using inner_dt = ekat::DataND<value_t,1>;
        auto v = m_src_field.get_view<inner_dt,HD>();
        return broadcast(v,m_extents);
      }
      case 2:
      {
        using inner_dt = ekat::DataND<value_t,1>;
        auto v = m_src_field.get_view<inner_dt,HD>();
        return broadcast(v,m_extents);
      }
    }
  }

protected:
  Field             m_src_field;
  std::vector<int>  m_extents;
};

} // namespace scream

#endif // EAMXX_FIELD_BROADCAST_HPP
