#include "share/field/field.hpp"

namespace scream {

template Field::get_view_type<float,Device> Field::get_view<float,Device> () const;
template Field::get_view_type<float*,Device> Field::get_view<float*,Device> () const;
template Field::get_view_type<float**,Device> Field::get_view<float**,Device> () const;
template Field::get_view_type<float***,Device> Field::get_view<float***,Device> () const;
template Field::get_view_type<float****,Device> Field::get_view<float****,Device> () const;
template Field::get_view_type<float*****,Device> Field::get_view<float*****,Device> () const;
template Field::get_view_type<float******,Device> Field::get_view<float******,Device> () const;

template Field::get_view_type<const float,Device> Field::get_view<const float,Device> () const;
template Field::get_view_type<const float*,Device> Field::get_view<const float*,Device> () const;
template Field::get_view_type<const float**,Device> Field::get_view<const float**,Device> () const;
template Field::get_view_type<const float***,Device> Field::get_view<const float***,Device> () const;
template Field::get_view_type<const float****,Device> Field::get_view<const float****,Device> () const;
template Field::get_view_type<const float*****,Device> Field::get_view<const float*****,Device> () const;
template Field::get_view_type<const float******,Device> Field::get_view<const float******,Device> () const;

} // namespace scream
