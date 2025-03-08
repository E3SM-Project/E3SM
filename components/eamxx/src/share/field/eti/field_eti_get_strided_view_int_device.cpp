#include "share/field/field.hpp"

namespace scream {

template Field::get_strided_view_type<int,Device> Field::get_strided_view<int,Device> () const;
template Field::get_strided_view_type<int*,Device> Field::get_strided_view<int*,Device> () const;
template Field::get_strided_view_type<int**,Device> Field::get_strided_view<int**,Device> () const;
template Field::get_strided_view_type<int***,Device> Field::get_strided_view<int***,Device> () const;
template Field::get_strided_view_type<int****,Device> Field::get_strided_view<int****,Device> () const;
template Field::get_strided_view_type<int*****,Device> Field::get_strided_view<int*****,Device> () const;
template Field::get_strided_view_type<int******,Device> Field::get_strided_view<int******,Device> () const;

template Field::get_strided_view_type<const int,Device> Field::get_strided_view<const int,Device> () const;
template Field::get_strided_view_type<const int*,Device> Field::get_strided_view<const int*,Device> () const;
template Field::get_strided_view_type<const int**,Device> Field::get_strided_view<const int**,Device> () const;
template Field::get_strided_view_type<const int***,Device> Field::get_strided_view<const int***,Device> () const;
template Field::get_strided_view_type<const int****,Device> Field::get_strided_view<const int****,Device> () const;
template Field::get_strided_view_type<const int*****,Device> Field::get_strided_view<const int*****,Device> () const;
template Field::get_strided_view_type<const int******,Device> Field::get_strided_view<const int******,Device> () const;

} // namespace scream
