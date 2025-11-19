#include "share/field/field.hpp"

namespace scream {

template Field::get_strided_view_type<float,Host> Field::get_strided_view<float,Host> () const;
template Field::get_strided_view_type<float*,Host> Field::get_strided_view<float*,Host> () const;
template Field::get_strided_view_type<float**,Host> Field::get_strided_view<float**,Host> () const;
template Field::get_strided_view_type<float***,Host> Field::get_strided_view<float***,Host> () const;
template Field::get_strided_view_type<float****,Host> Field::get_strided_view<float****,Host> () const;
template Field::get_strided_view_type<float*****,Host> Field::get_strided_view<float*****,Host> () const;
template Field::get_strided_view_type<float******,Host> Field::get_strided_view<float******,Host> () const;

template Field::get_strided_view_type<const float,Host> Field::get_strided_view<const float,Host> () const;
template Field::get_strided_view_type<const float*,Host> Field::get_strided_view<const float*,Host> () const;
template Field::get_strided_view_type<const float**,Host> Field::get_strided_view<const float**,Host> () const;
template Field::get_strided_view_type<const float***,Host> Field::get_strided_view<const float***,Host> () const;
template Field::get_strided_view_type<const float****,Host> Field::get_strided_view<const float****,Host> () const;
template Field::get_strided_view_type<const float*****,Host> Field::get_strided_view<const float*****,Host> () const;
template Field::get_strided_view_type<const float******,Host> Field::get_strided_view<const float******,Host> () const;

} // namespace scream
