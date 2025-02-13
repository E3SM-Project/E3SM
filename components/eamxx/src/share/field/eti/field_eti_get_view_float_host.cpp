#include "share/field/field.hpp"

namespace scream {

template Field::get_view_type<float,Host> Field::get_view<float,Host> () const;
template Field::get_view_type<float*,Host> Field::get_view<float*,Host> () const;
template Field::get_view_type<float**,Host> Field::get_view<float**,Host> () const;
template Field::get_view_type<float***,Host> Field::get_view<float***,Host> () const;
template Field::get_view_type<float****,Host> Field::get_view<float****,Host> () const;
template Field::get_view_type<float*****,Host> Field::get_view<float*****,Host> () const;
template Field::get_view_type<float******,Host> Field::get_view<float******,Host> () const;

template Field::get_view_type<const float,Host> Field::get_view<const float,Host> () const;
template Field::get_view_type<const float*,Host> Field::get_view<const float*,Host> () const;
template Field::get_view_type<const float**,Host> Field::get_view<const float**,Host> () const;
template Field::get_view_type<const float***,Host> Field::get_view<const float***,Host> () const;
template Field::get_view_type<const float****,Host> Field::get_view<const float****,Host> () const;
template Field::get_view_type<const float*****,Host> Field::get_view<const float*****,Host> () const;
template Field::get_view_type<const float******,Host> Field::get_view<const float******,Host> () const;

} // namespace scream
