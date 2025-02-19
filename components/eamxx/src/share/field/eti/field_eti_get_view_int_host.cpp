#include "share/field/field.hpp"

namespace scream {

template Field::get_view_type<int,Host> Field::get_view<int,Host> () const;
template Field::get_view_type<int*,Host> Field::get_view<int*,Host> () const;
template Field::get_view_type<int**,Host> Field::get_view<int**,Host> () const;
template Field::get_view_type<int***,Host> Field::get_view<int***,Host> () const;
template Field::get_view_type<int****,Host> Field::get_view<int****,Host> () const;
template Field::get_view_type<int*****,Host> Field::get_view<int*****,Host> () const;
template Field::get_view_type<int******,Host> Field::get_view<int******,Host> () const;

template Field::get_view_type<const int,Host> Field::get_view<const int,Host> () const;
template Field::get_view_type<const int*,Host> Field::get_view<const int*,Host> () const;
template Field::get_view_type<const int**,Host> Field::get_view<const int**,Host> () const;
template Field::get_view_type<const int***,Host> Field::get_view<const int***,Host> () const;
template Field::get_view_type<const int****,Host> Field::get_view<const int****,Host> () const;
template Field::get_view_type<const int*****,Host> Field::get_view<const int*****,Host> () const;
template Field::get_view_type<const int******,Host> Field::get_view<const int******,Host> () const;

} // namespace scream
