#include "share/field/field.hpp"

namespace scream {

template Field::get_strided_view_type<double,Host> Field::get_strided_view<double,Host> () const;
template Field::get_strided_view_type<double*,Host> Field::get_strided_view<double*,Host> () const;
template Field::get_strided_view_type<double**,Host> Field::get_strided_view<double**,Host> () const;
template Field::get_strided_view_type<double***,Host> Field::get_strided_view<double***,Host> () const;
template Field::get_strided_view_type<double****,Host> Field::get_strided_view<double****,Host> () const;
template Field::get_strided_view_type<double*****,Host> Field::get_strided_view<double*****,Host> () const;
template Field::get_strided_view_type<double******,Host> Field::get_strided_view<double******,Host> () const;

template Field::get_strided_view_type<const double,Host> Field::get_strided_view<const double,Host> () const;
template Field::get_strided_view_type<const double*,Host> Field::get_strided_view<const double*,Host> () const;
template Field::get_strided_view_type<const double**,Host> Field::get_strided_view<const double**,Host> () const;
template Field::get_strided_view_type<const double***,Host> Field::get_strided_view<const double***,Host> () const;
template Field::get_strided_view_type<const double****,Host> Field::get_strided_view<const double****,Host> () const;
template Field::get_strided_view_type<const double*****,Host> Field::get_strided_view<const double*****,Host> () const;
template Field::get_strided_view_type<const double******,Host> Field::get_strided_view<const double******,Host> () const;

} // namespace scream
