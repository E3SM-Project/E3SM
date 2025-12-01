#include "share/field/field.hpp"

namespace scream {

template Field::get_strided_view_type<double,Device> Field::get_strided_view<double,Device> () const;
template Field::get_strided_view_type<double*,Device> Field::get_strided_view<double*,Device> () const;
template Field::get_strided_view_type<double**,Device> Field::get_strided_view<double**,Device> () const;
template Field::get_strided_view_type<double***,Device> Field::get_strided_view<double***,Device> () const;
template Field::get_strided_view_type<double****,Device> Field::get_strided_view<double****,Device> () const;
template Field::get_strided_view_type<double*****,Device> Field::get_strided_view<double*****,Device> () const;
template Field::get_strided_view_type<double******,Device> Field::get_strided_view<double******,Device> () const;

template Field::get_strided_view_type<const double,Device> Field::get_strided_view<const double,Device> () const;
template Field::get_strided_view_type<const double*,Device> Field::get_strided_view<const double*,Device> () const;
template Field::get_strided_view_type<const double**,Device> Field::get_strided_view<const double**,Device> () const;
template Field::get_strided_view_type<const double***,Device> Field::get_strided_view<const double***,Device> () const;
template Field::get_strided_view_type<const double****,Device> Field::get_strided_view<const double****,Device> () const;
template Field::get_strided_view_type<const double*****,Device> Field::get_strided_view<const double*****,Device> () const;
template Field::get_strided_view_type<const double******,Device> Field::get_strided_view<const double******,Device> () const;

} // namespace scream
