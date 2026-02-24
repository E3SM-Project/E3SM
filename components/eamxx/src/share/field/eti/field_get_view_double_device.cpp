#include "share/field/field.hpp"

namespace scream {

template Field::get_view_type<double,Device> Field::get_view<double,Device> () const;
template Field::get_view_type<double*,Device> Field::get_view<double*,Device> () const;
template Field::get_view_type<double**,Device> Field::get_view<double**,Device> () const;
template Field::get_view_type<double***,Device> Field::get_view<double***,Device> () const;
template Field::get_view_type<double****,Device> Field::get_view<double****,Device> () const;
template Field::get_view_type<double*****,Device> Field::get_view<double*****,Device> () const;
template Field::get_view_type<double******,Device> Field::get_view<double******,Device> () const;

template Field::get_view_type<const double,Device> Field::get_view<const double,Device> () const;
template Field::get_view_type<const double*,Device> Field::get_view<const double*,Device> () const;
template Field::get_view_type<const double**,Device> Field::get_view<const double**,Device> () const;
template Field::get_view_type<const double***,Device> Field::get_view<const double***,Device> () const;
template Field::get_view_type<const double****,Device> Field::get_view<const double****,Device> () const;
template Field::get_view_type<const double*****,Device> Field::get_view<const double*****,Device> () const;
template Field::get_view_type<const double******,Device> Field::get_view<const double******,Device> () const;

} // namespace scream
