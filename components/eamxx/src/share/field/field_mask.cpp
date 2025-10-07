#include "field_mask.hpp"

namespace scream {

namespace detail {

template<int M>
using MDRange = Kokkos::MDRangePolicy<
                  typename Field::device_t::execution_space,
                  Kokkos::Rank<M,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                >;
} // namespace detail

FieldMask::
FieldMask(const std::string& name, const FieldLayout& layout, const std::string& grid_name, bool allocate)
 : Field(FieldIdentifier(name,layout,ekat::units::Units::nondimensional(),grid_name,DataType::IntType))
{
  if (allocate)
    allocate_view();
}

FieldMask (const Field& f, const std::string& mname, bool allocate)
 : FieldMask(mname,
             f.get_header().get_identifier().get_layout(),
             f.get_header().get_identifier().get_grid_name(),
             allocate)
{
  // Nothing to do here
}

void FieldMask::
update(const FieldMask& x, LogicalOp op, bool negate_lhs, bool negate_rhs)
{
  const auto& dims = get_header().get_identifier().get_layout().dims();
  int not_lhs = static_cast<int>(negate_lhs);
  int not_rhs = static_cast<int>(negate_rhs);
  bool and_op = op==LogicalOp::And;
  switch (rank()) {
    case 1:
    {
      auto xv = get_view<int*,Device>();
      auto yv = get_view<const int*,Device>();
      auto policy = Kokkos::RangePolicy<device_t>(0,dims[0]);
      auto lambda = KOKKOS_LAMBDA (const int i) {
        auto a = xv(i) ^ not_lhs;
        auto b = yv(i) ^ not_rhs;
        xv(i) = and_op ? a & b : a | b;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
    case 2:
    {
      auto xv = get_view<int**,Device>();
      auto yv = get_view<const int**,Device>();
      auto policy = detail::MDRange<2>({0,0},{dims[0],dims[1]});
      auto lambda = KOKKOS_LAMBDA (const int i, const int j) {
        auto a = xv(i,j) ^ not_lhs;
        auto b = yv(i,j) ^ not_rhs;
        xv(i,j) = and_op ? a & b : a | b;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
    case 3:
    {
      auto xv = get_view<int***,Device>();
      auto yv = get_view<const int***,Device>();
      auto policy = detail::MDRange<3>({0,0,0},{dims[0],dims[1],dims[2]});
      auto lambda = KOKKOS_LAMBDA (const int i, const int j, const int k) {
        auto a = xv(i,j.k) ^ not_lhs;
        auto b = yv(i,j.k) ^ not_rhs;
        xv(i,j,k) = and_op ? a & b : a | b;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
  }
}

FieldMask& FieldMask::operator&= (const FieldMask& x)
{
  update(x,true,false,false);
  return *this;
}


FieldMask& FieldMask::operator|= (const FieldMask& x)
{
  update(x,false,false,false);
  return this;
}

FieldMask& FieldMask::flip ()
{
  const auto& dims = get_header().get_identifier().get_layout().dims();
  switch (rank()) {
    case 1:
    {
      auto xv = get_view<int*,Device>();
      auto policy = Kokkos::RangePolicy<device_t>(0,dims[0]);
      auto lambda = KOKKOS_LAMBDA (const int i) {
        xv(i) ^= 1;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
    case 2:
    {
      auto xv = get_view<int**,Device>();
      auto yv = get_view<const int**,Device>();
      auto policy = detail::MDRange<2>({0,0},{dims[0],dims[1]});
      auto lambda = KOKKOS_LAMBDA (const int i, const int j) {
        xv(i,j) ^= 1;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
    case 3:
    {
      auto xv = get_view<int***,Device>();
      auto policy = detail::MDRange<3>({0,0,0},{dims[0],dims[1],dims[2]});
      auto lambda = KOKKOS_LAMBDA (const int i, const int j, const int k) {
        xv(i,j,k) ^= 1;
      };
      Kokkos::parallel_for(policy,lambda);
    } break;
  }
}

FieldMask FieldMask::operator!() const
{
  auto m = clone_mask();
  m.flip();
  return m;
}

FieldMask FieldMask::clone_mask(const std::string& name)
{
  FieldMask m2(*this,name,true);
  m2 |= *this;
  return m2;
}

} // namespace scream
