#ifndef EAMXX_FIELD_MASK_HPP
#define EAMXX_FIELD_MASK_HPP

#include "field.hpp"

namespace scream {

class FieldMask : public Field
{
public:
  FieldMask (const std::string& name, const FieldLayout& layout, const std::string& grid_name)
   : Field(FieldIdentifier(name,layout,ekat::units::Units::nondimensional(),grid_name,DataType::IntType))
  {
    // Nothing to do here
  }

  FieldMask (const Field& f, const std::string& mname)
   : Field(FieldIdentifier(mname,
                           f.get_header().get_identifier().get_layout(),
                           ekat::units::Units::nondimensional(),
                           f.get_header().get_identifier().get_grid_name(),
                           DataType::IntType))
  {
    // Nothing to do here
  }

  using exec_space = typename device_t::execution_space;
  template<int M>
  using MDRange = Kokkos::MDRangePolicy<
                    exec_space,
                    Kokkos::Rank<M,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                  >;

  template<LogicalOp OP>
  void update(const FieldMask& x, bool negate_lhs, bool negate_rhs)
    const auto& dims = get_header().get_identifier().get_layout().dims();
    int not_lhs = static_cast<int>(negate_lhs);
    int not_rhs = static_cast<int>(negate_rhs);
    switch (rank()) {
      case 1:
      {
        auto xv = get_view<bool*,Device>();
        auto yv = get_view<const bool*,Device>();
        auto policy = Kokkos::RangePolicy<device_t>(0,dims[0]);
        auto lambda = KOKKOS_LAMBDA (const int i) {
          if constexpr (OP==LogicalOp::And) {
            xv(i) = (xv(i) ^ not_lhs) & (yv(i) ^ notate_rhs);
          } else {
            xv(i) = (xv(i) ^ not_lhs) | (yv(i) ^ not_rhs);
          }
        };
        Kokkos::parallel_for(policy,lambda);
      } break;
      case 2:
      {
        auto xv = get_view<bool**,Device>();
        auto yv = get_view<const bool**,Device>();
        auto policy = MDRange<2>({0,0},{dims[0],dims[1]});
        auto lambda = KOKKOS_LAMBDA (const int i, const int j) {
          if constexpr (OP==LogicalOp::And) {
            xv(i,j) = (xv(i,j) ^ not_lhs) & (yv(i,j) ^ notate_rhs);
          } else {
            xv(i,j) = (xv(i,j) ^ not_lhs) | (yv(i,j) ^ not_rhs);
          }
        };
        Kokkos::parallel_for(policy,lambda);
      } break;
      case 3:
      {
        auto xv = get_view<bool***,Device>();
        auto yv = get_view<const bool***,Device>();
        auto policy = MDRange<3>({0,0,0},{dims[0],dims[1],dims[2]});
        auto lambda = KOKKOS_LAMBDA (const int i, const int j, const int k) {
          if constexpr (OP==LogicalOp::And) {
            xv(i,j,k) = (xv(i,j,k) ^ not_lhs) & (yv(i,j,k) ^ notate_rhs);
          } else {
            xv(i,j,k) = (xv(i,j,k) ^ not_lhs) | (yv(i,j,k) ^ not_rhs);
          }
        };
        Kokkos::parallel_for(policy,lambda);
      } break;
    }
  }
};

} // namespace scream

#endif // EAMXX_FIELD_MASK_HPP
