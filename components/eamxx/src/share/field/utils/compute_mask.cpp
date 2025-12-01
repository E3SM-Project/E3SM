#include "share/field/field_utils.hpp"

#include "share/field/field.hpp"

#include <type_traits>

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
namespace impl {

template<typename ST>
KOKKOS_INLINE_FUNCTION
bool compare (const ST lhs, const ST rhs, Comparison CMP)
{
  if (CMP==Comparison::EQ) {
    return lhs==rhs;
  } else if (CMP==Comparison::NE) {
    return lhs!=rhs;
  } else if (CMP==Comparison::GT) {
    return lhs>rhs;
  } else if (CMP==Comparison::GE) {
    return lhs>=rhs;
  } else if (CMP==Comparison::LT) {
    return lhs<rhs;
  } else if (CMP==Comparison::LE) {
    return lhs<=rhs;
  } else {
    EKAT_KERNEL_ERROR_MSG ("Error! Unsupported CMP\n");
  }
}

template<int N, typename ST>
void compute_mask (const Field& f, const ST value, Comparison CMP, Field& mask)
{
  using scalar_ND = typename ekat::DataND<const ST,N>::type;
  using int_ND    = typename ekat::DataND<int,N>::type;
  using KT        = Field::kt_dev;
  using exec_space = typename KT::ExeSpace;
  using mdrange_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<N>>;
  using range_t   = typename KT::RangePolicy;

  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto contiguous = f.get_header().get_alloc_properties().contiguous();

  int beg[N] = {};
  int end[N];
  for (int i=0; i<N; ++i) {
    end[i] = layout.dims()[i];
  }

  auto mv = mask.get_view<int_ND>();
  if (contiguous) {
    auto v = f.get_view<scalar_ND>();
    if constexpr (N==1) {
      range_t policy(0,end[0]);
      auto lambda = KOKKOS_LAMBDA (int i) {
        mv(i) = compare(v(i),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==2) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j) {
        mv(i,j) = compare(v(i,j),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==3) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k) {
        mv(i,j,k) = compare(v(i,j,k),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==4) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l) {
        mv(i,j,k,l) = compare(v(i,j,k,l),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==5) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m) {
        mv(i,j,k,l,m) = compare(v(i,j,k,l,m),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==6) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m, int n) {
        mv(i,j,k,l,m,n) = compare(v(i,j,k,l,m,n),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    }
  } else {
    auto v = f.get_strided_view<scalar_ND>();
    if constexpr (N==1) {
      range_t policy(0,end[0]);
      auto lambda = KOKKOS_LAMBDA (int i) {
        mv(i) = compare(v(i),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==2) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j) {
        mv(i,j) = compare(v(i,j),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==3) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k) {
        mv(i,j,k) = compare(v(i,j,k),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==4) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l) {
        mv(i,j,k,l) = compare(v(i,j,k,l),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==5) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m) {
        mv(i,j,k,l,m) = compare(v(i,j,k,l,m),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==6) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m, int n) {
        mv(i,j,k,l,m,n) = compare(v(i,j,k,l,m,n),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    }
  }
}

} // namespace impl

void compute_mask (const Field& f, const ScalarWrapper value, Comparison CMP, Field& mask)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! Input field was not yet allocated.\n");
  EKAT_REQUIRE_MSG (f.rank()>0 and f.rank()<=6,
      "Error! Input field rank not supported.\n");
  EKAT_REQUIRE_MSG (mask.is_allocated(),
      "Error! Mask field was not yet allocated.\n");
  EKAT_REQUIRE_MSG (not mask.is_read_only(),
      "Error! Cannot update mask field, as it is read-only.\n"
      " - mask name: " + mask.name() + "\n");
  EKAT_REQUIRE_MSG (mask.data_type()==DataType::IntType,
      "Error! The data type of the mask field must be 'int'.\n"
      " - mask field name: " << mask.name() << "\n"
      " - mask field data type: " << etoi(mask.data_type()) << "\n");

  const auto& f_layout = f.get_header().get_identifier().get_layout();
  const auto& m_layout = mask.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (m_layout.congruent(f_layout),
      "Error! Mask field layout is incompatible with this field.\n"
      " - field name  : " + f.name() + "\n"
      " - mask name   : " + mask.name() + "\n"
      " - field layout: " + f_layout.to_string() + "\n"
      " - mask layout : " + m_layout.to_string() + "\n");

  const auto f_dt   = f.data_type();
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(value.type,f_dt),
      "[compute_mask] Error! Target value may be narrowed when converted to field data type.\n"
      " - field data type: " + e2str(f_dt) + "\n"
      " - value data type: " + e2str(value.type) + "\n");

  switch (f_dt) {
    case DataType::IntType:
      switch(f.rank()) {
        case 1: impl::compute_mask<1>(f,value.as<int>(),CMP,mask); break;
        case 2: impl::compute_mask<2>(f,value.as<int>(),CMP,mask); break;
        case 3: impl::compute_mask<3>(f,value.as<int>(),CMP,mask); break;
        case 4: impl::compute_mask<4>(f,value.as<int>(),CMP,mask); break;
        case 5: impl::compute_mask<5>(f,value.as<int>(),CMP,mask); break;
        case 6: impl::compute_mask<6>(f,value.as<int>(),CMP,mask); break;
      } break;
    case DataType::FloatType:
      switch(f.rank()) {
        case 1: impl::compute_mask<1>(f,value.as<float>(),CMP,mask); break;
        case 2: impl::compute_mask<2>(f,value.as<float>(),CMP,mask); break;
        case 3: impl::compute_mask<3>(f,value.as<float>(),CMP,mask); break;
        case 4: impl::compute_mask<4>(f,value.as<float>(),CMP,mask); break;
        case 5: impl::compute_mask<5>(f,value.as<float>(),CMP,mask); break;
        case 6: impl::compute_mask<6>(f,value.as<float>(),CMP,mask); break;
      } break;
    case DataType::DoubleType:
      switch(f.rank()) {
        case 1: impl::compute_mask<1>(f,value.as<double>(),CMP,mask); break;
        case 2: impl::compute_mask<2>(f,value.as<double>(),CMP,mask); break;
        case 3: impl::compute_mask<3>(f,value.as<double>(),CMP,mask); break;
        case 4: impl::compute_mask<4>(f,value.as<double>(),CMP,mask); break;
        case 5: impl::compute_mask<5>(f,value.as<double>(),CMP,mask); break;
        case 6: impl::compute_mask<6>(f,value.as<double>(),CMP,mask); break;
      } break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported data type.\n");
  }
}

Field compute_mask (const Field& x, const ScalarWrapper value, Comparison CMP) {
  const auto& fid_x = x.get_header().get_identifier();
  const auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier fid(x.name()+"_mask",fid_x.get_layout(),nondim,fid_x.get_grid_name(),DataType::IntType);
  Field mask(fid,true);

  compute_mask(x,value,CMP,mask);
  return mask;
}

} // namespace scream
