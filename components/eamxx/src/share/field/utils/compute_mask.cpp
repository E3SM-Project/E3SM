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
  // Silence NVCC warning "missing return statement at end of non-void function"
  return false;
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
  const auto lr_ok = f.get_header().get_alloc_properties().allows_layout_right() and
                     mask.get_header().get_alloc_properties().allows_layout_right();

  int beg[N] = {};
  int end[N];
  for (int i=0; i<N; ++i) {
    end[i] = layout.dims()[i];
  }

  if (lr_ok) {
    auto mv = mask.get_view<int_ND>();
    auto v = f.get_view<scalar_ND>();
    if constexpr (N==0) {
      range_t policy(0,1);
      auto lambda = KOKKOS_LAMBDA (int) {
        mv() = compare(v(),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==1) {
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
    auto mv = mask.get_strided_view<int_ND>();
    if constexpr (N==0) {
      range_t policy(0,1);
      auto lambda = KOKKOS_LAMBDA (int) {
        mv() = compare(v(),value,CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==1) {
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

template<int N, typename ST>
void compute_mask (const Field& lhs, const Field& rhs, Comparison CMP, Field& mask)
{
  using scalar_ND = typename ekat::DataND<const ST,N>::type;
  using int_ND    = typename ekat::DataND<int,N>::type;

  using KT         = Field::kt_dev;
  using exec_space = typename KT::ExeSpace;
  using mdrange_t  = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<N>>;
  using range_t    = typename KT::RangePolicy;

  const auto& layout = lhs.get_header().get_identifier().get_layout();
  const auto lr_ok = lhs.get_header().get_alloc_properties().allows_layout_right() and
                     rhs.get_header().get_alloc_properties().allows_layout_right() and
                     mask.get_header().get_alloc_properties().allows_layout_right();

  int beg[N] = {};
  int end[N];
  for (int i=0; i<N; ++i) {
    end[i] = layout.dims()[i];
  }

  if (lr_ok) {
    auto lhs_v = lhs.get_view<scalar_ND>();
    auto rhs_v = rhs.get_view<scalar_ND>();
    auto mv = mask.get_view<int_ND>();
    if constexpr (N==0) {
      range_t policy(0,1);
      auto lambda = KOKKOS_LAMBDA (int) {
        mv() = compare(lhs_v(),rhs_v(),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==1) {
      range_t policy(0,end[0]);
      auto lambda = KOKKOS_LAMBDA (int i) {
        mv(i) = compare(lhs_v(i),rhs_v(i),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==2) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j) {
        mv(i,j) = compare(lhs_v(i,j),rhs_v(i,j),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==3) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k) {
        mv(i,j,k) = compare(lhs_v(i,j,k),rhs_v(i,j,k),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==4) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l) {
        mv(i,j,k,l) = compare(lhs_v(i,j,k,l),rhs_v(i,j,k,l),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==5) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m) {
        mv(i,j,k,l,m) = compare(lhs_v(i,j,k,l,m),rhs_v(i,j,k,l,m),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==6) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m, int n) {
        mv(i,j,k,l,m,n) = compare(lhs_v(i,j,k,l,m,n),rhs_v(i,j,k,l,m,n),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    }
  } else {
    auto lhs_v = lhs.get_strided_view<scalar_ND>();
    auto rhs_v = rhs.get_strided_view<scalar_ND>();
    auto mv = mask.get_strided_view<int_ND>();
    if constexpr (N==0) {
      range_t policy(0,1);
      auto lambda = KOKKOS_LAMBDA (int) {
        mv() = compare(lhs_v(),rhs_v(),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==1) {
      range_t policy(0,end[0]);
      auto lambda = KOKKOS_LAMBDA (int i) {
        mv(i) = compare(lhs_v(i),rhs_v(i),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==2) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j) {
        mv(i,j) = compare(lhs_v(i,j),rhs_v(i,j),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==3) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k) {
        mv(i,j,k) = compare(lhs_v(i,j,k),rhs_v(i,j,k),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==4) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l) {
        mv(i,j,k,l) = compare(lhs_v(i,j,k,l),rhs_v(i,j,k,l),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==5) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m) {
        mv(i,j,k,l,m) = compare(lhs_v(i,j,k,l,m),rhs_v(i,j,k,l,m),CMP);
      };
      Kokkos::parallel_for(policy,lambda);
    } else if constexpr (N==6) {
      mdrange_t policy(beg,end);
      auto lambda = KOKKOS_LAMBDA (int i, int j, int k, int l, int m, int n) {
        mv(i,j,k,l,m,n) = compare(lhs_v(i,j,k,l,m,n),rhs_v(i,j,k,l,m,n),CMP);
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
      "Error! Input field was not yet allocated.\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (f.rank()<=6,
      "Error! Input field rank not supported.\n"
      " - field name: " + f.name() + "\n"
      " - field rank: " + std::to_string(f.rank()) + "\n");
  EKAT_REQUIRE_MSG (mask.is_allocated(),
      "Error! Mask field was not yet allocated.\n"
      " - mask field name: " + mask.name() + "\n");
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
        case 0: impl::compute_mask<0>(f,value.as<int>(),CMP,mask); break;
        case 1: impl::compute_mask<1>(f,value.as<int>(),CMP,mask); break;
        case 2: impl::compute_mask<2>(f,value.as<int>(),CMP,mask); break;
        case 3: impl::compute_mask<3>(f,value.as<int>(),CMP,mask); break;
        case 4: impl::compute_mask<4>(f,value.as<int>(),CMP,mask); break;
        case 5: impl::compute_mask<5>(f,value.as<int>(),CMP,mask); break;
        case 6: impl::compute_mask<6>(f,value.as<int>(),CMP,mask); break;
      } break;
    case DataType::FloatType:
      switch(f.rank()) {
        case 0: impl::compute_mask<0>(f,value.as<float>(),CMP,mask); break;
        case 1: impl::compute_mask<1>(f,value.as<float>(),CMP,mask); break;
        case 2: impl::compute_mask<2>(f,value.as<float>(),CMP,mask); break;
        case 3: impl::compute_mask<3>(f,value.as<float>(),CMP,mask); break;
        case 4: impl::compute_mask<4>(f,value.as<float>(),CMP,mask); break;
        case 5: impl::compute_mask<5>(f,value.as<float>(),CMP,mask); break;
        case 6: impl::compute_mask<6>(f,value.as<float>(),CMP,mask); break;
      } break;
    case DataType::DoubleType:
      switch(f.rank()) {
        case 0: impl::compute_mask<0>(f,value.as<double>(),CMP,mask); break;
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

void compute_mask (const Field& lhs, const Field& rhs, Comparison CMP, Field& mask)
{
  // Sanity checks
  const auto& lhs_layout = lhs.get_header().get_identifier().get_layout();
  const auto& rhs_layout = rhs.get_header().get_identifier().get_layout();
  const auto& m_layout = mask.get_header().get_identifier().get_layout();

  // Layout
  EKAT_REQUIRE_MSG (rhs_layout.congruent(lhs_layout),
      "Error! RHS field layout is incompatible with the lhs field.\n"
      " - lhs name  : " + lhs.name() + "\n"
      " - rhs name  : " + rhs.name() + "\n"
      " - lhs layout: " + lhs_layout.to_string() + "\n"
      " - rhs layout: " + rhs_layout.to_string() + "\n");

  EKAT_REQUIRE_MSG (m_layout.congruent(lhs_layout),
      "Error! Mask field layout is incompatible with the lhs/rhs fields.\n"
      " - lhs name      : " + lhs.name() + "\n"
      " - rhs name      : " + rhs.name() + "\n"
      " - mask name     : " + mask.name() + "\n"
      " - lhs/rhs layout: " + rhs_layout.to_string() + "\n"
      " - mask layout   : " + m_layout.to_string() + "\n");
  EKAT_REQUIRE_MSG (lhs.rank()<=6,
      "Error! LHS field rank not supported.\n"
      " - LHS name: " + lhs.name() + "\n"
      " - LHS rank: " + std::to_string(lhs.rank()) + "\n");

  // Allocation
  EKAT_REQUIRE_MSG (lhs.is_allocated(),
      "Error! LHS field was not yet allocated.\n"
      " - LHS name: " + lhs.name() + "\n");
  EKAT_REQUIRE_MSG (rhs.is_allocated(),
      "Error! RHS field was not yet allocated.\n"
      " - RHS name: " + rhs.name() + "\n");
  EKAT_REQUIRE_MSG (mask.is_allocated(),
      "Error! Mask field was not yet allocated.\n"
      " - mask field name: " + mask.name() + "\n");
  EKAT_REQUIRE_MSG (not mask.is_read_only(),
      "Error! Cannot update mask field, as it is read-only.\n"
      " - mask name: " + mask.name() + "\n");

  // Data type
  EKAT_REQUIRE_MSG (mask.data_type()==DataType::IntType,
      "Error! The data type of the mask field must be 'int'.\n"
      " - mask field name: " << mask.name() << "\n"
      " - mask field data type: " << etoi(mask.data_type()) << "\n");
  EKAT_REQUIRE_MSG (rhs.data_type()==lhs.data_type(),
      "Error! RHS and LHS fields have different data type.\n"
      " - lhs name : " + lhs.name() + "\n"
      " - rhs name : " + rhs.name() + "\n"
      " - lhs dtype: " + e2str(lhs.data_type()) + "\n"
      " - rhs dtype: " + e2str(rhs.data_type()) + "\n");

  switch (lhs.data_type()) {
    case DataType::IntType:
      switch(lhs.rank()) {
        case 0: impl::compute_mask<0,int>(lhs,rhs,CMP,mask); break;
        case 1: impl::compute_mask<1,int>(lhs,rhs,CMP,mask); break;
        case 2: impl::compute_mask<2,int>(lhs,rhs,CMP,mask); break;
        case 3: impl::compute_mask<3,int>(lhs,rhs,CMP,mask); break;
        case 4: impl::compute_mask<4,int>(lhs,rhs,CMP,mask); break;
        case 5: impl::compute_mask<5,int>(lhs,rhs,CMP,mask); break;
        case 6: impl::compute_mask<6,int>(lhs,rhs,CMP,mask); break;
      } break;
    case DataType::FloatType:
      switch(lhs.rank()) {
        case 0: impl::compute_mask<0,float>(lhs,rhs,CMP,mask); break;
        case 1: impl::compute_mask<1,float>(lhs,rhs,CMP,mask); break;
        case 2: impl::compute_mask<2,float>(lhs,rhs,CMP,mask); break;
        case 3: impl::compute_mask<3,float>(lhs,rhs,CMP,mask); break;
        case 4: impl::compute_mask<4,float>(lhs,rhs,CMP,mask); break;
        case 5: impl::compute_mask<5,float>(lhs,rhs,CMP,mask); break;
        case 6: impl::compute_mask<6,float>(lhs,rhs,CMP,mask); break;
      } break;
    case DataType::DoubleType:
      switch(lhs.rank()) {
        case 0: impl::compute_mask<0,double>(lhs,rhs,CMP,mask); break;
        case 1: impl::compute_mask<1,double>(lhs,rhs,CMP,mask); break;
        case 2: impl::compute_mask<2,double>(lhs,rhs,CMP,mask); break;
        case 3: impl::compute_mask<3,double>(lhs,rhs,CMP,mask); break;
        case 4: impl::compute_mask<4,double>(lhs,rhs,CMP,mask); break;
        case 5: impl::compute_mask<5,double>(lhs,rhs,CMP,mask); break;
        case 6: impl::compute_mask<6,double>(lhs,rhs,CMP,mask); break;
      } break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported data type.\n");
  }
}

Field compute_mask (const Field& x, const ScalarWrapper value, Comparison CMP) {
  const auto& fid_x = x.get_header().get_identifier();
  auto fid = fid_x.clone(x.name()+"_mask")
                  .reset_units(ekat::units::none)
                  .reset_dtype(DataType::IntType);
  Field mask(fid,true);

  compute_mask(x,value,CMP,mask);
  return mask;
}

} // namespace scream
