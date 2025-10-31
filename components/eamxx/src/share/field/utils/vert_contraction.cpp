#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"

#include <ekat_comm.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {

namespace impl {

template <typename ST>
void vert_contraction(const Field &f_out, const Field &f_in, const Field &weight, bool AVG)
{
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto l_out = f_out.get_header().get_identifier().get_layout();
  auto l_in  = f_in.get_header().get_identifier().get_layout();
  auto l_w   = weight.get_header().get_identifier().get_layout();

  bool is_masked = f_in.get_header().has_extra_data("mask_data");
  bool is_avg_masked = AVG && is_masked && f_out.get_header().has_extra_data("mask_data");

  const auto fill_value = is_masked ? f_in.get_header().get_extra_data<Real>("mask_value") : 0;

  const int nlevs = l_in.dim(l_in.rank() - 1);

  // To avoid duplicating code for the 1d and 2d weight cases,
  // we use a view to access the weight ahead of time
  typename Field::get_view_type<const ST *, Device> w1d;
  typename Field::get_view_type<const ST **, Device> w2d;
  auto w_is_1d = l_w.rank() == 1;
  if(w_is_1d) {
    w1d = weight.get_view<const ST *>();
  } else {
    w2d = weight.get_view<const ST **>();
  }

  switch(l_in.rank()) {
    case 1: {
      auto v_in  = f_in.get_view<const ST *>();
      auto v_m   = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST *>() : v_in;
      auto v_out = f_out.get_view<ST>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST>() : v_out;
      ST n = 0, d = 0;
      Kokkos::parallel_reduce(
          f_out.name(), RangePolicy(0, nlevs),
          KOKKOS_LAMBDA(const int i, ST &n_acc, ST &d_acc) {
            auto mask = is_masked ? v_m(i) : ST(1.0);
            auto w = w1d(i);
            n_acc += w * v_in(i) * mask;
            d_acc += w * mask;
          },
          Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
      if (is_avg_masked) {
        ST tmp = d != 0 ? n / d : fill_value;
        Kokkos::deep_copy(v_out, tmp);
        ST mask_val = d != 0 ? 1 : 0;
        Kokkos::deep_copy(v_m_out, mask_val);
      } else {
        Kokkos::deep_copy(v_out, n);
      }
    } break;
    case 2: {
      auto v_in    = f_in.get_view<const ST **>();
      auto v_m     = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST **>() : v_in;
      auto v_out   = f_out.get_view<ST *>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST *>() : v_out;
      const int d0 = l_in.dim(0);
      auto p       = TPF::get_default_team_policy(d0, nlevs);
      Kokkos::parallel_for(
          f_out.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int i = tm.league_rank();
            ST n = 0, d = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, nlevs),
                [&](int j, ST &n_acc, ST &d_acc) {
                  auto mask = is_masked ? v_m(i, j) : ST(1.0);
                  auto w = w_is_1d ? w1d(j) : w2d(i, j); 
                  n_acc += w * v_in(i, j) * mask;
                  d_acc += w * mask;
                },
                Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
            if (is_avg_masked) {
              v_out(i) = d != 0 ? n / d : fill_value;
              v_m_out(i) = d != 0 ? 1 : 0;
            } else {
              v_out(i) = n;
            }
          });
    } break;
    case 3: {
      auto v_in    = f_in.get_view<const ST ***>();
      auto v_m     = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST ***>() : v_in;
      auto v_out   = f_out.get_view<ST **>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST **>() : v_out;
      const int d0 = l_in.dim(0);
      const int d1 = l_in.dim(1);
      auto p       = TPF::get_default_team_policy(d0 * d1, nlevs);
      Kokkos::parallel_for(
          f_out.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int idx = tm.league_rank();
            const int i   = idx / d1;
            const int j   = idx % d1;
            ST n = 0, d = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, nlevs),
                [&](int k, ST &n_acc, ST &d_acc) {
                  auto mask = is_masked ? v_m(i, j, k) : ST(1.0);
                  auto w = w_is_1d ? w1d(k) : w2d(i, k); 
                  n_acc += w * v_in(i, j, k) * mask;
                  d_acc += w * mask;
                },
                Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
            if (is_avg_masked) {
              v_out(i, j) = d != 0 ? n / d : fill_value;
              v_m_out(i, j) = d != 0 ? 1 : 0;
            } else {
              v_out(i, j) = n;
            }
          });
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank in vert_contraction.\n");
  }
}

} // namespace impl

void vert_contraction (const Field& f_out, const Field& f_in, const Field& weight, bool AVG)
{
  using namespace ShortFieldTagsNames;

  // Sanity checks before handing off to the implementation
  EKAT_REQUIRE_MSG (f_out.is_allocated() and f_in.is_allocated() and weight.is_allocated(),
    "[vert_contraction] Error! One (or more) between output, input, and weight field was not yet allocated.\n");

  const auto& l_out = f_out.get_header().get_identifier().get_layout();
  const auto& l_in  = f_in.get_header().get_identifier().get_layout();
  const auto& l_w   = weight.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (l_w.rank() == 1 or l_w.rank() == 2,
    "[vert_contraction] Error! The weight field must be at least rank-1 and at most rank-2.\n"
    " - weight field rank: " << l_w.rank() << "\n");
  EKAT_REQUIRE_MSG (l_w.tags().back() == LEV or l_w.tags().back() == ILEV,
    "[vert_contraction] Error! The weight field must have LEV (or ILEV) as its last dimension.\n"
    " - weight field layout: " << l_w.to_string() << "\n");
  EKAT_REQUIRE_MSG (l_in.rank() <= 3,
    "[vert_contraction] Error! The input field must be at most rank-3.\n"
    " - input field rank: " << l_in.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.rank() >= l_w.rank(),
    "[vert_contraction] Error! The input field must have at least as many dimensions as the weight field.\n"
    " - input field rank : " << l_in.rank() << "\n"
    " - weight field rank: " << l_w.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.tags().back() == LEV or l_in.tags().back() == ILEV,
    "[vert_contraction] Error! The input field layout must have a LEV (or ILEV) as its last dimension.\n"
    " - input field layout: " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_in.dims().back() == l_w.dims().back(),
    "[vert_contraction] Error! The input and weight fields must have the same level dimension tag.\n"
    " - weight field layout: " << l_w.to_string() + "\n"
    " - input field layout : " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_in.dims().back() > 0,
    "[vert_contraction] Error! The input field must have a non-zero level dimension.\n"
    " - input field layout: " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_w.rank()==1 or l_w.congruent(l_in.clone().strip_dim(CMP, false)),
    "[vert_contraction] Error! Input and weight fields have incompatible layouts.\n"
    " - input field layout : " + l_in.to_string() + "\n"
    " - weight field layout: " + l_w.to_string() + "\n");
  EKAT_REQUIRE_MSG (l_out == l_in.clone().strip_dim(l_in.rank() - 1),
    "[vert_contraction] Error! The output field layout must match the input field layout (without level tag)\n"
    " - input field layout : " + l_in.to_string() + "\n"
    " - output field layout: " + l_out.to_string() + ".\n");
  
  const auto dt = f_in.data_type();
  EKAT_REQUIRE_MSG (f_out.data_type()==dt,
    "[vert_contraction] Error! The input and output fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - output field data type: " + e2str(f_out.data_type()) + "\n");
  EKAT_REQUIRE_MSG (weight.data_type()==dt,
    "[vert_contraction] Error! The input and weight fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - weight field data type: " + e2str(weight.data_type()) + "\n");
  EKAT_REQUIRE_MSG (not AVG or dt==DataType::FloatType or dt==DataType::DoubleType,
    "[vert_contraction] Error! When AVG=true, the fields data type MUST be a floating point type.\n"
    " - input field data type : " + e2str(dt) + "\n");

  switch(dt) {
    case DataType::IntType:
      impl::vert_contraction<int>(f_out,f_in,weight,AVG);
      break;
    case DataType::FloatType:
      impl::vert_contraction<float>(f_out,f_in,weight,AVG);
      break;
    case DataType::DoubleType:
      impl::vert_contraction<double>(f_out,f_in,weight,AVG);
      break;
    default:
      EKAT_ERROR_MSG ("[vert_contraction] Error! Unsupported data type.\n");
  }
}

} // namespace scream
