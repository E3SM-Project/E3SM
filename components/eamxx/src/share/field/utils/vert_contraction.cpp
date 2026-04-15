#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"

#include <ekat_comm.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {

namespace impl {

template <typename ST>
void vert_contraction(const Field &f_out, const Field &f_in, const Field &weight)
{
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto l_in  = f_in.get_header().get_identifier().get_layout();
  auto l_w   = weight.get_header().get_identifier().get_layout();

  bool is_masked = f_in.has_valid_mask();

  const int nlevs = l_in.dim(l_in.rank() - 1);

  // To avoid duplicating code for the 1d and 2d weight cases,
  // we use a view to access the weight ahead of time
  Field::view_dev_t<const ST*> w1d;
  Field::view_dev_t<const ST**> w2d;
  auto w_is_1d = l_w.rank() == 1;
  if(w_is_1d) {
    w1d = weight.get_view<const ST *>();
  } else {
    w2d = weight.get_view<const ST **>();
  }

  switch(l_in.rank()) {
    case 1: {
      using mask_t = Field::view_dev_t<const int*>;

      auto mask = is_masked ? f_in.get_valid_mask().get_view<const int *>() : mask_t {};
      auto v_in  = f_in.get_view<const ST *>();
      auto v_out = f_out.get_view<ST>();

      auto policy = Kokkos::RangePolicy<Field::device_t::execution_space>(0,nlevs);
      ST n = 0;
      auto reducer = Kokkos::Sum<ST>(n);
      auto accum = KOKKOS_LAMBDA(const int i, ST &acc) {
        if (not is_masked or mask(i))
          acc += w1d(i) * v_in(i);
      };
      Kokkos::parallel_reduce(f_out.name(), policy, accum, reducer);
      Kokkos::deep_copy(v_out, n);
    } break;
    case 2: {
      using mask_t = Field::view_dev_t<const int**>;

      auto mask  = is_masked ? f_in.get_valid_mask().get_view<const int **>() : mask_t{};
      auto v_in  = f_in.get_view<const ST **>();
      auto v_out = f_out.get_view<ST *>();
      auto d0    = l_in.dim(0);

      auto policy = TPF::get_default_team_policy(d0, nlevs);
      auto outer  =  KOKKOS_LAMBDA(const TeamMember &tm) {
        const int i = tm.league_rank();
        ST n = 0;
        auto inner = [&](int j, ST &acc) {
          if (not is_masked or mask(i,j)) {
            auto w = w_is_1d ? w1d(j) : w2d(i, j);
            acc += w * v_in(i, j);
          }
        };
        auto tvr = Kokkos::TeamVectorRange(tm, nlevs);

        Kokkos::parallel_reduce(tvr, inner, Kokkos::Sum<ST>(n));
        v_out(i) = n;
      };
      Kokkos::parallel_for(f_out.name(), policy, outer);
    } break;
    case 3: {
      using mask_t = Field::view_dev_t<const int***>;

      auto mask  = is_masked ? f_in.get_valid_mask().get_view<const int ***>() : mask_t{};
      auto v_in  = f_in.get_view<const ST ***>();
      auto v_out = f_out.get_view<ST **>();
      auto d0    = l_in.dim(0);
      auto d1    = l_in.dim(1);

      auto policy = TPF::get_default_team_policy(d0*d1, nlevs);
      auto outer =  KOKKOS_LAMBDA(const TeamMember &tm) {
        const int idx = tm.league_rank();
        const int i   = idx / d1;
        const int j   = idx % d1;
        ST n = 0;
        auto inner = [&](int k, ST &acc) {
          if (not is_masked or mask(i,j,k)) {
            auto w = w_is_1d ? w1d(k) : w2d(i, k);
            acc += w * v_in(i, j, k);
          }
        };
        auto tvr = Kokkos::TeamVectorRange(tm, nlevs);

        Kokkos::parallel_reduce(tvr, inner, Kokkos::Sum<ST>(n));
        v_out(i,j) = n;
      };
      Kokkos::parallel_for(f_out.name(), policy, outer);
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank in vert_contraction.\n");
  }
}

} // namespace impl

void vert_contraction (const Field& f_out, const Field& f_in, const Field& weight)
{
  using namespace ShortFieldTagsNames;

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
  EKAT_REQUIRE_MSG (l_in.rank() >= 1 and l_in.rank() <= 3,
    "[vert_contraction] Error! The input field must be rank-1, rank-2, or rank-3.\n"
    " - input field rank: " << l_in.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.rank() >= l_w.rank(),
    "[vert_contraction] Error! The input field must have at least as many dimensions as the weight field.\n"
    " - input field rank : " << l_in.rank() << "\n"
    " - weight field rank: " << l_w.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.tags().back() == LEV or l_in.tags().back() == ILEV,
    "[vert_contraction] Error! The input field layout must have a LEV (or ILEV) as its last dimension.\n"
    " - input field layout: " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_in.dims().back() == l_w.dims().back(),
    "[vert_contraction] Error! The input and weight fields must have the same level dimension size.\n"
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
    "[vert_contraction] Error! The output field layout must match the input field layout (without level tag).\n"
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

  switch(dt) {
    case DataType::IntType:
      impl::vert_contraction<int>(f_out,f_in,weight);
      break;
    case DataType::FloatType:
      impl::vert_contraction<float>(f_out,f_in,weight);
      break;
    case DataType::DoubleType:
      impl::vert_contraction<double>(f_out,f_in,weight);
      break;
    default:
      EKAT_ERROR_MSG ("[vert_contraction] Error! Unsupported data type.\n");
  }
}

} // namespace scream

