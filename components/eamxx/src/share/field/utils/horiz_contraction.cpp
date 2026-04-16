#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"

#include <ekat_comm.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {

namespace impl {

template <typename ST>
void horiz_contraction(const Field &f_out, const Field &f_in, const Field &weight)
{
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto l_in  = f_in.get_header().get_identifier().get_layout();

  auto v_w = weight.get_view<const ST *>();

  const int ncols = l_in.dim(0);

  bool is_masked = f_in.has_valid_mask();

  switch(l_in.rank()) {
    case 1: {
      using mask_t = Field::view_dev_t<const int*>;

      auto v_in  = f_in.get_view<const ST *>();
      auto v_out = f_out.get_view<ST>();
      auto mask  = is_masked ? f_in.get_valid_mask().get_view<const int*>() : mask_t{};

      ST n = 0;
      auto policy = Kokkos::RangePolicy<Field::device_t::execution_space>(0,ncols);
      auto reducer = Kokkos::Sum<ST>(n);
      auto accum = KOKKOS_LAMBDA(const int i, ST &acc) {
        if (not is_masked or mask(i))
          acc += v_w(i) * v_in(i);
      };
      Kokkos::parallel_reduce(f_out.name(), policy, accum, reducer);
      Kokkos::deep_copy(v_out, n);
    } break;
    case 2: {
      using mask_t = Field::view_dev_t<const int**>;

      auto v_in   = f_in.get_view<const ST **>();
      auto v_out  = f_out.get_view<ST *>();
      auto d1     = l_in.dim(1);

      auto mask  = is_masked ? f_in.get_valid_mask().get_view<const int**>() : mask_t{};
      auto policy = TPF::get_default_team_policy(d1, ncols);
      auto outer =  KOKKOS_LAMBDA(const TeamMember &tm) {
        const int j = tm.league_rank();
        ST n = 0;
        auto inner = [&](int i, ST &acc) {
          if (not is_masked or mask(i,j)) {
            acc += v_w(i) * v_in(i, j);
          }
        };
        auto tvr = Kokkos::TeamVectorRange(tm, ncols);

        Kokkos::parallel_reduce(tvr, inner, Kokkos::Sum<ST>(n));
        v_out(j) = n;
      };
      Kokkos::parallel_for(f_out.name(), policy, outer);
    } break;
    case 3: {
      using mask_t = Field::view_dev_t<const int***>;

      auto v_in   = f_in.get_view<const ST ***>();
      auto v_out  = f_out.get_view<ST **>();
      auto d1     = l_in.dim(1);
      auto d2     = l_in.dim(2);

      auto mask  = is_masked ? f_in.get_valid_mask().get_view<const int***>() : mask_t{};
      auto policy = TPF::get_default_team_policy(d1*d2, ncols);
      auto outer =  KOKKOS_LAMBDA(const TeamMember &tm) {
        const int idx = tm.league_rank();
        const int j   = idx / d2;
        const int k   = idx % d2;
        ST n = 0;
        auto inner = [&](int i, ST &acc) {
          if (not is_masked or mask(i,j,k)) {
            acc += v_w(i) * v_in(i, j, k);
          }
        };
        auto tvr = Kokkos::TeamVectorRange(tm, ncols);

        Kokkos::parallel_reduce(tvr, inner, Kokkos::Sum<ST>(n));
        v_out(j,k) = n;
      };
      Kokkos::parallel_for(f_out.name(), policy, outer);
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank in horiz_contraction.\n");
  }
}

} // namespace impl

void horiz_contraction(const Field& f_out, const Field& f_in, const Field& weight)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (f_out.is_allocated() and f_in.is_allocated() and weight.is_allocated(),
    "[horiz_contraction] Error! One (or more) between output, input, and weight field was not yet allocated.\n");

  const auto &l_out = f_out.get_header().get_identifier().get_layout();
  const auto &l_in  = f_in.get_header().get_identifier().get_layout();
  const auto &l_w   = weight.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (l_w.rank() == 1,
    "[horiz_contraction] Error! The weight field must have rank=1.\n"
    " - weight field rank: " << l_w.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_w.tags()[0] == COL,
    "[horiz_contraction] Error! The weight field must have the COL dimension.\n"
    " - weight Field layout: " << l_w.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_in.rank() >= 1 and l_in.rank() <= 3,
    "[horiz_contraction] Error! The input field must be rank-1, rank-2, or rank-3.\n"
    " - input field rank: " << l_in.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.tags()[0] == COL,
    "[horiz_contraction] Error! The input field must have the COL tag as first dimension.\n"
    " - input field layout: " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_w.dim(0) == l_in.dim(0),
    "[horiz_contraction] Error! The input and weight fields must have the same extent along the COL dimension.\n"
    " - weight field layout: " << l_w.to_string() << "\n"
    " - input field layout : " << l_in.to_string() << "\n");
  EKAT_REQUIRE_MSG (l_out == l_in.clone().strip_dim(0),
    "[horiz_contraction] Error! The output field layout must match the input field layout (without COL).\n"
    " - input field layout : " << l_in.to_string() << "\n"
    " - output field layout: " << l_out.to_string() << ".\n");

  const auto dt = f_in.data_type();
  EKAT_REQUIRE_MSG (f_out.data_type() == dt,
    "[horiz_contraction] Error! Input and output fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - output field data type: " + e2str(f_out.data_type()) + "\n");
  EKAT_REQUIRE_MSG (weight.data_type() == dt,
    "[horiz_contraction] Error! The input and weight fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - weight field data type: " + e2str(weight.data_type()) + "\n");

  switch(dt) {
    case DataType::IntType:
      impl::horiz_contraction<int>(f_out, f_in, weight);
      break;
    case DataType::FloatType:
      impl::horiz_contraction<float>(f_out, f_in, weight);
      break;
    case DataType::DoubleType:
      impl::horiz_contraction<double>(f_out, f_in, weight);
      break;
    default:
      EKAT_ERROR_MSG ("[horiz_contraction] Error! Unsupported data type.\n");
  }
}

void horiz_contraction(const Field& f_out, const Field& f_in,
                       const Field& weight, const ekat::Comm& comm)
{
  horiz_contraction(f_out, f_in, weight);

  const auto& l_out = f_out.get_header().get_identifier().get_layout();
  // TODO: use device-aware MPI calls
  // TODO: doing cuda-aware MPI allreduce would be ~10% faster
  Kokkos::fence();
  f_out.sync_to_host();
  switch (f_out.data_type()) {
    case DataType::IntType:
      comm.all_reduce(f_out.get_internal_view_data<int, Host>(), l_out.size(), MPI_SUM);
      break;
    case DataType::FloatType:
      comm.all_reduce(f_out.get_internal_view_data<float, Host>(), l_out.size(), MPI_SUM);
      break;
    case DataType::DoubleType:
      comm.all_reduce(f_out.get_internal_view_data<double, Host>(), l_out.size(), MPI_SUM);
      break;
    default:
      EKAT_ERROR_MSG ("[horiz_contraction] Error! Unsupported data type.\n");
  }
  f_out.sync_to_dev();
}

} // namespace scream

