#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"

#include <ekat_comm.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {

namespace impl {

template <typename ST>
void horiz_contraction(const Field &f_out, const Field &f_in, const Field &weight,
                       bool AVG, const ekat::Comm *comm, const Field &f_tmp)
{
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto l_out = f_out.get_header().get_identifier().get_layout();
  auto l_in  = f_in.get_header().get_identifier().get_layout();

  auto v_w = weight.get_view<const ST *>();

  const int ncols = l_in.dim(0);

  bool is_masked = f_in.get_header().has_extra_data("mask_data");
  bool is_avg_masked = AVG && is_masked && f_out.get_header().has_extra_data("mask_data");
  bool is_comm_avg_masked = comm && is_avg_masked;

  const auto fill_value = is_masked ? f_in.get_header().get_extra_data<Real>("mask_value") : 0;

  switch(l_in.rank()) {
    case 1: {
      auto v_in  = f_in.get_view<const ST *>();
      auto v_m   = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST *>() : v_in;
      auto v_out = f_out.get_view<ST>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST>() : v_out;
      auto v_tmp = is_comm_avg_masked ? f_tmp.get_view<ST>() : v_out;
      ST n = 0, d = 0;
      Kokkos::parallel_reduce(
          f_out.name(), RangePolicy(0, ncols),
          KOKKOS_LAMBDA(const int i, ST &n_acc, ST &d_acc) {
            auto mask = is_masked ? v_m(i) : ST(1.0);
            n_acc += v_w(i) * v_in(i) * mask;
            d_acc += v_w(i) * mask;
          },
          Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
      Kokkos::deep_copy(v_out, n);
      if (is_comm_avg_masked) {
        Kokkos::deep_copy(v_tmp, d);
      } else if (is_avg_masked) {
        ST tmp = d != 0 ? n / d : fill_value;
        Kokkos::deep_copy(v_out, tmp);
        ST mask_val = d != 0 ? 1 : 0;
        Kokkos::deep_copy(v_m_out, mask_val);
      }
    } break;
    case 2: {
      auto v_in    = f_in.get_view<const ST **>();
      auto v_m     = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST **>() : v_in;
      auto v_out   = f_out.get_view<ST *>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST *>() : v_out;
      auto v_tmp   = is_comm_avg_masked ? f_tmp.get_view<ST *>() : v_out;
      const int d1 = l_in.dim(1);
      auto p       = TPF::get_default_team_policy(d1, ncols);
      Kokkos::parallel_for(
          f_out.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int j = tm.league_rank();
            ST n = 0, d = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, ncols),
                [&](int i, ST &n_acc, ST &d_acc) {
                  auto mask = is_masked ? v_m(i, j) : ST(1.0);
                  n_acc += v_w(i) * v_in(i, j) * mask;
                  d_acc += v_w(i) * mask;
                },
                Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
            v_out(j) = n;
            if (is_comm_avg_masked) {
              v_tmp(j) = d;
            } else if (is_avg_masked) {
              v_out(j) = d != 0 ? n / d : fill_value;
              v_m_out(j) = d != 0 ? 1 : 0;
            }
          });
    } break;
    case 3: {
      auto v_in    = f_in.get_view<const ST ***>();
      auto v_m     = is_masked ? f_in.get_header().get_extra_data<Field>("mask_data").get_view<const ST ***>() : v_in;
      auto v_out   = f_out.get_view<ST **>();
      auto v_m_out = is_avg_masked ? f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST **>() : v_out;
      auto v_tmp   = is_comm_avg_masked ? f_tmp.get_view<ST **>() : v_out;
      const int d1 = l_in.dim(1);
      const int d2 = l_in.dim(2);
      auto p       = TPF::get_default_team_policy(d1 * d2, ncols);
      Kokkos::parallel_for(
          f_out.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int idx = tm.league_rank();
            const int j   = idx / d2;
            const int k   = idx % d2;
            ST n = 0, d = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, ncols),
                [&](int i, ST &n_acc, ST &d_acc) {
                  auto mask = is_masked ? v_m(i, j, k) : ST(1.0);
                  n_acc += v_w(i) * v_in(i, j, k) * mask;
                  d_acc += v_w(i) * mask;
                },
                Kokkos::Sum<ST>(n), Kokkos::Sum<ST>(d));
            v_out(j, k) = n;
            if (is_comm_avg_masked) {
              v_tmp(j, k) = d;
            } else if (is_avg_masked) {
              v_out(j, k) = d != 0 ? n / d : fill_value;
              v_m_out(j, k) = d != 0 ? 1 : 0;
            }
          });
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank.\n");
  }

  if(comm) {
    // TODO: use device-side MPI calls
    // TODO: the dev ptr causes problems; revisit this later
    // TODO: doing cuda-aware MPI allreduce would be ~10% faster
    Kokkos::fence();
    f_out.sync_to_host();
    comm->all_reduce(f_out.template get_internal_view_data<ST, Host>(),
                     l_out.size(), MPI_SUM);
    f_out.sync_to_dev();
    if (is_comm_avg_masked) {
      f_tmp.sync_to_host();
      comm->all_reduce(f_tmp.template get_internal_view_data<ST, Host>(),
                      l_out.size(), MPI_SUM);
      f_tmp.sync_to_dev();

      // update f_out by dividing it with f_tmp
      switch(l_out.rank()) {
        case 0: {
          auto v_out = f_out.get_view<ST>();
          auto v_m_out = f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST>();
          auto v_tmp = f_tmp.get_view<const ST>();
          Kokkos::parallel_for(
              f_out.name(), RangePolicy(0, 1),
              KOKKOS_LAMBDA(const int idx) {
                v_out() = v_tmp() != 0 ? v_out() / v_tmp() : fill_value;
                v_m_out() = v_tmp() != 0 ? 1 : 0;
              });
        } break;
        case 1: {
          auto v_out = f_out.get_view<ST *>();
          auto v_m_out = f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST *>();
          auto v_tmp = f_tmp.get_view<const ST *>();
          Kokkos::parallel_for(
              f_out.name(), RangePolicy(0, l_out.dim(0)),
              KOKKOS_LAMBDA(const int i) {
                v_out(i) = v_tmp(i) != 0 ? v_out(i) / v_tmp(i) : fill_value;
                v_m_out(i) = v_tmp(i) != 0 ? 1 : 0;
              });
        } break;
        case 2: {
          auto v_out = f_out.get_view<ST **>();
          auto v_m_out = f_out.get_header().get_extra_data<Field>("mask_data").get_view<ST **>();
          auto v_tmp = f_tmp.get_view<const ST **>();
          const int d0 = l_out.dim(0);
          const int d1 = l_out.dim(1);
          auto p       = TPF::get_default_team_policy(d0*d1, 1);
          Kokkos::parallel_for(
              f_out.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
                const int i = tm.league_rank() / d1;
                const int j = tm.league_rank() % d1;
                v_out(i, j) = v_tmp(i, j) != 0 ? v_out(i, j) / v_tmp(i, j) : fill_value;
                v_m_out(i, j) = v_tmp(i, j) != 0 ? 1 : 0;
              });
        } break;
        default:
          EKAT_ERROR_MSG("Error! Unsupported field rank.\n");
      }
    }
  }
}

} // namespace impl

void horiz_contraction(const Field& f_out, const Field& f_in, const Field& weight,
                       bool AVG, const ekat::Comm* comm, const Field& f_tmp)
{
  using namespace ShortFieldTagsNames;

  // Sanity checks before handing off to the implementation
  EKAT_REQUIRE_MSG (f_out.is_allocated() and f_in.is_allocated() and weight.is_allocated(),
    "[vert_contraction] Error! One (or more) between output, input, and weight field was not yet allocated.\n");

  const auto &l_out = f_out.get_header().get_identifier().get_layout();
  const auto &l_in = f_in.get_header().get_identifier().get_layout();
  const auto &l_w = weight.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (l_w.rank() == 1,
    "[horiz_contraction] Error! The weight field must have rank=1.\n"
    " - weight field rank: " << l_w.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_w.tags()[0] == COL,
    "[horiz_contraction] Error! The weight field must have the COL dimension.\n"
    " - weight Field layout: " << l_w.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_in.rank() <= 3,
    "[horiz_contraction] Error! The input field must be at most rank-3.\n"
    " - input field rank: " << l_in.rank() << ".\n");
  EKAT_REQUIRE_MSG (l_in.tags()[0] == COL,
    "[horiz_contraction] Error! The input field must have the COL tag as first dimension.\n"
    " - input field layout: " << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_w.dim(0) == l_in.dim(0),
    "[horiz_contraction] Error! The input and weight fields must have the same extent along the COL dimension\n"
    " - weight field layout: " << l_w.to_string() << "\n"
    " - input field layout : " << l_in.to_string() << "\n");
  // EKAT_REQUIRE_MSG (l_in.dim(0) > 0,
  //     "Error! The input field must have a non-zero column dimension.\n"
  //     "The input field's layout is "
  //         << l_in.to_string() << ".\n");
  EKAT_REQUIRE_MSG (l_out == l_in.clone().strip_dim(0),
    "[horiz_contraction] Error! The output field layout must match the input field layout (without COL)\n"
    " - input field layout : " << l_in.to_string() << "\n"
    " - output field layout: " << l_out.to_string() << ".\n");
  const auto dt = f_in.data_type();
  EKAT_REQUIRE_MSG (f_out.data_type() == dt,
    "[horiz_contraction] Error! Input and output fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - output field data type: " + e2str(f_out.data_type()) + "\n");
  EKAT_REQUIRE_MSG (weight.data_type()==dt,
    "[horiz_contraction] Error! The input and weight fields must have the same data type.\n"
    " - input field data type : " + e2str(dt) + "\n"
    " - weight field data type: " + e2str(weight.data_type()) + "\n");

  bool is_masked = f_in.get_header().has_extra_data("mask_data");
  bool is_avg_masked = AVG && is_masked && f_out.get_header().has_extra_data("mask_data");
  bool is_comm_avg_masked = comm && is_avg_masked;
  if (is_comm_avg_masked) {
    // make sure f_tmp is allocated correctly
    EKAT_REQUIRE_MSG (f_tmp.is_allocated(),
      "[horiz_contraction] Error! Temporary field was not yet allocated.");
    EKAT_REQUIRE_MSG (f_tmp.data_type() == f_in.data_type(),
      "[horiz_contraction] Error! Temporary field must have the same data type as input field."
      " - input field data type: " + e2str(dt) + "\n"
      " - temp field data type : " + e2str(f_tmp.data_type()) + "\n");
    const auto l_tmp = f_tmp.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG ( l_tmp == l_out,
      "[horiz_contraction] Error! Temporary field must have the same layout as output field."
      " - input field layout: " << l_in.to_string() << "\n"
      " - temp field layout : " << l_tmp.to_string() << "\n");
  }

  switch(dt) {
    case DataType::IntType:
      impl::horiz_contraction<int>(f_out,f_in,weight,AVG,comm,f_tmp);
      break;
    case DataType::FloatType:
      impl::horiz_contraction<float>(f_out,f_in,weight,AVG,comm,f_tmp);
      break;
    case DataType::DoubleType:
      impl::horiz_contraction<double>(f_out,f_in,weight,AVG,comm,f_tmp);
      break;
    default:
      EKAT_ERROR_MSG ("[vert_contraction] Error! Unsupported data type.\n");
  }
}

} // namespace scream
