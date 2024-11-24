#ifndef SCREAM_FIELD_UTILS_IMPL_COLRED_HPP
#define SCREAM_FIELD_UTILS_IMPL_COLRED_HPP

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "share/field/field.hpp"

namespace scream {
namespace impl {

// Utility to compute the reduction of a field along its column dimension.
// This is equivalent to einsum('i,i...k->...k', f1, f2); i is the column.
// The layouts are such that:
// - The first dimension is for the columns (COL)
// - There can be only up to 3 dimensions

template <typename ST>
Field column_reduction(const Field &f1, const Field &f2, const ekat::Comm *co) {
  using KT          = ekat::KokkosTypes<DefaultDevice>;
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using namespace ShortFieldTagsNames;

  const auto &l1 = f1.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG(l1.rank() == 1,
                   "Error! First field f1 must be rank-1.\n"
                   "The input has rank "
                       << l1.rank() << ".\n");
  EKAT_REQUIRE_MSG(l1.tags() == std::vector<FieldTag>({COL}),
                   "Error! First field f1 must have a column dimension.\n"
                   "The input f1 layout is "
                       << l1.tags() << ".\n");

  const auto &n2 = f2.get_header().get_identifier().name();
  const auto &l2 = f2.get_header().get_identifier().get_layout();
  const auto &u2 = f2.get_header().get_identifier().get_units();
  const auto &g2 = f2.get_header().get_identifier().get_grid_name();

  EKAT_REQUIRE_MSG(l2.rank() <= 3,
                   "Error! Second field f2 must be at most rank-3.\n"
                   "The input f2 rank is "
                       << l2.rank() << ".\n");
  EKAT_REQUIRE_MSG(l2.tags()[0] == COL,
                   "Error! Second field f2 must have a column dimension.\n"
                   "The input f2 layout is "
                       << l2.tags() << ".\n");
  EKAT_REQUIRE_MSG(
      l1.dim(0) == l2.dim(0),
      "Error! The two input fields must have the same dimension along "
      "which we are taking the reducing the field.\n"
      "The first field f1 has dimension "
          << l1.dim(0)
          << " while "
             "the second field f2 has dimension "
          << l2.dim(0) << ".\n");

  auto v1 = f1.get_view<const ST *>();

  FieldIdentifier fo_id(n2 + "_colred", l2.clone().strip_dim(0), u2, g2);
  Field fo(fo_id);
  fo.allocate_view();
  fo.deep_copy(0);

  const int d0 = l2.dim(0);

  switch(l2.rank()) {
    case 1: {
      auto v2 = f2.get_view<const ST *>();
      auto vo = fo.get_view<ST>();
      Kokkos::parallel_reduce(
          fo.name(), RangePolicy(0, d0),
          KOKKOS_LAMBDA(const int i, ST &ls) { ls += v1(i) * v2(i); }, vo);
    } break;
    case 2: {
      auto v2      = f2.get_view<const ST **>();
      auto vo      = fo.get_view<ST *>();
      const int d1 = l2.dim(1);
      auto p       = ESU::get_default_team_policy(d1, d0);
      Kokkos::parallel_for(
          fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int j = tm.league_rank();
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, d0),
                [&](int i, ST &ac) { ac += v1(i) * v2(i, j); }, vo(j));
          });
    } break;
    case 3: {
      auto v2      = f2.get_view<const ST ***>();
      auto vo      = fo.get_view<ST **>();
      const int d1 = l2.dim(1);
      const int d2 = l2.dim(2);
      auto p       = ESU::get_default_team_policy(d1 * d2, d0);
      Kokkos::parallel_for(
          fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
            const int idx = tm.league_rank();
            const int j   = idx / d2;
            const int k   = idx % d2;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(tm, d0),
                [&](int i, ST &ac) { ac += v1(i) * v2(i, j, k); }, vo(j, k));
          });
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank.\n");
  }

  if(co) {
    Kokkos::fence();
    fo.sync_to_host();
    co->all_reduce(fo.template get_internal_view_data<ST, Host>(),
                   l2.size() / l2.dim(0), MPI_SUM);
    fo.sync_to_dev();
  }
  return fo;
}

}  // namespace impl
}  // namespace scream

#endif  // SCREAM_FIELD_UTILS_IMPL_COLRED_HPP
