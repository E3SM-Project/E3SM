#ifndef SCREAM_FIELD_UTILS_IMPL_HPP
#define SCREAM_FIELD_UTILS_IMPL_HPP

#include "share/field/field.hpp"

#include <ekat_comm.hpp>
#include <ekat_team_policy_utils.hpp>

#include <limits>
#include <type_traits>

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
namespace impl {

template<typename ST>
bool views_are_equal(const Field& f1, const Field& f2, const ekat::Comm* comm)
{
  // Get physical layout (shoudl be the same for both fields)
  const auto& l1 = f1.get_header().get_identifier().get_layout();
  const auto& l2 = f2.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (l1==l2,
      "Error! Input fields have different layouts.\n");

  // For simplicity, we perform the check on Host only. This is not a big
  // limitation, since this code is likely used only in testing.
  f1.sync_to_host();
  f2.sync_to_host();

  // Reshape based on the rank, then loop over all entries.
  // NOTE: because views_are_equal() is only used for testing, we generalize by
  // always calling get_strided_view(), even if it could be contiguous.
  // This allows us to call this function on any type of field,
  // including multi-slice subfields
  bool same_locally = true;
  const auto& dims = l1.dims();
  switch (l1.rank()) {
    case 0:
      {
        auto v1 = f1.template get_strided_view<ST,Host>();
        auto v2 = f2.template get_strided_view<ST,Host>();
        if (v1() != v2()) {
          same_locally = false;
          break;
        }
        break;
      }
    case 1:
      {
        auto v1 = f1.template get_strided_view<ST*,Host>();
        auto v2 = f2.template get_strided_view<ST*,Host>();
        for (int i=0; i<dims[0]; ++i) {
          if (v1(i) != v2(i)) {
            same_locally = false;
            break;
          }
        }
      }
      break;
    case 2:
      {
        auto v1 = f1.template get_strided_view<ST**,Host>();
        auto v2 = f2.template get_strided_view<ST**,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            if (v1(i,j) != v2(i,j)) {
              same_locally = false;
              break;
            }
        }}
      }
      break;
    case 3:
      {
        auto v1 = f1.template get_strided_view<ST***,Host>();
        auto v2 = f2.template get_strided_view<ST***,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              if (v1(i,j,k) != v2(i,j,k)) {
                same_locally = false;
                break;
              }
        }}}
      }
      break;
    case 4:
      {
        auto v1 = f1.template get_strided_view<ST****,Host>();
        auto v2 = f2.template get_strided_view<ST****,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                if (v1(i,j,k,l) != v2(i,j,k,l)) {
                  same_locally = false;
                  break;
                }
        }}}}
      }
      break;
    case 5:
      {
        auto v1 = f1.template get_strided_view<ST*****,Host>();
        auto v2 = f2.template get_strided_view<ST*****,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; same_locally && l<dims[3]; ++l) {
                for (int m=0; m<dims[4]; ++m) {
                  if (v1(i,j,k,l,m) != v2(i,j,k,l,m)) {
                    same_locally = false;
                    break;
                  }
        }}}}}
      }
      break;
    case 6:
      {
        auto v1 = f1.template get_strided_view<ST******,Host>();
        auto v2 = f2.template get_strided_view<ST******,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; same_locally && l<dims[3]; ++l) {
                for (int m=0; same_locally && m<dims[4]; ++m) {
                  for (int n=0; n<dims[5]; ++n) {
                    if (v1(i,j,k,l,m,n) != v2(i,j,k,l,m,n)) {
                      same_locally = false;
                      break;
                    }
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    bool same_globally;
    comm->all_reduce(&same_locally,&same_globally,1,MPI_LAND);
    return same_globally;
  } else {
    return same_locally;
  }
}

template<typename ST, typename Engine, typename PDF>
void randomize (const Field& f, Engine& engine, PDF&& pdf)
{
  const auto& fl = f.get_header().get_identifier().get_layout();
  switch (fl.rank()) {
    case 0:
      {
        auto v = f.template get_strided_view<ST,Host>();
        v() = pdf(engine);
      }
      break;
    case 1:
      {
        auto v = f.template get_strided_view<ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          v(i) = pdf(engine);
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            v(i,j) = pdf(engine);
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              v(i,j,k) = pdf(engine);
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<ST****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                v(i,j,k,l) = pdf(engine);
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<ST*****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  v(i,j,k,l,m) = pdf(engine);
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<ST******,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  for (int n=0; n<fl.dim(5); ++n) {
                    v(i,j,k,l,m,n) = pdf(engine);
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  // Sync the dev view with the host view.
  f.sync_to_dev();
}

template<typename ST, typename Engine, typename PDF, typename MaskType>
void perturb (Field& f,
              Engine& engine,
              PDF&& pdf,
              const unsigned int base_seed,
              const MaskType& level_mask,
              const Field& dof_gids)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  // Check to see if field has a column dimension
  using namespace ShortFieldTagsNames;
  const bool has_column_dim = fl.has_tag(COL);
  const bool has_lev_dim = fl.has_tag(LEV);

  if (has_column_dim) {
    // Because Column is the partitioned dimension, we must reset the
    // RNG seed to be the same on every column so that a column will
    // have the same value no matter where it exists in an MPI rank's
    // set of local columns.
    const auto gids = dof_gids.get_strided_view<const int*, Host>();

    // Create a field to store perturbation values with layout
    // the same as f, but stripped of column and level dimension.
    auto perturb_fl = fl.clone().strip_dims({COL,LEV});
    FieldIdentifier perturb_fid("perturb_field", perturb_fl, ekat::units::Units::nondimensional(), "");
    Field perturb_f(perturb_fid);
    perturb_f.allocate_view();

    // Loop through columns as reset RNG seed based on GID of column
    for (auto icol=0; icol<fl.dims().front(); ++icol) {
      const auto new_seed = base_seed+gids(icol);
      engine.seed(new_seed);

      if (has_lev_dim) {
        // Loop through levels. For each that satisfy the level_mask,
        // apply a random perturbation to f.
        for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
          if (level_mask(ilev)) {
            randomize(perturb_f, engine, pdf);
            f.subfield(COL, icol).subfield(LEV, ilev).scale(perturb_f);
          }
        }
      } else {
        randomize(perturb_f, engine, pdf);
        f.subfield(COL, icol).scale(perturb_f);
      }
    }
  } else {
    // If no Column tag exists, this field is not partitioned.
    // Set engine to base_seed to ensure computation is reproducible.
    engine.seed(base_seed);

    // Create a field to store perturbation values with layout
    // the same as f, but stripped of level dimension.
    auto perturb_fl = fl.clone().strip_dim(LEV);
    FieldIdentifier perturb_fid("perturb_field", perturb_fl, ekat::units::Units::nondimensional(), "");
    Field perturb_f(perturb_fid);
    perturb_f.allocate_view();

    if (has_lev_dim) {
      // Loop through levels. For each that satisfy the level_mask,
      // apply a random perturbation to f.
      for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
        if (level_mask(ilev)) {
          randomize(perturb_f, engine, pdf);
          f.subfield(LEV, ilev).scale(perturb_f);
        }
      }
    } else {
      randomize(perturb_f, engine, pdf);
      f.scale(perturb_f);
    }
  }
}

template <typename ST, int AVG>
void horiz_contraction(const Field &f_out, const Field &f_in,
                       const Field &weight, const ekat::Comm *comm, const Field &f_tmp) {
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

  if (is_comm_avg_masked) {
    // make sure f_tmp is allocated correctly
    EKAT_REQUIRE_MSG(f_tmp.is_allocated(),
                     "Error! Temporary field must be allocated.");
    EKAT_REQUIRE_MSG(f_tmp.data_type() == f_in.data_type(),
                     "Error! Temporary field must have the same data type as input fields.");
    EKAT_REQUIRE_MSG(f_tmp.get_header().get_identifier().get_layout() == l_out,
                     "Error! Temporary field must have the same layout as output field.");
  }

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

template <typename ST, int AVG>
void vert_contraction(const Field &f_out, const Field &f_in, const Field &weight) {
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

template<typename ST>
ST frobenius_norm(const Field& f, const ekat::Comm* comm)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  // TODO: compute directly on device
  f.sync_to_host();

  // Note: use Kahan algorithm to increase accuracy
  ST norm = 0;
  ST c = 0;
  ST temp,y;
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_strided_view<const ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          y = std::pow(v(i),2) - c;
          temp = norm + y;
          c = (temp - norm) - y;
          norm = temp;
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<const ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            y = std::pow(v(i,j),2) - c;
            temp = norm + y;
            c = (temp - norm) - y;
            norm = temp;
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<const ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              y = std::pow(v(i,j,k),2) - c;
              temp = norm + y;
              c = (temp - norm) - y;
              norm = temp;
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<const ST****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                y = std::pow(v(i,j,k,l),2) - c;
                temp = norm + y;
                c = (temp - norm) - y;
                norm = temp;
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<const ST*****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  y = std::pow(v(i,j,k,l,m),2) - c;
                  temp = norm + y;
                  c = (temp - norm) - y;
                  norm = temp;
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<const ST******,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  for (int n=0; n<fl.dim(5); ++n) {
                    y = std::pow(v(i,j,k,l,m,n),2) - c;
                    temp = norm + y;
                    c = (temp - norm) - y;
                    norm = temp;
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    ST global_norm;
    comm->all_reduce(&norm,&global_norm,1,MPI_SUM);
    return std::sqrt(global_norm);
  } else {
    return std::sqrt(norm);
  }
}

template<typename ST>
ST field_sum(const Field& f, const ekat::Comm* comm)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  // TODO: compute directly on device
  f.sync_to_host();

  // Note: use Kahan algorithm to increase accuracy
  ST sum = 0;
  ST c = 0;
  ST temp,y;
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_strided_view<const ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          y = v(i) - c;
          temp = sum + y;
          c = (temp - sum) - y;
          sum = temp;
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<const ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            y = v(i,j) - c;
            temp = sum + y;
            c = (temp - sum) - y;
            sum = temp;
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<const ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              y = v(i,j,k) - c;
              temp = sum + y;
              c = (temp - sum) - y;
              sum = temp;
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<const ST****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                y = v(i,j,k,l) - c;
                temp = sum + y;
                c = (temp - sum) - y;
                sum = temp;
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<const ST*****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  y = v(i,j,k,l,m) - c;
                  temp = sum + y;
                  c = (temp - sum) - y;
                  sum = temp;
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<const ST******,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  for (int n=0; n<fl.dim(5); ++n) {
                    y = v(i,j,k,l,m,n) - c;
                    temp = sum + y;
                    c = (temp - sum) - y;
                    sum = temp;
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    ST global_sum;
    comm->all_reduce(&sum,&global_sum,1,MPI_SUM);
    return global_sum;
  } else {
    return sum;
  }
}

template<typename ST>
ST field_max(const Field& f, const ekat::Comm* comm)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  // TODO: compute directly on device
  f.sync_to_host();

  ST max = std::numeric_limits<ST>::lowest();
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_strided_view<const ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          max = std::max(max,v(i));
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<const ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            max = std::max(max,v(i,j));
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<const ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              max = std::max(max,v(i,j,k));
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<const ST****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                max = std::max(max,v(i,j,k,l));
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<const ST*****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  max = std::max(max,v(i,j,k,l,m));
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<const ST******,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  for (int n=0; n<fl.dim(5); ++n) {
                    max = std::max(max,v(i,j,k,l,m,n));
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    ST global_max;
    comm->all_reduce(&max,&global_max,1,MPI_MAX);
    return global_max;
  } else {
    return max;
  }
}

template<typename ST>
ST field_min(const Field& f, const ekat::Comm* comm)
{
  const auto& fl = f.get_header().get_identifier().get_layout();

  // TODO: compute directly on device
  f.sync_to_host();

  ST min = std::numeric_limits<ST>::max();
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_strided_view<const ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          min = std::min(min,v(i));
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<const ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            min = std::min(min,v(i,j));
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<const ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              min = std::min(min,v(i,j,k));
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<const ST****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                min = std::min(min,v(i,j,k,l));
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<const ST*****,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  min = std::min(min,v(i,j,k,l,m));
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<const ST******,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              for (int l=0; l<fl.dim(3); ++l) {
                for (int m=0; m<fl.dim(4); ++m) {
                  for (int n=0; n<fl.dim(5); ++n) {
                    min = std::min(min,v(i,j,k,l,m,n));
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    ST global_min;
    comm->all_reduce(&min,&global_min,1,MPI_MIN);
    return global_min;
  } else {
    return min;
  }
}

template<Comparison CMP, typename ViewT, typename MaskT>
struct SetMaskHelper {
  using exec_space = typename ViewT::traits::execution_space;
  using ST     = typename ViewT::traits::non_const_value_type;
  using MaskST = typename MaskT::traits::non_const_value_type;

  template<int M>
  using MDRange = Kokkos::MDRangePolicy<
                    exec_space,
                    Kokkos::Rank<M,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                  >;

  static constexpr int N = ViewT::rank;

  void run (const std::vector<int>& dims) const {
    if constexpr (N==0) {
      Kokkos::RangePolicy<exec_space> policy(0,1);
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==1) {
      Kokkos::RangePolicy<exec_space> policy(0,dims[0]);
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==2) {
      MDRange<2> policy({0,0},{dims[0],dims[1]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==3) {
      MDRange<3> policy({0,0,0},{dims[0],dims[1],dims[2]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==4) {
      MDRange<4> policy({0,0,0,0},{dims[0],dims[1],dims[2],dims[3]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==5) {
      MDRange<5> policy({0,0,0,0,0},{dims[0],dims[1],dims[2],dims[3],dims[4]});
      Kokkos::parallel_for(policy,*this);
    } else if constexpr (N==6) {
      MDRange<6> policy({0,0,0,0,0,0},{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
      Kokkos::parallel_for(policy,*this);
    } else {
      EKAT_ERROR_MSG ("Unsupported rank! Should be in [2,6].\n");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void set_mask (const ST& x_val, MaskST& m_val) const {
    if constexpr (CMP==Comparison::EQ)
      m_val = static_cast<int>(x_val==val);
    else if constexpr (CMP==Comparison::NE)
      m_val = static_cast<int>(x_val!=val);
    else if constexpr (CMP==Comparison::GT)
      m_val = static_cast<int>(x_val>val);
    else if constexpr (CMP==Comparison::GE)
      m_val = static_cast<int>(x_val>=val);
    else if constexpr (CMP==Comparison::LT)
      m_val = static_cast<int>(x_val<val);
    else if constexpr (CMP==Comparison::LE)
      m_val = static_cast<int>(x_val<=val);
    else
      EKAT_KERNEL_ERROR_MSG ("Unsupported Comparison value.\n");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if constexpr (N==0)
      set_mask (x(),m());
    else
      set_mask (x(i),m(i));
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j) const {
    set_mask (x(i,j),m(i,j));
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k) const {
    set_mask (x(i,j,k),m(i,j,k));
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l) const {
    set_mask (x(i,j,k,l),m(i,j,k,l));
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int n) const {
    set_mask (x(i,j,k,l,n),m(i,j,k,l,n));
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (int i, int j, int k, int l, int n, int p) const {
    set_mask (x(i,j,k,l,n,p),m(i,j,k,l,n,p));
  }

  ViewT x;
  MaskT m;
  ST    val;
};

template<Comparison CMP, typename ViewT, typename MaskT>
void setMaskHelper (const ViewT& x, const MaskT& m,
                    const FieldHeader& mh,
                    typename ViewT::traits::value_type val,
                    const std::vector<int>& dims)
{
  SetMaskHelper<CMP,ViewT,MaskT> helper;
  helper.m = m;
  helper.x = x;
  helper.val = val;
  helper.run(dims);
}

template<Comparison CMP, typename ST>
void compute_mask (const Field& x, const ST value, Field& m)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();
  const auto contiguous = x.get_header().get_alloc_properties().contiguous();
  const auto& mh = m.get_header();

  switch (layout.rank()) {
    case 0:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST>(),m.get_view<int>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST>(),m.get_view<int>(),mh,value,dims);
      break;
    case 1:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST*>(),m.get_view<int*>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST*>(),m.get_view<int*>(),mh,value,dims);
      break;
    case 2:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST**>(),m.get_view<int**>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST**>(),m.get_view<int**>(),mh,value,dims);
      break;
    case 3:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST***>(),m.get_view<int***>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST***>(),m.get_view<int***>(),mh,value,dims);
      break;
    case 4:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST****>(),m.get_view<int****>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST****>(),m.get_view<int****>(),mh,value,dims);
      break;
    case 5:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST*****>(),m.get_view<int*****>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST*****>(),m.get_view<int*****>(),mh,value,dims);
      break;
    case 6:
      if (contiguous)
        setMaskHelper<CMP>(x.get_view<const ST******>(),m.get_view<int******>(),mh,value,dims);
      else
        setMaskHelper<CMP>(x.get_strided_view<const ST******>(),m.get_view<int******>(),mh,value,dims);
      break;
    default:
      EKAT_ERROR_MSG ("Unsupported field rank in compute_mask.\n"
          " - field name: " << x.name() << "\n"
          " - field rank: " << x.rank() << "\n");
  }
}

#define EAMXX_FIELD_UTILS_ETI_DECL_COMPUTE_MASK(C) \
extern template void impl::compute_mask<C, int>(const Field&, const int, Field&); \
extern template void impl::compute_mask<C, float>(const Field&, const float, Field&); \
extern template void impl::compute_mask<C, double>(const Field&, const double, Field&)

// Only these two, since they are the only ones used so far (in IO)
EAMXX_FIELD_UTILS_ETI_DECL_COMPUTE_MASK(Comparison::NE);
EAMXX_FIELD_UTILS_ETI_DECL_COMPUTE_MASK(Comparison::LE);

} // namespace impl

} // namespace scream

#endif // SCREAM_FIELD_UTILS_IMPL_HPP
