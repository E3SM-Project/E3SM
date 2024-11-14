#ifndef SCREAM_FIELD_UTILS_IMPL_HPP
#define SCREAM_FIELD_UTILS_IMPL_HPP

#include "share/field/field.hpp"

#include "ekat/mpi/ekat_comm.hpp"

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
void perturb (const Field& f,
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

  if (has_column_dim) {
    // Because Column is the partitioned dimension, we must reset the
    // RNG seed to be the same on every column so that a column will
    // have the same value no matter where it exists in an MPI rank's
    // set of local columns.
    const auto gids = dof_gids.get_strided_view<const int*, Host>();

    // Create a field to store perturbation values with layout
    // the same as f, but stripped of column and level dimension.
    auto perturb_fl = fl.clone().strip_dim(COL).strip_dim(LEV);
    FieldIdentifier perturb_fid("perturb_field", perturb_fl, ekat::units::Units::nondimensional(), "");
    Field perturb_f(perturb_fid);
    perturb_f.allocate_view();

    // Loop through columns as reset RNG seed based on GID of column
    for (auto icol=0; icol<fl.dims().front(); ++icol) {
      const auto new_seed = base_seed+gids(icol);
      engine.seed(new_seed);

      // Loop through levels. For each that satisfy the level_mask,
      // apply a random perturbation to f.
      for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
        if (level_mask(ilev)) {
          randomize(perturb_f, engine, pdf);
          f.subfield(0, icol).subfield(f.rank()-2, ilev).scale(perturb_f);
        }
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

    // Loop through levels. For each that satisfy the level_mask,
    // apply a random perturbation to f.
    for (auto ilev=0; ilev<fl.dims().back(); ++ilev) {
      if (level_mask(ilev)) {
        randomize(perturb_f, engine, pdf);
        f.subfield(f.rank()-1, ilev).scale(perturb_f);
      }
    }
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

template<typename T>
void print_field_hyperslab (const Field& f,
                            std::vector<FieldTag> tags,
                            std::vector<int> indices,
                            std::ostream& out,
                            const int orig_rank,
                            const size_t curr_idx)
{
  // General idea: call f.subfield with the proper index, and recurse
  // until all indices are exausted, then print the field that is left.
  //
  // We keep all the tags/indices we slice away, since we need them at
  // the end of recursion when we print the info of the field location.
  // E.g., if f has tags/dims <COL,CMP,LEV>/(2,3,4) and we subview at
  // tags/indices <COL,LEV>/(0,1), we want to print something like
  //   f(0,:,1):
  //     0.123, 0.456, 0.789

  EKAT_REQUIRE_MSG (tags.size()==indices.size(),
      "Error! Tags vector size differs from indices vector size.\n");

  const auto& layout = f.get_header().get_identifier().get_layout();

  // Get layout of original field (before all the slicing happened)
  auto get_orig_header = [&]() -> std::shared_ptr<const FieldHeader> {
    auto fh = f.get_header_ptr();
    while (fh->get_identifier().get_layout().rank()<orig_rank) {
      fh = fh->get_parent().lock();
    }
    return fh;
  };

  // Get a vector of strings, which we'll use to print the info of the field
  // slice/entry currently printed. We can already fill the entries corresponding
  // to dimensions we sliced away.
  auto get_dims_str = [&](const FieldLayout& orig_layout) -> std::vector<std::string> {
    const int orig_rank = orig_layout.rank();
    const int num_tags  = tags.size();
    std::vector<std::string> dims_str (orig_rank,"");
    for (int ii=0, jj=0; ii<orig_rank && jj<num_tags; ++ii) {
      if (orig_layout.tag(ii)==tags[jj]) {
        // Was sliced. Store slice idx as a stirng
        dims_str[ii] = std::to_string(indices[jj]);
        ++jj;
      }
    }
    return dims_str;
  };

  // Get the dimensions (in the original field) that were *NOT* sliced away,
  // since we'll have to loop over them. We'll use these indices to add the
  // missing info in the dims_str vector<string> as we loop and print.
  auto get_dims_left = [&](const FieldLayout& orig_layout) -> std::vector<int> {
    const int orig_rank = orig_layout.rank();
    const int num_tags  = tags.size();
    std::vector<int> dims_left;
    for (int ii=0, jj=0; ii<orig_rank; ++ii) {
      if (jj<num_tags && orig_layout.tag(ii)==tags[jj]) {
        // Was sliced. Skip
        ++jj;
      } else {
        // Was not sliced.
        dims_left.push_back(ii);
      }
    }
    return dims_left;
  };

  constexpr int max_per_line = 5;
  if (curr_idx==tags.size()) {
    // All slices have been taken. Print the whole input field.
    const auto& orig_layout = get_orig_header()->get_identifier().get_layout();
    const auto dims_left = get_dims_left(orig_layout);
    auto dims_str = get_dims_str(orig_layout);

    // NOTE: because print_field_hyperslab() is only used for testing, we
    // generalize by always calling get_strided_view(), even if it could be
    // contiguous.
    // This allows us to call this function on any type of field,
    // including multi-slice subfields
    f.sync_to_host();
    const int rank = layout.rank();
    out << "     " << f.name() << orig_layout.to_string() << "\n\n";
    switch (rank) {
      case 0:
      {
        out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
        // NOTE: add ", " at the end, to make rank0 behave the same as other ranks,
        //       for the sake of any script trying to manipulate output
        out << "\n    " << f.get_strided_view<const T,Host>()() << ", \n";
        break;
      }
      case 1:
      {
        dims_str[dims_left[0]] = ":";
        out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
        auto v = f.get_strided_view<const T*,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          if (i%max_per_line==0) {
            out << "\n    ";
          }
          out << v(i) << ", ";
        }
        out << "\n";
        break;
      }
      case 2:
      {
        dims_str[dims_left[1]] = ":";
        auto v = f.get_strided_view<const T**,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
          for (int j=0; j<layout.dim(1); ++j) {
            if (j%max_per_line==0) {
              out << "\n    ";
            }
            out << v(i,j) << ", ";
          }
          out << "\n";
        }
        break;
      }
      case 3:
      {
        dims_str[dims_left[2]] = ":";
        auto v = f.get_strided_view<const T***,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          for (int j=0; j<layout.dim(1); ++j) {
            dims_str[dims_left[1]] = std::to_string(j);
            out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
            for (int k=0; k<layout.dim(2); ++k) {
              if (k%max_per_line==0) {
                out << "\n    ";
              }
              out << v(i,j,k) << ", ";
            }
            out << "\n";
          }
        }
        break;
      }
      case 4:
      {
        dims_str[dims_left[3]] = ":";
        auto v = f.get_strided_view<const T****,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          for (int j=0; j<layout.dim(1); ++j) {
            dims_str[dims_left[1]] = std::to_string(j);
            for (int k=0; k<layout.dim(2); ++k) {
              dims_str[dims_left[2]] = std::to_string(k);
              out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
              for (int l=0; l<layout.dim(3); ++l) {
                if (l%max_per_line==0) {
                  out << "\n    ";
                }
                out << v(i,j,k,l) << ", ";
              }
              out << "\n";
            }
          }
        }
        break;
      }
      default:
        EKAT_ERROR_MSG (
            "Unsupported rank in print_field_hyperslab.\n"
            "  - field name  : " + f.name() + "\n"
            "  - field layout (upon slicing): " + layout.to_string() + "\n");
    }
  } else {
    auto tag = tags[curr_idx];
    auto idx = indices[curr_idx];

    auto it = ekat::find(layout.tags(),tag);
    EKAT_REQUIRE_MSG (it!=layout.tags().end(),
        "Error! Something went wrong while slicing field.\n"
        "  - field name  : " + f.name() + "\n"
        "  - field layout: " + layout.to_string() + "\n"
        "  - curr tag    : " + e2str(tag) + "\n");
    auto idim = std::distance(layout.tags().begin(),it);

    EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Cannot subview field for printing.\n"
        "  - field name  : " + f.name() + "\n"
        "  - field layout: " + layout.to_string() + "\n"
        "  - loc tags    : <" + ekat::join(tags,",") + ">\n"
        "  - loc indices : (" + ekat::join(indices,",") + ")\n");

    auto sub_f = f.subfield(idim,idx);
    return print_field_hyperslab<T>(sub_f,tags,indices,out,orig_rank,curr_idx+1);
  }
}

} // namespace impl

} // namespace scream

#endif // SCREAM_FIELD_UTILS_IMPL_HPP
