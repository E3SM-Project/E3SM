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
  bool same_locally = true;
  const auto& dims = l1.dims();
  switch (l1.rank()) {
    case 1:
      {
        auto v1 = f1.template get_view<ST*,Host>();
        auto v2 = f2.template get_view<ST*,Host>();
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
        auto v1 = f1.template get_view<ST**,Host>();
        auto v2 = f2.template get_view<ST**,Host>();
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
        auto v1 = f1.template get_view<ST***,Host>();
        auto v2 = f2.template get_view<ST***,Host>();
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
        auto v1 = f1.template get_view<ST****,Host>();
        auto v2 = f2.template get_view<ST****,Host>();
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
        auto v1 = f1.template get_view<ST*****,Host>();
        auto v2 = f2.template get_view<ST*****,Host>();
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
        auto v1 = f1.template get_view<ST******,Host>();
        auto v2 = f2.template get_view<ST******,Host>();
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
    case 1:
      {
        auto v = f.template get_view<ST*,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          v(i) = pdf(engine);
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<ST**,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            v(i,j) = pdf(engine);
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<ST***,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              v(i,j,k) = pdf(engine);
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<ST****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                v(i,j,k,l) = pdf(engine);
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_view<ST*****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  v(i,j,k,l,m) = pdf(engine);
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_view<ST******,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  for (int n=0; n<v.extent_int(5); ++n) {
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
        auto v = f.template get_view<ST*,Host>();
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
        auto v = f.template get_view<ST**,Host>();
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
        auto v = f.template get_view<ST***,Host>();
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
        auto v = f.template get_view<ST****,Host>();
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
        auto v = f.template get_view<ST*****,Host>();
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
        auto v = f.template get_view<ST******,Host>();
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
        auto v = f.template get_view<ST*,Host>();
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
        auto v = f.template get_view<ST**,Host>();
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
        auto v = f.template get_view<ST***,Host>();
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
        auto v = f.template get_view<ST****,Host>();
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
        auto v = f.template get_view<ST*****,Host>();
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
        auto v = f.template get_view<ST******,Host>();
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
        auto v = f.template get_view<ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          max = std::max(max,v(i));
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            max = std::max(max,v(i,j));
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              max = std::max(max,v(i,j,k));
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<ST****,Host>();
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
        auto v = f.template get_view<ST*****,Host>();
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
        auto v = f.template get_view<ST******,Host>();
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
        auto v = f.template get_view<ST*,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          min = std::min(min,v(i));
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<ST**,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            min = std::min(min,v(i,j));
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<ST***,Host>();
        for (int i=0; i<fl.dim(0); ++i) {
          for (int j=0; j<fl.dim(1); ++j) {
            for (int k=0; k<fl.dim(2); ++k) {
              min = std::min(min,v(i,j,k));
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<ST****,Host>();
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
        auto v = f.template get_view<ST*****,Host>();
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
        auto v = f.template get_view<ST******,Host>();
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

} // namespace impl

} // namespace scream

#endif // SCREAM_FIELD_UTILS_IMPL_HPP
