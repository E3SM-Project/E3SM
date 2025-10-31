#include "share/field/field_utils.hpp"

#include <ekat_comm.hpp>

namespace scream {

namespace impl {

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

} // namespace impl

ScalarWrapper frobenius_norm(const Field& f, const ekat::Comm* comm)
{
  EKAT_REQUIRE_MSG (f.is_allocated(),
    "[frobenius_norm] Error! Input field was not yet allocated.\n");

  ScalarWrapper norm;
  switch (f.data_type()) {
    case DataType::IntType:
      norm.set(impl::frobenius_norm<int>(f,comm));
      break;
    case DataType::FloatType:
      norm.set(impl::frobenius_norm<float>(f,comm));
      break;
    case DataType::DoubleType:
      norm.set(impl::frobenius_norm<double>(f,comm));
      break;
    default:
      EKAT_ERROR_MSG ("[print_field] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
  return norm;
}

ScalarWrapper field_sum(const Field& f, const ekat::Comm* comm)
{
  EKAT_REQUIRE_MSG (f.is_allocated(),
    "[field_sum] Error! Input field was not yet allocated.\n");

  ScalarWrapper norm;
  switch (f.data_type()) {
    case DataType::IntType:
      norm.set(impl::field_sum<int>(f,comm));
      break;
    case DataType::FloatType:
      norm.set(impl::field_sum<float>(f,comm));
      break;
    case DataType::DoubleType:
      norm.set(impl::field_sum<double>(f,comm));
      break;
    default:
      EKAT_ERROR_MSG ("[print_field] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
  return norm;
}

ScalarWrapper field_max(const Field& f, const ekat::Comm* comm)
{
  EKAT_REQUIRE_MSG (f.is_allocated(),
    "[field_max] Error! Input field was not yet allocated.\n");

  ScalarWrapper norm;
  switch (f.data_type()) {
    case DataType::IntType:
      norm.set(impl::field_max<int>(f,comm));
      break;
    case DataType::FloatType:
      norm.set(impl::field_max<float>(f,comm));
      break;
    case DataType::DoubleType:
      norm.set(impl::field_max<double>(f,comm));
      break;
    default:
      EKAT_ERROR_MSG ("[print_field] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
  return norm;
}

ScalarWrapper field_min(const Field& f, const ekat::Comm* comm)
{
  EKAT_REQUIRE_MSG (f.is_allocated(),
    "[field_min] Error! Input field was not yet allocated.\n");

  ScalarWrapper norm;
  switch (f.data_type()) {
    case DataType::IntType:
      norm.set(impl::field_min<int>(f,comm));
      break;
    case DataType::FloatType:
      norm.set(impl::field_min<float>(f,comm));
      break;
    case DataType::DoubleType:
      norm.set(impl::field_min<double>(f,comm));
      break;
    default:
      EKAT_ERROR_MSG ("[print_field] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
  return norm;
}

} // namespace scream
