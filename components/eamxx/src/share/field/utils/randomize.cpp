#include "share/field/field_utils.hpp"

#include <random>

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
namespace impl {

template<typename Engine, typename PDF>
void randomize (const Field& f, Engine& engine, PDF&& pdf)
{
  // Deduce scalar type from pdf
  using ST = decltype(pdf(engine));

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

} // namespace impl

void randomize_normal (const Field& f, int seed, const ScalarWrapper mean, const ScalarWrapper sigma)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  std::mt19937_64 engine(seed);

  // Check compatibility between PDF and field data type
  const auto dt = f.data_type();

  switch (dt) {
    case DataType::FloatType:
      {
        std::normal_distribution pdf(mean.as<float>(),sigma.as<float>());
        impl::randomize(f,engine,pdf);
        break;
      }
    case DataType::DoubleType:
      {
        std::normal_distribution pdf(mean.as<double>(),sigma.as<double>());
        impl::randomize(f,engine,pdf);
        break;
      }
    default:
      EKAT_ERROR_MSG ("[randomize_normal] Error! Unsupported field data type.\n");
  }
}

void randomize_uniform(const Field& f, int seed, const ScalarWrapper lb, const ScalarWrapper ub)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  std::mt19937_64 engine(seed);

  // Check compatibility between PDF and field data type
  const auto dt = f.data_type();
  switch (dt) {
    case DataType::IntType:
      {
        std::uniform_int_distribution pdf(lb.as<int>(),ub.as<int>());
        impl::randomize(f,engine,pdf);
        break;
      }
    case DataType::FloatType:
      {
        std::uniform_real_distribution pdf(lb.as<float>(),ub.as<float>());
        impl::randomize(f,engine,pdf);
        break;
      }
    case DataType::DoubleType:
      {
        std::uniform_real_distribution pdf(lb.as<double>(),ub.as<double>());
        impl::randomize(f,engine,pdf);
        break;
      }
    default:
      EKAT_ERROR_MSG ("[randomize] Error! Unsupported field data type.\n");
  }
}

void randomize_discrete(const Field& f, int seed,
                        const std::vector<ScalarWrapper>& values)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  std::mt19937_64 engine(seed);
  std::uniform_int_distribution<int> idx_pdf (0,values.size()-1);

  // Check compatibility between PDF and field data type
  const auto dt = f.data_type();
  switch (dt) {
    case DataType::IntType:
      {
        auto pdf = [&](std::mt19937_64& engine) -> int {
          int idx = idx_pdf(engine);
          return values[idx].as<int>();
        };
        impl::randomize(f,engine,pdf);
        break;
      }
    case DataType::FloatType:
      {
        auto pdf = [&](std::mt19937_64& engine) -> float {
          int idx = idx_pdf(engine);
          return values[idx].as<float>();
        };
        impl::randomize(f,engine,pdf);
        break;
      }
    case DataType::DoubleType:
      {
        auto pdf = [&](std::mt19937_64& engine) -> double {
          int idx = idx_pdf(engine);
          return values[idx].as<double>();
        };
        impl::randomize(f,engine,pdf);
        break;
      }
    default:
      EKAT_ERROR_MSG ("[randomize] Error! Unsupported field data type.\n");
  }
}


} // namespace scream
