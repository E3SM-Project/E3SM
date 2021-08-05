#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "field.hpp"

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
template<typename RT1, typename RT2>
bool views_are_equal(const Field<RT1>& f1, const Field<RT2>& f2) {
  static_assert(
      std::is_same<typename std::remove_cv<RT1>::type,
                   typename std::remove_cv<RT2>::type>::value,
      "Error! Real types must be the same (except possibly for cv qualifiers).\n");

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
  const auto& dims = l1.dims();
  switch (l1.rank()) {
    case 1:
      {
        auto v1 = f1.template get_view<RT1*,Host>();
        auto v2 = f2.template get_view<RT2*,Host>();
        for (int i=0; i<dims[0]; ++i) {
          if (v1(i) != v2(i)) {
            return false;
          }
        }
      }
      break;
    case 2:
      {
        auto v1 = f1.template get_view<RT1**,Host>();
        auto v2 = f2.template get_view<RT2**,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            if (v1(i,j) != v2(i,j)) {
              return false;
            }
        }}
      }
      break;
    case 3:
      {
        auto v1 = f1.template get_view<RT1***,Host>();
        auto v2 = f2.template get_view<RT2***,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              if (v1(i,j,k) != v2(i,j,k)) {
                return false;
              }
        }}}
      }
      break;
    case 4:
      {
        auto v1 = f1.template get_view<RT1****,Host>();
        auto v2 = f2.template get_view<RT2****,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                if (v1(i,j,k,l) != v2(i,j,k,l)) {
                  return false;
                }
        }}}}
      }
      break;
    case 5:
      {
        auto v1 = f1.template get_view<RT1*****,Host>();
        auto v2 = f2.template get_view<RT2*****,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                for (int m=0; m<dims[4]; ++m) {
                  if (v1(i,j,k,l,m) != v2(i,j,k,l,m)) {
                    return false;
                  }
        }}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  // If we get here, then all entries matched.
  return true;
}

template<typename RT, typename Engine, typename PDF>
void randomize (const Field<RT>& f, Engine& engine, PDF&& pdf)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  const auto& fl = f.get_header().get_identifier().get_layout();
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_view<RT*,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          v(i) = pdf(engine);
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<RT**,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            v(i,j) = pdf(engine);
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<RT***,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              v(i,j,k) = pdf(engine);
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<RT****,Host>();
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
        auto v = f.template get_view<RT*****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  v(i,j,k,l,m) = pdf(engine);
        }}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  // Sync the dev view with the host view.
  f.sync_to_dev();
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
