#ifndef EAMXX_DIAGNOSTIC_TEST_UTIL_HPP
#define EAMXX_DIAGNOSTIC_TEST_UTIL_HPP

#include "share/grid/mesh_free_grids_manager.hpp"
namespace scream{

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_local_elems = 4;
  const int np = 4;
  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<std::string>("Reference Grid", "Point Grid");
  gm_params.set<int>("Number of Global Columns", num_global_cols);
  gm_params.set<int>("Number of Local Elements", num_local_elems);
  gm_params.set<int>("Number of Vertical Levels", nlevs);
  gm_params.set<int>("Number of Gauss Points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_all_grids();

  return gm;
}

template<typename ScalarT, int NumLevels>
struct ChecksHelpers {

  static bool is_non_negative (const ScalarT& s, const int k) {
    return not ( k<NumLevels && (s<0 || std::isnan(s)) );
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return lhs==rhs;
  }
  static bool approx_equal (const ScalarT lhs, const ScalarT rhs,
                            const int k, const ScalarT tol) {
    using std::abs;
    return not ( k<NumLevels && abs(lhs-rhs)>=tol );
  }
  static bool approx_equal (const ScalarT computed, const ScalarT expected, const ScalarT tol) {
    using std::abs;
    return abs(computed-expected)/abs(expected) < tol;
  }
};

template<typename T, int N, int NumLevels>
struct ChecksHelpers<ekat::Pack<T,N>,NumLevels> {
  using ScalarT = ekat::Pack<T,N>;

  static bool is_non_negative (const ScalarT& s, const int k) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return ( range_mask && (s<0 || isnan(s) ) ).none();
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return (lhs==rhs).all();
  }
  static bool approx_equal (const ScalarT& lhs, const ScalarT& rhs,
                            const int k, const T tol) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return (range_mask && abs(lhs-rhs)>=tol).none();
  }
  static bool approx_equal (const ScalarT& computed, const ScalarT& expected, const T tol) {
    return (abs(computed-expected)/abs(expected) < tol).all();
  }
};

// Helper function. Create Mirror View and Deep-Copy (CMVDC)
template<typename ViewT>
auto cmvdc (const ViewT& v_d) -> typename ViewT::HostMirror {
  auto v_h = Kokkos::create_mirror_view(v_d);
  Kokkos::deep_copy(v_h,v_d);
  return v_h;
}

} // namespace scream

#endif // EAMXX_DIAGNOSTIC_TEST_UTIL_HPP
