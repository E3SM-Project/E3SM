#ifndef SCREAM_DEEP_COPY_HPP
#define SCREAM_DEEP_COPY_HPP

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_assert.hpp"
#include <vector>

namespace scream {

/*
 * This class provides functions which deep_copy between host and
 * device for 1d scalar views.
 */

struct ScreamDeepCopy {

  template<typename ViewT>
  static void copy_to_device(const std::vector<typename ViewT::value_type const*>& data,
                      const std::vector<int>&                               sizes,
                      std::vector<ViewT>&                                   views)
  {
    EKAT_ASSERT(ViewT::rank == 1);
    EKAT_ASSERT(data.size() == views.size());
    EKAT_ASSERT(data.size() == sizes.size());

    const int num_views = data.size();
    for (int n=0; n<num_views; ++n) {
      const int size = sizes[n];
      views[n] = ViewT("", size);
      const auto host_view = Kokkos::create_mirror_view(views[n]);
      for (int i=0; i<size; ++i) {
        host_view(i) = data[n][i];
      }
      Kokkos::deep_copy(views[n],host_view);
    }
  }

  template<typename ViewT>
  static void copy_to_device(const std::vector<typename ViewT::value_type const*>& data,
                      const int&                                            size,
                      std::vector<ViewT>&                                   views)
  {
    std::vector<int> sizes(data.size(),size);
    copy_to_device(data,sizes,views);
  }

  template<typename ViewT>
  static void copy_to_host(const std::vector<typename ViewT::value_type*>& data,
                    const std::vector<int>&                         sizes,
                    std::vector<ViewT>&                             views)
  {
    EKAT_ASSERT(ViewT::rank == 1);
    EKAT_ASSERT(data.size() == views.size());
    EKAT_ASSERT(data.size() == sizes.size());

    const int num_views = data.size();
    for (int n=0; n<num_views; ++n) {
      const int size = sizes[n];
      const auto host_view = Kokkos::create_mirror_view(views[n]);
      Kokkos::deep_copy(host_view,views[n]);
      for (int i=0; i<size; ++i) {
        data[n][i] = host_view(i);
      }
    }
  }

  template<typename ViewT>
  static void copy_to_host(const std::vector<typename ViewT::value_type*>& data,
                    const int&                                      size,
                    std::vector<ViewT>&                             views)
  {
    std::vector<int> sizes(data.size(),size);
    copy_to_host(data,sizes,views);
  }

};

} // namespace screa

#endif // SCREAM_DEEP_COPY_HPP
