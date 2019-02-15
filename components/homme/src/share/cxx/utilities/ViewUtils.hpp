/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VIEW_UTILS_HPP
#define HOMMEXX_VIEW_UTILS_HPP

#include "Types.hpp"

namespace Homme {

// Structure to define the type of a view that has a const data type,
// given the type of an input view
template<typename ViewT>
struct ViewConst{};

template<typename DataType, typename... Properties>
struct ViewConst<ViewType<DataType,Properties...>> {
  using type = ViewType<const DataType,Properties...>;
};

template<typename ViewT>
typename ViewConst<ViewT>::type
viewConst(const ViewT&& v) {
  return typename ViewConst<ViewT>::type(v);
}

template<typename ViewT>
typename ViewConst<ViewT>::type
viewConst(const ViewT& v) {
  return typename ViewConst<ViewT>::type(v);
}

} // namespace Homme

#endif // HOMMEXX_VIEW_UTILS_HPP
