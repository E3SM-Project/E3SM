#ifndef SCREAM_ARRAY_UTILS_HPP
#define SCREAM_ARRAY_UTILS_HPP

#include "share/scream_types.hpp"

#include <ekat/ekat_assert.hpp>

namespace scream {

// Given array dims, and a 1d index of flattened array,
// return vector<int> of the indices of idx in the original Nd array.
// NOTE: right-most dims are assumed to stride faster.
inline std::vector<int> unflatten_idx (const std::vector<int>& dims, const int idx) {
  const int r = dims.size();

  // Get dims from the back fastest to slowest
  const auto rbeg = dims.rbegin();
  const int dr0 = *rbeg;
  const int dr1 = r>1 ? *(rbeg+1) : 0;
  const int dr2 = r>2 ? *(rbeg+2) : 0;
  const int dr3 = r>3 ? *(rbeg+3) : 0;
  const int dr4 = r>4 ? *(rbeg+4) : 0;

  std::vector<int> indices(r);

  // Access the indices array backwards, so that ind(0) is the fastest
  // striding index, and ind(r-1) is the slowest
  auto ind = [&] (const int i) -> int& {
    return *(indices.rbegin()+i);
  };

  switch (r) {
    case 6: ind(5) = ((((idx/dr0)/dr1)/dr2)/dr3)/dr4; // Fallthrough
    case 5: ind(4) = ((((idx/dr0)/dr1)/dr2)/dr3)%dr4; // Fallthrough
    case 4: ind(3) =  (((idx/dr0)/dr1)/dr2)%dr3;      // Fallthrough
    case 3: ind(2) =   ((idx/dr0)/dr1)%dr2;           // Fallthrough
    case 2: ind(1) =    (idx/dr0)%dr1;                // Fallthrough
    case 1: ind(0) =     idx%dr0; break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank (" + std::to_string(r) + ").\n");
  }

  return indices;
}

// Kokkos-friendly versions of the above function, taking a Kokkos::View
// for the array dimensions, and unpacking directly into N integers.
template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==1,"Error! Wrong overload of unflatten_idx called.\n");
  i = idx;
  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i, int& j)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==2,"Error! Wrong overload of unflatten_idx called.\n");

  i = idx / dims[1];
  j = idx % dims[1];

  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (j>=0 && j<dims[1],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i, int& j, int& k)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==3,"Error! Wrong overload of unflatten_idx called.\n");

  i = (idx / dims[2]) / dims[1];
  j = (idx / dims[2]) % dims[1];
  k =  idx % dims[2];

  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (j>=0 && j<dims[1],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (k>=0 && k<dims[2],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i, int& j, int& k, int& l)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==4,"Error! Wrong overload of unflatten_idx called.\n");

  i = ((idx / dims[3]) / dims[2]) / dims[1];
  j = ((idx / dims[3]) / dims[2]) % dims[1];
  k =  (idx / dims[3]) % dims[2];
  l =   idx % dims[3];

  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (j>=0 && j<dims[1],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (k>=0 && k<dims[2],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (l>=0 && l<dims[3],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i, int& j, int& k, int& l, int& m)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==5,"Error! Wrong overload of unflatten_idx called.\n");

  i = (((idx / dims[4]) / dims[3]) / dims[2]) / dims[1];
  j = (((idx / dims[4]) / dims[3]) / dims[2]) % dims[1];
  k =  ((idx / dims[4]) / dims[3]) % dims[2];
  l =   (idx / dims[4]) % dims[3];
  m =    idx % dims[4];

  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (j>=0 && j<dims[1],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (k>=0 && k<dims[2],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (l>=0 && l<dims[3],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (m>=0 && m<dims[4],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims,
    int& i, int& j, int& k, int& l, int& m, int& n)
{
  EKAT_KERNEL_ASSERT_MSG (dims.size()==6,"Error! Wrong overload of unflatten_idx called.\n");

  i = ((((idx / dims[5]) / dims[4]) / dims[3]) / dims[2]) / dims[1];
  j = ((((idx / dims[5]) / dims[4]) / dims[3]) / dims[2]) % dims[1];
  k =  (((idx / dims[5]) / dims[4]) / dims[3]) % dims[2];
  l =   ((idx / dims[5]) / dims[4]) % dims[3];
  m =    (idx / dims[5]) % dims[4];
  n =     idx % dims[5];

  EKAT_KERNEL_ASSERT_MSG (i>=0 && i<dims[0],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (j>=0 && j<dims[1],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (k>=0 && k<dims[2],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (l>=0 && l<dims[3],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (m>=0 && m<dims[4],"Error! Flatten index out of bounds.\n");
  EKAT_KERNEL_ASSERT_MSG (n>=0 && n<dims[5],"Error! Flatten index out of bounds.\n");
}

template<typename... Props>
KOKKOS_INLINE_FUNCTION
void unflatten_idx (const int idx, const Kokkos::View<int*,Kokkos::LayoutRight,Props...>& dims, int* indices) {
  const int rank = dims.size();
  switch(rank) {
    case 1:
      unflatten_idx(idx,dims,indices[0]);
      break;
    case 2:
      unflatten_idx(idx,dims,indices[0],indices[1]);
      break;
    case 3:
      unflatten_idx(idx,dims,indices[0],indices[1],indices[2]);
      break;
    case 4:
      unflatten_idx(idx,dims,indices[0],indices[1],indices[2],indices[3]);
      break;
    case 5:
      unflatten_idx(idx,dims,indices[0],indices[1],indices[2],indices[3],indices[4]);
      break;
    case 6:
      unflatten_idx(idx,dims,indices[0],indices[1],indices[2],indices[3],indices[4],indices[5]);
      break;
    default:
      EKAT_KERNEL_ERROR_MSG ("Error! Unsupported rank.\n");
  }
}

} // namespace scream

#endif // SCREAM_ARRAY_UTILS_HPP
