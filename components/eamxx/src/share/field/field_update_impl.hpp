#ifndef SCREAM_FIELD_UPDATE_IMPL_HPP
#define SCREAM_FIELD_UPDATE_IMPL_HPP

#include "share/field/field.hpp"

#include "share/field/field_impl_details.hpp"

namespace scream
{

template<CombineMode CM, bool FillAware, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_contig = x.get_header().get_alloc_properties().contiguous();
  const auto y_contig = get_header().get_alloc_properties().contiguous();
  switch (layout.rank()) {
    case 0:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      break;
    case 1:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      break;
    case 2:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      break;
    case 3:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      break;
    case 4:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      break;
    case 5:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      break;
    case 6:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<bool use_mask, typename ST>
void Field::deep_copy_impl (const ST value, const Field& mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto contig = get_header().get_alloc_properties().contiguous();

  switch (rank) {
    case 0:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims,
                                 mask.get_view<const int>());
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims,
                                 mask.get_view<const int>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims);
      }
      break;
    case 1:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims);
      }
      break;
    case 2:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims);
      }
      break;
    case 3:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims);
      }
      break;
    case 4:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims);
      }
      break;
    case 5:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims);
      }
      break;
    case 6:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims);
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
