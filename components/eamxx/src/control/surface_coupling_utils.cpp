#include "control/surface_coupling_utils.hpp"

namespace scream {

void get_col_info_for_surface_values(const std::shared_ptr<const FieldHeader>& fh,
                                     int vecComp, int& col_offset, int& col_stride)
{
  using namespace ShortFieldTagsNames;

  // Get this field layout info
  const auto& layout = fh->get_identifier().get_layout();
  const auto& dims = layout.dims();

  auto lt = get_layout_type(layout.tags());
  const bool scalar   = lt==LayoutType::Scalar2D || lt==LayoutType::Scalar3D;
  const bool vector   = lt==LayoutType::Vector2D || lt==LayoutType::Vector3D;
  const bool layout3d = lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D;
  EKAT_REQUIRE_MSG(scalar || vector,
      "Error! Support for tensor fields not yet implemented in surface coupling.\n");
  EKAT_REQUIRE_MSG( (vecComp<0) == scalar,
      "Error! You can and must specify a vector component only for vector fields.\n");

  vecComp = std::max(0,vecComp);

  // Compute initial offset. We use get_last_extent() as it
  // accounts for padding (if any).
  col_offset = vecComp;
  if (layout3d) {
    if (lt==LayoutType::Vector3D) {
      col_offset *= fh->get_alloc_properties().get_last_extent();
    }
    col_offset += dims.back()-1;
  }

  // If rank>1, there always is some stride
  col_stride = 1;
  if (layout.rank()>1) {
    col_stride = fh->get_alloc_properties().get_last_extent();
    if (layout.rank()>2) {
      col_stride *= dims[1];
    }
  }

  // If this field has a parent, then the underlying data includes all entires
  // in the larger field. We then must treat the child field as having col_stride
  // and col_offset of a vector field with component as its subview_idx.
  std::shared_ptr<const FieldHeader> parent = fh->get_parent().lock();
  if (parent != nullptr) {

    EKAT_REQUIRE_MSG(parent->get_parent().lock() == nullptr,
                     "Error! Currently support isn't added for fields with grandparents.\n");

    const auto parent_lt = get_layout_type(parent->get_identifier().get_layout().tags());

    EKAT_REQUIRE_MSG(parent_lt==LayoutType::Vector3D,
                     "Error! SurfaceCoupling expects all subfields to have parents "
                     "with LayoutType::Vector3D.\n");

    const auto& sv_info = fh->get_alloc_properties().get_subview_info();

    // Recall: idx = (idim,k) = (dimension where slice happened, index along said dimension).
    // Field class only allows idim=0,1. But we should never be in the case of idim=0, here.
    // If we have idim=0, it means that the parent field did not have COL as tag[0].
    EKAT_REQUIRE_MSG(sv_info.dim_idx==1, "Error! Bizarre scenario discovered. Contact developers.\n");

    // Additional col_offset
    col_offset += sv_info.slice_idx*parent->get_alloc_properties().get_last_extent();

    // Additional product for col_stride
    col_stride *= parent->get_identifier().get_layout().dim(1);
  }
}

} // namespace scream
