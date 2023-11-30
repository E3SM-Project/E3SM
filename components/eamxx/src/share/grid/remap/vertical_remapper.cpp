#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/scream_vertical_interpolation.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"
#include "share/util/scream_universal_constants.hpp"

#include "ekat/util/ekat_units.hpp"
#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <ekat/ekat_pack_kokkos.hpp>

#include <numeric>

namespace scream
{
VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file,
                  const Field& lev_prof,
                  const Field& ilev_prof)
  : VerticalRemapper(src_grid,map_file,lev_prof,ilev_prof,constants::DefaultFillValue<float>::value)
{
  // Nothing to do here
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file,
                  const Field& lev_prof,
                  const Field& ilev_prof,
                  const Real mask_val)
 : AbstractRemapper()
 , m_comm (src_grid->get_comm())
 , m_mask_val(mask_val)
{
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (src_grid->type()==GridType::Point,
      "Error! VerticalRemapper only works on PointGrid grids.\n"
      "  - src grid name: " + src_grid->name() + "\n"
      "  - src_grid_type: " + e2str(src_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (src_grid->is_unique(),
      "Error! VerticalRemapper requires a unique source grid.\n");

  // This is a vertical remapper. We only go in one direction
  m_bwd_allowed = false;

  // Create tgt_grid that is a clone of the src grid but with
  // the correct number of levels.  Note that when vertically
  // remapping the target field will be defined on the same DOFs
  // as the source field, but will have a different number of 
  // vertical levels.
  scorpio::register_file(map_file,scorpio::FileMode::Read);
  m_num_remap_levs = scorpio::get_dimlen(map_file,"lev");

  auto tgt_grid = src_grid->clone("vertical_remap_tgt_grid",true);
  tgt_grid->reset_num_vertical_lev(m_num_remap_levs);
  this->set_grids(src_grid,tgt_grid);

  // Set the LEV and ILEV vertical profiles for interpolation from
  register_vertical_source_field(lev_prof);
  register_vertical_source_field(ilev_prof);

  // Gather the pressure level data for vertical remapping
  set_pressure_levels(map_file);

  // Add tgt pressure levels to the tgt grid
  tgt_grid->set_geometry_data(m_remap_pres);
  scorpio::eam_pio_closefile(map_file);
}

FieldLayout VerticalRemapper::
create_src_layout (const FieldLayout& tgt_layout) const
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[VerticalRemapper] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + to_string(tgt_layout));

  return create_layout(tgt_layout,m_src_grid);
}

FieldLayout VerticalRemapper::
create_tgt_layout (const FieldLayout& src_layout) const
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[VerticalRemapper] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + to_string(src_layout));

  return create_layout(src_layout,m_tgt_grid);
}

FieldLayout VerticalRemapper::
create_layout (const FieldLayout& fl_in,
               const grid_ptr_type& grid_out) const
{
  using namespace ShortFieldTagsNames;

  // NOTE: for the vert remapper, it doesn't really make sense to distinguish
  //       between midpoints and interfaces: we're simply asking for a quantity
  //       at a given set of pressure levels. So we choose to have fl_out
  //       to *always* have LEV as vertical tag.
  auto fl_out = FieldLayout::invalid();
  switch (fl_in.type()) {
    case LayoutType::Scalar0D: [[ fallthrough ]];
    case LayoutType::Vector0D: [[ fallthrough ]];
    case LayoutType::Scalar2D: [[ fallthrough ]];
    case LayoutType::Vector2D: [[ fallthrough ]];
    case LayoutType::Tensor2D:
      // These layouts do not have vertical dim tags, so no change
      fl_out = fl_in;
      break;
    case LayoutType::Scalar1D:
      fl_out = grid_out->get_vertical_layout(true);
      break;
    case LayoutType::Vector2D:
      tgt = m_tgt_grid->get_2d_vector_layout(fl_in.dim(CMP));
      break;
    case LayoutType::Scalar3D:
      fl_out = grid_out->get_3d_scalar_layout(true);
      break;
    case LayoutType::Vector3D:
      fl_out = grid_out->get_3d_vector_layout(true,vec_dim,fl_in.dim(CMP));
      break;
    default:
      // NOTE: this also include Tensor3D. We don't really have any atm proc
      //       that needs to handle a tensor3d quantity, so no need to add it
      EKAT_ERROR_MSG (
        "[VerticalRemapper] Error! Layout not supported by VerticalRemapper.\n"
        " - input layout: " + to_string(fl_in) + "\n");
  }
  return fl_out;
}

void VerticalRemapper::
set_pressure_levels(const std::string& map_file)
{
  // Ensure each map file gets a different decomp name
  static std::map<std::string,int> file2idx;
  if (file2idx.find(map_file)==file2idx.end()) {
    file2idx[map_file] = file2idx.size();
  }

  using namespace ShortFieldTagsNames;
  std::vector<FieldTag> tags = {LEV};
  std::vector<int>      dims = {m_num_remap_levs};
  FieldLayout layout(tags,dims);
  FieldIdentifier fid("p_levs",layout,ekat::units::Pa,m_tgt_grid->name());
  m_remap_pres = Field(fid);
  m_remap_pres.get_header().get_alloc_properties().request_allocation(mPack::n);
  m_remap_pres.allocate_view();

  auto remap_pres_scal = m_remap_pres.get_view<Real*,Host>();

  std::vector<scorpio::offset_t> dofs_offsets(m_num_remap_levs);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),0);
  const std::string decomp_tag = "VR::spl,nlev=" + std::to_string(m_num_remap_levs) + ",file-idx=" + std::to_string(file2idx[map_file]);
  scorpio::register_variable(map_file, "p_levs", "p_levs", {"lev"}, "real", decomp_tag);
  scorpio::set_dof(map_file,"p_levs",m_num_remap_levs,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"p_levs",-1,remap_pres_scal.data(),remap_pres_scal.size());

  m_remap_pres.sync_to_dev();
}

void VerticalRemapper::
register_vertical_source_field(const Field& src)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG(src.is_allocated(),
      "Error! Vertical level source field is not yet allocated.\n"
      " - field name: " + src.name() + "\n");

  const auto& layout = src.get_header().get_identifier().get_layout();
  const auto vert_tag = layout.tags().back();
  EKAT_REQUIRE_MSG (vert_tag==LEV or vert_tag==ILEV,
      "Error! Input vertical level field does not have a vertical level tag at the end.\n"
      " - field name: " + src.name() + "\n"
      " - field layout: " + to_string(layout) + "\n");

  if (vert_tag==LEV) {
    m_src_mid = src;
    m_mid_set = true; 
   } else {
    m_src_int = src;
    m_int_set = true; 
  }
}

void VerticalRemapper::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_src_fields.push_back(field_type(src));
  field_type tgt_f(tgt);
  m_tgt_fields.push_back(tgt_f);
}

void VerticalRemapper::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  using namespace ShortFieldTagsNames;
  auto name = src.name();
  auto src_layout = src.get_header().get_identifier().get_layout();
  auto tgt_layout = tgt.get_header().get_identifier().get_layout();
  const bool has_ilev = src_layout.has_tag(ILEV);
  EKAT_REQUIRE_MSG(src_layout.rank()==tgt_layout.rank(),
      "ERROR! vert_remap:do_bind_field:" + name + ", tgt and src do not have the same rank");
  // Note, for vertical remapper we set all target fields as having LEV as the vertical dimension.  So we check that all other tags
  // between source and target match if source has ILEV
  if (has_ilev) {
    EKAT_REQUIRE_MSG(src_layout.strip_dim(ILEV).tags()==tgt_layout.strip_dim(LEV).tags(),
        "ERROR! vert_remap:do_bind_field:" + name + ", tgt and src do not have the same set of field tags");
  } else {
    EKAT_REQUIRE_MSG(src_layout.tags()==tgt_layout.tags(),
        "ERROR! vert_remap:do_bind_field:" + name + ", tgt and src do not have the same set of field tags");
  } 

  EKAT_REQUIRE_MSG (
      src_layout.rank()>1 ||
      src.get_header().get_alloc_properties().get_padding()==0,
      "Error! vert_remap:do_bind_field:check_src:" + name + ", We don't support 2d scalar fields that are padded.\n");
  EKAT_REQUIRE_MSG (
      tgt_layout.rank()>1 ||
      tgt.get_header().get_alloc_properties().get_padding()==0,
      "Error! vert_remap:do_bind_field:check_tgt:" + name + ", We don't support 2d scalar fields that are padded.\n");

  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

  // Add mask tracking to the target field
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  auto nondim = Units::nondimensional();
  if (src_layout.has_tag(LEV) || src_layout.has_tag(ILEV)) {
    auto& f_tgt = m_tgt_fields[ifield];
    // NOTE: for now we assume that masking is determined only by the COL,LEV location in space
    //       and that fields with multiple components will have the same masking for each component
    //       at a specific COL,LEV
    auto src_lay = src_layout;
    auto tags = src_lay.tags();
    for (auto tag : tags) {
      if (tag != COL && tag != LEV && tag != ILEV) {
        src_lay = src_lay.strip_dim(tag);
      }
    }
    const auto  lname  = src.get_header().get_identifier().get_id_string()+"_mask";
    bool found = false;
    // Check if a field with lname has already been created:
    for (unsigned ii=0; ii<m_src_masks.size(); ii++) {
      const auto src_fld = m_src_masks[ii];
      if (lname == src_fld.name()) {
        auto& mask_tgt_fld = m_tgt_masks[ii];
        f_tgt.get_header().set_extra_data("mask_data",mask_tgt_fld);
        found = true;
        break;
      }
    }
    if (!found) {
      // We have to create this mask field and add it to the map so we can assign it to this tgt field as an extra data
      FieldIdentifier mask_src_fid (lname, src_lay, nondim, m_src_grid->name() );
      Field           mask_src_fld (mask_src_fid);
      mask_src_fld.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      mask_src_fld.allocate_view();
      const auto& tgt_lay = create_tgt_layout(src_lay);
      FieldIdentifier mask_tgt_fid (lname, tgt_lay, nondim, m_tgt_grid->name() );
      Field           mask_tgt_fld (mask_tgt_fid);
      mask_tgt_fld.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      mask_tgt_fld.allocate_view();
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_data"),
          "ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_data",mask_tgt_fld);
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
          "ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_value",m_mask_val);
      m_src_masks.push_back(mask_src_fld);
      m_tgt_masks.push_back(mask_tgt_fld);
    }
  } else {
    // If a field does not have LEV or ILEV it may still have mask tracking assigned from somewhere else.
    // In those cases we want to copy that mask tracking to the target field.
    // Note, we still make a new field to ensure it is defined on the target grid.
    if (src.get_header().has_extra_data("mask_data")) {
      auto f_src_mask = src.get_header().get_extra_data<Field>("mask_data");
      FieldIdentifier mask_tgt_fid (f_src_mask.name(), f_src_mask.get_header().get_identifier().get_layout(), nondim, m_tgt_grid->name() );
      Field           mask_tgt_fld (mask_tgt_fid);
      mask_tgt_fld.allocate_view();
      mask_tgt_fld.deep_copy(f_src_mask);

      auto& f_tgt    = m_tgt_fields[ifield];
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_data"),
          "ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_data",mask_tgt_fld);
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
          "ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_value",m_mask_val);
    }
  }
}

void VerticalRemapper::do_registration_ends ()
{
  // Check that the vertical profiles for the source data have been set
  EKAT_REQUIRE_MSG(m_mid_set,"Error::VerticalRemapper:registration_ends,\n"
    "Field for vertical profile of the source data for layout LEV has not been set.\n");
  EKAT_REQUIRE_MSG(m_int_set,"Error::VerticalRemapper:registration_ends,\n"
    "Field for vertical profile of the source data for layout ILEV has not been set.\n");
}

void VerticalRemapper::do_remap_fwd ()
{
  using namespace ShortFieldTagsNames;
  // Loop over each field
  const auto& tgt_pres_ap = m_remap_pres.get_header().get_alloc_properties();
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto  f_tgt    = m_tgt_fields[i];
    const auto& layout   = f_src.get_header().get_identifier().get_layout();
    const auto  src_tag  = layout.tags().back();
    const bool  do_remap = ekat::contains(std::vector<FieldTag>{ILEV,LEV},src_tag);
    if (do_remap) {
      // Dispatch kernel with the largest possible pack size
      const auto& src_ap = f_src.get_header().get_alloc_properties();
      const auto& tgt_ap = f_tgt.get_header().get_alloc_properties();
      const auto& src_pres_ap = src_tag == LEV ? m_src_mid.get_header().get_alloc_properties() : m_src_int.get_header().get_alloc_properties();
      if (src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          tgt_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          src_pres_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          tgt_pres_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>()) {
        apply_vertical_interpolation<SCREAM_PACK_SIZE>(f_src,f_tgt); 
      } else {
        apply_vertical_interpolation<1>(f_src,f_tgt); 
      }
    } else {
      // There is nothing to do, this field cannot be vertically interpolated,
      // so just copy it over.  Note, if this field has its own mask data make
      // sure that is copied too.
      if (f_tgt.get_header().has_extra_data("mask_data")) {
        auto f_tgt_mask = f_tgt.get_header().get_extra_data<Field>("mask_data");
        auto f_src_mask = f_src.get_header().get_extra_data<Field>("mask_data");
        f_tgt_mask.deep_copy(f_src_mask);
      }
      f_tgt.deep_copy(f_src);
    }
  }
  for (unsigned i=0; i<m_tgt_masks.size(); ++i) {
          auto& f_src    = m_src_masks[i];
          auto& f_tgt    = m_tgt_masks[i];
    const auto& layout   = f_src.get_header().get_identifier().get_layout();
    const auto  src_tag  = layout.tags().back();
    const bool  do_remap = ekat::contains(std::vector<FieldTag>{ILEV,LEV},src_tag);
    if (do_remap) {
      // If we are remapping then we need to initialize the mask source values to 1.0
      f_src.deep_copy(1.0);
      // Dispatch kernel with the largest possible pack size
      const auto& src_ap = f_src.get_header().get_alloc_properties();
      const auto& tgt_ap = f_tgt.get_header().get_alloc_properties();
      const auto& src_pres_ap = src_tag == LEV ? m_src_mid.get_header().get_alloc_properties() : m_src_int.get_header().get_alloc_properties();
      if (src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          tgt_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          src_pres_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
          tgt_pres_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>()) {
        apply_vertical_interpolation<SCREAM_PACK_SIZE>(f_src,f_tgt,true); 
      } else {
        apply_vertical_interpolation<1>(f_src,f_tgt,true); 
      }
    } else {
      // There is nothing to do, this field cannot be vertically interpolated,
      // so just copy it over.
      f_tgt.deep_copy(f_src);
    }
  }
}

template<int Packsize>
void VerticalRemapper::
apply_vertical_interpolation(const Field& f_src, const Field& f_tgt, const bool mask_interp) const
{
  using Pack = ekat::Pack<Real,Packsize>;
  using namespace ShortFieldTagsNames;
  using namespace scream::vinterp;
  const auto& layout = f_src.get_header().get_identifier().get_layout();
  const auto  rank   = f_src.rank();
  const auto src_tag = layout.tags().back();
  const auto src_num_levs = layout.dims().back();
  // ARG mask_interp checks if this is a vertical interpolation of the mask array that tracks masked 0.0 or not 1.0
  Real mask_val = mask_interp ? 0.0 : m_mask_val;

  Field    src_lev_f;
  if (src_tag == ILEV) {
    src_lev_f = m_src_int;
  } else {
    src_lev_f = m_src_mid;
  }
  auto src_lev  = src_lev_f.get_view<const Pack**>();
  auto remap_pres_view = m_remap_pres.get_view<Pack*>();
  switch(rank) {
    case 2:
    {
      auto src_view = f_src.get_view<const Pack**>();
      auto tgt_view = f_tgt.get_view<      Pack**>();
      perform_vertical_interpolation<Real,Packsize,2>(src_lev,remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs,mask_val);
      break;
    }
    case 3:
    {
      auto src_view = f_src.get_view<const Pack***>();
      auto tgt_view = f_tgt.get_view<      Pack***>();
      perform_vertical_interpolation<Real,Packsize,3>(src_lev,remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs,mask_val);
      break;
    }
    default:
      EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by VerticalRemapper.\n");
  }
}

} // namespace scream
