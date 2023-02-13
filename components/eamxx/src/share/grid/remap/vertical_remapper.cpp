#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/scream_vertical_interpolation.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"

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
                  const Field& ilev_prof,
                  const Real mask_val)
  : VerticalRemapper(src_grid,map_file,lev_prof,ilev_prof)
{
  m_mask_val = mask_val;
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file,
                  const Field& lev_prof,
                  const Field& ilev_prof)
 : AbstractRemapper()
 , m_comm (src_grid->get_comm())
 , m_mask_val(std::numeric_limits<float>::max()/10.0)
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
  m_num_remap_levs = scorpio::get_dimlen_c2f(map_file.c_str(),"nlevs");
  scorpio::eam_pio_closefile(map_file);

  auto tgt_grid_gids = src_grid->get_unique_gids();
  const int ngids = tgt_grid_gids.size();
  auto tgt_grid = std::make_shared<PointGrid>("vertical_remap_tgt_grid",ngids,m_num_remap_levs,m_comm);
  auto tgt_grid_gids_h = tgt_grid->get_dofs_gids().get_view<gid_t*,Host>();
  std::memcpy(tgt_grid_gids_h.data(),tgt_grid_gids.data(),ngids*sizeof(gid_t));
  tgt_grid->get_dofs_gids().sync_to_dev();
  this->set_grids(src_grid,tgt_grid);

  // Replicate the src grid geo data in the tgt grid.
  const auto& src_geo_data_names = src_grid->get_geometry_data_names();
  for (const auto& name : src_geo_data_names) {
    const auto& src_data = src_grid->get_geometry_data(name);
    const auto& src_data_fid = src_data.get_header().get_identifier();
    const auto& layout = src_data_fid.get_layout();
    // We only add geo data that is horizontal in nature,  vertical geo data won't be added because the vertical structure
    // has a rigid definition already.
    if (layout.tags().back()!=LEV && layout.tags().back()!=ILEV) {
      // Simply copy it in the tgt grid, but we still need to assign the new grid name.
      FieldIdentifier tgt_data_fid(src_data_fid.name(),src_data_fid.get_layout(),src_data_fid.get_units(),m_tgt_grid->name());
      auto tgt_data = tgt_grid->create_geometry_data(tgt_data_fid);
      tgt_data.deep_copy(src_data);
    } 
  }


  // Set the LEV and ILEV vertical profiles for interpolation from
  register_vertical_source_field(lev_prof,"mid");
  register_vertical_source_field(ilev_prof,"int");

  // Gather the pressure level data for vertical remapping
  set_pressure_levels(map_file);
}

FieldLayout VerticalRemapper::
create_src_layout (const FieldLayout& tgt_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(tgt_layout.tags());
  auto src = FieldLayout::invalid();
  const bool midpoints = tgt_layout.has_tag(LEV);
  const int vec_dim = tgt_layout.is_vector_layout() ? tgt_layout.dim(CMP) : -1;
  switch (lt) {
    case LayoutType::Scalar2D:
      src = m_src_grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D:
      src = m_src_grid->get_2d_vector_layout(CMP,vec_dim);
      break;
    case LayoutType::Scalar3D:
      src = m_src_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
      src = m_src_grid->get_3d_vector_layout(midpoints,CMP,vec_dim);
      break;
    default:
      EKAT_ERROR_MSG ("Layout not supported by VerticalRemapper: " + e2str(lt) + "\n");
  }
  return src;
}
FieldLayout VerticalRemapper::
create_tgt_layout (const FieldLayout& src_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(src_layout.tags());
  auto tgt = FieldLayout::invalid();
  const bool midpoints = true; //src_layout.has_tag(LEV);
  const int vec_dim = src_layout.is_vector_layout() ? src_layout.dim(CMP) : -1;
  switch (lt) {
    case LayoutType::Scalar2D:
      tgt = m_tgt_grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D:
      tgt = m_tgt_grid->get_2d_vector_layout(CMP,vec_dim);
      break;
    case LayoutType::Scalar3D:
      tgt = m_tgt_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
      tgt = m_tgt_grid->get_3d_vector_layout(midpoints,CMP,vec_dim);
      break;
    default:
      EKAT_ERROR_MSG ("Layout not supported by VerticalRemapper: " + e2str(lt) + "\n");
  }
  return tgt;
}

void VerticalRemapper::
set_pressure_levels(const std::string& map_file) {
  scorpio::register_file(map_file,scorpio::FileMode::Read);

  using namespace ShortFieldTagsNames;
  std::vector<FieldTag> tags = {LEV};
  std::vector<int>      dims = {m_num_remap_levs};
  FieldLayout layout(tags,dims);
  FieldIdentifier fid("p_remap",layout,ekat::units::Pa,m_tgt_grid->name());
  m_remap_pres = Field(fid);
  m_remap_pres.get_header().get_alloc_properties().request_allocation(mPack::n);
  m_remap_pres.allocate_view();

  auto remap_pres_scal = m_remap_pres.get_view<Real*>();

  std::vector<scorpio::offset_t> dofs_offsets(m_num_remap_levs);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),0);
  const std::string idx_decomp_tag = "vertical_remapper::" + std::to_string(m_num_remap_levs);
  scorpio::get_variable(map_file, "p_levs", "p_levs", {"nlevs"}, "real", idx_decomp_tag);
  scorpio::set_dof(map_file,"p_levs",m_num_remap_levs,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"p_levs",-1,remap_pres_scal.data(),remap_pres_scal.size());
  scorpio::eam_pio_closefile(map_file);

}

void VerticalRemapper::
register_vertical_source_field(const Field& src, const std::string& mode)
{
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG(mode=="mid" || mode=="int","Error: VerticalRemapper::register_vertical_source_field,"
    "mode arg must be 'mid' or 'int'\n");

  auto src_fid = src.get_header().get_identifier();
  if (mode=="mid") {
    auto layout = src_fid.get_layout();
    auto name   = src_fid.name();
    EKAT_REQUIRE_MSG(ekat::contains(std::vector<FieldTag>{LEV},layout.tags().back()),
      "Error::VerticalRemapper::register_vertical_source_field,\n"
      "mode = 'mid' expects a layour ending with LEV tag.\n"
      " - field name  : " + name + "\n"
      " - field layout: " + to_string(layout) + "\n");
    EKAT_REQUIRE_MSG(src.is_allocated(), "Error! LEV source field is not yet allocated.\n");
    m_src_mid = src;
    m_mid_set = true; 
   } else {  // mode=="int"
    auto layout = src_fid.get_layout();
    auto name   = src_fid.name();
    EKAT_REQUIRE_MSG(ekat::contains(std::vector<FieldTag>{ILEV},layout.tags().back()),
      "Error::VerticalRemapper::register_vertical_source_field,\n"
      "mode = 'int' expects a layour ending with ILEV tag.\n"
      " - field name  : " + name + "\n"
      " - field layout: " + to_string(layout) + "\n");
    EKAT_REQUIRE_MSG(src.is_allocated(), "Error! ILEV source field is not yet allocated.\n");
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
    const auto  lt     = get_layout_type(src_lay.tags());
    const auto  lname  = src.get_header().get_identifier().get_id_string()+"_mask";
    bool found = false;
    // Check if a field with lname has already been created:
    for (int ii=0; ii<m_src_masks.size(); ii++) {
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
      auto tgt_extra = tgt.get_header().get_extra_data();
      EKAT_REQUIRE_MSG(!tgt_extra.count("mask_data"),"ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_data",mask_tgt_fld);
      EKAT_REQUIRE_MSG(!tgt_extra.count("mask_value"),"ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_value",m_mask_val);
      m_src_masks.push_back(mask_src_fld);
      m_tgt_masks.push_back(mask_tgt_fld);
    }
  } else {
    // If a field does not have LEV or ILEV it may still have mask tracking assigned from somewhere else.
    // In those cases we want to copy that mask tracking to the target field.
    // Note, we still make a new field to ensure it is defined on the target grid.
    const auto src_extra = src.get_header().get_extra_data();
    if (src_extra.count("mask_data")) {
      auto f_src_mask = ekat::any_cast<Field>(src_extra.at("mask_data"));
      FieldIdentifier mask_tgt_fid (f_src_mask.name(), f_src_mask.get_header().get_identifier().get_layout(), nondim, m_tgt_grid->name() );
      Field           mask_tgt_fld (mask_tgt_fid);
      mask_tgt_fld.allocate_view();
      mask_tgt_fld.deep_copy(f_src_mask);

      auto& f_tgt    = m_tgt_fields[ifield];
      auto tgt_extra = tgt.get_header().get_extra_data();
      EKAT_REQUIRE_MSG(!tgt_extra.count("mask_data"),"ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
      f_tgt.get_header().set_extra_data("mask_data",mask_tgt_fld);
      EKAT_REQUIRE_MSG(!tgt_extra.count("mask_value"),"ERROR VerticalRemapper::do_bind_field " + src.name() + " already has mask_data assigned!");
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
  constexpr auto can_pack = SCREAM_PACK_SIZE>1;
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
      if (can_pack && src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
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
      auto f_tgt_extra = f_tgt.get_header().get_extra_data();
      if (f_tgt_extra.count("mask_data")) {
        auto f_src_extra = f_src.get_header().get_extra_data();
        auto f_tgt_mask = ekat::any_cast<Field>(f_tgt_extra.at("mask_data"));
        auto f_src_mask = ekat::any_cast<Field>(f_src_extra.at("mask_data"));
        f_tgt_mask.deep_copy(f_src_mask);
      }
      f_tgt.deep_copy(f_src);
    }
  }
  for (int i=0; i<m_tgt_masks.size(); ++i) {
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
      if (can_pack && src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
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
