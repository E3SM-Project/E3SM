#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/util/scream_vertical_interpolation.hpp"
#include "share/io/scorpio_input.hpp"

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
{
  VerticalRemapper(src_grid,map_file,lev_prof,ilev_prof,-999999.0);
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

  // Gather the pressure level data for vertical remapping
  set_pressure_levels(map_file);

  // Create tgt_grid that is a clone of the src grid but with
  // the correct number of level
  auto tgt_grid = src_grid->clone("Point Grid",false);  //<-- ?? For Luca, why does the 1st arg have to be PointGrid?  Shouldn't the tgt_grid be of the same type as the src_grid if being cloned?  Note, choosing any other string threw an error for an "unsupported grid".
  tgt_grid->reset_num_vertical_lev(m_num_remap_levs);
  printf("ASD - tgt_grid num of levs = %d\n",tgt_grid->get_num_vertical_levels());
  this->set_grids(src_grid,tgt_grid);

  // Set the LEV and ILEV vertical profiles for interpolation from
  register_vertical_source_field(lev_prof,"mid");
  register_vertical_source_field(ilev_prof,"int");
}

VerticalRemapper::
~VerticalRemapper ()
{
  // Nothing to do
}

void VerticalRemapper::
do_print() {
  // Helper function to print the setup of the vertical remapper:
  const auto src_grid = get_src_grid();
  const auto tgt_grid = get_tgt_grid();
  if (m_comm.am_i_root()) {
    printf("-------------------------------------\n");
    printf("  Setup for vertical remapper\n");
    printf("  ---------------------------\n");
    // Print the remap layout
    printf("  Layout of src->tgt:  (%d, %d) -> (%d, %d)\n",src_grid->get_num_global_dofs(),src_grid->get_num_vertical_levels(),
                  tgt_grid->get_num_global_dofs(),tgt_grid->get_num_vertical_levels());
    // Print the set of fields to be mapped
    printf("  ---------------------------\n");
    printf("  Map pairs for %d fields:\n",m_num_fields);
    for (int i=0; i<m_num_fields; ++i) {
      const auto& f_src  = m_src_fields[i];
      const auto& f_tgt  = m_tgt_fields[i];
      printf("      %s -> %s\n",f_src.name().c_str(),f_tgt.name());
    }
    printf("-------------------------------------\n");
  }
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
  const bool midpoints = src_layout.has_tag(LEV);
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
  m_num_remap_levs = scorpio::get_dimlen_c2f(map_file.c_str(),"nlevs");
  std::vector<Real> remap_pres_levs(m_num_remap_levs);
  std::vector<scorpio::offset_t> dofs_offsets(m_num_remap_levs);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),0);
  const std::string idx_decomp_tag = "vertical_remapper::" + std::to_string(m_num_remap_levs);
  scorpio::get_variable(map_file, "p_levs", "p_levs", {"nlevs"}, "real", idx_decomp_tag);
  scorpio::set_dof(map_file,"p_levs",m_num_remap_levs,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"p_levs",-1,remap_pres_levs.data(),remap_pres_levs.size());
  scorpio::eam_pio_closefile(map_file);

  auto npacks = ekat::PackInfo<SCREAM_PACK_SIZE>::num_packs(m_num_remap_levs);
  m_remap_pres_view = view_1d<Pack>("",npacks);
  auto remap_pres_host = Kokkos::create_mirror_view(m_remap_pres_view);
  auto remap_pres_scal = ekat::scalarize(remap_pres_host);
  for (int ii=0;ii<m_num_remap_levs;ii++) {
    remap_pres_scal(ii) = remap_pres_levs[ii];
  }
  Kokkos::deep_copy(m_remap_pres_view,remap_pres_host);  
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
    src_mid = src;
    mid_set = true; 
   } else {  // mode=="int"
    auto layout = src_fid.get_layout();
    auto name   = src_fid.name();
    EKAT_REQUIRE_MSG(ekat::contains(std::vector<FieldTag>{ILEV},layout.tags().back()),
      "Error::VerticalRemapper::register_vertical_source_field,\n"
      "mode = 'int' expects a layour ending with ILEV tag.\n"
      " - field name  : " + name + "\n"
      " - field layout: " + to_string(layout) + "\n");
    EKAT_REQUIRE_MSG(src.is_allocated(), "Error! ILEV source field is not yet allocated.\n");
    src_int = src;
    int_set = true; 
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
  EKAT_REQUIRE_MSG (
      src.get_header().get_identifier().get_layout().rank()>1 ||
      src.get_header().get_alloc_properties().get_padding()==0,
      "Error! We don't support 2d scalar fields that are padded.\n");
  EKAT_REQUIRE_MSG (
      tgt.get_header().get_identifier().get_layout().rank()>1 ||
      tgt.get_header().get_alloc_properties().get_padding()==0,
      "Error! We don't support 2d scalar fields that are padded.\n");
  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

}

void VerticalRemapper::do_registration_ends ()
{
  // Check that the vertical profiles for the source data have been set
  EKAT_REQUIRE_MSG(mid_set,"Error::VerticalRemapper:registration_ends,\n"
    "Field for vertical profile of the source data for layout LEV has not been set.\n");
  EKAT_REQUIRE_MSG(int_set,"Error::VerticalRemapper:registration_ends,\n"
    "Field for vertical profile of the source data for layout ILEV has not been set.\n");

  do_print();
}

void VerticalRemapper::do_remap_fwd ()
{

  using namespace ShortFieldTagsNames;
  using namespace scream::vinterp;
  FieldTag src_tag;
  Field    src_lev_f;
  int      src_num_levs;
  // Loop over each field
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src  = m_src_fields[i];
          auto  f_tgt  = m_tgt_fields[i];
    const auto& layout = f_src.get_header().get_identifier().get_layout();
    const auto  rank   = f_src.rank();
    src_tag = layout.tags().back();
    src_num_levs = layout.dim(src_tag);
    const bool do_remap = ekat::contains(std::vector<FieldTag>{ILEV,LEV},src_tag);
    if (src_tag == ILEV) {
      src_lev_f = src_int;
    } else {
      src_lev_f = src_mid;
    }
    auto src_lev  = src_lev_f.get_view<const Pack**>();
    printf("ASD - do remap: %s, SCREAM_PACK_SIZE = %d\n",f_src.name().c_str(),SCREAM_PACK_SIZE);
    auto fap = f_tgt.get_header().get_alloc_properties();
    if (fap.is_compatible<ekat::Pack<Real,16>>()) {
      printf("    -    Compatible w/ packsize=16    \n");
    }
    else if (fap.is_compatible<ekat::Pack<Real,8>>()) {
      printf("    -    Compatible w/ packsize=8    \n");
    }
    else if (fap.is_compatible<ekat::Pack<Real,4>>()) {
      printf("    -    Compatible w/ packsize=4    \n");
    } else {
      printf("    -    Compatible w/ REAL    \n");
    } 
    if (do_remap) { 
      switch(rank) {
        case 1:
        {
          // This must be a single vertical profile
          break;
        }
        case 2:
        {
          auto src_view = f_src.get_view<const Pack**>();
          auto tgt_view = f_tgt.get_view<      Pack**>();
          perform_vertical_interpolation(src_lev,m_remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs,m_mask_val);
          break;
        }
//ASD TODO; perform vertical interpolation needs to be extended to views of size >2
//ASD        case 3:
//ASD        {
//ASD          auto src_view = f_src.get_view<const Pack***>();
//ASD          auto tgt_view = f_src.get_view<      Pack***>();
//ASD          perform_vertical_interpolation(src_lev,m_remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs);
//ASD          break;
//ASD        }
//ASD        case 4:
//ASD        {
//ASD          auto src_view = f_src.get_view<const Pack****>();
//ASD          auto tgt_view = f_src.get_view<      Pack****>();
//ASD          perform_vertical_interpolation(src_lev,m_remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs);
//ASD          break;
//ASD        }
//ASD        case 5:
//ASD        {
//ASD          auto src_view = f_src.get_view<const Pack*****>();
//ASD          auto tgt_view = f_src.get_view<      Pack*****>();
//ASD          perform_vertical_interpolation(src_lev,m_remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs);
//ASD          break;
//ASD        }
//ASD        case 6:
//ASD        {
//ASD          auto src_view = f_src.get_view<const Pack******>();
//ASD          auto tgt_view = f_src.get_view<      Pack******>();
//ASD          perform_vertical_interpolation(src_lev,m_remap_pres_view,src_view,tgt_view,src_num_levs,m_num_remap_levs);
//ASD          break;
//ASD        }
        default:
          EKAT_ERROR_MSG ("Error! Field rank (" + std::to_string(rank) + ") not supported by VerticalRemapper.\n");
      }
    } else {
      // There is nothing to do, this field cannot be vertically interpolated,
      // so just copy it over.
      f_tgt.deep_copy(f_src);
    }
  }

}

} // namespace scream
