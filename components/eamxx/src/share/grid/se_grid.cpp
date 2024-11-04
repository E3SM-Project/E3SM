#include "share/grid/se_grid.hpp"
#include "share/field/field_utils.hpp"

#include <ekat/kokkos/ekat_subview_utils.hpp>

namespace scream {

SEGrid::
SEGrid (const std::string& grid_name,
        const int num_my_elements,
        const int num_gauss_pts,
        const int num_vertical_levels,
        const ekat::Comm& comm)
 : AbstractGrid (grid_name,GridType::SE,num_my_elements*num_gauss_pts*num_gauss_pts,num_vertical_levels,comm)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (num_my_elements>=0, "Error! Number of local elements must be non-negative.\n");
  EKAT_REQUIRE_MSG (num_gauss_pts>=2, "Error! Number of gauss points must be at least 2.\n");
  EKAT_REQUIRE_MSG (num_vertical_levels>=2, "Error! Number of vertical levels must be at least 2.\n");

  m_num_local_elem = num_my_elements;
  m_num_gp         = num_gauss_pts;
  get_comm().all_reduce(&m_num_local_elem,&m_num_global_elem,1,MPI_SUM);

  // Create dof_gids and lid2idx fields
  create_dof_fields (get_2d_scalar_layout().rank());

  // Create the cg dofs field
  using namespace ShortFieldTagsNames;
  const auto units = ekat::units::Units::nondimensional();
  m_cg_dofs_gids = Field(FieldIdentifier("cg_gids",FieldLayout({CMP},{get_num_local_dofs()}),units,this->name(),DataType::IntType));
  m_cg_dofs_gids.allocate_view();

  m_partitioned_dim_gids = Field(FieldIdentifier("el_gids",FieldLayout({EL},{m_num_local_elem}),units,this->name(),DataType::IntType));
  m_partitioned_dim_gids.allocate_view();
}

FieldLayout
SEGrid::get_2d_scalar_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({EL,GP,GP},{m_num_local_elem,m_num_gp,m_num_gp}).rename_dims(m_special_tag_names);
}

FieldLayout
SEGrid::get_2d_vector_layout (const int vector_dim, const std::string& vec_dim_name) const
{
  using namespace ShortFieldTagsNames;

  FieldLayout fl({EL,CMP,GP,GP},{m_num_local_elem,vector_dim,m_num_gp,m_num_gp});
  fl.rename_dim(1,vec_dim_name);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
SEGrid::get_2d_tensor_layout (const std::vector<int>& cmp_dims,
                              const std::vector<std::string>& cmp_names) const
{
  EKAT_REQUIRE_MSG (cmp_names.size()==cmp_dims.size(),
      "[SEGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names,",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims,",") + "\n");

  using namespace ShortFieldTagsNames;

  FieldLayout fl;

  fl = fl.append_dim(EL,m_num_local_elem);

  for (size_t i=0; i<cmp_dims.size(); ++i) {
    fl.append_dim(CMP,cmp_dims[i],cmp_names[i]);
  }
  fl.append_dim(GP,m_num_gp);

  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
SEGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({EL,GP,GP,VL},{m_num_local_elem,m_num_gp,m_num_gp,nvl}).rename_dims(m_special_tag_names);
}

FieldLayout
SEGrid::get_3d_vector_layout (const bool midpoints, const int vector_dim,
                              const std::string& vec_dim_name) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  FieldLayout fl({EL,CMP,GP,GP,VL},{m_num_local_elem,vector_dim,m_num_gp,m_num_gp,nvl});
  fl.rename_dim(1,vec_dim_name);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
SEGrid::get_3d_tensor_layout (const bool midpoints,
                              const std::vector<int>& cmp_dims,
                              const std::vector<std::string>& cmp_names) const
{
  EKAT_REQUIRE_MSG (cmp_names.size()==cmp_dims.size(),
      "[SEGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names,",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims,",") + "\n");

  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  FieldLayout fl;

  fl.append_dim(EL,m_num_local_elem);

  for (size_t i=0; i<cmp_dims.size(); ++i) {
    fl.append_dim(CMP,cmp_dims[i],cmp_names[i]);
  }
  fl.append_dim(GP,m_num_gp);
  fl.append_dim(VL,nvl);

  return fl.rename_dims(m_special_tag_names);
}

std::shared_ptr<AbstractGrid> SEGrid::clone (const std::string& clone_name, const bool shallow) const
{
  auto grid = std::make_shared<SEGrid>(clone_name,m_num_local_elem,m_num_gp,get_num_vertical_levels(),get_comm());

  grid->copy_data(*this,shallow);

  if (m_cg_dofs_gids.is_allocated()) {
    if (shallow) {
      grid->m_cg_dofs_gids = m_cg_dofs_gids;
    } else {
      grid->m_cg_dofs_gids = m_cg_dofs_gids.clone();
    }
  }

  return grid;
}

bool SEGrid::check_valid_dofs () const
{
  return is_unique();
}

bool SEGrid::check_valid_lid_to_idx () const
{
  auto h_l2i = get_lid_to_idx_map().get_view<const int**,Host>();
  for (int idof=0; idof<h_l2i.extent_int(0); ++idof) {
    auto elgpgp = ekat::subview(h_l2i,idof);
    const int el = elgpgp(0);
    const int ip = elgpgp(1);
    const int jp = elgpgp(2);

    if (el<0 || el>=m_num_local_elem) return false;
    if (ip<0 || ip>=m_num_gp) return false;
    if (jp<0 || jp>=m_num_gp) return false;
  }
  return true;
}

} // namespace scream
