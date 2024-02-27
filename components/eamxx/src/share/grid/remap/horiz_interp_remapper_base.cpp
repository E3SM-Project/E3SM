#include "horiz_interp_remapper_base.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <numeric>

namespace scream
{

HorizInterpRemapperBase::
HorizInterpRemapperBase (const grid_ptr_type& fine_grid,
                         const std::string& map_file,
                         const InterpType type)
 : m_fine_grid(fine_grid)
 , m_map_file (map_file)
 , m_type (type)
 , m_comm (fine_grid->get_comm())
{
  // Sanity checks
  EKAT_REQUIRE_MSG (fine_grid->type()==GridType::Point,
      "Error! Horizontal interpolatory remap only works on PointGrid grids.\n"
      "  - fine grid name: " + fine_grid->name() + "\n"
      "  - fine_grid_type: " + e2str(fine_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (fine_grid->is_unique(),
      "Error! HorizInterpRemapperBase requires a unique fine grid.\n");

  // This is a special remapper. We only go in one direction
  m_bwd_allowed = false;

  // Get the remap data (if not already present, it will be built)
  auto& data = s_remapper_data[m_map_file];
  if (data.num_customers==0) {
    data.build(m_map_file,m_fine_grid,m_comm,m_type);
  }
  ++data.num_customers;

  m_row_offsets = data.row_offsets;
  m_col_lids = data.col_lids;
  m_weights = data.weights;

  // The grids really only matter for the horiz part. We may have 2+ remappers with
  // fine grids that only differ in terms of number of levs. Such remappers cannot
  // store the same coarse grid. So we soft-clone the grid, and reset the number of levels
  auto coarse_grid = data.coarse_grid->clone(data.coarse_grid->name(),true);
  auto ov_coarse_grid = data.ov_coarse_grid->clone(data.ov_coarse_grid->name(),true);

  // Reset num levs, and remove any geo data that depends on levs
  using namespace ShortFieldTagsNames;
  for (std::shared_ptr<AbstractGrid> grid : {coarse_grid,ov_coarse_grid}) {
    grid->reset_num_vertical_lev(fine_grid->get_num_vertical_levels());
    for (const auto& name : grid->get_geometry_data_names()) {
      const auto& f = grid->get_geometry_data(name);
      const auto& fl = f.get_header().get_identifier().get_layout();
      if (fl.has_tag(LEV) or fl.has_tag(ILEV)) {
        grid->delete_geometry_data(name);
      }
    }
  }
  m_coarse_grid = coarse_grid;
  m_ov_coarse_grid = ov_coarse_grid;

  if (m_type==InterpType::Refine) {
    set_grids(m_coarse_grid,m_fine_grid);
  } else {
    set_grids(m_fine_grid,m_coarse_grid);
  }
}

HorizInterpRemapperBase::
~HorizInterpRemapperBase ()
{
  auto it = s_remapper_data.find(m_map_file);
  if (it==s_remapper_data.end()) {
    // This would be very suspicious. But since the error is "benign",
    // and since we want to avoid throwing inside a destructor, just issue a warning.
    std::cerr << "WARNING! Remapper data for this map file was already deleted!\n"
                 " - map file: " << m_map_file << "\n";
    return;
  }

  --it->second.num_customers;
  if (it->second.num_customers==0) {
    s_remapper_data.erase(it);
  }
}

FieldLayout HorizInterpRemapperBase::
create_src_layout (const FieldLayout& tgt_layout) const
{
  EKAT_REQUIRE_MSG (m_src_grid!=nullptr,
      "Error! Cannot create source layout until the source grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[HorizInterpRemapperBase] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + to_string(tgt_layout));

  return create_layout (tgt_layout, m_src_grid);
}

FieldLayout HorizInterpRemapperBase::
create_tgt_layout (const FieldLayout& src_layout) const
{
  EKAT_REQUIRE_MSG (m_tgt_grid!=nullptr,
      "Error! Cannot create target layout until the target grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[HorizInterpRemapperBase] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + to_string(src_layout));

  return create_layout (src_layout, m_tgt_grid);
}

FieldLayout HorizInterpRemapperBase::
create_layout (const FieldLayout& fl_in,
               const grid_ptr_type& grid) const
{
  using namespace ShortFieldTagsNames;
  const auto type = get_layout_type(fl_in.tags());
        auto fl_out = FieldLayout::invalid();
  const bool midpoints = fl_in.has_tag(LEV);
  const bool is3d = fl_in.has_tag(LEV) or fl_in.has_tag(ILEV);
  switch (type) {
    case LayoutType::Scalar2D: [[ fallthrough ]];
    case LayoutType::Scalar3D:
      fl_out = is3d
             ? grid->get_3d_scalar_layout(midpoints)
             : grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D: [[ fallthrough ]];
    case LayoutType::Vector3D:
    {
      auto vtag = fl_in.get_vector_tag();
      auto vdim = fl_in.dim(vtag);
      fl_out = is3d
             ? grid->get_3d_vector_layout(midpoints,vtag,vdim)
             : grid->get_2d_vector_layout(vtag,vdim);
      break;
    }

    case LayoutType::Tensor2D: [[ fallthrough ]];
    case LayoutType::Tensor3D:
    {
      auto ttags = fl_in.get_tensor_tags();
      std::vector<int> tdims;
      for (auto idx : fl_in.get_tensor_dims()) {
        tdims.push_back(fl_in.dim(idx));
      }
      fl_out = is3d
             ? grid->get_3d_tensor_layout(midpoints,ttags,tdims)
             : grid->get_2d_tensor_layout(ttags,tdims);
      break;
    }

    default:
      EKAT_ERROR_MSG ("Layout not supported by HorizInterpRemapperBase:\n"
                      " - layout: " + to_string(fl_in) + "\n");
  }
  return fl_out;
}

void HorizInterpRemapperBase::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_fields ();
    setup_mpi_data_structures ();
  }
}

void HorizInterpRemapperBase::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  constexpr auto COL = ShortFieldTagsNames::COL;
  EKAT_REQUIRE_MSG (src.get_layout().has_tag(COL),
      "Error! Cannot register a field without COL tag in RefiningRemapperP2P.\n"
      "  - field name: " + src.name() + "\n"
      "  - field layout: " + to_string(src.get_layout()) + "\n");
  m_src_fields.push_back(field_type(src));
  m_tgt_fields.push_back(field_type(tgt));
}

void HorizInterpRemapperBase::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  EKAT_REQUIRE_MSG (src.data_type()==DataType::RealType,
      "Error! RefiningRemapperRMA only allows fields with RealType data.\n"
      "  - src field name: " + src.name() + "\n"
      "  - src field type: " + e2str(src.data_type()) + "\n");
  EKAT_REQUIRE_MSG (tgt.data_type()==DataType::RealType,
      "Error! RefiningRemapperRMA only allows fields with RealType data.\n"
      "  - tgt field name: " + tgt.name() + "\n"
      "  - tgt field type: " + e2str(tgt.data_type()) + "\n");

  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

  // If this was the last field to be bound, we can setup the MPI schedule
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    create_ov_fields ();
    setup_mpi_data_structures ();
  }
}

void HorizInterpRemapperBase::create_ov_fields ()
{
  m_ov_fields.reserve(m_num_fields);
  const auto num_ov_gids = m_ov_coarse_grid->get_num_local_dofs();
  const auto ov_gn = m_ov_coarse_grid->name();
  const auto dt = DataType::RealType;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_type==InterpType::Refine ? m_tgt_fields[i] : m_src_fields[i];
    const auto& fid = f.get_header().get_identifier();
    const auto layout = fid.get_layout().clone_with_different_extent(0,num_ov_gids);
    FieldIdentifier ov_fid (fid.name(),layout,fid.get_units(),ov_gn,dt);

    auto& ov_f = m_ov_fields.emplace_back(ov_fid);

    // Use same alloc props as fine fields, to allow packing in local_mat_vec
    const auto pack_size = f.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

template<int PackSize>
void HorizInterpRemapperBase::
local_mat_vec (const Field& x, const Field& y) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto row_grid = m_type==InterpType::Refine ? m_fine_grid : m_ov_coarse_grid;
  const int  nrows    = row_grid->get_num_local_dofs();

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int   rank       = src_layout.rank();

  auto row_offsets = m_row_offsets;
  auto col_lids    = m_col_lids;
  auto weights     = m_weights;

  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      // Unlike get_view, get_strided_view returns a LayoutStride view,
      // therefore allowing the 1d field to be a subfield of a 2d field
      // along the 2nd dimension.
      auto x_view = x.get_strided_view<const Real*>();
      auto y_view = y.get_strided_view<      Real*>();
      Kokkos::parallel_for(RangePolicy(0,nrows),
                           KOKKOS_LAMBDA(const int& row) {
        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        y_view(row) = weights(beg)*x_view(col_lids(beg));
        for (int icol=beg+1; icol<end; ++icol) {
          y_view(row) += weights(icol)*x_view(col_lids(icol));
        }
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = ESU::get_default_team_policy(nrows,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                            [&](const int j){
          y_view(row,j) = weights(beg)*x_view(col_lids(beg),j);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j) += weights(icol)*x_view(col_lids(icol),j);
          }
        });
      });
      break;
    }
    case 3:
    {
      auto x_view = x.get_view<const Pack***>();
      auto y_view = y.get_view<      Pack***>();
      const int dim1 = src_layout.dim(1);
      const int dim2 = PackInfo::num_packs(src_layout.dim(2));
      auto policy = ESU::get_default_team_policy(nrows,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k);
          }
        });
      });
      break;
    }
    case 4:
    {
      auto x_view = x.get_view<const Pack****>();
      auto y_view = y.get_view<      Pack****>();
      const int dim1 = src_layout.dim(1);
      const int dim2 = src_layout.dim(2);
      const int dim3 = PackInfo::num_packs(src_layout.dim(3));
      auto policy = ESU::get_default_team_policy(nrows,dim1*dim2*dim3);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                            [&](const int idx){
          const int j = (idx / dim3) / dim2;
          const int k = (idx / dim3) % dim2;
          const int l =  idx % dim3;
          y_view(row,j,k,l) = weights(beg)*x_view(col_lids(beg),j,k,l);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k,l) += weights(icol)*x_view(col_lids(icol),j,k,l);
          }
        });
      });
      break;
    }
    default:
    {
      EKAT_ERROR_MSG("Error::refining_remapper::local_mat_vec doesn't support fields of rank 4 or greater");
    }
  }
}

void HorizInterpRemapperBase::clean_up ()
{
  // Clear all fields
  m_src_fields.clear();
  m_tgt_fields.clear();
  m_ov_fields.clear();

  // Reset the state of the base class
  m_state = RepoState::Clean;
  m_num_fields = 0;
  m_num_registered_fields = 0;
  m_fields_are_bound.clear();
  m_num_bound_fields = 0;
}

std::map<std::string,HorizRemapperData> HorizInterpRemapperBase::s_remapper_data;

// ETI, so derived classes can call this method
template
void HorizInterpRemapperBase::
local_mat_vec<1>(const Field&, const Field&) const;

#if SCREAM_PACK_SIZE>1
template
void HorizInterpRemapperBase::
local_mat_vec<SCREAM_PACK_SIZE>(const Field&, const Field&) const;
#endif

} // namespace scream
