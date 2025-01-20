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

void HorizInterpRemapperBase::registration_ends_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    const auto& src_dt = m_src_fields[i].get_header().get_identifier().data_type();
    const auto& tgt_dt = m_tgt_fields[i].get_header().get_identifier().data_type();
    EKAT_REQUIRE_MSG (src_dt==DataType::RealType and tgt_dt==DataType::RealType,
        "Error! HorizInterpRmapperBase requires src/tgt fields to have Real data type.\n"
        "  - src field name: " + m_src_fields[i].name() + "\n"
        "  - tgt field name: " + m_tgt_fields[i].name() + "\n"
        "  - src data type : " + e2str(src_dt) + "\n"
        "  - tgt data type : " + e2str(tgt_dt) + "\n");
  }

  create_ov_fields ();
  setup_mpi_data_structures ();
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
    const auto layout = fid.get_layout().clone().reset_dim(0,num_ov_gids);
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
      EKAT_ERROR_MSG("[HorizInterpRemapperBase::local_mat_vec] Error! Fields of rank 4 or greater are not supported.\n");
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
