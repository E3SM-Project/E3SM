#include "coarsening_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

CoarseningRemapper::
CoarseningRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file,
                    const std::vector<Field>& mask_fields,
                    const std::map<std::string,int>& mask_map)
  : CoarseningRemapper(src_grid,map_file)
{
  m_track_mask  = true;
  m_mask_fields_src = mask_fields;
  m_mask_map_src    = mask_map;
}

CoarseningRemapper::
CoarseningRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file)
 : AbstractRemapper()
 , m_comm (src_grid->get_comm())
{
  m_track_mask = false;
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (src_grid->type()==GridType::Point,
      "Error! CoarseningRemapper only works on PointGrid grids.\n"
      "  - src grid name: " + src_grid->name() + "\n"
      "  - src_grid_type: " + e2str(src_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (src_grid->is_unique(),
      "Error! CoarseningRemapper requires a unique source grid.\n");

  // This is a coarsening remapper. We only go in one direction
  m_bwd_allowed = false;

  // Create io_grid, containing the indices of the triplets
  // in the map file that this rank has to read
  auto my_gids = get_my_triplets_gids (map_file,src_grid);
  auto io_grid = std::make_shared<PointGrid>("",my_gids.size(),0,m_comm);

  auto dofs_gids = io_grid->get_dofs_gids();
  auto dofs_gids_h = dofs_gids.get_view<gid_t*,Host>();
  std::memcpy(dofs_gids_h.data(),my_gids.data(),my_gids.size()*sizeof(gid_t));
  dofs_gids.sync_to_dev();

  // Create CRS matrix views

  // Read in triplets.
  const int nlweights = io_grid->get_num_local_dofs();
  std::vector<gid_t> row_gids_h(nlweights);
  std::vector<gid_t> col_gids_h(nlweights);
  std::vector<Real>  S_h (nlweights);

  // scream's gids are of type int, while scorpio wants long int as offsets.
  std::vector<scorpio::offset_t> dofs_offsets(nlweights);
  for (int i=0; i<nlweights; ++i) {
    dofs_offsets[i] = dofs_gids_h[i];
  }
  const std::string idx_decomp_tag = "coarsening_remapper::constructor_int_nnz" + std::to_string(nlweights);
  const std::string val_decomp_tag = "coarsening_remapper::constructor_real_nnz" + std::to_string(nlweights);

  scorpio::register_file(map_file,scorpio::FileMode::Read);
  scorpio::get_variable(map_file, "row", "row", {"n_s"}, "int", idx_decomp_tag);
  scorpio::get_variable(map_file, "col", "col", {"n_s"}, "int", idx_decomp_tag);
  scorpio::get_variable(map_file, "S",   "S",   {"n_s"}, "real", val_decomp_tag);
  scorpio::set_dof(map_file,"row",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"S",nlweights,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"row",-1,row_gids_h.data(),nlweights);
  scorpio::grid_read_data_array(map_file,"col",-1,col_gids_h.data(),nlweights);
  scorpio::grid_read_data_array(map_file,"S",  -1,S_h.data(),       nlweights);
  scorpio::eam_pio_closefile(map_file);

  // Offset the cols ids to match the source grid.
  // Determine the min id among the cols array, we
  // also add the min_dof for the grid.
  // Note: The cols field is not guaranteed to have the
  // min offset value, but the rows will be numbered
  // by 1 from least to greatest.  So we use row_gids to calculate
  // the remap min.
  int remap_min_dof = std::numeric_limits<int>::max();  // Really big INT
  for (int id=0; id<nlweights; id++) {
    remap_min_dof = std::min(row_gids_h[id],remap_min_dof);
  }
  int global_remap_min_dof;
  m_comm.all_reduce(&remap_min_dof,&global_remap_min_dof,1,MPI_MIN);

  gid_t col_offset = global_remap_min_dof - src_grid->get_global_min_dof_gid();
  for (int ii=0; ii<nlweights; ii++) {
    col_gids_h[ii] -= col_offset;
  }

  // Create an "overlapped" tgt grid, that is, a grid where each rank
  // owns all tgt rows that are affected by at least one of the cols
  // in its src_grid
  std::set<gid_t> ov_tgt_gids;
  for (int i=0; i<nlweights; ++i) {
    ov_tgt_gids.insert(row_gids_h[i]-1);
  }
  const int num_ov_tgt_gids = ov_tgt_gids.size();
  auto ov_tgt_grid = std::make_shared<PointGrid>("ov_tgt_grid",num_ov_tgt_gids,0,m_comm);
  auto ov_tgt_gids_h = ov_tgt_grid->get_dofs_gids().get_view<gid_t*,Host>();
  auto it = ov_tgt_gids.begin();
  for (int i=0; i<num_ov_tgt_gids; ++i, ++it) {
    ov_tgt_gids_h[i] = *it;
  }
  ov_tgt_grid->get_dofs_gids().sync_to_dev();

  m_ov_tgt_grid = ov_tgt_grid;
  const int num_ov_row_gids = m_ov_tgt_grid->get_num_local_dofs();

  // Now we have to create the weights CRS matrix
  m_row_offsets = view_1d<int>("",num_ov_row_gids+1);
  m_col_lids    = view_1d<int>("",nlweights);
  m_weights     = view_1d<Real>("",nlweights);

  // Sort col_gids_h and row_gids_h by row gid. It is easier to sort
  // the array [0,...,n), and use it later to index the row/col/weight
  // views in the correct order.
  std::vector<int> id (nlweights);
  std::iota(id.begin(),id.end(),0);
  auto compare = [&] (const int i, const int j) -> bool {
    return row_gids_h[i] < row_gids_h[j];
  };
  std::sort(id.begin(),id.end(),compare);

  // Create mirror views
  auto row_offsets_h = Kokkos::create_mirror_view(m_row_offsets);
  auto col_lids_h    = Kokkos::create_mirror_view(m_col_lids);
  auto weights_h     = Kokkos::create_mirror_view(m_weights);

  for (int i=0; i<nlweights; ++i) {
    col_lids_h(i) = gid2lid(col_gids_h[id[i]],src_grid);
    weights_h(i)  = S_h[id[i]];
  }

  Kokkos::deep_copy(m_weights,weights_h);
  Kokkos::deep_copy(m_col_lids,col_lids_h);

  // Compute row offsets
  std::vector<int> row_counts(num_ov_row_gids);
  for (int i=0; i<nlweights; ++i) {
    ++row_counts[gid2lid(row_gids_h[i]-1,m_ov_tgt_grid)];
  }
  std::partial_sum(row_counts.begin(),row_counts.end(),row_offsets_h.data()+1);
  EKAT_REQUIRE_MSG (
      row_offsets_h(num_ov_row_gids)==nlweights,
      "Error! Something went wrong while computing row offsets.\n"
      "  - local nnz       : " + std::to_string(nlweights) + "\n"
      "  - row_offsets(end): " + std::to_string(row_offsets_h(num_ov_row_gids)) + "\n");

  Kokkos::deep_copy(m_row_offsets,row_offsets_h);

  const int nlevs  = src_grid->get_num_vertical_levels();

  auto tgt_grid_gids = m_ov_tgt_grid->get_unique_gids ();
  const int ngids = tgt_grid_gids.size();

  auto tgt_grid = std::make_shared<PointGrid>("horiz_remap_tgt_grid",ngids,nlevs,m_comm);

  auto tgt_grid_gids_h = tgt_grid->get_dofs_gids().get_view<gid_t*,Host>();
  std::memcpy(tgt_grid_gids_h.data(),tgt_grid_gids.data(),ngids*sizeof(gid_t));
  tgt_grid->get_dofs_gids().sync_to_dev();

  this->set_grids(src_grid,tgt_grid);

  // Replicate the src grid geo data in the tgt grid. We use this remapper to do
  // the remapping (if needed), and clean it up afterwards.
  const auto& src_geo_data_names = src_grid->get_geometry_data_names();
  registration_begins();
  for (const auto& name : src_geo_data_names) {
    const auto& src_data = src_grid->get_geometry_data(name);
    const auto& src_data_fid = src_data.get_header().get_identifier();
    const auto& layout = src_data_fid.get_layout();
    if (layout.tags()[0]!=COL) {
      // Not a field to be coarsened (perhaps a vertical coordinate field).
      // Simply copy it in the tgt grid, but we still need to assign the new grid name.
      FieldIdentifier tgt_data_fid(src_data_fid.name(),src_data_fid.get_layout(),src_data_fid.get_units(),m_tgt_grid->name());
      auto tgt_data = tgt_grid->create_geometry_data(tgt_data_fid);
      tgt_data.deep_copy(src_data);
    } else {
      // This field needs to be remapped
      auto tgt_data_fid = create_tgt_fid(src_data_fid);
      auto tgt_data = tgt_grid->create_geometry_data(tgt_data_fid);
      register_field(src_data,tgt_data);
    }
  }
  registration_ends();
  if (get_num_fields()>0) {
    remap(true);

    // The remap phase only alters the fields on device. We need to sync them to host as well
    for (int i=0; i<get_num_fields(); ++i) {
      auto tgt_data = get_tgt_field(i);
      tgt_data.sync_to_host();
    }
  }
  clean_up();

}

CoarseningRemapper::
~CoarseningRemapper ()
{
  // We need to free MPI requests
  for (size_t i=0; i<m_send_req.size(); ++i) {
    MPI_Request_free(&m_send_req[i]);
  }
  for (size_t i=0; i<m_recv_req.size(); ++i) {
    MPI_Request_free(&m_recv_req[i]);
  }
}

FieldLayout CoarseningRemapper::
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
      EKAT_ERROR_MSG ("Layout not supported by CoarseningRemapper: " + e2str(lt) + "\n");
  }
  return src;
}
FieldLayout CoarseningRemapper::
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
      EKAT_ERROR_MSG ("Layout not supported by CoarseningRemapper: " + e2str(lt) + "\n");
  }
  return tgt;
}

void CoarseningRemapper::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_src_fields.push_back(field_type(src));
  m_tgt_fields.push_back(field_type(tgt));
}

void CoarseningRemapper::
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

  // If this was the last field to be bound, we can setup the MPI schedule
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    create_ov_tgt_fields ();
    setup_mpi_data_structures ();
  }
}

void CoarseningRemapper::do_registration_ends ()
{
  // Before finishing the registration we need to also register the set of masks, if needed.
  // NOTE: Every time `register_field` is called the `m_num_registered_fields` value is
  //       updated.  So we can use that to point to the appropriate mask field in the
  //       m_mask_map_tgt map by taking `m_num_registered_fields-1`.
  if (m_track_mask) {
    for (int idx=0; idx<m_mask_fields_src.size(); ++idx) {
      auto f   = m_mask_fields_src[idx];
      auto fid = f.get_header().get_identifier();
      register_field_from_src(fid);
      m_mask_map_tgt.emplace(idx,m_num_registered_fields-1);
      auto f_mask_extra = f.get_header().get_extra_data();
      // Make sure fields representing masks are not themselves meant to be masked.
      EKAT_REQUIRE(!f_mask_extra.count("mask_data"));
    }
    // Update the number of fields.
    // NOTE: We need to do this again because the `registration_ends()` call in the abstract
    //       remapper base class set the `m_num_fields` value before calling
    //       `do_registration_ends()`.
    m_num_fields = m_num_registered_fields;

    // Make sure that all fields registered so far are represented in the mask map, if missing
    // flag to not be masked at all.
    for (auto f : m_src_fields) {
      auto name = f.name();
      if (!m_mask_map_src.count(name)) {
        m_mask_map_src.emplace(name,-1);
      } else {
        // If there is a mask attached to this add a check that it is compatible with the field layout,
        // i.e. if layout = COL then mask should equal COL.  If layout is COL,LEV then mask is COL,LEV and anything
        // with higher dimensionality will use a mask of COL,LEV (for now).
        using namespace ShortFieldTagsNames;
        auto m_idx = m_mask_map_src.at(name);
        if (m_idx >-1) {
        auto m_fld = m_mask_fields_src[m_idx];
        auto m_lt = m_fld.get_header().get_identifier().get_layout();
        auto f_lt = f.get_header().get_identifier().get_layout();
        EKAT_REQUIRE(f_lt.has_tag(COL) == m_lt.has_tag(COL));
        EKAT_REQUIRE(f_lt.has_tag(LEV) == m_lt.has_tag(LEV));
        }
      }
    }

  }

  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_tgt_fields ();
    setup_mpi_data_structures ();
  }
}

void CoarseningRemapper::do_remap_fwd ()
{
  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    int ierr = MPI_Startall(m_recv_req.size(),m_recv_req.data());
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while starting persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // TODO: Add check that if there are mask values they are either 1's or 0's for unmasked/masked.

  // Loop over each field
  constexpr auto can_pack = SCREAM_PACK_SIZE>1;
  for (int i=0; i<m_num_fields; ++i) {
    // First, perform the local mat-vec. Recall that in these y=Ax products,
    // x is the src field, and y is the overlapped tgt field.
    const auto& f_src    = m_src_fields[i];
    const auto& f_ov_tgt = m_ov_tgt_fields[i];

    int mask_idx = -1;
    if (m_track_mask) {
      mask_idx = m_mask_map_src.at(f_src.name());
    }

    if (mask_idx != -1) {
      // Pass the mask to the local_mat_vec routine
      auto mask = m_mask_fields_src[mask_idx];
      // Dispatch kernel with the largest possible pack size
      const auto& src_ap = f_src.get_header().get_alloc_properties();
      const auto& ov_tgt_ap = f_ov_tgt.get_header().get_alloc_properties();
      if (can_pack && src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
                      ov_tgt_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>()) {
        local_mat_vec<SCREAM_PACK_SIZE>(f_src,f_ov_tgt,&mask);
      } else {
        local_mat_vec<1>(f_src,f_ov_tgt,&mask);
      }
    } else {
      // Dispatch kernel with the largest possible pack size
      const auto& src_ap = f_src.get_header().get_alloc_properties();
      const auto& ov_tgt_ap = f_ov_tgt.get_header().get_alloc_properties();
      if (can_pack && src_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>() &&
                      ov_tgt_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>()) {
        local_mat_vec<SCREAM_PACK_SIZE>(f_src,f_ov_tgt);
      } else {
        local_mat_vec<1>(f_src,f_ov_tgt);
      }
    }
  }

  // Pack, then fire off the sends
  pack_and_send ();

  // Wait for all data to be received, then unpack
  recv_and_unpack ();

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    int ierr = MPI_Waitall(m_send_req.size(),m_send_req.data(), MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // Rescale any fields that had the mask applied.
  if (m_track_mask) {
    for (int i=0; i<m_num_fields; ++i) {
      const auto& f_tgt = m_tgt_fields[i];
      const int   mask_idx = m_mask_map_src.at(f_tgt.name());
      if (mask_idx != -1) {
        // Then this field did use a mask
        const int mask_f_idx = m_mask_map_tgt.at(mask_idx);
        const auto& mask = m_tgt_fields[mask_f_idx];
        const auto  rank = f_tgt.get_header().get_identifier().get_layout().rank();
        const auto& tgt_ap = f_tgt.get_header().get_alloc_properties();
        if (can_pack && tgt_ap.is_compatible<RPack<SCREAM_PACK_SIZE>>()) {
          rescale_masked_fields<SCREAM_PACK_SIZE>(f_tgt,mask);
        } else {
          rescale_masked_fields<1>(f_tgt,mask);
        }
      }
    }
  }
}

template<int PackSize>
void CoarseningRemapper::
rescale_masked_fields (const Field& x, const Field& mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& layout = x.get_header().get_identifier().get_layout();
  const int rank = layout.rank();
  const int ncols = m_tgt_grid->get_num_local_dofs();
  const auto x_extra  = x.get_header().get_extra_data();
  Real mask_val = std::numeric_limits<float>::max()/10.0;
  if (x_extra.count("mask_value")) {
    mask_val = ekat::any_cast<Real>(x_extra.at("mask_value"));
  }
  const Real mask_threshold = std::numeric_limits<Real>::epsilon();  // TODO: Should we not hardcode the threshold for simply masking out the column.
  switch (rank) {
    case 1:
    {
      auto x_view =    x.get_view<      Real*>();
      auto m_view = mask.get_view<const Real*>();
      Kokkos::parallel_for(RangePolicy(0,ncols),
                           KOKKOS_LAMBDA(const int& icol) {
        if (m_view(icol)>mask_threshold) {
          x_view(icol) /= m_view(icol);
        } else {
          x_view(icol) = mask_val;
        }
      });
      break;
    }
    case 2:
    {
      auto x_view =    x.get_view<      Pack**>();
      auto m_view = mask.get_view<const Pack**>();
      const int dim1 = PackInfo::num_packs(layout.dim(1));
      auto policy = ESU::get_default_team_policy(ncols,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto icol = team.league_rank();
        auto x_sub = ekat::subview(x_view,icol);
        auto m_sub = ekat::subview(m_view,icol);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                            [&](const int j){
          auto masked = m_sub(j) > mask_threshold;
          if (masked.any()) {
            x_sub(j).set(masked,x_sub(j)/m_sub(j));
          }
          x_sub(j).set(!masked,mask_val);
        });
      });
      break;
    }
    case 3:
    {
      auto x_view =    x.get_view<      Pack***>();
      auto m_view = mask.get_view<const Pack**>();
      const int dim1 = layout.dim(1);
      const int dim2 = PackInfo::num_packs(layout.dim(2));
      auto policy = ESU::get_default_team_policy(ncols,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto icol = team.league_rank();
        auto m_sub      = ekat::subview(m_view,icol);

        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          auto x_sub = ekat::subview(x_view,icol,j);
          auto masked = m_sub(k) > mask_threshold;

          if (masked.any()) {
            x_sub(k).set(masked,x_sub(k)/m_sub(k));
          }
          x_sub(k).set(!masked,mask_val);
        });
      });
      break;
    }
  }


}

template<int PackSize>
void CoarseningRemapper::
local_mat_vec (const Field& x, const Field& y, const Field* mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int rank = src_layout.rank();
  const int nrows = m_ov_tgt_grid->get_num_local_dofs();
  auto row_offsets = m_row_offsets;
  auto col_lids = m_col_lids;
  auto weights = m_weights;
  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      auto x_view = x.get_view<const Real*>();
      auto y_view = y.get_view<      Real*>();
      view_1d<Real> mask_view;
      if (mask != nullptr) {
        mask_view = mask->get_view<Real*>();
      }
      Kokkos::parallel_for(RangePolicy(0,nrows),
                           KOKKOS_LAMBDA(const int& row) {
        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        if (mask != nullptr) {
          y_view(row) = weights(beg)*x_view(col_lids(beg))*mask_view(col_lids(beg));
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row) += weights(icol)*x_view(col_lids(icol))*mask_view(col_lids(icol));
          }
        } else {
          y_view(row) = weights(beg)*x_view(col_lids(beg));
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row) += weights(icol)*x_view(col_lids(icol));
          }
        }
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      view_2d<Pack> mask_view;
      if (mask != nullptr) {
        mask_view = mask->get_view<Pack**>();
      }
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = ESU::get_default_team_policy(nrows,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                            [&](const int j){
          if (mask != nullptr) {
            y_view(row,j) = weights(beg)*x_view(col_lids(beg),j)*mask_view(col_lids(beg),j);
            for (int icol=beg+1; icol<end; ++icol) {
              y_view(row,j) += weights(icol)*x_view(col_lids(icol),j)*mask_view(col_lids(icol),j);
            }
          } else {
            y_view(row,j) = weights(beg)*x_view(col_lids(beg),j);
            for (int icol=beg+1; icol<end; ++icol) {
              y_view(row,j) += weights(icol)*x_view(col_lids(icol),j);
            }
          }
        });
      });
      break;
    }
    case 3:
    {
      auto x_view = x.get_view<const Pack***>();
      auto y_view = y.get_view<      Pack***>();
      // Note, the mask is still assumed to be defined on COLxLEV so still only 2D for case 3.
      view_2d<Pack> mask_view;
      if (mask != nullptr) {
        mask_view = mask->get_view<Pack**>();
      }
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
          if (mask != nullptr) {
            y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k)*mask_view(col_lids(beg),k);
            for (int icol=beg+1; icol<end; ++icol) {
              y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k)*mask_view(col_lids(icol),k);
            }
          } else {
            y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k);
            for (int icol=beg+1; icol<end; ++icol) {
              y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k);
            }
          }
        });
      });
      break;
    }
    default:
    {
      EKAT_ERROR_MSG("Error::coarsening_remapper::local_mat_vec doesn't support fields of rank 4 or greater");
    }
  }
}

void CoarseningRemapper::pack_and_send ()
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const int num_send_gids = m_ov_tgt_grid->get_num_local_dofs();
  const auto pid_lid_start = m_send_pid_lids_start;
  const auto lids_pids = m_send_lids_pids;
  const auto buf = m_send_buffer;

  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f  = m_ov_tgt_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(fl.tags());
    const auto f_pid_offsets = ekat::subview(m_send_f_pid_offsets,ifield);

    switch (lt) {
      case LayoutType::Scalar2D:
      {
        auto v = f.get_view<const Real*>();
        Kokkos::parallel_for(RangePolicy(0,num_send_gids),
                             KOKKOS_LAMBDA(const int& i){
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          buf (offset + lidpos) = v(lid);
        });
      } break;
      case LayoutType::Vector2D:
      {
        auto v = f.get_view<const Real**>();
        const int ndims = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_send_gids,ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims),
                               [&](const int idim) {
            buf(offset + lidpos*ndims + idim) = v(lid,idim);
          });
        });
      } break;
      case LayoutType::Scalar3D:
      {
        auto v = f.get_view<const Real**>();
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_send_gids,nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs),
                               [&](const int ilev) {
            buf(offset + lidpos*nlevs + ilev) = v(lid,ilev);
          });
        });
      } break;
      case LayoutType::Vector3D:
      {
        auto v = f.get_view<const Real***>();
        const int ndims = fl.dim(1);
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_send_gids,ndims*nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims*nlevs),
                               [&](const int idx) {
            const int idim = idx / nlevs;
            const int ilev = idx % nlevs;
            buf(offset + lidpos*ndims*nlevs + idim*nlevs + ilev) = v(lid,idim,ilev);
          });
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_mpi_send_buffer,m_send_buffer);
  }

  if (not m_send_req.empty()) {
    int ierr = MPI_Startall(m_send_req.size(),m_send_req.data());
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while starting persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }
}

void CoarseningRemapper::recv_and_unpack ()
{
  if (not m_recv_req.empty()) {
    int ierr = MPI_Waitall(m_recv_req.size(),m_recv_req.data(), MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const int num_tgt_dofs = m_tgt_grid->get_num_local_dofs();

  const auto buf = m_recv_buffer;
  const auto recv_lids_beg = m_recv_lids_beg;
  const auto recv_lids_end = m_recv_lids_end;
  const auto recv_lids_pidpos = m_recv_lids_pidpos;
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
          auto& f  = m_tgt_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(fl.tags());
    const auto f_pid_offsets = ekat::subview(m_recv_f_pid_offsets,ifield);


    f.deep_copy(0);
    switch (lt) {
      case LayoutType::Scalar2D:
      {
        auto v = f.get_view<Real*>();
        Kokkos::parallel_for(RangePolicy(0,num_tgt_dofs),
                             KOKKOS_LAMBDA(const int& lid){
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos;
            v(lid) += buf (offset);
          }
        });
      } break;
      case LayoutType::Vector2D:
      {
        auto v = f.get_view<Real**>();
        const int ndims = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid)+lidpos*ndims;
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims),
                                 [&](const int idim) {
              v(lid,idim) += buf (offset + idim);
            });
          }
        });
      } break;
      case LayoutType::Scalar3D:
      {
        auto v = f.get_view<Real**>();
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos*nlevs;

            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs),
                                 [&](const int ilev) {
              v(lid,ilev) += buf (offset + ilev);
            });
          }
        });
      } break;
      case LayoutType::Vector3D:
      {
        auto v = f.get_view<Real***>();
        const int ndims = fl.dim(1);
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,nlevs*ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos*ndims*nlevs;

            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs*ndims),
                                 [&](const int idx) {
              const int idim = idx / nlevs;
              const int ilev = idx % nlevs;
              v(lid,idim,ilev) += buf (offset + idim*nlevs + ilev);
            });
          }
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}


std::vector<CoarseningRemapper::gid_t>
CoarseningRemapper::
get_my_triplets_gids (const std::string& map_file,
                      const grid_ptr_type& src_grid) const
{
  using namespace ShortFieldTagsNames;

  scorpio::register_file(map_file,scorpio::FileMode::Read);
  // 1. Create a "helper" grid, with as many dofs as the number
  //    of triplets in the map file, and divided linearly across ranks
  const int ngweights = scorpio::get_dimlen_c2f(map_file.c_str(),"n_s");
  const auto io_grid_linear = create_point_grid ("helper",ngweights,1,m_comm);
  const int nlweights = io_grid_linear->get_num_local_dofs();

  gid_t offset = nlweights;
  m_comm.scan(&offset,1,MPI_SUM);
  offset -= nlweights; // scan is inclusive, but we need exclusive

  // 2. Read a chunk of triplets col indices
  std::vector<gid_t> cols(nlweights);
  std::vector<gid_t> rows(nlweights); // Needed to calculate min_dof
  const std::string idx_decomp_tag = "coarsening_remapper::get_my_triplet_gids_int_dim" + std::to_string(nlweights);
  scorpio::get_variable(map_file, "col", "col", {"n_s"}, "int", idx_decomp_tag);
  scorpio::get_variable(map_file, "row", "row", {"n_s"}, "int", idx_decomp_tag);
  std::vector<scorpio::offset_t> dofs_offsets(nlweights);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),offset);
  scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"row",nlweights,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"col",-1,cols.data(),cols.size());
  scorpio::grid_read_data_array(map_file,"row",-1,rows.data(),rows.size());
  scorpio::eam_pio_closefile(map_file);

  // Offset the cols ids to match the source grid.
  // Determine the min id among the cols array, we
  // also add the min_dof for the grid.
  int remap_min_dof = std::numeric_limits<int>::max();
  for (int id=0; id<nlweights; id++) {
    remap_min_dof = std::min(rows[id],remap_min_dof);
  }
  int global_remap_min_dof;
  m_comm.all_reduce(&remap_min_dof,&global_remap_min_dof,1,MPI_MIN);

  gid_t col_offset = global_remap_min_dof - src_grid->get_global_min_dof_gid();
  for (auto& id : cols) {
    id -= col_offset;
  }

  // 3. Get the owners of the cols gids we read in, according to the src grid
  auto owners = src_grid->get_owners(cols);

  // 4. Group gids we read by the pid we need to send them to
  std::map<int,std::vector<int>> pid2gids_send;
  for (int i=0; i<nlweights; ++i) {
    const auto pid = owners[i];
    pid2gids_send[pid].push_back(i+offset);
  }

  // 5. Obtain the dual map of the one above: a list of gids we need to
  //    receive, grouped by the pid we recv them from
  auto pid2gids_recv = recv_gids_from_pids(pid2gids_send);

  // 6. Cat all the list of gids in one
  int num_my_triplets = 0;
  for (const auto& it : pid2gids_recv) {
    num_my_triplets += it.second.size();
  }
  std::vector<gid_t> my_triplets_gids(num_my_triplets);
  int num_copied = 0;
  for (const auto& it : pid2gids_recv) {
    auto dst = my_triplets_gids.data()+num_copied;
    auto src = it.second.data();
    std::memcpy (dst,src,it.second.size()*sizeof(int));
    num_copied += it.second.size();
  }

  return my_triplets_gids;
}

std::vector<int>
CoarseningRemapper::get_pids_for_recv (const std::vector<int>& send_to_pids) const
{
  // Figure out how many sends each PID is doing
  std::vector<int> num_sends(m_comm.size(),0);
  num_sends[m_comm.rank()] = send_to_pids.size();
  m_comm.all_gather(num_sends.data(),1);

  // Offsets for send_pids coming from each pid
  // NOTE: the extra entry at the end if for ease of use later
  std::vector<int> sends_offsets (m_comm.size()+1,0);
  for (int pid=1; pid<=m_comm.size(); ++pid) {
    sends_offsets[pid] = sends_offsets[pid-1] + num_sends[pid-1];
  }

  // Gather all the pids that each rank is sending data to
  auto nglobal_sends = std::accumulate(num_sends.begin(),num_sends.end(),0);
  std::vector<int> global_send_pids(nglobal_sends,-1);
  MPI_Allgatherv (send_to_pids.data(),send_to_pids.size(),MPI_INT,
                  global_send_pids.data(),num_sends.data(),sends_offsets.data(),
                  MPI_INT, m_comm.mpi_comm());

  // Loop over all the ranks, and all the pids they send to, and look for my pid
  std::vector<int> recv_from_pids;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int beg = sends_offsets[pid];
    const int end = sends_offsets[pid+1];

    for (int i=beg; i<end; ++i) {
      if (global_send_pids[i]==m_comm.rank()) {
        recv_from_pids.push_back(pid);
        break;
      }
    }
  }

  return recv_from_pids;
}

std::map<int,std::vector<int>>
CoarseningRemapper::
recv_gids_from_pids (const std::map<int,std::vector<int>>& pid2gids_send) const
{
  const auto comm = m_comm.mpi_comm();

  // First, figure out which PIDs I need to recv from
  std::vector<int> send_to;
  for (const auto& it : pid2gids_send) {
    send_to.push_back(it.first);
  }
  auto recv_from = get_pids_for_recv (send_to);

  // Send num of lids send to each send_pid, and recv num of lids
  // recv from each recv_pid
  std::map<int,int> nsends, nrecvs;
  std::vector<MPI_Request> send_req, recv_req;
  for (const auto& it : pid2gids_send) {
    const int pid = it.first;
    nsends[pid] = it.second.size();
    send_req.emplace_back();
    MPI_Isend (&nsends[pid],1,MPI_INT,pid,0,comm,&send_req.back());
  }
  for (auto pid : recv_from) {
    recv_req.emplace_back();
    MPI_Irecv(&nrecvs[pid],1,MPI_INT,pid,0,comm,&recv_req.back());
  }
  MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE);

  send_req.resize(0);
  recv_req.resize(0);

  // Now that we know how many gids we will get from each pid, we can send/recv them
  for (const auto& it : pid2gids_send) {
    send_req.emplace_back();
    MPI_Isend(it.second.data(),it.second.size(),MPI_INT,
              it.first,0,comm,&send_req.back());
  }
  std::map<int,std::vector<int>> pid2gids_recv;
  for (const auto& it : nrecvs) {
    const int pid = it.first;
    const int n   = it.second;
    pid2gids_recv[pid].resize(n);
    recv_req.emplace_back();
    MPI_Irecv(pid2gids_recv[pid].data(),n,MPI_INT,
              pid,0,comm,&recv_req.back());
  }
  MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE);

  return pid2gids_recv;
}

void CoarseningRemapper::create_ov_tgt_fields ()
{
  using FL = FieldLayout;
  m_ov_tgt_fields.reserve(m_num_fields);
  const int num_ov_cols = m_ov_tgt_grid->get_num_local_dofs();
  const auto ov_gn = m_ov_tgt_grid->name();
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src = m_src_fields[i];
    const auto& f_tgt = m_tgt_fields[i];
    const auto& fid = f_tgt.get_header().get_identifier();
    auto tags = fid.get_layout().tags();
    auto dims = fid.get_layout().dims();
    dims[0] = num_ov_cols;
    FieldIdentifier ov_fid (fid.name(),FL(tags,dims),fid.get_units(),ov_gn,fid.data_type());

    // Note: with C++17, emplace_back already returns a ref
    m_ov_tgt_fields.emplace_back(ov_fid);
    auto& ov_f = m_ov_tgt_fields.back();

    // Use same alloc props as src fields, to allow packing in local_mat_vec
    const auto pack_size = f_src.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

void CoarseningRemapper::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  const auto mpi_comm  = m_comm.mpi_comm();
  const auto mpi_real  = ekat::get_mpi_type<Real>();

  const int last_rank = m_comm.size()-1;

  // Pre-compute the amount of data stored in each field on each dof
  std::vector<int> field_col_size (m_num_fields);
  int sum_fields_col_sizes = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f  = m_ov_tgt_fields[i];  // Doesn't matter if tgt or ov_tgt
    const auto& fl = f.get_header().get_identifier().get_layout();
    if (fl.dim(0)>0) {
      field_col_size[i] = fl.size() / fl.dim(0);
    } else {
      field_col_size[i] = fl.size();
    }
    sum_fields_col_sizes += field_col_size[i];
  }

  // --------------------------------------------------------- //
  //                   Setup SEND structures                   //
  // --------------------------------------------------------- //

  // 1. Retrieve pid (and associated lid) of all ov_tgt gids
  //    on the tgt grid
  const auto ov_gids = m_ov_tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  auto gids_owners = m_tgt_grid->get_owners (ov_gids);

  // 2. Group dofs to send by remote pid
  const int num_ov_gids = ov_gids.size();
  std::map<int,std::vector<int>> pid2lids_send;
  std::map<int,std::vector<gid_t>> pid2gids_send;
  for (int i=0; i<num_ov_gids; ++i) {
    const int pid = gids_owners[i];
    pid2lids_send[pid].push_back(i);
    pid2gids_send[pid].push_back(ov_gids(i));
  }
  const int num_send_pids = pid2lids_send.size();
  m_send_lids_pids = view_2d<int>("",num_ov_gids,2);
  m_send_pid_lids_start = view_1d<int>("",m_comm.size());
  auto send_lids_pids_h = Kokkos::create_mirror_view(m_send_lids_pids);
  auto send_pid_lids_start_h = Kokkos::create_mirror_view(m_send_pid_lids_start);
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    send_pid_lids_start_h(pid) = pos;
    for (auto lid : pid2lids_send[pid]) {
      send_lids_pids_h(pos,0) = lid;
      send_lids_pids_h(pos++,1) = pid;
    }
  }
  Kokkos::deep_copy(m_send_lids_pids,send_lids_pids_h);
  Kokkos::deep_copy(m_send_pid_lids_start,send_pid_lids_start_h);

  // 3. Compute offsets in send buffer for each pid/field pair
  m_send_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  auto send_f_pid_offsets_h = Kokkos::create_mirror_view(m_send_f_pid_offsets);
  std::vector<int> send_pid_offsets(m_comm.size());
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    send_pid_offsets[pid] = pos;
    for (int i=0; i<m_num_fields; ++i) {
      send_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*pid2lids_send[pid].size();
    }

    // At the end, pos must match the total amount of data in the overlapped fields
    if (pid==last_rank) {
      EKAT_REQUIRE_MSG (pos==num_ov_gids*sum_fields_col_sizes,
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
    }
  }
  Kokkos::deep_copy (m_send_f_pid_offsets,send_f_pid_offsets_h);

  // 4. Allocate send buffers
  m_send_buffer = view_1d<Real>("",sum_fields_col_sizes*num_ov_gids);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // 5. Setup send requests
  m_send_req.reserve(num_send_pids);
  for (const auto& it : pid2lids_send) {
    const int n = it.second.size()*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const int pid = it.first;
    const auto send_ptr = m_mpi_send_buffer.data() + send_pid_offsets[pid];

    m_send_req.emplace_back();
    auto& req = m_send_req.back();
    MPI_Send_init (send_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }

  // --------------------------------------------------------- //
  //                   Setup RECV structures                   //
  // --------------------------------------------------------- //

  // 1. Obtain the dual map of send_gids: a list of gids we need to
  //    receive, grouped by the pid we recv them from
  const int num_tgt_dofs = m_tgt_grid->get_num_local_dofs();
  auto pid2gids_recv = recv_gids_from_pids(pid2gids_send);
  const int num_recv_pids = pid2gids_recv.size();

  // 2. Convert the gids to lids, and arrange them by lid
  std::vector<std::vector<int>> lid2pids_recv(num_tgt_dofs);
  int num_total_recv_gids = 0;
  for (const auto& it : pid2gids_recv) {
    const int pid = it.first;
    for (auto gid : it.second) {
      const int lid = gid2lid(gid,m_tgt_grid);
      lid2pids_recv[lid].push_back(pid);
    }
    num_total_recv_gids += it.second.size();
  }

  // 3. Splice the vector-of-vectors above in a 1d view,
  //    keeping track of where each lid starts/ends
  m_recv_lids_pidpos = view_2d<int>("",num_total_recv_gids,2);
  m_recv_lids_beg = view_1d<int>("",num_tgt_dofs);
  m_recv_lids_end = view_1d<int>("",num_tgt_dofs);
  auto recv_lids_pidpos_h = Kokkos::create_mirror_view(m_recv_lids_pidpos);
  auto recv_lids_beg_h  = Kokkos::create_mirror_view(m_recv_lids_beg);
  auto recv_lids_end_h  = Kokkos::create_mirror_view(m_recv_lids_end);

  auto tgt_dofs_h = m_tgt_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  for (int i=0,pos=0; i<num_tgt_dofs; ++i) {
    recv_lids_beg_h(i) = pos;
    const int gid = tgt_dofs_h[i];
    for (auto pid : lid2pids_recv[i]) {
      auto it = std::find(pid2gids_recv.at(pid).begin(),pid2gids_recv.at(pid).end(),gid);
      EKAT_REQUIRE_MSG (it!=pid2gids_recv.at(pid).end(),
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
      recv_lids_pidpos_h(pos,0) = pid;
      recv_lids_pidpos_h(pos++,1) = std::distance(pid2gids_recv.at(pid).begin(),it);
    }
    recv_lids_end_h(i) = pos;
  }
  Kokkos::deep_copy(m_recv_lids_pidpos,recv_lids_pidpos_h);
  Kokkos::deep_copy(m_recv_lids_beg,recv_lids_beg_h);
  Kokkos::deep_copy(m_recv_lids_end,recv_lids_end_h);

  // 3. Splice gids from all pids, and convert to lid, and compute the offset of each pid
  //    in the spliced list of gids.
  std::vector<int> recv_pid_start(m_comm.size()+1);
  for (int pid=0,pos=0; pid<=m_comm.size(); ++pid) {
    recv_pid_start[pid] = pos;
    pos += pid2gids_recv[pid].size();
  }

  // 4. Compute offsets in recv buffer for each pid/field pair
  m_recv_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  auto recv_f_pid_offsets_h = Kokkos::create_mirror_view(m_recv_f_pid_offsets);
  std::vector<int> recv_pid_offsets(m_comm.size());
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    recv_pid_offsets[pid] = pos;
    const int num_recv_gids = recv_pid_start[pid+1] - recv_pid_start[pid];
    for (int i=0; i<m_num_fields; ++i) {
      recv_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*num_recv_gids;
    }

    // // At the end, pos must match the total amount of data received
    if (pid==last_rank) {
      EKAT_REQUIRE_MSG (pos==num_total_recv_gids*sum_fields_col_sizes,
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
    }
  }
  Kokkos::deep_copy (m_recv_f_pid_offsets,recv_f_pid_offsets_h);

  // 5. Allocate recv buffers
  m_recv_buffer = view_1d<Real>("",sum_fields_col_sizes*num_total_recv_gids);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // 6. Setup recv requests
  m_recv_req.reserve(num_recv_pids);
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int num_recv_gids = recv_pid_start[pid+1] - recv_pid_start[pid];
    const int n = num_recv_gids*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const auto recv_ptr = m_mpi_recv_buffer.data() + recv_pid_offsets[pid];

    m_recv_req.emplace_back();
    auto& req = m_recv_req.back();
    MPI_Recv_init (recv_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }
}

void CoarseningRemapper::clean_up ()
{
  // Clear all MPI related structures
  m_send_buffer         = view_1d<Real>();
  m_recv_buffer         = view_1d<Real>();
  m_mpi_send_buffer     = mpi_view_1d<Real>();
  m_mpi_recv_buffer     = mpi_view_1d<Real>();
  m_send_f_pid_offsets  = view_2d<int>();
  m_recv_f_pid_offsets  = view_2d<int>();
  m_send_lids_pids      = view_2d<int>();
  m_send_pid_lids_start = view_1d<int>();
  m_recv_lids_pidpos    = view_2d<int>();
  m_recv_lids_beg       = view_1d<int>();
  m_recv_lids_end       = view_1d<int>();
  m_send_req.clear();
  m_recv_req.clear();

  // Clear all fields
  m_src_fields.clear();
  m_tgt_fields.clear();
  m_ov_tgt_fields.clear();

  // Reset the state of the base class
  m_state = RepoState::Clean;
  m_num_fields = 0;
  m_num_registered_fields = 0;
  m_fields_are_bound.clear();
  m_num_bound_fields = 0;
}

} // namespace scream
