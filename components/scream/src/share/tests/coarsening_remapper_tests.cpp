#include <catch2/catch.hpp>

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace scream {

class CoarseningRemapperTester : public CoarseningRemapper {
public:
  CoarseningRemapperTester (const grid_ptr_type& src_grid,
                            const std::string& map_file)
   : CoarseningRemapper(src_grid,map_file)
  {
    // Nothing to do
  }
  view_1d<gid_t>::HostMirror
  test_triplet_gids (const std::string& map_file) const {
    return CoarseningRemapper::get_my_triplets_gids (map_file);
  }

  view_2d<int> get_row_col_lids () const {
    return m_row_col_lids;
  }
  view_1d<Real> get_weights () const {
    return m_weights;
  }

  grid_ptr_type get_ov_tgt_grid () const {
    return m_ov_tgt_grid;
  }
};

template<typename ViewT>
typename ViewT::HostMirror
cmvc (const ViewT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

TEST_CASE ("coarsening_remap") {

  // -------------------------------------- //
  //           Init MPI and PIO             //
  // -------------------------------------- //

  ekat::Comm comm(MPI_COMM_WORLD);

  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  // -------------------------------------- //
  //           Set grid/map sizes           //
  // -------------------------------------- //

  const int nldofs_src = 10;
  const int nldofs_tgt =  5;
  const int ngdofs_src = nldofs_src*comm.size();
  const int ngdofs_tgt = nldofs_tgt*comm.size();
  const int nnz = ngdofs_src;

  // -------------------------------------- //
  //           Create a map file            //
  // -------------------------------------- //

  printf ("creating map file ...\n");

  std::string filename = "coarsening_map_file_np" + std::to_string(comm.size()) + ".nc";
  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::register_dimension(filename,"n_a", "n_a", ngdofs_src);
  scorpio::register_dimension(filename,"n_b", "n_b", ngdofs_tgt);
  scorpio::register_dimension(filename,"n_s", "n_s", nnz);

  scorpio::register_variable(filename,"col","col","none",{"n_s"},"real","int","Real-nnz");
  scorpio::register_variable(filename,"row","row","none",{"n_s"},"real","int","Real-nnz");
  scorpio::register_variable(filename,"S","S","none",{"n_s"},"real","real","Real-nnz");

  std::vector<std::int64_t> dofs (nldofs_src);
  std::iota(dofs.begin(),dofs.end(),comm.rank()*nldofs_src);
  scorpio::set_dof(filename,"col",nldofs_src,dofs.data());
  scorpio::set_dof(filename,"row",nldofs_src,dofs.data());
  scorpio::set_dof(filename,"S",nldofs_src,dofs.data());
  
  scorpio::eam_pio_enddef(filename);

  std::vector<Real> col,row,S;
  for (int i=0; i<nldofs_tgt; ++i) {
    // Tgt entry K is the avg of src entries K and K+ngdofs_tgt
    row.push_back(i+nldofs_tgt*comm.rank());
    col.push_back(i+nldofs_tgt*comm.rank());
    S.push_back(0.5);

    row.push_back(i+nldofs_tgt*comm.rank());
    col.push_back(i+nldofs_tgt*comm.rank() + ngdofs_tgt);
    S.push_back(0.5);
  }
  scorpio::grid_write_data_array(filename,"row",row.data(),row.size());
  scorpio::grid_write_data_array(filename,"col",col.data(),row.size());
  scorpio::grid_write_data_array(filename,"S",S.data(),row.size());

  scorpio::eam_pio_closefile(filename);
  printf ("creating map file ... done!\n");

  // -------------------------------------- //
  //      Build src grid and remapper       //
  // -------------------------------------- //
  
  scorpio::register_file(filename, scorpio::FileMode::Read);

  printf ("creating grid and remapper ...\n");

  AbstractGrid::dofs_list_type src_dofs("",nldofs_src);
  auto src_dofs_h = cmvc(src_dofs);
  std::iota(src_dofs_h.data(),src_dofs_h.data()+nldofs_src,nldofs_src*comm.rank());
  Kokkos::deep_copy(src_dofs,src_dofs_h);

  auto src_grid = std::make_shared<PointGrid>("src",nldofs_src,20,comm);
  src_grid->set_dofs(src_dofs);

  auto remap = std::make_shared<CoarseningRemapperTester>(src_grid,filename);
  printf ("creating grid and remapper ... done!\n");

  // Test tgt grid
  auto tgt_grid = remap->get_tgt_grid();
  REQUIRE (tgt_grid->get_num_local_dofs()==nldofs_tgt);
  REQUIRE (tgt_grid->get_num_global_dofs()==ngdofs_tgt);

  // Test sparse matrix
  auto my_triplets = remap->test_triplet_gids (filename);
  const int num_triplets = my_triplets.size();
  REQUIRE (num_triplets==(nldofs_tgt*2));
  auto rowcol_h = cmvc(remap->get_row_col_lids());
  auto weights_h = cmvc(remap->get_weights());
  for (int i=0; i<num_triplets; ++i) {
    REQUIRE (my_triplets[i]==(comm.rank()*num_triplets + i));
  }
  for (int i=0; i<nldofs_tgt; ++i) {
    const auto tgt_gid = comm.rank()*nldofs_tgt+i;
    REQUIRE (rowcol_h(2*i  ,0)==tgt_gid);
    REQUIRE (rowcol_h(2*i+1,0)==tgt_gid);

    REQUIRE (rowcol_h(2*i  ,1)==tgt_gid);
    REQUIRE (rowcol_h(2*i+1,1)==tgt_gid+ngdofs_tgt);

    REQUIRE (weights_h(2*i  )==0.5);
    REQUIRE (weights_h(2*i+1)==0.5);
  }

  // Test overlapped tgt grid
  auto ov_tgt_grid = remap->get_ov_tgt_grid ();
  const int num_loc_ov_tgt_gids = ov_tgt_grid->get_num_local_dofs();
  REQUIRE (num_loc_ov_tgt_gids==nldofs_tgt);
  auto ov_gids = ov_tgt_grid->get_dofs_gids_host();
  for (int i=0; i<nldofs_tgt; ++i) {
    const auto tgt_gid = comm.rank()*nldofs_tgt+i;
    REQUIRE (ov_gids[i]==tgt_gid);
  }

  // -------------------------------------- //
  //      Create src/tgt grid fields        //
  // -------------------------------------- //

  printf ("creating fields ...\n");
  constexpr int vec_dim = 3;
  auto create_field = [&] (const std::string& name,
                           const std::shared_ptr<const AbstractGrid>& grid,
                           const bool twod, const bool vec, const bool mid = false, const int ps = 1) -> Field
  {
    constexpr auto CMP = FieldTag::Component;
    constexpr auto units = ekat::units::Units::nondimensional();
    auto fl = twod
            ? (vec ? grid->get_2d_vector_layout (CMP,vec_dim)
                   : grid->get_2d_scalar_layout ())
            : (vec ? grid->get_3d_vector_layout (mid,CMP,vec_dim)
                   : grid->get_3d_scalar_layout (mid));
    FieldIdentifier fid(name,fl,units,grid->name());
    Field f(fid);
    f.get_header().get_alloc_properties().request_allocation(ps);
    f.allocate_view();
    return f;
  };

  auto src_s2d   = create_field("s2d",  src_grid,true,false);
  auto src_v2d   = create_field("v2d",  src_grid,true,true);
  auto src_s3d_m = create_field("s3d_m",src_grid,false,false,true, 1);
  auto src_s3d_i = create_field("s3d_i",src_grid,false,false,false,std::min(SCREAM_PACK_SIZE,4));
  auto src_v3d_m = create_field("v3d_m",src_grid,false,true ,true, std::min(SCREAM_PACK_SIZE,8));
  auto src_v3d_i = create_field("v3d_i",src_grid,false,true ,false,std::min(SCREAM_PACK_SIZE,16));

  auto tgt_s2d   = create_field("s2d",  tgt_grid,true,false);
  auto tgt_v2d   = create_field("v2d",  tgt_grid,true,true);
  auto tgt_s3d_m = create_field("s3d_m",tgt_grid,false,false,true, 1);
  auto tgt_s3d_i = create_field("s3d_i",tgt_grid,false,false,false,std::min(SCREAM_PACK_SIZE,4));
  auto tgt_v3d_m = create_field("v3d_m",tgt_grid,false,true ,true, std::min(SCREAM_PACK_SIZE,8));
  auto tgt_v3d_i = create_field("v3d_i",tgt_grid,false,true ,false,std::min(SCREAM_PACK_SIZE,16));

  std::vector<Field> src_f = {src_s2d,src_v2d,src_s3d_m,src_s3d_i,src_v3d_m,src_v3d_i};
  std::vector<Field> tgt_f = {tgt_s2d,tgt_v2d,tgt_s3d_m,tgt_s3d_i,tgt_v3d_m,tgt_v3d_i};

  printf ("creating fields ... done!\n");

  // -------------------------------------- //
  //     Register fields in the remapper    //
  // -------------------------------------- //

  printf ("registering fields ...\n");
  remap->registration_begins();
  remap->register_field(src_s2d,  tgt_s2d);
  remap->register_field(src_v2d,  tgt_v2d);
  remap->register_field(src_s3d_m,tgt_s3d_m);
  remap->register_field(src_s3d_i,tgt_s3d_i);
  remap->register_field(src_v3d_m,tgt_v3d_m);
  remap->register_field(src_v3d_i,tgt_v3d_i);
  remap->registration_ends();

  // No bwd remap
  REQUIRE_THROWS(remap->remap(false));
  printf ("registering fields ... done!\n");

  // -------------------------------------- //
  //       Generate data for src fields     //
  // -------------------------------------- //

  // Generate data in a deterministic way, so that when we check results,
  // we know a priori what the input data that generated the tgt field's
  // values was, even if that data was off rank.
  auto src_gids = src_grid->get_dofs_gids_host();
  for (const auto& f : src_f) {
    const auto& l = f.get_header().get_identifier().get_layout();
    switch (get_layout_type(l.tags())) {
      case LayoutType::Scalar2D:
      {
        const auto v_src = f.get_view<Real*,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          v_src(i) = src_gids(i);
        }
      } break;
      case LayoutType::Vector2D:
      {
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            v_src(i,j) = src_gids(i)*vec_dim + j;
        }}
      } break;
      case LayoutType::Scalar3D:
      {
        const int nlevs = l.dims().back();
        const auto v_src = f.get_view<Real**,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<nlevs; ++j) {
            v_src(i,j) = src_gids(i)*nlevs + j;
        }}
      } break;
      case LayoutType::Vector3D:
      {
        const int nlevs = l.dims().back();
        const auto v_src = f.get_view<Real***,Host>();
        for (int i=0; i<nldofs_src; ++i) {
          for (int j=0; j<vec_dim; ++j) {
            for (int k=0; k<nlevs; ++k) {
              v_src(i,j,k) = src_gids(i)*vec_dim*nlevs + j*nlevs + k;
        }}}
      } break;
      default:
        EKAT_ERROR_MSG ("Unexpected layout.\n");
    }
    f.sync_to_dev();
  }

  remap->remap(true);

  // -------------------------------------- //
  //          Check remapped fields         //
  // -------------------------------------- //

  // Recall, tgt gid K should be the avg of src gids K and K+ngdofs_tgt
  auto tgt_gids = tgt_grid->get_dofs_gids_host();
  for (size_t ifield=0; ifield<tgt_f.size(); ++ifield) {
    const auto& f = tgt_f[ifield];

    f.sync_to_host();

    const auto& l = f.get_header().get_identifier().get_layout();
    switch (get_layout_type(l.tags())) {
      case LayoutType::Scalar2D:
      {
        const auto v_tgt = f.get_view<const Real*,Host>();
        for (int i=0; i<nldofs_tgt; ++i) {
          const auto gid = tgt_gids(i);
          const auto term1 = gid;
          const auto term2 = gid+ngdofs_tgt;
          REQUIRE ( v_tgt(i)== (term1 + term2)/2.0 );
        }
      } break;
      case LayoutType::Vector2D:
      {
        const auto v_tgt = f.get_view<const Real**,Host>();
        for (int i=0; i<nldofs_tgt; ++i) {
          const auto gid = tgt_gids(i);
          for (int j=0; j<vec_dim; ++j) {
            const auto term1 = gid*vec_dim+j;
            const auto term2 = (gid+ngdofs_tgt)*vec_dim+j;
            REQUIRE ( v_tgt(i,j)== (term1 + term2)/2.0 );
        }}
      } break;
      case LayoutType::Scalar3D:
      {
        const int nlevs = l.dims().back();
        const auto v_tgt = f.get_view<const Real**,Host>();
        for (int i=0; i<nldofs_tgt; ++i) {
          const auto gid = tgt_gids(i);
          for (int j=0; j<nlevs; ++j) {
            const auto term1 = gid*nlevs+j;
            const auto term2 = (gid+ngdofs_tgt)*nlevs+j;
            REQUIRE ( v_tgt(i,j)== (term1 + term2)/2.0 );
        }}
      } break;
      case LayoutType::Vector3D:
      {
        const int nlevs = l.dims().back();
        const auto v_tgt = f.get_view<const Real***,Host>();
        for (int i=0; i<nldofs_tgt; ++i) {
          const auto gid = tgt_gids(i);
          for (int j=0; j<vec_dim; ++j) {
            for (int k=0; k<nlevs; ++k) {
              const auto term1 = gid*vec_dim*nlevs+j*nlevs+k;
              const auto term2 = (gid+ngdofs_tgt)*vec_dim*nlevs+j*nlevs+k;

              printf("tgt(%d,%d,%d) = %f\n",i,j,k,v_tgt(i,j,k));
              printf("  -> exp = %f\n",(term1+term2)/2.0);
              REQUIRE ( v_tgt(i,j,k)== (term1 + term2)/2.0 );
        }}}
      } break;
      default:
        EKAT_ERROR_MSG ("Unexpected layout.\n");
    }
  }

  // Clean up scorpio stuff
  scorpio::eam_pio_closefile(filename);
  scorpio::eam_pio_finalize();
}

} // namespace scream
