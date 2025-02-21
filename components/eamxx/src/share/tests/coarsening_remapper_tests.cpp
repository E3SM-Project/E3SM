#include <catch2/catch.hpp>

#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

class CoarseningRemapperTester : public CoarseningRemapper {
public:
  using gid_type = AbstractGrid::gid_type;

  CoarseningRemapperTester (const grid_ptr_type& src_grid,
                            const std::string& map_file)
   : CoarseningRemapper(src_grid,map_file)
  {
    // Nothing to do
  }

  // Note: we use this instead of get_tgt_grid, b/c the nonconst grid
  //       will give use a not read-only gids field, so we can pass
  //       pointers to MPI_Bcast (which needs pointer to nonconst)
  std::shared_ptr<AbstractGrid> get_coarse_grid () const {
    return m_coarse_grid;
  }

  view_1d<int> get_row_offsets () const {
    return m_row_offsets;
  }
  view_1d<int> get_col_lids () const {
    return m_col_lids;
  }
  view_1d<Real> get_weights () const {
    return m_weights;
  }

  grid_ptr_type get_ov_tgt_grid () const {
    return m_ov_coarse_grid;
  }

  view_2d<int>::HostMirror get_send_f_pid_offsets () const {
    return cmvdc(m_send_f_pid_offsets);
  }
  view_2d<int>::HostMirror get_recv_f_pid_offsets () const {
    return cmvdc(m_recv_f_pid_offsets);
  }

  view_1d<int>::HostMirror get_recv_lids_beg () const {
    return cmvdc(m_recv_lids_beg);
  }
  view_1d<int>::HostMirror get_recv_lids_end () const {
    return cmvdc(m_recv_lids_end);
  }

  view_2d<int>::HostMirror get_send_lids_pids () const {
    return cmvdc(m_send_lids_pids );
  }
  view_2d<int>::HostMirror get_recv_lids_pidpos () const {
    return cmvdc(m_recv_lids_pidpos);
  }

  view_1d<int>::HostMirror get_send_pid_lids_start () const {
    return cmvdc(m_send_pid_lids_start);
  }
};

void root_print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// Create a source grid given number of global dofs.
// Dofs are scattered around randomly
template<typename Engine>
std::shared_ptr<AbstractGrid>
build_src_grid(const ekat::Comm& comm, const int ngdofs, Engine& engine) 
{
  using gid_type = AbstractGrid::gid_type;
  const int nlevs = 20;

  std::vector<gid_type> all_dofs (ngdofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs.data(),all_dofs.data()+all_dofs.size(),0);
    std::shuffle(all_dofs.data(),all_dofs.data()+ngdofs,engine);
  }
  comm.broadcast(all_dofs.data(),ngdofs,comm.root_rank());

  int nldofs = ngdofs / comm.size();
  int remainder = ngdofs % comm.size();
  int offset = nldofs * comm.rank() + std::min(comm.rank(),remainder);
  if (comm.rank()<remainder) {
    ++nldofs;
  }
  auto src_grid = std::make_shared<PointGrid>("src",nldofs,nlevs,comm);

  auto src_dofs = src_grid->get_dofs_gids();
  auto src_dofs_h = src_dofs.get_view<gid_type*,Host>();
  std::copy_n(all_dofs.data()+offset,nldofs,src_dofs_h.data());
  src_dofs.sync_to_dev();

  return src_grid;
}

constexpr int vec_dim = 2;
constexpr int tens_dim1 = 3;
constexpr int tens_dim2 = 4;
Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid, const bool midpoints)
{
  const auto u = ekat::units::Units::nondimensional();
  const auto& gn = grid.name();
  Field f;
  switch (lt) {
    case LayoutType::Scalar2D:
      f = Field(FieldIdentifier(name,grid.get_2d_scalar_layout(),u,gn));  break;
    case LayoutType::Vector2D:
      f = Field(FieldIdentifier(name,grid.get_2d_vector_layout(vec_dim),u,gn));  break;
    case LayoutType::Tensor2D:
      f = Field(FieldIdentifier(name,grid.get_2d_tensor_layout({tens_dim1,tens_dim2}),u,gn));  break;
    case LayoutType::Scalar3D:
      f = Field(FieldIdentifier(name,grid.get_3d_scalar_layout(midpoints),u,gn));
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      break;
    case LayoutType::Vector3D:
      f = Field(FieldIdentifier(name,grid.get_3d_vector_layout(midpoints,vec_dim),u,gn));
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      break;
    case LayoutType::Tensor3D:
      f = Field(FieldIdentifier(name,grid.get_3d_tensor_layout(midpoints,{tens_dim1,tens_dim2}),u,gn));
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      break;
    default:
      EKAT_ERROR_MSG ("Invalid layout type for this unit test.\n");
  }
  f.allocate_view();

  return f;
}

template<typename Engine>
Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid, const bool midpoints, Engine& engine) {
  auto f = create_field(name,lt,grid,midpoints);

  // Use discrete_distribution to get an integer, then use that as exponent for 2^-n.
  // This guarantees numbers that are exactly represented as FP numbers, which ensures
  // the test will produce the expected answer, regardless of how math ops are performed.
  using IPDF = std::discrete_distribution<int>;
  IPDF ipdf ({1,1,1,1,1,1,1,1,1,1});
  auto pdf = [&](Engine& e) {
    return Real(std::pow(2,ipdf(e)));
  };
  randomize(f,engine,pdf);

  return f;
}

template<typename T>
Field all_gather_field_impl (const Field& f, const ekat::Comm& comm) {
  constexpr auto COL = ShortFieldTagsNames::COL;
  const auto& fid = f.get_header().get_identifier();
  const auto& fl  = fid.get_layout();
  int col_size = fl.clone().strip_dim(COL).size();
  auto tags = fl.tags();
  auto dims = fl.dims();
  int my_cols = dims[0];;
  comm.all_reduce(&my_cols, &dims.front(), 1, MPI_SUM );
  FieldLayout gfl(tags,dims);
  FieldIdentifier gfid("g" + f.name(),gfl,fid.get_units(),fid.get_grid_name(),fid.data_type());
  Field gf(gfid);
  gf.allocate_view();
  std::vector<T> data_vec(col_size);
  f.sync_to_host();
  for (int pid=0,offset=0; pid<comm.size(); ++pid) {
    T* data;
    int ncols = fl.dims()[0];
    comm.broadcast(&ncols,1,pid);
    for (int icol=0; icol<ncols; ++icol,offset+=col_size) {
      switch (fl.rank()) {
        case 1:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<T*,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        case 2:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<T**,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        case 3:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<T***,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        case 4:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<T****,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        default:
          EKAT_ERROR_MSG (
              "Unexpected rank in RefiningRemapperRMA unit test.\n"
              "  - field name: " + f.name() + "\n");
      }
      comm.broadcast(data,col_size,pid);
      auto gdata = gf.get_internal_view_data<T,Host>()+offset;
      std::copy(data,data+col_size,gdata);
    }
  }
  return gf;
}

Field all_gather_field (const Field& f, const ekat::Comm& comm) {
  const auto dt = f.data_type();
  if (dt==DataType::RealType) {
    return all_gather_field_impl<Real>(f,comm);
  } else {
    return all_gather_field_impl<int>(f,comm);
  }
}

// Helper function to create a remap file
void create_remap_file(const std::string& filename, const int ngdofs_tgt)
{
  const int ngdofs_src = ngdofs_tgt + 1;
  const int nnz = 2*ngdofs_tgt;

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::define_dim(filename,"n_a", ngdofs_src);
  scorpio::define_dim(filename,"n_b", ngdofs_tgt);
  scorpio::define_dim(filename,"n_s", nnz);

  scorpio::define_var(filename,"col",{"n_s"},"int");
  scorpio::define_var(filename,"row",{"n_s"},"int");
  scorpio::define_var(filename,"S"  ,{"n_s"},"double");

  scorpio::enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz,0.5);
  for (int i=0; i<ngdofs_tgt; ++i) {
    row[2*i] = i;
    row[2*i+1] = i;
    col[2*i] = i;
    col[2*i+1] = i+1;
  }

  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",    S.data());

  scorpio::release_file(filename);
}

TEST_CASE("coarsening_remap")
{
  auto& catch_capture = Catch::getResultCapture();

  // This is a simple test to just make sure the coarsening remapper works
  // when the map itself has more remap triplets than the size of the 
  // source and target grid.  This is typical in monotone remappers from
  // fine to coarse meshes.

  // -------------------------------------- //
  //           Init MPI and PIO             //
  // -------------------------------------- //

  ekat::Comm comm(MPI_COMM_WORLD);

  root_print ("\n +---------------------------------+\n",comm);
  root_print (" |   Testing coarsening remapper   |\n",comm);
  root_print (" +---------------------------------+\n\n",comm);

  scorpio::init_subsystem(comm);
  auto engine = setup_random_test (&comm);

  // -------------------------------------- //
  //           Create a map file            //
  // -------------------------------------- //

  std::string filename = "cr_tests_map." + std::to_string(comm.size()) + ".nc";

  const int nldofs_tgt = 2;
  const int ngdofs_tgt = nldofs_tgt*comm.size();
  create_remap_file(filename, ngdofs_tgt);

  // -------------------------------------- //
  //      Build src grid and remapper       //
  // -------------------------------------- //

  const int ngdofs_src = ngdofs_tgt+1;
  auto src_grid = build_src_grid(comm, ngdofs_src, engine);
  auto remap = std::make_shared<CoarseningRemapperTester>(src_grid,filename);

  // -------------------------------------- //
  //      Create src/tgt grid fields        //
  // -------------------------------------- //

  // The other test checks remapping for fields of multiple dimensions.
  // Here we will simplify and just remap a simple 2D horizontal field.
  auto tgt_grid = remap->get_coarse_grid();

  auto src_s2d   = create_field("s2d",  LayoutType::Scalar2D, *src_grid, false, engine);
  auto src_v2d   = create_field("v2d",  LayoutType::Vector2D, *src_grid, false, engine);
  auto src_t2d   = create_field("t2d",  LayoutType::Tensor2D, *src_grid, false, engine);
  auto src_s3d_m = create_field("s3d_m",LayoutType::Scalar3D, *src_grid, true,  engine);
  auto src_s3d_i = create_field("s3d_i",LayoutType::Scalar3D, *src_grid, false, engine);
  auto src_v3d_m = create_field("v3d_m",LayoutType::Vector3D, *src_grid, true,  engine);
  auto src_v3d_i = create_field("v3d_i",LayoutType::Vector3D, *src_grid, false, engine);
  auto src_t3d_m = create_field("t3d_m",LayoutType::Tensor3D, *src_grid, true,  engine);
  auto src_t3d_i = create_field("t3d_i",LayoutType::Tensor3D, *src_grid, false, engine);

  auto tgt_s2d   = create_field("s2d",  LayoutType::Scalar2D, *tgt_grid, false);
  auto tgt_v2d   = create_field("v2d",  LayoutType::Vector2D, *tgt_grid, false);
  auto tgt_t2d   = create_field("t2d",  LayoutType::Tensor2D, *tgt_grid, false);
  auto tgt_s3d_m = create_field("s3d_m",LayoutType::Scalar3D, *tgt_grid, true );
  auto tgt_s3d_i = create_field("s3d_i",LayoutType::Scalar3D, *tgt_grid, false);
  auto tgt_v3d_m = create_field("v3d_m",LayoutType::Vector3D, *tgt_grid, true );
  auto tgt_v3d_i = create_field("v3d_i",LayoutType::Vector3D, *tgt_grid, false);
  auto tgt_t3d_m = create_field("t3d_m",LayoutType::Tensor3D, *tgt_grid, true );
  auto tgt_t3d_i = create_field("t3d_i",LayoutType::Tensor3D, *tgt_grid, false);

  std::vector<Field> src_f = {src_s2d,src_v2d,src_t2d,src_s3d_m,src_s3d_i,src_v3d_m,src_v3d_i,src_t3d_m,src_t3d_i};
  std::vector<Field> tgt_f = {tgt_s2d,tgt_v2d,tgt_t2d,tgt_s3d_m,tgt_s3d_i,tgt_v3d_m,tgt_v3d_i,tgt_t3d_m,tgt_t3d_i};

  // -------------------------------------- //
  //     Register fields in the remapper    //
  // -------------------------------------- //

  remap->registration_begins();
  for (size_t i=0; i<tgt_f.size(); ++i) {
    remap->register_field(src_f[i],tgt_f[i]);
  }
  remap->registration_ends();

  // -------------------------------------- //
  //          Check remapped fields         //
  // -------------------------------------- //

  Real w = 0.5;
  auto gids_tgt = all_gather_field(tgt_grid->get_dofs_gids(),comm);
  auto gids_src = all_gather_field(src_grid->get_dofs_gids(),comm);
  auto gids_src_v = gids_src.get_view<const AbstractGrid::gid_type*,Host>();
  auto gids_tgt_v = gids_tgt.get_view<const AbstractGrid::gid_type*,Host>();

  auto gid2lid = [&](const int gid, const auto gids_v) {
    auto data = gids_v.data();
    auto it = std::find(data,data+gids_v.size(),gid);
    return std::distance(data,it);
  };
  for (int irun=0; irun<5; ++irun) {
    root_print (" -> Run " + std::to_string(irun) + "\n",comm);
    remap->remap_fwd();

    // Recall, tgt gid K should be the avg of local src_gids
    for (size_t ifield=0; ifield<tgt_f.size(); ++ifield) {
      auto gsrc = all_gather_field(src_f[ifield],comm);
      auto gtgt = all_gather_field(tgt_f[ifield],comm);

      const auto& l = gsrc.get_header().get_identifier().get_layout();
      const auto ls = l.to_string();
      std::string dots (30-ls.size(),'.');
      auto msg = "   -> Checking field with layout " + ls + " " + dots;
      root_print (msg + "\n",comm);
      bool ok = true;
      switch (l.type()) {
        case LayoutType::Scalar2D:
        {
          const auto v_src = gsrc.get_view<const Real*,Host>();
          const auto v_tgt = gtgt.get_view<const Real*,Host>();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            Real expected = 0;
            auto gdof = gids_tgt_v(idof);
            for (int j=0; j<2; ++j) {
              auto src_gcol = gdof + j;
              auto src_lcol = gid2lid(src_gcol,gids_src_v);
              expected += w*v_src(src_lcol);
            }
            CHECK ( v_tgt(idof)== expected );
            ok &= catch_capture.lastAssertionPassed();
          }
        } break;
        case LayoutType::Vector2D:
        {
          const auto v_src = gsrc.get_view<const Real**,Host>();
          const auto v_tgt = gtgt.get_view<const Real**,Host>();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            for (int icmp=0; icmp<vec_dim; ++icmp) {
              Real expected = 0;
              auto gdof = gids_tgt_v(idof);
              for (int j=0; j<2; ++j) {
                auto src_gcol = gdof + j;
                auto src_lcol = gid2lid(src_gcol,gids_src_v);
                expected += w*v_src(src_lcol,icmp);
              }
              CHECK ( v_tgt(idof,icmp)== expected );
              ok &= catch_capture.lastAssertionPassed();
            }
          }
        } break;
        case LayoutType::Tensor2D:
        {
          const auto v_src = gsrc.get_view<const Real***,Host>();
          const auto v_tgt = gtgt.get_view<const Real***,Host>();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            for (int icmp=0; icmp<vec_dim; ++icmp) {
              for (int jcmp=0; jcmp<vec_dim; ++jcmp) {
                Real expected = 0;
                auto gdof = gids_tgt_v(idof);
                for (int j=0; j<2; ++j) {
                  auto src_gcol = gdof + j;
                  auto src_lcol = gid2lid(src_gcol,gids_src_v);
                  expected += w*v_src(src_lcol,icmp,jcmp);
                }
                CHECK ( v_tgt(idof,icmp,jcmp)== expected );
                ok &= catch_capture.lastAssertionPassed();
              }
            }
          }
        } break;
        case LayoutType::Scalar3D:
        {
          const auto v_src = gsrc.get_view<const Real**,Host>();
          const auto v_tgt = gtgt.get_view<const Real**,Host>();
          auto f_nlevs = gsrc.get_header().get_identifier().get_layout().dims().back();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            for (int ilev=0; ilev<f_nlevs; ++ilev) {
              Real expected = 0;
              auto gdof = gids_tgt_v(idof);
              for (int j=0; j<2; ++j) {
                auto src_gcol = gdof + j;
                auto src_lcol = gid2lid(src_gcol,gids_src_v);
                expected += w*v_src(src_lcol,ilev);
              }
              CHECK ( v_tgt(idof,ilev)== expected );
              ok &= catch_capture.lastAssertionPassed();
            }
          }
        } break;
        case LayoutType::Vector3D:
        {
          const auto v_src = gsrc.get_view<const Real***,Host>();
          const auto v_tgt = gtgt.get_view<const Real***,Host>();
          auto f_nlevs = gsrc.get_header().get_identifier().get_layout().dims().back();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            for (int icmp=0; icmp<vec_dim; ++icmp) {
              for (int ilev=0; ilev<f_nlevs; ++ilev) {
                Real expected = 0;
                auto gdof = gids_tgt_v(idof);
                for (int j=0; j<2; ++j) {
                  auto src_gcol = gdof + j;
                  auto src_lcol = gid2lid(src_gcol,gids_src_v);
                  expected += w*v_src(src_lcol,icmp,ilev);
                }
                CHECK ( v_tgt(idof,icmp,ilev)== expected );
                ok &= catch_capture.lastAssertionPassed();
              }
            }
          }
        } break;
        case LayoutType::Tensor3D:
        {
          const auto v_src = gsrc.get_view<const Real****,Host>();
          const auto v_tgt = gtgt.get_view<const Real****,Host>();
          auto f_nlevs = gsrc.get_header().get_identifier().get_layout().dims().back();
          for (int idof=0; idof<ngdofs_tgt; ++idof) {
            for (int icmp=0; icmp<tens_dim1; ++icmp) {
              for (int jcmp=0; jcmp<tens_dim2; ++jcmp) {
                for (int ilev=0; ilev<f_nlevs; ++ilev) {
                  Real expected = 0;
                  auto gdof = gids_tgt_v(idof);
                  for (int j=0; j<2; ++j) {
                    auto src_gcol = gdof + j;
                    auto src_lcol = gid2lid(src_gcol,gids_src_v);
                    expected += w*v_src(src_lcol,icmp,jcmp,ilev);
                  }
                  CHECK ( v_tgt(idof,icmp,jcmp,ilev)== expected );
                  ok &= catch_capture.lastAssertionPassed();
                }
              }
            }
          }
        } break;
        default:
          EKAT_ERROR_MSG ("Unexpected layout.\n");
      }
      root_print (msg + (ok ? "PASS" : "FAIL") + "\n",comm);
    }
  }

  // Clean up scorpio stuff
  scorpio::finalize_subsystem();
}

} // namespace scream
