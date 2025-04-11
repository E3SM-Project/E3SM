#include <catch2/catch.hpp>

#include "share/grid/remap/refining_remapper_rma.hpp"
#include "share/grid/point_grid.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

class RefiningRemapperRMATester : public RefiningRemapperRMA {
public:
  RefiningRemapperRMATester (const grid_ptr_type& tgt_grid,
                          const std::string& map_file)
   : RefiningRemapperRMA(tgt_grid,map_file) {}

  ~RefiningRemapperRMATester () = default;

  void test_internals () {
    // Test arrays size
    const size_t n = m_num_fields;
    REQUIRE (m_src_fields.size()==n);
    REQUIRE (m_ov_fields.size()==n);
    REQUIRE (m_tgt_fields.size()==n);
    REQUIRE (m_col_size.size()==n);
    REQUIRE (m_col_stride.size()==n);
    REQUIRE (m_col_offset.size()==n);
    REQUIRE (m_mpi_win.size()==n);
    REQUIRE (m_remote_lids.size()==static_cast<size_t>(m_ov_coarse_grid->get_num_local_dofs()));
    REQUIRE (m_remote_pids.size()==static_cast<size_t>(m_ov_coarse_grid->get_num_local_dofs()));

    // Test field specs
    constexpr auto COL = ShortFieldTagsNames::COL;
    for (int i=0; i<m_num_fields; ++i) {
      const auto& fh = m_src_fields[i].get_header();
      const auto& fl = fh.get_identifier().get_layout();
      const auto& fap = fh.get_alloc_properties();
      const auto col_alloc_size = fap.get_num_scalars() / fl.dim(COL);
      REQUIRE (m_col_size[i]==fl.clone().strip_dim(COL).size());
      if (fh.get_parent()) {
        REQUIRE (m_col_stride[i]==col_alloc_size*fap.get_subview_info().dim_extent);
        REQUIRE (m_col_offset[i]==col_alloc_size*fap.get_subview_info().slice_idx);
      } else {
        REQUIRE (m_col_stride[i]==col_alloc_size);
        REQUIRE (m_col_offset[i]==0);
      }
    }

    // Test CRS matirx
    int nldofs_tgt = m_tgt_grid->get_num_local_dofs();
    int ngdofs_src = m_src_grid->get_num_global_dofs();
    REQUIRE (m_row_offsets.extent_int(0)==nldofs_tgt+1);
    auto row_offsets_h = cmvdc(m_row_offsets);
    auto col_lids_h = cmvdc(m_col_lids);
    auto weights_h = cmvdc(m_weights);
    auto col_gids_h = m_ov_coarse_grid->get_dofs_gids().get_view<const AbstractGrid::gid_type*,Host>();

    auto row_gids_h = m_tgt_grid->get_dofs_gids().get_view<const AbstractGrid::gid_type*,Host>();
    for (int i=0; i<nldofs_tgt; ++i) {
      auto row_gid = row_gids_h[i];
      auto beg = row_offsets_h(i);
      auto end = row_offsets_h(i+1);
      REQUIRE (end>=beg);
      for (int j=beg; j<end; ++j) {
        if (row_gid<ngdofs_src) {
          REQUIRE (weights_h(j)==1);
        } else {
          REQUIRE (weights_h(j)==0.5);
        }
      }
    }
  }
};

Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid)
{
  const auto u = ekat::units::Units::nondimensional();
  const auto& gn = grid.name();
  const auto  ndims = 2;
  Field f;
  switch (lt) {
    case LayoutType::Scalar2D:
      f = Field(FieldIdentifier(name,grid.get_2d_scalar_layout(),u,gn));  break;
    case LayoutType::Vector2D:
      f = Field(FieldIdentifier(name,grid.get_2d_vector_layout(ndims),u,gn));  break;
    case LayoutType::Scalar3D:
      f = Field(FieldIdentifier(name,grid.get_3d_scalar_layout(true),u,gn));  break;
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    case LayoutType::Vector3D:
      f = Field(FieldIdentifier(name,grid.get_3d_vector_layout(false,ndims),u,gn));  break;
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    default:
      EKAT_ERROR_MSG ("Invalid layout type for this unit test.\n");
  }
  f.allocate_view();

  return f;
}

template<typename Engine>
Field create_field (const std::string& name, const LayoutType lt, const AbstractGrid& grid, Engine& engine) {
  auto f = create_field(name,lt,grid);

  // Use discrete_distribution to get an integer, then use that as exponent for 2^-n.
  // This guarantees numbers that are exactly represented as FP numbers, which ensures
  // the test will produce the expected answer, regardless of how math ops are performed.
  using IPDF = std::discrete_distribution<int>;
  IPDF ipdf ({1,1,1,1,1,1,1,1,1,1});
  auto pdf = [&](Engine& e) {
    return std::pow(2,ipdf(e));
  };
  randomize(f,engine,pdf);

  return f;
}

Field all_gather_field (const Field& f, const ekat::Comm& comm) {
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
  std::vector<Real> data_vec(col_size);
  for (int pid=0,offset=0; pid<comm.size(); ++pid) {
    Real* data;
    int ncols = fl.dims()[0];
    comm.broadcast(&ncols,1,pid);
    for (int icol=0; icol<ncols; ++icol,offset+=col_size) {
      switch (fl.rank()) {
        case 1:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<Real*,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        case 2:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<Real**,Host>(),icol).data();
          } else {
            data = data_vec.data();
          }
          break;
        case 3:
          if (pid==comm.rank()) {
            data = ekat::subview(f.get_view<Real***,Host>(),icol).data();
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
      auto gdata = gf.get_internal_view_data<Real,Host>()+offset;
      std::copy(data,data+col_size,gdata);
    }
  }
  return gf;
}

void write_map_file (const std::string& filename, const int ngdofs_src) {
  // Add a dof in the middle of two coarse dofs
  const int ngdofs_tgt = 2*ngdofs_src-1;

  // Existing dofs are "copied", added dofs are averaged from neighbors
  const int nnz = ngdofs_src + 2*(ngdofs_src-1);

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::define_dim(filename, "n_a", ngdofs_src);
  scorpio::define_dim(filename, "n_b", ngdofs_tgt);
  scorpio::define_dim(filename, "n_s", nnz);

  scorpio::define_var(filename, "col", {"n_s"}, "int");
  scorpio::define_var(filename, "row", {"n_s"}, "int");
  scorpio::define_var(filename, "S",   {"n_s"}, "double");

  scorpio::enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz);
  for (int i=0; i<ngdofs_src; ++i) {
    col[i] = i;
    row[i] = i;
      S[i] = 1.0;
  }
  for (int i=0; i<ngdofs_src-1; ++i) {
    col[ngdofs_src+2*i] = i;
    row[ngdofs_src+2*i] = ngdofs_src+i;
      S[ngdofs_src+2*i] = 0.5;

    col[ngdofs_src+2*i+1] = i+1;
    row[ngdofs_src+2*i+1] = ngdofs_src+i;
      S[ngdofs_src+2*i+1] = 0.5;
  }

  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",  S.data());

  scorpio::release_file(filename);
}

TEST_CASE ("refining_remapper") {
  using gid_type = AbstractGrid::gid_type;

  auto& catch_capture = Catch::getResultCapture();

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test (&comm);

  scorpio::init_subsystem(comm);

  // Create a map file
  const int ngdofs_src = 4*comm.size();
  const int ngdofs_tgt = 2*ngdofs_src-1;
  auto filename = "rr_rma_tests_map.np" + std::to_string(comm.size()) + ".nc";
  write_map_file(filename,ngdofs_src);

  // Create target grid. Ensure gids are numbered like in map file
  const int nlevs = std::max(SCREAM_PACK_SIZE,16);
  auto tgt_grid = create_point_grid("tgt",ngdofs_tgt,nlevs,comm);
  auto dofs_h = tgt_grid->get_dofs_gids().get_view<gid_type*,Host>();
  for (int i=0; i<tgt_grid->get_num_local_dofs(); ++i) {
    int q = dofs_h[i] / 2;
    if (dofs_h[i] % 2 == 0) {
      dofs_h[i] = q;
    } else {
      dofs_h[i] = ngdofs_src + q;
    }
  }
  tgt_grid->get_dofs_gids().sync_to_dev();

  // Test bad registrations separately, since they corrupt the remapper state for later
  {
    auto r = std::make_shared<RefiningRemapperRMATester>(tgt_grid,filename);
    auto src_grid = r->get_src_grid();
    r->registration_begins();
    Field bad_src(FieldIdentifier("",src_grid->get_2d_scalar_layout(),ekat::units::m,src_grid->name(),DataType::IntType));
    Field bad_tgt(FieldIdentifier("",tgt_grid->get_2d_scalar_layout(),ekat::units::m,tgt_grid->name(),DataType::IntType));
    CHECK_THROWS (r->register_field(bad_src,bad_tgt)); // not allocated
    bad_src.allocate_view();
    bad_tgt.allocate_view();
    r->register_field(bad_src,bad_tgt);
    CHECK_THROWS (r->registration_ends()); // bad data type (must be real)
  }

  auto r = std::make_shared<RefiningRemapperRMATester>(tgt_grid,filename);
  auto src_grid = r->get_src_grid();

  auto bundle_src = create_field("bundle3d_src",LayoutType::Vector3D,*src_grid,engine);
  auto s2d_src   = create_field("s2d_src",LayoutType::Scalar2D,*src_grid,engine);
  auto v2d_src   = create_field("v2d_src",LayoutType::Vector2D,*src_grid,engine);
  auto s3d_src   = create_field("s3d_src",LayoutType::Scalar3D,*src_grid,engine);
  auto v3d_src   = create_field("v3d_src",LayoutType::Vector3D,*src_grid,engine);

  auto bundle_tgt = create_field("bundle3d_tgt",LayoutType::Vector3D,*tgt_grid);
  auto s2d_tgt   = create_field("s2d_tgt",LayoutType::Scalar2D,*tgt_grid);
  auto v2d_tgt   = create_field("v2d_tgt",LayoutType::Vector2D,*tgt_grid);
  auto s3d_tgt   = create_field("s3d_tgt",LayoutType::Scalar3D,*tgt_grid);
  auto v3d_tgt   = create_field("v3d_tgt",LayoutType::Vector3D,*tgt_grid);

  r->registration_begins();
  r->register_field(s2d_src,s2d_tgt);
  r->register_field(v2d_src,v2d_tgt);
  r->register_field(s3d_src,s3d_tgt);
  r->register_field(v3d_src,v3d_tgt);
  r->register_field(bundle_src.get_component(0),bundle_tgt.get_component(0));
  r->register_field(bundle_src.get_component(1),bundle_tgt.get_component(1));
  r->registration_ends();

  // Test remapper internal state
  r->test_internals();

  // Run remap
  CHECK_THROWS (r->remap(false)); // No backward remap
  r->remap(true);

  // Gather global copies (to make checks easier) and check src/tgt fields
  auto gs2d_src = all_gather_field(s2d_src,comm);
  auto gv2d_src = all_gather_field(v2d_src,comm);
  auto gs3d_src = all_gather_field(s3d_src,comm);
  auto gv3d_src = all_gather_field(v3d_src,comm);
  auto gbundle_src = all_gather_field(bundle_src,comm);

  auto gs2d_tgt = all_gather_field(s2d_tgt,comm);
  auto gv2d_tgt = all_gather_field(v2d_tgt,comm);
  auto gs3d_tgt = all_gather_field(s3d_tgt,comm);
  auto gv3d_tgt = all_gather_field(v3d_tgt,comm);
  auto gbundle_tgt = all_gather_field(bundle_tgt,comm);

  Real avg;
  // Scalar 2D
  {
    if (comm.am_i_root()) {
      printf(" -> Checking 2d scalars .........\n");
    }
    bool ok = true;
    gs2d_src.sync_to_host();
    gs2d_tgt.sync_to_host();

    auto src_v = gs2d_src.get_view<const Real*,Host>();
    auto tgt_v = gs2d_tgt.get_view<const Real*,Host>();

    // Coarse grid cols are just copied
    for (int icol=0; icol<ngdofs_src; ++icol) {
      CHECK (tgt_v[2*icol]==src_v[icol]);
      ok &= catch_capture.lastAssertionPassed();
    }
    // Fine cols are an average of the two cols nearby
    for (int icol=0; icol<ngdofs_src-1; ++icol) {
      avg = (src_v[icol] + src_v[icol+1]) / 2;
      CHECK (tgt_v[2*icol+1]==avg);
      ok &= catch_capture.lastAssertionPassed();
    }
    if (comm.am_i_root()) {
      printf(" -> Checking 2d scalars ......... %s\n",ok ? "PASS" : "FAIL");
    }
  }

  // Vector 2D
  {
    if (comm.am_i_root()) {
      printf(" -> Checking 2d vectors .........\n");
    }
    bool ok = true;
    gv2d_src.sync_to_host();
    gv2d_tgt.sync_to_host();

    auto src_v = gv2d_src.get_view<const Real**,Host>();
    auto tgt_v = gv2d_tgt.get_view<const Real**,Host>();

    // Coarse grid cols are just copied
    for (int icol=0; icol<ngdofs_src; ++icol) {
      for (int icmp=0; icmp<2; ++icmp) {
        CHECK (tgt_v(2*icol,icmp)==src_v(icol,icmp));
        ok &= catch_capture.lastAssertionPassed();
      }
    }
    // Fine cols are an average of the two cols nearby
    for (int icol=0; icol<ngdofs_src-1; ++icol) {
      for (int icmp=0; icmp<2; ++icmp) {
        avg = (src_v(icol,icmp) + src_v(icol+1,icmp)) / 2;
        CHECK (tgt_v(2*icol+1,icmp)==avg);
        ok &= catch_capture.lastAssertionPassed();
      }
    }
    if (comm.am_i_root()) {
      printf(" -> Checking 2d vectors ......... %s\n",ok ? "PASS" : "FAIL");
    }  
  }

  // Scalar 3D
  {
    if (comm.am_i_root()) {
      printf(" -> Checking 3d scalars .........\n");
    }
    bool ok = true;
    gs3d_src.sync_to_host();
    gs3d_tgt.sync_to_host();

    auto src_v = gs3d_src.get_view<const Real**,Host>();
    auto tgt_v = gs3d_tgt.get_view<const Real**,Host>();

    // Coarse grid cols are just copied
    for (int icol=0; icol<ngdofs_src; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        CHECK (tgt_v(2*icol,ilev)==src_v(icol,ilev));
        ok &= catch_capture.lastAssertionPassed();
      }
    }
    // Fine cols are an average of the two cols nearby
    for (int icol=0; icol<ngdofs_src-1; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        avg = (src_v(icol,ilev) + src_v(icol+1,ilev)) / 2;
        CHECK (tgt_v(2*icol+1,ilev)==avg);
        ok &= catch_capture.lastAssertionPassed();
      }
    }
    if (comm.am_i_root()) {
      printf(" -> Checking 3d scalars ......... %s\n",ok ? "PASS" : "FAIL");
    }
  }

  // Vector 3D
  {
    if (comm.am_i_root()) {
      printf(" -> Checking 3d vectors .........\n");
    }
    bool ok = true;
    gv3d_src.sync_to_host();
    gv3d_tgt.sync_to_host();

    auto src_v = gv3d_src.get_view<const Real***,Host>();
    auto tgt_v = gv3d_tgt.get_view<const Real***,Host>();

    // Coarse grid cols are just copied
    for (int icol=0; icol<ngdofs_src; ++icol) {
      for (int icmp=0; icmp<2; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          CHECK (tgt_v(2*icol,icmp,ilev)==src_v(icol,icmp,ilev));
          ok &= catch_capture.lastAssertionPassed();
        }
      }
    }
    // Fine cols are an average of the two cols nearby
    for (int icol=0; icol<ngdofs_src-1; ++icol) {
      for (int icmp=0; icmp<2; ++icmp) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          avg = (src_v(icol,icmp,ilev) + src_v(icol+1,icmp,ilev)) / 2;
          CHECK (tgt_v(2*icol+1,icmp,ilev)==avg);
          ok &= catch_capture.lastAssertionPassed();
        }
      }
    }
    if (comm.am_i_root()) {
      printf(" -> Checking 3d vectors ......... %s\n",ok ? "PASS" : "FAIL");
    }
  }

  // Subfields
  {
    if (comm.am_i_root()) {
      printf(" -> Checking 3d subfields .......\n");
    }
    bool ok = true;
    gbundle_src.sync_to_host();
    gbundle_tgt.sync_to_host();

    for (int icmp=0; icmp<2; ++icmp) {
      auto sf_src = gbundle_src.get_component(icmp);
      auto sf_tgt = gbundle_tgt.get_component(icmp);

      auto src_v = sf_src.get_view<const Real**,Host>();
      auto tgt_v = sf_tgt.get_view<const Real**,Host>();

      // Coarse grid cols are just copied
      for (int icol=0; icol<ngdofs_src; ++icol) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          CHECK (tgt_v(2*icol,ilev)==src_v(icol,ilev));
          ok &= catch_capture.lastAssertionPassed();
        }
      }
      // Fine cols are an average of the two cols nearby
      for (int icol=0; icol<ngdofs_src-1; ++icol) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          avg = (src_v(icol,ilev) + src_v(icol+1,ilev)) / 2;
          CHECK (tgt_v(2*icol+1,ilev)==avg);
          ok &= catch_capture.lastAssertionPassed();
        }
      }
    }
    if (comm.am_i_root()) {
      printf(" -> Checking 3d subfields ....... %s\n",ok ? "PASS" : "FAIL");
    }
  }

  // Clean up
  r = nullptr;
  scorpio::finalize_subsystem();
}

} // namespace scream
