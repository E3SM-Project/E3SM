#include <catch2/catch.hpp>

#include "share/grid/remap/iop_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

void root_print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

constexpr int vec_dim = 2;
constexpr int tens_dim1 = 3;
constexpr int tens_dim2 = 4;

template<typename Engine, typename PDF>
Field create_field (const std::string& name,
                    const LayoutType lt,
                    const AbstractGrid& grid,
                    const bool midpoints,
                    Engine& engine, PDF& pdf)
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

  randomize(f,engine,pdf);

  return f;
}

TEST_CASE("iop_remap")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  root_print ("\n +---------------------------------+\n",comm);
  root_print (" |       Testing iop remapper      |\n",comm);
  root_print (" +---------------------------------+\n\n",comm);

  using IPDF = std::uniform_int_distribution<int>;
  using RPDF = std::uniform_real_distribution<Real>;

  auto engine = setup_random_test (&comm);
  IPDF ipdf (2*comm.size(),10*comm.size());
  RPDF rpdf (0,1);

  // -------------------------------------- //
  //      Build src grid and remapper       //
  // -------------------------------------- //

  int src_ngcols = ipdf(engine);
  int tgt_ngcols = ipdf(engine);
  comm.broadcast(&src_ngcols,1,0);
  comm.broadcast(&tgt_ngcols,1,0);

  auto src_grid = create_point_grid("src", src_ngcols,10,comm);
  auto tgt_grid = create_point_grid("tgt", tgt_ngcols,10,comm);
  int src_nlcols = src_grid->get_num_local_dofs();
  int tgt_nlcols = tgt_grid->get_num_local_dofs();

  REQUIRE_THROWS (std::make_shared<IOPRemapper>(src_grid,tgt_grid,0,0)); // not lat/lon in src

  auto bad_src = create_point_grid("src", src_ngcols,10,comm);
  bad_src->create_geometry_data("lat",src_grid->get_3d_scalar_layout(true));
  bad_src->create_geometry_data("lon",src_grid->get_3d_scalar_layout(true));
  REQUIRE_THROWS (std::make_shared<IOPRemapper>(bad_src,tgt_grid,0,0)); // lat/lon bad layout

  auto src_lat = src_grid->create_geometry_data("lat",src_grid->get_2d_scalar_layout());
  auto src_lon = src_grid->create_geometry_data("lon",src_grid->get_2d_scalar_layout());
  randomize(src_lat,engine,rpdf);
  randomize(src_lon,engine,rpdf);

  auto se_grid = std::make_shared<SEGrid>("se",5,4,10,comm);
  se_grid->create_geometry_data("lat",se_grid->get_2d_scalar_layout());
  se_grid->create_geometry_data("lon",se_grid->get_2d_scalar_layout());
  REQUIRE_THROWS (std::make_shared<IOPRemapper>(src_grid,se_grid,0,0)); // tgt!=PointGrid
  REQUIRE_THROWS (std::make_shared<IOPRemapper>(se_grid,tgt_grid,0,0)); // src!=PointGrid

  auto bad_tgt = tgt_grid->clone("bad_tgt",true);
  bad_tgt->reset_num_vertical_lev(12);
  REQUIRE_THROWS (std::make_shared<IOPRemapper>(src_grid,bad_tgt,0,0)); // nlevs doesn't match

  REQUIRE_THROWS (std::make_shared<IOPRemapper>(src_grid,bad_tgt,91,0)); // lat OOB
  REQUIRE_THROWS (std::make_shared<IOPRemapper>(src_grid,bad_tgt,0,-1)); // lat OOB

  // Finally a working remapper
  int closest_rank = IPDF(0,comm.size()-1)(engine);
  comm.broadcast(&closest_rank,1,0);
  int closest_lid = IPDF(0,src_nlcols-1)(engine);
  Real lat,lon;
  if (closest_rank==comm.rank()) {
    lat = src_lat.get_view<const Real*,Host>()[closest_lid];
    lon = src_lon.get_view<const Real*,Host>()[closest_lid];
  }
  comm.broadcast(&lat,1,closest_rank);
  comm.broadcast(&lon,1,closest_rank);

  auto remap = std::make_shared<IOPRemapper>(src_grid,tgt_grid,lat,lon);

  // -------------------------------------- //
  //     Create fields and register them    //
  // -------------------------------------- //

  auto layouts = {
    LayoutType::Scalar2D,
    LayoutType::Vector2D,
    LayoutType::Tensor2D,
    LayoutType::Scalar3D,
    LayoutType::Scalar3D,
    LayoutType::Vector3D,
    LayoutType::Vector3D,
    LayoutType::Tensor3D,
    LayoutType::Tensor3D
  };

  remap->registration_begins();

  bool midpoints = false; // midpoints is unused for 2d layouts
  for (auto l : layouts) {
    auto n = e2str(l);
    auto src = create_field(n+"_src",l,*src_grid,midpoints,engine,rpdf);
    auto tgt = create_field(n+"_tgt",l,*tgt_grid,midpoints,engine,rpdf);
    remap->register_field(src,tgt);

    midpoints = not midpoints; // ensure we create both mid and int fields
  }
  remap->registration_ends();

  // -------------------------------------- //
  //        Remap fields and check          //
  // -------------------------------------- //

  REQUIRE_THROWS (remap->remap_bwd()); // NO bwd remap
  remap->remap_fwd();

  using namespace ShortFieldTagsNames;
  for (int i=0; i<remap->get_num_fields(); ++i) {
    const auto& src = remap->get_src_field(i);
    const auto& tgt = remap->get_tgt_field(i);

    auto col = src.subfield(COL,0).clone("col");
    int col_size = col.get_header().get_identifier().get_layout().size();
    if (comm.rank()==closest_rank) {
      col.deep_copy(src.subfield(COL,closest_lid));
    }
#if SCREAM_MPI_ON_DEVICE
    comm.broadcast(col.get_internal_view_data<Real>(),col_size,closest_rank);
#else
    col.sync_to_host();
    comm.broadcast(col.get_internal_view_data<Real,Host>(),col_size,closest_rank);
    col.sync_to_dev();
#endif

    for (int icol=0; icol<tgt_nlcols; ++icol) {
      REQUIRE (views_are_equal(col,tgt.subfield(COL,icol)));
    }
  }
}

} // namespace scream
