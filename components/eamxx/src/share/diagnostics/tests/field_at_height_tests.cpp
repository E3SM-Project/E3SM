#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"

#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

void f_z_src(const Real y0, const Real m, const Field& z_data, Field& out_data);
// out_data must already have a valid mask (see Field::create_valid_mask):
// this sets both its data and its mask, matching HeightLevelIndex/
// FieldAtHeight's rule (see height_level_index.hpp): invalid above the
// top entry always, invalid below the bottom entry unless extrapolating.
void f_z_tgt(const Real y0, const Real m, const Real z_target, const Field& z_data,
            Field& out_data, const bool extrapolate_bottom);

TEST_CASE("field_at_height")
{
  using namespace ShortFieldTagsNames;

  register_diagnostics();

  // Get an MPI comm group for test
  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int nruns = 100;
  constexpr Real tol  = std::numeric_limits<Real>::epsilon()*1e5;

  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a pont grid
  int ncols = 3;
  int ndims = 4;
  int nlevs = 10;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  const auto m = ekat::units::m;
  // Create input data test fields
  FieldIdentifier s_mid_fid ("s_mid",FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier s_int_fid ("s_int",FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());
  FieldIdentifier v_mid_fid ("v_mid",FieldLayout({COL,CMP, LEV},{ncols,ndims,nlevs  }),m,grid->name());
  FieldIdentifier v_int_fid ("v_int",FieldLayout({COL,CMP,ILEV},{ncols,ndims,nlevs+1}),m,grid->name());
  // Create vertical fields z and geo on both midpoints and interfaces
  FieldIdentifier z_surf_fid ("z_surf",          FieldLayout({COL         },{ncols              }),m,grid->name());
  FieldIdentifier z_mid_fid   ("z_mid",           FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier z_int_fid   ("z_int",           FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());
  FieldIdentifier h_mid_fid ("height_mid",FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier h_int_fid ("height_int",FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());
  // Keep track of reference fields for comparison
  FieldIdentifier s_tgt_fid ("scalar_target",FieldLayout({COL    },{ncols      }),m,grid->name());
  FieldIdentifier v_tgt_fid ("vector_target",FieldLayout({COL,CMP},{ncols,ndims}),m,grid->name());

  Field s_mid   (s_mid_fid);
  Field s_int   (s_int_fid);
  Field v_mid   (v_mid_fid);
  Field v_int   (v_int_fid);
  Field z_surf  (z_surf_fid);
  Field z_mid   (z_mid_fid);
  Field z_int   (z_int_fid);
  Field h_mid (h_mid_fid);
  Field h_int (h_int_fid);
  Field s_tgt   (s_tgt_fid);
  Field v_tgt   (v_tgt_fid);

  s_mid.allocate_view();
  s_int.allocate_view();
  v_mid.allocate_view();
  v_int.allocate_view();
  z_surf.allocate_view();
  z_mid.allocate_view();
  z_int.allocate_view();
  h_mid.allocate_view();
  h_int.allocate_view();
  s_tgt.allocate_view();
  v_tgt.allocate_view();

  // FieldAtHeight's output always has a valid mask (see
  // field_at_height.cpp), so these hand-computed reference fields need
  // one too, for views_are_equal to compare them.
  s_tgt.create_valid_mask();
  v_tgt.create_valid_mask();

  s_mid.get_header().get_tracking().update_time_stamp(t0);
  s_int.get_header().get_tracking().update_time_stamp(t0);
  v_mid.get_header().get_tracking().update_time_stamp(t0);
  v_int.get_header().get_tracking().update_time_stamp(t0);
  z_surf.get_header().get_tracking().update_time_stamp(t0);
  z_mid.get_header().get_tracking().update_time_stamp(t0);
  z_int.get_header().get_tracking().update_time_stamp(t0);
  h_mid.get_header().get_tracking().update_time_stamp(t0);
  h_int.get_header().get_tracking().update_time_stamp(t0);
  s_tgt.get_header().get_tracking().update_time_stamp(t0);
  v_tgt.get_header().get_tracking().update_time_stamp(t0);

  auto print = [&](const std::string& msg) {
    if (comm.am_i_root()) {
      std::cout << msg;
    }
  };

  auto engine = scream::setup_random_test(&comm);
  using IPDF = std::uniform_int_distribution<int>;
  using RPDF = std::uniform_real_distribution<Real>;

  IPDF pdf_fields (0,1000);
  RPDF pdf_m  (1,10);
  RPDF pdf_y0 (0,5);

  // Lambda to create and run a diag, and return output
  auto run_diag = [&](const Field& f, const Field& z,
                      const double h, const std::string& surf_ref) {
    util::TimeStamp t0 ({2022,1,1},{0,0,0});
    auto& factory = DiagnosticFactory::instance();
    ekat::ParameterList pl;
    pl.set("surface_reference",surf_ref);
    pl.set("height_value",std::to_string(h));
    pl.set("height_units",std::string("m"));
    pl.set("field_name",f.name());
    pl.set("grid_name",grid->name());
    auto diag = factory.create("FieldAtheight",comm,pl, grid);
    diag->set_input_field(f);
    diag->set_input_field(z);

    // FieldAtHeight now depends on the (mid/int) bracket-index field
    // computed by HeightLevelIndex. Compute the variant matching z here,
    // and feed it in, since this test builds the diag by hand (rather
    // than through the io-stream machinery, which resolves such
    // dependencies automatically).
    const std::string layer = z.name().substr(z.name().size()-3)=="mid" ? "mid" : "int";
    ekat::ParameterList idx_pl;
    idx_pl.set("surface_reference",surf_ref);
    idx_pl.set("height_value",std::to_string(h));
    idx_pl.set("height_units",std::string("m"));
    idx_pl.set("vertical_layer",layer);
    auto idx_diag = factory.create("HeightLevelIndex",comm,idx_pl,grid);
    idx_diag->set_input_field(z);
    idx_diag->initialize();
    idx_diag->compute(t0);
    diag->set_input_field(idx_diag->get());

    diag->initialize();
    diag->compute(t0);
    diag->get().sync_to_host();
    return diag->get();
  };

  // Set up vertical structure for the tests.  Note,
  //   z_mid/int represents the height in m above sealevel
  //   h_mid/int represente the hegith in m above the surface
  // So we first construct z_mid/int using z_surf as reference, and
  // then can build h_mid/int from z_mid/int
  // Furthermore, z_mid is just the midpoint between two adjacent z_int
  // points, so we back z_mid out of z_int.
  //
  // To simplify the surface contribution we set the surface height to equal
  // the local column index.
  const int z_top = 1000;
  const Real surf_slope = z_top/10.0/ncols;
  const auto& zint_v   = z_int.get_view<Real**,Host>();
  const auto& zmid_v   = z_mid.get_view<Real**,Host>();
  const auto& zsurf_v  = z_surf.get_view<Real*,Host>();
  const auto& geoint_v = h_int.get_view<Real**,Host>();
  const auto& geomid_v = h_mid.get_view<Real**,Host>();
  int         min_col_thickness = z_top;
  int         max_surf = 0;
  for (int ii=0; ii<ncols; ++ii) {
    zsurf_v(ii) = ii*surf_slope;
    max_surf = zsurf_v(ii) > max_surf ? zsurf_v(ii) : max_surf;
    const Real col_thickness = z_top - zsurf_v(ii);
    min_col_thickness = min_col_thickness < col_thickness ? col_thickness : min_col_thickness;
    const Real dz = (z_top - zsurf_v(ii))/nlevs;
    zint_v(ii,0) = z_top;
    geoint_v(ii,0) = z_top - zsurf_v(ii); // Note, the distance above surface needs to consider the surface height.
    for (int jj=0; jj<nlevs; ++jj) {
      zint_v(ii,jj+1)   = zint_v(ii,jj)-dz;
      zmid_v(ii,jj)     = 0.5*(zint_v(ii,jj) + zint_v(ii,jj+1));
      geoint_v(ii,jj+1) = zint_v(ii,jj+1)- zsurf_v(ii); 
      geomid_v(ii,jj)   = zmid_v(ii,jj)  - zsurf_v(ii);
    }
  }
  z_mid.sync_to_dev();
  z_int.sync_to_dev();
  z_surf.sync_to_dev();
  h_int.sync_to_dev();
  h_mid.sync_to_dev();
  // Set the PDF for target height in the test to always be within the shortest column.
  // This ensures that we don't havea target z that extrapolates everywhere.
  // We test this case individually.
  IPDF pdf_levs (max_surf,min_col_thickness);
  // Sanity check that the geo and z vertical structures are in fact different,
  // so we know we are testing above_surface and above_sealevel as different cases.
  REQUIRE(! views_are_equal(z_int,h_int));
  REQUIRE(! views_are_equal(z_mid,h_mid));

  // Make sure that an unsupported reference height throws an error.
  print(" -> Testing throws error with unsupported reference height...\n");
  {
    REQUIRE_THROWS(run_diag (s_mid,h_mid,1.0,"foobar"));
  }
  print(" -> Testing throws error with unsupported reference height... OK\n");

  // Run many times
  int z_tgt;
  for (std::string surf_ref : {"sealevel","surface"}) {
    printf(" -> Testing for a reference height above %s...\n",surf_ref.c_str());
    const auto mid_src = surf_ref == "sealevel" ? z_mid : h_mid;
    const auto int_src = surf_ref == "sealevel" ? z_int : h_int;
    const int  max_surf_4test = surf_ref == "sealevel" ? max_surf : 0;
    // "above surface" extrapolates below the bottom entry; "above sealevel"
    // does not (see HeightLevelIndex). Either way, above the top entry is
    // always invalid, which f_z_tgt accounts for on its own.
    const bool extrapolate_bottom = surf_ref == "surface";
    for (int irun=0; irun<nruns; ++irun) {

      // Randomize fields using f_z_src function defined above:
      auto slope = pdf_m(engine); 
      auto inter = pdf_y0(engine);
      f_z_src(inter, slope, mid_src, s_mid);
      f_z_src(inter, slope, mid_src, v_mid);
      f_z_src(inter, slope, int_src, s_int);
      f_z_src(inter, slope, int_src, v_int);

      // Set target z-slice for testing to a random value.
      z_tgt = pdf_levs(engine)+max_surf_4test;
      printf("  -> test at height of %dm............\n",z_tgt);
      {
        print("    -> scalar midpoint field...............\n");
        auto d = run_diag(s_mid,mid_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,mid_src,s_tgt,extrapolate_bottom);
        REQUIRE (views_are_equal(d,s_tgt,tol));
        print("    -> scalar midpoint field............... OK!\n");
      }
      {
        print("    -> scalar interface field...............\n");
        auto d = run_diag (s_int,int_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,s_tgt,extrapolate_bottom);
        REQUIRE (views_are_equal(d,s_tgt,tol));
        print("    -> scalar interface field............... OK!\n");
      }
      {
        print("    -> vector midpoint field...............\n");
        auto d = run_diag (v_mid,mid_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,mid_src,v_tgt,extrapolate_bottom);
        REQUIRE (views_are_equal(d,v_tgt,tol));
        print("    -> vector midpoint field............... OK!\n");
      }
      {
        print("    -> vector interface field...............\n");
        auto d = run_diag (v_int,int_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,v_tgt,extrapolate_bottom);
        REQUIRE (views_are_equal(d,v_tgt,tol));
        print("    -> vector interface field............... OK!\n");
      }
      {
        print("    -> Forced fail, give incorrect location...............\n");
        const int z_tgt_adj = (z_tgt+max_surf_4test)/2;
        auto d = run_diag(s_int,int_src,z_tgt_adj,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,s_tgt,extrapolate_bottom);
        REQUIRE (!views_are_equal(d,s_tgt,tol));
        print("    -> Forced fail, give incorrect location............... OK!\n");
      }
    }
    // A target above the top entry is always flagged as invalid, regardless
    // of the surface reference (matches FieldAtPressureLevel at the model
    // top). A target below the bottom entry extrapolates for "above
    // surface" heights (always well defined near the surface), but is
    // flagged as invalid for "above sealevel" heights (e.g. a target
    // elevation below a mountain's surface is not well defined).
    auto check_invalid = [&](const Real z, const char* what) {
      print(std::string("    -> Forced out-of-range at ")+what+"...............\n");
      auto d = run_diag(s_int,int_src,z,surf_ref);
      REQUIRE (d.has_valid_mask());
      auto mask = d.get_valid_mask();
      mask.sync_to_host();
      auto mask_v = mask.get_view<const int*,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        REQUIRE (mask_v(icol)==0);
      }
      print(std::string("    -> Forced out-of-range at ")+what+"............... OK!\n");
    };
    check_invalid(2*z_top,"top");
    if (surf_ref=="surface") {
      print("    -> Forced extrapolation at bot...............\n");
      auto slope = pdf_m(engine);
      auto inter = pdf_y0(engine);
      f_z_src(inter, slope, int_src, s_int);
      z_tgt = 0;
      auto dbot = run_diag(s_int,int_src,z_tgt,surf_ref);
      f_z_tgt(inter,slope,z_tgt,int_src,s_tgt,extrapolate_bottom);
      REQUIRE (views_are_equal(dbot,s_tgt,tol));
      print("    -> Forced extrapolation at bot............... OK!\n");
    } else {
      // Strictly below every column's surface height (zsurf_v is in [0,(ncols-1)*surf_slope]).
      check_invalid(-1,"bot");
    }
    printf(" -> Testing for a reference height above %s... OK!\n",surf_ref.c_str());
  }
}

//-------------------------------
// Set up the inpute data.  To make the test simple we assume a linear distribution of the data
// with height.  That way we can exactly calculate what a linear interpolation to a random
// height would be.
void f_z_src(const Real y0, const Real m, const Field& z_data, Field& out_data) {
  using namespace ShortFieldTagsNames;
  const auto layout = out_data.get_header().get_identifier().get_layout();
  if (layout.has_tag(CMP)) { // Is a vector layout, meaning different dims than z_data.
    const auto& dims = layout.dims();
    const auto& z_view = z_data.get_view<const Real**,Host>();
    const auto& out_view = out_data.get_view<Real***,Host>();
    for (int ii=0; ii<dims[0]; ++ii) {
      for (int nd=0; nd<dims[1]; ++nd) {
        for (int jj=0; jj<dims[2]; ++jj) {
          out_view(ii,nd,jj) = y0 + m*(nd+1)*z_view(ii,jj);
        }
      }
    }
  } else { // Not a vector output, easier to deal with
    const auto z_view = z_data.get_internal_view_data<const Real,Host>();
    const auto& size = z_data.get_header().get_identifier().get_layout().size();
    auto out_view = out_data.get_internal_view_data<Real,Host>();
    for (int ii=0; ii<size; ++ii) {
      out_view[ii] = y0 + m*z_view[ii];
    }
  }
  out_data.sync_to_dev();
}
//-------------------------------
// Calculate the target data, AND the target mask (out_data must already
// have a valid mask). Note expression here must match the f_z_src above,
// and the validity rule must match HeightLevelIndex/FieldAtHeight: above
// the top entry is always invalid; below the bottom entry is invalid
// unless extrapolate_bottom.
void f_z_tgt(const Real y0, const Real m, const Real z_target, const Field& z_data,
            Field& out_data, const bool extrapolate_bottom) {
  using namespace ShortFieldTagsNames;
  const auto layout = out_data.get_header().get_identifier().get_layout();
  const auto& z_view = z_data.get_view<const Real**,Host>();
  const auto& zdims = z_data.get_header().get_identifier().get_layout().dims();
  auto mask = out_data.get_valid_mask();
  if (layout.has_tag(CMP)) { // Is a vector layout, meaning different dims than z_target.
    const auto& dims = layout.dims();
    const auto& out_view = out_data.get_view<Real**,Host>();
    const auto& mask_view = mask.get_view<int**,Host>();
    for (int ii=0; ii<dims[0]; ++ii) {
      // Check if FieldAtHeight would have had to extrapolate:
      const bool top_oob = z_target > z_view(ii,0);
      const bool bot_oob = z_target < z_view(ii,zdims[1]-1);
      const int valid = (not top_oob and (not bot_oob or extrapolate_bottom)) ? 1 : 0;
      for (int nd=0; nd<dims[1]; ++nd) {
        if (top_oob) {
          out_view(ii,nd) = y0 + m*(nd+1)*z_view(ii,0);
        } else if (bot_oob) {
          out_view(ii,nd) = y0 + m*(nd+1)*z_view(ii,zdims[1]-1);
        } else {
          out_view(ii,nd) = y0 + m*(nd+1)*z_target;
        }
        mask_view(ii,nd) = valid;
      }
    }
  } else { // Not a vector output, easier to deal with
    const auto& dims = layout.dims();
    const auto& out_view = out_data.get_view<Real*,Host>();
    const auto& mask_view = mask.get_view<int*,Host>();
    for (int ii=0; ii<dims[0]; ++ii) {
      // Check if FieldAtHeight would have had to extrapolate:
      const bool top_oob = z_target > z_view(ii,0);
      const bool bot_oob = z_target < z_view(ii,zdims[1]-1);
      if (top_oob) {
        out_view(ii) = y0 + m*z_view(ii,0);
      } else if (bot_oob) {
        out_view(ii) = y0 + m*z_view(ii,zdims[1]-1);
      } else {
        out_view(ii) = y0 + m*z_target;
      }
      mask_view(ii) = (not top_oob and (not bot_oob or extrapolate_bottom)) ? 1 : 0;
    }
  }
  out_data.sync_to_dev();
  mask.sync_to_dev();
}

} // namespace scream
