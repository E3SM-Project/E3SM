#include "catch2/catch.hpp"

#include "diagnostics/field_at_height.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

namespace scream {

void f_z_src(const Real y0, const Real m, const Field& z_data, Field& out_data);
void f_z_tgt(const Real y0, const Real m, const Real z_target, const Field& z_data, Field& out_data);
bool views_are_approx_equal(const Field& f0, const Field& f1, const Real tol, const bool msg = true);

TEST_CASE("field_at_height")
{
  using namespace ShortFieldTagsNames;

  register_diagnostics();

  // Get an MPI comm group for test
  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int nruns = 100;

  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager w/ a point grid
  constexpr Real tol            = std::numeric_limits<Real>::epsilon()*1e5;
  int ncols = 3;
  int ndims = 4;
  int nlevs = 10;
  int num_global_cols = ncols*comm.size();
  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,num_global_cols);
  gm->build_grids();
  auto grid = gm->get_grid("Point Grid");

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
    auto& factory = AtmosphereDiagnosticFactory::instance();
    ekat::ParameterList pl;
    pl.set("surface_reference",surf_ref);
    pl.set("height_value",std::to_string(h));
    pl.set("height_units",std::string("m"));
    pl.set("field_name",f.name());
    pl.set("grid_name",grid->name());
    auto diag = factory.create("FieldAtheight",comm,pl);
    diag->set_grids(gm);
    diag->set_required_field(f);
    diag->set_required_field(z);
    diag->initialize(t0,RunType::Initial);
    diag->compute_diagnostic();
    diag->get_diagnostic().sync_to_host();
    return diag->get_diagnostic();
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
        f_z_tgt(inter,slope,z_tgt,mid_src,s_tgt);
        REQUIRE (views_are_approx_equal(d,s_tgt,tol));
        print("    -> scalar midpoint field............... OK!\n");
      }
      {
        print("    -> scalar interface field...............\n");
        auto d = run_diag (s_int,int_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,s_tgt);
        REQUIRE (views_are_approx_equal(d,s_tgt,tol));
        print("    -> scalar interface field............... OK!\n");
      }
      {
        print("    -> vector midpoint field...............\n");
        auto d = run_diag (v_mid,mid_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,mid_src,v_tgt);
        REQUIRE (views_are_approx_equal(d,v_tgt,tol));
        print("    -> vector midpoint field............... OK!\n");
      }
      {
        print("    -> vector interface field...............\n");
        auto d = run_diag (v_int,int_src,z_tgt,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,v_tgt);
        REQUIRE (views_are_approx_equal(d,v_tgt,tol));
        print("    -> vector interface field............... OK!\n");
      }
      {
        print("    -> Forced fail, give incorrect location...............\n");
        const int z_tgt_adj = (z_tgt+max_surf_4test)/2;
        auto d = run_diag(s_int,int_src,z_tgt_adj,surf_ref);
        f_z_tgt(inter,slope,z_tgt,int_src,s_tgt);
        REQUIRE (!views_are_approx_equal(d,s_tgt,tol,false));
        print("    -> Forced fail, give incorrect location............... OK!\n");
      }
    }
    {
      print("    -> Forced extrapolation at top...............\n");
      auto slope = pdf_m(engine); 
      auto inter = pdf_y0(engine);
      f_z_src(inter, slope, int_src, s_int);
      z_tgt = 2*z_top;
      auto dtop = run_diag(s_int,int_src,z_tgt,surf_ref);
      f_z_tgt(inter,slope,z_tgt,int_src,s_tgt);
      REQUIRE (views_are_approx_equal(dtop,s_tgt,tol));
      print("    -> Forced extrapolation at top............... OK!\n");
      print("    -> Forced extrapolation at bot...............\n");
      z_tgt = 0;
      auto dbot = run_diag(s_int,int_src,z_tgt,surf_ref);
      f_z_tgt(inter,slope,z_tgt,int_src,s_tgt);
      REQUIRE (views_are_approx_equal(dbot,s_tgt,tol));
      print("    -> Forced extrapolation at bot............... OK!\n");
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
// Calculate the target data.  Note expression here must match the f_z_src abovel
void f_z_tgt(const Real y0, const Real m, const Real z_target, const Field& z_data, Field& out_data) {
  using namespace ShortFieldTagsNames;
  const auto layout = out_data.get_header().get_identifier().get_layout();
  const auto& z_view = z_data.get_view<const Real**,Host>();
  const auto& zdims = z_data.get_header().get_identifier().get_layout().dims();
  if (layout.has_tag(CMP)) { // Is a vector layout, meaning different dims than z_target.
    const auto& dims = layout.dims();
    const auto& out_view = out_data.get_view<Real**,Host>();
    for (int ii=0; ii<dims[0]; ++ii) {
      for (int nd=0; nd<dims[1]; ++nd) {
        // Check if FieldAtHeight would have had to extrapolate:
	if (z_target > z_view(ii,0)) {
          out_view(ii,nd) = y0 + m*(nd+1)*z_view(ii,0);
	} else if ( z_target < z_view(ii,zdims[1]-1)) {
          out_view(ii,nd) = y0 + m*(nd+1)*z_view(ii,zdims[1]-1);
	} else {
          out_view(ii,nd) = y0 + m*(nd+1)*z_target;
	}
      }
    }
  } else { // Not a vector output, easier to deal with
    const auto& dims = layout.dims();
    const auto& out_view = out_data.get_view<Real*,Host>();
    for (int ii=0; ii<dims[0]; ++ii) {
      // Check if FieldAtHeight would have had to extrapolate:
      if (z_target > z_view(ii,0)) {
        out_view(ii) = y0 + m*z_view(ii,0);
      } else if ( z_target < z_view(ii,zdims[1]-1)) {
        out_view(ii) = y0 + m*z_view(ii,zdims[1]-1);
      } else {
        out_view(ii) = y0 + m*z_target;
      }
    }
  }
  out_data.sync_to_dev();
}
/*-----------------------------------------------------------------------------------------------*/
bool views_are_approx_equal(const Field& f0, const Field& f1, const Real tol, const bool msg)
{
  const auto& l0    = f0.get_header().get_identifier().get_layout();
  const auto& l1    = f1.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(l0==l1,"Error! views_are_approx_equal - the two fields don't have matching layouts.");
  // Take advantage of field utils update, min and max to assess the max difference between the two fields
  // simply.
  auto ft = f0.clone();
  ft.update(f1,1.0,-1.0);
  auto d_min = field_min<Real>(ft);
  auto d_max = field_max<Real>(ft);
  if (std::abs(d_min) > tol or std::abs(d_max) > tol) {
    if (msg) {
      printf("The two copies of (%16s) are NOT approx equal within a tolerance of %e.\n     The min and max errors are %e and %e respectively.\n",f0.name().c_str(),tol,d_min,d_max);
    }
    return false;
  } else {
    return true;
  }

}

} // namespace scream
