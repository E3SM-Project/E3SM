#include "catch2/catch.hpp"

#include "diagnostics/field_at_height.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {

Field get_ref_field(const int lev, const Field& src_data, const std::string& surf_ref);
Field get_ref_field(const int levm, const int levp, const Field& src_data, const std::string& surf_ref);

TEST_CASE("field_at_height")
{
  using namespace ShortFieldTagsNames;

  register_diagnostics();

  // Get an MPI comm group for test
  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int nruns = 10;

  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager w/ a point grid
  int ncols = 3;
  int ndims = 4;
  int nlevs = 10;
  int num_global_cols = ncols*comm.size();
  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,num_global_cols);
  gm->build_grids();
  auto grid = gm->get_grid("Point Grid");

  // Create input test fields, as well as z_mid/int fields
  const auto m = ekat::units::m;
  FieldIdentifier s_mid_fid ("s_mid",FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier s_int_fid ("s_int",FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());
  FieldIdentifier v_mid_fid ("v_mid",FieldLayout({COL,CMP, LEV},{ncols,ndims,nlevs  }),m,grid->name());
  FieldIdentifier v_int_fid ("v_int",FieldLayout({COL,CMP,ILEV},{ncols,ndims,nlevs+1}),m,grid->name());
  FieldIdentifier z_mid_fid ("z_mid",FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier z_int_fid ("z_int",FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());
  FieldIdentifier geo_mid_fid ("geopotential_mid",FieldLayout({COL,     LEV},{ncols,      nlevs  }),m,grid->name());
  FieldIdentifier geo_int_fid ("geopotential_int",FieldLayout({COL,    ILEV},{ncols,      nlevs+1}),m,grid->name());

  Field s_mid (s_mid_fid);
  Field s_int (s_int_fid);
  Field v_mid (v_mid_fid);
  Field v_int (v_int_fid);
  Field z_mid (z_mid_fid);
  Field z_int (z_int_fid);
  Field geo_mid (geo_mid_fid);
  Field geo_int (geo_int_fid);

  s_mid.allocate_view();
  s_int.allocate_view();
  v_mid.allocate_view();
  v_int.allocate_view();
  z_mid.allocate_view();
  z_int.allocate_view();
  geo_mid.allocate_view();
  geo_int.allocate_view();

  s_mid.get_header().get_tracking().update_time_stamp(t0);
  s_int.get_header().get_tracking().update_time_stamp(t0);
  v_mid.get_header().get_tracking().update_time_stamp(t0);
  v_int.get_header().get_tracking().update_time_stamp(t0);
  z_mid.get_header().get_tracking().update_time_stamp(t0);
  z_int.get_header().get_tracking().update_time_stamp(t0);
  geo_mid.get_header().get_tracking().update_time_stamp(t0);
  geo_int.get_header().get_tracking().update_time_stamp(t0);

  auto print = [&](const std::string& msg) {
    if (comm.am_i_root()) {
      std::cout << msg;
    }
  };

  auto engine = scream::setup_random_test(&comm);
  using IPDF = std::uniform_int_distribution<int>;

  IPDF pdf_fields (0,1000);
  IPDF pdf_levs (1,nlevs-1);

  // Lambda to create and run a diag, and return output
  auto run_diag = [&](const Field& f, const Field& z,
                      const std::string& loc, const std::string& surf_ref) {
    util::TimeStamp t0 ({2022,1,1},{0,0,0});
    auto& factory = AtmosphereDiagnosticFactory::instance();
    ekat::ParameterList pl;
    pl.set<std::string>("surface_reference",surf_ref);
    pl.set("vertical_location",loc);
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

  // Create z(i,j)=nlevs-j, which makes testing easier
  for (auto f : {z_mid, z_int}) {
    auto v = f.get_view<Real**,Host>();
    const auto& dims = f.get_header().get_identifier().get_layout().dims();
    for (int i=0; i<dims[0]; ++i) {
      for (int j=0; j<dims[1]; ++j) {
        v(i,j) = dims[1]-j+i;
      }
    };
    f.sync_to_dev();
  }
  // Create geo(i,j)=nlevs-j + i, which makes testing easier
  for (auto f : {geo_mid, geo_int}) {
    auto v = f.get_view<Real**,Host>();
    const auto& dims = f.get_header().get_identifier().get_layout().dims();
    for (int i=0; i<dims[0]; ++i) {
      for (int j=0; j<dims[1]; ++j) {
        v(i,j) = dims[1]-j;
      }
    };
    f.sync_to_dev();
  }

  // Run many times
  Real z_tgt,lev_tgt;
  std::string loc;
  for (int irun=0; irun<nruns; ++irun) {

    // Randomize fields
    for (auto f : {s_mid, s_int, v_mid, v_int}) {
      // We know f is not a subfield, so we can use this power-user method
      auto data = f.get_internal_view_data<Real,Host>();
      const auto& size = f.get_header().get_identifier().get_layout().size();
      for (int i=0; i<size; ++i) {
        data[i] = pdf_fields(engine);
      }
      f.sync_to_dev();
    }

    // Make sure that an unsupported reference height throws an error.
    print(" -> Testing throws error with unsupported reference height...\n");
    {
      REQUIRE_THROWS(run_diag (s_mid,geo_mid,loc,"foobar"));
    }
    print(" -> Testing throws error with unsupported reference height... OK\n");

    // Tests for different reference cases:
    for (std::string surf_ref : {"surface", "sealevel"}) {
      z_tgt = pdf_levs(engine);
      loc = std::to_string(z_tgt) + "m";
      // z[ilev] = nlevs-ilev, so the tgt slice is nlevs-z_tgt
      lev_tgt = nlevs-z_tgt;

      print(" -> Testing with z_tgt above " + surf_ref + " coinciding with a z level\n");
      const Field ref_mid = surf_ref == "surface" ? geo_mid : z_mid;
      const Field ref_int = surf_ref == "surface" ? geo_int : z_int;
      {
        print("    -> scalar midpoint field...............\n");
        auto d = run_diag (s_mid,ref_mid,loc,surf_ref);
        auto tgt = get_ref_field(lev_tgt,s_mid,surf_ref);
        REQUIRE (views_are_equal(d,tgt,&comm));
        print("    -> scalar midpoint field............... OK!\n");
      }
      {
        print("    -> scalar interface field...............\n");
        auto d = run_diag (s_int,ref_int,loc,surf_ref);
        // geo_int = nlevs+1-ilev, so the tgt slice is nlevs+1-z_tgt
        // z_int = nlevs+1-ilev+i, so the tgt slice is nlevs+1-z_tgt+i
        auto tgt = get_ref_field(lev_tgt+1,s_int,surf_ref);
        REQUIRE (views_are_equal(d,tgt,&comm));
        print("    -> scalar interface field............... OK!\n");
      }
      {
        print("    -> vector midpoint field...............\n");
        auto d = run_diag (v_mid,ref_mid,loc,surf_ref);
        // We can't subview over 3rd index and keep layout right,
        // so do all cols separately
        for (int i=0; i<ncols; ++i) {
          auto fi = v_mid.subfield(0,i);
          auto di = d.subfield(0,i);
	  // Because we are already cycling through cols we don't need
	  // get_ref_field to adjust for above sealevel, we do it locally.
	  auto lev_adj = surf_ref == "surface" ? lev_tgt :
		  ( i + lev_tgt >= nlevs ? nlevs-1 : i + lev_tgt );
          auto tgt = get_ref_field(lev_adj,fi,"");
          REQUIRE (views_are_equal(di,tgt,&comm));
        }
        print("    -> vector midpoint field............... OK!\n");
      }
      {
        print("    -> vector interface field...............\n");
        auto d = run_diag (v_int,ref_int,loc,surf_ref);
        // We can't subview over 3rd index and keep layout right,
        // so do all cols separately
        for (int i=0; i<ncols; ++i) {
          auto fi = v_int.subfield(0,i);
          auto di = d.subfield(0,i);
	  // Because we are already cycling through cols we don't need
	  // get_ref_field to adjust for above sealevel, we do it locally.
	  auto lev_adj = surf_ref == "surface" ? lev_tgt :
		  ( i + lev_tgt >= nlevs ? nlevs-1 : i + lev_tgt );
          auto tgt = get_ref_field(lev_adj+1,fi,"");
          REQUIRE (views_are_equal(di,tgt,&comm));
        }
        print("    -> vector interface field............... OK!\n");
      }

      z_tgt = pdf_levs(engine) + 0.5;
      lev_tgt = nlevs-z_tgt;
      loc = std::to_string(z_tgt) + "m";

      auto zp1 = static_cast<int>(std::round(lev_tgt+0.5));
      auto zm1 = static_cast<int>(std::round(lev_tgt-0.5));

      print(" -> Testing with z_tgt above " + surf_ref + " between levels\n");
      {
        print("    -> scalar midpoint field...............\n");
        auto d = run_diag (s_mid,ref_mid,loc,surf_ref);
        auto tgt = get_ref_field(zm1,zp1,s_mid,surf_ref);
        REQUIRE (views_are_equal(d,tgt,&comm));
        print("    -> scalar midpoint field............... OK!\n");
      }
      {
        print("    -> scalar interface field...............\n");
        auto d = run_diag (s_int,ref_int,loc,surf_ref);
        auto tgt = get_ref_field(zm1+1,zp1+1,s_int,surf_ref);
        REQUIRE (views_are_equal(d,tgt,&comm));
        print("    -> scalar interface field............... OK!\n");
      }
      {
        print("    -> vector midpoint field...............\n");
        auto d = run_diag (v_mid,ref_mid,loc,surf_ref);
        // We can't subview over 3rd index and keep layout right,
        // so do all cols separately
        for (int i=0; i<ncols; ++i) {
          auto fi = v_mid.subfield(0,i);
          auto di = d.subfield(0,i);
	  // Because we are already cycling through cols we don't need
	  // get_ref_field to adjust for above sealevel, we do it locally.
	  auto zm1_adj = surf_ref == "surface" ? zm1 :
		  ( i + zm1 >= nlevs ? nlevs-1 : i + zm1 );
	  auto zp1_adj = surf_ref == "surface" ? zp1 :
		  ( i + zp1 >= nlevs ? nlevs-1 : i + zp1 );
          auto tgt = get_ref_field(zm1_adj,zp1_adj,fi,"");
          REQUIRE (views_are_equal(di,tgt,&comm));
        }
        print("    -> vector midpoint field............... OK!\n");
      }
      {
        print("    -> vector interface field...............\n");
        auto d = run_diag (v_int,ref_int,loc,surf_ref);
        // We can't subview over 3rd index and keep layout right,
        // so do all cols separately
        for (int i=0; i<ncols; ++i) {
          auto fi = v_int.subfield(0,i);
          auto di = d.subfield(0,i);
	  // Because we are already cycling through cols we don't need
	  // get_ref_field to adjust for above sealevel, we do it locally.
	  auto zm1_adj = surf_ref == "surface" ? zm1 :
		  ( i + zm1 >= nlevs ? nlevs-1 : i + zm1 );
	  auto zp1_adj = surf_ref == "surface" ? zp1 :
		  ( i + zp1 >= nlevs ? nlevs-1 : i + zp1 );
          auto tgt = get_ref_field(zm1_adj+1,zp1_adj+1,fi,"");
          REQUIRE (views_are_equal(di,tgt,&comm));
        }
        print("    -> vector interface field............... OK!\n");
      }
    }
  }
}

//-----------------
Field get_ref_field(const int lev, const Field& src_data, const std::string& surf_ref)
{
  // A single level implies this is a test where the height matches a single level.
  if (surf_ref == "sealevel") {
    auto tgt = src_data.subfield(1,lev).clone();
    const auto& tgt_view = tgt.get_view<Real*,Host>();
    const auto dims = src_data.get_header().get_identifier().get_layout().dims();
    const auto src_view = src_data.get_view<const Real**,Host>();
    for (int i = 0; i < dims[0]; ++i) {
      int ilev = lev + i >= dims[1] ? dims[1]-1 : lev + i;
      tgt_view(i) = src_view(i,ilev);
    }
    return tgt;
  } else { 
    return src_data.subfield(1,lev);
  }
}
//-----------------
Field get_ref_field(const int levm, const int levp, const Field& src_data, const std::string& surf_ref)
{
  // Two levels implies this is a test where the height is midway between two levels
  if (surf_ref == "sealevel") {
    auto tgtp = src_data.subfield(1,levp).clone();
    auto tgtm = src_data.subfield(1,levm).clone();

    const auto& tgtp_view = tgtp.get_view<Real*,Host>();
    const auto& tgtm_view = tgtm.get_view<Real*,Host>();

    const auto dims = src_data.get_header().get_identifier().get_layout().dims();
    const auto src_view = src_data.get_view<const Real**,Host>();
    for (int i = 0; i < dims[0]; ++i) {
      int ilevp = levp + i >= dims[1] ? dims[1]-1 : levp + i;
      int ilevm = levm + i >= dims[1] ? dims[1]-1 : levm + i;
      tgtp_view(i) = src_view(i,ilevp);
      tgtm_view(i) = src_view(i,ilevm);
    }
    tgtp.update(tgtm,0.5,0.5);
    return tgtp;
  } else {
    auto tgt = src_data.subfield(1,levp).clone();
    tgt.update(src_data.subfield(1,levm),0.5,0.5);
    return tgt;
  }
}

} // namespace scream
