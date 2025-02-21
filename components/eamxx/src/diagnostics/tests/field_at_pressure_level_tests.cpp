#include "catch2/catch.hpp"

#include "ekat/ekat_pack_utils.hpp"

#include "diagnostics/field_at_pressure_level.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

namespace scream {

const int packsize = SCREAM_SMALL_PACK_SIZE;
using Pack         = ekat::Pack<Real,packsize>;

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,num_global_cols);
  gm->build_grids();

  return gm;
}

template<typename T>
bool approx(const T a, const T b) {
  constexpr Real eps = std::numeric_limits<Real>::epsilon();
  if (std::abs(a-b)>eps) {
    printf("approx violated with std::abs(%e - %e) = %e > %e\n",a,b,std::abs(a-b),eps);
  }
  return std::abs(a-b)<=eps; 
}

struct PressureBnds
{
  Real p_top  = 10000.0;  //  100mb
  Real p_surf = 100000.0; // 1000mb
};

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<FieldAtPressureLevel>
get_test_diag(const ekat::Comm& comm, std::shared_ptr<const FieldManager> fm, std::shared_ptr<const GridsManager> gm, const std::string& type, const Real plevel);

Real get_test_pres(const int col, const int lev, const int num_lev, const int num_cols);
Real get_test_data(const Real pres);

TEST_CASE("field_at_pressure_level_p2")
{
  // Get an MPI comm group for test
  ekat::Comm comm(MPI_COMM_WORLD);

  // For any random case, set a number of iterations of that test
  int num_checks = 10; 

  // Create a grids manager w/ a point grid
  int ncols = 3;
  int nlevs = 10;
  auto gm   = create_gm(comm,ncols,nlevs);

  // Create a field manager for testing
  auto grid = gm->get_grid("Point Grid");
  auto fm   = get_test_fm(grid);
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Random values to be used in test
  PressureBnds pressure_bounds;
  auto engine = scream::setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  Real p_mid_bnds_dz = pressure_bounds.p_surf/nlevs;  // Note, p_mid will actually never make it to p_top or p_surf in all columns, so we offset the bounds.
  RPDF pdf_pmid(pressure_bounds.p_top+p_mid_bnds_dz,pressure_bounds.p_surf-p_mid_bnds_dz);
  RPDF pdf_pint(pressure_bounds.p_top,pressure_bounds.p_surf);
 
  // Get a diagnostic factory
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("FieldAtPressureLevel",&create_atmosphere_diagnostic<FieldAtPressureLevel>);

  {
    // Test 1: Take a slice at a random value for variable defined at midpoint.
    for (int test_itr=0;test_itr<num_checks;test_itr++) {
      Real plevel = std::round(pdf_pmid(engine));
      auto diag = get_test_diag(comm, fm, gm, "mid", plevel);
      diag->initialize(t0,RunType::Initial);
      diag->compute_diagnostic();
      auto diag_f = diag->get_diagnostic();
      diag_f.sync_to_host();
      auto test1_diag_v = diag_f.get_view<const Real*, Host>();
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(approx(test1_diag_v(icol),get_test_data(plevel)));
      }
    }
  } 
  {
    // Test 2: Take a slice at a random value for variable defined at interface.
    for (int test_itr=0;test_itr<num_checks;test_itr++) {
      Real plevel = std::round(pdf_pint(engine));
      auto diag = get_test_diag(comm, fm, gm, "int", plevel);
      diag->initialize(t0,RunType::Initial);
      diag->compute_diagnostic();
      auto diag_f = diag->get_diagnostic();
      diag_f.sync_to_host();
      auto test2_diag_v = diag_f.get_view<const Real*, Host>();
      // Check the mask field inside the diag_f
      auto mask_f = diag_f.get_header().get_extra_data<Field>("mask_data");
      mask_f.sync_to_host();
      auto test2_mask_v = mask_f.get_view<const Real*, Host>();
      //
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(approx(test2_diag_v(icol),get_test_data(plevel)));
        REQUIRE(approx(test2_mask_v(icol),Real(1.0)));
      }
    }
  } 
  {
    // Test 3: Take a slice at a value outside the bounds, which should return the default masked value
    for (int test_itr=0;test_itr<num_checks;test_itr++) {
      Real plevel = pressure_bounds.p_surf*2;
      auto diag = get_test_diag(comm, fm, gm, "int", plevel);
      diag->initialize(t0,RunType::Initial);
      diag->compute_diagnostic();
      auto diag_f = diag->get_diagnostic();
      diag_f.sync_to_host();
      auto test2_diag_v = diag_f.get_view<const Real*, Host>();
      // Check the mask field inside the diag_f
      auto mask_f = diag_f.get_header().get_extra_data<Field>("mask_data");
      mask_f.sync_to_host();
      auto test2_mask_v = mask_f.get_view<const Real*, Host>();
      auto mask_val = diag_f.get_header().get_extra_data<Real>("mask_value");
      //
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(approx(test2_diag_v(icol),Real(mask_val)));
        REQUIRE(approx(test2_mask_v(icol),Real(0.0)));
      }
    }
  } 
  
} // TEST_CASE("field_at_pressure_level")
/*==========================================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs  = grid->get_num_vertical_levels();

  // Create some fields for this fm
  std::vector<FieldTag> tag_mid = {COL,LEV};
  std::vector<FieldTag> tag_int = {COL,ILEV};

  std::vector<int>     dims_mid = {num_lcols,num_levs};
  std::vector<int>     dims_int = {num_lcols,num_levs+1};

  const std::string& gn = grid->name();

  FieldIdentifier fid1("V_mid",FL{tag_mid,dims_mid},m,gn);
  FieldIdentifier fid2("V_int",FL{tag_int,dims_int},kg,gn);
  FieldIdentifier fid3("p_mid",FL{tag_mid,dims_mid},Pa,gn);
  FieldIdentifier fid4("p_int",FL{tag_int,dims_int},Pa,gn);

  // Register fields with fm
  // Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
  fm->registration_begins();
  fm->register_field(FR{fid1,Pack::n});
  fm->register_field(FR{fid2,Pack::n});
  fm->register_field(FR{fid3,Pack::n});
  fm->register_field(FR{fid4,Pack::n});
  fm->registration_ends();

  // Initialize these fields
  auto f1 = fm->get_field(fid1);
  auto f2 = fm->get_field(fid2);
  auto f3 = fm->get_field(fid3);
  auto f4 = fm->get_field(fid4);
  auto f1_host = f1.get_view<Pack**,Host>();
  auto f2_host = f2.get_view<Pack**,Host>();
  auto f3_host = f3.get_view<Pack**,Host>();
  auto f4_host = f4.get_view<Pack**,Host>();

  for (int ii=0;ii<num_lcols;++ii) {
    for (int jj=0;jj<num_levs;++jj) {
      int ipack = jj / packsize;
      int ivec  = jj % packsize;
      Real p_mid = (get_test_pres(ii,jj,num_levs,num_lcols) + get_test_pres(ii,jj+1,num_levs,num_lcols)) / 2;
      Real p_int = get_test_pres(ii,jj,num_levs,num_lcols);
      f1_host(ii,ipack)[ivec] = get_test_data(p_mid); 
      f2_host(ii,ipack)[ivec] = get_test_data(p_int); 
      f3_host(ii,ipack)[ivec] = p_mid; 
      f4_host(ii,ipack)[ivec] = p_int;
    }
    // The interface values have one more vertical level
    int ipack = num_levs / packsize;
    int ivec  = num_levs % packsize;
    Real p_int = get_test_pres(ii,num_levs,num_levs,num_lcols);
    f2_host(ii,ipack)[ivec] = get_test_data(p_int); 
    f4_host(ii,ipack)[ivec] = p_int;
  }
  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  // Sync back to device
  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
  f4.sync_to_dev();

  return fm;
}
/*===================================================================================================*/
std::shared_ptr<FieldAtPressureLevel>
get_test_diag(const ekat::Comm& comm, std::shared_ptr<const FieldManager> fm, std::shared_ptr<const GridsManager> gm, const std::string& type, const Real plevel)
{
    std::string fname = "V_"+type;
    auto field = fm->get_field(fname);
    auto fid = field.get_header().get_identifier();
    ekat::ParameterList params;
    params.set("field_name",field.name());
    params.set("grid_name",fm->get_grid()->name());
    params.set("pressure_value",std::to_string(plevel));
    params.set("pressure_units",std::string("Pa"));
    auto diag = std::make_shared<FieldAtPressureLevel>(comm,params);
    diag->set_grids(gm);
    for (const auto& req : diag->get_required_field_requests()) {
      auto req_field = fm->get_field(req.fid);
      diag->set_required_field(req_field);
    }
    return diag;
}
/*===================================================================================================*/
Real get_test_pres(const int col, const int lev, const int num_lev, const int num_cols)
{
  PressureBnds pressure_bounds;
  Real p_surf = pressure_bounds.p_surf; 
  Real p_top  = pressure_bounds.p_top;  
  Real dp_dx  = p_top/(num_cols);        // Make sure that 1 column hits the min pressure of 0 at tom
  Real dp_dz  = (p_surf-(p_top-(col+1)*dp_dx))/(num_lev); // Make sure the max pressure at surface is the surf pressure
  return (p_top - (col+1)*dp_dx) + lev*dp_dz;
}

Real get_test_data(const Real pres)
{
  // Have the data follow the linear curve
  //   y = 100 + p, where p is the pressure
  return 100.0 + pres;
}
/*===================================================================================================*/
} // namespace scream
