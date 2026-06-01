#include "catch2/catch.hpp"

#include "share/diagnostics/field_at_pressure_level.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>

namespace scream {

const int packsize = SCREAM_PACK_SIZE;
using Pack         = ekat::Pack<Real,packsize>;


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

std::map<std::string,Field>
create_fields(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<FieldAtPressureLevel>
get_test_diag(const std::map<std::string,Field>& fields,
              std::shared_ptr<const AbstractGrid> grid,
              const std::string& type, const Real plevel);

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
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Create a field manager for testing
  auto fmap   = create_fields(grid);
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Random values to be used in test
  PressureBnds pressure_bounds;
  auto engine = scream::setup_random_test(&comm);
  using RPDF = std::uniform_real_distribution<Real>;
  Real p_mid_bnds_dz = pressure_bounds.p_surf/nlevs;  // Note, p_mid will actually never make it to p_top or p_surf in all columns, so we offset the bounds.
  RPDF pdf_pmid(pressure_bounds.p_top+p_mid_bnds_dz,pressure_bounds.p_surf-p_mid_bnds_dz);
  RPDF pdf_pint(pressure_bounds.p_top,pressure_bounds.p_surf);
 
  // Get a diagnostic factory
  auto& diag_factory = DiagnosticFactory::instance();
  diag_factory.register_product("FieldAtPressureLevel",&create_diagnostic<FieldAtPressureLevel>);

  {
    // Test 1: Take a slice at a random value for variable defined at midpoint.
    for (int test_itr=0;test_itr<num_checks;test_itr++) {
      Real plevel = std::round(pdf_pmid(engine));
      auto diag = get_test_diag(fmap, grid, "mid", plevel);
      diag->initialize();
      diag->compute(t0);
      auto diag_f = diag->get();
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
      auto diag = get_test_diag(fmap, grid, "int", plevel);
      diag->initialize();
      diag->compute(t0);
      auto diag_f = diag->get();
      diag_f.sync_to_host();
      auto test2_diag_v = diag_f.get_view<const Real*, Host>();
      // Check the mask field inside the diag_f
      auto mask_f = diag_f.get_valid_mask();
      mask_f.sync_to_host();
      auto test2_mask_v = mask_f.get_view<const int*, Host>();
      //
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(approx(test2_diag_v(icol),get_test_data(plevel)));
        REQUIRE(test2_mask_v(icol)==1);
      }
    }
  } 
  {
    // Test 3: Take a slice at a value outside the bounds, which should return the default masked value
    for (int test_itr=0;test_itr<num_checks;test_itr++) {
      Real plevel = pressure_bounds.p_surf*2;
      auto diag = get_test_diag(fmap, grid, "int", plevel);
      diag->initialize();
      diag->compute(t0);
      auto diag_f = diag->get();
      diag_f.sync_to_host();
      auto test2_diag_v = diag_f.get_view<const Real*, Host>();
      // Check the mask field inside the diag_f
      auto mask_f = diag_f.get_valid_mask();
      mask_f.sync_to_host();
      auto test2_mask_v = mask_f.get_view<const int*, Host>();

      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(approx(test2_diag_v(icol),constants::fill_value<Real>));
        REQUIRE(test2_mask_v(icol)==0);
      }
    }
  } 
  
} // TEST_CASE("field_at_pressure_level")
/*==========================================================================================================*/
std::map<std::string,Field>
create_fields(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  std::map<std::string,Field> fields;

  // Create both midpoints and interface fields
  auto layout_mid = grid->get_3d_scalar_layout(LEV);
  auto layout_int = grid->get_3d_scalar_layout(ILEV);
  const std::string& gn = grid->name();

  FieldIdentifier fid1("V_mid",layout_mid,m,gn);
  FieldIdentifier fid2("V_int",layout_int,kg,gn);
  FieldIdentifier fid3("p_mid",layout_mid,Pa,gn);
  FieldIdentifier fid4("p_int",layout_int,Pa,gn);

  auto& f1 = fields[fid1.name()] = Field(fid1,true);
  auto& f2 = fields[fid2.name()] = Field(fid2,true);
  auto& f3 = fields[fid3.name()] = Field(fid3,true);
  auto& f4 = fields[fid4.name()] = Field(fid4,true);

  // Initialize these fields
  auto f1_host = f1.get_view<Real**,Host>();
  auto f2_host = f2.get_view<Real**,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Real**,Host>();

  const int ncols = grid->get_num_local_dofs();
  const int nlevs = grid->get_num_vertical_levels();
  for (int icol=0;icol<ncols;++icol) {
    for (int ilev=0;ilev<nlevs;++ilev) {
      Real p_mid = (get_test_pres(icol,ilev,nlevs,ncols) + get_test_pres(icol,ilev+1,nlevs,ncols)) / 2;
      Real p_int = get_test_pres(icol,ilev,nlevs,ncols);
      f1_host(icol,ilev) = get_test_data(p_mid); 
      f2_host(icol,ilev) = get_test_data(p_int); 
      f3_host(icol,ilev) = p_mid; 
      f4_host(icol,ilev) = p_int;
    }
    // The interface values have one more vertical level
    Real p_int = get_test_pres(icol,nlevs,nlevs,ncols);
    f2_host(icol,nlevs) = get_test_data(p_int); 
    f4_host(icol,nlevs) = p_int;
  }
  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  f1.get_header().get_tracking().update_time_stamp(time);
  f2.get_header().get_tracking().update_time_stamp(time);
  f3.get_header().get_tracking().update_time_stamp(time);
  f4.get_header().get_tracking().update_time_stamp(time);

  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
  f4.sync_to_dev();

  return fields;
}
/*===================================================================================================*/
std::shared_ptr<FieldAtPressureLevel>
get_test_diag(const std::map<std::string,Field>& fields,
              std::shared_ptr<const AbstractGrid> grid,
              const std::string& type, const Real plevel)
{
  std::string fname = "V_"+type;
  ekat::ParameterList params;
  params.set("field_name",fname);
  params.set("pressure_value",std::to_string(plevel));
  params.set("pressure_units",std::string("Pa"));
  const auto& comm = grid->get_comm();
  auto diag = std::make_shared<FieldAtPressureLevel>(comm,params,grid);
  for (const auto& fname : diag->get_input_fields_names()) {
    auto f = fields.at(fname);
    diag->set_input_field(f);
  }
  return diag;
}

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
