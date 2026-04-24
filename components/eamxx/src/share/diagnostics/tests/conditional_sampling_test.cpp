#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <string>

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols, const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"point_grid"});
  auto &pl = gm_params.sublist("point_grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

const util::TimeStamp t0() {
  return util::TimeStamp({2024, 1, 1}, {0, 0, 0});
}

Field create_field (const std::string& name, int ncol, int ncmp, int nlev, const std::string& gn)
{
  using namespace ShortFieldTagsNames;
  FieldLayout fl{};
  if (ncol>0)
    fl.append_dim(COL,ncol);
  if (ncmp)
    fl.append_dim(CMP,ncmp);
  if (nlev)
    fl.append_dim(LEV,nlev);
  FieldIdentifier fid(name,fl,ekat::units::kg,gn);
  Field f(fid,true);
  f.get_header().get_tracking().update_time_stamp(t0());
  return f;
}

template<typename ST>
void iota_z (const Field& x) {
  int n1 = x.get_header().get_identifier().get_layout().dim(0);
  if (x.rank()>1) {
    for (int i=0; i<n1; ++i) {
      iota_z<ST>(x.subfield(0,i));
    }
  } else {
    auto vh = x.get_view<ST*,Host>();
    for (int i=0; i<n1; ++i) {
      vh(i) = i;
    }
    x.sync_to_dev();
  }
}

TEST_CASE("conditional_sampling") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  constexpr auto fv = constants::fill_value<Real>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 10;
  constexpr int ncols = 8;

  auto gm   = create_gm(comm, ncols, nlevs);
  auto grid = gm->get_grid("physics");

  // Random number generator seed
  int seed = get_random_test_seed();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  params.set("grid_name", grid->name());

  SECTION("exceptions") {
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params)); // No 'field_name' parameter
    params.set<std::string>("field_name", "x");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params)); // No 'condition_lhs' parameter
    params.set<std::string>("condition_lhs", "y");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params)); // No 'condition_cmp' parameter
    params.set<std::string>("condition_cmp", "z");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params)); // No 'condition_rhs'
  }

  SECTION("field_where_field_cmp_val") {
    const auto val = 0.5;
    params.set<std::string>("field_name", "x");
    params.set<std::string>("condition_lhs", "y");
    params.set<std::string>("condition_rhs", std::to_string(val));

    for (std::string cmp_s : {"ne", "gt"}) {
      for (int nx : {0,ncols}) {
        for (int ncmp : {0,2}) {
          for (int nz : {0,nlevs}) {
            auto x = create_field("x",nx,ncmp,nz,grid->name());
            auto y = create_field("y",nx,ncmp,nz,grid->name());
            randomize_uniform(x,seed++);
            randomize_uniform(y,seed++);
            params.set<std::string>("condition_cmp", cmp_s);
            auto diag = diag_factory.create("ConditionalSampling", comm, params);
            diag->set_grid(grid);
            diag->set_required_field(x);
            diag->set_required_field(y);
            diag->initialize();
            diag->compute_diagnostic();
            auto d = diag->get_diagnostic();

            REQUIRE (d.has_valid_mask());
            Field mask = d.get_valid_mask().clone();

            auto cmp = cmp_s=="ne" ? Comparison::NE : Comparison::GT;
            compute_mask(y,val,cmp,mask);
            if (not views_are_equal(d.get_valid_mask(),mask))
            {
              print_field_hyperslab(x);
              print_field_hyperslab(y);
              print_field_hyperslab(mask.alias("manual"));
              print_field_hyperslab(d.get_valid_mask().alias("computed"));
            }
            REQUIRE (views_are_equal(d.get_valid_mask(),mask));

            x.deep_copy(fv,mask,true);
            REQUIRE (views_are_equal(x,d));
          }
        }
      }
    }
  }

  SECTION("field_where_field_cmp_field") {
    params.set<std::string>("field_name", "x");
    params.set<std::string>("condition_lhs", "y");
    params.set<std::string>("condition_rhs", "z");

    for (std::string cmp_s : {"eq", "le"}) {
      for (int nx : {0,ncols}) {
        for (int ncmp : {0,2}) {
          for (int nz : {0,nlevs}) {
            auto x = create_field("x",nx,ncmp,nz,grid->name());
            auto y = create_field("y",nx,ncmp,nz,grid->name());
            auto z = create_field("z",nx,ncmp,nz,grid->name());
            randomize_uniform(x,seed++);
            randomize_uniform(y,seed++);
            randomize_uniform(z,seed++);
            params.set<std::string>("condition_cmp", cmp_s);
            auto diag = diag_factory.create("ConditionalSampling", comm, params);
            diag->set_grid(grid);
            diag->set_required_field(x);
            diag->set_required_field(y);
            diag->set_required_field(z);
            diag->initialize();
            diag->compute_diagnostic();
            auto d = diag->get_diagnostic();

            REQUIRE (d.has_valid_mask());
            Field mask = d.get_valid_mask().clone();

            auto cmp = cmp_s=="eq" ? Comparison::EQ : Comparison::LE;
            compute_mask(y,z,cmp,mask);
            REQUIRE (views_are_equal(d.get_valid_mask(),mask));

            x.deep_copy(fv,mask,true);
            REQUIRE (views_are_equal(x,d));
          }
        }
      }
    }
  }

  SECTION("field_where_lev_cmp_val") {
    const auto val = nlevs / 2;
    params.set<std::string>("field_name", "x");
    params.set<std::string>("condition_lhs", "lev");
    params.set<std::string>("condition_rhs", std::to_string(val));

    for (std::string cmp_s : {"ge", "le"}) {
      for (int nx : {0,ncols}) {
        for (int ncmp : {0,2}) {
          auto x = create_field("x",nx,ncmp,nlevs,grid->name());
          randomize_uniform(x,seed++);
          params.set<std::string>("condition_cmp", cmp_s);
          auto diag = diag_factory.create("ConditionalSampling", comm, params);
          diag->set_grid(grid);
          diag->set_required_field(x);
          diag->initialize();
          diag->compute_diagnostic();
          auto d = diag->get_diagnostic();

          Field mask = x.create_valid_mask();
          Field lev = mask.clone();
          iota_z<int>(lev);
          auto cmp = cmp_s=="ge" ? Comparison::GE : Comparison::LE;
          compute_mask(lev,val,cmp,mask);
          if (not views_are_equal(d.get_valid_mask(),mask))
          {
            print_field_hyperslab(lev.alias("lev"));
            print_field_hyperslab(mask.alias("manual"));
            print_field_hyperslab(d.get_valid_mask().alias("computed"));
          }
          REQUIRE (views_are_equal(d.get_valid_mask(),mask));

          x.deep_copy(fv,mask,true);
          REQUIRE (views_are_equal(x,d));
        }
      }
    }
  }

  SECTION("mask_where_lev_cmp_val") {
    const auto val = nlevs / 2;
    params.set<std::string>("field_name", "mask");
    params.set<std::string>("condition_lhs", "lev");
    params.set<std::string>("condition_rhs", std::to_string(val));

    for (std::string cmp_s : {"ge", "le"}) {
      params.set<std::string>("condition_cmp", cmp_s);
      auto diag = diag_factory.create("ConditionalSampling", comm, params);
      diag->set_grid(grid);
      diag->initialize();
      diag->compute_diagnostic();
      auto d = diag->get_diagnostic();

      REQUIRE (d.data_type()==DataType::IntType);
      REQUIRE (d.rank()==1);

      Field mask = d.clone();
      Field lev = mask.clone();
      iota_z<int>(lev);
      auto cmp = cmp_s=="ge" ? Comparison::GE : Comparison::LE;
      compute_mask(lev,val,cmp,mask);
      REQUIRE (views_are_equal(d,mask));
    }
  }

  SECTION("mask_where_field_cmp_field") {
    params.set<std::string>("field_name", "mask");
    params.set<std::string>("condition_lhs", "y");
    params.set<std::string>("condition_rhs", "z");

    for (std::string cmp_s : {"eq", "le"}) {
      for (int nx : {0,ncols}) {
        for (int ncmp : {0,2}) {
          for (int nz : {0,nlevs}) {
            auto y = create_field("y",nx,ncmp,nz,grid->name());
            auto z = create_field("z",nx,ncmp,nz,grid->name());
            randomize_uniform(y,seed++);
            randomize_uniform(z,seed++);
            params.set<std::string>("condition_cmp", cmp_s);
            auto diag = diag_factory.create("ConditionalSampling", comm, params);
            diag->set_grid(grid);
            diag->set_required_field(y);
            diag->set_required_field(z);
            diag->initialize();
            diag->compute_diagnostic();
            auto d = diag->get_diagnostic();

            REQUIRE (d.data_type()==DataType::IntType);
            REQUIRE (d.rank()==y.rank());

            Field mask = d.clone();
            auto cmp = cmp_s=="eq" ? Comparison::EQ : Comparison::LE;
            compute_mask(y,z,cmp,mask);
            REQUIRE (views_are_equal(d,mask));
          }
        }
      }
    }
  }

}

} // namespace scream
