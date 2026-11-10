#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <string>

namespace scream {


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
  const int nlevs = 10;
  constexpr int ncols = 8;

  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Random number generator seed
  int seed = get_random_test_seed();

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  params.set("grid_name", grid->name());

  SECTION("exceptions") {
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params, grid)); // No 'field_name' parameter
    params.set<std::string>("field_name", "x");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params, grid)); // No 'condition_lhs' parameter
    params.set<std::string>("condition_lhs", "y");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params, grid)); // No 'condition_cmp' parameter
    params.set<std::string>("condition_cmp", "z");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm, params, grid)); // No 'condition_rhs'
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
            auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
            diag->set_input_field(x);
            diag->set_input_field(y);
            diag->initialize();
            diag->compute(t0());
            auto d = diag->get();

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
            auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
            diag->set_input_field(x);
            diag->set_input_field(y);
            diag->set_input_field(z);
            diag->initialize();
            diag->compute(t0());
            auto d = diag->get();

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
          auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
          diag->set_input_field(x);
          diag->initialize();
          diag->compute(t0());
          auto d = diag->get();

          Field mask = x.create_valid_mask();
          Field lev = mask.clone();
          iota_z<int>(lev);
          auto cmp = cmp_s=="ge" ? Comparison::GE : Comparison::LE;
          compute_mask(lev,val,cmp,mask);
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
      auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
      diag->initialize();
      diag->compute(t0());
      auto d = diag->get();

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
            auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
            diag->set_input_field(y);
            diag->set_input_field(z);
            diag->initialize();
            diag->compute(t0());
            auto d = diag->get();

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

  SECTION("chained") {
    // x.where(y>0.5).where(C): the mask is the intersection of the two conditions
    auto x = create_field("x",ncols,0,nlevs,grid->name());
    auto y = create_field("y",ncols,0,nlevs,grid->name());
    auto z = create_field("z",ncols,0,nlevs,grid->name());
    randomize_uniform(x,seed++);

    params.set<std::string>("field_name", "x");
    params.set<std::string>("condition_lhs", "y");
    params.set<std::string>("condition_cmp", "gt");
    params.set<std::string>("condition_rhs", "0.5");
    auto diag1 = diag_factory.create("ConditionalSampling", comm, params, grid);
    diag1->set_input_field(x);
    diag1->set_input_field(y);
    diag1->initialize();
    auto d1 = diag1->get();

    // The second condition is either on a field, or on the level index
    int step = 0;
    for (std::string lhs2 : {"z", "lev"}) {
      params.set<std::string>("field_name", d1.name());
      params.set<std::string>("condition_lhs", lhs2);
      params.set<std::string>("condition_cmp", "gt");
      params.set<std::string>("condition_rhs", lhs2=="lev" ? "3" : "0.5");
      auto diag2 = diag_factory.create("ConditionalSampling", comm, params, grid);
      diag2->set_input_field(d1);
      if (lhs2!="lev") {
        diag2->set_input_field(z);
      }
      diag2->initialize();

      // Compute twice, with different data: the mask must not remember the first time
      for (int n=0; n<2; ++n) {
        ++step;
        randomize_uniform(y,seed++);
        randomize_uniform(z,seed++);
        for (auto f : {x,y,z}) {
          f.get_header().get_tracking().update_time_stamp(t0()+step);
        }
        diag1->compute(t0()+step);
        diag2->compute(t0()+step);
        auto d2 = diag2->get();

        Field mask1 = d2.get_valid_mask().clone("mask1");
        Field mask2 = d2.get_valid_mask().clone("mask2");
        compute_mask(y,0.5,Comparison::GT,mask1);
        if (lhs2=="lev") {
          auto m2_h = mask2.get_view<int**,Host>();
          for (int i=0; i<ncols; ++i) {
            for (int k=0; k<nlevs; ++k) {
              m2_h(i,k) = k>3 ? 1 : 0;
            }
          }
          mask2.sync_to_dev();
        } else {
          compute_mask(z,0.5,Comparison::GT,mask2);
        }
        mask1.scale(mask2);
        REQUIRE (views_are_equal(d2.get_valid_mask(),mask1));

        // Where the mask is 0, the diag is fill_value. Compare entries one by one,
        // since views_are_equal skips the entries that d2 masks out
        d2.sync_to_host();
        x.sync_to_host();
        mask1.sync_to_host();
        auto d2_h = d2.get_view<const Real**,Host>();
        auto x_h  = x.get_view<const Real**,Host>();
        auto m1_h = mask1.get_view<const int**,Host>();
        for (int i=0; i<ncols; ++i) {
          for (int k=0; k<nlevs; ++k) {
            REQUIRE (d2_h(i,k)==(m1_h(i,k)==1 ? x_h(i,k) : fv));
          }
        }
      }
    }
  }

  SECTION("units") {
    // Comparing two fields with different units is an error...
    auto x = create_field("x",ncols,0,nlevs,grid->name());
    FieldIdentifier y_fid("y",x.get_header().get_identifier().get_layout(),K, grid->name());
    FieldIdentifier z_fid("z",x.get_header().get_identifier().get_layout(),Pa,grid->name());
    Field y(y_fid,true), z(z_fid,true);
    for (auto f : {y,z}) {
      f.get_header().get_tracking().update_time_stamp(t0());
    }
    params.set<std::string>("field_name", "x");
    params.set<std::string>("condition_lhs", "y");
    params.set<std::string>("condition_cmp", "gt");
    params.set<std::string>("condition_rhs", "z");
    auto diag = diag_factory.create("ConditionalSampling", comm, params, grid);
    diag->set_input_field(x);
    diag->set_input_field(y);
    diag->set_input_field(z);
    REQUIRE_THROWS (diag->initialize());

    // ...unless one of them never declared its units
    FieldIdentifier w_fid("z",x.get_header().get_identifier().get_layout(),Units::invalid(),grid->name());
    Field w(w_fid,true);
    w.get_header().get_tracking().update_time_stamp(t0());
    auto diag_ok = diag_factory.create("ConditionalSampling", comm, params, grid);
    diag_ok->set_input_field(x);
    diag_ok->set_input_field(y);
    diag_ok->set_input_field(w);
    REQUIRE_NOTHROW (diag_ok->initialize());
  }
}

} // namespace scream
