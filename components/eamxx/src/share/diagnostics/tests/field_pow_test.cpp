#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_fpe.hpp>

#include <cfenv>
#include <cmath>
#include <string>

namespace scream {

namespace {

using namespace ShortFieldTagsNames;

ekat::ParameterList make_params (const std::string& grid_name,
                                 const std::string& field_name,
                                 const std::string& exponent)
{
  ekat::ParameterList params;
  params.set<std::string>("grid_name",grid_name);
  params.set<std::string>("field_name",field_name);
  params.set<std::string>("exponent",exponent);
  return params;
}

std::shared_ptr<AbstractDiagnostic>
make_diag (const ekat::ParameterList& params, const Field& f,
           const std::shared_ptr<const AbstractGrid>& grid)
{
  auto diag = DiagnosticFactory::instance().create("FieldPow",grid->get_comm(),params,grid);
  diag->set_input_field(f);
  diag->initialize();
  return diag;
}

// The fields in these tests are not padded, so their data is contiguous
const Real* host_data (const Field& f) {
  return f.get_internal_view_data<const Real,Host>();
}
const int* host_mask (const Field& f) {
  return f.get_internal_view_data<const int,Host>();
}

// Whether x^p is defined (the exponent is never 0). Whether p is an integer is passed in,
// since it is a property of the exponent as written, not of its (maybe rounded) value
bool in_domain (const Real x, const double p, const bool p_is_int) {
  return x>0 or (x==0 and p>0) or (x<0 and p_is_int);
}

} // anonymous namespace

TEST_CASE("field_pow") {
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  const int nlevs = 17;
  const int ncols = 13;
  const int ncmps = 3;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);
  const auto gn = grid->name();

  register_diagnostics();
  auto& factory = DiagnosticFactory::instance();

  int seed = get_random_test_seed(&comm);

  // A dimensionless field and a field with units, so that the units checks can be exercised
  const auto scalar3d = grid->get_3d_scalar_layout(LEV);
  Field q (FieldIdentifier("q", scalar3d, kg/kg, gn), true);
  Field pm(FieldIdentifier("pm",scalar3d, Pa,    gn), true);
  q.get_header().get_tracking().update_time_stamp(t0);
  pm.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(q, seed++,-2,2);
  randomize_uniform(pm,seed++, 1,2);

  SECTION ("params") {
    ekat::ParameterList params;
    params.set<std::string>("grid_name",gn);
    REQUIRE_THROWS (factory.create("FieldPow",comm,params,grid)); // No field_name
    params.set<std::string>("field_name","q");
    REQUIRE_THROWS (factory.create("FieldPow",comm,params,grid)); // No exponent

    // Numbers and fractions are fine...
    for (std::string good : {"2", "-1", "0.5", "-1.5", "1e-3", "1e20", "1/3", "-2/3", "7/2", "4/2"}) {
      REQUIRE_NOTHROW (factory.create("FieldPow",comm,make_params(gn,"q",good),grid));
    }
    // ...and anything else is not. Note that 0 is not a valid exponent: 0^0 is not defined.
    for (std::string bad : {"", "two", "2x", "nan", "inf", "0", "0.0", "1e-400", "0/3", "1/0", "1/", "/2", "1/2/3", "0.5/2", "1.5/2",
                            "1/-2",   // the sign goes on the numerator
                            // Out of range: strtoll would saturate these, which would be a different exponent
                            "18446744073709551614/9223372036854775807", "9223372036854775808/1", "1/99999999999999999999"}) {
      REQUIRE_THROWS (factory.create("FieldPow",comm,make_params(gn,"q",bad),grid));
    }

    auto p_inv = make_params(gn,"q","2");
    p_inv.set<std::string>("invalid","zero");
    REQUIRE_THROWS (factory.create("FieldPow",comm,p_inv,grid));
    for (std::string ok : {"nan","mask"}) {
      auto p = make_params(gn,"q","2");
      p.set<std::string>("invalid",ok);
      REQUIRE_NOTHROW (factory.create("FieldPow",comm,p,grid));
    }

    // NOTE: 'units' is read when the diag is initialized
    auto p_units = make_params(gn,"q","2");
    p_units.set<std::string>("units","none");
    REQUIRE_THROWS (make_diag(p_units,q,grid));
    for (std::string ok : {"check","ignore"}) {
      auto p = make_params(gn,"q","2");
      p.set<std::string>("units",ok);
      REQUIRE_NOTHROW (make_diag(p,q,grid));
    }
  }

  SECTION ("units") {
    auto units_of = [&](const std::string& expo, const std::string& units) {
      auto p = make_params(gn,"pm",expo);
      if (units!="") {
        p.set<std::string>("units",units);
      }
      return make_diag(p,pm,grid)->get().get_header().get_identifier().get_units();
    };
    REQUIRE (units_of("2",  "")==Pa*Pa);
    REQUIRE (units_of("3",  "")==Pa*Pa*Pa);
    REQUIRE (units_of("-1", "")==1/Pa);
    REQUIRE (units_of("0.5","")==sqrt(Pa));
    REQUIRE (units_of("-1.5","")==pow(Pa,ekat::RationalConstant(-3,2)));
    // A fraction is exact, and so are its equivalents
    REQUIRE (units_of("1/3","")==pow(Pa,ekat::RationalConstant(1,3)));
    REQUIRE (units_of("2/6","")==pow(Pa,ekat::RationalConstant(1,3)));
    REQUIRE (units_of("-5/17","")==pow(Pa,ekat::RationalConstant(-5,17)));  // d>12 is fine as a fraction
    REQUIRE (units_of("0.2","")==pow(Pa,ekat::RationalConstant(1,5)));   // even if Real is a float
    // A decimal that is not a fraction with a small denominator has no units...
    REQUIRE_THROWS (units_of("0.2857",""));
    REQUIRE_THROWS (units_of("0.3333333",""));
    REQUIRE_THROWS (units_of("1e20",""));                  // too large for a numerator
    // ...unless we say so
    REQUIRE (units_of("0.2857","ignore")==none);
    REQUIRE (units_of("1e20",  "ignore")==none);
    REQUIRE (units_of("2",     "ignore")==none);   // 'ignore' always gives none

    // A dimensionless, unscaled input stays so, for any exponent
    for (std::string expo : {"2","0.2857","1/3"}) {
      auto p = make_params(gn,"q",expo);
      REQUIRE (make_diag(p,q,grid)->get().get_header().get_identifier().get_units()==none);
    }
  }

  SECTION ("values") {
    // A rank-2 and a rank-3 layout, so that field_pfor is exercised on different ranks
    const auto vector3d = grid->get_3d_vector_layout(LEV,ncmps,"cmp");
    Field v(FieldIdentifier("v",vector3d,kg/kg,gn), true);
    v.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(v,seed++,-2,2);

    // Some exact zeros, to hit the edges of the domains
    for (auto f : {q,v}) {
      auto data = f.get_internal_view_data<Real,Host>();
      const auto n = f.get_header().get_identifier().get_layout().size();
      for (int i=0; i<n; i+=7) {
        data[i] = 0;
      }
      f.sync_to_dev();
    }

    // With FPEs trapped (fpe builds), invalid entries can only be tested with invalid='mask'
    const auto fpes = ekat::get_enabled_fpes();
    const bool trap = (fpes & (FE_INVALID | FE_DIVBYZERO))!=0;

    constexpr auto fv = constants::fill_value<Real>;

    // The exponents with a fast path (2, 3, -1) and the general ones, integer or not. The last
    // one is not an integer, but rounds to 2 in a float: it must not take the square's path, nor be
    // seen as an integer. Then what a negative x^p gives with invalid='nan' depends on the precision.
    struct Expo {
      std::string s;
      double      p;
      bool        is_int;
      bool        nan_is_not_finite;
    };
    const std::vector<Expo> expos = {
      {"2",2,true,true}, {"3",3,true,true}, {"-1",-1,true,true}, {"-2",-2,true,true},
      {"4",4,true,true}, {"0.5",0.5,false,true}, {"1.5",1.5,false,true}, {"-1.5",-1.5,false,true},
      {"1/3",1.0/3,false,true}, {"-2/3",-2.0/3,false,true}, {"1",1,true,true},
      {"4/2",2,true,true}, {"33554433/16777216",2.0+1.0/16777216,false,false},
    };
    for (const auto& e : expos) {
      const Real p = static_cast<Real>(e.p);
      for (auto x : {q,v}) {
        const auto n = x.get_header().get_identifier().get_layout().size();
        const auto x_h = host_data(x);

        for (std::string invalid : {"nan","mask"}) {
          if (invalid=="nan" and trap) {
            bool any_invalid = false;
            for (int i=0; i<n; ++i) {
              any_invalid |= not in_domain(x_h[i],e.p,e.is_int);
            }
            if (any_invalid) {
              continue;
            }
          }
          auto params = make_params(gn,x.name(),e.s);
          params.set<std::string>("invalid",invalid);

          auto diag = make_diag(params,x,grid);
          diag->compute(t0);
          auto d = diag->get();
          d.sync_to_host();
          const auto d_h = host_data(d);

          const bool masked = invalid=="mask";
          REQUIRE (d.has_valid_mask()==masked);
          const int* m_h = nullptr;
          if (masked) {
            d.get_valid_mask().sync_to_host();
            m_h = host_mask(d.get_valid_mask());
          }

          for (int i=0; i<n; ++i) {
            const bool valid = in_domain(x_h[i],e.p,e.is_int);
            if (masked) {
              REQUIRE (m_h[i]==(valid ? 1 : 0));
            }
            if (valid) {
              REQUIRE (d_h[i]==Approx(std::pow(x_h[i],p)));
            } else if (masked) {
              REQUIRE (d_h[i]==fv);
            } else if (e.nan_is_not_finite) {
              // Whatever the math library gives, NaN or inf
              REQUIRE (not std::isfinite(d_h[i]));
            }
          }
        }
      }
    }
  }

  SECTION ("masked_input") {
    // Masked entries of X are not computed, so they may hold anything, e.g.
    // the fill value, whose cube overflows.
    Field x = q.clone("x",CloneFlags::CopyData);
    x.get_header().get_tracking().update_time_stamp(t0);
    x.create_valid_mask();
    randomize_uniform(x.get_valid_mask(),seed++);
    x.get_valid_mask().sync_to_host();

    constexpr auto fv = constants::fill_value<Real>;
    const auto n = x.get_header().get_identifier().get_layout().size();
    const auto xm_h = host_mask(x.get_valid_mask());
    {
      auto data = x.get_internal_view_data<Real,Host>();
      for (int i=0; i<n; ++i) {
        if (xm_h[i]==0) {
          data[i] = fv;
        }
      }
      x.sync_to_dev();
    }
    const auto x_h = host_data(x);

    // With invalid='nan', the output mask is X's; with invalid='mask', X's AND the domain
    for (std::string invalid : {"nan","mask"}) {
      auto params = make_params(gn,"x","3");
      params.set<std::string>("invalid",invalid);
      auto diag = make_diag(params,x,grid);
      diag->compute(t0);
      auto d = diag->get();
      REQUIRE (d.has_valid_mask());
      d.sync_to_host();
      d.get_valid_mask().sync_to_host();
      const auto d_h = host_data(d);
      const auto m_h = host_mask(d.get_valid_mask());
      for (int i=0; i<n; ++i) {
        REQUIRE (m_h[i]==xm_h[i]);
        if (xm_h[i]) {
          REQUIRE (d_h[i]==Approx(std::pow(x_h[i],Real(3))));
        } else {
          REQUIRE (d_h[i]==fv);
        }
      }
    }

    auto params = make_params(gn,"x","0.5");
    params.set<std::string>("invalid","mask");
    auto diag = make_diag(params,x,grid);
    diag->compute(t0);
    auto d = diag->get();
    d.sync_to_host();
    d.get_valid_mask().sync_to_host();
    const auto d_h = host_data(d);
    const auto m_h = host_mask(d.get_valid_mask());
    for (int i=0; i<n; ++i) {
      const bool valid = xm_h[i]!=0 and x_h[i]>=0;
      REQUIRE (m_h[i]==(valid ? 1 : 0));
      if (valid) {
        REQUIRE (d_h[i]==Approx(std::sqrt(x_h[i])));
      } else {
        REQUIRE (d_h[i]==fv);
      }
    }
  }
}

} // namespace scream
