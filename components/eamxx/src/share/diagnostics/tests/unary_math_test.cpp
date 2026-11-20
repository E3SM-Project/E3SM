#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_fpe.hpp>

#include <cfenv>
#include <cmath>
#include <functional>
#include <string>

namespace scream {

namespace {

using namespace ShortFieldTagsNames;

ekat::ParameterList make_params (const std::string& grid_name,
                                 const std::string& field_name,
                                 const std::string& op)
{
  ekat::ParameterList params;
  params.set<std::string>("grid_name",grid_name);
  params.set<std::string>("field_name",field_name);
  params.set<std::string>("op",op);
  return params;
}

std::shared_ptr<AbstractDiagnostic>
make_diag (const ekat::ParameterList& params, const Field& f,
           const std::shared_ptr<const AbstractGrid>& grid)
{
  auto diag = DiagnosticFactory::instance().create("UnaryMath",grid->get_comm(),params,grid);
  diag->set_input_field(f);
  diag->initialize();
  return diag;
}

// The fields in these tests are not padded, so their data is contiguous
// and has layout.size() entries
const Real* host_data (const Field& f) {
  return f.get_internal_view_data<const Real,Host>();
}
const int* host_mask (const Field& f) {
  return f.get_internal_view_data<const int,Host>();
}

} // anonymous namespace

TEST_CASE("unary_math") {
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
    REQUIRE_THROWS (factory.create("UnaryMath",comm,params,grid)); // No field_name
    params.set<std::string>("field_name","q");
    REQUIRE_THROWS (factory.create("UnaryMath",comm,params,grid)); // No op

    // NOTE: ekat::ParameterList::operator= merges into the existing entries
    //       (it does not drop them), so each case below uses a fresh list.
    auto p_sin = make_params(gn,"q","sin");
    REQUIRE_THROWS (factory.create("UnaryMath",comm,p_sin,grid)); // Unknown op

    auto p_pow = make_params(gn,"q","pow");
    REQUIRE_THROWS (factory.create("UnaryMath",comm,p_pow,grid)); // pow needs an exponent
    for (std::string bad : {"", "two", "2x", "nan", "inf"}) {
      p_pow.set<std::string>("exponent",bad);
      REQUIRE_THROWS (factory.create("UnaryMath",comm,p_pow,grid));
    }
    p_pow.set<std::string>("exponent","2");
    REQUIRE_NOTHROW (factory.create("UnaryMath",comm,p_pow,grid));

    auto p_log = make_params(gn,"q","log");
    p_log.set<std::string>("exponent","2");
    REQUIRE_THROWS (factory.create("UnaryMath",comm,p_log,grid)); // exponent is only for pow

    for (std::string op : {"log","exp","sqrt","abs","pow"}) {
      const bool has_invalid = op=="log" or op=="sqrt" or op=="pow";
      const bool has_units   = op=="log" or op=="exp"  or op=="pow";

      auto make = [&]() {
        auto pp = make_params(gn,"q",op);
        if (op=="pow") {
          pp.set<std::string>("exponent","2");
        }
        return pp;
      };
      auto p_inv = make();
      p_inv.set<std::string>("invalid","zero");
      REQUIRE_THROWS (factory.create("UnaryMath",comm,p_inv,grid));
      p_inv.set<std::string>("invalid","mask");
      if (has_invalid) {
        REQUIRE_NOTHROW (factory.create("UnaryMath",comm,p_inv,grid));
      } else {
        REQUIRE_THROWS (factory.create("UnaryMath",comm,p_inv,grid));
      }

      auto p_units = make();
      p_units.set<std::string>("units","none");
      REQUIRE_THROWS (factory.create("UnaryMath",comm,p_units,grid));
      p_units.set<std::string>("units","ignore");
      if (has_units) {
        REQUIRE_NOTHROW (factory.create("UnaryMath",comm,p_units,grid));
      } else {
        REQUIRE_THROWS (factory.create("UnaryMath",comm,p_units,grid));
      }
    }
  }

  SECTION ("units") {
    // log/exp need a dimensionless input, unless units='ignore'
    for (std::string op : {"log","exp"}) {
      auto p_q = make_params(gn,"q",op);
      REQUIRE (make_diag(p_q,q,grid)->get().get_header().get_identifier().get_units()==none);

      auto p_pm = make_params(gn,"pm",op);
      REQUIRE_THROWS (make_diag(p_pm,pm,grid));
      p_pm.set<std::string>("units","ignore");
      REQUIRE (make_diag(p_pm,pm,grid)->get().get_header().get_identifier().get_units()==none);
    }

    auto units_of = [&](const std::string& op, const std::string& expo,
                        const std::string& units) {
      auto p = make_params(gn,"pm",op);
      if (expo!="") {
        p.set<std::string>("exponent",expo);
      }
      if (units!="") {
        p.set<std::string>("units",units);
      }
      return make_diag(p,pm,grid)->get().get_header().get_identifier().get_units();
    };
    REQUIRE (units_of("sqrt","",  "")==sqrt(Pa));
    REQUIRE (units_of("abs", "",  "")==Pa);
    REQUIRE (units_of("pow", "2", "")==Pa*Pa);
    REQUIRE (units_of("pow", "-1","")==1/Pa);
    REQUIRE (units_of("pow", "0.5","")==sqrt(Pa));
    REQUIRE (units_of("pow", "1.5","")==pow(Pa,ekat::RationalConstant(3,2)));
    REQUIRE_THROWS (units_of("pow","0.2857",""));        // Not n/d with d<=12
    REQUIRE (units_of("pow","0.2857","ignore")==none);
    REQUIRE (units_of("pow","2",     "ignore")==none);   // 'ignore' always gives none

    // A dimensionless, unscaled input stays so, for any exponent
    auto p = make_params(gn,"q","pow");
    p.set<std::string>("exponent","0.2857");
    REQUIRE (make_diag(p,q,grid)->get().get_header().get_identifier().get_units()==none);
  }

  SECTION ("values") {
    // A rank-2 and a rank-3 layout, so that field_pfor is exercised on different ranks
    const auto vector3d = grid->get_3d_vector_layout(LEV,ncmps,"cmp");
    Field v(FieldIdentifier("v",vector3d,kg/kg,gn), true);
    v.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(v,seed++,-2,2);

    // Some exact zeros (and, below, integers), to hit the edges of the domains
    for (auto f : {q,v}) {
      auto data = f.get_internal_view_data<Real,Host>();
      const auto n = f.get_header().get_identifier().get_layout().size();
      for (int i=0; i<n; i+=7) {
        data[i] = 0;
      }
      f.sync_to_dev();
    }

    // With FE_INVALID trapped (fpe builds), invalid entries can only be tested with invalid='mask'
    const bool trap_invalid = (ekat::get_enabled_fpes() & FE_INVALID)!=0;
    const bool trap_divbyzero = (ekat::get_enabled_fpes() & FE_DIVBYZERO)!=0;

    struct Case {
      std::string op;
      std::string expo;
      std::function<Real(Real)> ref;
      std::function<bool(Real)> valid;
    };
    const std::vector<Case> cases = {
      {"log", "",    [](Real x){ return std::log(x);    }, [](Real x){ return x>0; }},
      {"exp", "",    [](Real x){ return std::exp(x);    }, [](Real)  { return true; }},
      {"sqrt","",    [](Real x){ return std::sqrt(x);   }, [](Real x){ return x>=0; }},
      {"abs", "",    [](Real x){ return std::abs(x);    }, [](Real)  { return true; }},
      {"pow", "3",   [](Real x){ return std::pow(x,Real(3));   }, [](Real)  { return true; }},
      {"pow", "-2",  [](Real x){ return std::pow(x,Real(-2));  }, [](Real x){ return x!=0; }},
      {"pow", "0.5", [](Real x){ return std::pow(x,Real(0.5)); }, [](Real x){ return x>=0; }},
      {"pow", "-1.5",[](Real x){ return std::pow(x,Real(-1.5));}, [](Real x){ return x>0; }},
    };

    constexpr auto fv = constants::fill_value<Real>;
    for (auto x : {q,v}) {
      const auto n = x.get_header().get_identifier().get_layout().size();
      const auto x_h = host_data(x);
      for (const auto& c : cases) {
        const bool has_invalid = c.op=="log" or c.op=="sqrt" or c.op=="pow";
        for (std::string invalid : {"nan","mask"}) {
          if (invalid=="mask" and not has_invalid) {
            continue;
          }
          // log(0) and 0^(-p) raise FE_DIVBYZERO, the others raise FE_INVALID
          if (invalid=="nan" and (trap_invalid or trap_divbyzero)) {
            bool any_invalid = false;
            for (int i=0; i<n; ++i) {
              any_invalid |= not c.valid(x_h[i]);
            }
            if (any_invalid) {
              continue;
            }
          }
          auto p = make_params(gn,x.name(),c.op);
          if (c.expo!="") {
            p.set<std::string>("exponent",c.expo);
          }
          p.set<std::string>("invalid",invalid);

          auto diag = make_diag(p,x,grid);
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
            const bool valid = c.valid(x_h[i]);
            if (masked) {
              REQUIRE (m_h[i]==(valid ? 1 : 0));
            }
            if (valid) {
              REQUIRE (d_h[i]==Approx(c.ref(x_h[i])));
            } else if (masked) {
              REQUIRE (d_h[i]==fv);
            } else {
              // Whatever the math library gives, NaN or +/-inf
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

    for (std::string invalid : {"nan","mask"}) {
      auto p = make_params(gn,"x","pow");
      p.set<std::string>("exponent","3");
      p.set<std::string>("invalid",invalid);
      auto diag = make_diag(p,x,grid);
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

    // With invalid='mask', the output mask is X's mask AND the domain
    auto p = make_params(gn,"x","log");
    p.set<std::string>("invalid","mask");
    auto diag = make_diag(p,x,grid);
    diag->compute(t0);
    auto d = diag->get();
    d.sync_to_host();
    d.get_valid_mask().sync_to_host();
    const auto d_h = host_data(d);
    const auto m_h = host_mask(d.get_valid_mask());
    for (int i=0; i<n; ++i) {
      const bool valid = xm_h[i]!=0 and x_h[i]>0;
      REQUIRE (m_h[i]==(valid ? 1 : 0));
      if (valid) {
        REQUIRE (d_h[i]==Approx(std::log(x_h[i])));
      } else {
        REQUIRE (d_h[i]==fv);
      }
    }
  }
}

} // namespace scream
