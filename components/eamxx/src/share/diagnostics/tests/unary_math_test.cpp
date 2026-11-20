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

    // pow is not a unary op
    auto p_pow = make_params(gn,"q","pow");
    p_pow.set<std::string>("exponent","2");
    REQUIRE_THROWS (factory.create("UnaryMath",comm,p_pow,grid));

    for (std::string op : {"log","exp","sqrt","abs"}) {
      const bool has_invalid = op=="log" or op=="sqrt";
      const bool has_units   = op=="log" or op=="exp";

      auto p_inv = make_params(gn,"q",op);
      p_inv.set<std::string>("invalid","zero");
      REQUIRE_THROWS (factory.create("UnaryMath",comm,p_inv,grid));
      p_inv.set<std::string>("invalid","mask");
      if (has_invalid) {
        REQUIRE_NOTHROW (factory.create("UnaryMath",comm,p_inv,grid));
      } else {
        REQUIRE_THROWS (factory.create("UnaryMath",comm,p_inv,grid));
      }

      // NOTE: 'units' is read when the diag is initialized
      auto p_units = make_params(gn,"q",op);
      p_units.set<std::string>("units","none");
      REQUIRE_THROWS (make_diag(p_units,q,grid));
      auto p_ignore = make_params(gn,"q",op);
      p_ignore.set<std::string>("units","ignore");
      if (has_units) {
        REQUIRE_NOTHROW (make_diag(p_ignore,q,grid));
      } else {
        REQUIRE_THROWS (make_diag(p_ignore,q,grid));
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

    auto units_of = [&](const std::string& op) {
      auto p = make_params(gn,"pm",op);
      return make_diag(p,pm,grid)->get().get_header().get_identifier().get_units();
    };
    REQUIRE (units_of("sqrt")==sqrt(Pa));
    REQUIRE (units_of("abs")==Pa);
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
      std::function<Real(Real)> ref;
      std::function<bool(Real)> valid;
    };
    const std::vector<Case> cases = {
      {"log", [](Real x){ return std::log(x);    }, [](Real x){ return x>0; }},
      {"exp", [](Real x){ return std::exp(x);    }, [](Real)  { return true; }},
      {"sqrt", [](Real x){ return std::sqrt(x);   }, [](Real x){ return x>=0; }},
      {"abs", [](Real x){ return std::abs(x);    }, [](Real)  { return true; }},
    };

    constexpr auto fv = constants::fill_value<Real>;
    for (auto x : {q,v}) {
      const auto n = x.get_header().get_identifier().get_layout().size();
      const auto x_h = host_data(x);
      for (const auto& c : cases) {
        const bool has_invalid = c.op=="log" or c.op=="sqrt";
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

  SECTION ("wind_speed") {
    // sqrt(U*U+V*V), with U and V the components of horiz_winds, is what the wind_speed
    // diag computes: an early check of sqrt against a diag that does it by hand
    const auto uv_layout = grid->get_3d_vector_layout(LEV,2,"cmp");
    Field winds (FieldIdentifier("horiz_winds",uv_layout,m/s,gn), true);
    winds.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(winds,seed++,-10,10);
    const auto U = winds.get_component(0);
    const auto V = winds.get_component(1);

    const auto speed2 = (m/s)*(m/s);
    Field s2 (FieldIdentifier("s2",scalar3d,speed2,gn), true);
    Field v2 (FieldIdentifier("v2",scalar3d,speed2,gn), true);
    s2.get_header().get_tracking().update_time_stamp(t0);
    s2.deep_copy(U);
    s2.scale(U);
    v2.deep_copy(V);
    v2.scale(V);
    s2.update(v2,1,1);

    auto sqrt_diag = make_diag(make_params(gn,"s2","sqrt"),s2,grid);
    sqrt_diag->compute(t0);

    ekat::ParameterList ws_params;
    ws_params.set<std::string>("grid_name",gn);
    auto ws_diag = factory.create("wind_speed",comm,ws_params,grid);
    ws_diag->set_input_field(winds);
    ws_diag->initialize();
    ws_diag->compute(t0);

    auto sq = sqrt_diag->get();
    auto ws = ws_diag->get();
    REQUIRE (sq.get_header().get_identifier().get_units()==m/s);
    REQUIRE (sq.get_header().get_identifier().get_units()==ws.get_header().get_identifier().get_units());
    sq.sync_to_host();
    ws.sync_to_host();
    const auto sq_h = sq.get_view<const Real**,Host>();
    const auto ws_h = ws.get_view<const Real**,Host>();
    for (int icol=0; icol<grid->get_num_local_dofs(); ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        REQUIRE (sq_h(icol,ilev)==Approx(ws_h(icol,ilev)));
      }
    }
  }

  SECTION ("masked_input") {
    // Masked entries of X are not computed, so they may hold anything, e.g.
    // the fill value, whose exp overflows.
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
      auto p = make_params(gn,"x","log");
      p.set<std::string>("invalid",invalid);
      if (invalid=="nan") {
        p = make_params(gn,"x","exp");
      }
      auto diag = make_diag(p,x,grid);
      diag->compute(t0);
      auto d = diag->get();
      REQUIRE (d.has_valid_mask());
      d.sync_to_host();
      d.get_valid_mask().sync_to_host();
      const auto d_h = host_data(d);
      const auto m_h = host_mask(d.get_valid_mask());
      for (int i=0; i<n; ++i) {
        const bool valid = xm_h[i]!=0 and (invalid=="nan" or x_h[i]>0);
        REQUIRE (m_h[i]==(valid ? 1 : 0));
        if (valid) {
          const auto ref = invalid=="nan" ? std::exp(x_h[i]) : std::log(x_h[i]);
          REQUIRE (d_h[i]==Approx(ref));
        } else {
          REQUIRE (d_h[i]==fv);
        }
      }
    }

    // sqrt, with invalid='mask'
    auto p = make_params(gn,"x","sqrt");
    p.set<std::string>("invalid","mask");
    auto diag = make_diag(p,x,grid);
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
