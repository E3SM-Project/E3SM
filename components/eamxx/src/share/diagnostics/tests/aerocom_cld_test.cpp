#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace scream {

namespace {

using HostView1 = Kokkos::View<Real*,Kokkos::HostSpace>;
using HostView2 = Kokkos::View<Real**,Kokkos::HostSpace>;
using fields_t  = std::map<std::string,Field>;

const std::vector<std::string> phases = {"none","liq","ice"};
const util::TimeStamp t0({2024,1,1},{0,0,0});

HostView2 host_view (const Field& f) {
  f.sync_to_host();
  return f.get_view<Real**,Host>();
}

// A host, serial reference of the weights. Top looks at each layer from the
// one above it, starting at the second; Bot is its mirror.
HostView2 weights_reference (const std::string& position, const std::string& phase,
                             const fields_t& fields)
{
  const auto cld = host_view(fields.at("cldfrac_tot"));
  const auto qc  = host_view(fields.at("qc"));
  const auto qi  = host_view(fields.at("qi"));
  const int ncols = cld.extent(0);
  const int nlevs = cld.extent(1);
  const Real q_thresh   = AeroComWeights::q_threshold;
  const Real cld_thresh = AeroComWeights::cldfrac_tot_threshold;

  HostView2 w("w",ncols,nlevs);
  for (int icol=0; icol<ncols; ++icol) {
    Real clr = 1;
    auto visit = [&](const int k, const int nb) {
      const Real q = qc(icol,k) + qi(icol,k);
      if (q>q_thresh and cld(icol,k)>cld_thresh) {
        const Real clr_k = clr * (1 - std::max(cld(icol,nb),cld(icol,k))) /
                                 (1 - std::min(cld(icol,nb),1-cld_thresh));
        const Real phi = qc(icol,k)/q;
        const Real Phi = phase=="liq" ? phi : (phase=="ice" ? 1-phi : 1);
        w(icol,k) = (clr - clr_k)*Phi;
        clr = clr_k;
      }
    };
    if (position=="top") {
      for (int k=1; k<nlevs; ++k) visit(k,k-1);
    } else {
      for (int k=nlevs-2; k>=0; --k) visit(k,k+1);
    }
  }
  return w;
}

// The reference of X.aerocom(position,phase=phase): the sum of X*W over the
// levels with nonzero weight (and valid X, if X is masked). X=1 for 'mask'.
HostView1 sample_reference (const std::string& X, const std::string& position,
                            const std::string& phase, const fields_t& fields)
{
  const auto w = weights_reference(position,phase,fields);
  const int ncols = w.extent(0);
  const int nlevs = w.extent(1);
  HostView2 x;
  Kokkos::View<int**,Kokkos::HostSpace> m;
  if (X!="mask") {
    x = host_view(fields.at(X));
    if (fields.at(X).has_valid_mask()) {
      fields.at(X).get_valid_mask().sync_to_host();
      m = fields.at(X).get_valid_mask().get_view<int**,Host>();
    }
  }
  HostView1 out("out",ncols);
  for (int icol=0; icol<ncols; ++icol) {
    for (int k=0; k<nlevs; ++k) {
      if (w(icol,k)!=0 and (m.size()==0 or m(icol,k)!=0)) {
        out(icol) += (X=="mask" ? 1 : x(icol,k))*w(icol,k);
      }
    }
  }
  return out;
}

std::shared_ptr<AbstractDiagnostic>
create_weights (const std::string& position, const std::string& phase,
                const fields_t& fields, const std::shared_ptr<const AbstractGrid>& grid)
{
  ekat::ParameterList params;
  params.set<std::string>("position",position);
  params.set<std::string>("phase",phase);
  auto w = DiagnosticFactory::instance().create("AeroComWeights",grid->get_comm(),params,grid);
  for (const std::string n : {"cldfrac_tot","qc","qi"}) {
    w->set_input_field(fields.at(n));
  }
  w->initialize();
  return w;
}

std::shared_ptr<AbstractDiagnostic>
create_sampler (const std::string& X, const std::shared_ptr<AbstractDiagnostic>& w,
                const fields_t& fields, const std::shared_ptr<const AbstractGrid>& grid)
{
  ekat::ParameterList params;
  params.set<std::string>("field_name",X);
  params.set<std::string>("position",w->get_params().get<std::string>("position"));
  params.set<std::string>("phase",w->get_params().get<std::string>("phase"));
  auto d = DiagnosticFactory::instance().create("AeroComCld",grid->get_comm(),params,grid);
  if (X!="mask") {
    d->set_input_field(fields.at(X));
  }
  d->set_input_field(w->get());
  d->initialize();
  return d;
}

// X.aerocom(position,phase=phase), from new diags
HostView1 sample_diag (const std::string& X, const std::string& position, const std::string& phase,
                  const fields_t& fields, const std::shared_ptr<const AbstractGrid>& grid)
{
  auto w = create_weights(position,phase,fields,grid);
  auto d = create_sampler(X,w,fields,grid);
  w->compute(t0);
  d->compute(t0);
  auto f = d->get();
  f.sync_to_host();
  // A copy: the field's views do not keep its memory alive
  HostView1 x("x",f.get_header().get_identifier().get_layout().dim(0));
  Kokkos::deep_copy(x,f.get_view<const Real*,Host>());
  return x;
}

// Every value finite and equal to the reference, up to rounding (the diag may
// run on device, where the compiler may contract to fma). Prints the first few
// mismatches, so that a failure says what kind of wrong it is.
template<typename V, typename R>
void check (const V& v, const R& ref, const std::string& label)
{
  const Real tol = std::sqrt(std::numeric_limits<Real>::epsilon());
  REQUIRE (v.size()==ref.size());
  int nbad = 0, nonfinite = 0;
  for (size_t i=0; i<v.size(); ++i) {
    const Real d = v.data()[i];
    const Real r = ref.data()[i];
    const bool ok = std::isfinite(d) and d==Approx(r).epsilon(tol).margin(tol);
    if (not ok) {
      if (nbad<5) {
        std::cout << " AeroCom mismatch (" << label << ") - entry: " << i
                  << " - diag: " << d << " - reference: " << r << "\n";
      }
      ++nbad;
      nonfinite += not std::isfinite(d);
    }
  }
  INFO (label + ": " + std::to_string(nbad) + " wrong values, " +
        std::to_string(nonfinite) + " of them not finite");
  REQUIRE (nbad==0);
}

// The weights, and the samples of two fields and of 'mask', at both positions
// and in every phase, against the reference
void check_all (const fields_t& fields, const std::shared_ptr<const AbstractGrid>& grid,
                const std::string& label)
{
  for (const std::string pos : {"top","bot"}) {
    for (const auto& ph : phases) {
      const auto lbl = label + ", " + pos + ", phase=" + ph;
      auto w = create_weights(pos,ph,fields,grid);
      w->compute(t0);
      check(host_view(w->get()),weights_reference(pos,ph,fields),lbl + ", weights");
      for (const std::string X : {"T_mid","nc","mask"}) {
        check(sample_diag(X,pos,ph,fields,grid),sample_reference(X,pos,ph,fields),lbl + ", " + X);
      }
    }
  }
}

fields_t create_inputs (const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto micron = m / 1000000;
  const auto layout = grid->get_3d_scalar_layout(LEV);
  const std::vector<std::pair<std::string,ekat::units::Units>> inputs = {
    {"T_mid", K}, {"p_mid", Pa}, {"qc", kg/kg}, {"qi", kg/kg},
    {"eff_radius_qc", micron}, {"eff_radius_qi", micron},
    {"cldfrac_tot", none}, {"nc", 1/kg}, {"ni", 1/kg}
  };
  fields_t fields;
  for (const auto& [n,u] : inputs) {
    Field f(FieldIdentifier(n,layout,u,grid->name()));
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    fields.emplace(n,f);
  }
  return fields;
}

// Set inputs that look like an atmosphere: pressure growing downward, and clear,
// liquid-only, ice-only and fully cloudy layers mixed in every column.
void set_realistic_inputs (const fields_t& fields, int& seed)
{
  auto set = [&](const std::string& n, const Real lb, const Real ub) {
    randomize_uniform(fields.at(n),seed++,lb,ub);
    fields.at(n).sync_to_host();
  };
  set("T_mid",200,300);
  set("p_mid",1000,2000);
  set("qc",-5e-4,1e-3);
  set("qi",-5e-4,1e-3);
  set("eff_radius_qc",5,20);
  set("eff_radius_qi",10,50);
  set("cldfrac_tot",-0.3,1.2);
  set("nc",1e7,1e9);
  set("ni",1e3,1e6);

  auto pmid = fields.at("p_mid").get_view<Real**,Host>();
  auto qc   = fields.at("qc").get_view<Real**,Host>();
  auto qi   = fields.at("qi").get_view<Real**,Host>();
  auto cld  = fields.at("cldfrac_tot").get_view<Real**,Host>();
  for (int icol=0; icol<static_cast<int>(pmid.extent(0)); ++icol) {
    Real p = 100;
    for (int k=0; k<static_cast<int>(pmid.extent(1)); ++k) {
      p += pmid(icol,k);
      pmid(icol,k) = p;
      qc(icol,k)  = std::max(qc(icol,k),Real(0));
      qi(icol,k)  = std::max(qi(icol,k),Real(0));
      cld(icol,k) = std::min(std::max(cld(icol,k),Real(0)),Real(1));
    }
  }
  for (const std::string n : {"p_mid","qc","qi","cldfrac_tot"}) {
    fields.at(n).sync_to_dev();
  }
}

void hand_cases (const int nlevs)
{
  ekat::Comm comm(MPI_COMM_WORLD);
  auto grid = create_point_grid("physics",comm.size(),nlevs,comm);
  auto fields = create_inputs(grid);
  const std::string nl = "nlevs=" + std::to_string(nlevs);

  auto& T  = fields.at("T_mid");
  auto& p  = fields.at("p_mid");
  auto& qc = fields.at("qc");
  auto& qi = fields.at("qi");
  auto& cd = fields.at("cldfrac_tot");
  auto& nc = fields.at("nc");
  auto& ni = fields.at("ni");
  auto cd_v = cd.get_view<Real**,Host>();
  auto qc_v = qc.get_view<Real**,Host>();
  auto qi_v = qi.get_view<Real**,Host>();
  auto T_v  = T.get_view<Real**,Host>();

  // The value at the (only) column
  auto at = [&](const std::string& X, const std::string& pos, const std::string& ph) {
    return sample_diag(X,pos,ph,fields,grid)(0);
  };

  int seed = get_random_test_seed();
  for (auto& it : fields) {
    randomize_uniform(it.second,seed++,0,0.05);
  }

  // Case 1: no cloud, no weight
  cd.deep_copy(0.0);
  for (const std::string pos : {"top","bot"}) {
    for (const auto& ph : phases) {
      REQUIRE (at("T_mid",pos,ph)==0);
      REQUIRE (at("mask",pos,ph)==0);
    }
  }
  check_all(fields,grid,"case 1, " + nl);

  // Case 2: overcast, with as much liquid as ice: the first layer seen, from
  // either side, is the cloud
  cd.deep_copy(1.0);
  T.deep_copy(300.0);
  p.deep_copy(100.0);
  qc.deep_copy(1.0);
  qi.deep_copy(1.0);
  nc.deep_copy(5.0);
  ni.deep_copy(1.0);
  for (const std::string pos : {"top","bot"}) {
    REQUIRE (at("T_mid",pos,"none")==300.0);
    REQUIRE (at("p_mid",pos,"none")==100.0);
    REQUIRE (at("mask",pos,"none")==1.0);
    REQUIRE (at("mask",pos,"liq")==0.5);
    REQUIRE (at("mask",pos,"ice")==0.5);
    REQUIRE (at("nc",pos,"liq")==2.5);
    REQUIRE (at("ni",pos,"ice")==0.5);
  }
  check_all(fields,grid,"case 2, " + nl);

  // Case 3: contiguous cloudy layers overlap maximally
  Kokkos::deep_copy(cd_v,0);
  cd_v(0,1) = 0.5;
  cd_v(0,2) = 0.7;
  cd_v(0,3) = 0.3;
  cd_v(0,4) = 0.2;
  cd.sync_to_dev();
  REQUIRE (at("mask","top","none")==Approx(0.7));
  REQUIRE (at("mask","bot","none")==Approx(0.7));
  check_all(fields,grid,"case 3, " + nl);

  // Case 3xtra: a local minimum (0.1 between 0.2 and 0.4) enhances it slightly
  cd_v(0,5) = 0.1;
  cd_v(0,6) = 0.4;
  cd_v(0,7) = 0.2;
  cd.sync_to_dev();
  REQUIRE (at("mask","top","none")>0.7);
  REQUIRE (at("mask","bot","none")>0.7);
  check_all(fields,grid,"case 3xtra, " + nl);

  // Case 4: non-contiguous cloudy layers overlap randomly
  cd_v(0,6) = 0.0;
  cd_v(0,7) = 0.1;
  cd.sync_to_dev();
  REQUIRE (at("mask","top","none")>0.7);
  check_all(fields,grid,"case 4, " + nl);

  // Case 5: liquid only, then ice only
  cd_v(0,3) = 1.0;
  cd_v(0,7) = 1.0;
  cd_v(0,8) = 0.2;
  cd.sync_to_dev();
  qc.deep_copy(1.0);
  qi.deep_copy(0.0);
  for (const std::string pos : {"top","bot"}) {
    REQUIRE (at("mask",pos,"none")==Approx(1.0));
    REQUIRE (at("mask",pos,"liq")==Approx(1.0));
    REQUIRE (at("mask",pos,"ice")==0.0);
    REQUIRE (at("ni",pos,"ice")==0.0);
  }
  check_all(fields,grid,"case 5a, " + nl);
  qc.deep_copy(0.0);
  qi.deep_copy(1.0);
  for (const std::string pos : {"top","bot"}) {
    REQUIRE (at("mask",pos,"liq")==0.0);
    REQUIRE (at("mask",pos,"ice")==Approx(1.0));
    REQUIRE (at("nc",pos,"liq")==0.0);
  }
  check_all(fields,grid,"case 5b, " + nl);

  // Case 6: ice clouds above liquid ones hide them from the top, but not
  // from the bottom
  Kokkos::deep_copy(cd_v,0);
  Kokkos::deep_copy(qc_v,0);
  Kokkos::deep_copy(qi_v,0);
  cd_v(0,1) = 0.5;  qi_v(0,1) = 100;
  cd_v(0,2) = 0.7;  qi_v(0,2) = 200;
  cd_v(0,3) = 0.3;  qi_v(0,3) = 50;
  cd_v(0,5) = 0.2;  qc_v(0,5) = 20;
  cd_v(0,6) = 0.5;  qc_v(0,6) = 50;
  cd_v(0,7) = 0.1;  qc_v(0,7) = 10;
  cd.sync_to_dev();
  qc.sync_to_dev();
  qi.sync_to_dev();
  REQUIRE (at("mask","top","ice")==Approx(0.7));
  REQUIRE (at("mask","top","liq")<0.5);
  REQUIRE (at("mask","bot","liq")==Approx(0.5));
  REQUIRE (at("mask","bot","ice")<0.7);
  check_all(fields,grid,"case 6, " + nl);

  // Case 7: a single cloudy layer is both the cloud top and the cloud bottom
  Kokkos::deep_copy(cd_v,0);
  cd_v(0,nlevs/2) = 0.5;
  cd.sync_to_dev();
  qc.deep_copy(1.0);
  qi.deep_copy(3.0);
  T.deep_copy(300.0);
  REQUIRE (at("T_mid","top","none")==150.0);
  REQUIRE (at("mask","top","liq")==0.125);
  REQUIRE (at("mask","top","none")==0.5);
  for (const auto& ph : phases) {
    for (const std::string X : {"T_mid","p_mid","nc","ni","mask"}) {
      REQUIRE (at(X,"top",ph)==at(X,"bot",ph));
    }
  }
  check_all(fields,grid,"case 7, " + nl);

  // Case 8: two adjacent layers, the lower one cloudier. From the top, both
  // are seen; from the bottom, the lower one hides the upper one.
  const int ka = nlevs/2;
  const int kb = ka+1;
  Kokkos::deep_copy(cd_v,0);
  cd_v(0,ka) = 0.2;
  cd_v(0,kb) = 0.6;
  T_v(0,ka) = 250;
  T_v(0,kb) = 280;
  cd.sync_to_dev();
  T.sync_to_dev();
  auto w_top = create_weights("top","none",fields,grid);
  auto w_bot = create_weights("bot","none",fields,grid);
  w_top->compute(t0);
  w_bot->compute(t0);
  const auto wt = host_view(w_top->get());
  const auto wb = host_view(w_bot->get());
  REQUIRE (wt(0,ka)==Approx(0.2));
  REQUIRE (wt(0,kb)==Approx(0.4));
  REQUIRE (wb(0,kb)==Approx(0.6));
  REQUIRE (wb(0,ka)==Approx(0).margin(1e-12));
  REQUIRE (at("T_mid","top","none")==Approx(0.2*250 + 0.4*280));
  REQUIRE (at("T_mid","bot","none")==Approx(0.6*280));
  REQUIRE (at("mask","top","none")==Approx(0.6));
  REQUIRE (at("mask","bot","none")==Approx(0.6));
  check_all(fields,grid,"case 8, " + nl);
}

} // anonymous namespace

TEST_CASE("aerocom_cld") {
  register_diagnostics();

  // One level count below a GPU warp, and two above it
  for (const int nlevs : {9, 72, 128}) {
    hand_cases(nlevs);
  }
}

TEST_CASE("aerocom_cld_setup") {
  register_diagnostics();
  ekat::Comm comm(MPI_COMM_WORLD);
  auto grid = create_point_grid("physics",3*comm.size(),9,comm);
  auto& factory = DiagnosticFactory::instance();

  ekat::ParameterList params;
  params.set<std::string>("field_name","T_mid");
  params.set<std::string>("position","middle");
  REQUIRE_THROWS (factory.create("AeroComWeights",comm,params,grid));
  REQUIRE_THROWS (factory.create("AeroComCld",comm,params,grid));
  params.set<std::string>("position","top");
  params.set<std::string>("phase","water");
  REQUIRE_THROWS (factory.create("AeroComWeights",comm,params,grid));
  REQUIRE_THROWS (factory.create("AeroComCld",comm,params,grid));

  // The sampler asks for the weights by the name they are built with
  auto fields = create_inputs(grid);
  auto w = create_weights("top","liq",fields,grid);
  REQUIRE (w->get().name()=="cldfrac_tot.aerocom_weights('top',phase='liq')");
  auto d = create_sampler("T_mid",w,fields,grid);
  const auto& d_in = d->get_input_fields_names();
  REQUIRE (std::find(d_in.begin(),d_in.end(),w->get().name())!=d_in.end());
  REQUIRE (d->get().get_header().get_identifier().get_units()==ekat::units::K);

  // X must be a Real (COL,LEV) field
  Field ps(FieldIdentifier("ps",grid->get_2d_scalar_layout(),ekat::units::Pa,grid->name()));
  ps.allocate_view();
  fields.emplace("ps",ps);
  REQUIRE_THROWS (create_sampler("ps",w,fields,grid));
}

TEST_CASE("aerocom_cld_skips") {
  register_diagnostics();
  ekat::Comm comm(MPI_COMM_WORLD);
  auto grid = create_point_grid("physics",64*comm.size(),72,comm);
  auto fields = create_inputs(grid);
  int seed = get_random_test_seed();
  set_realistic_inputs(fields,seed);

  constexpr auto nan = std::numeric_limits<Real>::quiet_NaN();
  for (const std::string pos : {"top","bot"}) {
    for (const auto& ph : phases) {
      const auto lbl = pos + ", phase=" + ph;
      const auto w = weights_reference(pos,ph,fields);
      fields_t xf = fields;

      // X may not be finite where the weight is zero (e.g., divided by a zero
      // cloud fraction); it must not reach the sum there
      auto X = fields.at("nc").clone("X",CloneFlags::CopyData);
      auto X_h = host_view(X);
      int nzero = 0;
      for (int i=0; i<static_cast<int>(X_h.extent(0)); ++i) {
        for (int k=0; k<static_cast<int>(X_h.extent(1)); ++k) {
          if (w(i,k)==0) {
            X_h(i,k) = nan;
            ++nzero;
          }
        }
      }
      REQUIRE (nzero>0);
      X.sync_to_dev();
      xf.emplace("X",X);
      check(sample_diag("X",pos,ph,xf,grid),sample_reference("X",pos,ph,xf),lbl + ", NaN where unweighted");

      // A masked X counts only where it is valid
      auto Y = fields.at("nc").clone("Y",CloneFlags::CopyData);
      auto& m = Y.create_valid_mask();
      auto m_h = m.get_view<int**,Host>();
      for (int i=0; i<static_cast<int>(m_h.extent(0)); ++i) {
        for (int k=0; k<static_cast<int>(m_h.extent(1)); ++k) {
          m_h(i,k) = (i+k)%3!=0;
        }
      }
      m.sync_to_dev();
      xf.emplace("Y",Y);
      const auto ref_all = sample_reference("nc",pos,ph,fields);
      const auto ref_y   = sample_reference("Y",pos,ph,xf);
      int ndiff = 0;
      for (int i=0; i<static_cast<int>(ref_y.extent(0)); ++i) {
        ndiff += ref_y(i)!=ref_all(i);
      }
      REQUIRE (ndiff>0);
      check(sample_diag("Y",pos,ph,xf,grid),ref_y,lbl + ", masked");
    }
  }
}

TEST_CASE("aerocom_cld_realistic") {
  ekat::Comm comm(MPI_COMM_WORLD);
  register_diagnostics();

  int seed = get_random_test_seed();

  // Enough columns to fill a GPU: with only a few, the threads working on one
  // column run in near lockstep, and a race between them can go unseen
  const int ncols = 4096;
  for (const int nlevs : {9, 72, 128}) {
    auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);
    auto fields = create_inputs(grid);
    for (const std::string pos : {"top","bot"}) {
      for (const auto& ph : phases) {
        const auto lbl = "nlevs=" + std::to_string(nlevs) + ", " + pos + ", phase=" + ph;
        set_realistic_inputs(fields,seed);
        auto w = create_weights(pos,ph,fields,grid);
        auto d = create_sampler("T_mid",w,fields,grid);
        auto t = t0;
        w->compute(t);
        d->compute(t);
        check(host_view(w->get()),weights_reference(pos,ph,fields),lbl + ", weights");
        auto first = d->get().clone(CloneFlags::CopyData);
        first.sync_to_host();
        check(first.get_view<const Real*,Host>(),sample_reference("T_mid",pos,ph,fields),lbl);

        // Same inputs, new time stamp: not a single bit may change
        t += 1;
        w->compute(t);
        d->compute(t);
        auto second = d->get();
        second.sync_to_host();
        const auto v1 = first.get_view<const Real*,Host>();
        const auto v2 = second.get_view<const Real*,Host>();
        int ndiff = 0;
        for (int i=0; i<static_cast<int>(v1.extent(0)); ++i) {
          ndiff += v1(i)!=v2(i);
        }
        INFO (lbl + ": " + std::to_string(ndiff) + " values changed between two computes");
        REQUIRE (ndiff==0);
      }
      const auto lbl = "nlevs=" + std::to_string(nlevs) + ", " + pos + ", mask";
      check(sample_diag("mask",pos,"none",fields,grid),sample_reference("mask",pos,"none",fields),lbl);
    }
  }
}

}  // namespace scream
