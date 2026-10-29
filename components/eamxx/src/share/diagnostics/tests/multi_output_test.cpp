#include "catch2/catch.hpp"

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_units.hpp>

namespace scream {

// A diag with three outputs on different layouts, one of them masked:
//  - MOT_prof   (COL,LEV): a+b
//  - MOT_col    (COL):     the vertical sum of a
//  - MOT_masked (COL):     a at the first level, valid where b>0
class MOTest : public AbstractDiagnostic
{
public:
  MOTest (const ekat::Comm& comm, const ekat::ParameterList& params,
          const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    m_field_in_names = {"a","b"};

    const auto& gn = m_grid->name();
    m_prof   = Field(FieldIdentifier("MOT_prof",  m_grid->get_3d_scalar_layout(LEV),K,gn),true);
    m_col    = Field(FieldIdentifier("MOT_col",   m_grid->get_2d_scalar_layout(),K,gn),true);
    m_masked = Field(FieldIdentifier("MOT_masked",m_grid->get_2d_scalar_layout(),K,gn),true);
    m_masked.create_valid_mask();
    m_masked.get_header().set_may_be_filled(true);

    add_output(m_prof);
    add_output(m_col);
    add_output(m_masked);
  }

  std::string name () const override { return "MOTest"; }

  int get_num_evaluations () const { return m_num_evaluations; }

protected:
  void compute_impl () override {
    const auto& a = m_fields_in.at("a");
    const auto& b = m_fields_in.at("b");
    a.sync_to_host();
    b.sync_to_host();
    auto a_h    = a.get_view<const Real**,Host>();
    auto b_h    = b.get_view<const Real*,Host>();
    auto prof_h = m_prof.get_view<Real**,Host>();
    auto col_h  = m_col.get_view<Real*,Host>();
    auto msk_h  = m_masked.get_view<Real*,Host>();
    auto mask_h = m_masked.get_valid_mask().get_view<int*,Host>();
    const int ncols = a_h.extent(0);
    const int nlevs = a_h.extent(1);
    for (int icol=0; icol<ncols; ++icol) {
      col_h(icol) = 0;
      for (int ilev=0; ilev<nlevs; ++ilev) {
        prof_h(icol,ilev) = a_h(icol,ilev) + b_h(icol);
        col_h(icol) += a_h(icol,ilev);
      }
      mask_h(icol) = b_h(icol)>0 ? 1 : 0;
      msk_h(icol)  = b_h(icol)>0 ? a_h(icol,0) : constants::fill_value<Real>;
    }
    m_prof.sync_to_dev();
    m_col.sync_to_dev();
    m_masked.sync_to_dev();
    m_masked.get_valid_mask().sync_to_dev();
    ++m_num_evaluations;
  }

  Field m_prof, m_col, m_masked;
  int m_num_evaluations = 0;
};

// A diag with one output, set the usual way (m_diagnostic_output in the ctor)
class SOTest : public AbstractDiagnostic
{
public:
  SOTest (const ekat::Comm& comm, const ekat::ParameterList& params,
          const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    m_field_in_names = {"b"};
    FieldIdentifier fid("SOT",m_grid->get_2d_scalar_layout(),ekat::units::K,m_grid->name());
    m_diagnostic_output = Field(fid,true);
    m_diagnostic_output.create_valid_mask(Field::MaskInit::Valid);
  }

  std::string name () const override { return "SOTest"; }

protected:
  void compute_impl () override {
    m_diagnostic_output.deep_copy(m_fields_in.at("b"));
  }
};

// The 'misuse' param picks a rule of the base class for this diag to break
class BadMOTest : public AbstractDiagnostic
{
public:
  BadMOTest (const ekat::Comm& comm, const ekat::ParameterList& params,
             const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    m_misuse = m_params.get<std::string>("misuse");
    if (m_misuse=="add_in_init") {
      add_output(make("x",m_grid));
    } else if (m_misuse=="both") {
      add_output(make("x",m_grid));
      m_diagnostic_output = make("y",m_grid);
    } else if (m_misuse=="duplicate") {
      add_output(make("x",m_grid));
      add_output(make("X",m_grid));
    } else if (m_misuse=="wrong_grid") {
      auto other = create_point_grid("other",m_grid->get_num_global_dofs(),
                                     m_grid->get_num_vertical_levels(),m_comm);
      add_output(make("x",other));
    } else if (m_misuse=="unallocated") {
      FieldIdentifier fid("x",m_grid->get_2d_scalar_layout(),ekat::units::none,m_grid->name());
      add_output(Field(fid));
    } else if (m_misuse=="unshared_mask" or m_misuse=="shared_mask") {
      auto x = make("x",m_grid);
      auto y = make("y",m_grid);
      x.create_valid_mask("region_mask");
      if (m_misuse=="shared_mask") {
        y.set_valid_mask(x.get_valid_mask());
      } else {
        y.create_valid_mask("region_mask");
      }
      add_output(x);
      add_output(y);
    }
  }

  std::string name () const override { return "BadMOTest"; }

protected:
  void initialize_impl () override {
    if (m_misuse=="add_in_init") {
      add_output(make("y",m_grid));
    } else if (m_misuse=="two_in_init") {
      add_output(make("x",m_grid));
      add_output(make("y",m_grid));
    }
  }

  void compute_impl () override {}

  static Field make (const std::string& name, const std::shared_ptr<const AbstractGrid>& grid) {
    FieldIdentifier fid(name,grid->get_2d_scalar_layout(),ekat::units::none,grid->name());
    return Field(fid,true);
  }

  std::string m_misuse;
};

TEST_CASE ("multi_output") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  int seed = get_random_test_seed(&comm);

  const int ncols = 6;
  const int nlevs = 7;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  util::TimeStamp t0({2022,1,1},{0,0,0});
  const util::TimeStamp t1 = t0 + 600;

  // Inputs. b alternates sign across columns, so the mask is valid in some of them
  Field a(FieldIdentifier("a",grid->get_3d_scalar_layout(LEV),K,grid->name()),true);
  Field b(FieldIdentifier("b",grid->get_2d_scalar_layout(),K,grid->name()),true);
  randomize_uniform(a,seed++,200,300);
  auto b_h = b.get_view<Real*,Host>();
  for (int icol=0; icol<grid->get_num_local_dofs(); ++icol) {
    b_h(icol) = (icol%2==0 ? 1 : -1)*(icol+1);
  }
  b.sync_to_dev();
  a.get_header().get_tracking().update_time_stamp(t0);
  b.get_header().get_tracking().update_time_stamp(t0);

  ekat::ParameterList params;

  SECTION ("several_outputs") {
    auto diag = std::make_shared<MOTest>(comm,params,grid);

    // All outputs are known right after construction, in declaration order
    REQUIRE (diag->get_outputs().size()==3);
    REQUIRE (diag->get_outputs()[0].name()=="MOT_prof");
    REQUIRE (diag->get_outputs()[1].name()=="MOT_col");
    REQUIRE (diag->get_outputs()[2].name()=="MOT_masked");

    // get() is for diags with one output, and the error lists the outputs
    REQUIRE_THROWS_WITH (diag->get(),
        Catch::Contains("MOT_prof, MOT_col, MOT_masked") and Catch::Contains("get(name)"));
    REQUIRE (diag->get("mot_col").name()=="MOT_col");
    REQUIRE_THROWS_WITH (diag->get("nope"),
        Catch::Contains("requested output: nope") and
        Catch::Contains("available outputs: MOT_prof, MOT_col, MOT_masked"));

    diag->set_input_field(a);
    diag->set_input_field(b);
    REQUIRE_THROWS (diag->compute(t0));  // Not yet initialized
    diag->initialize();
    REQUIRE (diag->get_outputs().size()==3);
    REQUIRE_THROWS (diag->initialize());  // Already initialized

    // One compute evaluates all outputs, and stamps them all
    diag->compute(t0);
    REQUIRE (diag->get_num_evaluations()==1);
    for (const auto& f : diag->get_outputs()) {
      REQUIRE (f.get_header().get_tracking().get_time_stamp()==t0);
    }

    // Same inputs and time: no evaluation. New input time: one more.
    diag->compute(t0);
    REQUIRE (diag->get_num_evaluations()==1);
    a.get_header().get_tracking().update_time_stamp(t1);
    b.get_header().get_tracking().update_time_stamp(t1);
    diag->compute(t1);
    REQUIRE (diag->get_num_evaluations()==2);
    for (const auto& f : diag->get_outputs()) {
      REQUIRE (f.get_header().get_tracking().get_time_stamp()==t1);
    }

    // Values
    auto prof = diag->get("MOT_prof");
    auto col  = diag->get("MOT_col");
    auto msk  = diag->get("MOT_masked");
    prof.sync_to_host();
    col.sync_to_host();
    msk.sync_to_host();
    msk.get_valid_mask().sync_to_host();
    REQUIRE (msk.get_valid_mask().name()=="MOT_masked_valid_mask");
    REQUIRE (msk.get_header().may_be_filled());
    auto a_h    = a.get_view<const Real**,Host>();
    auto prof_h = prof.get_view<const Real**,Host>();
    auto col_h  = col.get_view<const Real*,Host>();
    auto msk_h  = msk.get_view<const Real*,Host>();
    auto mask_h = msk.get_valid_mask().get_view<const int*,Host>();
    int nvalid = 0;
    for (int icol=0; icol<grid->get_num_local_dofs(); ++icol) {
      Real sum = 0;
      for (int ilev=0; ilev<nlevs; ++ilev) {
        REQUIRE (prof_h(icol,ilev)==Approx(a_h(icol,ilev)+b_h(icol)));
        sum += a_h(icol,ilev);
      }
      REQUIRE (col_h(icol)==Approx(sum));
      if (b_h(icol)>0) {
        REQUIRE (mask_h(icol)==1);
        REQUIRE (msk_h(icol)==Approx(a_h(icol,0)));
        ++nvalid;
      } else {
        REQUIRE (mask_h(icol)==0);
        REQUIRE (msk_h(icol)==constants::fill_value<Real>);
      }
    }
    // Both branches were exercised
    REQUIRE (nvalid>0);
    REQUIRE (nvalid<grid->get_num_local_dofs());
  }

  SECTION ("no_renaming") {
    // Several outputs cannot be renamed: one is selected by name, then aliased
    params.set<std::string>("output_field_name","foo");
    auto diag = std::make_shared<MOTest>(comm,params,grid);
    diag->set_input_field(a);
    diag->set_input_field(b);
    REQUIRE_THROWS_WITH (diag->initialize(),
        Catch::Contains("requested name: foo") and Catch::Contains("<name> := <output>"));
  }

  SECTION ("misuse") {
    auto create = [&](const std::string& misuse) {
      params.set("misuse",misuse);
      return std::make_shared<BadMOTest>(comm,params,grid);
    };

    // Outputs must be valid and uniquely named
    REQUIRE_THROWS (create("duplicate"));
    REQUIRE_THROWS (create("wrong_grid"));
    REQUIRE_THROWS (create("unallocated"));

    // Several outputs must all be declared in the ctor
    REQUIRE_THROWS (create("add_in_init")->initialize());
    REQUIRE_THROWS (create("two_in_init")->initialize());

    // Use add_output or m_diagnostic_output, not both
    REQUIRE_THROWS (create("both")->initialize());

    // Masks with the same name must be the same mask
    REQUIRE_THROWS (create("unshared_mask")->initialize());
    REQUIRE_NOTHROW (create("shared_mask")->initialize());
  }

  SECTION ("one_output") {
    // A diag with one output works as it always did
    auto diag = std::make_shared<SOTest>(comm,params,grid);
    REQUIRE (diag->get_outputs().empty());
    REQUIRE (diag->get().name()=="SOT");
    REQUIRE (diag->get("sot").name()=="SOT");
    diag->set_input_field(b);
    diag->initialize();
    REQUIRE (diag->get_outputs().size()==1);
    REQUIRE (diag->get_outputs()[0].is_aliasing(diag->get()));

    diag->compute(t0);
    REQUIRE (diag->get().get_header().get_tracking().get_time_stamp()==t0);
    REQUIRE (views_are_equal(diag->get(),b));

    // It can be renamed, and the alias keeps the mask
    params.set<std::string>("output_field_name","foo");
    auto renamed = std::make_shared<SOTest>(comm,params,grid);
    renamed->set_input_field(b);
    renamed->initialize();
    REQUIRE (renamed->get_outputs().size()==1);
    REQUIRE (renamed->get_outputs()[0].name()=="foo");
    REQUIRE (renamed->get().name()=="foo");
    REQUIRE (renamed->get("FOO").name()=="foo");
    REQUIRE (renamed->get().get_valid_mask().name()=="SOT_valid_mask");
    REQUIRE_THROWS (renamed->get("SOT"));
  }
}

} // namespace scream
