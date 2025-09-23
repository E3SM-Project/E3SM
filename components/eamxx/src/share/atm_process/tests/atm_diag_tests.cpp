#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_parameter_list.hpp>

namespace scream {

// ========================= A dummy diagnostic ========================= //

class DiagSum : public AtmosphereDiagnostic {

public:
  DiagSum (const ekat::Comm& comm, const ekat::ParameterList&params)
    : AtmosphereDiagnostic(comm, params)
  {
    m_grid_name = params.get<std::string> ("grid_name");
  }

  // Return some sort of name, linked to PType
  std::string name () const { return "Sum diagnostic"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;

    const auto grid = gm->get_grid(m_grid_name);
    const auto lt = grid->get_2d_scalar_layout ();

    add_field<Required>("Field A",lt,K,m_grid_name);
    add_field<Required>("Field B",lt,K,m_grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, K, m_grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
  }

protected:

  void initialize_impl (const RunType /* run_type */ ) {}

  void finalize_impl () {}

  void compute_diagnostic_impl () {
    auto f_A = get_field_in("Field A", m_grid_name);
    auto f_B = get_field_in("Field B", m_grid_name);
    m_diagnostic_output.deep_copy(f_A);
    m_diagnostic_output.update(f_B,1,1);
  }

  std::string m_grid_name;
};

// ================================ TESTS ============================== //

TEST_CASE ("diagnostics") {

  //TODO: This test needs a field manager so that changes in Field A are seen everywhere.
  using namespace scream;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager
  const int nlcols = 3;
  const int nlevs = 10;
  auto grid = create_point_grid ("point_grid",nlcols*comm.size(),nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>(grid);

  // Create the sum diagnostic
  ekat::ParameterList params_sum("DiagSum");
  params_sum.set<std::string>("grid_name", "point_grid");
  auto diag_sum = std::make_shared<DiagSum>(comm,params_sum);
  diag_sum->set_grids(gm);

  std::map<std::string,Field> input_fields;
  for (const auto& req : diag_sum->get_required_field_requests()) {
    Field f(req.fid);
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    input_fields[f.name()] = f;
  }
  auto f_A = input_fields["Field A"];
  auto f_B = input_fields["Field B"];
  f_A.deep_copy(1.0);
  f_B.deep_copy(1.0);

  diag_sum->set_required_field(f_A.get_const());
  diag_sum->set_required_field(f_B.get_const());
  // diag fields are created INSIDE diags, not set from outside
  REQUIRE_THROWS(diag_sum->set_computed_field(f_A));

  diag_sum->initialize(t0,RunType::Initial);

  // Run the diagnostic
  diag_sum->compute_diagnostic();

  // Get diagnostics outputs
  const auto& f_sum      = diag_sum->get_diagnostic();

  f_B.update(f_A,1,1);
  REQUIRE (views_are_equal(f_B,f_sum));
}

} // empty namespace
