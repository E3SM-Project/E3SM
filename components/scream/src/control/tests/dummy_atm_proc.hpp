#include "control/atmosphere_driver.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"

namespace scream {

// === A dummy atm process, on Physics grid === //

class DummyProcess : public scream::AtmosphereProcess {
public:
  using exec_space  = typename DefaultDevice::execution_space;

  enum DummyType {
    A2G,
    G2G,
    G2A
  };

  DummyProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
   : m_comm(comm)
  {
    m_params = params;
    m_name = m_params.get<std::string>("Sub Name");
    if (m_name=="Group to Group") {
      m_params.set<std::string>("Grid Name","Point Grid B");
      m_dummy_type = G2G;
    } else {
      m_params.set<std::string>("Grid Name","Point Grid A");
      m_dummy_type = (m_name=="A to BC") ? A2G : G2A;
    }
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::set<std::string> get_required_grids () const {
    return std::set<std::string> {m_params.get<std::string>("Grid Name")};
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
    using namespace ShortFieldTagsNames;

    m_grid = grids_manager->get_grid(m_params.get<std::string>("Grid Name"));

    const auto num_cols = m_grid->get_num_local_dofs();
    const auto num_levs = m_grid->get_num_vertical_levels();

    std::vector<FieldTag> tags;
    std::vector<int> dims;
    if (m_grid->name()=="Point Grid A") {
      tags = {COL,LEV};
      dims = {num_cols, num_levs};
    } else {
      tags = {LEV,COL};
      dims = {num_levs,num_cols};
    }
    FieldLayout layout (tags,dims);
    // To test vector input
    FieldLayout layout_vec ( {COL,CMP,LEV}, {num_cols,2,num_levs} );

    if (m_dummy_type==A2G) {
      add_required_field("A",layout,ekat::units::m,m_grid->name());
      add_computed_field("B",layout,ekat::units::m,m_grid->name());
      add_computed_field("C",layout,ekat::units::m,m_grid->name());
      add_required_field("D",layout_vec,ekat::units::m,m_grid->name());
      add_required_field("E",layout,ekat::units::m,m_grid->name());
    } else if (m_dummy_type == G2A) {
      add_computed_field("A",layout,ekat::units::m,m_grid->name());
      add_required_group("The Group",m_grid->name());
    } else {
      add_updated_group("The Group",m_grid->name());
    }
  }

  // Register all fields in the proper fm(s).
  // Note: field_mgrs[grid_name] is the FM on grid $grid_name
  void register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const {
    auto fm = field_mgrs.at(m_grid->name());
    if (m_dummy_type==A2G) {
      for (const auto& fid : get_required_fields()) {
        fm->register_field(fid);
      }
      for (const auto& fid : get_computed_fields()) {
        fm->register_field(fid,"The Group");
      }
    } else if (m_dummy_type == G2A) {
      for (const auto& fid : get_computed_fields()) {
        fm->register_field(fid);
      }
    }
  }

  void set_required_group (const FieldGroup<const Real>& field_group) {
    EKAT_REQUIRE_MSG (m_dummy_type==G2A,
                      "Error! This atmosphere process does not require a group of fields.\n");

    for (const auto& it : field_group.m_fields) {
      const auto& f = it.second;
      const auto& fid = f->get_header().get_identifier();
      m_inputs.emplace(fid.name(),*f);
      add_required_field(fid);
    }
  }
  void set_updated_group (const FieldGroup<Real>& field_group) {
    EKAT_REQUIRE_MSG (m_dummy_type==G2G,
                      "Error! This atmosphere process does not require a group of fields.\n");

    for (const auto& it : field_group.m_fields) {
      const auto& f = it.second;
      const auto& fid = f->get_header().get_identifier();
      add_required_field(fid);
      add_computed_field(fid);

      m_inputs.emplace(fid.name(),f->get_const());
      m_outputs.emplace(fid.name(),*f);
    }
  }

protected:

  void initialize_impl (const util::TimeStamp&) {
    // Do nothing
  }

// CUDA needs top level lambdas to be enclosed by a method that is public.
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void run_impl (const Real /* dt */) {
    const int ncols = m_grid->get_num_local_dofs();
    const int nlevs = m_grid->get_num_vertical_levels();
    auto policy = KokkosTypes<exec_space>::RangePolicy(0,ncols*nlevs);
    if (m_name=="A to BC") {
      const auto view_A = m_inputs["A"].get_reshaped_view<const Real**>();
      const auto view_B = m_outputs["B"].get_reshaped_view<Real**>();
      const auto view_C = m_outputs["C"].get_reshaped_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
        // A to BC is on grid A: (COL,LEV)
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        view_B(icol,ilev) = view_A(icol,ilev) / 2;
        view_C(icol,ilev) = view_A(icol,ilev) / 2;
      });
    } else if (m_name=="Group to Group") {
      const auto& fB = m_outputs["B"];
      const auto view_B = fB.get_reshaped_view<Real**>();
      const auto view_C = m_outputs["C"].get_reshaped_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
        // Group to Group is on grid B: (LEV,COL)
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        view_B(ilev,icol) = view_B(ilev,icol) / 2;
        view_C(ilev,icol) = view_C(ilev,icol) / 2;
      });
    } else {
      const auto view_B = m_inputs["B"].get_reshaped_view<const Real**>();
      const auto view_C = m_inputs["C"].get_reshaped_view<const Real**>();
      const auto view_A = m_outputs["A"].get_reshaped_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
        // Group to A is on grid A: (COL,LEV)
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        view_A(icol,ilev) = view_B(icol,ilev) + view_C(icol,ilev);
      });
    }
  }

protected:

  void finalize_impl () {
    // Do nothing
  }

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real>& f) {
    m_inputs[f.get_header().get_identifier().name()] = f;
  }
  void set_computed_field_impl (const Field<Real>& f) {
    m_outputs[f.get_header().get_identifier().name()] = f;
  }

  std::map<std::string,Field<const Real>>   m_inputs;
  std::map<std::string,Field<      Real>>   m_outputs;

  std::shared_ptr<const AbstractGrid>   m_grid;

  std::string m_name;

  ekat::ParameterList m_params;
  DummyType     m_dummy_type; 

  ekat::Comm    m_comm;
};

} // namespace scream
