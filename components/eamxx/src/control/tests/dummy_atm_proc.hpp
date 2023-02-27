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
    : AtmosphereProcess(comm, params)
  {
    m_name = m_params.get<std::string>("Sub Name");
    if (m_name=="Group to Group") {
      m_dummy_type = G2G;
    } else {
      m_dummy_type = (m_name=="A to Group") ? A2G : G2A;
    }
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
    using namespace ShortFieldTagsNames;

    m_grid = grids_manager->get_grid(m_params.get<std::string>("Grid Name"));

    const auto num_cols = m_grid->get_num_local_dofs();
    const auto num_levs = m_grid->get_num_vertical_levels();

    FieldLayout layout ({COL,LEV},{num_cols,num_levs});
    // To test vector input
    FieldLayout layout_vec ( {COL,CMP,LEV}, {num_cols,2,num_levs} );

    if (m_dummy_type==A2G) {
      add_field<Required>("A",layout,ekat::units::m,m_grid->name());
      add_field<Computed>("B",layout,ekat::units::m,m_grid->name(),"The Group");
      add_field<Computed>("C",layout,ekat::units::m,m_grid->name(),"The Group");
      // These are not used at run time, but we use them to test
      // the initialization of IC fields
      add_field<Required>("V",layout_vec,ekat::units::m,m_grid->name());
      add_field<Required>("Z",layout,ekat::units::m,m_grid->name());
    } else if (m_dummy_type == G2A) {
      add_field<Computed>("A",layout,ekat::units::m,m_grid->name());
      add_group<Required>("The Group",m_grid->name());
    } else {
      add_field<Required>("B",layout,ekat::units::m,m_grid->name());
      add_group<Updated>("The Group",m_grid->name());
    }
  }

protected:

  void initialize_impl (const RunType /* run_type */) {
    // Do nothing
  }

// CUDA needs top level lambdas to be enclosed by a method that is public.
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void run_impl (const double /* dt */) {
    const int ncols = m_grid->get_num_local_dofs();
    const int nlevs = m_grid->get_num_vertical_levels();
    auto policy = KokkosTypes<exec_space>::RangePolicy(0,ncols*nlevs);
    if (m_name=="A to Group") {
      const auto view_A = get_field_in("A").get_view<const Real**>();
      const auto view_B = get_field_out("B").get_view<Real**>();
      const auto view_C = get_field_out("C").get_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        view_B(icol,ilev) = view_A(icol,ilev) / 2;
        view_C(icol,ilev) = view_A(icol,ilev) / 2;
      });
    } else if (m_name=="Group to Group") {
      const auto& g = get_group_out("The Group");
      const auto view_B = g.m_fields.at("B")->get_view<Real**>();
      const auto view_C = g.m_fields.at("C")->get_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
        const int icol = idx / nlevs;
        const int ilev = idx % nlevs;

        view_B(icol,ilev) = view_B(icol,ilev) / 2;
        view_C(icol,ilev) = view_C(icol,ilev) / 2;
      });
    } else {
      const auto& g = get_group_in("The Group");
      const auto view_B = g.m_fields.at("B")->get_view<const Real**>();
      const auto view_C = g.m_fields.at("C")->get_view<const Real**>();
      const auto view_A = get_field_out("A").get_view<Real**>();

      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
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

  std::shared_ptr<const AbstractGrid>   m_grid;

  std::string m_name;

  DummyType     m_dummy_type; 
};

} // namespace scream
