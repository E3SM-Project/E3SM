#include "control/atmosphere_driver.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"

namespace scream {

// === A dummy atm process, on Physics grid === //

template<typename DeviceType, int PackSize, bool forward>
class DummyProcess : public scream::AtmosphereProcess {
public:
  using device_type = DeviceType;
  using exec_space  = typename device_type::execution_space;

  DummyProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
   : m_comm(comm)
  {
    m_iter = 0;
    m_params = params;
    m_id = comm.rank();

    if (forward) {
      m_name = "Physics_fwd";
    } else {
      m_name = "Physics_bwd";
    }
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::set<std::string> get_required_grids () const {
    // TODO: define what grid the coupling runs on. Check with MOAB folks.
    static std::set<std::string> s;
    s.insert(m_name);
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
    m_grid = grids_manager->get_grid(m_name);

    auto num_cols = m_grid->get_num_local_dofs();

    std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
    std::vector<int> dims = {num_cols, m_params.get<int>("Number of vector components")};
    FieldLayout layout (tags,dims);

    std::string in_name = "field_";
    std::string out_name = "field_";
    if (forward) {
      in_name  += "0";
      out_name += "1";
    } else {
      in_name  += "1";
      out_name += "0";
    }

    m_input_fids.emplace(in_name,layout,ekat::units::m,m_grid->name());
    m_output_fids.emplace(out_name,layout,ekat::units::m,m_grid->name());
  }

  void initialize_impl (const util::TimeStamp& t0) {
  }

  void run_impl (const Real dt) {
    auto in = m_input.get_view();
    auto out = m_output.get_view();
    auto iter = m_iter % 4;
    Kokkos::parallel_for(Kokkos::RangePolicy<exec_space>(0,in.size()),
      KOKKOS_LAMBDA(const int i) {
        switch (iter) {
          case 0: out(i) = in(i) + 2.0; break;
          case 1: out(i) = in(i) * 2.0; break;
          case 2: out(i) = in(i) - 2.0; break;
          case 3: out(i) = in(i) / 2.0; break;
          default:
            EKAT_KERNEL_ASSERT(false);
        }
    });
    Kokkos::fence();

    ++m_iter;
    auto ts = timestamp();
    ts += dt;
    m_output.get_header().get_tracking().update_time_stamp(ts);
  }

  // Clean up
  void finalize_impl ( ) {}

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const {
    using pack_type = ekat::Pack<Real,PackSize>;
    field_repo.template register_field<pack_type>(*m_input_fids.begin());
    field_repo.template register_field<pack_type>(*m_output_fids.begin());
  }

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_input_fids; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_output_fids; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& f) {
    m_input = f;
  }
  void set_computed_field_impl (const Field<      Real, device_type>& f) {
    m_output = f;
  }

  int m_iter;

  std::set<FieldIdentifier> m_input_fids;
  std::set<FieldIdentifier> m_output_fids;

  Field<const Real,device_type> m_input;
  Field<Real,device_type>       m_output;

  std::shared_ptr<const AbstractGrid>   m_grid;

  std::string m_name;

  ekat::ParameterList m_params;
  int           m_id;

  ekat::Comm    m_comm;
};

} // namespace scream
