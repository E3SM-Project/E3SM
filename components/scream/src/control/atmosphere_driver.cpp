#include "atmosphere_driver.hpp"

#include <control/atmosphere_process_group.hpp>
#include <share/scream_assert.hpp>
#include <share/util/string_utils.hpp>

#include "control/surface_coupling.hpp"

#include "dynamics/atmosphere_dynamics.hpp"

namespace scream {
namespace control {

void AtmosphereDriver::initialize (const Comm& atm_comm, const ParameterList& params /* inputs? */ ) {
  m_atm_comm = atm_comm;
  m_atm_params = params;

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_params.sublist("Atmosphere Processes"));

  // Initialize the processes
  m_atm_process_group->initialize(m_atm_comm);

  // Let the processes register their fields in the repo
  m_device_field_repo.registration_begins();
  m_atm_process_group->register_fields(m_device_field_repo);

  // TODO: this is the correct location where you should make sure all the fields
  //       dimensions have been set. Some may have not been known at construction
  //       time: e.g., a field may have been created following the request of an
  //       atmosphere process which did not have access to all the dimensions
  //       information. Now, before you call registration_complete, you MUST
  //       make sure all dimensions are set.

  // Prohibit further additions to the repo, and allocate fields.
  m_device_field_repo.registration_ends();

  // TODO: this is a good place where we can insert a DAG analysis, to make sure all
  //       processes have their dependency met.

  // Set all the fields in the processes needing them (before, they only had headers)
  // Input fields will be handed to the processes as const
  const auto& inputs  = m_atm_process_group->get_required_fields();
  const auto& outputs = m_atm_process_group->get_computed_fields();
  for (const auto& id : inputs) {
    m_atm_process_group->set_required_field(m_device_field_repo.get_field(id));
  }
  // Output fields are handed to the processes as writable
  for (const auto& id : outputs) {
    m_atm_process_group->set_computed_field(m_device_field_repo.get_field(id));
  }
}

void AtmosphereDriver::run ( /* inputs ? */ ) {
  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run( /* inputs ? */ );
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );
}

// ==================== STUB for testing ================= //

int driver_stub() {
  // Test that we can at least create an AD...
  AtmosphereDriver ad;
  (void) ad;

  // Return the answer to life, universe, everything...
  return 42;
}

}  // namespace control
}  // namespace scream
