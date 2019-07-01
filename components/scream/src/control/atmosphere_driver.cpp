#include "atmosphere_driver.hpp"

#include "share/atmosphere_process_group.hpp"
#include "share/scream_assert.hpp"
#include "share/util/string_utils.hpp"

namespace scream {

namespace control {

void AtmosphereDriver::initialize (const Comm& atm_comm, const ParameterList& params) {
  m_atm_comm = atm_comm;
  m_atm_params = params;

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_comm,m_atm_params.sublist("Atmosphere Processes"));

  // Create the grids manager
  auto& gm_params = m_atm_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  m_grids_manager = GridsManagerFactory::instance().create(gm_type,m_atm_comm,gm_params);

  // Tell the grid manager to build all the grids required by the atm processes
  m_grids_manager->build_grids(m_atm_process_group->get_required_grids());

  // Initialize the processes
  m_atm_process_group->set_grid(m_grids_manager);

  // By now, the processes should have fully built the ids of their
  // required/computed fields. Let them register them in the repo
  m_device_field_repo.registration_begins();
  m_atm_process_group->register_fields(m_device_field_repo);
  m_device_field_repo.registration_ends();

  // TODO: this is a good place where we can insert a DAG analysis, to make sure all
  //       processes have their dependency met.

  // Set all the fields in the processes needing them (before, they only had ids)
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

  // Initialize the processes
  m_atm_process_group->initialize();
}

void AtmosphereDriver::run ( /* inputs ? */ ) {
  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run( /* inputs ? */ );
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );

  m_device_field_repo.clean_up();
}

}  // namespace control
}  // namespace scream
