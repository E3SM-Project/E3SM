#include "atmosphere_driver.hpp"

#include "share/atmosphere_process_group.hpp"
#include "share/scream_assert.hpp"
#include "share/util/string_utils.hpp"

namespace scream {

namespace control {

void AtmosphereDriver::initialize (const Comm& atm_comm,
                                   const ParameterList& params,
                                   const util::TimeStamp& t0)
{
  m_atm_comm = atm_comm;
  m_atm_params = params;
  m_current_ts = t0;

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_comm,m_atm_params.sublist("Atmosphere Processes"));

  // Create the grids manager
  auto& gm_params = m_atm_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  m_grids_manager = GridsManagerFactory::instance().create(gm_type,m_atm_comm,gm_params);

  // Tell the grid manager to build all the grids required
  // by the atm processes, as well as the reference grid
  m_grids_manager->build_grids(m_atm_process_group->get_required_grids());

  // Set the grids in the processes. Do this by passing the grids manager.
  // Each process will grab what they need
  m_atm_process_group->set_grids(m_grids_manager);

  // By now, the processes should have fully built the ids of their
  // required/computed fields. Let them register them in the repo
  m_device_field_repo.registration_begins();
  m_atm_process_group->register_fields(m_device_field_repo);
  m_device_field_repo.registration_ends();

  // TODO: this is a good place where we can insert a DAG analysis, to make sure all
  //       processes have their dependency met.

  // Note: remappers can be setup anytime after the field repo has allocated
  //       the fields. We put it here, but there's no dependency with the
  //       other atm process setup stages
  m_atm_process_group->setup_remappers(m_device_field_repo);

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
  m_atm_process_group->initialize(t0);

  // Set time steamp t0 to all fields
  for (auto& field_map_it : m_device_field_repo) {
    for (auto& f_it : field_map_it.second) {
      f_it.second.get_header().get_tracking().update_time_stamp(t0);
    }
  }

  // ===  Fields remapping structures === //
  //
  // Fields that are used by more than one process may need to
  // be remapped. Here remap could mean a change in the layout
  // (e.g., swapping dimensions), a change in the MPI distribution
  // (e.g., going from a 2 ranks distribution to a 4 ranks distribution),
  // or a change in grid (e.g. going from physics to dynamics).
  // In any case, a field with one layout has to be copied/interpolated
  // to a field with different layout.
  // The driver does this by using a 'default grid'. All fields are
  // kept up to date on the default grid, and remapped on the particular
  // process grid on a per-need basis. In other words, the AD makes sure
  // that upon return from an atm process run call, all its outputs
  // are up to date on the default grid; viceversa, before passing
  // the control to the process, it makes sure that its inputs are
  // up to date on the atm process grid. These two checks involve
  // a remap procedure in case the atm process does not use the
  // default grid. In particular, outputs are copied from the process
  // grid to the default grid, while inputs are copied from the
  // default grid to the process grid.

  // For fields that need to be used by more than one process,
  // we need to make sure that one copy is stored on the so-called
  // "default grid". The driver will use the 
}

void AtmosphereDriver::run (const double dt) {
  // Make sure the end of the time step is after the current start_time
  scream_require_msg (dt>0, "Error! Input time step must be positive.\n");

  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run(dt);

  // Update old and current time stamps
  m_old_ts = m_current_ts;
  m_current_ts += dt;

  // Check that the computed fields have un updated time stamp
  for (const auto& id : m_atm_process_group->get_computed_fields()) {
    const auto& f = m_device_field_repo.get_field(id);
    const auto& ts = f.get_header().get_tracking().get_time_stamp();
    scream_require_msg(ts==m_current_ts,
                       "Error! The time stamp of field '" + id.get_id_string() + "' has not been updated.\n"
                       "       Expected time stamp : " + m_current_ts.to_string() + "\n"
                       "       Field time stamp    : " + ts.to_string() + "\n");

  }
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );

  m_device_field_repo.clean_up();
}

}  // namespace control
}  // namespace scream
