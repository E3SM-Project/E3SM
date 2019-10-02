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
  m_grids_manager->build_grid(gm_params.get<std::string>("Reference Grid"));

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

  // Set all the fields in the processes needing them (before, they only had ids)
  // Input fields will be handed to the processes as const
  const auto& inputs  = m_atm_process_group->get_required_fields();
  const auto& outputs = m_atm_process_group->get_computed_fields();
  for (const auto& id : inputs) {
    m_atm_process_group->set_required_field(m_device_field_repo.get_field(id).get_const());
  }
  // Internal fields are fields that the atm proc group both computes and requires
  // (in that order). These are present only for sequential splitting
  for (const auto& id : m_atm_process_group->get_internal_fields()) {
    m_atm_process_group->set_internal_field(m_device_field_repo.get_field(id));
  }
  // Output fields are handed to the processes as writable
  for (const auto& id : outputs) {
    m_atm_process_group->set_computed_field(m_device_field_repo.get_field(id));
  }

  // Note: remappers should be setup *after* fields have been set in the atm processes,
  //       just in case some atm proc sets some extra data in the field header,
  //       that some remappers may actually need.
  m_atm_process_group->setup_remappers(m_device_field_repo);

  // Initialize the processes
  m_atm_process_group->initialize(t0);

  // Set time steamp t0 to all fields
  for (auto& field_map_it : m_device_field_repo) {
    for (auto& f_it : field_map_it.second) {
      f_it.second.get_header().get_tracking().update_time_stamp(t0);
    }
  }

#ifdef SCREAM_DEBUG
  create_bkp_device_field_repo();
  m_atm_process_group->set_field_repos(m_device_field_repo,m_bkp_device_field_repo);
#endif
}

void AtmosphereDriver::run (const Real dt) {
  // Make sure the end of the time step is after the current start_time
  scream_require_msg (dt>0, "Error! Input time step must be positive.\n");

  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run(dt);

  // Update current time stamps
  m_current_ts += dt;
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );

  m_device_field_repo.clean_up();
#ifdef SCREAM_DEBUG
  m_bkp_device_field_repo.clean_up();
#endif
}

#ifdef SCREAM_DEBUG
void AtmosphereDriver::create_bkp_device_field_repo () {
  m_bkp_device_field_repo.registration_begins();
  for (const auto& it : m_device_field_repo) {
    for (const auto& id_field : it.second) {
      const auto& id = id_field.first;
      const auto& f = id_field.second;
      const auto& groups = f.get_header().get_tracking().get_groups_names();
      // Unfortunately, set<string> and set<CaseInsensitiveString>
      // are unrelated types for the compiler
      std::set<std::string> grps;
      for (const auto& group : groups) {
        grps.insert(group);
      }
      m_bkp_device_field_repo.register_field(id,grps);
    }
  }
  m_bkp_device_field_repo.registration_ends();

  // Deep copy the fields
  for (const auto& it : m_device_field_repo) {
    for (const auto& id_field : it.second) {
      const auto& id = id_field.first;
      const auto& f  = id_field.second;
      auto src = f.get_view();
      auto tgt = m_bkp_device_field_repo.get_field(id).get_view();

      Kokkos::deep_copy(tgt,src);
    }
  }
}
#endif

}  // namespace control
}  // namespace scream
