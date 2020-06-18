#include "control/atmosphere_driver.hpp"

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/field/field_initializer.hpp"
#include "share/field/field_utils.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/string_utils.hpp"

namespace scream {

namespace control {

/*
 * IMPORTANT: read carefully this banner before attempting any change to the initialize method!
 *
 * The order in which the AD initializes all its internal stuff matters. Here's the order in
 * which operation currently happen, and why. If you alter the method, then a) make sure you
 * are not breaking any logic here explained (or else fix it!), and b) modify this banner to
 * update the explanation of the initialization sequence.
 *
 *  1) Create all atm processes. Each proc is allowed to start some sort of setup during creation,
 *     but will not be able to fully set up its required/computed fields, due to lack of grids info.
 *     However, and this is important, each process MUST establish what grid it needs.
 *  2) Create the grid manager, and query the atm procs for the grids they need. The GM will then
 *     proceed to build those grids (and only those grids).
 *  3) The GM is passed back to the atm procs, which can grab the needed grids, from which they can
 *     get the information needed to complete the setup of the FieldIdentifiers of their fields
 *     (both required and computed). Their field identifiers MUST be completed upon return from
 *     the 'set_grids' method.
 *  4) Register all fields from all atm procs inside the field manager (or field repo, whatever you
 *     want to call it).
 *  5) Set all the fields into the atm procs. Before this point, all the atm procs had were the
 *     FieldIdentifiers for their input/output fields. Now, we pass actual Field objects to them,
 *     where both the data (Kokkos::View) and metadata (FieldHeader) inside will be shared across
 *     all processes using the field. This allow data and metadata to be always in sync.
 *     Note: output fields are passed to an atm proc as read-write (i.e., non-const data type),
 *           while input fields are passed as read-only (i.e., const data type). Yes, the atm proc
 *           could cheat, and cast away the const, but we can't prevent that. However, in debug builds,
 *           we store 2 copies of each field, and use the extra copy to check, at run time, that
 *           no process alters the values of any of its input fields.
 *     Note: when setting input fields in an atm proc, the atm proc is also allowed to specify how
 *           the field should be initialized. We only allow ONE process to attempt to do that,
 *           so atm procs should ALWAYS check that no initialization is already set for that field.
 *           Also, notice that the atm proc has no idea whether an initializer will be needed for
 *           this field or not. Checking whether there are providers already registered for it
 *           seems plausible, but it relies on atm procs being parsed in the correct order. Although
 *           this IS currently the case, we don't want to rely on it. Therefore, an atm proc should
 *           be aware that, even if it sets a certain type of initialization for a field, that
 *           initialization may not be needed. E.g., proc A needs field F, and sets it to be init-ed
 *           to 0. However, upon completion of the atm setup, proc B computes field F, and proc B
 *           is evaluated *before* proc A in the atm DAG. Therefore, no initialization is needed
 *           for field F.
 *  6) The atm procs group builds the remappers. Remappers are objects that map fields from one grid
 *     to another. Currently, we perform remappings only to/from a reference grid (set in the GM).
 *     Before calling each atm proc in the group, the AtmosphereProcessGroup (APG) class will
 *     remap all its inputs from the ref grid to the grid(s) needed by the atm proc. Similarly,
 *     upon completion of the atm proc run call, all the outputs will be remapped back to the
 *     ref grid.
 *  7) All the atm process are initialized. During this call, atm process are able to set up
 *     all the internal structures that they were not able to init previously.
 *  8) All the atm inputs are initialized, according to the initialization type that was
 *     specified during step 5. This MUST happen AFTER step 7, since it MIGHT depend on
 *     an atm proc being fully init-ed. On the other hand, we cannot make this happen
 *     inside step 7 (assuming a field is marked as to be init-ed by a specific atm proc),
 *     since the atm proc does not yet know if the field needs to be inited or not (we haven't
 *     yet built/parsed the atm DAG).
 *  9) We build the atm DAG, and inspect its consistency. Each input of each atm proc must be
 *     computed by some other atm proc (possibly at the previous time step) or have a FieldInitializer
 *     set (could be the case for some geometry-related field, constant during the whole simulation).
 *     Additionally, every field that comes from the previous time step (i.e., an input of the
 *     atmosphere as a whole), MUST have an initializer structure set. If some dependency is not
 *     met, the AD will error out, and store the atm DAG in a file.
 * 10) Finally, set the initial time stamp on all fields, and perform some debug structure setup.
 *
 */


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
  m_grids_manager->build_grids(m_atm_process_group->get_required_grids(),
                               gm_params.get<std::string>("Reference Grid"));

  // Set the grids in the processes. Do this by passing the grids manager.
  // Each process will grab what they need
  m_atm_process_group->set_grids(m_grids_manager);

  // By now, the processes should have fully built the ids of their
  // required/computed fields. Let them register them in the repo
  m_device_field_repo.registration_begins();
  m_atm_process_group->register_fields(m_device_field_repo);
  m_device_field_repo.registration_ends();

  // Set all the fields in the processes needing them (before, they only had ids)
  // Input fields will be handed to the processes as const
  const auto& inputs  = m_atm_process_group->get_required_fields();
  const auto& outputs = m_atm_process_group->get_computed_fields();
  for (const auto& id : inputs) {
    m_atm_process_group->set_required_field(m_device_field_repo.get_field(id).get_const());
  }
  // Internal fields are fields that the atm proc group both computes and requires
  // (in that order). These are present only in case of sequential splitting
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

  // Initialize atm inputs
  init_atm_inputs ();

  // Now we can inspect the dag, including also checking that
  // all fields are correctly inited. Any unmet dependency in
  // the dag at this point has to be treated as an error.
  inspect_atm_dag ();

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

void AtmosphereDriver::init_atm_inputs () {
  const auto& atm_inputs = m_atm_process_group->get_required_fields();
  for (const auto& id : atm_inputs) {
    auto& f = m_device_field_repo.get_field(id);
    auto init_type = f.get_header_ptr()->get_tracking().get_init_type();
    if (init_type!=InitType::None) {
      auto initializer = f.get_header_ptr()->get_tracking().get_initializer().lock();
      scream_require_msg (static_cast<bool>(initializer),
                          "Error! Field '" + f.get_header().get_identifier().name() + "' has initialization type '" + e2str(init_type) + "',\n" +
                          "       but its initializer pointer is not valid.\n");

      initializer->add_field(f);
      m_field_initializers.insert(initializer);
    }
  }

  // Now loop over all the initializers, and make them init their fields.
  for (auto it : m_field_initializers) {
    auto initializer = it.lock();
    if (initializer->get_inited_fields().size()>0) {
      initializer->initialize_fields();
    }
  }
}

void AtmosphereDriver::inspect_atm_dag () {

  // First, process the dag
  AtmProcDAG dag;
  dag.create_dag(*m_atm_process_group);

  for (const auto& it : m_field_initializers) {
    dag.add_field_initializer(*it.lock());
  }

  auto& deb_pl = m_atm_params.sublist("Debug");
  if (dag.has_unmet_dependencies()) {
    const int err_verb_lev = deb_pl.get<int>("Atmosphere DAG Verbosity Level",int(AtmProcDAG::VERB_MAX));
    dag.write_dag("error_atm_dag.dot",err_verb_lev);
    scream_error_msg("Error! There are unmet dependencies in the atmosphere internal dag.\n"
                     "       Use the graphviz package to inspect the dependency graph:\n"
                     "    \n"
                     "    $ dot -Tjpg -o error_atm_dag.jpg error_atm_dag.dot\n"
                     "    $ eog error_atm_dag.dot\n");
  }

  // If requested, write a dot file for visualization
  const int verb_lev = deb_pl.get<int>("Atmosphere DAG Verbosity Level",0);
  dag.write_dag("scream_atm_dag.dot",verb_lev);
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
