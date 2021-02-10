#include "control/atmosphere_driver.hpp"

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/field/field_initializer.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"

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


void AtmosphereDriver::initialize (const ekat::Comm& atm_comm,
                                   const ekat::ParameterList& params,
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
  m_field_repo = std::make_shared<FieldRepository<Real> >();
  m_field_repo->registration_begins();
  m_atm_process_group->register_fields(*m_field_repo);
  register_groups();
  m_field_repo->registration_ends();

  // Set all the fields in the processes needing them (before, they only had ids)
  // Input fields will be handed to the processes as const
  const auto& inputs  = m_atm_process_group->get_required_fields();
  const auto& outputs = m_atm_process_group->get_computed_fields();
  for (const auto& id : inputs) {
    m_atm_process_group->set_required_field(m_field_repo->get_field(id).get_const());
  }
  // Output fields are handed to the processes as writable
  for (const auto& id : outputs) {
    m_atm_process_group->set_computed_field(m_field_repo->get_field(id));
  }
  // Set all groups of fields
  for (const auto& it : m_atm_process_group->get_required_groups()) {
    auto group = m_field_repo->get_const_field_group(it.first,it.second);
    m_atm_process_group->set_required_group(group);
  }
  for (const auto& it : m_atm_process_group->get_updated_groups()) {
    auto group = m_field_repo->get_field_group(it.first,it.second);
    m_atm_process_group->set_updated_group(group);
  }

  // Initialize the processes
  m_atm_process_group->initialize(t0);

  // Note 1: remappers should be setup *after* fields have been set in the atm processes,
  //       just in case some atm proc sets some extra data in the field header,
  //       that some remappers may actually need.
  // Note 2: setup the remapper *after* atm procs are fully init-ed. This allows
  //         dynamics to complete its initialization, which is needed by the
  //         phys-dyn remapper (which uses homme's mpi infrastructure).
  m_atm_process_group->setup_remappers(*m_field_repo);

  // Initialize atm inputs
  init_atm_inputs ();

  // Set time steamp t0 to all fields
  for (auto& field_map_it : *m_field_repo) {
    for (auto& f_it : field_map_it.second) {
      f_it.second->get_header().get_tracking().update_time_stamp(t0);
    }
  }

  // Now that all fields have been set and checked we can establish an output manager
  // for all output streams.
  if (m_atm_params.isSublist("Output Manager"))
  {
    auto& out_params = m_atm_params.sublist("Output Manager");
    m_output_manager.set_params(out_params);
    m_output_manager.set_comm(atm_comm);
    m_output_manager.set_grids(m_grids_manager);
    m_output_manager.set_repo(m_field_repo);
    // If not true then m_output_manager.init() won't do anything, leaving no_output flag as true and skipping output throughout simulation.
  }
  m_output_manager.init();

#if defined(SCREAM_CIME_BUILD)
  // If this is a CIME build, we need to prepare the surface coupling
  m_surface_coupling = std::make_shared<SurfaceCoupling>(m_grids_manager->get_reference_grid(),*m_field_repo);
#endif

  if (!m_surface_coupling) {
    // Standalone runs *may* require initialization of atm inputs.
    init_atm_inputs ();

    // In standalone runs, we already have everything we need for the dag,
    // so we can proceed to do its inspection. Any unmet dependency in
    // the dag at this point has to be treated as an error.
    inspect_atm_dag ();
  }

#ifdef SCREAM_DEBUG
  create_bkp_field_repo();
  m_atm_process_group->set_field_repos(*m_field_repo,m_bkp_field_repo);
#endif
}

void AtmosphereDriver::run (const Real dt) {
  // Make sure the end of the time step is after the current start_time
  EKAT_REQUIRE_MSG (dt>0, "Error! Input time step must be positive.\n");

  if (m_surface_coupling) {
    // Import fluxes from the component coupler (if any)
    m_surface_coupling->do_import();
  }

  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run(dt);

  // Update current time stamps
  m_current_ts += dt;

  // Update output streams
  m_output_manager.run(m_current_ts);

  if (m_surface_coupling) {
    // Export fluxes from the component coupler (if any)
    m_surface_coupling->do_export();
  }
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );

  // Finalize output streams, make sure files are closed
  m_output_manager.finalize();

  m_field_repo->clean_up();
#ifdef SCREAM_DEBUG
  m_bkp_field_repo.clean_up();
#endif
}

void AtmosphereDriver::register_groups () {
  using ci_string = ekat::CaseInsensitiveString;
  using ci_string_pair = std::pair<ci_string,ci_string>;

  // Given a list of group-grid pairs (A,B), make sure there is a copy
  // of each field in group A on grid B registered in the repo.
  auto lambda = [&](const std::set<ci_string_pair>& groups_grids) {
    const auto& groups_info = m_field_repo->get_groups_info();

    for (const auto& gg : groups_grids) {
      const auto& group = gg.first;
      const auto& grid = gg.second;

      // Lambda helper fcn, that register field $name with group $group on grid $grid
      // if not yet already registered
      auto register_if_not_there = [&](const std::string& name) {
        const auto& aliases = m_field_repo->get_aliases(name);
        EKAT_REQUIRE_MSG(aliases.size()>0,
          "Error! Something went wrong while looking for field '" << name << "' in the repo.\n");

        // Check if a copy of this field on the right grid is already registered.
        bool found = false;
        for (const auto& alias : aliases) {
          if (alias.first.get_grid_name()==grid) {
            found = true;
            break;
          }
        }

        if (!found) {
          // Field $name in group $group has no copy on grid $grid.
          // Lets take any fid in the repo for this field, and register
          // a copy of it on grid $grid. We can do this by creating
          // a remapper and using its capabilities.
          const auto& alias = *aliases.begin();
          const auto& fid = alias.first;
          auto r = m_grids_manager->create_remapper(fid.get_grid_name(),grid);
          auto f_units = fid.get_units();
          auto src_layout = fid.get_layout();
          auto tgt_layout = r->create_tgt_layout(src_layout);
          FieldIdentifier new_fid(name,tgt_layout,f_units,grid);
          m_field_repo->register_field(new_fid,group);
        }
      };

      auto group_it = groups_info.find(group);
      EKAT_REQUIRE_MSG(group_it!=groups_info.end(),
        "Error! Group '" << group << "' not found in the repo.\n");

      const auto& fnames = group_it->second->m_fields_names;
      for (const auto& name : fnames) {
        register_if_not_there(name);
      }

      if (group_it->second->m_bundled) {
        // The group was allocated as a single bundled field, with each
        // field in the group later subviewing the bundle.
        // We need to ensure the bundle also exists on $grid
        const auto& name = *fnames.begin();
        const auto& f = m_field_repo->get_field(name,grid);
        const auto& bundle_name = f.get_header().get_parent().lock()->get_identifier().name();
        register_if_not_there(bundle_name);
      }
    }
  };
  lambda( m_atm_process_group->get_required_groups() );
  lambda( m_atm_process_group->get_updated_groups() );
}

void AtmosphereDriver::init_atm_inputs () {
  const auto& atm_inputs = m_atm_process_group->get_required_fields();
  for (const auto& input : atm_inputs) {
    // We loop on all ids with the same name, cause it might be that some
    // initializer is set to init this field on a non-ref grid. In that case,
    // such non-ref field would *not* appear as an atm input (since it's the
    // output of a remapper).
    auto& aliases = m_field_repo->get_aliases(input.name());
    for (auto& it : aliases) {
      const auto& f = *it.second;
      const auto& hdr = f.get_header();
      const auto& id = hdr.get_identifier();
      auto init_type = hdr.get_tracking().get_init_type();
      if (init_type!=InitType::None && init_type!=InitType::Inherited) {
        auto initializer = hdr.get_tracking().get_initializer().lock();
        EKAT_REQUIRE_MSG (static_cast<bool>(initializer),
                            "Error! Field '" + id.name() + "' has initialization type '" + e2str(init_type) + "',\n" +
                            "       but its initializer pointer is not valid.\n");

        // If f is not on the ref grid, the initializer should init also the field
        // on the ref grid.
        const auto& ref_grid = m_grids_manager->get_reference_grid();
        const auto&     grid = m_grids_manager->get_grid(id.get_grid_name());
        if (grid->name()!=ref_grid->name()) {
          auto& f_ref = m_field_repo->get_field(id.name(),ref_grid->name());
          const auto& remapper = m_grids_manager->create_remapper(grid,ref_grid);
          initializer->add_field(f,f_ref,remapper);
        } else {
          initializer->add_field(f);
        }
        m_field_initializers.insert(initializer);
      }
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

void AtmosphereDriver::inspect_atm_dag () const {

  // Sanity checks
  EKAT_REQUIRE_MSG (!(m_surface_coupling && (m_field_initializers.size()>0)),
      "Error! You cannot have surface coupling *and* field initializers.\n"
      "       In CIME runs, you only have the former, while in standalone runs you can only have the latter.\n");

  EKAT_REQUIRE_MSG ( !m_surface_coupling || m_surface_coupling->get_repo_state()==RepoState::Closed,
      "Error! You cannot inspect the dag if the surface coupling has been fully initialized.\n");

  AtmProcDAG dag;

  // First, add all atm processes
  dag.create_dag(*m_atm_process_group,m_field_repo);

  // Then, add all field initializers (if any) and surface coupling (if any).
  // Recall that at most one of these two steps can be non-trivial.
  for (const auto& it : m_field_initializers) {
    dag.add_field_initializer(*it.lock());
  }

  // Then, add all surface coupling dependencies, if any
  if (m_surface_coupling) {
    dag.add_surface_coupling(m_surface_coupling->get_import_fids(),
                             m_surface_coupling->get_export_fids());
  }

  // Finally, we can proceed with checking the dag
  auto& deb_pl = m_atm_params.sublist("Debug");
  int verb_lvl = -1;
  if (deb_pl.isParameter("Atmosphere DAG Verbosity Level")) {
    verb_lvl = deb_pl.get<int>("Atmosphere DAG Verbosity Level");
  }

  if (dag.has_unmet_dependencies()) {
    const int err_verb_lev = verb_lvl >=0 ? verb_lvl : int(AtmProcDAG::VERB_MAX);
    dag.write_dag("error_atm_dag.dot",err_verb_lev);
    EKAT_ERROR_MSG("Error! There are unmet dependencies in the atmosphere internal dag.\n"
                   "       Use the graphviz package to inspect the dependency graph:\n"
                   "    \n"
                   "    $ dot -Tjpg -o error_atm_dag.jpg error_atm_dag.dot\n"
                   "    $ eog error_atm_dag.dot\n");
  }

  // If requested, write a dot file for visualization
  dag.write_dag("scream_atm_dag.dot",std::max(verb_lvl,0));
}

#ifdef SCREAM_DEBUG
void AtmosphereDriver::create_bkp_field_repo () {
  m_bkp_field_repo.registration_begins();
  for (const auto& it : *m_field_repo) {
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
      m_bkp_field_repo.register_field(id,grps);
    }
  }
  m_bkp_field_repo.registration_ends();

  // Deep copy the fields
  for (const auto& it : *m_field_repo) {
    for (const auto& id_field : it.second) {
      const auto& id = id_field.first;
      const auto& f  = id_field.second;
      auto src = f.get_view();
      auto tgt = m_bkp_field_repo.get_field(id).get_view();

      Kokkos::deep_copy(tgt,src);
    }
  }
}
#endif

}  // namespace control
}  // namespace scream
