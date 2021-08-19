#include "control/atmosphere_driver.hpp"
#include <memory>

#include "ekat/ekat_parameter_list.hpp"
#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
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
 *     Note: at this stage, atm procs that act on non-ref grid(s) should be able to create their
 *           remappers. The AD will *not* take care of remapping inputs/outputs of the process.
 *  4) Register all fields from all atm procs inside the field manager
 *  5) Set all the fields into the atm procs. Before this point, all the atm procs had were the
 *     FieldIdentifiers for their input/output fields. Now, we pass actual Field objects to them,
 *     where both the data (Kokkos::View) and metadata (FieldHeader) inside will be shared across
 *     all processes using the field. This allow data and metadata to be always in sync.
 *     Note: output fields are passed to an atm proc as read-write (i.e., non-const data type),
 *           while input fields are passed as read-only (i.e., const data type). Yes, the atm proc
 *           could cheat, and cast away the const, but we can't prevent that. However, in debug builds,
 *           we store 2 copies of each field, and use the extra copy to check, at run time, that
 *           no process alters the values of any of its input fields.
 *  6) All the atm inputs (that the AD can deduce by asking the atm proc group for the required fiedls)
 *     are initialized, by reading values from an initial conditions netcdf file.
 *     If an atm input is not found in the IC file, we'll error out, saving a DAG of the
 *     atm processes, which the user can inspect (to see what's missing in the IC file).
 *  7) All the atm process are initialized. During this call, atm process are able to set up
 *     all the internal structures that they were not able to init previously. They can also
 *     utilize their input fields to perform initialization of some internal data structure.
 *  8) Finally, set the initial time stamp on all fields, and perform some debug structure setup.
 *
 */

AtmosphereDriver::
AtmosphereDriver(const ekat::Comm& atm_comm,
                 const ekat::ParameterList& params)
{
  set_comm(atm_comm);
  set_params(params);
}

void AtmosphereDriver::
set_comm(const ekat::Comm& atm_comm)
{
  // I can't think of a scenario where changing the comm is a good idea,
  // so let's forbid it, for now.
  check_ad_status (s_comm_set, false);

  m_atm_comm = atm_comm;

  m_ad_status |= s_comm_set;
}

void AtmosphereDriver::
set_params(const ekat::ParameterList& atm_params)
{
  // I can't think of a scenario where changing the params is useful,
  // so let's forbid it, for now.
  check_ad_status (~s_params_set);

  m_atm_params = atm_params;

  m_ad_status |= s_params_set;
}

void AtmosphereDriver::create_atm_processes()
{
  // At this point, must have comm and params set.
  check_ad_status(s_comm_set | s_params_set);

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_comm,m_atm_params.sublist("Atmosphere Processes"));

  m_ad_status |= s_procs_created;
}

void AtmosphereDriver::create_grids()
{
  // Must have procs created by now (and comm/params set)
  check_ad_status (s_procs_created | s_comm_set | s_params_set);

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

  m_ad_status |= s_grids_created;
}

void AtmosphereDriver::create_fields()
{
  // Must have grids and procs at this point
  check_ad_status (s_procs_created | s_grids_created);

  // By now, the processes should have fully built the ids of their
  // required/computed fields and groups. Let them register them in the FM
  for (auto it : m_grids_manager->get_repo()) {
    auto grid = it.second;
    m_field_mgrs[grid->name()] = std::make_shared<field_mgr_type>(grid);
    m_field_mgrs[grid->name()]->registration_begins();
  }

  // Register required/computed fields
  for (const auto& req : m_atm_process_group->get_required_fields()) {
    m_field_mgrs.at(req.fid.get_grid_name())->register_field(req);
  }
  for (const auto& req : m_atm_process_group->get_computed_fields()) {
    m_field_mgrs.at(req.fid.get_grid_name())->register_field(req);
  }

  // Register required/updated groups
  // IMPORTANT: Some FMs on non-ref grids *might* depend on the ref-grid FM.
  //            E.g., there could be a GroupRequest gr1 that requires to be an
  //            alias of gr2, with gr2 defined on ref grid. Most obvious case:
  //            dyn needs a copy of group "tracers" on dyn grid. But nobody
  //            register tracers on dyn grid, so at first sight, it would appear
  //            that the tracers group on dyn grid is empty.
  //            To overcome this issue, dyn can register the tracers group on the
  //            dyn grid as an "alias" of the tracers group on the ref grid,
  //            and the AD must take care of ensuring that the two contain the
  //            same fields (effectively "adding" fields on the dyn grid).
  //            To do this, we loop over the requests, and if we find a request
  //            gr1 that has a 'relative' request gr2 (see field_request.hpp for a
  //            definition of 'relative'), and if gr2 is on a different grid, we
  //            make sure that:
  //              a) if the relative group is a 'child' or 'alias' (see field_request.hpp
  //                 for an explanation of these terms), we ensure the relative
  //                 group is also registered in this grid, with fields in the same order.
  //              b) all the fields of the relative request are registered
  //            in the FM on gr1's grid.

  // Helper lambda to reduce code duplication
  auto import_relative_group = [&](const GroupRequest& greq) {
    // Given request on grid g1, with relative (alias or child) on grid g2,
    // imports all fields of the relative from grid g2 to grid g1.

    const auto& fm = m_field_mgrs.at(greq.grid);
    const auto& rel = *greq.relative;
    const auto& rel_fm   = m_field_mgrs.at(rel.grid);
    const auto& rel_info = rel_fm->get_groups_info().at(rel.name);

    auto r = m_grids_manager->create_remapper(rel.grid,greq.grid);
    // Loop over all fields in group rel_name on grid rel_grid.
    for (const auto& fname : rel_info->m_fields_names) {
      // Get field on rel_grid
      auto f = rel_fm->get_field_ptr(fname);

      // Build a FieldRequest for the same field on greq's grid
      auto fid = r->create_tgt_fid(f->get_header().get_identifier());
      FieldRequest freq(fid,rel.name,rel.pack_size);
      fm->register_field(freq);
    }

    // Register also the relative group on this grid
    GroupRequest rel_on_my_grid(rel.name,greq.grid,greq.pack_size,greq.bundling);
    fm->register_group(rel_on_my_grid);
  };

  for (const auto& req : m_atm_process_group->get_required_groups()) {
    if (req.relative!=nullptr && req.grid!=req.relative->grid) {
      EKAT_REQUIRE_MSG (req.relative_type==Relationship::Alias || req.relative_type==Relationship::Child,
          "Error! Linking group requests on different grids is only allowed for 'Alias' or 'Child' relatives.\n"
          "       See field_request.hpp for more info on what that means.\n");
      import_relative_group(req);
    }
    m_field_mgrs.at(req.grid)->register_group(req);
  }
  for (const auto& req : m_atm_process_group->get_updated_groups()) {
    if (req.relative!=nullptr && req.grid!=req.relative->grid) {
      EKAT_REQUIRE_MSG (req.relative_type==Relationship::Alias || req.relative_type==Relationship::Child,
          "Error! Linking group requests on different grids is only allowed for 'Alias' or 'Child' relatives.\n"
          "       See field_request.hpp for more info on what that means.\n");
      import_relative_group(req);
    }
    m_field_mgrs.at(req.grid)->register_group(req);
  }

  // Close the FM's, allocate all fields
  for (auto it : m_grids_manager->get_repo()) {
    auto grid = it.second;
    m_field_mgrs[grid->name()]->registration_ends();
  }

  m_ad_status |= s_fields_created;
}

void AtmosphereDriver::initialize_output_managers (const bool restarted_run) {
  check_ad_status (s_comm_set | s_params_set | s_grids_created | s_fields_created);

  // For each grid in the grids manager, grab a homonymous sublist of the
  // "Output Managers" atm params sublist (might be empty), and create
  // an output manager for that grid.
  auto& io_params = m_atm_params.sublist("SCORPIO");
  for (auto it : m_grids_manager->get_repo()) {
    auto fm = m_field_mgrs.at(it.first);
    m_output_managers[it.first].setup(m_atm_comm,io_params,fm,restarted_run);
  }

  m_ad_status |= s_output_inited;
}

void AtmosphereDriver::
initialize_fields (const util::TimeStamp& t0)
{
  // See if we need to print a DAG. We do this first, cause if any input
  // field is missing from the initial condition file, an error will be thrown.
  // By printing the DAG first, we give the user the possibility of seeing
  // what fields are inputs to the atm time step, so he/she can fix the i.c. file.

  auto& deb_pl = m_atm_params.sublist("Debug");
  const int verb_lvl = deb_pl.get<int>("Atmosphere DAG Verbosity Level",-1);
  if (verb_lvl>0) {
    // Check the atm DAG for missing stuff
    AtmProcDAG dag;

    // First, add all atm processes
    dag.create_dag(*m_atm_process_group,m_field_mgrs);

    // Then, add all surface coupling dependencies, if any
    if (m_surface_coupling) {
      dag.add_surface_coupling(m_surface_coupling->get_import_fids(),
                               m_surface_coupling->get_export_fids());
    }

    // Write a dot file for visualization
    dag.write_dag("scream_atm_dag.dot",std::max(verb_lvl,0));
  }

  // Figure out the list of inputs for the atmosphere.
  const auto& fields_in = m_atm_process_group->get_required_fields();

  const auto& ic_pl = m_atm_params.sublist("Initial Conditions");
  const auto& ref_grid_name = m_grids_manager->get_reference_grid()->name();

  // Create parameter list for AtmosphereInput
  ekat::ParameterList ic_reader_params;
  ic_reader_params.set("GRID",ref_grid_name);
  std::vector<std::string> ic_fields_names;
  int ifield=0;
  std::vector<FieldIdentifier> ic_fields_to_copy;
  for (const auto& req : fields_in) {
    const auto& fid = req.fid;
    const auto& name = fid.name();
    const auto& fm = get_field_mgr(fid.get_grid_name());

    auto f = fm->get_field(fid);
    // First, check if the input file contains constant values for some of the fields
    if (ic_pl.isParameter(name)) {
      // The user provided a constant value for this field. Simply use that.
      if (ic_pl.isType<double>(name) or ic_pl.isType<std::vector<double>>(name)) {
        initialize_constant_field(req, ic_pl);
      } else if (ic_pl.isType<std::string>(name)) {
        ic_fields_to_copy.push_back(req.fid);
      } else {
        EKAT_REQUIRE_MSG (false, "ERROR: invalid assignment for variable " + name + ", only scalar double or string, or vector double arguments are allowed");
      }
    } else {
      // The field does not have a constant value, so we expect to find it in the nc file
      // A requirement for this, is that the field is on the ref grid
      EKAT_REQUIRE_MSG (req.fid.get_grid_name()==ref_grid_name,
          "Error! So far, only reference grid fields can be inited via Scorpio input.\n"
          "       Ref grid name:    " + ref_grid_name + "\n"
          "       Input field grid: " + req.fid.get_grid_name() + "\n");

      ic_fields_names.push_back(name);
      ++ifield;

    }
    f.get_header().get_tracking().update_time_stamp(t0);
  }
  ic_reader_params.set("FIELDS",ic_fields_names);

  // Check whether we need to load latitude/longitude of reference grid dofs.
  // This option allows the user to set lat or lon in their own
  // test or run setup code rather than by file.
  bool load_latitude  = false;
  bool load_longitude = false;
  if (ic_pl.isParameter("I")) {
    load_latitude = ic_pl.get<bool>("Load Latitude");
  }
  if (ic_pl.isParameter("load_longitude")) {
    load_longitude = ic_pl.get<bool>("Load Longitude");
  }

  if (ifield>0 || load_longitude || load_latitude) {
    // There are fields to read from the nc file. We must have a valid nc file then.
    ic_reader_params.set("FILENAME",ic_pl.get<std::string>("Initial Conditions File"));

    MPI_Fint fcomm = MPI_Comm_c2f(m_atm_comm.mpi_comm());
    if (!scorpio::is_eam_pio_subsystem_inited()) {
      scorpio::eam_init_pio_subsystem(fcomm);
    } else {
      EKAT_REQUIRE_MSG (fcomm==scorpio::eam_pio_subsystem_comm(),
          "Error! EAM subsystem was inited with a comm different from the current atm comm.\n");
    }

    // Case where there are fields to load from initial condition file.
    if (ifield>0) {
      AtmosphereInput ic_reader(m_atm_comm,ic_reader_params,get_ref_grid_field_mgr());
      ic_reader.read_variables();
      ic_reader.finalize();
    }

    // Case where lat and/or lon pulled from initial condition file
    if ( load_latitude || load_longitude) {
      using namespace ShortFieldTagsNames;
      using view_d = AbstractGrid::geo_view_type;
      using view_h = view_d::HostMirror;
      auto ref_grid = m_grids_manager->get_reference_grid();
      int ncol  = ref_grid->get_num_local_dofs();
      FieldLayout layout ({COL},{ncol});

      std::vector<std::string> fnames;
      std::map<std::string,FieldLayout> layouts;
      std::map<std::string,view_h> host_views;
      std::map<std::string,view_d> dev_views;
      if (load_latitude) {
        dev_views["lat"] = ref_grid->get_geometry_data("lat");
        host_views["lat"] = Kokkos::create_mirror_view(dev_views["lat"]);
        layouts.emplace("lat",layout);
        fnames.push_back("lat");
      }
      if (load_longitude) {
        dev_views["lon"] = ref_grid->get_geometry_data("lon");
        host_views["lon"] = Kokkos::create_mirror_view(dev_views["lon"]);
        layouts.emplace("lon",layout);
        fnames.push_back("lon");
      }

      ekat::ParameterList lat_lon_params;
      lat_lon_params.set("FIELDS",fnames);
      lat_lon_params.set("GRID",ref_grid->name());
      lat_lon_params.set("FILENAME",ic_pl.get<std::string>("Initial Conditions File"));

      AtmosphereInput lat_lon_reader(m_atm_comm,lat_lon_params);
      lat_lon_reader.init(ref_grid,host_views,layouts);
      lat_lon_reader.read_variables();
      lat_lon_reader.finalize();
    }
  }

  // If there were any fields that needed to be copied per the input yaml file, now we copy them.
  for (const auto& tgt_fid : ic_fields_to_copy) {
    auto fm = get_field_mgr(tgt_fid.get_grid_name());

    const auto& src_name = ic_pl.get<std::string>(tgt_fid.name());
    // The target field must exist in the fm on the input field's grid
    EKAT_REQUIRE_MSG (fm->has_field(src_name),
        "Error! Source field for initial condition not found in the field manager.\n"
        "       Grid name:     " + tgt_fid.get_grid_name() + "\n"
        "       Field to init: " + tgt_fid.name() + "\n"
        "       Source field:  " + src_name + " (NOT FOUND)\n");

    // Get the two fields, and copy src to tgt
    auto f_tgt = fm->get_field(tgt_fid.name());
    auto f_src = fm->get_field(src_name);
    f_tgt.deep_copy(f_src);
  }

  m_current_ts = t0;

  m_ad_status |= s_fields_inited;
}

void AtmosphereDriver::initialize_constant_field(const FieldRequest& freq, const ekat::ParameterList& ic_pl)
{
  const auto& name = freq.fid.name();
  const auto& grid = freq.fid.get_grid_name();
  auto f = get_field_mgr(grid)->get_field(name);
  // The user provided a constant value for this field. Simply use that.
  const auto& layout = f.get_header().get_identifier().get_layout();

  // For vector fields, we expect something like "fname: [val0,...,valN],
  // where the field dim is N+1. For scalars, "fname: val". So check the
  // field layout first, so we know what to get from the parameter list.
  if (layout.is_vector_layout()) {
    const auto idim = layout.get_vector_dim();
    const auto vec_dim = layout.dim(idim);
    const auto& values = ic_pl.get<std::vector<double>>(name);
    EKAT_REQUIRE_MSG (values.size()==static_cast<size_t>(vec_dim),
        "Error! Initial condition values array for '" + name + "' has the wrong dimension.\n"
        "       Field dimension: " + std::to_string(vec_dim) + "\n"
        "       Array dimenions: " + std::to_string(values.size()) + "\n");

    // Extract a subfield for each component. This is not "too" expensive, expecially
    // considering that this code is executed during initialization only.
    for (int comp=0; comp<vec_dim; ++comp) {
      auto f_i = f.get_component(comp);
      f_i.deep_copy(values[comp]);
    }
  } else {
    const auto& value = ic_pl.get<double>(name);
    f.deep_copy(value);
  }
}

void AtmosphereDriver::initialize_atm_procs ()
{
  // Set all the fields in the processes needing them (before, they only had ids)
  // Input fields will be handed to the processes as const
  const auto& inputs  = m_atm_process_group->get_required_fields();
  const auto& outputs = m_atm_process_group->get_computed_fields();
  for (const auto& req : inputs) {
    const auto& fid = req.fid;
    auto fm = get_field_mgr(fid.get_grid_name());
    m_atm_process_group->set_required_field(fm->get_field(fid).get_const());
  }
  // Output fields are handed to the processes as writable
  for (const auto& req : outputs) {
    const auto& fid = req.fid;
    auto fm = get_field_mgr(fid.get_grid_name());
    m_atm_process_group->set_computed_field(fm->get_field(fid));
  }
  // Set all groups of fields
  for (const auto& it : m_atm_process_group->get_required_groups()) {
    auto fm = get_field_mgr(it.grid);
    auto group = fm->get_const_field_group(it.name);
    m_atm_process_group->set_required_group(group);
  }
  for (const auto& it : m_atm_process_group->get_updated_groups()) {
    auto fm = get_field_mgr(it.grid);
    auto group = fm->get_field_group(it.name);
    m_atm_process_group->set_updated_group(group);
  }

  // Initialize memory buffer for all atm processes
  m_atm_process_group->initialize_atm_memory_buffer(m_memory_buffer);

  // Initialize the processes
  m_atm_process_group->initialize(m_current_ts);

  m_ad_status |= s_procs_inited;
}

void AtmosphereDriver::
initialize (const ekat::Comm& atm_comm,
            const ekat::ParameterList& params,
            const util::TimeStamp& t0,
            const bool restarted_run)
{
  set_comm(atm_comm);
  set_params(params);

  create_atm_processes ();

  create_grids ();

  create_fields ();

  initialize_fields (t0);

  initialize_output_managers (restarted_run);

  initialize_atm_procs ();
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
  for (auto& out_mgr : m_output_managers) {
    out_mgr.second.run(m_current_ts);
  }

  if (m_surface_coupling) {
    // Export fluxes from the component coupler (if any)
    m_surface_coupling->do_export();
  }
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  m_atm_process_group->finalize( /* inputs ? */ );

  // Finalize output streams, make sure files are closed
  for (auto& out_mgr : m_output_managers) {
    out_mgr.second.finalize();
  }

  for (auto it : m_field_mgrs) {
    it.second->clean_up();
  }

  if (scorpio::is_eam_pio_subsystem_inited()) {
    scorpio::eam_pio_finalize();
  }
}

AtmosphereDriver::field_mgr_ptr
AtmosphereDriver::get_ref_grid_field_mgr () const {
  EKAT_REQUIRE_MSG (m_ad_status & s_grids_created,
      "Error! Field manager(s) are created *after* the grids.\n");

  auto ref_grid = m_grids_manager->get_reference_grid();
  return get_field_mgr(ref_grid->name());
}

AtmosphereDriver::field_mgr_ptr
AtmosphereDriver::get_field_mgr (const std::string& grid_name) const {
  EKAT_REQUIRE_MSG (m_ad_status & s_grids_created,
      "Error! Field manager(s) are created *after* the grids.\n");
  // map::at would throw, but you won't know which map threw.
  // With our own msg, we can tell you where the throw happened.
  EKAT_REQUIRE_MSG(m_field_mgrs.find(grid_name)!=m_field_mgrs.end(),
      "Error! Request for field manager on a non-existing grid '" + grid_name + "'.\n");

  return m_field_mgrs.at(grid_name);
}

void AtmosphereDriver::
check_ad_status (const int flag, const bool must_be_set)
{
  if (must_be_set) {
    EKAT_REQUIRE_MSG (m_ad_status & flag,
        "Error! Failed AD status check:\n"
        "        expected flag:  " + std::to_string(flag) + "\n"
        "        ad status flag: " + std::to_string(m_ad_status) + "\n");
  } else {
    EKAT_REQUIRE_MSG (~m_ad_status & flag,
        "Error! Failed AD status check:\n"
        "        not expected flag:  " + std::to_string(flag) + "\n"
        "        ad status flag: " + std::to_string(m_ad_status) + "\n");
  }
}

}  // namespace control
}  // namespace scream
