#include "control/atmosphere_driver.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
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
 *     get the information needed to complete the setup of the FieldRequest and GroupRequest for
 *     the required/computed fields/groups. Their requests MUST be completed upon return from
 *     the 'set_grids' method.
 *     Note: at this stage, atm procs that act on non-ref grid(s) should be able to create their
 *           remappers. The AD will *not* take care of remapping inputs/outputs of the process.
 *  4) Register all fields and all groups from all atm procs inside the field managers, and proceed
 *     to allocate fields. Each field manager (there is one FM per grid) will take care of
 *     accommodating all requests for packing as well as (if possible) bundling of groups.
 *     For more details, see the documentation in the share/field/field_request.hpp header.
 *  5) Set all the fields into the atm procs. Before this point, all the atm procs had were the
 *     FieldIdentifiers for their input/output fields and FieldGroupInfo for their input/output
 *     field groups. Now, we pass actual Field and FieldGroup objects to them, where both the
 *     data (Kokkos::View) and metadata (FieldHeader) inside will be shared across all copies
 *     of the field. This allow data and metadata to be always in sync.
 *     Note: output fields/groups are passed to an atm proc as read-write (i.e., non-const data type),
 *           while input ones are passed as read-only (i.e., const data type). Yes, the atm proc
 *           could cheat, and cast away the const, but we can't prevent that.
 *  6) All the atm inputs (that the AD can deduce by asking the atm proc group for the required fiedls)
 *     are initialized. For restart runs, all fields are read from a netcdf file (to allow BFB
 *     restarts), while for initial runs we offer a few more options (e.g., init a field to
 *     a constant, or as a copy of another field). During this process, we also set the initial
 *     time stamp on all the atm input fields.
 *     If an atm input is not found in the IC file, we'll error out, saving a DAG of the
 *     atm processes, which the user can inspect (to see what's missing in the IC file).
 *  7) All the atm process are initialized. During this call, atm process are able to set up
 *     all the internal structures that they were not able to init previously. For instance,
 *     they can set up remappers from the reference grid to the grid they operate on. They can
 *     also utilize their input fields to perform initialization of some internal data structure.
 *
 * For more info see header comments in the proper files:
 *  - for field                -> src/share/field/field.hpp
 *  - for field manager        -> src/share/field/field_manager.hpp
 *  - for field groups         -> src/share/field/field_group.hpp
 *  - for field/group requests -> src/share/field/field_request.hpp
 *  - for grid                 -> src/share/grid/abstract_grid.hpp
 *  - for grid manager         -> src/share/grid/grids_manager.hpp
 *  - for atm proc             -> src/share/atm_process/atmosphere_process.hpp
 *  - for atm proc group       -> src/share/atm_process/atmosphere_process_group.hpp
 *  - for scorpio input/output -> src/share/io/scorpio_[input|output].hpp
 *  - for output manager       -> src/share/io/scream_output_manager.hpp
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
  check_ad_status (s_params_set, false);

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
  m_grids_manager->build_grids(m_atm_process_group->get_required_grids());

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
  for (const auto& req : m_atm_process_group->get_required_field_requests()) {
    m_field_mgrs.at(req.fid.get_grid_name())->register_field(req);
  }
  for (const auto& req : m_atm_process_group->get_computed_field_requests()) {
    m_field_mgrs.at(req.fid.get_grid_name())->register_field(req);
  }

  // Register required/updated groups
  // IMPORTANT: Some group requests might be formulated in terms of another
  //   group request (see field_requests.hpp). The FieldManager class is able
  //   to handle most of these. However, a GroupRequest that requests to Import
  //   a group from another grid cannot be handled internally by the FieldManager.
  //   In particular, the FM needs the fields in the imported group to be
  //   registered in the FM before it can handle the request.
  //   For this reason, we scan required/updated group requests from all atm procs,
  //   and if we find an Import request, for each field in the "source" group,
  //   we register a corresponding version of the field on the "target" FM.

  // Helper lambda to reduce code duplication
  auto process_imported_groups = [&](const std::set<GroupRequest>& group_requests) {
    for (auto req : group_requests) {
      if (req.derived_type==DerivationType::Import) {
        EKAT_REQUIRE_MSG (req.grid!=req.src_grid,
            "Error! A group request with 'Import' derivation type is meant to import\n"
            "       a group of fields from a grid to another. Found same grid name instead.\n"
            "   group name: " + req.name + "\n"
            "   group to import: " + req.src_name + "\n"
            "   grid name: " + req.grid + "\n");
        // Given request for group A on grid g1 to be an Import of
        // group B on grid g2, register each field in group B in the 
        // field manager on grid g1.
        const auto& fm = m_field_mgrs.at(req.grid);
        const auto& rel_fm   = m_field_mgrs.at(req.src_grid);
        const auto& rel_info = rel_fm->get_groups_info().at(req.src_name);

        auto r = m_grids_manager->create_remapper(req.src_grid,req.grid);
        // Loop over all fields in group src_name on grid src_grid.
        for (const auto& fname : rel_info->m_fields_names) {
          // Get field on src_grid
          auto f = rel_fm->get_field_ptr(fname);

          // Build a FieldRequest for the same field on greq's grid
          auto fid = r->create_tgt_fid(f->get_header().get_identifier());
          FieldRequest freq(fid,req.src_name,req.pack_size);
          fm->register_field(freq);
        }

        // Now that the fields have been imported on this grid's GM,
        // this group request is a "normal" group request, and we can
        // reset the derivation type to "None"
        req.derived_type=DerivationType::None;
      }

      m_field_mgrs.at(req.grid)->register_group(req);
    }
  };

  process_imported_groups (m_atm_process_group->get_required_group_requests());
  process_imported_groups (m_atm_process_group->get_computed_group_requests());

  // Close the FM's, allocate all fields
  for (auto it : m_grids_manager->get_repo()) {
    auto grid = it.second;
    m_field_mgrs[grid->name()]->registration_ends();
  }

  // Set all the fields/groups in the processes. Input fields/groups will be handed
  // to the processes with const scalar type (const Real), to prevent them from
  // overwriting them (though, they can always cast away const...).
  // IMPORTANT: set all computed fields/groups first, since the AtmProcGroup class
  // needs to inspect those before deciding whether a required group is indeed
  // required or not. E.g., in AtmProcGroup [A, B], if A computes group "blah" (or all
  // the fields contained in group "blah"), then group "blah" is not a required
  // group for the AtmProcGroup, even if it is a required group for B.
  for (const auto& req : m_atm_process_group->get_computed_field_requests()) {
    const auto& fid = req.fid;
    auto fm = get_field_mgr(fid.get_grid_name());
    m_atm_process_group->set_computed_field(fm->get_field(fid));
  }
  for (const auto& it : m_atm_process_group->get_computed_group_requests()) {
    auto fm = get_field_mgr(it.grid);
    auto group = fm->get_field_group(it.name);
    m_atm_process_group->set_computed_group(group);
  }
  for (const auto& it : m_atm_process_group->get_required_group_requests()) {
    auto fm = get_field_mgr(it.grid);
    auto group = fm->get_const_field_group(it.name);
    m_atm_process_group->set_required_group(group);
  }
  for (const auto& req : m_atm_process_group->get_required_field_requests()) {
    const auto& fid = req.fid;
    auto fm = get_field_mgr(fid.get_grid_name());
    m_atm_process_group->set_required_field(fm->get_field(fid).get_const());
  }

  m_ad_status |= s_fields_created;
}

void AtmosphereDriver::initialize_output_managers (const bool restarted_run) {
  check_ad_status (s_comm_set | s_params_set | s_grids_created | s_fields_created);

  auto& io_params = m_atm_params.sublist("Scorpio");

  // Build one manager per output yaml file
  using vos_t = std::vector<std::string>;
  const auto& output_yaml_files = io_params.get<vos_t>("Output YAML Files",vos_t{});
  for (const auto& fname : output_yaml_files) {
    ekat::ParameterList params;
    ekat::parse_yaml_file(fname,params);
    m_output_managers.emplace_back();
    auto& om = m_output_managers.back();
    om.setup(m_atm_comm,params,m_field_mgrs,m_grids_manager,m_current_ts,false,restarted_run);
  }

  // Check for model restart output
  if (io_params.isSublist("Model Restart")) {
    auto restart_pl = io_params.sublist("Model Restart");
    // Signal that this is not a normal output, but the model restart one
    m_output_managers.emplace_back();
    auto& om = m_output_managers.back();
    om.setup(m_atm_comm,restart_pl,m_field_mgrs,m_grids_manager,m_current_ts,true,restarted_run);
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
  // TODO: would be nice to do the IC input first, and mark the fields in the
  //       DAG node "Begin of atm time step" in red if there's no initialization
  //       mechanism set for them. That is, allow field XYZ to not be found in
  //       the IC file, and throw an error when the dag is created.

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

  auto& ic_pl = m_atm_params.sublist("Initial Conditions");

  // When processing groups and fields separately, we might end up processing the same
  // field twice. E.g., say we have the required field group G=(f1,f2). If f1 is also 
  // listed as a required field on itself, we might process it twice. To prevent that,
  // we store the processed fields
  // std::set<FieldIdentifier> inited_fields;

  // Check which fields need to have an initial condition.
  std::map<std::string,std::vector<std::string>> ic_fields_names;
  std::vector<FieldIdentifier> ic_fields_to_copy;
  std::map<std::string,std::vector<std::string>> fields_inited;

  // Helper lambda, to reduce code duplication
  auto process_ic_field = [&](const Field<const Real>& f) {
    const auto& fid = f.get_header().get_identifier();
    const auto& grid_name = fid.get_grid_name();
    const auto& fname = fid.name();

    const auto& ic_pl_grid = ic_pl.sublist(grid_name);

    // First, check if the input file contains constant values for some of the fields
    if (ic_pl_grid.isParameter(fname)) {
      // The user provided a constant value for this field. Simply use that.
      if (ic_pl_grid.isType<double>(fname) or ic_pl_grid.isType<std::vector<double>>(fname)) {
        initialize_constant_field(fid, ic_pl_grid);
        fields_inited[grid_name].push_back(fname);

        // Note: f is const, so we can't modify the tracking. So get the same field from the fm
        auto f_nonconst = m_field_mgrs.at(grid_name)->get_field(fid.name());
        f_nonconst.get_header().get_tracking().update_time_stamp(t0);
      } else if (ic_pl_grid.isType<std::string>(fname)) {
        ic_fields_to_copy.push_back(fid);
        fields_inited[grid_name].push_back(fname);
      } else {
        EKAT_REQUIRE_MSG (false, "ERROR: invalid assignment for variable " + fname + ", only scalar double or string, or vector double arguments are allowed");
      }
    } else {
      // If this field is the parent of other subfields, we only read from file the subfields.
      auto c = f.get_header().get_children();
      if (c.size()==0) {
        auto& this_grid_ic_fnames = ic_fields_names[fid.get_grid_name()];
        if (not ekat::contains(this_grid_ic_fnames,fname)) {
          this_grid_ic_fnames.push_back(fname);
        }
      }
    }
  };

  // First the individual input fields...
  for (const auto& f : m_atm_process_group->get_fields_in()) {
    process_ic_field (f);
  }
  // ...then the input groups
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    if (g.m_bundle) {
      process_ic_field(*g.m_bundle);
    }
    for (auto it : g.m_fields) {
      process_ic_field(*it.second);
    }
  }

  // Some fields might be the subfield of a group's bundled field. In that case,
  // we only need to init one: either the bundled field, or all the individual subfields.
  // So loop over the fields that appear to require loading from file, and remove
  // them from the list if they are the subfield of a bundled field already inited
  // (perhaps via initialize_constant_field, or copied from another field).
  for (auto& it1 : ic_fields_names) {
    const auto& grid_name =  it1.first;
    auto fm = m_field_mgrs.at(grid_name);

    // Note: every time we erase an entry in the vector, all iterators are
    //       invalidated, so we need to re-start the for loop.
    bool run_again = true;
    while (run_again) {
      run_again = false;
      auto& names = it1.second;
      for (auto it2=names.begin(); it2!=names.end(); ++it2) {
        const auto& fname = *it2;
        auto f = fm->get_field(fname);
        auto p = f.get_header().get_parent().lock();
        if (p) {
          const auto& pname = p->get_identifier().name();
          if (ekat::contains(fields_inited[grid_name],pname)) {
            // The parent is already inited. No need to init this field as well.
            names.erase(it2);
            run_again = true;
            break;
          }
        }

      }
    }
  }

  // Now loop over all grids, and load from file the needed fields on each grid (if any).
  for (const auto& it : m_field_mgrs) {
    const auto& grid_name = it.first;
    const auto& ic_pl_grid = ic_pl.sublist(grid_name);

    // Check whether we need to load latitude/longitude of reference grid dofs.
    // This option allows the user to set lat or lon in their own
    // test or run setup code rather than by file.
    bool load_latitude  = false;
    bool load_longitude = false;
    if (ic_pl_grid.isParameter("Load Latitude")) {
      load_latitude = ic_pl_grid.get<bool>("Load Latitude");
    }
    if (ic_pl_grid.isParameter("Load Longitude")) {
      load_longitude = ic_pl_grid.get<bool>("Load Longitude");
    }

    if (ic_fields_names[grid_name].size()>0 || load_longitude || load_latitude) {
      MPI_Fint fcomm = MPI_Comm_c2f(m_atm_comm.mpi_comm());
      if (!scorpio::is_eam_pio_subsystem_inited()) {
        scorpio::eam_init_pio_subsystem(fcomm);
      } else {
        EKAT_REQUIRE_MSG (fcomm==scorpio::eam_pio_subsystem_comm(),
            "Error! EAM subsystem was inited with a comm different from the current atm comm.\n");
      }

      // Case where there are fields to load from initial condition file.
      if (ic_fields_names[grid_name].size()>0) {
        // There are fields to read from the nc file. We must have a valid nc file then.
        ekat::ParameterList ic_reader_params;
        ic_reader_params.set("Fields",ic_fields_names[grid_name]);
        ic_reader_params.set("Filename",ic_pl_grid.get<std::string>("Filename"));

        const auto& field_mgr = it.second;
        AtmosphereInput ic_reader(m_atm_comm,ic_reader_params,field_mgr);
        ic_reader.read_variables();
        ic_reader.finalize();

        for (const auto& fname : ic_fields_names[grid_name]) {
          // Set the initial time stamp
          auto f = m_field_mgrs.at(grid_name)->get_field(fname);
          f.get_header().get_tracking().update_time_stamp(t0);
        }
      }

      // Case where lat and/or lon are pulled from initial condition file
      if ( load_latitude || load_longitude) {
        using namespace ShortFieldTagsNames;
        using view_d = AbstractGrid::geo_view_type;
        using view_h = view_d::HostMirror;
        auto grid = m_grids_manager->get_grid(grid_name);
        int ncol  = grid->get_num_local_dofs();
        FieldLayout layout ({COL},{ncol});

        std::vector<std::string> fnames;
        std::map<std::string,FieldLayout> layouts;
        std::map<std::string,view_h> host_views;
        std::map<std::string,view_d> dev_views;
        if (load_latitude) {
          dev_views["lat"] = grid->get_geometry_data("lat");
          host_views["lat"] = Kokkos::create_mirror_view(dev_views["lat"]);
          layouts.emplace("lat",layout);
          fnames.push_back("lat");
        }
        if (load_longitude) {
          dev_views["lon"] = grid->get_geometry_data("lon");
          host_views["lon"] = Kokkos::create_mirror_view(dev_views["lon"]);
          layouts.emplace("lon",layout);
          fnames.push_back("lon");
        }

        ekat::ParameterList lat_lon_params;
        lat_lon_params.set("Fields",fnames);
        lat_lon_params.set("Filename",ic_pl_grid.get<std::string>("Filename"));

        AtmosphereInput lat_lon_reader(m_atm_comm,lat_lon_params,grid,host_views,layouts);
        lat_lon_reader.read_variables();
        lat_lon_reader.finalize();
        for (auto& fname : fnames) {
          Kokkos::deep_copy(dev_views[fname],host_views[fname]);
        }
      }
    }
  }

  // If there were any fields that needed to be copied per the input yaml file, now we copy them.
  for (const auto& tgt_fid : ic_fields_to_copy) {
    const auto& tgt_fname = tgt_fid.name();
    const auto& tgt_gname = tgt_fid.get_grid_name();

    auto tgt_fm = get_field_mgr(tgt_gname);

    // The user can request to init a field from its copy on another grid
    std::string src_fname, src_gname;

    const auto& param_value = ic_pl.sublist(tgt_gname).get<std::string>(tgt_fname);
    const auto tokens = ekat::split(param_value,',');
    EKAT_REQUIRE_MSG (tokens.size()==1 || tokens.size()==2,
        "Error! To copy an initial condition for a field from another, use one of the following ways:\n"
        "    - field_1: field_2\n"
        "    - field_1: field_2, grid_2\n"
        "The first assumes field_2 is on the same grid as field_1, while the second allows\n"
        "cross-grid imports.\n\n"
        "Parameter list entry:\n"
        "   " + tgt_fname + ": " + param_value + "\n");

    src_fname = ekat::trim(tokens[0]);
    if (tokens.size()==2) {
      src_gname = ekat::trim(tokens[1]);
    } else {
      src_gname = tgt_gname;
    }

    auto src_fm = get_field_mgr(src_gname);

    // The field must exist in the fm on the input field's grid
    EKAT_REQUIRE_MSG (src_fm->has_field(src_fname),
        "Error! Source field for initial condition not found in the field manager.\n"
        "       Grid name:     " + tgt_gname + "\n"
        "       Field to init: " + tgt_fname + "\n"
        "       Source field:  " + src_fname + " (NOT FOUND)\n");

    // Get the two fields, and copy src to tgt
    auto f_tgt = tgt_fm->get_field(tgt_fname);
    auto f_src = src_fm->get_field(src_fname);
    if (src_gname==tgt_gname) {
      // Same grid: simply copy the field
      f_tgt.deep_copy(f_src);
    } else {
      // Different grid: create a remapper on the fly
      auto r = m_grids_manager->create_remapper(src_gname,tgt_gname);
      r->registration_begins();
      r->register_field(f_src,f_tgt);
      r->registration_ends();
      r->remap(true);
    }
    // Set the initial time stamp
    f_tgt.get_header().get_tracking().update_time_stamp(t0);
  }

  // Final step: it is possible to have a bundled group G1=(f1,f2,f3),
  // where the IC are read from file for f1, f2, and f3. In that case,
  // the time stamp for the bundled G1 has not be inited, but the data
  // is valid (all entries have been inited). Let's fix that.
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    if (g.m_bundle) {
      auto& track = g.m_bundle->get_header().get_tracking();
      if (not track.get_time_stamp().is_valid()) {
        // The bundled field has not been inited. Check if all the subfields
        // have been inited. If so, init the timestamp of the bundled field too.
        const auto& children = track.get_children();
        bool all_good = children.size()>0; // If no children, then something is off, so mark as not good
        for (auto wp : children) {
          auto sp = wp.lock();
          if (not sp->get_time_stamp().is_valid()) {
            all_good = false;
            break;
          }
        }
        if (all_good) {
          track.update_time_stamp(t0);
        }
      }
    }
  }

  m_current_ts = t0;

  m_ad_status |= s_fields_inited;
}

void AtmosphereDriver::
initialize_constant_field(const FieldIdentifier& fid,
                          const ekat::ParameterList& ic_pl)
{
  const auto& name = fid.name();
  const auto& grid = fid.get_grid_name();
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
  // Initialize memory buffer for all atm processes
  m_memory_buffer = std::make_shared<ATMBufferManager>();
  m_atm_process_group->initialize_atm_memory_buffer(*m_memory_buffer);

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

void AtmosphereDriver::run (const int dt) {
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
    out_mgr.run(m_current_ts);
  }

  if (m_surface_coupling) {
    // Export fluxes from the component coupler (if any)
    m_surface_coupling->do_export();
  }
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {

  // Finalize and destroy output streams, make sure files are closed
  for (auto& out_mgr : m_output_managers) {
    out_mgr.finalize();
  }
  m_output_managers.clear();

  // Finalize, and then destroy all atmosphere processes
  m_atm_process_group->finalize( /* inputs ? */ );
  m_atm_process_group = nullptr;

  // Destroy the buffer manager
  m_memory_buffer = nullptr;

  // Destroy the surface coupling (if any)
  m_surface_coupling = nullptr;

  // Destroy the grids manager
  m_grids_manager = nullptr;

  // Destroy all the fields manager
  for (auto it : m_field_mgrs) {
    it.second->clean_up();
  }

  // Finalize scorpio
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
