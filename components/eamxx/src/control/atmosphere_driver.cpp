#include "control/atmosphere_driver.hpp"
#include "control/atmosphere_surface_coupling_importer.hpp"
#include "control/atmosphere_surface_coupling_exporter.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/util/scream_timing.hpp"
#include "share/util/scream_utils.hpp"
#include "share/io/scream_io_utils.hpp"
#include "share/property_checks/mass_and_energy_column_conservation_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

// The global variable fvphyshack is used to help the initial pgN implementation
// work around some current AD constraints. Search the code for "fvphyshack" to
// find blocks that eventually should be removed in favor of a design that
// accounts for pg2. Some blocks may turn out to be unnecessary, and I simply
// didn't realize I could do without the workaround.
#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"

#ifndef SCREAM_CIME_BUILD
#include <unistd.h>
#endif

#include <fstream>
#include <random>

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
 *     they can set up remappers from the physics grid to the grid they operate on. They can
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

AtmosphereDriver::~AtmosphereDriver ()
{
  finalize();
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

  create_logger ();

  m_ad_status |= s_params_set;
}

void AtmosphereDriver::
init_scorpio(const int atm_id)
{
  check_ad_status (s_comm_set, true);

  // Init scorpio right away, in case some class (atm procs, grids,...)
  // needs to source some data from NC files during construction,
  // before we start processing IC files.
  EKAT_REQUIRE_MSG (!scorpio::is_subsystem_inited(),
      "Error! The PIO subsystem was alreday inited before the driver was constructed.\n"
      "       This is an unexpected behavior. Please, contact developers.\n");
  scorpio::init_subsystem(m_atm_comm,atm_id);

  // In CIME runs, gptl is already inited. In standalone runs, it might
  // not be, depending on what scorpio does.
  init_gptl(m_gptl_externally_handled);

  m_ad_status |= s_scorpio_inited;
}

void AtmosphereDriver::
init_time_stamps (const util::TimeStamp& run_t0, const util::TimeStamp& case_t0, int run_type)
{
  m_atm_logger->info("  [EAMxx] Run  start time stamp: " + run_t0.to_string());
  m_atm_logger->info("  [EAMxx] Case start time stamp: " + case_t0.to_string());

  EKAT_REQUIRE_MSG (case_t0<=run_t0,
      "Error! Case t0 time stamp must precede the run t0 time stamp.\n"
      "  - case t0: " + case_t0.to_string() + "\n"
      "  - run  t0: " + run_t0.to_string() + "\n");

  // Initialize time stamps
  m_run_t0 = m_current_ts = run_t0;
  m_case_t0 = case_t0;

  switch (run_type) {
    case 0:
      m_run_type = RunType::Initial; break;
    case 1:
      m_run_type = RunType::Restart; break;
    case -1:
      m_run_type = case_t0==run_t0 ? RunType::Initial : RunType::Restart; break;
    default:
      EKAT_ERROR_MSG ("Unsupported/unrecognized run_type: " + std::to_string(run_type) + "\n");
  }
}

void AtmosphereDriver::
setup_iop ()
{
  // At this point, must have comm, params, initialized timestamps created.
  check_ad_status(s_comm_set | s_params_set | s_ts_inited);

  // Check to make sure iop is not already initialized
  EKAT_REQUIRE_MSG(not m_iop, "Error! setup_iop() is called, but IOP already set up.\n");

  // This function should only be called if we are enabling IOP
  const bool enable_iop =
    m_atm_params.sublist("driver_options").get("enable_iop", false);
  EKAT_REQUIRE_MSG(enable_iop, "Error! setup_iop() is called, but enable_iop=false "
                               "in driver_options parameters.\n");

  // Params must include iop_options sublist.
  const auto iop_sublist_exists = m_atm_params.isSublist("iop_options");
  EKAT_REQUIRE_MSG(iop_sublist_exists,
                   "Error! setup_iop() is called, but no iop_options "
                   "defined in parameters.\n");

  const auto iop_params = m_atm_params.sublist("iop_options");
  const auto phys_grid = m_grids_manager->get_grid("Physics");
  const auto nlevs = phys_grid->get_num_vertical_levels();
  const auto hyam = phys_grid->get_geometry_data("hyam");
  const auto hybm = phys_grid->get_geometry_data("hybm");

  m_iop = std::make_shared<IntensiveObservationPeriod>(m_atm_comm,
                                                       iop_params,
                                                       m_run_t0,
                                                       nlevs,
                                                       hyam,
                                                       hybm);

  // Set IOP object in atm processes
  m_atm_process_group->set_iop(m_iop);
}

void AtmosphereDriver::create_atm_processes()
{
  m_atm_logger->info("[EAMxx] create_atm_processes  ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_atm_processes");

  // At this point, must have comm and params set.
  check_ad_status(s_comm_set | s_params_set);

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  auto& atm_proc_params = m_atm_params.sublist("atmosphere_processes");
  atm_proc_params.rename("EAMxx");
  atm_proc_params.set("Logger",m_atm_logger);
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_comm,atm_proc_params);

  m_ad_status |= s_procs_created;
  stop_timer("EAMxx::create_atm_processes");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_atm_processes  ... done!");
}

void AtmosphereDriver::create_grids()
{
  m_atm_logger->info("[EAMxx] create_grids ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_grids");

  // Must have procs created by now (and comm/params set)
  check_ad_status (s_procs_created | s_comm_set | s_params_set | s_ts_inited);

  // Create the grids manager
  auto& gm_params = m_atm_params.sublist("grids_manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");

  // The GridsManager might load some geometric data from IC file.
  // To avoid having to pass the same data twice in the input file,
  // we have the AD add the IC file name to the GM params
  const auto& ic_pl = m_atm_params.sublist("initial_conditions");
  if (m_run_type==RunType::Restart) {
    // Restarted run -> read geo data from restart file
    const auto& provenance = m_atm_params.sublist("provenance");
    const auto& casename = provenance.get<std::string>("rest_caseid");
    auto filename = find_filename_in_rpointer (casename+".scream",true,m_atm_comm,m_run_t0);
    gm_params.set("ic_filename", filename);
    m_atm_params.sublist("provenance").set("initial_conditions_file",filename);
  } else if (ic_pl.isParameter("Filename")) {
    // Initial run, if an IC file is present, pass it.
    auto filename = ic_pl.get<std::string>("Filename");
    gm_params.set("ic_filename", filename);
    m_atm_params.sublist("provenance").set("initial_conditions_file",filename);
  }

  m_atm_logger->debug("  [EAMxx] Creating grid manager '" + gm_type + "' ...");
  m_grids_manager = GridsManagerFactory::instance().create(gm_type,m_atm_comm,gm_params);

  m_atm_logger->debug("  [EAMxx] Creating grid manager '" + gm_type + "' ... done!");

  // Tell the grid manager to build all the grids required
  // by the atm processes
  m_grids_manager->build_grids();

  m_atm_logger->debug("  [EAMxx] Grids created.");

  // If TMS process is enabled, SHOC needs to know to request tms' surface drag coefficient
  // as a required field during the set_grid() call below, but SHOC does not have knowledge
  // of other processes. The driver needs propgate this information to SHOC.
  if(m_atm_process_group->has_process("tms") &&
     m_atm_process_group->has_process("shoc")) {
    setup_shoc_tms_links();
  }

  // IOP object needs the grids_manager to have been created, but is then needed in set_grids()
  // implementation of some processes, so setup here.
  const bool enable_iop =
    m_atm_params.sublist("driver_options").get("enable_iop", false);
  if (enable_iop) {
    setup_iop ();
  }

  // Set the grids in the processes. Do this by passing the grids manager.
  // Each process will grab what they need
  m_atm_process_group->set_grids(m_grids_manager);

  // Also make each atm proc build requests for tendency fields, if needed
  m_atm_process_group->setup_tendencies_requests();

  m_ad_status |= s_grids_created;

  stop_timer("EAMxx::create_grids");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_grids ... done!");
}

void AtmosphereDriver::setup_surface_coupling_data_manager(SurfaceCouplingTransferType transfer_type,
                                                           const int num_cpl_fields, const int num_scream_fields,
                                                           const int field_size, Real* data_ptr,
#ifdef HAVE_MOAB
                                                           Real* data_ptr_moab,
#endif
                                                           char* names_ptr, int* cpl_indices_ptr, int* vec_comps_ptr,
                                                           Real* constant_multiple_ptr, bool* do_transfer_during_init_ptr)
{
  std::shared_ptr<SCDataManager> sc_data_mgr;

  if (transfer_type==SurfaceCouplingTransferType::Import) {

    m_surface_coupling_import_data_manager = std::make_shared<SCDataManager>();
    sc_data_mgr = m_surface_coupling_import_data_manager;

  } else if (transfer_type==SurfaceCouplingTransferType::Export) {

    m_surface_coupling_export_data_manager = std::make_shared<SCDataManager>();
    sc_data_mgr= m_surface_coupling_export_data_manager;

  } else EKAT_ERROR_MSG("Error! Unexpected SurfaceCouplingTransferType.");

  sc_data_mgr->setup_internals(num_cpl_fields, num_scream_fields, field_size, data_ptr,
#ifdef HAVE_MOAB
                               data_ptr_moab,
#endif
                               names_ptr, cpl_indices_ptr, vec_comps_ptr,
                               constant_multiple_ptr, do_transfer_during_init_ptr);
}

void AtmosphereDriver::setup_surface_coupling_processes () const
{
  // Loop through atmosphere processes and look for importer/exporter. If one is
  // found, cast to derived class type and call setup_surface_coupling_data()
  bool importer_found = false;
  bool exporter_found = false;

  for (int proc=0; proc<m_atm_process_group->get_num_processes(); ++proc) {

    const auto atm_proc = m_atm_process_group->get_process_nonconst(proc);
    if (atm_proc->type() == AtmosphereProcessType::SurfaceCouplingImporter) {
      importer_found = true;

      EKAT_REQUIRE_MSG(m_surface_coupling_import_data_manager != nullptr,
                       "Error! SurfaceCouplingImporter atm process found, "
                       "but m_surface_coupling_import_data_manager was not "
                       "setup.\n");

      std::shared_ptr<SurfaceCouplingImporter> importer = std::dynamic_pointer_cast<SurfaceCouplingImporter>(atm_proc);
      importer->setup_surface_coupling_data(*m_surface_coupling_import_data_manager);
    }
    if (atm_proc->type() == AtmosphereProcessType::SurfaceCouplingExporter) {
      exporter_found = true;

      EKAT_REQUIRE_MSG(m_surface_coupling_export_data_manager != nullptr,
                       "Error! SurfaceCouplingExporter atm process found, "
                       "but m_surface_coupling_export_data_manager was not "
                       "setup.\n");

      std::shared_ptr<SurfaceCouplingExporter> exporter = std::dynamic_pointer_cast<SurfaceCouplingExporter>(atm_proc);
      exporter->setup_surface_coupling_data(*m_surface_coupling_export_data_manager);
    }
  }

  // If import or export data manager is defined,
  // ensure corresponding atm process was found.
  if (m_surface_coupling_import_data_manager) {
    EKAT_REQUIRE_MSG(importer_found, "Error! SurfaceCoupling importer data was setup, but no atm process "
                                     "of type AtmosphereProcessType::SurfaceCouplingImporter exists.\n");
  }
  if (m_surface_coupling_export_data_manager) {
    EKAT_REQUIRE_MSG(exporter_found, "Error! SurfaceCoupling exporter data was setup, but no atm process "
                                     "of type AtmosphereProcessType::SurfaceCouplingExporter exists.\n");
  }
}

void AtmosphereDriver::reset_accumulated_fields ()
{
  constexpr Real zero = 0;
  for (auto fm_it : m_field_mgrs) {
    const auto& fm = fm_it.second;
    if (not fm->has_group("ACCUMULATED")) {
      continue;
    }

    auto accum_group = fm->get_field_group("ACCUMULATED");
    for (auto f_it : accum_group.m_fields) {
      auto& track = f_it.second->get_header().get_tracking();
      f_it.second->deep_copy(zero);
      track.set_accum_start_time(m_current_ts);
    }
  }
}

void AtmosphereDriver::setup_column_conservation_checks ()
{
  // Query m_atm_process_group if any process enables the conservation check,
  // and if not, return before creating and passing the check.
  if (not m_atm_process_group->are_column_conservation_checks_enabled()) {
    return;
  }

  auto phys_grid = m_grids_manager->get_grid("Physics");
  const auto phys_grid_name = phys_grid->name();
  auto phys_field_mgr = m_field_mgrs[phys_grid_name];

  // Get fields needed to run the mass and energy conservation checks. Require that
  // all fields exist.
  EKAT_REQUIRE_MSG (
    phys_field_mgr->has_field("pseudo_density") and
    phys_field_mgr->has_field("ps") and
    phys_field_mgr->has_field("phis") and
    phys_field_mgr->has_field("horiz_winds") and
    phys_field_mgr->has_field("T_mid") and
    phys_field_mgr->has_field("qv") and
    phys_field_mgr->has_field("qc") and
    phys_field_mgr->has_field("qr") and
    phys_field_mgr->has_field("qi") and
    phys_field_mgr->has_field("vapor_flux") and
    phys_field_mgr->has_field("water_flux") and
    phys_field_mgr->has_field("ice_flux") and
    phys_field_mgr->has_field("heat_flux"),
                   "Error! enable_column_conservation_checks=true for some atm process, "
                   "but not all fields needed for this check exist in the FieldManager.\n");

  // Get tolerances for mass and energy checks from driver_option parameters.
  auto& driver_options_pl = m_atm_params.sublist("driver_options");
  const Real mass_error_tol   = driver_options_pl.get<double>("mass_column_conservation_error_tolerance",   1e-10);
  const Real energy_error_tol = driver_options_pl.get<double>("energy_column_conservation_error_tolerance", 1e-14);

  // Create energy checker
  const auto pseudo_density = phys_field_mgr->get_field("pseudo_density");
  const auto ps             = phys_field_mgr->get_field("ps");
  const auto phis           = phys_field_mgr->get_field("phis");
  const auto horiz_winds    = phys_field_mgr->get_field("horiz_winds");
  const auto T_mid          = phys_field_mgr->get_field("T_mid");
  const auto qv             = phys_field_mgr->get_field("qv");
  const auto qc             = phys_field_mgr->get_field("qc");
  const auto qr             = phys_field_mgr->get_field("qr");
  const auto qi             = phys_field_mgr->get_field("qi");
  const auto vapor_flux     = phys_field_mgr->get_field("vapor_flux");
  const auto water_flux     = phys_field_mgr->get_field("water_flux");
  const auto ice_flux       = phys_field_mgr->get_field("ice_flux");
  const auto heat_flux      = phys_field_mgr->get_field("heat_flux");

  auto conservation_check =
    std::make_shared<MassAndEnergyColumnConservationCheck>(phys_grid,
                                                           mass_error_tol, energy_error_tol,
                                                           pseudo_density, ps, phis,
                                                           horiz_winds, T_mid, qv,
                                                           qc, qr, qi,
                                                           vapor_flux, water_flux,
                                                           ice_flux, heat_flux);

  //Get fail handling type from driver_option parameters.
  const std::string fail_handling_type_str =
      driver_options_pl.get<std::string>("column_conservation_checks_fail_handling_type", "Warning");

  CheckFailHandling fail_handling_type;
  if (fail_handling_type_str == "Warning") {
    fail_handling_type = CheckFailHandling::Warning;
  } else if (fail_handling_type_str == "Fatal") {
    fail_handling_type = CheckFailHandling::Fatal;
  } else {
    EKAT_ERROR_MSG("Error! Unknown column_conservation_checks_fail_handling_type parameter. "
                   "Acceptable types are \"Warning\" and \"Fatal\".\n");
  }

  // Pass energy checker to the process group to be added
  // to postcondition checks of appropriate processes.
  m_atm_process_group->setup_column_conservation_checks(conservation_check, fail_handling_type);
}

void AtmosphereDriver::setup_shoc_tms_links ()
{
  EKAT_REQUIRE_MSG(m_atm_process_group->has_process("tms"),
                   "Error! Attempting to setup link between "
                   "SHOC and TMS, but TMS is not defined.\n");
  EKAT_REQUIRE_MSG(m_atm_process_group->has_process("shoc"),
                   "Error! Attempting to setup link between "
                   "SHOC and TMS, but SHOC is not defined.\n");

  auto shoc_process = m_atm_process_group->get_process_nonconst("shoc");
  shoc_process->get_params().set<bool>("apply_tms", true);
}

void AtmosphereDriver::add_additional_column_data_to_property_checks () {
  // Get list of additional data fields from driver_options parameters.
  // If no fields given, return.
  using vos_t = std::vector<std::string>;
  auto additional_data_fields = m_atm_params.sublist("driver_options").get<vos_t>("property_check_data_fields",
                                                                                  {"NONE"});
  if (additional_data_fields == vos_t{"NONE"}) return;

  // Add requested fields to property checks
  auto phys_field_mgr = m_field_mgrs[m_grids_manager->get_grid("Physics")->name()];
  for (auto fname : additional_data_fields) {
    EKAT_REQUIRE_MSG(phys_field_mgr->has_field(fname), "Error! The field "+fname+" is requested for property check output "
                                                       "but does not exist in the physics field manager.\n");

    m_atm_process_group->add_additional_data_fields_to_property_checks(phys_field_mgr->get_field(fname));
  }
}

void AtmosphereDriver::create_fields()
{
  m_atm_logger->info("[EAMxx] create_fields ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_fields");

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

        // In the case of pg2, can't use create_remapper. But the remapper is
        // used only for a convenience function, anyway, create_tgt_fid.
        auto r = fvphyshack ? nullptr : m_grids_manager->create_remapper(req.src_grid,req.grid);
        // Loop over all fields in group src_name on grid src_grid.
        for (const auto& fname : rel_info->m_fields_names) {
          // Get field on src_grid
          const auto& rel_fid = rel_fm->get_field_id(fname);

          // Build a FieldRequest for the same field on greq's grid,
          // and add it to the group of this request
          if (fvphyshack) {
            auto dims = rel_fid.get_layout().dims();
            dims[0] = fm->get_grid()->get_num_local_dofs();
            FieldLayout fl(rel_fid.get_layout().tags(), dims);
            FieldIdentifier fid(rel_fid.name(), fl, rel_fid.get_units(), req.grid);
            FieldRequest freq(fid,req.name,req.pack_size);
            fm->register_field(freq);
          } else {
            const auto fid = r->create_tgt_fid(rel_fid);
            FieldRequest freq(fid,req.name,req.pack_size);
            fm->register_field(freq);
          }
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
    auto fm = m_field_mgrs.at(grid->name());
    fm->registration_ends();
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
    auto group = fm->get_field_group(it.name).get_const();
    m_atm_process_group->set_required_group(group);
  }
  for (const auto& req : m_atm_process_group->get_required_field_requests()) {
    const auto& fid = req.fid;
    auto fm = get_field_mgr(fid.get_grid_name());
    m_atm_process_group->set_required_field(fm->get_field(fid).get_const());
  }

  // Now that all processes have all the required/computed fields/groups, they
  // have also created any possible internal field (if needed). Notice that some
  // atm proc might have created internal fields already during the set_grids
  // call. However, some atm proc might need to create an internal field based
  // on the dimension of a group, or create an internal field as a subfield
  // of another one. Therefore, we had to wait till this point to query the
  // atm proc group for any internal field, and add it to the field manager
  // Besides, the field manager(s) can accept pre-built fields only after
  // the registration phase has ended.
  m_atm_process_group->gather_internal_fields();
  for (const auto& f : m_atm_process_group->get_internal_fields()) {
    const auto& fid = f.get_header().get_identifier();
    auto fm = get_field_mgr(fid.get_grid_name());
    fm->add_field(f);
  }

  // Now go through the input fields/groups to the atm proc group, as well as
  // the internal fields/groups, and mark them as part of the RESTART group.
  for (const auto& f : m_atm_process_group->get_fields_in()) {
    const auto& fid = f.get_header().get_identifier();
    const auto& gn = fid.get_grid_name();
    auto fm = get_field_mgr(gn);
    fm->add_to_group(fid.name(),"RESTART");
  }
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    const auto& gn = g.grid_name();
    auto fm = get_field_mgr(gn);
    if (g.m_bundle) {
      fm->add_to_group(g.m_bundle->get_header().get_identifier().name(),"RESTART");
    } else {
      for (const auto& fn : g.m_info->m_fields_names) {
        fm->add_to_group(fn,"RESTART");
      }
    }
  }
  for (const auto& f : m_atm_process_group->get_internal_fields()) {
    const auto& fid = f.get_header().get_identifier();
    const auto& gn = fid.get_grid_name();
    auto fm = get_field_mgr(gn);
    fm->add_to_group(fid.name(),"RESTART");
  }

  m_ad_status |= s_fields_created;

  stop_timer("EAMxx::create_fields");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_fields ... done!");
}

void AtmosphereDriver::create_output_managers () {
  m_atm_logger->info("[EAMxx] create_output_managers ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_output_managers");

  check_ad_status (s_comm_set | s_params_set | s_ts_inited);

  auto& io_params = m_atm_params.sublist("Scorpio");

  ekat::ParameterList checkpoint_params;
  checkpoint_params.set("frequency_units",std::string("never"));
  checkpoint_params.set("Frequency",-1);

  // Create model restart OutputManager first. This OM will be in charge
  // of creating rpointer.atm, while other OM's will simply append to it.
  // If this assumption is not verified, we must always append to rpointer, which
  // can make the rpointer file a bit confusing.
  if (io_params.isSublist("model_restart")) {
    // Create model restart manager
    auto params = io_params.sublist("model_restart");
    params.set<std::string>("filename_prefix",m_casename+".scream");
    params.set<std::string>("Averaging Type","Instant");
    params.sublist("provenance") = m_atm_params.sublist("provenance");

    m_restart_output_manager = std::make_shared<OutputManager>();
    m_restart_output_manager->initialize(m_atm_comm,
                                         params,
                                         m_run_t0,
                                         m_case_t0,
                                         /*is_model_restart_output*/ true);

    // Store the "Output Control" pl of the model restart as the "Checkpoint Control" for all other output streams
    checkpoint_params.set<std::string>("frequency_units",params.sublist("output_control").get<std::string>("frequency_units"));
    checkpoint_params.set("Frequency",params.sublist("output_control").get<int>("Frequency"));
  }

  // Create one output manager per output yaml file
  using vos_t = std::vector<std::string>;
  const auto& output_yaml_files = io_params.get<vos_t>("output_yaml_files",vos_t{});
  for (const auto& fname : output_yaml_files) {
    ekat::ParameterList params;
    ekat::parse_yaml_file(fname,params);
    params.rename(ekat::split(fname,"/").back());
    auto& checkpoint_pl = params.sublist("Checkpoint Control");
    checkpoint_pl.set("frequency_units",checkpoint_params.get<std::string>("frequency_units"));
    checkpoint_pl.set("Frequency",checkpoint_params.get<int>("Frequency"));

    // Check if the filename prefix for this file has already been set.  If not, use the simulation casename.
    if (not params.isParameter("filename_prefix")) {
      params.set<std::string>("filename_prefix",m_casename+".scream.h");
    }
    params.sublist("provenance") = m_atm_params.sublist("provenance");

    auto& om = m_output_managers.emplace_back();
    om.initialize(m_atm_comm,
                  params,
                  m_run_t0,
                  m_case_t0,
                  /*is_model_restart_output*/ false);
  }

  m_ad_status |= s_output_created;

  stop_timer("EAMxx::create_output_managers");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_output_managers ... done!");
}

void AtmosphereDriver::initialize_output_managers () {
  m_atm_logger->info("[EAMxx] initialize_output_managers ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::initialize_output_managers");

  check_ad_status (s_output_created | s_grids_created | s_fields_created);

  // Check for model restart output manager and setup if it exists.
  if (m_restart_output_manager) {
    if (fvphyshack) {
      // Don't save CGLL fields from ICs to the restart file.
      std::map<std::string,field_mgr_ptr> fms;
      for (auto& it : m_field_mgrs) {
        if (it.first == "Physics GLL") continue;
        fms[it.first] = it.second;
      }
      m_restart_output_manager->setup(fms, m_grids_manager);
    } else {
      m_restart_output_manager->setup(m_field_mgrs,m_grids_manager);
    }
    m_restart_output_manager->set_logger(m_atm_logger);
    for (const auto& it : m_atm_process_group->get_restart_extra_data()) {
      m_restart_output_manager->add_global(it.first,it.second);
    }
  }

  // Setup output managers
  for (auto& om : m_output_managers) {
    EKAT_REQUIRE_MSG(not om.is_restart(),
                     "Error! No restart output should be in m_output_managers. Model restart "
                     "output should be setup in m_restart_output_manager./n");

    om.set_logger(m_atm_logger);
    om.setup(m_field_mgrs,m_grids_manager);
  }

  m_ad_status |= s_output_inited;

  stop_timer("EAMxx::initialize_output_managers");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] initialize_output_managers ... done!");
}

void AtmosphereDriver::
set_provenance_data (std::string caseid,
                     std::string rest_caseid,
                     std::string hostname,
                     std::string username,
		     std::string versionid)
{
#ifdef SCREAM_CIME_BUILD
  // Check the inputs are valid
  EKAT_REQUIRE_MSG (caseid!="", "Error! Invalid case id: " + caseid + "\n");
  EKAT_REQUIRE_MSG (m_run_type==RunType::Initial or rest_caseid!="",
      "Error! Invalid restart case id: " + rest_caseid + "\n");
  EKAT_REQUIRE_MSG (hostname!="", "Error! Invalid hostname: " + hostname + "\n");
  EKAT_REQUIRE_MSG (username!="", "Error! Invalid username: " + username + "\n");
  EKAT_REQUIRE_MSG (versionid!="", "Error! Invalid version: " + versionid + "\n");
#else
  caseid = rest_caseid = m_casename;
  char* user = new char[32];
  char* host = new char[256];
  int err;
  err = gethostname(host,255);
  if (err==0) {
    hostname = std::string(host);
  } else {
    hostname = "UNKNOWN";
  }
  err = getlogin_r(user,31);
  if (err==0) {
    username = std::string(user);
  } else {
    username = "UNKNOWN";
  }
  delete[] user;
  delete[] host;
  versionid = EAMXX_GIT_VERSION;
#endif
  auto& provenance = m_atm_params.sublist("provenance");
  provenance.set("caseid",caseid);
  provenance.set("rest_caseid",rest_caseid);
  provenance.set("hostname",hostname);
  provenance.set("username",username);
  provenance.set("git_version",versionid);
}

void AtmosphereDriver::
initialize_fields ()
{
  check_ad_status (s_fields_created | s_ts_inited);

  m_atm_logger->info("[EAMxx] initialize_fields ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::initialize_fields");

  // See the [rrtmgp active gases] note in share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp
  if (fvphyshack) {
    TraceGasesWorkaround::singleton().run_type = m_run_type;
  }

  // See if we need to print a DAG. We do this first, cause if any input
  // field is missing from the initial condition file, an error will be thrown.
  // By printing the DAG first, we give the user the possibility of seeing
  // what fields are inputs to the atm time step, so he/she can fix the i.c. file.
  // TODO: would be nice to do the IC input first, and mark the fields in the
  //       DAG node "Begin of atm time step" in red if there's no initialization
  //       mechanism set for them. That is, allow field XYZ to not be found in
  //       the IC file, and throw an error when the dag is created.

  auto& driver_options_pl = m_atm_params.sublist("driver_options");
  const int verb_lvl = driver_options_pl.get<int>("atmosphere_dag_verbosity_level",-1);
  if (verb_lvl>0) {
    // Check the atm DAG for missing stuff
    AtmProcDAG dag;

    // First, add all atm processes
    dag.create_dag(*m_atm_process_group);

    // Write a dot file for visualization
    dag.write_dag("scream_atm_dag.dot",std::max(verb_lvl,0));
  }

  // Initialize fields
  if (m_run_type==RunType::Restart) {
    restart_model ();
  } else {
    set_initial_conditions ();
  }

  // Now that IC have been read, add U/V subfields of horiz_winds,
  // as well as U/V component of surf_mom_flux
  // NOTE: if you add them _before_ the IC read, set_initial_conditions
  //       will skip horiz_winds, and only process U/V, which, being
  //       missing in the IC file, would cause horiz_winds=0.
  for (auto it : m_grids_manager->get_repo()) {
    auto grid = it.second;
    auto fm = m_field_mgrs.at(grid->name());
    if (fm->has_field("horiz_winds")) {
      using namespace ShortFieldTagsNames;
      auto hw = fm->get_field("horiz_winds");
      const auto& fid = hw.get_header().get_identifier();
      const auto& layout = fid.get_layout();
      const int vec_dim = layout.get_vector_component_idx();
      const auto& units = fid.get_units();
      auto U = hw.subfield("U",units,vec_dim,0);
      auto V = hw.subfield("V",units,vec_dim,1);
      if (not fm->has_field("U")) {
        fm->add_field(U);
      }
      if (not fm->has_field("V")) {
        fm->add_field(V);
      }
    }
    if (fm->has_field("surf_mom_flux")) {
      using namespace ShortFieldTagsNames;
      auto hw = fm->get_field("surf_mom_flux");
      const auto& fid = hw.get_header().get_identifier();
      const auto& layout = fid.get_layout();
      const int vec_dim = layout.get_vector_component_idx();
      const auto& units = fid.get_units();
      auto surf_mom_flux_U = hw.subfield("surf_mom_flux_U",units,vec_dim,0);
      auto surf_mom_flux_V = hw.subfield("surf_mom_flux_V",units,vec_dim,1);
      if (not fm->has_field("surf_mom_flux_U")) {
        fm->add_field(surf_mom_flux_U);
      }
      if (not fm->has_field("surf_mom_flux_V")) {
        fm->add_field(surf_mom_flux_V);
      }
    }
  }

#ifdef SCREAM_HAS_MEMORY_USAGE
  long long my_mem_usage = get_mem_usage(MB);
  long long max_mem_usage;
  m_atm_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
  m_atm_logger->debug("[EAMxx::init::initialize_fields] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif
  stop_timer("EAMxx::initialize_fields");
  stop_timer("EAMxx::init");
  m_ad_status |= s_fields_inited;
  m_atm_logger->info("[EAMxx] initialize_fields ... done!");
}

void AtmosphereDriver::restart_model ()
{
  m_atm_logger->info("  [EAMxx] restart_model ...");

  // First, figure out the name of the netcdf file containing the restart data
  const auto& provenance = m_atm_params.sublist("provenance");
  const auto& casename = provenance.get<std::string>("rest_caseid");
  auto filename = find_filename_in_rpointer (casename+".scream",true,m_atm_comm,m_run_t0);

  m_atm_logger->info("    [EAMxx] Restart filename: " + filename);

  for (auto& it : m_field_mgrs) {
    if (fvphyshack and it.second->get_grid()->name() == "Physics GLL") continue;
    if (not it.second->has_group("RESTART")) {
      // No field needs to be restarted on this grid.
      continue;
    }
    const auto& restart_group = it.second->get_groups_info().at("RESTART");
    std::vector<std::string> fnames;
    for (const auto& fn : restart_group->m_fields_names) {
      fnames.push_back(fn);
    }
    read_fields_from_file (fnames,it.second->get_grid(),filename,m_current_ts);
  }

  // Restart the num steps counter in the atm time stamp
  int nsteps = scorpio::get_attribute<int>(filename,"GLOBAL","nsteps");
  m_current_ts.set_num_steps(nsteps);
  m_run_t0.set_num_steps(nsteps);

  for (auto& it : m_atm_process_group->get_restart_extra_data()) {
    const auto& name = it.first;
          auto& any  = it.second;

    if (any.isType<int>()) {
      ekat::any_cast<int>(any) = scorpio::get_attribute<int>(filename,"GLOBAL",name);
    } else if (any.isType<std::int64_t>()) {
      ekat::any_cast<std::int64_t>(any) = scorpio::get_attribute<std::int64_t>(filename,"GLOBAL",name);
    } else if (any.isType<float>()) {
      ekat::any_cast<float>(any) = scorpio::get_attribute<float>(filename,"GLOBAL",name);
    } else if (any.isType<double>()) {
      ekat::any_cast<double>(any) = scorpio::get_attribute<double>(filename,"GLOBAL",name);
    } else if (any.isType<std::string>()) {
      ekat::any_cast<std::string>(any) = scorpio::get_attribute<std::string>(filename,"GLOBAL",name);
    } else {
      EKAT_ERROR_MSG (
          "Error! Unrecognized/unsupported concrete type for restart extra data.\n"
          " - extra data name  : " + name + "\n"
          " - extra data typeid: " + any.content().type().name() + "\n");
    }
  }

  m_atm_logger->info("  [EAMxx] restart_model ... done!");
}

void AtmosphereDriver::create_logger () {
  using namespace ekat::logger;
  using ci_string = ekat::CaseInsensitiveString;

  auto& driver_options_pl = m_atm_params.sublist("driver_options");

  ci_string log_fname = driver_options_pl.get<std::string>("Atm Log File","atm.log");
  ci_string log_level = driver_options_pl.get<std::string>("atm_log_level","info");
  ci_string flush_level = driver_options_pl.get<std::string>("atm_flush_level","warn");
  EKAT_REQUIRE_MSG (log_fname!="",
      "Invalid string for 'Atm Log File': '" + log_fname + "'.\n");

  auto str2lev = [](const std::string& s, const std::string& name) {
    LogLevel lev;
    if (s=="trace") {
      lev = LogLevel::trace;
    } else if (s=="debug") {
      lev = LogLevel::debug;
    } else if (s=="info") {
      lev = LogLevel::info;
    } else if (s=="warn") {
      lev = LogLevel::warn;
    } else if (s=="err") {
      lev = LogLevel::err;
    } else if (s=="off") {
      lev = LogLevel::off;
    } else {
      EKAT_ERROR_MSG ("Invalid choice for '" + name + "': " + s + "\n");
    }
    return lev;
  };

  using logger_t = Logger<LogBasicFile,LogRootRank>;
  m_atm_logger = std::make_shared<logger_t>(log_fname,str2lev(log_level,"atm_log_level"),m_atm_comm,"");
  m_atm_logger->flush_on(str2lev(flush_level,"atm_flush_level"));
  m_atm_logger->set_no_format();

  // In CIME runs, this is already set to false, so atm log does not pollute e3sm.loc.
  // In standalone, we default to true, so we see output to screen.
  if (not driver_options_pl.get<bool>("output_to_screen",true)) {
     m_atm_logger->set_console_level(LogLevel::off);
  }

  // Record the CASENAME for this run, set default to EAMxx
  if (m_atm_params.isSublist("e3sm_parameters")) {
    auto e3sm_params = m_atm_params.sublist("e3sm_parameters");
    m_casename = e3sm_params.get<std::string>("e3sm_casename","EAMxx");
  } else {
    m_casename = "EAMxx";
  }
}

void AtmosphereDriver::set_initial_conditions ()
{
  m_atm_logger->info("  [EAMxx] set_initial_conditions ...");

  auto& ic_pl = m_atm_params.sublist("initial_conditions");

  // Check which fields need to have an initial condition.
  std::map<std::string,std::vector<std::string>> ic_fields_names;
  std::vector<FieldIdentifier> ic_fields_to_copy;
  std::map<std::string,std::vector<std::string>> fields_inited;

  // Check which fields should be loaded from the topography file
  std::map<std::string,std::vector<std::string>> topography_file_fields_names;
  std::map<std::string,std::vector<std::string>> topography_eamxx_fields_names;

  // Helper lambda, to reduce code duplication
  auto process_ic_field = [&](const Field& f) {
    const auto& fid = f.get_header().get_identifier();
    const auto& fname = fid.name();
    const auto& grid_name = fid.get_grid_name();

    if (ic_pl.isParameter(fname)) {
      // This is the case that the user provided an initialization
      // for this field in the parameter file.
      if (ic_pl.isType<double>(fname) or ic_pl.isType<std::vector<double>>(fname)) {
        // Initial condition is a constant
        initialize_constant_field(fid, ic_pl);

        // Note: f is const, so we can't modify the tracking. So get the same field from the fm
        auto f_nonconst = m_field_mgrs.at(grid_name)->get_field(fid.name());
        f_nonconst.get_header().get_tracking().update_time_stamp(m_current_ts);
      } else if (ic_pl.isType<std::string>(fname)) {
        // Initial condition is a string
        ic_fields_to_copy.push_back(fid);
      } else {
        EKAT_ERROR_MSG ("ERROR: invalid assignment for variable " + fname + ", only scalar "
                        "double or string, or vector double arguments are allowed");
      }
      fields_inited[grid_name].push_back(fname);
    } else if (fname == "phis" or fname == "sgh30") {
      // Both phis and sgh30 need to be loaded from the topography file
      auto& this_grid_topo_file_fnames = topography_file_fields_names[grid_name];
      auto& this_grid_topo_eamxx_fnames = topography_eamxx_fields_names[grid_name];

      if (fname == "phis") {
        // For GLL points, phis corresponds to "PHIS_d" in the
        // topography file. On PG2 grid, dynamics will take care
        // of computing phis, so do not add to initialized fields.
        if (grid_name == "Physics PG2") {
          // Skip
        } else if (grid_name == "Physics GLL" ||
                   grid_name == "Point Grid") {
          this_grid_topo_file_fnames.push_back("PHIS_d");
          this_grid_topo_eamxx_fnames.push_back(fname);
          fields_inited[grid_name].push_back(fname);
        } else {
          EKAT_ERROR_MSG ("Error! Requesting phis on an unknown grid: " + grid_name + ".\n");
        }
      } else if (fname == "sgh30") {
        // The eamxx field "sgh30" is called "SGH30" in the
        // topography file and is only available on the PG2 grid.
        EKAT_ASSERT_MSG(grid_name == "Physics PG2",
                        "Error! Requesting sgh30 field on " + grid_name +
                        " topo file only has sgh30 for Physics PG2.\n");
        topography_file_fields_names[grid_name].push_back("SGH30");
        topography_eamxx_fields_names[grid_name].push_back(fname);
	fields_inited[grid_name].push_back(fname);
      }
    } else if (not (fvphyshack and grid_name == "Physics PG2")) {
      // The IC file is written for the GLL grid, so we only load
      // fields from there. Any other input fields on the PG2 grid
      // will be properly computed in the dynamics interface.
      auto& this_grid_ic_fnames = ic_fields_names[grid_name];
      auto c = f.get_header().get_children();
      if (c.size()==0) {
        // If this field is the parent of other subfields, we only read from file the subfields.
        if (not ekat::contains(this_grid_ic_fnames,fname)) {
          this_grid_ic_fnames.push_back(fname);
	  fields_inited[grid_name].push_back(fname);
        }
      } else if (fvphyshack and grid_name == "Physics GLL") {
        // [CGLL ICs in pg2] I tried doing something like this in
        // HommeDynamics::set_grids, but I couldn't find the means to get the
        // list of fields. I think the issue is that you can't access group
        // objects until some registration period ends. So instead do it here,
        // where the list is definitely available.
        for (const auto& e : c) {
          const auto f = e.lock();
          const auto& fid = f->get_identifier();
          const auto& fname = fid.name();
          if (ic_pl.isParameter(fname) and ic_pl.isType<double>(fname)) {
            initialize_constant_field(fid, ic_pl);
          } else {
            this_grid_ic_fnames.push_back(fname);
          }
	  fields_inited[grid_name].push_back(fname);
        }
      }
    }
  };

  // First the individual input fields...
  m_atm_logger->debug("    [EAMxx] Processing input fields ...");
  for (const auto& f : m_atm_process_group->get_fields_in()) {
    process_ic_field (f);
  }
  m_atm_logger->debug("    [EAMxx] Processing input fields ... done!");

  // ...then the input groups
  m_atm_logger->debug("    [EAMxx] Processing input groups ...");
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    if (g.m_bundle) {
      process_ic_field(*g.m_bundle);
    }
    for (auto it : g.m_fields) {
      process_ic_field(*it.second);
    }
  }
  m_atm_logger->debug("    [EAMxx] Processing input groups ... done!");

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

  if (m_iop) {
    // For runs with IOP, call to setup io grids and lat
    // lon information needed for reading from file
    // We use a single topo file for both GLL and PG2 runs. All
    // lat/lon data in topo file is defined in terms of PG2 grid,
    // so if we have a topo field on GLL grid, we need to setup
    // io info using the IC file (which is always GLL).
    for (const auto& it : m_field_mgrs) {
      const auto& grid_name = it.first;
      if (ic_fields_names[grid_name].size() > 0 or
	  topography_eamxx_fields_names[grid_name].size() > 0) {
        const auto& file_name = grid_name == "Physics GLL"
                                ?
                                ic_pl.get<std::string>("Filename")
                                :
                                ic_pl.get<std::string>("topography_filename");
        m_iop->setup_io_info(file_name, it.second->get_grid());
      }
    }
  }

  // If a filename is specified, use it to load inputs on all grids
  if (ic_pl.isParameter("Filename")) {
    // Now loop over all grids, and load from file the needed fields on each grid (if any).
    const auto& file_name = ic_pl.get<std::string>("Filename");
    m_atm_logger->info("    [EAMxx] IC filename: " + file_name);
    for (const auto& it : m_field_mgrs) {
      const auto& grid_name = it.first;
      if (not m_iop) {
        read_fields_from_file (ic_fields_names[grid_name],it.second->get_grid(),file_name,m_current_ts);
      } else {
        // For IOP enabled, we load from file and copy data from the closest
        // lat/lon column to every other column
        m_iop->read_fields_from_file_for_iop(file_name,
                                             ic_fields_names[grid_name],
                                             m_current_ts,
                                             it.second);
      }
    }
  }

  // If there were any fields that needed to be copied per the input yaml file, now we copy them.
  m_atm_logger->debug("    [EAMxx] Processing fields to copy ...");
  for (const auto& tgt_fid : ic_fields_to_copy) {
    const auto& tgt_fname = tgt_fid.name();
    const auto& gname = tgt_fid.get_grid_name();

    const auto& src_fname = ic_pl.get<std::string>(tgt_fname);

    auto fm = get_field_mgr(gname);

    // The field must exist in the fm on the input field's grid
    EKAT_REQUIRE_MSG (fm->has_field(src_fname),
        "Error! Source field for initial condition not found in the field manager.\n"
        "       Grid name:     " + gname + "\n"
        "       Field to init: " + tgt_fname + "\n"
        "       Source field:  " + src_fname + " (NOT FOUND)\n");

    // Get the two fields, and copy src to tgt
    auto f_tgt = fm->get_field(tgt_fname);
    auto f_src = fm->get_field(src_fname);
    f_tgt.deep_copy(f_src);

    // Set the initial time stamp
    f_tgt.get_header().get_tracking().update_time_stamp(m_current_ts);
  }
  m_atm_logger->debug("    [EAMxx] Processing fields to copy ... done!");

  // It is possible to have a bundled group G1=(f1,f2,f3),
  // where the IC are read from file for f1, f2, and f3. In that case,
  // the time stamp for the bundled G1 has not be inited, but the data
  // is valid (all entries have been inited). Let's fix that.
  m_atm_logger->debug("    [EAMxx] Processing subfields ...");
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    if (g.m_bundle) {
      auto& track = g.m_bundle->get_header().get_tracking();
      if (not track.get_time_stamp().is_valid()) {
        // The bundled field has not been inited. Check if all the subfields
        // have been inited. If so, init the timestamp of the bundled field too.
        const auto& children = track.get_children();
        bool all_inited = children.size()>0; // If no children, then something is off, so mark as not good
        for (auto wp : children) {
          auto sp = wp.lock();
          if (not sp->get_time_stamp().is_valid()) {
            all_inited = false;
            break;
          }
        }
        if (all_inited) {
          track.update_time_stamp(m_current_ts);
        }
      }
    }
  }
  m_atm_logger->debug("    [EAMxx] Processing subfields ... done!");

  // Load topography from file if topography file is given.
  if (ic_pl.isParameter("topography_filename")) {
    m_atm_logger->info("    [EAMxx] Reading topography from file ...");
    const auto& file_name = ic_pl.get<std::string>("topography_filename");
    m_atm_logger->info("        filename: " + file_name);
    for (const auto& it : m_field_mgrs) {
      const auto& grid_name = it.first;
      if (not m_iop) {
        // Topography files always use "ncol_d" for the GLL grid value of ncol.
        // To ensure we read in the correct value, we must change the name for that dimension
        auto io_grid = it.second->get_grid();
        if (grid_name=="Physics GLL") {
          using namespace ShortFieldTagsNames;
          auto grid = io_grid->clone(io_grid->name(),true);
          grid->reset_field_tag_name(COL,"ncol_d");
          io_grid = grid;
        }
        read_fields_from_file (topography_file_fields_names[grid_name],
                               topography_eamxx_fields_names[grid_name],
                               io_grid,file_name,m_current_ts);
      } else {
        // For IOP enabled, we load from file and copy data from the closest
        // lat/lon column to every other column
        m_iop->read_fields_from_file_for_iop(file_name,
                                             topography_file_fields_names[grid_name],
                                             topography_eamxx_fields_names[grid_name],
                                             m_current_ts,
                                             it.second);
      }
    }
    // Store in provenance list, for later usage in output file metadata
    m_atm_params.sublist("provenance").set("topography_file",file_name);
    m_atm_logger->debug("    [EAMxx] Processing topography from file ... done!");
  } else {
    // Ensure that, if no topography_filename is given, no
    // processes is asking for topography data (assuming a
    // separate IC param entry isn't given for the field).
    for (const auto& it : m_field_mgrs) {
      const auto& grid_name = it.first;
      EKAT_REQUIRE_MSG(topography_file_fields_names[grid_name].size()==0,
                      "Error! Topography data was requested in the FM, but no "
                      "topography_filename or entry matching the field name "
                      "was given in IC parameters.\n");
    }

    m_atm_params.sublist("provenance").set<std::string>("topography_file","NONE");
  }

  if (m_iop) {
    // Load IOP data file data for initial time stamp
    m_iop->read_iop_file_data(m_current_ts);

    // Now that ICs are processed, set appropriate fields using IOP file data.
    // Since ICs are loaded on GLL grid, we set those fields only and dynamics
    // will take care of the rest (for PG2 case).
    if (m_field_mgrs.count("Physics GLL") > 0) {
      const auto& fm = m_field_mgrs.at("Physics GLL");
      m_iop->set_fields_from_iop_data(fm);
    }
  }

  // Compute IC perturbations of GLL fields (if requested)
  using vos = std::vector<std::string>;
  const auto perturbed_fields = ic_pl.get<vos>("perturbed_fields", {});
  const auto num_perturb_fields = perturbed_fields.size();
  if (num_perturb_fields > 0) {
    m_atm_logger->info("    [EAMxx] Adding random perturbation to ICs ...");

    EKAT_REQUIRE_MSG(m_field_mgrs.count("Physics GLL") > 0,
                     "Error! Random perturbation can only be applied to fields on "
                     "the GLL grid, but no Physics GLL FieldManager was defined.\n");
    const auto& fm = m_field_mgrs.at("Physics GLL");

    // Setup RNG. There are two relevant params: generate_perturbation_random_seed and
    // perturbation_random_seed. We have 3 cases:
    //   1. Parameter generate_perturbation_random_seed is set true, assert perturbation_random_seed
    //      is not given and generate a random seed using std::rand() to get an integer random value.
    //   2. Parameter perturbation_random_seed is given, use this value for the seed.
    //   3. Parameter perturbation_random_seed is not given and generate_perturbation_random_seed is
    //      not given, use 0 as the random seed.
    // Case 3 is considered the default (using seed=0).
    int seed;
    if (ic_pl.get<bool>("generate_perturbation_random_seed", false)) {
      EKAT_REQUIRE_MSG(not ic_pl.isParameter("perturbation_random_seed"),
                       "Error! Param generate_perturbation_random_seed=true, and "
                       "a perturbation_random_seed is given. Only one of these can "
                       "be defined for a simulation.\n");
      std::srand(std::time(nullptr));
      seed = std::rand();
    } else {
      seed = ic_pl.get<int>("perturbation_random_seed", 0);
    }
    m_atm_logger->info("      For IC perturbation, random seed: "+std::to_string(seed));
    std::mt19937_64 engine(seed);

    // Get perturbation limit. Defines a range [1-perturbation_limit, 1+perturbation_limit]
    // for which the perturbation value will be randomly generated from. Create a uniform
    // distribution for this range.
    const auto perturbation_limit = ic_pl.get<Real>("perturbation_limit", 0.001);
    std::uniform_real_distribution<Real> pdf(1-perturbation_limit, 1+perturbation_limit);

    // Define a level mask using reference pressure and the perturbation_minimum_pressure parameter.
    // This mask dictates which levels we apply a perturbation.
    const auto gll_grid = m_grids_manager->get_grid("Physics GLL");
    const auto hyam_h = gll_grid->get_geometry_data("hyam").get_view<const Real*, Host>();
    const auto hybm_h = gll_grid->get_geometry_data("hybm").get_view<const Real*, Host>();
    constexpr auto ps0 = physics::Constants<Real>::P0;
    const auto min_pressure = ic_pl.get<Real>("perturbation_minimum_pressure", 1050.0);
    auto pressure_mask = [&] (const int ilev) {
      const auto pref = (hyam_h(ilev)*ps0 + hybm_h(ilev)*ps0)/100; // Reference pressure ps0 is in Pa, convert to millibar
      return pref > min_pressure;
    };

    // Loop through fields and apply perturbation.
    for (size_t f=0; f<perturbed_fields.size(); ++f) {
      const auto fname = perturbed_fields[f];
      EKAT_REQUIRE_MSG(ekat::contains(fields_inited[fm->get_grid()->name()], fname),
                       "Error! Attempting to apply perturbation to field not in initial_conditions.\n"
                       "  - Field: "+fname+"\n"
                       "  - Grid:  "+fm->get_grid()->name()+"\n");

      auto field = fm->get_field(fname);
      perturb(field, engine, pdf, seed, pressure_mask, fm->get_grid()->get_dofs_gids());
    }

    m_atm_logger->info("    [EAMxx] Adding random perturbation to ICs ... done!");
  }

  m_atm_logger->info("  [EAMxx] set_initial_conditions ... done!");
}

void AtmosphereDriver::
read_fields_from_file (const std::vector<std::string>& field_names_nc,
                       const std::vector<std::string>& field_names_eamxx,
                       const std::shared_ptr<const AbstractGrid>& grid,
                       const std::string& file_name,
                       const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(field_names_nc.size()==field_names_eamxx.size(),
                   "Error! Field name arrays must have same size.\n");

  if (field_names_nc.size()==0) {
    return;
  }

  // NOTE: we cannot pass the field_mgr and m_grids_mgr, since the input
  //       grid may not be in the grids_manager and may not be the grid
  //       of the field mgr. This sounds weird, but there is a precise
  //       use case: when grid is a shallow clone of the fm grid, where
  //       we changed the name of some field tags (e.g., we set the name
  //       of COL to ncol_d). This is used when reading the topography,
  //       since the topo file *always* uses ncol_d for GLL points data,
  //       while a non-PG2 run would have the tag name be "ncol".
  const auto& field_mgr = m_field_mgrs.at(grid->name());
  std::vector<Field> fields;
  for (size_t i=0; i<field_names_nc.size(); ++i) {
    const auto& eamxx_name = field_names_eamxx[i];
    const auto& nc_name    = field_names_nc[i];
    fields.push_back(field_mgr->get_field(eamxx_name).alias(nc_name));
  }

  AtmosphereInput ic_reader(file_name,grid,fields);
  ic_reader.set_logger(m_atm_logger);
  ic_reader.read_variables();
  ic_reader.finalize();

  for (auto& f : fields) {
    // Set the initial time stamp
    // NOTE: f is an alias of the field from field_mgr, so it shares all
    //       pointers to the metadata (except for the FieldIdentifier),
    //       so changing its timestamp will also change the timestamp
    //       of the field in field_mgr
    f.get_header().get_tracking().update_time_stamp(t0);
  }
}

void AtmosphereDriver::
read_fields_from_file (const std::vector<std::string>& field_names,
                       const std::shared_ptr<const AbstractGrid>& grid,
                       const std::string& file_name,
                       const util::TimeStamp& t0)
{
  if (field_names.size()==0) {
    return;
  }

  // NOTE: we cannot pass the field_mgr and m_grids_mgr, since the input
  //       grid may not be in the grids_manager and may not be the grid
  //       of the field mgr. This sounds weird, but there is a precise
  //       use case: when grid is a shallow clone of the fm grid, where
  //       we changed the name of some field tags (e.g., we set the name
  //       of COL to ncol_d). This is used when reading the topography,
  //       since the topo file *always* uses ncol_d for GLL points data,
  //       while a non-PG2 run would have the tag name be "ncol".
  const auto& field_mgr = m_field_mgrs.at(grid->name());
  std::vector<Field> fields;
  for (const auto& fn : field_names) {
    fields.push_back(field_mgr->get_field(fn));
  }

  AtmosphereInput ic_reader(file_name,grid,fields);
  ic_reader.set_logger(m_atm_logger);
  ic_reader.read_variables();
  ic_reader.finalize();

  for (auto& f : fields) {
    // Set the initial time stamp
    f.get_header().get_tracking().update_time_stamp(t0);
  }
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

  // For vector fields, we allow either single value init or vector value init.
  // That is, both these are ok
  //   fname: val
  //   fname: [val1,...,valN]
  // In the first case, all entries of the field are inited to val, while in the latter,
  // each component is inited to the corresponding entry of the array.
  if (layout.is_vector_layout() and ic_pl.isType<std::vector<double>>(name)) {
    const auto idim = layout.get_vector_component_idx();
    const auto vec_dim = layout.get_vector_dim();
    const auto& values = ic_pl.get<std::vector<double>>(name);
    EKAT_REQUIRE_MSG (values.size()==static_cast<size_t>(vec_dim),
        "Error! Initial condition values array for '" + name + "' has the wrong dimension.\n"
        "       Field dimension: " + std::to_string(vec_dim) + "\n"
        "       Array dimenions: " + std::to_string(values.size()) + "\n");

    if (layout.rank()==2 && idim==1) {
      // We cannot use 'get_component' for views of rank 2 with vector dimension
      // striding fastest, since we would not get a LayoutRight view. For these views,
      // simply do a manual loop
      using kt = Field::kt_dev;
      typename kt::view_1d<double> data("data",vec_dim);
      auto data_h = Kokkos::create_mirror_view(data);
      for (int i=0; i<vec_dim; ++i) {
        data_h(i) = values[i];
      }
      Kokkos::deep_copy(data,data_h);

      const int n = layout.dim(0);
      auto v = f.get_view<double**>();
      Kokkos::parallel_for(typename kt::RangePolicy(0,n),
                           KOKKOS_LAMBDA(const int i) {
        for (int j=0; j<vec_dim; ++j) {
          v(i,j) = data(j);
        }
      });
    } else {
      // Extract a subfield for each component. This is not "too" expensive, expecially
      // considering that this code is executed during initialization only.
      for (int comp=0; comp<vec_dim; ++comp) {
        auto f_i = f.get_component(comp);
        f_i.deep_copy(values[comp]);
      }
    }
  } else {
    const auto& value = ic_pl.get<double>(name);
    f.deep_copy(value);
  }
}

void AtmosphereDriver::initialize_atm_procs ()
{
  m_atm_logger->info("[EAMxx] initialize_atm_procs ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::initialize_atm_procs");

  // Initialize memory buffer for all atm processes
  m_memory_buffer = std::make_shared<ATMBufferManager>();
  m_memory_buffer->request_bytes(m_atm_process_group->requested_buffer_size_in_bytes());
  m_memory_buffer->allocate();
  m_atm_process_group->init_buffers(*m_memory_buffer);

  // Setup SurfaceCoupling import and export (if they exist)
  if (m_surface_coupling_import_data_manager || m_surface_coupling_export_data_manager) {
    setup_surface_coupling_processes();
  }

  // Initialize the processes
  m_atm_process_group->initialize(m_current_ts, m_run_type);

  // Create and add energy and mass conservation check to appropriate atm procs
  setup_column_conservation_checks();

  // If user requests it, we set up NaN checks for all computed fields after each atm proc run
  if (m_atm_params.sublist("driver_options").get("check_all_computed_fields_for_nans",true)) {
    m_atm_process_group->add_postcondition_nan_checks();
  }

  // Add additional column data fields to pre/postcondition checks (if they exist)
  add_additional_column_data_to_property_checks();

  if (fvphyshack) {
    // [CGLL ICs in pg2] See related notes in atmosphere_dynamics.cpp.
    const auto gn = "Physics GLL";
    m_field_mgrs[gn]->clean_up();
    m_field_mgrs.erase(gn);
  }

  m_ad_status |= s_procs_inited;

  stop_timer("EAMxx::initialize_atm_procs");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] initialize_atm_procs ... done!");

  report_res_dep_memory_footprint ();
}

void AtmosphereDriver::
initialize (const ekat::Comm& atm_comm,
            const ekat::ParameterList& params,
            const util::TimeStamp& run_t0,
            const util::TimeStamp& case_t0)
{
  set_comm(atm_comm);
  set_params(params);
  set_provenance_data ();

  init_scorpio ();

  init_time_stamps (run_t0, case_t0);

  create_output_managers ();

  create_atm_processes ();

  create_grids ();

  create_fields ();

  initialize_fields ();

  initialize_atm_procs ();

  // Do this before init-ing the output managers,
  // so the fields are valid if outputing at t=0
  reset_accumulated_fields();

  initialize_output_managers ();
}

void AtmosphereDriver::run (const int dt) {
  start_timer("EAMxx::run");

  // DEBUG option: Check if user has set the run to fail at a specific timestep.
  auto& debug = m_atm_params.sublist("driver_debug_options");
  auto fail_step = debug.get<int>("force_crash_nsteps",-1);
  if (fail_step==m_current_ts.get_num_steps()) {
    std::abort();
  }

  // Make sure the end of the time step is after the current start_time
  EKAT_REQUIRE_MSG (dt>0, "Error! Input time step must be positive.\n");

  // Print current timestamp information
  m_atm_logger->log(ekat::logger::LogLevel::info,
    "Atmosphere step = " + std::to_string(m_current_ts.get_num_steps()) + "\n" +
    "  model start-of-step time = " + m_current_ts.get_date_string() + " " + m_current_ts.get_time_string() + "\n");

  // Reset accum fields to 0
  // Note: at the 1st timestep this is redundant, since we did it at init,
  //       to ensure t=0 INSTANT output was correct. However, it's not a
  //       very expensive operation, so it's not worth the effort of the
  //       nano-opt of removing the call for the 1st timestep.
  reset_accumulated_fields();

  // Tell the output managers that we're starting a timestep. This is usually
  // a no-op, but some diags *may* require to do something. E.g., a diag that
  // computes tendency of an arbitrary quantity may want to store a copy of
  // that quantity at the beginning of the timestep. Or they may need to store
  // the timestamp at the beginning of the timestep, so that we can compute
  // dt at the end.
  if (m_restart_output_manager) m_restart_output_manager->init_timestep(m_current_ts, dt);
  for (auto& it : m_output_managers) {
    it.init_timestep(m_current_ts,dt);
  }

  // The class AtmosphereProcessGroup will take care of dispatching arguments to
  // the individual processes, which will be called in the correct order.
  m_atm_process_group->run(dt);

  // Some accumulated fields need to be divided by dt at the end of the atm step
  for (auto fm_it : m_field_mgrs) {
    const auto& fm = fm_it.second;
    if (not fm->has_group("DIVIDE_BY_DT")) {
      continue;
    }

    auto rescale_group = fm->get_field_group("DIVIDE_BY_DT");
    for (auto f_it : rescale_group.m_fields) {
      f_it.second->scale(Real(1) / dt);
    }
  }

  // Update current time stamps
  m_current_ts += dt;

  // Update output streams
  m_atm_logger->debug("[EAMxx::run] running output managers...");
  if (m_restart_output_manager) m_restart_output_manager->run(m_current_ts);
  for (auto& out_mgr : m_output_managers) {
    out_mgr.run(m_current_ts);
  }

#ifdef SCREAM_HAS_MEMORY_USAGE
  long long my_mem_usage = get_mem_usage(MB);
  long long max_mem_usage;
  m_atm_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
  m_atm_logger->info("[EAMxx::run] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif

  // Flush the logger at least once per time step.
  // Without this flush, depending on how much output we are loggin,
  // it might be several time steps before the file is updated.
  // This way, we give the user a chance to follow the log more real-time.
  m_atm_logger->flush();

  stop_timer("EAMxx::run");
}

void AtmosphereDriver::finalize ( /* inputs? */ ) {
  start_timer("EAMxx::finalize");

  if (m_ad_status==0) {
    return;
  }

  m_atm_logger->info("[EAMxx] Finalize ...");

  // Finalize and destroy output streams, make sure files are closed
  if (m_restart_output_manager) {
    m_restart_output_manager->finalize();
    m_restart_output_manager = nullptr;
  }
  for (auto& out_mgr : m_output_managers) {
    out_mgr.finalize();
  }
  m_output_managers.clear();

  // Finalize, and then destroy all atmosphere processes
  if (m_atm_process_group.get()) {
    m_atm_process_group->finalize( /* inputs ? */ );
    m_atm_process_group = nullptr;
  }

  // Destroy iop
  m_iop = nullptr;

  // Destroy the buffer manager
  m_memory_buffer = nullptr;

  // Destroy the surface coupling data managers
  m_surface_coupling_import_data_manager = nullptr;
  m_surface_coupling_export_data_manager = nullptr;

  // Destroy the grids manager
  m_grids_manager = nullptr;

  // Destroy all the fields manager
  for (auto it : m_field_mgrs) {
    it.second->clean_up();
  }

  // Write all timers to file, and possibly finalize gptl
  if (not m_gptl_externally_handled) {
    write_timers_to_file (m_atm_comm,"scream_timing.txt");
    finalize_gptl();
  }

  // Finalize scorpio. Check, just in case we're calling finalize after
  // an exception, thrown before the AD (and scorpio) was inited
  if (scorpio::is_subsystem_inited()) {
    scorpio::finalize_subsystem();
  }

  m_atm_logger->info("[EAMxx] Finalize ... done!");

#ifdef SCREAM_HAS_MEMORY_USAGE
  long long my_mem_usage = get_mem_usage(MB);
  long long max_mem_usage;
  m_atm_comm.all_reduce(&my_mem_usage,&max_mem_usage,1,MPI_MAX);
  m_atm_logger->debug("[EAMxx::finalize] memory usage: " + std::to_string(max_mem_usage) + "MB");
#endif
  m_atm_logger->flush();

  m_ad_status = 0;

  stop_timer("EAMxx::finalize");
}

AtmosphereDriver::field_mgr_ptr
AtmosphereDriver::get_field_mgr (const std::string& grid_name) const {
  EKAT_REQUIRE_MSG (m_ad_status & s_grids_created,
      "Error! Field manager(s) are created *after* the grids.\n");
  // map::at would throw, but you won't know which map threw.
  // With our own msg, we can tell you where the throw happened.
  EKAT_REQUIRE_MSG(m_grids_manager->has_grid(grid_name),
      "Error! Request for field manager on a non-existing grid '" + grid_name + "'.\n");

  // Do not use $grid_name, since it might be an alias. Fetch grid,
  // then get the actual name from the grid itself
  auto grid = m_grids_manager->get_grid(grid_name);

  return m_field_mgrs.at(grid->name());
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

void AtmosphereDriver::report_res_dep_memory_footprint () const {
  // Log the amount of memory used that is linked to the grid(s) sizes
  long long my_dev_mem_usage = 0;
  long long my_host_mem_usage = 0;
  long long max_dev_mem_usage, max_host_mem_usage;

  // The first report includes memory used by 1) fields (metadata excluded),
  // 2) grids data (dofs, maps, geo views), 3) atm buff manager, and 4) IO.

  // Fields
  for (const auto& fm_it : m_field_mgrs) {
    for (const auto& it : *fm_it.second) {
      const auto& fap = it.second->get_header().get_alloc_properties();
      if (fap.is_subfield()) {
        continue;
      }
      my_dev_mem_usage += fap.get_alloc_size();
      my_host_mem_usage += fap.get_alloc_size();
    }
  }
  // Grids
  for (const auto& it : m_grids_manager->get_repo()) {
    const auto& grid = it.second;
    const int nldofs = grid->get_num_local_dofs();

    my_dev_mem_usage += sizeof(AbstractGrid::gid_type)*nldofs;

    my_dev_mem_usage += sizeof(int)*grid->get_lid_to_idx_map().get_header().get_identifier().get_layout().size();

    const auto& geo_names = grid->get_geometry_data_names();
    my_dev_mem_usage += sizeof(Real)*geo_names.size()*nldofs;
  }
  // Atm buffer
  my_dev_mem_usage += m_memory_buffer->allocated_bytes();
  // Output
  if (m_restart_output_manager) {
    const auto om_footprint = m_restart_output_manager->res_dep_memory_footprint();
    my_dev_mem_usage += om_footprint;
    my_host_mem_usage += om_footprint;
  }
  for (const auto& om : m_output_managers) {
    const auto om_footprint = om.res_dep_memory_footprint ();
    my_dev_mem_usage += om_footprint;
    my_host_mem_usage += om_footprint;
  }

  m_atm_comm.all_reduce(&my_dev_mem_usage,&max_dev_mem_usage,1,MPI_MAX);
  m_atm_logger->info("[EAMxx::init] resolution-dependent device memory footprint: " + std::to_string(max_dev_mem_usage/1e6) + "MB");

  if (not std::is_same<HostDevice,DefaultDevice>::value) {
    m_atm_comm.all_reduce(&my_host_mem_usage,&max_host_mem_usage,1,MPI_MAX);
    m_atm_logger->info("[EAMxx::init] resolution-dependent host memory footprint: " + std::to_string(max_host_mem_usage/1e6) + "MB");
  }

  // The following is a memory usage based on probing some OS tools
#ifdef SCREAM_HAS_MEMORY_USAGE
  long long my_mem_usage_from_os = get_mem_usage(MB);
  long long max_mem_usage_from_os;
  m_atm_comm.all_reduce(&my_mem_usage_from_os,&max_mem_usage_from_os,1,MPI_MAX);
  m_atm_logger->info("[EAMxx::init] memory usage from OS probing tools: " + std::to_string(max_mem_usage_from_os) + "MB");
#endif
}

}  // namespace control
}  // namespace scream
