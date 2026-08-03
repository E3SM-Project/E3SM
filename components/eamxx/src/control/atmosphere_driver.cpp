#include "control/atmosphere_driver.hpp"
#include "control/atmosphere_surface_coupling_importer.hpp"
#include "control/atmosphere_surface_coupling_exporter.hpp"
#include "control/default_model_init.hpp"

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/atm_process/atmosphere_process_dag.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/property_checks/mass_and_energy_conservation_check.hpp"
#include "share/core/eamxx_config.hpp"
#include "eamxx_version.h"

#include <ekat_assert.hpp>
#include <ekat_string_utils.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_yaml.hpp>
#include <ekat_std_utils.hpp>

// The global variable fvphyshack is used to help the initial pgN implementation
// work around some current AD constraints. Search the code for "fvphyshack" to
// find blocks that eventually should be removed in favor of a design that
// accounts for pg2. Some blocks may turn out to be unnecessary, and I simply
// didn't realize I could do without the workaround.
#include "share/algorithm/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"

#ifndef SCREAM_CIME_BUILD
#include <unistd.h>
#endif



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
 *     to allocate fields. For more details, see the documentation in the share/data_managers/field_request.hpp header.
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
 *  - for field manager        -> src/share/data_managers/field_manager.hpp
 *  - for field groups         -> src/share/field/field_group.hpp
 *  - for field/group requests -> src/share/data_managers/field_request.hpp
 *  - for grid                 -> src/share/grid/abstract_grid.hpp
 *  - for grid manager         -> src/share/data_managers/grids_manager.hpp
 *  - for atm proc             -> src/share/atm_process/atmosphere_process.hpp
 *  - for atm proc group       -> src/share/atm_process/atmosphere_process_group.hpp
 *  - for scorpio input/output -> src/share/io/scorpio_[input|output].hpp
 *  - for output manager       -> src/share/io/eamxx_output_manager.hpp
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
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

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
    case 2:
      m_run_type = RunType::Restart;
      m_branch_run = true;
      break;
    case -1:
      m_run_type = case_t0==run_t0 ? RunType::Initial : RunType::Restart; break;
    default:
      EKAT_ERROR_MSG ("Unsupported/unrecognized run_type: " + std::to_string(run_type) + "\n");
  }

  // If it is a restarted run, make sure we have the correct num steps.
  // If num_steps is left at default (0), it messes up output managers logic
  if (m_run_type==RunType::Restart) {
    // Figure out the name of the netcdf file containing the restart data
    const auto& provenance = m_atm_params.sublist("provenance");
    const auto& casename = provenance.get<std::string>("rest_caseid");
    auto filename = find_filename_in_rpointer (casename+".scream",true,m_atm_comm,m_run_t0);

    // Restart the num steps counter in the atm time stamp
    int nsteps = scorpio::get_attribute<int>(filename,"GLOBAL","nsteps");
    m_run_t0.set_num_steps(nsteps);
    m_current_ts.set_num_steps(nsteps);
  }
}

void AtmosphereDriver::
setup_iop_data_manager ()
{
  // At this point, must have comm, params, initialized timestamps created.
  check_ad_status(s_comm_set | s_params_set | s_ts_inited);

  // Check to make sure iop is not already initialized
  EKAT_REQUIRE_MSG(not m_iop_data_manager, "Error! setup_iop_data_manager() is called, but IOP already set up.\n");

  // This function should only be called if we are enabling IOP
  const bool enable_iop =
    m_atm_params.sublist("driver_options").get("enable_iop", false);
  EKAT_REQUIRE_MSG(enable_iop, "Error! setup_iop_data_manager() is called, but enable_iop=false "
                               "in driver_options parameters.\n");

  // Params must include iop_options sublist.
  const auto iop_sublist_exists = m_atm_params.isSublist("iop_options");
  EKAT_REQUIRE_MSG(iop_sublist_exists,
                   "Error! setup_iop_data_manager() is called, but no iop_options "
                   "defined in parameters.\n");

  const auto iop_params = m_atm_params.sublist("iop_options");
  const auto phys_grid = m_grids_manager->get_grid("physics");
  const auto nlevs = phys_grid->get_num_vertical_levels();
  const auto hyam = phys_grid->get_geometry_data("hyam");
  const auto hybm = phys_grid->get_geometry_data("hybm");

  m_iop_data_manager = std::make_shared<IOPDataManager>(m_atm_comm,
                                                        iop_params,
                                                        m_run_t0,
                                                        nlevs,
                                                        hyam,
                                                        hybm);

  // Set IOP object in atm processes
  m_atm_process_group->set_iop_data_manager(m_iop_data_manager);
}

void AtmosphereDriver::create_atm_processes()
{
  m_atm_logger->info("[EAMxx] create_atm_processes  ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_atm_processes");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

  // At this point, must have comm and params set.
  check_ad_status(s_comm_set | s_params_set);

  // Create the group of processes. This will recursively create the processes
  // tree, storing also the information regarding parallel execution (if needed).
  // See AtmosphereProcessGroup class documentation for more details.
  auto& atm_proc_params = m_atm_params.sublist("eamxx");
  atm_proc_params.rename("EAMxx");
  atm_proc_params.set("logger",m_atm_logger);
  m_atm_process_group = std::make_shared<AtmosphereProcessGroup>(m_atm_comm,atm_proc_params);

  m_ad_status |= s_procs_created;
  stop_timer("EAMxx::create_atm_processes");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_atm_processes  ... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
}

void AtmosphereDriver::create_grids()
{
  m_atm_logger->info("[EAMxx] create_grids ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_grids");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

  // Must have procs created by now (and comm/params set)
  check_ad_status (s_procs_created | s_comm_set | s_params_set | s_ts_inited);

  // Create the grids manager
  auto& gm_params = m_atm_params.sublist("grids_manager");
  const std::string& gm_type = gm_params.get<std::string>("type");

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
  } else if (ic_pl.isParameter("filename")) {
    // Initial run, if an IC file is present, pass it.
    auto filename = ic_pl.get<std::string>("filename");
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
    setup_iop_data_manager ();
  }

  // Set the grids in the processes. Do this by passing the grids manager.
  // Each process will grab what they need
  m_atm_process_group->set_grids(m_grids_manager);

  m_ad_status |= s_grids_created;

  stop_timer("EAMxx::create_grids");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_grids ... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
}

void AtmosphereDriver::setup_surface_coupling_data_manager(SurfaceCouplingTransferType transfer_type,
                                                           const int num_cpl_fields, const int num_scream_fields,
                                                           const int field_size, Real* data_ptr,
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
  for (auto grid_it : m_grids_manager->get_repo()) {
    const auto& grid_name = grid_it.second->name();
    if (not m_field_mgr->has_group("ACCUMULATED", grid_name)) {
      continue;
    }

    auto accum_group = m_field_mgr->get_field_group("ACCUMULATED", grid_name);
    for (auto f_it : accum_group.m_individual_fields) {
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
  if (not m_atm_process_group->are_conservation_checks_enabled()) {
    return;
  }

  auto phys_grid = m_grids_manager->get_grid("physics");
  const auto phys_grid_name = phys_grid->name();
  // Get fields needed to run the mass and energy conservation checks. Require that
  // all fields exist (except h2otemp which is optional).
  EKAT_REQUIRE_MSG (
    m_field_mgr->has_field("pseudo_density", phys_grid_name) and
    m_field_mgr->has_field("ps",             phys_grid_name) and
    m_field_mgr->has_field("phis",           phys_grid_name) and
    m_field_mgr->has_field("horiz_winds",    phys_grid_name) and
    m_field_mgr->has_field("T_mid",          phys_grid_name) and
    m_field_mgr->has_field("qv",             phys_grid_name) and
    m_field_mgr->has_field("qc",             phys_grid_name) and
    m_field_mgr->has_field("qr",             phys_grid_name) and
    m_field_mgr->has_field("qi",             phys_grid_name) and
    m_field_mgr->has_field("vapor_flux",     phys_grid_name) and
    m_field_mgr->has_field("water_flux",     phys_grid_name) and
    m_field_mgr->has_field("ice_flux",       phys_grid_name) and
    m_field_mgr->has_field("heat_flux",      phys_grid_name),
    "Error! enable_column_conservation_checks=true for some atm process, "
    "but not all fields needed for this check exist in the FieldManager.\n");

  // Get tolerances for mass and energy checks from driver_option parameters.
  auto& driver_options_pl = m_atm_params.sublist("driver_options");
  const Real mass_error_tol   = driver_options_pl.get<double>("mass_column_conservation_error_tolerance",   1e-10);
  const Real energy_error_tol = driver_options_pl.get<double>("energy_column_conservation_error_tolerance", 1e-14);

  // Create energy checker
  const auto pseudo_density = m_field_mgr->get_field("pseudo_density", phys_grid_name);
  const auto ps             = m_field_mgr->get_field("ps",             phys_grid_name);
  const auto phis           = m_field_mgr->get_field("phis",           phys_grid_name);
  const auto horiz_winds    = m_field_mgr->get_field("horiz_winds",    phys_grid_name);
  const auto T_mid          = m_field_mgr->get_field("T_mid",          phys_grid_name);
  const auto qv             = m_field_mgr->get_field("qv",             phys_grid_name);
  const auto qc             = m_field_mgr->get_field("qc",             phys_grid_name);
  const auto qr             = m_field_mgr->get_field("qr",             phys_grid_name);
  const auto qi             = m_field_mgr->get_field("qi",             phys_grid_name);
  const auto vapor_flux     = m_field_mgr->get_field("vapor_flux",     phys_grid_name);
  const auto water_flux     = m_field_mgr->get_field("water_flux",     phys_grid_name);
  const auto ice_flux       = m_field_mgr->get_field("ice_flux",       phys_grid_name);
  const auto heat_flux      = m_field_mgr->get_field("heat_flux",      phys_grid_name);
  
  // h2otemp is optional - only available when SurfaceCouplingImporter is active
  Field h2otemp;
  if (m_field_mgr->has_field("h2otemp", phys_grid_name)) {
    h2otemp = m_field_mgr->get_field("h2otemp", phys_grid_name);
  }

  auto conservation_check =
    std::make_shared<MassAndEnergyConservationCheck>(m_atm_comm,phys_grid,
                                                           mass_error_tol, energy_error_tol,
                                                           pseudo_density, ps, phis,
                                                           horiz_winds, T_mid, qv,
                                                           qc, qr, qi,
                                                           vapor_flux, water_flux,
                                                           ice_flux, heat_flux,
                                                           h2otemp);

  //Get fail handling type from driver_option parameters.
  const std::string fail_handling_type_str =
      driver_options_pl.get<std::string>("column_conservation_checks_fail_handling_type", "warning");

  CheckFailHandling fail_handling_type = CheckFailHandling::Warning;
  if (fail_handling_type_str == "warning") {
    fail_handling_type = CheckFailHandling::Warning;
  } else if (fail_handling_type_str == "fatal") {
    fail_handling_type = CheckFailHandling::Fatal;
  } else {
    EKAT_ERROR_MSG("Error! Unknown column_conservation_checks_fail_handling_type parameter. "
                   "Acceptable types are \"warning\" and \"fatal\".\n");
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
  auto additional_data_fields = m_atm_params.sublist("driver_options").get<strvec_t>("property_check_data_fields",
                                                                                  {"NONE"});
  if (additional_data_fields == strvec_t{"NONE"}) return;

  // Add requested fields to property checks
  const auto& grid_name = m_grids_manager->get_grid("physics")->name();
  for (auto fname : additional_data_fields) {
    EKAT_REQUIRE_MSG(m_field_mgr->has_field(fname, grid_name),
      "Error! The field " + fname + " is requested for property "
      "check output but does not exist on the physics grid.\n");
    m_atm_process_group->add_additional_data_fields_to_property_checks(m_field_mgr->get_field(fname, grid_name));
  }
}

void AtmosphereDriver::create_fields()
{
  m_atm_logger->info("[EAMxx] create_fields ...");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_fields");

  // Must have grids and procs at this point
  check_ad_status (s_procs_created | s_grids_created);

  // Create FM
  m_field_mgr = std::make_shared<field_mgr_type>(m_grids_manager);

  // Before registering fields, check that Field Requests for tracers are compatible
  // and store the correct type of turbulence advection for each tracer
  m_atm_process_group->pre_process_tracer_requests();

  // Register required/computed fields. By now, the processes should have
  // fully built the ids of their required/computed fields and groups
  for (const auto& req : m_atm_process_group->get_field_requests()) {
    m_field_mgr->register_field(req);
  }
  for (const auto& greq : m_atm_process_group->get_group_requests()) {
    m_field_mgr->register_group(greq);
  }

  // Closes the FM, allocate all fields
  m_field_mgr->registration_ends();

  // Set all the fields/groups in the processes. Input fields/groups will be handed
  // to the processes with const scalar type (const Real), to prevent them from
  // overwriting them (though, they can always cast away const...).
  for (const auto& req : m_atm_process_group->get_field_requests()) {
    const auto& fid = req.fid;
    const auto& f = m_field_mgr->get_field(fid);
    if (req.usage & Required)
      m_atm_process_group->set_required_field(f.get_const());
    if (req.usage & Computed)
      m_atm_process_group->set_computed_field(f);
  }
  for (const auto& req : m_atm_process_group->get_group_requests()) {
    auto group = m_field_mgr->get_field_group(req.name, req.grid);
    if (req.usage & Required)
      m_atm_process_group->set_required_group(group.get_const());
    if (req.usage & Computed)
      m_atm_process_group->set_computed_group(group);
  }

  // Make atm procs create the proc-level tendency fields (if requested)
  m_atm_process_group->setup_step_tendencies(m_grids_manager->get_grid("physics")->name());

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
    m_field_mgr->add_field(f);

    // Internal fields have their group names set by the processes that create them.
    // Hence, simply add them to all the groups they are marked as part of
    const auto& ftrack = f.get_header().get_tracking();
    const auto& fid    = f.get_header().get_identifier();
    for (const auto& gn : ftrack.get_groups_names()) {
      m_field_mgr->add_to_group(fid, gn);
    }
  }

  // Now go through the input fields/groups to the atm proc group,
  // and mark them as part of the RESTART/STARTUP groups.
  // Skip fields in the ACCUMULATED group, since those are reset to 0
  // at the beginning of each atm step, so there is no need to read
  // them from the IC or restart file.
  auto is_topography_field = [] (const std::string& name) {
    return name=="phis" or name=="sgh" or name=="sgh30";
  };
  for (const auto& f : m_atm_process_group->get_fields_in()) {
    const auto& fid = f.get_header().get_identifier();
    const auto& fname = fid.name();
    const auto& fgroups = f.get_header().get_tracking().get_groups_names();
    if (not ekat::contains(fgroups, "ACCUMULATED")) {
      m_field_mgr->add_to_group(fid, "RESTART");
      m_field_mgr->add_to_group(fid, "STARTUP");
      if (is_topography_field(fname)) {
        m_field_mgr->add_to_group(fid, "TOPOGRAPHY");
      }
    }
  }
  for (const auto& g : m_atm_process_group->get_groups_in()) {
    if (g.m_monolithic_field) {
      const auto& mf = *g.m_monolithic_field;
      const auto& mfn = mf.get_header().get_identifier().name();
      const auto& mfgroups = mf.get_header().get_tracking().get_groups_names();
      if (not ekat::contains(mfgroups, "ACCUMULATED")) {
        m_field_mgr->add_to_group(mf.get_header().get_identifier(), "RESTART");
        m_field_mgr->add_to_group(mf.get_header().get_identifier(), "STARTUP");
        if (is_topography_field(mfn)) {
          m_field_mgr->add_to_group(mf.get_header().get_identifier(), "TOPOGRAPHY");
        }
      }
    } else {
      for (const auto& fn : g.m_info->m_fields_names) {
        auto field = m_field_mgr->get_field(fn, g.grid_name());
        const auto& fgroups = field.get_header().get_tracking().get_groups_names();
        if (not ekat::contains(fgroups, "ACCUMULATED")) {
          m_field_mgr->add_to_group(fn, g.grid_name(), "RESTART");
          m_field_mgr->add_to_group(fn, g.grid_name(), "STARTUP");
          if (is_topography_field(fn)) {
            m_field_mgr->add_to_group(fn, g.grid_name(), "TOPOGRAPHY");
          }
        }
      }
    }
  }

  auto& driver_options_pl = m_atm_params.sublist("driver_options");
  const int verb_lvl = driver_options_pl.get<int>("atmosphere_dag_verbosity_level",-1);
  if (verb_lvl>0) {
    // now that we've got fields, generate a DAG with fields and dependencies
    // NOTE: at this point, fields provided by initial conditions may (will)
    // appear as unmet dependencies
    AtmProcDAG dag;
    // First, add all atm processes
    dag.create_dag(*m_atm_process_group);
    // Write a dot file for visualizing the DAG
    if (m_atm_comm.am_i_root()) {
      std::string filename = "scream_atm_createField_dag";
      if (is_scream_standalone()) {
        filename += ".np" + std::to_string(m_atm_comm.size());
      }
      filename += ".dot";
      dag.write_dag(filename, verb_lvl);
    }
  }

  m_ad_status |= s_fields_created;

  // If the user requested it, we can save a dictionary of the FM fields to file
  if (driver_options_pl.get("save_field_manager_content",false)) {
    auto grid_name = m_grids_manager->get_grid("physics")->name();
    auto& phys_fields = m_field_mgr->get_repo(grid_name);
    ekat::ParameterList pl_out("field_manager_content");
    pl_out.sublist("provenance") = m_atm_params.sublist("provenance");
    DefaultMetadata std_names;
    std::string desc;
    desc = "content of the EAMxx FieldManager corresponding to the 'physics' grid.\n"
           "The dict keys are the field names as used in EAMxx.\n"
           "For each field, we add the following entries:\n"
           "  - standard_name: the name commonly used to refer to this field in atm sciences (if applicable)\n"
           "  - units: the units for this field used in EAMxx\n"
           "  - layout: the names of the dimensions for this field (time excluded)\n"
           "  - providers: the atm processes that update/compute this field\n"
           "  - customers: the atm processes that require this field as an input\n";
    pl_out.set("description", desc);
    auto& dict = pl_out.sublist("fields");
    for (const auto& it : phys_fields) {
      const auto& fid = it.second->get_header().get_identifier();
      auto& pl = dict.sublist(fid.name());

      pl.set("units",fid.get_units().to_string());
      pl.set("layout",fid.get_layout().names());
      pl.set("standard_name",std_names.get_standardname(fid.name()));
      const auto& track = it.second->get_header().get_tracking();
      const auto& p = track.get_providers();
      const auto& c = track.get_customers();
      pl.set("providers",strvec_t(p.begin(),p.end()));
      pl.set("customers",strvec_t(c.begin(),c.end()));
    }

    ekat::write_yaml_file("eamxx_field_manager_content.yaml",pl_out);
  }

  stop_timer("EAMxx::create_fields");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] create_fields ... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
}

void AtmosphereDriver::create_output_managers () {
  m_atm_logger->info("[EAMxx] create_output_managers ...");
  start_timer("EAMxx::init");
  start_timer("EAMxx::create_output_managers");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

  check_ad_status (s_comm_set | s_params_set | s_ts_inited);

  auto& io_params = m_atm_params.sublist("scorpio");

  ekat::ParameterList checkpoint_params;
  checkpoint_params.set("frequency_units",std::string("never"));
  checkpoint_params.set("frequency",-1);

  // Create model restart OutputManager first. This OM will be in charge
  // of creating rpointer.atm, while other OM's will simply append to it.
  // If this assumption is not verified, we must always append to rpointer, which
  // can make the rpointer file a bit confusing.
  if (io_params.isSublist("model_restart")) {
    // Create model restart manager
    auto params = io_params.sublist("model_restart");
    params.set<std::string>("filename_prefix",m_casename+".scream");
    params.set<std::string>("averaging_type","instant");
    params.sublist("provenance") = m_atm_params.sublist("provenance");

    m_restart_output_manager = std::make_shared<OutputManager>();
    m_restart_output_manager->initialize(m_atm_comm,
                                         params,
                                         m_run_t0,
                                         m_case_t0,
                                         /*is_model_restart_output*/ true);

    // Store the "Output Control" pl of the model restart as the "checkpoint_control" for all other output streams
    checkpoint_params.set<std::string>("frequency_units",params.sublist("output_control").get<std::string>("frequency_units"));
    checkpoint_params.set("frequency",params.sublist("output_control").get<int>("frequency"));
  }

  // Create one output manager per output yaml file
  const auto& output_yaml_files = io_params.get<strvec_t>("output_yaml_files",strvec_t{});
  for (const auto& fname : output_yaml_files) {
    ekat::ParameterList params;
    ekat::parse_yaml_file(fname,params);
    params.rename(ekat::split(fname,"/").back());
    auto& checkpoint_pl = params.sublist("checkpoint_control");
    checkpoint_pl.set("frequency_units",checkpoint_params.get<std::string>("frequency_units"));
    checkpoint_pl.set("frequency",checkpoint_params.get<int>("frequency"));

    // Check if the filename prefix for this file has already been set.  If not, use the simulation casename.
    if (not params.isParameter("filename_prefix")) {
      params.set<std::string>("filename_prefix",m_casename+".scream.h");
    }
    params.sublist("provenance") = m_atm_params.sublist("provenance");
    params.sublist("restart").set("branch_run",m_branch_run);

    if (not params.isParameter("enable_fine_grain_timers")) {
      params.set("enable_fine_grain_timers",io_params.get("enable_fine_grain_timers",false));
    }

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
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
}

void AtmosphereDriver::initialize_output_managers () {
  m_atm_logger->info("[EAMxx] initialize_output_managers ...");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
  start_timer("EAMxx::init");
  start_timer("EAMxx::initialize_output_managers");

  check_ad_status (s_output_created | s_grids_created | s_fields_created);

  // Check for model restart output manager and setup if it exists.
  if (m_restart_output_manager) {
    auto output_grids = m_grids_manager->get_grid_names();

    // Don't save CGLL fields from ICs to the restart file if we are running with PG2.
    if (fvphyshack and output_grids.find("physics_gll")!=output_grids.end()) {
      output_grids.erase("physics_gll");
    }

    m_restart_output_manager->set_logger(m_atm_logger);
    m_restart_output_manager->setup(m_field_mgr, output_grids);
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
    om.setup(m_field_mgr,m_grids_manager->get_grid_names());
  }

  m_ad_status |= s_output_inited;

  stop_timer("EAMxx::initialize_output_managers");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] initialize_output_managers ... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)
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

  // Create the model init object (if not already set externally).
  if (m_model_init == nullptr) {
    register_default_model_init();
    auto model_init_type = std::string("default");
    const auto gm_type = m_atm_params.sublist("grids_manager").get<std::string>("type");
    if (gm_type == "homme") {
      model_init_type = "homme";
    }
    if (m_atm_params.isSublist("model_init")) {
      auto& mi_pl = m_atm_params.sublist("model_init");
      model_init_type = mi_pl.get<std::string>("type", model_init_type);
    }
    m_model_init = ModelInitFactory::instance().create(model_init_type, m_atm_params);
  }

  // Provide all data the model init needs, then run it.
  m_model_init->setup(m_atm_comm, m_atm_params, m_field_mgr, m_grids_manager,
                      m_atm_process_group, m_iop_data_manager, m_atm_logger,
                      m_current_ts, m_run_type);
  m_model_init->run();

  // Capture which fields were initialized (used later for DAG output).
  m_fields_inited = m_model_init->get_fields_inited();

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

void AtmosphereDriver::create_logger () {
  using namespace ekat::logger;
  using ci_string = ekat::CaseInsensitiveString;

  auto& driver_options_pl = m_atm_params.sublist("driver_options");

  ci_string log_fname = driver_options_pl.get<std::string>("Atm Log File","atm.log");
  if (is_scream_standalone()) {
    // Useful for concurrent unit tests in the same folder
    log_fname += ".np" + std::to_string (m_atm_comm.size());
  }
  ci_string log_level = driver_options_pl.get<std::string>("atm_log_level","info");
  ci_string flush_level = driver_options_pl.get<std::string>("atm_flush_level","warn");
  EKAT_REQUIRE_MSG (log_fname!="",
      "Invalid string for 'Atm Log File': '" + log_fname + "'.\n");

  auto str2lev = [](const std::string& s, const std::string& name) {
    LogLevel lev = LogLevel::info;
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

void AtmosphereDriver::initialize_atm_procs ()
{
  m_atm_logger->info("[EAMxx] initialize_atm_procs ...");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

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
    const auto gn = "physics_gll";
    m_field_mgr->clean_up(gn);
  }

  m_ad_status |= s_procs_inited;

  stop_timer("EAMxx::initialize_atm_procs");
  stop_timer("EAMxx::init");
  m_atm_logger->info("[EAMxx] initialize_atm_procs ... done!");
  m_atm_logger->flush(); // During init, flush often (to help debug crashes)

  report_res_dep_memory_footprint ();

  auto& driver_options_pl = m_atm_params.sublist("driver_options");
  const int verb_lvl = driver_options_pl.get<int>("atmosphere_dag_verbosity_level",-1);
  if (verb_lvl>0) {
    // now that we've got fields, generate a DAG with fields and dependencies
    // NOTE: at this point, fields provided by initial conditions may (will)
    // appear as unmet dependencies
    AtmProcDAG dag;
    // First, add all atm processes
    dag.create_dag(*m_atm_process_group);
    // process the initial conditions to maybe fulfill unmet dependencies
    dag.process_initial_conditions(m_fields_inited);
    // Write a dot file for visualizing the DAG
    if (m_atm_comm.am_i_root()) {
      std::string filename = "scream_atm_initProc_dag";
      if (is_scream_standalone()) {
        filename += ".np" + std::to_string(m_atm_comm.size());
      }
      filename += ".dot";
      dag.write_dag(filename, verb_lvl);
    }
  }
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
  auto end_of_step = m_current_ts + dt;
  m_atm_logger->log(ekat::logger::LogLevel::info,
    "Atmosphere step = " + std::to_string(end_of_step.get_num_steps()) + "\n" +
    "  model beg-of-step timestamp: " + m_current_ts.get_date_string() + " " + m_current_ts.get_time_string() + "\n"
    "  model end-of-step timestamp: " + end_of_step.get_date_string() + " " + end_of_step.get_time_string() + "\n");

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
  for (auto gname : m_grids_manager->get_grid_names()) {
    if (not m_field_mgr->has_group("DIVIDE_BY_DT", gname)) {
      continue;
    }

    auto rescale_group = m_field_mgr->get_field_group("DIVIDE_BY_DT", gname);
    for (auto f_it : rescale_group.m_individual_fields) {
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
  m_iop_data_manager = nullptr;
  m_model_init = nullptr;

  // Destroy the buffer manager
  m_memory_buffer = nullptr;

  // Destroy the surface coupling data managers
  m_surface_coupling_import_data_manager = nullptr;
  m_surface_coupling_export_data_manager = nullptr;

  // Destroy the grids manager
  m_grids_manager = nullptr;

  // Destroy all the fields manager
  m_field_mgr = nullptr;

  // Write all timers to file, and possibly finalize gptl
  if (not m_gptl_externally_handled) {
    write_timers_to_file (m_atm_comm,"eamxx_timing.txt");
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
  for (auto gname : m_field_mgr->get_grids_manager()->get_grid_names()) {
    for (const auto& it : m_field_mgr->get_repo(gname)) {
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
