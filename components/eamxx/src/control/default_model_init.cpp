#include "control/default_model_init.hpp"

#include "share/atm_process/atmosphere_process_group.hpp"
#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/grids_manager.hpp"
#include "share/data_managers/IOPDataManager.hpp"
#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/algorithm/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"

#include <ekat_assert.hpp>
#include <ekat_std_utils.hpp>
#include <ekat_logger.hpp>

#include <fstream>
#include <random>

namespace scream {

DefaultModelInit::
DefaultModelInit (const ekat::ParameterList& /* params */)
{
  // Data is provided via setup(); nothing to do here.
}

void DefaultModelInit::
run ()
{
  if (m_run_type == RunType::Restart) {
    m_logger->info("  [EAMxx] restart_model ...");
    m_logger->flush();
    restart_model();
    m_logger->info("  [EAMxx] restart_model ... done!");
  } else {
    m_logger->info("  [EAMxx] set_initial_conditions ...");
    m_logger->flush();
    set_initial_conditions();
    m_logger->info("  [EAMxx] set_initial_conditions ... done!");
  }

  // Now that IC/restart data have been read, add U/V subfields of
  // horiz_winds and surf_mom_flux.
  // NOTE: adding them before the IC read would cause set_initial_conditions
  //       to skip horiz_winds and only process U/V, which — being absent from
  //       the IC file — would leave horiz_winds zeroed out.
  for (auto it : m_grids_mgr->get_repo()) {
    auto grid = it.second;
    auto gn   = grid->name();
    if (m_field_mgr->has_field("horiz_winds", gn)) {
      using namespace ShortFieldTagsNames;
      auto hw = m_field_mgr->get_field("horiz_winds", gn);
      const auto& fid     = hw.get_header().get_identifier();
      const auto& layout  = fid.get_layout();
      const int   vec_dim = layout.get_vector_component_idx();
      const auto& units   = fid.get_units();
      auto U = hw.subfield("U", units, vec_dim, 0);
      auto V = hw.subfield("V", units, vec_dim, 1);
      if (not m_field_mgr->has_field("U", gn)) {
        m_field_mgr->add_field(U);
      }
      if (not m_field_mgr->has_field("V", gn)) {
        m_field_mgr->add_field(V);
      }
    }
    if (m_field_mgr->has_field("surf_mom_flux", gn)) {
      using namespace ShortFieldTagsNames;
      auto smf    = m_field_mgr->get_field("surf_mom_flux", gn);
      const auto& fid     = smf.get_header().get_identifier();
      const auto& layout  = fid.get_layout();
      const int   vec_dim = layout.get_vector_component_idx();
      const auto& units   = fid.get_units();
      auto surf_mom_flux_U = smf.subfield("surf_mom_flux_U", units, vec_dim, 0);
      auto surf_mom_flux_V = smf.subfield("surf_mom_flux_V", units, vec_dim, 1);
      if (not m_field_mgr->has_field("surf_mom_flux_U", gn)) {
        m_field_mgr->add_field(surf_mom_flux_U);
      }
      if (not m_field_mgr->has_field("surf_mom_flux_V", gn)) {
        m_field_mgr->add_field(surf_mom_flux_V);
      }
    }
  }
}

void DefaultModelInit::
restart_model ()
{
  // The restart filename was already resolved by the driver during
  // create_grids() and stored in provenance for reuse here.
  const auto& provenance = m_atm_params.sublist("provenance");
  const auto& filename   = provenance.get<std::string>("initial_conditions_file");

  m_logger->info("    [EAMxx] Restart filename: " + filename);

  for (auto& gn : m_grids_mgr->get_grid_names()) {
    if (fvphyshack and gn == "physics_gll") continue;
    if (not m_field_mgr->has_group("RESTART", gn)) {
      // No field needs to be restarted on this grid.
      continue;
    }
    const auto& restart_fnames =
        m_field_mgr->get_group_info("RESTART", gn).m_fields_names;
    std::vector<Field> fields;
    for (const auto& fn : restart_fnames) {
      // If the field has a parent that is also in the RESTART group,
      // skip it — restarting the parent will restart the child too.
      auto f = m_field_mgr->get_field(fn, gn);
      auto p = f.get_header().get_parent();
      if (p and ekat::contains(p->get_tracking().get_groups_names(), "RESTART")) {
        continue;
      }
      fields.push_back(m_field_mgr->get_field(fn, gn));
    }
    auto grid = m_grids_mgr->get_grid(gn);
    read_fields(filename, fields, grid->get_partitioned_dim_gids(), m_comm);
    for (auto& f : fields) {
      f.get_header().get_tracking().update_time_stamp(m_current_ts);
    }
  }

  for (auto& it : m_atm_procs->get_restart_extra_data()) {
    const auto& name =  it.first;
          auto& any  = *it.second;

    if (any.type() == typeid(int)) {
      std::any_cast<int&>(any) =
          scorpio::get_attribute<int>(filename, "GLOBAL", name);
    } else if (any.type() == typeid(std::int64_t)) {
      std::any_cast<std::int64_t&>(any) =
          scorpio::get_attribute<std::int64_t>(filename, "GLOBAL", name);
    } else if (any.type() == typeid(float)) {
      std::any_cast<float&>(any) =
          scorpio::get_attribute<float>(filename, "GLOBAL", name);
    } else if (any.type() == typeid(double)) {
      std::any_cast<double&>(any) =
          scorpio::get_attribute<double>(filename, "GLOBAL", name);
    } else if (any.type() == typeid(std::string)) {
      std::any_cast<std::string&>(any) =
          scorpio::get_attribute<std::string>(filename, "GLOBAL", name);
    } else {
      EKAT_ERROR_MSG (
          "Error! Unrecognized/unsupported concrete type for restart extra data.\n"
          " - extra data name  : " + name + "\n"
          " - extra data typeid: " + std::string(any.type().name()) + "\n");
    }
  }
}

void DefaultModelInit::
set_initial_conditions ()
{
  m_logger->flush(); // During init, flush often (to help debug crashes)

  // See the [rrtmgp active gases] note in
  // share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp
  if (fvphyshack) {
    TraceGasesWorkaround::singleton().run_type = m_run_type;
  }

  auto& ic_pl = m_atm_params.sublist("initial_conditions");

  // Check which fields need an initial condition.
  strmap_t ic_fields_names;
  std::vector<FieldIdentifier> ic_fields_to_copy;

  // Check which fields should be loaded from the topography file.
  strmap_t topography_file_fields_names;
  strmap_t topography_eamxx_fields_names;

  // Helper lambda to reduce code duplication
  auto process_ic_field = [&](const Field& f) {
    const auto& fid       = f.get_header().get_identifier();
    const auto& fname     = fid.name();
    const auto& grid_name = fid.get_grid_name();

    if (ic_pl.isParameter(fname)) {
      // User provided a constant or copy-field initialization in the yaml.
      if (ic_pl.isType<int>(fname) or ic_pl.isType<double>(fname) or
          ic_pl.isType<std::vector<double>>(fname)) {
        initialize_constant_field(fid, ic_pl);

        auto f_nonconst = m_field_mgr->get_field(fid);
        f_nonconst.get_header().get_tracking().update_time_stamp(m_current_ts);
      } else if (ic_pl.isType<std::string>(fname)) {
        ic_fields_to_copy.push_back(fid);
      } else {
        EKAT_ERROR_MSG (
            "ERROR: invalid assignment for variable " + fname +
            ", only scalar double or string, or vector double arguments are allowed");
      }
      m_fields_inited[grid_name].push_back(fname);
    } else if (fname == "phis" or fname == "sgh30" or fname == "sgh") {
      // These fields are loaded from the topography file.
      auto& this_grid_topo_file_fnames  = topography_file_fields_names[grid_name];
      auto& this_grid_topo_eamxx_fnames = topography_eamxx_fields_names[grid_name];

      if (fname == "phis") {
        // For GLL points, phis corresponds to "PHIS_d" in the topography file.
        // On PG2, dynamics computes phis, so skip it.
        if (grid_name == "physics_pg2") {
          // Skip
        } else if (grid_name == "physics_gll" || grid_name == "point_grid") {
          this_grid_topo_file_fnames.push_back("PHIS_d");
          this_grid_topo_eamxx_fnames.push_back(fname);
          m_fields_inited[grid_name].push_back(fname);
        } else {
          EKAT_ERROR_MSG (
              "Error! Requesting phis on an unknown grid: " + grid_name + ".\n");
        }
      } else if (fname == "sgh30") {
        EKAT_ASSERT_MSG(grid_name == "physics_pg2",
                        "Error! Requesting sgh30 field on " + grid_name +
                        " topo file only has sgh30 for physics_pg2.\n");
        topography_file_fields_names[grid_name].push_back("SGH30");
        topography_eamxx_fields_names[grid_name].push_back(fname);
        m_fields_inited[grid_name].push_back(fname);
      } else if (fname == "sgh") {
        EKAT_ASSERT_MSG(grid_name == "physics_pg2",
                        "Error! Requesting sgh field on " + grid_name +
                        " topo file only has sgh for physics_pg2.\n");
        topography_file_fields_names[grid_name].push_back("SGH");
        topography_eamxx_fields_names[grid_name].push_back(fname);
        m_fields_inited[grid_name].push_back(fname);
      }
    } else if (not (fvphyshack and grid_name == "physics_pg2")) {
      // The IC file is written for the GLL grid, so we only load fields from
      // there. Any other input fields on the PG2 grid will be properly
      // computed in the dynamics interface.
      auto& this_grid_ic_fnames = ic_fields_names[grid_name];
      auto  c = f.get_header().get_children();
      if (c.size() == 0) {
        if (not ekat::contains(this_grid_ic_fnames, fname)) {
          this_grid_ic_fnames.push_back(fname);
          m_fields_inited[grid_name].push_back(fname);
        }
      } else if (fvphyshack and grid_name == "physics_gll") {
        // [CGLL ICs in pg2] I tried doing something like this in
        // HommeDynamics::set_grids, but I couldn't find the means to get the
        // list of fields. I think the issue is that you can't access group
        // objects until some registration period ends. So instead do it here,
        // where the list is definitely available.
        for (const auto& e : c) {
          const auto sf = e.lock();
          const auto& sfid  = sf->get_identifier();
          const auto& sfname = sfid.name();
          if (ic_pl.isParameter(sfname) and ic_pl.isType<double>(sfname)) {
            initialize_constant_field(sfid, ic_pl);
          } else {
            this_grid_ic_fnames.push_back(sfname);
          }
          m_fields_inited[grid_name].push_back(sfname);
        }
      }
    }
  };

  // Process all startup fields from the FM startup groups.
  m_logger->debug("    [EAMxx] Processing startup groups ...");
  for (const auto& it : m_grids_mgr->get_repo()) {
    const auto& grid      = it.second;
    const auto& grid_name = grid->name();
    if (not m_field_mgr->has_group("STARTUP", grid_name)) {
      continue;
    }
    const auto& startup_fnames =
        m_field_mgr->get_group_info("STARTUP", grid_name).m_fields_names;
    for (const auto& fn : startup_fnames) {
      process_ic_field(m_field_mgr->get_field(fn, grid_name));
    }
  }
  m_logger->debug("    [EAMxx] Processing startup groups ... done!");

  // If a field is a subfield of a group monolithic field that is already
  // inited (e.g. via initialize_constant_field or copied), remove it from
  // the list of fields to load from file.
  for (auto& it1 : ic_fields_names) {
    const auto& grid_name = it1.first;

    bool run_again = true;
    while (run_again) {
      run_again = false;
      auto& names = it1.second;
      for (auto it2 = names.begin(); it2 != names.end(); ++it2) {
        const auto& fn = *it2;
        auto f = m_field_mgr->get_field(fn, grid_name);
        auto p = f.get_header().get_parent();
        if (p) {
          const auto& pname = p->get_identifier().name();
          if (ekat::contains(m_fields_inited[grid_name], pname)) {
            names.erase(it2);
            run_again = true;
            break;
          }
        }
      }
    }
  }

  if (m_iop_data_mgr) {
    // For IOP runs, set up io grids and lat/lon info needed for reading
    // from file.  The topo file always uses PG2 lat/lon, so if we have
    // topo fields on the GLL grid we use the IC file for io info.
    for (const auto& it : m_grids_mgr->get_repo()) {
      const auto& grid      = it.second;
      const auto& grid_name = grid->name();
      if (ic_fields_names[grid_name].size() > 0 or
          topography_eamxx_fields_names[grid_name].size() > 0) {
        const auto& file_name =
            grid_name == "physics_gll"
            ? ic_pl.get<std::string>("filename")
            : ic_pl.get<std::string>("topography_filename");
        m_iop_data_mgr->setup_io_info(file_name, grid);
      }
    }
  }

  // If a filename is specified, load IC fields from it on all grids.
  if (ic_pl.isParameter("filename")) {
    const auto& file_name = ic_pl.get<std::string>("filename");
    m_logger->info("    [EAMxx] IC filename: " + file_name);

    for (const auto& it : m_grids_mgr->get_repo()) {
      const auto& grid      = it.second;
      const auto& grid_name = grid->name();

      if (ic_fields_names[grid_name].size() == 0) continue;

      std::vector<Field> ic_fields;
      for (const auto& fn : ic_fields_names[grid_name]) {
        ic_fields.push_back(m_field_mgr->get_field(fn, grid_name));
      }
      if (not m_iop_data_mgr) {
        read_fields(file_name, ic_fields,
                    grid->get_partitioned_dim_gids(), m_comm);
      } else {
        m_iop_data_mgr->read_fields_from_file_for_iop(
            file_name, ic_fields, grid);
      }
      for (auto& f : ic_fields) {
        f.get_header().get_tracking().update_time_stamp(m_current_ts);
      }
    }
  }

  // Copy fields that were specified as string references in the IC yaml.
  m_logger->debug("    [EAMxx] Processing fields to copy ...");
  for (const auto& tgt_fid : ic_fields_to_copy) {
    const auto& tgt_fname = tgt_fid.name();
    const auto& gname     = tgt_fid.get_grid_name();
    const auto& src_fname = ic_pl.get<std::string>(tgt_fname);

    EKAT_REQUIRE_MSG (m_field_mgr->has_field(src_fname, gname),
        "Error! Source field for initial condition not found in the field manager.\n"
        "       Grid name:     " + gname + "\n"
        "       Field to init: " + tgt_fname + "\n"
        "       Source field:  " + src_fname + " (NOT FOUND)\n");

    auto f_tgt = m_field_mgr->get_field(tgt_fname, gname);
    auto f_src = m_field_mgr->get_field(src_fname, gname);
    f_tgt.deep_copy(f_src);
    f_tgt.get_header().get_tracking().update_time_stamp(m_current_ts);
  }
  m_logger->debug("    [EAMxx] Processing fields to copy ... done!");

  // If all subfields of a monolithically-allocated group have been inited,
  // stamp the monolithic field's time stamp too.
  m_logger->debug("    [EAMxx] Processing subfields ...");
  for (const auto& g : m_atm_procs->get_groups_in()) {
    if (g.m_monolithic_field) {
      auto& track = g.m_monolithic_field->get_header().get_tracking();
      if (not track.get_time_stamp().is_valid()) {
        const auto& children = track.get_children();
        bool all_inited = children.size() > 0;
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
  m_logger->debug("    [EAMxx] Processing subfields ... done!");

  // Load topography from file.
  if (ic_pl.isParameter("topography_filename")) {
    m_logger->info("    [EAMxx] Reading topography from file ...");
    const auto& file_name = ic_pl.get<std::string>("topography_filename");
    m_logger->info("        filename: " + file_name);

    for (const auto& it : m_grids_mgr->get_repo()) {
      const auto& grid      = it.second;
      const auto& grid_name = grid->name();

      int nfields = topography_eamxx_fields_names[grid_name].size();
      if (nfields == 0) continue;

      std::vector<Field> topo_fields;
      for (int i = 0; i < nfields; ++i) {
        auto eamxx_fname = topography_eamxx_fields_names[grid_name][i];
        auto file_fname  = topography_file_fields_names[grid_name][i];
        topo_fields.push_back(
            m_field_mgr->get_field(eamxx_fname, grid_name).alias(file_fname));
      }

      if (not m_iop_data_mgr) {
        // Topography files always use "ncol_d" for the GLL ncol.
        strmap_t tag_rename;
        if (grid_name == "physics_gll") {
          tag_rename["ncol"] = "ncol_d";
        }
        FieldReader reader;
        reader.set_file_specs(file_name, tag_rename);
        reader.set_dim_decomp(grid->get_partitioned_dim_gids(), m_comm);
        reader.set_fields(topo_fields);
        reader.read();
      } else {
        m_iop_data_mgr->read_fields_from_file_for_iop(
            file_name, topo_fields, grid);
      }
      for (auto& f : topo_fields) {
        f.get_header().get_tracking().update_time_stamp(m_current_ts);
      }
    }
    m_atm_params.sublist("provenance").set("topography_file", file_name);
    m_logger->debug("    [EAMxx] Processing topography from file ... done!");
  } else {
    for (const auto& grid_name : m_grids_mgr->get_grid_names()) {
      EKAT_REQUIRE_MSG(topography_file_fields_names[grid_name].size() == 0,
                       "Error! Topography data was requested in the FM, but no "
                       "topography_filename or entry matching the field name "
                       "was given in IC parameters.\n");
    }
    m_atm_params.sublist("provenance").set<std::string>("topography_file", "NONE");
  }

  if (m_iop_data_mgr) {
    // Load IOP data for the initial time stamp.
    m_iop_data_mgr->read_iop_file_data(m_current_ts);

    // Set fields from IOP data (ICs are on GLL grid; dynamics handles PG2).
    if (m_field_mgr->get_grids_manager()->get_grid_names().count("physics_gll") > 0) {
      m_iop_data_mgr->set_fields_from_iop_data(m_field_mgr, "physics_gll");
    }
  }

  // Compute IC perturbations of GLL fields (if requested).
  const auto perturbed_fields =
      ic_pl.get<strvec_t>("perturbed_fields", {});
  const auto num_perturb_fields = perturbed_fields.size();
  if (num_perturb_fields > 0) {
    m_logger->info("    [EAMxx] Adding random perturbation to ICs ...");

    EKAT_REQUIRE_MSG(
        m_field_mgr->get_grids_manager()->get_grid_names().count("physics_gll") > 0,
        "Error! Random perturbation can only be applied to fields on "
        "the GLL grid, but no physics GLL grid was defined in FieldManager.\n");

    int seed;
    if (ic_pl.get<bool>("generate_perturbation_random_seed", false)) {
      EKAT_REQUIRE_MSG(
          not ic_pl.isParameter("perturbation_random_seed"),
          "Error! Param generate_perturbation_random_seed=true, and "
          "a perturbation_random_seed is given. Only one of these can "
          "be defined for a simulation.\n");
      std::srand(std::time(nullptr));
      seed = std::rand();
    } else {
      seed = ic_pl.get<int>("perturbation_random_seed", 0);
    }
    m_logger->info("      For IC perturbation, random seed: " +
                   std::to_string(seed));

    const auto perturbation_limit =
        ic_pl.get<Real>("perturbation_limit", 0.001);

    const auto gll_grid = m_grids_mgr->get_grid("physics_gll");
    const auto hyam_h =
        gll_grid->get_geometry_data("hyam").get_view<const Real*, Host>();
    const auto hybm_h =
        gll_grid->get_geometry_data("hybm").get_view<const Real*, Host>();
    constexpr auto ps0 = physics::Constants<Real>::P0.value;
    const auto min_pressure =
        ic_pl.get<Real>("perturbation_minimum_pressure", 1050.0);

    using namespace ShortFieldTagsNames;
    const auto& pmask_lt = gll_grid->get_vertical_layout(LEV);
    const auto  nondim   = ekat::units::none;
    FieldIdentifier pmask_fid("lev_mask", pmask_lt, nondim,
                               gll_grid->name(), DataType::IntType);
    Field pressure_mask(pmask_fid, true);
    auto  pmask_h = pressure_mask.get_view<int*, Host>();
    for (int ilev = 0; ilev < pmask_lt.dim(0); ++ilev) {
      const auto pref =
          (hyam_h(ilev) * ps0 + hybm_h(ilev) * ps0) / 100;
      pmask_h(ilev) = static_cast<int>(pref > min_pressure);
    }
    pressure_mask.sync_to_dev();

    const std::string gll_grid_name = gll_grid->name();
    auto dofs_gids = gll_grid->get_dofs_gids();
    for (size_t f = 0; f < perturbed_fields.size(); ++f) {
      const auto fname = perturbed_fields[f];
      EKAT_REQUIRE_MSG(
          ekat::contains(m_fields_inited[gll_grid_name], fname),
          "Error! Attempting to apply perturbation to field not in "
          "initial_conditions.\n"
          "  - Field: " + fname + "\n"
          "  - Grid:  " + gll_grid_name + "\n");

      auto field = m_field_mgr->get_field(fname, gll_grid_name);
      perturb(field, perturbation_limit, seed, pressure_mask, dofs_gids);
    }

    m_logger->info("    [EAMxx] Adding random perturbation to ICs ... done!");
  }

  m_logger->flush();
}

void DefaultModelInit::
initialize_constant_field (const FieldIdentifier& fid,
                            const ekat::ParameterList& ic_pl)
{
  auto f = m_field_mgr->get_field(fid);
  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto& name   = fid.name();

  if (layout.is_vector_layout() and
      ic_pl.isType<std::vector<double>>(name)) {
    const int idim    = layout.get_vector_component_idx();
    const int vec_dim = layout.get_vector_dim();
    const auto& values = ic_pl.get<std::vector<double>>(name);
    EKAT_REQUIRE_MSG (values.size() == static_cast<size_t>(vec_dim),
        "Error! Initial condition values array for '" + name +
        "' has the wrong dimension.\n"
        "       Field dimension: " + std::to_string(vec_dim) + "\n"
        "       Array dimenions: " + std::to_string(values.size()) + "\n");

    if (layout.rank() == 2 and idim == 1) {
      using kt = Field::kt_dev;
      typename kt::view_1d<double> data("data", vec_dim);
      auto data_h = Kokkos::create_mirror_view(data);
      for (int i = 0; i < vec_dim; ++i) {
        data_h(i) = values[i];
      }
      Kokkos::deep_copy(data, data_h);

      const int n = layout.dim(0);
      auto v = f.get_view<double**>();
      Kokkos::parallel_for(
          typename kt::RangePolicy(0, n),
          KOKKOS_LAMBDA(const int i) {
            for (int j = 0; j < vec_dim; ++j) {
              v(i, j) = data(j);
            }
          });
    } else {
      for (int comp = 0; comp < vec_dim; ++comp) {
        auto f_i = f.get_component(comp);
        f_i.deep_copy(values[comp]);
      }
    }
  } else {
    if (ic_pl.isType<int>(name)) {
      f.deep_copy(ic_pl.get<int>(name));
    } else {
      f.deep_copy(ic_pl.get<double>(name));
    }
  }
}

std::shared_ptr<ModelInit>
create_default_model_init (const ekat::ParameterList& params)
{
  return std::make_shared<DefaultModelInit>(params);
}

void register_default_model_init ()
{
  static bool registered = false;
  if (not registered) {
    auto& f = ModelInitFactory::instance();
    f.register_product("default", &create_default_model_init);
    // "homme" uses DefaultModelInit until a dedicated HommeModelInit exists.
    f.register_product("homme",   &create_default_model_init);
    registered = true;
  }
}

} // namespace scream
