#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/p3_main_impl.hpp"

#include "ekat/ekat_assert.hpp"

#include <array>

namespace scream
{
/*
 * P3 Microphysics routines
*/

  using namespace p3;
  using P3F          = Functions<Real, DefaultDevice>;
  using Spack        = typename P3F::Spack;
  using Pack         = ekat::Pack<Real,Spack::n>;

// =========================================================================================
P3Microphysics::P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_p3_comm (comm)
 , m_p3_params (params)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, p3 options.
*/
  m_initializer = create_field_initializer<P3InputsInitializer>();
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  auto nondim = m/m;
  auto mm = m/1000;

  // Retrieve physical dimension extents
  m_num_levs = SCREAM_NUM_VERTICAL_LEV;  // Number of levels per column

  const auto& grid_name = m_p3_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this ranks

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  FieldLayout scalar3d_layout_mid { {COL,VL}, {m_num_cols,m_num_levs} };   // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level 
  FieldLayout scalar3d_layout_int { {COL,VL}, {m_num_cols,m_num_levs+1} }; // Layout for 3D (2d horiz X 1d vertical) variable defined at the interfaces

  // Define fields needed in P3.
  // Note: p3_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  // These variables are needed, but not actually passed to p3_main.  For organization
  // purposes they are listed here:
  m_required_fields.emplace("ast",   scalar3d_layout_mid, nondim, grid_name);
  m_required_fields.emplace("pmid",  scalar3d_layout_mid, Pa,     grid_name);
  m_required_fields.emplace("zi",    scalar3d_layout_int, m,      grid_name);
  m_required_fields.emplace("T_atm", scalar3d_layout_mid, K,      grid_name);
  m_computed_fields.emplace("T_atm", scalar3d_layout_mid, K,      grid_name);  // T_atm is the only one of these variables that is also updated.

  // Prognostic State:  (all fields are both input and output)
  m_required_fields.emplace("qv",     scalar3d_layout_mid, Q,    grid_name);
  m_required_fields.emplace("qc",     scalar3d_layout_mid, Q,    grid_name);
  m_required_fields.emplace("qr",     scalar3d_layout_mid, Q,    grid_name);
  m_required_fields.emplace("qi",     scalar3d_layout_mid, Q,    grid_name);
  m_required_fields.emplace("qm",     scalar3d_layout_mid, Q,    grid_name);
  m_required_fields.emplace("nc",     scalar3d_layout_mid, 1/kg, grid_name);
  m_required_fields.emplace("nr",     scalar3d_layout_mid, 1/kg, grid_name);
  m_required_fields.emplace("ni",     scalar3d_layout_mid, 1/kg, grid_name);
  m_required_fields.emplace("bm",     scalar3d_layout_mid, 1/kg, grid_name);
  m_required_fields.emplace("th_atm", scalar3d_layout_mid, K,    grid_name);  //TODO: Delete, don't acutally need this as required.
  //
  m_computed_fields.emplace("qv",     scalar3d_layout_mid, Q,    grid_name);
  m_computed_fields.emplace("qc",     scalar3d_layout_mid, Q,    grid_name);
  m_computed_fields.emplace("qr",     scalar3d_layout_mid, Q,    grid_name);
  m_computed_fields.emplace("qi",     scalar3d_layout_mid, Q,    grid_name);
  m_computed_fields.emplace("qm",     scalar3d_layout_mid, Q,    grid_name);
  m_computed_fields.emplace("nc",     scalar3d_layout_mid, 1/kg, grid_name);
  m_computed_fields.emplace("nr",     scalar3d_layout_mid, 1/kg, grid_name);
  m_computed_fields.emplace("ni",     scalar3d_layout_mid, 1/kg, grid_name);
  m_computed_fields.emplace("bm",     scalar3d_layout_mid, 1/kg, grid_name);
  m_computed_fields.emplace("th_atm", scalar3d_layout_mid, K,    grid_name);
  // Diagnostic Inputs: (only the X_prev fields are both input and output, all others are just inputs)
  m_required_fields.emplace("nc_nuceat_tend",  scalar3d_layout_mid, 1/(kg*s), grid_name);
  m_required_fields.emplace("nccn_prescribed", scalar3d_layout_mid, nondim,   grid_name);
  m_required_fields.emplace("ni_activated",    scalar3d_layout_mid, 1/kg,     grid_name);
  m_required_fields.emplace("inv_qc_relvar",   scalar3d_layout_mid, nondim,   grid_name);
  m_required_fields.emplace("dp",              scalar3d_layout_mid, Pa,       grid_name);
  m_required_fields.emplace("qv_prev",         scalar3d_layout_mid, Q,        grid_name);
  m_required_fields.emplace("T_prev",          scalar3d_layout_mid, K,        grid_name); 
  //
  m_computed_fields.emplace("qv_prev",         scalar3d_layout_mid, Q,        grid_name);
  m_computed_fields.emplace("T_prev",          scalar3d_layout_mid, K,        grid_name);
  // Diagnostic Outputs: (all fields are just outputs w.r.t. P3)
  m_computed_fields.emplace("mu_c",               scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("lamc",               scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("diag_eff_radius_qc", scalar3d_layout_mid, m,      grid_name);
  m_computed_fields.emplace("diag_eff_radius_qi", scalar3d_layout_mid, m,      grid_name);
  m_computed_fields.emplace("precip_total_tend",  scalar3d_layout_mid, mm,     grid_name);
  m_computed_fields.emplace("nevapr",             scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("qr_evap_tend",       scalar3d_layout_mid, mm/s,   grid_name);
  // History Only: (all fields are just outputs and are really only meant for I/O purposes)
  m_computed_fields.emplace("liq_ice_exchange", scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("vap_liq_exchange", scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("vap_ice_exchange", scalar3d_layout_mid, nondim, grid_name);
  // TODO: AaronDonahue - The following should actually be set locally and therefore don't
  //       need to be fields registered with the field manager.  Currently we are use them
  //       as fields so they can be initialized by the p3_inputs_initialzer for testing
  //       purposes.  A future task is to remove these and define them locally as 2d views
  //       in run_impl.
  m_required_fields.emplace("dz", scalar3d_layout_mid, K,         grid_name);
  m_required_fields.emplace("exner", scalar3d_layout_mid, K,      grid_name);
  m_required_fields.emplace("cld_frac_l", scalar3d_layout_mid, K, grid_name);
  m_required_fields.emplace("cld_frac_i", scalar3d_layout_mid, K, grid_name);
  m_required_fields.emplace("cld_frac_r", scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("dz", scalar3d_layout_mid, K,         grid_name);
  m_computed_fields.emplace("exner", scalar3d_layout_mid, K,      grid_name);
  m_computed_fields.emplace("cld_frac_l", scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("cld_frac_i", scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("cld_frac_r", scalar3d_layout_mid, K, grid_name);

}

// =========================================================================================
void P3Microphysics::initialize_impl (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // We may have to init some fields from within P3. This can be the case in a P3 standalone run.
  // Some options:
  //  - we can tell P3 it can init all inputs or specify which ones it can init. We call the
  //    resulting list of inputs the 'initializaable' (or initable) inputs. The default is
  //    that no inputs can be inited.
  //  - we can request that P3 either inits no inputs or all of the initable ones (as specified
  //    at the previous point). The default is that P3 must be in charge of init ing ALL or NONE
  //    of its initable inputs.
  // Recall that:
  //  - initable fields may not need initialization (e.g., some other atm proc that
  //    appears earlier in the atm dag might provide them).

  std::vector<std::string> p3_inputs = {"T_atm","ast","ni_activated","nc_nuceat_tend","pmid","dp","zi","qv_prev","T_prev",
                                        "qv", "qc", "qr", "qi", "qm", "nc", "nr", "ni", "bm","nccn_prescribed","inv_qc_relvar"
                                       ,"th_atm","exner","dz","cld_frac_l","cld_frac_r","cld_frac_i"  // TODO: Delete these, should be local.
                                       };
  using strvec = std::vector<std::string>;
  const strvec& allowed_to_init = m_p3_params.get<strvec>("Initializable Inputs",strvec(0));
  const bool can_init_all = m_p3_params.get<bool>("Can Initialize All Inputs", false);
  const bool init_all_or_none = m_p3_params.get<bool>("Must Init All Inputs Or None", true);

  const strvec& initable = can_init_all ? p3_inputs : allowed_to_init;
  if (initable.size()>0) {
    bool all_inited = true, all_uninited = true;
    for (const auto& name : initable) {
      const auto& f = m_p3_fields_in.at(name);
      const auto& track = f.get_header().get_tracking();
      if (track.get_init_type()==InitType::None) {
        // Nobody claimed to init this field. P3InputsInitializer will take care of it
        m_initializer->add_me_as_initializer(f);
        all_uninited &= true;
        all_inited &= false;
      } else {
        all_uninited &= false;
        all_inited &= true;
      }
    }

    // In order to gurantee some consistency between inputs, it is best if P3
    // initializes either none or all of the inputs.
    EKAT_REQUIRE_MSG (!init_all_or_none || all_inited || all_uninited,
                      "Error! Some p3 inputs were marked to be inited by P3, while others weren't.\n"
                      "       P3 was requested to init either all or none of the inputs.\n");
  }
}

// =========================================================================================
void P3Microphysics::run_impl (const Real dt)
{
  // std::array<const char*, num_views> view_names = {"q", "FQ", "T", "zi", "pmid", "dpres", "ast", "ni_activated", "nc_nuceat_tend"};

  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_p3_fields_in) {
    Kokkos::deep_copy(m_p3_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(m_p3_host_views_out.at(it.first),it.second.get_view());
  }

  // Copy outputs back to device
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(it.second.get_view(),m_p3_host_views_out.at(it.first));
  }

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it, updating the p3 fields.
  auto ts = timestamp();
  ts += dt;
//  m_p3_fields_out.at("q").get_header().get_tracking().update_time_stamp(ts);
//  m_p3_fields_out.at("FQ").get_header().get_tracking().update_time_stamp(ts);
//  m_p3_fields_out.at("T").get_header().get_tracking().update_time_stamp(ts);
}

// =========================================================================================
void P3Microphysics::finalize_impl()
{
  // Do nothing
}

// =========================================================================================
void P3Microphysics::register_fields (FieldRepository<Real>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field<Pack>(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field<Pack>(fid);
  }
}

void P3Microphysics::set_required_field_impl (const Field<const Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_in.emplace(name,f);
  m_p3_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_p3_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void P3Microphysics::set_computed_field_impl (const Field<      Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_out.emplace(name,f);
  m_p3_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_p3_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
