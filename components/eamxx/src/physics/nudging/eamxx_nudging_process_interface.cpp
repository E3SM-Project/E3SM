#include "eamxx_nudging_process_interface.hpp"

namespace scream
{

  //using namespace spa;
// =========================================================================================
Nudging::Nudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_datafiles  = m_params.get<std::vector<std::string>>("nudging_filename");
  m_timescale = m_params.get<int>("nudging_timescale",0);
  m_fields_nudge = m_params.get<std::vector<std::string>>("nudging_fields");
  // TODO: Add some warning messages here.
  // 1. if m_timescale is <= 0 we will do direct replacement.
  // 2. if m_fields_nudge is empty or =NONE then we skip nudging altogether.
}

// =========================================================================================
void Nudging::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };

  constexpr int ps = 1;
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_nudge variable
   */
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("qv",    scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  add_field<Updated>("u",     scalar3d_layout_mid, m/s, grid_name, ps);
  add_field<Updated>("v",     scalar3d_layout_mid, m/s, grid_name, ps);
  /* ----------------------- WARNING --------------------------------*/

  //Now need to read in the file
  scorpio::register_file(m_datafiles[0],scorpio::Read);
  m_num_src_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
  scorpio::eam_pio_closefile(m_datafiles[0]);
}
// =========================================================================================
void Nudging::apply_tendency(Field& base, const Field& next, const int dt)
{
  // Calculate the weight to apply the tendency
  const Real dtend = Real(dt)/Real(m_timescale);
  EKAT_REQUIRE_MSG(dtend>=0,"Error! Nudging::apply_tendency - timescale tendency of " << std::to_string(dt) 
		  << " / " << std::to_string(m_timescale) << " = " << std::to_string(dtend) 
		  << " is invalid.  Please check the timescale and/or dt");
  // Now apply the tendency.
  Field tend = base.clone();
  // Use update internal to set tendency, will be (1.0*next - 1.0*base), note tend=base at this point.
  tend.update(next,1.0,-1.0);
  base.update(tend,dtend,1.0);
}
// =========================================================================================
void Nudging::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;

  // Initialize the time interpolator
  auto grid_ext = m_grid->clone("Point Grid", false);
  grid_ext->reset_num_vertical_lev(m_num_src_levs);
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  m_time_interp = util::TimeInterpolation(grid_ext, m_datafiles);

  constexpr int ps = 1;  // TODO: I think this could be the regular packsize, right?
  const auto& grid_name = m_grid->name();
  create_helper_field("p_mid_ext", scalar3d_layout_mid, grid_name, ps);
  auto pmid_ext = get_helper_field("p_mid_ext");
  m_time_interp.add_field(pmid_ext,"p_mid",true);
  for (auto name : m_fields_nudge) {
    std::string name_ext = name + "_ext";
    // Helper fields that will temporarily store the target state, which can then
    // be used to back out a nudging tendency
    auto field  = get_field_out(name);
    auto layout = field.get_header().get_identifier().get_layout();
    create_helper_field(name,     layout,              grid_name, ps);
    create_helper_field(name_ext, scalar3d_layout_mid, grid_name, ps);
    auto field_ext = get_helper_field(name_ext);
    m_time_interp.add_field(field_ext.alias(name),true);
  }
  m_time_interp.initialize_data_from_files();

}

// =========================================================================================
void Nudging::run_impl (const double dt)
{
  using namespace scream::vinterp;

  //Have to add dt because first time iteration is at 0 seconds where you will
  //not have any data from the field. The timestamp is only iterated at the
  //end of the full step in scream.
  auto ts = timestamp()+dt;

  //perform time interpolation
  m_time_interp.perform_time_interpolation(ts);

  // Process data and nudge the atmosphere state
  const auto& p_mid_v     = get_field_in("p_mid").get_view<const mPack**>();
  const auto& p_mid_ext   = get_helper_field("p_mid_ext").get_view<mPack**>();
  const view_Nd<mPack,2> p_mid_ext_p(reinterpret_cast<mPack*>(p_mid_ext.data()),
				     m_num_cols,m_num_src_levs);
  for (auto name : m_fields_nudge) {
    auto atm_state_field = get_field_out(name);
    auto int_state_field = get_helper_field(name);
    auto ext_state_field = get_helper_field(name+"_ext").get_view<mPack**>();
    auto atm_state_view  = atm_state_field.get_view<mPack**>();  // TODO: Right now assume whatever field is defined on COLxLEV
    auto int_state_view  = int_state_field.get_view<mPack**>();
    const view_Nd<mPack,2> ext_state_view(reinterpret_cast<mPack*>(ext_state_field.data()),
                                          m_num_cols,m_num_src_levs);
    // Vertical Interpolation onto atmosphere state pressure levels
    perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                             p_mid_v,
                                             ext_state_view,
                                             int_state_view,
                                             m_num_src_levs,
                                             m_num_levs);
    if (m_timescale <= 0) {
      // We do direct replacement
      Kokkos::deep_copy(atm_state_view,int_state_view);
    } else {
      // Back out a tendency and apply it.
      apply_tendency(atm_state_field, int_state_field, dt);
    }

    // Remove mask values, setting them to the nearest unmasked value in the column.
    // TODO: This chunk of code has a bug in it where we are assigning values by PACK.
    // Note, that in our case the packsize is hard-coded to 1, so it shouldn't actually
    // be a problem, but the way the code is set up here could lead to confusion if we
    // ever change that hard-coded value.  We will need to replace this with something
    // better, in a subsequent commit.
    const int num_cols           = atm_state_view.extent(0);
    const int num_vert_packs     = atm_state_view.extent(1);
    const int num_vert_packs_ext = ext_state_view.extent(1);
    const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
        const int icol = team.league_rank();
        auto atm_state_view_1d = ekat::subview(atm_state_view,icol);
        auto ext_state_view_1d = ekat::subview(ext_state_view,icol);
        auto p_mid_1d = ekat::subview(p_mid_v,icol);
        auto p_ext_1d = ekat::subview(p_mid_ext,icol);
        const auto range = Kokkos::TeamThreadRange(team, num_vert_packs);
        Kokkos::parallel_for(range, [&] (const Int & k) {
          const auto above_max = p_mid_1d(k) > p_ext_1d(num_vert_packs_ext-1);
          const auto below_min = p_mid_1d(k) < p_ext_1d(0);
  	  if (above_max.any()){
  	    atm_state_view_1d(k).set(above_max, ext_state_view_1d(num_vert_packs_ext-1));
	  }
	  if (below_min.any()){
	    atm_state_view_1d(k).set(below_min,ext_state_view_1d(0));
  	  }
        });
        team.team_barrier();
    });
    Kokkos::fence();
  }
}

// =========================================================================================
void Nudging::finalize_impl()
{
  m_time_interp.finalize();
}
// =========================================================================================
void Nudging::create_helper_field (const std::string& name,
                                             const FieldLayout& layout,
                                             const std::string& grid_name,
                                             const int ps)
{
  using namespace ekat::units;
  // For helper fields we don't bother w/ units, so we set them to non-dimensional
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  if (ps>=0) {
    f.get_header().get_alloc_properties().request_allocation(ps);
  } else {
    f.get_header().get_alloc_properties().request_allocation();
  }
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}

} // namespace scream
