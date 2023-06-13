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
  double time_value_1= scorpio::read_time_at_index_c2f(m_datafiles[0].c_str(),1);
  double time_value_2= scorpio::read_time_at_index_c2f(m_datafiles[0].c_str(),2);

  // Here we are assuming that the time in the netcdf file is in days
  // Internally we want this in seconds so need to convert
  // Only consider integer time steps in seconds to resolve any roundoff error
  m_time_step_file=int((time_value_2-time_value_1)*86400);
  scorpio::eam_pio_closefile(m_datafiles[0]);
}
// =========================================================================================
void Nudging::apply_tendency(Field& base, const Field& next, const int dt)
{
  // Calculate the weight to apply the tendency
  const Real dtend = Real(dt)/Real(m_timescale);
  EKAT_REQUIRE_MSG(dtend>=0,"Error! Nudging::apply_tendency - timescale tendency of " << std::to_string(dt) << " / " << std::to_string(m_timescale) << " = " << std::to_string(dtend) << " is invalid.  Please check the timescale and/or dt");
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

  std::map<std::string,view_1d_host<Real>> host_views;
  std::map<std::string,FieldLayout>  layouts;

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_src_levs} };

  constexpr int ps = 1;  // TODO: I think this could be the regular packsize, right?
  const auto& grid_name = m_grid->name();
  for (auto name : m_fields_nudge) {
    // Helper fields that will temporarily store the target state, which can then
    // be used to back out a nudging tendency
    auto field  = get_field_out(name);
    auto layout = field.get_header().get_identifier().get_layout();
    create_helper_field(name, layout, grid_name, ps);
    // Extended fields for handling data from file.
    m_fields_ext[name]   = view_2d<Real>(name,m_num_cols,m_num_src_levs);
    m_fields_ext_h[name] = Kokkos::create_mirror_view(m_fields_ext[name]);
    auto view_h          = m_fields_ext_h[name];
    host_views[name]     = view_1d_host<Real>(view_h.data(),view_h.size());
    layouts.emplace(name, scalar3d_layout_mid);
  }

  m_fields_ext["p_mid"] = view_2d<Real>("p_mid",m_num_cols,m_num_src_levs);
  m_fields_ext_h["p_mid"]       = Kokkos::create_mirror_view(m_fields_ext["p_mid"]);
  auto p_mid_h=m_fields_ext_h["p_mid"];
  host_views["p_mid"] = view_1d_host<Real>(p_mid_h.data(),p_mid_h.size());
  layouts.emplace("p_mid", scalar3d_layout_mid);

  auto grid_l = m_grid->clone("Point Grid", false);
  grid_l->reset_num_vertical_lev(m_num_src_levs);
  ekat::ParameterList data_in_params;
  std::vector<std::string> fnames = {"T_mid","p_mid","qv","u","v"};
  data_in_params.set("Field Names",fnames);
  data_in_params.set("Filename",m_datafiles[0]);
  // We need to skip grid checks because multiple ranks
  // may want the same column of source data.
  data_in_params.set("Skip_Grid_Checks",true);
  m_data_input.init(data_in_params,grid_l,host_views,layouts);

  m_ts0=timestamp();

  //Check that internal timestamp starts at same point as time in external file
  auto case_t0 = scorpio::read_timestamp(m_datafiles[0],"case_t0");

  EKAT_REQUIRE_MSG(case_t0.get_year()==m_ts0.get_year(),
		   "ERROR: The start year from the nudging file is "\
		   "different than the internal simulation start year\n");
  EKAT_REQUIRE_MSG(case_t0.get_month()==m_ts0.get_month(),
		   "ERROR: The start month from the nudging file is "\
		   "different than the internal simulation start month\n");
  EKAT_REQUIRE_MSG(case_t0.get_day()==m_ts0.get_day(),
		   "ERROR: The start day from the nudging file is "\
		   "different than the internal simulation start day\n");
  EKAT_REQUIRE_MSG(case_t0.get_hours()==m_ts0.get_hours(),
		   "ERROR: The start hour from the nudging file is "\
		   "different than the internal simulation start hour\n");
  EKAT_REQUIRE_MSG(case_t0.get_minutes()==m_ts0.get_minutes(),
		   "ERROR: The start minute from the nudging file is "\
		   "different than the internal simulation start minute\n");
  EKAT_REQUIRE_MSG(case_t0.get_seconds()==m_ts0.get_seconds(),
		   "ERROR: The start second from the nudging file is "\
		   "different than the internal simulation start second\n");

  //Initialize before and after data
  m_NudgingData_bef.init(m_num_cols,m_num_src_levs,true);
  m_NudgingData_bef.time = -999;
  m_NudgingData_aft.init(m_num_cols,m_num_src_levs,true);
  m_NudgingData_aft.time = -999;

  //Read in the first time step
  m_data_input.read_variables(0);
  Kokkos::deep_copy(m_NudgingData_bef.T_mid,m_fields_ext_h["T_mid"]);
  Kokkos::deep_copy(m_NudgingData_bef.p_mid,m_fields_ext_h["p_mid"]);
  Kokkos::deep_copy(m_NudgingData_bef.u,m_fields_ext_h["u"]);
  Kokkos::deep_copy(m_NudgingData_bef.v,m_fields_ext_h["v"]);
  Kokkos::deep_copy(m_NudgingData_bef.qv,m_fields_ext_h["qv"]);
  m_NudgingData_bef.time = 0.;

  //Read in the first time step
  m_data_input.read_variables(1);
  Kokkos::deep_copy(m_NudgingData_aft.T_mid,m_fields_ext_h["T_mid"]);
  Kokkos::deep_copy(m_NudgingData_aft.p_mid,m_fields_ext_h["p_mid"]);
  Kokkos::deep_copy(m_NudgingData_aft.u,m_fields_ext_h["u"]);
  Kokkos::deep_copy(m_NudgingData_aft.v,m_fields_ext_h["v"]);
  Kokkos::deep_copy(m_NudgingData_aft.qv,m_fields_ext_h["qv"]);
  m_NudgingData_aft.time = m_time_step_file;

}

void Nudging::time_interpolation (const int time_s) {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using MemberType = typename KT::MemberType;

  const int time_index = time_s/m_time_step_file;
  double time_step_file_d = m_time_step_file;
  double w_bef = ((time_index+1)*m_time_step_file-time_s) / time_step_file_d;
  double w_aft = (time_s-(time_index)*m_time_step_file) / time_step_file_d;
  const int num_cols = m_NudgingData_aft.T_mid.extent(0);
  const int num_vert_packs = m_NudgingData_aft.T_mid.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);

  auto T_mid_ext = m_fields_ext["T_mid"];
  auto p_mid_ext = m_fields_ext["p_mid"];
  auto qv_ext = m_fields_ext["qv"];
  auto u_ext = m_fields_ext["u"];
  auto v_ext = m_fields_ext["v"];

  auto T_mid_bef = m_NudgingData_bef.T_mid;
  auto T_mid_aft = m_NudgingData_aft.T_mid;
  auto u_bef = m_NudgingData_bef.u;
  auto u_aft = m_NudgingData_aft.u;
  auto v_bef = m_NudgingData_bef.v;
  auto v_aft = m_NudgingData_aft.v;
  auto p_mid_bef = m_NudgingData_bef.p_mid;
  auto p_mid_aft = m_NudgingData_aft.p_mid;
  auto qv_bef = m_NudgingData_bef.qv;
  auto qv_aft = m_NudgingData_aft.qv;

  Kokkos::parallel_for("nudging_time_interpolation", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();

      auto T_mid_1d = ekat::subview(T_mid_ext,icol);
      auto T_mid_bef_1d = ekat::subview(T_mid_bef,icol);
      auto T_mid_aft_1d = ekat::subview(T_mid_aft,icol);
      auto u_1d = ekat::subview(u_ext,icol);
      auto u_bef_1d = ekat::subview(u_bef,icol);
      auto u_aft_1d = ekat::subview(u_aft,icol);
      auto v_1d = ekat::subview(v_ext,icol);
      auto v_bef_1d = ekat::subview(v_bef,icol);
      auto v_aft_1d = ekat::subview(v_aft,icol);
      auto p_mid_1d = ekat::subview(p_mid_ext,icol);
      auto p_mid_bef_1d = ekat::subview(p_mid_bef,icol);
      auto p_mid_aft_1d = ekat::subview(p_mid_aft,icol);
      auto qv_1d = ekat::subview(qv_ext,icol);
      auto qv_bef_1d = ekat::subview(qv_bef,icol);
      auto qv_aft_1d = ekat::subview(qv_aft,icol);

      const auto range = Kokkos::TeamThreadRange(team, num_vert_packs);
      Kokkos::parallel_for(range, [&] (const Int & k) {
        T_mid_1d(k)=w_bef*T_mid_bef_1d(k) + w_aft*T_mid_aft_1d(k);
        p_mid_1d(k)=w_bef*p_mid_bef_1d(k) + w_aft*p_mid_aft_1d(k);
        u_1d(k)=w_bef*u_bef_1d(k) + w_aft*u_aft_1d(k);
        v_1d(k)=w_bef*v_bef_1d(k) + w_aft*v_aft_1d(k);
        qv_1d(k)=w_bef*qv_bef_1d(k) + w_aft*qv_aft_1d(k);
      });
      team.team_barrier();
  });
  Kokkos::fence();
}

// =========================================================================================
void Nudging::update_time_step (const int time_s)
{
  //Check to see in time state needs to be updated
  if (time_s >= m_NudgingData_aft.time)
    {
      const int time_index = time_s/m_time_step_file;
      std::swap (m_NudgingData_bef,m_NudgingData_aft);
      m_NudgingData_bef.time = m_NudgingData_aft.time;

      m_data_input.read_variables(time_index+1);
      Kokkos::deep_copy(m_NudgingData_aft.T_mid,m_fields_ext_h["T_mid"]);
      Kokkos::deep_copy(m_NudgingData_aft.p_mid,m_fields_ext_h["p_mid"]);
      Kokkos::deep_copy(m_NudgingData_aft.u,m_fields_ext_h["u"]);
      Kokkos::deep_copy(m_NudgingData_aft.v,m_fields_ext_h["v"]);
      Kokkos::deep_copy(m_NudgingData_aft.qv,m_fields_ext_h["qv"]);
      m_NudgingData_aft.time = m_time_step_file*(time_index+1);
    }

}


// =========================================================================================
void Nudging::run_impl (const double dt)
{
  using namespace scream::vinterp;

  //Have to add dt because first time iteration is at 0 seconds where you will
  //not have any data from the field. The timestamp is only iterated at the
  //end of the full step in scream.
  auto ts = timestamp()+dt;
  auto time_since_zero = ts.seconds_from(m_ts0);

  //update time state information and check whether need to update data
  update_time_step(time_since_zero);

  //perform time interpolation
  time_interpolation(time_since_zero);

  // Process data and nudge the atmosphere state
  const auto& p_mid_v     = get_field_in("p_mid").get_view<const mPack**>();
  const auto& p_mid_ext   = m_fields_ext.at("p_mid");
  const view_Nd<mPack,2> p_mid_ext_p(reinterpret_cast<mPack*>(p_mid_ext.data()),
				     m_num_cols,m_num_src_levs);
  for (auto name : m_fields_nudge) {
    auto atm_state_field = get_field_out(name);
    auto int_state_field = get_helper_field(name);
    auto ext_state_field = m_fields_ext.at(name);
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
  m_data_input.finalize();
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
