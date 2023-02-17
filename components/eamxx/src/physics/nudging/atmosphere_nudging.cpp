#include "atmosphere_nudging.hpp"

namespace scream
{

  //using namespace spa;
// =========================================================================================
Nudging::Nudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  datafile=m_params.get<std::string>("Nudging_Filename");
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
  // Layout for horiz_wind field
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };

  constexpr int ps = 1;
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("qv", scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  add_field<Updated>("horiz_winds",   horiz_wind_layout,   m/s,     grid_name, ps);
  
  //Now need to read in the file
  scorpio::register_file(datafile,scorpio::Read);
  m_num_src_levs = scorpio::get_dimlen_c2f(datafile.c_str(),"lev");
  double time_value_1= scorpio::read_time_at_index_c2f(datafile.c_str(),1);
  double time_value_2= scorpio::read_time_at_index_c2f(datafile.c_str(),2);
    
  //Here we are assuming that the time in the netcdf file is in days
  //Internally we want this in seconds so need to convert
  //Only consider integer time steps in seconds to resolve any roundoff error
  time_step_file=int((time_value_2-time_value_1)*86400);
  scorpio::eam_pio_closefile(datafile);
}

// =========================================================================================
void Nudging::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_src_levs} };
  fields_ext["T_mid"] = view_2d<Real>("T_mid",m_num_cols,m_num_src_levs);
  auto T_mid_h       = Kokkos::create_mirror_view(fields_ext["T_mid"]);
  host_views["T_mid"] = view_1d_host<Real>(T_mid_h.data(),T_mid_h.size());
  layouts.emplace("T_mid", scalar3d_layout_mid);

  fields_ext["p_mid"] = view_2d<Real>("p_mid",m_num_cols,m_num_src_levs);
  auto p_mid_h       = Kokkos::create_mirror_view(fields_ext["p_mid"]);
  host_views["p_mid"] = view_1d_host<Real>(p_mid_h.data(),p_mid_h.size());
  layouts.emplace("p_mid", scalar3d_layout_mid);
  
  fields_ext["qv"] = view_2d<Real>("qv",m_num_cols,m_num_src_levs);
  auto qv_h       = Kokkos::create_mirror_view(fields_ext["qv"]);
  host_views["qv"] = view_1d_host<Real>(qv_h.data(),qv_h.size());
  layouts.emplace("qv", scalar3d_layout_mid);

  fields_ext_3d["horiz_winds"] = view_3d<Real>("horiz_winds",m_num_cols,2,m_num_src_levs);
  auto horiz_winds_h       = Kokkos::create_mirror_view(fields_ext_3d["horiz_winds"]);
  host_views["horiz_winds"] = view_1d_host<Real>(horiz_winds_h.data(),horiz_winds_h.size());
  layouts.emplace("horiz_winds", horiz_wind_layout);
  
  auto grid_l = m_grid->clone("Point Grid", false);
  grid_l->reset_num_vertical_lev(m_num_src_levs);
  ekat::ParameterList data_in_params;
  m_fnames = {"T_mid","p_mid","qv","horiz_winds"};
  data_in_params.set("Field Names",m_fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input0(m_comm,data_in_params);
  data_input0.init(grid_l,host_views,layouts);
  data_input=std::make_shared<AtmosphereInput>(data_input0);

  T_mid_ext = fields_ext["T_mid"];
  p_mid_ext = fields_ext["p_mid"];
  horiz_winds_ext = fields_ext_3d["horiz_winds"];
  qv_ext = fields_ext["qv"];
  ts0=timestamp();

  //Check that internal timestamp starts at same point as time in external file
  int start_date=scorpio::get_int_attribute_c2f(datafile.c_str(),"start_date");
  int start_time=scorpio::get_int_attribute_c2f(datafile.c_str(),"start_time");
  int start_year=int(start_date/10000);
  int start_month=int((start_date-start_year*10000)/100);
  int start_day=int(start_date-start_year*10000-start_month*100);
  int start_hour=int(start_time/10000);
  int start_min=int((start_time-start_hour*10000)/100);
  int start_sec=int(start_time-start_hour*10000-start_min*100);

  EKAT_REQUIRE_MSG(start_year==ts0.get_date()[0],
		   "ERROR: The start year from the nudging file is "\ 
		   "different than the internal simulation start year\n");
  EKAT_REQUIRE_MSG(start_month==ts0.get_date()[1],
		   "ERROR: The start month from the nudging file is "\ 
		   "different than the internal simulation start month\n");
  EKAT_REQUIRE_MSG(start_day==ts0.get_date()[2],
		   "ERROR: The start day from the nudging file is "\
		   "different than the internal simulation start day\n");
  EKAT_REQUIRE_MSG(start_hour==ts0.get_time()[0],
		   "ERROR: The start hour from the nudging file is "\
		   "different than the internal simulation start hour\n");
  EKAT_REQUIRE_MSG(start_min==ts0.get_time()[1],
		   "ERROR: The start minute from the nudging file is "\
		   "different than the internal simulation start minute\n");
  EKAT_REQUIRE_MSG(start_sec==ts0.get_time()[2],
		   "ERROR: The start second from the nudging file is "\
		   "different than the internal simulation start second\n");
		   
  //Initialize before and after data
  NudgingData_bef.init(m_num_cols,m_num_src_levs,true);
  NudgingData_bef.time = -999;
  NudgingData_aft.init(m_num_cols,m_num_src_levs,true);
  NudgingData_aft.time = -999;
}

void Nudging::time_interpolation (const int time_s) {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using MemberType = typename KT::MemberType;
  
  const int time_index = time_s/time_step_file;
  double time_step_file_d = time_step_file;
  double w_bef = ((time_index+1)*time_step_file-time_s) / time_step_file_d;
  double w_aft = (time_s-(time_index)*time_step_file) / time_step_file_d;
  
  const int num_cols = NudgingData_aft.T_mid.extent(0);
  const int num_vert_packs = NudgingData_aft.T_mid.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);

  auto T_mid_bef = view_2d<Real>("",NudgingData_bef.T_mid.extent_int(0),
				 NudgingData_bef.T_mid.extent_int(1));
  Kokkos::deep_copy(T_mid_bef,NudgingData_bef.T_mid);
  auto T_mid_aft = view_2d<Real>("",NudgingData_aft.T_mid.extent_int(0),
				 NudgingData_aft.T_mid.extent_int(1));
  Kokkos::deep_copy(T_mid_aft,NudgingData_aft.T_mid);
  auto hw_bef = view_3d<Real>("",NudgingData_bef.hw.extent_int(0),2,
				 NudgingData_bef.hw.extent_int(2));
  Kokkos::deep_copy(hw_bef,NudgingData_bef.hw);
  auto hw_aft = view_3d<Real>("",NudgingData_aft.hw.extent_int(0),2,
				 NudgingData_aft.hw.extent_int(2));
  Kokkos::deep_copy(hw_aft,NudgingData_aft.hw);
  auto p_mid_bef = view_2d<Real>("",NudgingData_bef.p_mid.extent_int(0),
				 NudgingData_bef.p_mid.extent_int(1));
  Kokkos::deep_copy(p_mid_bef,NudgingData_bef.p_mid);
  auto p_mid_aft = view_2d<Real>("",NudgingData_aft.p_mid.extent_int(0),
				 NudgingData_aft.p_mid.extent_int(1));
  Kokkos::deep_copy(p_mid_aft,NudgingData_aft.p_mid);
  auto qv_bef = view_2d<Real>("",NudgingData_bef.qv.extent_int(0),
				 NudgingData_bef.qv.extent_int(1));
  Kokkos::deep_copy(qv_bef,NudgingData_bef.qv);
  auto qv_aft = view_2d<Real>("",NudgingData_aft.qv.extent_int(0),
				 NudgingData_aft.qv.extent_int(1));
  Kokkos::deep_copy(qv_aft,NudgingData_aft.qv);
  
  Kokkos::parallel_for("nudging_time_interpolation", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto T_mid_1d = ekat::subview(T_mid_ext,icol);
      auto T_mid_bef_1d = ekat::subview(T_mid_bef,icol);
      auto T_mid_aft_1d = ekat::subview(T_mid_aft,icol);
      auto hw_2d = ekat::subview(horiz_winds_ext, icol);
      auto hw_bef_2d = ekat::subview(hw_bef, icol);
      auto hw_aft_2d = ekat::subview(hw_aft, icol);
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
	  hw_2d(0,k)=w_bef*hw_bef_2d(0,k) + w_aft*hw_aft_2d(0,k);
	  hw_2d(1,k)=w_bef*hw_bef_2d(1,k) + w_aft*hw_aft_2d(1,k);
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
  if (time_s >= NudgingData_aft.time)
    {
      const int time_index = time_s/time_step_file;
      data_input->read_variables(time_index-1);
      Kokkos::deep_copy(NudgingData_bef.T_mid,T_mid_ext);
      Kokkos::deep_copy(NudgingData_bef.p_mid,p_mid_ext);
      Kokkos::deep_copy(NudgingData_bef.hw,horiz_winds_ext);
      Kokkos::deep_copy(NudgingData_bef.qv,qv_ext);
      NudgingData_bef.time = time_step_file*(time_index);

      data_input->read_variables(time_index);
      Kokkos::deep_copy(NudgingData_aft.T_mid,T_mid_ext);
      Kokkos::deep_copy(NudgingData_aft.p_mid,p_mid_ext);
      Kokkos::deep_copy(NudgingData_aft.hw,horiz_winds_ext);
      Kokkos::deep_copy(NudgingData_aft.qv,qv_ext);
      NudgingData_aft.time = time_step_file*(time_index+1);
    }
  
}

  
// =========================================================================================
void Nudging::run_impl (const int dt)
{
  using namespace scream::vinterp;

  //Have to add dt because first time iteration is at 0 seconds where you will not have
  //any data from the field. The timestamp is only iterated at the end of the
  //full step in scream.
  auto ts = timestamp()+dt;
  auto time_since_zero = ts.seconds_from(ts0);

  //Don't nudge for times less than first time step in nudged file to avoid
  //extrapolation
  if (time_since_zero < time_step_file){ return;}

  //update time state information and check whether need to update data
  update_time_step(time_since_zero);
  
  //perform time interpolation
  time_interpolation(time_since_zero);

  auto T_mid       = get_field_out("T_mid").get_view<mPack**>();
  auto qv          = get_field_out("qv").get_view<mPack**>();
  auto hw          = get_field_out("horiz_winds").get_view<mPack***>();

  const auto& p_mid    = get_field_in("p_mid").get_view<const mPack**>();
  auto p_mid_v = view_Nd<mPack,2>("",p_mid.extent_int(0),p_mid.extent_int(1));
  Kokkos::deep_copy(p_mid_v,p_mid);
  
  const view_Nd<mPack,2> T_mid_ext_p(reinterpret_cast<mPack*>(T_mid_ext.data()),
				     m_num_cols,m_num_src_levs); 
  const view_Nd<mPack,2> p_mid_ext_p(reinterpret_cast<mPack*>(p_mid_ext.data()),
				     m_num_cols,m_num_src_levs);
  const view_Nd<mPack,2> qv_ext_p(reinterpret_cast<mPack*>(qv_ext.data()),
				  m_num_cols,m_num_src_levs);
  const view_Nd<mPack,3> hw_ext_p(reinterpret_cast<mPack*>(horiz_winds_ext.data()),
				  m_num_cols,2,m_num_src_levs);

  //Now perform the vertical interpolation from the external file to the internal
  //field levels
  perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                           p_mid_v,
                                           T_mid_ext_p,
                                           T_mid,
                                           m_num_src_levs,
                                           m_num_levs);

  perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                           p_mid_v,
                                           qv_ext_p,
                                           qv,
                                           m_num_src_levs,
                                           m_num_levs);

  perform_vertical_interpolation<Real,1,3>(p_mid_ext_p,
                                           p_mid_v,
                                           hw_ext_p,
                                           hw,
                                           m_num_src_levs,
                                           m_num_levs);

  //This code removes any masked values and sets them to the corresponding
  //values at the lowest/highest pressure levels from the external file
  const int num_cols = T_mid.extent(0);
  const int num_vert_packs = T_mid.extent(1);
  const int num_vert_packs_ext = T_mid_ext.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
  Kokkos::parallel_for("correct_for_masked_values", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto T_mid_1d = ekat::subview(T_mid,icol);
      auto T_mid_ext_1d = ekat::subview(T_mid_ext,icol);
      auto hw_2d = ekat::subview(hw, icol);
      auto hw_ext_2d = ekat::subview(horiz_winds_ext, icol);
      auto p_mid_1d = ekat::subview(p_mid,icol);
      auto p_mid_ext_1d = ekat::subview(p_mid_ext,icol);
      auto qv_1d = ekat::subview(qv,icol);
      auto qv_ext_1d = ekat::subview(qv_ext,icol);
      const auto range = Kokkos::TeamThreadRange(team, num_vert_packs);
      Kokkos::parallel_for(range, [&] (const Int & k) {
        const auto above_max = p_mid_1d(k) > p_mid_ext_1d(num_vert_packs_ext-1);
        const auto below_min = p_mid_1d(k) < p_mid_ext_1d(0);
        T_mid_1d(k).set(above_max,T_mid_ext_1d(num_vert_packs_ext-1));
        T_mid_1d(k).set(below_min,T_mid_ext_1d(0));
        qv_1d(k).set(above_max,qv_ext_1d(num_vert_packs_ext-1));
        qv_1d(k).set(below_min,qv_ext_1d(0));
        hw_2d(0,k).set(above_max,hw_ext_2d(0,num_vert_packs_ext-1));
        hw_2d(0,k).set(below_min,hw_ext_2d(0,0));
        hw_2d(1,k).set(above_max,hw_ext_2d(1,num_vert_packs_ext-1));
        hw_2d(1,k).set(below_min,hw_ext_2d(1,0));
      });
      team.team_barrier();
  });
  Kokkos::fence();   
    
}

// =========================================================================================
void Nudging::finalize_impl()
{
  data_input->finalize();
}

} // namespace scream
