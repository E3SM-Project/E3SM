#include "atmosphere_nudging.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream
{

  //using namespace spa;
// =========================================================================================
NUDGING::NUDGING (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_fnames=m_params.get<std::vector<std::string>>("Field_Names");
  datafile=m_params.get<std::string>("Nudging_Filename");
  time_step_file=m_params.get<Int>("Time_Step_File",-9999);
}

// =========================================================================================
void NUDGING::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  std::cout<<"Grid name is : "<<grid_name<<std::endl;
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column
  std::cout<<"m_num_levs: "<<m_num_levs<<std::endl;
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
  scorpio::eam_pio_closefile(datafile);
  std::cout<<"m_num_src_levs: "<<m_num_src_levs<<std::endl;

  
}

// =========================================================================================
void NUDGING::init_buffers(const ATMBufferManager &buffer_manager)
{

}

// =========================================================================================
void NUDGING::initialize_impl (const RunType /* run_type */)
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
  data_in_params.set("Field Names",m_fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input0(m_comm,data_in_params);
  data_input0.init(grid_l,host_views,layouts);
  data_input=data_input0;

  T_mid_ext = fields_ext["T_mid"];
  p_mid_ext = fields_ext["p_mid"];
  horiz_winds_ext = fields_ext_3d["horiz_winds"];
  qv_ext = fields_ext["qv"];
  ts0=timestamp();
}



void NUDGING::time_interpolation (const int time_s) {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using MemberType = typename KT::MemberType;
  
  const int time_index = time_s/time_step_file;
  std::cout<<"time_index: "<<time_index<<std::endl;

  data_input.read_variables(time_index-1);
  //return;
  view_2d<Real> T_mid_bef("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  view_2d<Real> p_mid_bef("",p_mid_ext.extent_int(0),p_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_bef,T_mid_ext);
  Kokkos::deep_copy(p_mid_bef,p_mid_ext);

  view_3d<Real> hw_bef("",horiz_winds_ext.extent_int(0),2,horiz_winds_ext.extent_int(2));
  view_2d<Real> qv_bef("",qv_ext.extent_int(0),qv_ext.extent_int(1));
  Kokkos::deep_copy(hw_bef,horiz_winds_ext);
  Kokkos::deep_copy(qv_bef,qv_ext);

  data_input.read_variables(time_index);
  view_2d<Real> T_mid_aft("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  view_2d<Real> p_mid_aft("",p_mid_ext.extent_int(0),p_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_aft,T_mid_ext);
  Kokkos::deep_copy(p_mid_aft,p_mid_ext);
  
  view_3d<Real> hw_aft("",horiz_winds_ext.extent_int(0),2,horiz_winds_ext.extent_int(2));
  view_2d<Real> qv_aft("",qv_ext.extent_int(0),qv_ext.extent_int(1));
  Kokkos::deep_copy(hw_aft,horiz_winds_ext);
  Kokkos::deep_copy(qv_aft,qv_ext);
  double time_step_file_d = time_step_file;
  double w_bef = ((time_index+1)*time_step_file-time_s) / time_step_file_d;
  double w_aft = (time_s-(time_index)*time_step_file) / time_step_file_d;
  
  const int num_cols = T_mid_ext.extent(0);
  const int num_vert_packs = T_mid_ext.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
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
void NUDGING::run_impl (const int dt)
{
  using namespace scream::vinterp;
  std::cout<<"I get in run_impl of atmosphere nudging"<<std::endl;

  std::cout<<"dt: "<<dt<<std::endl;
  auto ts = timestamp()+dt;
  //std::cout<<"time since zero: "<<ts.seconds_from(ts0)<<std::endl;
  //std::cout<<"Time_Step_File: "<<time_step_file<<std::endl;
  auto time_since_zero = ts.seconds_from(ts0);
  //Don't nudge for times less than first time step in nudged file to avoid
  //extrapolation
  if (time_since_zero < time_step_file){ return;}
  //if (time_since_zero > 7200){ return;}

  //Have to add dt because first time iteration is at 0 seconds where you will not have
  //any data from the field. The timestamp is only iterated at the end of the
  //full step in scream.
  time_interpolation(time_since_zero);

  auto T_mid       = get_field_out("T_mid").get_view<mPack**>();
  auto qv          = get_field_out("qv").get_view<mPack**>();
  auto hw          = get_field_out("horiz_winds").get_view<mPack***>();

  const auto& p_mid    = get_field_in("p_mid").get_view<const mPack**>();
  auto p_mid_v = view_Nd<mPack,2>("",p_mid.extent_int(0),p_mid.extent_int(1));
  Kokkos::deep_copy(p_mid_v,p_mid);
  
  view_Nd<mPack,2> T_mid_out("T_mid_out",m_num_cols,T_mid.extent_int(1));
  view_Nd<mPack,2> T_mid_ext_p("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_ext_p,T_mid_ext);

  view_Nd<mPack,2> qv_out("qv_out",m_num_cols,qv.extent_int(1));
  view_Nd<mPack,2> qv_ext_p("qv_ext_p",qv_ext.extent_int(0),qv_ext.extent_int(1));
  Kokkos::deep_copy(qv_ext_p,qv_ext);

  view_Nd<mPack,3> hw_out("hw_out",m_num_cols,hw.extent_int(1),hw.extent(2));
  view_Nd<mPack,3> hw_ext_p("hw_ext_p",
  			    horiz_winds_ext.extent_int(0),
  			    horiz_winds_ext.extent_int(1),
  			    horiz_winds_ext.extent_int(2));
  Kokkos::deep_copy(hw_ext_p,horiz_winds_ext);

  view_Nd<mPack,2> p_mid_ext_p("",p_mid_ext.extent_int(0),p_mid_ext.extent_int(1));
  Kokkos::deep_copy(p_mid_ext_p,p_mid_ext);
  
  perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                           p_mid_v,
                                           T_mid_ext_p,
                                           T_mid_out,
                                           m_num_src_levs,
                                           m_num_levs);
  Kokkos::deep_copy(T_mid,T_mid_out);

  perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                           p_mid_v,
                                           qv_ext_p,
                                           qv_out,
                                           m_num_src_levs,
                                           m_num_levs);
  Kokkos::deep_copy(qv,qv_out);

  perform_vertical_interpolation<Real,1,3>(p_mid_ext_p,
                                           p_mid_v,
                                           hw_ext_p,
                                           hw_out,
                                           m_num_src_levs,
                                           m_num_levs);
  Kokkos::deep_copy(hw,hw_out);

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
void NUDGING::finalize_impl()
{
  // Do nothing
  data_input.finalize();
}

} // namespace scream
