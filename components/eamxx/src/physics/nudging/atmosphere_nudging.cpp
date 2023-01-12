#include "atmosphere_nudging.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream
{

  //using namespace spa;
// =========================================================================================
NUDGING::NUDGING (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
  //m_comm=comm;
  m_fnames=m_params.get<std::vector<std::string>>("Field_Names");
  datafile=m_params.get<std::string>("Nudging_Filename");
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
  //constexpr int ps = Pack::n;
  constexpr int ps = 1;
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  /*
  for (auto &s: m_fnames) {
  //  std::cout<<"s: "<<s<<std::endl;
    if (s == "p_mid"){continue;}
      add_field<Updated>(s, scalar3d_layout_mid, K, grid_name, ps);
    //add_field<Required>(s, scalar3d_layout_mid, K, grid_name, ps);
    //add_field<Computed>(s, scalar3d_layout_mid, K, grid_name, ps);
  }
  */
  //add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);

  //Now need to read in the file
  //datafile = "io_output_test.INSTANT.nsteps_x1.np1.2000-01-01-00250.nc";
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

  //Need to make this more general so it is not ju
  //This defines the fields we are going to get from the nudging file
  using namespace ShortFieldTagsNames;
  for (auto &s: m_fnames) {
    FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
    fields_ext[s] = view_2d<Real>(s,m_num_cols,m_num_src_levs);
    auto T_mid_r_v_h       = Kokkos::create_mirror_view(fields_ext[s]);      
    host_views[s] = view_1d_host<Real>(T_mid_r_v_h.data(),T_mid_r_v_h.size());
    layouts.emplace(s, scalar3d_layout_mid);
  }

  auto grid_l = m_grid->clone("Point Grid", false);
  grid_l->reset_num_vertical_lev(m_num_src_levs);
  
  ekat::ParameterList data_in_params;
  data_in_params.set("Field Names",m_fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input0(m_comm,data_in_params);
  data_input0.init(grid_l,host_views,layouts);
  data_input=data_input0;

  //T_mid_ext = T_mid_r_v;
  T_mid_ext = fields_ext["T_mid"];
  //p_mid_ext = p_mid_r_v;
  p_mid_ext = fields_ext["p_mid"];
  ts0 = timestamp();
  //time_step_file=250;
  time_step_file=300;
  //time_step_internal=100;
  //std::cout<<"ts at start time is: "<<ts.get_seconds()<<std::endl;
  //data_input.read_variables(0);
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
  Kokkos::deep_copy(T_mid_bef,T_mid_ext);
  data_input.read_variables(time_index);
  view_2d<Real> T_mid_aft("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_aft,T_mid_ext);

  double w_bef = (time_index+1)*time_step_file-time_s;
  double w_aft = time_s-(time_index)*time_step_file;
  
  const int num_cols = T_mid_ext.extent(0);
  const int num_vert_packs = T_mid_ext.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
  Kokkos::parallel_for("nudging_time_interpolation", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto T_mid_1d = ekat::subview(T_mid_ext,icol);
      auto T_mid_bef_1d = ekat::subview(T_mid_bef,icol);
      auto T_mid_aft_1d = ekat::subview(T_mid_aft,icol);
      //auto T_mid_r_m_1d = ekat::subview(T_mid_ext,icol);
      const auto range = Kokkos::TeamThreadRange(team, num_vert_packs);
      Kokkos::parallel_for(range, [&] (const Int & k) {
	  //T_mid_1d(k)=T_mid_r_m_1d(k);
	  T_mid_1d(k)=( w_bef*T_mid_bef_1d(k) + w_aft*T_mid_aft_1d(k) )
	               / time_step_file; 
      });
      team.team_barrier();
  });
  Kokkos::fence();   
  
  
  
  //return T_mid_bef;  
}
  
// =========================================================================================
void NUDGING::run_impl (const int dt)
{
  using namespace scream::vinterp;
  std::cout<<"I get in run_impl of atmosphere nudging"<<std::endl;

  std::cout<<"dt: "<<dt<<std::endl;
  auto ts = timestamp();
  /*
  std::cout<<"time_stamp year: "<<ts.get_year()<<std::endl;
  std::cout<<"time_stamp month: "<<ts.get_month()<<std::endl;
  std::cout<<"time_stamp day: "<<ts.get_day()<<std::endl;
  std::cout<<"time_stamp hours: "<<ts.get_hours()<<std::endl;
  std::cout<<"time_stamp minutes: "<<ts.get_minutes()<<std::endl;
  std::cout<<"time_stamp seconds: "<<ts.get_seconds()<<std::endl;
  */
  std::cout<<"time since zero: "<<ts.seconds_from(ts0)<<std::endl;
  //std::cout<<"time_stamp seconds: "<<ts.get_seconds()<<std::endl;
  auto time_since_zero = ts.seconds_from(ts0);
  //if (dt < time_step_file){ return;}
  if (time_step_internal<0){
    std::cout<<"time_step_internal: "<<time_since_zero<<std::endl;
    time_step_internal=time_since_zero;
  }  
  if (time_since_zero < time_step_file){ return;}
  std::cout<<"Get before time interpolation"<<std::endl;
  time_interpolation(time_since_zero);
  std::cout<<"Get after time interpolation"<<std::endl;
  
  //TO ADD: need to loop over fields instead of just define them below.
  //probably want to define p_mid/p_in differently because will always them
  //using aPack = ekat::Pack<Real,1>;
  
  std::cout<<"Get before fields"<<std::endl;
  //These are fields before modifications
  const auto& T_mid          = get_field_out("T_mid").get_view<mPack**>();
  std::cout<<"I get here 0"<<std::endl;
  const auto& p_mid    = get_field_in("p_mid").get_view<const ekat::Pack<Real,1>**>();
  //const auto& p_mid    = get_field_in("p_mid").get_view<const ekat::Pack<Real,1>**>();
  std::cout<<"I get here 1"<<std::endl;
  std::cout<<"p_mid.extent_int(0): "<<p_mid.extent_int(0)<<std::endl;
  std::cout<<"p_mid.extent_int(1): "<<p_mid.extent_int(1)<<std::endl;
  std::cout<<"T_mid.extent_int(0): "<<T_mid.extent_int(0)<<std::endl;
  std::cout<<"T_mid.extent_int(1): "<<T_mid.extent_int(1)<<std::endl;
  std::cout<<"I get here 1"<<std::endl;
  auto p_mid_v = view_Nd<mPack,2>("",p_mid.extent_int(0),p_mid.extent_int(1));
  std::cout<<"Get after fields"<<std::endl;
  Kokkos::deep_copy(p_mid_v,p_mid);
  
  //This is field out of vertical interpolation
  //view_Nd<mPack,2> T_mid_out("T_mid_r_m_out",m_num_cols,m_num_levs);
  view_Nd<mPack,2> T_mid_out("T_mid_r_m_out",m_num_cols,T_mid.extent(1));
  //These are field values from external file
  const view_Nd<mPack,2> T_mid_ext_p(reinterpret_cast<mPack*>(T_mid_ext.data()),
  				m_num_cols,m_num_src_levs); 
  const view_Nd<mPack,2> p_mid_ext_p(reinterpret_cast<mPack*>(p_mid_ext.data()),
				m_num_cols,m_num_src_levs);
  std::cout<<"p_mid_v.extent(1): "<<p_mid_v.extent(1)<<std::endl;
  std::cout<<"T_mid_out.extent(1): "<<T_mid_out.extent(1)<<std::endl;
  std::cout<<"Get after define new fields"<<std::endl;
  
  perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                           p_mid_v,
                                           T_mid_ext_p,
                                           T_mid_out,
                                           m_num_src_levs,
                                           m_num_levs);
  std::cout<<"Get after vertical interpolation"<<std::endl;
  Kokkos::deep_copy(T_mid,T_mid_out);  

}

// =========================================================================================
void NUDGING::finalize_impl()
{
  // Do nothing
  data_input.finalize();
}

} // namespace scream
