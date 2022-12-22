#include "atmosphere_nudging.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream
{


std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}


  //using namespace spa;
// =========================================================================================
NUDGING::NUDGING (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
  //m_comm=comm;
  m_fnames=m_params.get<std::vector<std::string>>("Field Names");
  datafile=m_params.get<std::string>("Nudging Filename");
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

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  //constexpr int ps = Pack::n;
  constexpr int ps = 1;
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);
  for (auto &s: m_fnames) {
  //  std::cout<<"s: "<<s<<std::endl;
    if (s == "p_mid"){continue;}
    add_field<Updated>(s, scalar3d_layout_mid, K, grid_name, ps);
  }
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

  //DO I NEED TO ADD FIELDS FOR BELOW/ABOVE SO THAT CAN INTERPOLATE?
  //WHAT ABOUT MAKING THE T_mid_r same name, so T_mid?
  using namespace ShortFieldTagsNames;
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  view_2d<Real> T_mid_r_v("T_mid",m_num_cols,m_num_src_levs);
  view_2d<Real> p_mid_r_v("p_mid",m_num_cols,m_num_src_levs);
  auto T_mid_r_v_h       = Kokkos::create_mirror_view(T_mid_r_v);      
  auto p_mid_r_v_h       = Kokkos::create_mirror_view(p_mid_r_v);      
  host_views["T_mid"] = view_1d_host<Real>(T_mid_r_v_h.data(),T_mid_r_v_h.size());
  layouts.emplace("T_mid", scalar3d_layout_mid);
  host_views["p_mid"] = view_1d_host<Real>(p_mid_r_v_h.data(),p_mid_r_v_h.size());
  layouts.emplace("p_mid", scalar3d_layout_mid);

  //auto grid_l = std::make_shared<PointGrid>("grid_l",m_num_cols,m_num_src_levs,m_comm);
  //grid_l->set_dofs(m_num_cols);
  auto gm_l = create_gm(m_comm,m_num_cols,m_num_src_levs);
  auto grid_l = gm_l->get_grid("Point Grid");


  ekat::ParameterList data_in_params;
  data_in_params.set("Field Names",m_fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input0(m_comm,data_in_params);
  data_input0.init(grid_l,host_views,layouts);
  data_input=data_input0;

  T_mid_ext = T_mid_r_v;
  p_mid_ext = p_mid_r_v;
  auto ts = timestamp();
  time_step_file=250;
  time_step_internal=100;
  //std::cout<<"ts at start time is: "<<time.seconds_from(ts)<<std::endl;
  //data_input.read_variables(0);
}



void NUDGING::time_interpolation (const int dt) {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using MemberType = typename KT::MemberType;
  
  const int time_index = dt/time_step_file;
  std::cout<<"time_index: "<<time_index<<std::endl;
  data_input.read_variables(time_index-1);
  //return;
  view_2d<Real> T_mid_bef("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_bef,T_mid_ext);
  data_input.read_variables(time_index);
  view_2d<Real> T_mid_aft("",T_mid_ext.extent_int(0),T_mid_ext.extent_int(1));
  Kokkos::deep_copy(T_mid_aft,T_mid_ext);

  double w_bef = (time_index+1)*time_step_file-dt;
  double w_aft = dt-(time_index)*time_step_file;
  
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
  
  //auto ts = timestamp()+dt;
  //std::cout<<"time_stamp seconds: "<<ts.get_seconds()<<std::endl;

  if (dt < time_step_file){ return;}
  time_interpolation(dt);
  
  //These are fields before modifications
  auto T_mid          = get_field_in("T_mid").get_view<Pack**>();
  const auto p_mid    = get_field_in("p_mid").get_view<Pack**>();
 
  //This is field out of vertical interpolation
  view_2d<Pack> T_mid_out("T_mid_r_m_out",m_num_cols,m_num_levs);
  //These are field values from external file
  const view_2d<Pack> T_mid_ext_p(reinterpret_cast<Pack*>(T_mid_ext.data()),
  				m_num_cols,m_num_src_levs); 
  const view_2d<Pack> p_mid_ext_p(reinterpret_cast<Pack*>(p_mid_ext.data()),
				m_num_cols,m_num_src_levs);

  perform_vertical_interpolation(p_mid_ext_p,
                                 p_mid,
                                 T_mid_ext_p,
                                 T_mid_out,
                                 m_num_src_levs,
                                 m_num_levs);
  
  std::cout<<"T_mid.extent(0): "<<T_mid.extent(0)<<std::endl;
  std::cout<<"T_mid.extent(1): "<<T_mid.extent(1)<<std::endl;
  Kokkos::deep_copy(T_mid,T_mid_out);  

}

// =========================================================================================
void NUDGING::finalize_impl()
{
  // Do nothing
  data_input.finalize();
}

} // namespace scream
