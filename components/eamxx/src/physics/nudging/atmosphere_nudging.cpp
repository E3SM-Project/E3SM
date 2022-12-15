#include "atmosphere_nudging.hpp"

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
  //std::cout<<"fnames_ size: "<<fnames_.size()<<std::endl;
  //std::cout<<fnames_[0]<<std::endl;
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
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);

  //Now need to read in the file
  datafile = "io_output_test.INSTANT.nsteps_x1.np1.2000-01-01-00001.nc";
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
  view_2d<Real> T_mid_r_v("T_mid_r",m_num_cols,m_num_src_levs);
  view_2d<Real> p_mid_r_v("p_mid_r",m_num_cols,m_num_src_levs);
  auto T_mid_r_v_h       = Kokkos::create_mirror_view(T_mid_r_v);      
  auto p_mid_r_v_h       = Kokkos::create_mirror_view(p_mid_r_v);      
  //auto T_mid_r_v_h_s       = ekat::scalarize(T_mid_r_v);      
  //auto p_mid_r_v_h_s       = ekat::scalarize(p_mid_r_v);      

  host_views["T_mid_r"] = view_1d_host<Real>(T_mid_r_v_h.data(),
                                             T_mid_r_v_h.size());
  layouts.emplace("T_mid_r", scalar3d_layout_mid);

  host_views["p_mid_r"] = view_1d_host<Real>(p_mid_r_v_h.data(),
                                             p_mid_r_v_h.size());
  layouts.emplace("p_mid_r", scalar3d_layout_mid);

  //auto grid_l = std::make_shared<PointGrid>("grid_l",m_num_cols,m_num_src_levs,m_comm);
  //grid_l->set_dofs(m_num_cols);
  auto gm_l = create_gm(m_comm,m_num_cols,m_num_src_levs);
  auto grid_l = gm_l->get_grid("Point Grid");

  //This needs to be input parameter
  ekat::ParameterList data_in_params;
  data_in_params.set("Field Names",m_fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input(m_comm,data_in_params);
  //data_input.init(m_grid,host_views,layouts);
  data_input.init(grid_l,host_views,layouts);
  //This is where I read in time index
  data_input.read_variables(0);
  data_input.finalize();

  T_mid_r_v_g = T_mid_r_v;
  p_mid_r_v_g = p_mid_r_v;
  
}

// =========================================================================================
void NUDGING::run_impl (const int dt)
{
  using namespace scream::vinterp;
  std::cout<<"I get in run_impl of atmosphere nudging"<<std::endl;

  //auto ts = timestamp()+dt;
  auto ts = timestamp();
  std::cout<<"time_stamp year: "<<ts.get_year()<<std::endl;
  std::cout<<"time_stamp month: "<<ts.get_month()<<std::endl;
  std::cout<<"time_stamp day: "<<ts.get_day()<<std::endl;
  std::cout<<"time_stamp hours: "<<ts.get_hours()<<std::endl;
  std::cout<<"time_stamp minutes: "<<ts.get_minutes()<<std::endl;
  std::cout<<"time_stamp seconds: "<<ts.get_seconds()<<std::endl;

  
  //These are fields before modifications
  auto T_mid          = get_field_in("T_mid").get_view<Pack**>();
  const auto p_mid    = get_field_in("p_mid").get_view<Pack**>();

  		       
  for (int i=0; i<T_mid.extent(0); ++i) { 
    for (int k=0; k<T_mid.extent(1); ++k) {
      //T_mid(i,k)=T_mid_r_m(i,k);
      std::cout<<"Before T_mid("<<i<<","<<k<<"): "<<T_mid(i,k)<<std::endl;
    }
  }
  

  
  /*
  auto ft = T_mid.get_header_ptr()->get_tracking();
  auto pt = p_mid.get_header_ptr()->get_tracking();
  auto ts = ft.get_time_stamp();
  auto ts_p = pt.get_time_stamp();
  ft.update_time_stamp(ts+dt);
  pt.update_time_stamp(ts_p+dt);
  */
  
  //These are field values to be modified to
  view_2d<Pack> T_mid_r_m_out("T_mid_r_m_out",m_num_cols,m_num_levs);

  const view_2d<Pack> T_mid_r_m(reinterpret_cast<Pack*>(T_mid_r_v_g.data()),m_num_cols,m_num_src_levs); 
  const view_2d<Pack> p_mid_r_m(reinterpret_cast<Pack*>(p_mid_r_v_g.data()),m_num_cols,m_num_src_levs);

  perform_vertical_interpolation(p_mid_r_m,
                                 p_mid,
                                 T_mid_r_m,
                                 T_mid_r_m_out,
                                 m_num_src_levs,
                                 m_num_levs);

  std::cout<<"T_mid.extent(0): "<<T_mid.extent(0)<<std::endl;
  std::cout<<"T_mid.extent(1): "<<T_mid.extent(1)<<std::endl;
  
  const int num_cols = T_mid.extent(0);
  const int num_vert_packs = T_mid.extent(1);
  const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto T_mid_1d = ekat::subview(T_mid,icol);
      auto T_mid_r_m_1d = ekat::subview(T_mid_r_m,icol);
      const auto range = Kokkos::TeamThreadRange(team, num_vert_packs);
      Kokkos::parallel_for(range, [&] (const Int & k) {
        T_mid_1d(k)=T_mid_r_m_1d(k);
      });
      team.team_barrier();
  });
  Kokkos::fence();   
  

  /*		       
  for (int i=0; i<T_mid.extent(0); ++i) { 
    for (int k=0; k<T_mid.extent(1); ++k) {
      T_mid(i,k)=T_mid_r_m(i,k);
      std::cout<<"T_mid("<<i<<","<<k<<"): "<<T_mid(i,k)<<std::endl;
    }
  }
  */

}

// =========================================================================================
void NUDGING::finalize_impl()
{
  // Do nothing
}

} // namespace scream
