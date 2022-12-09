#include "atmosphere_nudging.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

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
  //m_dofs_gids = m_grid->get_dofs_gids();
  //m_min_global_dof    = m_grid->get_global_min_dof_gid();


  //Here get all the fields that are defined in parameter list
  //Then need to find the corresponding field layout for each of them
  //Grid should be the same

  //FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  //scalar3d_layout_mid_=scalar3d_layout_mid;
  constexpr int ps = Pack::n;
  //constexpr int ps = 16;
  add_field<Updated>("T_mid"      , scalar3d_layout_mid, K,     grid_name, ps);
  add_field<Required>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);

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

  //In here make

  //need to make grid for input file
  //auto grid = std::make_shared<PointGrid>("grid",num_local_cols,source_data_nlevs,comm);

  //Is there a way to get list of fields
  //Is there a way to get   

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
  data_input.read_variables(0);
  data_input.finalize();

  //T_mid_r_m = T_mid_r_v;
  //view_2d<Real> T_mid_r_p("T_mid_r",m_num_cols,m_num_src_levs);
  //T_mid_r_m = T_mid_r_p;
  //T_mid_r_m = T_mid_r_v;
  T_mid_r_m=  view_2d<Pack>(reinterpret_cast<Pack*>(T_mid_r_v.data()),m_num_cols,m_num_src_levs);
  p_mid_r_m=  view_2d<Pack>(reinterpret_cast<Pack*>(p_mid_r_v.data()),m_num_cols,m_num_src_levs);



  /*  
  //This is a way to get change from real to pack with pack size 16
  auto T_mid_r_m_s = ekat::scalarize(T_mid_r_m);
  int pp = Pack::n;
  //Need to parallize this
  for (int jj=0; jj<m_num_cols; jj++){
    for (int kk=0; kk<m_num_levs; kk++) {
      T_mid_r_m_s(jj,kk) = T_mid_r_v(jj,kk);
    }
  }
  */


  //p_mid_r_m = p_mid_r_v;


  //Get proper data from file, i.e.
  //Is this an initialization that happens at each time step or at the beginning before
  //everything?
  //Use AtmosphereInput from scorpio_input to read in file
  //    Is there a file that I can start with? 
  //Is there where we need to vertically interpolate to the new grid?

  
  
}

// =========================================================================================
void NUDGING::run_impl (const int dt)
{


  std::cout<<"I get in run_impl of atmosphere nudging"<<std::endl;

  //Eventually loop over all varialbes not just T_mid
  //auto T_mid          = get_field_in("T_mid").get_view<Spack**>();
  //auto T_mid          = get_field_in("T_mid").get_view<Pack**>();
  auto T_mid          = get_field_in("T_mid").get_view<Pack**>();
  //auto T_mid          = get_field_in("T_mid").get_view<Real**>();
  //Need to parallelize this

  //First do vertical interpolation
  // Need T_mid_r_m in Pack form
  // Need p_mid source in Pack form
  // Need p_mid target in Pack form
  //auto T_mid_r_m_s = ekat::scalarize(T_mid_r_m);
  std::cout<<"T_mid.extent(0): "<<T_mid.extent(0)<<std::endl;
  std::cout<<"T_mid.extent(1): "<<T_mid.extent(1)<<std::endl;
  for (int i=0; i<T_mid.extent(0); ++i) { 
    for (int k=0; k<T_mid.extent(1); ++k) {
    //for (int k=0; k<3; ++k) {
      //T_mid(i,k)=T_mid(i,k)+1;
      //T_mid(i,k)=T_mid_r_m(i,k);
      //T_mid(i,k)=T_mid_r_m_s(i,k);
      T_mid(i,k)=T_mid_r_m(i,k);
      std::cout<<"T_mid("<<i<<","<<k<<"): "<<T_mid(i,k)<<std::endl;
    }
  }
 
  //Take T_mid field from model and T_mid from file (vertically interpolated)
  //and apply the nudging so that the new T_mid output has been nudged 
  //Then generalize to any variable
  //Question of how to do this, need to get values from yaml file
  //Then put into a vector of string names?
  //Need ability for user to specify variables.
  //Similar to 
  //Also need to check whether the particular cell should be nudged
}

// =========================================================================================
void NUDGING::finalize_impl()
{
  // Do nothing
}

} // namespace scream
