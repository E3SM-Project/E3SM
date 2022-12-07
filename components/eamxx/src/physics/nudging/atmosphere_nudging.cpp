#include "atmosphere_nudging.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"


namespace scream
{
  //using namespace spa;
// =========================================================================================
NUDGING::NUDGING (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
  m_comm=comm;
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

  //FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  //scalar3d_layout_mid_=scalar3d_layout_mid;
  constexpr int ps = Pack::n;
  //constexpr int ps = 16;
  add_field<Updated>("T_mid"      , scalar3d_layout_mid, K,     grid_name, ps);
  //add_field<Required>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);

  // Initialize the pio_subsystem for this test:
  //scorpio::eam_init_pio_subsystem(m_comm);   // Gather the initial PIO subsystem data creater by component coupler

  //Now need to read in the file
  datafile = "io_output_test.INSTANT.nsteps_x1.np1.2000-01-01-00001.nc";
  scorpio::register_file(datafile,scorpio::Read);
  m_num_src_levs = scorpio::get_dimlen_c2f(datafile.c_str(),"lev");
  scorpio::eam_pio_closefile(datafile);
  std::cout<<"m_num_src_levs: "<<m_num_src_levs<<std::endl;
  //need to make grid for input file
  //auto grid = std::make_shared<PointGrid>("grid",num_local_cols,source_data_nlevs,comm);

  //Is there a way to get list of fields
  //Is there a way to get   



  std::vector<std::string> fnames = {"T_mid_r"};
  ekat::ParameterList data_in_params;
  data_in_params.set("Field Names",fnames);
  data_in_params.set("Filename",datafile);
  data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput data_input(m_comm,data_in_params);
  view_2d<Real> T_mid_r_v("T_mid_r",m_num_cols,m_num_src_levs);
  auto T_mid_r_v_h       = Kokkos::create_mirror_view(T_mid_r_v);      

  std::map<std::string,view_1d_host<Real>> host_views;
  std::map<std::string,FieldLayout>  layouts;
  //FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  //FieldLayout scalar3d_layout_mid_t { {COL,LEV}, {3, 33} };
  //FieldLayout scalar3d_layout_mid { {COL, LEV}, {m_num_cols, m_num_levs} };
  //FieldLayout scalar3d_layout_mid { {COL, LEV}, {m_num_cols, m_num_levs} };
  host_views["T_mid_r"] = view_1d_host<Real>(T_mid_r_v_h.data(),T_mid_r_v_h.size());
  layouts.emplace("T_mid_r", scalar3d_layout_mid);

  data_input.init(m_grid,host_views,layouts);
  data_input.read_variables(0);
  data_input.finalize();
  T_mid_r_m = T_mid_r_v;
  
}
// =========================================================================================
/*
size_t NUDGING::requested_buffer_size_in_bytes() const
{

}
*/

// =========================================================================================
void NUDGING::init_buffers(const ATMBufferManager &buffer_manager)
{

}

// =========================================================================================
void NUDGING::initialize_impl (const RunType /* run_type */)
{
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
  //auto& T_mid          = get_field_out("T_mid").get_view<Spack**>();
  //const Field& f = get_field_in("T_mid");
  //auto T_mid          = get_field_in("T_mid").get_view<Spack**>();
  auto T_mid          = get_field_in("T_mid").get_view<Real**>();
  //Need to parallelize this
  for (int i=0; i<m_num_cols; ++i) { 
    for (int k=0; k<m_num_levs; ++k) {
      //T_mid(i,k)=T_mid(i,k)+1;
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
