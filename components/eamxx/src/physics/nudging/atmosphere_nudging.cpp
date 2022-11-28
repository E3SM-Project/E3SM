#include "atmosphere_nudging.hpp"

namespace scream
{
  //using namespace spa;
// =========================================================================================
NUDGING::NUDGING (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
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

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  constexpr int ps = Pack::n;
  //constexpr int ps = 16;
  add_field<Updated>("T_mid"      , scalar3d_layout_mid, K,     grid_name, ps);
  //add_field<Required>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);
  
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
  for (int i=0; i<m_num_cols; ++i) { 
    for (int k=0; k<m_num_levs; ++k) {
      T_mid(i,k)=T_mid(i,k)+1;
      std::cout<<"T_mid("<<i<<","<<k<<"): "<<T_mid(i,k)<<std::endl;
    }
  }
 
  //Take T_mid field from model and T_mid from file (vertically interpolated)
  //and apply the nudging so that the new T_mid output has been nudged 
  //Then generalize to any variable
  //Need ability for user to specify variables.
  //Similar to 
}

// =========================================================================================
void NUDGING::finalize_impl()
{
  // Do nothing
}

} // namespace scream
