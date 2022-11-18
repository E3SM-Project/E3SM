#include "atmosphere_nudging.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

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
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column
  m_dofs_gids = m_grid->get_dofs_gids();
  m_min_global_dof    = m_grid->get_global_min_dof_gid();

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  constexpr int ps = Pack::n;
  add_field<Updated>("T_mid"      , scalar3d_layout_mid, K,     grid_name, ps);
  add_field<Updated>("p_mid"      , scalar3d_layout_mid, Pa,     grid_name, ps);
  
}
// =========================================================================================
size_t NUDGING::requested_buffer_size_in_bytes() const
{

}

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
