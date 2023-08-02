#ifndef SCREAM_NUDGING_HPP
#define SCREAM_NUDGING_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_lin_interp.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/util/scream_vertical_interpolation.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "physics/nudging/nudging_functions.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the nudging of variables
*/

class Nudging : public AtmosphereProcess
{
public:
  using NudgingFunc = nudging::NudgingFunctions;
  using mPack = ekat::Pack<Real,1>;
  using KT = KokkosTypes<DefaultDevice>;
  
  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_Nd_host = typename KT::template view_ND<S,N>::HostMirror;

  template <typename S>
  using view_1d_host = view_Nd_host<S,1>;

  template <typename S>
  using view_2d_host = view_Nd_host<S,2>;

  // Constructors
  Nudging (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Nudging"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  //Update the time step
  void update_time_step(const int time_s);
  
  //Time interpolation function
  void time_interpolation(const int time_s);

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void run_impl        (const double dt);
  
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Keep track of field dimensions and the iteration count
  int m_num_cols; 
  int m_num_levs;
  int m_num_src_levs;
  int time_step_file;
  std::string datafile;
  std::map<std::string,view_1d_host<Real>> host_views;
  std::map<std::string,FieldLayout>  layouts;
  std::vector<std::string> m_fnames;
  std::map<std::string,view_2d<Real>> fields_ext;
  std::map<std::string,view_2d_host<Real>> fields_ext_h;
  view_2d<Real> T_mid_ext;
  view_2d<Real> p_mid_ext;
  view_2d<Real> qv_ext;
  view_2d<Real> u_ext;
  view_2d<Real> v_ext;
  TimeStamp ts0;
  NudgingFunc::NudgingData NudgingData_bef;
  NudgingFunc::NudgingData NudgingData_aft;
  AtmosphereInput data_input;
}; // class Nudging

} // namespace scream

#endif // SCREAM_NUDGING_HPP
