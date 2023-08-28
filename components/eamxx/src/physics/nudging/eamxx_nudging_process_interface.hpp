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

  // The two other main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const FieldLayout& layout,
                            const std::string& grid_name,
                            const int ps=0);

  // Query if a local field exists
  bool has_helper_field (const std::string& name) const { return m_helper_fields.find(name)!=m_helper_fields.end(); }
  // Retrieve a helper field
  Field get_helper_field (const std::string& name) const { return m_helper_fields.at(name); }
  // Internal function to apply nudging at specific timescale
  void apply_tendency(Field& base, const Field& next, const int dt);

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Keep track of field dimensions and the iteration count
  int m_num_cols;
  int m_num_levs;
  int m_num_src_levs;
  int m_time_step_file;
  int m_timescale;
  std::string m_datafile;

  // Some helper fields.
  std::map<std::string,Field> m_helper_fields;

  std::map<std::string,view_2d<Real>> m_fields_ext;
  std::map<std::string,view_2d_host<Real>> m_fields_ext_h;

  TimeStamp m_ts0;
  NudgingFunc::NudgingData m_NudgingData_bef;
  NudgingFunc::NudgingData m_NudgingData_aft;
  AtmosphereInput m_data_input;
}; // class Nudging

} // namespace scream

#endif // SCREAM_NUDGING_HPP
