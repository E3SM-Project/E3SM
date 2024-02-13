#ifndef SCREAM_NUDGING_HPP
#define SCREAM_NUDGING_HPP

#include "share/atm_process/atmosphere_process.hpp"

#include "share/util/eamxx_time_interpolation.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <string>

namespace scream
{

/*
 * The class responsible to handle the nudging of variables
*/
class Nudging : public AtmosphereProcess
{
public:
  // enum to track how the source pressure levels are defined
  enum SourcePresType {
    // DEFAULT - source data should include time/spatially varying p_mid with dimensions (time, col, lev)
    TIME_DEPENDENT_3D_PROFILE,
    // source data includes p_levs which is a static set of levels in both space and time, with dimensions (lev),
    STATIC_1D_VERTICAL_PROFILE,
    // hybrid data includes hyai(ilev),hybi(ilev),hyam(lev),hybm(lev),ilev(ilev),lev(lev),P0,PS(time, ncol)
    TIME_DEPENDENT_3D_HYBRID
  };

  // Constructors
  Nudging (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "Nudging"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) override;

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void run_impl (const double dt) override;

  // Internal function to apply nudging at specific timescale
  // NOTE: this method will handle weighted and cutoff cases as well
  void apply_tendency (Field &state, const Field &nudge, const Real dt) const;

protected:

  Field get_field_out_wrap(const std::string& field_name);

  // The two other main overrides for the subcomponent
  void initialize_impl (const RunType run_type) override;
  void finalize_impl   () override;

  // Creates an helper field, not to be shared with the AD's FieldManager
  Field create_helper_field (const std::string& name,
                            const FieldLayout& layout,
                            const std::string& grid_name,
                            const int ps = 1);

  // Retrieve a helper field
  Field get_helper_field (const std::string& name) const { return m_helper_fields.at(name); }

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Keep track of field dimensions and the iteration count
  int m_num_cols;
  int m_num_levs;
  int m_num_src_levs;
  int m_timescale;
  bool m_use_weights;
  std::vector<std::string> m_datafiles;
  std::string              m_static_vertical_pressure_file;
  // add nudging weights for regional nudging update
  std::string              m_weights_file;

  SourcePresType m_src_pres_type;
  
  std::map<std::string,Field> m_helper_fields;

  std::vector<std::string> m_fields_nudge;

  /* Nudge from coarse data */
  // if true, remap coarse data to fine grid
  bool m_refine_remap;
  // file containing coarse data mapping
  std::string m_refine_remap_file;
  // (refining) remapper object
  std::shared_ptr<scream::AbstractRemapper> m_horiz_remapper;
  // (refining) remapper vertical cutoff
  Real m_refine_remap_vert_cutoff;

  util::TimeInterpolation m_time_interp;
}; // class Nudging

} // namespace scream

#endif // SCREAM_NUDGING_HPP
