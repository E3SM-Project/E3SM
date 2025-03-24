#ifndef SCREAM_ML_CORRECTION_HPP
#define SCREAM_ML_CORRECTION_HPP

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <array>
#include <string>
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_lin_interp.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

/*
 * The class responsible to handle the calculation of the subgrid cloud
 * fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 */

class MLCorrection : public AtmosphereProcess {
 public:
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  // Constructors
  MLCorrection(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The type of subcomponent
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name() const { return "MLCorrection"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
  // The three main overrides for the subcomponent
  void initialize_impl(const RunType run_type);
  void run_impl(const double dt);
  void finalize_impl();
  void apply_tendency(Field& base, const Field& next, const int dt);

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Keep track of field dimensions and the iteration count
  Int m_num_cols;
  Int m_num_levs;
  Field m_lat;
  Field m_lon;
  std::string m_ML_model_path_tq;
  std::string m_ML_model_path_uv;
  std::string m_ML_model_path_sfc_fluxes;
  std::vector<std::string> m_fields_ml_output_variables;
  bool m_ML_correction_unit_test;
  pybind11::module py_correction;
  pybind11::object ML_model_tq;
  pybind11::object ML_model_uv;
  pybind11::object ML_model_sfc_fluxes;
  int fpe_mask;
};  // class MLCorrection

}  // namespace scream

#endif  // SCREAM_ML_CORRECTION_HPP
