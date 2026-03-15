#ifndef SCREAM_COSP_HPP
#define SCREAM_COSP_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{

/*
 * COSP (CFMIP Observation Simulator Package) diagnostic.
 * Simulates satellite observations from model state for comparison
 * with real satellite data. Produces multiple output fields:
 *   isccp_cldtot, isccp_ctptau, modis_ctptau, misr_cthtau
 *
 * This is a multi-output diagnostic: the output manager creates it
 * when any of its output fields are requested in the output YAML.
 * COSP runs at the frequency determined by the output stream.
*/

class Cosp : public AtmosphereDiagnostic
{

public:
  // Constructors
  Cosp (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return "Cosp"; }

  // Set the grid and declare fields
  void create_requests ();

  // Return all diagnostic output field names
  std::vector<std::string> get_diagnostic_names () const override;

protected:

  void initialize_impl (const RunType run_type);
#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void finalize_impl   ();

  // Keep track of field dimensions
  Int m_num_cols;
  Int m_num_subcols;
  Int m_num_levs;
  Int m_num_tau = 7;
  Int m_num_ctp = 7;
  Int m_num_cth = 16;

  std::shared_ptr<const AbstractGrid> m_grid;

  // Scratch fields for height computation
  Field m_z_mid;
  Field m_z_int;
}; // class Cosp

} // namespace scream

#endif // SCREAM_COSP_HPP
