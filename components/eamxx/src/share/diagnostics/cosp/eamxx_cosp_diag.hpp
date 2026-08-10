#ifndef SCREAM_COSP_DIAGNOSTIC_HPP
#define SCREAM_COSP_DIAGNOSTIC_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * The COSP satellite simulator, as a diagnostic.
 *
 * COSP produces a whole suite of fields from a single (expensive) pass, so this
 * is a multi-output diagnostic: it declares all of output_names() via
 * add_output(), and a request for any one of them builds this one diag. Two
 * streams asking for two different COSP fields therefore run COSP once.
 *
 * This is a port of the Cosp atmosphere process, and computes exactly what it
 * computed, including deriving z_mid internally rather than taking it as an
 * input. The two must NOT both be active in one run: COSP keeps global state in
 * its Fortran modules, which both would initialize.
 *
 * COSP is expensive, so -- like the process before it -- it throttles itself,
 * via 'cosp_frequency' and 'cosp_frequency_units'. Between runs the output
 * fields simply keep the values of the last one, which is what the process did
 * too. A diagnostic has no yaml section of its own, so these come from the
 * 'diag_params: Cosp:' sublist of the output yaml (see build_diag_tree).
 */
class CospDiag : public AbstractDiagnostic
{
public:
  CospDiag (const ekat::Comm& comm, const ekat::ParameterList& params,
            const std::shared_ptr<const AbstractGrid>& grid);

  ~CospDiag ();

  // The name of the diagnostic CLASS (not the computed fields)
  std::string name () const override { return "Cosp"; }

  // The fields this diag computes. Static, because the factory registration has
  // to advertise them before any instance exists.
  static std::vector<std::string> output_names ();

protected:
  void initialize_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif
  void compute_impl () override;

protected:
  // Create the coordinate vars for the tau/prs/cth dimensions, unless the grid
  // already carries them.
  void create_cosp_geometry_data () const;

  // Whether COSP should actually run at the timestamp we are being computed at.
  // Always true the first time, so that the outputs are never left unset.
  bool cosp_do () const;

  // How often to run COSP, and in what units ('steps' or 'hours')
  Int                         m_cosp_frequency;
  ekat::CaseInsensitiveString m_cosp_frequency_units;

  // Timestamp of the last COSP run, and whether there was one at all
  util::TimeStamp m_last_cosp_ts;
  bool            m_cosp_ran = false;

  // Field dimensions
  Int m_num_cols;
  Int m_num_subcols;
  Int m_num_levs;
  Int m_num_tau = 7;
  Int m_num_ctp = 7;
  Int m_num_cth = 16;

  // Temporaries
  // TODO: use atm buffer instead
  Field m_z_mid;
  Field m_z_int;
  Field m_sunlit_real;
}; // class CospDiag

} // namespace scream

#endif // SCREAM_COSP_DIAGNOSTIC_HPP
