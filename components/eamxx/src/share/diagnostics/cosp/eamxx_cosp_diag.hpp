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
 * NOTE: the process throttled itself with 'cosp_frequency'. A diagnostic has no
 *       step counter, and the output stream decides when it is evaluated, so
 *       that knob is gone: COSP now runs whenever the stream asks for it.
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
