// ============================================================
//  atm_osc_intermittency.hpp
//
//  Oscillation intermittency: outputs 1.0 at each (col, lev) if a
//  2Δt oscillation is currently active, 0.0 otherwise.  Time-
//  averaging to produce a fractional occurrence rate is handled
//  at the higher output level (e.g. avg over the output interval).
//
//  State machine per (col, lev):
//   - A step is a "sign-alternation" if
//       sign(BT(t)) ≠ sign(BT(t-1))  AND  |BT(t)| > noise_floor
//                                     AND  |BT(t-1)| > noise_floor
//   - m_n_alt is the count of consecutive sign-alternations.
//     Increments when the condition holds; resets to 0 otherwise.
//   - The step is flagged "oscillating" when m_n_alt >= min_alt_count (K).
//     Requiring K consecutive alternations filters isolated spikes,
//     which produce at most ~2 alternations before the pattern breaks.
//
//  Carry-forward fields (never reset): m_f_pre, m_bt_pre, m_n_alt
// ============================================================
#ifndef SCREAM_ATM_OSC_INTERMITTENCY_HPP
#define SCREAM_ATM_OSC_INTERMITTENCY_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

class AtmOscIntermittency : public AtmosphereDiagnostic {
public:
  AtmOscIntermittency(const ekat::Comm& comm,
                      const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "AtmOscIntermittency"; }

  void create_requests ();

protected:
  void initialize_impl(const RunType run_type) override;
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl() override;

private:
  std::string m_name;
  std::string m_diag_name;

  Real m_noise_floor;   // minimum |BT| to count as a sign-alternation
  int  m_min_alt_count; // K: consecutive alternations required to flag

  // Carry-forward state (persist across output intervals)
  Field m_f_pre;    // φ(t-1)
  Field m_bt_pre;   // BT(t-1)
  // m_n_alt stored as Real so it resides on device; values are small exact ints.
  Field m_n_alt;
};

} // namespace scream
#endif // SCREAM_ATM_OSC_INTERMITTENCY_HPP
