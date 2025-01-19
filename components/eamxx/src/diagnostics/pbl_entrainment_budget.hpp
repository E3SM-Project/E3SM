#ifndef EAMXX_ENTRAINMENT_BUDGET_DIAG
#define EAMXX_ENTRAINMENT_BUDGET_DIAG

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"

namespace scream {

/*
 * This diagnostic will compute pbl entrainment budget terms.
 */

class PBLEntrainmentBudget : public AtmosphereDiagnostic {
 public:
  using PF      = scream::PhysicsFunctions<DefaultDevice>;
  using KT      = KokkosTypes<DefaultDevice>;
  using view_2d = typename KT::template view_2d<Real>;

  // Constructors
  PBLEntrainmentBudget(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic
  std::string name() const { return "PBLEntrainmentBudget"; };

  // Set the grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

  // Let's override the init time step method
  void init_timestep(const util::TimeStamp &start_of_step) override;

  // Let's override the initialize method to set the fields below
  void initialize_impl(const RunType /*run_type*/) override;

  // A function to calculate the tracked fields of interest
  void calc_tl_qt(const view_2d &tm_v, const view_2d &pm_v, const view_2d &qv_v,
                  const view_2d &qc_v, const view_2d &tl_v,
                  const view_2d &qt_v);

  // Store fields at init_timestep
  // That is before anything takes place in atm
  Field m_prev_qt;
  Field m_prev_tl;
  util::TimeStamp m_start_t;

  // Grid info
  int m_ncols;
  int m_nlevs;

  // How many diags we have
  int m_ndiag;

  // Attribute maps for self-documentation
  std::map<std::string, int> m_index_map;
  std::map<std::string, std::string> m_units_map;

  // Which PBL inversion algorithm to use
  int m_pblinvalg;
};

}  // namespace scream

#endif  // EAMXX_ENTRAINMENT_BUDGET_DIAG
