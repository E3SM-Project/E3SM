#ifndef EAMXX_MAM_WETSCAV_HPP
#define EAMXX_MAM_WETSCAV_HPP

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For declaring wetscav class derived from atm process class
#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

// For component name
#include <string>

namespace scream {

/*
 * The class responsible to handle the aerosol wetscavenging
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 */

class MAMWetscav : public MAMGenericInterface {
  using KT      = ekat::KokkosTypes<DefaultDevice>;
  using view_2d = typename KT::template view_2d<Real>;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

 public:
  // Constructor
  MAMWetscav(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the subcomponent
  std::string name() const override { return "mam4_wetscav"; }

  // Set the grid and input output variables
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  // ON HOST, returns the number of bytes of device memory needed by the above
  // Buffer type given the number of columns and vertical levels
  size_t requested_buffer_size_in_bytes() const override {
    return mam_coupling::buffer_size(ncol_, nlev_, num_2d_scratch_, work_len_);
  }
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl() override{/*Do nothing*/};

 private:
  // -----------------------------------------------
  // Local variables
  // ------------------------------------------------

  // Number of aerosol modes
  static constexpr int ntot_amode_ = mam4::AeroConfig::num_modes();

  // Work arrays
  view_2d work_;

  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now

  view_2d sh_frac_;

  // Deep convective cloud fraction [fraction]
  view_2d dp_frac_;

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  view_2d evapcsh_;

  view_2d evapcdp_;

  // Rain production, shallow convection [kg/kg/s]
  view_2d rprdsh_;

  // Rain production, deep convection [kg/kg/s]
  view_2d rprddp_;

  // In cloud water mixing ratio, deep convection
  view_2d icwmrdp_;

  // In cloud water mixing ratio, shallow convection
  view_2d icwmrsh_;

  // Detraining cld H20 from deep convection [kg/kg/s]
  view_2d dlf_;

  int num_2d_scratch_= 48;

  int work_len_=0;

  void set_work_len();

  // Aerosol states
  mam_coupling::AerosolState dry_aero_tends_;

  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

};  // class MAMWetscav

}  // namespace scream

#endif  // EAMXX_MAM_WETSCAV_HPP
