#ifndef EAMXX_MAM_OPTICS_HPP
#define EAMXX_MAM_OPTICS_HPP

#include <physics/mam/mam_coupling.hpp>
#include <physics/mam/mam_aerosol_optics_read_tables.hpp>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>
#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosol optical properties. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMOptics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

public:

  // Constructor
  MAMOptics(const ekat::Comm& comm, const ekat::ParameterList& params);

protected_except_cuda:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

private_except_cuda:
  // state variable
  mam_coupling::view_3d state_q_,  qqcw_;// odap_aer_,
  mam_coupling::view_2d ext_cmip6_lw_;

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of shortwave and longwave radiation bands
  int nswbands_, nlwbands_;



  mam_coupling::view_2d mass_, radsurf_, logradsurf_  ;
  mam_coupling::view_3d cheb_, dgnumwet_m_, dgnumdry_m_;
  mam_coupling::complex_view_3d specrefindex_; // work array
  mam_coupling::view_3d qaerwat_m_, ext_cmip6_lw_inv_m_;
  // FIXME: move this values to mam_coupling
  mam_coupling::const_view_2d z_mid_, z_iface_, p_int_, p_del_;

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // atmospheric and aerosol state variables
  // mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;
  mam_coupling::AerosolState  wet_aero_;//,

  // inputs:
  mam_coupling::view_3d ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_;

  // These inputs maybe are from a netCDF file:
  mam_coupling::complex_view_2d specrefndxsw_; // complex refractive index for water visible
  mam_coupling::complex_view_2d specrefndxlw_; // complex refractive index for water infrared

  // FIXME: we need to save these values in a different file.
  // set complex representation of refractive indices as module data
  // I need netcdf files : read_water_refindex(modal_aer_opt.F90)
  Kokkos::complex<Real> crefwlw_[mam4::modal_aer_opt::nlwbands];
  Kokkos::complex<Real> crefwsw_[mam4::modal_aer_opt::nswbands];

  // Inputs from netCDF files. I already got files and code to read them.
  // ready means: I added code to read this table from netcdf. I will delete this comment.
  mam_coupling::view_3d  abspsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  mam_coupling::view_3d  extpsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  mam_coupling::view_3d  asmpsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready

  mam_coupling::view_1d refrtabsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  mam_coupling::view_1d refitabsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready

  mam_coupling::view_3d absplw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready
  mam_coupling::view_1d refrtablw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready
  mam_coupling::view_1d refitablw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready

  // work arrays
  mam_coupling::view_2d air_density_;
  mam_coupling::view_3d ext_cmip6_sw_inv_m_;

  // aerosol processes
  //std::unique_ptr<mam4::OpticsProcess> optics_;
  // std::unique_ptr<mam4::CalcSizeProcess> calcsize_process_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMOptics

} // namespace scream

#endif // EAMXX_MAM_OPTICS_HPP
