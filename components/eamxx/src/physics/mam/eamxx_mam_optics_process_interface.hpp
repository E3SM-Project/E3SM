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

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of shortwave and longwave radiation bands
  int nswbands_, nlwbands_;

  // FIXME: move these values to mam_coupling
  mam_coupling::const_view_2d z_mid_, z_iface_, p_int_, p_del_;

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // atmospheric and aerosol state variables
  // mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;
  mam_coupling::AerosolState  wet_aero_;//,

  // inputs:
  // I do not know how to get these inputs.
  // components/eam/src/physics/rrtmg/radiation.F90 in Line 1330
  // Obtain read in values for ssa and asymmetry factor (af) from the
  //volcanic input file
  // call pbuf_get_field(pbuf, idx_ssa_sw, ssa_cmip6_sw)
  // pbuf_get_field(pbuf, idx_af_sw,  af_cmip6_sw)
  // Get extinction so as to supply to modal_aero_sw routine for computing EXTINCT variable
  // ext_cmip6_sw => null()
  // call pbuf_get_field(pbuf, idx_ext_sw, ext_cmip6_sw)
  // is_cmip6_volc true if cmip6 style volcanic file is read otherwise false
  // these inputs are prescribed.
  mam_coupling::view_3d ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_;
  //long wave extinction in the units of [1/km]
  mam_coupling::view_3d ext_cmip6_lw_;

  mam4::modal_aer_opt::AerosolOpticsDeviceData aerosol_optics_device_data_;
  mam4::modal_aer_opt::DiagnosticsAerosolOpticsSW diagnostics_aerosol_optics_sw_;

  // // These inputs maybe are from a netCDF file:
  // // complex refractive index for aersol species
  // mam_coupling::complex_view_2d specrefndxsw_; // ready
  // // complex refractive index for aersol species
  // mam_coupling::complex_view_2d specrefndxlw_;// ready

  // // FIXME: Maybe use a 1D view instead of arrays of complex.
  // Kokkos::complex<Real> crefwlw_[mam4::modal_aer_opt::nlwbands]; // ready
  // Kokkos::complex<Real> crefwsw_[mam4::modal_aer_opt::nswbands];// ready

  // // Inputs from netCDF files. I already got files and code to read them.
  // // ready means: I added code to read this table from netcdf. I will delete this comment.
  // mam_coupling::view_3d  abspsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  // mam_coupling::view_3d  extpsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  // mam_coupling::view_3d  asmpsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready

  // mam_coupling::view_1d refrtabsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready
  // mam_coupling::view_1d refitabsw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands]; // ready

  // mam_coupling::view_3d absplw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready
  // mam_coupling::view_1d refrtablw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready
  // mam_coupling::view_1d refitablw_[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands]; // ready

    // aerosol processes
  //std::unique_ptr<mam4::OpticsProcess> optics_;
  // std::unique_ptr<mam4::CalcSizeProcess> calcsize_process_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

  // work arrays
  // FIXME: try to remove this work arrays:I may need to reorganized data in specrefndxsw_ and specrefndxlw_
  mam_coupling::complex_view_3d specrefindex_; // work array
  mam_coupling::view_2d work_;
}; // MAMOptics

} // namespace scream

#endif // EAMXX_MAM_OPTICS_HPP
