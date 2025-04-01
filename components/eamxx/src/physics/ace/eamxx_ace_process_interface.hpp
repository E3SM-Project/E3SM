#ifndef EAMXX_ACE_PROCESS_INTERFACE_HPP
#define EAMXX_ACE_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

namespace scream {

struct TorchData;

class ACE : public AtmosphereProcess {
 public:
  ACE(const ekat::Comm &comm, const ekat::ParameterList &params);

  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  std::string name() const override { return "ace"; }

  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

#ifndef KOKKOS_ENABLE_CUDA
 protected:
#endif
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // need three grids (lat-lon, point-grid, surface-coupling)
  std::shared_ptr<const AbstractGrid> m_ll_grid, m_pt_grid, m_sc_grid;

  // a pointer to TorchData struct
  std::shared_ptr<TorchData> m_torch_data;

  // tensor geometries (hardcoded for now)
  int m_batch  = 1;
  int m_in_ch  = 39;
  int m_out_ch = 44;
  int m_height = 180;
  int m_width  = 360;
  // levels
  int m_levels = 8;

  // checkpoint path from the namelist
  std::string m_checkpoint_path;

  // a bool to print some weak verification
  bool m_ace_print_verification = true;

  // mapfiles for remapping
  // sc2ace_map: surface coupling to ACE
  // ace2sc_map: ACE to surface coupling
  std::string m_sc2ace_map, m_ace2sc_map;
  // Need two horiz remappers corresponding to above
  std::shared_ptr<scream::AbstractRemapper> m_sc2ace_remapper,
      m_ace2sc_remapper;

  // helper functions to organize data
  void preprocess();
  void postprocess();

  // store fields to seamlessly move between ll and pt grids
  Field m_land_fraction_ll, m_ocean_fraction_ll, m_sea_ice_fraction_ll;
  Field m_surf_radiative_T_ll, m_T_mid_ll, m_qv_ll, m_wind_ll;
  Field m_p_mid_ll, m_p_int_ll, m_dp_mid_ll;
  Field m_phis_ll;
  Field m_sfc_flux_dir_nir, m_sfc_flux_dir_vis, m_sfc_flux_dif_nir,
      m_sfc_flux_dif_vis, m_sfc_flux_sw_net, m_sfc_flux_lw_dn;
  Field m_precip_liq_surf_mass, m_precip_ice_surf_mass;

  // store two fields for big in/out fields
  Field m_input_field, m_output_field;

  // a map to store additional fields away from FM
  std::map<std::string, Field> m_helper_fields;

  // create helper field, not to be shared with the FM
  Field create_helper_field(const std::string &name, const FieldLayout &layout,
                            const std::string &grid_name);

  // Retrieve a helper field
  // HACK: return different names
  Field get_helper_field(const std::string &name) const {
    return m_helper_fields.at(name + "_helper");
  }
};

}  // namespace scream

#endif  // EAMXX_ACE_PROCESS_INTERFACE_HPP
