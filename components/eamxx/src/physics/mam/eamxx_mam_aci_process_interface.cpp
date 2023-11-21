#include <physics/mam/eamxx_mam_aci_process_interface.hpp>

namespace scream
{
MAMAci::MAMAci(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params){
}


//Return type of the process
AtmosphereProcessType MAMAci::type() const {
  return AtmosphereProcessType::Physics;
}

//return name of the process
std::string MAMAci::name() const{
  return "mam4_aci";
  }

//set grid for all the inputs and outputs
void MAMAci::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI set grid = ");

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} }; // mid points
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {ncol_, nlev_+1} }; //interfaces

  using namespace ekat::units;
  const auto q_unit = kg/kg; // units of mass mixing ratios of tracers
  const auto n_unit = 1/kg; // units of number mixing ratios of tracers

  add_field<Required>("qc",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("qi",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud ice mass mixing ratio [kg/kg]
  add_field<Required>("ni",             scalar3d_layout_mid, n_unit, grid_name, "tracers");// cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,      grid_name); // Temperature[K] at midpoints
  add_field<Required>("omega",          scalar3d_layout_mid, Pa/s,   grid_name); // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,     grid_name); // Total pressure [Pa] at midpoints
  add_field<Required>("p_int",          scalar3d_layout_int, Pa,     grid_name); // Total pressure [Pa] at interfaces
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,     grid_name); // Layer thickness(pdel) [Pa] at midpoints

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit, grid_name, "tracers");
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
      }
    }
    for (int c = 0; c < mam_coupling::num_aero_species(); ++c) {
      const char* cld_mmr_field_name = mam_coupling::cld_aero_mmr_field_name(m, c);
      std::cout<<cld_mmr_field_name<<std::endl;
      if (strlen(cld_mmr_field_name) > 0) {
        //NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are NOT advected
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name);
      }
    }
  }

  //input for aci codes that existed in PBUF in EAM
  //These outputs should come from the cloud macrophysics process (e.g., SHOC)
  const auto m2 = m*m;
  const auto s2 = s*s;
  //FIXME BALLI: w_sec,  is at OLD time step; strat_cld_frac and liq_strat_cld_frac may also need OLD time
  add_field<Required>("w_sec",              scalar3d_layout_mid, m2/s2,  grid_name); // Vertical velocity variance (wp2) at midpoints

  auto nondim = Units::nondimensional();
  add_field<Required>("strat_cld_frac",     scalar3d_layout_mid, nondim, grid_name); // Stratiform cloud fraction at midpoints
  add_field<Required>("liq_strat_cld_frac", scalar3d_layout_mid, nondim, grid_name); // Liquid stratiform cloud fraction  at midpoints
  add_field<Required>("kvh",                scalar3d_layout_mid, m2/s, grid_name); // Eddy diffusivity for heat
  
  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  FieldLayout scalar4d_layout_mid{ {COL, LEV, NUM_MODES}, {ncol_, nlev_, num_modes_} }; // mid points
  add_field<Required>("dgnum", scalar4d_layout_mid, m, grid_name); // dry diameter of aerosols
  

  /*NOTE on other inputs for the aci process:
  1. reciprocal of pseudo_density (rpdel): computed from the pseudo_density
  2. geopotential height at midpoints: computed geopotential height at interfaces, which inturn is computed using
     pseudo_density, p_mid, T_mid and qv_mid (see dry_static_energy.cpp's "compute_diagnostic_impl" function).
     qv_mid can be obtained from "get_field_in" call*/

}

void MAMAci::initialize_impl(const RunType run_type) {
   m_atm_logger->log(ekat::logger::LogLevel::info,"ACI init = ");
   std::cout<<"BALLI==============================="<<std::endl;
}

void MAMAci::run_impl(const double dt) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI run = ");
}

void MAMAci::finalize_impl(){
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI final = ");
}

} // namespace scream
