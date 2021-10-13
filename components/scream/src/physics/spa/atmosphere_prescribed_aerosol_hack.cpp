#include "atmosphere_prescribed_aerosol.hpp"

namespace scream {

// Had to move this code into it's own compilation unit to avoid
// single-prec nvcc crashes in atmosphere_prescribed_aerosol.cpp
// on weaver.

void SPA::initialize_spa_impl ()
{
  SPAPressureState.ncols         = m_num_cols;
  SPAPressureState.nlevs         = m_num_levs;
  SPAPressureState.hyam          = get_field_in("hyam").get_view<const Pack*>();
  SPAPressureState.hybm          = get_field_in("hybm").get_view<const Pack*>();
//  SPAPressureState.ps_this_month = get_field_out("PS_beg").get_view<Real*>();
//  SPAPressureState.ps_next_month = get_field_out("PS_end").get_view<Real*>();
  SPAPressureState.pmid          = get_field_in("p_mid").get_view<const Pack**>();

//  SPAData_start.CCN3             = get_field_out("CCN3_beg").get_view<Pack**>();
//  SPAData_end.CCN3               = get_field_out("CCN3_end").get_view<Pack**>();
  SPAData_out.CCN3               = get_field_out("nc_activated").get_view<Pack**>();

//  SPAData_start.AER_G_SW         = get_field_out("AER_G_SW_beg").get_view<Pack***>();
//  SPAData_end.AER_G_SW           = get_field_out("AER_G_SW_end").get_view<Pack***>();
  SPAData_out.AER_G_SW           = get_field_out("aero_g_sw").get_view<Pack***>();

//  SPAData_start.AER_SSA_SW       = get_field_out("AER_SSA_SW_beg").get_view<Pack***>();
//  SPAData_end.AER_SSA_SW         = get_field_out("AER_SSA_SW_end").get_view<Pack***>();
  SPAData_out.AER_SSA_SW         = get_field_out("aero_ssa_sw").get_view<Pack***>();

//  SPAData_start.AER_TAU_SW       = get_field_out("AER_TAU_SW_beg").get_view<Pack***>();
//  SPAData_end.AER_TAU_SW         = get_field_out("AER_TAU_SW_end").get_view<Pack***>();
  SPAData_out.AER_TAU_SW         = get_field_out("aero_tau_sw").get_view<Pack***>();

//  SPAData_start.AER_TAU_LW       = get_field_out("AER_TAU_LW_beg").get_view<Pack***>();
//  SPAData_end.AER_TAU_LW         = get_field_out("AER_TAU_LW_end").get_view<Pack***>();
  SPAData_out.AER_TAU_LW         = get_field_out("aero_tau_lw").get_view<Pack***>();
}

} // namespace scream
