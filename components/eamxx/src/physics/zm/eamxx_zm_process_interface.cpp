#include "eamxx_config.h"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_zm_process_interface.hpp"
#include "physics/share/physics_constants.hpp"

#include "zm_eamxx_bridge.hpp"

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
// Constructor for the zm_deep_convection interface
zm_deep_convection::zm_deep_convection(const ekat::Comm& comm,
                                       const ekat::ParameterList& params)
 : AtmosphereProcess(comm,params)
{
  // params holds all runtime options - what do we need for ZM?
}

/*------------------------------------------------------------------------------------------------*/
void zm_deep_convection::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  const auto layout = m_grid->get_3d_scalar_layout(true);
  const auto comm = m_grid->get_comm();

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();

  // get max ncol value across ranks to mimic how pcols is used on the fortran side
  m_pcol = m_ncol;
  comm.all_reduce(&m_pcol,1,MPI_MAX);

  constexpr int pack_size = Spack::n;

  const auto nondim = Units::nondimensional();
  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();        // Layout for 2D variable
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);    // Layout for 3D variable at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);   // Layout for 3D variable at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);  // Layout for horiz_wind field

  // Input variables
  add_field<Required>("p_mid",                scalar3d_mid, Pa,    grid_name, pack_size);
  add_field<Required>("p_int",                scalar3d_int, Pa,    grid_name, pack_size);
  add_field<Required>("pseudo_density",       scalar3d_mid, Pa,    grid_name, pack_size);
  add_field<Required>("phis",                 scalar2d    , m2/s2, grid_name);
  add_field<Required>("omega",                scalar3d_mid, Pa/s,  grid_name, pack_size);
  add_field<Required>("pbl_height",           scalar2d    , m,     grid_name);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,     grid_name, pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,   grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,            pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,            pack_size);

  // Output variables
  add_field<Updated>("precip_liq_surf_mass",  scalar2d,     kg/m2, grid_name, "ACCUMULATED");
  add_field<Updated>("precip_ice_surf_mass",  scalar2d,     kg/m2, grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  // ???

}

/*------------------------------------------------------------------------------------------------*/
void zm_deep_convection::initialize_impl (const RunType)
{
  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);

  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_liq_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_ice_surf_mass"),m_grid,0.0,false);

  // initialize variables on the fortran side
  zm::zm_eamxx_bridge_init( m_pcol, m_nlev );
}

/*------------------------------------------------------------------------------------------------*/
void zm_deep_convection::run_impl (const double dt)
{

  auto ts_start = start_of_step_ts();

  // get fields
  const auto& phis     = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid    = get_field_in("p_mid")         .get_view<const Spack**, Host>();
  const auto& p_int    = get_field_in("p_int")         .get_view<const Spack**, Host>();
  const auto& p_del    = get_field_in("pseudo_density").get_view<const Spack**, Host>();
  
  const auto& T_mid    = get_field_out("T_mid")         .get_view<Spack**, Host>();
  const auto& qv       = get_field_out("qv")            .get_view<Spack**, Host>();
  const auto& qc       = get_field_out("qc")            .get_view<Spack**, Host>();
  const auto& omega    = get_field_in("omega")          .get_view<const Spack**, Host>();
  const auto& pblh     = get_field_in("pbl_height")     .get_view<const Real*>();

  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();

  // view_2d<Spack>  z_del;
  // view_2d<Scalar> z_surf;
  // // calculate_z_int contains a team-level parallel_scan, which requires a special policy
  // const auto scan_policy = ekat::ExeSpaceUtils<ZMF::KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol, nlev_packs);
  // Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const TMSFunctions::KT::MemberType& team) {
  //   const int i = team.league_rank();

  //   const auto z_del_i = ekat::subview(z_del, i);
  //   const auto z_int_i = ekat::subview(z_int, i);
  //   const auto z_mid_i = ekat::subview(z_mid, i);

  //   z_surf(i) = phis(i) / PC::gravit;

  //   // Calculate z_mid
  //   PF::calculate_dz(team, pseudo_density_i, p_mid_i, T_mid_i, qv_i, z_del_i);
  //   const Real z_surf = 0.0; // For now, set z_int(i,nlevs) = z_surf = 0
  //   team.team_barrier();
  //   PF::calculate_z_int(team, nlevs, z_del_i, z_surf(i), z_int_i);
  //   team.team_barrier();
  //   PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
  // });


  // prepare inputs
  zm_input.ncol           = m_ncol;
  zm_input.pcol           = m_pcol;
  zm_input.dtime          = dt;
  zm_input.is_first_step  = (ts_start.get_num_steps()==0);
  zm_input.phis           = phis;
  // zm_input.z_mid          = ???;
  // zm_input.z_int          = ???;
  zm_input.p_mid          = p_mid;
  zm_input.p_int          = p_int;
  zm_input.p_del          = p_del;
  zm_input.T_mid          = T_mid;
  zm_input.qv             = qv;
  zm_input.qc             = qc;
  zm_input.omega          = omega;
  zm_input.pblh           = pblh;

  // initial buffer variables
  zm_buff.init(m_pcol,m_nlev);
  
  // prepare outputs
  zm_output.ncol          = m_ncol;
  zm_output.pcol          = m_pcol;
  zm_output.precip        = zm_buff.precip;
  zm_output.tend_s        = zm_buff.tend_s;
  zm_output.tend_q        = zm_buff.tend_q;
  zm_output.prec_flux     = zm_buff.prec_flux;
  zm_output.mass_flux     = zm_buff.mass_flux;

  zm_output.f_tend_s      = zm_buff.f_tend_s;
  zm_output.f_tend_q      = zm_buff.f_tend_q;

  // Run ZM
  zm_eamxx_bridge_run( m_nlev, zm_input, zm_output, zm_buff );

  // Update output fields
  // ???
}

/*------------------------------------------------------------------------------------------------*/
void zm_deep_convection::finalize_impl ()
{
  // placeholder for final cleanup
}

/*------------------------------------------------------------------------------------------------*/

size_t zm_deep_convection::requested_buffer_size_in_bytes() const
{
  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  const int nlevi_packs = ekat::npack<Spack>(m_nlev+1);
  size_t zm_buffer_size = 0;
  zm_buffer_size+= ZMF::zm_buffer_data::num_1d_scalr_views * sizeof(Real)  * m_pcol;
  zm_buffer_size+= ZMF::zm_buffer_data::num_2d_midlv_views * sizeof(Real)  * m_pcol * m_nlev; // "f_" views
  zm_buffer_size+= ZMF::zm_buffer_data::num_2d_midlv_views * sizeof(Spack) * m_pcol * nlevm_packs;
  zm_buffer_size+= ZMF::zm_buffer_data::num_2d_intfc_views * sizeof(Spack) * m_pcol * nlevi_packs;
  return zm_buffer_size;
}

/*------------------------------------------------------------------------------------------------*/

void zm_deep_convection::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  const int nlevi_packs = ekat::npack<Spack>(m_nlev+1);

  auto num_1d_scalr_views = ZMF::zm_buffer_data::num_1d_scalr_views;
  auto num_2d_midlv_views = ZMF::zm_buffer_data::num_2d_midlv_views;
  auto num_2d_intfc_views = ZMF::zm_buffer_data::num_2d_intfc_views;

  //----------------------------------------------------------------------------

  using scalar_1d_view_t = decltype(zm_buff.precip);

  // 1D scalar variables
  scalar_1d_view_t* scl_ptrs[num_1d_scalr_views]  = { &zm_buff.precip
                                                    };
  for (int i=0; i<num_1d_scalr_views; ++i) {
    *scl_ptrs[i] = scalar_1d_view_t(mem, m_pcol);
    mem += scl_ptrs[i]->size();
  }

  // 2D "f_" views on mid-point levels
  uview_2dl* midlv_f_ptrs[num_2d_midlv_views] = { &zm_buff.f_tend_s,
                                                  &zm_buff.f_tend_q
                                                };
  for (int i=0; i<num_2d_midlv_views; ++i) {
    *midlv_f_ptrs[i] = uview_2dl(mem, m_pcol, m_nlev);
    mem += midlv_f_ptrs[i]->size();
  }

  //----------------------------------------------------------------------------

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2D views on mid-point levels
  uview_2d* midlv_ptrs[num_2d_midlv_views]  = { &zm_buff.tend_s,
                                                &zm_buff.tend_q
                                              };
  for (int i=0; i<num_2d_midlv_views; ++i) {
    *midlv_ptrs[i] = uview_2d(s_mem, m_pcol, nlevm_packs);
    s_mem += midlv_ptrs[i]->size();
  }


  // 2D variables on interface levels
  uview_2d* intfc_ptrs[num_2d_intfc_views]  = { &zm_buff.prec_flux,
                                                &zm_buff.mass_flux
                                              };
  for (int i=0; i<num_2d_intfc_views; ++i) {
    *intfc_ptrs[i] = uview_2d(s_mem, m_pcol, nlevi_packs);
    s_mem += intfc_ptrs[i]->size();
  }

  //----------------------------------------------------------------------------

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for zm_deep_convection.");
}

/*------------------------------------------------------------------------------------------------*/


} // namespace scream
