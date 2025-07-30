#include "eamxx_config.h"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_zm_process_interface.hpp"
#include "physics/share/physics_constants.hpp"

#include "zm_eamxx_bridge.hpp"

#include <mpi.h> // Include the MPI header for special print statement diagnostics

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
  add_field<Required>("cldfrac_tot",          scalar3d_mid, nondim,grid_name, pack_size);
  add_field<Required>("pbl_height",           scalar2d    , m,     grid_name);
  add_field<Required>("landfrac",             scalar2d    , nondim,grid_name);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,     grid_name, pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,   grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,            pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,            pack_size);

  // Output variables
  add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2, grid_name, "ACCUMULATED");
  add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2, grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  add_field<Computed>("zm_prec",              scalar2d,     m/s,   grid_name);
  add_field<Computed>("zm_cape",              scalar2d,     m/s,   grid_name);

  add_field<Computed>("zm_T_mid_tend",        scalar3d_mid, K,     grid_name, pack_size);
  add_field<Computed>("zm_qv_tend",           scalar3d_mid, K,     grid_name, pack_size);

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
  const auto& cldfrac  = get_field_in("cldfrac_tot")    .get_view<const Spack**, Host>();
  const auto& pblh     = get_field_in("pbl_height")     .get_view<const Real*>();
  const auto& landfrac = get_field_in("landfrac")       .get_view<const Real*>();

  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();

  // prepare inputs
  zm_input.ncol           = m_ncol;
  zm_input.pcol           = m_pcol;
  zm_input.dtime          = dt;
  zm_input.is_first_step  = (ts_start.get_num_steps()==0);
  zm_input.phis           = phis;
  // zm_input.z_mid          = z_mid;
  // zm_input.z_int          = z_int;
  zm_input.p_mid          = p_mid;
  zm_input.p_int          = p_int;
  zm_input.p_del          = p_del;
  zm_input.T_mid          = T_mid;
  zm_input.qv             = qv;
  zm_input.qc             = qc;
  zm_input.omega          = omega;
  zm_input.cldfrac        = cldfrac;
  zm_input.pblh           = pblh;
  zm_input.landfrac       = landfrac;

  // initialize buffer variables
  zm_output.init(m_pcol,m_nlev);
  
  // prepare outputs
  zm_output.ncol          = m_ncol;
  zm_output.pcol          = m_pcol;

  //----------------------------------------------------------------------------
  // calculate interface and mid-point altitudes
  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  // calculate_z_int contains a team-level parallel_scan, which requires a special policy
  const auto scan_policy = ekat::ExeSpaceUtils<ZMF::KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_ncol, nlevm_packs);
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const ZMF::KT::MemberType& team) {
    const int i = team.league_rank();
    const auto z_mid_i = ekat::subview(zm_input.z_mid, i);
    const auto z_del_i = ekat::subview(zm_input.z_del, i);
    const auto z_int_i = ekat::subview(zm_input.z_int, i);
    const auto p_mid_i = ekat::subview(zm_input.p_mid, i);
    const auto p_del_i = ekat::subview(zm_input.p_del, i);
    const auto T_mid_i = ekat::subview(zm_input.T_mid, i);
    const auto qv_i    = ekat::subview(zm_input.qv,    i);
    auto z_surf = phis(i) / PC::gravit;
    // const Real z_surf = 0.0; // For now, set z_int(i,m_nlev) = z_surf = 0
    // Calculate z_mid
    PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
    team.team_barrier();
    PF::calculate_z_int(team, m_nlev, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, m_nlev, z_int_i, z_mid_i);
  });
  //----------------------------------------------------------------------------

  // Run ZM
  zm_eamxx_bridge_run( m_nlev, zm_input, zm_output );

  //----------------------------------------------------------------------------
  // // print max precip and cape
  // const auto comm = m_grid->get_comm();
  // auto max_prec = zm_output.precip.data();
  // auto max_cape = zm_output.cape  .data();
  // comm.all_reduce(&max_prec,1,MPI_MAX);
  // comm.all_reduce(&max_cape,1,MPI_MAX);
  // int my_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // if (my_rank==0) {
  //   std::cout << "max prec / cape : " << max_prec(0) << " / " << max_cape(0) << std::endl;
  // }

  // int world_size;
  // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // // Each process has a local value
  // double local_value = (double)world_rank + 1.0; // Example: process 0 has 1.0, process 1 has 2.0, etc.

  // double global_sum;
  // // Perform an Allreduce to sum all local_values into global_sum
  // MPI_Allreduce(&local_value, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // // Calculate the global average
  // double global_average = global_sum / world_size;

  //----------------------------------------------------------------------------
  // Update output fields

  const Real cpair = PC::Cpair;

  const auto& zm_prec = get_field_out("zm_prec").get_view<Real*>();
  const auto& zm_cape = get_field_out("zm_cape").get_view<Real*>();

  // const typename KT::RangePolicy policy (0,m_ncol);
  // Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int& i) {
  //   zm_prec(i) = zm_output.precip (i);
  //   zm_cape(i) = zm_output.cape   (i);
  // });

  const auto& zm_T_mid_tend = get_field_out("zm_T_mid_tend").get_view<Spack**, Host>();
  const auto& zm_qv_tend    = get_field_out("zm_qv_tend").get_view<Spack**, Host>();

  // Kokkos::deep_copy(zm_T_mid_tend,  zm_output.tend_s/cpair);
  // Kokkos::deep_copy(zm_qv_tend,     zm_output.tend_q);
  
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_ncol, nlevm_packs);
  Kokkos::parallel_for("zm_diag_outputs",policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const auto i = team.league_rank();
    // 2D output
    zm_prec(i) = zm_output.precip (i);
    zm_cape(i) = zm_output.cape   (i);
    // 3D output (vertically resolved)
    const auto zm_T_mid_tend_i    = ekat::subview( zm_T_mid_tend,     i);
    const auto zm_qv_tend_i       = ekat::subview( zm_qv_tend,        i);
    const auto zm_output_tend_s_i = ekat::subview( zm_output.tend_s,  i);
    const auto zm_output_tend_q_i = ekat::subview( zm_output.tend_q,  i);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevm_packs), [&](const int& k) {
      zm_T_mid_tend_i(k) = zm_output_tend_s_i(k)/cpair;
      zm_qv_tend_i   (k) = zm_output_tend_q_i(k);
    });
  });


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

  zm_buffer_size+= ZMF::zm_input_state::num_2d_midlv_c_views * sizeof(Spack) * m_pcol * nlevm_packs;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_intfc_c_views * sizeof(Spack) * m_pcol * nlevi_packs;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_midlv_f_views * sizeof(Real)  * m_pcol * m_nlev;
  zm_buffer_size+= ZMF::zm_input_state::num_2d_intfc_f_views * sizeof(Real)  * m_pcol * (m_nlev+1);

  zm_buffer_size+= ZMF::zm_output_tend::num_1d_scalr_views * sizeof(Scalar)  * m_pcol;

  zm_buffer_size+= ZMF::zm_output_tend::num_2d_midlv_c_views * sizeof(Spack) * m_pcol * nlevm_packs;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_intfc_c_views * sizeof(Spack) * m_pcol * nlevi_packs;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_midlv_f_views * sizeof(Real)  * m_pcol * m_nlev;
  zm_buffer_size+= ZMF::zm_output_tend::num_2d_intfc_f_views * sizeof(Real)  * m_pcol * (m_nlev+1);

  return zm_buffer_size;
}

/*------------------------------------------------------------------------------------------------*/

void zm_deep_convection::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");

  const int nlevm_packs = ekat::npack<Spack>(m_nlev);
  const int nlevi_packs = ekat::npack<Spack>(m_nlev+1);

  auto num_1d_scalr_views   = ZMF::zm_output_tend::num_1d_scalr_views;
  auto num_2d_midlv_c_views = ZMF::zm_input_state::num_2d_midlv_c_views + ZMF::zm_output_tend::num_2d_midlv_c_views;
  auto num_2d_intfc_c_views = ZMF::zm_input_state::num_2d_intfc_c_views + ZMF::zm_output_tend::num_2d_intfc_c_views;
  auto num_2d_midlv_f_views = ZMF::zm_input_state::num_2d_midlv_f_views + ZMF::zm_output_tend::num_2d_midlv_f_views;
  auto num_2d_intfc_f_views = ZMF::zm_input_state::num_2d_intfc_f_views + ZMF::zm_output_tend::num_2d_intfc_f_views;

  //----------------------------------------------------------------------------
  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // 1D scalar variables
  ZMF::uview_1d<Scalar>* scl_ptrs[num_1d_scalr_views] = { &zm_output.precip,
                                                          &zm_output.cape,
                                                        };
  for (int i=0; i<num_1d_scalr_views; ++i) {
    *scl_ptrs[i] = ZMF::uview_1d<Scalar>(mem, m_pcol);
    mem += scl_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  // 2D "f_" views on mid-point levels
  ZMF::uview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f_views]  = { &zm_input.f_z_mid,
                                                                &zm_input.f_p_mid,
                                                                &zm_input.f_p_del,
                                                                &zm_input.f_T_mid,
                                                                &zm_input.f_qv,
                                                                &zm_input.f_qc,
                                                                &zm_input.f_omega,
                                                                &zm_input.f_cldfrac,
                                                                &zm_output.f_tend_s,
                                                                &zm_output.f_tend_q
                                                              };
  for (int i=0; i<num_2d_midlv_f_views; ++i) {
    *midlv_f_ptrs[i] = ZMF::uview_2dl<Real>(mem, m_pcol, m_nlev);
    mem += midlv_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  // 2D "f_" views on interface levels
  ZMF::uview_2dl<Real>* intfc_f_ptrs[num_2d_intfc_f_views]  = { &zm_input.f_z_int,
                                                                &zm_input.f_p_int,
                                                                &zm_output.f_prec_flux,
                                                                &zm_output.f_mass_flux
                                                              };
  for (int i=0; i<num_2d_intfc_f_views; ++i) {
    *intfc_f_ptrs[i] = ZMF::uview_2dl<Real>(mem, m_pcol, (m_nlev+1));
    mem += intfc_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Spack* s_mem = reinterpret_cast<Spack*>(mem);
  //----------------------------------------------------------------------------
  // 2D views on mid-point levels
  ZMF::uview_2d<Spack>* midlv_c_ptrs[num_2d_midlv_c_views]  = { &zm_input.z_mid,
                                                                &zm_input.z_del,
                                                                &zm_output.tend_s,
                                                                &zm_output.tend_q
                                                              };
  for (int i=0; i<num_2d_midlv_c_views; ++i) {
    *midlv_c_ptrs[i] = ZMF::uview_2d<Spack>(s_mem, m_pcol, nlevm_packs);
    s_mem += midlv_c_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  // 2D variables on interface levels
  ZMF::uview_2d<Spack>* intfc_c_ptrs[num_2d_intfc_c_views]  = { &zm_input.z_int,
                                                                &zm_output.prec_flux,
                                                                &zm_output.mass_flux
                                                              };
  for (int i=0; i<num_2d_intfc_c_views; ++i) {
    *intfc_c_ptrs[i] = ZMF::uview_2d<Spack>(s_mem, m_pcol, nlevi_packs);
    s_mem += intfc_c_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for zm_deep_convection.");
}

/*------------------------------------------------------------------------------------------------*/


} // namespace scream
