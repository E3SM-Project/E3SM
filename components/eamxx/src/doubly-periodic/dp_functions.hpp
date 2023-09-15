#ifndef DP_FUNCTIONS_HPP
#define DP_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"
#include "dp_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

#include "Elements.hpp"
#include "Tracers.hpp"

namespace scream {
namespace dp {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for DP. We use the ETI pattern for
 * these functions.
 *
 * DP assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

using element_t = Homme::Elements;
using tracer_t  = Homme::Tracers;
struct hvcoord_t{};
struct timelevel_t{};
struct hybrid_t{};

template <typename ScalarT, typename DeviceT>
struct Functions
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask  = ekat::Mask<Pack::n>;
  using Smask = ekat::Mask<Spack::n>;

  using KT       = ekat::KokkosTypes<Device>;
  using ExeSpace = typename KT::ExeSpace;

  using C  = physics::Constants<Scalar>;
  using DPC = dp::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using WorkspaceMgr = typename ekat::WorkspaceManager<Spack, Device>;
  using Workspace    = typename WorkspaceMgr::Workspace;

  //
  // --------- Functions ---------
  //

  // ---------------------------------------------------------------------
  // Define the pressures of the interfaces and midpoints from the
  // coordinate definitions and the surface pressure.
  // ---------------------------------------------------------------------
  KOKKOS_FUNCTION
  static void plevs0(
    // Input arguments
    const Int& nver,                   // vertical dimension
    const Scalar& ps,                  // Surface pressure (pascals)
    const uview_1d<const Spack>& hyai, // ps0 component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hyam, // ps0 component of hybrid coordinate - midpoints
    const uview_1d<const Spack>& hybi, // ps component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hybm, // ps component of hybrid coordinate - midpoints
    // Kokkos stuff
    const MemberType& team,
    // Output arguments
    const uview_1d<Spack>& pint,     // Pressure at model interfaces
    const uview_1d<Spack>& pmid,     // Pressure at model levels
    const uview_1d<Spack>& pdel);    // Layer thickness (pint(k+1) - pint(k))

  //-----------------------------------------------------------------------
  // advance_iop_forcing
  // Purpose:
  // Apply large scale forcing for t, q, u, and v as provided by the
  //   case IOP forcing file.
  //
  // Author:
  // Original version: Adopted from CAM3.5/CAM5
  // Updated version for E3SM: Peter Bogenschutz (bogenschutz1@llnl.gov)
  //  and replaces the forecast.F90 routine in CAM3.5/CAM5/CAM6/E3SMv1/E3SMv2
  // CXX version: James Foucar (jgfouca@sandia.gov)
  //
  //-----------------------------------------------------------------------
  KOKKOS_FUNCTION
  static void advance_iop_forcing(
    // Input arguments
    const Int& plev,                         // number of vertical levels
    const Int& pcnst,                        // number of advected constituents including cloud water
    const bool& have_u,                      // dataset contains u
    const bool& have_v,                      // dataset contains v
    const bool& dp_crm,                      // use 3d forcing
    const bool& use_3dfrc,                   // use 3d forcing
    const Scalar& scm_dt,                    // model time step [s]
    const Scalar& ps_in,                     // surface pressure [Pa]
    const uview_1d<const Spack>& u_in,       // zonal wind [m/s]
    const uview_1d<const Spack>& v_in,       // meridional wind [m/s]
    const uview_1d<const Spack>& t_in,       // temperature [K]
    const uview_2d<const Spack>& q_in,       // q tracer array [units vary]
    const uview_1d<const Spack>& t_phys_frc, // temperature forcing from physics [K/s]
    const uview_1d<const Spack>& divt3d,     // 3D T advection
    const uview_2d<const Spack>& divq3d,     // 3D q advection
    const uview_1d<const Spack>& divt,       // Divergence of temperature
    const uview_2d<const Spack>& divq,       // Divergence of moisture
    const uview_1d<const Spack>& wfld,       // Vertical motion (slt)
    const uview_1d<const Spack>& uobs,       // actual u wind
    const uview_1d<const Spack>& vobs,       // actual v wind
    const uview_1d<const Spack>& hyai,       // ps0 component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hyam,       // ps0 component of hybrid coordinate - midpoints
    const uview_1d<const Spack>& hybi,       // ps component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hybm,       // ps component of hybrid coordinate - midpoints
    // Kokkos stuff
    const MemberType& team,
    const Workspace& workspace,
    // Output arguments
    const uview_1d<Spack>& u_update,         // updated temperature [K]
    const uview_1d<Spack>& v_update,         // updated q tracer array [units vary]
    const uview_1d<Spack>& t_update,         // updated zonal wind [m/s]
    const uview_2d<Spack>& q_update);        // updated meridional wind [m/s]

  //-----------------------------------------------------------------------
  // advance_iop_nudging
  // Purpose:
  // Option to nudge t and q to observations as specified by the IOP file
  //
  // Author:
  // Original version: Adopted from CAM3.5/CAM5
  // Updated version for E3SM: Peter Bogenschutz (bogenschutz1@llnl.gov)
  // CXX version: Conrad Clevenger (tccleve@sandia.gov)
  //
  //-----------------------------------------------------------------------
  KOKKOS_FUNCTION
  static void advance_iop_nudging(
    // Input arguments
    const Int& plev,                   // number of vertical levels
    const Scalar& scm_dt,              // model time step [s]
    const Scalar& ps_in,               // surface pressure [Pa]
    const uview_1d<const Spack>& t_in, // temperature [K]
    const uview_1d<const Spack>& q_in, // water vapor mixing ratio [kg/kg]
    const uview_1d<const Spack>& tobs, // observed temperature [K]
    const uview_1d<const Spack>& qobs, // observed vapor mixing ratio [kg/kg]
    const uview_1d<const Spack>& hyai, // ps0 component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hyam, // ps0 component of hybrid coordinate - midpoints
    const uview_1d<const Spack>& hybi, // ps component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hybm, // ps component of hybrid coordinate - midpoints
    // Kokkos stuff
    const MemberType& team,
    const Workspace& workspace,
    // Output arguments
    const uview_1d<Spack>& t_update,   // updated temperature [K]
    const uview_1d<Spack>& q_update,   // updated water vapor [kg/kg]
    const uview_1d<Spack>& relaxt,     // relaxation of temperature [K/s]
    const uview_1d<Spack>& relaxq);    // relaxation of vapor [kg/kg/s]

  KOKKOS_INLINE_FUNCTION
  static void do_advance_iop_subsidence_update(
    const Int& k,
    const Int& plev,
    const Spack& fac,
    const Spack& swfldint,
    const Spack& swfldint_p1,
    const uview_1d<const Spack>& in,
    const uview_1d<const Scalar>& in_s,
    const uview_1d<Spack>& update);

  //-----------------------------------------------------------------------
  //
  // Purpose:
  // Option to compute effects of large scale subsidence on T, q, u, and v.
  // Code originated from CAM3.5/CAM5 Eulerian subsidence computation for SCM
  // in the old forecast.f90 routine.
  //-----------------------------------------------------------------------
  KOKKOS_FUNCTION
  static void advance_iop_subsidence(
    // Input arguments
    const Int& plev,                   // number of vertical levels
    const Int& pcnst,                  // number of advected constituents including cloud water
    const Scalar& scm_dt,              // model time step [s]
    const Scalar& ps_in,               // surface pressure [Pa]
    const uview_1d<const Spack>& u_in, // zonal wind [m/s]
    const uview_1d<const Spack>& v_in, // meridional wind [m/s]
    const uview_1d<const Spack>& t_in, // temperature [K]
    const uview_2d<const Spack>& q_in, // tracer [vary]
    const uview_1d<const Spack>& hyai, // ps0 component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hyam, // ps0 component of hybrid coordinate - midpoints
    const uview_1d<const Spack>& hybi, // ps component of hybrid coordinate - interfaces
    const uview_1d<const Spack>& hybm, // ps component of hybrid coordinate - midpoints
    const uview_1d<const Spack>& wfld, // Vertical motion (slt)
    // Kokkos stuff
    const MemberType& team,
    const Workspace& workspace,
    // Output arguments
    const uview_1d<Spack>& u_update,   // zonal wind [m/s]
    const uview_1d<Spack>& v_update,   // meridional wind [m/s]
    const uview_1d<Spack>& t_update,   // temperature [m/s]
    const uview_2d<Spack>& q_update);  // tracer [vary]

  //---------------------------------------------------------
  // Purpose: Set initial values from IOP files (where available)
  //   when running SCM or DP-CRM.
  //----------------------------------------------------------
  static void iop_setinitial(
    // Input arguments
    const Int& plev,                        // number of vertical levels
    const Int& pcnst,                       // number of advected constituents including cloud water
    const Int& nelemd,                      // number of elements per MPI task
    const Int& np,                          // NP
    const Int& nstep,                       // the timestep number
    const bool& use_replay,                 // use e3sm generated forcing
    const bool& dynproc,                    // Designation of a dynamics processor - AaronDonahue
    const bool& have_t,                     // dataset contains t
    const bool& have_q,                     // dataset contains q
    const bool& have_ps,                    // dataset contains ps
    const bool& have_u,                     // dataset contains u
    const bool& have_v,                     // dataset contains v
    const bool& have_numliq,                // dataset contains numliq
    const bool& have_cldliq,                // dataset contains cldliq
    const bool& have_numice,                // dataset contains numice
    const bool& have_cldice,                // dataset contains cldice
    const bool& scm_zero_non_iop_tracers,   // Ignore all tracers from initial conditions file, and default all tracers not specified in IOP to minimum value (usually zero)
    const bool& is_first_restart_step,      // is first restart step
    const uview_1d<const Scalar>& qmin,     // minimum permitted constituent concentration (kg/kg) (pcnst)
    const uview_1d<const Spack>& uobs,      // actual u wind
    const uview_1d<const Spack>& vobs,      // actual v wind
    const uview_1d<const Spack>& numliqobs, // actual ???
    const uview_1d<const Spack>& numiceobs, // actual ???
    const uview_1d<const Spack>& cldliqobs, // actual ???
    const uview_1d<const Spack>& cldiceobs, // actual ???
    const Scalar& psobs,                    // ???
    const uview_1d<const Scalar>& dx_short, // short length scale in km (nelemd)
    // Input/Output arguments
    Scalar& dyn_dx_size,                    // for use in doubly periodic CRM mode
    tracer_t& tracers,                      // tracers
    element_t& elem,                        // elements
    const uview_1d<Spack>& tobs,            // actual temperature, dims=(plev)
    const uview_1d<Spack>& qobs);           // actual W.V. Mixing ratio

  KOKKOS_FUNCTION
  static void iop_broadcast();

  KOKKOS_FUNCTION
  static void apply_iop_forcing(const Int& nelemd, const uview_1d<element_t>& elem, hvcoord_t& hvcoord, const hybrid_t& hybrid, const timelevel_t& tl, const Int& n, const bool& t_before_advance, const Int& nets, const Int& nete);

  KOKKOS_FUNCTION
  static void iop_domain_relaxation(const Int& nelemd, const Int& np, const Int& nlev, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const uview_1d<Spack>& dp, const Int& nelemd_todo, const Int& np_todo, const Spack& dt);

  KOKKOS_FUNCTION
  static void crm_resolved_turb(const Int& nelemd, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const Int& nelemd_todo, const Int& np_todo);

  static void iop_default_opts(Spack& scmlat_out, Spack& scmlon_out, std::string& iopfile_out, bool& single_column_out, bool& scm_iop_srf_prop_out, bool& iop_nudge_tq_out, bool& iop_nudge_uv_out, Spack& iop_nudge_tq_low_out, Spack& iop_nudge_tq_high_out, Spack& iop_nudge_tscale_out, bool& scm_observed_aero_out, bool& iop_dosubsidence_out, bool& scm_multcols_out, bool& dp_crm_out, Spack& iop_perturb_high_out, bool& precip_off_out, bool& scm_zero_non_iop_tracers_out);

  static void iop_setopts(const Spack& scmlat_in, const Spack& scmlon_in, const std::string& iopfile_in, const bool& single_column_in, const bool& scm_iop_srf_prop_in, const bool& iop_nudge_tq_in, const bool& iop_nudge_uv_in, const Spack& iop_nudge_tq_low_in, const Spack& iop_nudge_tq_high_in, const Spack& iop_nudge_tscale_in, const bool& scm_observed_aero_in, const bool& iop_dosubsidence_in, const bool& scm_multcols_in, const bool& dp_crm_in, const Spack& iop_perturb_high_in, const bool& precip_off_in, const bool& scm_zero_non_iop_tracers_in);

  KOKKOS_FUNCTION
  static void setiopupdate_init();

  KOKKOS_FUNCTION
  static void setiopupdate();

  KOKKOS_FUNCTION
  static void readiopdata(const Int& plev, const bool& iop_update_phase1, const uview_1d<const Spack>& hyam, const uview_1d<const Spack>& hybm);

  KOKKOS_FUNCTION
  static void iop_intht();
}; // struct Functions

} // namespace dp
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)  \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

# include "impl/dp_advance_iop_forcing_impl.hpp"
# include "impl/dp_advance_iop_nudging_impl.hpp"
# include "impl/dp_advance_iop_subsidence_impl.hpp"
# include "impl/dp_iop_setinitial_impl.hpp"
# include "impl/dp_iop_broadcast_impl.hpp"
# include "impl/dp_apply_iop_forcing_impl.hpp"
# include "impl/dp_iop_domain_relaxation_impl.hpp"
# include "impl/dp_crm_resolved_turb_impl.hpp"
# include "impl/dp_iop_default_opts_impl.hpp"
# include "impl/dp_iop_setopts_impl.hpp"
# include "impl/dp_setiopupdate_init_impl.hpp"
# include "impl/dp_setiopupdate_impl.hpp"
# include "impl/dp_readiopdata_impl.hpp"
# include "impl/dp_iop_intht_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

#endif // DP_FUNCTIONS_HPP
