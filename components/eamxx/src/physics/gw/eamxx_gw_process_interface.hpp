#ifndef SCREAM_GW_DRAG_HPP
#define SCREAM_GW_DRAG_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "physics/gw/gw_functions.hpp"

#include <ekat_parameter_list.hpp>
#include <string>

namespace scream
{

/* Gravity Wave Drag Parameterization Suite

This suite of parameterizations can represent the drag from these sources:
  orographic
  frontogenesis
  deep convection
*/

class GWDrag : public AtmosphereProcess
{

  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using GWF = gw::Functions<Real, DefaultDevice>;
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;
  using WSM = typename GWF::WorkspaceManager;

  using Scalar   = typename GWF::Scalar;
  using Pack     = typename GWF::Pack;
  using IntPack  = typename GWF::IntPack;

  using uview_2d = GWF::uview_2d<Pack>;
  using uview_3d = GWF::uview_3d<Pack>;

  public:
    // Constructors
    GWDrag (const ekat::Comm& comm, const ekat::ParameterList& params);

    // The type of subcomponent
    AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

    // The name of the subcomponent
    std::string name () const { return "gw"; }

    // Set the grid
    void create_requests ();

    // Structure for storing local variables initialized using the ATMBufferManager
    struct Buffer {
      static constexpr int pcnst               = 3; // number of constituents (qv, qc, qi)
      static constexpr int num_3d_int_views    = 1; // tau uses interface levels (pver+1)
      static constexpr int num_3d_pcnst_views  = 3;
      static constexpr int num_3d_cd_int_views = 1; // taucd uses interface levels (pver+1)
      static constexpr int num_3d_pgw_views    = 1;
      static constexpr int num_2d_mid_views    = 14;
      static constexpr int num_2d_int_views    = 8;
      static constexpr int num_2d_pgw_views    = 1;
      uview_2d z_mid;       // mid-point altitude
      uview_2d z_int;       // interface altitude
      uview_2d z_del;       // thickness of layer altitudes
      uview_2d p_del_rcp;   // reciprcal of p_del
      uview_2d p_int_log;   // natural log of p_int
      uview_2d T_int;       // interface absolute temperature (dimension must equal T_mid)
      uview_2d N_mid;       // mid-point Brunt-Vaisalla frequency
      uview_2d N_int;       // interface Brunt-Vaisalla frequency
      uview_2d rho_int;     // interface density
      uview_3d q_combined;  // combined qv, qc, qi [ncol][pver][3]

      uview_3d tau;         // gravity wave Reynolds stress
      uview_2d ubm;         // mid-point projection of wind
      uview_2d ubi;         // interface projection of wind
      uview_2d c;           // calculated gravity wave phase speeds

      uview_2d kvtt;        // molecular diffusivity
      uview_2d dse;         // dry static energy

      uview_2d utgw;        // intermediate zonal wind tendency [m/s/s]
      uview_2d vtgw;        // intermediate meridional wind tendency [m/s/s]
      uview_2d ttgw;        // intermediate temperature tendency [J/s]
      uview_3d qtgw;        // intermediate constituents tendencies [kg/kg/s]

      uview_2d gw_tend_u;   // aggregated output zonal wind tendency [m/s/s]
      uview_2d gw_tend_v;   // aggregated output meridional wind tendency [m/s/s]
      uview_2d gw_tend_t;   // aggregated output dry static energy tendency [J/s]
      uview_3d gw_tend_q;   // aggregated output constituents tendencies [kg/kg/s]

      uview_3d taucd;       // reynolds stress for waves propagating in each cardinal direction

      uview_2d egwdffi;     // effective gw diffusivity at interfaces needed for output
      uview_3d gwut;        // gravity wave wind tendency for each wave
      uview_2d dttdf;       // temperature tendencies from diffusion
      uview_2d dttke;       // temperature tendencies from kinetic energy

    };

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
  protected:
#endif

    void initialize_impl (const RunType run_type);
    void run_impl        (const double dt);

  protected:

    void finalize_impl   ();

    // Computes total number of bytes needed for local variables
    size_t requested_buffer_size_in_bytes() const;

    // Set local variables using memory provided by the ATMBufferManager
    void init_buffers(const ATMBufferManager &buffer_manager);

    // Struct which contains local variables
    Buffer m_buffer;

    std::shared_ptr<const AbstractGrid> m_grid;
    int m_ncol;
    int m_nlev;
    int m_npgw;

    Field m_lat;

}; // class GWDrag

} // namespace scream

#endif // SCREAM_GW_DRAG_HPP
