#include "eamxx_rayleigh_friction_process_interface.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <sstream>

namespace scream
{

// =========================================================================================
RayleighFriction::RayleighFriction (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_rayk0     = m_params.get<int>("rayleigh_friction_vertical_level", 2);
  m_raykrange = m_params.get<double>("rayleigh_friction_range", 0.0);
  m_raytau0   = m_params.get<double>("rayleigh_friction_decay_time", 5.0);
}

// =========================================================================================
void RayleighFriction::create_requests()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // Initialize grid from grids manager
  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();

  m_ncols = m_grid->get_num_local_dofs();       // Number of columns on this rank
  m_nlevs = m_grid->get_num_vertical_levels();  // Number of levels per column

  EKAT_REQUIRE_MSG(m_rayk0 >= 1 && m_rayk0 <= m_nlevs,
                   "Error! rayleigh_friction_vertical_level must be in [1, nlevs].\n"
                   "  - rayleigh_friction_vertical_level: " + std::to_string(m_rayk0) + "\n"
                   "  - nlevs: " + std::to_string(m_nlevs) + "\n");

  auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);
  auto vector3d_mid = m_grid->get_3d_vector_layout(LEV,2);

  add_field<Updated>("horiz_winds", vector3d_mid, m/s, grid_name);
  add_field<Updated>("T_mid",       scalar3d_mid, K,   grid_name);
}

// =========================================================================================
void RayleighFriction::initialize_impl (const RunType /* run_type */)
{
  m_otau = view_1d("otau", m_nlevs);
  RFFunctions::rayleigh_friction_init(m_nlevs, m_rayk0, m_raykrange, m_raytau0, m_otau);

  if (m_comm.am_i_root()) {
    auto otau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_otau);
    std::stringstream ss;
    ss << "Rayleigh friction parameters:\n"
       << "  rayleigh_friction_vertical_level: " << m_rayk0     << "\n"
       << "  rayleigh_friction_range:          " << m_raykrange << "\n"
       << "  rayleigh_friction_decay_time:     " << m_raytau0   << " [days]\n"
       << "Rayleigh friction decay rate profile:\n";
    for (int k=0; k<m_nlevs; ++k) {
      ss << "  k = " << k+1 << "   otau = " << otau_h(k) << "\n";
    }
    m_atm_logger->info(ss.str());
  }
}

// =========================================================================================
void RayleighFriction::run_impl (const double dt)
{
  // If decay time is 0, no Rayleigh friction is applied
  if (m_raytau0 == 0) return;

  const auto horiz_winds = get_field_out("horiz_winds").get_view<Real***>();
  const auto T_mid       = get_field_out("T_mid").get_view<Real**>();

  RFFunctions::rayleigh_friction_tend(m_ncols, m_nlevs, dt, m_otau, horiz_winds, T_mid);
}

// =========================================================================================
void RayleighFriction::finalize_impl()
{
  // Do nothing
}

} // namespace scream
