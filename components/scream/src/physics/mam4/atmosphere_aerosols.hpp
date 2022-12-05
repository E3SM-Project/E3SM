#ifndef SCREAM_MAM4_AEROSOLS_HPP
#define SCREAM_MAM4_AEROSOLS_HPP

#include <share/atm_process/atmosphere_process.hpp>
//#include <share/util/scream_common_physics_functions.hpp>
//#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected public
#define private public
#endif

namespace scream
{

/*
 * The class responsible for handling MAM4 aerosols
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
*/

class MAM4Aerosols final : public scream::AtmosphereProcess
{
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;

  using ColumnView   = mam4::ColumnView;
  using ThreadTeam   = mam4::ThreadTeam;

public:

  // Constructor
  MAM4Aerosols(const ekat::Comm& comm, const ekat::ParameterList& params);

protected:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const int dt) override;
  void finalize_impl() override;

//  void set_computed_group_impl(const FieldGroup& group) override;

private:

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // WSM for internal local variables
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
};

} // namespace scream

#endif // SCREAM_MAM4_AEROSOLS_HPP
