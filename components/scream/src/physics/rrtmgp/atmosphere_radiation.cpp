#include "share/scream_assert.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"

namespace scream {
    RRTMGPRadiation::RRTMGPRadiation (const Comm& comm, const ParameterList& params) : m_rrtmgp_comm (comm), m_rrtmgp_params (params) {
        /*
         * Anything that can be initialized without grid information can be initialized here.
         * I.e., universal constants, options, etc.
         */
    }  // RRTMGPRadiation::RRTMGPRadiation

    void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

    }  // RRTMGPRadiation::set_grids

    void RRTMGPRadiation::initialize(const util::TimeStamp& t0) {
        // Call RRTMGP initialize routine
        //rrtmgp_initialize()
    }
    void RRTMGPRadiation::run       (const Real dt) {}
    void RRTMGPRadiation::finalize  () {}

    // Register all required and computed field in the field repository
    void RRTMGPRadiation::register_fields(FieldRepository<Real, device_type>& field_repo) const {
        for (auto& fid: m_required_fields) { field_repo.register_field(fid); }
        for (auto& fid: m_computed_fields) { field_repo.register_field(fid); }
    }

}  // namespace scream
