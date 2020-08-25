#ifndef SCREAM_RRTMGP_RADIATION_HPP
#define SCREAM_RRTMGP_RADIATION_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include <string>

namespace scream {
    /*
     * Class responsible for atmosphere radiative transfer. The AD should store
     * exactly ONE instance of this class in its list of subcomponents.
     */

    class RRTMGPRadiation : public AtmosphereProcess {
        public:
            using field_type       = Field<      Real, device_type>;
            using const_field_type = Field<const Real, device_type>;

            // Constructors
            RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params);

            // The type of the subcomponent
            AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

            // The name of the subcomponent
            std::string name () const { return "Radiation"; }

            // The communicator used by the subcomponent
            const ekat::Comm& get_comm () const { return m_rrtmgp_comm; }

            // Required grid for the subcomponent (??)
            std::set<std::string> get_required_grids () const {
                static std::set<std::string> s;
                s.insert(m_rrtmgp_params.get<std::string>("Grid"));
                return s;
            }

            // Set the grid
            void set_grids (const std::shared_ptr<const GridsManager> grid_manager);

            // The three main interfaces for the subcomponent
            void initialize (const util::TimeStamp& t0);
            void run        (const Real dt);
            void finalize   ();

            // Register all fields in the given repo
            void register_fields (FieldRepository<Real, device_type>& field_repo) const;

            // Get the set of required/computed fields
            const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
            const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

        protected:
            // Set fields in the atmosphere process
            void set_required_field_impl (const Field<const Real, device_type>& f);
            void set_computed_field_impl (const Field<      Real, device_type>& f);

            std::set<FieldIdentifier> m_required_fields;
            std::set<FieldIdentifier> m_computed_fields;

            std::map<std::string,const_field_type> m_rrtmgp_fields_in;
            std::map<std::string,field_type>       m_rrtmgp_fields_out;

            std::shared_ptr<FieldInitializer> m_initializer;

            util::TimeStamp m_current_ts;
            ekat::Comm            m_rrtmgp_comm;
            ekat::ParameterList   m_rrtmgp_params;
    };  // class RRTMGPRadiation
}  // namespace scream

#endif  // SCREAM_RRTMGP_RADIATION_HPP
