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
            using field_type       = Field<      Real>;
            using const_field_type = Field<const Real>;

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

            // Register all fields in the given repo
            void register_fields (FieldRepository<Real>& field_repo) const;

            // Get the set of required/computed fields
            const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
            const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

        protected:
            // The three main interfaces for the subcomponent
            void initialize_impl (const util::TimeStamp& t0);
            void run_impl        (const Real dt);
            void finalize_impl   ();

            // Set fields in the atmosphere process
            void set_required_field_impl (const Field<const Real>& f);
            void set_computed_field_impl (const Field<      Real>& f);

            std::set<FieldIdentifier> m_required_fields;
            std::set<FieldIdentifier> m_computed_fields;

            // Input and input/output fields
            std::map<std::string,const_field_type> m_rrtmgp_fields_in;
            std::map<std::string,field_type>       m_rrtmgp_fields_out;

            std::map<std::string,const_field_type::view_type::HostMirror>  m_rrtmgp_host_views_in;
            std::map<std::string,field_type::view_type::HostMirror>        m_rrtmgp_host_views_out;

            std::map<std::string,const Real*>  m_raw_ptrs_in;
            std::map<std::string,Real*>        m_raw_ptrs_out;

            std::shared_ptr<FieldInitializer> m_initializer;

            util::TimeStamp m_current_ts;
            ekat::Comm            m_rrtmgp_comm;
            ekat::ParameterList   m_rrtmgp_params;

            // Keep track of number of columns and levels
            int ncol;
            int nlay;

            // Need to hard-code some dimension sizes for now. 
            // TODO: find a better way of configuring this
            const int nswbands = 14;
            const int nlwbands = 16;

            // These are the gases that we keep track of
            const int ngas = 8;
            const std::string gas_names[8] = {
                "h2o", "co2", "o3", "n2o",
                "co" , "ch4", "o2", "n2"
            };

        private: 
            void require_unpadded(const Field<const Real>& f);

    };  // class RRTMGPRadiation
}  // namespace scream

#endif  // SCREAM_RRTMGP_RADIATION_HPP
