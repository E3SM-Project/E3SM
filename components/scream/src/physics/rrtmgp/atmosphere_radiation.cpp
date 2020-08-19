#include "ekat/ekat_assert.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/rrtmgp_inputs_initializer.hpp"
#include "mo_gas_concentrations.h"

namespace scream {
    RRTMGPRadiation::RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params) 
        : AtmosphereProcess::AtmosphereProcess(), m_rrtmgp_comm (comm), m_rrtmgp_params (params) {
        /*
         * Anything that can be initialized without grid information can be initialized here.
         * I.e., universal constants, options, etc.
         */
          if (!m_rrtmgp_params.isParameter("Grid")) {
              m_rrtmgp_params.set("Grid",std::string("SE Physics"));
          }

          m_initializer = create_field_initializer<RRTMGPInputsInitializer>();
    }  // RRTMGPRadiation::RRTMGPRadiation

    void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

        using namespace ekat::units;

        auto kgkg = kg/kg;
        kgkg.set_string("kg/kg");
        auto m3 = m * m * m;
        m3.set_string("m3");
        auto Wm2 = W / m / m;
        Wm2.set_string("Wm2");
        auto nondim = m/m;  // dummy unit for non-dimensional fields
        auto micron = m / 1000000;

        auto VL  = FieldTag::VerticalLevel;
        auto COL = FieldTag::Column;
        auto VAR = FieldTag::Variable;
        auto CMP = FieldTag::Component;
        constexpr int NVL = 42;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

        auto grid = grids_manager->get_grid("Physics");
        const int ncol = grid->get_num_local_dofs();

        // Need to hard-code some dimension sizes for now. TODO: find a better way of configuring this
        int nswbands = 14;
        int nlwbands = 16;
        int ngas = 8;

        // Set up dimension layouts
        FieldLayout scalar2d_layout     { {COL   }, {ncol    } };
        FieldLayout scalar3d_layout_mid { {COL,VL}, {ncol,NVL} };
        FieldLayout scalar3d_layout_int { {COL,VL}, {ncol,NVL+1} };
        // Use VAR field tag for gases for now; consider adding a tag?
        FieldLayout gas_layout          { {VAR,COL,VL}, {ngas,ncol,NVL} };
        FieldLayout scalar2d_swband_layout { {CMP,COL}, {nswbands,ncol} };

        // Set required (input) fields here
        m_required_fields.emplace("pmid" , scalar3d_layout_mid, Pa, grid->name());
        m_required_fields.emplace("pint", scalar3d_layout_int, Pa, grid->name());
        m_required_fields.emplace("tmid" , scalar3d_layout_mid, K , grid->name());
        m_required_fields.emplace("tint" , scalar3d_layout_int, K , grid->name());
        m_required_fields.emplace("col_dry", scalar3d_layout_mid, kgkg, grid->name());
        //m_required_fields.emplace("gas_names", gas_names_layout, nondim, grid->name());
        m_required_fields.emplace("gas_vmr", gas_layout, kgkg, grid->name());
        m_required_fields.emplace("sfc_alb_dir", scalar2d_swband_layout, nondim, grid->name());
        m_required_fields.emplace("sfc_alb_dif", scalar2d_swband_layout, nondim, grid->name());
        m_required_fields.emplace("mu0", scalar2d_layout, nondim, grid->name());
        m_required_fields.emplace("lwp", scalar3d_layout_mid, kg/m3, grid->name());
        m_required_fields.emplace("iwp", scalar3d_layout_mid, kg/m3, grid->name());
        m_required_fields.emplace("rel", scalar3d_layout_mid, micron, grid->name());
        m_required_fields.emplace("rei", scalar3d_layout_mid, micron, grid->name());

        // Outputs; these needed to be added to both required and computed fields?
        m_required_fields.emplace("sw_flux_dn", scalar3d_layout_int, Wm2, grid->name());
        m_required_fields.emplace("sw_flux_up", scalar3d_layout_int, Wm2, grid->name());
        m_required_fields.emplace("sw_flux_dn_dir", scalar3d_layout_int, Wm2, grid->name());
        m_required_fields.emplace("lw_flux_up", scalar3d_layout_int, Wm2, grid->name());
        m_required_fields.emplace("lw_flux_dn", scalar3d_layout_int, Wm2, grid->name());

        // Set computed (output) fields
        m_computed_fields.emplace("sw_flux_dn", scalar3d_layout_int, Wm2, grid->name());
        m_computed_fields.emplace("sw_flux_up", scalar3d_layout_int, Wm2, grid->name());
        m_computed_fields.emplace("sw_flux_dn_dir", scalar3d_layout_int, Wm2, grid->name());
        m_computed_fields.emplace("lw_flux_up", scalar3d_layout_int, Wm2, grid->name());
        m_computed_fields.emplace("lw_flux_dn", scalar3d_layout_int, Wm2, grid->name());

    }  // RRTMGPRadiation::set_grids

    void RRTMGPRadiation::initialize_impl(const util::TimeStamp& t0) {
        // We may have to init some fields from within RRTMGP. This can be the case in a RRTMGP standalone run.
        // Some options:
        //  - we can tell RRTMGP it can init all inputs or specify which ones it can init. We call the
        //    resulting list of inputs the 'initializable' (or initable) inputs. The default is
        //    that no inputs can be inited.
        //  - we can request that RRTMGP either inits no inputs or all of the initable ones (as specified
        //    at the previous point). The default is that RRTMGP must be in charge of init ing ALL or NONE
        //    of its initable inputs.
        // Recall that:
        //  - initable fields may not need initialization (e.g., some other atm proc that
        //    appears earlier in the atm dag might provide them).
        std::vector<std::string> rrtmgp_inputs = {
            "pmid","pint","tmid","tint","col_dry","gas_vmr","sfc_alb_dir","sfc_alb_dif",
            "mu0","lwp","iwp","rel","rei",
            "sw_flux_up", "sw_flux_dn", "sw_flux_dn_dir", "lw_flux_up", "lw_flux_dn"
        };
        using strvec = std::vector<std::string>;
        const strvec& allowed_to_init = m_rrtmgp_params.get<strvec>("Initializable Inputs",strvec(0));
        const bool can_init_all = m_rrtmgp_params.get<bool>("Can Initialize All Inputs", false);
        const bool init_all_or_none = m_rrtmgp_params.get<bool>("Must Init All Inputs Or None", true);
        const strvec& initable = can_init_all ? rrtmgp_inputs : allowed_to_init;
        //const strvec& initable = rrtmgp_inputs;
        if (initable.size()>0) {
            bool all_inited = true, all_uninited = true;
            for (const auto& name : initable) {
                const auto& f = m_rrtmgp_fields_in.at(name);
                const auto& track = f.get_header().get_tracking();
                if (track.get_init_type()==InitType::None) {
                    // Nobody claimed to init this field. RRTMGPInputsInitializer will take care of it
                    m_initializer->add_me_as_initializer(f);
                    all_uninited &= true;
                    all_inited &= false;
                } else {
                    all_uninited &= false;
                    all_inited &= true;
                }
            }
    
            // In order to gurantee some consistency between inputs, it is best if RRTMGP
            // initializes either none or all of the inputs.
            EKAT_REQUIRE_MSG (
                !init_all_or_none || all_inited || all_uninited,
                "Error! Some rrtmgp inputs were marked to be inited by RRTMGP, while others weren't.\n"
                "       RRTMGP was requested to init either all or none of the inputs.\n"
            );
        }
        rrtmgp::rrtmgp_initialize();
    }

    void RRTMGPRadiation::run_impl      (const Real dt) {
        // Get data from AD; RRTMGP wants YAKL views
        // TODO: how can I just keep these around without having to create every time?
        // They are just pointers, so should be able to keep them somewhere else and just associate them once?

        // Get device views
        auto d_pmid = m_rrtmgp_fields_in.at("pmid").get_view();
        auto d_pint = m_rrtmgp_fields_in.at("pint").get_view();
        auto d_tmid = m_rrtmgp_fields_in.at("tmid").get_view();
        auto d_tint = m_rrtmgp_fields_in.at("tint").get_view();
        auto d_col_dry = m_rrtmgp_fields_in.at("col_dry").get_view();
        auto d_gas_vmr = m_rrtmgp_fields_in.at("gas_vmr").get_view();
        auto d_sfc_alb_dir = m_rrtmgp_fields_in.at("sfc_alb_dir").get_view();
        auto d_sfc_alb_dif = m_rrtmgp_fields_in.at("sfc_alb_dif").get_view();
        auto d_mu0 = m_rrtmgp_fields_in.at("mu0").get_view();
        auto d_lwp = m_rrtmgp_fields_in.at("lwp").get_view();
        auto d_iwp = m_rrtmgp_fields_in.at("iwp").get_view();
        auto d_rel = m_rrtmgp_fields_in.at("rel").get_view();
        auto d_rei = m_rrtmgp_fields_in.at("rei").get_view();
        auto d_sw_flux_up = m_rrtmgp_fields_out.at("sw_flux_up").get_view();
        auto d_sw_flux_dn = m_rrtmgp_fields_out.at("sw_flux_dn").get_view();
        auto d_sw_flux_dn_dir = m_rrtmgp_fields_out.at("sw_flux_dn_dir").get_view();
        auto d_lw_flux_up = m_rrtmgp_fields_out.at("lw_flux_up").get_view();
        auto d_lw_flux_dn = m_rrtmgp_fields_out.at("lw_flux_dn").get_view();
 
        // Map to YAKL
        int ngas =  8;
        int ncol = 128;
        int nlay = 42;
        int nswbands = 14;
        yakl::Array<double,2,memDevice,yakl::styleFortran> p_lay  ("p_lay", const_cast<Real*>(d_pmid.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> t_lay  ("t_lay", const_cast<Real*>(d_tmid.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> p_lev  ("p_lev",const_cast<Real*>(d_pint.data()), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> t_lev  ("t_lev",const_cast<Real*>(d_tint.data()), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> col_dry("col_dry",const_cast<Real*>(d_col_dry.data()), ncol, nlay);
        yakl::Array<double,3,memDevice,yakl::styleFortran> gas_vmr("gas_vmr",const_cast<Real*>(d_gas_vmr.data()), ngas, ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sfc_alb_dir("sfc_alb_dir",const_cast<Real*>(d_sfc_alb_dir.data()), nswbands, ncol);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sfc_alb_dif("sfc_alb_dif",const_cast<Real*>(d_sfc_alb_dif.data()), nswbands, ncol);
        yakl::Array<double,1,memDevice,yakl::styleFortran> mu0("mu0",const_cast<Real*>(d_mu0.data()), ncol);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lwp("lwp",const_cast<Real*>(d_lwp.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> iwp("iwp",const_cast<Real*>(d_iwp.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> rel("rel",const_cast<Real*>(d_rel.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> rei("rei",const_cast<Real*>(d_rei.data()), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_up("sw_flux_up", d_sw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn("sw_flux_dn", d_sw_flux_dn.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn_dir("sw_flux_dn_dir", d_sw_flux_dn_dir.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_up("lw_flux_up", d_lw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_dn("lw_flux_dn", d_lw_flux_dn.data(), ncol, nlay+1);

        // Make GasConcs from gas_vmr and gas_names
        string1d gas_names("gas_names",ngas);
        gas_names(1) = std::string("h2o");
        gas_names(2) = std::string("co2");
        gas_names(3) = std::string("o3" );
        gas_names(4) = std::string("n2o");
        gas_names(5) = std::string("co" );
        gas_names(6) = std::string("ch4");
        gas_names(7) = std::string("o2" );
        gas_names(8) = std::string("n2" );

        // Initialize GasConcs object with an "ncol" given from the calling program
        GasConcs gas_concs;
        gas_concs.init(gas_names,ncol,nlay);
        real2d tmp2d;
        tmp2d = real2d("tmp", ncol, nlay);
        for (int igas = 1; igas <= ngas; igas++) {
            for (int icol = 1; icol <= ncol; icol++) {
                for (int ilay = 1; ilay <= nlay; ilay++) {
                    tmp2d(icol,ilay) = gas_vmr(igas,icol,ilay);
                }
            }
            gas_concs.set_vmr(gas_names(igas), tmp2d);
        } 

        // Run RRTMGP driver
        rrtmgp::rrtmgp_main( 
          p_lay, t_lay, p_lev, t_lev,
          gas_concs, col_dry,
          sfc_alb_dir, sfc_alb_dif, mu0,
          lwp, iwp, rel, rei,
          sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
          lw_flux_up, lw_flux_dn
        );
    }

    void RRTMGPRadiation::finalize_impl  () {
        rrtmgp::rrtmgp_finalize();
    }

    // Register all required and computed field in the field repository
    void RRTMGPRadiation::register_fields(FieldRepository<Real>& field_repo) const {
        for (auto& fid: m_required_fields) { field_repo.register_field(fid); }
        for (auto& fid: m_computed_fields) { field_repo.register_field(fid); }
    }

    // Private function to check that fields are not padded
    void RRTMGPRadiation::require_unpadded(const Field<const Real>& f) {
        const auto& fid = f.get_header().get_identifier();
        const auto& layout = fid.get_layout();
        auto v = f.get_view();
        EKAT_REQUIRE_MSG (
            layout.size() == v.size(), 
            "ERROR: field " << fid.name() << " was padded (to allow packing), but currently RRTMGP does not work with padded views."
        );
    }

    void RRTMGPRadiation::set_required_field_impl(const Field<const Real>& f) {
        const auto& name = f.get_header().get_identifier().name();
        m_rrtmgp_fields_in.emplace(name,f);
        m_rrtmgp_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
        m_raw_ptrs_in[name] = m_rrtmgp_host_views_in[name].data();

        // Add myself as customer to the field
        add_me_as_customer(f);

        // Check to make sure that fields are not padded because we are not 
        // sure how to handle that with RRTMGP yet
        require_unpadded(f);
    }

    void RRTMGPRadiation::set_computed_field_impl(const Field<      Real>& f) {
        const auto& name = f.get_header().get_identifier().name();
        m_rrtmgp_fields_out.emplace(name,f);
        m_rrtmgp_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
        m_raw_ptrs_out[name] = m_rrtmgp_host_views_out[name].data();

        // Add myself as provider for the field
        add_me_as_provider(f);

        // Check to make sure that fields are not padded because we are not 
        // sure how to handle that with RRTMGP yet
        require_unpadded(f);
    }

}  // namespace scream
