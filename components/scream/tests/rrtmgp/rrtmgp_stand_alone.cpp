#include <iostream>
#include <cmath>
#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/share/physics_only_grids_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "Intrinsics.h"
#include "rrtmgp_test_utils.hpp"

/*
 * Run standalone test problem for RRTMGP and compare with baseline
 */

namespace scream {
       
    /* 
     * Run standalone test through SCREAM driver this time
     */
    TEST_CASE("rrtmgp_scream_stand_alone", "") {
        using namespace scream;
        using namespace scream::control;

        // Get baseline name (needs to be passed as an arg)
        std::string inputfile = ekat::TestSession::get().params.at("rrtmgp_inputfile");
        std::string baseline = ekat::TestSession::get().params.at("rrtmgp_baseline");

        // Check if files exists
        REQUIRE(rrtmgpTest::file_exists(inputfile.c_str()));
        REQUIRE(rrtmgpTest::file_exists(baseline.c_str()));

        // Initialize yakl
        if(!yakl::isInitialized()) { yakl::init(); }

        // Read reference fluxes from baseline file
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        rrtmgpTest::read_fluxes(baseline, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

        // Load ad parameter list
        std::string fname = "input.yaml";
        ekat::ParameterList ad_params("Atmosphere Driver");
        REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );
        // Create a MPI communicator
        ekat::Comm atm_comm (MPI_COMM_WORLD);

        // Need to register products in the factory *before* we create any atm process or grids manager.,
        auto& proc_factory = AtmosphereProcessFactory::instance();
        auto& gm_factory = GridsManagerFactory::instance();
        proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);
        gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);

        // Create the grids manager
        auto& gm_params = ad_params.sublist("Grids Manager");
        const std::string& gm_type = gm_params.get<std::string>("Type");
        auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

        // Create the driver
        AtmosphereDriver ad;

        // Dummy timestamp
        util::TimeStamp time (0,0,0,0);

        // Initialize the driver
        ad.initialize(atm_comm, ad_params, time);

        /*
         * Setup the dummy problem and overwrite default initial conditions
         */

        // Get dimension sizes from the field manager
        auto& field_repo = ad.get_field_repo();
        const auto& grid = ad.get_grids_manager()->get_grid("Physics");
        int ncol = grid->get_num_local_dofs();
        int nlay = grid->get_num_vertical_levels();

        // Get number of shortwave bands and number of gases from RRTMGP
        int nswbands = scream::rrtmgp::k_dist_sw.get_nband();
        int ngas     =   8;  // TODO: get this intelligently

        // Make sure we have the right dimension sizes
        REQUIRE(nlay == sw_flux_up_ref.dimension[1]-1);

        // Grab reshaped views from the field manager and wrap pointers in yakl arrays
        auto& d_pmid= field_repo.get_field("pmid", "Physics").get_view();
        auto& d_tmid= field_repo.get_field("tmid", "Physics").get_view();
        auto& d_pint= field_repo.get_field("pint", "Physics").get_view();
        auto& d_tint= field_repo.get_field("tint", "Physics").get_view();
        auto& d_sfc_alb_dir= field_repo.get_field("sfc_alb_dir", "Physics").get_view();
        auto& d_sfc_alb_dif= field_repo.get_field("sfc_alb_dif", "Physics").get_view();
        auto& d_lwp= field_repo.get_field("lwp", "Physics").get_view();
        auto& d_iwp= field_repo.get_field("iwp", "Physics").get_view();
        auto& d_rel= field_repo.get_field("rel", "Physics").get_view();
        auto& d_rei= field_repo.get_field("rei", "Physics").get_view();
        auto& d_mu0= field_repo.get_field("mu0", "Physics").get_view();
        auto& d_gas_vmr = field_repo.get_field("gas_vmr", "Physics").get_view();
        yakl::Array<double,2,memDevice,yakl::styleFortran> p_lay("p_lay", d_pmid.data(), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> t_lay("t_lay", d_tmid.data(), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> p_lev("p_lev", d_pint.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> t_lev("t_lev", d_tint.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sfc_alb_dir("sfc_alb_dir", d_sfc_alb_dir.data(), ncol, nswbands);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sfc_alb_dif("sfc_alb_dif", d_sfc_alb_dif.data(), ncol, nswbands);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lwp("lwp", d_lwp.data(), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> iwp("iwp", d_iwp.data(), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> rel("rel", d_rel.data(), ncol, nlay);
        yakl::Array<double,2,memDevice,yakl::styleFortran> rei("rei", d_rei.data(), ncol, nlay);
        yakl::Array<double,1,memDevice,yakl::styleFortran> mu0("mu0", d_mu0.data(), ncol);
        yakl::Array<double,3,memDevice,yakl::styleFortran> gas_vmr("gas_vmr", d_gas_vmr.data(), ncol, nlay, ngas);

        // Read in dummy Garand atmosphere; if this were an actual model simulation,
        // these would be passed as inputs to the driver
        // NOTE: set ncol to size of col_flx dimension in the input file. This is so
        // that we can compare to the reference data provided in that file. Note that
        // this will copy the first column of the input data (the first profile) ncol
        // times. We will then fill some fraction of these columns with clouds for
        // the test problem.
        GasConcs gas_concs;
        read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, ncol);

        // Setup dummy problem
        rrtmgpTest::dummy_atmos(
            inputfile, ncol, p_lay, t_lay,
            sfc_alb_dir, sfc_alb_dif, mu0,
            lwp, iwp, rel, rei
        );

        // Copy gases from gas_concs to gas_vmr array
        parallel_for(Bounds<3>(ncol,nlay,ngas), YAKL_LAMBDA(int icol, int ilay, int igas) {
            gas_vmr(icol,ilay,igas) = gas_concs.concs(icol,ilay,igas);
        });
        gas_concs.reset();

        // Run driver
        ad.run(300.0);

        // Check values; need to get fluxes from field manager first
        // The AD should have called RRTMGP to calculate these values in the ad.run() call
        auto& d_sw_flux_up = field_repo.get_field("sw_flux_up", "Physics").get_view();
        auto& d_sw_flux_dn = field_repo.get_field("sw_flux_dn", "Physics").get_view();
        auto& d_sw_flux_dn_dir = field_repo.get_field("sw_flux_dn_dir", "Physics").get_view();
        auto& d_lw_flux_up = field_repo.get_field("lw_flux_up", "Physics").get_view();
        auto& d_lw_flux_dn = field_repo.get_field("lw_flux_dn", "Physics").get_view();
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_up_test("sw_flux_up_test", d_sw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn_test("sw_flux_dn_test", d_sw_flux_dn.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> sw_flux_dn_dir_test("sw_flux_dn_dir_test", d_sw_flux_dn_dir.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_up_test("lw_flux_up_test", d_lw_flux_up.data(), ncol, nlay+1);
        yakl::Array<double,2,memDevice,yakl::styleFortran> lw_flux_dn_test("lw_flux_dn_test", d_lw_flux_dn.data(), ncol, nlay+1);

        // Make sure fluxes from field manager that were calculated in AD call of RRTMGP match reference fluxes from input file
        REQUIRE(rrtmgpTest::all_equals(sw_flux_up_ref    , sw_flux_up_test  ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_ref    , sw_flux_dn_test    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_dir_ref, sw_flux_dn_dir_test));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_up_ref    , lw_flux_up_test    ));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_dn_ref    , lw_flux_dn_test    ));

        // Deallocate YAKL arrays
        sw_flux_up_ref.deallocate();
        sw_flux_dn_ref.deallocate();
        sw_flux_dn_dir_ref.deallocate();
        lw_flux_up_ref.deallocate();
        lw_flux_dn_ref.deallocate();

        // Finalize the driver; needs to come before yakl::finalize because
        // rrtmgp::finalize() frees YAKL arrays
        ad.finalize();

        // Finalize YAKL
        yakl::finalize();

        // If we got this far, we were able to run the code through the AD
        REQUIRE(true);
    }
}
