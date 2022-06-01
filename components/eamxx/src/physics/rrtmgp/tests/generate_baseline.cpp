#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/mo_garand_atmos_io.h"
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"
#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

#include "cpp/rrtmgp/mo_gas_concentrations.h"

#include "YAKL.h"

#include <ekat/logging/ekat_logger.hpp>

#include <iostream>
#include <cmath>

using namespace scream;

// Names of input files we will need.
// NOTE: use full spectral resolution absorption data for consistency with reference problem
std::string coefficients_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-data-sw-g224-2018-12-04.nc";
std::string coefficients_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-data-lw-g256-2018-12-04.nc";
std::string cloud_optics_file_sw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-sw.nc";
std::string cloud_optics_file_lw = SCREAM_DATA_DIR "/init/rrtmgp-cloud-optics-coeffs-lw.nc";

int main (int argc, char** argv) {
    MPI_Init(&argc,&argv);

    using namespace ekat::logger;
    using logger_t = Logger<LogNoFile,LogRootRank>;

    ekat::Comm comm(MPI_COMM_WORLD);
    auto logger = std::make_shared<logger_t>("",LogLevel::info,comm);

    // Get filenames from command line
    if (argc != 3) {
      std::string msg = "Missing required inputs. Usage:\n";
      msg += argv[0];
      msg += " inputfile baseline\n";
      logger->error(msg);
      return 1;
    }
    std::string inputfile(argv[argc-2]);
    std::string baseline(argv[argc-1]);

    // Initialize yakl
    yakl::init();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    real2d sw_flux_up_ref;
    real2d sw_flux_dn_ref;
    real2d sw_flux_dn_dir_ref;
    real2d lw_flux_up_ref;
    real2d lw_flux_dn_ref;
    logger->info("read_fluxes...");
    rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

    // Get dimension sizes
    int ncol = sw_flux_up_ref.dimension[0];
    int nlev = sw_flux_up_ref.dimension[1];
    int nlay = nlev - 1;

    // Read in dummy Garand atmosphere; if this were an actual model simulation,
    // these would be passed as inputs to the driver
    // NOTE: set ncol to size of col_flx dimension in the input file. This is so
    // that we can compare to the reference data provided in that file. Note that
    // this will copy the first column of the input data (the first profile) ncol
    // times. We will then fill some fraction of these columns with clouds for
    // the test problem.
    real2d p_lay ("p_lay", ncol, nlay);
    real2d t_lay ("t_lay", ncol, nlay);
    real2d p_lev ("p_lev", ncol, nlay+1);
    real2d t_lev ("t_lev", ncol, nlay+1);
    GasConcs gas_concs;
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, ncol);

    // Initialize the RRTMGP interface; this will read in the k-distribution
    // data that contains information about absorption coefficients for gases
    logger->info("rrtmgp_initialize...");
    rrtmgp::rrtmgp_initialize(gas_concs, coefficients_file_sw, coefficients_file_lw, cloud_optics_file_sw, cloud_optics_file_lw, logger);

    // Setup dummy all-sky problem
    real1d sfc_alb_dir_vis ("sfc_alb_dir_vis", ncol);
    real1d sfc_alb_dir_nir ("sfc_alb_dir_nir", ncol);
    real1d sfc_alb_dif_vis ("sfc_alb_dif_vis", ncol);
    real1d sfc_alb_dif_nir ("sfc_alb_dif_nir", ncol);
    real1d mu0 ("mu0", ncol);
    real2d lwp ("lwp", ncol, nlay);
    real2d iwp ("iwp", ncol, nlay);
    real2d rel ("rel", ncol, nlay);
    real2d rei ("rei", ncol, nlay);
    real2d cld ("cld", ncol, nlay);
    rrtmgpTest::dummy_atmos(
        inputfile, ncol, p_lay, t_lay,
        sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir,
        mu0,
        lwp, iwp, rel, rei, cld
    );

    // Setup flux outputs; In a real model run, the fluxes would be
    // input/outputs into the driver (persisting between calls), and
    // we would just have to setup the pointers to them in the
    // FluxesBroadband object
    const auto nswbands = scream::rrtmgp::k_dist_sw.get_nband();
    const auto nlwbands = scream::rrtmgp::k_dist_lw.get_nband();
    real2d sw_flux_up ("sw_flux_up" , ncol, nlay+1);
    real2d sw_flux_dn ("sw_flux_dn" , ncol, nlay+1);
    real2d sw_flux_dn_dir("sw_flux_dn_dir", ncol, nlay+1);
    real2d lw_flux_up ("lw_flux_up" , ncol, nlay+1);
    real2d lw_flux_dn ("lw_flux_dn" , ncol, nlay+1);
    real2d sw_clrsky_flux_up ("sw_clrsky_flux_up" , ncol, nlay+1);
    real2d sw_clrsky_flux_dn ("sw_clrsky_flux_dn" , ncol, nlay+1);
    real2d sw_clrsky_flux_dn_dir("sw_clrsky_flux_dn_dir", ncol, nlay+1);
    real2d lw_clrsky_flux_up ("lw_clrsky_flux_up" , ncol, nlay+1);
    real2d lw_clrsky_flux_dn ("lw_clrsky_flux_dn" , ncol, nlay+1);
    real3d sw_bnd_flux_up ("sw_bnd_flux_up" , ncol, nlay+1, nswbands);
    real3d sw_bnd_flux_dn ("sw_bnd_flux_dn" , ncol, nlay+1, nswbands);
    real3d sw_bnd_flux_dir("sw_bnd_flux_dir", ncol, nlay+1, nswbands);
    real3d lw_bnd_flux_up ("lw_bnd_flux_up" ,ncol, nlay+1, nlwbands);
    real3d lw_bnd_flux_dn ("lw_bnd_flux_dn" ,ncol, nlay+1, nlwbands);

    // Compute band-by-band surface_albedos.
    real2d sfc_alb_dir("sfc_alb_dir", ncol, nswbands);
    real2d sfc_alb_dif("sfc_alb_dif", ncol, nswbands);
    rrtmgp::compute_band_by_band_surface_albedos(
      ncol, nswbands,
      sfc_alb_dir_vis, sfc_alb_dir_nir,
      sfc_alb_dif_vis, sfc_alb_dif_nir,
      sfc_alb_dir, sfc_alb_dif);

    // Setup some dummy aerosol optical properties
    auto aer_tau_sw = real3d("aer_tau_sw", ncol, nlay, nswbands);
    auto aer_ssa_sw = real3d("aer_ssa_sw", ncol, nlay, nswbands);
    auto aer_asm_sw = real3d("aer_asm_sw", ncol, nlay, nswbands);
    auto aer_tau_lw = real3d("aer_tau_lw", ncol, nlay, nlwbands);
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nswbands,nlay,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        aer_tau_sw(icol,ilay,ibnd) = 0;
        aer_ssa_sw(icol,ilay,ibnd) = 0;
        aer_asm_sw(icol,ilay,ibnd) = 0;
    });
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(nlwbands,nlay,ncol), YAKL_LAMBDA(int ibnd, int ilay, int icol) {
        aer_tau_lw(icol,ilay,ibnd) = 0;
    });

    // These are returned as outputs now from rrtmgp_main
    // TODO: provide as inputs consistent with how aerosol is treated?
    const auto nswgpts = scream::rrtmgp::k_dist_sw.get_ngpt();
    const auto nlwgpts = scream::rrtmgp::k_dist_lw.get_ngpt();
    auto cld_tau_sw = real3d("cld_tau_sw", ncol, nlay, nswgpts);
    auto cld_tau_lw = real3d("cld_tau_lw", ncol, nlay, nlwgpts);

    // Run RRTMGP standalone codes and compare with AD run
    // Do something interesting here...
    // NOTE: these will get replaced with AD stuff that handles these
    logger->info("rrtmgp_main...");
    const Real tsi_scaling = 1;
    rrtmgp::rrtmgp_main(
        ncol, nlay,
        p_lay, t_lay, p_lev, t_lev, gas_concs,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei, cld,
        aer_tau_sw, aer_ssa_sw, aer_asm_sw, aer_tau_lw,
        cld_tau_sw, cld_tau_lw,  // outputs
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn,
        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
        lw_clrsky_flux_up, lw_clrsky_flux_dn,
        sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
        lw_bnd_flux_up, lw_bnd_flux_dn, tsi_scaling,
        logger
    );

    // Write fluxes
    logger->info("write_fluxes...");
    rrtmgpTest::write_fluxes(
        baseline, 
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    logger->info("cleaning up...");
    // Clean up from test; this is probably not necessary, these things
    // should be deallocated when they fall out of scope, but we should be
    // good citizens and clean up our mess.
    p_lay.deallocate();
    t_lay.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();
    sfc_alb_dir_vis.deallocate();
    sfc_alb_dir_nir.deallocate();
    sfc_alb_dif_vis.deallocate();
    sfc_alb_dif_nir.deallocate();
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();
    mu0.deallocate();
    lwp.deallocate();
    iwp.deallocate();
    rel.deallocate();
    rei.deallocate();
    cld.deallocate();
    aer_tau_sw.deallocate();
    aer_ssa_sw.deallocate();
    aer_asm_sw.deallocate();
    aer_tau_lw.deallocate();
    cld_tau_sw.deallocate();
    cld_tau_lw.deallocate();
    sw_flux_up_ref.deallocate();
    sw_flux_dn_ref.deallocate();
    sw_flux_dn_dir_ref.deallocate();
    lw_flux_up_ref.deallocate();
    lw_flux_dn_ref.deallocate();
    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dn_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();
    sw_clrsky_flux_up.deallocate();
    sw_clrsky_flux_dn.deallocate();
    sw_clrsky_flux_dn_dir.deallocate();
    lw_clrsky_flux_up.deallocate();
    lw_clrsky_flux_dn.deallocate();
    sw_bnd_flux_up.deallocate();
    sw_bnd_flux_dn.deallocate();
    sw_bnd_flux_dir.deallocate();
    lw_bnd_flux_up.deallocate();
    lw_bnd_flux_dn.deallocate();

    gas_concs.reset();
    rrtmgp::rrtmgp_finalize();
    yakl::finalize();

    MPI_Finalize();

    return 0;
}
